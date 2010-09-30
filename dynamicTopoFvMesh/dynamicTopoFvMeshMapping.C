/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    dynamicTopoFvMesh

Description
    Functions specific to conservative mapping

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*----------------------------------------------------------------------------*/

#include "dynamicTopoFvMesh.H"

#include "meshOps.H"
#include "IOmanip.H"
#include "triFace.H"
#include "objectMap.H"
#include "StaticHashTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Compute mapping weights for modified entities
void dynamicTopoFvMesh::computeMapping
(
    const scalar matchTol,
    const bool skipMapping,
    const label faceStart,
    const label faceSize,
    const label cellStart,
    const label cellSize
)
{
    // Convex-set algorithm for cells
    cellSetAlgorithm cellAlgorithm
    (
        (*this),
        oldPoints_,
        cells_,
        faces_,
        owner_,
        neighbour_,
        pointsFromPoints_,
        modPoints_
    );

    // Compute cell mapping
    for (label cellI = cellStart; cellI < (cellStart + cellSize); cellI++)
    {
        label precisionAttempts = 0;
        label cIndex = cellsFromCells_[cellI].index();

        if (skipMapping)
        {
            // Set empty mapping parameters
            const labelList& mo = cellParents_[cIndex];

            cellsFromCells_[cellI].masterObjects() = mo;
            cellWeights_[cellI].setSize(mo.size(), (1.0/(mo.size() + VSMALL)));
            cellCentres_[cellI].setSize(mo.size(), vector::zero);
        }
        else
        {
            // Obtain weighting factors for this cell.
            computeWeights
            (
                cIndex,
                cellParents_[cIndex],
                polyMesh::cellCells(),
                matchTol,
                cellAlgorithm,
                precisionAttempts,
                cellsFromCells_[cellI].masterObjects(),
                cellWeights_[cellI],
                cellCentres_[cellI]
            );
        }
    }

    // Convex-set algorithm for faces
    faceSetAlgorithm faceAlgorithm
    (
        (*this),
        oldPoints_,
        cells_,
        faces_,
        owner_,
        neighbour_,
        pointsFromPoints_,
        modPoints_
    );

    // Compute face mapping
    for (label faceI = faceStart; faceI < (faceStart + faceSize); faceI++)
    {
        label precisionAttempts = 0;
        label fIndex = facesFromFaces_[faceI].index();
        label patchIndex = whichPatch(fIndex);

        // Skip mapping for internal faces.
        if (patchIndex == -1)
        {
            // Set dummy masters, so that the conventional
            // faceMapper doesn't incur a seg-fault.
            facesFromFaces_[faceI].masterObjects() = labelList(1, 0);
            continue;
        }

        if (skipMapping)
        {
            // Set empty mapping parameters
            const labelList& mo = faceParents_[fIndex];

            facesFromFaces_[faceI].masterObjects() = mo;
            faceWeights_[faceI].setSize(mo.size(), (1.0/(mo.size() + VSMALL)));
            faceCentres_[faceI].setSize(mo.size(), vector::zero);
        }
        else
        {
            // Obtain weighting factors for this face.
            computeWeights
            (
                fIndex,
                faceParents_[fIndex],
                boundaryMesh()[patchIndex].faceFaces(),
                matchTol,
                faceAlgorithm,
                precisionAttempts,
                facesFromFaces_[faceI].masterObjects(),
                faceWeights_[faceI],
                faceCentres_[faceI]
            );
        }
    }
}


// Static equivalent for multiThreading
void dynamicTopoFvMesh::computeMappingThread(void *argument)
{
    // Recast the argument
    meshHandler *thread = static_cast<meshHandler*>(argument);

    if (thread->slave())
    {
        thread->sendSignal(meshHandler::START);
    }

    dynamicTopoFvMesh& mesh = thread->reference();

    // Recast the pointers for the argument
    scalar& matchTol  = *(static_cast<scalar*>(thread->operator()(0)));
    bool& skipMapping = *(static_cast<bool*>(thread->operator()(1)));
    label& faceStart  = *(static_cast<label*>(thread->operator()(2)));
    label& faceSize   = *(static_cast<label*>(thread->operator()(3)));
    label& cellStart  = *(static_cast<label*>(thread->operator()(4)));
    label& cellSize   = *(static_cast<label*>(thread->operator()(5)));

    // Now calculate addressing
    mesh.computeMapping
    (
        matchTol,
        skipMapping,
        faceStart, faceSize,
        cellStart, cellSize
    );

    if (thread->slave())
    {
        thread->sendSignal(meshHandler::STOP);
    }
}


// Routine to invoke threaded mapping
void dynamicTopoFvMesh::threadedMapping
(
    scalar matchTol,
    bool skipMapping
)
{
    label nThreads = threader_->getNumThreads();

    // If mapping is being skipped, issue a warning.
    if (skipMapping)
    {
        Info << " *** Mapping is being skipped *** " << endl;
    }

    // Check if single-threaded
    if (nThreads == 1)
    {
        computeMapping
        (
            matchTol,
            skipMapping,
            0, facesFromFaces_.size(),
            0, cellsFromCells_.size()
        );

        return;
    }

    // Set one handler per thread
    PtrList<meshHandler> hdl(nThreads);

    forAll(hdl, i)
    {
        hdl.set(i, new meshHandler(*this, threader()));
    }

    // Simple load-balancing scheme
    FixedList<label, 2> index(-1);
    FixedList<labelList, 2> tStarts(labelList(nThreads, 0));
    FixedList<labelList, 2> tSizes(labelList(nThreads, 0));

    index[0] = facesFromFaces_.size();
    index[1] = cellsFromCells_.size();

    if (debug > 2)
    {
        Info << " Mapping Faces: " << index[0] << endl;
        Info << " Mapping Cells: " << index[1] << endl;
    }

    forAll(index, indexI)
    {
        label j = 0, total = 0;

        while (index[indexI]--)
        {
            tSizes[indexI][(j = tSizes[indexI].fcIndex(j))]++;
        }

        for (label i = 1; i < tStarts[indexI].size(); i++)
        {
            tStarts[indexI][i] = tSizes[indexI][i-1] + total;

            total += tSizes[indexI][i-1];
        }

        if (debug > 2)
        {
            Info << " Load starts: " << tStarts[indexI] << endl;
            Info << " Load sizes: " << tSizes[indexI] << endl;
        }
    }

    // Set the argument list for each thread
    forAll(hdl, i)
    {
        // Size up the argument list
        hdl[i].setSize(6);

        // Set match tolerance
        hdl[i].set(0, &matchTol);

        // Set the skipMapping flag
        hdl[i].set(1, &skipMapping);

        // Set the start/size indices
        hdl[i].set(2, &(tStarts[0][i]));
        hdl[i].set(3, &(tSizes[0][i]));
        hdl[i].set(4, &(tStarts[1][i]));
        hdl[i].set(5, &(tSizes[1][i]));
    }

    // Prior to multi-threaded operation,
    // force calculation of demand-driven data.
    polyMesh::cells();
    primitiveMesh::cellCells();

    const polyBoundaryMesh& boundary = boundaryMesh();

    forAll(boundary, patchI)
    {
        boundary[patchI].faceFaces();
    }

    // Execute threads in linear sequence
    executeThreads(identity(nThreads), hdl, &computeMappingThread);
}


// Obtain map weighting factors
bool dynamicTopoFvMesh::computeWeights
(
    const label index,
    const labelList& mapCandidates,
    const labelListList& oldNeighbourList,
    const scalar mTol,
    const convexSetAlgorithm& algorithm,
    label& precisionAttempts,
    labelList& parents,
    scalarField& weights,
    vectorField& centres
) const
{
    if (parents.size() || weights.size() || centres.size())
    {
        FatalErrorIn
        (
            "\n\n"
            "void dynamicTopoFvMesh::computeWeights\n"
            "(\n"
            "    const label index,\n"
            "    const labelList& mapCandidates,\n"
            "    const labelListList& oldNeighbourList,\n"
            "    const scalar mTol,\n"
            "    const convexSetAlgorithm<T>& algorithm,\n"
            "    label& precisionAttempts,\n"
            "    labelList& parents,\n"
            "    scalarField& weights,\n"
            "    vectorField& centres\n"
            ") const\n"
        )
            << " Addressing has already been calculated." << nl
            << " Index: " << index << nl
            << " Type: "
            << (algorithm.dimension() == 2 ? "Face" : "Cell") << nl
            << " mapCandidates: " << mapCandidates << nl
            << " Parents: " << parents << nl
            << " Weights: " << weights << nl
            << " Centres: " << centres << nl
            << abort(FatalError);
    }

    bool changed;
    scalar matchTol = mTol;
    label nAttempts = 0, nIntersects = 0;

    // Figure out the patch offset
    label offset = -1;

    if (algorithm.dimension() == 2)
    {
        offset = boundaryMesh()[whichPatch(index)].start();
    }
    else
    if (algorithm.dimension() == 3)
    {
        offset = 0;
    }

    // Calculate the algorithm normFactor
    algorithm.computeNormFactor(index);

    // Maintain a check-list
    StaticHashTable<empty, label, Hash<label> > checked, skipped;

    // Loop and add intersections until nothing changes
    do
    {
        // Reset flag
        changed = false;

        // Fetch the set of candidates
        labelList checkList;

        if (nAttempts == 0)
        {
            checkList = mapCandidates;
        }
        else
        {
            checkList = checked.toc();
        }

        forAll(checkList, indexI)
        {
            labelList checkEntities;

            if (nAttempts == 0)
            {
                checkEntities = labelList(1, checkList[indexI] - offset);
            }
            else
            {
                checkEntities = oldNeighbourList[checkList[indexI]];
            }

            forAll(checkEntities, entityI)
            {
                label checkEntity = checkEntities[entityI];

                // Skip if this is already
                // on the checked / skipped list
                if
                (
                    (checked.found(checkEntity)) ||
                    (skipped.found(checkEntity))
                )
                {
                    continue;
                }

                bool intersect =
                (
                    algorithm.computeInsersection
                    (
                        index,
                        checkEntity + offset,
                        matchTol,
                        false
                    )
                );

                if (intersect)
                {
                    nIntersects++;

                    checked.insert(checkEntity, empty());

                    changed = true;
                }
                else
                {
                    // Add to the skipped list
                    skipped.insert(checkEntity, empty());
                }
            }
        }

        if (nAttempts == 0 && !changed)
        {
            // Need to setup a rescue mechanism.
            StaticHashTable<empty, label, Hash<label> > rescue;

            forAll(mapCandidates, cI)
            {
                rescue.insert(mapCandidates[cI] - offset, empty());
            }

            for (label level = 0; level < 10; level++)
            {
                labelList initList = rescue.toc();

                forAll(initList, fI)
                {
                    const labelList& ff = oldNeighbourList[initList[fI]];

                    forAll(ff, entityI)
                    {
                        rescue.insert(ff[entityI], empty());
                    }
                }
            }

            labelList finalList = rescue.toc();

            forAll(finalList, entityI)
            {
                label checkEntity = finalList[entityI];

                bool intersect =
                (
                    algorithm.computeInsersection
                    (
                        index,
                        checkEntity + offset,
                        matchTol,
                        false
                    )
                );

                if (intersect)
                {
                    nIntersects++;

                    checked.insert(checkEntity, empty());

                    changed = true;
                    break;
                }
            }

            if (!changed)
            {
                // No point in continuing further...
                break;
            }
        }

        nAttempts++;

        // Break out if we're taking too long
        if (nAttempts > 20)
        {
            break;
        }

    } while (changed);

    // Test weights for consistency
    bool consistent = algorithm.consistent(1e-13);
    bool normByWeights = false;

    if (!consistent)
    {
        // Inconsistent weights.
        switch (algorithm.dimension())
        {
            case 2:
            {
                // Check whether any edges lie on bounding curves.
                // These faces can have relaxed weights to account
                // for addressing into patches on the other side
                // of the curve.
                const labelList& fEdges = faceEdges_[index];

                forAll(fEdges, edgeI)
                {
                    if (checkBoundingCurve(fEdges[edgeI]))
                    {
                        consistent = true;
                    }
                }

                break;
            }

            case 3:
            {
                // Check whether any edges lie on boundary patches.
                // These cells can have relaxed weights to account
                // for mild convexity.
                const cell& cellToCheck = cells_[index];

                if (twoDMesh_)
                {
                    forAll(parents, cellI)
                    {
                        const cell& pCell = polyMesh::cells()[parents[cellI]];

                        forAll(pCell, faceI)
                        {
                            const face& pFace = polyMesh::faces()[pCell[faceI]];

                            if (pFace.size() == 3)
                            {
                                continue;
                            }

                            if (boundaryMesh().whichPatch(pCell[faceI]) > -1)
                            {
                                consistent = true;
                                break;
                            }
                        }

                        if (consistent)
                        {
                            break;
                        }
                    }
                }
                else
                {
                    forAll(cellToCheck, faceI)
                    {
                        const labelList& fE = faceEdges_[cellToCheck[faceI]];

                        forAll(fE, edgeI)
                        {
                            if (whichEdgePatch(fE[edgeI]) > -1)
                            {
                                consistent = true;
                                break;
                            }
                        }

                        if (consistent)
                        {
                            break;
                        }
                    }
                }

                break;
            }
        }

        if (consistent)
        {
            // Optionally output for post-processing
            if (debug > 4)
            {
                labelList uList = skipped.toc();

                // Renumber to global indices (for faces)
                forAll(uList, entityI)
                {
                    uList[entityI] += offset;
                }

                label pT = algorithm.dimension();

                // Populate lists
                algorithm.populateLists(parents, centres, weights);

                writeVTK("n_" + Foam::name(index), index, pT, false, true);
                writeVTK("m_" + Foam::name(index), parents, pT, true, true);
                writeVTK("u_" + Foam::name(index), uList, pT, true, true);

                // Write out intersections for post-processing
                forAll(parents, indexI)
                {
                    algorithm.computeInsersection
                    (
                        index,
                        parents[indexI],
                        matchTol,
                        true
                    );
                }
            }

            // Normalize by sum of weights instead
            normByWeights = true;
        }
        else
        if (precisionAttempts < 2)
        {
            // Could be a precision problem.
            // Recurse until consistency is obtained.
            matchTol *= 0.1;

            // Toggle higher precision
            algorithm.setHighPrecision();

            consistent =
            (
                computeWeights
                (
                    index,
                    mapCandidates,
                    oldNeighbourList,
                    matchTol,
                    algorithm,
                    ++precisionAttempts,
                    parents,
                    weights,
                    centres
                )
            );
        }
    }

    if (!consistent)
    {
        // Write out for post-processing
        labelList uList = skipped.toc();

        // Renumber to global indices (for faces)
        forAll(uList, entityI)
        {
            uList[entityI] += offset;
        }

        label pT = algorithm.dimension();

        // Normalize weights
        algorithm.normalize(normByWeights);

        // Populate lists
        algorithm.populateLists(parents, centres, weights);

        writeVTK("nE_" + Foam::name(index), index, pT, false, true);
        writeVTK("oE_" + Foam::name(index), mapCandidates, pT, true, true);
        writeVTK("mE_" + Foam::name(index), parents, pT, true, true);
        writeVTK("uE_" + Foam::name(index), uList, pT, true, true);

        // Write out intersections for post-processing
        forAll(parents, indexI)
        {
            algorithm.computeInsersection
            (
                index,
                parents[indexI],
                matchTol,
                true
            );
        }

        FatalErrorIn
        (
            "\n\n"
            "void dynamicTopoFvMesh::computeWeights\n"
            "(\n"
            "    const label index,\n"
            "    const labelList& mapCandidates,\n"
            "    const labelListList& oldNeighbourList,\n"
            "    const scalar mTol,\n"
            "    const convexSetAlgorithm<T>& algorithm,\n"
            "    label& precisionAttempts,\n"
            "    labelList& parents,\n"
            "    scalarField& weights,\n"
            "    vectorField& centres\n"
            ") const\n"
        )
            << "Encountered non-conservative weighting factors." << nl
            << " Index: " << index << nl
            << " Type: "
            << (algorithm.dimension() == 2 ? "Face" : "Cell") << nl
            << " mapCandidates: " << mapCandidates << nl
            << " nParents: " << parents.size() << nl
            << " nAttempts: " << nAttempts << nl
            << " precisionAttempts: " << precisionAttempts << nl
            << " matchTolerance: " << matchTol << nl
            << setprecision(16)
            << " Magnitude: " << algorithm.normFactor() << nl
            << " Norm Sum(Weights): " << sum(weights) << nl
            << " Norm Error: " << mag(1.0 - sum(weights)) << nl
            << " Parents: " << parents << nl
            << " Weights: " << weights
            << abort(FatalError);
    }

    // Return normalized weights,
    // but only if we're at the top
    // of the recursion stack,
    if (precisionAttempts)
    {
        precisionAttempts--;
    }
    else
    {
        // Normalize weights
        algorithm.normalize(normByWeights);

        // Revert precision, if necessary
        if (algorithm.highPrecision())
        {
            algorithm.unsetHighPrecision();
        }

        // Populate lists
        algorithm.populateLists(parents, centres, weights);

        if (debug > 2)
        {
            if (debug > 3)
            {
                label pT = algorithm.dimension();

                writeVTK("nE_" + Foam::name(index), index, pT, false, true);
                writeVTK("mE_" + Foam::name(index), parents, pT, true, true);

                // Write out intersections for post-processing
                forAll(parents, indexI)
                {
                    algorithm.computeInsersection
                    (
                        index,
                        parents[indexI],
                        matchTol,
                        true
                    );
                }
            }

            Info << " Index: " << index << nl
                 << setprecision(16)
                 << " Magnitude: " << algorithm.normFactor() << nl
                 << " Norm Error: " << mag(1.0 - sum(weights)) << nl
                 << " Sum(Weights): " << sum(weights) << nl
                 << " Weights: " << weights << nl
                 << endl;
        }
    }

    return consistent;
}


// Set fill-in mapping information for a particular cell
void dynamicTopoFvMesh::setCellMapping
(
    const label cIndex,
    const labelList& mapCells,
    bool addEntry
)
{
    if (addEntry)
    {
        if (debug > 3)
        {
            Info << "Inserting mapping cell: " << cIndex << nl
                 << " mapCells: " << mapCells
                 << endl;
        }

        // Insert index into the list, and overwrite if necessary
        label index = -1;

        forAll(cellsFromCells_, indexI)
        {
            if (cellsFromCells_[indexI].index() == cIndex)
            {
                index = indexI;
                break;
            }
        }

        if (index == -1)
        {
            meshOps::sizeUpList
            (
                objectMap(cIndex, labelList(0)),
                cellsFromCells_
            );
        }
        else
        {
            cellsFromCells_[index].masterObjects() = labelList(0);
        }
    }

    // Update cell-parents information
    labelHashSet masterCells;

    forAll(mapCells, cellI)
    {
        if (mapCells[cellI] < 0)
        {
            continue;
        }

        if (mapCells[cellI] < nOldCells_)
        {
            masterCells.insert(mapCells[cellI]);
        }
        else
        if (cellParents_.found(mapCells[cellI]))
        {
            const labelList& nParents = cellParents_[mapCells[cellI]];

            forAll(nParents, cI)
            {
                masterCells.insert(nParents[cI]);
            }
        }
    }

    cellParents_.set(cIndex, masterCells.toc());
}


// Set fill-in mapping information for a particular face
void dynamicTopoFvMesh::setFaceMapping
(
    const label fIndex,
    const labelList& mapFaces
)
{
    label patch = whichPatch(fIndex);

    if (debug > 3)
    {
        Info << "Inserting mapping face: " << fIndex << nl
             << " patch: " << patch << nl
             << " mapFaces: " << mapFaces
             << endl;
    }

    bool foundError = false;

    // Check to ensure that internal faces are not mapped
    // from any faces in the mesh
    if (patch == -1 && mapFaces.size())
    {
        foundError = true;
    }

    // Check to ensure that boundary faces map
    // only from other faces on the same patch
    if (patch > -1 && mapFaces.empty())
    {
        foundError = true;
    }

    if (foundError)
    {
        FatalErrorIn
        (
            "\n"
            "void dynamicTopoFvMesh::setFaceMapping\n"
            "(\n"
            "    const label fIndex,\n"
            "    const labelList& mapFaces\n"
            ")"
        )
            << nl << " Incompatible mapping. " << nl
            << "  Possible reasons: " << nl
            << "    1. No mapping specified for a boundary face; " << nl
            << "    2. Mapping specified for an internal face, " << nl
            << "       when none was expected." << nl << nl
            << " Face: " << fIndex << nl
            << " Patch: " << patch << nl
            << " Owner: " << owner_[fIndex] << nl
            << " Neighbour: " << neighbour_[fIndex] << nl
            << " mapFaces: " << mapFaces << nl
            << abort(FatalError);
    }

    // Insert addressing into the list, and overwrite if necessary
    label index = -1;

    forAll(facesFromFaces_, indexI)
    {
        if (facesFromFaces_[indexI].index() == fIndex)
        {
            index = indexI;
            break;
        }
    }

    if (index == -1)
    {
        meshOps::sizeUpList
        (
            objectMap(fIndex, labelList(0)),
            facesFromFaces_
        );
    }
    else
    {
        facesFromFaces_[index].masterObjects() = labelList(0);
    }

    // For internal faces, set dummy maps / weights, and bail out
    if (patch == -1)
    {
        return;
    }

    // Update face-parents information
    labelHashSet masterFaces;

    forAll(mapFaces, faceI)
    {
        if (mapFaces[faceI] < 0)
        {
            continue;
        }

        if (mapFaces[faceI] < nOldFaces_)
        {
            masterFaces.insert(mapFaces[faceI]);
        }
        else
        if (faceParents_.found(mapFaces[faceI]))
        {
            const labelList& nParents = faceParents_[mapFaces[faceI]];

            forAll(nParents, fI)
            {
                masterFaces.insert(nParents[fI]);
            }
        }
    }

    faceParents_.set(fIndex, masterFaces.toc());
}


} // End namespace Foam

// ************************************************************************* //
