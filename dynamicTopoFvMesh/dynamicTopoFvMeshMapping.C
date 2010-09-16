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
            computeCellWeights
            (
                cIndex,
                cellParents_[cIndex],
                matchTol,
                precisionAttempts,
                cellsFromCells_[cellI].masterObjects(),
                cellWeights_[cellI],
                cellCentres_[cellI]
            );
        }
    }

    // Compute face mapping
    for (label faceI = faceStart; faceI < (faceStart + faceSize); faceI++)
    {
        label precisionAttempts = 0;
        label fIndex = facesFromFaces_[faceI].index();

        // Skip mapping for internal faces.
        if (whichPatch(fIndex) == -1)
        {
            // Set dummy masters, so that the conventional
            // faceMapper doesn't incur a seg-fault.
            facesFromFaces_[faceI].masterObjects() = labelList(1, 0);
            continue;
        }

        // Obtain weighting factors for this face.
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
            computeFaceWeights
            (
                fIndex,
                faceParents_[fIndex],
                matchTol,
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


// Obtain map weighting factors for a face
bool dynamicTopoFvMesh::computeFaceWeights
(
    const label fIndex,
    const labelList& mapCandidates,
    const scalar mTol,
    label& precisionAttempts,
    labelList& parents,
    scalarField& weights,
    vectorField& centres
) const
{
    // Determine which patch this is...
    label patchIndex = whichPatch(fIndex);
    const polyBoundaryMesh& boundary = boundaryMesh();

    if (parents.size() || weights.size() || centres.size() || patchIndex == -1)
    {
        FatalErrorIn
        (
            "\n\n"
            "void dynamicTopoFvMesh::computeFaceWeights\n"
            "(\n"
            "    const label fIndex,\n"
            "    const labelList& mapCandidates,\n"
            "    const scalar mTol,\n"
            "    label& precisionAttempts,\n"
            "    labelList& parents,\n"
            "    scalarField& weights,\n"
            "    vectorField& centres\n"
            ") const\n"
        )
            << " Addressing has already been calculated." << nl
            << " Face: " << fIndex << nl
            << " mapCandidates: " << mapCandidates << nl
            << " Patch: " << patchIndex << nl
            << " Parents: " << parents << nl
            << " Weights: " << weights << nl
            << " Centres: " << centres << nl
            << abort(FatalError);
    }

    bool changed, output = false;
    label nAttempts = 0, nIntersects = 0;
    scalar matchTol = mTol, sumAreas = 0.0;

    // Fetch the area / centre of the cell
    vector fNormal =
    (
        meshOps::faceNormal
        (
            faces_[fIndex],
            oldPoints_
        )
    );

    scalar faceArea = mag(fNormal);

    // Fetch old patch start
    label patchStart = boundary[patchIndex].start();

    StaticHashTable<empty, label, Hash<label> > checkedFaces, skippedFaces;

    // VectorField of intersection points
    vectorField tP(0);

    // Fetch connectivity from the old mesh.
    const labelListList& oldFaceFaces = boundary[patchIndex].faceFaces();

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
            checkList = checkedFaces.toc();
        }

        forAll(checkList, indexI)
        {
            labelList checkFaces;

            if (nAttempts == 0)
            {
                checkFaces = labelList(1, checkList[indexI] - patchStart);
            }
            else
            {
                checkFaces = oldFaceFaces[checkList[indexI]];
            }

            forAll(checkFaces, faceI)
            {
                label checkFace = checkFaces[faceI];

                // Skip if this is already
                // on the checked / skipped list
                if
                (
                    (checkedFaces.found(checkFace)) ||
                    (skippedFaces.found(checkFace))
                )
                {
                    continue;
                }

                bool intersect =
                (
                    faceIntersection
                    (
                        fIndex,
                        checkFace + patchStart,
                        matchTol,
                        tP
                    )
                );

                if (intersect)
                {
                    scalar cvxArea = 0.0;
                    vector cvxCentre = vector::zero;

                    // We need a reference normal. Use the new face.
                    vector refNorm =
                    (
                        meshOps::faceNormal
                        (
                            faces_[fIndex],
                            oldPoints_
                        )
                    );

                    refNorm /= mag(refNorm) + VSMALL;

                    // Compute weights
                    meshOps::convexSetArea
                    (
                        fIndex,
                        checkFace + patchStart,
                        tP,
                        refNorm,
                        cvxArea,
                        cvxCentre,
                        output
                    );

                    // Size-up lists
                    meshOps::sizeUpList(checkFace + patchStart, parents);
                    meshOps::sizeUpList(cvxArea, weights);
                    meshOps::sizeUpList(cvxCentre, centres);

                    // Accumulate area
                    sumAreas += cvxArea;
                    nIntersects++;

                    checkedFaces.insert(checkFace, empty());

                    changed = true;
                }
                else
                {
                    // Add to the skipped list
                    skippedFaces.insert(checkFace, empty());
                }
            }
        }

        if (nAttempts == 0 && !changed)
        {
            // Need to setup a rescue mechanism.
            StaticHashTable<empty, label, Hash<label> > rescue;

            forAll(mapCandidates, cI)
            {
                rescue.insert(mapCandidates[cI] - patchStart, empty());
            }

            for (label level = 0; level < 10; level++)
            {
                labelList initList = rescue.toc();

                forAll(initList, fI)
                {
                    const labelList& ff = oldFaceFaces[initList[fI]];

                    forAll(ff, faceI)
                    {
                        rescue.insert(ff[faceI], empty());
                    }
                }
            }

            labelList finalList = rescue.toc();

            forAll(finalList, faceI)
            {
                label checkFace = finalList[faceI];

                bool intersect =
                (
                    faceIntersection
                    (
                        fIndex,
                        checkFace + patchStart,
                        matchTol,
                        tP
                    )
                );

                if (intersect)
                {
                    scalar cvxArea = 0.0;
                    vector cvxCentre = vector::zero;

                    // We need a reference normal. Use the new face.
                    vector refNorm =
                    (
                        meshOps::faceNormal
                        (
                            faces_[fIndex],
                            oldPoints_
                        )
                    );

                    refNorm /= mag(refNorm) + VSMALL;

                    // Compute weights
                    meshOps::convexSetArea
                    (
                        fIndex,
                        checkFace + patchStart,
                        tP,
                        refNorm,
                        cvxArea,
                        cvxCentre,
                        output
                    );

                    // Size-up lists
                    meshOps::sizeUpList(checkFace + patchStart, parents);
                    meshOps::sizeUpList(cvxArea, weights);
                    meshOps::sizeUpList(cvxCentre, centres);

                    // Accumulate area
                    sumAreas += cvxArea;
                    nIntersects++;

                    checkedFaces.insert(checkFace, empty());

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
    bool consistent = false;

    scalar consistency = mag(1.0 - sum(weights/faceArea));

    // Test weights for consistency
    if (consistency > 1e-10)
    {
        // Inconsistent weights. Check whether any edges
        // lie on bounding curves. These faces can have
        // relaxed weights to account for addressing into
        // patches on the other side of the curve.
        if (weights.size())
        {
            const labelList& fEdges = faceEdges_[fIndex];

            forAll(fEdges, edgeI)
            {
                if (checkBoundingCurve(fEdges[edgeI]))
                {
                    consistent = true;
                }
            }
        }

        if (consistent)
        {
            // Normalize by sum of weights instead
            faceArea = sum(weights);
        }
        else
        if (precisionAttempts < 12)
        {
            // Could be a precision problem.
            // Recurse until consistency is obtained.
            matchTol *= 0.1;

            parents.clear();
            weights.clear();
            centres.clear();

            consistent =
            (
                computeFaceWeights
                (
                    fIndex,
                    mapCandidates,
                    matchTol,
                    ++precisionAttempts,
                    parents,
                    weights,
                    centres
                )
            );
        }
    }
    else
    {
        // Weights are consistent
        consistent = true;
    }

    if (!consistent)
    {
        // Write out for post-processing
        labelList uFaces = skippedFaces.toc();

        // Renumber to global indices
        forAll(uFaces, faceI)
        {
            uFaces[faceI] += patchStart;
        }

        writeVTK("nFace_" + Foam::name(fIndex), fIndex, 2, false, true);
        writeVTK("oFace_" + Foam::name(fIndex), mapCandidates, 2, true, true);
        writeVTK("mFace_" + Foam::name(fIndex), parents, 2, true, true);
        writeVTK("uFace_" + Foam::name(fIndex), uFaces, 2, true, true);

        // Write out intersections for post-processing
        forAll(parents, faceI)
        {
            faceIntersection(fIndex, parents[faceI], matchTol, tP, true);
        }

        FatalErrorIn
        (
            "\n\n"
            "bool dynamicTopoFvMesh::computeFaceWeights\n"
            "(\n"
            "    const label fIndex,\n"
            "    const labelList& mapCandidates,\n"
            "    const scalar mTol,\n"
            "    label& precisionAttempts,\n"
            "    labelList& parents,\n"
            "    scalarField& weights,\n"
            "    vectorField& centres\n"
            ") const\n"
        )
            << "Encountered non-conservative weighting factors." << nl
            << " Face: " << fIndex << ":: " << faces_[fIndex] << nl
            << " Patch: " << boundaryMesh()[patchIndex].name() << nl
            << " mapCandidates: " << mapCandidates << nl
            << " nParents: " << parents.size() << nl
            << " nAttempts: " << nAttempts << nl
            << " precisionAttempts: " << precisionAttempts << nl
            << " matchTolerance: " << matchTol << nl
            << " nFaces: " << nFaces_ << nl
            << " nOldFaces: " << nOldFaces_ << nl
            << " nInternalFaces: " << nInternalFaces_ << nl
            << " nOldInternalFaces: " << nOldInternalFaces_ << nl
            << setprecision(16)
            << " Face area: " << faceArea << nl
            << " Sum(Weights): " << sum(weights) << nl
            << " Error: " << (faceArea - sum(weights)) << nl
            << " Norm Sum(Weights): " << sum(weights/faceArea) << nl
            << " Norm Error: " << consistency << nl
            << " Parents: " << parents << nl
            << " Weights: " << (weights/faceArea)
            << abort(FatalError);
    }

    if (debug > 2)
    {
        Info << " Face: " << fIndex << nl
             << setprecision(16)
             << " Face area: " << faceArea << nl
             << " Sum(Weights): " << sum(weights) << nl
             << " Error: " << (faceArea - sum(weights)) << nl
             << " Norm Sum(Weights): " << sum(weights/faceArea) << nl
             << " Norm Error: " << consistency
             << endl;
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
        weights /= faceArea;
    }

    return consistent;
}


// Obtain map weighting factors for a cell
bool dynamicTopoFvMesh::computeCellWeights
(
    const label cIndex,
    const labelList& mapCandidates,
    const scalar mTol,
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
            "void dynamicTopoFvMesh::computeCellWeights\n"
            "(\n"
            "    const label cIndex,\n"
            "    const labelList& mapCandidates,\n"
            "    const scalar mTol,\n"
            "    label& precisionAttempts,\n"
            "    labelList& parents,\n"
            "    scalarField& weights,\n"
            "    vectorField& centres\n"
            ") const\n"
        )
            << " Addressing has already been calculated." << nl
            << " Cell: " << cIndex << nl
            << " mapCandidates: " << mapCandidates << nl
            << " Parents: " << parents << nl
            << " Weights: " << weights << nl
            << " Centres: " << centres << nl
            << abort(FatalError);
    }

    bool changed, output = false;
    label nAttempts = 0, nIntersects = 0;
    scalar matchTol = mTol, sumVols = 0.0;

    // Fetch the volume / centre of the cell
    scalar cellVolume = 0.0;
    vector cellCentre = vector::zero;

    meshOps::cellCentreAndVolume
    (
        cIndex,
        oldPoints_,
        faces_,
        cells_,
        owner_,
        cellCentre,
        cellVolume
    );

    // Maintain a check-list
    StaticHashTable<empty, label, Hash<label> > checkedCells, skippedCells;

    // VectorField of intersection points
    vectorField tP(0);

    // Fetch connectivity from the old mesh.
    const labelListList& oldCellCells = polyMesh::cellCells();

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
            checkList = checkedCells.toc();
        }

        forAll(checkList, indexI)
        {
            labelList checkCells;

            if (nAttempts == 0)
            {
                checkCells = labelList(1, checkList[indexI]);
            }
            else
            {
                checkCells = oldCellCells[checkList[indexI]];
            }

            forAll(checkCells, cellI)
            {
                label checkCell = checkCells[cellI];

                // Skip if this is already
                // on the checked / skipped list
                if
                (
                    (checkedCells.found(checkCell)) ||
                    (skippedCells.found(checkCell))
                )
                {
                    continue;
                }

                bool intersect =
                (
                    cellIntersection
                    (
                        cIndex,
                        checkCell,
                        matchTol,
                        tP
                    )
                );

                if (intersect)
                {
                    scalar cvxVolume = 0.0;
                    vector cvxCentre = vector::zero;

                    // Compute weights
                    meshOps::convexSetVolume
                    (
                        cIndex,
                        checkCell,
                        tP,
                        cvxVolume,
                        cvxCentre,
                        output
                    );

                    // Size-up lists
                    meshOps::sizeUpList(checkCell, parents);
                    meshOps::sizeUpList(cvxVolume, weights);
                    meshOps::sizeUpList(cvxCentre, centres);

                    // Accumulate volume
                    sumVols += cvxVolume;
                    nIntersects++;

                    checkedCells.insert(checkCell, empty());

                    changed = true;
                }
                else
                {
                    // Add to the skipped list
                    skippedCells.insert(checkCell, empty());
                }
            }
        }

        if (nAttempts == 0 && !changed)
        {
            // Need to setup a rescue mechanism.
            StaticHashTable<empty, label, Hash<label> > rescue;

            forAll(mapCandidates, cI)
            {
                rescue.insert(mapCandidates[cI], empty());
            }

            // Agglomerate a larger set of cells
            for (label level = 0; level < 10; level++)
            {
                labelList initList = rescue.toc();

                forAll(initList, cI)
                {
                    const labelList& cc = oldCellCells[initList[cI]];

                    forAll(cc, cellI)
                    {
                        rescue.insert(cc[cellI], empty());
                    }
                }
            }

            labelList finalList = rescue.toc();

            forAll(finalList, cellI)
            {
                label checkCell = finalList[cellI];

                bool intersect =
                (
                    cellIntersection
                    (
                        cIndex,
                        checkCell,
                        matchTol,
                        tP
                    )
                );

                if (intersect)
                {
                    scalar cvxVolume = 0.0;
                    vector cvxCentre = vector::zero;

                    // Compute weights
                    meshOps::convexSetVolume
                    (
                        cIndex,
                        checkCell,
                        tP,
                        cvxVolume,
                        cvxCentre,
                        output
                    );

                    // Size-up lists
                    meshOps::sizeUpList(checkCell, parents);
                    meshOps::sizeUpList(cvxVolume, weights);
                    meshOps::sizeUpList(cvxCentre, centres);

                    // Accumulate volume
                    sumVols += cvxVolume;
                    nIntersects++;

                    checkedCells.insert(checkCell, empty());

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
    bool consistent = false;

    scalar consistency = mag(1.0 - sum(weights/cellVolume));

    if (consistency > 1e-10)
    {
        // Inconsistent weights. Check whether any edges
        // lie on boundary patches. These cells can have
        // relaxed weights to account for mild convexity.
        if (weights.size())
        {
            const cell& cellToCheck = cells_[cIndex];

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
                    const labelList& fEdges = faceEdges_[cellToCheck[faceI]];

                    forAll(fEdges, edgeI)
                    {
                        if (whichEdgePatch(fEdges[edgeI]) > -1)
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
        }

        if (consistent)
        {
            // Normalize by sum of weights instead
            cellVolume = sum(weights);
        }
        else
        if (precisionAttempts < 12)
        {
            // Could be a precision problem.
            // Recurse until consistency is obtained.
            matchTol *= 0.1;

            parents.clear();
            weights.clear();
            centres.clear();

            consistent =
            (
                computeCellWeights
                (
                    cIndex,
                    mapCandidates,
                    matchTol,
                    ++precisionAttempts,
                    parents,
                    weights,
                    centres
                )
            );
        }
    }
    else
    {
        // Weights are consistent
        consistent = true;
    }

    if (!consistent)
    {
        // Write out for post-processing
        labelList uCells = skippedCells.toc();

        writeVTK("nCell_" + Foam::name(cIndex), cIndex, 3, false, true);
        writeVTK("oCell_" + Foam::name(cIndex), mapCandidates, 3, true, true);
        writeVTK("mCell_" + Foam::name(cIndex), parents, 3, true, true);
        writeVTK("uCell_" + Foam::name(cIndex), uCells, 3, true, true);

        // Write out intersections for post-processing
        forAll(parents, cellI)
        {
            cellIntersection(cIndex, parents[cellI], matchTol, tP, true);
        }

        FatalErrorIn
        (
            "\n\n"
            "bool dynamicTopoFvMesh::computeCellWeights\n"
            "(\n"
            "    const label cIndex,\n"
            "    const labelList& mapCandidates,\n"
            "    const scalar mTol,\n"
            "    label& precisionAttempts,\n"
            "    labelList& parents,\n"
            "    scalarField& weights,\n"
            "    vectorField& centres\n"
            ") const\n"
        )
            << "Encountered non-conservative weighting factors." << nl
            << " Cell: " << cIndex << nl
            << " mapCandidates: " << mapCandidates << nl
            << " nParents: " << parents.size() << nl
            << " nIntersects: " << nIntersects << nl
            << " nAttempts: " << nAttempts << nl
            << " precisionAttempts: " << precisionAttempts << nl
            << " matchTolerance: " << matchTol << nl
            << " nCells: " << nCells_ << nl
            << " nOldCells: " << nOldCells_ << nl
            << setprecision(16)
            << " Cell volume: " << cellVolume << nl
            << " Sum(Weights): " << sum(weights) << nl
            << " Error: " << (cellVolume - sum(weights)) << nl
            << " Norm Sum(Weights): " << sum(weights/cellVolume) << nl
            << " Norm Error: " << consistency << nl
            << " Parents: " << parents << nl
            << " Weights: " << (weights/cellVolume)
            << abort(FatalError);
    }

    if (debug > 2)
    {
        Info << " Cell: " << cIndex << nl
             << setprecision(16)
             << " Cell volume: " << cellVolume << nl
             << " Sum(Weights): " << sum(weights) << nl
             << " Error: " << (cellVolume - sum(weights)) << nl
             << " Norm Sum(Weights): " << sum(weights/cellVolume) << nl
             << " Norm Error: " << consistency
             << endl;
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
        weights /= cellVolume;
    }

    return consistent;
}


// Return the intersection points between faces in old/new meshes
bool dynamicTopoFvMesh::faceIntersection
(
    const label newFaceIndex,
    const label oldFaceIndex,
    const scalar matchTol,
    vectorField& tP,
    bool output
) const
{
    // Reset inputs
    tP.clear();

    // Fetch references for each mesh
    const face& fromFace = polyMesh::faces()[oldFaceIndex];
    const face& toFace = faces_[newFaceIndex];

    // Obtain face centre and projection normal
    vector xf = meshOps::faceCentre(toFace, oldPoints_);
    vector nf = meshOps::faceNormal(toFace, oldPoints_);

    nf /= mag(nf) + VSMALL;

    // Track all possible intersections from here on.
    label nInts = 0;
    Map<vector> intersections;
    vector intPoint = vector::zero;

    // Topologically check for common points,
    // and project uniques ones on to the face plane
    Map<label> projPoints;
    Map<labelList> commonPoints;
    vectorField projections(fromFace.size(), vector::zero);

    forAll(fromFace, pointI)
    {
        label fromPoint = fromFace[pointI];
        label pIndex = findIndex(toFace, fromPoint);

        vector r = polyMesh::points()[fromPoint];

        if (pIndex == -1)
        {
            // Project this point on to the toFace plane.
            projections[pointI] = xf + ((r - xf) - ((r - xf) & nf)*nf);

            projPoints.insert(fromPoint, pointI);
        }
        else
        {
            // If this point was modified by a collapse
            // to an edge mid-point, it can't be a common point.
            if (modPoints_.found(toFace[pIndex]))
            {
                // Project this point on to the toFace plane.
                projections[pointI] = xf + ((r - xf) - ((r - xf) & nf)*nf);

                projPoints.insert(fromPoint, pointI);
            }
            else
            {
                commonPoints.insert(toFace[pIndex], labelList(0));

                projections[pointI] = r;

                intersections.set(++nInts, r);
            }
        }
    }

    // Add points if they resulted from
    // bisections of old face edges.
    forAll(toFace, pointI)
    {
        label pIndex = toFace[pointI];

        if (pIndex >= nOldPoints_)
        {
            // Check pointsFromPoints info
            label index = -1;

            forAll(pointsFromPoints_, indexI)
            {
                if (pointsFromPoints_[indexI].index() == pIndex)
                {
                    index = indexI;
                    break;
                }
            }

            const labelList& mObj = pointsFromPoints_[index].masterObjects();

            // Check if the old face contains all master points
            bool allMaster = true;

            forAll(mObj, pointJ)
            {
                if (findIndex(fromFace, mObj[pointJ]) == -1)
                {
                    allMaster = false;
                    break;
                }
            }

            if (allMaster)
            {
                commonPoints.insert(toFace[pointI], mObj);

                intersections.set(++nInts, oldPoints_[toFace[pointI]]);
            }
        }
    }

    // If all points are common, this is identical to the old face.
    if (nInts == fromFace.size())
    {
        // Copy intersections
        tP.setSize(nInts, vector::zero);

        nInts = 0;

        forAllConstIter(Map<vector>, intersections, pI)
        {
            tP[nInts++] = pI();
        }

        if (debug > 3 || output)
        {
            if (meshOps::checkPointNearness(tP, 1e-20))
            {
                writeVTK(Foam::name(newFaceIndex),newFaceIndex,2,false,true);
                writeVTK(Foam::name(oldFaceIndex),oldFaceIndex,2,true,true);
            }

            meshOps::writeVTK
            (
                (*this),
                "ccSet_"
              + Foam::name(newFaceIndex)
              + '<' + Foam::name(oldFaceIndex) + '>',
                tP.size(),
                tP.size(),
                tP.size(),
                tP
            );
        }

        return true;
    }

    // Check for point-segment intersections
    bool pointIntersections = false;

    forAll(fromFace, pointI)
    {
        label fromPoint = fromFace[pointI];

        if (commonPoints.found(fromPoint))
        {
            continue;
        }

        const point& checkPoint = projections[pointI];

        // Loop through all new edges, and find possible intersections
        // with (projections of) old face points,
        forAll(toFace, pointJ)
        {
            edge toEdge = toFace.faceEdge(pointJ);

            if
            (
                meshOps::pointSegmentIntersection
                (
                    linePointRef
                    (
                        oldPoints_[toEdge.start()],
                        oldPoints_[toEdge.end()]
                    ),
                    checkPoint,
                    matchTol
                )
            )
            {
                commonPoints.insert
                (
                    fromPoint,
                    labelList(toEdge)
                );

                intersections.set(++nInts, checkPoint);

                pointIntersections = true;
            }
        }
    }

    forAll(toFace, pointI)
    {
        label toPoint = toFace[pointI];

        if (commonPoints.found(toPoint))
        {
            continue;
        }

        const point& checkPoint = oldPoints_[toPoint];

        forAll(fromFace, pointJ)
        {
            label nextJ = fromFace.fcIndex(pointJ);
            edge fromEdge = fromFace.faceEdge(pointJ);

            if
            (
                meshOps::pointSegmentIntersection
                (
                    linePointRef
                    (
                        projections[pointJ],
                        projections[nextJ]
                    ),
                    checkPoint,
                    matchTol
                )
             || meshOps::pointSegmentIntersection
                (
                    linePointRef
                    (
                        polyMesh::points()[fromEdge.start()],
                        polyMesh::points()[fromEdge.end()]
                    ),
                    checkPoint,
                    matchTol
                )
            )
            {
                commonPoints.insert
                (
                    toPoint,
                    labelList(fromEdge)
                );

                intersections.set(++nInts, checkPoint);

                pointIntersections = true;
            }
        }
    }

    if (pointIntersections && debug > 3)
    {
        Info << "Point Intersections exist: " << nl
             << " newFaceIndex: " << newFaceIndex
             << " oldFaceIndex: " << oldFaceIndex
             << endl;
    }

    if (fromFace.size() == 3 && toFace.size() == 3)
    {
        // Perform tests specific to triangular faces

        // Check whether any old projections are within
        // the new face. Count these as 'intersections'.
        forAll(fromFace, pointI)
        {
            label fromPoint = fromFace[pointI];

            if (commonPoints.found(fromPoint))
            {
                // Only skip for shared-points.
                // If the point-position was modified
                // due to a collapse, then this point
                // could be inside the new face.
                if (commonPoints[fromPoint].empty())
                {
                    continue;
                }
            }

            const point& checkPoint = projections[pointI];

            if
            (
                meshOps::pointInTriFace
                (
                    triPointRef
                    (
                        oldPoints_[toFace[0]],
                        oldPoints_[toFace[1]],
                        oldPoints_[toFace[2]]
                    ),
                    checkPoint
                )
            )
            {
                intersections.set(++nInts, checkPoint);
            }
        }

        // Check whether and new points are within
        // projected old faces. Count these as 'intersections'.
        forAll(toFace, pointI)
        {
            label toPoint = toFace[pointI];

            if (commonPoints.found(toPoint))
            {
                continue;
            }

            const point& checkPoint = oldPoints_[toPoint];

            if
            (
                meshOps::pointInTriFace
                (
                    triPointRef
                    (
                        projections[0],
                        projections[1],
                        projections[2]
                    ),
                    checkPoint
                )
            )
            {
                intersections.set(++nInts, checkPoint);
            }
        }

        // Loop through all new edges, and find possible intersections
        // with (projections of) old face edges,
        forAll(toFace, pointI)
        {
            edge toEdge = toFace.faceEdge(pointI);
            label nextLabel = toFace.nextLabel(pointI);

            forAll(fromFace, pointJ)
            {
                label nextJ = fromFace.fcIndex(pointJ);
                edge fromEdge = fromFace.faceEdge(pointJ);

                // Form an edge-pair
                Pair<edge> edgePair(fromEdge, toEdge);

                bool disableCheck = false;

                // Check edges topologically
                if (edgePair.first() == edgePair.second())
                {
                    const edge& checkEdge = edgePair.first();

                    // Check if points were modified by a collapse.
                    // If both were modified, continue with check.
                    if
                    (
                        !modPoints_.found(checkEdge.start()) &&
                        !modPoints_.found(checkEdge.end())
                    )
                    {
                        disableCheck = true;
                    }

                    // Skip shared points
                    if (commonPoints.found(checkEdge.start()))
                    {
                        if (commonPoints[checkEdge.start()].empty())
                        {
                            disableCheck = true;
                        }
                    }

                    if (commonPoints.found(checkEdge.end()))
                    {
                        if (commonPoints[checkEdge.end()].empty())
                        {
                            disableCheck = true;
                        }
                    }
                }
                else
                {
                    // Check for common vertices
                    label cV = edgePair.first().commonVertex(edgePair.second());

                    if (cV > -1)
                    {
                        // If this point was modified by a collapse
                        // to an edge mid-point, it can't be a common point.
                        // So, allow the check to continue.
                        if (!modPoints_.found(cV))
                        {
                            disableCheck = true;
                        }

                        // Skip shared points
                        if (commonPoints.found(cV))
                        {
                            if (commonPoints[cV].empty())
                            {
                                disableCheck = true;
                            }
                        }
                    }
                }

                if (disableCheck)
                {
                    continue;
                }

                // Also check for bisection / point-on-edge cases
                bool foundPointOnEdge = false;

                forAll(edgePair, indexI)
                {
                    const edge thisEdge = edgePair[indexI];

                    const edge otherEdge =
                    (
                        (thisEdge == edgePair.first()) ?
                        edgePair.second() : edgePair.first()
                    );

                    forAll(otherEdge, pointI)
                    {
                        label pIndex = otherEdge[pointI];

                        if (commonPoints.found(pIndex))
                        {
                            // Fetch masterObjects
                            const labelList& mObj = commonPoints[pIndex];

                            // Skip shared-points.
                            if (mObj.size())
                            {
                                // Check if the old edge
                                // contains all master points
                                bool allMaster = true;

                                forAll(mObj, pointJ)
                                {
                                    if (findIndex(thisEdge, mObj[pointJ]) == -1)
                                    {
                                        allMaster = false;
                                        break;
                                    }
                                }

                                if (allMaster)
                                {
                                    foundPointOnEdge = true;
                                }
                            }
                        }

                        if (foundPointOnEdge)
                        {
                            break;
                        }
                    }

                    if (foundPointOnEdge)
                    {
                        break;
                    }
                }

                if (foundPointOnEdge)
                {
                    continue;
                }

                point p1 = projections[pointJ];
                point p2 = projections[nextJ];
                point p3 = oldPoints_[toFace[pointI]];
                point p4 = oldPoints_[nextLabel];

                // Compute edge normal and tangent-to-edge
                vector te = (p4 - p3);
                vector n = (nf ^ te);
                n /= mag(n) + VSMALL;

                // Compute uValues
                scalar numOld = n & (p3 - p1);
                scalar denOld = n & (p2 - p1);

                // Check if the edges are parallel
                if (mag(denOld) < VSMALL)
                {
                    continue;
                }

                scalar tolerance = (matchTol * mag(p2 - p1));

                scalar u = (numOld / denOld);
                vector checkPoint = p1 + u*(p2 - p1);
                scalar v = (te & (checkPoint - p3)) / (te & te);

                // Check for intersection along lines.
                if
                (
                    ((u > tolerance) && (u < (1.0 - tolerance))) &&
                    ((v > tolerance) && (v < (1.0 - tolerance)))
                )
                {
                    intersections.set(++nInts, checkPoint);
                }
            }
        }
    }
    else
    if (fromFace.size() == 4 && toFace.size() == 4)
    {
        // Perform tests specific to quad faces
        notImplemented
        (
            "\n\n"
            "bool dynamicTopoFvMesh::faceIntersection\n"
            "(\n"
            "    const label newFaceIndex,\n"
            "    const label oldFaceIndex,\n"
            "    const scalar matchTol,\n"
            "    vectorField& tP\n"
            ") const\n"
        );
    }
    else
    {
        FatalErrorIn
        (
            "\n\n"
            "bool dynamicTopoFvMesh::faceIntersection\n"
            "(\n"
            "    const label newFaceIndex,\n"
            "    const label oldFaceIndex,\n"
            "    const scalar matchTol,\n"
            "    vectorField& tP\n"
            ") const\n"
        )
            << " Invalid face pair: " << nl
            << " Old face: "
            << oldFaceIndex << "::" << fromFace << nl
            << " New face: "
            << newFaceIndex << "::" << toFace << nl
            << abort(FatalError);
    }

    // Copy intersections
    tP.setSize(nInts, vector::zero);

    nInts = 0;

    forAllConstIter(Map<vector>, intersections, pI)
    {
        tP[nInts++] = pI();
    }

    // Check for concurrent points.
    if (debug > 3 || output)
    {
        if (meshOps::checkPointNearness(tP, 1e-20))
        {
            writeVTK(Foam::name(newFaceIndex),newFaceIndex,2,false,true);
            writeVTK(Foam::name(oldFaceIndex),oldFaceIndex,2,true,true);
        }

        meshOps::writeVTK
        (
            (*this),
            "ccSet_"
          + Foam::name(newFaceIndex)
          + '<' + Foam::name(oldFaceIndex) + '>',
            tP.size(),
            tP.size(),
            tP.size(),
            tP
        );
    }

    // Found a convex set of points.
    if (nInts >= 3)
    {
        return true;
    }

    // Does not intersect
    return false;
}


// Return the intersection points between cells in old/new meshes
bool dynamicTopoFvMesh::cellIntersection
(
    const label newCellIndex,
    const label oldCellIndex,
    const scalar matchTol,
    vectorField& tP,
    bool output
) const
{
    // Reset inputs
    tP.clear();

    // Assume XY plane here for 2D meshes
    vector planeNormal = vector(0,0,1);

    // Fetch references for each mesh
    const cell& fromCell = polyMesh::cells()[oldCellIndex];
    const edgeList fromCellEdges = fromCell.edges(polyMesh::faces());
    const labelList fromCellPoints = fromCell.labels(polyMesh::faces());

    const cell& toCell = cells_[newCellIndex];
    const edgeList toCellEdges = toCell.edges(faces_);
    const labelList toCellPoints = toCell.labels(faces_);

    // Track all possible intersections from here on.
    label nInts = 0;
    Map<vector> intersections;
    vector intPoint = vector::zero;

    // Topologically check for common points
    Map<labelList> commonPoints;

    forAll(fromCellPoints, pointI)
    {
        label fromPoint = fromCellPoints[pointI];
        label pIndex = findIndex(toCellPoints, fromPoint);

        if (pIndex > -1)
        {
            // If this point was modified by a collapse
            // to an edge mid-point, it can't be a common point.
            if (!modPoints_.found(toCellPoints[pIndex]))
            {
                commonPoints.insert(toCellPoints[pIndex], labelList(0));

                intersections.set(++nInts, oldPoints_[toCellPoints[pIndex]]);
            }
        }
    }

    // Add points if they resulted from
    // bisections of old cell edges.
    forAll(toCellPoints, pointI)
    {
        label pIndex = toCellPoints[pointI];

        if (pIndex >= nOldPoints_)
        {
            // Check pointsFromPoints info
            label index = -1;

            forAll(pointsFromPoints_, indexI)
            {
                if (pointsFromPoints_[indexI].index() == pIndex)
                {
                    index = indexI;
                    break;
                }
            }

            const labelList& mObj = pointsFromPoints_[index].masterObjects();

            // Check if the old cell contains all master points
            bool allMaster = true;

            forAll(mObj, pointJ)
            {
                if (findIndex(fromCellPoints, mObj[pointJ]) == -1)
                {
                    allMaster = false;
                    break;
                }
            }

            if (allMaster)
            {
                commonPoints.insert(toCellPoints[pointI], mObj);

                intersections.set(++nInts, oldPoints_[toCellPoints[pointI]]);
            }
        }
    }

    // If all points are common, this is identical to the old cell.
    if (nInts == fromCellPoints.size())
    {
        // Copy intersections
        tP.setSize(nInts, vector::zero);

        nInts = 0;

        forAllConstIter(Map<vector>, intersections, pI)
        {
            tP[nInts++] = pI();
        }

        if (debug > 3 || output)
        {
            if (meshOps::checkPointNearness(tP, 1e-20))
            {
                writeVTK(Foam::name(newCellIndex),newCellIndex,3,false,true);
                writeVTK(Foam::name(oldCellIndex),oldCellIndex,3,true,true);
            }

            meshOps::writeVTK
            (
                (*this),
                "ccSet_"
              + Foam::name(newCellIndex)
              + '<' + Foam::name(oldCellIndex) + '>',
                tP.size(),
                tP.size(),
                tP.size(),
                tP
            );
        }

        return true;
    }

    // Check for point-segment intersections
    bool pointIntersections = false;

    forAll(fromCellPoints, pointI)
    {
        label fromPoint = fromCellPoints[pointI];

        if (commonPoints.found(fromPoint))
        {
            continue;
        }

        const point& checkPoint = polyMesh::points()[fromPoint];

        forAll(toCellEdges, edgeI)
        {
            const edge edgeToCheck = toCellEdges[edgeI];

            if
            (
                meshOps::pointSegmentIntersection
                (
                    linePointRef
                    (
                        oldPoints_[edgeToCheck.start()],
                        oldPoints_[edgeToCheck.end()]
                    ),
                    checkPoint,
                    matchTol
                )
            )
            {
                commonPoints.insert
                (
                    fromPoint,
                    labelList(edgeToCheck)
                );

                intersections.set(++nInts, checkPoint);

                pointIntersections = true;
            }
        }
    }

    forAll(toCellPoints, pointI)
    {
        label toPoint = toCellPoints[pointI];

        if (commonPoints.found(toPoint))
        {
            continue;
        }

        const point& checkPoint = oldPoints_[toPoint];

        forAll(fromCellEdges, edgeI)
        {
            const edge edgeToCheck = fromCellEdges[edgeI];

            if
            (
                meshOps::pointSegmentIntersection
                (
                    linePointRef
                    (
                        polyMesh::points()[edgeToCheck.start()],
                        polyMesh::points()[edgeToCheck.end()]
                    ),
                    checkPoint,
                    matchTol
                )
            )
            {
                commonPoints.insert
                (
                    toPoint,
                    labelList(edgeToCheck)
                );

                intersections.set(++nInts, checkPoint);

                pointIntersections = true;
            }
        }
    }

    if (pointIntersections && debug > 3)
    {
        Info << "Point Intersections exist: " << nl
             << " newCellIndex: " << newCellIndex
             << " oldCellIndex: " << oldCellIndex
             << endl;
    }

    if (twoDMesh_)
    {
        // Check if edge mid-points are clearly within the cell.
        // If so, add edge points as 'intersections'.
        forAll(fromCellEdges, edgeI)
        {
            const edge edgeToCheck = fromCellEdges[edgeI];

            if
            (
                commonPoints.found(edgeToCheck.start()) &&
                commonPoints.found(edgeToCheck.end())
            )
            {
                continue;
            }

            vector edgeVec =
            (
                polyMesh::points()[edgeToCheck.start()] -
                polyMesh::points()[edgeToCheck.end()]
            );

            edgeVec /= mag(edgeVec) + VSMALL;

            if (mag(edgeVec & planeNormal) < 0.5)
            {
                continue;
            }

            vector checkPoint =
            (
                0.5 *
                (
                    polyMesh::points()[edgeToCheck.start()] +
                    polyMesh::points()[edgeToCheck.end()]
                )
            );

            if
            (
                meshOps::pointInCell
                (
                    newCellIndex,
                    toCell,
                    faces_,
                    owner_,
                    oldPoints_,
                    checkPoint
                )
            )
            {
                intersections.set
                (
                    ++nInts,
                    polyMesh::points()[edgeToCheck.start()]
                );

                intersections.set
                (
                    ++nInts,
                    polyMesh::points()[edgeToCheck.end()]
                );
            }
        }

        forAll(toCellEdges, edgeI)
        {
            const edge edgeToCheck = toCellEdges[edgeI];

            if
            (
                commonPoints.found(edgeToCheck.start()) &&
                commonPoints.found(edgeToCheck.end())
            )
            {
                continue;
            }

            vector edgeVec =
            (
                oldPoints_[edgeToCheck.start()] -
                oldPoints_[edgeToCheck.end()]
            );

            edgeVec /= mag(edgeVec) + VSMALL;

            if (mag(edgeVec & planeNormal) < 0.5)
            {
                continue;
            }

            vector checkPoint =
            (
                0.5 *
                (
                    oldPoints_[edgeToCheck.start()] +
                    oldPoints_[edgeToCheck.end()]
                )
            );

            if
            (
                meshOps::pointInCell
                (
                    oldCellIndex,
                    fromCell,
                    polyMesh::faces(),
                    polyMesh::faceOwner(),
                    polyMesh::points(),
                    checkPoint
                )
            )
            {
                intersections.set
                (
                    ++nInts,
                    oldPoints_[edgeToCheck.start()]
                );

                intersections.set
                (
                    ++nInts,
                    oldPoints_[edgeToCheck.end()]
                );
            }
        }
    }
    else
    {
        // Check whether any old points are within
        // the new cell. Count these as 'intersections'.
        forAll(fromCellPoints, pointI)
        {
            label fromPoint = fromCellPoints[pointI];

            if (commonPoints.found(fromPoint))
            {
                // Only skip for shared-points.
                // If the point-position was modified
                // due to a collapse, then this point
                // could be inside the new cell.
                if (commonPoints[fromPoint].empty())
                {
                    continue;
                }
            }

            const point& checkPoint = polyMesh::points()[fromPoint];

            if
            (
                meshOps::pointInCell
                (
                    newCellIndex,
                    toCell,
                    faces_,
                    owner_,
                    oldPoints_,
                    checkPoint
                )
            )
            {
                intersections.set(++nInts, checkPoint);
            }
        }

        // Check whether any new points are within
        // the old cell. Count these as 'intersections'.
        forAll(toCellPoints, pointI)
        {
            label toPoint = toCellPoints[pointI];

            if (commonPoints.found(toPoint))
            {
                continue;
            }

            const point& checkPoint = oldPoints_[toPoint];

            if
            (
                meshOps::pointInCell
                (
                    oldCellIndex,
                    fromCell,
                    polyMesh::faces(),
                    polyMesh::faceOwner(),
                    polyMesh::points(),
                    checkPoint
                )
            )
            {
                intersections.set(++nInts, checkPoint);
            }
        }
    }

    bool foundIntersection = false, edgeIntersections = false;

    // Loop through edges from each cell, and check whether they intersect.
    List<Pair<edge> > FeToTe, TeToFe;

    // Define edge-vectors in 2D
    vector fromVec(vector::zero), toVec(vector::zero);

    forAll(fromCellEdges, edgeI)
    {
        // For 2D meshes, only select edges on wedge/empty planes
        if (twoDMesh_)
        {
            fromVec =
            (
                polyMesh::points()[fromCellEdges[edgeI].start()] -
                polyMesh::points()[fromCellEdges[edgeI].end()]
            );

            fromVec /= mag(fromVec) + VSMALL;

            if (mag(fromVec & planeNormal) > 0.5)
            {
                continue;
            }
        }

        forAll(toCellEdges, edgeJ)
        {
            // For 2D meshes, only select edges on wedge/empty planes
            if (twoDMesh_)
            {
                toVec =
                (
                    oldPoints_[toCellEdges[edgeJ].start()] -
                    oldPoints_[toCellEdges[edgeJ].end()]
                );

                toVec /= mag(toVec) + VSMALL;

                if (mag(toVec & planeNormal) > 0.5)
                {
                    continue;
                }
            }

            // Form an edge-pair
            Pair<edge> edgePair(fromCellEdges[edgeI], toCellEdges[edgeJ]);

            bool disableCheck = false;

            // Check edges topologically
            if (edgePair.first() == edgePair.second())
            {
                const edge& checkEdge = edgePair.first();

                // Check if points were modified by a collapse.
                // If both were modified, continue with check.
                if
                (
                    !modPoints_.found(checkEdge.start()) &&
                    !modPoints_.found(checkEdge.end())
                )
                {
                    disableCheck = true;
                }

                // Skip shared points
                if (commonPoints.found(checkEdge.start()))
                {
                    if (commonPoints[checkEdge.start()].empty())
                    {
                        disableCheck = true;
                    }
                }

                if (commonPoints.found(checkEdge.end()))
                {
                    if (commonPoints[checkEdge.end()].empty())
                    {
                        disableCheck = true;
                    }
                }
            }
            else
            {
                // Check for common vertices
                label cV = edgePair.first().commonVertex(edgePair.second());

                if (cV > -1)
                {
                    // If this point was modified by a collapse
                    // to an edge mid-point, it can't be a common point.
                    // So, allow the check to continue.
                    if (!modPoints_.found(cV))
                    {
                        disableCheck = true;
                    }

                    // Skip shared points
                    if (commonPoints.found(cV))
                    {
                        if (commonPoints[cV].empty())
                        {
                            disableCheck = true;
                        }
                    }
                }
            }

            if (disableCheck)
            {
                continue;
            }

            // Deal with edge-bisection / point-on-edge cases
            bool foundPointOnEdge = false;

            forAll(edgePair, indexI)
            {
                const edge thisEdge = edgePair[indexI];

                const edge otherEdge =
                (
                    (thisEdge == edgePair.first()) ?
                    edgePair.second() : edgePair.first()
                );

                forAll(otherEdge, pointI)
                {
                    label pIndex = otherEdge[pointI];

                    if (commonPoints.found(pIndex))
                    {
                        // Fetch masterObjects
                        const labelList& mObj = commonPoints[pIndex];

                        // Skip shared-points.
                        if (mObj.size())
                        {
                            // Check if the old edge
                            // contains all master points
                            bool allMaster = true;

                            forAll(mObj, pointJ)
                            {
                                if (findIndex(thisEdge, mObj[pointJ]) == -1)
                                {
                                    allMaster = false;
                                    break;
                                }
                            }

                            if (allMaster)
                            {
                                foundPointOnEdge = true;
                            }
                        }
                    }

                    if (foundPointOnEdge)
                    {
                        break;
                    }
                }

                if (foundPointOnEdge)
                {
                    break;
                }
            }

            if (foundPointOnEdge)
            {
                continue;
            }

            foundIntersection = false;

            foundIntersection =
            (
                meshOps::segmentSegmentIntersection
                (
                    linePointRef
                    (
                        polyMesh::points()[edgePair.first().start()],
                        polyMesh::points()[edgePair.first().end()]
                    ),
                    linePointRef
                    (
                        oldPoints_[edgePair.second().start()],
                        oldPoints_[edgePair.second().end()]
                    ),
                    matchTol,
                    intPoint
                )
            );

            if (foundIntersection)
            {
                FeToTe.setSize
                (
                    FeToTe.size() + 1,
                    edgePair
                );

                TeToFe.setSize
                (
                    TeToFe.size() + 1,
                    edgePair.reversePair()
                );

                intersections.set(++nInts, intPoint);

                // Note for later that edge-intersections exist.
                edgeIntersections = true;
            }
        }
    }

    // If this is a 2D mesh, we're done.
    if (twoDMesh_)
    {
        // Does not intersect.
        if (nInts < 6)
        {
            return false;
        }

        // Copy intersections
        tP.setSize(nInts, vector::zero);

        nInts = 0;

        forAllConstIter(Map<vector>, intersections, pI)
        {
            tP[nInts++] = pI();
        }

        if (debug > 3 || output)
        {
            if (meshOps::checkPointNearness(tP, 1e-20))
            {
                writeVTK(Foam::name(newCellIndex),newCellIndex,3,false,true);
                writeVTK(Foam::name(oldCellIndex),oldCellIndex,3,true,true);
            }

            meshOps::writeVTK
            (
                (*this),
                "ccSet_"
              + Foam::name(newCellIndex)
              + '<' + Foam::name(oldCellIndex) + '>',
                tP.size(),
                tP.size(),
                tP.size(),
                tP
            );
        }

        // Found a convex set of points
        return true;
    }

    if (edgeIntersections && debug > 1)
    {
        Info << "Edge Intersections exist: " << nl
             << " newCellIndex: " << newCellIndex
             << " oldCellIndex: " << oldCellIndex
             << endl;
    }

    // Loop through all old edges, and find possible
    // intersections with faces of the new cell.
    forAll(fromCellEdges, edgeI)
    {
        const edge edgeToCheck(fromCellEdges[edgeI]);

        if
        (
            commonPoints.found(edgeToCheck.start()) &&
            commonPoints.found(edgeToCheck.end())
        )
        {
            // Are both points only shared?
            if
            (
                commonPoints[edgeToCheck.start()].empty() &&
                commonPoints[edgeToCheck.end()].empty()
            )
            {
                continue;
            }
        }

        forAll(toCell, faceI)
        {
            const triFace faceToCheck(faces_[toCell[faceI]]);

            // Avoid point-edge / edge-edge intersections, if any.
            if (edgeIntersections)
            {
                // Is edgeToCheck in the list?
                bool foundEdge = false;
                const edgeList fEdges = faceToCheck.edges();

                forAll(FeToTe, indexI)
                {
                    if (FeToTe[indexI].first() == edgeToCheck)
                    {
                        // Check whether the intersecting edge
                        // exists on this face.
                        forAll(fEdges, edgeJ)
                        {
                            if (fEdges[edgeJ] == FeToTe[indexI].second())
                            {
                                foundEdge = true;
                                break;
                            }
                        }

                        if (foundEdge)
                        {
                            break;
                        }
                    }
                }

                if (foundEdge)
                {
                    continue;
                }
            }

            bool foundCommon = false;

            forAllConstIter(Map<labelList>, commonPoints, pIter)
            {
                if (findIndex(faceToCheck, pIter.key()) > -1)
                {
                    // Avoid shared points, since this implies that
                    // the edge intersects at a face point
                    if (edgeToCheck[0] == pIter.key())
                    {
                        if (pIter().empty())
                        {
                            foundCommon = true;
                            break;
                        }
                    }

                    if (edgeToCheck[1] == pIter.key())
                    {
                        if (pIter().empty())
                        {
                            foundCommon = true;
                            break;
                        }
                    }

                    // Also check for bisection points.
                    // This accounts for successive bisections.
                    if
                    (
                        (findIndex(pIter(), edgeToCheck[0]) > -1) &&
                        (findIndex(pIter(), edgeToCheck[1]) > -1)
                    )
                    {
                        foundCommon = true;
                        break;
                    }
                }
            }

            if (foundCommon)
            {
                continue;
            }

            foundIntersection = false;

            foundIntersection =
            (
                meshOps::segmentTriFaceIntersection
                (
                    triPointRef
                    (
                        oldPoints_[faceToCheck[0]],
                        oldPoints_[faceToCheck[1]],
                        oldPoints_[faceToCheck[2]]
                    ),
                    linePointRef
                    (
                        polyMesh::points()[edgeToCheck.start()],
                        polyMesh::points()[edgeToCheck.end()]
                    ),
                    matchTol,
                    intPoint
                )
            );

            if (foundIntersection)
            {
                // Add to the list.
                intersections.set(++nInts, intPoint);
            }
        }
    }

    // Loop through all new edges, and find possible
    // intersections with faces of the old cell.
    forAll(toCellEdges, edgeI)
    {
        const edge edgeToCheck(toCellEdges[edgeI]);

        if
        (
            commonPoints.found(edgeToCheck.start()) &&
            commonPoints.found(edgeToCheck.end())
        )
        {
            // Are both points only shared?
            if
            (
                commonPoints[edgeToCheck.start()].empty() &&
                commonPoints[edgeToCheck.end()].empty()
            )
            {
                continue;
            }
        }

        forAll(fromCell, faceI)
        {
            const triFace faceToCheck(polyMesh::faces()[fromCell[faceI]]);

            // Avoid point-edge / edge-edge intersections, if any.
            if (edgeIntersections)
            {
                // Is edgeToCheck in the list?
                bool foundEdge = false;
                const edgeList fEdges = faceToCheck.edges();

                forAll(TeToFe, indexI)
                {
                    if (TeToFe[indexI].first() == edgeToCheck)
                    {
                        // Check whether the intersecting edge
                        // exists on this face.
                        forAll(fEdges, edgeJ)
                        {
                            if (fEdges[edgeJ] == TeToFe[indexI].second())
                            {
                                foundEdge = true;
                                break;
                            }
                        }

                        if (foundEdge)
                        {
                            break;
                        }
                    }
                }

                if (foundEdge)
                {
                    continue;
                }
            }

            bool foundCommon = false;

            forAllConstIter(Map<labelList>, commonPoints, pIter)
            {
                if (findIndex(faceToCheck, pIter.key()) > -1)
                {
                    // Avoid shared points, since this implies that
                    // the edge intersects at a face point
                    if (edgeToCheck[0] == pIter.key())
                    {
                        if (pIter().empty())
                        {
                            foundCommon = true;
                            break;
                        }
                    }

                    if (edgeToCheck[1] == pIter.key())
                    {
                        if (pIter().empty())
                        {
                            foundCommon = true;
                            break;
                        }
                    }

                    // Also check for bisection points
                    if
                    (
                        (findIndex(pIter(), edgeToCheck[0]) > -1) &&
                        (findIndex(pIter(), edgeToCheck[1]) > -1)
                    )
                    {
                        foundCommon = true;
                        break;
                    }

                    // Check for point-on-edge cases
                    bool foundPointOnEdge = false;

                    if (pIter().size())
                    {
                        bool allMaster = true;

                        const labelList& mObj = pIter();

                        forAll(mObj, pointJ)
                        {
                            if (findIndex(faceToCheck, mObj[pointJ]) == -1)
                            {
                                allMaster = false;
                                break;
                            }
                        }

                        if (allMaster)
                        {
                            foundPointOnEdge = true;
                        }
                    }

                    if (foundPointOnEdge)
                    {
                        foundCommon = true;
                        break;
                    }
                }
            }

            if (foundCommon)
            {
                continue;
            }

            foundIntersection = false;

            foundIntersection =
            (
                meshOps::segmentTriFaceIntersection
                (
                    triPointRef
                    (
                        polyMesh::points()[faceToCheck[0]],
                        polyMesh::points()[faceToCheck[1]],
                        polyMesh::points()[faceToCheck[2]]
                    ),
                    linePointRef
                    (
                        oldPoints_[edgeToCheck.start()],
                        oldPoints_[edgeToCheck.end()]
                    ),
                    matchTol,
                    intPoint
                )
            );

            if (foundIntersection)
            {
                // Add to the list.
                intersections.set(++nInts, intPoint);
            }
        }
    }

    if (nInts < 4)
    {
        // Does not intersect.
        return false;
    }

    // Copy intersections
    tP.setSize(nInts, vector::zero);

    nInts = 0;

    forAllConstIter(Map<vector>, intersections, pI)
    {
        tP[nInts++] = pI();
    }

    if (debug > 3 || output)
    {
        if (meshOps::checkPointNearness(tP, 1e-20))
        {
            writeVTK(Foam::name(newCellIndex),newCellIndex,3,false,true);
            writeVTK(Foam::name(oldCellIndex),oldCellIndex,3,true,true);
        }

        meshOps::writeVTK
        (
            (*this),
            "ccSet_"
          + Foam::name(newCellIndex)
          + '<' + Foam::name(oldCellIndex) + '>',
            tP.size(),
            tP.size(),
            tP.size(),
            tP
        );
    }

    // Found a convex set of points.
    return true;
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
