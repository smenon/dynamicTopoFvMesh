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
#include "objectMap.H"
#include "StaticHashTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Compute mapping weights for modified entities
void dynamicTopoFvMesh::computeMapping
(
    const scalar matchTol,
    const label faceStart,
    const label faceSize,
    const label cellStart,
    const label cellSize
)
{
    // Compute cell mapping
    for (label cellI = cellStart; cellI < (cellStart + cellSize); cellI++)
    {
        label cIndex = cellsFromCells_[cellI].index();

        // Obtain weighting factors for this cell.
        computeCellWeights
        (
            cIndex,
            cellParents_[cIndex],
            matchTol,
            cellsFromCells_[cellI].masterObjects(),
            cellWeights_[cellI],
            cellCentres_[cellI]
        );
    }

    // Compute face mapping
    for (label faceI = faceStart; faceI < (faceStart + faceSize); faceI++)
    {
        label fIndex = facesFromFaces_[faceI].index();

        // Skip mapping for internal faces.
        if (whichPatch(fIndex) == -1)
        {
            continue;
        }

        // Obtain weighting factors for this face.
        computeFaceWeights
        (
            fIndex,
            faceParents_[fIndex],
            matchTol,
            facesFromFaces_[faceI].masterObjects(),
            faceWeights_[faceI],
            faceCentres_[faceI]
        );
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
    scalar& matchTol = *(static_cast<scalar*>(thread->operator()(0)));
    label& faceStart = *(static_cast<label*>(thread->operator()(1)));
    label& faceSize  = *(static_cast<label*>(thread->operator()(2)));
    label& cellStart = *(static_cast<label*>(thread->operator()(3)));
    label& cellSize  = *(static_cast<label*>(thread->operator()(4)));

    // Now calculate addressing
    mesh.computeMapping
    (
        matchTol,
        faceStart, faceSize,
        cellStart, cellSize
    );

    if (thread->slave())
    {
        thread->sendSignal(meshHandler::STOP);
    }
}


// Routine to invoke threaded mapping
void dynamicTopoFvMesh::threadedMapping(scalar matchTol)
{
    label nThreads = threader_->getNumThreads();

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
        hdl[i].setSize(5);

        // Set match tolerance
        hdl[i].set(0, &matchTol);

        // Set the start/size indices
        hdl[i].set(1, &(tStarts[0][0]));
        hdl[i].set(2, &(tSizes[0][1]));
        hdl[i].set(3, &(tStarts[1][0]));
        hdl[i].set(4, &(tSizes[1][1]));
    }

    // Prior to multi-threaded operation,
    // force calculation of demand-driven data.
    polyMesh::cells();
    primitiveMesh::cellCells();
    primitiveMesh::faceAreas();
    primitiveMesh::cellCentres();
    primitiveMesh::faceCentres();

    const polyBoundaryMesh& boundary = boundaryMesh();

    forAll(boundary, patchI)
    {
        boundary[patchI].faceFaces();
    }

    // Execute threads in linear sequence
    executeThreads(identity(nThreads), hdl, &computeMappingThread);
}


// Obtain map weighting factors for a face
void dynamicTopoFvMesh::computeFaceWeights
(
    const label fIndex,
    const labelList& mapCandidates,
    const scalar mTol,
    labelList& Parents,
    scalarField& Weights,
    vectorField& Centres
) const
{
    if (Parents.size() || Weights.size() || Centres.size())
    {
        FatalErrorIn
        (
            "\n\n"
            "void dynamicTopoFvMesh::computeFaceWeights\n"
            "(\n"
            "    const label fIndex,\n"
            "    const labelList& mapCandidates,\n"
            "    labelList& parents,\n"
            "    scalarField& weights,\n"
            "    vectorField& centres\n"
            ") const\n"
        )
            << " Addressing has already been calculated." << nl
            << " Face: " << fIndex << nl
            << " mapCandidates: " << mapCandidates << nl
            << abort(FatalError);
    }

    // Determine which patch this is...
    label patchIndex = whichPatch(fIndex);

    scalar searchFactor = 1.0, matchTol = mTol;
    label nOldIntersects = -1, nIntersects = 0, nAttempts = 0;

    // Maintain a list of candidates and intersection points
    boolList oldIntersects, intersects;
    labelList oldCandidates, candidates;

    // Output option for the convex set algorithm
    bool output = false;

    // Fetch the area / centre of the cell
    vector fNormal =
    (
        meshOps::faceNormal
        (
            faces_[fIndex],
            oldPoints_
        )
    );

    scalar fArea = mag(fNormal);

    vector fCentre =
    (
        meshOps::faceCentre
        (
            faces_[fIndex],
            oldPoints_
        )
    );

    // Local variables
    labelList parents;
    scalarField weights;
    vectorField centres;

    while (nAttempts < 10)
    {
        // Reset counter first
        nIntersects = 0;

        // Obtain candidate parents for this face
        candidates =
        (
            faceParents
            (
                fIndex,
                fNormal,
                fCentre,
                searchFactor,
                mapCandidates
            )
        );

        // For empty patches, skip calculations
        if (patchIndex > -1)
        {
            if (boundaryMesh()[patchIndex].type() == "empty")
            {
                // Set sizes
                parents.setSize(1, candidates[0]);
                weights.setSize(1, fArea);
                centres.setSize(1, vector::zero);

                break;
            }
        }

        // Set sizes
        intersects.setSize(candidates.size());

        // Test for intersections
        forAll(candidates, indexI)
        {
            vectorField tP;
            intersects[indexI] = false;

            intersects[indexI] =
            (
                faceIntersection
                (
                    fIndex,
                    candidates[indexI],
                    matchTol,
                    tP
                )
            );

            if (intersects[indexI])
            {
                nIntersects++;
            }
        }

        if ((nIntersects == nOldIntersects) && (nIntersects != 0))
        {
            if (debug > 3)
            {
                Info << " Face: " << fIndex << nl
                     << " nCandidates: " << candidates.size() << nl
                     << " nOldCandidates: " << oldCandidates.size() << nl
                     << " nIntersects: " << nIntersects
                     << endl;

                // Specify output option as well
                output = true;
            }

            // Set sizes
            parents.setSize(nIntersects, -1);
            weights.setSize(nIntersects, 0.0);
            centres.setSize(nIntersects, vector::zero);

            // Reset counter
            nIntersects = 0;

            // Compute actual intersections
            forAll(oldIntersects, indexI)
            {
                if (oldIntersects[indexI])
                {
                    vectorField tP(0);

                    oldIntersects[indexI] =
                    (
                        faceIntersection
                        (
                            fIndex,
                            oldCandidates[indexI],
                            matchTol,
                            tP
                        )
                    );

                    // Skip false positives
                    if (oldIntersects[indexI])
                    {
                        parents[nIntersects] = oldCandidates[indexI];

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
                            parents[nIntersects],
                            tP,
                            refNorm,
                            weights[nIntersects],
                            centres[nIntersects],
                            output
                        );

                        nIntersects++;
                    }
                }
            }

            // Shorten to actual sizes
            parents.setSize(nIntersects);
            weights.setSize(nIntersects);
            centres.setSize(nIntersects);

            break;
        }
        else
        {
            nAttempts++;

            // Copy / reset parameters
            oldIntersects = intersects;
            oldCandidates = candidates;
            nOldIntersects = nIntersects;

            // Expand the search radius and try again.
            searchFactor *= 1.6;
        }
    }

    // Test weights for consistency
    if (mag(1.0 - sum(weights/fArea)) > 1e-10)
    {
        bool inconsistent = true;

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
                    // Normalize by sum of weights instead
                    fArea = sum(weights);
                    inconsistent = false;
                }
            }
        }

        if (inconsistent)
        {
            // Write out for post-processing
            label uIdx = 0;
            labelList unMatched(oldCandidates.size() - parents.size(), -1);

            forAll(oldCandidates, cI)
            {
                if (findIndex(parents, oldCandidates[cI]) == -1)
                {
                    unMatched[uIdx++] = oldCandidates[cI];
                }
            }

            writeVTK("nFace_" + Foam::name(fIndex), fIndex, 2, false, true);
            writeVTK("oFace_" + Foam::name(fIndex), oldCandidates, 2, true, true);
            writeVTK("mFace_" + Foam::name(fIndex), parents, 2, true, true);
            writeVTK("uFace_" + Foam::name(fIndex), unMatched, 2, true, true);

            // Write out weights
            forAll(parents, indexI)
            {
                Info << parents[indexI] << ": "
                     << setprecision(16)
                     << weights[indexI]
                     << endl;
            }

            FatalErrorIn
            (
                "\n\n"
                "void dynamicTopoFvMesh::computeFaceWeights\n"
                "(\n"
                "    const label fIndex,\n"
                "    const labelList& mapCandidates,\n"
                "    labelList& parents,\n"
                "    scalarField& weights,\n"
                "    vectorField& centres\n"
                ") const\n"
            )
                << "Encountered non-conservative weighting factors." << nl
                << " Face: " << fIndex << ":: " << faces_[fIndex] << nl
                << " Patch: " << boundaryMesh()[patchIndex].name() << nl
                << " mapCandidates: " << mapCandidates << nl
                << " nCandidates: " << candidates.size() << nl
                << " nOldCandidates: " << oldCandidates.size() << nl
                << " nIntersects: " << nIntersects << nl
                << " nOldIntersects: " << nOldIntersects << nl
                << " nParents: " << parents.size() << nl
                << " nAttempts: " << nAttempts << nl
                << " nFaces: " << nFaces_ << nl
                << " nOldFaces: " << nOldFaces_ << nl
                << " nInternalFaces: " << nInternalFaces_ << nl
                << " nOldInternalFaces: " << nOldInternalFaces_ << nl
                << setprecision(16)
                << " Face area: " << fArea << nl
                << " Sum(Weights): " << sum(weights) << nl
                << " Error: " << (fArea - sum(weights)) << nl
                << " Norm Sum(Weights): " << sum(weights/fArea) << nl
                << " Norm Error: " << mag(1.0 - sum(weights/fArea))
                << abort(FatalError);
        }
    }
    else
    if (debug > 2)
    {
        Info << " Face: " << fIndex << nl
             << setprecision(16)
             << " Face area: " << fArea << nl
             << " Sum(Weights): " << sum(weights) << nl
             << " Error: " << (fArea - sum(weights)) << nl
             << " Norm Sum(Weights): " << sum(weights/fArea) << nl
             << " Norm Error: " << mag(1.0 - sum(weights/fArea))
             << endl;
    }

    // Return normalized weights
    weights /= fArea;

    // Set inputs
    Parents.setSize(parents.size(), -1);
    Weights.setSize(weights.size(), 0.0);
    Centres.setSize(centres.size(), vector::zero);

    // Copy fields
    forAll(Parents, indexI)
    {
        Parents[indexI] = parents[indexI];
        Weights[indexI] = weights[indexI];
        Centres[indexI] = centres[indexI];
    }
}


// Obtain map weighting factors for a cell
void dynamicTopoFvMesh::computeCellWeights
(
    const label cIndex,
    const labelList& mapCandidates,
    const scalar mTol,
    labelList& Parents,
    scalarField& Weights,
    vectorField& Centres
) const
{
    if (Parents.size() || Weights.size() || Centres.size())
    {
        FatalErrorIn
        (
            "\n\n"
            "void dynamicTopoFvMesh::computeCellWeights\n"
            "(\n"
            "    const label cIndex,\n"
            "    const labelList& mapCandidates,\n"
            "    labelList& parents,\n"
            "    scalarField& weights,\n"
            "    vectorField& centres\n"
            ") const\n"
        )
            << " Addressing has already been calculated." << nl
            << " Cell: " << cIndex << nl
            << " mapCandidates: " << mapCandidates << nl
            << abort(FatalError);
    }

    label nAttempts = 0;
    scalar searchFactor = 1.0, matchTol = mTol;
    label nOldIntersects = -1, nIntersects = 0, realIntersects = 0;

    // Maintain a list of candidates and intersection points
    boolList oldIntersects, intersects;
    labelList oldCandidates, candidates;

    // Output option for the convex set algorithm
    bool output = false;

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

    // Local variables
    labelList parents;
    scalarField weights;
    vectorField centres;

    while (nAttempts < 10)
    {
        // Reset counter first
        nIntersects = 0;

        // Obtain candidate parents for this cell
        candidates =
        (
            cellParents
            (
                cIndex,
                cellCentre,
                searchFactor,
                mapCandidates
            )
        );

        // Set sizes and reset
        intersects.setSize(candidates.size());

        // Test for intersections
        forAll(candidates, indexI)
        {
            intersects[indexI] = false;

            intersects[indexI] =
            (
                testCellIntersection
                (
                    cIndex,
                    candidates[indexI]
                )
            );

            if (intersects[indexI])
            {
                nIntersects++;
            }
        }

        label nInnerAttempts = 0;
        bool attainedAccuracy = false;

        if ((nIntersects == nOldIntersects) && (nIntersects != 0))
        {
            if (debug > 3)
            {
                Info << " Cell: " << cIndex << nl
                     << " nCandidates: " << candidates.size() << nl
                     << " nOldCandidates: " << oldCandidates.size() << nl
                     << " nIntersects: " << nIntersects
                     << endl;

                // Specify output option as well
                output = true;
            }

            // Set sizes
            parents.setSize(nIntersects, -1);
            weights.setSize(nIntersects, 0.0);
            centres.setSize(nIntersects, vector::zero);

            // Compute actual intersections
            while (nInnerAttempts < 5)
            {
                // Reset counter
                realIntersects = 0;
                scalar sumVols = 0.0;

                forAll(oldIntersects, indexI)
                {
                    if (oldIntersects[indexI])
                    {
                        vectorField tP(0);

                        bool realIntersect =
                        (
                            cellIntersection
                            (
                                cIndex,
                                oldCandidates[indexI],
                                matchTol,
                                tP
                            )
                        );

                        // Skip false positives
                        if (realIntersect)
                        {
                            parents[realIntersects] = oldCandidates[indexI];

                            // Compute weights
                            meshOps::convexSetVolume
                            (
                                cIndex,
                                parents[realIntersects],
                                tP,
                                weights[realIntersects],
                                centres[realIntersects],
                                output
                            );

                            // Accumulate volume
                            sumVols += weights[realIntersects];

                            realIntersects++;
                        }
                    }
                }

                scalar mismatch = (1.0 - (sumVols/cellVolume));

                // Check for consistency
                if (mag(mismatch) > 1e-10)
                {
                    if (mismatch < 0.0)
                    {
                        // Reduce geometric tolerance, and try again.
                        matchTol *= 0.1;
                        nInnerAttempts++;
                    }
                    else
                    {
                        // Expand search radius
                        break;
                    }
                }
                else
                {
                    // Attained sufficient accuracy
                    attainedAccuracy = true;
                    break;
                }
            }
        }

        // Shorten to actual sizes
        parents.setSize(realIntersects);
        weights.setSize(realIntersects);
        centres.setSize(realIntersects);

        if (attainedAccuracy)
        {
            break;
        }
        else
        {
            // Reset match tolerance, if necessary
            if (nInnerAttempts)
            {
                matchTol *= Foam::pow(10, nInnerAttempts);
                nInnerAttempts = 0;
            }

            nAttempts++;

            // Copy parameters
            oldIntersects = intersects;
            oldCandidates = candidates;
            nOldIntersects = nIntersects;

            // Expand the search radius and try again.
            searchFactor *= 1.4;
        }
    }

    // Test weights for consistency
    if (mag(1.0 - sum(weights/cellVolume)) > 1e-10)
    {
        // Inconsistent weights. Check whether any edges
        // lie on boundary patches. These cells can have
        // relaxed weights to account for mild convexity.
        bool inconsistent = true;

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
                        inconsistent = false;
                        break;
                    }
                }

                if (!inconsistent)
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
                        inconsistent = false;
                        break;
                    }
                }

                if (!inconsistent)
                {
                    break;
                }
            }
        }

        if (!inconsistent && (sum(weights) > VSMALL))
        {
            // Normalize by sum of weights instead
            cellVolume = sum(weights);
        }
        else
        {
            // Write out for post-processing
            label uIdx = 0;
            labelList unMatched(oldCandidates.size() - parents.size(), -1);

            forAll(oldCandidates, cI)
            {
                if (findIndex(parents, oldCandidates[cI]) == -1)
                {
                    unMatched[uIdx++] = oldCandidates[cI];
                }
            }

            writeVTK("nCell_" + Foam::name(cIndex), cIndex, 3, false, true);
            writeVTK("oCell_" + Foam::name(cIndex), oldCandidates, 3, true, true);
            writeVTK("mCell_" + Foam::name(cIndex), parents, 3, true, true);
            writeVTK("uCell_" + Foam::name(cIndex), unMatched, 3, true, true);

            // Write out weights
            forAll(parents, indexI)
            {
                Info << parents[indexI] << ": "
                     << setprecision(16)
                     << weights[indexI]
                     << endl;
            }

            // Write out intersection points
            forAll(oldCandidates, indexI)
            {
                vectorField tP(0);

                cellIntersection
                (
                    cIndex,
                    oldCandidates[indexI],
                    matchTol,
                    tP
                );

                if (tP.size() >= 4)
                {
                    // Write out intersection points to VTK
                    meshOps::writeVTK
                    (
                        (*this),
                        "cvxSet_"
                      + Foam::name(cIndex)
                      + '<' + Foam::name(oldCandidates[indexI]) + '>',
                        tP.size(),
                        tP.size(),
                        tP.size(),
                        tP
                    );

                    // Write out convex set info to screen
                    scalar dummyVolume = 0.0;
                    vector dummyCentre = vector::zero;

                    meshOps::convexSetVolume
                    (
                        cIndex,
                        oldCandidates[indexI],
                        tP,
                        dummyVolume,
                        dummyCentre,
                        true
                    );
                }
            }

            FatalErrorIn
            (
                "\n\n"
                "void dynamicTopoFvMesh::computeCellWeights\n"
                "(\n"
                "    const label cIndex,\n"
                "    const labelList& mapCandidates,\n"
                "    labelList& parents,\n"
                "    scalarField& weights,\n"
                "    vectorField& centres\n"
                ") const\n"
            )
                << "Encountered non-conservative weighting factors." << nl
                << " Cell: " << cIndex << nl
                << " mapCandidates: " << mapCandidates << nl
                << " nCandidates: " << candidates.size() << nl
                << " nOldCandidates: " << oldCandidates.size() << nl
                << " nIntersects: " << nIntersects << nl
                << " nOldIntersects: " << nOldIntersects << nl
                << " nRealIntersects: " << realIntersects << nl
                << " nParents: " << parents.size() << nl
                << " nAttempts: " << nAttempts << nl
                << " matchTolerance: " << matchTol << nl
                << " nCells: " << nCells_ << nl
                << " nOldCells: " << nOldCells_ << nl
                << setprecision(16)
                << " Cell volume: " << cellVolume << nl
                << " Sum(Weights): " << sum(weights) << nl
                << " Error: " << (cellVolume - sum(weights)) << nl
                << " Norm Sum(Weights): " << sum(weights/cellVolume) << nl
                << " Norm Error: " << mag(1.0 - sum(weights/cellVolume))
                << abort(FatalError);
        }
    }
    else
    if (debug > 2)
    {
        Info << " Cell: " << cIndex << nl
             << setprecision(16)
             << " Cell volume: " << cellVolume << nl
             << " Sum(Weights): " << sum(weights) << nl
             << " Error: " << (cellVolume - sum(weights)) << nl
             << " Norm Sum(Weights): " << sum(weights/cellVolume) << nl
             << " Norm Error: " << mag(1.0 - sum(weights/cellVolume))
             << endl;
    }

    // Return normalized weights
    weights /= cellVolume;

    // Set inputs
    Parents.setSize(parents.size(), -1);
    Weights.setSize(weights.size(), 0.0);
    Centres.setSize(centres.size(), vector::zero);

    // Copy fields
    forAll(Parents, indexI)
    {
        Parents[indexI] = parents[indexI];
        Weights[indexI] = weights[indexI];
        Centres[indexI] = centres[indexI];
    }
}


// Test for intersection between cells in old/new meshes
//   - Uses the static separating axis test for polyhedra,
//     outlined in work by David Eberly
//     'Intersection of Convex Objects: The Method of Separating Axes'
//     http://www.geometrictools.com/
bool dynamicTopoFvMesh::testCellIntersection
(
    const label newCellIndex,
    const label oldCellIndex
) const
{
    // Direction vector
    vector dir(vector::zero);

    // Check indices first
    if (newCellIndex >= cells_.size() || newCellIndex < 0)
    {
        FatalErrorIn
        (
            "\n"
            "bool dynamicTopoFvMesh::testCellIntersection\n"
            "(\n"
            "    const label newCellIndex,\n"
            "    const label oldCellIndex\n"
            ") const"
        )
            << " Wrong newCellIndex: " << newCellIndex << nl
            << " nCells: " << nCells_
            << abort(FatalError);
    }

    if (oldCellIndex >= nOldCells_ || oldCellIndex < 0)
    {
        FatalErrorIn
        (
            "\n"
            "bool dynamicTopoFvMesh::testCellIntersection\n"
            "(\n"
            "    const label newCellIndex,\n"
            "    const label oldCellIndex\n"
            ") const"
        )
            << " Wrong oldCellIndex: " << oldCellIndex << nl
            << " nOldCells: " << nOldCells_
            << abort(FatalError);
    }

    // Fetch references for each mesh
    const faceList& fromFaces = polyMesh::faces();
    const labelList& fromOwner = polyMesh::faceOwner();
    const cell& fromCell = polyMesh::cells()[oldCellIndex];
    const edgeList fromCellEdges = fromCell.edges(polyMesh::faces());
    const labelList fromCellPoints = fromCell.labels(polyMesh::faces());

    const cell& toCell = cells_[newCellIndex];
    const edgeList toCellEdges = toCell.edges(faces_);
    const labelList toCellPoints = toCell.labels(faces_);

    // Test faces of oldCell for separation
    forAll(fromCell, faceI)
    {
        const label fIndex = fromCell[faceI];

        dir = meshOps::faceNormal(fromFaces[fIndex], polyMesh::points());

        // Reverse normal if necessary
        if (fromOwner[fIndex] != oldCellIndex)
        {
            dir *= -1.0;
        }

        if
        (
            meshOps::whichSide
            (
                toCellPoints,
                oldPoints_,
                dir,
                polyMesh::points()[fromFaces[fIndex][0]]
            ) > 0
        )
        {
            return false;
        }
    }

    // Test faces of newCell for separation
    forAll(toCell, faceI)
    {
        const label fIndex = toCell[faceI];

        dir = meshOps::faceNormal(faces_[fIndex], oldPoints_);

        // Reverse normal if necessary
        if (owner_[fIndex] != newCellIndex)
        {
            dir *= -1.0;
        }

        if
        (
            meshOps::whichSide
            (
                fromCellPoints,
                polyMesh::points(),
                dir,
                oldPoints_[faces_[fIndex][0]]
            ) > 0
        )
        {
            return false;
        }
    }

    // Test edges of both cells for separation
    forAll(fromCellEdges, edgeI)
    {
        const edge& fromEdge = fromCellEdges[edgeI];

        vector fromVec =
        (
            polyMesh::points()[fromEdge[1]]
          - polyMesh::points()[fromEdge[0]]
        );

        forAll(toCellEdges, edgeJ)
        {
            const edge& toEdge = toCellEdges[edgeJ];

            vector toVec =
            (
                oldPoints_[toEdge[1]] - oldPoints_[toEdge[0]]
            );

            dir = (fromVec ^ toVec);

            label firstSide =
            (
                meshOps::whichSide
                (
                    fromCellPoints,
                    polyMesh::points(),
                    dir,
                    polyMesh::points()[fromEdge[0]]
                )
            );

            if (firstSide == 0)
            {
                continue;
            }

            label secondSide =
            (
                meshOps::whichSide
                (
                    toCellPoints,
                    oldPoints_,
                    dir,
                    polyMesh::points()[fromEdge[0]]
                )
            );

            if (secondSide == 0)
            {
                continue;
            }

            if ((firstSide * secondSide) < 0)
            {
                return false;
            }
        }
    }

    return true;
}


// Return the intersection points between faces in old/new meshes
bool dynamicTopoFvMesh::faceIntersection
(
    const label newFaceIndex,
    const label oldFaceIndex,
    const scalar matchTol,
    vectorField& tP
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

    // Add all new points, if they resulted
    // from bisections of old face edges.
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

    forAll(fromFace, pointI)
    {
        label pIndex = findIndex(toFace, fromFace[pointI]);

        vector r = polyMesh::points()[fromFace[pointI]];

        if (pIndex == -1)
        {
            // Project this point on to the toFace plane.
            projections[pointI] = xf + ((r - xf) - ((r - xf) & nf)*nf);

            projPoints.insert(fromFace[pointI], pointI);
        }
        else
        {
            // If this point was modified by a collapse
            // to an edge mid-point, it can't be a common point.
            if (modPoints_.found(toFace[pIndex]))
            {
                // Project this point on to the toFace plane.
                projections[pointI] = xf + ((r - xf) - ((r - xf) & nf)*nf);

                projPoints.insert(fromFace[pointI], pointI);
            }
            else
            {
                commonPoints.insert(toFace[pIndex], labelList(0));

                projections[pointI] = r;

                intersections.set(++nInts, r);
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

        if (debug)
        {
            if (meshOps::checkPointNearness(tP, 1e-20))
            {
                writeVTK(Foam::name(newFaceIndex),newFaceIndex,2,false,true);
                writeVTK(Foam::name(oldFaceIndex),oldFaceIndex,2,true,true);

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
        }

        return true;
    }

    // Check whether any old projections are within
    // the new face. Count these as 'intersections'.
    forAll(fromFace, pointI)
    {
        if (commonPoints.found(fromFace[pointI]))
        {
            continue;
        }

        const point& checkPoint = projections[pointI];

        if (meshOps::pointInFace(toFace, oldPoints_, checkPoint))
        {
            intersections.set(++nInts, checkPoint);
        }
    }

    // Check whether and new points are within
    // projected old faces. Count these as 'intersections'.
    face ifFace(identity(fromFace.size()));

    forAll(toFace, pointI)
    {
        if (commonPoints.found(toFace[pointI]))
        {
            continue;
        }

        const point& checkPoint = oldPoints_[toFace[pointI]];

        if (meshOps::pointInFace(ifFace, projections, checkPoint))
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

            // Avoid common points
            label cv = toEdge.commonVertex(fromEdge);

            if (cv > -1 && !modPoints_.found(cv))
            {
                continue;
            }

            // Also check for bisection points

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

    // Copy intersections
    tP.setSize(nInts, vector::zero);

    nInts = 0;

    forAllConstIter(Map<vector>, intersections, pI)
    {
        tP[nInts++] = pI();
    }

    // Check for concurrent points.
    if (debug)
    {
        if (meshOps::checkPointNearness(tP, 1e-20))
        {
            writeVTK(Foam::name(newFaceIndex),newFaceIndex,2,false,true);
            writeVTK(Foam::name(oldFaceIndex),oldFaceIndex,2,true,true);

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
    vectorField& tP
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
        label pIndex = findIndex(toCellPoints, fromCellPoints[pointI]);

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

    // Add all new points as well, if they resulted
    // from bisections of old cell edges.
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

        if (debug > 3)
        {
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
        if (commonPoints.found(fromCellPoints[pointI]))
        {
            continue;
        }

        const point& checkPoint = polyMesh::points()[fromCellPoints[pointI]];

        forAll(toCellEdges, edgeI)
        {
            const edge& edgeToCheck = toCellEdges[edgeI];

            if
            (
                meshOps::pointSegmentIntersection
                (
                    edgeToCheck,
                    oldPoints_,
                    checkPoint
                )
            )
            {
                commonPoints.insert
                (
                    fromCellPoints[pointI],
                    labelList(edgeToCheck)
                );

                intersections.set(++nInts, checkPoint);

                pointIntersections = true;
            }
        }
    }

    forAll(toCellPoints, pointI)
    {
        if (commonPoints.found(toCellPoints[pointI]))
        {
            continue;
        }

        const point& checkPoint = oldPoints_[toCellPoints[pointI]];

        forAll(fromCellEdges, edgeI)
        {
            const edge& edgeToCheck = fromCellEdges[edgeI];

            if
            (
                meshOps::pointSegmentIntersection
                (
                    edgeToCheck,
                    polyMesh::points(),
                    checkPoint
                )
            )
            {
                commonPoints.insert
                (
                    toCellPoints[pointI],
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
            const edge& edgeToCheck = fromCellEdges[edgeI];

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
            const edge& edgeToCheck = toCellEdges[edgeI];

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
            if (commonPoints.found(fromCellPoints[pointI]))
            {
                continue;
            }

            const point& checkPoint =
            (
                polyMesh::points()[fromCellPoints[pointI]]
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
                intersections.set(++nInts, checkPoint);
            }
        }

        // Check whether any new points are within
        // the old cell. Count these as 'intersections'.
        forAll(toCellPoints, pointI)
        {
            if (commonPoints.found(toCellPoints[pointI]))
            {
                continue;
            }

            const point& checkPoint =
            (
                oldPoints_[toCellPoints[pointI]]
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

                if
                (
                    commonPoints.found(checkEdge.start()) ||
                    commonPoints.found(checkEdge.end())
                )
                {
                    disableCheck = true;
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

                    if (commonPoints.found(cV))
                    {
                        disableCheck = true;
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
                    edgePair.first(),
                    edgePair.second(),
                    polyMesh::points(),
                    oldPoints_,
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
        const edge& edgeToCheck = fromCellEdges[edgeI];

        if
        (
            commonPoints.found(edgeToCheck.start()) &&
            commonPoints.found(edgeToCheck.end())
        )
        {
            continue;
        }

        forAll(toCell, faceI)
        {
            const face& faceToCheck = faces_[toCell[faceI]];

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
                if (faceToCheck.which(pIter.key()) > -1)
                {
                    // Avoid common points, since this implies that
                    // the edge intersects at a face point
                    if
                    (
                        (edgeToCheck[0] == pIter.key()) ||
                        (edgeToCheck[1] == pIter.key())
                    )
                    {
                        foundCommon = true;
                        break;
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
                }
            }

            if (foundCommon)
            {
                continue;
            }

            foundIntersection = false;

            foundIntersection =
            (
                meshOps::segmentFaceIntersection
                (
                    edgeToCheck,
                    faceToCheck,
                    polyMesh::points(),
                    oldPoints_,
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
        const edge& edgeToCheck = toCellEdges[edgeI];

        if
        (
            commonPoints.found(edgeToCheck.start()) &&
            commonPoints.found(edgeToCheck.end())
        )
        {
            continue;
        }

        forAll(fromCell, faceI)
        {
            const face& faceToCheck = polyMesh::faces()[fromCell[faceI]];

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
                // Avoid common points, since this implies that
                // the edge intersects at a face point
                if (faceToCheck.which(pIter.key()) > -1)
                {
                    if
                    (
                        (edgeToCheck[0] == pIter.key()) ||
                        (edgeToCheck[1] == pIter.key())
                    )
                    {
                        foundCommon = true;
                        break;
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
                }
            }

            if (foundCommon)
            {
                continue;
            }

            foundIntersection = false;

            foundIntersection =
            (
                meshOps::segmentFaceIntersection
                (
                    edgeToCheck,
                    faceToCheck,
                    oldPoints_,
                    polyMesh::points(),
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

    if (debug > 3)
    {
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


// Obtain a list of possible parent faces from the old mesh
labelList dynamicTopoFvMesh::faceParents
(
    const label fIndex,
    const vector& fNormal,
    const vector& fCentre,
    const scalar searchFactor,
    const labelList& oldCandidates
) const
{
    const face& faceToCheck = faces_[fIndex];
    const polyBoundaryMesh& boundary = boundaryMesh();

    // Determine the patch that this face belongs to..
    label patchIndex = whichPatch(fIndex);

    if (patchIndex == -1 || faceToCheck.empty())
    {
        FatalErrorIn
        (
            "\n"
            "labelList dynamicTopoFvMesh::faceParents\n"
            "(\n"
            "    const label fIndex,\n"
            "    const vector& fCentre,\n"
            "    const scalar searchFactor,\n"
            "    const labelList& oldCandidates\n"
            ") const"
        )
            << nl << " Illegal request for face: "
            << fIndex << ":: " << faces_[fIndex] << nl
            << " oldCandidates: " << oldCandidates
            << abort(FatalError);
    }

    // Determine the bounds of the new face
    scalar maxDist = 0.0;
    vector maxDistPos = vector::zero;

    forAll(faceToCheck, pointI)
    {
        scalar dist = magSqr(oldPoints_[faceToCheck[pointI]] - fCentre);

        if (dist > maxDist)
        {
            maxDist = dist;
            maxDistPos = oldPoints_[faceToCheck[pointI]];
        }
    }

    // Define a search radius
    vector bMax = searchFactor * (maxDistPos - fCentre);

    StaticHashTable<empty, label, Hash<label> > masterFaces;

    // Fetch old patch start
    label patchStart = boundary[patchIndex].start();

    // Insert the old candidates first
    // Assume that candidates / parents are global face indices,
    // and all addressing into the same patch.
    forAll(oldCandidates, faceI)
    {
        if (whichPatch(oldCandidates[faceI]) == patchIndex)
        {
            masterFaces.insert
            (
                (oldCandidates[faceI] - patchStart),
                empty()
            );
        }
    }

    // Fetch geometry / connectivity from the old mesh.
    const vectorField& faceNormals = boundary[patchIndex].faceAreas();
    const vectorField& faceCentres = boundary[patchIndex].faceCentres();
    const labelListList& oldFaceFaces = boundary[patchIndex].faceFaces();

    label nAttempts = 0;
    bool changed;

    do
    {
        // Reset flag
        changed = false;

        // Fetch the initial set of candidates
        labelList initList = masterFaces.toc();

        // Accumulate a larger stencil of face neighbours
        forAll(initList, indexI)
        {
            const labelList& ff = oldFaceFaces[initList[indexI]];

            forAll(ff, faceI)
            {
                scalar dP = (faceNormals[ff[faceI]] & fNormal);
                vector xC = (faceCentres[ff[faceI]] - fCentre);

                if (((xC & xC) < (bMax & bMax)) && (dP > 0.0))
                {
                    if (!masterFaces.found(ff[faceI]))
                    {
                        masterFaces.insert(ff[faceI], empty());
                        changed = true;
                    }
                }
            }
        }

        nAttempts++;

    } while (changed);

    if (debug > 3)
    {
        Info << " Face: " << fIndex
             << " No. of parent candidates: "
             << masterFaces.size()
             << " searchFactor: "
             << searchFactor
             << endl;
    }

    // Now renumber all masterFaces back to global indices
    labelList finalList = masterFaces.toc();

    forAll(finalList, faceI)
    {
        finalList[faceI] += patchStart;
    }

    return finalList;
}


// Obtain a list of possible parent cells from the old mesh.
labelList dynamicTopoFvMesh::cellParents
(
    const label cIndex,
    const vector& cCentre,
    const scalar searchFactor,
    const labelList& oldCandidates
) const
{
    // Fetch the new cell, and determine its bounds.
    const cell& newCell = cells_[cIndex];

    scalar maxDist = 0.0;
    vector maxDistPos = vector::zero;

    forAll(newCell, faceI)
    {
        const face& faceToCheck = faces_[newCell[faceI]];

        forAll(faceToCheck, pointI)
        {
            scalar dist = magSqr(oldPoints_[faceToCheck[pointI]] - cCentre);

            if (dist > maxDist)
            {
                maxDist = dist;
                maxDistPos = oldPoints_[faceToCheck[pointI]];
            }
        }
    }

    // Define a search radius
    vector bMax = searchFactor * (maxDistPos - cCentre);

    StaticHashTable<empty, label, Hash<label> > masterCells;

    // Insert the old candidates first
    forAll(oldCandidates, cellI)
    {
        masterCells.insert(oldCandidates[cellI], empty());
    }

    // Fetch connectivity from the old mesh.
    const vectorField& cellCentres = polyMesh::cellCentres();
    const labelListList& oldCellCells = polyMesh::cellCells();

    label nAttempts = 0;
    bool changed;

    do
    {
        // Reset flag
        changed = false;

        // Fetch the initial set of candidates
        labelList initList = masterCells.toc();

        // Accumulate a larger stencil of cell neighbours
        forAll(initList, indexI)
        {
            const labelList& cc = oldCellCells[initList[indexI]];

            forAll(cc, cellI)
            {
                vector xC = (cellCentres[cc[cellI]] - cCentre);

                if ((xC & xC) < (bMax & bMax))
                {
                    if (!masterCells.found(cc[cellI]))
                    {
                        masterCells.insert(cc[cellI], empty());
                        changed = true;
                    }
                }
            }
        }

        nAttempts++;

    } while (changed);

    if (debug > 3)
    {
        Info << " Cell: " << cIndex
             << " No. of parent candidates: "
             << masterCells.size()
             << " searchFactor: "
             << searchFactor
             << endl;
    }

    return masterCells.toc();
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
