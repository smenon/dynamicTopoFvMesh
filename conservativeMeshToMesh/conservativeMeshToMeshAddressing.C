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

Description
    Private member of conservativeMeshToMesh.

    Calculates mesh to mesh addressing pattern (for each cell from one mesh
    find the closest cells in the other mesh).

\*---------------------------------------------------------------------------*/

#include "conservativeMeshToMesh.H"
#include "IOmanip.H"
#include "ListOps.H"
#include "clockTime.H"
#include "tetPointRef.H"
#include "SortableList.H"
#include "StaticHashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void conservativeMeshToMesh::calcAddressingAndWeights
(
    const label cellStart,
    const label cellSize,
    bool report
)
{
    if (debug)
    {
        Info << "conservativeMeshToMesh::calculateIntersectionAddressing() : "
             << "calculating mesh-to-mesh cell addressing" << endl;
    }

    // Fetch nearest-cell addressing from meshToMesh
    const labelList& cAddr = meshToMesh::cellAddressing();

    clockTime sTimer;
    bool reported = false;
    label count = 0, oldStart = cellStart;
    scalar interval = 2.0, oIndex = 0.0, nIndex = 0.0;

    oIndex = ::floor(sTimer.elapsedTime() / interval);

    for (label cellI = cellStart; cellI < (cellStart + cellSize); cellI++)
    {
        count++;

        // Update the index, if its changed
        nIndex = ::floor(sTimer.elapsedTime() / interval);

        if ((nIndex - oIndex) > VSMALL)
        {
            oIndex = nIndex;

            // Report progress
            if (report)
            {
                Info << "\r  Progress: "
                     << 100.0 * (double(count) / (cellSize + VSMALL))
                     << "% : Cells processed: " << count
                     << " out of " << cellSize << " total."
                     << "             "
                     << flush;

                reported = true;
            }
        }

        // Fetch references
        labelList& parents = addressing_[cellI];
        scalarField& weights = weights_[cellI];
        vectorField& centres = centres_[cellI];

        // Obtain weighting factors for this cell.
        computeCellWeights
        (
            cellI,
            cAddr[cellI],
            parents,
            weights,
            centres
        );

        if ((cellI - oldStart) > 50)
        {
            ctrMutex_.lock();

            counter_ += (cellI - oldStart);

            ctrMutex_.unlock();

            // Reset start index
            oldStart = cellI;
        }
    }

    // Final addition to counter
    ctrMutex_.lock();

    counter_ += ((cellStart + cellSize) - oldStart);

    ctrMutex_.unlock();

    if (reported && report)
    {
        Info << "\r  Progress: "
             << 100.0 * (double(count) / (cellSize + VSMALL))
             << "% : Entities processed: " << count
             << " out of " << cellSize << " total."
             << "             "
             << endl;
    }
}


// Invert addressing from source to target
bool conservativeMeshToMesh::invertAddressing()
{
    // Read in addressing arrays from the source mesh
    IOList<labelList> srcAddressing
    (
        IOobject
        (
            "addressing",
            fromMesh().time().timeName(),
            fromMesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Check for compatibility
    bool compatible = true;
    label targetCells = toMesh().nCells();
    labelList nCellsPerCell(targetCells, 0);

    forAll(srcAddressing, cellI)
    {
        const labelList& srcAddr = srcAddressing[cellI];

        label maxIndex = findMax(srcAddr);

        if (srcAddressing[cellI][maxIndex] >= targetCells)
        {
            Info << " Addressing is not compatible. " << endl;

            compatible = false;
            break;
        }

        forAll(srcAddr, j)
        {
            nCellsPerCell[srcAddr[j]]++;
        }
    }

    if (!compatible)
    {
        addressing_.clear();
        addressing_.setSize(targetCells);

        return false;
    }

    // Read weights
    IOList<scalarField> srcWeights
    (
        IOobject
        (
            "weights",
            fromMesh().time().timeName(),
            fromMesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read centres
    IOList<vectorField> srcCentres
    (
        IOobject
        (
            "centres",
            fromMesh().time().timeName(),
            fromMesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Set sizes
    forAll(nCellsPerCell, cellI)
    {
        addressing_[cellI].setSize(nCellsPerCell[cellI]);
        weights_[cellI].setSize(nCellsPerCell[cellI]);
        centres_[cellI].setSize(nCellsPerCell[cellI]);
    }

    nCellsPerCell = 0;

    // Invert addressing
    forAll(srcAddressing, cellI)
    {
        const labelList& srcAddr = srcAddressing[cellI];
        const scalarField& srcWt = srcWeights[cellI];
        const vectorField& srcCt = srcCentres[cellI];

        forAll(srcAddr, j)
        {
            label cellJ = srcAddr[j];

            addressing_[cellJ][nCellsPerCell[cellJ]] = cellI;
            weights_[cellJ][nCellsPerCell[cellJ]] = srcWt[j];
            centres_[cellJ][nCellsPerCell[cellJ]] = srcCt[j];

            nCellsPerCell[cellJ]++;
        }
    }

    // Check weights for consistency
    const scalarField& V = toMesh().cellVolumes();

    forAll(V, cellI)
    {
        if (mag(V[cellI] - sum(weights_[cellI])) > 1e-16)
        {
            Info << " Weights are not compatible. " << nl
                 << " Cell: " << cellI << nl
                 << " Volume: " << V[cellI] << nl
                 << " Sum(weights): " << sum(weights_[cellI]) << nl
                 << " Error: " << mag(V[cellI] - sum(weights_[cellI]))
                 << endl;

            compatible = false;
            break;
        }
    }

    if (!compatible)
    {
        addressing_.clear();
        weights_.clear();
        centres_.clear();

        addressing_.setSize(targetCells);
        weights_.setSize(targetCells);
        centres_.setSize(targetCells);

        return false;
    }

    // Return a successful inversion
    return true;
}


// Compute cell weighting factors for a particular cell
void conservativeMeshToMesh::computeCellWeights
(
    const label newCellIndex,
    const label oldCandidate,
    labelList& parents,
    scalarField& weights,
    vectorField& centres
) const
{
    scalar searchFactor = 1.5;

    label nOldIntersects = -1, nIntersects = 0;
    label nAttempts = 0, nInnerAttempts = 0;

    // Maintain a list of candidates and intersection points
    boolList oldIntersects, intersects;
    labelList oldCandidates, candidates;

    // Fetch the volume of the new cell
    scalar newCellVolume = toMesh().cellVolumes()[newCellIndex];

    while (nAttempts < 10)
    {
        // Reset counter first
        nIntersects = 0;

        // Obtain candidate parents for this cell
        candidates =
        (
            cellParents
            (
                newCellIndex,
                searchFactor,
                oldCandidate
            )
        );

        // Set sizes
        intersects.setSize(candidates.size());

        // Test for intersections
        forAll(candidates, indexI)
        {
            intersects[indexI] = false;

            intersects[indexI] =
            (
                testIntersection
                (
                    newCellIndex,
                    candidates[indexI]
                )
            );

            if (intersects[indexI])
            {
                nIntersects++;
            }
        }

        if ((nIntersects == nOldIntersects) && (nIntersects != 0))
        {
            if (debug)
            {
                Info << " nCandidates: " << candidates.size() << nl
                     << " nIntersects: " << nIntersects
                     << endl;
            }

            // Set sizes
            parents.setSize(nIntersects, -1);
            weights.setSize(nIntersects, 0.0);
            centres.setSize(nIntersects, vector::zero);

            // Compute actual intersections
            nInnerAttempts = 0;
            bool attainedAccuracy = false;

            while (nInnerAttempts < 10)
            {
                // Reset counter
                nIntersects = 0;
                scalar sumVols = 0.0;

                // Loop through old candidates
                // to minimize intersection calculations
                forAll(oldIntersects, indexI)
                {
                    if (oldIntersects[indexI])
                    {
                        vectorField tP(0);

                        bool realIntersect =
                        (
                            cellIntersection
                            (
                                newCellIndex,
                                oldCandidates[indexI],
                                tP
                            )
                        );

                        // Skip false positives
                        if (realIntersect)
                        {
                            parents[nIntersects] = oldCandidates[indexI];

                            // Compute weights
                            convexSetVolume
                            (
                                newCellIndex,
                                parents[nIntersects],
                                tP,
                                weights[nIntersects],
                                centres[nIntersects]
                            );

                            // Accumulate volume
                            sumVols += weights[nIntersects];

                            nIntersects++;
                        }
                    }
                }

                // Check for consistency
                if (mag(1.0 - (sumVols/newCellVolume)) > 1e-10)
                {
                    // Reduce geometric tolerance, and try again.
                    matchTol_ *= 0.1;
                    nInnerAttempts++;
                }
                else
                {
                    // Attained sufficient accuracy
                    attainedAccuracy = true;
                    break;
                }
            }

            // Shorten to actual sizes
            parents.setSize(nIntersects);
            weights.setSize(nIntersects);
            centres.setSize(nIntersects);

            if (attainedAccuracy)
            {
                break;
            }
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
    if (mag(1.0 - sum(weights/newCellVolume)) > 1e-10)
    {
        // Write out for post-processing
        label uIdx = 0, cellI = newCellIndex;
        labelList unMatch(oldCandidates.size() - parents.size(), -1);

        forAll(oldCandidates, cI)
        {
            if (findIndex(parents, oldCandidates[cI]) == -1)
            {
                unMatch[uIdx++] = oldCandidates[cI];
            }
        }

        writeVTK("nCell_" + Foam::name(cellI), toMesh(), cellI, 3);
        writeVTK("oCell_" + Foam::name(cellI), fromMesh(), oldCandidates, 3);
        writeVTK("mCell_" + Foam::name(cellI), fromMesh(), parents, 3);
        writeVTK("uCell_" + Foam::name(cellI), fromMesh(), unMatch, 3);

        // Write out intersection points
        forAll(oldCandidates, indexI)
        {
            vectorField tP(0);

            cellIntersection
            (
                newCellIndex,
                oldCandidates[indexI],
                tP
            );

            if (tP.size() >= 4)
            {
                writeVTK(cellI, oldCandidates[indexI], tP);

                // Write out to screen
                scalar dummyWeight = 0.0;
                vector dummyCentre = vector::zero;

                convexSetVolume
                (
                    cellI,
                    oldCandidates[indexI],
                    tP,
                    dummyWeight,
                    dummyCentre,
                    true
                );
            }
        }

        // Sort list for easier post-processing
        SortableList<label> sortedParents(parents);

        // In-place reorder weights
        const labelList& indices = sortedParents.indices();

        forAll(indices, indexI)
        {
            Info << sortedParents[indexI] << ": "
                 << setprecision(16)
                 << weights[indices[indexI]]
                 << endl;
        }

        FatalErrorIn
        (
            "\n\n"
            "void conservativeMeshToMesh::computeCellWeights\n"
            "(\n"
            "    const label cIndex,\n"
            "    const labelList& mapCandidate,\n"
            "    labelList& parents,\n"
            "    scalarField& weights,\n"
            "    vectorField& centres\n"
            ") const\n"
        )
            << "Encountered non-conservative weighting factors." << nl
            << " Cell: " << newCellIndex << nl
            << " Candidate parent: " << oldCandidate << nl
            << " nCandidates: " << candidates.size() << nl
            << " nOldCandidates: " << oldCandidates.size() << nl
            << " nIntersects: " << nIntersects << nl
            << " nOldIntersects: " << nOldIntersects << nl
            << " nParents: " << parents.size() << nl
            << " nAttempts: " << nAttempts << nl
            << " nInnerAttempts: " << nInnerAttempts << nl
            << setprecision(16)
            << " New cell volume: " << newCellVolume << nl
            << " Sum(Weights): " << sum(weights) << nl
            << " Error: " << (newCellVolume - sum(weights)) << nl
            << " Norm Sum(Weights): " << sum(weights/newCellVolume) << nl
            << " Norm Error: " << mag(1.0 - sum(weights/newCellVolume))
            << abort(FatalError);
    }

    // Check if tolerances were adjusted.
    if (nInnerAttempts)
    {
        /*
        InfoIn
        (
            "\n\n"
            "void conservativeMeshToMesh::computeCellWeights\n"
            "(\n"
            "    const label cIndex,\n"
            "    const labelList& mapCandidate,\n"
            "    labelList& parents,\n"
            "    scalarField& weights,\n"
            "    vectorField& centres\n"
            ") const\n"
        )
            << " matchTol_ was reduced to: "
            << matchTol_ << nl
            << endl;
        */

        // Reset match tolerance, if necessary
        matchTol_ *= Foam::pow(10, nInnerAttempts);
    }
}


// Obtain a list of possible parent cells from the old mesh.
labelList conservativeMeshToMesh::cellParents
(
    const label newCellIndex,
    const scalar searchFactor,
    const label oldCandidate
) const
{
    // Fetch the new cell, and determine its bounds.
    const pointField& toPoints = toMesh().points();
    const cell& toCell = toMesh().cells()[newCellIndex];
    const vectorField& cellCentres = fromMesh().cellCentres();

    scalar maxDist = 0.0;
    vector maxDistPos = vector::zero;
    vector cCentre = toMesh().cellCentres()[newCellIndex];

    forAll(toCell, faceI)
    {
        const face& faceToCheck = toMesh().faces()[toCell[faceI]];

        forAll(faceToCheck, pointI)
        {
            scalar dist = magSqr(toPoints[faceToCheck[pointI]] - cCentre);

            if (dist > maxDist)
            {
                maxDist = dist;
                maxDistPos = toPoints[faceToCheck[pointI]];
            }
        }
    }

    // Define a search radius
    vector bMax = searchFactor * (maxDistPos - cCentre);

    typedef StaticHashTable<empty, label, Hash<label> > labelStaticHashSet;

    labelStaticHashSet masterCells;

    // Fetch connectivity from the old mesh.
    const labelListList& fromCellCells = fromMesh().cellCells();

    // Insert the old candidate first
    masterCells.insert(oldCandidate, empty());

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
            const labelList& cc = fromCellCells[initList[indexI]];

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

    if (debug)
    {
        Info << " Cell: " << newCellIndex
             << " No. of parent candidates: "
             << masterCells.size()
             << " searchFactor: "
             << searchFactor
             << endl;
    }

    return masterCells.toc();
}


label conservativeMeshToMesh::whichSide
(
    const labelList& cellPoints,
    const pointField& points,
    const point& dir,
    const point& p
) const
{
    label nP = 0, nN = 0;

    forAll(cellPoints, pointI)
    {
        scalar t = (dir & (points[cellPoints[pointI]] - p));

        if (t > 0.0)
        {
            nP++;
        }
        else
        if (t < 0.0)
        {
            nN++;
        }

        if (nP && nN)
        {
            return 0;
        }
    }

    return (nP ? 1 : -1);
}

// Test for intersection between two cells
//   - Uses the static separating axis test for polyhedra,
//     outlined in work by David Eberly
//     'Intersection of Convex Objects: The Method of Separating Axes'
//     http://www.geometrictools.com/
bool conservativeMeshToMesh::testIntersection
(
    const label newCellIndex,
    const label oldCellIndex
) const
{
    // Direction vector
    vector dir(vector::zero);

    // Fetch references for each mesh
    const edgeList& fromEdges = fromMesh().edges();
    const faceList& fromFaces = fromMesh().faces();
    const pointField& fromPoints = fromMesh().points();
    const labelList& fromOwner = fromMesh().faceOwner();
    const vectorField& fromAreas = fromMesh().faceAreas();
    const cell& fromCell = fromMesh().cells()[oldCellIndex];
    const labelList& fromCellEdges = fromMesh().cellEdges()[oldCellIndex];
    const labelList& fromCellPoints = fromMesh().cellPoints()[oldCellIndex];

    const edgeList& toEdges = toMesh().edges();
    const faceList& toFaces = toMesh().faces();
    const pointField& toPoints = toMesh().points();
    const labelList& toOwner = toMesh().faceOwner();
    const vectorField& toAreas = toMesh().faceAreas();
    const cell& toCell = toMesh().cells()[newCellIndex];
    const labelList& toCellEdges = toMesh().cellEdges()[newCellIndex];
    const labelList& toCellPoints = toMesh().cellPoints()[newCellIndex];

    // Test faces of oldCell for separation
    forAll(fromCell, faceI)
    {
        const label fIndex = fromCell[faceI];

        dir = fromAreas[fIndex];

        // Reverse normal if necessary
        if (fromOwner[fIndex] != oldCellIndex)
        {
            dir *= -1.0;
        }

        if
        (
            whichSide
            (
                toCellPoints,
                toPoints,
                dir,
                fromPoints[fromFaces[fIndex][0]]
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

        dir = toAreas[fIndex];

        // Reverse normal if necessary
        if (toOwner[fIndex] != newCellIndex)
        {
            dir *= -1.0;
        }

        if
        (
            whichSide
            (
                fromCellPoints,
                fromPoints,
                dir,
                toPoints[toFaces[fIndex][0]]
            ) > 0
        )
        {
            return false;
        }
    }

    // Test edges of oldCell for separation
    forAll(fromCellEdges, edgeI)
    {
        const edge& fromEdge = fromEdges[fromCellEdges[edgeI]];

        vector fromVec = fromEdge.vec(fromPoints);

        forAll(toCellEdges, edgeJ)
        {
            const edge& toEdge = toEdges[toCellEdges[edgeJ]];

            vector toVec = toEdge.vec(toPoints);

            dir = (fromVec ^ toVec);

            label firstSide =
            (
                whichSide
                (
                    fromCellPoints,
                    fromPoints,
                    dir,
                    fromPoints[fromEdge[0]]
                )
            );

            if (firstSide == 0)
            {
                continue;
            }

            label secondSide =
            (
                whichSide
                (
                    toCellPoints,
                    toPoints,
                    dir,
                    fromPoints[fromEdge[0]]
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


// Return the intersection points between cells in old/new meshes
bool conservativeMeshToMesh::cellIntersection
(
    const label newCellIndex,
    const label oldCellIndex,
    vectorField& tP
) const
{
    // Reset inputs
    tP.clear();

    // Fetch references for each mesh
    const edgeList& fromEdges = fromMesh().edges();
    const faceList& fromFaces = fromMesh().faces();
    const pointField& fromPoints = fromMesh().points();
    const cell& fromCell = fromMesh().cells()[oldCellIndex];
    const labelListList& fromFaceEdges = fromMesh().faceEdges();
    const labelList& fromCellEdges = fromMesh().cellEdges()[oldCellIndex];
    const labelList& fromCellPoints = fromMesh().cellPoints()[oldCellIndex];

    const edgeList& toEdges = toMesh().edges();
    const faceList& toFaces = toMesh().faces();
    const pointField& toPoints = toMesh().points();
    const cell& toCell = toMesh().cells()[newCellIndex];
    const labelListList& toFaceEdges = toMesh().faceEdges();
    const labelList& toCellEdges = toMesh().cellEdges()[newCellIndex];
    const labelList& toCellPoints = toMesh().cellPoints()[newCellIndex];

    // Track all possible intersections from here on.
    label nInts = 0;
    Map<vector> intersections;
    vector intPoint = vector::zero;
    FixedList<vector,2> segment(vector::zero);

    scalar minEdgeMag = GREAT;

    forAll(fromCellEdges, edgeI)
    {
        scalar edgeMag = magSqr(fromEdges[fromCellEdges[edgeI]].vec(fromPoints));

        minEdgeMag = Foam::min(edgeMag, minEdgeMag);
    }

    forAll(toCellEdges, edgeI)
    {
        scalar edgeMag = magSqr(toEdges[toCellEdges[edgeI]].vec(toPoints));

        minEdgeMag = Foam::min(edgeMag, minEdgeMag);
    }

    scalar pointMergeTol = matchTol_ * minEdgeMag;

    // Check if any points are coincident.
    Map<labelList> FtoT, TtoF;

    forAll(fromCellPoints, pointI)
    {
        forAll(toCellPoints, pointJ)
        {
            if
            (
                magSqr
                (
                    fromPoints[fromCellPoints[pointI]]
                  - toPoints[toCellPoints[pointJ]]
                ) < pointMergeTol
            )
            {
                FtoT.insert
                (
                    fromCellPoints[pointI],
                    labelList(1, toCellPoints[pointJ])
                );

                TtoF.insert
                (
                    toCellPoints[pointJ],
                    labelList(1, fromCellPoints[pointI])
                );

                intersections.set(++nInts, toPoints[toCellPoints[pointJ]]);
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

        return true;
    }

    // Check for point-segment intersections
    bool pointIntersections = false;

    forAll(fromCellPoints, pointI)
    {
        if (FtoT.found(fromCellPoints[pointI]))
        {
            continue;
        }

        const point& checkPoint = fromPoints[fromCellPoints[pointI]];

        forAll(toCellEdges, edgeI)
        {
            const edge& edgeToCheck = toEdges[toCellEdges[edgeI]];

            if
            (
                pointSegmentIntersection
                (
                    edgeToCheck,
                    toPoints,
                    checkPoint
                )
            )
            {
                FtoT.insert(fromCellPoints[pointI], labelList(edgeToCheck));

                intersections.set(++nInts, checkPoint);

                pointIntersections = true;
            }
        }
    }

    forAll(toCellPoints, pointI)
    {
        if (TtoF.found(toCellPoints[pointI]))
        {
            continue;
        }

        const point& checkPoint = toPoints[toCellPoints[pointI]];

        forAll(fromCellEdges, edgeI)
        {
            const edge& edgeToCheck = fromEdges[fromCellEdges[edgeI]];

            if
            (
                pointSegmentIntersection
                (
                    edgeToCheck,
                    fromPoints,
                    checkPoint
                )
            )
            {
                TtoF.insert(toCellPoints[pointI], labelList(edgeToCheck));

                intersections.set(++nInts, checkPoint);

                pointIntersections = true;
            }
        }
    }

    if (twoDMesh_)
    {
        // Check if edge mid-points are clearly within the cell.
        // If so, add edge points as 'intersections'.
        forAll(fromCellEdges, edgeI)
        {
            const edge& edgeToCheck = fromEdges[fromCellEdges[edgeI]];

            if (FtoT.found(edgeToCheck.start()) && FtoT.found(edgeToCheck.end()))
            {
                continue;
            }

            vector checkPoint =
            (
                0.5 * (fromPoints[edgeToCheck.start()] + fromPoints[edgeToCheck.end()])
            );

            if
            (
                pointInCell
                (
                    newCellIndex,
                    checkPoint,
                    false
                )
            )
            {
                intersections.set(++nInts, fromPoints[edgeToCheck.start()]);
                intersections.set(++nInts, fromPoints[edgeToCheck.end()]);
            }
        }

        forAll(toCellEdges, edgeI)
        {
            const edge& edgeToCheck = toEdges[toCellEdges[edgeI]];

            if (TtoF.found(edgeToCheck.start()) && TtoF.found(edgeToCheck.end()))
            {
                continue;
            }

            vector checkPoint =
            (
                0.5 * (toPoints[edgeToCheck.start()] + toPoints[edgeToCheck.end()])
            );

            if
            (
                pointInCell
                (
                    oldCellIndex,
                    checkPoint,
                    true
                )
            )
            {
                intersections.set(++nInts, toPoints[edgeToCheck.start()]);
                intersections.set(++nInts, toPoints[edgeToCheck.end()]);
            }
        }
    }
    else
    {
        // Check whether any old points are within
        // the new cell. Count these as 'intersections'.
        forAll(fromCellPoints, pointI)
        {
            if (FtoT.found(fromCellPoints[pointI]))
            {
                continue;
            }

            const point& checkPoint = fromPoints[fromCellPoints[pointI]];

            if
            (
                pointInCell
                (
                    newCellIndex,
                    checkPoint,
                    false
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
            if (TtoF.found(toCellPoints[pointI]))
            {
                continue;
            }

            const point& checkPoint = toPoints[toCellPoints[pointI]];

            if
            (
                pointInCell
                (
                    oldCellIndex,
                    checkPoint,
                    true
                )
            )
            {
                intersections.set(++nInts, checkPoint);
            }
        }
    }

    bool foundIntersection = false, edgeIntersections = false;

    // Loop through edges from each cell, and check whether they intersect.
    labelListList FeToTe(fromCellEdges.size(), labelList(0));
    labelListList TeToFe(toCellEdges.size(), labelList(0));

    forAll(fromCellEdges, edgeI)
    {
        const edge& fromEdge = fromEdges[fromCellEdges[edgeI]];

        forAll(toCellEdges, edgeJ)
        {
            const edge& toEdge = toEdges[toCellEdges[edgeJ]];

            // If an edge-point was found in common-points or
            // point-segment intersections, skip this pair
            if (FtoT.found(fromEdge.start()))
            {
                if (findInEdge(FtoT[fromEdge.start()], toEdge))
                {
                    continue;
                }
            }

            if (FtoT.found(fromEdge.end()))
            {
                if (findInEdge(FtoT[fromEdge.end()], toEdge))
                {
                    continue;
                }
            }

            if (TtoF.found(toEdge.start()))
            {
                if (findInEdge(TtoF[toEdge.start()], fromEdge))
                {
                    continue;
                }
            }

            if (TtoF.found(toEdge.end()))
            {
                if (findInEdge(TtoF[toEdge.end()], fromEdge))
                {
                    continue;
                }
            }

            foundIntersection = false;

            foundIntersection =
            (
                segmentSegmentIntersection
                (
                    fromEdge,
                    toEdge,
                    intPoint
                )
            );

            if (foundIntersection)
            {
                FeToTe[edgeI].setSize
                (
                    FeToTe[edgeI].size() + 1,
                    toCellEdges[edgeJ]
                );

                TeToFe[edgeJ].setSize
                (
                    TeToFe[edgeJ].size() + 1,
                    fromCellEdges[edgeI]
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

    // Loop through all old edges, and find possible
    // intersections with faces of the new cell.
    forAll(fromCellEdges, edgeI)
    {
        const edge& edgeToCheck = fromEdges[fromCellEdges[edgeI]];

        if (FtoT.found(edgeToCheck.start()) && FtoT.found(edgeToCheck.end()))
        {
            continue;
        }

        forAll(toCell, faceI)
        {
            // Avoid edge-edge intersections, if any.
            if (edgeIntersections)
            {
                const labelList& fEdges = toFaceEdges[toCell[faceI]];

                bool foundEdge = false;

                forAll(FeToTe[edgeI], edgeJ)
                {
                    if (findIndex(fEdges, FeToTe[edgeI][edgeJ]) > -1)
                    {
                        foundEdge = true;
                        break;
                    }
                }

                if (foundEdge)
                {
                    continue;
                }
            }

            // Avoid common points / points-on-edge,
            // since this implies that the edge
            // intersects at a face point
            const face& faceToCheck = toFaces[toCell[faceI]];

            bool foundCommon = false;

            forAllConstIter(Map<labelList>, TtoF, pIter)
            {
                if (faceToCheck.which(pIter.key()) > -1)
                {
                    if (findInEdge(pIter(), edgeToCheck))
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

            segment[0] = fromPoints[edgeToCheck[0]];
            segment[1] = fromPoints[edgeToCheck[1]];

            foundIntersection = false;

            foundIntersection =
            (
                segmentFaceIntersection
                (
                    segment,
                    edgeToCheck,
                    toCell[faceI],
                    intPoint,
                    false
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
        const edge& edgeToCheck = toEdges[toCellEdges[edgeI]];

        if (TtoF.found(edgeToCheck.start()) && TtoF.found(edgeToCheck.end()))
        {
            continue;
        }

        forAll(fromCell, faceI)
        {
            // Avoid edge-edge intersections, if any.
            if (edgeIntersections)
            {
                const labelList& fEdges = fromFaceEdges[fromCell[faceI]];

                bool foundEdge = false;

                forAll(TeToFe[edgeI], edgeJ)
                {
                    if (findIndex(fEdges, TeToFe[edgeI][edgeJ]) > -1)
                    {
                        foundEdge = true;
                        break;
                    }
                }

                if (foundEdge)
                {
                    continue;
                }
            }

            // Avoid common points / points-on-edge,
            // since this implies that the edge
            // intersects at a face point
            const face& faceToCheck = fromFaces[fromCell[faceI]];

            bool foundCommon = false;

            forAllConstIter(Map<labelList>, FtoT, pIter)
            {
                if (faceToCheck.which(pIter.key()) > -1)
                {
                    if (findInEdge(pIter(), edgeToCheck))
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

            segment[0] = toPoints[edgeToCheck[0]];
            segment[1] = toPoints[edgeToCheck[1]];

            foundIntersection = false;

            foundIntersection =
            (
                segmentFaceIntersection
                (
                    segment,
                    edgeToCheck,
                    fromCell[faceI],
                    intPoint,
                    true
                )
            );

            if (foundIntersection)
            {
                // Add to the list.
                intersections.set(++nInts, intPoint);
            }
        }
    }

    // Merge degenerate points
    forAllConstIter(Map<vector>, intersections, pI)
    {
        forAllIter(Map<vector>, intersections, pJ)
        {
            if (pI.key() == pJ.key())
            {
                continue;
            }

            if (magSqr(pI() - pJ()) < pointMergeTol)
            {
                intersections.erase(pJ);

                nInts--;
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

    // Found a convex set of points.
    if (nInts >= 4)
    {
        return true;
    }

    // Does not intersect
    return false;
}


// Find all indices in edge
bool conservativeMeshToMesh::findInEdge
(
    const labelList& cList,
    const edge& edgeToCheck
) const
{
    bool foundAll = true;

    forAll(cList, pI)
    {
        if (findIndex(edgeToCheck, cList[pI]) == -1)
        {
            foundAll = false;
            break;
        }
    }

    return foundAll;
}


// Determine whether the particular point lies
// inside the given cell
bool conservativeMeshToMesh::pointInCell
(
    const label cellIndex,
    const point& checkPoint,
    const bool useFromMesh
) const
{
    const cell* cellPtr = NULL;

    if (useFromMesh)
    {
        cellPtr = &(fromMesh().cells()[cellIndex]);
    }
    else
    {
        cellPtr = &(toMesh().cells()[cellIndex]);
    }

    // Fetch the reference
    const cell& cellToCheck = *cellPtr;

    label own = -1;
    vector xf(vector::zero), nf(vector::zero);

    forAll(cellToCheck, faceI)
    {
        if (useFromMesh)
        {
            own = fromMesh().faceOwner()[cellToCheck[faceI]];
            xf = fromMesh().faceCentres()[cellToCheck[faceI]];
            nf = fromMesh().faceAreas()[cellToCheck[faceI]];
        }
        else
        {
            own = toMesh().faceOwner()[cellToCheck[faceI]];
            xf = toMesh().faceCentres()[cellToCheck[faceI]];
            nf = toMesh().faceAreas()[cellToCheck[faceI]];
        }

        if (((xf - checkPoint) & nf) > 0.0)
        {
            if (own != cellIndex)
            {
                return false;
            }
        }
        else
        {
            if (own == cellIndex)
            {
                return false;
            }
        }
    }

    // Passed test with all faces
    return true;
}


// Determine whether a point lies on a segment
inline bool conservativeMeshToMesh::pointSegmentIntersection
(
    const edge& segment,
    const pointField& points,
    const point& checkPoint
) const
{
    vector u = points[segment.end()] - points[segment.start()];
    vector v = checkPoint - points[segment.start()];

    scalar magU = mag(u) + VSMALL;
    scalar magV = mag(v) + VSMALL;

    scalar tolerance = (matchTol_ * magU);

    // Compare dot-products
    if ( 1.0 - ((u/magU) & (v/magV)) > tolerance )
    {
        return false;
    }

    // Compute uValue
    scalar uValue = (u & v) / ((u & u) + VSMALL);

    // Does point fall-off ends?
    // Note: Does not check if end-points coincide
    //       with checkPoint, so be careful here.
    if (uValue < 0.0 || uValue > 1.0)
    {
        return false;
    }

    // Intersects segment
    return true;
}


// Determine whether a two line segments intersect
bool conservativeMeshToMesh::segmentSegmentIntersection
(
    const edge& fromSegment,
    const edge& toSegment,
    point& intPoint
) const
{
    vector u = fromSegment.vec(fromMesh().points());
    vector v = toSegment.vec(toMesh().points());

    vector p = fromMesh().points()[fromSegment.start()];
    vector q = toMesh().points()[toSegment.start()];

    vector w = (p - q);

    // Get parameters
    scalar a = (u & u), b = (u & v), c = (v & v), d = (u & w), e = (v & w);

    // Parallel / collinear check
    scalar denom = (a * c) - (b * b);

    if (mag(denom) < VSMALL)
    {
        return false;
    }

    scalar s = ((b * e) - (c * d))/(denom);
    scalar t = ((a * e) - (b * d))/(denom);

    // Out-of-bounds check
    if (s < 0.0 || t < 0.0 || s > 1.0 || t > 1.0)
    {
        return false;
    }

    // Proximity check
    scalar dist = mag( w + (s * u) - (t * v) );

    scalar tolerance = (matchTol_ * Foam::min(mag(u), mag(v)));

    if (dist > tolerance)
    {
        return false;
    }

    // Calculate the intersection point
    intPoint = fromMesh().points()[fromSegment.start()] + (s * u);

    return true;
}


// Determine whether a particular line segment
// intersects with a given face
bool conservativeMeshToMesh::segmentFaceIntersection
(
    const FixedList<vector,2>& segment,
    const edge& segmentIndices,
    const label faceIndex,
    point& intPoint,
    const bool useFromMesh
) const
{
    vector nf(vector::zero), fp(vector::zero);

    if (useFromMesh)
    {
        nf = fromMesh().faceAreas()[faceIndex];
        fp = fromMesh().points()[fromMesh().faces()[faceIndex][0]];
    }
    else
    {
        nf = toMesh().faceAreas()[faceIndex];
        fp = toMesh().points()[toMesh().faces()[faceIndex][0]];
    }

    // Normalize
    nf /= mag(nf) + VSMALL;

    // Compute uValue
    scalar numerator = nf & (fp - segment[0]);
    scalar denominator = nf & (segment[1] - segment[0]);

    // Check if the edge is parallel to the face
    if (mag(denominator) < VSMALL)
    {
        return false;
    }

    scalar u = (numerator / denominator);

    scalar tolerance = (matchTol_ * mag(segment[1] - segment[0]));

    // Check for intersection along line.
    if ((u > tolerance) && (u < (1.0 - tolerance)))
    {
        // Compute point of intersection
        intPoint = segment[0] + u*(segment[1] - segment[0]);

        // Also make sure that intPoint lies within triFace
        if (pointInFace(faceIndex, intPoint, useFromMesh))
        {
            return true;
        }
    }

    // Failed to fall within edge-bounds, or within face
    return false;
}


// Determine whether the particular point lies
// inside the given face
bool conservativeMeshToMesh::pointInFace
(
    const label faceIndex,
    const point& checkPoint,
    const bool useFromMesh
) const
{
    FixedList<vector, 2> segment(vector::zero);

    vector nf(vector::zero);
    const face* facePtr = NULL;

    if (useFromMesh)
    {
        nf = fromMesh().faceAreas()[faceIndex];
        facePtr = &(fromMesh().faces()[faceIndex]);
    }
    else
    {
        nf = toMesh().faceAreas()[faceIndex];
        facePtr = &(toMesh().faces()[faceIndex]);
    }

    // Normalize
    nf /= mag(nf) + VSMALL;

    const face& faceToCheck = *facePtr;

    forAll(faceToCheck, pI)
    {
        vector p, q;

        label nI = faceToCheck.fcIndex(pI);

        if (useFromMesh)
        {
            segment[0] = fromMesh().points()[faceToCheck[pI]];
            segment[1] = fromMesh().points()[faceToCheck[nI]];
        }
        else
        {
            segment[0] = toMesh().points()[faceToCheck[pI]];
            segment[1] = toMesh().points()[faceToCheck[nI]];
        }

        // Identify vectors
        p = (segment[1] - segment[0]);
        q = (checkPoint - segment[1]);

        scalar area = 0.5 * ((p ^ q) & nf);

        if (area < 0.0)
        {
            return false;
        }
    }

    // Passed test with all edges
    return true;
}


// Compute the volume of a polyhedron
// formed by a convex set of points.
void conservativeMeshToMesh::convexSetVolume
(
    const label newCellIndex,
    const label oldCellIndex,
    const vectorField& cvxSet,
    scalar& cVolume,
    vector& cCentre,
    bool output
) const
{
    // Reset inputs
    cVolume = 0.0;
    cCentre = vector::zero;

    // Try the trivial case for a tetrahedron.
    // No checking for orientation here.
    if (cvxSet.size() == 4)
    {
        cCentre = average(cvxSet);

        cVolume =
        (
            mag
            (
                tetPointRef
                (
                    cvxSet[0],
                    cvxSet[1],
                    cvxSet[2],
                    cvxSet[3]
                ).mag()
            )
        );

        if (output)
        {
            Info << " newCellIndex: " << newCellIndex
                 << " oldCellIndex: " << oldCellIndex << nl
                 << " Volume: " << cVolume << nl
                 << " Centre: " << cCentre << nl
                 << endl;
        }

        return;
    }

    // Track faces
    face tmpFace(3);
    label nFaces = 0;
    faceList testFaces(0);
    labelHashSet uniquePts;

    // Loop through all points, and build faces with every
    // other point in the set
    forAll(cvxSet, i)
    {
        forAll(cvxSet, j)
        {
            // Skip duplicates.
            if (j == i)
            {
                continue;
            }

            forAll(cvxSet, k)
            {
                // Skip duplicates.
                if (k == i || k == j)
                {
                    continue;
                }

                // Configure the face.
                tmpFace[0] = i;
                tmpFace[1] = j;
                tmpFace[2] = k;

                // Quick-reject test:
                //   If this is a subset of an existing face, skip it.
                bool foundSubSet = false;

                forAll(testFaces, faceI)
                {
                    const face& checkFace = testFaces[faceI];

                    if (checkFace.size() >= tmpFace.size())
                    {
                        bool foundUniquePoint = false;

                        forAll(tmpFace, pI)
                        {
                            if (findIndex(checkFace, tmpFace[pI]) == -1)
                            {
                                foundUniquePoint = true;
                                break;
                            }
                        }

                        if (!foundUniquePoint)
                        {
                            foundSubSet = true;
                            break;
                        }
                    }
                }

                if (foundSubSet)
                {
                    continue;
                }

                // Specify a tolerance for planarity
                scalar tolerance = 1e-14;

                // Compute the normal to this face
                vector n = tmpFace.normal(cvxSet);

                n /= mag(n) + VSMALL;

                label curFaceSign = 0;
                bool foundInternalFace = false;

                // Quick-reject test:
                //   Check all other points in the set,
                //   and decide if all points lie on one side.
                forAll(cvxSet, l)
                {
                    // Skip duplicates.
                    if (tmpFace[0] == l || tmpFace[1] == l || tmpFace[2] == l)
                    {
                        continue;
                    }

                    vector rfVec = (cvxSet[l] - cvxSet[i]);
                    scalar dotProd = (rfVec/(mag(rfVec) + VSMALL)) & n;

                    // Skip nearly co-planar points.
                    if (mag(dotProd) < tolerance)
                    {
                        continue;
                    }

                    // Obtain the sign of this point.
                    label fSign = Foam::sign(dotProd);

                    // Update the current sign if necessary.
                    if (curFaceSign == 0)
                    {
                        curFaceSign = fSign;
                    }
                    else
                    if (curFaceSign != fSign)
                    {
                        // Interior face. Bail out.
                        foundInternalFace = true;
                        break;
                    }
                }

                if (foundInternalFace)
                {
                    continue;
                }

                // Looks like we found a face on the boundary.
                // Check its sign to ensure that it points outward.
                if (curFaceSign == 1)
                {
                    n *= -1.0;
                    tmpFace = tmpFace.reverseFace();
                }

                // Ensure that the face wasn't checked in.
                bool alreadyCheckedIn = false;

                forAll(testFaces, faceI)
                {
                    // Fetch a non-const reference, since this face
                    // might be modified in this loop.
                    face& checkFace = testFaces[faceI];

                    label nCommon = 0;

                    uniquePts.clear();

                    forAll(tmpFace, pI)
                    {
                        if (findIndex(checkFace, tmpFace[pI]) > -1)
                        {
                            nCommon++;
                        }
                        else
                        {
                            uniquePts.insert(tmpFace[pI]);
                        }
                    }

                    if (nCommon >= 2)
                    {
                        if (checkFace.size() >= tmpFace.size())
                        {
                            // Check for unique points
                            if (uniquePts.size() > 0)
                            {
                                // Compute the existing normal
                                vector eNorm = checkFace.normal(cvxSet);

                                scalar dotProd =
                                (
                                    n & (eNorm/(mag(eNorm) + VSMALL))
                                );

                                if
                                (
                                    (mag(1.0 - dotProd) < tolerance) &&
                                    (dotProd > 0.0)
                                )
                                {
                                    // Add all unique points to checkFace
                                    insertPointLabels
                                    (
                                        n,
                                        cvxSet,
                                        uniquePts,
                                        checkFace
                                    );

                                    alreadyCheckedIn = true;
                                    break;
                                }
                            }
                            else
                            {
                                // Subset face
                                alreadyCheckedIn = true;
                                break;
                            }
                        }
                        else
                        {
                            // checkFace is a subset. Replace it.
                            checkFace = tmpFace;

                            alreadyCheckedIn = true;
                            break;
                        }
                    }
                }

                // Add this face to the list of faces.
                if (!alreadyCheckedIn)
                {
                    testFaces.setSize(++nFaces, tmpFace);
                }

                // Reset the face size.
                tmpFace.setSize(3, -1);
            }
        }
    }

    // Account for planarity test failure.
    //  - Check for subsets.
    //  - Loop until no more merges are made.
    bool changed;

    do
    {
        // Reset flag
        changed = false;

        forAll(testFaces, faceI)
        {
            // Fetch a non-const reference, since this face
            // might be modified in this loop.
            face& checkFace = testFaces[faceI];

            // Account for deleted testFaces
            if (checkFace.empty())
            {
                continue;
            }

            // Compute the normal to this face
            vector n = checkFace.normal(cvxSet);

            forAll(testFaces, faceJ)
            {
                if (faceI == faceJ)
                {
                    continue;
                }

                // Fetch a non-const reference, since this face
                // might be modified in this loop.
                face& testFace = testFaces[faceJ];

                label nCommon = 0;
                uniquePts.clear();

                if (checkFace.size() >= testFace.size())
                {
                    forAll(testFace, pI)
                    {
                        if (findIndex(checkFace, testFace[pI]) > -1)
                        {
                            nCommon++;
                        }
                        else
                        {
                            uniquePts.insert(testFace[pI]);
                        }
                    }

                    if (nCommon >= 3)
                    {
                        // Delete the test face
                        testFace.clear();

                        // Add all unique points to checkFace
                        // Failed the tolerance test before,
                        // so don't check for it now
                        if (uniquePts.size())
                        {
                            insertPointLabels
                            (
                                n,
                                cvxSet,
                                uniquePts,
                                checkFace
                            );
                        }

                        // Note that changes were made
                        changed = true;
                    }
                }
                else
                {
                    // Check if this is a subset
                    forAll(checkFace, pI)
                    {
                        if (findIndex(testFace, checkFace[pI]) > -1)
                        {
                            nCommon++;
                        }
                        else
                        {
                            uniquePts.insert(checkFace[pI]);
                        }
                    }

                    if (nCommon >= 3)
                    {
                        // This is a subset. Delete it.
                        checkFace.clear();

                        // Add all unique points to checkFace
                        // Failed the tolerance test before,
                        // so don't check for it now
                        if (uniquePts.size())
                        {
                            insertPointLabels
                            (
                                n,
                                cvxSet,
                                uniquePts,
                                testFace
                            );
                        }

                        // Note that changes were made
                        changed = true;

                        break;
                    }
                }
            }
        }
    } while (changed);

    // Find an approximate cell-centroid
    vector xC = average(cvxSet);

    // Calculate volume from all accumulated faces.
    forAll(testFaces, faceI)
    {
        const face& checkFace = testFaces[faceI];

        if (checkFace.empty())
        {
            continue;
        }

        vector xF = checkFace.centre(cvxSet);
        vector Sf = checkFace.normal(cvxSet);

        // Calculate 3*face-pyramid volume
        scalar pyr3Vol = Sf & (xF - xC);

        // Calculate face-pyramid centre
        vector pc = (3.0/4.0)*xF + (1.0/4.0)*xC;

        // Accumulate volume-weighted face-pyramid centre
        cCentre += pyr3Vol*pc;

        // Accumulate face-pyramid volume
        cVolume += pyr3Vol;
    }

    cCentre /= cVolume + VSMALL;
    cVolume *= (1.0/3.0);

    if (output)
    {
        Info << " newCellIndex: " << newCellIndex
             << " oldCellIndex: " << oldCellIndex << nl
             << " Faces: " << testFaces << nl
             << " Volume: " << cVolume << nl
             << " Centre: " << cCentre << nl
             << endl;
    }
}


// Method to insert labels in a face, so that
// right-handedness is preserved.
void conservativeMeshToMesh::insertPointLabels
(
    const vector& refNorm,
    const vectorField& points,
    const labelHashSet& pLabels,
    face& modFace
) const
{
    // Need to configure a new face.
    face newFace(modFace);

    forAllConstIter(labelHashSet, pLabels, pIter)
    {
        forAll(newFace, pI)
        {
            label nI = newFace.fcIndex(pI);

            // Compute the normal.
            vector newNorm =
            (
                triPointRef
                (
                    points[newFace[pI]],
                    points[pIter.key()],
                    points[newFace[nI]]
                ).normal()
            );

            if ((refNorm & newNorm) > 0.0)
            {
                // Insert the point.
                insertLabel
                (
                    pIter.key(),
                    newFace[pI],
                    newFace[nI],
                    newFace
                );

                break;
            }
        }
    }

    // Take over storage
    modFace.transfer(newFace);
}


// Method to insert a label between two labels in a face
// Assumes that all labels are unique.
void conservativeMeshToMesh::insertLabel
(
    const label newLabel,
    const label labelA,
    const label labelB,
    labelList& list
) const
{
    // Create a new list
    bool found = false;
    label origSize = list.size();
    labelList newList(origSize + 1);

    label index = 0, nextI = -1;

    // Start a linear search
    forAll(list, itemI)
    {
        newList[index++] = list[itemI];

        nextI = list.fcIndex(itemI);

        if
        (
            (
                (list[itemI] == labelA && list[nextI] == labelB) ||
                (list[itemI] == labelB && list[nextI] == labelA)
            ) &&
           !found
        )
        {
            found = true;
            newList[index++] = newLabel;
        }
    }

    if (!found)
    {
        FatalErrorIn("label conservativeMeshToMesh::insertLabel()")
            << "Cannot insert " << newLabel << " in list: " << list << endl
            << " Labels: "
            << labelA << " and " << labelB << " were not found in sequence."
            << abort(FatalError);
    }

    // Transfer the list
    list.transfer(newList);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
