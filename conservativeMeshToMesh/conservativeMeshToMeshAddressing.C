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

#include "triFace.H"
#include "IOmanip.H"
#include "ListOps.H"
#include "clockTime.H"
#include "tetPointRef.H"
#include "SortableList.H"

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
    scalar interval = 0.5, oIndex = 0.0, nIndex = 0.0;

    oIndex = ::floor(sTimer.elapsedTime() / interval);

    for (label cellI = cellStart; cellI < (cellStart + cellSize); cellI++)
    {
        count++;

        // Update the index, if its changed
        nIndex = ::floor(sTimer.elapsedTime() / interval);

        if ((nIndex - oIndex) > VSMALL)
        {
            oIndex = nIndex;

            scalar percent = 100.0 * (double(count) / (cellSize + VSMALL));

            // Report progress
            if (report)
            {
                Info << "  Progress: " << percent << "% : "
                     << "  Cells processed: " << count
                     << "  out of " << cellSize << " total."
                     << "             \r"
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
        Info << "  Progress: 100%"
             << "  Entities processed: " << count
             << "  out of " << cellSize << " total."
             << "             \r"
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
    //- Geometric match tolerance
    bool changed;
    scalar matchTol = 1e-4;;
    label nAttempts = 0, nIntersects = 0;

    label mapCandidate = -1;

    if (oldCandidate < 0)
    {
        // Loop through all cells, and find minimum
        scalar minDist = GREAT;
        label minCell = -1;

        const vectorField& oldCentres = fromMesh().cellCentres();
        const vector newCentre = toMesh().cellCentres()[newCellIndex];

        forAll(oldCentres, cellI)
        {
            if (magSqr(newCentre - oldCentres[cellI]) < minDist)
            {
                minDist = magSqr(newCentre - oldCentres[cellI]);
                minCell = cellI;
            }
        }

        // Reassign old candidate
        mapCandidate = minCell;
    }
    else
    {
        mapCandidate = oldCandidate;
    }

    // Fetch the volume of the new cell
    scalar newCellVolume = toMesh().cellVolumes()[newCellIndex];

    // Fetch cellCells addressing
    const labelListList& oldNeighbourList = fromMesh().cellCells();

    // Maintain a check-list
    labelHashSet checked, skipped;

    // Loop and add intersections until nothing changes
    do
    {
        // Reset flag
        changed = false;

        // Fetch the set of candidates
        labelList checkList;

        if (nAttempts == 0)
        {
            checkList = labelList(1, mapCandidate);
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
                checkEntities = labelList(1, checkList[indexI]);
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

                vectorField intPoints(0);

                bool intersect =
                (
                    cellIntersection
                    (
                        newCellIndex,
                        checkEntity,
                        matchTol,
                        intPoints
                    )
                );

                if (intersect)
                {
                    scalar weight = 0.0;
                    vector centre = vector::zero;

                    // Compute weights
                    convexSetVolume
                    (
                        newCellIndex,
                        checkEntity,
                        intPoints,
                        weight,
                        centre
                    );

                    label oldSize = parents.size();

                    parents.setSize(oldSize + 1, checkEntity);
                    weights.setSize(oldSize + 1, weight);
                    centres.setSize(oldSize + 1, centre);

                    nIntersects++;

                    if (!checked.found(checkEntity))
                    {
                        checked.insert(checkEntity);
                    }

                    changed = true;
                }
                else
                if (nAttempts == 0)
                {
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
                        << endl;
                }
                else
                {
                    // Add to the skipped list
                    if (!skipped.found(checkEntity))
                    {
                        skipped.insert(checkEntity);
                    }
                }
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
    if (mag(newCellVolume - sum(weights)) > 1e-17)
    {
        // Write out for post-processing
        labelList uList = skipped.toc();

        writeVTK("nE_" + Foam::name(newCellIndex), newCellIndex);
        writeVTK("oE_" + Foam::name(newCellIndex), mapCandidate, 3, true);
        writeVTK("mE_" + Foam::name(newCellIndex), parents, 3, true);
        writeVTK("uE_" + Foam::name(newCellIndex), uList, 3, true);
    }

    /*
    // Check if tolerances were adjusted.
    if (nInnerAttempts)
    {
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

        // Reset match tolerance, if necessary
        matchTol_ *= Foam::pow(10, nInnerAttempts);
    }
    */
}


// Return the intersection points between cells in old/new meshes
template <class T>
bool conservativeMeshToMesh::cellIntersection
(
    const label newIndex,
    const label oldIndex,
    const T& matchTol,
    Field<Vector<T> >& intPoints,
    bool output
) const
{
    // Reset inputs
    intPoints.clear();

    // Assume XY plane here for 2D meshes
    Vector<T> planeNormal = convert<T>(vector(0,0,1));

    // Fetch references for each mesh
    const cell& oldCell = fromMesh().cells()[oldIndex];
    const faceList& oldFaces = fromMesh().faces();
    const pointField& oldPoints = fromMesh().points();
    const labelList& oldOwner = fromMesh().faceOwner();
    const edgeList oldCellEdges = oldCell.edges(oldFaces);
    const labelList oldCellPoints = oldCell.labels(oldFaces);

    const cell& newCell = toMesh().cells()[newIndex];
    const faceList& newFaces = toMesh().faces();
    const pointField& newPoints = toMesh().points();
    const labelList& newOwner = toMesh().faceOwner();
    const edgeList newCellEdges = newCell.edges(newFaces);
    const labelList newCellPoints = newCell.labels(newFaces);

    // Track all possible intersections from here on.
    label nInts = 0;
    Map<Vector<T> > intersections;
    Vector<T> intPoint = Vector<T>::zero;
    Vector<T> checkPoint = Vector<T>::zero;

    FixedList<vector,2> segment(vector::zero);

    scalar minEdgeMag = GREAT;

    forAll(oldCellEdges, edgeI)
    {
        scalar edgeMag = mag(oldCellEdges[edgeI].vec(oldPoints));

        minEdgeMag = Foam::min(edgeMag, minEdgeMag);
    }

    forAll(newCellEdges, edgeI)
    {
        scalar edgeMag = mag(newCellEdges[edgeI].vec(newPoints));

        minEdgeMag = Foam::min(edgeMag, minEdgeMag);
    }

    scalar pointMergeTol = magSqr(matchTol) * minEdgeMag;

    // Check if any points are coincident.
    Map<labelList> OtoN, NtoO;

    forAll(oldCellPoints, pointI)
    {
        forAll(newCellPoints, pointJ)
        {
            if
            (
                mag
                (
                    oldPoints[oldCellPoints[pointI]]
                  - newPoints[newCellPoints[pointJ]]
                ) < pointMergeTol
            )
            {
                OtoN.insert
                (
                    oldCellPoints[pointI],
                    labelList(1, newCellPoints[pointJ])
                );

                NtoO.insert
                (
                    newCellPoints[pointJ],
                    labelList(1, oldCellPoints[pointI])
                );

                intersections.set(++nInts, newPoints[newCellPoints[pointJ]]);
            }
        }
    }

    // If all points are common, this is identical to the old cell.
    if (nInts == oldCellPoints.size())
    {
        // Copy intersections
        intPoints.setSize(nInts, vector::zero);

        nInts = 0;

        forAllConstIter(Map<vector>, intersections, pI)
        {
            intPoints[nInts++] = pI();
        }

        return true;
    }

    // Check for point-segment intersections
    bool pointIntersections = false;

    forAll(oldCellPoints, pointI)
    {
        label oldPoint = oldCellPoints[pointI];

        if (OtoN.found(oldPoint))
        {
            continue;
        }

        checkPoint = convert<T>(oldPoints[oldPoint]);

        forAll(newCellEdges, edgeI)
        {
            const edge& edgeToCheck = newCellEdges[edgeI];

            if
            (
                meshOps::pointSegmentIntersection
                (
                    line<Vector<T>, const Vector<T>&>
                    (
                        convert<T>(newPoints[edgeToCheck.start()]),
                        convert<T>(newPoints[edgeToCheck.end()])
                    ),
                    checkPoint,
                    matchTol
                )
            )
            {
                OtoN.insert
                (
                    oldPoint,
                    labelList(edgeToCheck)
                );

                intersections.set(++nInts, checkPoint);

                pointIntersections = true;
            }
        }
    }

    forAll(newCellPoints, pointI)
    {
        label newPoint = newCellPoints[pointI];

        if (NtoO.found(newPoint))
        {
            continue;
        }

        checkPoint = convert<T>(newPoints[newPoint]);

        forAll(oldCellEdges, edgeI)
        {
            const edge& edgeToCheck = oldCellEdges[edgeI];

            if
            (
                meshOps::pointSegmentIntersection
                (
                    line<Vector<T>, const Vector<T>&>
                    (
                        convert<T>(oldPoints[edgeToCheck.start()]),
                        convert<T>(oldPoints[edgeToCheck.end()])
                    ),
                    checkPoint,
                    matchTol
                )
            )
            {
                NtoO.insert
                (
                    newPoint,
                    labelList(edgeToCheck)
                );

                intersections.set(++nInts, checkPoint);

                pointIntersections = true;
            }
        }
    }

    if (pointIntersections && output)
    {
        Info << "Point Intersections exist: " << nl
             << " newCellIndex: " << newIndex
             << " oldCellIndex: " << oldIndex
             << endl;
    }

    if (twoDMesh_)
    {
        // Check if edge mid-points are clearly within the cell.
        // If so, add edge points as 'intersections'.
        forAll(oldCellEdges, edgeI)
        {
            const edge& edgeToCheck = oldCellEdges[edgeI];

            if
            (
                OtoN.found(edgeToCheck.start()) &&
                OtoN.found(edgeToCheck.end())
            )
            {
                continue;
            }

            vector edgeVec = edgeToCheck.vec(oldPoints);

            edgeVec /= mag(edgeVec) + VSMALL;

            if (mag(edgeVec & planeNormal) < 0.5)
            {
                continue;
            }

            vector checkPoint =
            (
                0.5 *
                (
                    oldPoints[edgeToCheck.start()] +
                    oldPoints[edgeToCheck.end()]
                )
            );

            if
            (
                meshOps::pointInCell
                (
                    newIndex,
                    newCell,
                    newFaces,
                    newOwner,
                    newPoints,
                    0.0,
                    checkPoint
                )
            )
            {
                intersections.set(++nInts, oldPoints[edgeToCheck.start()]);
                intersections.set(++nInts, oldPoints[edgeToCheck.end()]);
            }
        }

        forAll(newCellEdges, edgeI)
        {
            const edge& edgeToCheck = newCellEdges[edgeI];

            if
            (
                NtoO.found(edgeToCheck.start()) &&
                NtoO.found(edgeToCheck.end())
            )
            {
                continue;
            }

            vector edgeVec = edgeToCheck.vec(newPoints);

            edgeVec /= mag(edgeVec) + VSMALL;

            if (mag(edgeVec & planeNormal) < 0.5)
            {
                continue;
            }

            vector checkPoint =
            (
                0.5 *
                (
                    newPoints[edgeToCheck.start()] +
                    newPoints[edgeToCheck.end()]
                )
            );

            if
            (
                meshOps::pointInCell
                (
                    oldIndex,
                    oldCell,
                    oldFaces,
                    oldOwner,
                    oldPoints,
                    0.0,
                    checkPoint
                )
            )
            {
                intersections.set(++nInts, newPoints[edgeToCheck.start()]);
                intersections.set(++nInts, newPoints[edgeToCheck.end()]);
            }
        }
    }
    else
    {
        // Check whether any old points are within
        // the new cell. Count these as 'intersections'.
        forAll(oldCellPoints, pointI)
        {
            if (OtoN.found(oldCellPoints[pointI]))
            {
                continue;
            }

            checkPoint = oldPoints[oldCellPoints[pointI]];

            if
            (
                meshOps::pointInCell
                (
                    newIndex,
                    newCell,
                    newFaces,
                    newOwner,
                    newPoints,
                    0.0,
                    checkPoint
                )
            )
            {
                intersections.set(++nInts, checkPoint);
            }
        }

        // Check whether any new points are within
        // the old cell. Count these as 'intersections'.
        forAll(newCellPoints, pointI)
        {
            if (NtoO.found(newCellPoints[pointI]))
            {
                continue;
            }

            checkPoint = newPoints[newCellPoints[pointI]];

            if
            (
                meshOps::pointInCell
                (
                    oldIndex,
                    oldCell,
                    oldFaces,
                    oldOwner,
                    oldPoints,
                    0.0,
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
    List<Pair<edge> > OeToNe, NeToOe;

    // Define edge-vectors in 2D
    Vector<T> oldVec(Vector<T>::zero), newVec(Vector<T>::zero);

    forAll(oldCellEdges, edgeI)
    {
        // For 2D meshes, only select edges on wedge/empty planes
        if (twoDMesh_)
        {
            oldVec =
            (
                convert<T>(oldPoints[oldCellEdges[edgeI].start()])
              - convert<T>(oldPoints[oldCellEdges[edgeI].end()])
            );

            oldVec /= mag(oldVec) + VSMALL;

            if (mag(oldVec & planeNormal) > 0.5)
            {
                continue;
            }
        }

        const edge& oldEdge = oldCellEdges[edgeI];

        forAll(newCellEdges, edgeJ)
        {
            // For 2D meshes, only select edges on wedge/empty planes
            if (twoDMesh_)
            {
                newVec =
                (
                    convert<T>(newPoints[newCellEdges[edgeJ].start()])
                  - convert<T>(newPoints[newCellEdges[edgeJ].end()])
                );

                newVec /= mag(newVec) + VSMALL;

                if (mag(newVec & planeNormal) > 0.5)
                {
                    continue;
                }
            }

            const edge& newEdge = newCellEdges[edgeJ];

            // Form an edge-pair
            Pair<edge> edgePair(oldEdge, newEdge);

            // If an edge-point was found in common-points or
            // point-segment intersections, skip this pair
            if (OtoN.found(oldEdge.start()))
            {
                if (findInEdge(OtoN[oldEdge.start()], newEdge))
                {
                    continue;
                }
            }

            if (OtoN.found(oldEdge.end()))
            {
                if (findInEdge(OtoN[oldEdge.end()], newEdge))
                {
                    continue;
                }
            }

            if (NtoO.found(newEdge.start()))
            {
                if (findInEdge(NtoO[newEdge.start()], oldEdge))
                {
                    continue;
                }
            }

            if (NtoO.found(newEdge.end()))
            {
                if (findInEdge(NtoO[newEdge.end()], oldEdge))
                {
                    continue;
                }
            }

            foundIntersection = false;

            foundIntersection =
            (
                meshOps::segmentSegmentIntersection
                (
                    line<Vector<T>, const Vector<T>&>
                    (
                        convert<T>(oldPoints[oldEdge.start()]),
                        convert<T>(oldPoints[oldEdge.end()])
                    ),
                    line<Vector<T>, const Vector<T>&>
                    (
                        convert<T>(newPoints[newEdge.start()]),
                        convert<T>(newPoints[newEdge.end()])
                    ),
                    matchTol,
                    intPoint
                )
            );

            if (foundIntersection)
            {
                OeToNe.setSize
                (
                    OeToNe[edgeI].size() + 1,
                    edgePair
                );

                NeToOe.setSize
                (
                    NeToOe[edgeJ].size() + 1,
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
        intPoints.setSize(nInts, Vector<T>::zero);

        nInts = 0;

        forAllConstIter(typename Map<Vector<T> >, intersections, pI)
        {
            intPoints[nInts++] = pI();
        }

        if (output)
        {
            meshOps::checkPointNearness(intPoints, T(1e-20));

            meshOps::writeVTK
            (
                toMesh(),
                "ccSet_" + Foam::name(newIndex)
              + '<' + Foam::name(oldIndex) + '>',
                intPoints.size(),
                intPoints.size(),
                intPoints.size(),
                convert<scalar>(intPoints)
            );
        }

        // Found a convex set of points
        return true;
    }

    if (edgeIntersections && output)
    {
        Info << "Edge Intersections exist: " << nl
             << " newCellIndex: " << newIndex
             << " oldCellIndex: " << oldIndex
             << endl;
    }

    // Loop through all old edges, and find possible
    // intersections with faces of the new cell.
    forAll(oldCellEdges, edgeI)
    {
        const edge& edgeToCheck = oldCellEdges[edgeI];

        if
        (
            OtoN.found(edgeToCheck.start()) &&
            OtoN.found(edgeToCheck.end())
        )
        {
            continue;
        }

        forAll(newCell, faceI)
        {
            const triFace faceToCheck(newFaces[newCell[faceI]]);

            // Avoid point-edge / edge-edge intersections, if any.
            if (edgeIntersections)
            {
                // Is edgeToCheck in the list?
                bool foundEdge = false;
                const edgeList fEdges = faceToCheck.edges();

                forAll(OeToNe, indexI)
                {
                    if (OeToNe[indexI].first() == edgeToCheck)
                    {
                        // Check whether the intersecting edge
                        // exists on this face.
                        forAll(fEdges, edgeJ)
                        {
                            if (fEdges[edgeJ] == OeToNe[indexI].second())
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

            // Avoid common points / points-on-edge,
            // since this implies that the edge
            // intersects at a face point
            bool foundCommon = false;

            forAllConstIter(Map<labelList>, NtoO, pIter)
            {
                if (findIndex(faceToCheck, pIter.key()) > -1)
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

            foundIntersection = false;

            foundIntersection =
            (
                meshOps::segmentTriFaceIntersection
                (
                    triangle<Vector<T>, const Vector<T>&>
                    (
                        convert<T>(newPoints[faceToCheck[0]]),
                        convert<T>(newPoints[faceToCheck[1]]),
                        convert<T>(newPoints[faceToCheck[2]])
                    ),
                    line<Vector<T>, const Vector<T>&>
                    (
                        convert<T>(oldPoints[edgeToCheck.start()]),
                        convert<T>(oldPoints[edgeToCheck.end()])
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
    forAll(newCellEdges, edgeI)
    {
        const edge& edgeToCheck = newCellEdges[edgeI];

        if
        (
            NtoO.found(edgeToCheck.start()) &&
            NtoO.found(edgeToCheck.end())
        )
        {
            continue;
        }

        forAll(oldCell, faceI)
        {
            const triFace faceToCheck(oldFaces[oldCell[faceI]]);

            // Avoid point-edge / edge-edge intersections, if any.
            if (edgeIntersections)
            {
                // Is edgeToCheck in the list?
                bool foundEdge = false;
                const edgeList fEdges = faceToCheck.edges();

                forAll(NeToOe, indexI)
                {
                    if (NeToOe[indexI].first() == edgeToCheck)
                    {
                        // Check whether the intersecting edge
                        // exists on this face.
                        forAll(fEdges, edgeJ)
                        {
                            if (fEdges[edgeJ] == NeToOe[indexI].second())
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

            // Avoid common points / points-on-edge,
            // since this implies that the edge
            // intersects at a face point
            bool foundCommon = false;

            forAllConstIter(Map<labelList>, OtoN, pIter)
            {
                if (findIndex(faceToCheck, pIter.key()) > -1)
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

            foundIntersection = false;

            foundIntersection =
            (
                meshOps::segmentTriFaceIntersection
                (
                    triangle<Vector<T>, const Vector<T>&>
                    (
                        convert<T>(oldPoints[faceToCheck[0]]),
                        convert<T>(oldPoints[faceToCheck[1]]),
                        convert<T>(oldPoints[faceToCheck[2]])
                    ),
                    line<Vector<T>, const Vector<T>&>
                    (
                        convert<T>(newPoints[edgeToCheck.start()]),
                        convert<T>(newPoints[edgeToCheck.end()])
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
    intPoints.setSize(nInts, vector::zero);

    nInts = 0;

    forAllConstIter(typename Map<Vector<T> >, intersections, pI)
    {
        intPoints[nInts++] = pI();
    }

    if (output)
    {
        meshOps::checkPointNearness(intPoints, T(1e-20));

        meshOps::writeVTK
        (
            toMesh(),
            "ccSet_" + Foam::name(newIndex)
          + '<' + Foam::name(oldIndex) + '>',
            intPoints.size(),
            intPoints.size(),
            intPoints.size(),
            convert<scalar>(intPoints)
        );
    }

    // Found a convex set of points.
    return true;
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
    DynamicList<label> uniquePts(10);

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
                            uniquePts.append(tmpFace[pI]);
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
                                    meshOps::insertPointLabels
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
                            uniquePts.append(testFace[pI]);
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
                            meshOps::insertPointLabels
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
                            uniquePts.append(checkFace[pI]);
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
                            meshOps::insertPointLabels
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
