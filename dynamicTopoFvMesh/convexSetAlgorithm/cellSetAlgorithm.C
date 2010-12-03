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
    cellSetAlgorithm

Description
    Class for convexSetAlgorithms pertaining to cells

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "cellSetAlgorithm.H"

#include "meshOps.H"
#include "triFace.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
cellSetAlgorithm::cellSetAlgorithm
(
    const polyMesh& mesh,
    const pointField& newPoints,
    const UList<edge>& newEdges,
    const UList<face>& newFaces,
    const UList<cell>& newCells,
    const UList<label>& newOwner,
    const UList<label>& newNeighbour,
    const List<objectMap>& pointsFromPoints,
    const Map<labelList>& modPoints
)
:
    convexSetAlgorithm
    (
        mesh,
        newPoints,
        newEdges,
        newFaces,
        newCells,
        newOwner,
        newNeighbour,
        pointsFromPoints,
        modPoints
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void cellSetAlgorithm::computeNormFactor(const label index) const
{
    // Clear existing fields
    parents_.clear();

    centres_.clear();
    weights_.clear();

#   if USE_MPFR

    mpCentres_.clear();
    mpWeights_.clear();

    if (highPrecision())
    {
        meshOps::cellCentreAndVolumeT
        (
            index,
            newPoints_,
            newFaces_,
            newCells_,
            newOwner_,
            mpRefNorm_,
            mpNormFactor_
        );
    }
    else
#   endif
    {
        meshOps::cellCentreAndVolume
        (
            index,
            newPoints_,
            newFaces_,
            newCells_,
            newOwner_,
            refNorm_,
            normFactor_
        );
    }
}


bool cellSetAlgorithm::computeInsersection
(
    const label newIndex,
    const label oldIndex,
    const scalar& matchTol,
    bool output
) const
{
    bool intersects = false;

#   if USE_MPFR
    if (convexSetAlgorithm::highPrecision())
    {
        Field<mpVector> mpIntPoints;
        const mpScalar mpMatchTol(matchTol);

        // Invoke the high-precision variant
        intersects =
        (
            cellIntersection
            (
                newIndex,
                oldIndex,
                mpMatchTol,
                mpIntPoints,
                output
            )
        );

        if (intersects)
        {
            mpScalar volume;
            mpVector centre;

            convexSetVolume
            (
                newIndex,
                oldIndex,
                mpIntPoints,
                volume,
                centre,
                output
            );

            // Size-up the internal lists
            if (!output)
            {
                meshOps::sizeUpList(oldIndex, parents_);
                meshOps::sizeUpList(volume, mpWeights_);
                meshOps::sizeUpList(centre, mpCentres_);
            }
        }
    }
    else
#   endif
    {
        vectorField intPoints;

        // Invoke the conventional variant
        intersects =
        (
            cellIntersection
            (
                newIndex,
                oldIndex,
                matchTol,
                intPoints,
                output
            )
        );

        if (intersects)
        {
            scalar volume;
            vector centre;

            convexSetVolume
            (
                newIndex,
                oldIndex,
                intPoints,
                volume,
                centre,
                output
            );

            // Size-up the internal lists
            if (!output)
            {
                meshOps::sizeUpList(oldIndex, parents_);
                meshOps::sizeUpList(volume, weights_);
                meshOps::sizeUpList(centre, centres_);
            }
        }
    }

    return intersects;
}


template <class T>
bool cellSetAlgorithm::cellIntersection
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
    const cell& oldCell = this->mesh_.cells()[oldIndex];
    const UList<face>& oldFaces = this->mesh_.faces();
    const UList<label>& oldOwner = this->mesh_.faceOwner();
    const edgeList oldCellEdges = oldCell.edges(oldFaces);
    const labelList oldCellPoints = oldCell.labels(oldFaces);

    const cell& newCell = this->newCells_[newIndex];
    const UList<face>& newFaces = this->newFaces_;
    const UList<label>& newOwner = this->newOwner_;
    const edgeList newCellEdges = newCell.edges(newFaces);
    const labelList newCellPoints = newCell.labels(newFaces);

    // Alias references
    const UList<point>& newPoints = this->newPoints_;
    const UList<point>& oldPoints = this->mesh_.points();

    const Map<labelList>& modPoints = this->modPoints_;
    const List<objectMap>& pfp = this->pointsFromPoints_;

    // Track all possible intersections from here on.
    label nInts = 0;
    Map<Vector<T> > intersections;
    Vector<T> intPoint = Vector<T>::zero;
    Vector<T> checkPoint = Vector<T>::zero;

    // Topologically check for common points
    Map<labelList> commonPoints;

    forAll(oldCellPoints, pointI)
    {
        label oldPoint = oldCellPoints[pointI];
        label pIndex = findIndex(newCellPoints, oldPoint);

        if (pIndex > -1)
        {
            // If this point was modified by a collapse
            // to an edge mid-point, it can't be a common point.
            if (!modPoints.found(newCellPoints[pIndex]))
            {
                commonPoints.insert(newCellPoints[pIndex], labelList(0));

                intersections.set
                (
                    ++nInts,
                    convert<T>(oldPoints[newCellPoints[pIndex]])
                );
            }
        }
    }

    // Add points if they resulted from
    // bisections of old cell edges.
    forAll(newCellPoints, pointI)
    {
        label pIndex = newCellPoints[pointI];

        if (pIndex >= nOldPoints_ && !modPoints.found(pIndex))
        {
            // Check pointsFromPoints info
            label index = -1;

            forAll(pfp, indexI)
            {
                if (pfp[indexI].index() == pIndex)
                {
                    index = indexI;
                    break;
                }
            }

            const labelList& mObj = pfp[index].masterObjects();

            // Check if the old cell contains all master points
            bool allMaster = true;

            forAll(mObj, pointJ)
            {
                if (findIndex(oldCellPoints, mObj[pointJ]) == -1)
                {
                    allMaster = false;
                    break;
                }
            }

            if (allMaster)
            {
                commonPoints.insert(pIndex, mObj);

                intersections.set
                (
                    ++nInts,
                    convert<T>(newPoints[pIndex])
                );
            }
        }
    }

    // If all points are common, this is identical to the old cell.
    if (nInts == oldCellPoints.size())
    {
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
                this->mesh_,
                "ccSet_" + Foam::name(newIndex)
              + '<' + Foam::name(oldIndex) + '>',
                intPoints.size(),
                intPoints.size(),
                intPoints.size(),
                convert<scalar>(intPoints)
            );
        }

        return true;
    }

    // Compute the merge-tolerance
    scalar mergeTol = -1.0;

    forAll(oldCellEdges, edgeI)
    {
        const edge& e = oldCellEdges[edgeI];
        scalar edgeMag = mag(oldPoints[e[1]] - oldPoints[e[0]]);
        mergeTol = Foam::max(edgeMag, mergeTol);
    }

    forAll(newCellEdges, edgeI)
    {
        const edge& e = newCellEdges[edgeI];
        scalar edgeMag = mag(newPoints[e[1]] - newPoints[e[0]]);
        mergeTol = Foam::max(edgeMag, mergeTol);
    }

    mergeTol *= 1e-06;

    // Check for point-segment intersections
    bool pointIntersections = false;

    forAll(oldCellPoints, pointI)
    {
        label oldPoint = oldCellPoints[pointI];

        if (commonPoints.found(oldPoint))
        {
            continue;
        }

        checkPoint = convert<T>(oldPoints[oldPoint]);

        forAll(newCellEdges, edgeI)
        {
            const edge edgeToCheck = newCellEdges[edgeI];

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
                commonPoints.insert
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

        if (commonPoints.found(newPoint))
        {
            continue;
        }

        checkPoint = convert<T>(newPoints[newPoint]);

        forAll(oldCellEdges, edgeI)
        {
            const edge edgeToCheck = oldCellEdges[edgeI];

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
                commonPoints.insert
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
            const edge edgeToCheck = oldCellEdges[edgeI];

            if
            (
                commonPoints.found(edgeToCheck.start()) &&
                commonPoints.found(edgeToCheck.end())
            )
            {
                continue;
            }

            Vector<T> oldS = convert<T>(oldPoints[edgeToCheck.start()]);
            Vector<T> oldE = convert<T>(oldPoints[edgeToCheck.end()]);

            Vector<T> edgeVec = (oldS - oldE);
            edgeVec /= mag(edgeVec) + VSMALL;

            if (mag(edgeVec & planeNormal) < 0.5)
            {
                continue;
            }

            checkPoint = (0.5 * (oldS + oldE));

            if
            (
                meshOps::pointInCell
                (
                    newIndex,
                    newCell,
                    newFaces,
                    newOwner,
                    newPoints,
                    T(0.0),
                    checkPoint
                )
            )
            {
                intersections.set(++nInts, oldS);
                intersections.set(++nInts, oldE);
            }
        }

        forAll(newCellEdges, edgeI)
        {
            const edge edgeToCheck = newCellEdges[edgeI];

            if
            (
                commonPoints.found(edgeToCheck.start()) &&
                commonPoints.found(edgeToCheck.end())
            )
            {
                continue;
            }

            Vector<T> newS = convert<T>(newPoints[edgeToCheck.start()]);
            Vector<T> newE = convert<T>(newPoints[edgeToCheck.end()]);

            Vector<T> edgeVec = (newS - newE);
            edgeVec /= mag(edgeVec) + VSMALL;

            if (mag(edgeVec & planeNormal) < 0.5)
            {
                continue;
            }

            checkPoint = (0.5 * (newS + newE));

            if
            (
                meshOps::pointInCell
                (
                    oldIndex,
                    oldCell,
                    oldFaces,
                    oldOwner,
                    oldPoints,
                    T(0.0),
                    checkPoint
                )
            )
            {
                intersections.set(++nInts, newS);
                intersections.set(++nInts, newE);
            }
        }
    }
    else
    {
        // Check whether any old points are within
        // the new cell. Count these as 'intersections'.
        forAll(oldCellPoints, pointI)
        {
            label oldPoint = oldCellPoints[pointI];

            if (commonPoints.found(oldPoint))
            {
                // Only skip for shared-points.
                // If the point-position was modified
                // due to a collapse, then this point
                // could be inside the new cell.
                if (commonPoints[oldPoint].empty())
                {
                    continue;
                }
            }

            checkPoint = convert<T>(oldPoints[oldPoint]);

            if
            (
                meshOps::pointInCell
                (
                    newIndex,
                    newCell,
                    newFaces,
                    newOwner,
                    newPoints,
                    T(0.0),
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
            label newPoint = newCellPoints[pointI];

            if (commonPoints.found(newPoint))
            {
                continue;
            }

            checkPoint = convert<T>(newPoints[newPoint]);

            if
            (
                meshOps::pointInCell
                (
                    oldIndex,
                    oldCell,
                    oldFaces,
                    oldOwner,
                    oldPoints,
                    T(0.0),
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

            bool disableCheck = false;

            // Check edges topologically
            if (edgePair.first() == edgePair.second())
            {
                const edge& checkEdge = edgePair.first();

                // Check if points were modified by a collapse.
                // If both were modified, continue with check.
                if
                (
                    !modPoints.found(checkEdge.start()) &&
                    !modPoints.found(checkEdge.end())
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
                    if (!modPoints.found(cV))
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
                    line<Vector<T>, const Vector<T>&>
                    (
                        convert<T>(oldPoints[edgePair.first().start()]),
                        convert<T>(oldPoints[edgePair.first().end()])
                    ),
                    line<Vector<T>, const Vector<T>&>
                    (
                        convert<T>(newPoints[edgePair.second().start()]),
                        convert<T>(newPoints[edgePair.second().end()])
                    ),
                    matchTol,
                    intPoint
                )
            );

            if (foundIntersection)
            {
                OeToNe.setSize
                (
                    OeToNe.size() + 1,
                    edgePair
                );

                NeToOe.setSize
                (
                    NeToOe.size() + 1,
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
                this->mesh_,
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

    // Check for point-face intersections
    forAll(newCellPoints, pointI)
    {
        label newPoint = newCellPoints[pointI];

        if (commonPoints.found(newPoint))
        {
            continue;
        }

        checkPoint = convert<T>(newPoints[newPoint]);

        forAll(oldCell, faceI)
        {
            const triFace faceToCheck(oldFaces[oldCell[faceI]]);

            foundIntersection = false;

            foundIntersection =
            (
                meshOps::pointInTriFace
                (
                    triangle<Vector<T>, const Vector<T>&>
                    (
                        convert<T>(oldPoints[faceToCheck[0]]),
                        convert<T>(oldPoints[faceToCheck[1]]),
                        convert<T>(oldPoints[faceToCheck[2]])
                    ),
                    checkPoint,
                    matchTol,
                    true
                )
            );

            if (foundIntersection)
            {
                // Add to the list.
                intersections.set(++nInts, checkPoint);
            }
        }
    }

    // Loop through all old edges, and find possible
    // intersections with faces of the new cell.
    forAll(oldCellEdges, edgeI)
    {
        const edge edgeToCheck(oldCellEdges[edgeI]);

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
        const edge edgeToCheck(newCellEdges[edgeI]);

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
                }

                if
                (
                    pIter().size() &&
                    (
                        (edgeToCheck[0] == pIter.key()) ||
                        (edgeToCheck[1] == pIter.key())
                    )
                )
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

    // Points may become coincident due to round-off errors.
    // Merge as necessary.
    {
        forAllConstIter(typename Map<Vector<T> >, intersections, pI)
        {
            forAllIter(typename Map<Vector<T> >, intersections, pJ)
            {
                if (pI.key() != pJ.key())
                {
                    if (mag(pI() - pJ()) < mergeTol)
                    {
                        intersections.erase(pJ);
                    }
                }
            }
        }

        // Update actual intersections
        nInts = intersections.size();
    }

    if (nInts < 4)
    {
        // Does not intersect.
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
            this->mesh_,
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


// Compute the volume / centre of a polyhedron
// formed by a convex set of points.
template <class T>
void cellSetAlgorithm::convexSetVolume
(
    const label newCellIndex,
    const label oldCellIndex,
    const Field<Vector<T> >& cvxSet,
    T& cVolume,
    Vector<T>& cCentre,
    bool output
) const
{
    // Reset inputs
    cVolume = pTraits<T>::zero;
    cCentre = Vector<T>::zero;

    // Try the trivial case for a tetrahedron.
    // No checking for orientation here.
    if (cvxSet.size() == 4)
    {
        const Vector<T>& a = cvxSet[0];
        const Vector<T>& b = cvxSet[1];
        const Vector<T>& c = cvxSet[2];
        const Vector<T>& d = cvxSet[3];

        cCentre = ( T(0.25) * (a + b + c + d) );

        cVolume =
        (
            mag
            (
                (pTraits<T>::one / T(6.0)) *
                (
                    ((b - a) ^ (c - a)) & (d - a)
                )
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
    DynamicList<label> uniquePts;

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
                T tolerance(1e-14);

                // Compute the normal to this face
                Vector<T> n;

                meshOps::faceNormal(tmpFace, cvxSet, n);
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

                    Vector<T> rfVec = (cvxSet[l] - cvxSet[i]);
                    T dotProd = (rfVec/(mag(rfVec) + VSMALL)) & n;

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
                                Vector<T> eNorm;

                                meshOps::faceNormal(checkFace, cvxSet, eNorm);

                                T dotProd =
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
            Vector<T> n;
            meshOps::faceNormal(checkFace, cvxSet, n);

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

    // Prepare temporary connectivity
    // for volume / centre computation.
    labelList owner(testFaces.size(), 0);
    cellList cells(1, cell(identity(testFaces.size())));

    bool validVolume =
    (
        meshOps::cellCentreAndVolumeT
        (
            0,
            cvxSet,
            testFaces,
            cells,
            owner,
            cCentre,
            cVolume
        )
    );

    // Check faces for consistency
    label nValidFaces = 0;

    forAll(testFaces, faceI)
    {
        if (testFaces[faceI].size())
        {
            nValidFaces++;
        }
    }

    if (nValidFaces <= 3 || !validVolume)
    {
        meshOps::checkPointNearness(cvxSet, T(1e-20));

        // Write out cells
        writeVTK("newCell_" + Foam::name(newCellIndex), newCellIndex, 3, false);
        writeVTK("oldCell_" + Foam::name(oldCellIndex), oldCellIndex, 3, true);

        meshOps::writeVTK
        (
            this->mesh_,
            "tfSet_" + Foam::name(newCellIndex)
          + '<' + Foam::name(oldCellIndex) + '>',
            cvxSet.size(),
            cvxSet.size(),
            cvxSet.size(),
            convert<scalar>(cvxSet)
        );

        FatalErrorIn("void cellSetAlgorithm::convexSetVolume() const")
            << " Incorrect number of valid faces." << nl
            << "   newCellIndex: " << newCellIndex << nl
            << "   oldCellIndex: " << oldCellIndex << nl
            << "   nFaces: " << nValidFaces << nl
            << "   Volume: " << cVolume << nl
            << "   testFaces: " << nl << testFaces << nl
            << "   Point set: " << nl << cvxSet << nl
            << abort(FatalError);
    }

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

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
