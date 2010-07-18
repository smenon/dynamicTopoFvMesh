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
    meshOps

Description
    Various utility functions that perform mesh-related operations.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*----------------------------------------------------------------------------*/

#include "meshOps.H"
#include "ListOps.H"
#include "triPointRef.H"
#include "tetPointRef.H"
#include "labelHashSet.H"

namespace Foam
{

namespace meshOps
{

// Compute the area / centre of a polygon
// formed by a convex set of points.
void convexSetArea
(
    const label newFaceIndex,
    const label oldFaceIndex,
    const vectorField& cvxSet,
    const vector& refNorm,
    scalar& fArea,
    vector& fCentre,
    bool output
)
{
    // Reset inputs
    fArea = 0.0;
    fCentre = vector::zero;

    // Try the trivial case for a triangle.
    if (cvxSet.size() == 3)
    {
        triPointRef tpr(cvxSet[0], cvxSet[1], cvxSet[2]);

        fArea = tpr.mag();
        fCentre = tpr.centre();

        return;
    }

    // Track edges
    label nEdges = 0;
    edgeList testEdges(0);

    // Loop through all points, and build edges with every
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

            // Define the edge
            edge tmpEdge(i, j);

            // If this is an existing edge, skip it.
            bool foundExisting = false;

            forAll(testEdges, edgeI)
            {
                if (testEdges[edgeI] == tmpEdge)
                {
                    foundExisting = true;
                    break;
                }
            }

            if (foundExisting)
            {
                continue;
            }

            // Specify a tolerance for collinearity
            scalar tolerance = 1e-14;

            // Compute the normal to this edge
            vector n = (tmpEdge.vec(cvxSet) ^ refNorm);

            n /= mag(n) + VSMALL;

            label curEdgeSign = 0;
            bool foundInternalEdge = false;

            // Quick-reject test:
            //   Check all other points in the set,
            //   and decide if all points lie on one side.
            forAll(cvxSet, k)
            {
                // Skip duplicates.
                if (tmpEdge[0] == k || tmpEdge[1] == k)
                {
                    continue;
                }

                vector rfVec = (cvxSet[k] - cvxSet[i]);
                scalar dotProd = (rfVec/(mag(rfVec) + VSMALL)) & n;

                // Skip nearly collinear points.
                if (mag(dotProd) < tolerance)
                {
                    continue;
                }

                // Obtain the sign of this point.
                label eSign = Foam::sign(dotProd);

                // Update the current sign if necessary.
                if (curEdgeSign == 0)
                {
                    curEdgeSign = eSign;
                }
                else
                if (curEdgeSign != eSign)
                {
                    // Interior edge. Bail out.
                    foundInternalEdge = true;
                    break;
                }
            }

            if (foundInternalEdge)
            {
                continue;
            }

            // Looks like we found an edge on the boundary.
            // Check its sign to ensure that it points outward.
            if (curEdgeSign == 1)
            {
                n *= -1.0;
                tmpEdge = tmpEdge.reverseEdge();
            }

            // Add to the list of edges.
            testEdges.setSize(++nEdges, tmpEdge);
        }
    }

    // Sanity check - do points match edges?
    if (testEdges.size() != cvxSet.size())
    {
        FatalErrorIn
        (
            "\n"
            "void convexSetArea\n"
            "(\n"
            "    const label newFaceIndex,\n"
            "    const label oldFaceIndex,\n"
            "    const vectorField& cvxSet,\n"
            "    const vector& refNorm,\n"
            "    scalar& fArea,\n"
            "    vector& fCentre,\n"
            "    bool output\n"
            ") const"
        )
            << " Points do not match edges. " << nl
            << " nPoints: " << cvxSet.size() << nl
            << " nEdges: " << testEdges.size() << nl
            << " Edge list: " << testEdges
            << abort(FatalError);
    }

    // Find an approximate face-centroid
    scalar sumA = 0.0;
    vector sumAc = vector::zero;
    vector xC = average(cvxSet);

    forAll(testEdges, edgeI)
    {
        const edge& e = testEdges[edgeI];

        vector c = cvxSet[e[0]] + cvxSet[e[1]] + xC;
        scalar a = mag(e.vec(cvxSet) ^ (xC - cvxSet[e[0]]));

        sumA += a;
        sumAc += a*c;
    }

    fCentre = (1.0/3.0)*sumAc/(sumA + VSMALL);
    fArea = 0.5*sumA;

    if (output)
    {
        Info << " newFaceIndex: " << newFaceIndex
             << " oldFaceIndex: " << oldFaceIndex << nl
             << " Edges: " << testEdges << nl
             << " Area: " << fArea << nl
             << " Centre: " << fCentre << nl
             << endl;
    }
}


// Compute the volume / centre of a polyhedron
// formed by a convex set of points.
void convexSetVolume
(
    const label newCellIndex,
    const label oldCellIndex,
    const vectorField& cvxSet,
    scalar& cVolume,
    vector& cCentre,
    bool output
)
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

            if (checkFace.size() >= testFace.size())
            {
                label nCommon = 0;

                uniquePts.clear();

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
                        meshOps::insertPointLabels
                        (
                            n,
                            cvxSet,
                            uniquePts,
                            checkFace
                        );
                    }
                }
            }
        }
    }

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


// Given a set of points and edges, find the shortest path
// between the start and end point, using Dijkstra's algorithm.
//  - Takes a Map of points and edges that use those points.
//  - Edge weights are currently edge-lengths, but can easily be adapted.
//  - Returns true if the endPoint was found by the algorithm.
//  - The Map 'pi' returns a preceding point for every point in 'points'.
//
//  Algorithm is inspired by:
//    Renaud Waldura
//    Dijkstra's Shortest Path Algorithm in Java
//    http://renaud.waldura.com/
bool Dijkstra
(
    const Map<point>& points,
    const Map<edge>& edges,
    const label startPoint,
    const label endPoint,
    Map<label>& pi
)
{
    bool foundEndPoint = false;

    // Set of unvisited (Q) / visited (S) points and distances (d)
    labelHashSet Q, S;
    Map<scalar> d;

    // Initialize distances to large values
    forAllConstIter(Map<point>, points, pIter)
    {
        d.insert(pIter.key(), GREAT);
    }

    // Invert edges to make a local pointEdges list
    Map<labelList> localPointEdges;

    forAllConstIter(Map<edge>, edges, eIter)
    {
        const edge& edgeToCheck = eIter();

        forAll(edgeToCheck, pointI)
        {
            if (!localPointEdges.found(edgeToCheck[pointI]))
            {
                localPointEdges.insert(edgeToCheck[pointI], labelList(0));
            }

            meshOps::sizeUpList
            (
                eIter.key(),
                localPointEdges[edgeToCheck[pointI]]
            );
        }
    }

    // Mark the startPoint as having the smallest distance
    d[startPoint] = 0.0;

    // Add the startPoint to the list of unvisited points
    Q.insert(startPoint);

    while (Q.size())
    {
        // Step 1: Find the node with the smallest distance from the start.
        labelHashSet::iterator smallest = Q.begin();

        for
        (
            labelHashSet::iterator iter = ++Q.begin();
            iter != Q.end();
            iter++
        )
        {
            if (d[iter.key()] < d[smallest.key()])
            {
                smallest = iter;
            }
        }

        label pointIndex = smallest.key();
        scalar smallestDistance = d[pointIndex];

        // Move to the visited points list
        S.insert(pointIndex);
        Q.erase(pointIndex);

        // Step 2: Build a list of points adjacent to pointIndex
        //         but not in the visited list
        DynamicList<label> adjacentPoints(10);

        const labelList& pEdges = localPointEdges[pointIndex];

        forAll(pEdges, edgeI)
        {
            const edge& edgeToCheck = edges[pEdges[edgeI]];

            label otherPoint = edgeToCheck.otherVertex(pointIndex);

            if (!S.found(otherPoint))
            {
                adjacentPoints.append(otherPoint);
            }
        }

        // Step 3: Perform distance-based checks for adjacent points
        forAll(adjacentPoints, pointI)
        {
            label adjPoint = adjacentPoints[pointI];

            scalar distance =
            (
                mag(points[adjPoint] - points[pointIndex])
              + smallestDistance
            );

            // Check if the end-point has been touched.
            if (adjPoint == endPoint)
            {
                foundEndPoint = true;
            }

            if (distance < d[adjPoint])
            {
                // Update to the shorter distance
                d[adjPoint] = distance;

                // Update the predecessor
                if (pi.found(adjPoint))
                {
                    pi[adjPoint] = pointIndex;
                }
                else
                {
                    pi.insert(adjPoint, pointIndex);
                }

                // Add to the list of unvisited points
                Q.insert(adjPoint);
            }
        }
    }

    /*
    // Write out the path
    if (debug > 3)
    {
        if (foundEndPoint)
        {
            DynamicList<label> pathNodes(50);

            label currentPoint = endPoint;

            while (currentPoint != startPoint)
            {
                pathNodes.append(currentPoint);

                currentPoint = pi[currentPoint];
            }

            pathNodes.append(startPoint);

            pathNodes.shrink();

            writeVTK
            (
                "DijkstraPath_"
              + Foam::name(startPoint)
              + '_'
              + Foam::name(endPoint),
                pathNodes,
                0
            );
        }
    }
    */

    return foundEndPoint;
}


// Utility method to check for concurrent points.
bool checkPointNearness
(
    const pointField& points,
    const scalar magSqrTol
)
{
    forAll(points, pI)
    {
        forAll(points, pJ)
        {
            if (pI == pJ)
            {
                continue;
            }

            scalar magSqrDist = magSqr(points[pI] - points[pJ]);

            if (magSqrDist < magSqrTol)
            {
                SeriousErrorIn
                (
                    "void meshOps::checkPointNearness"
                    "(const pointField&, const scalar)"
                )
                    << " Found concurrent points: " << nl
                    << " pI: " << pI << " pJ: " << pJ << nl
                    << " point: " << points[pI] << nl
                    << " Points: " << points << nl
                    << " distance: " << magSqrDist << nl
                    << " tolerance: " << magSqrTol
                    << endl;

                return true;
            }
        }
    }

    // No errors detected.
    return false;
}

} // End namespace meshOps


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
