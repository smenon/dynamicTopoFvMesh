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

#include "Time.H"
#include "meshOps.H"
#include "ListOps.H"
#include "Pstream.H"
#include "triFace.H"
#include "IOmanip.H"
#include "polyMesh.H"
#include "triPointRef.H"
#include "tetPointRef.H"
#include "labelHashSet.H"

#if USE_CGAL
#include <CGAL/Gmpz.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#endif

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{

namespace meshOps
{

// Utility method to build a hull of cells
// connected to the edge (for 2D simplical meshes)
void constructPrismHull
(
    const label eIndex,
    const UList<face>& faces,
    const UList<cell>& cells,
    const UList<label>& owner,
    const UList<label>& neighbour,
    const UList<labelList>& edgeFaces,
    labelList& hullTriFaces,
    labelList& hullCells
)
{
    labelHashSet cellSet, triFaceSet;

    // Obtain references
    const labelList& eFaces = edgeFaces[eIndex];

    // Loop through edgeFaces and add cells
    forAll(eFaces, faceI)
    {
        label c0 = owner[eFaces[faceI]];
        label c1 = neighbour[eFaces[faceI]];

        if (!cellSet.found(c0))
        {
            // Add this cell
            cellSet.insert(c0);

            // Find associated triFaces and add them too
            const cell& cC = cells[c0];

            forAll(cC, faceJ)
            {
                const face& cF = faces[cC[faceJ]];

                if ((cF.size() == 3) && !(triFaceSet.found(cC[faceJ])))
                {
                    triFaceSet.insert(cC[faceJ]);
                }
            }
        }

        if (!cellSet.found(c1) && (c1 != -1))
        {
            // Add this cell
            cellSet.insert(c1);

            // Find associated triFaces and add them too
            const cell& cC = cells[c1];

            forAll(cC, faceJ)
            {
                const face& cF = faces[cC[faceJ]];

                if ((cF.size() == 3) && !(triFaceSet.found(cC[faceJ])))
                {
                    triFaceSet.insert(cC[faceJ]);
                }
            }
        }
    }

    // Obtain lists from hashSets
    hullCells = cellSet.toc();
    hullTriFaces = triFaceSet.toc();
}


// Utility method to build a hull of cells (and faces)
// around an edge (for 3D simplical meshes)
void constructHull
(
    const label eIndex,
    const UList<face>& faces,
    const UList<edge>& edges,
    const UList<cell>& cells,
    const UList<label>& owner,
    const UList<label>& neighbour,
    const UList<labelList>& faceEdges,
    const UList<labelList>& edgeFaces,
    const UList<labelList>& edgePoints,
    labelList& hullEdges,
    labelList& hullFaces,
    labelList& hullCells,
    labelListList& ringEntities
)
{
    // [1] hullEdges is an ordered list of edge-labels around eIndex,
    //     but not connected to it.
    //      - Ordering is in the same manner as edgePoints.
    // [2] hullFaces is an ordered list of face-labels connected to eIndex.
    //      - Ordering is in the same manner as edgePoints.
    // [3] hullCells is an ordered list of cell-labels connected to eIndex.
    //      - For boundary hulls, the last cell label is -1
    // [4] ringEntities are edges and faces connected to eIndex[0] and eIndex[1]
    //      - ringEntities[0]: edges connected to eIndex[0]
    //      - ringEntities[1]: faces connected to eIndex[0]
    //      - ringEntities[2]: edges connected to eIndex[1]
    //      - ringEntities[3]: faces connected to eIndex[1]

    bool found;
    label otherPoint = -1, nextPoint = -1;

    // Obtain a reference to this edge, and its edgeFaces
    const edge& edgeToCheck = edges[eIndex];
    const labelList& eFaces = edgeFaces[eIndex];
    const labelList& hullVertices = edgePoints[eIndex];

    // Loop through all faces of this edge and add them to hullFaces
    forAll(eFaces, faceI)
    {
        const face& faceToCheck = faces[eFaces[faceI]];

        // Find the isolated point on this face,
        // and compare it with hullVertices
        meshOps::findIsolatedPoint
        (
            faceToCheck,
            edgeToCheck,
            otherPoint,
            nextPoint
        );

        found = false;

        forAll(hullVertices, indexI)
        {
            if (hullVertices[indexI] == otherPoint)
            {
                // Fill in the position of this face on the hull
                hullFaces[indexI] = eFaces[faceI];

                // Obtain edges connected to top and bottom
                // vertices of edgeToCheck
                const labelList& fEdges = faceEdges[hullFaces[indexI]];

                forAll(fEdges, edgeI)
                {
                    if
                    (
                        edges[fEdges[edgeI]]
                     == edge(edgeToCheck[0], otherPoint)
                    )
                    {
                        ringEntities[0][indexI] = fEdges[edgeI];
                    }

                    if
                    (
                        edges[fEdges[edgeI]]
                     == edge(edgeToCheck[1], otherPoint)
                    )
                    {
                        ringEntities[2][indexI] = fEdges[edgeI];
                    }
                }

                // Depending on the orientation of this face,
                // fill in hull cell indices as well
                if (nextPoint == edgeToCheck[0])
                {
                    hullCells[indexI] = owner[eFaces[faceI]];
                }
                else
                if (nextPoint == edgeToCheck[1])
                {
                    hullCells[indexI] = neighbour[eFaces[faceI]];
                }
                else
                {
                    // Something's terribly wrong
                    FatalErrorIn("void meshOps::constructHull()")
                        << nl << " Failed to construct hull. "
                        << nl << " Possibly not a tetrahedral mesh. "
                        << abort(FatalError);
                }

                if (hullCells[indexI] != -1)
                {
                    label nextI = hullVertices.fcIndex(indexI);
                    label nextHullPoint = hullVertices[nextI];
                    const cell& currCell = cells[hullCells[indexI]];

                    // Look for the ring-faces
                    forAll(currCell, faceI)
                    {
                        const face& cFace = faces[currCell[faceI]];

                        // Check if this face contains edgeToCheck[0]
                        if
                        (
                            triFace::compare
                            (
                                triFace(cFace),
                                triFace
                                (
                                    edgeToCheck[0],
                                    otherPoint,
                                    nextHullPoint
                                )
                            )
                        )
                        {
                            ringEntities[1][indexI] = currCell[faceI];
                        }

                        // Check if this face contains edgeToCheck[1]
                        if
                        (
                            triFace::compare
                            (
                                triFace(cFace),
                                triFace
                                (
                                    edgeToCheck[1],
                                    nextHullPoint,
                                    otherPoint
                                )
                            )
                        )
                        {
                            ringEntities[3][indexI] = currCell[faceI];
                        }
                    }

                    // Scan one the faces for the ring-edge
                    const labelList& rFaceEdges =
                    (
                        faceEdges[ringEntities[1][indexI]]
                    );

                    forAll(rFaceEdges, edgeI)
                    {
                        if
                        (
                            edges[rFaceEdges[edgeI]]
                         == edge(otherPoint,nextHullPoint)
                        )
                        {
                            hullEdges[indexI] = rFaceEdges[edgeI];
                            break;
                        }
                    }
                }

                // Done with this index. Break out.
                found = true;
                break;
            }
        }

        // Throw an error if the point wasn't found
        if (!found)
        {
            // Something's terribly wrong
            FatalErrorIn("void meshOps::constructHull()")
                << " Failed to construct hull. " << nl
                << " edgeFaces connectivity is inconsistent. " << nl
                << " Edge: " << eIndex << ":: " << edgeToCheck << nl
                << " edgeFaces: " << eFaces << nl
                << " edgePoints: " << hullVertices
                << abort(FatalError);
        }
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
template <class T>
bool checkPointNearness
(
    const Field<Vector<T> >& points,
    const T& magSqrTol
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

            T magSqrDist = magSqr(points[pI] - points[pJ]);

            if (magSqrDist < magSqrTol)
            {
                SeriousErrorIn
                (
                    "void meshOps::checkPointNearness"
                    "(const Field<vector<T> >&, const T&)"
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


// Parallel blocking send
void pWrite
(
    const label toID,
    const label& data
)
{
    OPstream::write
    (
        Pstream::blocking,
        toID,
        reinterpret_cast<const char*>
        (
            &data
        ),
        sizeof(label)
    );
}


// Parallel blocking receive
void pRead
(
    const label fromID,
    label& data
)
{
    IPstream::read
    (
        Pstream::blocking,
        fromID,
        reinterpret_cast<char*>
        (
            &data
        ),
        sizeof(label)
    );
}


// Parallel non-blocking send for fixed lists
template <class Type, label Size>
void pWrite
(
    const label toID,
    const FixedList<Type, Size>& data
)
{
    OPstream::write
    (
        Pstream::blocking,
        toID,
        reinterpret_cast<const char*>(&data[0]),
        Size*sizeof(Type)
    );
}


// Parallel non-blocking receive for fixed lists
template <class Type, label Size>
void pRead
(
    const label fromID,
    FixedList<Type, Size>& data
)
{
    IPstream::read
    (
        Pstream::blocking,
        fromID,
        reinterpret_cast<char*>(data.begin()),
        Size*sizeof(Type)
    );
}


// Parallel non-blocking send for lists
template <class Type>
void pWrite
(
    const label toID,
    const UList<Type>& data
)
{
    OPstream::write
    (
        Pstream::nonBlocking,
        toID,
        reinterpret_cast<const char*>(&data[0]),
        data.size()*sizeof(Type)
    );
}


// Parallel non-blocking receive for lists
template <class Type>
void pRead
(
    const label fromID,
    UList<Type>& data
)
{
    IPstream::read
    (
        Pstream::nonBlocking,
        fromID,
        reinterpret_cast<char*>(&data[0]),
        data.size()*sizeof(Type)
    );
}


// Wait for buffer transfer completion.
void waitForBuffers()
{
    if (Pstream::parRun())
    {
        OPstream::waitRequests();
        IPstream::waitRequests();
    }
}


// Select a list of elements from connectivity,
// and output to a VTK format
void writeVTK
(
    const polyMesh& mesh,
    const word& name,
    const labelList& cList,
    const label primitiveType,
    const UList<point>& meshPoints,
    const UList<edge>& edges,
    const UList<face>& faces,
    const UList<cell>& cells,
    const UList<label>& owner
)
{
    label nTotalCells = 0;
    label nPoints = 0, nCells = 0;

    // Estimate a size for points and cellPoints
    pointField points(6*cList.size());

    // Connectivity lists
    labelListList cpList(cList.size());

    // Create a map for local points
    Map<label> pointMap, reversePointMap, reverseCellMap;

    forAll(cList, cellI)
    {
        if (cList[cellI] < 0)
        {
            continue;
        }

        // Are we looking at points?
        if (primitiveType == 0)
        {
            // Size the list
            cpList[nCells].setSize(1);

            cpList[nCells] = cList[cellI];

            nTotalCells++;
        }

        // Are we looking at edges?
        if (primitiveType == 1)
        {
            // Size the list
            cpList[nCells].setSize(2);

            const edge& tEdge = edges[cList[cellI]];

            cpList[nCells][0] = tEdge[0];
            cpList[nCells][1] = tEdge[1];

            nTotalCells += 2;
        }

        // Are we looking at faces?
        if (primitiveType == 2)
        {
            const face& tFace = faces[cList[cellI]];

            if (tFace.size() == 3)
            {
                // Size the list
                cpList[nCells].setSize(3);

                // Write out in order
                cpList[nCells][0] = tFace[0];
                cpList[nCells][1] = tFace[1];
                cpList[nCells][2] = tFace[2];

                nTotalCells += 3;
            }
            else
            if (tFace.size() == 4)
            {
                // Size the list
                cpList[nCells].setSize(4);

                // Write out in order
                cpList[nCells][0] = tFace[0];
                cpList[nCells][1] = tFace[1];
                cpList[nCells][2] = tFace[2];
                cpList[nCells][3] = tFace[3];

                nTotalCells += 4;
            }
        }

        // Are we looking at cells?
        if (primitiveType == 3)
        {
            const cell& tCell = cells[cList[cellI]];

            if (tCell.size() == 4)
            {
                // Point-ordering for tetrahedra
                const face& baseFace = faces[tCell[0]];
                const face& checkFace = faces[tCell[1]];

                // Size the list
                cpList[nCells].setSize(4);

                // Get the fourth point
                label apexPoint =
                (
                    meshOps::findIsolatedPoint(baseFace, checkFace)
                );

                // Something's wrong with connectivity.
                if (apexPoint == -1)
                {
                    FatalErrorIn
                    (
                        "void writeVTK\n"
                        "(\n"
                        "    const polyMesh& mesh,\n"
                        "    const word& name,\n"
                        "    const labelList& cList,\n"
                        "    const label primitiveType,\n"
                        "    const UList<point>& points,\n"
                        "    const UList<edge>& edges,\n"
                        "    const UList<face>& faces,\n"
                        "    const UList<cell>& cells,\n"
                        "    const UList<label>& owner\n"
                        ") const\n"
                    )
                        << "Cell: " << cList[cellI]
                        << ":: " << tCell
                        << " has inconsistent connectivity."
                        << abort(FatalError);
                }

                // Write-out in order
                label ownCell = owner[tCell[0]];

                if (ownCell == cList[cellI])
                {
                    cpList[nCells][0] = baseFace[2];
                    cpList[nCells][1] = baseFace[1];
                    cpList[nCells][2] = baseFace[0];
                    cpList[nCells][3] = apexPoint;
                }
                else
                {
                    cpList[nCells][0] = baseFace[0];
                    cpList[nCells][1] = baseFace[1];
                    cpList[nCells][2] = baseFace[2];
                    cpList[nCells][3] = apexPoint;
                }

                nTotalCells += 4;
            }
            else
            if (tCell.size() == 5)
            {
                // Point-ordering for wedge cells
                label firstTriFace = -1;

                // Size the list
                cpList[nCells].setSize(6);

                // Figure out triangle faces
                forAll(tCell, faceI)
                {
                    const face& cFace = faces[tCell[faceI]];

                    if (cFace.size() == 3)
                    {
                        if (firstTriFace == -1)
                        {
                            firstTriFace = tCell[faceI];

                            // Right-handedness is assumed here.
                            // Tri-faces are always on the boundary.
                            cpList[nCells][0] = cFace[0];
                            cpList[nCells][1] = cFace[1];
                            cpList[nCells][2] = cFace[2];
                        }
                        else
                        {
                            // Detect the three other points.
                            forAll(tCell, faceJ)
                            {
                                const face& nFace = faces[tCell[faceJ]];

                                if (nFace.size() == 4)
                                {
                                    // Search for vertices on cFace
                                    // in this face.
                                    forAll(cFace, I)
                                    {
                                        label i = nFace.which(cFace[I]);

                                        if (i != -1)
                                        {
                                            label p = nFace.prevLabel(i);
                                            label n = nFace.nextLabel(i);

                                            if (p == cpList[nCells][0])
                                            {
                                                cpList[nCells][3] = cFace[I];
                                            }

                                            if (p == cpList[nCells][1])
                                            {
                                                cpList[nCells][4] = cFace[I];
                                            }

                                            if (p == cpList[nCells][2])
                                            {
                                                cpList[nCells][5] = cFace[I];
                                            }

                                            if (n == cpList[nCells][0])
                                            {
                                                cpList[nCells][3] = cFace[I];
                                            }

                                            if (n == cpList[nCells][1])
                                            {
                                                cpList[nCells][4] = cFace[I];
                                            }

                                            if (n == cpList[nCells][2])
                                            {
                                                cpList[nCells][5] = cFace[I];
                                            }
                                        }
                                    }
                                }
                            }

                            break;
                        }
                    }
                }

                nTotalCells += 6;
            }
        }

        // Renumber to local ordering
        forAll(cpList[nCells], pointI)
        {
            // Check if this point was added to the map
            if (!pointMap.found(cpList[nCells][pointI]))
            {
                // Point was not found, so add it
                points[nPoints] = meshPoints[cpList[nCells][pointI]];

                // Update the map
                pointMap.insert(cpList[nCells][pointI], nPoints);
                reversePointMap.insert(nPoints, cpList[nCells][pointI]);

                // Increment the number of points
                nPoints++;
            }

            // Renumber it.
            cpList[nCells][pointI] = pointMap[cpList[nCells][pointI]];
        }

        // Update the cell map.
        reverseCellMap.insert(nCells, cList[cellI]);

        nCells++;
    }

    // Finally write it out
    meshOps::writeVTK
    (
        mesh,
        name,
        nPoints,
        nCells,
        nTotalCells,
        points,
        cpList,
        primitiveType,
        reversePointMap,
        reverseCellMap
    );
}


// Actual routine to write out the VTK file
void writeVTK
(
    const polyMesh& mesh,
    const word& name,
    const label nPoints,
    const label nCells,
    const label nTotalCells,
    const vectorField& points,
    const labelListList& cpList,
    const label primitiveType,
    const Map<label>& reversePointMap,
    const Map<label>& reverseCellMap
)
{
    // Make the directory
    fileName dirName(mesh.time().path()/"VTK"/mesh.time().timeName());

    mkDir(dirName);

    // Open stream for output
    OFstream file(dirName/name+".vtk");

    // Write out the header
    file << "# vtk DataFile Version 2.0" << nl
         << name << ".vtk" << nl
         << "ASCII" << nl
         << "DATASET UNSTRUCTURED_GRID" << nl
         << "POINTS " << nPoints << " double" << nl;

    for (label i = 0; i < nPoints; i++)
    {
        file << setprecision(15)
             << points[i].x() << ' '
             << points[i].y() << ' '
             << points[i].z() << ' '
             << nl;
    }

    file << "CELLS " << nCells << " " << nTotalCells + nCells << endl;

    if (cpList.size())
    {
        forAll(cpList, i)
        {
            if (cpList[i].size())
            {
                file << cpList[i].size() << ' ';

                forAll(cpList[i], j)
                {
                    file << cpList[i][j] << ' ';
                }

                file << nl;
            }
        }
    }
    else
    {
        // List of points
        for (label i = 0; i < nPoints; i++)
        {
            file << 1 << ' ' << i << nl;
        }
    }

    file << "CELL_TYPES " << nCells << endl;

    if (cpList.size())
    {
        forAll(cpList, i)
        {
            if (cpList[i].size() == 1)
            {
                // Vertex
                file << "1" << nl;
            }

            if (cpList[i].size() == 2)
            {
                // Edge
                file << "3" << nl;
            }

            if (cpList[i].size() == 3)
            {
                // Triangle face
                file << "5" << nl;
            }

            if
            (
                (cpList[i].size() == 4) &&
                (primitiveType == 2)
            )
            {
                // Quad face
                file << "9" << nl;
            }

            if
            (
                (cpList[i].size() == 4) &&
                (primitiveType == 3)
            )
            {
                // Tetrahedron
                file << "10" << nl;
            }

            if (cpList[i].size() == 6)
            {
                // Wedge
                file << "13" << nl;
            }
        }
    }
    else
    {
        // List of points
        for (label i = 0; i < nPoints; i++)
        {
            // Vertex
            file << '1' << nl;
        }
    }

    // Write out indices for visualization.
    if (reverseCellMap.size())
    {
        file << "CELL_DATA " << nCells << endl;

        file << "FIELD CellFields 1" << endl;

        file << "CellIds 1 " << nCells << " int" << endl;

        for (label i = 0; i < nCells; i++)
        {
            file << reverseCellMap[i] << ' ';
        }

        file << endl;
    }

    // Write out indices for visualization.
    if (reversePointMap.size())
    {
        file << "POINT_DATA " << nPoints << endl;

        file << "FIELD PointFields 1" << endl;

        file << "PointIds 1 " << nPoints << " int" << endl;

        for (label i = 0; i < nPoints; i++)
        {
            file << reversePointMap[i] << ' ';
        }

        file << endl;
    }
}


} // End namespace meshOps


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
convexSetAlgorithm::convexSetAlgorithm
(
    const polyMesh& mesh,
    const UList<point>& newPoints,
    const UList<edge>& newEdges,
    const UList<face>& newFaces,
    const UList<cell>& newCells,
    const UList<label>& newOwner,
    const UList<label>& newNeighbour,
    const List<objectMap>& pointsFromPoints,
    const Map<labelList>& modPoints
)
:
    twoDMesh_(mesh.nGeometricD() == 2),
    nOldPoints_(mesh.nPoints()),
    mesh_(mesh),
    newPoints_(newPoints),
    newEdges_(newEdges),
    newFaces_(newFaces),
    newCells_(newCells),
    newOwner_(newOwner),
    newNeighbour_(newNeighbour),
    pointsFromPoints_(pointsFromPoints),
    modPoints_(modPoints),
    highPrecision_(false)
{}


faceSetAlgorithm::faceSetAlgorithm
(
    const polyMesh& mesh,
    const UList<point>& newPoints,
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


cellSetAlgorithm::cellSetAlgorithm
(
    const polyMesh& mesh,
    const UList<point>& newPoints,
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

// Output an entity as a VTK file
void convexSetAlgorithm::writeVTK
(
    const word& name,
    const label entity,
    const label primitiveType,
    const bool useOldConnectivity
) const
{
    writeVTK
    (
        name,
        labelList(1, entity),
        primitiveType,
        useOldConnectivity
    );
}


// Output a list of entities as a VTK file
void convexSetAlgorithm::writeVTK
(
    const word& name,
    const labelList& cList,
    const label primitiveType,
    const bool useOldConnectivity
) const
{
    if (useOldConnectivity)
    {
        const polyMesh& mesh = this->mesh_;

        meshOps::writeVTK
        (
            mesh,
            name,
            cList,
            primitiveType,
            mesh.points(),
            mesh.edges(),
            mesh.faces(),
            mesh.cells(),
            mesh.faceOwner()
        );
    }
    else
    {
        meshOps::writeVTK
        (
            this->mesh_,
            name,
            cList,
            primitiveType,
            newPoints_,
            newEdges_,
            newFaces_,
            newCells_,
            newOwner_
        );
    }
}


template <class T>
bool faceSetAlgorithm::faceIntersection
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

    // Fetch face references for each mesh
    const face& newFace = this->newFaces_[newIndex];
    const face& oldFace = this->mesh_.faces()[oldIndex];

    // Alias references
    const UList<point>& newPoints = this->newPoints_;
    const UList<point>& oldPoints = this->mesh_.points();

    const Map<labelList>& modPoints = this->modPoints_;
    const List<objectMap>& pfp = this->pointsFromPoints_;

    // Obtain face centre and projection normal
    Vector<T> xf, nf;

    meshOps::faceCentre(newFace, newPoints, xf);
    meshOps::faceNormal(newFace, newPoints, nf);

    nf /= mag(nf) + VSMALL;

    // Track all possible intersections from here on.
    label nInts = 0;
    Map<Vector<T> > intersections;
    Vector<T> intPoint = Vector<T>::zero;

    // Topologically check for common points,
    // and project uniques ones on to the face plane
    Map<label> projPoints;
    Map<labelList> commonPoints;
    Field<Vector<T> > projections(oldFace.size(), Vector<T>::zero);

    forAll(oldFace, pointI)
    {
        label oldPoint = oldFace[pointI];
        label pIndex = findIndex(newFace, oldPoint);

        Vector<T> r = convert<T>(oldPoints[oldPoint]);

        if (pIndex == -1)
        {
            // Project this point on to the newFace plane.
            projections[pointI] = xf + ((r - xf) - ((r - xf) & nf)*nf);

            projPoints.insert(oldPoint, pointI);
        }
        else
        {
            // If this point was modified by a collapse
            // to an edge mid-point, it can't be a common point.
            if (modPoints.found(newFace[pIndex]))
            {
                // Project this point on to the newFace plane.
                projections[pointI] = xf + ((r - xf) - ((r - xf) & nf)*nf);

                projPoints.insert(oldPoint, pointI);
            }
            else
            {
                commonPoints.insert(newFace[pIndex], labelList(0));

                projections[pointI] = r;

                intersections.set(++nInts, r);
            }
        }
    }

    // Add points if they resulted from
    // bisections of old face edges.
    forAll(newFace, pointI)
    {
        label pIndex = newFace[pointI];

        if (pIndex >= nOldPoints_)
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

            // Check if the old face contains all master points
            bool allMaster = true;

            forAll(mObj, pointJ)
            {
                if (findIndex(oldFace, mObj[pointJ]) == -1)
                {
                    allMaster = false;
                    break;
                }
            }

            if (allMaster)
            {
                commonPoints.insert(newFace[pointI], mObj);

                intersections.set
                (
                    ++nInts,
                    convert<T>(newPoints[newFace[pointI]])
                );
            }
        }
    }

    // If all points are common, this is identical to the old face.
    if (nInts == oldFace.size())
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
            // Check for point nearness
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

    // Check for point-segment intersections
    bool pointIntersections = false;

    forAll(oldFace, pointI)
    {
        label oldPoint = oldFace[pointI];

        if (commonPoints.found(oldPoint))
        {
            continue;
        }

        const Vector<T>& checkPoint = projections[pointI];

        // Loop through all new edges, and find possible intersections
        // with (projections of) old face points,
        forAll(newFace, pointJ)
        {
            edge newEdge = newFace.faceEdge(pointJ);

            if
            (
                meshOps::pointSegmentIntersection
                (
                    line<Vector<T>, const Vector<T>&>
                    (
                        convert<T>(newPoints[newEdge.start()]),
                        convert<T>(newPoints[newEdge.end()])
                    ),
                    checkPoint,
                    matchTol
                )
            )
            {
                commonPoints.insert
                (
                    oldPoint,
                    labelList(newEdge)
                );

                intersections.set(++nInts, checkPoint);

                pointIntersections = true;
            }
        }
    }

    forAll(newFace, pointI)
    {
        label newPoint = newFace[pointI];

        if (commonPoints.found(newPoint))
        {
            continue;
        }

        // Check if this point was modified during a collapse.
        // We might be able to avoid a calculation if it is.
        if (modPoints.found(newPoint))
        {
            // Fetch the two points that it originated from
            const labelList& mP = modPoints[newPoint];

            bool foundPoint = false;

            forAll(oldFace, pointJ)
            {
                const edge edgeToCheck = oldFace.faceEdge(pointJ);

                if (edgeToCheck == edge(mP[0], mP[1]))
                {
                    commonPoints.insert
                    (
                        newPoint,
                        labelList(edgeToCheck)
                    );

                    intersections.set
                    (
                        ++nInts,
                        convert<T>(newPoints[newPoint])
                    );

                    foundPoint = true;
                    pointIntersections = true;

                    break;
                }
            }

            if (foundPoint)
            {
                continue;
            }
        }

        const Vector<T> checkPoint = convert<T>(newPoints[newPoint]);

        forAll(oldFace, pointJ)
        {
            label nextJ = oldFace.fcIndex(pointJ);
            edge oldEdge = oldFace.faceEdge(pointJ);

            if
            (
                meshOps::pointSegmentIntersection
                (
                    line<Vector<T>, const Vector<T>&>
                    (
                        projections[pointJ],
                        projections[nextJ]
                    ),
                    checkPoint,
                    matchTol
                )
            )
            {
                commonPoints.insert
                (
                    newPoint,
                    labelList(oldEdge)
                );

                intersections.set(++nInts, checkPoint);

                pointIntersections = true;
            }
        }
    }

    if (pointIntersections && output)
    {
        Info << "Point Intersections exist: " << nl
             << " newFaceIndex: " << newIndex
             << " oldFaceIndex: " << oldIndex
             << endl;
    }

    if (oldFace.size() == 3 && newFace.size() == 3)
    {
        // Perform tests specific to triangular faces

        // Check whether any old projections are within
        // the new face. Count these as 'intersections'.
        forAll(oldFace, pointI)
        {
            label oldPoint = oldFace[pointI];

            if (commonPoints.found(oldPoint))
            {
                // Only skip for shared-points.
                // If the point-position was modified
                // due to a collapse, then this point
                // could be inside the new face.
                if (commonPoints[oldPoint].empty())
                {
                    continue;
                }
                else
                {
                    // Fetch master objects
                    const labelList& mObj = commonPoints[oldPoint];

                    // Check if the new face
                    // contains all master points
                    bool allMaster = true;

                    forAll(mObj, pointJ)
                    {
                        if (findIndex(newFace, mObj[pointJ]) == -1)
                        {
                            allMaster = false;
                            break;
                        }
                    }

                    if (allMaster)
                    {
                        continue;
                    }
                }
            }

            const Vector<T>& checkPoint = projections[pointI];

            if
            (
                meshOps::pointInTriFace
                (
                    triangle<Vector<T>, const Vector<T>&>
                    (
                        convert<T>(newPoints[newFace[0]]),
                        convert<T>(newPoints[newFace[1]]),
                        convert<T>(newPoints[newFace[2]])
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
        forAll(newFace, pointI)
        {
            label newPoint = newFace[pointI];

            if (commonPoints.found(newPoint))
            {
                continue;
            }

            const Vector<T> checkPoint = convert<T>(newPoints[newPoint]);

            if
            (
                meshOps::pointInTriFace
                (
                    triangle<Vector<T>, const Vector<T>&>
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
        forAll(newFace, pointI)
        {
            edge newEdge = newFace.faceEdge(pointI);
            label nextLabel = newFace.nextLabel(pointI);

            forAll(oldFace, pointJ)
            {
                label nextJ = oldFace.fcIndex(pointJ);
                edge oldEdge = oldFace.faceEdge(pointJ);

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

                bool foundIntersection = false;

                foundIntersection =
                (
                    meshOps::segmentSegmentIntersection
                    (
                        line<Vector<T>, const Vector<T>&>
                        (
                            projections[pointJ],
                            projections[nextJ]
                        ),
                        line<Vector<T>, const Vector<T>&>
                        (
                            convert<T>(newPoints[newFace[pointI]]),
                            convert<T>(newPoints[nextLabel])
                        ),
                        matchTol,
                        intPoint
                    )
                );

                if (foundIntersection)
                {
                    intersections.set(++nInts, intPoint);
                }
            }
        }
    }
    else
    if (oldFace.size() == 4 && newFace.size() == 4)
    {
        // Perform tests specific to quad faces
        notImplemented
        (
            "\n\n"
            "bool faceSetAlgorithm::faceIntersection\n"
            "(\n"
            "    const label,\n"
            "    const label,\n"
            "    const T&,\n"
            "    Field<Vector<T> >&\n"
            ") const\n"
        );
    }
    else
    {
        FatalErrorIn
        (
            "\n\n"
            "bool faceSetAlgorithm::faceIntersection\n"
            "(\n"
            "    const label,\n"
            "    const label,\n"
            "    const T&,\n"
            "    Field<Vector<T> >&\n"
            ") const\n"
        )
            << " Invalid face pair: " << nl
            << " Old face: " << oldIndex << "::" << oldFace << nl
            << " New face: " << newIndex << "::" << newFace << nl
            << abort(FatalError);
    }

    // Copy intersections
    intPoints.setSize(nInts, Vector<T>::zero);

    nInts = 0;

    forAllConstIter(typename Map<Vector<T> >, intersections, pI)
    {
        intPoints[nInts++] = pI();
    }

    // Check for concurrent points.
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
    if (nInts >= 3)
    {
        return true;
    }

    // Does not intersect
    return false;
}


// Compute the area / centre of a polygon
// formed by a convex set of points.
template <class T>
void faceSetAlgorithm::convexSetArea
(
    const label newFaceIndex,
    const label oldFaceIndex,
    const Field<Vector<T> >& cvxSet,
    const Vector<T>& refNorm,
    T& fArea,
    Vector<T>& fCentre,
    bool output
) const
{
    // Reset inputs
    fArea = pTraits<T>::zero;
    fCentre = Vector<T>::zero;

    // Try the trivial case for a triangle.
    if (cvxSet.size() == 3)
    {
        const Vector<T>& a = cvxSet[0];
        const Vector<T>& b = cvxSet[1];
        const Vector<T>& c = cvxSet[2];

        fArea = mag(0.5 * ((b - a)^(c - a)));
        fCentre = (pTraits<T>::one / T(3.0)) * (a + b + c);

        if (output)
        {
            Info << " newFaceIndex: " << newFaceIndex
                 << " oldFaceIndex: " << oldFaceIndex << nl
                 << " Area: " << fArea << nl
                 << " Centre: " << fCentre << nl
                 << endl;
        }

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
            T tolerance(1e-14);

            // Compute the normal to this edge
            Vector<T> n;

            n = ((cvxSet[tmpEdge.end()] - cvxSet[tmpEdge.start()]) ^ refNorm);
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

                Vector<T> rfVec = (cvxSet[k] - cvxSet[i]);
                T dotProd = (rfVec/(mag(rfVec) + VSMALL)) & n;

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
        WarningIn
        (
            "\n"
            "void meshOps::convexSetArea\n"
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
            << " newFaceIndex: " << newFaceIndex
            << " oldFaceIndex: " << oldFaceIndex
            << " nPoints: " << cvxSet.size() << nl
            << " nEdges: " << testEdges.size() << nl
            << " Edge list: " << testEdges << nl
            << " Set: " << cvxSet << nl
            << endl;
    }

    // Find an approximate face-centroid
    T sumA = 0.0;
    Vector<T> sumAc = Vector<T>::zero;
    Vector<T> xC = average(cvxSet);

    forAll(testEdges, edgeI)
    {
        const edge& e = testEdges[edgeI];

        Vector<T> c = cvxSet[e[0]] + cvxSet[e[1]] + xC;
        T a = mag((cvxSet[e[1]] - cvxSet[e[0]]) ^ (xC - cvxSet[e[0]]));

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
    const edgeList newCellEdges = newCell.edges(this->newFaces_);
    const labelList newCellPoints = newCell.labels(this->newFaces_);

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

        if (pIndex >= nOldPoints_)
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
                commonPoints.insert(newCellPoints[pointI], mObj);

                intersections.set
                (
                    ++nInts,
                    convert<T>(newPoints[newCellPoints[pointI]])
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

        // Check if this point was modified during a collapse.
        // We might be able to avoid a calculation if it is.
        if (modPoints.found(newPoint))
        {
            // Fetch the two points that it originated from
            const labelList& mP = modPoints[newPoint];

            bool foundPoint = false;

            forAll(oldCellEdges, edgeI)
            {
                const edge edgeToCheck = oldCellEdges[edgeI];

                if (edgeToCheck == edge(mP[0], mP[1]))
                {
                    commonPoints.insert
                    (
                        newPoint,
                        labelList(edgeToCheck)
                    );

                    intersections.set
                    (
                        ++nInts,
                        convert<T>(newPoints[newPoint])
                    );

                    foundPoint = true;
                    pointIntersections = true;

                    break;
                }
            }

            if (foundPoint)
            {
                continue;
            }
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

            // Form an edge-pair
            Pair<edge> edgePair(oldCellEdges[edgeI], newCellEdges[edgeJ]);

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

#ifdef USE_CGAL
bool cellSetAlgorithm::cellIntersection_CGAL
(
    const label newIndex,
    const label oldIndex,
    scalar& cVolume,
    vector& cCentre,
    bool output
) const
{
    // Fetch references for each mesh
    const cell& oldCell = this->mesh_.cells()[oldIndex];
    const labelList oldCellPoints = oldCell.labels(this->mesh_.faces());

    const cell& newCell = this->newCells_[newIndex];
    const labelList newCellPoints = newCell.labels(this->newFaces_);

    // Alias references
    const UList<point>& newPoints = this->newPoints_;
    const UList<point>& oldPoints = this->mesh_.points();

    // Define typedefs for convenience
    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

    // Create two Polyhedron_3 objects
    CGAL::Polyhedron_3<Kernel> Po, Pn;

    // If both cells are tets, generate them.
    if (oldCell.size() == 4 && newCell.size() == 4)
    {
        // Fetch point-positions
        const point& Oa = oldPoints[oldCellPoints[0]];
        const point& Ob = oldPoints[oldCellPoints[1]];
        const point& Oc = oldPoints[oldCellPoints[2]];
        const point& Od = oldPoints[oldCellPoints[3]];

        Po.make_tetrahedron
        (
            Kernel::Point_3(Oa.x(), Oa.y(), Oa.z()),
            Kernel::Point_3(Ob.x(), Ob.y(), Ob.z()),
            Kernel::Point_3(Oc.x(), Oc.y(), Oc.z()),
            Kernel::Point_3(Od.x(), Od.y(), Od.z())
        );

        // Fetch point-positions
        const point& Na = newPoints[newCellPoints[0]];
        const point& Nb = newPoints[newCellPoints[1]];
        const point& Nc = newPoints[newCellPoints[2]];
        const point& Nd = newPoints[newCellPoints[3]];

        Pn.make_tetrahedron
        (
            Kernel::Point_3(Na.x(), Na.y(), Na.z()),
            Kernel::Point_3(Nb.x(), Nb.y(), Nb.z()),
            Kernel::Point_3(Nc.x(), Nc.y(), Nc.z()),
            Kernel::Point_3(Nd.x(), Nd.y(), Nd.z())
        );
    }

    // Create Nef_Polyhedron equivalents
    CGAL::Nef_polyhedron_3<Kernel> nPo(Po), nPn(Pn);

    // Compute the intersection
    nPo *= nPn;

    if (nPo.is_simple())
    {
        // Convert back to polyhedron
        nPo.convert_to_polyhedron(Po);

        label nFaces = Po.size_of_facets();
        label nPoints = Po.size_of_vertices();

        if (nPoints >= 4 && nFaces >= 4)
        {
            // Prepare temporary connectivity
            // for volume / centre computation.
            labelList owner(nFaces, 0);
            faceList polyFaces(nFaces);
            vectorField cvxSet(nPoints, vector::zero);
            cellList cells(1, cell(identity(nFaces)));

            CGAL::Polyhedron_3<Kernel>::Vertex_iterator v;

            nPoints = 0;

            for (v = Po.vertices_begin(); v != Po.vertices_end(); ++v)
            {
                cvxSet[nPoints].x() = to_double(v->point().x());
                cvxSet[nPoints].y() = to_double(v->point().y());
                cvxSet[nPoints].z() = to_double(v->point().z());

                nPoints++;
            }

            CGAL::Polyhedron_3<Kernel>::Facet_iterator f;
            CGAL::Polyhedron_3<Kernel>::Halfedge_around_facet_circulator hfc;

            nFaces = 0;

            for (f = Po.facets_begin(); f != Po.facets_end(); ++f)
            {
                hfc = f->facet_begin();

                label j = 0;
                face& pF = polyFaces[nFaces];

                // Size the face
                pF.setSize(CGAL::circulator_size(hfc));

                do
                {
                    pF[j++] =
                    (
                        std::distance
                        (
                            Po.vertices_begin(),
                            hfc->vertex()
                        )
                    );

                } while (++hfc != f->facet_begin());

                nFaces++;
            }

            meshOps::cellCentreAndVolume
            (
                0,
                cvxSet,
                polyFaces,
                cells,
                owner,
                cCentre,
                cVolume
            );

            if (output)
            {
                std::cout << "   newIndex: " << newIndex << nl
                          << "   oldIndex: " << oldIndex << nl
                          << "   Volume: " << cVolume << nl
                          << "   Polyhedron info: " << nl
                          << Po << std::endl;
            }

            // Found an intersection
            return true;
        }
    }

    // Failed to find an intersection
    return false;
}
#endif

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

    meshOps::cellCentreAndVolume
    (
        0,
        cvxSet,
        testFaces,
        cells,
        owner,
        cCentre,
        cVolume
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

    if (nValidFaces <= 3)
    {
        meshOps::checkPointNearness(cvxSet, T(1e-20));

        // Write out cells
        writeVTK("newCell_" + Foam::name(newCellIndex), newCellIndex, 3, false);
        writeVTK("oldCell_" + Foam::name(oldCellIndex),oldCellIndex, 3, true);

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
