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

\*---------------------------------------------------------------------------*/

#include "objectMap.H"
#include "dynamicTopoFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Reorder points after a topology change
void dynamicTopoFvMesh::reOrderPoints
(
    pointField& points
)
{
    // *** Point renumbering *** //
    // If points were deleted during topology change, the numerical order ceases
    // to be continuous. Loop through all points and renumber sequentially.

    // Allocate for the mapping information
    pointMap_.setSize(this->nPoints_, -1);

    label pointInOrder = 0;

    addedPointRenumbering_.clear();

    for (label pointI = 0; pointI < nOldPoints_; pointI++)
    {
        // Check if this is a deleted point
        if (reversePointMap_[pointI] == -1)
        {
            continue;
        }

        // Update the point info
        points[pointInOrder] = this->points_[pointI];

        // Update maps
        pointMap_[pointInOrder]  = pointI;
        reversePointMap_[pointI] = pointInOrder;

        // Update the counter
        pointInOrder++;
    }

    for (label pointI = nOldPoints_; pointI < points_.size(); pointI++)
    {
        // Was this point removed after addition?
        if (deletedPoints_.found(pointI))
        {
            continue;
        }

        // Update the point info
        points[pointInOrder] = points_[pointI];

        // Put inserted points in a seperate hashSet
        addedPointRenumbering_.insert(pointI, pointInOrder);

        // Update the counter
        pointInOrder++;
    }

    // Final check to ensure everything went okay
    if (pointInOrder != nPoints_)
    {
        FatalErrorIn("dynamicTopoFvMesh::reOrderPoints()") << nl
                << " Algorithm did not visit every point in the mesh."
                << " Something's messed up." << nl
                << abort(FatalError);
    }

    // Update the local copy
    points_.setSize(nPoints_);

    points_ = points;

    // Clear the deleted entity map
    deletedPoints_.clear();
}

// Reorder edges after a topology change
void dynamicTopoFvMesh::reOrderEdges
(
    edgeList& edges,
    labelListList& edgeFaces,
    labelListList& faceEdges
)
{
    // *** Edge renumbering *** //
    // If edges were deleted during topology change, the numerical order ceases
    // to be continuous. Edges are added to respective internal/boundary patches

    // Allocate for mapping information
    edgeMap_.setSize(nEdges_, -1);

    label edgeInOrder = 0, allEdges = edges_.size();
    edgeList oldEdges(allEdges);
    labelListList oldEdgeFaces(allEdges);
    labelListList oldEdgePoints(allEdges);

    addedEdgeRenumbering_.clear();
    Map<label> addedEdgeReverseRenumbering;

    // Transfer old edge-based lists, and clear them
    forAll(edges_, edgeI)
    {
        oldEdges[edgeI] = edges_[edgeI];
        oldEdgeFaces[edgeI].transfer(edgeFaces_[edgeI]);
    }

    edges_.setSize(nEdges_); edgeFaces_.setSize(nEdges_);

    if (!twoDMesh_)
    {
        forAll(edgePoints_, edgeI)
        {
            oldEdgePoints[edgeI].transfer(edgePoints_[edgeI]);
        }

        edgePoints_.setSize(nEdges_);
    }

    // Keep track of inserted boundary edge indices
    labelList boundaryPatchIndices(edgePatchStarts_);

    // Loop through all edges and add internal ones first
    forAll(oldEdges, edgeI)
    {
        // Ensure that we're adding valid edges
        if (oldEdgeFaces[edgeI].empty())
        {
            continue;
        }

        // Determine which patch this edge belongs to
        label patch = whichEdgePatch(edgeI);

        // Obtain references
        edge& thisEdge = oldEdges[edgeI];
        labelList& thisEF = oldEdgeFaces[edgeI];

        // Renumber edges
        if (thisEdge[0] < nOldPoints_)
        {
            thisEdge[0] = reversePointMap_[thisEdge[0]];
        }
        else
        {
            thisEdge[0] = addedPointRenumbering_[thisEdge[0]];
        }

        if (thisEdge[1] < nOldPoints_)
        {
            thisEdge[1] = reversePointMap_[thisEdge[1]];
        }
        else
        {
            thisEdge[1] = addedPointRenumbering_[thisEdge[1]];
        }

        // Renumber edgeFaces
        forAll(thisEF,faceI)
        {
            if (thisEF[faceI] < nOldFaces_)
            {
                thisEF[faceI] = reverseFaceMap_[thisEF[faceI]];
            }
            else
            {
                thisEF[faceI] = addedFaceRenumbering_[thisEF[faceI]];
            }
        }

        // Renumber edgePoints
        if (!twoDMesh_)
        {
            labelList& thisEP = oldEdgePoints[edgeI];

            forAll(thisEP,pointI)
            {
                if (thisEP[pointI] < nOldPoints_)
                {
                    thisEP[pointI] = reversePointMap_[thisEP[pointI]];
                }
                else
                {
                    thisEP[pointI] = addedPointRenumbering_[thisEP[pointI]];
                }
            }
        }

        // Update maps for boundary edges. Edge insertion for
        // boundaries will be done after internal edges.
        if (patch >= 0)
        {
            label bEdgeIndex = boundaryPatchIndices[patch]++;

            // Update the maps
            if (edgeI < nOldEdges_)
            {
                edgeMap_[bEdgeIndex] = edgeI;
                reverseEdgeMap_[edgeI] = bEdgeIndex;
            }
            else
            {
                addedEdgeRenumbering_.insert(edgeI, bEdgeIndex);
                addedEdgeReverseRenumbering.insert(bEdgeIndex, edgeI);

                edgeMap_[bEdgeIndex] = -1;
            }
        }
        else
        {
            // Renumber internal edges and add normally.
            if (edgeI < nOldEdges_)
            {
                edgeMap_[edgeInOrder] = edgeI;
                reverseEdgeMap_[edgeI] = edgeInOrder;
            }
            else
            {
                addedEdgeRenumbering_.insert(edgeI, edgeInOrder);
            }

            // Insert entities into local lists...
            edges_[edgeInOrder] = thisEdge;
            edgeFaces_[edgeInOrder] = thisEF;

            // Insert entities into mesh-reset lists
            edges[edgeInOrder] = thisEdge;
            edgeFaces[edgeInOrder].transfer(thisEF);

            if (!twoDMesh_)
            {
                edgePoints_[edgeInOrder] = oldEdgePoints[edgeI];
            }

            edgeInOrder++;
        }
    }

    // All internal edges have been inserted. Now insert boundary edges.
    label oldIndex;
    for(label i = nInternalEdges_; i < nEdges_; i++)
    {
        if (edgeMap_[i] == -1)
        {
            // This boundary edge was added during the topology change
            oldIndex = addedEdgeReverseRenumbering[i];
        }
        else
        {
            oldIndex = edgeMap_[i];
        }

        // Insert entities into local Lists...
        edges_[edgeInOrder] = oldEdges[oldIndex];
        edgeFaces_[edgeInOrder] = oldEdgeFaces[oldIndex];

        // Insert entities into mesh-reset lists
        edges[edgeInOrder] = oldEdges[oldIndex];
        edgeFaces[edgeInOrder].transfer(oldEdgeFaces[oldIndex]);

        if (!twoDMesh_)
        {
            edgePoints_[edgeInOrder] = oldEdgePoints[oldIndex];
        }

        edgeInOrder++;
    }

    // Final check to ensure everything went okay
    if (edgeInOrder != nEdges_)
    {
        FatalErrorIn("dynamicTopoFvMesh::reOrderEdges()") << nl
                << " Algorithm did not visit every edge in the mesh."
                << " Something's messed up." << nl
                << abort(FatalError);
    }

    // Renumber all faceEdges
    forAll(faceEdges_, faceI)
    {
        // Obtain references
        labelList& fEdges = faceEdges_[faceI];
        labelList& rfEdges = faceEdges[faceI];

        forAll(fEdges, edgeI)
        {
            if (fEdges[edgeI] < nOldEdges_)
            {
                fEdges[edgeI] = reverseEdgeMap_[fEdges[edgeI]];
                rfEdges[edgeI] = reverseEdgeMap_[rfEdges[edgeI]];
            }
            else
            {
                fEdges[edgeI] = addedEdgeRenumbering_[fEdges[edgeI]];
                rfEdges[edgeI] = addedEdgeRenumbering_[rfEdges[edgeI]];
            }
        }
    }

    // Invert edges to obtain pointEdges
    if (!twoDMesh_)
    {
        // Number of edges per point
        labelList nEdgesPerPoint(nPoints_, 0);

        forAll(edges_, edgeI)
        {
            nEdgesPerPoint[edges_[edgeI][0]]++;
            nEdgesPerPoint[edges_[edgeI][1]]++;
        }

        // Size pointEdges
        pointEdges_.setSize(nPoints_);

        forAll(nEdgesPerPoint, pointI)
        {
            pointEdges_[pointI].setSize(nEdgesPerPoint[pointI]);
        }

        nEdgesPerPoint = 0;

        // Fill pointEdges
        forAll(edges_, edgeI)
        {
            const edge& thisEdge = edges_[edgeI];

            pointEdges_[thisEdge[0]][nEdgesPerPoint[thisEdge[0]]++] = edgeI;
            pointEdges_[thisEdge[1]][nEdgesPerPoint[thisEdge[1]]++] = edgeI;
        }
    }

    // Clear the deleted entity map
    this->deletedEdges_.clear();
}

// Reorder faces in upper-triangular order after a topology change
void dynamicTopoFvMesh::reOrderFaces
(
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    labelListList& faceEdges
)
{
    // *** Face renumbering *** //
    // Faces have to be renumbered if any were added/deleted/modified
    // Boundary faces are added to respective patches.
    // Internal faces, however, have to be added in upper-triangular ordering;
    // i.e., in the increasing order of neighbours

    // Allocate for mapping information
    faceMap_.setSize(nFaces_, -1);

    label faceInOrder = 0, allFaces = faces_.size();
    faceList oldFaces(allFaces);
    labelList oldOwner(allFaces), oldNeighbour(allFaces), visited(allFaces,0);
    labelListList oldFaceEdges(allFaces);

    addedFaceRenumbering_.clear();
    Map<label> addedFaceReverseRenumbering;

    // Make a copy of the old face-based lists, and clear them
    forAll(faces_, faceI)
    {
        oldFaces[faceI].transfer(faces_[faceI]);
        oldOwner[faceI] = owner_[faceI];
        oldNeighbour[faceI] = neighbour_[faceI];
        oldFaceEdges[faceI].transfer(faceEdges_[faceI]);
    }

    faces_.setSize(nFaces_);
    owner_.setSize(nFaces_);
    neighbour_.setSize(nFaces_);
    faceEdges_.setSize(nFaces_);

    // Mark the internal faces with -2 so that they are inserted first
    forAll(cells_, cellI)
    {
        const cell& curFaces = cells_[cellI];

        forAll(curFaces, faceI)
        {
            visited[curFaces[faceI]]--;
        }
    }

    // Upper-triangular ordering of faces:

    // Keep track of inserted boundary face indices
    labelList boundaryPatchIndices(patchStarts_);

    // Insertion cannot be done in one go as the faces need to be
    // added into the list in the increasing order of neighbour
    // cells.  Therefore, all neighbours will be detected first
    // and then added in the correct order.
    forAll(cells_, cellI)
    {
        // Record the neighbour cell
        const cell& curFaces = cells_[cellI];

        labelList neiCells(curFaces.size(), -1);

        label nNeighbours = 0;

        forAll(curFaces, faceI)
        {
            if (visited[curFaces[faceI]] == -2)
            {
                // Face is internal and gets reordered
                label own =
                (
                    oldOwner[curFaces[faceI]] < nOldCells_
                  ? reverseCellMap_[oldOwner[curFaces[faceI]]]
                  : addedCellRenumbering_[oldOwner[curFaces[faceI]]]
                );

                label nei =
                (
                    oldNeighbour[curFaces[faceI]] < nOldCells_
                  ? reverseCellMap_[oldNeighbour[curFaces[faceI]]]
                  : addedCellRenumbering_[oldNeighbour[curFaces[faceI]]]
                );

                label smallerIndex = own < nei ? own : nei;
                label largerIndex  = own > nei ? own : nei;

                if (cellI == smallerIndex)
                {
                    neiCells[faceI] = largerIndex;
                    nNeighbours++;
                }
            }

            // Boundary faces are inserted normally. Update maps for now.
            // Face insertion for boundaries will be done after internal faces.
            if (visited[curFaces[faceI]] == -1)
            {
                label patchID = whichPatch(curFaces[faceI]);
                label bFaceIndex = boundaryPatchIndices[patchID]++;

                // Renumber the point-labels for this boundary-face
                face& faceRenumber = oldFaces[curFaces[faceI]];

                forAll(faceRenumber,pointI)
                {
                    if (faceRenumber[pointI] < nOldPoints_)
                    {
                        faceRenumber[pointI] =
                        (
                            reversePointMap_[faceRenumber[pointI]]
                        );
                    }
                    else
                    {
                        faceRenumber[pointI] =
                        (
                            addedPointRenumbering_[faceRenumber[pointI]]
                        );
                    }
                }

                // Update the maps
                if (curFaces[faceI] < nOldFaces_)
                {
                    faceMap_[bFaceIndex] = curFaces[faceI];
                    reverseFaceMap_[curFaces[faceI]] = bFaceIndex;
                }
                else
                {
                    addedFaceRenumbering_.insert
                    (
                        curFaces[faceI],
                        bFaceIndex
                    );

                    addedFaceReverseRenumbering.insert
                    (
                        bFaceIndex,
                        curFaces[faceI]
                    );

                    faceMap_[bFaceIndex] = -1;
                }

                // Mark this face as visited
                visited[curFaces[faceI]] = 0;
            }
        }

        // Add internal faces in the increasing order of neighbours
        for (label neiSearch = 0; neiSearch < nNeighbours; neiSearch++)
        {
            // Find the lowest neighbour which is still valid
            label nextNei = -1;
            label minNei = nCells_;

            forAll(neiCells, ncI)
            {
                if (neiCells[ncI] > -1 && neiCells[ncI] < minNei)
                {
                    nextNei = ncI;
                    minNei = neiCells[ncI];
                }
            }

            if (nextNei > -1)
            {
                // Face is internal and gets reordered
                if (curFaces[nextNei] < nOldFaces_)
                {
                    faceMap_[faceInOrder] = curFaces[nextNei];
                    reverseFaceMap_[curFaces[nextNei]] = faceInOrder;
                }
                else
                {
                    addedFaceRenumbering_.insert
                    (
                        curFaces[nextNei],
                        faceInOrder
                    );
                }

                // Renumber the point labels in this face
                face& faceRenumber = oldFaces[curFaces[nextNei]];

                forAll(faceRenumber, pointI)
                {
                    if (faceRenumber[pointI] < nOldPoints_)
                    {
                        faceRenumber[pointI] =
                        (
                            reversePointMap_[faceRenumber[pointI]]
                        );
                    }
                    else
                    {
                        faceRenumber[pointI] =
                        (
                            addedPointRenumbering_[faceRenumber[pointI]]
                        );
                    }
                }

                // Renumber owner/neighbour
                label oldOwn = oldOwner[curFaces[nextNei]];
                label oldNei = oldNeighbour[curFaces[nextNei]];

                label ownerRenumber =
                (
                    oldOwn < nOldCells_
                  ? reverseCellMap_[oldOwn] : addedCellRenumbering_[oldOwn]
                );

                label neighbourRenumber =
                (
                    oldNei < nOldCells_
                  ? reverseCellMap_[oldNei] : addedCellRenumbering_[oldNei]
                );

                // Cell-reordering may cause flipped faces.. Correct them.
                if (neighbourRenumber < ownerRenumber)
                {
                    faceRenumber = faceRenumber.reverseFace();
                }

                // Insert entities into local lists...
                faces_[faceInOrder] = faceRenumber;
                owner_[faceInOrder] = cellI;
                neighbour_[faceInOrder] = minNei;
                faceEdges_[faceInOrder] =
                (
                    oldFaceEdges[curFaces[nextNei]]
                );

                // Insert entities into mesh-reset lists
                faces[faceInOrder].transfer(faceRenumber);
                owner[faceInOrder] = cellI;
                neighbour[faceInOrder] = minNei;
                faceEdges[faceInOrder].transfer
                (
                    oldFaceEdges[curFaces[nextNei]]
                );

                // Stop the neighbour from being used again
                neiCells[nextNei] = -1;

                // Mark this face as visited
                visited[curFaces[nextNei]] = 0;

                faceInOrder++;
            }
            else
            {
                FatalErrorIn("dynamicTopoFvMesh::reOrderFaces()") << nl
                    << "Error in internal face insertion" << nl
                    << abort(FatalError);
            }
        }
    }

    // All internal faces have been inserted. Now insert boundary faces.
    label oldIndex;

    for(label i=nInternalFaces_; i<nFaces_; i++)
    {
        if (faceMap_[i] == -1)
        {
            // This boundary face was added during the topology change
            oldIndex = addedFaceReverseRenumbering[i];
        }
        else
        {
            oldIndex = faceMap_[i];
        }

        // Renumber owner/neighbour
        label ownerRenumber =   oldOwner[oldIndex] < nOldCells_
                              ? reverseCellMap_[oldOwner[oldIndex]]
                              : addedCellRenumbering_[oldOwner[oldIndex]];

        // Insert entities into local listsLists...
        faces_[faceInOrder] = oldFaces[oldIndex];
        owner_[faceInOrder] = ownerRenumber;
        neighbour_[faceInOrder] = -1;
        faceEdges_[faceInOrder] = oldFaceEdges[oldIndex];

        // Insert entities into mesh-reset lists
        // NOTE: From OF-1.5 onwards, neighbour array
        //       does not store -1 for boundary faces
        faces[faceInOrder].transfer(oldFaces[oldIndex]);
        owner[faceInOrder] = ownerRenumber;
        faceEdges[faceInOrder].transfer(oldFaceEdges[oldIndex]);

        faceInOrder++;
    }

    // Renumber all cells
    forAll(cells_, cellI)
    {
        cell& cellFaces = cells_[cellI];

        forAll(cellFaces, faceI)
        {
            if (cellFaces[faceI] < nOldFaces_)
            {
                cellFaces[faceI] = reverseFaceMap_[cellFaces[faceI]];
            }
            else
            {
                cellFaces[faceI] = addedFaceRenumbering_[cellFaces[faceI]];
            }
        }
    }

    // Final check to ensure everything went okay
    if (debug > 1)
    {
        if (sum(visited) != 0)
        {
            FatalErrorIn("dynamicTopoFvMesh::reOrderFaces()") << nl
                    << " Algorithm did not visit every face in the mesh."
                    << " Something's messed up." << nl
                    << abort(FatalError);
        }
    }

    // Clear the deleted entity map
    deletedFaces_.clear();
}

// Reorder & renumber cells with bandwidth reduction after a topology change
void dynamicTopoFvMesh::reOrderCells()
{
    // *** Cell renumbering *** //
    // If cells were deleted during topology change, the numerical order ceases
    // to be continuous. Also, cells are always added at the end of the list by
    // virtue of the append method. Thus, cells would now have to be
    // reordered so that bandwidth is reduced and renumbered to be sequential.

    // Allocate for mapping information
    cellMap_.setSize(nCells_, -1);

    label currentCell, cellInOrder = 0, allCells = cells_.size();
    SLList<label> nextCell;
    labelList ncc(allCells, 0);
    labelList visited(allCells, 0);
    labelListList cellCellAddr(allCells);
    cellList oldCells(allCells);

    addedCellRenumbering_.clear();

    // Make a copy of the old cell-based lists, and clear them
    forAll(cells_, cellI)
    {
        oldCells[cellI].transfer(cells_[cellI]);
    }

    cells_.setSize(nCells_);

    // Build a cell-cell addressing list
    forAll(owner_, faceI)
    {
        if ((neighbour_[faceI] > -1) && (owner_[faceI] > -1))
        {
            ncc[owner_[faceI]]++;
            ncc[neighbour_[faceI]]++;
        }
    }

    forAll(cellCellAddr, cellI)
    {
        cellCellAddr[cellI].setSize(ncc[cellI]);

        // Mark off deleted cells as "visited"
        if (ncc[cellI] == 0)
        {
            visited[cellI] = 1;
        }
    }

    ncc = 0;
    forAll(owner_, faceI)
    {
        if ((owner_[faceI] > -1) && (neighbour_[faceI] > -1))
        {
            cellCellAddr[owner_[faceI]][ncc[owner_[faceI]]++] =
            (
                neighbour_[faceI]
            );

            cellCellAddr[neighbour_[faceI]][ncc[neighbour_[faceI]]++] =
            (
                owner_[faceI]
            );
        }
    }

    // Let's get to the "business bit" of the band-compression
    forAll(visited, cellI)
    {
        // Find the first cell that has not been visited yet
        if (visited[cellI] == 0)
        {
            // Use this cell as a start
            currentCell = cellI;

            nextCell.append(currentCell);

            // Loop through the nextCell list. Add the first cell into the
            // cell order if it has not already been visited and ask for its
            // neighbours. If the neighbour in question has not been visited,
            // add it to the end of the nextCell list
            while (nextCell.size() > 0)
            {
                currentCell = nextCell.removeHead();

                if (visited[currentCell] == 0)
                {
                    // Mark as visited and update cell mapping info
                    visited[currentCell] = 1;

                    if (currentCell < nOldCells_)
                    {
                        cellMap_[cellInOrder] = currentCell;
                        reverseCellMap_[currentCell] = cellInOrder;
                    }
                    else
                    {
                        addedCellRenumbering_.insert(currentCell,cellInOrder);
                    }

                    // Insert entities into local lists...
                    cells_[cellInOrder].transfer(oldCells[currentCell]);

                    cellInOrder++;

                    // Find if the neighbours have been visited
                    const labelList& neighbours = cellCellAddr[currentCell];

                    forAll(neighbours, nI)
                    {
                        if (visited[neighbours[nI]] == 0)
                        {
                            // Not visited, add to the list
                            nextCell.append(neighbours[nI]);
                        }
                    }
                }
            }
        }
    }

    // Loop through the cellsFromCells list, and renumber the map indices
    // HashTable keys, however, are not altered.
    forAllIter(Map<objectMap>, cellsFromCells_, cellI)
    {
        objectMap& thisMap = cellI();

        if (thisMap.index() < nOldCells_)
        {
            thisMap.index() = reverseCellMap_[thisMap.index()];
        }
        else
        {
            thisMap.index() = addedCellRenumbering_[thisMap.index()];
        }
    }

    if (debug > 1)
    {
        if (sum(visited) != allCells)
        {
            FatalErrorIn("dynamicTopoFvMesh::reOrderCells()") << nl
                    << " Algorithm did not visit every cell in the mesh."
                    << " Something's messed up." << nl
                    << abort(FatalError);
        }
    }

    // Clear the deleted entity map
    deletedCells_.clear();
}

// Reorder the faces in upper-triangular order, and generate mapping information
void dynamicTopoFvMesh::reOrderMesh
(
    pointField& points,
    edgeList& edges,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    labelListList& faceEdges,
    labelListList& edgeFaces
)
{
    if (debug)
    {
        Info << endl;
        Info << "=================" << endl;
        Info << " Mesh reOrdering " << endl;
        Info << "=================" << endl;
        Info << "Mesh Info [n]:" << endl;
        Info << "Points: " << nOldPoints_ << endl;
        Info << "Edges: " << nOldEdges_ << endl;
        Info << "Faces: " << nOldFaces_ << endl;
        Info << "Cells: " << nOldCells_ << endl;
        Info << "Internal Edges: " << nOldInternalEdges_ << endl;
        Info << "Internal Faces: " << nOldInternalFaces_ << endl;
        Info << "Patch Starts [Face]: " << oldPatchStarts_ << endl;
        Info << "Patch Sizes  [Face]: " << oldPatchSizes_ << endl;
        Info << "Patch Starts [Edge]: " << oldEdgePatchStarts_ << endl;
        Info << "Patch Sizes  [Edge]: " << oldEdgePatchSizes_ << endl;
        Info << "=================" << endl;
        Info << "Mesh Info [n+1]:" << endl;
        Info << "Points: " << nPoints_ << endl;
        Info << "Edges: " << nEdges_ << endl;
        Info << "Faces: " << nFaces_ << endl;
        Info << "Cells: " << nCells_ << endl;
        Info << "Internal Edges: " << nInternalEdges_ << endl;
        Info << "Internal Faces: " << nInternalFaces_ << endl;
        Info << "Patch Starts [Face]: " << patchStarts_ << endl;
        Info << "Patch Sizes: [Face]: " << patchSizes_ << endl;
        Info << "Patch Starts [Edge]: " << edgePatchStarts_ << endl;
        Info << "Patch Sizes: [Edge]: " << edgePatchSizes_ << endl;
        Info << "=================" << endl;

        // Check connectivity structures for consistency
        // before entering the reOrdering phase.
        checkConnectivity();
    }

    // Reorder the points
    if (debug) Info << "ReOrdering points..." << endl;
    reOrderPoints(points);

    // Reorder the cells
    if (debug) Info << "ReOrdering cells..." << endl;
    reOrderCells();

    // Reorder the faces
    if (debug) Info << "ReOrdering faces..." << endl;
    reOrderFaces(faces, owner, neighbour, faceEdges);

    // Reorder the edges
    if (debug) Info << "ReOrdering edges..." << endl;
    reOrderEdges(edges, edgeFaces, faceEdges);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
