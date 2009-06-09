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
    pointMap_.setSize(nPoints_, -1);

    label pointRenum = 0;

    addedPointRenumbering_.clear();

    forAllIter(HashList<point>::iterator, meshPoints_, ptIter)
    {
        // Obtain the index for this point
        label pIndex = ptIter.index();

        // Update the point info
        points[pointRenum] = ptIter();

        // Added points are always numbered after nOldPoints_
        // (by virtue of the HashList append method)
        if (pIndex < nOldPoints_)
        {
            pointMap_[pointRenum]    = pIndex;
            reversePointMap_[pIndex] = pointRenum;
        }
        else
        {
            addedPointRenumbering_.insert(pIndex,pointRenum);
        }

        // Update the counter
        pointRenum++;
    }

    // Clear the HashList and reset.
    meshPoints_.clear();
    forAll(points, pointI)
    {
        meshPoints_.append(points[pointI]);
    }
}

// Reorder edges after a topology change
void dynamicTopoFvMesh::reOrderEdges
(
    edgeList& edges,
    labelListList& edgeFaces
)
{
    // *** Edge renumbering *** //
    // If edges were deleted during topology change, the numerical order ceases
    // to be continuous. Edges are added to respective internal/boundary patches

    // Allocate for mapping information
    edgeMap_.setSize(nEdges_, -1);

    label edgeInOrder = 0, allEdges = edges_.lastIndex() + 1;
    edgeList oldEdges(allEdges);
    labelListList oldEdgeFaces(allEdges);
    labelListList oldEdgePoints(allEdges);

    addedEdgeRenumbering_.clear();
    Map<label> addedEdgeReverseRenumbering;

    // Transfer old edge-based HashLists, and clear them
    HashList<edge>::iterator eIter = edges_.begin();
    HashList<labelList>::iterator efIter = edgeFaces_.begin();

    while (eIter != edges_.end())
    {
        oldEdges[eIter.index()] = eIter();
        oldEdgeFaces[efIter.index()].transfer(efIter());
        eIter++; efIter++;
    }

    edges_.clear(); edgeFaces_.clear();

    if (!twoDMesh_)
    {
        HashList<labelList>::iterator epIter = edgePoints_.begin();

        while (epIter != edgePoints_.end())
        {
            oldEdgePoints[epIter.index()].transfer(epIter());
            epIter++;
        }

        edgePoints_.clear();
    }

    // Keep track of inserted boundary edge indices
    labelList boundaryPatchIndices(edgePatchStarts_);

    // Loop through all edges and add internal ones first
    forAll(oldEdges, edgeI)
    {
        // Ensure that we're adding valid edges
        if (oldEdgeFaces[edgeI].size() == 0)
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
                addedEdgeRenumbering_.insert(edgeI,bEdgeIndex);
                addedEdgeReverseRenumbering.insert(bEdgeIndex,edgeI);
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
                addedEdgeRenumbering_.insert(edgeI,edgeInOrder);
            }

            // Insert entities into HashLists...
            edges_.append(thisEdge);
            edgeFaces_.append(thisEF);

            // Insert entities into mesh-reset lists
            edges[edgeInOrder] = thisEdge;
            edgeFaces[edgeInOrder].transfer(thisEF);

            if (!twoDMesh_)
            {
                edgePoints_.append(oldEdgePoints[edgeI]);
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

        // Insert entities into HashLists...
        edges_.append(oldEdges[oldIndex]);
        edgeFaces_.append(oldEdgeFaces[oldIndex]);

        // Insert entities into mesh-reset lists
        edges[edgeInOrder] = oldEdges[oldIndex];
        edgeFaces[edgeInOrder].transfer(oldEdgeFaces[oldIndex]);

        if (!twoDMesh_)
        {
            edgePoints_.append(oldEdgePoints[oldIndex]);
        }

        edgeInOrder++;
    }
}

// Reorder faces in upper-triangular order after a topology change
void dynamicTopoFvMesh::reOrderFaces
(
    faceList& faces,
    labelList& owner,
    labelList& neighbour
)
{
    // *** Face renumbering *** //
    // Faces have to be renumbered if any were added/deleted/modified
    // Boundary faces are added to respective patches.
    // Internal faces, however, have to be added in upper-triangular ordering;
    // i.e., in the increasing order of neighbours

    // Allocate for mapping information
    faceMap_.setSize(nFaces_, -1);

    label faceInOrder = 0, allFaces = faces_.lastIndex() + 1;
    faceList oldFaces(allFaces);
    labelList oldOwner(allFaces), oldNeighbour(allFaces), visited(allFaces,0);

    addedFaceRenumbering_.clear();
    Map<label> addedFaceReverseRenumbering;

    // Make a copy of the old face-based HashLists, and clear them
    HashList<face>::iterator fIter = faces_.begin();
    HashList<label>::iterator oIter = owner_.begin();
    HashList<label>::iterator nIter = neighbour_.begin();

    while(fIter != faces_.end())
    {
        oldFaces[fIter.index()].transfer(fIter());
        oldOwner[oIter.index()] = oIter();
        oldNeighbour[nIter.index()] = nIter();
        fIter++; oIter++; nIter++;
    }

    faces_.clear(); owner_.clear(); neighbour_.clear();

    // Mark the internal faces with -2 so that they are inserted first
    forAllIter(HashList<cell>::iterator, cells_, cIter)
    {
        const cell& curFaces = cIter();
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
    forAllIter(HashList<cell>::iterator, cells_, cIter)
    {
        // Record the neighbour cell
        label cellI = cIter.index();
        const cell& curFaces = cIter();
        labelList neiCells(curFaces.size(), -1);

        label nNeighbours = 0;

        forAll(curFaces, faceI)
        {
            if (visited[curFaces[faceI]] == -2)
            {
                // Face is internal and gets reordered
                label own =   oldOwner[curFaces[faceI]] < nOldCells_
                            ? reverseCellMap_[oldOwner[curFaces[faceI]]]
                            : addedCellRenumbering_[oldOwner[curFaces[faceI]]];
                label nei =   oldNeighbour[curFaces[faceI]] < nOldCells_
                            ? reverseCellMap_[oldNeighbour[curFaces[faceI]]]
                            : addedCellRenumbering_[oldNeighbour[curFaces[faceI]]];

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
                            reversePointMap_[faceRenumber[pointI]];
                    }
                    else
                    {
                        faceRenumber[pointI] =
                            addedPointRenumbering_[faceRenumber[pointI]];
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
                    addedFaceRenumbering_.insert(curFaces[faceI],bFaceIndex);
                    addedFaceReverseRenumbering.insert(bFaceIndex,curFaces[faceI]);
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
                    addedFaceRenumbering_.insert(curFaces[nextNei],faceInOrder);
                }

                // Renumber the point labels in this face
                face& faceRenumber = oldFaces[curFaces[nextNei]];
                forAll(faceRenumber, pointI)
                {
                    if (faceRenumber[pointI] < nOldPoints_)
                    {
                        faceRenumber[pointI] =
                            reversePointMap_[faceRenumber[pointI]];
                    }
                    else
                    {
                        faceRenumber[pointI] =
                            addedPointRenumbering_[faceRenumber[pointI]];
                    }
                }

                // Renumber owner/neighbour
                label oldOwn = oldOwner[curFaces[nextNei]];
                label oldNei = oldNeighbour[curFaces[nextNei]];

                label ownerRenumber =
                    oldOwner[curFaces[nextNei]] < nOldCells_
                  ? reverseCellMap_[oldOwn] : addedCellRenumbering_[oldOwn];

                label neighbourRenumber =
                    oldNeighbour[curFaces[nextNei]] < nOldCells_
                  ? reverseCellMap_[oldNei] : addedCellRenumbering_[oldNei];

                // Cell-reordering may cause flipped faces.. Correct them.
                if (neighbourRenumber < ownerRenumber)
                {
                    faceRenumber = faceRenumber.reverseFace();
                }

                // Insert entities into HashLists...
                faces_.append(faceRenumber);
                owner_.append(cellI);
                neighbour_.append(minNei);

                // Insert entities into mesh-reset lists
                faces[faceInOrder].transfer(faceRenumber);
                owner[faceInOrder] = cellI;
                neighbour[faceInOrder] = minNei;

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

        // Insert entities into HashLists...
        faces_.append(oldFaces[oldIndex]);
        owner_.append(ownerRenumber);
        neighbour_.append(-1);

        // Insert entities into mesh-reset lists
        // NOTE: From OF-1.5 onwards, neighbour array
        //       does not store -1 for boundary faces
        faces[faceInOrder].transfer(oldFaces[oldIndex]);
        owner[faceInOrder] = ownerRenumber;

        faceInOrder++;
    }

    // Renumber all cells
    forAllIter(HashList<cell>::iterator, cells_, cIter)
    {
        cell& cellFaces = cIter();
        forAll(cellFaces,faceI)
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
}

// Reorder & renumber cells with bandwidth reduction after a topology change
void dynamicTopoFvMesh::reOrderCells()
{
    // *** Cell renumbering *** //
    // If cells were deleted during topology change, the numerical order ceases
    // to be continuous. Also, cells are always added at the end of the list by
    // virtue of the HashList append method. Thus, cells would now have to be
    // reordered so that bandwidth is reduced and renumbered to be sequential.

    // Allocate for mapping information
    cellMap_.setSize(nCells_, -1);

    label currentCell, cellInOrder = 0, allCells = cells_.lastIndex() + 1;
    SLList<label> nextCell;
    labelList ncc(allCells, 0);
    labelList visited(allCells, 0);
    labelListList cellCellAddr(allCells);
    cellList oldCells(allCells);
    scalarField oldLengthScale(0);

    addedCellRenumbering_.clear();

    // Make a copy of the old cell-based HashLists, and clear them
    forAllIter(HashList<cell>::iterator, cells_, cIter)
    {
        oldCells[cIter.index()].transfer(cIter());
    }
    cells_.clear();

    if (edgeModification_)
    {
        oldLengthScale.setSize(allCells);
        forAllIter(HashList<scalar>::iterator, lengthScale_, cIter)
        {
            oldLengthScale[cIter.index()] = cIter();
        }
        lengthScale_.clear();
    }

    // Build a cell-cell addressing list
    HashList<label>::iterator ownIter = owner_.begin();
    HashList<label>::iterator neiIter = neighbour_.begin();

    while(ownIter != owner_.end())
    {
        if (neiIter() != -1)
        {
            ncc[ownIter()]++;
            ncc[neiIter()]++;
        }
        ownIter++; neiIter++;
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
    ownIter = owner_.begin(); neiIter = neighbour_.begin();
    while(ownIter != owner_.end())
    {
        if (neiIter() != -1)
        {
            cellCellAddr[ownIter()][ncc[ownIter()]++] = neiIter();
            cellCellAddr[neiIter()][ncc[neiIter()]++] = ownIter();
        }
        ownIter++; neiIter++;
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

                    // Insert entities into HashLists...
                    cells_.append(oldCells[currentCell]);

                    if (edgeModification_)
                    {
                        lengthScale_.append(oldLengthScale[currentCell]);
                    }

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
}

// Reorder the faces in upper-triangular order, and generate mapping information
void dynamicTopoFvMesh::reOrderMesh
(
    pointField& points,
    edgeList& edges,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
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
    }

    // Reorder the points
    if (debug) Info << "ReOrdering points..." << endl;
    reOrderPoints(points);

    // Reorder the cells
    if (debug) Info << "ReOrdering cells..." << endl;
    reOrderCells();

    // Reorder the faces
    if (debug) Info << "ReOrdering faces..." << endl;
    reOrderFaces(faces, owner, neighbour);

    // Reorder the edges
    if (debug) Info << "ReOrdering edges..." << endl;
    reOrderEdges(edges, edgeFaces);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
