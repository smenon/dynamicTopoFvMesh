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
#include "multiThreader.H"
#include "dynamicTopoFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Method for the swapping of a quad-face in 2D
void dynamicTopoFvMesh::swapQuadFace
(
    const label fIndex
)
{
    face f;
    bool found = false;
    label commonIndex = -1;
    FixedList<label,4> otherPointIndex(-1), nextToOtherPoint(-1);
    FixedList<label,2> c0BdyIndex, c0IntIndex, c1BdyIndex, c1IntIndex;
    FixedList<face,2>  c0BdyFace,  c0IntFace,  c1BdyFace,  c1IntFace;
    FixedList<label,2> commonEdgeIndex(-1);
    FixedList<edge,2>  commonEdges;
    FixedList<label,4> otherEdgeIndex(-1);
    FixedList<label,4> commonFaceIndex(-1), modifiedEdgeIndex(-1);
    FixedList<face,4>  commonFaces(face(3)), commonIntFaces(face(4));
    FixedList<label,4> commonIntFaceIndex(-1);
    FixedList<bool,2> foundTriFace0(false), foundTriFace1(false);
    FixedList<face,2> triFaces0(face(3)), triFaces1(face(3));

    // Get the two cells on either side...
    label c0 = owner_[fIndex];
    label c1 = neighbour_[fIndex];

    // Get cell references
    cell& cell_0 = cells_[c0];
    cell& cell_1 = cells_[c1];

    // Need to find common faces and edges...
    // At the end of this loop, commonFaces [0] & [1] share commonEdge[0]
    // and commonFaces [2] & [3] share commonEdge[1]
    // Also, commonFaces[0] & [2] are connected to cell[0],
    // and commonFaces[1] & [3] are connected to cell[1]
    labelList& fEdges = faceEdges_[fIndex];

    forAll(fEdges, edgeI)
    {
        // Break out if all triangular faces are found
        if
        (
            foundTriFace0[0] && foundTriFace0[1]
         && foundTriFace1[0] && foundTriFace1[1]
        )
        {
            break;
        }

        // Obtain edgeFaces for this edge
        labelList& eFaces = edgeFaces_[fEdges[edgeI]];

        forAll(eFaces, faceI)
        {
            face& eFace = faces_[eFaces[faceI]];

            if (eFace.size() == 3)
            {
                // Found a triangular face. Determine which cell it belongs to.
                if (owner_[eFaces[faceI]] == c0)
                {
                    if (foundTriFace0[0])
                    {
                        // Update the second face on cell[0].
                        commonIndex = 2;
                        foundTriFace0[1] = true;

                        if (foundTriFace1[1])
                        {
                            commonEdgeIndex[1] = fEdges[edgeI];
                            commonEdges[1] = edges_[fEdges[edgeI]];
                        }
                    }
                    else
                    {
                        // Update the first face on cell[0].
                        commonIndex = 0;
                        foundTriFace0[0] = true;

                        if (foundTriFace1[0])
                        {
                            commonEdgeIndex[0] = fEdges[edgeI];
                            commonEdges[0] = edges_[fEdges[edgeI]];
                        }
                    }
                }
                else
                {
                    if (foundTriFace1[0])
                    {
                        // Update the second face on cell[1].
                        commonIndex = 3;
                        foundTriFace1[1] = true;

                        if (foundTriFace0[1])
                        {
                            commonEdgeIndex[1] = fEdges[edgeI];
                            commonEdges[1] = edges_[fEdges[edgeI]];
                        }
                    }
                    else
                    {
                        // Update the first face on cell[1].
                        commonIndex = 1;
                        foundTriFace1[0] = true;

                        if (foundTriFace0[0])
                        {
                            commonEdgeIndex[0] = fEdges[edgeI];
                            commonEdges[0] = edges_[fEdges[edgeI]];
                        }
                    }
                }

                // Store the face and index
                commonFaces[commonIndex][0] = eFace[0];
                commonFaces[commonIndex][1] = eFace[1];
                commonFaces[commonIndex][2] = eFace[2];

                commonFaceIndex[commonIndex] = eFaces[faceI];
            }
        }
    }

    if (debug > 1)
    {
        Info << nl << nl << "Face: " << fIndex
             << " needs to be flipped. " << endl;

        Info << "Cell[0]: " << c0 << ": " << cell_0 << endl;
        Info << "Cell[1]: " << c1 << ": " << cell_1 << endl;

        if (debug > 2)
        {
            Info << "Common Faces: Set 1: "
                 << commonFaceIndex[0] << ": " << commonFaces[0] << ", "
                 << commonFaceIndex[1] << ": " << commonFaces[1] << endl;

            Info << "Common Faces: Set 2: "
                 << commonFaceIndex[2] << ": " << commonFaces[2] << ", "
                 << commonFaceIndex[3] << ": " << commonFaces[3] << endl;

            Info << "Old face: " << faces_[fIndex] << endl;
        }
    }

    // Find the interior/boundary faces.
    findPrismFaces
    (
        fIndex,
        c0,
        c0BdyFace,
        c0BdyIndex,
        c0IntFace,
        c0IntIndex
    );

    findPrismFaces
    (
        fIndex,
        c1,
        c1BdyFace,
        c1BdyIndex,
        c1IntFace,
        c1IntIndex
    );

    // Find the points that don't lie on shared edges
    // and the points next to them (for orientation)
    findIsolatedPoint
    (
        commonFaces[1],
        commonEdges[0],
        otherPointIndex[1],
        nextToOtherPoint[1]
    );

    findIsolatedPoint
    (
        commonFaces[0],
        commonEdges[0],
        otherPointIndex[0],
        nextToOtherPoint[0]
    );

    findIsolatedPoint
    (
        commonFaces[2],
        commonEdges[1],
        otherPointIndex[2],
        nextToOtherPoint[2]
    );

    findIsolatedPoint
    (
        commonFaces[3],
        commonEdges[1],
        otherPointIndex[3],
        nextToOtherPoint[3]
    );

    // Find the other two edges on the face being flipped
    forAll(fEdges, edgeI)
    {
        if
        (
            fEdges[edgeI] != commonEdgeIndex[0]
         && fEdges[edgeI] != commonEdgeIndex[1]
        )
        {
            // Obtain a reference to this edge
            edge& eThis = edges_[fEdges[edgeI]];

            if
            (
                eThis[0] == nextToOtherPoint[0]
             || eThis[1] == nextToOtherPoint[0]
            )
            {
                otherEdgeIndex[0] = fEdges[edgeI];
            }
            else
            {
                otherEdgeIndex[1] = fEdges[edgeI];
            }
        }
    }

    // At the end of this loop, commonIntFaces [0] & [1] share otherEdges[0]
    // and commonIntFaces [2] & [3] share the otherEdges[1],
    // where [0],[2] lie on cell[0] and [1],[3] lie on cell[1]
    found = false;

    labelList& e1 = faceEdges_[c0IntIndex[0]];

    forAll(e1,edgeI)
    {
        if (e1[edgeI] == otherEdgeIndex[0])
        {
            commonIntFaces[0] = c0IntFace[0];
            commonIntFaces[2] = c0IntFace[1];
            commonIntFaceIndex[0] = c0IntIndex[0];
            commonIntFaceIndex[2] = c0IntIndex[1];
            found = true; break;
        }
    }

    if (!found)
    {
        // The edge was not found before
        commonIntFaces[0] = c0IntFace[1];
        commonIntFaces[2] = c0IntFace[0];
        commonIntFaceIndex[0] = c0IntIndex[1];
        commonIntFaceIndex[2] = c0IntIndex[0];
    }

    found = false;

    labelList& e3 = faceEdges_[c1IntIndex[0]];

    forAll(e3,edgeI)
    {
        if (e3[edgeI] == otherEdgeIndex[0])
        {
            commonIntFaces[1] = c1IntFace[0];
            commonIntFaces[3] = c1IntFace[1];
            commonIntFaceIndex[1] = c1IntIndex[0];
            commonIntFaceIndex[3] = c1IntIndex[1];
            found = true; break;
        }
    }

    if (!found)
    {
        // The edge was not found before
        commonIntFaces[1] = c1IntFace[1];
        commonIntFaces[3] = c1IntFace[0];
        commonIntFaceIndex[1] = c1IntIndex[1];
        commonIntFaceIndex[3] = c1IntIndex[0];
    }

    // Find two common edges between interior/interior faces
    findCommonEdge
    (
        c0IntIndex[0],
        c0IntIndex[1],
        otherEdgeIndex[2]
    );

    findCommonEdge
    (
        c1IntIndex[0],
        c1IntIndex[1],
        otherEdgeIndex[3]
    );

    // Find four common edges between interior/boundary faces
    findCommonEdge
    (
        commonFaceIndex[1],
        commonIntFaceIndex[1],
        modifiedEdgeIndex[0]
    );

    findCommonEdge
    (
        commonFaceIndex[3],
        commonIntFaceIndex[1],
        modifiedEdgeIndex[1]
    );

    findCommonEdge
    (
        commonFaceIndex[0],
        commonIntFaceIndex[2],
        modifiedEdgeIndex[2]
    );

    findCommonEdge
    (
        commonFaceIndex[2],
        commonIntFaceIndex[2],
        modifiedEdgeIndex[3]
    );

    // Modify the five faces belonging to this hull
    face& newFace = faces_[fIndex];
    face& newBdyFace0 = faces_[commonFaceIndex[0]];
    face& newBdyFace1 = faces_[commonFaceIndex[1]];
    face& newBdyFace2 = faces_[commonFaceIndex[2]];
    face& newBdyFace3 = faces_[commonFaceIndex[3]];
    labelList& newFEdges = faceEdges_[fIndex];
    label c0count=0, c1count=0;

    // Size down edgeFaces for the original face
    sizeDownList
    (
        fIndex,
        edgeFaces_[otherEdgeIndex[0]]
    );

    sizeDownList
    (
        fIndex,
        edgeFaces_[otherEdgeIndex[1]]
    );

    // Size up edgeFaces for the face after flipping
    sizeUpList
    (
        fIndex,
        edgeFaces_[otherEdgeIndex[2]]
    );

    sizeUpList
    (
        fIndex,
        edgeFaces_[otherEdgeIndex[3]]
    );

    // Replace edgeFaces and faceEdges
    replaceLabel
    (
        modifiedEdgeIndex[0],
        modifiedEdgeIndex[1],
        faceEdges_[commonFaceIndex[1]]
    );

    replaceLabel
    (
        commonFaceIndex[1],
        commonFaceIndex[0],
        edgeFaces_[modifiedEdgeIndex[0]]
    );

    replaceLabel
    (
        modifiedEdgeIndex[1],
        modifiedEdgeIndex[0],
        faceEdges_[commonFaceIndex[3]]
    );

    replaceLabel
    (
        commonFaceIndex[3],
        commonFaceIndex[2],
        edgeFaces_[modifiedEdgeIndex[1]]
    );

    replaceLabel
    (
        modifiedEdgeIndex[2],
        modifiedEdgeIndex[3],
        faceEdges_[commonFaceIndex[0]]
    );

    replaceLabel
    (
        commonFaceIndex[0],
        commonFaceIndex[1],
        edgeFaces_[modifiedEdgeIndex[2]]
    );

    replaceLabel
    (
        modifiedEdgeIndex[3],
        modifiedEdgeIndex[2],
        faceEdges_[commonFaceIndex[2]]
    );

    replaceLabel
    (
        commonFaceIndex[2],
        commonFaceIndex[3],
        edgeFaces_[modifiedEdgeIndex[3]]
    );

    // Define parameters for the new flipped face
    newFace[0] = otherPointIndex[0];
    newFace[1] = otherPointIndex[1];
    newFace[2] = otherPointIndex[3];
    newFace[3] = otherPointIndex[2];
    newFEdges[0] = otherEdgeIndex[2];
    newFEdges[1] = commonEdgeIndex[0];
    newFEdges[2] = otherEdgeIndex[3];
    newFEdges[3] = commonEdgeIndex[1];
    cell_0[c0count++] = fIndex;
    cell_1[c1count++] = fIndex;
    owner_[fIndex] = c0;
    neighbour_[fIndex] = c1;

    // Four modified boundary faces need to be constructed,
    // but right-handedness is also important.
    // Take a cue from the existing boundary-face orientation

    // Zeroth boundary face - Owner c[0], Neighbour -1
    newBdyFace0[0] = otherPointIndex[0];
    newBdyFace0[1] = nextToOtherPoint[0];
    newBdyFace0[2] = otherPointIndex[1];
    cell_0[c0count++] = commonFaceIndex[0];
    owner_[commonFaceIndex[0]] = c0;
    neighbour_[commonFaceIndex[0]] = -1;

    // First boundary face - Owner c[1], Neighbour -1
    newBdyFace1[0] = otherPointIndex[1];
    newBdyFace1[1] = nextToOtherPoint[1];
    newBdyFace1[2] = otherPointIndex[0];
    cell_1[c1count++] = commonFaceIndex[1];
    owner_[commonFaceIndex[1]] = c1;
    neighbour_[commonFaceIndex[1]] = -1;

    // Second boundary face - Owner c[0], Neighbour -1
    newBdyFace2[0] = otherPointIndex[3];
    newBdyFace2[1] = nextToOtherPoint[3];
    newBdyFace2[2] = otherPointIndex[2];
    cell_0[c0count++] = commonFaceIndex[2];
    owner_[commonFaceIndex[2]] = c0;
    neighbour_[commonFaceIndex[2]] = -1;

    // Third boundary face - Owner c[1], Neighbour -1
    newBdyFace3[0] = otherPointIndex[2];
    newBdyFace3[1] = nextToOtherPoint[2];
    newBdyFace3[2] = otherPointIndex[3];
    cell_1[c1count++] = commonFaceIndex[3];
    owner_[commonFaceIndex[3]] = c1;
    neighbour_[commonFaceIndex[3]] = -1;

    if (debug > 1)
    {
        Info << "New flipped face: " << newFace << endl;

        if (debug > 2)
        {
            Info << "New boundary face[0]" << commonFaceIndex[0]
                 << ": " << newBdyFace0 << endl;
            Info << "New boundary face[1]" << commonFaceIndex[1]
                 << ": " << newBdyFace1 << endl;
            Info << "New boundary face[2]" << commonFaceIndex[2]
                 << ": " << newBdyFace2 << endl;
            Info << "New boundary face[3]" << commonFaceIndex[3]
                 << ": " << newBdyFace3 << endl;
        }
    }

    // Check the orientation of the two quad faces, and modify as necessary
    label newOwn=0, newNei=0;

    // The quad face belonging to cell[1] now becomes a part of cell[0]
    if (neighbour_[commonIntFaceIndex[1]] == -1)
    {
        // Boundary face
        // Face doesn't need to be flipped, just update the owner
        f          = commonIntFaces[1];
        newOwn     = c0;
        newNei     = -1;
    }
    else
    if (owner_[commonIntFaceIndex[1]] == c1)
    {
        // This face is on the interior, check for previous owner
        // Upper-triangular ordering has to be maintained, however...
        if (c0 > neighbour_[commonIntFaceIndex[1]])
        {
            // Flip is necessary
            f          = commonIntFaces[1].reverseFace();
            newOwn     = neighbour_[commonIntFaceIndex[1]];
            newNei     = c0;
        }
        else
        {
            // Flip isn't necessary, just change the owner
            f          = commonIntFaces[1];
            newOwn     = c0;
            newNei     = neighbour_[commonIntFaceIndex[1]];
        }
    }
    else
    if (neighbour_[commonIntFaceIndex[1]] == c1)
    {
        // This face is on the interior, check for previous neighbour
        // Upper-triangular ordering has to be maintained, however...
        if (c0 < owner_[commonIntFaceIndex[1]])
        {
            // Flip is necessary
            f          = commonIntFaces[1].reverseFace();
            newOwn     = c0;
            newNei     = owner_[commonIntFaceIndex[1]];
        }
        else
        {
            // Flip isn't necessary, just change the neighbour
            f          = commonIntFaces[1];
            newOwn     = owner_[commonIntFaceIndex[1]];
            newNei     = c0;
        }
    }

    faces_[commonIntFaceIndex[1]] = f;
    cell_0[c0count++] = commonIntFaceIndex[0];
    cell_0[c0count++] = commonIntFaceIndex[1];
    owner_[commonIntFaceIndex[1]] = newOwn;
    neighbour_[commonIntFaceIndex[1]] = newNei;

    // The quad face belonging to cell[0] now becomes a part of cell[1]
    if (neighbour_[commonIntFaceIndex[2]] == -1)
    {
        // Boundary face
        // Face doesn't need to be flipped, just update the owner
        f          = commonIntFaces[2];
        newOwn     = c1;
        newNei     = -1;
    }
    else
    if (owner_[commonIntFaceIndex[2]] == c0)
    {
        // This face is on the interior, check for previous owner
        // Upper-triangular ordering has to be maintained, however...
        if (c1 > neighbour_[commonIntFaceIndex[2]])
        {
            // Flip is necessary
            f          = commonIntFaces[2].reverseFace();
            newOwn     = neighbour_[commonIntFaceIndex[2]];
            newNei     = c1;
        }
        else
        {
            // Flip isn't necessary, just change the owner
            f          = commonIntFaces[2];
            newOwn     = c1;
            newNei     = neighbour_[commonIntFaceIndex[2]];
        }
    }
    else
    if (neighbour_[commonIntFaceIndex[2]] == c0)
    {
        // This face is on the interior, check for previous neighbour
        // Upper-triangular ordering has to be maintained, however...
        if (c1 < owner_[commonIntFaceIndex[2]])
        {
            // Flip is necessary
            f          = commonIntFaces[2].reverseFace();
            newOwn     = c1;
            newNei     = owner_[commonIntFaceIndex[2]];
        }
        else
        {
            // Flip isn't necessary, just change the neighbour
            f          = commonIntFaces[2];
            newOwn     = owner_[commonIntFaceIndex[2]];
            newNei     = c1;
        }
    }

    faces_[commonIntFaceIndex[2]] = f;
    cell_1[c1count++] = commonIntFaceIndex[2];
    cell_1[c1count++] = commonIntFaceIndex[3];
    owner_[commonIntFaceIndex[2]] = newOwn;
    neighbour_[commonIntFaceIndex[2]] = newNei;

    // Generate mapping information for both cells
    label firstParent, secondParent;
    const labelListList& cc = cellCells();
    labelHashSet c0MasterObjects(6);
    labelHashSet c1MasterObjects(6);

    if (c0 < nOldCells_)
    {
        firstParent = c0;
    }
    else
    {
        firstParent = cellParents_[c0];
    }

    if (c1 < nOldCells_)
    {
        secondParent = c1;
    }
    else
    {
        secondParent = cellParents_[c1];
    }

    // Find the cell's neighbours in the old mesh
    c0MasterObjects.insert(firstParent);
    c1MasterObjects.insert(firstParent);

    forAll(cc[firstParent],cellI)
    {
        if (!c0MasterObjects.found(cc[firstParent][cellI]))
        {
            c0MasterObjects.insert(cc[firstParent][cellI]);
        }
        if (!c1MasterObjects.found(cc[firstParent][cellI]))
        {
            c1MasterObjects.insert(cc[firstParent][cellI]);
        }
    }

    c0MasterObjects.insert(secondParent);
    c1MasterObjects.insert(secondParent);

    forAll(cc[secondParent],cellI)
    {
        if (!c0MasterObjects.found(cc[secondParent][cellI]))
        {
            c0MasterObjects.insert(cc[secondParent][cellI]);
        }
        if (!c1MasterObjects.found(cc[secondParent][cellI]))
        {
            c1MasterObjects.insert(cc[secondParent][cellI]);
        }
    }

    // Insert mapping info into the HashTable
    cellsFromCells_.insert(c0,objectMap(c0,c0MasterObjects.toc()));
    cellsFromCells_.insert(c1,objectMap(c1,c1MasterObjects.toc()));

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    nSwaps_++;
}

// Routine to perform 2-3 swaps
void dynamicTopoFvMesh::swap23
(
    const label isolatedVertex,
    const label eIndex,
    const label triangulationIndex,
    const label numTriangulations,
    const labelListList& triangulations,
    const labelList& hullVertices,
    const labelList& hullFaces,
    const labelList& hullCells
)
{
    // A 2-3 swap performs the following operations:
    //      [1] Remove face: [ edge[0] edge[1] isolatedVertex ]
    //      [2] Remove two cells on either side of removed face
    //      [3] Add one edge
    //      [4] Add three new faces
    //      [5] Add three new cells
    //      Update faceEdges, edgeFaces and edgePoints information

    // Figure out which edge this is...
    label tIndex = self();

    // Obtain edge reference
    edge& edgeToCheck = edges_[eIndex];

    label faceForRemoval = hullFaces[isolatedVertex];
    label vertexForRemoval = hullVertices[isolatedVertex];

    // Determine the two cells to be removed
    FixedList<label,2> cellsForRemoval;
    cellsForRemoval[0] = owner_[faceForRemoval];
    cellsForRemoval[1] = neighbour_[faceForRemoval];

    if (debug > 2)
    {
        // Print out arguments
        Info << endl;
        Info << "== Swapping 2-3 ==" << endl;
        Info << "Edge: " << eIndex << ": " << edgeToCheck << endl;
        Info << "Ring: " << hullVertices << endl;
        Info << "Faces: " << hullFaces << endl;
        Info << "Cells: " << hullCells << endl;
        Info << "Triangulation: "
             << triangulations[0][triangulationIndex] << " "
             << triangulations[1][triangulationIndex] << " "
             << triangulations[2][triangulationIndex] << " "
             << endl;
        Info << "Isolated vertex: " << isolatedVertex << endl;

        if (debug > 3)
        {
            writeVTK
            (
                Foam::name(eIndex)
              + "_beforeSwap_"
              + Foam::name(numTriangulations) + "_"
              + Foam::name(triangulationIndex),
                cellsForRemoval
            );
        }
    }

    // Add three new cells to the end of the cell list
    FixedList<label,3> newCellIndex(-1);
    newCellIndex[0] = cells_.append(cell(4));
    newCellIndex[1] = cells_.append(cell(4));
    newCellIndex[2] = cells_.append(cell(4));

    cell &newTetCell0 = cells_[newCellIndex[0]];
    cell &newTetCell1 = cells_[newCellIndex[1]];
    cell &newTetCell2 = cells_[newCellIndex[2]];

    // Update length-scale info
    if (edgeModification_)
    {
        scalar avgScale =
        (
             lengthScale_[cellsForRemoval[0]]
           + lengthScale_[cellsForRemoval[1]]
        )/2.0;

        for (label i = 0; i < 3; i++)
        {
            lengthScale_.append(avgScale);
        }
    }

    // Obtain point-ordering for the other vertices
    // otherVertices[0] is the point before isolatedVertex
    // otherVertices[1] is the point after isolatedVertex
    FixedList<label,2> otherVertices;

    if (triangulations[0][triangulationIndex] == isolatedVertex)
    {
        otherVertices[0] = hullVertices[triangulations[2][triangulationIndex]];
        otherVertices[1] = hullVertices[triangulations[1][triangulationIndex]];
    }
    else
    if (triangulations[1][triangulationIndex] == isolatedVertex)
    {
        otherVertices[0] = hullVertices[triangulations[0][triangulationIndex]];
        otherVertices[1] = hullVertices[triangulations[2][triangulationIndex]];
    }
    else
    {
        otherVertices[0] = hullVertices[triangulations[1][triangulationIndex]];
        otherVertices[1] = hullVertices[triangulations[0][triangulationIndex]];
    }

    // Insert three new internal faces
    FixedList<label,3> newFaceIndex;
    face tmpTriFace(3);

    // First face: The actual triangulation
    tmpTriFace[0] = otherVertices[0];
    tmpTriFace[1] = vertexForRemoval;
    tmpTriFace[2] = otherVertices[1];

    newFaceIndex[0] = insertFace
                      (
                          -1,
                          tmpTriFace,
                          newCellIndex[0],
                          newCellIndex[1]
                      );

    // Second face: Triangle involving edgeToCheck[0]
    tmpTriFace[0] = otherVertices[0];
    tmpTriFace[1] = edgeToCheck[0];
    tmpTriFace[2] = otherVertices[1];

    newFaceIndex[1] = insertFace
                      (
                          -1,
                          tmpTriFace,
                          newCellIndex[1],
                          newCellIndex[2]
                      );

    // Third face: Triangle involving edgeToCheck[1]
    tmpTriFace[0] = otherVertices[1];
    tmpTriFace[1] = edgeToCheck[1];
    tmpTriFace[2] = otherVertices[0];

    newFaceIndex[2] = insertFace
                      (
                          -1,
                          tmpTriFace,
                          newCellIndex[0],
                          newCellIndex[2]
                      );

    // Append three dummy faceEdges entries.
    for (label i = 0; i < 3; i++)
    {
        faceEdges_.append();
    }

    // Add an entry to edgeFaces
    labelList newEdgeFaces(3);
    newEdgeFaces[0] = newFaceIndex[0];
    newEdgeFaces[1] = newFaceIndex[1];
    newEdgeFaces[2] = newFaceIndex[2];

    // Add an entry for edgePoints as well
    labelList newEdgePoints(3);
    newEdgePoints[0] = vertexForRemoval;
    newEdgePoints[1] = edgeToCheck[0];
    newEdgePoints[2] = edgeToCheck[1];

    // Add a new internal edge to the mesh
    label newEdgeIndex = insertEdge
                         (
                             -1,
                             edge
                             (
                                 otherVertices[0],
                                 otherVertices[1]
                             ),
                             newEdgeFaces,
                             newEdgePoints
                         );

    // Define the six edges to check while building faceEdges:
    FixedList<edge,6> check;

    check[0][0] = vertexForRemoval; check[0][1] = otherVertices[0];
    check[1][0] = vertexForRemoval; check[1][1] = otherVertices[1];
    check[2][0] = edgeToCheck[0];   check[2][1] = otherVertices[0];
    check[3][0] = edgeToCheck[1];   check[3][1] = otherVertices[0];
    check[4][0] = edgeToCheck[0];   check[4][1] = otherVertices[1];
    check[5][0] = edgeToCheck[1];   check[5][1] = otherVertices[1];

    // Add three new entries to faceEdges
    label nE0 = 0, nE1 = 0, nE2 = 0;
    labelListList newFaceEdges(3,labelList(3));

    newFaceEdges[0][nE0++] = newEdgeIndex;
    newFaceEdges[1][nE1++] = newEdgeIndex;
    newFaceEdges[2][nE2++] = newEdgeIndex;

    // Fill-in information for the three new cells,
    // and correct info on existing neighbouring cells
    label nF0 = 0, nF1 = 0, nF2 = 0;
    FixedList<bool,2> foundEdge;

    // Add the newly created faces to cells
    newTetCell0[nF0++] = newFaceIndex[0];
    newTetCell0[nF0++] = newFaceIndex[2];
    newTetCell1[nF1++] = newFaceIndex[0];
    newTetCell1[nF1++] = newFaceIndex[1];
    newTetCell2[nF2++] = newFaceIndex[1];
    newTetCell2[nF2++] = newFaceIndex[2];

    forAll(cellsForRemoval, cellI)
    {
        label cellIndex = cellsForRemoval[cellI];
        cell& cellToCheck = cells_[cellIndex];

        forAll(cellToCheck, faceI)
        {
            label faceIndex = cellToCheck[faceI];
            face& faceToCheck = faces_[faceIndex];

            foundEdge[0] = false; foundEdge[1] = false;

            // Check if face contains edgeToCheck[0]
            if
            (
                (faceToCheck[0] == edgeToCheck[0])
             || (faceToCheck[1] == edgeToCheck[0])
             || (faceToCheck[2] == edgeToCheck[0])
            )
            {
                foundEdge[0] = true;
            }

            // Check if face contains edgeToCheck[1]
            if
            (
                (faceToCheck[0] == edgeToCheck[1])
             || (faceToCheck[1] == edgeToCheck[1])
             || (faceToCheck[2] == edgeToCheck[1])
            )
            {
                foundEdge[1] = true;
            }

            // Face is connected to edgeToCheck[0]
            if (foundEdge[0] && !foundEdge[1])
            {
                // Check if a face-flip is necessary
                if (owner_[faceIndex] == cellIndex)
                {
                    if (neighbour_[faceIndex] == -1)
                    {
                        // Change the owner
                        owner_[faceIndex] = newCellIndex[1];
                    }
                    else
                    {
                        // Flip this face
                        faces_[faceIndex] = faceToCheck.reverseFace();
                        owner_[faceIndex] = neighbour_[faceIndex];
                        neighbour_[faceIndex] = newCellIndex[1];
                    }
                }
                else
                {
                    // Flip is unnecessary. Just update neighbour
                    neighbour_[faceIndex] = newCellIndex[1];
                }

                // Add this face to the cell
                newTetCell1[nF1++] = faceIndex;

                // Update faceEdges, edgeFaces, and edgePoints.
                // Add them to the stack as well.
                const labelList& fEdges = faceEdges_[faceIndex];

                forAll(fEdges, edgeI)
                {
                    if (edges_[fEdges[edgeI]] == check[0])
                    {
                        newFaceEdges[0][nE0++] = fEdges[edgeI];

                        sizeUpList
                        (
                            newFaceIndex[0],
                            edgeFaces_[fEdges[edgeI]]
                        );

                        insertLabel
                        (
                            otherVertices[1],
                            edgeToCheck[0],
                            edgeToCheck[1],
                            edgePoints_[fEdges[edgeI]]
                        );

                        edgeStack(tIndex).push(fEdges[edgeI]);
                    }

                    if (edges_[fEdges[edgeI]] == check[1])
                    {
                        newFaceEdges[0][nE0++] = fEdges[edgeI];

                        sizeUpList
                        (
                            newFaceIndex[0],
                            edgeFaces_[fEdges[edgeI]]
                        );

                        insertLabel
                        (
                            otherVertices[0],
                            edgeToCheck[0],
                            edgeToCheck[1],
                            edgePoints_[fEdges[edgeI]]
                        );

                        edgeStack(tIndex).push(fEdges[edgeI]);
                    }

                    if (edges_[fEdges[edgeI]] == check[2])
                    {
                        newFaceEdges[1][nE1++] = fEdges[edgeI];

                        sizeUpList
                        (
                            newFaceIndex[1],
                            edgeFaces_[fEdges[edgeI]]
                        );

                        insertLabel
                        (
                            otherVertices[1],
                            vertexForRemoval,
                            edgeToCheck[1],
                            edgePoints_[fEdges[edgeI]]
                        );

                        edgeStack(tIndex).push(fEdges[edgeI]);
                    }

                    if (edges_[fEdges[edgeI]] == check[4])
                    {
                        newFaceEdges[1][nE1++] = fEdges[edgeI];

                        sizeUpList
                        (
                            newFaceIndex[1],
                            edgeFaces_[fEdges[edgeI]]
                        );

                        insertLabel
                        (
                            otherVertices[0],
                            vertexForRemoval,
                            edgeToCheck[1],
                            edgePoints_[fEdges[edgeI]]
                        );

                        edgeStack(tIndex).push(fEdges[edgeI]);
                    }
                }
            }

            // Face is connected to edgeToCheck[1]
            if (!foundEdge[0] && foundEdge[1])
            {
                // Check if a face-flip is necessary
                if (owner_[faceIndex] == cellIndex)
                {
                    if (neighbour_[faceIndex] == -1)
                    {
                        // Change the owner
                        owner_[faceIndex] = newCellIndex[0];
                    }
                    else
                    {
                        // Flip this face
                        faces_[faceIndex] = faceToCheck.reverseFace();
                        owner_[faceIndex] = neighbour_[faceIndex];
                        neighbour_[faceIndex] = newCellIndex[0];
                    }
                }
                else
                {
                    // Flip is unnecessary. Just update neighbour
                    neighbour_[faceIndex] = newCellIndex[0];
                }

                // Add this face to the cell
                newTetCell0[nF0++] = faceIndex;

                // Update faceEdges, edgeFaces, and edgePoints.
                // Add them to the stack as well.
                const labelList& fEdges = faceEdges_[faceIndex];

                forAll(fEdges, edgeI)
                {
                    if (edges_[fEdges[edgeI]] == check[3])
                    {
                        newFaceEdges[2][nE2++] = fEdges[edgeI];

                        sizeUpList
                        (
                            newFaceIndex[2],
                            edgeFaces_[fEdges[edgeI]]
                        );

                        insertLabel
                        (
                            otherVertices[1],
                            vertexForRemoval,
                            edgeToCheck[0],
                            edgePoints_[fEdges[edgeI]]
                        );

                        edgeStack(tIndex).push(fEdges[edgeI]);
                    }

                    if (edges_[fEdges[edgeI]] == check[5])
                    {
                        newFaceEdges[2][nE2++] = fEdges[edgeI];

                        sizeUpList
                        (
                            newFaceIndex[2],
                            edgeFaces_[fEdges[edgeI]]
                        );

                        insertLabel
                        (
                            otherVertices[0],
                            vertexForRemoval,
                            edgeToCheck[0],
                            edgePoints_[fEdges[edgeI]]
                        );

                        edgeStack(tIndex).push(fEdges[edgeI]);
                    }
                }
            }

            // Face is connected to both edgeToCheck [0] and [1]
            if
            (
                (foundEdge[0] && foundEdge[1])
             && (faceIndex != faceForRemoval)
            )
            {
                // Check if a face-flip is necessary
                if (owner_[faceIndex] == cellIndex)
                {
                    if (neighbour_[faceIndex] == -1)
                    {
                        // Change the owner
                        owner_[faceIndex] = newCellIndex[2];
                    }
                    else
                    {
                        // Flip this face
                        faces_[faceIndex] = faceToCheck.reverseFace();
                        owner_[faceIndex] = neighbour_[faceIndex];
                        neighbour_[faceIndex] = newCellIndex[2];
                    }
                }
                else
                {
                    // Flip is unnecessary. Just update neighbour
                    neighbour_[faceIndex] = newCellIndex[2];
                }

                // Add this face to the cell
                newTetCell2[nF2++] = faceIndex;
            }
        }
    }

    // Now update faceEdges for the three new faces
    forAll(newFaceEdges, faceI)
    {
        faceEdges_[newFaceIndex[faceI]] = newFaceEdges[faceI];
    }

    // Generate mapping information for the three new cells
    // Prepare a list of master-objects to map from
    const labelListList& cc = cellCells();
    labelHashSet masterObjects;
    FixedList<label,2> parents(-1);

    forAll(cellsForRemoval, indexI)
    {
        // Determine an appropriate parent cell
        if (cellsForRemoval[indexI] < nOldCells_)
        {
            parents[indexI] = cellsForRemoval[indexI];
        }
        else
        {
            parents[indexI] = cellParents_[cellsForRemoval[indexI]];
        }

        // Find the cell's neighbours in the old mesh
        masterObjects.insert(parents[indexI]);
        forAll(cc[parents[indexI]],cellI)
        {
            if (!masterObjects.found(cc[parents[indexI]][cellI]))
            {
                masterObjects.insert(cc[parents[indexI]][cellI]);
            }
        }
    }

    forAll(newCellIndex, cellI)
    {
        // Insert the parent cell [from first by default]
        cellParents_.insert(newCellIndex[cellI], parents[0]);

        // Insert mapping info into the HashTable
        cellsFromCells_.insert
        (
            newCellIndex[cellI],
            objectMap
            (
                newCellIndex[cellI],
                masterObjects.toc()
            )
        );
    }

    // Update edgeFaces and edgePoints for edges of the removed face
    label otherPoint = -1, nextPoint = -1;
    face& checkFace = faces_[faceForRemoval];
    labelList& fEdges = faceEdges_[faceForRemoval];

    forAll(fEdges, edgeI)
    {
        sizeDownList
        (
            faceForRemoval,
            edgeFaces_[fEdges[edgeI]]
        );

        // Find the isolated point and remove it
        findIsolatedPoint
        (
            checkFace,
            edges_[fEdges[edgeI]],
            otherPoint,
            nextPoint
        );

        sizeDownList
        (
            otherPoint,
            edgePoints_[fEdges[edgeI]]
        );

        edgeStack(tIndex).push(fEdges[edgeI]);
    }

    // Remove the face
    removeFace(faceForRemoval);

    // Update the number of cells, and the reverse cell map
    nCells_++;

    forAll(cellsForRemoval, cellI)
    {
        label cIndex = cellsForRemoval[cellI];

        if (debug > 2)
        {
            Info << "Removing cell: "
                 << cIndex << ": "
                 << cells_[cIndex]
                 << endl;
        }

        cells_.remove(cIndex);

        if (edgeModification_)
        {
            lengthScale_.remove(cIndex);
        }

        if (cIndex < nOldCells_)
        {
            reverseCellMap_[cIndex] = -1;
        }

        // Check if the cell was added in the current morph, and delete
        if (cellsFromCells_.found(cIndex))
        {
            cellsFromCells_.erase(cIndex);
        }
    }

    if (debug > 2)
    {
        Info << "Added edge: " << endl;

        Info << newEdgeIndex << ":: "
             << edges_[newEdgeIndex]
             << " edgeFaces: "
             << edgeFaces_[newEdgeIndex]
             << endl;

        Info << "Added faces: " << endl;

        forAll(newFaceIndex, faceI)
        {
            Info << newFaceIndex[faceI] << ":: "
                 << faces_[newFaceIndex[faceI]]
                 << " faceEdges: "
                 << faceEdges_[newFaceIndex[faceI]]
                 << endl;
        }

        Info << "Added cells: " << endl;

        forAll(newCellIndex, cellI)
        {
            Info << newCellIndex[cellI] << ":: "
                 << cells_[newCellIndex[cellI]]
                 << endl;
        }

        if (debug > 3)
        {
            writeVTK
            (
                Foam::name(eIndex)
              + "_afterSwap_"
              + Foam::name(numTriangulations) + "_"
              + Foam::name(triangulationIndex),
                newCellIndex
            );
        }
    }
}

// Routine to perform 3-2 or 2-2 swaps
void dynamicTopoFvMesh::swap32
(
    const label eIndex,
    const label triangulationIndex,
    const label numTriangulations,
    const labelListList& triangulations,
    const labelList& hullVertices,
    const labelList& hullFaces,
    const labelList& hullCells
)
{
    // A 2-2 / 3-2 swap performs the following operations:
    //      [1] Remove three faces surrounding edgeToCheck
    //      [2] Remove two (2-2 swap) or three(3-2 swap)
    //          cells surrounding edgeToCheck
    //      [3] Add one internal face
    //      [4] Add two new cells
    //      [5] If edgeToCheck is on a boundary,
    //          add two boundary faces and a boundary edge (2-2 swap)
    //      edgeToCheck is removed later by swap3DEdges
    //      Update faceEdges, edgeFaces and edgePoints information

    // Figure out which edge this is...
    label tIndex = self();

    // Obtain edge reference
    edge& edgeToCheck = edges_[eIndex];

    // Determine the patch this edge belongs to
    label edgePatch = whichEdgePatch(eIndex);

    // Determine the three faces to be removed
    FixedList<label,3> facesForRemoval;
    labelHashSet cellsForRemoval(3);

    forAll(facesForRemoval, faceI)
    {
        facesForRemoval[faceI]
            = hullFaces[triangulations[faceI][triangulationIndex]];

        label own = owner_[facesForRemoval[faceI]];
        label nei = neighbour_[facesForRemoval[faceI]];

        // Check and add cells as well
        if (!cellsForRemoval.found(own))
        {
            cellsForRemoval.insert(own);
        }

        if (!cellsForRemoval.found(nei) && nei != -1)
        {
            cellsForRemoval.insert(nei);
        }
    }

    labelList cellRemovalList = cellsForRemoval.toc();

    if (debug > 2)
    {
        // Print out arguments
        Info << endl;

        if (edgePatch < 0)
        {
            Info << "== Swapping 3-2 ==" << endl;
        }
        else
        {
            Info << "== Swapping 2-2 ==" << endl;
        }

        Info << "Edge: " << eIndex << ": " << edgeToCheck << endl;
        Info << "Ring: " << hullVertices << endl;
        Info << "Faces: " << hullFaces << endl;
        Info << "Cells: " << hullCells << endl;
        Info << "Triangulation: "
             << triangulations[0][triangulationIndex] << " "
             << triangulations[1][triangulationIndex] << " "
             << triangulations[2][triangulationIndex] << " "
             << endl;

        if (debug > 3)
        {
            writeVTK
            (
                Foam::name(eIndex)
              + "_beforeSwap_"
              + Foam::name(numTriangulations) + "_"
              + Foam::name(triangulationIndex),
                cellRemovalList
            );
        }
    }

    // Add two new cells to the end of the cell list
    FixedList<label,2> newCellIndex(-1);
    newCellIndex[0] = cells_.append(cell(4));
    newCellIndex[1] = cells_.append(cell(4));

    cell &newTetCell0 = cells_[newCellIndex[0]];
    cell &newTetCell1 = cells_[newCellIndex[1]];

    // Update length-scale info
    if (edgeModification_)
    {
        scalar avgScale = 0.0;

        forAll(cellRemovalList, cellI)
        {
            avgScale += lengthScale_[cellRemovalList[cellI]];
        }

        avgScale /= cellRemovalList.size();

        for (label i = 0; i < 2; i++)
        {
            lengthScale_.append(avgScale);
        }
    }

    // Insert a new internal face
    face newTriFace(3);

    newTriFace[0] = hullVertices[triangulations[0][triangulationIndex]];
    newTriFace[1] = hullVertices[triangulations[1][triangulationIndex]];
    newTriFace[2] = hullVertices[triangulations[2][triangulationIndex]];

    label newFaceIndex = insertFace
                         (
                             -1,
                             newTriFace,
                             newCellIndex[0],
                             newCellIndex[1]
                         );

    // Add faceEdges for the new face as well.
    faceEdges_.append(labelList(3));

    // Define the three edges to check while building faceEdges:
    FixedList<edge,3> check;

    check[0][0] = newTriFace[0]; check[0][1] = newTriFace[1];
    check[1][0] = newTriFace[1]; check[1][1] = newTriFace[2];
    check[2][0] = newTriFace[2]; check[2][1] = newTriFace[0];

    // New faceEdge entry for the interior face
    label nE = 0;
    labelList& newFaceEdges = faceEdges_[newFaceIndex];

    // For 2-2 swaps, two faces are introduced
    FixedList<label,2> nBE(0);
    labelListList bdyFaceEdges(2, labelList(3, -1));

    // Fill-in information for the two new cells,
    // and correct info on existing neighbouring cells
    label nF0 = 0, nF1 = 0;
    label otherPoint = -1, nextPoint = -1;
    FixedList<bool,2> foundEdge;

    // For a 2-2 swap on a boundary edge,
    // add two boundary faces and an edge
    FixedList<label,2> newBdyFaceIndex(-1);
    label newEdgeIndex = -1;

    if (edgePatch > -1)
    {
        // Temporary local variables
        label facePatch = -1;
        edge newEdge(-1, -1);
        FixedList<label,2> nBEdge(0);
        FixedList<FixedList<label,2>,2> bdyEdges;
        FixedList<face,2> newBdyTriFace(face(3));

        // Get a cue for face orientation from existing faces
        forAll(facesForRemoval, faceI)
        {
            if (neighbour_[facesForRemoval[faceI]] == -1)
            {
                facePatch = whichPatch(facesForRemoval[faceI]);
                face& faceToCheck = faces_[facesForRemoval[faceI]];

                findIsolatedPoint
                (
                    faceToCheck,
                    edgeToCheck,
                    otherPoint,
                    nextPoint
                );

                if (nextPoint == edgeToCheck[0])
                {
                    newEdge[1] = otherPoint;
                    newBdyTriFace[0][0] = otherPoint;
                    newBdyTriFace[0][1] = edgeToCheck[0];
                    newBdyTriFace[1][2] = otherPoint;
                }
                else
                {
                    newEdge[0] = otherPoint;
                    newBdyTriFace[1][0] = otherPoint;
                    newBdyTriFace[1][1] = edgeToCheck[1];
                    newBdyTriFace[0][2] = otherPoint;
                }

                // Also update faceEdges for the new boundary faces
                labelList& fEdges = faceEdges_[facesForRemoval[faceI]];

                forAll(fEdges, edgeI)
                {
                    if
                    (
                        edges_[fEdges[edgeI]]
                     == edge(edgeToCheck[0], otherPoint)
                    )
                    {
                        bdyFaceEdges[0][nBE[0]++] = fEdges[edgeI];
                        bdyEdges[0][nBEdge[0]++] = fEdges[edgeI];

                        edgeStack(tIndex).push(fEdges[edgeI]);
                    }

                    if
                    (
                        edges_[fEdges[edgeI]]
                     == edge(edgeToCheck[1], otherPoint)
                    )
                    {
                        bdyFaceEdges[1][nBE[1]++] = fEdges[edgeI];
                        bdyEdges[1][nBEdge[1]++] = fEdges[edgeI];

                        edgeStack(tIndex).push(fEdges[edgeI]);
                    }
                }
            }
        }

        // Insert the two new faces
        newBdyFaceIndex[0] = insertFace
                             (
                                 facePatch,
                                 newBdyTriFace[0],
                                 newCellIndex[1],
                                 -1
                             );

        newBdyFaceIndex[1] = insertFace
                             (
                                 facePatch,
                                 newBdyTriFace[1],
                                 newCellIndex[0],
                                 -1
                             );

        // Update the new cells
        newTetCell0[nF0++] = newBdyFaceIndex[1];
        newTetCell1[nF1++] = newBdyFaceIndex[0];

        // Add an edgeFaces entry
        labelList newBdyEdgeFaces(3, -1);
        newBdyEdgeFaces[0] = newBdyFaceIndex[0];
        newBdyEdgeFaces[1] = newFaceIndex;
        newBdyEdgeFaces[2] = newBdyFaceIndex[1];

        // Find the point other than the new edge
        // on the new triangular face
        findIsolatedPoint
        (
            newTriFace,
            newEdge,
            otherPoint,
            nextPoint
        );

        // Add an edgePoints entry
        labelList newBdyEdgePoints(3, -1);
        newBdyEdgePoints[0] = edgeToCheck[0];
        newBdyEdgePoints[1] = otherPoint;
        newBdyEdgePoints[2] = edgeToCheck[1];

        // Insert the edge
        newEdgeIndex = insertEdge
                       (
                           edgePatch,
                           newEdge,
                           newBdyEdgeFaces,
                           newBdyEdgePoints
                       );

        // Update faceEdges with the new edge
        newFaceEdges[nE++] = newEdgeIndex;
        bdyFaceEdges[0][nBE[0]++] = newEdgeIndex;
        bdyFaceEdges[1][nBE[1]++] = newEdgeIndex;

        // Update edgeFaces and edgePoints with the two new faces
        forAll(bdyEdges[0], edgeI)
        {
            sizeUpList(newBdyFaceIndex[0], edgeFaces_[bdyEdges[0][edgeI]]);
            sizeUpList(newBdyFaceIndex[1], edgeFaces_[bdyEdges[1][edgeI]]);

            // Replace the edgePoints label, and preserve position on the list
            findIsolatedPoint
            (
                newBdyTriFace[0],
                edges_[bdyEdges[0][edgeI]],
                otherPoint,
                nextPoint
            );

            replaceLabel
            (
                edgeToCheck[1],
                otherPoint,
                edgePoints_[bdyEdges[0][edgeI]]
            );

            // Size up edgePoints again, so that it is sized down later
            sizeUpList(edgeToCheck[1], edgePoints_[bdyEdges[0][edgeI]]);

            // Replace the edgePoints label, and preserve position on the list
            findIsolatedPoint
            (
                newBdyTriFace[1],
                edges_[bdyEdges[1][edgeI]],
                otherPoint,
                nextPoint
            );

            replaceLabel
            (
                edgeToCheck[0],
                otherPoint,
                edgePoints_[bdyEdges[1][edgeI]]
            );

            // Size up edgePoints again, so that it is sized down later
            sizeUpList(edgeToCheck[0], edgePoints_[bdyEdges[1][edgeI]]);
        }

        // Add faceEdges for the two new boundary faces
        faceEdges_.append(bdyFaceEdges[0]);
        faceEdges_.append(bdyFaceEdges[1]);
    }

    newTetCell0[nF0++] = newFaceIndex;
    newTetCell1[nF1++] = newFaceIndex;

    forAll(cellRemovalList, cellI)
    {
        label cellIndex = cellRemovalList[cellI];
        cell& cellToCheck = cells_[cellIndex];

        forAll(cellToCheck, faceI)
        {
            label faceIndex = cellToCheck[faceI];
            face& faceToCheck = faces_[faceIndex];

            foundEdge[0] = false; foundEdge[1] = false;

            // Find the face that contains only
            // edgeToCheck[0] or edgeToCheck[1]
            forAll(faceToCheck, pointI)
            {
                if (faceToCheck[pointI] == edgeToCheck[0])
                {
                    foundEdge[0] = true;
                }

                if (faceToCheck[pointI] == edgeToCheck[1])
                {
                    foundEdge[1] = true;
                }
            }

            // Face is connected to edgeToCheck[0]
            if (foundEdge[0] && !foundEdge[1])
            {
                // Check if a face-flip is necessary
                if (owner_[faceIndex] == cellIndex)
                {
                    if (neighbour_[faceIndex] == -1)
                    {
                        // Change the owner
                        owner_[faceIndex] = newCellIndex[1];
                    }
                    else
                    {
                        // Flip this face
                        faces_[faceIndex] = faceToCheck.reverseFace();
                        owner_[faceIndex] = neighbour_[faceIndex];
                        neighbour_[faceIndex] = newCellIndex[1];
                    }
                }
                else
                {
                    // Flip is unnecessary. Just update neighbour
                    neighbour_[faceIndex] = newCellIndex[1];
                }

                // Add this face to the cell
                newTetCell1[nF1++] = faceIndex;

                // Update faceEdges, edgeFaces and edgePoints
                labelList& fEdges = faceEdges_[faceIndex];

                forAll(fEdges, edgeI)
                {
                    edge& checkEdge = edges_[fEdges[edgeI]];

                    if
                    (
                        (checkEdge == check[0])
                     || (checkEdge == check[1])
                     || (checkEdge == check[2])
                    )
                    {
                        newFaceEdges[nE++] = fEdges[edgeI];

                        sizeUpList
                        (
                            newFaceIndex,
                            edgeFaces_[fEdges[edgeI]]
                        );

                        // Find the isolated point and insert it
                        findIsolatedPoint
                        (
                            newTriFace,
                            checkEdge,
                            otherPoint,
                            nextPoint
                        );

                        insertLabel
                        (
                            otherPoint,
                            edgeToCheck[0],
                            edgeToCheck[1],
                            edgePoints_[fEdges[edgeI]]
                        );

                        edgeStack(tIndex).push(fEdges[edgeI]);

                        break;
                    }
                }
            }

            // Face is connected to edgeToCheck[1]
            if (!foundEdge[0] && foundEdge[1])
            {
                // Check if a face-flip is necessary
                if (owner_[faceIndex] == cellIndex)
                {
                    if (neighbour_[faceIndex] == -1)
                    {
                        // Change the owner
                        owner_[faceIndex] = newCellIndex[0];
                    }
                    else
                    {
                        // Flip this face
                        faces_[faceIndex] = faceToCheck.reverseFace();
                        owner_[faceIndex] = neighbour_[faceIndex];
                        neighbour_[faceIndex] = newCellIndex[0];
                    }
                }
                else
                {
                    // Flip is unnecessary. Just update neighbour
                    neighbour_[faceIndex] = newCellIndex[0];
                }

                // Add this face to the cell
                newTetCell0[nF0++] = faceIndex;
            }
        }
    }

    // Generate mapping information for the two new cells
    const labelListList& cc = cellCells();
    labelHashSet masterObjects;
    FixedList<label,3> parents(-1);

    forAll(cellRemovalList, indexI)
    {
        // Determine an appropriate parent cell
        if (cellRemovalList[indexI] < nOldCells_)
        {
            parents[indexI] = cellRemovalList[indexI];
        }
        else
        {
            parents[indexI] = cellParents_[cellRemovalList[indexI]];
        }

        // Find the cell's neighbours in the old mesh
        masterObjects.insert(parents[indexI]);
        forAll(cc[parents[indexI]],cellI)
        {
            if (!masterObjects.found(cc[parents[indexI]][cellI]))
            {
                masterObjects.insert(cc[parents[indexI]][cellI]);
            }
        }
    }

    forAll(newCellIndex, cellI)
    {
        // Insert the parent cell [from first by default]
        cellParents_.insert(newCellIndex[cellI], parents[0]);

        // Insert mapping info into the HashTable
        cellsFromCells_.insert
        (
            newCellIndex[cellI],
            objectMap
            (
                newCellIndex[cellI],
                masterObjects.toc()
            )
        );
    }

    // Remove the faces and update associated edges
    forAll(facesForRemoval, faceI)
    {
        // Update edgeFaces and edgePoints
        face& checkFace = faces_[facesForRemoval[faceI]];
        labelList& fEdges = faceEdges_[facesForRemoval[faceI]];

        forAll(fEdges, edgeI)
        {
            label edgeIndex = fEdges[edgeI];

            if (edgeIndex != eIndex)
            {
                sizeDownList
                (
                    facesForRemoval[faceI],
                    edgeFaces_[edgeIndex]
                );

                // Find the isolated point and remove it
                findIsolatedPoint
                (
                    checkFace,
                    edges_[edgeIndex],
                    otherPoint,
                    nextPoint
                );

                sizeDownList
                (
                    otherPoint,
                    edgePoints_[edgeIndex]
                );

                edgeStack(tIndex).push(edgeIndex);
            }
        }

        // Now remove the face...
        removeFace(facesForRemoval[faceI]);
    }

    if (edgePatch < 0)
    {
        // Update the number of cells only for 3-2 swaps
        nCells_--;
    }

    forAll(cellRemovalList, cellI)
    {
        label cIndex = cellRemovalList[cellI];

        if (debug > 2)
        {
            Info << "Removing cell: "
                 << cIndex << ": "
                 << cells_[cIndex]
                 << endl;
        }

        cells_.remove(cIndex);

        if (edgeModification_)
        {
            lengthScale_.remove(cIndex);
        }

        if (cIndex < nOldCells_)
        {
            reverseCellMap_[cIndex] = -1;
        }

        // Check if the cell was added in the current morph, and delete
        if (cellsFromCells_.found(cIndex))
        {
            cellsFromCells_.erase(cIndex);
        }
    }

    if (debug > 2)
    {
        if (edgePatch > -1)
        {
            Info << "Added edge: " << endl;

            Info << newEdgeIndex << ":: "
                 << edges_[newEdgeIndex]
                 << " edgeFaces: "
                 << edgeFaces_[newEdgeIndex]
                 << endl;
        }

        Info << "Added face(s): " << endl;

        Info << newFaceIndex << ":: "
             << faces_[newFaceIndex];

        Info << " faceEdges: "
             << faceEdges_[newFaceIndex]
             << endl;

        if (edgePatch > -1)
        {
            forAll(newBdyFaceIndex, faceI)
            {
                Info << newBdyFaceIndex[faceI] << ":: "
                     << faces_[newBdyFaceIndex[faceI]]
                     << " faceEdges: "
                     << faceEdges_[newBdyFaceIndex[faceI]]
                     << endl;
            }
        }

        Info << "Added cells: " << endl;

        forAll(newCellIndex, cellI)
        {
            Info << newCellIndex[cellI] << ":: "
                 << cells_[newCellIndex[cellI]]
                 << endl;
        }

        if (debug > 3)
        {
            writeVTK
            (
                Foam::name(eIndex)
              + "_afterSwap_"
              + Foam::name(numTriangulations) + "_"
              + Foam::name(triangulationIndex),
                newCellIndex
            );
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
