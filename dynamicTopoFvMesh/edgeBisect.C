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

// Method for the bisection of a quad-face in 2D
void dynamicTopoFvMesh::bisectQuadFace
(
    const label fIndex
)
{
    // Quad-face bisection performs the following operations:
    //      [1] Add two points at the middle of the face
    //      [2] Create a new internal face for each bisected cell
    //      [3] Modify existing face and create a new half-face
    //      [4] Modify triangular boundary faces, and create new ones as well
    //      [5] Create edges for new faces
    //      Update faceEdges and edgeFaces information

    bool found;
    label replaceFace = -1;
    face tmpQuadFace(4), tmpTriFace(3);
    FixedList<label,7> newFaceIndex(-1), newEdgeIndex(-1);
    FixedList<edge,2> commonEdges;
    FixedList<label,4> modifiedEdgeIndex(-1);
    FixedList<label,2> commonEdgeIndex(-1), commonFaceIndex(-1);
    FixedList<label,2> newPtIndex(-1), newCellIndex(-1), otherEdgePoint(-1);
    FixedList<label,4> otherEdgeIndex(-1);
    FixedList<label,4> otherPointIndex(-1), nextToOtherPoint(-1);
    FixedList<label,2> c0BdyIndex, c0IntIndex, c1BdyIndex, c1IntIndex;
    FixedList<face,2>  c0BdyFace,  c0IntFace,  c1BdyFace,  c1IntFace;

    // Obtain a reference for this face...
    face& thisFace = faces_[fIndex];

    // Get the two cells on either side...
    label c0 = owner_[fIndex], c1 = neighbour_[fIndex];

    // Find the prism faces for cell[0].
    cell &cell_0 = cells_[c0];

    findPrismFaces
    (
        fIndex,
        c0,
        c0BdyFace,
        c0BdyIndex,
        c0IntFace,
        c0IntIndex
    );

    if (debug > 1)
    {
        Info << nl << nl << "Face: " << fIndex
             << ": " << thisFace << " is to be bisected. " << endl;

        if (debug > 2)
        {
            Info << "Cell[0]: " << c0 << ": " << cell_0 << endl;

            forAll(cell_0, faceI)
            {
                Info << cell_0[faceI] << ": "
                     << faces_[cell_0[faceI]] << endl;
            }
        }
    }

    // Find the common-edge between the triangular boundary faces
    // and the face under consideration.
    findCommonEdge(c0BdyIndex[0], fIndex, commonEdgeIndex[0]);
    findCommonEdge(c0BdyIndex[1], fIndex, commonEdgeIndex[1]);

    commonEdges[0] = edges_[commonEdgeIndex[0]];
    commonEdges[1] = edges_[commonEdgeIndex[1]];

    // Find the isolated point on both boundary faces of cell[0]
    findIsolatedPoint
    (
        c0BdyFace[0],
        commonEdges[0],
        otherPointIndex[0],
        nextToOtherPoint[0]
    );

    findIsolatedPoint
    (
        c0BdyFace[1],
        commonEdges[1],
        otherPointIndex[1],
        nextToOtherPoint[1]
    );

    // For convenience...
    otherEdgePoint[0] = commonEdges[0].otherVertex(nextToOtherPoint[0]);
    otherEdgePoint[1] = commonEdges[1].otherVertex(nextToOtherPoint[1]);

    // Add two new points to the end of the list
    newPtIndex[0] = points_.append
                    (
                        0.5*
                        (
                            points_[commonEdges[0][0]]
                          + points_[commonEdges[0][1]]
                        )
                    );

    newPtIndex[1] = points_.append
                    (
                        0.5*
                        (
                            points_[commonEdges[1][0]]
                          + points_[commonEdges[1][1]]
                        )
                    );
    nPoints_ += 2;

    // Add a new prism cell to the end of the list
    newCellIndex[0] = cells_.append(cell(5));
    cell &newCell0 = cells_[newCellIndex[0]];

    // Generate mapping information for this new cell
    label firstParent;
    const labelListList& cc = cellCells();
    labelHashSet c0MasterObjects;

    if (c0 < nOldCells_)
    {
        firstParent = c0;
    }
    else
    {
        firstParent = cellParents_[c0];
    }

    // Insert the parent cell
    cellParents_.insert(newCellIndex[0],firstParent);

    // Find the cell's neighbours in the old mesh
    c0MasterObjects.insert(firstParent);
    forAll(cc[firstParent],cellI)
    {
        if (!c0MasterObjects.found(cc[firstParent][cellI]))
        {
            c0MasterObjects.insert(cc[firstParent][cellI]);
        }
    }

    // Insert mapping info into the HashTable
    cellsFromCells_.insert
    (
        newCellIndex[0],
        objectMap
        (
            newCellIndex[0],
            c0MasterObjects.toc()
        )
    );

    // Add a new element to the lengthScale field
    // (Currently the same as cell[0])
    lengthScale_.append(lengthScale_[c0]);

    // Modify the two existing triangle boundary faces

    // Zeroth boundary face - Owner = c[0] & Neighbour [-1] (unchanged)
    replaceLabel
    (
        otherEdgePoint[0],
        newPtIndex[0],
        c0BdyFace[0]
    );

    // First boundary face - Owner = newCell[0], Neighbour = -1
    replaceLabel
    (
        otherEdgePoint[1],
        newPtIndex[1],
        c0BdyFace[1]
    );

    owner_[c0BdyIndex[1]] = newCellIndex[0];
    replaceLabel(c0BdyIndex[1],-1,cell_0);

    // Detect edges other than commonEdges
    labelList& fEdges = faceEdges_[fIndex];

    forAll(fEdges, edgeI)
    {
        if
        (
            fEdges[edgeI] != commonEdgeIndex[0] &&
            fEdges[edgeI] != commonEdgeIndex[1]
        )
        {
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

    // Modify point-labels on the quad face under consideration
    replaceLabel
    (
        otherEdgePoint[0],
        newPtIndex[0],
        thisFace
    );

    replaceLabel
    (
        nextToOtherPoint[1],
        newPtIndex[1],
        thisFace
    );

    if (debug > 1)
    {
        Info << "Modified thisFace: " << fIndex
             << ": " << thisFace << endl;
    }

    // Find the interior face that contains otherEdgeIndex[1]
    found = false;

    labelList& e1 = faceEdges_[c0IntIndex[0]];

    forAll(e1, edgeI)
    {
        if (e1[edgeI] == otherEdgeIndex[1])
        {
            replaceLabel(c0IntIndex[0], -1, cell_0);
            replaceFace = c0IntIndex[0];
            found = true; break;
        }
    }

    if (!found)
    {
        // The edge was not found before
        replaceLabel(c0IntIndex[1], -1, cell_0);
        replaceFace = c0IntIndex[1];
    }

    // Check if face reversal is necessary for the replacement
    if (owner_[replaceFace] == c0)
    {
        if (neighbour_[replaceFace] == -1)
        {
            // Change the owner
            owner_[replaceFace] = newCellIndex[0];
        }
        else
        {
            // This face has to be reversed
            faces_[replaceFace] = faces_[replaceFace].reverseFace();
            owner_[replaceFace] = neighbour_[replaceFace];
            neighbour_[replaceFace] = newCellIndex[0];
        }
    }
    else
    {
        // Keep owner, but change neighbour
        neighbour_[replaceFace] = newCellIndex[0];
    }

    // Define the faces for the new cell
    newCell0[0] = c0BdyIndex[1];
    newCell0[1] = replaceFace;

    // Define the set of new faces and insert...

    // New interior face; Owner = cell[0] & Neighbour = newCell[0]
    tmpQuadFace[0] = otherPointIndex[0];
    tmpQuadFace[1] = newPtIndex[0];
    tmpQuadFace[2] = newPtIndex[1];
    tmpQuadFace[3] = otherPointIndex[1];

    newFaceIndex[0] = insertFace
                      (
                          -1,
                          tmpQuadFace,
                          c0,
                          newCellIndex[0]
                      );

    replaceLabel(-1, newFaceIndex[0], newCell0);
    replaceLabel(-1, newFaceIndex[0], cell_0);

    // remove2DSliver requires this face index for removal
    bisectInterior_ = newFaceIndex[0];

    // Second boundary face; Owner = newCell[0] & Neighbour = [-1]
    tmpTriFace[0] = otherPointIndex[0];
    tmpTriFace[1] = newPtIndex[0];
    tmpTriFace[2] = otherEdgePoint[0];

    newFaceIndex[1] = insertFace
                      (
                          whichPatch(c0BdyIndex[0]),
                          tmpTriFace,
                          newCellIndex[0],
                          -1
                      );

    replaceLabel(-1, newFaceIndex[1], newCell0);

    // Third boundary face; Owner = c[0] & Neighbour = [-1]
    tmpTriFace[0] = otherPointIndex[1];
    tmpTriFace[1] = newPtIndex[1];
    tmpTriFace[2] = otherEdgePoint[1];

    newFaceIndex[2] = insertFace
                      (
                          whichPatch(c0BdyIndex[1]),
                          tmpTriFace,
                          c0,
                          -1
                      );

    replaceLabel(-1, newFaceIndex[2], cell_0);

    // Create / modify edges...
    labelList tmpTriEdgeFaces(3, -1);

    // The edge bisecting the zeroth boundary triangular face
    tmpTriEdgeFaces[0] = c0BdyIndex[0];
    tmpTriEdgeFaces[1] = newFaceIndex[0];
    tmpTriEdgeFaces[2] = newFaceIndex[1];

    newEdgeIndex[1] = insertEdge
                      (
                          whichEdgePatch(commonEdgeIndex[0]),
                          edge(newPtIndex[0], otherPointIndex[0]),
                          tmpTriEdgeFaces
                      );

    // The edge bisecting the first boundary triangular face
    tmpTriEdgeFaces[0] = c0BdyIndex[1];
    tmpTriEdgeFaces[1] = newFaceIndex[0];
    tmpTriEdgeFaces[2] = newFaceIndex[2];

    newEdgeIndex[2] = insertEdge
                      (
                          whichEdgePatch(commonEdgeIndex[1]),
                          edge(newPtIndex[1], otherPointIndex[1]),
                          tmpTriEdgeFaces
                      );

    if (c1 == -1)
    {
        // The quad boundary face resulting from bisection;
        // Owner = newCell[0] & Neighbour = [-1]
        tmpQuadFace[0] = newPtIndex[1];
        tmpQuadFace[1] = nextToOtherPoint[1];
        tmpQuadFace[2] = otherEdgePoint[0];
        tmpQuadFace[3] = newPtIndex[0];

        newFaceIndex[3] = insertFace
                          (
                              whichPatch(fIndex),
                              tmpQuadFace,
                              newCellIndex[0],
                              -1
                          );

        // Add a faceEdges entry as well
        faceEdges_.append(labelList(4, -1));

        replaceLabel(-1, newFaceIndex[3], newCell0);

        labelList tmpBiEdgeFaces(2, -1);

        // The edge bisecting the face
        tmpTriEdgeFaces[0] = fIndex;
        tmpTriEdgeFaces[1] = newFaceIndex[0];
        tmpTriEdgeFaces[2] = newFaceIndex[3];

        newEdgeIndex[0] = insertEdge
                          (
                              whichEdgePatch(otherEdgeIndex[0]),
                              edge(newPtIndex[0], newPtIndex[1]),
                              tmpTriEdgeFaces
                          );

        // Create / replace side edges created from face bisection
        tmpBiEdgeFaces[0] = newFaceIndex[1];
        tmpBiEdgeFaces[1] = newFaceIndex[3];

        newEdgeIndex[3] = insertEdge
                          (
                              whichEdgePatch(commonEdgeIndex[0]),
                              edge(newPtIndex[0], otherEdgePoint[0]),
                              tmpBiEdgeFaces
                          );

        tmpBiEdgeFaces[0] = c0BdyIndex[1];
        tmpBiEdgeFaces[1] = newFaceIndex[3];

        newEdgeIndex[4] = insertEdge
                          (
                              whichEdgePatch(commonEdgeIndex[1]),
                              edge(newPtIndex[1], nextToOtherPoint[1]),
                              tmpBiEdgeFaces
                          );

        // Replace an edge on the bisected face
        replaceLabel
        (
            otherEdgeIndex[1],
            newEdgeIndex[0],
            faceEdges_[fIndex]
        );

        // Now that edges are defined, configure faceEdges
        replaceLabel
        (
            commonEdgeIndex[1],
            newEdgeIndex[4],
            faceEdges_[c0BdyIndex[1]]
        );

        if (debug > 2)
        {
            Info << "Modified Cell[0]: "
                 << c0 << ": " << cell_0 << endl;

            forAll(cell_0, faceI)
            {
                Info << cell_0[faceI]
                     << ": " << faces_[cell_0[faceI]] << endl;
            }

            Info << "New Cell[0]: " << newCellIndex[0]
                 << ": " << newCell0 << endl;

            forAll(newCell0, faceI)
            {
                Info << newCell0[faceI]
                     << ": " << faces_[newCell0[faceI]] << endl;
            }
        }
    }
    else
    {
        cell &cell_1 = cells_[c1];

        // Find the prism faces for cell[1].
        findPrismFaces
        (
            fIndex,
            c1,
            c1BdyFace,
            c1BdyIndex,
            c1IntFace,
            c1IntIndex
        );

        // Add a new prism cell to the end of the list
        newCellIndex[1] = cells_.append(cell(5));
        cell &newCell1 = cells_[newCellIndex[1]];

        // Generate mapping information for this new cell
        label secondParent;
        labelHashSet c1MasterObjects;

        if (c1 < nOldCells_)
        {
            secondParent = c1;
        }
        else
        {
            secondParent = cellParents_[c1];
        }

        // Insert the parent cell
        cellParents_.insert(newCellIndex[1],secondParent);

        // Find the cell's neighbours in the old mesh
        c1MasterObjects.insert(secondParent);
        forAll(cc[secondParent],cellI)
        {
            if (!c1MasterObjects.found(cc[secondParent][cellI]))
            {
                c1MasterObjects.insert(cc[secondParent][cellI]);
            }
        }

        // Insert mapping info into the HashTable
        cellsFromCells_.insert
        (
            newCellIndex[1],
            objectMap
            (
                newCellIndex[1],
                c1MasterObjects.toc()
            )
        );

        // Add a new element to the lengthScale field
        // (Currently the same as cell[1])
        lengthScale_.append(lengthScale_[c1]);

        if (debug > 2)
        {
            Info << "Cell[1]: " << c1 << ": " << cell_1 << endl;

            forAll(cell_1, faceI)
            {
                Info << cell_1[faceI] << ": "
                     << faces_[cell_1[faceI]] << endl;
            }
        }

        // Find the interior face that contains otherEdgeIndex[1]
        found = false;

        labelList& e2 = faceEdges_[c1IntIndex[0]];

        forAll(e2, edgeI)
        {
            if (e2[edgeI] == otherEdgeIndex[1])
            {
                replaceLabel(c1IntIndex[0], -1, cell_1);
                replaceFace = c1IntIndex[0];
                found = true; break;
            }
        }

        if (!found)
        {
            // The edge was not found before
            replaceLabel(c1IntIndex[1], -1, cell_1);
            replaceFace = c1IntIndex[1];
        }

        // Check if face reversal is necessary for the replacement
        if (owner_[replaceFace] == c1)
        {
            if (neighbour_[replaceFace] == -1)
            {
                // Change the owner
                owner_[replaceFace] = newCellIndex[1];
            }
            else
            {
                // This face has to be reversed
                faces_[replaceFace] = faces_[replaceFace].reverseFace();
                owner_[replaceFace] = neighbour_[replaceFace];
                neighbour_[replaceFace] = newCellIndex[1];
            }
        }
        else
        {
            // Keep owner, but change neighbour
            neighbour_[replaceFace] = newCellIndex[1];
        }

        // Define attributes for the new prism cell
        newCell1[0] = replaceFace;

        // The interior face resulting from bisection;
        // Owner = newCell[0] & Neighbour = newCell[1]
        tmpQuadFace[0] = newPtIndex[1];
        tmpQuadFace[1] = nextToOtherPoint[1];
        tmpQuadFace[2] = otherEdgePoint[0];
        tmpQuadFace[3] = newPtIndex[0];

        newFaceIndex[3] = insertFace
                          (
                              -1,
                              tmpQuadFace,
                              newCellIndex[0],
                              newCellIndex[1]
                          );

        replaceLabel(-1, newFaceIndex[3], newCell0);
        replaceLabel(-1, newFaceIndex[3], newCell1);
        newCell1[1] = newFaceIndex[3];

        // Check for common edges among the two boundary faces
        // Find the isolated point on both boundary faces of cell[1]
        if
        (
            findCommonEdge(c1BdyIndex[0], c0BdyIndex[0], commonEdgeIndex[0])
        )
        {
            findCommonEdge(c1BdyIndex[1], c0BdyIndex[1], commonEdgeIndex[1]);

            commonFaceIndex[0] = c1BdyIndex[0];
            commonFaceIndex[1] = c1BdyIndex[1];
        }
        else
        {
            findCommonEdge(c1BdyIndex[0], c0BdyIndex[1], commonEdgeIndex[1]);
            findCommonEdge(c1BdyIndex[1], c0BdyIndex[0], commonEdgeIndex[0]);

            commonFaceIndex[0] = c1BdyIndex[1];
            commonFaceIndex[1] = c1BdyIndex[0];
        }

        commonEdges[0] = edges_[commonEdgeIndex[0]];
        commonEdges[1] = edges_[commonEdgeIndex[1]];

        findIsolatedPoint
        (
            faces_[commonFaceIndex[0]],
            commonEdges[0],
            otherPointIndex[2],
            nextToOtherPoint[2]
        );

        findIsolatedPoint
        (
            faces_[commonFaceIndex[1]],
            commonEdges[1],
            otherPointIndex[3],
            nextToOtherPoint[3]
        );

        // For convenience...
        otherEdgePoint[0] = commonEdges[0].otherVertex(nextToOtherPoint[2]);
        otherEdgePoint[1] = commonEdges[1].otherVertex(nextToOtherPoint[3]);

        // Modify the two existing triangle boundary faces
        // Zeroth boundary face - Owner = newCell[1], Neighbour = -1
        replaceLabel
        (
            otherEdgePoint[0],
            newPtIndex[0],
            faces_[commonFaceIndex[0]]
        );

        owner_[commonFaceIndex[0]] = newCellIndex[1];
        replaceLabel(commonFaceIndex[0], -1, cell_1);
        newCell1[2] = commonFaceIndex[0];

        // First boundary face - Owner = c[1] & Neighbour [-1] (unchanged)
        replaceLabel
        (
            otherEdgePoint[1],
            newPtIndex[1],
            faces_[commonFaceIndex[1]]
        );

        // New interior face; Owner = cell[1] & Neighbour = newCell[1]
        tmpQuadFace[0] = newPtIndex[0];
        tmpQuadFace[1] = otherPointIndex[2];
        tmpQuadFace[2] = otherPointIndex[3];
        tmpQuadFace[3] = newPtIndex[1];

        newFaceIndex[4] = insertFace
                          (
                              -1,
                              tmpQuadFace,
                              c1,
                              newCellIndex[1]
                          );

        replaceLabel(-1, newFaceIndex[4], newCell1);
        replaceLabel(-1, newFaceIndex[4], cell_1);

        // Second boundary face; Owner = cell[1] & Neighbour [-1]
        tmpTriFace[0] = otherPointIndex[2];
        tmpTriFace[1] = newPtIndex[0];
        tmpTriFace[2] = otherEdgePoint[0];

        newFaceIndex[5] = insertFace
                          (
                              whichPatch(commonFaceIndex[0]),
                              tmpTriFace,
                              c1,
                              -1
                          );

        replaceLabel(-1, newFaceIndex[5], cell_1);

        // Third boundary face; Owner = newCell[1] & Neighbour [-1]
        tmpTriFace[0] = otherPointIndex[3];
        tmpTriFace[1] = newPtIndex[1];
        tmpTriFace[2] = otherEdgePoint[1];

        newFaceIndex[6] = insertFace
                          (
                              whichPatch(commonFaceIndex[1]),
                              tmpTriFace,
                              newCellIndex[1],
                              -1
                          );

        replaceLabel(-1, newFaceIndex[6], newCell1);

        // Create/modify edges...
        labelList tmpQuadEdgeFaces(4, -1);

        // The internal edge bisecting the face
        tmpQuadEdgeFaces[0] = fIndex;
        tmpQuadEdgeFaces[1] = newFaceIndex[0];
        tmpQuadEdgeFaces[2] = newFaceIndex[3];
        tmpQuadEdgeFaces[3] = newFaceIndex[4];

        newEdgeIndex[0] = insertEdge
                          (
                              -1,
                              edge(newPtIndex[0], newPtIndex[1]),
                              tmpQuadEdgeFaces
                          );

        // Create / replace side edges created from face bisection
        tmpTriEdgeFaces[0] = commonFaceIndex[0];
        tmpTriEdgeFaces[1] = newFaceIndex[3];
        tmpTriEdgeFaces[2] = newFaceIndex[1];

        newEdgeIndex[3] = insertEdge
                          (
                              whichEdgePatch(commonEdgeIndex[0]),
                              edge(newPtIndex[0], nextToOtherPoint[0]),
                              tmpTriEdgeFaces
                          );

        tmpTriEdgeFaces[0] = c0BdyIndex[1];
        tmpTriEdgeFaces[1] = newFaceIndex[3];
        tmpTriEdgeFaces[2] = newFaceIndex[6];

        newEdgeIndex[4] = insertEdge
                          (
                              whichEdgePatch(commonEdgeIndex[1]),
                              edge(newPtIndex[1], otherEdgePoint[1]),
                              tmpTriEdgeFaces
                          );

        // The edge bisecting the second boundary triangular face
        tmpTriEdgeFaces[0] = commonFaceIndex[0];
        tmpTriEdgeFaces[1] = newFaceIndex[4];
        tmpTriEdgeFaces[2] = newFaceIndex[5];

        newEdgeIndex[5] = insertEdge
                          (
                              whichEdgePatch(commonEdgeIndex[0]),
                              edge(newPtIndex[0], otherPointIndex[2]),
                              tmpTriEdgeFaces
                          );

        // The edge bisecting the third boundary triangular face
        tmpTriEdgeFaces[0] = commonFaceIndex[1];
        tmpTriEdgeFaces[1] = newFaceIndex[4];
        tmpTriEdgeFaces[2] = newFaceIndex[6];

        newEdgeIndex[6] = insertEdge
                          (
                              whichEdgePatch(commonEdgeIndex[1]),
                              edge(newPtIndex[1], otherPointIndex[3]),
                              tmpTriEdgeFaces
                          );

        if (debug > 2)
        {
            Info << nl << "Modified Cell[0]: "
                 << c0 << ": " << cell_0 << endl;

            forAll(cell_0, faceI)
            {
                Info << cell_0[faceI]
                     << ": " << faces_[cell_0[faceI]] << endl;
            }

            Info << "New Cell[0]: "
                 << newCellIndex[0] << ": " << newCell0 << endl;

            forAll(newCell0, faceI)
            {
                Info << newCell0[faceI] << ": "
                     << faces_[newCell0[faceI]] << endl;
            }

            Info << nl << "Modified Cell[1]: "
                 << c1 << ": " << cell_1 << endl;

            forAll(cell_1, faceI)
            {
                Info << cell_1[faceI] << ": "
                     << faces_[cell_1[faceI]] << endl;
            }

            Info << "New Cell[1]: "
                 << newCellIndex[1] << ": " << newCell1 << endl;

            forAll(newCell1, faceI)
            {
                Info << newCell1[faceI] << ": "
                     << faces_[newCell1[faceI]] << endl;
            }
        }
    }

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    nBisections_++;

    // Increment the number of modifications
    nModifications_++;

    // Update the number of cells
    nCells_++;

    if (c1 != -1)
    {
        nCells_++;
    }
}

// Method for the bisection of an edge in 3D
// - Returns a changeMap with a type specifying:
//    -1: Bisection failed since max number of topo-changes was reached.
// - AddedPoint is the index of the newly added point.
const dynamicTopoFvMesh::changeMap
dynamicTopoFvMesh::bisectEdge
(
    const label eIndex
)
{
    // Edge bisection performs the following operations:
    //      [1] Add a point at middle of the edge
    //      [2] Bisect all faces surrounding this edge
    //      [3] Bisect all cells surrounding this edge
    //      [4] Create internal/external edges for each bisected face
    //      [5] Create internal faces for each bisected cell
    //      Update faceEdges, edgeFaces and edgePoints information

    // Figure out which thread this is...
    label tIndex = self(), pIndex = -1;

    // Prepare the changeMaps
    changeMap map, coupleMap;
    bool bisectingSlave = false;

    if
    (
        (nModifications_ > maxModifications_) &&
        (maxModifications_ > -1)
    )
    {
        // Reached the max allowable topo-changes.
        edgeStack(tIndex).clear();

        return map;
    }

    if (coupledModification_)
    {
        // Is this a locally coupled edge?
        if (locallyCoupledEdge(eIndex))
        {
            // Temporarily turn off coupledModification
            coupledModification_ = false;

            // Loop through master/slave maps
            // and determine the coupled edge index.
            forAll(masterToSlave_, indexI)
            {
                if (masterToSlave_[indexI].found(eIndex))
                {
                    // Bisect the slave.
                    coupleMap = bisectEdge(masterToSlave_[indexI][eIndex]);

                    // Keep this index for master/slave mapping.
                    pIndex = indexI;
                }

                // The following bit happens only during the sliver
                // exudation process, since slave edges are
                // usually not added to the coupled edge-stack.
                if (slaveToMaster_[indexI].found(eIndex))
                {
                    // Bisect the master.
                    coupleMap = bisectEdge(slaveToMaster_[indexI][eIndex]);

                    // Keep this index for master/slave mapping.
                    pIndex = indexI;

                    // Notice that we are bisecting a slave edge.
                    bisectingSlave = true;
                }
            }

            // Turn it back on.
            coupledModification_ = true;
        }
        else
        {
            // Bisect edge on the patchSubMesh.

        }
    }

    // Hull variables
    face tmpTriFace(3);
    labelList tmpEdgeFaces(3,-1);
    labelList tmpIntEdgeFaces(4,-1);
    labelList tmpEdgePoints(3,-1);
    labelList tmpIntEdgePoints(4,-1);
    labelList tmpFaceEdges(3,-1);
    edge& thisEdge = edges_[eIndex];
    labelList& vertexHull = edgePoints_[eIndex];
    label m = vertexHull.size();

    // Size up the hull lists
    labelList cellHull(m, -1);
    labelList faceHull(m, -1);
    labelList edgeHull(m, -1);
    labelListList ringEntities(4, labelList(m, -1));

    // Construct a hull around this edge
    constructHull
    (
        eIndex,
        edgeHull,
        faceHull,
        cellHull,
        ringEntities
    );

    if (debug > 1)
    {
        Info << nl << nl << "Edge: " << eIndex
             << ": " << thisEdge << " is to be bisected. " << endl;

        // Write out VTK files prior to change
        if (debug > 3)
        {
            writeVTK
            (
                Foam::name(eIndex)
              + "Bisect_0",
                cellHull
            );
        }
    }

    // Add a new point to the end of the list
    label newPointIndex =
    (
        points_.append
        (
            0.5*
            (
                points_[thisEdge[0]]
              + points_[thisEdge[1]]
            )
        )
    );

    // Add an entry to pointEdges as well
    pointEdges_.append(labelList(0));

    nPoints_++;

    // Add this point to the map.
    map.addPoint(newPointIndex);

    // Add a new edge to the end of the list
    label newEdgeIndex =
    (
        insertEdge
        (
            whichEdgePatch(eIndex),
            edge(newPointIndex,thisEdge[1]),
            labelList(faceHull.size(),-1),
            vertexHull
        )
    );

    // Add this edge to the map.
    map.addEdge(newEdgeIndex);

    // Remove the existing edge from the pointEdges list
    // of the modified point, and add it to the new point
    sizeDownList(eIndex, pointEdges_[thisEdge[1]]);
    sizeUpList(eIndex, pointEdges_[newPointIndex]);

    // Modify the existing edge
    thisEdge[1] = newPointIndex;

    // Obtain new references
    edge& newEdge = edges_[newEdgeIndex];
    labelList& newEdgeFaces = edgeFaces_[newEdgeIndex];

    // Keep track of added entities
    labelList addedCellIndices(cellHull.size(),-1);
    labelList addedFaceIndices(faceHull.size(),-1);
    labelList addedEdgeIndices(faceHull.size(),-1);
    labelList addedIntFaceIndices(faceHull.size(),-1);

    // Obtain cellCells for mapping information
    const labelListList& cc = cellCells();

    // Now loop through the hull and bisect individual entities
    forAll(vertexHull, indexI)
    {
        // Fetch the existing face
        face& currFace = faces_[faceHull[indexI]];

        // Modify the existing face
        replaceLabel
        (
            newEdge[1],
            newPointIndex,
            currFace
        );

        // Modify edgePoints for the edge connected to thisEdge[0]
        replaceLabel
        (
            newEdge[1],
            newPointIndex,
            edgePoints_[ringEntities[0][indexI]]
        );

        // Obtain circular indices
        label nextI = vertexHull.fcIndex(indexI);
        label prevI = vertexHull.rcIndex(indexI);

        // Check if this is an interior/boundary face
        if (cellHull[indexI] != -1)
        {
            cell& currCell = cells_[cellHull[indexI]];

            // Create a new cell
            addedCellIndices[indexI] = cells_.append(cell(4));
            cell& newCell = cells_[addedCellIndices[indexI]];
            nCells_++;

            // Add this cell to the map.
            map.addCell(addedCellIndices[indexI]);

            // Generate mapping information for this new cell
            label parent;
            labelHashSet masterObjects;

            if (cellHull[indexI] < nOldCells_)
            {
                parent = cellHull[indexI];
            }
            else
            {
                parent = cellParents_[cellHull[indexI]];
            }

            // Insert the parent cell
            cellParents_.insert(addedCellIndices[indexI], parent);

            // Find the cell's neighbours in the old mesh
            masterObjects.insert(parent);
            forAll(cc[parent], cellI)
            {
                if (!masterObjects.found(cc[parent][cellI]))
                {
                    masterObjects.insert(cc[parent][cellI]);
                }
            }

            // Insert mapping info into the HashTable
            cellsFromCells_.insert
            (
                addedCellIndices[indexI],
                objectMap
                (
                    addedCellIndices[indexI],
                    masterObjects.toc()
                )
            );

            // Add a new element to the lengthScale field
            lengthScale_.append(lengthScale_[cellHull[indexI]]);

            // Configure the interior face
            tmpTriFace[0] = vertexHull[nextI];
            tmpTriFace[1] = vertexHull[indexI];
            tmpTriFace[2] = newPointIndex;

            // Insert the face
            addedIntFaceIndices[indexI] =
            (
                insertFace
                (
                    -1,
                    tmpTriFace,
                    cellHull[indexI],
                    addedCellIndices[indexI]
                )
            );

            // Add a faceEdges entry as well
            faceEdges_.append(tmpFaceEdges);

            // Add this face to the map.
            map.addFace(addedIntFaceIndices[indexI]);

            // Add to the new cell
            newCell[0] = addedIntFaceIndices[indexI];

            // Modify the existing ring face connected to newEdge[1]
            label replaceFace = ringEntities[3][indexI];

            // Check if face reversal is necessary
            if (owner_[replaceFace] == cellHull[indexI])
            {
                if (neighbour_[replaceFace] == -1)
                {
                    // Change the owner
                    owner_[replaceFace] = addedCellIndices[indexI];
                }
                else
                {
                    // This face has to be reversed
                    faces_[replaceFace] = faces_[replaceFace].reverseFace();
                    owner_[replaceFace] = neighbour_[replaceFace];
                    neighbour_[replaceFace] = addedCellIndices[indexI];
                }
            }
            else
            {
                // Keep owner, but change neighbour
                neighbour_[replaceFace] = addedCellIndices[indexI];
            }

            // Modify the edge on the ring.
            // Add the new interior face to edgeFaces.
            sizeUpList
            (
                addedIntFaceIndices[indexI],
                edgeFaces_[edgeHull[indexI]]
            );

            // Insert the new point to edgePoints for the ring edge
            insertLabel
            (
                newPointIndex,
                thisEdge[0],
                newEdge[1],
                edgePoints_[edgeHull[indexI]]
            );

            // Add this edge to faceEdges for the new interior face
            faceEdges_[addedIntFaceIndices[indexI]][0] = edgeHull[indexI];

            // Replace face labels
            replaceLabel
            (
                replaceFace,
                addedIntFaceIndices[indexI],
                currCell
            );

            // Add to the new cell
            newCell[1] = replaceFace;

            // Check if this is a boundary face
            if (cellHull[prevI] == -1)
            {
                // Configure the boundary face
                tmpTriFace[0] = newPointIndex;
                tmpTriFace[1] = newEdge[1];
                tmpTriFace[2] = vertexHull[indexI];

                // Insert the face
                addedFaceIndices[indexI] =
                (
                    insertFace
                    (
                        whichPatch(faceHull[indexI]),
                        tmpTriFace,
                        addedCellIndices[indexI],
                        -1
                    )
                );

                // Add this face to the map.
                map.addFace(addedFaceIndices[indexI]);

                // Generate mapping information for this new face
                if (faceHull[indexI] < nOldFaces_)
                {
                    parent = faceHull[indexI];
                }
                else
                {
                    parent = faceParents_[faceHull[indexI]];
                }

                // Insert the parent face
                faceParents_.insert(addedFaceIndices[indexI], parent);

                // Configure edgeFaces
                tmpEdgeFaces[0] = faceHull[indexI];
                tmpEdgeFaces[1] = addedIntFaceIndices[indexI];
                tmpEdgeFaces[2] = addedFaceIndices[indexI];

                // Configure edgePoints
                tmpEdgePoints[0] = thisEdge[0];
                tmpEdgePoints[1] = vertexHull[nextI];
                tmpEdgePoints[2] = newEdge[1];

                // Add an edge
                addedEdgeIndices[indexI] =
                (
                    insertEdge
                    (
                        whichPatch(faceHull[indexI]),
                        edge(newPointIndex,vertexHull[indexI]),
                        tmpEdgeFaces,
                        tmpEdgePoints
                    )
                );

                // Add this edge to the map.
                map.addEdge(addedEdgeIndices[indexI]);

                // Add this edge to the interior-face faceEdges entry
                faceEdges_[addedIntFaceIndices[indexI]][1] =
                (
                    addedEdgeIndices[indexI]
                );

                // Configure faceEdges for this boundary face
                tmpFaceEdges[0] = addedEdgeIndices[indexI];
                tmpFaceEdges[1] = newEdgeIndex;
                tmpFaceEdges[2] = ringEntities[2][indexI];

                // Modify faceEdges for the hull face
                replaceLabel
                (
                    ringEntities[2][indexI],
                    addedEdgeIndices[indexI],
                    faceEdges_[faceHull[indexI]]
                );

                // Modify edgeFaces for the edge connected to newEdge[1]
                replaceLabel
                (
                    faceHull[indexI],
                    addedFaceIndices[indexI],
                    edgeFaces_[ringEntities[2][indexI]]
                );

                // Modify edgePoints for the edge connected to newEdge[1]
                replaceLabel
                (
                    thisEdge[0],
                    newPointIndex,
                    edgePoints_[ringEntities[2][indexI]]
                );

                // Add the faceEdges entry
                faceEdges_.append(tmpFaceEdges);

                // Add an entry to newEdgeFaces
                newEdgeFaces[indexI] = addedFaceIndices[indexI];

                // Add an entry for this cell
                newCell[2] = addedFaceIndices[indexI];
            }
            else
            // Check if a cell was added before this
            if (addedCellIndices[prevI] != -1)
            {
                // Configure the interior face
                tmpTriFace[0] = vertexHull[indexI];
                tmpTriFace[1] = newEdge[1];
                tmpTriFace[2] = newPointIndex;

                // Insert the face
                addedFaceIndices[indexI] =
                (
                    insertFace
                    (
                        -1,
                        tmpTriFace,
                        addedCellIndices[prevI],
                        addedCellIndices[indexI]
                    )
                );

                // Add this face to the map.
                map.addFace(addedFaceIndices[indexI]);

                // Configure edgeFaces
                tmpIntEdgeFaces[0] = faceHull[indexI];
                tmpIntEdgeFaces[1] = addedIntFaceIndices[indexI];
                tmpIntEdgeFaces[2] = addedFaceIndices[indexI];
                tmpIntEdgeFaces[3] = addedIntFaceIndices[prevI];

                // Configure edgePoints
                tmpIntEdgePoints[0] = thisEdge[0];
                tmpIntEdgePoints[1] = vertexHull[nextI];
                tmpIntEdgePoints[2] = newEdge[1];
                tmpIntEdgePoints[3] = vertexHull[prevI];

                // Add an internal edge
                addedEdgeIndices[indexI] =
                (
                    insertEdge
                    (
                        -1,
                        edge(newPointIndex,vertexHull[indexI]),
                        tmpIntEdgeFaces,
                        tmpIntEdgePoints
                    )
                );

                // Add this edge to the map.
                map.addEdge(addedEdgeIndices[indexI]);

                // Add this edge to the interior-face faceEdges entry..
                faceEdges_[addedIntFaceIndices[indexI]][1] =
                (
                    addedEdgeIndices[indexI]
                );

                // ... and to the previous interior face as well
                faceEdges_[addedIntFaceIndices[prevI]][2] =
                (
                    addedEdgeIndices[indexI]
                );

                // Configure faceEdges for this split interior face
                tmpFaceEdges[0] = addedEdgeIndices[indexI];
                tmpFaceEdges[1] = newEdgeIndex;
                tmpFaceEdges[2] = ringEntities[2][indexI];

                // Modify faceEdges for the hull face
                replaceLabel
                (
                    ringEntities[2][indexI],
                    addedEdgeIndices[indexI],
                    faceEdges_[faceHull[indexI]]
                );

                // Modify edgeFaces for the edge connected to newEdge[1]
                replaceLabel
                (
                    faceHull[indexI],
                    addedFaceIndices[indexI],
                    edgeFaces_[ringEntities[2][indexI]]
                );

                // Modify edgePoints for the edge connected to newEdge[1]
                replaceLabel
                (
                    thisEdge[0],
                    newPointIndex,
                    edgePoints_[ringEntities[2][indexI]]
                );

                // Add the faceEdges entry
                faceEdges_.append(tmpFaceEdges);

                // Add an entry to newEdgeFaces
                newEdgeFaces[indexI] = addedFaceIndices[indexI];

                // Add an entry for this cell
                newCell[2] = addedFaceIndices[indexI];

                // Make the final entry for the previous cell
                cells_[addedCellIndices[prevI]][3] = addedFaceIndices[indexI];
            }

            // Do the first interior face at the end
            if (indexI == vertexHull.size() - 1)
            {
                // Configure the interior face
                tmpTriFace[0] = newPointIndex;
                tmpTriFace[1] = newEdge[1];
                tmpTriFace[2] = vertexHull[0];

                // Insert the face
                addedFaceIndices[0] =
                (
                    insertFace
                    (
                        -1,
                        tmpTriFace,
                        addedCellIndices[0],
                        addedCellIndices[indexI]
                    )
                );

                // Add this face to the map.
                map.addFace(addedFaceIndices[0]);

                // Configure edgeFaces
                tmpIntEdgeFaces[0] = faceHull[0];
                tmpIntEdgeFaces[1] = addedIntFaceIndices[0];
                tmpIntEdgeFaces[2] = addedFaceIndices[0];
                tmpIntEdgeFaces[3] = addedIntFaceIndices[indexI];

                // Configure edgePoints
                tmpIntEdgePoints[0] = thisEdge[0];
                tmpIntEdgePoints[1] = vertexHull[1];
                tmpIntEdgePoints[2] = newEdge[1];
                tmpIntEdgePoints[3] = vertexHull[indexI];

                // Add an internal edge
                addedEdgeIndices[0] =
                (
                    insertEdge
                    (
                        -1,
                        edge(newPointIndex,vertexHull[0]),
                        tmpIntEdgeFaces,
                        tmpIntEdgePoints
                    )
                );

                // Add this edge to the map.
                map.addEdge(addedEdgeIndices[0]);

                // Add this edge to the interior-face faceEdges entry..
                faceEdges_[addedIntFaceIndices[0]][1] =
                (
                    addedEdgeIndices[0]
                );

                // ... and to the previous interior face as well
                faceEdges_[addedIntFaceIndices[indexI]][2] =
                (
                    addedEdgeIndices[0]
                );

                // Configure faceEdges for the first split face
                tmpFaceEdges[0] = addedEdgeIndices[0];
                tmpFaceEdges[1] = newEdgeIndex;
                tmpFaceEdges[2] = ringEntities[2][0];

                // Modify faceEdges for the hull face
                replaceLabel
                (
                    ringEntities[2][0],
                    addedEdgeIndices[0],
                    faceEdges_[faceHull[0]]
                );

                // Modify edgeFaces for the edge connected to newEdge[1]
                replaceLabel
                (
                    faceHull[0],
                    addedFaceIndices[0],
                    edgeFaces_[ringEntities[2][0]]
                );

                // Modify edgePoints for the edge connected to newEdge[1]
                replaceLabel
                (
                    thisEdge[0],
                    newPointIndex,
                    edgePoints_[ringEntities[2][0]]
                );

                // Add the faceEdges entry
                faceEdges_.append(tmpFaceEdges);

                // Add an entry to newEdgeFaces
                newEdgeFaces[0] = addedFaceIndices[0];

                // Add an entry for this cell
                newCell[3] = addedFaceIndices[0];

                // Make the final entry for the first cell
                cells_[addedCellIndices[0]][2] = addedFaceIndices[0];
            }
        }
        else
        {
            // Configure the final boundary face
            tmpTriFace[0] = vertexHull[indexI];
            tmpTriFace[1] = newEdge[1];
            tmpTriFace[2] = newPointIndex;

            // Insert the face
            addedFaceIndices[indexI] =
            (
                insertFace
                (
                    whichPatch(faceHull[indexI]),
                    tmpTriFace,
                    addedCellIndices[prevI],
                    -1
                )
            );

            // Add this face to the map.
            map.addFace(addedFaceIndices[indexI]);

            // Generate mapping information for this new face
            label parent;

            if (faceHull[indexI] < nOldFaces_)
            {
                parent = faceHull[indexI];
            }
            else
            {
                parent = faceParents_[faceHull[indexI]];
            }

            // Insert the parent face
            faceParents_.insert(addedFaceIndices[indexI], parent);

            // Configure edgeFaces
            tmpEdgeFaces[0] = addedFaceIndices[indexI];
            tmpEdgeFaces[1] = addedIntFaceIndices[prevI];
            tmpEdgeFaces[2] = faceHull[indexI];

            // Configure edgePoints
            tmpEdgePoints[0] = newEdge[1];
            tmpEdgePoints[1] = vertexHull[prevI];
            tmpEdgePoints[2] = thisEdge[0];

            // Add an edge
            addedEdgeIndices[indexI] =
            (
                insertEdge
                (
                    whichPatch(faceHull[indexI]),
                    edge(newPointIndex,vertexHull[indexI]),
                    tmpEdgeFaces,
                    tmpEdgePoints
                )
            );

            // Add this edge to the map.
            map.addEdge(addedEdgeIndices[indexI]);

            // Add a faceEdges entry to the previous interior face
            faceEdges_[addedIntFaceIndices[prevI]][2] =
            (
                addedEdgeIndices[indexI]
            );

            // Configure faceEdges for the final boundary face
            tmpFaceEdges[0] = addedEdgeIndices[indexI];
            tmpFaceEdges[1] = newEdgeIndex;
            tmpFaceEdges[2] = ringEntities[2][indexI];

            // Modify faceEdges for the hull face
            replaceLabel
            (
                ringEntities[2][indexI],
                addedEdgeIndices[indexI],
                faceEdges_[faceHull[indexI]]
            );

            // Modify edgeFaces for the edge connected to newEdge[1]
            replaceLabel
            (
                faceHull[indexI],
                addedFaceIndices[indexI],
                edgeFaces_[ringEntities[2][indexI]]
            );

            // Modify edgePoints for the edge connected to newEdge[1]
            replaceLabel
            (
                thisEdge[0],
                newPointIndex,
                edgePoints_[ringEntities[2][indexI]]
            );

            // Add the faceEdges entry
            faceEdges_.append(tmpFaceEdges);

            // Add an entry to newEdgeFaces
            newEdgeFaces[indexI] = addedFaceIndices[indexI];

            // Make the final entry for the previous cell
            cells_[addedCellIndices[prevI]][3] = addedFaceIndices[indexI];
        }
    }

    if (coupledModification_)
    {
        // Create a master/slave entry for the new edges on the patch.
        if (locallyCoupledEdge(eIndex))
        {
            // Add the new edge to the coupled stack.
            // Although slave edges might be added here,
            // they are removed during stack initialization anyway.
            edgeStack(tIndex).push(newEdgeIndex);

            // Since we don't know the corresponding edges on the couple patch,
            // loop through recently added edges and perform a geometric match.
            FixedList<bool, 3> foundMatch(false);
            FixedList<label, 3> checkList(-1);
            FixedList<point, 3> cCentres(vector::zero);

            // Fill in the edges to be check for...
            label eCounter = 0;

            // The new edge...
            edge& newEdge = edges_[newEdgeIndex];

            checkList[eCounter] = newEdgeIndex;
            cCentres[eCounter++] =
            (
                0.5*(points_[newEdge[0]] + points_[newEdge[1]])
            );

            // ... and two new boundary edges.
            forAll(addedEdgeIndices, edgeI)
            {
                if (whichEdgePatch(addedEdgeIndices[edgeI]) != -1)
                {
                    edge& bEdge = edges_[addedEdgeIndices[edgeI]];

                    checkList[eCounter] = addedEdgeIndices[edgeI];

                    cCentres[eCounter++] =
                    (
                        0.5*(points_[bEdge[0]] + points_[bEdge[1]])
                    );
                }
            }

            const labelList aeList = coupleMap.addedEdgeList();

            forAll(aeList, edgeI)
            {
                // Get an edge reference.
                edge& check = edges_[aeList[edgeI]];

                // Get the centre of the edge.
                vector centre = 0.5*(points_[check[0]] + points_[check[1]]);

                // Compare with all three entries.
                forAll (checkList, indexI)
                {
                    if (mag(cCentres[indexI] - centre) < gTol_)
                    {
                        if (bisectingSlave)
                        {
                            slaveToMaster_[pIndex].insert
                            (
                                checkList[indexI], aeList[edgeI]
                            );

                            masterToSlave_[pIndex].insert
                            (
                                aeList[edgeI], checkList[indexI]
                            );
                        }
                        else
                        {
                            masterToSlave_[pIndex].insert
                            (
                                checkList[indexI], aeList[edgeI]
                            );

                            slaveToMaster_[pIndex].insert
                            (
                                aeList[edgeI], checkList[indexI]
                            );
                        }

                        foundMatch[indexI] = true;

                        break;
                    }
                }

                // Are we done checking?
                if (foundMatch[0] && foundMatch[1] && foundMatch[2])
                {
                    break;
                }
            }

            if (!(foundMatch[0] && foundMatch[1] && foundMatch[2]))
            {
                FatalErrorIn
                (
                    "dynamicTopoFvMesh::bisectEdge"
                ) << "Failed to build coupled maps." << abort(FatalError);
            }
        }
        else
        if (processorCoupledEdge(eIndex))
        {
            // Look for matching slave edges on the patchSubMesh.

        }
    }

    if (debug > 2)
    {
        Info << "Vertices: " << vertexHull << endl;
        Info << "Edges: " << edgeHull << endl;
        Info << "Faces: " << faceHull << endl;
        Info << "Cells: " << cellHull << endl;

        Info << "Modified cells: " << endl;
        forAll(cellHull, cellI)
        {
            if (cellHull[cellI] == -1)
            {
                continue;
            }

            Info << cellHull[cellI] << ":: "
                 << cells_[cellHull[cellI]]
                 << endl;
        }

        Info << "Added cells: " << endl;
        forAll(addedCellIndices, cellI)
        {
            if (addedCellIndices[cellI] == -1)
            {
                continue;
            }

            Info << addedCellIndices[cellI] << ":: "
                 << cells_[addedCellIndices[cellI]]
                 << endl;
        }

        Info << "Modified faces: " << endl;
        forAll(faceHull, faceI)
        {
            Info << faceHull[faceI] << ":: "
                 << faces_[faceHull[faceI]] << ": "
                 << owner_[faceHull[faceI]] << ": "
                 << neighbour_[faceHull[faceI]] << " "
                 << "faceEdges:: " << faceEdges_[faceHull[faceI]]
                 << endl;
        }

        Info << "Added faces: " << endl;
        forAll(addedFaceIndices, faceI)
        {
            Info << addedFaceIndices[faceI] << ":: "
                 << faces_[addedFaceIndices[faceI]] << ": "
                 << owner_[addedFaceIndices[faceI]] << ": "
                 << neighbour_[addedFaceIndices[faceI]] << " "
                 << "faceEdges:: " << faceEdges_[addedFaceIndices[faceI]]
                 << endl;
        }
        forAll(addedIntFaceIndices, faceI)
        {
            if (addedIntFaceIndices[faceI] == -1)
            {
                continue;
            }

            Info << addedIntFaceIndices[faceI] << ":: "
                 << faces_[addedIntFaceIndices[faceI]] << ": "
                 << owner_[addedIntFaceIndices[faceI]] << ": "
                 << neighbour_[addedIntFaceIndices[faceI]] << " "
                 << "faceEdges:: "
                 << faceEdges_[addedIntFaceIndices[faceI]]
                 << endl;
        }

        Info << "New edge:: " << newEdgeIndex
             << ": " << edges_[newEdgeIndex] << nl
             << " edgeFaces:: " << edgeFaces_[newEdgeIndex] << nl
             << " edgePoints:: " << edgePoints_[newEdgeIndex]
             << endl;

        Info << "Added edges: " << endl;
        forAll(addedEdgeIndices, edgeI)
        {
            Info << addedEdgeIndices[edgeI]
                 << ":: " << edges_[addedEdgeIndices[edgeI]] << nl
                 << " edgeFaces:: " << edgeFaces_[addedEdgeIndices[edgeI]] << nl
                 << " edgePoints:: " << edgePoints_[addedEdgeIndices[edgeI]]
                 << endl;
        }

        Info << "New Point:: " << newPointIndex << endl;
        Info << "pointEdges:: " << pointEdges_[newPointIndex] << endl;

        // Write out VTK files after change
        if (debug > 3)
        {
            labelList newHull(cellHull.size() + addedCellIndices.size(), 0);

            // Combine both lists into one.
            forAll(cellHull, i)
            {
                newHull[i] = cellHull[i];
            }

            label start = cellHull.size();

            for(label i = start; i < newHull.size(); i++)
            {
                newHull[i] = addedCellIndices[i - start];
            }

            writeVTK
            (
                Foam::name(eIndex)
              + "Bisect_1",
                newHull
            );
        }
    }

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    nBisections_++;

    // Increment the number of modifications
    nModifications_++;

    // Return the changeMap
    return map;
}

// Method for the trisection of a face in 3D
// - Returns a changeMap with a type specifying:
//    -1: Bisection failed since max number of topo-changes was reached.
// - AddedPoint is the index of the newly added point.
const dynamicTopoFvMesh::changeMap
dynamicTopoFvMesh::trisectFace
(
    const label fIndex
)
{
    // Face trisection performs the following operations:
    //      [1] Add a point at middle of the face
    //      [2] Remove the face and add three new faces in place.
    //      [3] Add three cells for each trisected cell (remove the originals).
    //      [4] Create one internal edge for each trisected cell.
    //      [5] Create three edges for the trisected face.
    //      [6] Create three internal faces for each trisected cell.
    //      Update faceEdges, edgeFaces and edgePoints information.

    // Figure out which thread this is...
    label tIndex = self(), pIndex = -1;

    // Prepare the changeMaps
    changeMap map, coupleMap;
    bool bisectingSlave = false;

    if
    (
        (nModifications_ > maxModifications_) &&
        (maxModifications_ > -1)
    )
    {
        // Reached the max allowable topo-changes.
        edgeStack(tIndex).clear();

        return map;
    }

    if (coupledModification_)
    {
        // Are all edges of this face locally coupled?
        bool locallyCoupledFace = false;

        labelList& fEdges = faceEdges_[fIndex];

        if
        (
            locallyCoupledEdge(fEdges[0]) &&
            locallyCoupledEdge(fEdges[1]) &&
            locallyCoupledEdge(fEdges[2])
        )
        {
            locallyCoupledFace = true;
        }

        if (locallyCoupledFace)
        {
            // Temporarily turn off coupledModification
            coupledModification_ = false;

            // Loop through master/slave maps
            // and determine the coupled edge index.
            forAll(masterToSlave_, indexI)
            {
                if
                (
                    masterToSlave_[indexI].found(fEdges[0]) &&
                    masterToSlave_[indexI].found(fEdges[1]) &&
                    masterToSlave_[indexI].found(fEdges[2])
                )
                {
                    bool foundFace = false;
                    label slaveFace = -1;

                    // Search for a boundary face shared
                    // by two of the edges.
                    labelList& eFirstFaces = edgeFaces_[fEdges[0]];
                    labelList& eSecondFaces = edgeFaces_[fEdges[1]];

                    forAll(eFirstFaces, edgeI)
                    {
                        if (neighbour_[eFirstFaces[edgeI]] == -1)
                        {
                            forAll(eSecondFaces, edgeJ)
                            {
                                if
                                (
                                    eFirstFaces[edgeI]
                                 == eSecondFaces[edgeJ]
                                )
                                {
                                    slaveFace = eFirstFaces[edgeI];
                                    foundFace = true;
                                    break;
                                }
                            }
                        }

                        if (foundFace)
                        {
                            break;
                        }
                    }

                    // Trisect the slave.
                    coupleMap = trisectFace(slaveFace);

                    // Keep this index for master/slave mapping.
                    pIndex = indexI;
                }

                // The following bit happens only during the sliver
                // exudation process.
                if
                (
                    slaveToMaster_[indexI].found(fEdges[0]) &&
                    slaveToMaster_[indexI].found(fEdges[1]) &&
                    slaveToMaster_[indexI].found(fEdges[2])
                )
                {
                    bool foundFace = false;
                    label masterFace = -1;

                    // Search for a boundary face shared
                    // by two of the edges.
                    labelList& eFirstFaces = edgeFaces_[fEdges[0]];
                    labelList& eSecondFaces = edgeFaces_[fEdges[1]];

                    forAll(eFirstFaces, edgeI)
                    {
                        if (neighbour_[eFirstFaces[edgeI]] == -1)
                        {
                            forAll(eSecondFaces, edgeJ)
                            {
                                if
                                (
                                    eFirstFaces[edgeI]
                                 == eSecondFaces[edgeJ]
                                )
                                {
                                    masterFace = eFirstFaces[edgeI];
                                    foundFace = true;
                                    break;
                                }
                            }
                        }

                        if (foundFace)
                        {
                            break;
                        }
                    }

                    // Trisect the master.
                    coupleMap = trisectFace(masterFace);

                    // Keep this index for master/slave mapping.
                    pIndex = indexI;

                    // Notice that we are bisecting a slave edge.
                    bisectingSlave = true;
                }
            }
        }
    }

    // References
    face& thisFace = faces_[fIndex];
    cell& firstCell = cells_[owner_[fIndex]];

    // Hull variables
    bool foundApex;
    face tmpTriFace(3);
    labelList newTriEdgeFaces(3), newTriEdgePoints(3);
    labelList newQuadEdgeFaces(4), newQuadEdgePoints(4);

    FixedList<label,2> apexPoint(-1);
    FixedList<face,3>  checkFace(face(3));
    FixedList<label,5> newEdgeIndex(-1);
    FixedList<label,9> newFaceIndex(-1);
    FixedList<label,6> newCellIndex(-1);
    FixedList<labelList, 9> newFaceEdges(labelList(3));

    // Counters for entities
    FixedList<label, 9> nE(0);
    FixedList<label, 6> nF(0);

    // Determine the two cells to be removed
    FixedList<label,2> cellsForRemoval;
    cellsForRemoval[0] = owner_[fIndex];
    cellsForRemoval[1] = neighbour_[fIndex];

    if (debug > 1)
    {
        Info << nl << nl << "Face: " << fIndex
             << ": " << thisFace << " is to be trisected. " << endl;

        // Write out VTK files prior to change
        if (debug > 3)
        {
            labelList vtkCells;

            if (neighbour_[fIndex] == -1)
            {
                vtkCells.setSize(1);
                vtkCells[0] = owner_[fIndex];
            }
            else
            {
                vtkCells.setSize(2);
                vtkCells[0] = owner_[fIndex];
                vtkCells[1] = neighbour_[fIndex];
            }

            writeVTK
            (
                Foam::name(fIndex)
              + "Trisect_0",
                vtkCells
            );
        }
    }

    // Add a new point to the end of the list
    label newPointIndex =
    (
        points_.append
        (
            (
                points_[thisFace[0]]
              + points_[thisFace[1]]
              + points_[thisFace[2]]
            )
            /3.0
        )
    );

    // Add an entry to pointEdges as well
    pointEdges_.append(labelList(0));

    nPoints_++;

    // Add this point to the map.
    map.addPoint(newPointIndex);

    // Add three new cells to the end of the cell list
    newCellIndex[0] = cells_.append(cell(4));
    newCellIndex[1] = cells_.append(cell(4));
    newCellIndex[2] = cells_.append(cell(4));

    cell &newTetCell0 = cells_[newCellIndex[0]];
    cell &newTetCell1 = cells_[newCellIndex[1]];
    cell &newTetCell2 = cells_[newCellIndex[2]];

    // Update length-scale info
    if (edgeModification_)
    {
        scalar parentScale = lengthScale_[cellsForRemoval[0]];

        for (label i = 0; i < 3; i++)
        {
            lengthScale_.append(parentScale);
        }
    }

    // Find the apex point for this cell
    foundApex = false;

    forAll(firstCell, faceI)
    {
        face& faceToCheck = faces_[firstCell[faceI]];

        forAll(faceToCheck, pointI)
        {
            if
            (
                faceToCheck[pointI] != thisFace[0] &&
                faceToCheck[pointI] != thisFace[1] &&
                faceToCheck[pointI] != thisFace[2]
            )
            {
                apexPoint[0] = faceToCheck[pointI];

                foundApex = true;
                break;
            }
        }

        if (foundApex)
        {
            break;
        }
    }

    // Insert three new internal faces

    // First face: Owner: newCellIndex[0], Neighbour: newCellIndex[1]
    tmpTriFace[0] = newPointIndex;
    tmpTriFace[1] = thisFace[0];
    tmpTriFace[2] = apexPoint[0];

    newFaceIndex[0] =
    (
        insertFace
        (
            -1,
            tmpTriFace,
            newCellIndex[0],
            newCellIndex[1]
        )
    );

    // Second face: Owner: newCellIndex[1], Neighbour: newCellIndex[2]
    tmpTriFace[0] = newPointIndex;
    tmpTriFace[1] = thisFace[1];
    tmpTriFace[2] = apexPoint[0];

    newFaceIndex[1] =
    (
        insertFace
        (
            -1,
            tmpTriFace,
            newCellIndex[1],
            newCellIndex[2]
        )
    );

    // Third face: Owner: newCellIndex[0], Neighbour: newCellIndex[2]
    tmpTriFace[0] = newPointIndex;
    tmpTriFace[1] = apexPoint[0];
    tmpTriFace[2] = thisFace[2];

    newFaceIndex[2] =
    (
        insertFace
        (
            -1,
            tmpTriFace,
            newCellIndex[0],
            newCellIndex[2]
        )
    );

    // Add an entry to edgeFaces
    newTriEdgeFaces[0] = newFaceIndex[0];
    newTriEdgeFaces[1] = newFaceIndex[1];
    newTriEdgeFaces[2] = newFaceIndex[2];

    // Add an entry for edgePoints as well
    newTriEdgePoints[0] = thisFace[0];
    newTriEdgePoints[1] = thisFace[1];
    newTriEdgePoints[2] = thisFace[2];

    // Add a new internal edge to the mesh
    newEdgeIndex[0] =
    (
        insertEdge
        (
            -1,
            edge
            (
               newPointIndex,
               apexPoint[0]
            ),
            newTriEdgeFaces,
            newTriEdgePoints
        )
    );

    // Configure faceEdges with the new internal edge
    newFaceEdges[0][nE[0]++] = newEdgeIndex[0];
    newFaceEdges[1][nE[1]++] = newEdgeIndex[0];
    newFaceEdges[2][nE[2]++] = newEdgeIndex[0];

    // Add the newly created faces to cells
    newTetCell0[nF[0]++] = newFaceIndex[0];
    newTetCell0[nF[0]++] = newFaceIndex[2];
    newTetCell1[nF[1]++] = newFaceIndex[0];
    newTetCell1[nF[1]++] = newFaceIndex[1];
    newTetCell2[nF[2]++] = newFaceIndex[1];
    newTetCell2[nF[2]++] = newFaceIndex[2];

    // Define the three faces to check for orientation:
    checkFace[0][0] = thisFace[2];
    checkFace[0][1] = apexPoint[0];
    checkFace[0][2] = thisFace[0];

    checkFace[1][0] = thisFace[0];
    checkFace[1][1] = apexPoint[0];
    checkFace[1][2] = thisFace[1];

    checkFace[2][0] = thisFace[1];
    checkFace[2][1] = apexPoint[0];
    checkFace[2][2] = thisFace[2];

    // Check the orientation of faces on the first cell.
    forAll(firstCell, faceI)
    {
        label faceIndex = firstCell[faceI];

        if (faceIndex == fIndex)
        {
            continue;
        }

        face& faceToCheck = faces_[faceIndex];
        label cellIndex = cellsForRemoval[0];
        label newIndex = -1;

        // Check against faces.
        if (compare(faceToCheck, checkFace[0]) != 0)
        {
            newIndex = newCellIndex[0];
            newTetCell0[nF[0]++] = faceIndex;
        }
        else
        if (compare(faceToCheck, checkFace[1]) != 0)
        {
            newIndex = newCellIndex[1];
            newTetCell1[nF[1]++] = faceIndex;
        }
        else
        if (compare(faceToCheck, checkFace[2]) != 0)
        {
            newIndex = newCellIndex[2];
            newTetCell2[nF[2]++] = faceIndex;
        }
        else
        {
            // Something's terribly wrong.
            FatalErrorIn("dynamicTopoFvMesh::trisectFace()")
                << "Failed to determine a face match."
                << abort(FatalError);
        }

        // Check if a face-flip is necessary
        if (owner_[faceIndex] == cellIndex)
        {
            if (neighbour_[faceIndex] == -1)
            {
                // Change the owner
                owner_[faceIndex] = newIndex;
            }
            else
            {
                // Flip this face
                faces_[faceIndex] = faceToCheck.reverseFace();
                owner_[faceIndex] = neighbour_[faceIndex];
                neighbour_[faceIndex] = newIndex;
            }
        }
        else
        {
            // Flip is unnecessary. Just update neighbour
            neighbour_[faceIndex] = newIndex;
        }
    }

    if (cellsForRemoval[1] == -1)
    {
        // Boundary face. Determine its patch.
        label facePatch = whichPatch(fIndex);

        // Add three new boundary faces.

        // Fourth face: Owner: newCellIndex[0], Neighbour: -1
        tmpTriFace[0] = newPointIndex;
        tmpTriFace[1] = thisFace[2];
        tmpTriFace[2] = thisFace[0];

        newFaceIndex[3] =
        (
            insertFace
            (
                facePatch,
                tmpTriFace,
                newCellIndex[0],
                -1
            )
        );

        // Fifth face: Owner: newCellIndex[1], Neighbour: -1
        tmpTriFace[0] = newPointIndex;
        tmpTriFace[1] = thisFace[0];
        tmpTriFace[2] = thisFace[1];

        newFaceIndex[4] =
        (
            insertFace
            (
                facePatch,
                tmpTriFace,
                newCellIndex[1],
                -1
            )
        );

        // Sixth face: Owner: newCellIndex[2], Neighbour: -1
        tmpTriFace[0] = newPointIndex;
        tmpTriFace[1] = thisFace[1];
        tmpTriFace[2] = thisFace[2];

        newFaceIndex[5] =
        (
            insertFace
            (
                facePatch,
                tmpTriFace,
                newCellIndex[2],
                -1
            )
        );

        // Add the newly created faces to cells
        newTetCell0[nF[0]++] = newFaceIndex[3];
        newTetCell1[nF[1]++] = newFaceIndex[4];
        newTetCell2[nF[2]++] = newFaceIndex[5];

        // Configure edgeFaces and edgePoints for three new boundary edges.
        newTriEdgeFaces[0] = newFaceIndex[4];
        newTriEdgeFaces[1] = newFaceIndex[0];
        newTriEdgeFaces[2] = newFaceIndex[3];

        newTriEdgePoints[0] = thisFace[1];
        newTriEdgePoints[1] = apexPoint[0];
        newTriEdgePoints[2] = thisFace[2];

        newEdgeIndex[1] =
        (
            insertEdge
            (
                facePatch,
                edge
                (
                   newPointIndex,
                   thisFace[0]
                ),
                newTriEdgeFaces,
                newTriEdgePoints
            )
        );

        newTriEdgeFaces[0] = newFaceIndex[5];
        newTriEdgeFaces[1] = newFaceIndex[1];
        newTriEdgeFaces[2] = newFaceIndex[4];

        newTriEdgePoints[0] = thisFace[2];
        newTriEdgePoints[1] = apexPoint[0];
        newTriEdgePoints[2] = thisFace[0];

        newEdgeIndex[2] =
        (
            insertEdge
            (
                facePatch,
                edge
                (
                   newPointIndex,
                   thisFace[1]
                ),
                newTriEdgeFaces,
                newTriEdgePoints
            )
        );

        newTriEdgeFaces[0] = newFaceIndex[3];
        newTriEdgeFaces[1] = newFaceIndex[2];
        newTriEdgeFaces[2] = newFaceIndex[5];

        newTriEdgePoints[0] = thisFace[0];
        newTriEdgePoints[1] = apexPoint[0];
        newTriEdgePoints[2] = thisFace[1];

        newEdgeIndex[3] =
        (
            insertEdge
            (
                facePatch,
                edge
                (
                   newPointIndex,
                   thisFace[2]
                ),
                newTriEdgeFaces,
                newTriEdgePoints
            )
        );

        // Configure faceEdges with the three new edges.
        newFaceEdges[0][nE[0]++] = newEdgeIndex[1];
        newFaceEdges[1][nE[1]++] = newEdgeIndex[2];
        newFaceEdges[2][nE[2]++] = newEdgeIndex[3];

        newFaceEdges[3][nE[3]++] = newEdgeIndex[1];
        newFaceEdges[3][nE[3]++] = newEdgeIndex[3];
        newFaceEdges[4][nE[4]++] = newEdgeIndex[1];
        newFaceEdges[4][nE[4]++] = newEdgeIndex[2];
        newFaceEdges[5][nE[5]++] = newEdgeIndex[2];
        newFaceEdges[5][nE[5]++] = newEdgeIndex[3];

        // Define the six edges to check while building faceEdges:
        FixedList<edge,6> check;

        check[0][0] = apexPoint[0]; check[0][1] = thisFace[0];
        check[1][0] = apexPoint[0]; check[1][1] = thisFace[1];
        check[2][0] = apexPoint[0]; check[2][1] = thisFace[2];

        check[3][0] = thisFace[2]; check[3][1] = thisFace[0];
        check[4][0] = thisFace[0]; check[4][1] = thisFace[1];
        check[5][0] = thisFace[1]; check[5][1] = thisFace[2];

        // Build a list of cellEdges
        labelHashSet cellEdges;

        forAll(firstCell, faceI)
        {
            labelList fEdges = faceEdges_[firstCell[faceI]];

            forAll(fEdges, edgeI)
            {
                if (!cellEdges.found(fEdges[edgeI]))
                {
                    cellEdges.insert(fEdges[edgeI]);
                }
            }
        }

        // Loop through cellEdges, and perform appropriate actions.
        forAllIter(labelHashSet::iterator, cellEdges, eIter)
        {
            edge& edgeToCheck = edges_[eIter.key()];

            // Check against the specified edges.
            if (edgeToCheck == check[0])
            {
                insertLabel
                (
                    newPointIndex,
                    thisFace[1],
                    thisFace[2],
                    edgePoints_[eIter.key()]
                );

                sizeUpList(newFaceIndex[0], edgeFaces_[eIter.key()]);
                newFaceEdges[0][nE[0]++] = eIter.key();
            }

            if (edgeToCheck == check[1])
            {
                insertLabel
                (
                    newPointIndex,
                    thisFace[0],
                    thisFace[2],
                    edgePoints_[eIter.key()]
                );

                sizeUpList(newFaceIndex[1], edgeFaces_[eIter.key()]);
                newFaceEdges[1][nE[1]++] = eIter.key();
            }

            if (edgeToCheck == check[2])
            {
                insertLabel
                (
                    newPointIndex,
                    thisFace[0],
                    thisFace[1],
                    edgePoints_[eIter.key()]
                );

                sizeUpList(newFaceIndex[2], edgeFaces_[eIter.key()]);
                newFaceEdges[2][nE[2]++] = eIter.key();
            }

            if (edgeToCheck == check[3])
            {
                replaceLabel
                (
                    thisFace[1],
                    newPointIndex,
                    edgePoints_[eIter.key()]
                );

                replaceLabel
                (
                    fIndex,
                    newFaceIndex[3],
                    edgeFaces_[eIter.key()]
                );

                newFaceEdges[3][nE[3]++] = eIter.key();
            }

            if (edgeToCheck == check[4])
            {
                replaceLabel
                (
                    thisFace[2],
                    newPointIndex,
                    edgePoints_[eIter.key()]
                );

                replaceLabel
                (
                    fIndex,
                    newFaceIndex[4],
                    edgeFaces_[eIter.key()]
                );

                newFaceEdges[4][nE[4]++] = eIter.key();
            }

            if (edgeToCheck == check[5])
            {
                replaceLabel
                (
                    thisFace[0],
                    newPointIndex,
                    edgePoints_[eIter.key()]
                );

                replaceLabel
                (
                    fIndex,
                    newFaceIndex[5],
                    edgeFaces_[eIter.key()]
                );

                newFaceEdges[5][nE[5]++] = eIter.key();
            }
        }

        // Now that faceEdges has been configured, append them to the list.
        for (label i = 0; i < 6; i++)
        {
            faceEdges_.append(newFaceEdges[i]);
        }
    }
    else
    {
        cell& secondCell = cells_[neighbour_[fIndex]];

        // Add three new cells to the end of the cell list
        newCellIndex[3] = cells_.append(cell(4));
        newCellIndex[4] = cells_.append(cell(4));
        newCellIndex[5] = cells_.append(cell(4));

        cell &newTetCell3 = cells_[newCellIndex[3]];
        cell &newTetCell4 = cells_[newCellIndex[4]];
        cell &newTetCell5 = cells_[newCellIndex[5]];

        // Update length-scale info
        if (edgeModification_)
        {
            scalar parentScale = lengthScale_[cellsForRemoval[1]];

            for (label i = 0; i < 3; i++)
            {
                lengthScale_.append(parentScale);
            }
        }

        // Find the apex point for this cell
        foundApex = false;

        forAll(secondCell, faceI)
        {
            face& faceToCheck = faces_[secondCell[faceI]];

            forAll(faceToCheck, pointI)
            {
                if
                (
                    faceToCheck[pointI] != thisFace[0] &&
                    faceToCheck[pointI] != thisFace[1] &&
                    faceToCheck[pointI] != thisFace[2]
                )
                {
                    apexPoint[1] = faceToCheck[pointI];

                    foundApex = true;
                    break;
                }
            }

            if (foundApex)
            {
                break;
            }
        }

        // Add six new interior faces.

        // Fourth face: Owner: newCellIndex[0], Neighbour: newCellIndex[3]
        tmpTriFace[0] = newPointIndex;
        tmpTriFace[1] = thisFace[2];
        tmpTriFace[2] = thisFace[0];

        newFaceIndex[3] =
        (
            insertFace
            (
                -1,
                tmpTriFace,
                newCellIndex[0],
                newCellIndex[3]
            )
        );

        // Fifth face: Owner: newCellIndex[1], Neighbour: newCellIndex[4]
        tmpTriFace[0] = newPointIndex;
        tmpTriFace[1] = thisFace[0];
        tmpTriFace[2] = thisFace[1];

        newFaceIndex[4] =
        (
            insertFace
            (
                -1,
                tmpTriFace,
                newCellIndex[1],
                newCellIndex[4]
            )
        );

        // Sixth face: Owner: newCellIndex[2], Neighbour: newCellIndex[5]
        tmpTriFace[0] = newPointIndex;
        tmpTriFace[1] = thisFace[1];
        tmpTriFace[2] = thisFace[2];

        newFaceIndex[5] =
        (
            insertFace
            (
                -1,
                tmpTriFace,
                newCellIndex[2],
                newCellIndex[5]
            )
        );

        // Seventh face: Owner: newCellIndex[3], Neighbour: newCellIndex[4]
        tmpTriFace[0] = newPointIndex;
        tmpTriFace[1] = apexPoint[1];
        tmpTriFace[2] = thisFace[0];

        newFaceIndex[6] =
        (
            insertFace
            (
                -1,
                tmpTriFace,
                newCellIndex[3],
                newCellIndex[4]
            )
        );

        // Eighth face: Owner: newCellIndex[4], Neighbour: newCellIndex[5]
        tmpTriFace[0] = newPointIndex;
        tmpTriFace[1] = apexPoint[1];
        tmpTriFace[2] = thisFace[1];

        newFaceIndex[7] =
        (
            insertFace
            (
                -1,
                tmpTriFace,
                newCellIndex[4],
                newCellIndex[5]
            )
        );

        // Ninth face: Owner: newCellIndex[3], Neighbour: newCellIndex[5]
        tmpTriFace[0] = newPointIndex;
        tmpTriFace[1] = thisFace[2];
        tmpTriFace[2] = apexPoint[1];

        newFaceIndex[8] =
        (
            insertFace
            (
                -1,
                tmpTriFace,
                newCellIndex[3],
                newCellIndex[5]
            )
        );

        // Add the newly created faces to cells
        newTetCell3[nF[3]++] = newFaceIndex[6];
        newTetCell3[nF[3]++] = newFaceIndex[8];
        newTetCell4[nF[4]++] = newFaceIndex[6];
        newTetCell4[nF[4]++] = newFaceIndex[7];
        newTetCell5[nF[5]++] = newFaceIndex[7];
        newTetCell5[nF[5]++] = newFaceIndex[8];

        newTetCell0[nF[0]++] = newFaceIndex[3];
        newTetCell1[nF[1]++] = newFaceIndex[4];
        newTetCell2[nF[2]++] = newFaceIndex[5];

        newTetCell3[nF[3]++] = newFaceIndex[3];
        newTetCell4[nF[4]++] = newFaceIndex[4];
        newTetCell5[nF[5]++] = newFaceIndex[5];

        // Define the three faces to check for orientation:
        checkFace[0][0] = thisFace[2];
        checkFace[0][1] = apexPoint[1];
        checkFace[0][2] = thisFace[0];

        checkFace[1][0] = thisFace[0];
        checkFace[1][1] = apexPoint[1];
        checkFace[1][2] = thisFace[1];

        checkFace[2][0] = thisFace[1];
        checkFace[2][1] = apexPoint[1];
        checkFace[2][2] = thisFace[2];

        // Check the orientation of faces on the second cell.
        forAll(secondCell, faceI)
        {
            label faceIndex = secondCell[faceI];

            if (faceIndex == fIndex)
            {
                continue;
            }

            face& faceToCheck = faces_[faceIndex];
            label cellIndex = cellsForRemoval[1];
            label newIndex = -1;

            // Check against faces.
            if (compare(faceToCheck, checkFace[0]) != 0)
            {
                newIndex = newCellIndex[3];
                newTetCell3[nF[3]++] = faceIndex;
            }
            else
            if (compare(faceToCheck, checkFace[1]) != 0)
            {
                newIndex = newCellIndex[4];
                newTetCell4[nF[4]++] = faceIndex;
            }
            else
            if (compare(faceToCheck, checkFace[2]) != 0)
            {
                newIndex = newCellIndex[5];
                newTetCell5[nF[5]++] = faceIndex;
            }
            else
            {
                // Something's terribly wrong.
                FatalErrorIn("dynamicTopoFvMesh::trisectFace()")
                    << "Failed to determine a face match."
                    << abort(FatalError);
            }

            // Check if a face-flip is necessary
            if (owner_[faceIndex] == cellIndex)
            {
                if (neighbour_[faceIndex] == -1)
                {
                    // Change the owner
                    owner_[faceIndex] = newIndex;
                }
                else
                {
                    // Flip this face
                    faces_[faceIndex] = faceToCheck.reverseFace();
                    owner_[faceIndex] = neighbour_[faceIndex];
                    neighbour_[faceIndex] = newIndex;
                }
            }
            else
            {
                // Flip is unnecessary. Just update neighbour
                neighbour_[faceIndex] = newIndex;
            }
        }

        // Configure edgeFaces and edgePoints for four new interior edges.
        newQuadEdgeFaces[0] = newFaceIndex[4];
        newQuadEdgeFaces[1] = newFaceIndex[0];
        newQuadEdgeFaces[2] = newFaceIndex[3];
        newQuadEdgeFaces[3] = newFaceIndex[6];

        newQuadEdgePoints[0] = thisFace[1];
        newQuadEdgePoints[1] = apexPoint[0];
        newQuadEdgePoints[2] = thisFace[2];
        newQuadEdgePoints[3] = apexPoint[1];

        newEdgeIndex[1] =
        (
            insertEdge
            (
                -1,
                edge
                (
                   newPointIndex,
                   thisFace[0]
                ),
                newQuadEdgeFaces,
                newQuadEdgePoints
            )
        );

        newQuadEdgeFaces[0] = newFaceIndex[5];
        newQuadEdgeFaces[1] = newFaceIndex[1];
        newQuadEdgeFaces[2] = newFaceIndex[4];
        newQuadEdgeFaces[3] = newFaceIndex[7];

        newQuadEdgePoints[0] = thisFace[2];
        newQuadEdgePoints[1] = apexPoint[0];
        newQuadEdgePoints[2] = thisFace[0];
        newQuadEdgePoints[3] = apexPoint[1];

        newEdgeIndex[2] =
        (
            insertEdge
            (
                -1,
                edge
                (
                   newPointIndex,
                   thisFace[1]
                ),
                newQuadEdgeFaces,
                newQuadEdgePoints
            )
        );

        newQuadEdgeFaces[0] = newFaceIndex[3];
        newQuadEdgeFaces[1] = newFaceIndex[2];
        newQuadEdgeFaces[2] = newFaceIndex[5];
        newQuadEdgeFaces[3] = newFaceIndex[8];

        newQuadEdgePoints[0] = thisFace[0];
        newQuadEdgePoints[1] = apexPoint[0];
        newQuadEdgePoints[2] = thisFace[1];
        newQuadEdgePoints[3] = apexPoint[1];

        newEdgeIndex[3] =
        (
            insertEdge
            (
                -1,
                edge
                (
                   newPointIndex,
                   thisFace[2]
                ),
                newQuadEdgeFaces,
                newQuadEdgePoints
            )
        );

        newTriEdgeFaces[0] = newFaceIndex[6];
        newTriEdgeFaces[1] = newFaceIndex[7];
        newTriEdgeFaces[2] = newFaceIndex[8];

        newTriEdgePoints[0] = thisFace[0];
        newTriEdgePoints[1] = thisFace[1];
        newTriEdgePoints[2] = thisFace[2];

        newEdgeIndex[4] =
        (
            insertEdge
            (
                -1,
                edge
                (
                   apexPoint[1],
                   newPointIndex
                ),
                newTriEdgeFaces,
                newTriEdgePoints
            )
        );

        // Configure faceEdges with the new internal edges
        newFaceEdges[0][nE[0]++] = newEdgeIndex[1];
        newFaceEdges[1][nE[1]++] = newEdgeIndex[2];
        newFaceEdges[2][nE[2]++] = newEdgeIndex[3];

        newFaceEdges[3][nE[3]++] = newEdgeIndex[1];
        newFaceEdges[3][nE[3]++] = newEdgeIndex[3];
        newFaceEdges[4][nE[4]++] = newEdgeIndex[1];
        newFaceEdges[4][nE[4]++] = newEdgeIndex[2];
        newFaceEdges[5][nE[5]++] = newEdgeIndex[2];
        newFaceEdges[5][nE[5]++] = newEdgeIndex[3];

        newFaceEdges[6][nE[6]++] = newEdgeIndex[1];
        newFaceEdges[7][nE[7]++] = newEdgeIndex[2];
        newFaceEdges[8][nE[8]++] = newEdgeIndex[3];

        newFaceEdges[6][nE[6]++] = newEdgeIndex[4];
        newFaceEdges[7][nE[7]++] = newEdgeIndex[4];
        newFaceEdges[8][nE[8]++] = newEdgeIndex[4];

        // Define the nine edges to check while building faceEdges:
        FixedList<edge,9> check;

        check[0][0] = apexPoint[0]; check[0][1] = thisFace[0];
        check[1][0] = apexPoint[0]; check[1][1] = thisFace[1];
        check[2][0] = apexPoint[0]; check[2][1] = thisFace[2];

        check[3][0] = thisFace[2]; check[3][1] = thisFace[0];
        check[4][0] = thisFace[0]; check[4][1] = thisFace[1];
        check[5][0] = thisFace[1]; check[5][1] = thisFace[2];

        check[6][0] = apexPoint[1]; check[6][1] = thisFace[0];
        check[7][0] = apexPoint[1]; check[7][1] = thisFace[1];
        check[8][0] = apexPoint[1]; check[8][1] = thisFace[2];

        // Build a list of cellEdges
        labelHashSet cellEdges;

        forAll(firstCell, faceI)
        {
            labelList fEdges = faceEdges_[firstCell[faceI]];

            forAll(fEdges, edgeI)
            {
                if (!cellEdges.found(fEdges[edgeI]))
                {
                    cellEdges.insert(fEdges[edgeI]);
                }
            }
        }

        forAll(secondCell, faceI)
        {
            labelList fEdges = faceEdges_[secondCell[faceI]];

            forAll(fEdges, edgeI)
            {
                if (!cellEdges.found(fEdges[edgeI]))
                {
                    cellEdges.insert(fEdges[edgeI]);
                }
            }
        }

        // Loop through cellEdges, and perform appropriate actions.
        forAllIter(labelHashSet::iterator, cellEdges, eIter)
        {
            edge& edgeToCheck = edges_[eIter.key()];

            // Check against the specified edges.
            if (edgeToCheck == check[0])
            {
                insertLabel
                (
                    newPointIndex,
                    thisFace[1],
                    thisFace[2],
                    edgePoints_[eIter.key()]
                );

                sizeUpList(newFaceIndex[0], edgeFaces_[eIter.key()]);
                newFaceEdges[0][nE[0]++] = eIter.key();
            }

            if (edgeToCheck == check[1])
            {
                insertLabel
                (
                    newPointIndex,
                    thisFace[0],
                    thisFace[2],
                    edgePoints_[eIter.key()]
                );

                sizeUpList(newFaceIndex[1], edgeFaces_[eIter.key()]);
                newFaceEdges[1][nE[1]++] = eIter.key();
            }

            if (edgeToCheck == check[2])
            {
                insertLabel
                (
                    newPointIndex,
                    thisFace[0],
                    thisFace[1],
                    edgePoints_[eIter.key()]
                );

                sizeUpList(newFaceIndex[2], edgeFaces_[eIter.key()]);
                newFaceEdges[2][nE[2]++] = eIter.key();
            }

            if (edgeToCheck == check[3])
            {
                replaceLabel
                (
                    thisFace[1],
                    newPointIndex,
                    edgePoints_[eIter.key()]
                );

                replaceLabel
                (
                    fIndex,
                    newFaceIndex[3],
                    edgeFaces_[eIter.key()]
                );

                newFaceEdges[3][nE[3]++] = eIter.key();
            }

            if (edgeToCheck == check[4])
            {
                replaceLabel
                (
                    thisFace[2],
                    newPointIndex,
                    edgePoints_[eIter.key()]
                );

                replaceLabel
                (
                    fIndex,
                    newFaceIndex[4],
                    edgeFaces_[eIter.key()]
                );

                newFaceEdges[4][nE[4]++] = eIter.key();
            }

            if (edgeToCheck == check[5])
            {
                replaceLabel
                (
                    thisFace[0],
                    newPointIndex,
                    edgePoints_[eIter.key()]
                );

                replaceLabel
                (
                    fIndex,
                    newFaceIndex[5],
                    edgeFaces_[eIter.key()]
                );

                newFaceEdges[5][nE[5]++] = eIter.key();
            }

            if (edgeToCheck == check[6])
            {
                insertLabel
                (
                    newPointIndex,
                    thisFace[1],
                    thisFace[2],
                    edgePoints_[eIter.key()]
                );

                sizeUpList(newFaceIndex[6], edgeFaces_[eIter.key()]);
                newFaceEdges[6][nE[6]++] = eIter.key();
            }

            if (edgeToCheck == check[7])
            {
                insertLabel
                (
                    newPointIndex,
                    thisFace[0],
                    thisFace[2],
                    edgePoints_[eIter.key()]
                );

                sizeUpList(newFaceIndex[7], edgeFaces_[eIter.key()]);
                newFaceEdges[7][nE[7]++] = eIter.key();
            }

            if (edgeToCheck == check[8])
            {
                insertLabel
                (
                    newPointIndex,
                    thisFace[0],
                    thisFace[1],
                    edgePoints_[eIter.key()]
                );

                sizeUpList(newFaceIndex[8], edgeFaces_[eIter.key()]);
                newFaceEdges[8][nE[8]++] = eIter.key();
            }
        }

        // Now that faceEdges has been configured, append them to the list.
        for (label i = 0; i < 9; i++)
        {
            faceEdges_.append(newFaceEdges[i]);
        }
    }

    // Now generate mapping info and remove entities.
    forAll(cellsForRemoval, cellI)
    {
        label cIndex = cellsForRemoval[cellI];

        if (cIndex == -1)
        {
            continue;
        }

        // Determine an appropriate parent cell
        label parent = -1;

        if (cIndex < nOldCells_)
        {
            parent = cIndex;
        }
        else
        {
            parent = cellParents_[cIndex];
        }

        if (cellI == 0)
        {
            for (label i = 0; i < 3; i++)
            {
                // Insert the parent cell [from first by default]
                cellParents_.insert(newCellIndex[i], parent);

                // Insert mapping info into the HashTable
                cellsFromCells_.insert
                (
                    newCellIndex[i],
                    objectMap
                    (
                        newCellIndex[i],
                        labelList(1, parent)
                    )
                );
            }
        }
        else
        {
            for (label i = 3; i < 6; i++)
            {
                // Insert the parent cell [from first by default]
                cellParents_.insert(newCellIndex[i], parent);

                // Insert mapping info into the HashTable
                cellsFromCells_.insert
                (
                    newCellIndex[i],
                    objectMap
                    (
                        newCellIndex[i],
                        labelList(1, parent)
                    )
                );
            }
        }

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

        // Update the number of cells.
        nCells_--;
    }

    // Now finally remove the face...
    removeFace(fIndex);

    if (debug > 2)
    {
        Info << "New Point:: " << newPointIndex << endl;

        labelList& pEdges = pointEdges_[newPointIndex];

        Info << "pointEdges:: " << pEdges << endl;

        Info << "Added edges: " << endl;
        forAll(pEdges, edgeI)
        {
            Info << pEdges[edgeI]
                 << ":: " << edges_[pEdges[edgeI]] << nl
                 << " edgeFaces:: " << edgeFaces_[pEdges[edgeI]] << nl
                 << " edgePoints:: " << edgePoints_[pEdges[edgeI]]
                 << endl;
        }

        Info << "Added faces: " << endl;
        forAll(newFaceIndex, faceI)
        {
            if (newFaceIndex[faceI] == -1)
            {
                continue;
            }

            Info << newFaceIndex[faceI] << ":: "
                 << cells_[newFaceIndex[faceI]]
                 << endl;
        }

        Info << "Added cells: " << endl;
        forAll(newCellIndex, cellI)
        {
            if (newCellIndex[cellI] == -1)
            {
                continue;
            }

            Info << newCellIndex[cellI] << ":: "
                 << cells_[newCellIndex[cellI]]
                 << endl;
        }

        // Write out VTK files after change
        if (debug > 3)
        {
            labelList vtkCells;

            if (cellsForRemoval[1] == -1)
            {
                vtkCells.setSize(3);

                // Fill in cell indices
                vtkCells[0] = newCellIndex[0];
                vtkCells[1] = newCellIndex[1];
                vtkCells[2] = newCellIndex[2];
            }
            else
            {
                vtkCells.setSize(6);

                // Fill in cell indices
                forAll(newCellIndex, indexI)
                {
                    vtkCells[indexI] = newCellIndex[indexI];
                }
            }

            writeVTK
            (
                Foam::name(fIndex)
              + "Trisect_1",
                vtkCells
            );
        }
    }

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    nBisections_++;

    // Increment the number of modifications
    nModifications_++;

    // Return the changeMap
    return map;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
