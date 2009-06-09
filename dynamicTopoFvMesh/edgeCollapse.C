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

// Method for the collapse of a quad-face in 2D
void dynamicTopoFvMesh::collapseQuadFace
(
    const label fIndex
)
{
    // Obtain a reference for this face...
    face& thisFace = faces_[fIndex];

    // This face is to be collapsed...
    if (debug > 1)
    {
        Info << nl << nl
             << "Face: " << fIndex << ": " << thisFace
             << " is to be collapsed. " << endl;
    }

    // Local variables
    FixedList<label,2> c0BdyIndex, c0IntIndex, c1BdyIndex, c1IntIndex;
    FixedList<face,2> c0BdyFace, c0IntFace, c1BdyFace, c1IntFace;
    FixedList<edge,4> checkEdge(edge(-1,-1));
    FixedList<label,4> checkEdgeIndex;
    face tmpTriFace(3);

    // Define checkEdges
    checkEdgeIndex[0] = getTriBoundaryEdge(fIndex);
    checkEdge[0] = edges_[checkEdgeIndex[0]];

    labelList& fEdges = faceEdges_[fIndex];

    forAll(fEdges, edgeI)
    {
        if (checkEdgeIndex[0] != fEdges[edgeI])
        {
            edge& thisEdge = edges_[fEdges[edgeI]];

            if
            (
                (
                    checkEdge[0].start() == thisEdge[0]
                 || checkEdge[0].start() == thisEdge[1]
                )
            )
            {
                checkEdgeIndex[1] = fEdges[edgeI];
                checkEdge[1] = thisEdge;
            }
            else
            if
            (
                (
                    checkEdge[0].end() == thisEdge[0]
                 || checkEdge[0].end() == thisEdge[1]
                )
            )
            {
                checkEdgeIndex[2] = fEdges[edgeI];
                checkEdge[2] = thisEdge;
            }
            else
            {
                checkEdgeIndex[3] = fEdges[edgeI];
                checkEdge[3] = thisEdge;
            }
        }
    }

    // Determine if either edge belongs to a boundary
    bool firstEdgeBoundary  = (whichEdgePatch(checkEdgeIndex[1]) > -1);
    bool secondEdgeBoundary = (whichEdgePatch(checkEdgeIndex[2]) > -1);

    // Build a hull of cells and tri-faces that are connected to each edge
    labelHashSet firstHullCells, secondHullCells;
    labelHashSet firstHullTriFaces, secondHullTriFaces;

    constructPrismHull
    (
        checkEdgeIndex[1],
        firstHullTriFaces,
        firstHullCells
    );

    constructPrismHull
    (
        checkEdgeIndex[2],
        secondHullTriFaces,
        secondHullCells
    );

    // Obtain lists from hashSets
    labelList firstCells = firstHullCells.toc();
    labelList secondCells = secondHullCells.toc();
    labelList firstTriFaces = firstHullTriFaces.toc();
    labelList secondTriFaces = secondHullTriFaces.toc();

    // Obtain references to edgeFaces
    labelList& firstEdgeFaces = edgeFaces_[checkEdgeIndex[1]];
    labelList& secondEdgeFaces = edgeFaces_[checkEdgeIndex[2]];

    if (debug > 2)
    {
        Info << endl;
        Info << "-------------------------" << endl;
        Info << "Hulls before modification" << endl;
        Info << "-------------------------" << endl;

        Info << nl << "Cells belonging to first Edge Hull: "
             << firstCells << endl;

        forAll(firstCells,cellI)
        {
            cell &firstCurCell = cells_[firstCells[cellI]];

            Info << "Cell: " << firstCells[cellI]
                 << ": " << firstCurCell << endl;

            forAll(firstCurCell,faceI)
            {
                Info << firstCurCell[faceI]
                     << ": " << faces_[firstCurCell[faceI]] << endl;
            }
        }

        Info << nl << "First Edge Face Hull: " << firstEdgeFaces << endl;

        forAll(firstEdgeFaces,indexI)
        {
            Info << firstEdgeFaces[indexI]
                 << ": " << faces_[firstEdgeFaces[indexI]] << endl;
        }

        Info << nl << "Cells belonging to second Edge Hull: "
             << secondCells << endl;

        forAll(secondCells, cellI)
        {
            cell &secondCurCell = cells_[secondCells[cellI]];

            Info << "Cell: " << secondCells[cellI]
                 << ": " << secondCurCell << endl;

            forAll(secondCurCell, faceI)
            {
                Info << secondCurCell[faceI]
                     << ": " << faces_[secondCurCell[faceI]] << endl;
            }
        }

        Info << nl << "Second Edge Face Hull: " << secondEdgeFaces << endl;

        forAll(secondEdgeFaces, indexI)
        {
            Info << secondEdgeFaces[indexI]
                 << ": " << faces_[secondEdgeFaces[indexI]] << endl;
        }
    }

    // Determine the common vertices for the first and second edges
    label cv0 = checkEdge[1].commonVertex(checkEdge[0]);
    label cv1 = checkEdge[1].commonVertex(checkEdge[3]);
    label cv2 = checkEdge[2].commonVertex(checkEdge[0]);
    label cv3 = checkEdge[2].commonVertex(checkEdge[3]);

    // Determine the neighbouring cells
    label c0 = owner_[fIndex], c1 = neighbour_[fIndex];

    // Find the prism-faces
    FixedList<label,2> faceToKeep(0), faceToThrow(0);

    findPrismFaces
    (
        fIndex,
        c0,
        c0BdyFace,
        c0BdyIndex,
        c0IntFace,
        c0IntIndex
    );

    if (c1 != -1)
    {
        findPrismFaces
        (
            fIndex,
            c1,
            c1BdyFace,
            c1BdyIndex,
            c1IntFace,
            c1IntIndex
        );
    }

    // Collapse preferentially towards a symmetryPlane.
    if (firstEdgeBoundary && secondEdgeBoundary)
    {
        WarningIn
        (
            "dynamicTopoFvMesh::collapseQuadFace"
            "(const label, face&)"
        )   << "Collapsing a face that lies on two boundary patches. "
            << "Algorithm will look for a symmetryPlane and collapse "
            << "the face preferentially towards it.\n"
            << "Face: " << fIndex << ": " << thisFace << endl;

        if
        (
            boundaryMesh()[whichEdgePatch(checkEdgeIndex[1])].type()
            == "symmetryPlane"
        )
        {
            secondEdgeBoundary = false;
        }

        if (secondEdgeBoundary)
        {
            if
            (
                boundaryMesh()[whichEdgePatch(checkEdgeIndex[2])].type()
                == "symmetryPlane"
            )
            {
                firstEdgeBoundary = false;
            }
        }
    }

    if (!firstEdgeBoundary && secondEdgeBoundary)
    {
        // Check whether the collapse is possible.
        forAll(firstTriFaces, indexI)
        {
            if
            (
                (firstTriFaces[indexI] == c0BdyIndex[0])
             || (firstTriFaces[indexI] == c0BdyIndex[1])
            )
            {
                continue;
            }

            if (c1 != -1)
            {
                if
                (
                    (firstTriFaces[indexI] == c1BdyIndex[0])
                 || (firstTriFaces[indexI] == c1BdyIndex[1])
                )
                {
                    continue;
                }
            }

            face &triFace = faces_[firstTriFaces[indexI]];

            forAll(triFace, pointI)
            {
                tmpTriFace[pointI] = triFace[pointI];

                if (triFace[pointI] == cv0)
                {
                    tmpTriFace[pointI] = cv2;
                }

                if (triFace[pointI] == cv1)
                {
                    tmpTriFace[pointI] = cv3;
                }
            }

            // Compute the area and check if it's zero/negative
            scalar origArea = triFaceArea(triFace);
            scalar newArea  = triFaceArea(tmpTriFace);

            if
            (
                (Foam::sign(origArea) != Foam::sign(newArea))
             || (mag(newArea) < VSMALL)
            )
            {
                return;
            }
        }

        // Collapse to the second node...
        forAll(firstEdgeFaces,faceI)
        {
            face& replacementFace = faces_[firstEdgeFaces[faceI]];
            replaceLabel(cv0,cv2,replacementFace);
            replaceLabel(cv1,cv3,replacementFace);

            // Determine the quad-face in cell[0] & cell[1]
            // that belongs to firstEdgeFaces
            if (firstEdgeFaces[faceI] == c0IntIndex[0])
            {
                faceToKeep[0]  = c0IntIndex[1];
                faceToThrow[0] = c0IntIndex[0];
            }

            if (firstEdgeFaces[faceI] == c0IntIndex[1])
            {
                faceToKeep[0]  = c0IntIndex[0];
                faceToThrow[0] = c0IntIndex[1];
            }

            if (c1 != -1)
            {
                if (firstEdgeFaces[faceI] == c1IntIndex[0])
                {
                    faceToKeep[1]  = c1IntIndex[1];
                    faceToThrow[1] = c1IntIndex[0];
                }

                if (firstEdgeFaces[faceI] == c1IntIndex[1])
                {
                    faceToKeep[1]  = c1IntIndex[0];
                    faceToThrow[1] = c1IntIndex[1];
                }
            }
        }

        // All triangular boundary faces also need to have point labels replaced
        forAll(firstCells,cellI)
        {
            cell& cellToCheck = cells_[firstCells[cellI]];

            forAll(cellToCheck,faceI)
            {
                face& faceToCheck = faces_[cellToCheck[faceI]];

                if (faceToCheck.size() == 3)
                {
                    forAll(faceToCheck,pointI)
                    {
                        if (faceToCheck[pointI] == cv0)
                        {
                            faceToCheck[pointI] = cv2;
                        }

                        if (faceToCheck[pointI] == cv1)
                        {
                            faceToCheck[pointI] = cv3;
                        }
                    }
                }
            }
        }

        // Delete the two points...
        meshPoints_.remove(cv0);
        meshPoints_.remove(cv1);
        nPoints_ -= 2;

        // Update the reverse point map
        if (cv0 < nOldPoints_)
        {
            reversePointMap_[cv0] = -1;
        }

        if (cv1 < nOldPoints_)
        {
            reversePointMap_[cv1] = -1;
        }
    }
    else
    {
        // Check whether the collapse is possible.
        forAll(secondTriFaces, indexI)
        {
            if
            (
                (secondTriFaces[indexI] == c0BdyIndex[0])
             || (secondTriFaces[indexI] == c0BdyIndex[1])
            )
            {
                continue;
            }

            if (c1 != -1)
            {
                if
                (
                    (secondTriFaces[indexI] == c1BdyIndex[0])
                 || (secondTriFaces[indexI] == c1BdyIndex[1])
                )
                {
                    continue;
                }
            }

            face &triFace = faces_[secondTriFaces[indexI]];

            forAll(triFace, pointI)
            {
                tmpTriFace[pointI] = triFace[pointI];

                if (triFace[pointI] == cv2)
                {
                    tmpTriFace[pointI] = cv0;
                }

                if (triFace[pointI] == cv3)
                {
                    tmpTriFace[pointI] = cv1;
                }
            }

            // Compute the area and check if it's zero/negative
            scalar origArea = triFaceArea(triFace);
            scalar newArea  = triFaceArea(tmpTriFace);

            if
            (
                (Foam::sign(origArea) != Foam::sign(newArea))
             || (mag(newArea) < VSMALL)
            )
            {
                return;
            }
        }

        // Collapse to the first node by default...
        forAll(secondEdgeFaces,faceI)
        {
            face& replacementFace = faces_[secondEdgeFaces[faceI]];
            replaceLabel(cv2, cv0, replacementFace);
            replaceLabel(cv3, cv1, replacementFace);

            // Determine the quad-face(s) in cell[0] & cell[1]
            // that belongs to secondEdgeFaces
            if (secondEdgeFaces[faceI] == c0IntIndex[0])
            {
                faceToKeep[0]  = c0IntIndex[1];
                faceToThrow[0] = c0IntIndex[0];
            }

            if (secondEdgeFaces[faceI] == c0IntIndex[1])
            {
                faceToKeep[0]  = c0IntIndex[0];
                faceToThrow[0] = c0IntIndex[1];
            }

            if (c1 != -1)
            {
                if (secondEdgeFaces[faceI] == c1IntIndex[0])
                {
                    faceToKeep[1]  = c1IntIndex[1];
                    faceToThrow[1] = c1IntIndex[0];
                }

                if (secondEdgeFaces[faceI] == c1IntIndex[1])
                {
                    faceToKeep[1]  = c1IntIndex[0];
                    faceToThrow[1] = c1IntIndex[1];
                }
            }
        }

        // All triangular boundary faces also need to have point labels replaced
        forAll(secondCells, cellI)
        {
            cell& cellToCheck = cells_[secondCells[cellI]];
            forAll(cellToCheck, faceI)
            {
                face& faceToCheck = faces_[cellToCheck[faceI]];
                if (faceToCheck.size() == 3)
                {
                    forAll(faceToCheck, pointI)
                    {
                        if (faceToCheck[pointI] == cv2)
                        {
                            faceToCheck[pointI] = cv0;
                        }

                        if (faceToCheck[pointI] == cv3)
                        {
                            faceToCheck[pointI] = cv1;
                        }
                    }
                }
            }
        }

        // Delete the two points...
        meshPoints_.remove(cv2);
        meshPoints_.remove(cv3);
        nPoints_ -= 2;

        // Update the reverse point map
        if (cv2 < nOldPoints_) reversePointMap_[cv2] = -1;
        if (cv3 < nOldPoints_) reversePointMap_[cv3] = -1;
    }

    if (debug > 2)
    {
        Info << endl;
        Info << "------------------------" << endl;
        Info << "Hulls after modification" << endl;
        Info << "------------------------" << endl;

        Info << nl << "Cells belonging to first Edge Hull: "
             << firstCells << endl;

        forAll(firstCells, cellI)
        {
            cell &firstCurCell = cells_[firstCells[cellI]];

            Info << "Cell: " << firstCells[cellI]
                 << ": " << firstCurCell << endl;

            forAll(firstCurCell, faceI)
            {
                Info << firstCurCell[faceI]
                     << ": " << faces_[firstCurCell[faceI]] << endl;
            }
        }

        Info << nl << "First Edge Face Hull: " << firstEdgeFaces << endl;

        forAll(firstEdgeFaces, indexI)
        {
            Info << firstEdgeFaces[indexI]
                 << ": " << faces_[firstEdgeFaces[indexI]] << endl;
        }

        Info << nl << "Cells belonging to second Edge Hull: "
             << secondCells << endl;

        forAll(secondCells, cellI)
        {
            cell &secondCurCell = cells_[secondCells[cellI]];

            Info << "Cell: " << secondCells[cellI]
                 << ": " << secondCurCell << endl;

            forAll(secondCurCell, faceI)
            {
                Info << secondCurCell[faceI]
                     << ": " << faces_[secondCurCell[faceI]] << endl;
            }
        }

        Info << nl << "Second Edge Face Hull: " << secondEdgeFaces << endl;

        forAll(secondEdgeFaces, indexI)
        {
            Info << secondEdgeFaces[indexI]
                 << ": " << faces_[secondEdgeFaces[indexI]] << endl;
        }

        Info << endl;

        Info << "Retained face: "
             << faceToKeep[0] << ": "
             << " owner: " << owner_[faceToKeep[0]]
             << " neighbour: " << neighbour_[faceToKeep[0]]
             << endl;

        Info << "Discarded face: "
             << faceToThrow[0] << ": "
             << " owner: " << owner_[faceToThrow[0]]
             << " neighbour: " << neighbour_[faceToThrow[0]]
             << endl;

        if (c1 != -1)
        {
            Info << "Retained face: "
                 << faceToKeep[1] << ": "
                 << " owner: " << owner_[faceToKeep[1]]
                 << " neighbour: " << neighbour_[faceToKeep[1]]
                 << endl;

            Info << "Discarded face: "
                 << faceToThrow[1] << ": "
                 << " owner: " << owner_[faceToThrow[1]]
                 << " neighbour: " << neighbour_[faceToThrow[1]]
                 << endl;
        }
    }

    // Ensure proper orientation for the two retained faces
    FixedList<label,2> cellCheck(0);

    if (owner_[faceToThrow[0]] == c0)
    {
        cellCheck[0] = neighbour_[faceToThrow[0]];

        if (owner_[faceToKeep[0]] == c0)
        {
            if
            (
                (neighbour_[faceToThrow[0]] > neighbour_[faceToKeep[0]])
             && (neighbour_[faceToKeep[0]] != -1)
            )
            {
                // This face is to be flipped
                faces_[faceToKeep[0]] = faces_[faceToKeep[0]].reverseFace();
                owner_[faceToKeep[0]] = neighbour_[faceToKeep[0]];
                neighbour_[faceToKeep[0]] = neighbour_[faceToThrow[0]];
            }
            else
            {
                if (neighbour_[faceToThrow[0]] != -1)
                {
                    // Keep orientation intact, and update the owner
                    owner_[faceToKeep[0]] = neighbour_[faceToThrow[0]];
                }
                else
                {
                    // This face will need to be flipped and converted
                    // to a boundary face. Flip it now, so that conversion
                    // happens later.
                    faces_[faceToKeep[0]] = faces_[faceToKeep[0]].reverseFace();
                    owner_[faceToKeep[0]] = neighbour_[faceToKeep[0]];
                    neighbour_[faceToKeep[0]] = -1;
                }
            }
        }
        else
        {
            // Keep orientation intact, and update the neighbour
            neighbour_[faceToKeep[0]] = neighbour_[faceToThrow[0]];
        }
    }
    else
    {
        cellCheck[0] = owner_[faceToThrow[0]];

        if (neighbour_[faceToKeep[0]] == c0)
        {
            if (owner_[faceToThrow[0]] < owner_[faceToKeep[0]])
            {
                // This face is to be flipped
                faces_[faceToKeep[0]] = faces_[faceToKeep[0]].reverseFace();
                neighbour_[faceToKeep[0]] = owner_[faceToKeep[0]];
                owner_[faceToKeep[0]] = owner_[faceToThrow[0]];
            }
            else
            {
                // Keep orientation intact, and update the neighbour
                neighbour_[faceToKeep[0]] = owner_[faceToThrow[0]];
            }
        }
        else
        {
            // Keep orientation intact, and update the owner
            owner_[faceToKeep[0]] = owner_[faceToThrow[0]];
        }
    }

    if (c1 != -1)
    {
        if (owner_[faceToThrow[1]] == c1)
        {
            cellCheck[1] = neighbour_[faceToThrow[1]];

            if (owner_[faceToKeep[1]] == c1)
            {
                if
                (
                    (neighbour_[faceToThrow[1]] > neighbour_[faceToKeep[1]])
                 && (neighbour_[faceToKeep[1]] != -1)
                )
                {
                    // This face is to be flipped
                    faces_[faceToKeep[1]] = faces_[faceToKeep[1]].reverseFace();
                    owner_[faceToKeep[1]] = neighbour_[faceToKeep[1]];
                    neighbour_[faceToKeep[1]] = neighbour_[faceToThrow[1]];
                }
                else
                {
                    if (neighbour_[faceToThrow[1]] != -1)
                    {
                        // Keep orientation intact, and update the owner
                        owner_[faceToKeep[1]] = neighbour_[faceToThrow[1]];
                    }
                    else
                    {
                        // This face will need to be flipped and converted
                        // to a boundary face. Flip it now, so that conversion
                        // happens later.
                        faces_[faceToKeep[1]] =
                            faces_[faceToKeep[1]].reverseFace();
                        owner_[faceToKeep[1]] = neighbour_[faceToKeep[1]];
                        neighbour_[faceToKeep[1]] = -1;
                    }
                }
            }
            else
            {
                // Keep orientation intact, and update the neighbour
                neighbour_[faceToKeep[1]] = neighbour_[faceToThrow[1]];
            }
        }
        else
        {
            cellCheck[1] = owner_[faceToThrow[1]];
            if (neighbour_[faceToKeep[1]] == c1)
            {
                if (owner_[faceToThrow[1]] < owner_[faceToKeep[1]])
                {
                    // This face is to be flipped
                    faces_[faceToKeep[1]] = faces_[faceToKeep[1]].reverseFace();
                    neighbour_[faceToKeep[1]] = owner_[faceToKeep[1]];
                    owner_[faceToKeep[1]] = owner_[faceToThrow[1]];
                }
                else
                {
                    // Keep orientation intact, and update the neighbour
                    neighbour_[faceToKeep[1]] = owner_[faceToThrow[1]];
                }
            }
            else
            {
                // Keep orientation intact, and update the owner
                owner_[faceToKeep[1]] = owner_[faceToThrow[1]];
            }
        }
    }

    // Remove orphaned faces
    if (owner_[faceToKeep[0]] == -1)
    {
        removeFace(faceToKeep[0]);
    }
    else
    if
    (
        (neighbour_[faceToKeep[0]] == -1)
     && (whichPatch(faceToKeep[0]) < 0)
    )
    {
        // This face is being converted from interior to boundary. Remove
        // from the interior list and add as a boundary face to the end.
        label newFaceIndex = insertFace
                             (
                                 whichPatch(faceToThrow[0]),
                                 faces_[faceToKeep[0]],
                                 owner_[faceToKeep[0]],
                                 -1
                             );

        replaceLabel
        (
            faceToKeep[0],
            newFaceIndex,
            cells_[owner_[faceToKeep[0]]]
        );

        // Renumber the neighbour so that this face is removed correctly.
        neighbour_[faceToKeep[0]] = 0;
        removeFace(faceToKeep[0]);
    }

    // Remove the unwanted faces in the cell(s) adjacent to this face,
    // and correct the cells that contain discarded faces
    cell &cell_0 = cells_[c0];

    forAll(cell_0,faceI)
    {
        if (cell_0[faceI] != fIndex && cell_0[faceI] != faceToKeep[0])
        {
           removeFace(cell_0[faceI]);
        }
    }

    cells_.remove(c0);
    lengthScale_.remove(c0);

    if (cellCheck[0] != -1)
    {
        replaceLabel(faceToThrow[0], faceToKeep[0], cells_[cellCheck[0]]);
    }

    // Update the number of cells, and the reverse map
    nCells_--;
    if (c0 < nOldCells_)
    {
        reverseCellMap_[c0] = -1;
    }

    // Check if the cell was added in the current morph, and delete
    if (cellsFromCells_.found(c0))
    {
        cellsFromCells_.erase(c0);
    }

    if (c1 != -1)
    {
        // Remove orphaned faces
        if (owner_[faceToKeep[1]] == -1)
        {
            removeFace(faceToKeep[1]);
        }
        else
        if
        (
            (neighbour_[faceToKeep[1]] == -1)
         && (whichPatch(faceToKeep[1]) < 0)
        )
        {
            // This face is being converted from interior to boundary. Remove
            // from the interior list and add as a boundary face to the end.
            label newFaceIndex = insertFace
                                 (
                                     whichPatch(faceToThrow[1]),
                                     faces_[faceToKeep[1]],
                                     owner_[faceToKeep[1]],
                                     -1
                                 );

            replaceLabel
            (
                faceToKeep[1],
                newFaceIndex,
                cells_[owner_[faceToKeep[1]]]
            );

            // Renumber the neighbour so that this face is removed correctly.
            neighbour_[faceToKeep[1]] = 0;
            removeFace(faceToKeep[1]);
        }

        cell &cell_1 = cells_[c1];

        forAll(cell_1, faceI)
        {
            if (cell_1[faceI] != fIndex && cell_1[faceI] != faceToKeep[1])
            {
               removeFace(cell_1[faceI]);
            }
        }

        cells_.remove(c1);
        lengthScale_.remove(c1);

        if (cellCheck[1] != -1)
        {
            replaceLabel(faceToThrow[1], faceToKeep[1], cells_[cellCheck[1]]);
        }

        // Update the number of cells, and the reverse map
        nCells_--;
        if (c1 < nOldCells_)
        {
            reverseCellMap_[c1] = -1;
        }

        // Check if the cell was added in the current morph, and delete
        if (cellsFromCells_.found(c1))
        {
            cellsFromCells_.erase(c1);
        }
    }

    // Finally remove the face
    removeFace(fIndex);

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    nCollapses_++;

    // Increment the number of modifications
    nModifications_++;
}

// Method for the collapse of an edge in 3D
void dynamicTopoFvMesh::collapseEdge
(
    const label eIndex
)
{
    // Edge collapse performs the following operations:
    //      [1] Checks if either vertex of the edge is on a boundary
    //      [2] Checks whether cells attached to deleted vertices will be valid
    //          after the edge-collapse operation
    //      [3] Deletes all cells surrounding this edge
    //      [4] Deletes all faces surrounding this edge
    //      [5] Deletes all faces surrounding the deleted vertex attached
    //          to the cells in [3]
    //      [6] Checks the orientation of faces connected to the retained
    //          vertices
    //      [7] Remove one of the vertices of the edge
    //      Update faceEdges, edgeFaces and edgePoints information

    // Figure out which thread this is...
    label tIndex = self();

    // Try to write-lock this edge.
    if (tryEdgeLock(eIndex, rwMutex::WRITE_LOCK))
    {
        edgeStack(tIndex).push(eIndex);
        return;
    }

    if
    (
        (nModifications_ > maxModifications_)
     && (maxModifications_ > -1)
    )
    {
        // Reached the max allowable topo-changes.
        edgeStack(tIndex).clear();
        return;
    }

    // Hull variables
    bool found = false;
    edge& thisEdge = edges_[eIndex];
    labelList& vertexHull = edgePoints_[eIndex];
    label replaceIndex = -1, m = vertexHull.size();
    FixedList<bool,2> edgeBoundary(false);

    // Size up the hull lists
    labelList cellHull(m, -1);
    labelList faceHull(m, -1);
    labelList edgeHull(m, -1);
    labelListList ringEntities(4, labelList(m, -1));

    // Construct a hull around this edge, and write-lock entities
    if
    (
        constructHull
        (
            eIndex,
            edgeHull,
            faceHull,
            cellHull,
            ringEntities,
            rwMutex::WRITE_LOCK
        )
    )
    {
        // Put this edge back on the stack and bail out
        edgeStack(tIndex).push(eIndex);
        return;
    }

    if (debug > 1)
    {
        Info << nl << nl << "Edge: " << eIndex
             << ": " << thisEdge << " is to be collapsed. " << endl;

        // Write out VTK files prior to change
        if (debug > 3)
        {
            writeVTK
            (
                Foam::name(eIndex)+"Collapse_0",
                cellHull
            );
        }
    }

    // Check whether points of the edge lies on a boundary
    checkEdgeBoundary(eIndex, edgeBoundary);

    // Configure the new point-position
    point newPoint = vector::zero;

    // Decide which point to remove
    FixedList<label,2> checkPoints(-1);
    label collapseCase = -1;
    label collapsePoint = -1, replacePoint = -1;
    label removeEdgeIndex = -1, removeFaceIndex = -1;
    label replaceEdgeIndex = -1, replaceFaceIndex = -1;

    if (edgeBoundary[0] && !edgeBoundary[1])
    {
        collapseCase = 1;
    }
    else
    if (!edgeBoundary[0] && edgeBoundary[1])
    {
        collapseCase = 2;
    }
    else
    if (edgeBoundary[0] && edgeBoundary[1])
    {
        // Looks like both points are on the boundary.
        // Check if either point touches a boundary face, and retain that.
        FixedList<bool,2> faceCheck(false);

        forAll(ringEntities[1], faceI)
        {
            if (whichPatch(ringEntities[1][faceI]) > -1)
            {
                faceCheck[0] = true;
                break;
            }
        }

        forAll(ringEntities[3], faceI)
        {
            if (whichPatch(ringEntities[3][faceI]) > -1)
            {
                faceCheck[1] = true;
                break;
            }
        }

        if (faceCheck[0] && !faceCheck[1])
        {
            collapseCase = 1;
        }
        else
        if (!faceCheck[0] && faceCheck[1])
        {
            collapseCase = 2;
        }
        else
        {
            collapseCase = 2;
        }
    }
    else
    {
        // Looks like this is an interior edge.
        // Collapse case [2] by default
        collapseCase = 2;
    }

    switch (collapseCase)
    {
        case 1:

            // Collapse to the first node
            replacePoint = thisEdge[0];
            collapsePoint = thisEdge[1];
            replaceEdgeIndex = 0;
            replaceFaceIndex = 1;
            removeEdgeIndex = 2;
            removeFaceIndex = 3;
            checkPoints[0] = collapsePoint;
            newPoint = meshPoints_[thisEdge[0]];

            break;

        case 2:

            // Collapse to the second node
            replacePoint = thisEdge[1];
            collapsePoint = thisEdge[0];
            removeEdgeIndex = 0;
            removeFaceIndex = 1;
            replaceEdgeIndex = 2;
            replaceFaceIndex = 3;
            checkPoints[0] = collapsePoint;
            newPoint = meshPoints_[thisEdge[1]];

            break;

        default:

            // Don't think this will ever happen.
            FatalErrorIn("dynamicTopoFvMesh::collapseEdge()")
                << "Edge: " << eIndex << ": " << thisEdge
                << ". Couldn't decide on collapseCase."
                << abort(FatalError);

            break;
    }

    if (debug > 2)
    {
        Info << "Vertices: " << vertexHull << endl;
        Info << "Edges: " << edgeHull << endl;
        Info << "Faces: " << faceHull << endl;
        Info << "Cells: " << cellHull << endl;
        Info << "replacePoint: " << replacePoint << endl;
        Info << "collapsePoint: " << collapsePoint << endl;
        Info << "ringEntities (removed faces): " << endl;

        forAll(ringEntities[removeFaceIndex], faceI)
        {
            label fIndex = ringEntities[removeFaceIndex][faceI];

            if (fIndex == -1)
            {
                continue;
            }

            Info << fIndex << ": " << faces_[fIndex] << endl;
        }

        Info << "ringEntities (removed edges): " << endl;
        forAll(ringEntities[removeEdgeIndex], edgeI)
        {
            label ieIndex = ringEntities[removeEdgeIndex][edgeI];

            if (ieIndex == -1)
            {
                continue;
            }

            Info << ieIndex << ": " << edges_[ieIndex] << endl;
        }

        Info << "ringEntities (replacement faces): " << endl;
        forAll(ringEntities[replaceFaceIndex], faceI)
        {
            label fIndex = ringEntities[replaceFaceIndex][faceI];

            if (fIndex == -1)
            {
                continue;
            }

            Info << fIndex << ": " << faces_[fIndex] << endl;
        }

        Info << "ringEntities (replacement edges): " << endl;
        forAll(ringEntities[replaceEdgeIndex], edgeI)
        {
            label ieIndex = ringEntities[replaceEdgeIndex][edgeI];

            if (ieIndex == -1)
            {
                continue;
            }

            Info << ieIndex << ": " << edges_[ieIndex] << endl;
        }

        labelList& collapsePointEdges = pointEdges_[collapsePoint];

        Info << "pointEdges (collapsePoint): ";

        forAll(collapsePointEdges, edgeI)
        {
            Info << collapsePointEdges[edgeI] << " ";
        }

        Info << endl;
    }

    // Loop through edges and check for feasibility of collapse
    labelHashSet cellsChecked;

    // Add all hull cells as 'checked'
    forAll(cellHull, cellI)
    {
        if (cellHull[cellI] == -1)
        {
            continue;
        }

        cellsChecked.insert(cellHull[cellI]);
    }

    // Check collapsibility of cells around edges with the re-configured point
    forAll(checkPoints, pointI)
    {
        if (checkPoints[pointI] == -1)
        {
            continue;
        }

        labelList& checkPointEdges = pointEdges_[checkPoints[pointI]];

        forAll(checkPointEdges, edgeI)
        {
            labelList& eFaces = edgeFaces_[checkPointEdges[edgeI]];

            // Build a list of cells to check
            forAll(eFaces, faceI)
            {
                label own = owner_[eFaces[faceI]];
                label nei = neighbour_[eFaces[faceI]];

                // Check owner cell
                if (!cellsChecked.found(own))
                {
                    // Check if a collapse is feasible
                    if
                    (
                        checkCollapse
                        (
                            newPoint,
                            checkPoints[pointI],
                            own,
                            edgeBoundary,
                            cellsChecked
                        )
                    )
                    {
                        // Unlock all entities
                        unlockMutexLists(tIndex);

                        return;
                    }
                }

                // Check neighbour cell
                if (!cellsChecked.found(nei) && nei != -1)
                {
                    // Check if a collapse is feasible
                    if
                    (
                        checkCollapse
                        (
                            newPoint,
                            checkPoints[pointI],
                            nei,
                            edgeBoundary,
                            cellsChecked
                        )
                    )
                    {
                        // Unlock all entities
                        unlockMutexLists(tIndex);

                        return;
                    }
                }
            }
        }
    }

    // Renumber all hull faces and edges
    forAll(faceHull, indexI)
    {
        // Loop through all faces of the edge to be removed
        // and reassign them to the replacement edge
        label edgeToRemove = ringEntities[removeEdgeIndex][indexI];
        label faceToRemove = ringEntities[removeFaceIndex][indexI];
        label cellToRemove = cellHull[indexI];
        label replaceEdge = ringEntities[replaceEdgeIndex][indexI];
        label replaceFace = ringEntities[replaceFaceIndex][indexI];

        labelList& rmvEdgeFaces = edgeFaces_[edgeToRemove];
        labelList& rplEdgeFaces = edgeFaces_[replaceEdge];

        // Replace edgePoints for all edges emanating from hullVertices
        // except ring-edges; those are sized-down later
        labelList& hullPointEdges = pointEdges_[vertexHull[indexI]];

        forAll(hullPointEdges, edgeI)
        {
            labelList& hullEdgePoints = edgePoints_[hullPointEdges[edgeI]];

            if
            (
                 foundInList(collapsePoint, hullEdgePoints)
             && !foundInList(replacePoint, hullEdgePoints)
            )
            {
                replaceLabel
                (
                    collapsePoint,
                    replacePoint,
                    hullEdgePoints
                );
            }
        }

        forAll(rmvEdgeFaces, faceI)
        {
            // Replace edge labels for faces
            replaceLabel
            (
                edgeToRemove,
                replaceEdge,
                faceEdges_[rmvEdgeFaces[faceI]]
            );

            // Loop through faces associated with this edge,
            // and renumber them as well.
            face& faceToCheck = faces_[rmvEdgeFaces[faceI]];
            if ((replaceIndex = faceToCheck.which(collapsePoint)) > -1)
            {
                if (debug > 2)
                {
                    Info << "Renumbering face: "
                         << rmvEdgeFaces[faceI] << ": "
                         << faceToCheck << endl;
                }

                faceToCheck[replaceIndex] = replacePoint;
            }

            // Hull faces should be removed for the replacement edge
            if (rmvEdgeFaces[faceI] == faceHull[indexI])
            {
                sizeDownList
                (
                    faceHull[indexI],
                    rplEdgeFaces
                );

                continue;
            }

            found = false;

            // Need to avoid ring faces as well.
            forAll(ringEntities[removeFaceIndex], faceII)
            {
                if
                (
                    rmvEdgeFaces[faceI]
                 == ringEntities[removeFaceIndex][faceII]
                )
                {
                    found = true;
                    break;
                }
            }

            // Size-up the replacement edge list if the face hasn't been found.
            // These faces are connected to the edge slated for
            // removal, but do not belong to the hull.
            if (!found)
            {
                sizeUpList
                (
                    rmvEdgeFaces[faceI],
                    rplEdgeFaces
                );
            }
        }

        if (cellToRemove == -1)
        {
            continue;
        }

        // Size down edgeFaces for the ring edges
        sizeDownList
        (
            faceToRemove,
            edgeFaces_[edgeHull[indexI]]
        );

        // Size down edgePoints for the ring edges
        sizeDownList
        (
            collapsePoint,
            edgePoints_[edgeHull[indexI]]
        );

        // Ensure proper orientation of retained faces
        if (owner_[faceToRemove] == cellToRemove)
        {
            if (owner_[replaceFace] == cellToRemove)
            {
                if
                (
                    (neighbour_[faceToRemove] > neighbour_[replaceFace])
                 && (neighbour_[replaceFace] != -1)
                )
                {
                    // This face is to be flipped
                    faces_[replaceFace] = faces_[replaceFace].reverseFace();
                    owner_[replaceFace] = neighbour_[replaceFace];
                    neighbour_[replaceFace] = neighbour_[faceToRemove];
                }
                else
                {
                    // Keep orientation intact, and update the owner
                    owner_[replaceFace] = neighbour_[faceToRemove];
                }
            }
            else
            {
                // Keep orientation intact, and update the neighbour
                neighbour_[replaceFace] = neighbour_[faceToRemove];
            }

            // Update the cell
            if (neighbour_[faceToRemove] != -1)
            {
                replaceLabel
                (
                    faceToRemove,
                    replaceFace,
                    cells_[neighbour_[faceToRemove]]
                );
            }
        }
        else
        {
            if (neighbour_[replaceFace] == cellToRemove)
            {
                if (owner_[faceToRemove] < owner_[replaceFace])
                {
                    // This face is to be flipped
                    faces_[replaceFace] = faces_[replaceFace].reverseFace();
                    neighbour_[replaceFace] = owner_[replaceFace];
                    owner_[replaceFace] = owner_[faceToRemove];
                }
                else
                {
                    // Keep orientation intact, and update the neighbour
                    neighbour_[replaceFace] = owner_[faceToRemove];
                }
            }
            else
            {
                // Keep orientation intact, and update the owner
                owner_[replaceFace] = owner_[faceToRemove];
            }

            // Update the cell
            if (owner_[faceToRemove] != -1)
            {
                replaceLabel
                (
                    faceToRemove,
                    replaceFace,
                    cells_[owner_[faceToRemove]]
                );
            }
        }

        // Check orientation of faces
        if (owner_[replaceFace] == -1)
        {
            FatalErrorIn("dynamicTopoFvMesh::collapseEdge()")
                << "Face: " << replaceFace << ": " << faces_[replaceFace]
                << " is an orphan, i.e, no owner cell."
                << abort(FatalError);
        }
        else
        if
        (
            (neighbour_[replaceFace] == -1)
         && (whichPatch(replaceFace) < 0)
        )
        {
            FatalErrorIn("dynamicTopoFvMesh::collapseEdge()")
                << "Face: " << replaceFace << faces_[replaceFace]
                << " is being converted to a boundary."
                << abort(FatalError);
        }
    }

    // Write lock mutexes
    eMutex_.lock(rwMutex::WRITE_LOCK);
    fMutex_.lock(rwMutex::WRITE_LOCK);
    cMutex_.lock(rwMutex::WRITE_LOCK);

    // Remove all hull entities
    forAll(faceHull, indexI)
    {
        label edgeToRemove = ringEntities[removeEdgeIndex][indexI];
        label faceToRemove = ringEntities[removeFaceIndex][indexI];
        label cellToRemove = cellHull[indexI];

        if (cellToRemove != -1)
        {
            // Remove faceToRemove and associated faceEdges
            removeFace(faceToRemove);

            // Remove from list of locked faces
            removeFaceLock(faceToRemove);

            // Remove the hull cell
            cells_.remove(cellToRemove);
            lengthScale_.remove(cellToRemove);

            // Remove from list of locked cells
            removeCellLock(cellToRemove);

            // Remove the cell mutex
            if (threader_->multiThreaded())
            {
                // Unlock it first
                cellMutex_[cellToRemove].unlock();
                cellMutex_.remove(cellToRemove);
            }

            // Update the number of cells, and the reverse cell map
            nCells_--;

            if (cellToRemove < nOldCells_)
            {
                reverseCellMap_[cellToRemove] = -1;
            }

            // Check if the cell was added in the current morph, and delete
            if (cellsFromCells_.found(cellToRemove))
            {
                cellsFromCells_.erase(cellToRemove);
            }
        }

        // Remove the hull edge and associated edgeFaces
        removeEdge(edgeToRemove);

        // Remove from list of locked edges
        removeEdgeLock(edgeToRemove);

        // Remove the hull face
        removeFace(faceHull[indexI]);

        // Remove from list of locked faces
        removeFaceLock(faceHull[indexI]);
    }

    // Unlock mutexes from write lock
    cMutex_.unlock();
    fMutex_.unlock();
    eMutex_.unlock();

    // Loop through pointEdges for the collapsePoint,
    // and replace all occurrences with replacePoint.
    // Size-up pointEdges for the replacePoint as well.
    labelList& pEdges = pointEdges_[collapsePoint];

    forAll(pEdges, edgeI)
    {
        // Renumber edges
        edge& edgeToCheck = edges_[pEdges[edgeI]];
        labelList& eFaces = edgeFaces_[pEdges[edgeI]];

        if (pEdges[edgeI] != eIndex)
        {
            if (debug > 2)
            {
                Info << "Renumbering [edge]: "
                     << pEdges[edgeI] << ": "
                     << edgeToCheck << endl;
            }

            if (edgeToCheck[0] == collapsePoint)
            {
                edgeToCheck[0] = replacePoint;

                sizeUpList
                (
                    pEdges[edgeI],
                    pointEdges_[replacePoint]
                );
            }
            else
            if (edgeToCheck[1] == collapsePoint)
            {
                edgeToCheck[1] = replacePoint;

                sizeUpList
                (
                    pEdges[edgeI],
                    pointEdges_[replacePoint]
                );
            }
            else
            {
                // Looks like pointEdges is inconsistent
                FatalErrorIn("dynamicTopoFvMesh::collapseEdge()") << nl
                    << "pointEdges is inconsistent." << nl
                    << "Point: " << collapsePoint << nl
                    << "pointEdges: " << pEdges << nl
                    << abort(FatalError);
            }

            // Loop through faces associated with this edge,
            // and renumber them as well.
            forAll(eFaces, faceI)
            {
                face& faceToCheck = faces_[eFaces[faceI]];

                if ((replaceIndex = faceToCheck.which(collapsePoint)) > -1)
                {
                    if (debug > 2)
                    {
                        Info << "Renumbering face: "
                             << eFaces[faceI] << ": "
                             << faceToCheck << endl;
                    }

                    faceToCheck[replaceIndex] = replacePoint;

                    // Look for an edge on this face that doesn't
                    // contain collapsePoint or replacePoint.
                    label rplIndex = -1;
                    labelList& fEdges = faceEdges_[eFaces[faceI]];

                    forAll(fEdges, edgeI)
                    {
                        edge& eCheck = edges_[fEdges[edgeI]];

                        if
                        (
                            eCheck[0] != collapsePoint
                         && eCheck[1] != collapsePoint
                         && eCheck[0] != replacePoint
                         && eCheck[1] != replacePoint
                        )
                        {
                            rplIndex = fEdges[edgeI];
                            break;
                        }
                    }

                    // Modify edgePoints for this edge
                    replaceLabel
                    (
                        collapsePoint,
                        replacePoint,
                        edgePoints_[rplIndex]
                    );
                }
            }
        }
    }

    // At this point, edgePoints for the replacement edges are broken,
    // but edgeFaces are consistent. So use this information to re-build
    // edgePoints for all replacement edges.
    forAll(ringEntities[replaceEdgeIndex], edgeI)
    {
        if (debug > 2)
        {
            Info << "Building edgePoints for edge: "
                 << ringEntities[replaceEdgeIndex][edgeI] << ": "
                 << edges_[ringEntities[replaceEdgeIndex][edgeI]]
                 << endl;
        }

        buildEdgePoints(ringEntities[replaceEdgeIndex][edgeI]);
    }

    // Write lock the point mutex
    pMutex_.lock(rwMutex::WRITE_LOCK);

    // Move to the new point
    meshPoints_[replacePoint] = newPoint;

    // Remove the collapse point
    meshPoints_.remove(collapsePoint);
    nPoints_--;

    // Remove from list of locked points
    removePointLock(collapsePoint);

    // Remove the point mutex
    if (threader_->multiThreaded())
    {
        // Unlock it first
        pointMutex_[collapsePoint].unlock();
        pointMutex_.remove(collapsePoint);
    }

    // Null pointEdges so that removeEdge deletes it.
    pointEdges_[collapsePoint] = labelList(0);

    // Unlock the point mutex from write lock
    pMutex_.unlock();

    // Update the reverse point map
    if (collapsePoint < nOldPoints_)
    {
        reversePointMap_[collapsePoint] = -1;
    }

    // Write lock the edge mutex
    eMutex_.lock(rwMutex::WRITE_LOCK);

    // Remove the edge
    removeEdge(eIndex);

    // Remove from list of locked edges
    removeEdgeLock(eIndex);

    // Unlock the edge mutex from write lock
    eMutex_.unlock();

    // Unlock all entities (from write lock)
    unlockMutexLists(tIndex);

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    nCollapses_++;

    // Increment the number of modifications
    nModifications_++;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
