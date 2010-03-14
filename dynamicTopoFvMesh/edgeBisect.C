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
#include "interpolator.H"
#include "resizableList.H"
#include "multiThreader.H"
#include "dynamicTopoFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Method for the bisection of a quad-face in 2D
// - Returns a changeMap with a type specifying:
//     1: Bisection was successful
//    -1: Bisection failed since max number of topo-changes was reached.
const changeMap dynamicTopoFvMesh::bisectQuadFace
(
    const label fIndex,
    bool checkOnly,
    const changeMap& masterMap
)
{
    // Quad-face bisection performs the following operations:
    //      [1] Add two points at the middle of the face
    //      [2] Create a new internal face for each bisected cell
    //      [3] Modify existing face and create a new half-face
    //      [4] Modify triangular boundary faces, and create new ones as well
    //      [5] Create edges for new faces
    //      Update faceEdges and edgeFaces information

    // Figure out which thread this is...
    label tIndex = self(), pIndex = -1;

    // Prepare the changeMaps
    changeMap map, slaveMap;
    bool bisectingSlave = false;

    if
    (
        (nModifications_ > maxModifications_) &&
        (maxModifications_ > -1)
    )
    {
        // Reached the max allowable topo-changes.
        faceStack(tIndex).clear();

        return map;
    }

    // Sanity check: Is the index legitimate?
    if (fIndex < 0 || fIndex >= nFaces_)
    {
        FatalErrorIn("dynamicTopoFvMesh::bisectQuadFace()")
            << " Invalid index: " << fIndex
            << abort(FatalError);
    }

    bool found;
    label replaceFace = -1, retainFace = -1;
    face tmpQuadFace(4), tmpTriFace(3);
    labelList tmpQFEdges(4, -1), tmpTFEdges(3, -1);
    FixedList<label,7> newFaceIndex(-1), newEdgeIndex(-1);
    FixedList<edge,4> commonEdges;
    FixedList<label,4> cornerEdgeIndex(-1), commonEdgeIndex(-1);
    FixedList<label,4> commonFaceIndex(-1);
    FixedList<label,2> newPointIndex(-1), newCellIndex(-1);
    FixedList<label,4> otherEdgeIndex(-1), otherEdgePoint(-1);
    FixedList<label,4> otherPointIndex(-1), nextToOtherPoint(-1);
    FixedList<label,2> c0BdyIndex, c0IntIndex, c1BdyIndex, c1IntIndex;
    FixedList<face,2>  c0BdyFace,  c0IntFace,  c1BdyFace,  c1IntFace;

    // Get the two cells on either side...
    label c0 = owner_[fIndex], c1 = neighbour_[fIndex];

    // Keep track of old / new cells
    FixedList<cell, 2> oldCells(cell(5));
    FixedList<cell, 2> newCells(cell(5));

    // Find the prism faces for cell[0].
    oldCells[0] = cells_[c0];

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
             << ": " << faces_[fIndex] << " is to be bisected. " << endl;

        label epIndex = whichPatch(fIndex);

        Info << "Patch: ";

        if (epIndex == -1)
        {
            Info << "Internal" << endl;
        }
        else
        {
            Info << boundaryMesh()[epIndex].name() << endl;
        }

        if (debug > 2)
        {
            Info << "Cell[0]: " << c0 << ": " << oldCells[0] << endl;

            forAll(oldCells[0], faceI)
            {
                const labelList& fE = faceEdges_[oldCells[0][faceI]];

                Info << oldCells[0][faceI] << ": "
                     << faces_[oldCells[0][faceI]]
                     << " fE: " << fE
                     << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Info << '\t' << fE[edgeI]
                         << ": " << edges_[fE[edgeI]]
                         << " eF: " << eF
                         << endl;
                }
            }
        }

        // Write out VTK files prior to change
        if (debug > 3)
        {
            labelList cellHull(2, -1);

            cellHull[0] = owner_[fIndex];
            cellHull[1] = neighbour_[fIndex];

            writeVTK
            (
                Foam::name(fIndex)
              + "_Bisect_0",
                cellHull
            );
        }
    }

    // Find the common-edge between the triangular boundary faces
    // and the face under consideration.
    findCommonEdge(c0BdyIndex[0], fIndex, commonEdgeIndex[0]);
    findCommonEdge(c0BdyIndex[1], fIndex, commonEdgeIndex[1]);

    commonEdges[0] = edges_[commonEdgeIndex[0]];
    commonEdges[1] = edges_[commonEdgeIndex[1]];

    if (coupledModification_)
    {
        // Is this a locally coupled face?
        if (locallyCoupledFace(fIndex))
        {
            label sIndex = -1;

            // Determine the slave index.
            forAll(patchCoupling_, patchI)
            {
                if (patchCoupling_(patchI))
                {
                    const label faceEnum  = coupleMap::FACE;
                    const coupleMap& cMap = patchCoupling_[patchI].patchMap();

                    if ((sIndex = cMap.findSlaveIndex(faceEnum, fIndex)) > -1)
                    {
                        // Keep this index for master/slave mapping.
                        pIndex = patchI;

                        break;
                    }

                    // The following bit happens only during the sliver
                    // exudation process, since slave faces are
                    // usually not added to the coupled face-stack.
                    if ((sIndex = cMap.findMasterIndex(faceEnum, fIndex)) > -1)
                    {
                        // Keep this index for master/slave mapping.
                        pIndex = patchI;

                        // Notice that we are bisecting a slave face.
                        bisectingSlave = true;

                        break;
                    }
                }
            }

            if (sIndex == -1)
            {
                FatalErrorIn("dynamicTopoFvMesh::bisectQuadFace()")
                    << "Coupled maps were improperly specified." << nl
                    << " Slave index not found for: " << nl
                    << " Face: " << fIndex << nl
                    << abort(FatalError);
            }

            // Temporarily turn off coupledModification.
            unsetCoupledModification();
            setSlaveModification();

            // First check the slave for bisection feasibility.
            slaveMap = bisectQuadFace(sIndex, true);

            if (slaveMap.type() == 1)
            {
                // Can the master be bisected as well?
                changeMap masterMap = bisectQuadFace(fIndex, true);

                // Master couldn't perform bisection
                if (masterMap.type() != 1)
                {
                    setCoupledModification();
                    unsetSlaveModification();

                    return masterMap;
                }

                // Fill the masterMap with points that
                // we seek maps for...
                masterMap.addPoint(commonEdges[0].start());
                masterMap.addPoint(commonEdges[1].start());

                // Bisect the slave edge
                slaveMap = bisectQuadFace(sIndex, false, masterMap);
            }
            else
            {
                // Slave couldn't perform collapse.
                setCoupledModification();
                unsetSlaveModification();

                map.type() = -2;

                return map;
            }

            // Turn it back on.
            setCoupledModification();
            unsetSlaveModification();
        }
        else
        {
            // Bisect edge on the patchSubMesh.

        }
    }

    // Are we performing only checks?
    if (checkOnly)
    {
        map.type() = 1;
        return map;
    }

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

    labelList mP(2, -1);

    // Set mapping for this point
    mP[0] = commonEdges[0][0];
    mP[1] = commonEdges[0][1];

    // Add two new points to the end of the list
    newPointIndex[0] =
    (
        insertPoint
        (
            0.5 * (points_[mP[0]] + points_[mP[1]]),
            0.5 * (oldPoints_[mP[0]] + oldPoints_[mP[1]]),
            mP
        )
    );

    // Set mapping for this point
    mP[0] = commonEdges[1][0];
    mP[1] = commonEdges[1][1];

    newPointIndex[1] =
    (
        insertPoint
        (
            0.5 * (points_[mP[0]] + points_[mP[1]]),
            0.5 * (oldPoints_[mP[0]] + oldPoints_[mP[1]]),
            mP
        )
    );

    // Add the points to the map. Since this might require master mapping,
    // first check to see if a slave is being bisected.
    if (slaveModification_)
    {
        const Map<label>& pMap = masterMap.addedPointList();

        // Look through the reverse point map
        // to check which point it corresponds to.
        forAll(patchCoupling_, patchI)
        {
            if (!patchCoupling_(patchI))
            {
                continue;
            }

            const coupleMap& cMap = patchCoupling_[patchI].patchMap();

            Map<label>& rPointMap = cMap.reverseEntityMap(coupleMap::POINT);

            forAll(commonEdgeIndex, edgeI)
            {
                if (commonEdgeIndex[edgeI] == -1)
                {
                    continue;
                }

                forAll(commonEdges[edgeI], pointI)
                {
                    bool found = false;
                    label mPoint = -1;

                    if (rPointMap.found(commonEdges[edgeI][pointI]))
                    {
                        mPoint = rPointMap[commonEdges[edgeI][pointI]];

                        // Now check the master map.
                        if (pMap.found(mPoint))
                        {
                            found = true;
                        }
                    }

                    // Add the map entry for the point.
                    if (found)
                    {
                        map.addPoint
                        (
                            newPointIndex[edgeI],
                            mPoint
                        );

                        break;
                    }
                }
            }
        }
    }
    else
    {
        map.addPoint(newPointIndex[0]);
        map.addPoint(newPointIndex[1]);
    }

    // Fill-in mapping information
    labelList mC0(1, c0);
    scalarField mW0(1, 1.0);

    // Add a new prism cell to the end of the list.
    // Currently invalid, but will be updated later.
    newCellIndex[0] = insertCell(newCells[0], mC0, mW0, lengthScale_[c0]);

    // Modify the two existing triangle boundary faces

    // Zeroth boundary face - Owner = c[0] & Neighbour [-1] (unchanged)
    replaceLabel
    (
        otherEdgePoint[0],
        newPointIndex[0],
        c0BdyFace[0]
    );

    // First boundary face - Owner = newCell[0], Neighbour = -1
    replaceLabel
    (
        otherEdgePoint[1],
        newPointIndex[1],
        c0BdyFace[1]
    );

    // Update faces.
    faces_[c0BdyIndex[0]] = c0BdyFace[0];
    faces_[c0BdyIndex[1]] = c0BdyFace[1];

    owner_[c0BdyIndex[1]] = newCellIndex[0];
    replaceLabel(c0BdyIndex[1], -1, oldCells[0]);

    // Detect edges other than commonEdges
    const labelList& fEdges = faceEdges_[fIndex];

    forAll(fEdges, edgeI)
    {
        if
        (
            fEdges[edgeI] != commonEdgeIndex[0] &&
            fEdges[edgeI] != commonEdgeIndex[1]
        )
        {
            // Obtain a reference to this edge
            const edge& eThis = edges_[fEdges[edgeI]];

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
        newPointIndex[0],
        faces_[fIndex]
    );

    replaceLabel
    (
        nextToOtherPoint[1],
        newPointIndex[1],
        faces_[fIndex]
    );

    // Add this face to the map.
    // Although this face isn't technically 'added', it's
    // required for coupled patch mapping.
    map.addFace(fIndex);

    // Set face flux for the modified face.
    scalar oldPhi = iPtr_->getPhi(fIndex);

    iPtr_->setPhi(fIndex, 0.5 * oldPhi);

    if (debug > 1)
    {
        Info << "Modified face: " << fIndex
             << ": " << faces_[fIndex] << endl;

        if (debug > 2)
        {
            Info << "Common edges: " << nl
                 << commonEdgeIndex[0] << ": " << commonEdges[0] << nl
                 << commonEdgeIndex[1] << ": " << commonEdges[1]
                 << endl;
        }
    }

    // Find the quad face that contains otherEdgeIndex[1]
    found = false;

    const labelList& e1 = faceEdges_[c0IntIndex[0]];

    forAll(e1, edgeI)
    {
        if (e1[edgeI] == otherEdgeIndex[1])
        {
            replaceLabel(c0IntIndex[0], -1, oldCells[0]);
            replaceFace = c0IntIndex[0];
            retainFace = c0IntIndex[1];
            found = true; break;
        }
    }

    if (!found)
    {
        // The edge was not found before
        replaceLabel(c0IntIndex[1], -1, oldCells[0]);
        replaceFace = c0IntIndex[1];
        retainFace = c0IntIndex[0];
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

            iPtr_->setFlip(replaceFace);
        }
    }
    else
    {
        // Keep owner, but change neighbour
        neighbour_[replaceFace] = newCellIndex[0];
    }

    // Define the faces for the new cell
    newCells[0][0] = c0BdyIndex[1];
    newCells[0][1] = replaceFace;

    // Define the set of new faces and insert...

    // New interior face; Owner = cell[0] & Neighbour = newCell[0]
    tmpQuadFace[0] = otherPointIndex[0];
    tmpQuadFace[1] = newPointIndex[0];
    tmpQuadFace[2] = newPointIndex[1];
    tmpQuadFace[3] = otherPointIndex[1];

    newFaceIndex[0] =
    (
        insertFace
        (
            -1,
            tmpQuadFace,
            c0,
            newCellIndex[0]
        )
    );

    // Add a faceEdges entry as well
    faceEdges_.append(tmpQFEdges);

    // Find the common edge between quad/quad faces...
    findCommonEdge
    (
        c0IntIndex[0],
        c0IntIndex[1],
        otherEdgeIndex[2]
    );

    // ... and size-up edgeFaces for the edge
    sizeUpList
    (
        newFaceIndex[0],
        edgeFaces_[otherEdgeIndex[2]]
    );

    replaceLabel(-1, newFaceIndex[0], newCells[0]);
    replaceLabel(-1, newFaceIndex[0], oldCells[0]);

    // remove2DSliver requires this face index for removal
    map.addFace(newFaceIndex[0]);

    // Second boundary face; Owner = newCell[0] & Neighbour = [-1]
    tmpTriFace[0] = otherPointIndex[0];
    tmpTriFace[1] = newPointIndex[0];
    tmpTriFace[2] = otherEdgePoint[0];

    newFaceIndex[1] =
    (
        insertFace
        (
            whichPatch(c0BdyIndex[0]),
            tmpTriFace,
            newCellIndex[0],
            -1,
            labelList(1, c0BdyIndex[0]),
            scalarField(1, 1.0)
        )
    );

    // Add a faceEdges entry as well
    faceEdges_.append(tmpTFEdges);

    replaceLabel(-1, newFaceIndex[1], newCells[0]);

    // Third boundary face; Owner = c[0] & Neighbour = [-1]
    tmpTriFace[0] = otherPointIndex[1];
    tmpTriFace[1] = newPointIndex[1];
    tmpTriFace[2] = otherEdgePoint[1];

    newFaceIndex[2] =
    (
        insertFace
        (
            whichPatch(c0BdyIndex[1]),
            tmpTriFace,
            c0,
            -1,
            labelList(1, c0BdyIndex[1]),
            scalarField(1, 1.0)
        )
    );

    // Add a faceEdges entry as well
    faceEdges_.append(tmpTFEdges);

    replaceLabel(-1, newFaceIndex[2], oldCells[0]);

    // Create / modify edges...
    labelList tmpTriEdgeFaces(3, -1);

    // The edge bisecting the zeroth boundary triangular face
    tmpTriEdgeFaces[0] = c0BdyIndex[0];
    tmpTriEdgeFaces[1] = newFaceIndex[0];
    tmpTriEdgeFaces[2] = newFaceIndex[1];

    newEdgeIndex[1] =
    (
        insertEdge
        (
            whichEdgePatch(commonEdgeIndex[0]),
            edge(newPointIndex[0], otherPointIndex[0]),
            tmpTriEdgeFaces
        )
    );

    // Find the common edge between the quad/tri faces...
    findCommonEdge
    (
        c0BdyIndex[0],
        replaceFace,
        cornerEdgeIndex[0]
    );

    // ...and correct faceEdges / edgeFaces
    replaceLabel
    (
        cornerEdgeIndex[0],
        newEdgeIndex[1],
        faceEdges_[c0BdyIndex[0]]
    );

    replaceLabel
    (
        c0BdyIndex[0],
        newFaceIndex[1],
        edgeFaces_[cornerEdgeIndex[0]]
    );

    // The edge bisecting the first boundary triangular face
    tmpTriEdgeFaces[0] = c0BdyIndex[1];
    tmpTriEdgeFaces[1] = newFaceIndex[0];
    tmpTriEdgeFaces[2] = newFaceIndex[2];

    newEdgeIndex[2] =
    (
        insertEdge
        (
            whichEdgePatch(commonEdgeIndex[1]),
            edge(newPointIndex[1], otherPointIndex[1]),
            tmpTriEdgeFaces
        )
    );

    // Find the common edge between the quad/tri faces...
    findCommonEdge
    (
        c0BdyIndex[1],
        retainFace,
        cornerEdgeIndex[1]
    );

    // ...and correct faceEdges / edgeFaces
    replaceLabel
    (
        cornerEdgeIndex[1],
        newEdgeIndex[2],
        faceEdges_[c0BdyIndex[1]]
    );

    replaceLabel
    (
        c0BdyIndex[1],
        newFaceIndex[2],
        edgeFaces_[cornerEdgeIndex[1]]
    );

    if (c1 == -1)
    {
        // The quad boundary face resulting from bisection;
        // Owner = newCell[0] & Neighbour = [-1]
        tmpQuadFace[0] = newPointIndex[1];
        tmpQuadFace[1] = nextToOtherPoint[1];
        tmpQuadFace[2] = otherEdgePoint[0];
        tmpQuadFace[3] = newPointIndex[0];

        newFaceIndex[3] =
        (
            insertFace
            (
                whichPatch(fIndex),
                tmpQuadFace,
                newCellIndex[0],
                -1,
                labelList(1, fIndex),
                scalarField(1, 1.0)
            )
        );

        // Add this face to the map.
        map.addFace(newFaceIndex[3]);

        // Add a faceEdges entry as well
        faceEdges_.append(tmpQFEdges);

        // Set the flux
        iPtr_->setPhi(newFaceIndex[3], 0.5 * oldPhi);

        // Correct edgeFaces for otherEdgeIndex[1]
        replaceLabel
        (
            fIndex,
            newFaceIndex[3],
            edgeFaces_[otherEdgeIndex[1]]
        );

        replaceLabel(-1, newFaceIndex[3], newCells[0]);

        labelList tmpBiEdgeFaces(2, -1);

        // The edge bisecting the face
        tmpTriEdgeFaces[0] = newFaceIndex[3];
        tmpTriEdgeFaces[1] = newFaceIndex[0];
        tmpTriEdgeFaces[2] = fIndex;

        newEdgeIndex[0] =
        (
            insertEdge
            (
                whichEdgePatch(otherEdgeIndex[0]),
                edge(newPointIndex[0], newPointIndex[1]),
                tmpTriEdgeFaces
            )
        );

        // Replace an edge on the bisected face
        replaceLabel
        (
            otherEdgeIndex[1],
            newEdgeIndex[0],
            faceEdges_[fIndex]
        );

        // Create / replace side edges created from face bisection
        tmpBiEdgeFaces[0] = newFaceIndex[1];
        tmpBiEdgeFaces[1] = newFaceIndex[3];

        newEdgeIndex[3] =
        (
            insertEdge
            (
                whichEdgePatch(commonEdgeIndex[0]),
                edge(newPointIndex[0], otherEdgePoint[0]),
                tmpBiEdgeFaces
            )
        );

        tmpBiEdgeFaces[0] = c0BdyIndex[1];
        tmpBiEdgeFaces[1] = newFaceIndex[3];

        newEdgeIndex[4] =
        (
            insertEdge
            (
                whichEdgePatch(commonEdgeIndex[1]),
                edge(newPointIndex[1], nextToOtherPoint[1]),
                tmpBiEdgeFaces
            )
        );

        // Now that edges are defined, configure faceEdges
        // for all new faces

        // The quad interior face; Owner = cell[0] & Neighbour = newCell[0]
        tmpQFEdges[0] = newEdgeIndex[0];
        tmpQFEdges[1] = newEdgeIndex[1];
        tmpQFEdges[2] = otherEdgeIndex[2];
        tmpQFEdges[3] = newEdgeIndex[2];
        faceEdges_[newFaceIndex[0]] = tmpQFEdges;

        // Second boundary face; Owner = newCell[0] & Neighbour = [-1]
        tmpTFEdges[0] = newEdgeIndex[3];
        tmpTFEdges[1] = cornerEdgeIndex[0];
        tmpTFEdges[2] = newEdgeIndex[1];
        faceEdges_[newFaceIndex[1]] = tmpTFEdges;

        // Third boundary face; Owner = c[0] & Neighbour = [-1]
        tmpTFEdges[0] = newEdgeIndex[2];
        tmpTFEdges[1] = cornerEdgeIndex[1];
        tmpTFEdges[2] = commonEdgeIndex[1];
        faceEdges_[newFaceIndex[2]] = tmpTFEdges;

        // The quad face from bisection:
        tmpQFEdges[0] = otherEdgeIndex[1];
        tmpQFEdges[1] = newEdgeIndex[3];
        tmpQFEdges[2] = newEdgeIndex[0];
        tmpQFEdges[3] = newEdgeIndex[4];
        faceEdges_[newFaceIndex[3]] = tmpQFEdges;

        replaceLabel
        (
            commonEdgeIndex[1],
            newEdgeIndex[4],
            faceEdges_[c0BdyIndex[1]]
        );

        replaceLabel
        (
            c0BdyIndex[1],
            newFaceIndex[2],
            edgeFaces_[commonEdgeIndex[1]]
        );

        if (debug > 2)
        {
            Info << "Modified Cell[0]: "
                 << c0 << ": " << oldCells[0] << endl;

            forAll(oldCells[0], faceI)
            {
                const labelList& fE = faceEdges_[oldCells[0][faceI]];

                Info << oldCells[0][faceI]
                     << ": " << faces_[oldCells[0][faceI]]
                     << " fE: " << fE
                     << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Info << '\t' << fE[edgeI]
                         << ": " << edges_[fE[edgeI]]
                         << " eF: " << eF
                         << endl;
                }
            }

            Info << "New Cell[0]: " << newCellIndex[0]
                 << ": " << newCells[0] << endl;

            forAll(newCells[0], faceI)
            {
                const labelList& fE = faceEdges_[newCells[0][faceI]];

                Info << newCells[0][faceI]
                     << ": " << faces_[newCells[0][faceI]]
                     << " fE: " << fE
                     << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Info << '\t' << fE[edgeI]
                         << ": " << edges_[fE[edgeI]]
                         << " eF: " << eF
                         << endl;
                }
            }
        }
    }
    else
    {
        oldCells[1] = cells_[c1];

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

        // Fill-in mapping information
        labelList mC1(1, c1);
        scalarField mW1(1, 1.0);

        newCellIndex[1] = insertCell(newCells[1], mC1, mW1, lengthScale_[c1]);

        if (debug > 2)
        {
            Info << "Cell[1]: " << c1 << ": " << oldCells[1] << endl;

            forAll(oldCells[1], faceI)
            {
                const labelList& fE = faceEdges_[oldCells[1][faceI]];

                Info << oldCells[1][faceI] << ": "
                     << faces_[oldCells[1][faceI]]
                     << " fE: " << fE
                     << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Info << '\t' << fE[edgeI]
                         << ": " << edges_[fE[edgeI]]
                         << " eF: " << eF
                         << endl;
                }
            }
        }

        // Find the interior face that contains otherEdgeIndex[1]
        found = false;

        const labelList& e2 = faceEdges_[c1IntIndex[0]];

        forAll(e2, edgeI)
        {
            if (e2[edgeI] == otherEdgeIndex[1])
            {
                replaceLabel(c1IntIndex[0], -1, oldCells[1]);
                replaceFace = c1IntIndex[0];
                retainFace = c1IntIndex[1];
                found = true; break;
            }
        }

        if (!found)
        {
            // The edge was not found before
            replaceLabel(c1IntIndex[1], -1, oldCells[1]);
            replaceFace = c1IntIndex[1];
            retainFace = c1IntIndex[0];
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

                iPtr_->setFlip(replaceFace);
            }
        }
        else
        {
            // Keep owner, but change neighbour
            neighbour_[replaceFace] = newCellIndex[1];
        }

        // Define attributes for the new prism cell
        newCells[1][0] = replaceFace;

        // The quad interior face resulting from bisection;
        // Owner = newCell[0] & Neighbour = newCell[1]
        tmpQuadFace[0] = newPointIndex[1];
        tmpQuadFace[1] = nextToOtherPoint[1];
        tmpQuadFace[2] = otherEdgePoint[0];
        tmpQuadFace[3] = newPointIndex[0];

        newFaceIndex[3] =
        (
            insertFace
            (
                -1,
                tmpQuadFace,
                newCellIndex[0],
                newCellIndex[1]
            )
        );

        // Add a faceEdges entry as well
        faceEdges_.append(tmpQFEdges);

        // Set the flux
        iPtr_->setPhi(newFaceIndex[3], 0.5 * oldPhi);

        // Correct edgeFaces for otherEdgeIndex[1]
        replaceLabel
        (
            fIndex,
            newFaceIndex[3],
            edgeFaces_[otherEdgeIndex[1]]
        );

        replaceLabel(-1, newFaceIndex[3], newCells[0]);
        replaceLabel(-1, newFaceIndex[3], newCells[1]);
        newCells[1][1] = newFaceIndex[3];

        // Check for common edges among the two boundary faces
        // Find the isolated point on both boundary faces of cell[1]
        if
        (
            findCommonEdge(c1BdyIndex[0], c0BdyIndex[0], commonEdgeIndex[2])
        )
        {
            findCommonEdge(c1BdyIndex[1], c0BdyIndex[1], commonEdgeIndex[3]);

            commonFaceIndex[2] = c1BdyIndex[0];
            commonFaceIndex[3] = c1BdyIndex[1];
        }
        else
        {
            findCommonEdge(c1BdyIndex[0], c0BdyIndex[1], commonEdgeIndex[3]);
            findCommonEdge(c1BdyIndex[1], c0BdyIndex[0], commonEdgeIndex[2]);

            commonFaceIndex[2] = c1BdyIndex[1];
            commonFaceIndex[3] = c1BdyIndex[0];
        }

        commonEdges[2] = edges_[commonEdgeIndex[2]];
        commonEdges[3] = edges_[commonEdgeIndex[3]];

        if (debug > 2)
        {
            Info << "Common edges: " << nl
                 << commonEdgeIndex[2] << ": " << commonEdges[2] << nl
                 << commonEdgeIndex[3] << ": " << commonEdges[3]
                 << endl;
        }

        findIsolatedPoint
        (
            faces_[commonFaceIndex[2]],
            commonEdges[2],
            otherPointIndex[2],
            nextToOtherPoint[2]
        );

        findIsolatedPoint
        (
            faces_[commonFaceIndex[3]],
            commonEdges[3],
            otherPointIndex[3],
            nextToOtherPoint[3]
        );

        // For convenience...
        otherEdgePoint[2] = commonEdges[2].otherVertex(nextToOtherPoint[2]);
        otherEdgePoint[3] = commonEdges[3].otherVertex(nextToOtherPoint[3]);

        // Modify the two existing triangle boundary faces

        // Zeroth boundary face - Owner = newCell[1], Neighbour = -1
        replaceLabel
        (
            otherEdgePoint[2],
            newPointIndex[0],
            faces_[commonFaceIndex[2]]
        );

        owner_[commonFaceIndex[2]] = newCellIndex[1];
        replaceLabel(commonFaceIndex[2], -1, oldCells[1]);
        newCells[1][2] = commonFaceIndex[2];

        // First boundary face - Owner = c[1] & Neighbour [-1] (unchanged)
        replaceLabel
        (
            otherEdgePoint[3],
            newPointIndex[1],
            faces_[commonFaceIndex[3]]
        );

        // New interior face; Owner = cell[1] & Neighbour = newCell[1]
        tmpQuadFace[0] = newPointIndex[0];
        tmpQuadFace[1] = otherPointIndex[2];
        tmpQuadFace[2] = otherPointIndex[3];
        tmpQuadFace[3] = newPointIndex[1];

        newFaceIndex[4] =
        (
            insertFace
            (
                -1,
                tmpQuadFace,
                c1,
                newCellIndex[1]
            )
        );

        // Add a faceEdges entry as well
        faceEdges_.append(tmpQFEdges);

        // remove2DSliver requires this face index for removal
        map.addFace(newFaceIndex[4]);

        // Find the common edge between quad/quad faces...
        findCommonEdge
        (
            c1IntIndex[0],
            c1IntIndex[1],
            otherEdgeIndex[3]
        );

        // ... and size-up edgeFaces for the edge
        sizeUpList
        (
            newFaceIndex[4],
            edgeFaces_[otherEdgeIndex[3]]
        );

        replaceLabel(-1, newFaceIndex[4], newCells[1]);
        replaceLabel(-1, newFaceIndex[4], oldCells[1]);

        // Second boundary face; Owner = cell[1] & Neighbour [-1]
        tmpTriFace[0] = otherPointIndex[2];
        tmpTriFace[1] = newPointIndex[0];
        tmpTriFace[2] = otherEdgePoint[2];

        newFaceIndex[5] =
        (
            insertFace
            (
                whichPatch(commonFaceIndex[2]),
                tmpTriFace,
                c1,
                -1,
                labelList(1, commonFaceIndex[2]),
                scalarField(1, 1.0)
            )
        );

        // Add a faceEdges entry as well
        faceEdges_.append(tmpTFEdges);

        replaceLabel(-1, newFaceIndex[5], oldCells[1]);

        // Third boundary face; Owner = newCell[1] & Neighbour [-1]
        tmpTriFace[0] = otherPointIndex[3];
        tmpTriFace[1] = newPointIndex[1];
        tmpTriFace[2] = otherEdgePoint[3];

        newFaceIndex[6] =
        (
            insertFace
            (
                whichPatch(commonFaceIndex[3]),
                tmpTriFace,
                newCellIndex[1],
                -1,
                labelList(1, commonFaceIndex[3]),
                scalarField(1, 1.0)
            )
        );

        // Add a faceEdges entry as well
        faceEdges_.append(tmpTFEdges);

        replaceLabel(-1, newFaceIndex[6], newCells[1]);

        // Create / modify edges...
        labelList tmpQuadEdgeFaces(4, -1);

        // The internal edge bisecting the face
        tmpQuadEdgeFaces[0] = fIndex;
        tmpQuadEdgeFaces[1] = newFaceIndex[0];
        tmpQuadEdgeFaces[2] = newFaceIndex[3];
        tmpQuadEdgeFaces[3] = newFaceIndex[4];

        newEdgeIndex[0] =
        (
            insertEdge
            (
                -1,
                edge(newPointIndex[0], newPointIndex[1]),
                tmpQuadEdgeFaces
            )
        );

        // Replace an edge on the bisected face
        replaceLabel
        (
            otherEdgeIndex[1],
            newEdgeIndex[0],
            faceEdges_[fIndex]
        );

        // Create / replace side edges created from face bisection
        tmpTriEdgeFaces[0] = commonFaceIndex[2];
        tmpTriEdgeFaces[1] = newFaceIndex[3];
        tmpTriEdgeFaces[2] = newFaceIndex[1];

        newEdgeIndex[3] =
        (
            insertEdge
            (
                whichEdgePatch(commonEdgeIndex[2]),
                edge(newPointIndex[0], nextToOtherPoint[2]),
                tmpTriEdgeFaces
            )
        );

        tmpTriEdgeFaces[0] = c0BdyIndex[1];
        tmpTriEdgeFaces[1] = newFaceIndex[3];
        tmpTriEdgeFaces[2] = newFaceIndex[6];

        newEdgeIndex[4] =
        (
            insertEdge
            (
                whichEdgePatch(commonEdgeIndex[3]),
                edge(newPointIndex[1], otherEdgePoint[3]),
                tmpTriEdgeFaces
            )
        );

        // The edge bisecting the second boundary triangular face
        tmpTriEdgeFaces[0] = commonFaceIndex[2];
        tmpTriEdgeFaces[1] = newFaceIndex[4];
        tmpTriEdgeFaces[2] = newFaceIndex[5];

        newEdgeIndex[5] =
        (
            insertEdge
            (
                whichEdgePatch(commonEdgeIndex[2]),
                edge(newPointIndex[0], otherPointIndex[2]),
                tmpTriEdgeFaces
            )
        );

        // Find the common edge between the quad/tri faces...
        findCommonEdge
        (
            commonFaceIndex[2],
            retainFace,
            cornerEdgeIndex[2]
        );

        // ...and correct faceEdges / edgeFaces
        replaceLabel
        (
            cornerEdgeIndex[2],
            newEdgeIndex[5],
            faceEdges_[commonFaceIndex[2]]
        );

        replaceLabel
        (
            commonFaceIndex[2],
            newFaceIndex[5],
            edgeFaces_[cornerEdgeIndex[2]]
        );

        // The edge bisecting the third boundary triangular face
        tmpTriEdgeFaces[0] = commonFaceIndex[3];
        tmpTriEdgeFaces[1] = newFaceIndex[4];
        tmpTriEdgeFaces[2] = newFaceIndex[6];

        newEdgeIndex[6] =
        (
            insertEdge
            (
                whichEdgePatch(commonEdgeIndex[3]),
                edge(newPointIndex[1], otherPointIndex[3]),
                tmpTriEdgeFaces
            )
        );

        // Find the common edge between the quad/tri faces...
        findCommonEdge
        (
            commonFaceIndex[3],
            replaceFace,
            cornerEdgeIndex[3]
        );

        // ...and correct faceEdges / edgeFaces
        replaceLabel
        (
            cornerEdgeIndex[3],
            newEdgeIndex[6],
            faceEdges_[commonFaceIndex[3]]
        );

        replaceLabel
        (
            commonFaceIndex[3],
            newFaceIndex[6],
            edgeFaces_[cornerEdgeIndex[3]]
        );

        // Now that edges are defined, configure faceEdges
        // for all new faces

        // The quad interior face; Owner = c[0] & Neighbour = newCell[0]
        tmpQFEdges[0] = newEdgeIndex[0];
        tmpQFEdges[1] = newEdgeIndex[1];
        tmpQFEdges[2] = otherEdgeIndex[2];
        tmpQFEdges[3] = newEdgeIndex[2];
        faceEdges_[newFaceIndex[0]] = tmpQFEdges;

        // Second boundary face; Owner = newCell[0] & Neighbour = [-1]
        tmpTFEdges[0] = newEdgeIndex[3];
        tmpTFEdges[1] = cornerEdgeIndex[0];
        tmpTFEdges[2] = newEdgeIndex[1];
        faceEdges_[newFaceIndex[1]] = tmpTFEdges;

        // Third boundary face; Owner = c[0] & Neighbour = [-1]
        tmpTFEdges[0] = newEdgeIndex[2];
        tmpTFEdges[1] = cornerEdgeIndex[1];
        tmpTFEdges[2] = commonEdgeIndex[3];
        faceEdges_[newFaceIndex[2]] = tmpTFEdges;

        // The quad face from bisection:
        tmpQFEdges[0] = otherEdgeIndex[1];
        tmpQFEdges[1] = newEdgeIndex[3];
        tmpQFEdges[2] = newEdgeIndex[0];
        tmpQFEdges[3] = newEdgeIndex[4];
        faceEdges_[newFaceIndex[3]] = tmpQFEdges;

        // The quad interior face; Owner = c[1] & Neighbour = newCell[1]
        tmpQFEdges[0] = newEdgeIndex[0];
        tmpQFEdges[1] = newEdgeIndex[5];
        tmpQFEdges[2] = otherEdgeIndex[3];
        tmpQFEdges[3] = newEdgeIndex[6];
        faceEdges_[newFaceIndex[4]] = tmpQFEdges;

        // Second boundary face; Owner = c[1] & Neighbour = [-1]
        tmpTFEdges[0] = commonEdgeIndex[2];
        tmpTFEdges[1] = cornerEdgeIndex[2];
        tmpTFEdges[2] = newEdgeIndex[5];
        faceEdges_[newFaceIndex[5]] = tmpTFEdges;

        // Third boundary face; Owner = newCell[1] & Neighbour = [-1]
        tmpTFEdges[0] = newEdgeIndex[4];
        tmpTFEdges[1] = cornerEdgeIndex[3];
        tmpTFEdges[2] = newEdgeIndex[6];
        faceEdges_[newFaceIndex[6]] = tmpTFEdges;

        replaceLabel
        (
            commonEdgeIndex[1],
            newEdgeIndex[4],
            faceEdges_[c0BdyIndex[1]]
        );

        replaceLabel
        (
            c0BdyIndex[1],
            newFaceIndex[2],
            edgeFaces_[commonEdgeIndex[1]]
        );

        replaceLabel
        (
            commonEdgeIndex[2],
            newEdgeIndex[3],
            faceEdges_[commonFaceIndex[2]]
        );

        replaceLabel
        (
            commonFaceIndex[2],
            newFaceIndex[5],
            edgeFaces_[commonEdgeIndex[2]]
        );

        if (debug > 2)
        {
            Info << nl << "Modified Cell[0]: "
                 << c0 << ": " << oldCells[0] << endl;

            forAll(oldCells[0], faceI)
            {
                const labelList& fE = faceEdges_[oldCells[0][faceI]];

                Info << oldCells[0][faceI]
                     << ": " << faces_[oldCells[0][faceI]]
                     << " fE: " << fE
                     << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Info << '\t' << fE[edgeI]
                         << ": " << edges_[fE[edgeI]]
                         << " eF: " << eF
                         << endl;
                }
            }

            Info << "New Cell[0]: "
                 << newCellIndex[0] << ": " << newCells[0] << endl;

            forAll(newCells[0], faceI)
            {
                const labelList& fE = faceEdges_[newCells[0][faceI]];

                Info << newCells[0][faceI] << ": "
                     << faces_[newCells[0][faceI]]
                     << " fE: " << fE
                     << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Info << '\t' << fE[edgeI]
                         << ": " << edges_[fE[edgeI]]
                         << " eF: " << eF
                         << endl;
                }
            }

            Info << nl << "Modified Cell[1]: "
                 << c1 << ": " << oldCells[1] << endl;

            forAll(oldCells[1], faceI)
            {
                const labelList& fE = faceEdges_[oldCells[1][faceI]];

                Info << oldCells[1][faceI] << ": "
                     << faces_[oldCells[1][faceI]]
                     << " fE: " << fE
                     << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Info << '\t' << fE[edgeI]
                         << ": " << edges_[fE[edgeI]]
                         << " eF: " << eF
                         << endl;
                }
            }

            Info << "New Cell[1]: "
                 << newCellIndex[1] << ": " << newCells[1] << endl;

            forAll(newCells[1], faceI)
            {
                const labelList& fE = faceEdges_[newCells[1][faceI]];

                Info << newCells[1][faceI] << ": "
                     << faces_[newCells[1][faceI]]
                     << " fE: " << fE
                     << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Info << '\t' << fE[edgeI]
                         << ": " << edges_[fE[edgeI]]
                         << " eF: " << eF
                         << endl;
                }
            }
        }

        // Update the cell list.
        cells_[c1] = oldCells[1];
        cells_[newCellIndex[1]] = newCells[1];
    }

    // Update the cell list.
    cells_[c0] = oldCells[0];
    cells_[newCellIndex[0]] = newCells[0];

    // Modify point labels for common edges
    if (edges_[commonEdgeIndex[0]].start() == otherEdgePoint[0])
    {
        edges_[commonEdgeIndex[0]].start() = newPointIndex[0];
    }
    else
    {
        edges_[commonEdgeIndex[0]].end() = newPointIndex[0];
    }

    if (edges_[commonEdgeIndex[1]].start() == nextToOtherPoint[1])
    {
        edges_[commonEdgeIndex[1]].start() = newPointIndex[1];
    }
    else
    {
        edges_[commonEdgeIndex[1]].end() = newPointIndex[1];
    }

    if (coupledModification_)
    {
        // Create a master/slave entry for the new face on the patch.
        if (locallyCoupledFace(fIndex))
        {
            FixedList<bool, 2> foundMatch(false);
            FixedList<label, 2> checkFaces(-1);

            // Add the new point to the coupling map
            const label pointEnum = coupleMap::POINT;

            const Map<label>& apList = slaveMap.addedPointList();

            FixedList<label, 2> mapPointIndex(-1);

            forAllConstIter(Map<label>, apList, apIter)
            {
                if (apIter() == commonEdges[0].start())
                {
                    mapPointIndex[0] = apIter.key();
                }

                if (apIter() == commonEdges[1].start())
                {
                    mapPointIndex[1] = apIter.key();
                }
            }

            // Update pointMap
            Map<label>& pointMap =
            (
                patchCoupling_[pIndex].patchMap().entityMap
                (
                    pointEnum
                )
            );

            pointMap.insert(newPointIndex[0], mapPointIndex[0]);
            pointMap.insert(newPointIndex[1], mapPointIndex[1]);

            // Update reverse pointMap
            Map<label>& rPointMap =
            (
                patchCoupling_[pIndex].patchMap().reverseEntityMap
                (
                    pointEnum
                )
            );

            rPointMap.insert(mapPointIndex[0], newPointIndex[0]);
            rPointMap.insert(mapPointIndex[1], newPointIndex[1]);

            // Fill in the faces to be check for...
            checkFaces[0] = fIndex;
            checkFaces[1] = newFaceIndex[3];

            const Map<label>& afList = slaveMap.addedFaceList();

            forAll (checkFaces, indexI)
            {
                const face& mFace = faces_[checkFaces[indexI]];

                label sFaceIndex = -1;

                forAllConstIter(Map<label>, afList, sfIter)
                {
                    const face& tFace = faces_[sfIter.key()];

                    FixedList <bool, 4> cP(false);

                    forAll(mFace, pointI)
                    {
                        if (tFace.which(pointMap[mFace[pointI]]) > -1)
                        {
                            cP[pointI] = true;
                        }
                    }

                    if (cP[0] && cP[1] && cP[2] && cP[3])
                    {
                        sFaceIndex = sfIter.key();
                        foundMatch[indexI] = true;
                        break;
                    }
                }

                if (foundMatch[indexI])
                {
                    if (bisectingSlave)
                    {
                        patchCoupling_[pIndex].patchMap().mapMaster
                        (
                            coupleMap::FACE,
                            checkFaces[indexI],
                            sFaceIndex
                        );

                        patchCoupling_[pIndex].patchMap().mapSlave
                        (
                            coupleMap::FACE,
                            sFaceIndex,
                            checkFaces[indexI]
                        );
                    }
                    else
                    {
                        patchCoupling_[pIndex].patchMap().mapSlave
                        (
                            coupleMap::FACE,
                            checkFaces[indexI],
                            sFaceIndex
                        );

                        patchCoupling_[pIndex].patchMap().mapMaster
                        (
                            coupleMap::FACE,
                            sFaceIndex,
                            checkFaces[indexI]
                        );
                    }
                }
                else
                {
                    FatalErrorIn("dynamicTopoFvMesh::bisectQuadFace")
                        << "Failed to build coupled maps."
                        << abort(FatalError);
                }
            }
        }
        else
        if (processorCoupledFace(fIndex))
        {
            // Look for matching slave faces on the patchSubMesh.

        }
    }

    // Write out VTK files after change
    if (debug > 3)
    {
        labelList cellHull(4, -1);

        cellHull[0] = owner_[fIndex];
        cellHull[1] = neighbour_[fIndex];
        cellHull[2] = owner_[newFaceIndex[3]];
        cellHull[3] = neighbour_[newFaceIndex[3]];

        writeVTK
        (
            Foam::name(fIndex)
          + "_Bisect_1",
            cellHull
        );
    }

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    nBisections_++;

    // Increment the number of modifications
    nModifications_++;

    // Specify that the operation was successful
    map.type() = 1;

    // Return the changeMap
    return map;
}

// Method for the bisection of an edge in 3D
// - Returns a changeMap with a type specifying:
//     1: Bisection was successful
//    -1: Bisection failed since max number of topo-changes was reached.
//    -2: Bisection failed since resulting quality would be really bad.
// - AddedPoints contain the index of the newly added point.
const changeMap dynamicTopoFvMesh::bisectEdge
(
    const label eIndex,
    bool checkOnly,
    bool forceOp,
    const changeMap& masterMap
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
    changeMap map, slaveMap;
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

    // Sanity check: Is the index legitimate?
    if (eIndex < 0 || eIndex >= nEdges_)
    {
        FatalErrorIn("dynamicTopoFvMesh::bisectEdge()")
            << " Invalid index: " << eIndex
            << abort(FatalError);
    }

    if (coupledModification_)
    {
        // Is this a locally coupled edge?
        if (locallyCoupledEdge(eIndex))
        {
            label sIndex = -1;

            // Determine the slave index.
            forAll(patchCoupling_, patchI)
            {
                if (patchCoupling_(patchI))
                {
                    const label edgeEnum  = coupleMap::EDGE;
                    const coupleMap& cMap = patchCoupling_[patchI].patchMap();

                    if ((sIndex = cMap.findSlaveIndex(edgeEnum, eIndex)) > -1)
                    {
                        // Keep this index for master/slave mapping.
                        pIndex = patchI;

                        break;
                    }

                    // The following bit happens only during the sliver
                    // exudation process, since slave edges are
                    // usually not added to the coupled edge-stack.
                    if ((sIndex = cMap.findMasterIndex(edgeEnum, eIndex)) > -1)
                    {
                        // Keep this index for master/slave mapping.
                        pIndex = patchI;

                        // Notice that we are bisecting a slave edge first.
                        bisectingSlave = true;

                        break;
                    }
                }
            }

            // Temporarily turn off coupledModification.
            unsetCoupledModification();
            setSlaveModification();

            // First check the slave for bisection feasibility.
            slaveMap = bisectEdge(sIndex, true);

            if (slaveMap.type() == 1)
            {
                // Can the master be bisected as well?
                changeMap masterMap = bisectEdge(eIndex, true);

                // Master couldn't perform bisection
                if (masterMap.type() != 1)
                {
                    setCoupledModification();
                    unsetSlaveModification();

                    return masterMap;
                }

                const labelList& ePoints = edgePoints_[eIndex];

                // Fill the masterMap with points that
                // we seek edge-maps for...
                masterMap.addPoint(edges_[eIndex].end());
                masterMap.addPoint(ePoints[0]);
                masterMap.addPoint(ePoints[ePoints.size() - 1]);

                // Bisect the slave edge
                slaveMap = bisectEdge(sIndex, false, false, masterMap);
            }
            else
            {
                // Slave couldn't perform collapse.
                setCoupledModification();
                unsetSlaveModification();

                map.type() = -2;

                return map;
            }

            // Turn it back on.
            setCoupledModification();
            unsetSlaveModification();
        }
        else
        {
            // Bisect edge on the patchSubMesh.

        }
    }

    // Before we bisect this edge, check whether the operation will
    // yield an acceptable cell-quality.
    scalar minQ = 0.0;

    if ((minQ = computeBisectionQuality(eIndex)) < sliverThreshold_)
    {
        // Check if the quality is actually valid before forcing it.
        if (forceOp && (minQ < 0.0))
        {
            FatalErrorIn("dynamicTopoFvMesh::bisectEdge()")
                << " Forcing bisection on edge: " << eIndex
                << " will yield an invalid cell."
                << abort(FatalError);
        }
        else
        if (!forceOp)
        {
            map.type() = -2;
            return map;
        }
    }

    // Are we performing only checks?
    if (checkOnly)
    {
        map.type() = 1;
        return map;
    }

    // Hull variables
    face tmpTriFace(3);
    labelList tmpEdgeFaces(3,-1);
    labelList tmpIntEdgeFaces(4,-1);
    labelList tmpEdgePoints(3,-1);
    labelList tmpIntEdgePoints(4,-1);
    labelList tmpFaceEdges(3,-1);

    // Make a copy of existing entities
    const labelList vertexHull = edgePoints_[eIndex];
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
        Info << nl << nl
             << "Edge: " << eIndex
             << ": " << edges_[eIndex]
             << " is to be bisected. " << endl;

        label epIndex = whichEdgePatch(eIndex);

        Info << "Patch: ";

        if (epIndex == -1)
        {
            Info << "Internal" << endl;
        }
        else
        {
            Info << boundaryMesh()[epIndex].name() << endl;
        }

        // Write out VTK files prior to change
        if (debug > 3)
        {
            writeVTK
            (
                Foam::name(eIndex)
              + "_Bisect_0",
                cellHull
            );
        }
    }

    labelList mP(2, -1);

    // Set mapping for this point
    mP[0] = edges_[eIndex][0];
    mP[1] = edges_[eIndex][1];

    // Add a new point to the end of the list
    label newPointIndex =
    (
        insertPoint
        (
            0.5 * (points_[mP[0]] + points_[mP[1]]),
            0.5 * (oldPoints_[mP[0]] + oldPoints_[mP[1]]),
            mP
        )
    );

    // Add this point to the map.
    map.addPoint(newPointIndex);

    // Add a new edge to the end of the list
    label newEdgeIndex =
    (
        insertEdge
        (
            whichEdgePatch(eIndex),
            edge(newPointIndex,edges_[eIndex][1]),
            labelList(faceHull.size(),-1),
            vertexHull
        )
    );

    // Add this edge to the map. Since this might require master mapping,
    // first check to see if a slave is being bisected.
    if (slaveModification_)
    {
        const Map<label>& pMap = masterMap.addedPointList();

        // Look through the reverse point map
        // to check which point it corresponds to.
        forAll(patchCoupling_, patchI)
        {
            if (!patchCoupling_(patchI))
            {
                continue;
            }

            const coupleMap& cMap = patchCoupling_[patchI].patchMap();

            Map<label>& rPointMap = cMap.reverseEntityMap(coupleMap::POINT);

            forAll(edges_[eIndex], pointI)
            {
                bool found = false;
                label mPoint = -1;

                if (rPointMap.found(edges_[eIndex][pointI]))
                {
                    mPoint = rPointMap[edges_[eIndex][pointI]];

                    // Now check the master map.
                    if (pMap.found(mPoint))
                    {
                        found = true;
                    }
                }

                // Add the map entry for the point.
                if (found)
                {
                    map.addEdge
                    (
                        newEdgeIndex,
                        mPoint
                    );

                    break;
                }
            }
        }
    }
    else
    {
        map.addEdge(newEdgeIndex);
    }

    // Remove the existing edge from the pointEdges list
    // of the modified point, and add it to the new point
    sizeDownList(eIndex, pointEdges_[edges_[eIndex][1]]);
    sizeUpList(eIndex, pointEdges_[newPointIndex]);

    // Modify the existing edge
    edges_[eIndex][1] = newPointIndex;

    // Add this edge to the map.
    // Although this edge isn't technically 'added', it's
    // required for coupled patch mapping.
    map.addEdge(eIndex);

    // Keep track of added entities
    labelList addedCellIndices(cellHull.size(),-1);
    labelList addedFaceIndices(faceHull.size(),-1);
    labelList addedEdgeIndices(faceHull.size(),-1);
    labelList addedIntFaceIndices(faceHull.size(),-1);

    scalar oldPhi = 0.0;

    // Now loop through the hull and bisect individual entities
    forAll(vertexHull, indexI)
    {
        // Modify the existing face
        replaceLabel
        (
            edges_[newEdgeIndex][1],
            newPointIndex,
            faces_[faceHull[indexI]]
        );

        // Set face flux for the modified face.
        oldPhi = iPtr_->getPhi(faceHull[indexI]);

        iPtr_->setPhi(faceHull[indexI], 0.5 * oldPhi);

        // Prepare a mapping master face list (and corresponding weight)
        labelList mF(1, faceHull[indexI]);
        scalarField mW(1, 1.0);

        // Modify edgePoints for the edge
        replaceLabel
        (
            edges_[newEdgeIndex][1],
            newPointIndex,
            edgePoints_[ringEntities[0][indexI]]
        );

        // Obtain circular indices
        label nextI = vertexHull.fcIndex(indexI);
        label prevI = vertexHull.rcIndex(indexI);

        // Check if this is an interior/boundary face
        if (cellHull[indexI] != -1)
        {
            // Create a new cell. Add it for now, but update later.
            cell newCell(4);

            // Fill-in mapping information
            labelList mC(1, cellHull[indexI]);
            scalarField mW(1, 1.0);

            addedCellIndices[indexI] =
            (
                insertCell(newCell, mC, mW, lengthScale_[cellHull[indexI]])
            );

            // Add this cell to the map.
            map.addCell(addedCellIndices[indexI]);

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

                    iPtr_->setFlip(replaceFace);
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
                edges_[eIndex][0],
                edges_[newEdgeIndex][1],
                edgePoints_[edgeHull[indexI]]
            );

            // Add this edge to faceEdges for the new interior face
            faceEdges_[addedIntFaceIndices[indexI]][0] = edgeHull[indexI];

            // Replace face labels
            replaceLabel
            (
                replaceFace,
                addedIntFaceIndices[indexI],
                cells_[cellHull[indexI]]
            );

            // Add to the new cell
            newCell[1] = replaceFace;

            // Check if this is a boundary face
            if (cellHull[prevI] == -1)
            {
                // Configure the boundary face
                tmpTriFace[0] = newPointIndex;
                tmpTriFace[1] = edges_[newEdgeIndex][1];
                tmpTriFace[2] = vertexHull[indexI];

                // Insert the face
                addedFaceIndices[indexI] =
                (
                    insertFace
                    (
                        whichPatch(faceHull[indexI]),
                        tmpTriFace,
                        addedCellIndices[indexI],
                        -1,
                        mF,
                        mW
                    )
                );

                // Add this face to the map.
                map.addFace(addedFaceIndices[indexI]);

                // Set the flux
                iPtr_->setPhi(addedFaceIndices[indexI], 0.5 * oldPhi);

                // Configure edgeFaces
                tmpEdgeFaces[0] = faceHull[indexI];
                tmpEdgeFaces[1] = addedIntFaceIndices[indexI];
                tmpEdgeFaces[2] = addedFaceIndices[indexI];

                // Configure edgePoints
                tmpEdgePoints[0] = edges_[eIndex][0];
                tmpEdgePoints[1] = vertexHull[nextI];
                tmpEdgePoints[2] = edges_[newEdgeIndex][1];

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

                // Add this edge to the map. Since this might require
                // master mapping, first check to see if a slave
                // is being bisected.
                if (slaveModification_)
                {
                    const Map<label>& pMap = masterMap.addedPointList();

                    // Look through the reverse point map
                    // to check which point it corresponds to.
                    forAll(patchCoupling_, patchI)
                    {
                        if (!patchCoupling_(patchI))
                        {
                            continue;
                        }

                        const coupleMap& cMap =
                        (
                            patchCoupling_[patchI].patchMap()
                        );

                        const label pointEnum = coupleMap::POINT;

                        Map<label>& rPointMap =
                        (
                            cMap.reverseEntityMap(pointEnum)
                        );

                        if (rPointMap.found(vertexHull[indexI]))
                        {
                            label mPoint = rPointMap[vertexHull[indexI]];

                            // Now check the master map.
                            if (pMap.found(mPoint))
                            {
                                // Add the map entry for the point.
                                map.addEdge
                                (
                                    addedEdgeIndices[indexI],
                                    mPoint
                                );

                                break;
                            }
                        }
                    }
                }
                else
                {
                    map.addEdge(addedEdgeIndices[indexI]);
                }

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
                    edges_[eIndex][0],
                    newPointIndex,
                    edgePoints_[ringEntities[2][indexI]]
                );

                // Add the faceEdges entry
                faceEdges_.append(tmpFaceEdges);

                // Add an entry to edgeFaces
                edgeFaces_[newEdgeIndex][indexI] = addedFaceIndices[indexI];

                // Add an entry for this cell
                newCell[2] = addedFaceIndices[indexI];
            }
            else
            // Check if a cell was added before this
            if (addedCellIndices[prevI] != -1)
            {
                // Configure the interior face
                tmpTriFace[0] = vertexHull[indexI];
                tmpTriFace[1] = edges_[newEdgeIndex][1];
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

                // Set the flux, but check orientation first.
                scalar fSign = 1.0;

                if (cellHull[indexI] < cellHull[prevI])
                {
                    fSign = -1.0;
                }

                iPtr_->setPhi(addedFaceIndices[indexI], (0.5*fSign*oldPhi));

                // Configure edgeFaces
                tmpIntEdgeFaces[0] = faceHull[indexI];
                tmpIntEdgeFaces[1] = addedIntFaceIndices[indexI];
                tmpIntEdgeFaces[2] = addedFaceIndices[indexI];
                tmpIntEdgeFaces[3] = addedIntFaceIndices[prevI];

                // Configure edgePoints
                tmpIntEdgePoints[0] = edges_[eIndex][0];
                tmpIntEdgePoints[1] = vertexHull[nextI];
                tmpIntEdgePoints[2] = edges_[newEdgeIndex][1];
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
                    edges_[eIndex][0],
                    newPointIndex,
                    edgePoints_[ringEntities[2][indexI]]
                );

                // Add the faceEdges entry
                faceEdges_.append(tmpFaceEdges);

                // Add an entry to edgeFaces
                edgeFaces_[newEdgeIndex][indexI] = addedFaceIndices[indexI];

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
                tmpTriFace[1] = edges_[newEdgeIndex][1];
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

                // Set the flux, but check orientation first.
                scalar fSign = 1.0;

                if (cellHull[indexI] < cellHull[0])
                {
                    fSign = -1.0;
                }

                iPtr_->setPhi(addedFaceIndices[0], (0.5*fSign*oldPhi));

                // Configure edgeFaces
                tmpIntEdgeFaces[0] = faceHull[0];
                tmpIntEdgeFaces[1] = addedIntFaceIndices[0];
                tmpIntEdgeFaces[2] = addedFaceIndices[0];
                tmpIntEdgeFaces[3] = addedIntFaceIndices[indexI];

                // Configure edgePoints
                tmpIntEdgePoints[0] = edges_[eIndex][0];
                tmpIntEdgePoints[1] = vertexHull[1];
                tmpIntEdgePoints[2] = edges_[newEdgeIndex][1];
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
                    edges_[eIndex][0],
                    newPointIndex,
                    edgePoints_[ringEntities[2][0]]
                );

                // Add the faceEdges entry
                faceEdges_.append(tmpFaceEdges);

                // Add an entry to edgeFaces
                edgeFaces_[newEdgeIndex][0] = addedFaceIndices[0];

                // Add an entry for this cell
                newCell[3] = addedFaceIndices[0];

                // Make the final entry for the first cell
                cells_[addedCellIndices[0]][2] = addedFaceIndices[0];
            }

            // Update the cell list with the new cell.
            cells_[addedCellIndices[indexI]] = newCell;
        }
        else
        {
            // Configure the final boundary face
            tmpTriFace[0] = vertexHull[indexI];
            tmpTriFace[1] = edges_[newEdgeIndex][1];
            tmpTriFace[2] = newPointIndex;

            // Insert the face
            addedFaceIndices[indexI] =
            (
                insertFace
                (
                    whichPatch(faceHull[indexI]),
                    tmpTriFace,
                    addedCellIndices[prevI],
                    -1,
                    mF,
                    mW
                )
            );

            // Add this face to the map.
            map.addFace(addedFaceIndices[indexI]);

            // Set the flux
            iPtr_->setPhi(addedFaceIndices[indexI], 0.5 * oldPhi);

            // Configure edgeFaces
            tmpEdgeFaces[0] = addedFaceIndices[indexI];
            tmpEdgeFaces[1] = addedIntFaceIndices[prevI];
            tmpEdgeFaces[2] = faceHull[indexI];

            // Configure edgePoints
            tmpEdgePoints[0] = edges_[newEdgeIndex][1];
            tmpEdgePoints[1] = vertexHull[prevI];
            tmpEdgePoints[2] = edges_[eIndex][0];

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

            // Add this edge to the map. Since this might require
            // master mapping, first check to see if a slave
            // is being bisected.
            if (slaveModification_)
            {
                const Map<label>& pMap = masterMap.addedPointList();

                // Look through the reverse point map
                // to check which point it corresponds to.
                forAll(patchCoupling_, patchI)
                {
                    if (!patchCoupling_(patchI))
                    {
                        continue;
                    }

                    const label pEnum = coupleMap::POINT;
                    const coupleMap& cMap = patchCoupling_[patchI].patchMap();

                    Map<label>& rPointMap = cMap.reverseEntityMap(pEnum);

                    if (rPointMap.found(vertexHull[indexI]))
                    {
                        label mPoint = rPointMap[vertexHull[indexI]];

                        // Now check the master map.
                        if (pMap.found(mPoint))
                        {
                            // Add the map entry for the point.
                            map.addEdge
                            (
                                addedEdgeIndices[indexI],
                                mPoint
                            );

                            break;
                        }
                    }
                }
            }
            else
            {
                map.addEdge(addedEdgeIndices[indexI]);
            }

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
                edges_[eIndex][0],
                newPointIndex,
                edgePoints_[ringEntities[2][indexI]]
            );

            // Add the faceEdges entry
            faceEdges_.append(tmpFaceEdges);

            // Add an entry to edgeFaces
            edgeFaces_[newEdgeIndex][indexI] = addedFaceIndices[indexI];

            // Make the final entry for the previous cell
            cells_[addedCellIndices[prevI]][3] = addedFaceIndices[indexI];
        }
    }

    // Now that all old / new cells possess correct connectivity,
    // compute values for old cell volume.
    forAll(cellHull, indexI)
    {
        if (cellHull[indexI] == -1)
        {
            continue;
        }

        // Compute old volumes, using old point positions.
        scalar modOldVol = tetVolume(cellHull[indexI], true);
        scalar newOldVol = tetVolume(addedCellIndices[indexI], true);

        if (modOldVol < 0.0 || newOldVol < 0.0)
        {
            FatalErrorIn
            (
                "dynamicTopoFvMesh::bisectEdge()"
            )
                << "Negative old-volumes encountered." << nl
                << cellHull[indexI] << ": " << modOldVol
                << addedCellIndices[indexI] << ": " << newOldVol
                << abort(FatalError);
        }

        if (debug > 2)
        {
            Info << "Cell: " << cellHull[indexI]
                 << " Old volume: " << modOldVol
                 << " New volume: " << tetVolume(cellHull[indexI])
                 << endl;

            Info << "Cell: " << addedCellIndices[indexI]
                 << " Old volume: " << newOldVol
                 << " New volume: " << tetVolume(addedCellIndices[indexI])
                 << endl;
        }
    }

    if (coupledModification_)
    {
        // Create a master/slave entry for the new edges on the patch.
        if (locallyCoupledEdge(eIndex))
        {
            FixedList<bool, 3> foundMatch(false), reqCheck(false);
            FixedList<label, 3> checkPoints(-1), checkEdges(-1);

            // Fill in the edges to be check for...
            label eCounter = 0;

            // Add the check-point.
            // Note that new edges are of the form:
            //   [newPointIndex, existingPoint]
            checkEdges[eCounter] = newEdgeIndex;
            checkPoints[eCounter] = edges_[newEdgeIndex].end();

            if (patchCoupling_(pIndex))
            {
                // Add the new point to the coupling map
                const coupleMap& cMap = patchCoupling_[pIndex].patchMap();

                // Update pointMap
                Map<label>& pointMap = cMap.entityMap(coupleMap::POINT);

                pointMap.insert
                (
                    newPointIndex,
                    slaveMap.addedPointList().begin().key()
                );

                // Update reverse pointMap
                Map<label>& rPointMap = cMap.reverseEntityMap(coupleMap::POINT);

                rPointMap.insert
                (
                    slaveMap.addedPointList().begin().key(),
                    newPointIndex
                );

                reqCheck[eCounter++] = true;
            }

            // ... and two new boundary edges.
            forAll(addedEdgeIndices, edgeI)
            {
                label chIndex = whichEdgePatch(addedEdgeIndices[edgeI]);

                if ((chIndex != -1) && patchCoupling_(chIndex))
                {
                    // Add the check-point
                    checkEdges[eCounter] = addedEdgeIndices[edgeI];
                    checkPoints[eCounter] = edges_[checkEdges[eCounter]].end();

                    reqCheck[eCounter++] = true;
                }
            }

            // Mark off un-necessary searches
            forAll(reqCheck, indexI)
            {
                if (!reqCheck[indexI])
                {
                    foundMatch[indexI] = true;
                }
            }

            const Map<label>& aeList = slaveMap.addedEdgeList();

            // Compare with all check entries.
            forAll(reqCheck, indexI)
            {
                if (!reqCheck[indexI])
                {
                    continue;
                }

                label mapEdgeIndex = -1;

                forAllConstIter(Map<label>, aeList, eIter)
                {
                    if (eIter() == checkPoints[indexI])
                    {
                        mapEdgeIndex = eIter.key();
                        break;
                    }
                }

                if (bisectingSlave)
                {
                    patchCoupling_[pIndex].patchMap().mapMaster
                    (
                        coupleMap::EDGE,
                        checkEdges[indexI],
                        mapEdgeIndex
                    );

                    patchCoupling_[pIndex].patchMap().mapSlave
                    (
                        coupleMap::EDGE,
                        mapEdgeIndex,
                        checkEdges[indexI]
                    );

                    if (debug > 2)
                    {
                        Info << "Mapping: " << endl;
                        Info << " Slave: " << checkEdges[indexI]
                             << " Master: " << mapEdgeIndex
                             << endl;
                    }
                }
                else
                {
                    patchCoupling_[pIndex].patchMap().mapSlave
                    (
                        coupleMap::EDGE,
                        checkEdges[indexI],
                        mapEdgeIndex
                    );

                    patchCoupling_[pIndex].patchMap().mapMaster
                    (
                        coupleMap::EDGE,
                        mapEdgeIndex,
                        checkEdges[indexI]
                    );

                    if (debug > 2)
                    {
                        Info << "Mapping: " << endl;
                        Info << " Master: " << checkEdges[indexI]
                             << " Slave: " << mapEdgeIndex
                             << endl;
                    }
                }

                foundMatch[indexI] = true;
            }

            if (!(foundMatch[0] && foundMatch[1] && foundMatch[2]))
            {
                forAll(checkEdges, edgeI)
                {
                    Info << "checkEdges: " << endl;
                    Info << checkEdges[edgeI] << ": "
                         << edges_[checkEdges[edgeI]]
                         << endl;
                }

                FatalErrorIn("dynamicTopoFvMesh::bisectEdge")
                    << "Failed to build coupled maps." << nl
                    << " foundMatch: " << foundMatch << nl
                    << " eCounter: " << eCounter << nl
                    << abort(FatalError);
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
        label bPatch = whichEdgePatch(eIndex);

        if (bPatch == -1)
        {
            Info << "Patch: Internal" << endl;
        }
        else
        {
            Info << "Patch: " << boundaryMesh()[bPatch].name() << endl;
        }

        Info << "EdgePoints: " << vertexHull << endl;
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
                 << cells_[addedCellIndices[cellI]] << nl
                 << "lengthScale: " << lengthScale_[addedCellIndices[cellI]]
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
              + "_Bisect_1",
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

    // Specify that the operation was successful
    map.type() = 1;

    // Return the changeMap
    return map;
}

// Method for the trisection of a face in 3D
// - Returns a changeMap with a type specifying:
//     1: Trisection was successful
//    -1: Trisection failed since max number of topo-changes was reached.
//    -2: Trisection failed since resulting quality would be really bad.
// - AddedPoint is the index of the newly added point.
const changeMap dynamicTopoFvMesh::trisectFace
(
    const label fIndex,
    bool checkOnly,
    bool forceOp,
    const changeMap& masterMap
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
    changeMap map, slaveMap;
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

    // Sanity check: Is the index legitimate?
    if (fIndex < 0 || fIndex >= nFaces_)
    {
        FatalErrorIn("dynamicTopoFvMesh::trisectFace()")
            << " Invalid index: " << fIndex
            << abort(FatalError);
    }

    if (coupledModification_)
    {
        // Are all edges of this face locally coupled?
        bool locallyCoupledFace = false;

        const labelList& fEdges = faceEdges_[fIndex];

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
            label slaveIndex = -1;

            // Loop through master/slave maps
            // and determine the coupled edge index.
            forAll(patchCoupling_, patchI)
            {
                if (!patchCoupling_(patchI))
                {
                    continue;
                }

                const label edgeEnum  = coupleMap::EDGE;
                const coupleMap& cMap = patchCoupling_[patchI].patchMap();

                if
                (
                    (cMap.findSlaveIndex(edgeEnum, fEdges[0]) > -1) &&
                    (cMap.findSlaveIndex(edgeEnum, fEdges[1]) > -1) &&
                    (cMap.findSlaveIndex(edgeEnum, fEdges[2]) > -1)
                )
                {
                    bool foundFace = false;

                    // Search for a boundary face shared
                    // by two of the edges.
                    const labelList& eFirstFaces = edgeFaces_[fEdges[0]];
                    const labelList& eSecondFaces = edgeFaces_[fEdges[1]];

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
                                    slaveIndex = eFirstFaces[edgeI];
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

                    // Keep this index for master/slave mapping.
                    pIndex = patchI;
                }

                // The following bit happens only during the sliver
                // exudation process.
                if
                (
                    (cMap.findMasterIndex(edgeEnum, fEdges[0]) > -1) &&
                    (cMap.findMasterIndex(edgeEnum, fEdges[1]) > -1) &&
                    (cMap.findMasterIndex(edgeEnum, fEdges[2]) > -1)
                )
                {
                    bool foundFace = false;

                    // Search for a boundary face shared
                    // by two of the edges.
                    const labelList& eFirstFaces = edgeFaces_[fEdges[0]];
                    const labelList& eSecondFaces = edgeFaces_[fEdges[1]];

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
                                    slaveIndex = eFirstFaces[edgeI];
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

                    // Keep this index for master/slave mapping.
                    pIndex = patchI;

                    // Notice that we are trisecting a slave edge.
                    bisectingSlave = true;
                }
            }

            // Temporarily turn off coupledModification.
            unsetCoupledModification();
            setSlaveModification();

            // First check the slave for trisection feasibility.
            slaveMap = trisectFace(slaveIndex, true);

            if (slaveMap.type() == 1)
            {
                // Can the master be trisected as well?
                changeMap masterMap = trisectFace(fIndex, true);

                // Master couldn't perform bisection
                if (masterMap.type() != 1)
                {
                    setCoupledModification();
                    unsetSlaveModification();

                    return masterMap;
                }

                // Fill the masterMap with points that
                // we seek edge-maps for...
                const face& thisFace = faces_[fIndex];

                masterMap.addPoint(thisFace[0]);
                masterMap.addPoint(thisFace[1]);
                masterMap.addPoint(thisFace[2]);

                // Trisect the slave face
                slaveMap = trisectFace(slaveIndex, false, false, masterMap);
            }
            else
            {
                // Slave couldn't perform collapse.
                setCoupledModification();
                unsetSlaveModification();

                map.type() = -2;

                return map;
            }

            // Turn it back on.
            setCoupledModification();
            unsetSlaveModification();
        }
    }

    // Before we trisect this face, check whether the operation will
    // yield an acceptable cell-quality.
    scalar minQ = 0.0;

    if ((minQ = computeTrisectionQuality(fIndex)) < sliverThreshold_)
    {
        // Check if the quality is actually valid before forcing it.
        if (forceOp && (minQ < 0.0))
        {
            FatalErrorIn("dynamicTopoFvMesh::trisectFace()")
                << " Forcing trisection on face: " << fIndex
                << " will yield an invalid cell."
                << abort(FatalError);
        }
        else
        if (!forceOp)
        {
            map.type() = -2;
            return map;
        }
    }

    // Are we performing only checks?
    if (checkOnly)
    {
        map.type() = 1;
        return map;
    }

    // Hull variables
    face tmpTriFace(3);
    labelList newTriEdgeFaces(3), newTriEdgePoints(3);
    labelList newQuadEdgeFaces(4), newQuadEdgePoints(4);

    FixedList<label,2> apexPoint(-1);
    FixedList<face, 3> checkFace(face(3));
    FixedList<label,5> newEdgeIndex(-1);
    FixedList<label,9> newFaceIndex(-1);
    FixedList<label,6> newCellIndex(-1);
    FixedList<cell, 6> newTetCell(cell(4));
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
        Info << nl << nl
             << "Face: " << fIndex
             << ": " << faces_[fIndex]
             << " is to be trisected. " << endl;

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

    labelList mP(3, -1);

    // Fill in mapping information
    mP[0] = faces_[fIndex][0];
    mP[1] = faces_[fIndex][1];
    mP[2] = faces_[fIndex][2];

    // Add a new point to the end of the list
    scalar oT = (1.0/3.0);

    label newPointIndex =
    (
        insertPoint
        (
            oT * (points_[mP[0]] + points_[mP[1]] + points_[mP[2]]),
            oT * (oldPoints_[mP[0]] + oldPoints_[mP[1]] + oldPoints_[mP[2]]),
            mP
        )
    );

    // Add this point to the map.
    map.addPoint(newPointIndex);

    // Add three new cells to the end of the cell list
    for (label i = 0; i < 3; i++)
    {
        scalar parentScale = -1.0;

        if (edgeRefinement_)
        {
            parentScale = lengthScale_[cellsForRemoval[0]];
        }

        // Fill-in mapping information
        labelList mC(1, cellsForRemoval[0]);
        scalarField mW(1, 1.0);

        newCellIndex[i] = insertCell(newTetCell[i], mC, mW, parentScale);

        // Add cells to the map
        map.addCell(newCellIndex[i]);
    }

    // Find the apex point for this cell
    apexPoint[0] = tetApexPoint(owner_[fIndex], fIndex);

    // Insert three new internal faces

    // First face: Owner: newCellIndex[0], Neighbour: newCellIndex[1]
    tmpTriFace[0] = newPointIndex;
    tmpTriFace[1] = faces_[fIndex][0];
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
    tmpTriFace[1] = faces_[fIndex][1];
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
    tmpTriFace[2] = faces_[fIndex][2];

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
    newTriEdgePoints[0] = faces_[fIndex][0];
    newTriEdgePoints[1] = faces_[fIndex][1];
    newTriEdgePoints[2] = faces_[fIndex][2];

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
    newTetCell[0][nF[0]++] = newFaceIndex[0];
    newTetCell[0][nF[0]++] = newFaceIndex[2];
    newTetCell[1][nF[1]++] = newFaceIndex[0];
    newTetCell[1][nF[1]++] = newFaceIndex[1];
    newTetCell[2][nF[2]++] = newFaceIndex[1];
    newTetCell[2][nF[2]++] = newFaceIndex[2];

    // Define the three faces to check for orientation:
    checkFace[0][0] = faces_[fIndex][2];
    checkFace[0][1] = apexPoint[0];
    checkFace[0][2] = faces_[fIndex][0];

    checkFace[1][0] = faces_[fIndex][0];
    checkFace[1][1] = apexPoint[0];
    checkFace[1][2] = faces_[fIndex][1];

    checkFace[2][0] = faces_[fIndex][1];
    checkFace[2][1] = apexPoint[0];
    checkFace[2][2] = faces_[fIndex][2];

    // Check the orientation of faces on the first cell.
    forAll(cells_[owner_[fIndex]], faceI)
    {
        label faceIndex = cells_[owner_[fIndex]][faceI];

        if (faceIndex == fIndex)
        {
            continue;
        }

        const face& faceToCheck = faces_[faceIndex];
        label cellIndex = cellsForRemoval[0];
        label newIndex = -1;

        // Check against faces.
        if (compare(faceToCheck, checkFace[0]) != 0)
        {
            newIndex = newCellIndex[0];
            newTetCell[0][nF[0]++] = faceIndex;
        }
        else
        if (compare(faceToCheck, checkFace[1]) != 0)
        {
            newIndex = newCellIndex[1];
            newTetCell[1][nF[1]++] = faceIndex;
        }
        else
        if (compare(faceToCheck, checkFace[2]) != 0)
        {
            newIndex = newCellIndex[2];
            newTetCell[2][nF[2]++] = faceIndex;
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

                iPtr_->setFlip(faceIndex);
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
        tmpTriFace[1] = faces_[fIndex][2];
        tmpTriFace[2] = faces_[fIndex][0];

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
        tmpTriFace[1] = faces_[fIndex][0];
        tmpTriFace[2] = faces_[fIndex][1];

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
        tmpTriFace[1] = faces_[fIndex][1];
        tmpTriFace[2] = faces_[fIndex][2];

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
        newTetCell[0][nF[0]++] = newFaceIndex[3];
        newTetCell[1][nF[1]++] = newFaceIndex[4];
        newTetCell[2][nF[2]++] = newFaceIndex[5];

        // Configure edgeFaces and edgePoints for three new boundary edges.
        newTriEdgeFaces[0] = newFaceIndex[4];
        newTriEdgeFaces[1] = newFaceIndex[0];
        newTriEdgeFaces[2] = newFaceIndex[3];

        newTriEdgePoints[0] = faces_[fIndex][1];
        newTriEdgePoints[1] = apexPoint[0];
        newTriEdgePoints[2] = faces_[fIndex][2];

        newEdgeIndex[1] =
        (
            insertEdge
            (
                facePatch,
                edge
                (
                   newPointIndex,
                   faces_[fIndex][0]
                ),
                newTriEdgeFaces,
                newTriEdgePoints
            )
        );

        newTriEdgeFaces[0] = newFaceIndex[5];
        newTriEdgeFaces[1] = newFaceIndex[1];
        newTriEdgeFaces[2] = newFaceIndex[4];

        newTriEdgePoints[0] = faces_[fIndex][2];
        newTriEdgePoints[1] = apexPoint[0];
        newTriEdgePoints[2] = faces_[fIndex][0];

        newEdgeIndex[2] =
        (
            insertEdge
            (
                facePatch,
                edge
                (
                   newPointIndex,
                   faces_[fIndex][1]
                ),
                newTriEdgeFaces,
                newTriEdgePoints
            )
        );

        newTriEdgeFaces[0] = newFaceIndex[3];
        newTriEdgeFaces[1] = newFaceIndex[2];
        newTriEdgeFaces[2] = newFaceIndex[5];

        newTriEdgePoints[0] = faces_[fIndex][0];
        newTriEdgePoints[1] = apexPoint[0];
        newTriEdgePoints[2] = faces_[fIndex][1];

        newEdgeIndex[3] =
        (
            insertEdge
            (
                facePatch,
                edge
                (
                   newPointIndex,
                   faces_[fIndex][2]
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

        check[0][0] = apexPoint[0]; check[0][1] = faces_[fIndex][0];
        check[1][0] = apexPoint[0]; check[1][1] = faces_[fIndex][1];
        check[2][0] = apexPoint[0]; check[2][1] = faces_[fIndex][2];

        check[3][0] = faces_[fIndex][2]; check[3][1] = faces_[fIndex][0];
        check[4][0] = faces_[fIndex][0]; check[4][1] = faces_[fIndex][1];
        check[5][0] = faces_[fIndex][1]; check[5][1] = faces_[fIndex][2];

        // Build a list of cellEdges
        labelHashSet cellEdges;

        forAll(cells_[owner_[fIndex]], faceI)
        {
            const labelList& fEdges =
            (
                faceEdges_[cells_[owner_[fIndex]][faceI]]
            );

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
            const edge& edgeToCheck = edges_[eIter.key()];

            // Check against the specified edges.
            if (edgeToCheck == check[0])
            {
                insertLabel
                (
                    newPointIndex,
                    faces_[fIndex][1],
                    faces_[fIndex][2],
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
                    faces_[fIndex][0],
                    faces_[fIndex][2],
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
                    faces_[fIndex][0],
                    faces_[fIndex][1],
                    edgePoints_[eIter.key()]
                );

                sizeUpList(newFaceIndex[2], edgeFaces_[eIter.key()]);
                newFaceEdges[2][nE[2]++] = eIter.key();
            }

            if (edgeToCheck == check[3])
            {
                replaceLabel
                (
                    faces_[fIndex][1],
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
                    faces_[fIndex][2],
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
                    faces_[fIndex][0],
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

            // Add faces to the map.
            map.addFace(newFaceIndex[i]);
        }

        // If modification is coupled, generate mapping info.
        if (coupledModification_)
        {
            // Create a master/slave entry for the new edges on the patch.
            if (locallyCoupledFace(fIndex))
            {
                FixedList<bool, 3> foundMatch(false);
                FixedList<label, 3> checkPoints(-1), checkEdges(-1);

                // Fill in the edges to be check for...
                label eCounter = 0;

                for (label i = 1; i <= 3; i++)
                {
                    checkEdges[eCounter] = newEdgeIndex[i];
                    checkEdges[eCounter] = faces_[fIndex][i-1];

                    eCounter++;
                }

                // Add the new point to the coupling map
                const coupleMap& cMap = patchCoupling_[pIndex].patchMap();

                // Update pointMap
                Map<label>& pointMap = cMap.entityMap(coupleMap::POINT);

                pointMap.insert
                (
                    newPointIndex,
                    slaveMap.addedPointList().begin().key()
                );

                // Update reverse pointMap
                Map<label>& rPointMap = cMap.reverseEntityMap(coupleMap::POINT);

                rPointMap.insert
                (
                    slaveMap.addedPointList().begin().key(),
                    newPointIndex
                );

                const Map<label>& aeList = slaveMap.addedEdgeList();

                // Compare with all three entries.
                forAll(checkPoints, indexI)
                {
                    label mapEdgeIndex = -1;

                    forAllConstIter(Map<label>, aeList, eIter)
                    {
                        if (eIter() == checkPoints[indexI])
                        {
                            mapEdgeIndex = eIter.key();
                            break;
                        }
                    }

                    if (bisectingSlave)
                    {
                        cMap.mapMaster
                        (
                            coupleMap::EDGE,
                            checkEdges[indexI],
                            mapEdgeIndex
                        );

                        cMap.mapSlave
                        (
                            coupleMap::EDGE,
                            mapEdgeIndex,
                            checkEdges[indexI]
                        );
                    }
                    else
                    {
                        cMap.mapSlave
                        (
                            coupleMap::EDGE,
                            checkEdges[indexI],
                            mapEdgeIndex
                        );

                        cMap.mapMaster
                        (
                            coupleMap::EDGE,
                            mapEdgeIndex,
                            checkEdges[indexI]
                        );
                    }

                    foundMatch[indexI] = true;
                }

                if (!(foundMatch[0] && foundMatch[1] && foundMatch[2]))
                {
                    forAll(checkEdges, edgeI)
                    {
                        Info << "checkEdges: " << endl;
                        Info << checkEdges[edgeI] << ": "
                             << edges_[checkEdges[edgeI]]
                             << endl;
                    }

                    FatalErrorIn("dynamicTopoFvMesh::trisectFace()")
                        << "Failed to build coupled maps."
                        << abort(FatalError);
                }
            }
            else
            if (processorCoupledFace(fIndex))
            {
                // Look for matching slave edges on the patchSubMesh.

            }
        }
    }
    else
    {
        // Add three new cells to the end of the cell list
        for (label i = 3; i < 6; i++)
        {
            scalar parentScale = -1.0;

            if (edgeRefinement_)
            {
                parentScale = lengthScale_[cellsForRemoval[1]];
            }

            // Fill-in mapping information
            labelList mC(1, cellsForRemoval[1]);
            scalarField mW(1, 1.0);

            newCellIndex[i] = insertCell(newTetCell[i], mC, mW, parentScale);

            // Add to the map.
            map.addCell(newCellIndex[i]);
        }

        // Find the apex point for this cell
        apexPoint[1] = tetApexPoint(neighbour_[fIndex], fIndex);

        // Add six new interior faces.

        // Fourth face: Owner: newCellIndex[0], Neighbour: newCellIndex[3]
        tmpTriFace[0] = newPointIndex;
        tmpTriFace[1] = faces_[fIndex][2];
        tmpTriFace[2] = faces_[fIndex][0];

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
        tmpTriFace[1] = faces_[fIndex][0];
        tmpTriFace[2] = faces_[fIndex][1];

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
        tmpTriFace[1] = faces_[fIndex][1];
        tmpTriFace[2] = faces_[fIndex][2];

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
        tmpTriFace[2] = faces_[fIndex][0];

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
        tmpTriFace[2] = faces_[fIndex][1];

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
        tmpTriFace[1] = faces_[fIndex][2];
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
        newTetCell[3][nF[3]++] = newFaceIndex[6];
        newTetCell[3][nF[3]++] = newFaceIndex[8];
        newTetCell[4][nF[4]++] = newFaceIndex[6];
        newTetCell[4][nF[4]++] = newFaceIndex[7];
        newTetCell[5][nF[5]++] = newFaceIndex[7];
        newTetCell[5][nF[5]++] = newFaceIndex[8];

        newTetCell[0][nF[0]++] = newFaceIndex[3];
        newTetCell[1][nF[1]++] = newFaceIndex[4];
        newTetCell[2][nF[2]++] = newFaceIndex[5];

        newTetCell[3][nF[3]++] = newFaceIndex[3];
        newTetCell[4][nF[4]++] = newFaceIndex[4];
        newTetCell[5][nF[5]++] = newFaceIndex[5];

        // Define the three faces to check for orientation:
        checkFace[0][0] = faces_[fIndex][2];
        checkFace[0][1] = apexPoint[1];
        checkFace[0][2] = faces_[fIndex][0];

        checkFace[1][0] = faces_[fIndex][0];
        checkFace[1][1] = apexPoint[1];
        checkFace[1][2] = faces_[fIndex][1];

        checkFace[2][0] = faces_[fIndex][1];
        checkFace[2][1] = apexPoint[1];
        checkFace[2][2] = faces_[fIndex][2];

        // Check the orientation of faces on the second cell.
        forAll(cells_[neighbour_[fIndex]], faceI)
        {
            label faceIndex = cells_[neighbour_[fIndex]][faceI];

            if (faceIndex == fIndex)
            {
                continue;
            }

            const face& faceToCheck = faces_[faceIndex];
            label cellIndex = cellsForRemoval[1];
            label newIndex = -1;

            // Check against faces.
            if (compare(faceToCheck, checkFace[0]) != 0)
            {
                newIndex = newCellIndex[3];
                newTetCell[3][nF[3]++] = faceIndex;
            }
            else
            if (compare(faceToCheck, checkFace[1]) != 0)
            {
                newIndex = newCellIndex[4];
                newTetCell[4][nF[4]++] = faceIndex;
            }
            else
            if (compare(faceToCheck, checkFace[2]) != 0)
            {
                newIndex = newCellIndex[5];
                newTetCell[5][nF[5]++] = faceIndex;
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

                    iPtr_->setFlip(faceIndex);
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

        newQuadEdgePoints[0] = faces_[fIndex][1];
        newQuadEdgePoints[1] = apexPoint[0];
        newQuadEdgePoints[2] = faces_[fIndex][2];
        newQuadEdgePoints[3] = apexPoint[1];

        newEdgeIndex[1] =
        (
            insertEdge
            (
                -1,
                edge
                (
                   newPointIndex,
                   faces_[fIndex][0]
                ),
                newQuadEdgeFaces,
                newQuadEdgePoints
            )
        );

        newQuadEdgeFaces[0] = newFaceIndex[5];
        newQuadEdgeFaces[1] = newFaceIndex[1];
        newQuadEdgeFaces[2] = newFaceIndex[4];
        newQuadEdgeFaces[3] = newFaceIndex[7];

        newQuadEdgePoints[0] = faces_[fIndex][2];
        newQuadEdgePoints[1] = apexPoint[0];
        newQuadEdgePoints[2] = faces_[fIndex][0];
        newQuadEdgePoints[3] = apexPoint[1];

        newEdgeIndex[2] =
        (
            insertEdge
            (
                -1,
                edge
                (
                   newPointIndex,
                   faces_[fIndex][1]
                ),
                newQuadEdgeFaces,
                newQuadEdgePoints
            )
        );

        newQuadEdgeFaces[0] = newFaceIndex[3];
        newQuadEdgeFaces[1] = newFaceIndex[2];
        newQuadEdgeFaces[2] = newFaceIndex[5];
        newQuadEdgeFaces[3] = newFaceIndex[8];

        newQuadEdgePoints[0] = faces_[fIndex][0];
        newQuadEdgePoints[1] = apexPoint[0];
        newQuadEdgePoints[2] = faces_[fIndex][1];
        newQuadEdgePoints[3] = apexPoint[1];

        newEdgeIndex[3] =
        (
            insertEdge
            (
                -1,
                edge
                (
                   newPointIndex,
                   faces_[fIndex][2]
                ),
                newQuadEdgeFaces,
                newQuadEdgePoints
            )
        );

        newTriEdgeFaces[0] = newFaceIndex[6];
        newTriEdgeFaces[1] = newFaceIndex[7];
        newTriEdgeFaces[2] = newFaceIndex[8];

        newTriEdgePoints[0] = faces_[fIndex][0];
        newTriEdgePoints[1] = faces_[fIndex][1];
        newTriEdgePoints[2] = faces_[fIndex][2];

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

        check[0][0] = apexPoint[0]; check[0][1] = faces_[fIndex][0];
        check[1][0] = apexPoint[0]; check[1][1] = faces_[fIndex][1];
        check[2][0] = apexPoint[0]; check[2][1] = faces_[fIndex][2];

        check[3][0] = faces_[fIndex][2]; check[3][1] = faces_[fIndex][0];
        check[4][0] = faces_[fIndex][0]; check[4][1] = faces_[fIndex][1];
        check[5][0] = faces_[fIndex][1]; check[5][1] = faces_[fIndex][2];

        check[6][0] = apexPoint[1]; check[6][1] = faces_[fIndex][0];
        check[7][0] = apexPoint[1]; check[7][1] = faces_[fIndex][1];
        check[8][0] = apexPoint[1]; check[8][1] = faces_[fIndex][2];

        // Build a list of cellEdges
        labelHashSet cellEdges;

        forAll(cells_[owner_[fIndex]], faceI)
        {
            const labelList& fEdges =
            (
                faceEdges_[cells_[owner_[fIndex]][faceI]]
            );

            forAll(fEdges, edgeI)
            {
                if (!cellEdges.found(fEdges[edgeI]))
                {
                    cellEdges.insert(fEdges[edgeI]);
                }
            }
        }

        forAll(cells_[neighbour_[fIndex]], faceI)
        {
            const labelList& fEdges =
            (
                faceEdges_[cells_[neighbour_[fIndex]][faceI]]
            );

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
            const edge& edgeToCheck = edges_[eIter.key()];

            // Check against the specified edges.
            if (edgeToCheck == check[0])
            {
                insertLabel
                (
                    newPointIndex,
                    faces_[fIndex][1],
                    faces_[fIndex][2],
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
                    faces_[fIndex][0],
                    faces_[fIndex][2],
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
                    faces_[fIndex][0],
                    faces_[fIndex][1],
                    edgePoints_[eIter.key()]
                );

                sizeUpList(newFaceIndex[2], edgeFaces_[eIter.key()]);
                newFaceEdges[2][nE[2]++] = eIter.key();
            }

            if (edgeToCheck == check[3])
            {
                replaceLabel
                (
                    faces_[fIndex][1],
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
                    faces_[fIndex][2],
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
                    faces_[fIndex][0],
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
                    faces_[fIndex][1],
                    faces_[fIndex][2],
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
                    faces_[fIndex][0],
                    faces_[fIndex][2],
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
                    faces_[fIndex][0],
                    faces_[fIndex][1],
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

            // Add faces to the map.
            map.addFace(newFaceIndex[i]);
        }
    }

    // Added edges are those connected to the new point
    const labelList& pointEdges = pointEdges_[newPointIndex];

    forAll(pointEdges, edgeI)
    {
        if (slaveModification_)
        {
            const Map<label>& pMap = masterMap.addedPointList();

            // Look through the reverse point map
            // to check which point it corresponds to.
            forAll(patchCoupling_, patchI)
            {
                if (!patchCoupling_(patchI))
                {
                    continue;
                }

                const coupleMap& cMap = patchCoupling_[patchI].patchMap();

                Map<label>& rPointMap = cMap.reverseEntityMap(coupleMap::POINT);

                bool found = false;
                label mPoint = -1;

                if (rPointMap.found(edges_[pointEdges[edgeI]].end()))
                {
                    mPoint = rPointMap[edges_[pointEdges[edgeI]].end()];

                    // Now check the master map.
                    if (pMap.found(mPoint))
                    {
                        found = true;
                    }
                }

                // Add the map entry for the point.
                if (found)
                {
                    map.addEdge
                    (
                        pointEdges[edgeI],
                        mPoint
                    );
                }
            }
        }
        else
        {
            map.addEdge(pointEdges[edgeI]);
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

        if (cellI == 0)
        {
            for (label i = 0; i < 3; i++)
            {
                // Update the cell list with newly configured cells.
                cells_[newCellIndex[i]] = newTetCell[i];
            }
        }
        else
        {
            for (label i = 3; i < 6; i++)
            {
                // Update the cell list with newly configured cells.
                cells_[newCellIndex[i]] = newTetCell[i];
            }
        }

        removeCell(cIndex);
    }

    // Now finally remove the face...
    removeFace(fIndex);

    if (debug > 2)
    {
        Info << "New Point:: " << newPointIndex << endl;

        const labelList& pEdges = pointEdges_[newPointIndex];

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
                 << faces_[newFaceIndex[faceI]]
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

    // Specify that the operation was successful
    map.type() = 1;

    // Return the changeMap
    return map;
}

// Slice the mesh at a particular location
void dynamicTopoFvMesh::sliceMesh
(
    const labelPair& pointPair
)
{
    if (debug > 1)
    {
        Info << nl << nl
             << "Pair: " << pointPair
             << " is to be used for mesh slicing. " << endl;
    }

    label patchIndex = -1;
    scalar dx = 0.0;
    vector gCentre = vector::zero;

    if (twoDMesh_)
    {
        gCentre =
        (
            0.5 *
            (
                quadFaceCentre(pointPair.first())
              + quadFaceCentre(pointPair.second())
            )
        );

        // Specify a search distance
        dx =
        (
            mag
            (
                quadFaceCentre(pointPair.first())
              - quadFaceCentre(pointPair.second())
            )
        );

        patchIndex = whichPatch(pointPair.first());
    }
    else
    {
        // Find the patch that the edge-vertex is connected to.
        const labelList& pEdges = pointEdges_[pointPair.first()];

        forAll(pEdges, edgeI)
        {
            if ((patchIndex = whichEdgePatch(pEdges[edgeI])) > -1)
            {
                break;
            }
        }

        // Specify the centre.
        gCentre =
        (
            0.5 * (points_[pointPair.first()] + points_[pointPair.second()])
        );

        // Specify a search distance
        dx = mag(points_[pointPair.first()] - points_[pointPair.second()]);
    }

    // Is this edge in the vicinity of a previous slice-point?
    forAll(sliceBoxes_, boxI)
    {
        if (sliceBoxes_[boxI].contains(gCentre))
        {
            if (debug > 1)
            {
                Info << nl << nl
                     << "Pair: " << pointPair
                     << " is too close to another slice point. "
                     << endl;
            }

            // Too close to another slice-point. Bail out.
            return;
        }
    }

    // Choose a box around the centre and scan all
    // surface entities that fall into this region.
    boundBox bBox
    (
        gCentre - vector(dx, dx, dx),
        gCentre + vector(dx, dx, dx)
    );

    vector p = vector::zero, N = vector::zero;
    Map<vector> checkPoints, surfFaces;
    Map<edge> checkEdges;

    if (twoDMesh_)
    {
        // Assign plane point / normal
        p = gCentre;

        vector gNorm = quadFaceNormal(faces_[pointPair.first()]);

        gNorm /= (mag(gNorm) + VSMALL);

        // Since this is 2D, assume XY-plane here.
        N = (gNorm ^ vector(0.0, 0.0, 1.0));
    }
    else
    {
        // Prepare surface points / edges for Dijkstra's algorithm
        for (label edgeI = nOldInternalEdges_; edgeI < nEdges_; edgeI++)
        {
            if (edgeFaces_[edgeI].empty())
            {
                continue;
            }

            if (whichEdgePatch(edgeI) == patchIndex)
            {
                const edge& surfaceEdge = edges_[edgeI];

                if
                (
                    (bBox.contains(points_[surfaceEdge[0]])) &&
                    (bBox.contains(points_[surfaceEdge[1]]))
                )
                {
                    checkEdges.insert(edgeI, surfaceEdge);

                    if (!checkPoints.found(surfaceEdge[0]))
                    {
                        checkPoints.insert
                        (
                            surfaceEdge[0],
                            points_[surfaceEdge[0]]
                        );
                    }

                    if (!checkPoints.found(surfaceEdge[1]))
                    {
                        checkPoints.insert
                        (
                            surfaceEdge[1],
                            points_[surfaceEdge[1]]
                        );
                    }

                    // Add surface faces as well.
                    const labelList& eFaces = edgeFaces_[edgeI];

                    forAll(eFaces, faceI)
                    {
                        if
                        (
                            (neighbour_[eFaces[faceI]] == -1) &&
                            (!surfFaces.found(eFaces[faceI]))
                        )
                        {
                            surfFaces.insert
                            (
                                eFaces[faceI],
                                triFaceNormal(faces_[eFaces[faceI]])
                            );
                        }
                    }
                }
            }
        }

        if (debug > 1)
        {
            Info << nl << nl
                 << " Point [0]: " << points_[pointPair.first()] << nl
                 << " Point [1]: " << points_[pointPair.second()] << endl;

            if (debug > 3)
            {
                writeVTK("slicePoints", checkPoints.toc(), 0);
                writeVTK("sliceEdges", checkEdges.toc(), 1);
            }
        }

        // Find the shortest path using Dijkstra's algorithm.
        Map<label> shortestPath;

        bool foundPath =
        (
            Dijkstra
            (
                checkPoints,
                checkEdges,
                pointPair.first(),
                pointPair.second(),
                shortestPath
            )
        );

        // First normalize all face-normals
        forAllIter(Map<vector>, surfFaces, sIter)
        {
            sIter() /= (mag(sIter()) + VSMALL);
        }

        if (foundPath)
        {
            // Next, take cross-products with every other
            // vector in the list, and accumulate.
            forAllIter(Map<vector>, surfFaces, sIterI)
            {
                forAllIter(Map<vector>, surfFaces, sIterJ)
                {
                    if (sIterI.key() != sIterJ.key())
                    {
                        vector n = (sIterI() ^ sIterJ());

                        n /= (mag(n) + VSMALL);

                        // Reverse the vector if necessary
                        if ((N & n) < 0.0)
                        {
                            n = -n;
                        }

                        N += n;
                    }
                }
            }

            N /= (mag(N) + VSMALL);

            // Obtain point position
            p = gCentre;
        }
        else
        {
            // Probably a membrane-type configuration.
            labelHashSet checkCells;

            // Prepare a bounding cylinder with radius dx.
            forAllIter(Map<vector>, surfFaces, sIter)
            {
                const face& thisFace = faces_[sIter.key()];

                if (thisFace.which(pointPair.first()) > -1)
                {
                    N += sIter();
                }
            }

            // Normalize and reverse.
            N /= -(mag(N) + VSMALL);

            vector a0 = points_[pointPair.first()];
            vector a1 = points_[pointPair.second()];
            scalar dist = mag(a1 - a0);

            forAll(cells_, cellI)
            {
                if (cells_[cellI].empty())
                {
                    continue;
                }

                vector x = tetCellCentre(cellI);

                vector rx = (x - a0);
                vector ra = (rx & N)*N;

                // Check if point falls off cylinder ends.
                if (mag(ra) > dist || mag(ra) < 0.0)
                {
                    continue;
                }

                vector r = (rx - ra);

                // Check if the magnitude of 'r' is within radius.
                if (mag(r) < dx)
                {
                    checkCells.insert(cellI);
                }
            }

            labelList cList = checkCells.toc();

            if (debug > 1)
            {
                Info << "Dijkstra's algorithm could not find a path." << endl;

                if (debug > 3)
                {
                    writeVTK("checkCells", cList, 3);
                }
            }

            removeCells(cList, patchIndex);

            checkConnectivity(10);

            // Add an entry to sliceBoxes.
            label currentSize = sliceBoxes_.size();

            sliceBoxes_.setSize(currentSize + 1);

            sliceBoxes_[currentSize] = bBox;

            return;
        }
    }

    // if (debug > 1)
    {
        Info << nl << nl
             << " Plane point: " << p << nl
             << " Plane normal: " << N << endl;
    }

    // Mark cells and interior faces that fall
    // within the bounding box.
    labelHashSet checkCells, checkFaces, splitFaces;
    Map<bool> cellColors;

    forAll(faces_, faceI)
    {
        if (faces_[faceI].empty())
        {
            continue;
        }

        if (twoDMesh_ && faces_[faceI].size() == 3)
        {
            continue;
        }

        vector fCentre = vector::zero;

        if (twoDMesh_)
        {
            fCentre = quadFaceCentre(faceI);
        }
        else
        {
            fCentre = triFaceCentre(faceI);
        }

        FixedList<label, 2> cellsToCheck(-1);
        cellsToCheck[0] = owner_[faceI];
        cellsToCheck[1] = neighbour_[faceI];

        if (bBox.contains(fCentre) && cellsToCheck[1] != -1)
        {
            // Add this internal face to the list.
            checkFaces.insert(faceI);

            vector centre = vector::zero;

            forAll(cellsToCheck, cellI)
            {
                if (!checkCells.found(cellsToCheck[cellI]))
                {
                    if (twoDMesh_)
                    {
                        centre = prismCellCentre(cellsToCheck[cellI]);
                    }
                    else
                    {
                        centre = tetCellCentre(cellsToCheck[cellI]);
                    }

                    checkCells.insert(cellsToCheck[cellI]);

                    if (((centre - p) & N) > 0.0)
                    {
                        cellColors.insert(cellsToCheck[cellI], true);
                    }
                    else
                    {
                        cellColors.insert(cellsToCheck[cellI], false);
                    }
                }
            }
        }
    }

    // Prepare a list of internal faces for mesh splitting.
    forAllIter(labelHashSet, checkFaces, fIter)
    {
        if
        (
            cellColors[owner_[fIter.key()]]
         != cellColors[neighbour_[fIter.key()]]
        )
        {
            splitFaces.insert(fIter.key());
        }

        // Loop through all points (and associated pointEdges)
        // for this face, and check if connected cells are also
        // present in the checkCells/cellColors list
        if (twoDMesh_)
        {
            const labelList& fEdges = faceEdges_[fIter.key()];

            forAll(fEdges, edgeI)
            {
                const labelList& eFaces = edgeFaces_[fEdges[edgeI]];

                forAll(eFaces, faceI)
                {
                    label own = owner_[eFaces[faceI]];
                    label nei = neighbour_[eFaces[faceI]];

                    if (!checkCells.found(own))
                    {
                        vector centre = prismCellCentre(own);

                        checkCells.insert(own);

                        if (((centre - p) & N) > 0.0)
                        {
                            cellColors.insert(own, true);
                        }
                        else
                        {
                            cellColors.insert(own, false);
                        }
                    }

                    if (!checkCells.found(nei) && nei != -1)
                    {
                        vector centre = prismCellCentre(nei);

                        checkCells.insert(nei);

                        if (((centre - p) & N) > 0.0)
                        {
                            cellColors.insert(nei, true);
                        }
                        else
                        {
                            cellColors.insert(nei, false);
                        }
                    }
                }
            }
        }
        else
        {
            const face& faceToCheck = faces_[fIter.key()];

            forAll(faceToCheck, pointI)
            {
                const labelList& pEdges = pointEdges_[faceToCheck[pointI]];

                forAll(pEdges, edgeI)
                {
                    const labelList& eFaces = edgeFaces_[pEdges[edgeI]];

                    forAll(eFaces, faceI)
                    {
                        label own = owner_[eFaces[faceI]];
                        label nei = neighbour_[eFaces[faceI]];

                        if (!checkCells.found(own))
                        {
                            vector centre = tetCellCentre(own);

                            checkCells.insert(own);

                            if (((centre - p) & N) > 0.0)
                            {
                                cellColors.insert(own, true);
                            }
                            else
                            {
                                cellColors.insert(own, false);
                            }
                        }

                        if (!checkCells.found(nei) && nei != -1)
                        {
                            vector centre = tetCellCentre(nei);

                            checkCells.insert(nei);

                            if (((centre - p) & N) > 0.0)
                            {
                                cellColors.insert(nei, true);
                            }
                            else
                            {
                                cellColors.insert(nei, false);
                            }
                        }
                    }
                }
            }
        }
    }

    if (debug > 3)
    {
        writeVTK("splitFaces", splitFaces.toc(), 2);
        writeVTK("checkCells", checkCells.toc(), 3);
    }

    // Pass this info into the splitInternalFaces routine.
    splitInternalFaces
    (
        patchIndex,
        splitFaces.toc(),
        cellColors
    );

    // Add an entry to sliceBoxes.
    label currentSize = sliceBoxes_.size();

    sliceBoxes_.setSize(currentSize + 1);

    sliceBoxes_[currentSize] = bBox;
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
bool dynamicTopoFvMesh::Dijkstra
(
    const Map<point>& points,
    const Map<edge>& edges,
    const label startPoint,
    const label endPoint,
    Map<label>& pi
) const
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

            sizeUpList(eIter.key(), localPointEdges[edgeToCheck[pointI]]);
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

    // Write out the path
    // if (debug > 3)
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

    return foundEndPoint;
}

// Split a set of internal faces into boundary faces
//   - Add boundary faces and edges to the patch specified by 'patchIndex'
//   - Cell color should specify a binary value dictating either side
//     of the split face.
void dynamicTopoFvMesh::splitInternalFaces
(
    const label patchIndex,
    const labelList& internalFaces,
    const Map<bool>& cellColors
)
{
    Map<label> mirrorPointLabels;
    FixedList<Map<label>, 2> mirrorEdgeLabels, mirrorFaceLabels;

    // First loop through the list and accumulate a list of
    // points and edges that need to be duplicated.
    forAll(internalFaces, faceI)
    {
        const face& faceToCheck = faces_[internalFaces[faceI]];

        forAll(faceToCheck, pointI)
        {
            if (!mirrorPointLabels.found(faceToCheck[pointI]))
            {
                mirrorPointLabels.insert(faceToCheck[pointI], -1);
            }
        }

        const labelList& fEdges = faceEdges_[internalFaces[faceI]];

        forAll(fEdges, edgeI)
        {
            if (!mirrorEdgeLabels[0].found(fEdges[edgeI]))
            {
                mirrorEdgeLabels[0].insert(fEdges[edgeI], -1);
            }
        }
    }

    // Now for every point in the list, add a new one.
    // Add a mapping entry as well.
    forAllIter(Map<label>, mirrorPointLabels, pIter)
    {
        // Obtain a copy of the point before adding it,
        // since the reference might become invalid during list resizing.
        point newPoint = points_[pIter.key()];
        point oldPoint = oldPoints_[pIter.key()];

        pIter() = insertPoint(newPoint, oldPoint, labelList(1, pIter.key()));

        if (!twoDMesh_)
        {
            const labelList& pEdges = pointEdges_[pIter.key()];

            labelHashSet edgesToRemove;

            forAll(pEdges, edgeI)
            {
                const labelList& eFaces = edgeFaces_[pEdges[edgeI]];

                bool allTrue = true;

                forAll(eFaces, faceI)
                {
                    label own = owner_[eFaces[faceI]];
                    label nei = neighbour_[eFaces[faceI]];

                    // Check if an owner/neighbour cell is false
                    if (!cellColors[own])
                    {
                        allTrue = false;
                        break;
                    }

                    if (nei != -1)
                    {
                        if (!cellColors[nei])
                        {
                            allTrue = false;
                            break;
                        }
                    }
                }

                if (allTrue)
                {
                    // Mark this edge label to be discarded later
                    edgesToRemove.insert(pEdges[edgeI]);
                }
            }

            // It is dangerous to use the pointEdges references,
            // so call it using array-lookup instead.
            forAllIter(labelHashSet, edgesToRemove, hsIter)
            {
                // Add the edge to the mirror point list
                sizeUpList
                (
                    hsIter.key(),
                    pointEdges_[pIter()]
                );

                // Remove the edge from the original point list
                sizeDownList
                (
                    hsIter.key(),
                    pointEdges_[pIter.key()]
                );
            }
        }
    }

    if (debug > 3)
    {
        label i = 0;
        labelList mPoints(mirrorPointLabels.size());

        if (!twoDMesh_)
        {
            forAllIter(Map<label>, mirrorPointLabels, pIter)
            {
                writeVTK
                (
                    "pEdges_o_" + Foam::name(pIter.key()) + '_',
                    pointEdges_[pIter.key()],
                    1
                );

                writeVTK
                (
                    "pEdges_m_" + Foam::name(pIter()) + '_',
                    pointEdges_[pIter()],
                    1
                );

                mPoints[i++] = pIter();
            }

            writeVTK
            (
                "points_o_",
                mirrorPointLabels.toc(),
                0
            );

            writeVTK
            (
                "points_m_",
                mPoints,
                0
            );
        }
    }

    // For every internal face, add a new one.
    //  - Stick to the rule:
    //    [1] Every cell marked false keeps the existing entities.
    //    [2] Every cell marked true gets new points/edges/faces.
    //  - If faces are improperly oriented, reverse them.
    forAll(internalFaces, faceI)
    {
        FixedList<face, 2> newFace;
        FixedList<label, 2> newFaceIndex(-1);
        FixedList<label, 2> newOwner(-1);

        label oldOwn = owner_[internalFaces[faceI]];
        label oldNei = neighbour_[internalFaces[faceI]];

        if (cellColors[oldOwn] && !cellColors[oldNei])
        {
            // The owner gets a new boundary face.
            // Note that orientation is already correct.
            newFace[0] = faces_[internalFaces[faceI]];

            // The neighbour needs to have its face reversed
            // and moved to the boundary patch, thereby getting
            // deleted in the process.
            newFace[1] = newFace[0].reverseFace();

            newOwner[0] = oldOwn;
            newOwner[1] = oldNei;
        }
        else
        if (!cellColors[oldOwn] && cellColors[oldNei])
        {
            // The neighbour gets a new boundary face.
            // The face is oriented in the opposite sense, however.
            newFace[0] = faces_[internalFaces[faceI]].reverseFace();

            // The owner keeps the existing face and orientation.
            // But it also needs to be moved to the boundary.
            newFace[1] = faces_[internalFaces[faceI]];

            newOwner[0] = oldNei;
            newOwner[1] = oldOwn;
        }
        else
        {
            // Something's wrong here.
            FatalErrorIn
            (
                "dynamicTopoFvMesh::splitInternalFaces()"
            )
                << nl << " Face: "
                << internalFaces[faceI]
                << " has cells which are improperly marked: " << nl
                << oldOwn << ":: " << cellColors[oldOwn] << nl
                << oldNei << ":: " << cellColors[oldNei]
                << abort(FatalError);
        }

        // Renumber point labels for the first new face.
        forAll(newFace[0], pointI)
        {
            newFace[0][pointI] = mirrorPointLabels[newFace[0][pointI]];
        }

        // Insert the new boundary faces.
        forAll(newFace, indexI)
        {
            newFaceIndex[indexI] =
            (
                insertFace
                (
                    patchIndex,
                    newFace[indexI],
                    newOwner[indexI],
                    -1
                )
            );

            // Make an identical faceEdges entry.
            // This will be renumbered once new edges are added.
            labelList newFaceEdges(faceEdges_[internalFaces[faceI]]);

            faceEdges_.append(newFaceEdges);

            // Replace face labels on cells
            replaceLabel
            (
                internalFaces[faceI],
                newFaceIndex[indexI],
                cells_[newOwner[indexI]]
            );

            // Make face mapping entries for posterity.
            mirrorFaceLabels[indexI].insert
            (
                internalFaces[faceI],
                newFaceIndex[indexI]
            );
        }
    }

    // For every edge in the list, add a new one.
    // We'll deal with correcting edgeFaces and edgePoints later.
    forAllIter(Map<label>, mirrorEdgeLabels[0], eIter)
    {
        // Obtain copies for the append method
        edge origEdge = edges_[eIter.key()];
        labelList newEdgeFaces(edgeFaces_[eIter.key()]);

        eIter() =
        (
            insertEdge
            (
                patchIndex,
                edge
                (
                    mirrorPointLabels[origEdge[0]],
                    mirrorPointLabels[origEdge[1]]
                ),
                newEdgeFaces,
                labelList(0)
            )
        );

        // Is the original edge an internal one?
        // If it is, we need to move it to the boundary.
        if (whichEdgePatch(eIter.key()) == -1)
        {
            label rplEdgeIndex =
            (
                insertEdge
                (
                    patchIndex,
                    origEdge,
                    newEdgeFaces,
                    labelList(0)
                )
            );

            // Map the new entry.
            mirrorEdgeLabels[1].insert(eIter.key(), rplEdgeIndex);
        }
        else
        {
            // This is already a boundary edge.
            // Make an identical map.
            mirrorEdgeLabels[1].insert(eIter.key(), eIter.key());
        }
    }

    // Correct edgeFaces for all new edges
    forAll(mirrorEdgeLabels, indexI)
    {
        forAllIter(Map<label>, mirrorEdgeLabels[indexI], eIter)
        {
            labelList& eFaces = edgeFaces_[eIter()];

            labelHashSet facesToRemove;

            forAll(eFaces, faceI)
            {
                bool remove = false;

                label own = owner_[eFaces[faceI]];
                label nei = neighbour_[eFaces[faceI]];

                if
                (
                    (!cellColors[own] && !indexI) ||
                    ( cellColors[own] &&  indexI)
                )
                {
                    remove = true;
                }

                if (nei != -1)
                {
                    if
                    (
                        (!cellColors[nei] && !indexI) ||
                        ( cellColors[nei] &&  indexI)
                    )
                    {
                        remove = true;
                    }
                }

                if (mirrorFaceLabels[indexI].found(eFaces[faceI]))
                {
                    // Perform a replacement instead of a removal.
                    eFaces[faceI] = mirrorFaceLabels[indexI][eFaces[faceI]];

                    remove = false;
                }

                if (remove)
                {
                    facesToRemove.insert(eFaces[faceI]);
                }
            }

            // Now successively size down edgeFaces.
            // It is dangerous to use the eFaces reference,
            // so call it using array-lookup.
            forAllIter(labelHashSet, facesToRemove, hsIter)
            {
                sizeDownList(hsIter.key(), edgeFaces_[eIter()]);
            }
        }
    }

    // Renumber faceEdges for all faces connected to new edges
    forAll(mirrorEdgeLabels, indexI)
    {
        forAllIter(Map<label>, mirrorEdgeLabels[indexI], eIter)
        {
            const labelList& eFaces = edgeFaces_[eIter()];

            forAll(eFaces, faceI)
            {
                labelList& fEdges = faceEdges_[eFaces[faceI]];

                forAll(fEdges, edgeI)
                {
                    if (mirrorEdgeLabels[indexI].found(fEdges[edgeI]))
                    {
                        fEdges[edgeI] =
                        (
                            mirrorEdgeLabels[indexI][fEdges[edgeI]]
                        );
                    }
                }
            }
        }
    }

    if (twoDMesh_)
    {
        // Renumber edges and faces
        forAllIter(Map<label>, mirrorEdgeLabels[0], eIter)
        {
            const labelList& eFaces = edgeFaces_[eIter()];

            // Two levels of indirection to ensure
            // that all entities we renumbered.
            // A flip-side for the lack of a pointEdges list in 2D.
            forAll(eFaces, faceI)
            {
                const labelList& fEdges = faceEdges_[eFaces[faceI]];

                forAll(fEdges, edgeI)
                {
                    // Renumber this edge.
                    edge& edgeToCheck = edges_[fEdges[edgeI]];

                    forAll(edgeToCheck, pointI)
                    {
                        if (mirrorPointLabels.found(edgeToCheck[pointI]))
                        {
                            edgeToCheck[pointI] =
                            (
                                mirrorPointLabels[edgeToCheck[pointI]]
                            );
                        }
                    }

                    // Also renumber faces connected to this edge.
                    const labelList& efFaces = edgeFaces_[fEdges[edgeI]];

                    forAll(efFaces, faceJ)
                    {
                        face& faceToCheck = faces_[efFaces[faceJ]];

                        forAll(faceToCheck, pointI)
                        {
                            if
                            (
                                mirrorPointLabels.found(faceToCheck[pointI])
                            )
                            {
                                faceToCheck[pointI] =
                                (
                                    mirrorPointLabels[faceToCheck[pointI]]
                                );
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        // Point renumbering of entities connected to mirror points
        forAllIter(Map<label>, mirrorPointLabels, pIter)
        {
            const labelList& pEdges = pointEdges_[pIter()];

            forAll(pEdges, edgeI)
            {
                // Renumber this edge.
                edge& edgeToCheck = edges_[pEdges[edgeI]];

                forAll(edgeToCheck, pointI)
                {
                    if (edgeToCheck[pointI] == pIter.key())
                    {
                        edgeToCheck[pointI] = pIter();
                    }
                }

                // Also renumber faces connected to this edge.
                const labelList& eFaces = edgeFaces_[pEdges[edgeI]];

                forAll(eFaces, faceI)
                {
                    face& faceToCheck = faces_[eFaces[faceI]];

                    forAll(faceToCheck, pointI)
                    {
                        if (faceToCheck[pointI] == pIter.key())
                        {
                            faceToCheck[pointI] = pIter();
                        }
                    }
                }
            }
        }

        // Scan edges connected to mirror points,
        // and correct any edgePoints entries for
        // edges connected to the other vertex.
        forAllIter(Map<label>, mirrorPointLabels, pIter)
        {
            const labelList& pEdges = pointEdges_[pIter()];

            forAll(pEdges, edgeI)
            {
                label otherVertex = edges_[pEdges[edgeI]].otherVertex(pIter());

                // Scan edgePoints for edges connected to this point
                const labelList& opEdges = pointEdges_[otherVertex];

                forAll(opEdges, edgeJ)
                {
                    labelList& ePoints = edgePoints_[opEdges[edgeJ]];

                    forAll(ePoints, pointI)
                    {
                        if (mirrorPointLabels.found(ePoints[pointI]))
                        {
                            // Replace this point with the mirror point
                            ePoints[pointI] =
                            (
                                mirrorPointLabels[ePoints[pointI]]
                            );
                        }
                    }
                }
            }
        }

        // Build edgePoints for new edges
        forAll(mirrorEdgeLabels, indexI)
        {
            forAllIter(Map<label>, mirrorEdgeLabels[indexI], eIter)
            {
                buildEdgePoints(eIter());
            }
        }
    }

    // Now that we're done with the internal faces, remove them.
    forAll(internalFaces, faceI)
    {
        removeFace(internalFaces[faceI]);
    }

    // Remove old internal edges as well.
    forAllIter(Map<label>, mirrorEdgeLabels[1], eIter)
    {
        if (eIter.key() != eIter())
        {
            removeEdge(eIter.key());
        }
    }

    // Set the flag
    topoChangeFlag_ = true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
