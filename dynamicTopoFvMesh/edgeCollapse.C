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
#include "multiThreader.H"
#include "dynamicTopoFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Method for the collapse of a quad-face in 2D
// - Returns a changeMap with a type specifying:
//    -1: Collapse failed since max number of topo-changes was reached.
//     0: Collapse could not be performed.
//     1: Collapsed to first node.
//     2: Collapsed to second node.
// - overRideCase is used to force a certain collapse configuration.
//    -1: Use this value to let collapseQuadFace decide a case.
//     1: Force collapse to first node.
//     2: Force collapse to second node.
// - checkOnly performs a feasibility check and returns without modifications.
const changeMap dynamicTopoFvMesh::collapseQuadFace
(
    const label fIndex,
    label overRideCase,
    bool checkOnly
)
{
    // Figure out which thread this is...
    label tIndex = self();

    // Prepare the changeMap
    changeMap map;

    if
    (
        (nModifications_ > maxModifications_)
     && (maxModifications_ > -1)
    )
    {
        // Reached the max allowable topo-changes.
        faceStack(tIndex).clear();

        return map;
    }

    // Sanity check: Is the index legitimate?
    if (fIndex < 0 || fIndex >= nFaces_)
    {
        FatalErrorIn("dynamicTopoFvMesh::collapseQuadFace()")
            << " Invalid index: " << fIndex
            << abort(FatalError);
    }

    // Local variables
    FixedList<bool,2> edgeBoundary(false);
    FixedList<label,2> c0BdyIndex, c0IntIndex, c1BdyIndex, c1IntIndex;
    FixedList<face,2>  c0BdyFace,  c0IntFace,  c1BdyFace,  c1IntFace;
    FixedList<edge,4> checkEdge(edge(-1,-1));
    FixedList<label,4> checkEdgeIndex(-1);

    // Define checkEdges
    checkEdgeIndex[0] = getTriBoundaryEdge(fIndex);
    checkEdge[0] = edges_[checkEdgeIndex[0]];

    const labelList& fEdges = faceEdges_[fIndex];

    forAll(fEdges, edgeI)
    {
        if (checkEdgeIndex[0] != fEdges[edgeI])
        {
            const edge& thisEdge = edges_[fEdges[edgeI]];

            if
            (
                checkEdge[0].start() == thisEdge[0] ||
                checkEdge[0].start() == thisEdge[1]
            )
            {
                checkEdgeIndex[1] = fEdges[edgeI];
                checkEdge[1] = thisEdge;

                // Update the map
                map.firstEdge() = checkEdgeIndex[1];
            }
            else
            if
            (
                checkEdge[0].end() == thisEdge[0] ||
                checkEdge[0].end() == thisEdge[1]
            )
            {
                checkEdgeIndex[2] = fEdges[edgeI];
                checkEdge[2] = thisEdge;

                // Update the map
                map.secondEdge() = checkEdgeIndex[2];
            }
            else
            {
                checkEdgeIndex[3] = fEdges[edgeI];
                checkEdge[3] = thisEdge;
            }
        }
    }

    // If coupled modification is set, and this is a
    // master edge, collapse its slaves first.
    if (coupledModification_)
    {
        // Is this a locally coupled edge?
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
                        break;
                    }
                }
            }

            if (sIndex == -1)
            {
                FatalErrorIn("dynamicTopoFvMesh::collapseQuadFace()")
                    << "Coupled maps were improperly specified." << nl
                    << " Slave index not found for: " << nl
                    << " Face: " << fIndex << nl
                    << abort(FatalError);
            }

            // Temporarily turn off coupledModification.
            unsetCoupledModification();

            // First check the slave for collapse feasibility.
            changeMap slaveMap = collapseQuadFace(sIndex, -1, true);

            if (slaveMap.type() > 0)
            {
                const label pointEnum = coupleMap::POINT;

                // Obtain the pointMap
                Map<label>& pointMap =
                (
                    patchCoupling_[whichPatch(fIndex)].patchMap().entityMap
                    (
                        pointEnum
                    )
                );

                FixedList<edge, 2> mEdge(edge(-1, -1)), sEdge(edge(-1, -1));

                mEdge[0][0] = pointMap[checkEdge[1].start()];
                mEdge[0][1] = pointMap[checkEdge[1].end()];

                mEdge[1][0] = pointMap[checkEdge[2].start()];
                mEdge[1][1] = pointMap[checkEdge[2].end()];

                sEdge[0] = edges_[slaveMap.firstEdge()];
                sEdge[1] = edges_[slaveMap.secondEdge()];

                // Set the overRideCase for this edge
                changeMap masterMap;

                // Perform a geometric comparison.
                switch (slaveMap.type())
                {
                    case 1:

                        if (mEdge[0] == sEdge[0])
                        {
                            overRideCase = 1;
                        }
                        else
                        if (mEdge[1] == sEdge[0])
                        {
                            overRideCase = 2;
                        }
                        else
                        {
                            FatalErrorIn("dynamicTopoFvMesh::collapseEdge()")
                                << "Coupled collapse failed." << nl
                                << "Masters: " << nl
                                << checkEdgeIndex[1] << ": "
                                << checkEdge[1] << nl
                                << checkEdgeIndex[2] << ": "
                                << checkEdge[2] << nl
                                << "Slaves: " << nl
                                << slaveMap.firstEdge() << ": "
                                << sEdge[0] << nl
                                << slaveMap.secondEdge() << ": "
                                << sEdge[1] << nl
                                << abort(FatalError);
                        }

                        break;

                    case 2:

                        if (mEdge[1] == sEdge[1])
                        {
                            overRideCase = 2;
                        }
                        else
                        if (mEdge[0] == sEdge[1])
                        {
                            overRideCase = 1;
                        }
                        else
                        {
                            FatalErrorIn("dynamicTopoFvMesh::collapseEdge()")
                                << "Coupled collapse failed." << nl
                                << "Masters: " << nl
                                << checkEdgeIndex[1] << ": "
                                << checkEdge[1] << nl
                                << checkEdgeIndex[2] << ": "
                                << checkEdge[2] << nl
                                << "Slaves: " << nl
                                << slaveMap.firstEdge() << ": "
                                << sEdge[0] << nl
                                << slaveMap.secondEdge() << ": "
                                << sEdge[1] << nl
                                << abort(FatalError);
                        }

                        break;
                }

                // Can the overRideCase be used for this edge?
                masterMap = collapseQuadFace(fIndex, overRideCase, true);

                // Master couldn't perform collapse.
                if (masterMap.type() <= 0)
                {
                    // Turn coupledModification back on before bailing out.
                    setCoupledModification();

                    return masterMap;
                }

                // Collapse the slave.
                collapseQuadFace(sIndex);
            }
            else
            {
                // Slave couldn't perform collapse.
                setCoupledModification();

                map.type() = 0;

                return map;
            }

            // Turn it back on.
            setCoupledModification();
        }
        else
        if (processorCoupledFace(fIndex))
        {
            // Collapse face on the patchSubMesh.

        }
    }

    // Determine if either edge belongs to a boundary
    edgeBoundary[0] = (whichEdgePatch(checkEdgeIndex[1]) > -1);
    edgeBoundary[1] = (whichEdgePatch(checkEdgeIndex[2]) > -1);

    if (debug > 1)
    {
        const labelList& fE = faceEdges_[fIndex];

        Info << nl << nl
             << "Face: " << fIndex << ": " << faces_[fIndex] << nl
             << "faceEdges: " << fE
             << " is to be collapsed. " << endl;

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
            const cell& firstCurCell = cells_[firstCells[cellI]];

            Info << "Cell: " << firstCells[cellI]
                 << ": " << firstCurCell << endl;

            forAll(firstCurCell,faceI)
            {
                Info << firstCurCell[faceI]
                     << ": " << faces_[firstCurCell[faceI]] << endl;
            }
        }

        const labelList& firstEdgeFaces = edgeFaces_[checkEdgeIndex[1]];

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
            const cell& secondCurCell = cells_[secondCells[cellI]];

            Info << "Cell: " << secondCells[cellI]
                 << ": " << secondCurCell << endl;

            forAll(secondCurCell, faceI)
            {
                Info << secondCurCell[faceI]
                     << ": " << faces_[secondCurCell[faceI]] << endl;
            }
        }

        const labelList& secondEdgeFaces = edgeFaces_[checkEdgeIndex[2]];

        Info << nl << "Second Edge Face Hull: " << secondEdgeFaces << endl;

        forAll(secondEdgeFaces, indexI)
        {
            Info << secondEdgeFaces[indexI]
                 << ": " << faces_[secondEdgeFaces[indexI]] << endl;
        }

        // Write out VTK files prior to change
        if (debug > 3)
        {
            labelHashSet vtkCells;

            forAll(firstCells, cellI)
            {
                if (!vtkCells.found(firstCells[cellI]))
                {
                    vtkCells.insert(firstCells[cellI]);
                }
            }

            forAll(secondCells, cellI)
            {
                if (!vtkCells.found(secondCells[cellI]))
                {
                    vtkCells.insert(secondCells[cellI]);
                }
            }

            writeVTK
            (
                Foam::name(fIndex)
              + "_Collapse_0",
                vtkCells.toc()
            );
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

    // Lists for feasibility checks
    label collapseCase = -1;
    FixedList<bool, 2> nBoundCurves(false);
    FixedList<label,2> original(-1), replacement(-1), ends(-1);
    FixedList<label,2> faceToKeep(-1), faceToThrow(-1);
    FixedList<label,4> edgeToKeep(-1), edgeToThrow(-1);

    // Set the collapseCase
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
        if (c1 != -1)
        {
            if (debug > 2)
            {
                WarningIn
                (
                    "dynamicTopoFvMesh::collapseQuadFace"
                )   << "Collapsing an internal face that "
                    << "lies on two boundary patches. "
                    << "Face: " << fIndex << ": " << faces_[fIndex] << endl;
            }

            // Bail out for now. If proximity based refinement is
            // switched on, mesh may be sliced at this point.
            return map;
        }

        // Check if either edge lies on a bounding curve.
        if (checkBoundingCurve(checkEdgeIndex[1]))
        {
            nBoundCurves[0] = true;
        }

        if (checkBoundingCurve(checkEdgeIndex[2]))
        {
            nBoundCurves[1] = true;
        }

        // Collapse towards a bounding curve
        if (nBoundCurves[0] && !nBoundCurves[1])
        {
            collapseCase = 1;
        }
        else
        if (!nBoundCurves[0] && nBoundCurves[1])
        {
            collapseCase = 2;
        }
        else
        {
            // Two bounding curves? This might cause pinching.
            // Bail out for now.
            return map;
        }
    }

    // Perform an override if requested.
    if (overRideCase != -1)
    {
        collapseCase = overRideCase;
    }

    if (collapseCase == 2)
    {
        original[0] = cv0; original[1] = cv1;
        replacement[0] = cv2; replacement[1] = cv3;

        // Check whether the collapse is possible.
        if
        (
            checkCollapse
            (
                firstTriFaces,
                c0BdyIndex,
                c1BdyIndex,
                original,
                replacement,
                (c1 != -1)
            )
        )
        {
            map.type() = 0;
            return map;
        }

        // Are we only performing checks?
        if (checkOnly)
        {
            map.type() = collapseCase;
            return map;
        }

        const labelList& firstEdgeFaces = edgeFaces_[checkEdgeIndex[1]];

        // Collapse to the second node...
        forAll(firstEdgeFaces,faceI)
        {
            // Replace point indices on faces.
            replaceLabel(cv0, cv2, faces_[firstEdgeFaces[faceI]]);
            replaceLabel(cv1, cv3, faces_[firstEdgeFaces[faceI]]);

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

        // Find common edges between quad / tri faces...
        findCommonEdge(c0BdyIndex[0], faceToKeep[0], edgeToKeep[0]);
        findCommonEdge(c0BdyIndex[1], faceToKeep[0], edgeToKeep[1]);
        findCommonEdge(c0BdyIndex[0], faceToThrow[0], edgeToThrow[0]);
        findCommonEdge(c0BdyIndex[1], faceToThrow[0], edgeToThrow[1]);

        // Size down edgeFaces for the ends.
        findCommonEdge(faceToThrow[0], faceToKeep[0], ends[0]);
        sizeDownList(faceToThrow[0], edgeFaces_[ends[0]]);

        if (c1 != -1)
        {
            findCommonEdge(c1BdyIndex[0], faceToKeep[1], edgeToKeep[2]);
            findCommonEdge(c1BdyIndex[1], faceToKeep[1], edgeToKeep[3]);
            findCommonEdge(c1BdyIndex[0], faceToThrow[1], edgeToThrow[2]);
            findCommonEdge(c1BdyIndex[1], faceToThrow[1], edgeToThrow[3]);

            // Size down edgeFaces for the ends.
            findCommonEdge(faceToThrow[1], faceToKeep[1], ends[1]);
            sizeDownList(faceToThrow[1], edgeFaces_[ends[1]]);
        }

        // Correct edgeFaces for triangular faces...
        forAll(edgeToThrow, indexI)
        {
            if (edgeToThrow[indexI] == -1)
            {
                continue;
            }

            const labelList& eF = edgeFaces_[edgeToThrow[indexI]];

            label origTriFace = -1, retTriFace = -1;

            // Find the original triangular face index
            forAll(eF, faceI)
            {
                if (eF[faceI] == c0BdyIndex[0])
                {
                    origTriFace = c0BdyIndex[0];
                    break;
                }

                if (eF[faceI] == c0BdyIndex[1])
                {
                    origTriFace = c0BdyIndex[1];
                    break;
                }

                if (eF[faceI] == c1BdyIndex[0])
                {
                    origTriFace = c1BdyIndex[0];
                    break;
                }

                if (eF[faceI] == c1BdyIndex[1])
                {
                    origTriFace = c1BdyIndex[1];
                    break;
                }
            }

            // Find the retained triangular face index
            forAll(eF, faceI)
            {
                if
                (
                    (eF[faceI] != origTriFace) &&
                    (eF[faceI] != faceToThrow[0]) &&
                    (eF[faceI] != faceToThrow[1])
                )
                {
                    retTriFace = eF[faceI];
                    break;
                }
            }

            // Finally replace the face index
            if (retTriFace == -1)
            {
                // Couldn't find a retained face.
                // This must be a boundary edge, so size-down instead.
                sizeDownList(origTriFace, edgeFaces_[edgeToKeep[indexI]]);
            }
            else
            {
                replaceLabel
                (
                    origTriFace,
                    retTriFace,
                    edgeFaces_[edgeToKeep[indexI]]
                );

                replaceLabel
                (
                    edgeToThrow[indexI],
                    edgeToKeep[indexI],
                    faceEdges_[retTriFace]
                );
            }
        }

        // Correct faceEdges / edgeFaces for quad-faces...
        forAll(firstEdgeFaces,faceI)
        {
            replaceLabel
            (
                checkEdgeIndex[1],
                checkEdgeIndex[2],
                faceEdges_[firstEdgeFaces[faceI]]
            );

            // Renumber the edges on this face
            const labelList& fE = faceEdges_[firstEdgeFaces[faceI]];

            forAll(fE, edgeI)
            {
                if (edges_[fE[edgeI]].start() == cv0)
                {
                    edges_[fE[edgeI]].start() = cv2;
                }

                if (edges_[fE[edgeI]].end() == cv0)
                {
                    edges_[fE[edgeI]].end() = cv2;
                }

                if (edges_[fE[edgeI]].start() == cv1)
                {
                    edges_[fE[edgeI]].start() = cv3;
                }

                if (edges_[fE[edgeI]].end() == cv1)
                {
                    edges_[fE[edgeI]].end() = cv3;
                }
            }

            // Size-up edgeFaces for the replacement
            if
            (
                (firstEdgeFaces[faceI] != faceToThrow[0]) &&
                (firstEdgeFaces[faceI] != faceToThrow[1]) &&
                (firstEdgeFaces[faceI] != fIndex)
            )
            {
                sizeUpList
                (
                    firstEdgeFaces[faceI],
                    edgeFaces_[checkEdgeIndex[2]]
                );
            }
        }

        // Remove the current face from the replacement edge
        sizeDownList(fIndex, edgeFaces_[checkEdgeIndex[2]]);

        // Replace point labels on all triangular boundary faces.
        forAll(firstCells,cellI)
        {
            const cell& cellToCheck = cells_[firstCells[cellI]];

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

        // Now that we're done with all edges, remove them.
        removeEdge(checkEdgeIndex[0]);
        removeEdge(checkEdgeIndex[1]);
        removeEdge(checkEdgeIndex[3]);

        forAll(edgeToThrow, indexI)
        {
            if (edgeToThrow[indexI] == -1)
            {
                continue;
            }

            removeEdge(edgeToThrow[indexI]);
        }

        // Delete the two points...
        removePoint(cv0);
        removePoint(cv1);

        if (coupledModification_)
        {
            if (locallyCoupledFace(fIndex))
            {
                // Remove the point entries.
                const label pEnum = coupleMap::POINT;

                forAll(patchCoupling_, patchI)
                {
                    if (!patchCoupling_(patchI))
                    {
                        continue;
                    }

                    const coupleMap& cMap = patchCoupling_[patchI].patchMap();

                    // Obtain references
                    Map<label>& pointMap = cMap.entityMap(pEnum);
                    Map<label>& rPointMap = cMap.reverseEntityMap(pEnum);

                    if (pointMap.found(cv0))
                    {
                        // Erase the reverse map first
                        rPointMap.erase(pointMap[cv0]);

                        // Update pointMap
                        pointMap.erase(cv0);
                    }

                    if (pointMap.found(cv1))
                    {
                        // Erase the reverse map first
                        rPointMap.erase(pointMap[cv1]);

                        // Update pointMap
                        pointMap.erase(cv1);
                    }
                }
            }
        }
    }
    else
    {
        original[0] = cv2; original[1] = cv3;
        replacement[0] = cv0; replacement[1] = cv1;

        // Check whether the collapse is possible.
        if
        (
            checkCollapse
            (
                secondTriFaces,
                c0BdyIndex,
                c1BdyIndex,
                original,
                replacement,
                (c1 != -1)
            )
        )
        {
            map.type() = 0;
            return map;
        }

        // Are we only performing checks?
        if (checkOnly)
        {
            map.type() = collapseCase;
            return map;
        }

        const labelList& secondEdgeFaces = edgeFaces_[checkEdgeIndex[2]];

        // Collapse to the first node by default...
        forAll(secondEdgeFaces,faceI)
        {
            // Replace point indices on faces.
            replaceLabel(cv2, cv0, faces_[secondEdgeFaces[faceI]]);
            replaceLabel(cv3, cv1, faces_[secondEdgeFaces[faceI]]);

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

        // Find common edges between quad / tri faces...
        findCommonEdge(c0BdyIndex[0], faceToKeep[0], edgeToKeep[0]);
        findCommonEdge(c0BdyIndex[1], faceToKeep[0], edgeToKeep[1]);
        findCommonEdge(c0BdyIndex[0], faceToThrow[0], edgeToThrow[0]);
        findCommonEdge(c0BdyIndex[1], faceToThrow[0], edgeToThrow[1]);

        // Size down edgeFaces for the ends.
        findCommonEdge(faceToThrow[0], faceToKeep[0], ends[0]);
        sizeDownList(faceToThrow[0], edgeFaces_[ends[0]]);

        if (c1 != -1)
        {
            findCommonEdge(c1BdyIndex[0], faceToKeep[1], edgeToKeep[2]);
            findCommonEdge(c1BdyIndex[1], faceToKeep[1], edgeToKeep[3]);
            findCommonEdge(c1BdyIndex[0], faceToThrow[1], edgeToThrow[2]);
            findCommonEdge(c1BdyIndex[1], faceToThrow[1], edgeToThrow[3]);

            // Size down edgeFaces for the ends.
            findCommonEdge(faceToThrow[1], faceToKeep[1], ends[1]);
            sizeDownList(faceToThrow[1], edgeFaces_[ends[1]]);
        }

        // Correct edgeFaces for triangular faces...
        forAll(edgeToThrow, indexI)
        {
            if (edgeToThrow[indexI] == -1)
            {
                continue;
            }

            const labelList& eF = edgeFaces_[edgeToThrow[indexI]];

            label origTriFace = -1, retTriFace = -1;

            // Find the original triangular face index
            forAll(eF, faceI)
            {
                if (eF[faceI] == c0BdyIndex[0])
                {
                    origTriFace = c0BdyIndex[0];
                    break;
                }

                if (eF[faceI] == c0BdyIndex[1])
                {
                    origTriFace = c0BdyIndex[1];
                    break;
                }

                if (eF[faceI] == c1BdyIndex[0])
                {
                    origTriFace = c1BdyIndex[0];
                    break;
                }

                if (eF[faceI] == c1BdyIndex[1])
                {
                    origTriFace = c1BdyIndex[1];
                    break;
                }
            }

            // Find the retained triangular face index
            forAll(eF, faceI)
            {
                if
                (
                    (eF[faceI] != origTriFace) &&
                    (eF[faceI] != faceToThrow[0]) &&
                    (eF[faceI] != faceToThrow[1])
                )
                {
                    retTriFace = eF[faceI];
                    break;
                }
            }

            // Finally replace the face/edge indices
            if (retTriFace == -1)
            {
                // Couldn't find a retained face.
                // This must be a boundary edge, so size-down instead.
                sizeDownList(origTriFace, edgeFaces_[edgeToKeep[indexI]]);
            }
            else
            {
                replaceLabel
                (
                    origTriFace,
                    retTriFace,
                    edgeFaces_[edgeToKeep[indexI]]
                );

                replaceLabel
                (
                    edgeToThrow[indexI],
                    edgeToKeep[indexI],
                    faceEdges_[retTriFace]
                );
            }
        }

        // Correct faceEdges / edgeFaces for quad-faces...
        forAll(secondEdgeFaces,faceI)
        {
            replaceLabel
            (
                checkEdgeIndex[2],
                checkEdgeIndex[1],
                faceEdges_[secondEdgeFaces[faceI]]
            );

            // Renumber the edges on this face
            const labelList& fE = faceEdges_[secondEdgeFaces[faceI]];

            forAll(fE, edgeI)
            {
                if (edges_[fE[edgeI]].start() == cv2)
                {
                    edges_[fE[edgeI]].start() = cv0;
                }

                if (edges_[fE[edgeI]].end() == cv2)
                {
                    edges_[fE[edgeI]].end() = cv0;
                }

                if (edges_[fE[edgeI]].start() == cv3)
                {
                    edges_[fE[edgeI]].start() = cv1;
                }

                if (edges_[fE[edgeI]].end() == cv3)
                {
                    edges_[fE[edgeI]].end() = cv1;
                }
            }

            // Size-up edgeFaces for the replacement
            if
            (
                (secondEdgeFaces[faceI] != faceToThrow[0]) &&
                (secondEdgeFaces[faceI] != faceToThrow[1]) &&
                (secondEdgeFaces[faceI] != fIndex)
            )
            {
                sizeUpList
                (
                    secondEdgeFaces[faceI],
                    edgeFaces_[checkEdgeIndex[1]]
                );
            }
        }

        // Remove the current face from the replacement edge
        sizeDownList(fIndex, edgeFaces_[checkEdgeIndex[1]]);

        // Replace point labels on all triangular boundary faces.
        forAll(secondCells, cellI)
        {
            const cell& cellToCheck = cells_[secondCells[cellI]];

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

        // Now that we're done with all edges, remove them.
        removeEdge(checkEdgeIndex[0]);
        removeEdge(checkEdgeIndex[2]);
        removeEdge(checkEdgeIndex[3]);

        forAll(edgeToThrow, indexI)
        {
            if (edgeToThrow[indexI] == -1)
            {
                continue;
            }

            removeEdge(edgeToThrow[indexI]);
        }

        // Delete the two points...
        removePoint(cv2);
        removePoint(cv3);

        if (coupledModification_)
        {
            if (locallyCoupledFace(fIndex))
            {
                // Remove the point entries.
                const label pEnum = coupleMap::POINT;

                forAll(patchCoupling_, patchI)
                {
                    if (!patchCoupling_(patchI))
                    {
                        continue;
                    }

                    const coupleMap& cMap = patchCoupling_[patchI].patchMap();

                    // Obtain references
                    Map<label>& pointMap = cMap.entityMap(pEnum);
                    Map<label>& rPointMap = cMap.reverseEntityMap(pEnum);

                    if (pointMap.found(cv2))
                    {
                        // Erase the reverse map first
                        rPointMap.erase(pointMap[cv2]);

                        // Update pointMap
                        pointMap.erase(cv2);
                    }

                    if (pointMap.found(cv3))
                    {
                        // Erase the reverse map first
                        rPointMap.erase(pointMap[cv3]);

                        // Update pointMap
                        pointMap.erase(cv3);
                    }
                }
            }
        }
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
            const cell& firstCurCell = cells_[firstCells[cellI]];

            Info << "Cell: " << firstCells[cellI]
                 << ": " << firstCurCell << endl;

            forAll(firstCurCell, faceI)
            {
                Info << firstCurCell[faceI]
                     << ": " << faces_[firstCurCell[faceI]] << endl;
            }
        }

        const labelList& firstEdgeFaces = edgeFaces_[checkEdgeIndex[1]];

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
            const cell& secondCurCell = cells_[secondCells[cellI]];

            Info << "Cell: " << secondCells[cellI]
                 << ": " << secondCurCell << endl;

            forAll(secondCurCell, faceI)
            {
                Info << secondCurCell[faceI]
                     << ": " << faces_[secondCurCell[faceI]] << endl;
            }
        }

        const labelList& secondEdgeFaces = edgeFaces_[checkEdgeIndex[2]];

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
                faces_[faceToKeep[0]] =
                (
                    faces_[faceToKeep[0]].reverseFace()
                );

                owner_[faceToKeep[0]] = neighbour_[faceToKeep[0]];
                neighbour_[faceToKeep[0]] = neighbour_[faceToThrow[0]];

                iPtr_->setFlip(faceToKeep[0]);
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
                    faces_[faceToKeep[0]] =
                    (
                        faces_[faceToKeep[0]].reverseFace()
                    );

                    owner_[faceToKeep[0]] = neighbour_[faceToKeep[0]];
                    neighbour_[faceToKeep[0]] = -1;

                    iPtr_->setFlip(faceToKeep[0]);
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
                faces_[faceToKeep[0]] =
                (
                    faces_[faceToKeep[0]].reverseFace()
                );

                neighbour_[faceToKeep[0]] = owner_[faceToKeep[0]];
                owner_[faceToKeep[0]] = owner_[faceToThrow[0]];

                iPtr_->setFlip(faceToKeep[0]);
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
                    faces_[faceToKeep[1]] =
                    (
                        faces_[faceToKeep[1]].reverseFace()
                    );

                    owner_[faceToKeep[1]] = neighbour_[faceToKeep[1]];
                    neighbour_[faceToKeep[1]] = neighbour_[faceToThrow[1]];

                    iPtr_->setFlip(faceToKeep[1]);
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
                        (
                            faces_[faceToKeep[1]].reverseFace()
                        );

                        owner_[faceToKeep[1]] = neighbour_[faceToKeep[1]];
                        neighbour_[faceToKeep[1]] = -1;

                        iPtr_->setFlip(faceToKeep[1]);
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
                    faces_[faceToKeep[1]] =
                    (
                        faces_[faceToKeep[1]].reverseFace()
                    );

                    neighbour_[faceToKeep[1]] = owner_[faceToKeep[1]];
                    owner_[faceToKeep[1]] = owner_[faceToThrow[1]];

                    iPtr_->setFlip(faceToKeep[1]);
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
        // Obtain a copy before adding the new face,
        // since the reference might become invalid during list resizing.
        face newFace = faces_[faceToKeep[0]];
        label newOwn = owner_[faceToKeep[0]];
        labelList newFaceEdges = faceEdges_[faceToKeep[0]];

        // This face is being converted from interior to boundary. Remove
        // from the interior list and add as a boundary face to the end.
        label newFaceIndex = insertFace
                             (
                                 whichPatch(faceToThrow[0]),
                                 newFace,
                                 newOwn,
                                 -1
                             );

        // Add a faceEdges entry as well.
        // Edges don't have to change, since they're
        // all on the boundary anyway.
        faceEdges_.append(newFaceEdges);

        replaceLabel
        (
            faceToKeep[0],
            newFaceIndex,
            cells_[newOwn]
        );

        // Correct edgeFaces with the new face label.
        forAll(newFaceEdges, edgeI)
        {
            replaceLabel
            (
                faceToKeep[0],
                newFaceIndex,
                edgeFaces_[newFaceEdges[edgeI]]
            );
        }

        // Renumber the neighbour so that this face is removed correctly.
        neighbour_[faceToKeep[0]] = 0;
        removeFace(faceToKeep[0]);
    }

    // Remove the unwanted faces in the cell(s) adjacent to this face,
    // and correct the cells that contain discarded faces
    const cell& cell_0 = cells_[c0];

    forAll(cell_0,faceI)
    {
        if (cell_0[faceI] != fIndex && cell_0[faceI] != faceToKeep[0])
        {
           removeFace(cell_0[faceI]);
        }
    }

    if (cellCheck[0] != -1)
    {
        replaceLabel(faceToThrow[0], faceToKeep[0], cells_[cellCheck[0]]);
    }

    // Remove the cell
    removeCell(c0);

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
            // Obtain a copy before adding the new face,
            // since the reference might become invalid during list resizing.
            face newFace = faces_[faceToKeep[1]];
            label newOwn = owner_[faceToKeep[1]];
            labelList newFaceEdges = faceEdges_[faceToKeep[1]];

            // This face is being converted from interior to boundary. Remove
            // from the interior list and add as a boundary face to the end.
            label newFaceIndex = insertFace
                                 (
                                     whichPatch(faceToThrow[1]),
                                     newFace,
                                     newOwn,
                                     -1
                                 );

            // Add a faceEdges entry as well.
            // Edges don't have to change, since they're
            // all on the boundary anyway.
            faceEdges_.append(newFaceEdges);

            replaceLabel
            (
                faceToKeep[1],
                newFaceIndex,
                cells_[newOwn]
            );

            // Correct edgeFaces with the new face label.
            forAll(newFaceEdges, edgeI)
            {
                replaceLabel
                (
                    faceToKeep[1],
                    newFaceIndex,
                    edgeFaces_[newFaceEdges[edgeI]]
                );
            }

            // Renumber the neighbour so that this face is removed correctly.
            neighbour_[faceToKeep[1]] = 0;
            removeFace(faceToKeep[1]);
        }

        const cell& cell_1 = cells_[c1];

        forAll(cell_1, faceI)
        {
            if (cell_1[faceI] != fIndex && cell_1[faceI] != faceToKeep[1])
            {
               removeFace(cell_1[faceI]);
            }
        }

        if (cellCheck[1] != -1)
        {
            replaceLabel(faceToThrow[1], faceToKeep[1], cells_[cellCheck[1]]);
        }

        // Remove the cell
        removeCell(c1);
    }

    // Finally remove the face
    removeFace(fIndex);

    // Write out VTK files after change
    if (debug > 3)
    {
        labelHashSet vtkCells;

        forAll(firstCells, cellI)
        {
            if (firstCells[cellI] == c0 || firstCells[cellI] == c1)
            {
                continue;
            }

            if (!vtkCells.found(firstCells[cellI]))
            {
                vtkCells.insert(firstCells[cellI]);
            }
        }

        forAll(secondCells, cellI)
        {
            if (secondCells[cellI] == c0 || secondCells[cellI] == c1)
            {
                continue;
            }

            if (!vtkCells.found(secondCells[cellI]))
            {
                vtkCells.insert(secondCells[cellI]);
            }
        }

        writeVTK
        (
            Foam::name(fIndex)
          + "_Collapse_1",
            vtkCells.toc()
        );
    }

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    nCollapses_++;

    // Increment the number of modifications
    nModifications_++;

    // Return a succesful collapse
    map.type() = collapseCase;

    return map;
}

// Method for the collapse of an edge in 3D
// - Returns a changeMap with a type specifying:
//    -1: Collapse failed since max number of topo-changes was reached.
//     0: Collapse could not be performed.
//     1: Collapsed to first node.
//     2: Collapsed to second node.
//     3: Collapsed to mid-point (default)
// - overRideCase is used to force a certain collapse configuration.
//    -1: Use this value to let collapseEdge decide a case.
//     1: Force collapse to first node.
//     2: Force collapse to second node.
//     3: Force collapse to mid-point
// - checkOnly performs a feasibility check and returns without modifications.
// - forceOp to force the collapse.
const changeMap dynamicTopoFvMesh::collapseEdge
(
    const label eIndex,
    label overRideCase,
    bool checkOnly,
    bool forceOp
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

    // Prepare the changeMap
    changeMap map;

    if
    (
        (nModifications_ > maxModifications_)
     && (maxModifications_ > -1)
    )
    {
        // Reached the max allowable topo-changes.
        edgeStack(tIndex).clear();

        return map;
    }

    // Sanity check: Is the index legitimate?
    if (eIndex < 0 || eIndex >= nEdges_)
    {
        FatalErrorIn("dynamicTopoFvMesh::collapseEdge()")
            << " Invalid index: " << eIndex
            << abort(FatalError);
    }

    // If coupled modification is set, and this is a
    // master edge, collapse its slaves first.
    if (coupledModification_)
    {
        // Is this a locally coupled edge?
        if (locallyCoupledEdge(eIndex))
        {
            label sIndex = -1, pIndex = -1;

            // Determine the slave index.
            forAll(patchCoupling_, patchI)
            {
                if (patchCoupling_(patchI))
                {
                    const label edgeEnum  = coupleMap::EDGE;
                    const coupleMap& cMap = patchCoupling_[patchI].patchMap();

                    if ((sIndex = cMap.findSlaveIndex(edgeEnum, eIndex)) > -1)
                    {
                        pIndex = patchI;

                        break;
                    }
                }
            }

            if (sIndex == -1)
            {
                FatalErrorIn("dynamicTopoFvMesh::collapseEdge")
                    << "Coupled maps were improperly specified." << nl
                    << " Slave index not found for: " << nl
                    << " Edge: " << eIndex << nl
                    << abort(FatalError);
            }

            // Temporarily turn off coupledModification
            unsetCoupledModification();

            // First check the slave for collapse feasibility.
            changeMap slaveMap = collapseEdge(sIndex, -1, true);

            if (slaveMap.type() > 0)
            {
                edge mEdge = edges_[eIndex];
                edge sEdge = edges_[sIndex];

                // Set the overRideCase for this edge
                changeMap masterMap;

                const label pointEnum = coupleMap::POINT;

                // Obtain the pointMap
                Map<label>& pointMap =
                (
                    patchCoupling_[pIndex].patchMap().entityMap(pointEnum)
                );

                // Perform a topological comparison.
                switch (slaveMap.type())
                {
                    case 1:

                        if (pointMap[mEdge[0]] == sEdge[0])
                        {
                            overRideCase = 1;
                        }
                        else
                        if (pointMap[mEdge[1]] == sEdge[0])
                        {
                            overRideCase = 2;
                        }
                        else
                        {
                            FatalErrorIn("dynamicTopoFvMesh::collapseEdge()")
                                << "Coupled collapse failed." << nl
                                << "Master: " << mEdge << nl
                                << "Slave: " << sEdge << nl
                                << abort(FatalError);
                        }

                        break;

                    case 2:

                        if (pointMap[mEdge[1]] == sEdge[1])
                        {
                            overRideCase = 2;
                        }
                        else
                        if (pointMap[mEdge[0]] == sEdge[1])
                        {
                            overRideCase = 1;
                        }
                        else
                        {
                            FatalErrorIn("dynamicTopoFvMesh::collapseEdge()")
                                << "Coupled collapse failed." << nl
                                << "Master: " << mEdge << nl
                                << "Slave: " << sEdge << nl
                                << abort(FatalError);
                        }

                        break;

                    case 3:

                        overRideCase = 3;

                        break;
                }

                // Can the overRideCase be used for this edge?
                masterMap = collapseEdge(eIndex, overRideCase, true);

                // Master couldn't perform collapse.
                if (masterMap.type() <= 0)
                {
                    // Turn coupledModification back on before bailing out.
                    setCoupledModification();

                    return masterMap;
                }

                // Collapse the slave edge.
                collapseEdge(sIndex);
            }
            else
            {
                // Slave couldn't perform collapse.
                setCoupledModification();

                map.type() = 0;

                return map;
            }

            // Turn coupledModification back on.
            setCoupledModification();
        }
        else
        if (processorCoupledEdge(eIndex))
        {
            // Collapse edge on the patchSubMesh.

        }
    }

    // Hull variables
    bool found = false;
    label replaceIndex = -1, m = edgePoints_[eIndex].size();

    // Size up the hull lists
    labelList cellHull(m, -1);
    labelList faceHull(m, -1);
    labelList edgeHull(m, -1);
    labelListList ringEntities(4, labelList(m, -1));

    if (debug > 1)
    {
        Info << nl << nl
             << "Edge: " << eIndex
             << ": " << edges_[eIndex]
             << " is to be collapsed. " << endl;

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
    }

    // Construct a hull around this edge
    constructHull
    (
        eIndex,
        edgeHull,
        faceHull,
        cellHull,
        ringEntities
    );

    // Check whether points of the edge lies on a boundary
    FixedList<label, 2> nBoundCurves(0), checkPoints(-1);
    const FixedList<bool,2> edgeBoundary = checkEdgeBoundary(eIndex);

    // Configure the new point-position
    point newPoint = vector::zero;
    point oldPoint = vector::zero;

    // Decide which point to remove
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
        // If this is an interior edge with two boundary points.
        // Bail out for now. If proximity based refinement is
        // switched on, mesh may be sliced at this point.
        if (whichEdgePatch(eIndex) == -1)
        {
            return map;
        }

        // Check if either point lies on a bounding curve
        // Used to ensure that collapses happen towards bounding curves.
        // If the edge itself is on a bounding curve, collapse is valid.
        forAll(edges_[eIndex], pointI)
        {
            const labelList& pEdges = pointEdges_[edges_[eIndex][pointI]];

            forAll(pEdges, edgeI)
            {
                if (checkBoundingCurve(pEdges[edgeI]))
                {
                    nBoundCurves[pointI]++;
                }
            }
        }

        // Pick the point which is connected to more bounding curves
        if (nBoundCurves[0] > nBoundCurves[1])
        {
            collapseCase = 1;
        }
        else
        if (nBoundCurves[1] > nBoundCurves[0])
        {
            collapseCase = 2;
        }
        else
        {
            // Bounding edge: collapseEdge can collapse this edge
            collapseCase = 3;
        }

        // Override previous decision for rare cases.
        // Check if either point touches a hull boundary face, and retain that.
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

        // Check if we can continue...
        if (faceCheck[0] && !faceCheck[1])
        {
            if (collapseCase != 1)
            {
                return map;
            }
        }
        else
        if (!faceCheck[0] && faceCheck[1])
        {
            if (collapseCase != 2)
            {
                return map;
            }
        }
        else
        if (collapseCase != 3)
        {
            return map;
        }
    }
    else
    {
        // Looks like this is an interior edge.
        // Collapse case [3] by default
        collapseCase = 3;
    }

    // Perform an override if requested.
    if (overRideCase != -1)
    {
        collapseCase = overRideCase;
    }

    switch (collapseCase)
    {
        case 1:

            // Collapse to the first node
            replacePoint = edges_[eIndex][0];
            collapsePoint = edges_[eIndex][1];
            replaceEdgeIndex = 0;
            replaceFaceIndex = 1;
            removeEdgeIndex = 2;
            removeFaceIndex = 3;
            newPoint = points_[edges_[eIndex][0]];
            oldPoint = oldPoints_[edges_[eIndex][0]];
            checkPoints[0] = collapsePoint;

            break;

        case 2:

            // Collapse to the second node
            replacePoint = edges_[eIndex][1];
            collapsePoint = edges_[eIndex][0];
            removeEdgeIndex = 0;
            removeFaceIndex = 1;
            replaceEdgeIndex = 2;
            replaceFaceIndex = 3;
            newPoint = points_[edges_[eIndex][1]];
            oldPoint = oldPoints_[edges_[eIndex][1]];
            checkPoints[0] = collapsePoint;

            break;

        case 3:

            // Collapse to the mid-point
            // Uses connectivity of the second case.
            replacePoint = edges_[eIndex][1];
            collapsePoint = edges_[eIndex][0];
            removeEdgeIndex = 0;
            removeFaceIndex = 1;
            replaceEdgeIndex = 2;
            replaceFaceIndex = 3;
            newPoint = 0.5 *
            (
                points_[edges_[eIndex][0]]
              + points_[edges_[eIndex][1]]
            );
            oldPoint = oldPoints_[edges_[eIndex][1]];
            checkPoints[0] = collapsePoint;
            checkPoints[1] = replacePoint;

            break;

        default:

            // Don't think this will ever happen.
            FatalErrorIn("dynamicTopoFvMesh::collapseEdge()")
                << "Edge: " << eIndex << ": " << edges_[eIndex]
                << ". Couldn't decide on collapseCase."
                << abort(FatalError);

            break;
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

        Info << "Vertices: " << edgePoints_[eIndex] << endl;
        Info << "Edges: " << edgeHull << endl;
        Info << "Faces: " << faceHull << endl;
        Info << "Cells: " << cellHull << endl;
        Info << "replacePoint: " << replacePoint << endl;
        Info << "collapsePoint: " << collapsePoint << endl;
        Info << "checkPoints: " << checkPoints << endl;;
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

        const labelList& checkPointEdges = pointEdges_[checkPoints[pointI]];

        forAll(checkPointEdges, edgeI)
        {
            const labelList& eFaces = edgeFaces_[checkPointEdges[edgeI]];

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
                            oldPoint,
                            collapsePoint,
                            own,
                            cellsChecked,
                            forceOp
                        )
                    )
                    {
                        map.type() = 0;
                        return map;
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
                            oldPoint,
                            collapsePoint,
                            nei,
                            cellsChecked,
                            forceOp
                        )
                    )
                    {
                        map.type() = 0;
                        return map;
                    }
                }
            }
        }
    }

    // Are we only performing checks?
    if (checkOnly)
    {
        map.type() = collapseCase;
        return map;
    }

    // Write out VTK files prior to change
    if (debug > 3)
    {
        labelList vtkCells = cellsChecked.toc();

        writeVTK
        (
            Foam::name(eIndex)
          + '(' + Foam::name(edges_[eIndex][0])
          + ',' + Foam::name(edges_[eIndex][1]) + ')'
          + "_Collapse_0",
            vtkCells
        );
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

        const labelList& rmvEdgeFaces = edgeFaces_[edgeToRemove];

        // Replace edgePoints for all edges emanating from hullVertices
        // except ring-edges; those are sized-down later
        const labelList& hullPointEdges =
        (
            pointEdges_[edgePoints_[eIndex][indexI]]
        );

        forAll(hullPointEdges, edgeI)
        {
            if
            (
                findIndex
                (
                    edgePoints_[hullPointEdges[edgeI]],
                    collapsePoint
                ) != -1
             && findIndex
                (
                    edgePoints_[hullPointEdges[edgeI]],
                    replacePoint
                ) == -1
            )
            {
                replaceLabel
                (
                    collapsePoint,
                    replacePoint,
                    edgePoints_[hullPointEdges[edgeI]]
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
            const face& faceToCheck = faces_[rmvEdgeFaces[faceI]];

            if ((replaceIndex = faceToCheck.which(collapsePoint)) > -1)
            {
                if (debug > 2)
                {
                    Info << "Renumbering face: "
                         << rmvEdgeFaces[faceI] << ": "
                         << faceToCheck << endl;
                }

                // Renumber the face...
                faces_[rmvEdgeFaces[faceI]][replaceIndex] = replacePoint;
            }

            // Hull faces should be removed for the replacement edge
            if (rmvEdgeFaces[faceI] == faceHull[indexI])
            {
                sizeDownList
                (
                    faceHull[indexI],
                    edgeFaces_[replaceEdge]
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
                    edgeFaces_[replaceEdge]
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

                    iPtr_->setFlip(replaceFace);
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

                    iPtr_->setFlip(replaceFace);
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
                << " is an orphan, i.e, no owner cell." << nl
                << " Error occurred while collapsing edge: "
                << eIndex << " :: " << edges_[eIndex]
                << " Patch: " << whichEdgePatch(eIndex)
                << " edgeBoundary: " << edgeBoundary
                << " nBoundCurves: " << nBoundCurves
                << " collapseCase: " << collapseCase
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
                << " is being converted to a boundary." << nl
                << " Error occurred while collapsing edge: "
                << eIndex << " :: " << edges_[eIndex]
                << " Patch: " << whichEdgePatch(eIndex)
                << " edgeBoundary: " << edgeBoundary
                << " nBoundCurves: " << nBoundCurves
                << " collapseCase: " << collapseCase
                << abort(FatalError);
        }
    }

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

            // Remove the hull cell
            removeCell(cellToRemove);
        }

        // Remove the hull edge and associated edgeFaces
        removeEdge(edgeToRemove);

        // Remove the hull face
        removeFace(faceHull[indexI]);
    }

    // Loop through pointEdges for the collapsePoint,
    // and replace all occurrences with replacePoint.
    // Size-up pointEdges for the replacePoint as well.
    const labelList& pEdges = pointEdges_[collapsePoint];

    // While we're at it, compile a set of cells which need
    // their old-volumes to be modified. Once all faces are
    // renumbered, old-volumes can be safely computed.
    labelHashSet modCells;

    forAll(pEdges, edgeI)
    {
        // Renumber edges
        label edgeIndex = pEdges[edgeI];

        if (edgeIndex != eIndex)
        {
            if (debug > 2)
            {
                Info << "Renumbering [edge]: "
                     << edgeIndex << ": "
                     << edges_[edgeIndex] << endl;
            }

            if (edges_[edgeIndex][0] == collapsePoint)
            {
                edges_[edgeIndex][0] = replacePoint;

                sizeUpList
                (
                    edgeIndex,
                    pointEdges_[replacePoint]
                );
            }
            else
            if (edges_[edgeIndex][1] == collapsePoint)
            {
                edges_[edgeIndex][1] = replacePoint;

                sizeUpList
                (
                    edgeIndex,
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
            const labelList& eFaces = edgeFaces_[edgeIndex];

            forAll(eFaces, faceI)
            {
                // Check owner / neighbour for cells.
                label own = owner_[eFaces[faceI]];
                label nei = neighbour_[eFaces[faceI]];

                // Check owner cell
                if (!modCells.found(own))
                {
                    modCells.insert(own);
                }

                // Check neighbour cell
                if (!modCells.found(nei) && nei != -1)
                {
                    modCells.insert(nei);
                }

                const face& faceToCheck = faces_[eFaces[faceI]];

                if ((replaceIndex = faceToCheck.which(collapsePoint)) > -1)
                {
                    if (debug > 2)
                    {
                        Info << "Renumbering face: "
                             << eFaces[faceI] << ": "
                             << faceToCheck << endl;
                    }

                    // Renumber the face...
                    faces_[eFaces[faceI]][replaceIndex] = replacePoint;

                    // Look for an edge on this face that doesn't
                    // contain collapsePoint or replacePoint.
                    label rplIndex = -1;
                    const labelList& fEdges = faceEdges_[eFaces[faceI]];

                    forAll(fEdges, edgeI)
                    {
                        const edge& eCheck = edges_[fEdges[edgeI]];

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

    // Loop through accumulated cells and compute old-volumes.
    forAllConstIter(labelHashSet, modCells, cIter)
    {
        scalar modOldVol = tetVolume(cIter.key(), true);

        if (modOldVol < 0.0)
        {
            FatalErrorIn
            (
                "dynamicTopoFvMesh::collapseEdge()"
            )
                << "Negative old-volumes encountered." << nl
                << cIter.key() << ": " << modOldVol
                << abort(FatalError);
        }

        if (debug > 2)
        {
            Info << "Cell: " << cIter.key()
                 << " Old volume: " << modOldVol
                 << " New volume: " << tetVolume(cIter.key())
                 << endl;
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

    // Move to the new point
    points_[replacePoint] = newPoint;

    // Remove the collapse point
    removePoint(collapsePoint);

    // Write out VTK files after change
    if (debug > 3)
    {
        // Since cellsChecked is no longer used,
        // we'll use it for post-processing.
        forAll(cellHull, indexI)
        {
            if (cellsChecked.found(cellHull[indexI]))
            {
                cellsChecked.erase(cellHull[indexI]);
            }
        }

        labelList vtkCells = cellsChecked.toc();

        writeVTK
        (
            Foam::name(eIndex)
          + '(' + Foam::name(edges_[eIndex][0])
          + ',' + Foam::name(edges_[eIndex][1]) + ')'
          + "_Collapse_1",
            vtkCells
        );
    }

    // Remove the edge
    removeEdge(eIndex);

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    nCollapses_++;

    // Increment the number of modifications
    nModifications_++;

    // Return a succesful collapse
    map.type() = collapseCase;

    return map;
}

// Merge two triangular boundary faces
//  - If a separation distance exists, attempt to fill the space
//    with tetrahedra, but otherwise, perform an exact merge, removing one
//    of the faces.
void dynamicTopoFvMesh::mergeBoundaryFaces
(
    const label firstFace,
    const label secondFace
)
{
    if (debug > 2)
    {
        Info << "Merging faces: "
             << firstFace << " and "
             << secondFace << endl;
    }

    // Sanity check: Are these actually boundary faces?
    if (neighbour_[firstFace] != -1 || neighbour_[secondFace] != -1)
    {
        FatalErrorIn
        (
            "dynamicTopoFvMesh::mergeBoundaryFaces()"
        )
            << nl << " Faces: "
            << firstFace << " and " << secondFace
            << " are not on boundaries. "
            << abort(FatalError);
    }

    // Perform distance-based checks to determine corresponding points
    Map<label> mapPoints;
    const face& firstPolyFace = faces_[firstFace];
    const face& secondPolyFace = faces_[secondFace];

    bool matchGeomTol = true;

    forAll(firstPolyFace, pointI)
    {
        scalarList pointDistance(secondPolyFace.size(), 0.0);

        forAll(secondPolyFace, pointJ)
        {
            pointDistance[pointJ] =
            (
                magSqr
                (
                    points_[firstPolyFace[pointI]]
                  - points_[secondPolyFace[pointJ]]
                )
            );
        }

        bool matchedPoint = false;

        while (!matchedPoint)
        {
            label minIndex = findMin(pointDistance);

            // Check if this point was mapped for another point.
            bool foundPoint = false;

            forAllIter(Map<label>::iterator, mapPoints, pIter)
            {
                if (pIter() == minIndex)
                {
                    // Point was mapped before. Repeat search.
                    pointDistance[minIndex] = GREAT;
                    foundPoint = true;
                    break;
                }
            }

            if (!foundPoint)
            {
                // Good. Point wasn't matched before.
                // Does it satisfy the geometric match tolerance?
                if (pointDistance[minIndex] > gTol_)
                {
                    matchGeomTol = false;
                }

                mapPoints.insert(pointI, minIndex);

                matchedPoint = true;
                break;
            }
        }
    }

    // Sanity check: Were all points matched up?
    if (mapPoints.size() != 3)
    {
        FatalErrorIn
        (
            "dynamicTopoFvMesh::mergeBoundaryFaces()"
        )
            << nl << " Faces: "
            << firstFace << " and " << secondFace
            << " do not match up. "
            << abort(FatalError);
    }

    if (matchGeomTol)
    {
        // Obtain owners for both faces, and compare their labels
        face newFace;
        label removedFace = -1, retainedFace = -1;
        label newOwner = -1, newNeighbour = -1;

        if (owner_[firstFace] < owner_[secondFace])
        {
            // Retain the first face
            newFace = firstPolyFace;
            newOwner = owner_[firstFace];
            newNeighbour = owner_[secondFace];
        }
        else
        {
            // Retain the second face
            newFace = secondPolyFace;
            newOwner = owner_[secondFace];
            newNeighbour = owner_[firstFace];
        }

        const labelList& faceEdges = faceEdges_[retainedFace];

        // Replace cell with the new face label
        replaceLabel
        (
            removedFace,
            retainedFace,
            cells_[newNeighbour]
        );

        // Remove the boundary face
        removeFace(removedFace);

        // Insert a new interior face
        insertFace(-1, newFace, newOwner, newNeighbour);

        // Add the faceEdges entry
        faceEdges_.append(faceEdges);

        // Done with work here. Bail out.
        return;
    }

    // Looks like there's a separation distance.
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
