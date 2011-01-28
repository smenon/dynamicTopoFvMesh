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

#include "Stack.H"
#include "triFace.H"
#include "objectMap.H"
#include "changeMap.H"
#include "multiThreader.H"
#include "coupledPatchInfo.H"
#include "dynamicTopoFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Method to collapse a quad-face in 2D
// - Returns a changeMap with a type specifying:
//    -1: Collapse failed since max number of topo-changes was reached.
//     0: Collapse could not be performed.
//     1: Collapsed to first node.
//     2: Collapsed to second node.
//     3: Collapse to mid-point.
// - overRideCase is used to force a certain collapse configuration.
//    -1: Use this value to let collapseQuadFace decide a case.
//     1: Force collapse to first node.
//     2: Force collapse to second node.
//     3: Force collapse to mid-point.
// - checkOnly performs a feasibility check and returns without modifications.
const changeMap dynamicTopoFvMesh::collapseQuadFace
(
    const label fIndex,
    label overRideCase,
    bool checkOnly,
    bool forceOp
)
{
    // Figure out which thread this is...
    label tIndex = self();

    // Prepare the changeMaps
    changeMap map;
    List<changeMap> slaveMaps;
    bool collapsingSlave = false;

    if
    (
        (statistics_[0] > maxModifications_)
     && (maxModifications_ > -1)
    )
    {
        // Reached the max allowable topo-changes.
        stack(tIndex).clear();

        return map;
    }

    // Check if edgeRefinements are to be avoided on patch.
    if (!isSubMesh_)
    {
        if (lengthEstimator().checkRefinementPatch(whichPatch(fIndex)))
        {
            return map;
        }
    }

    // Sanity check: Is the index legitimate?
    if (fIndex < 0)
    {
        FatalErrorIn
        (
            "\n"
            "const changeMap "
            "dynamicTopoFvMesh::collapseQuadFace\n"
            "(\n"
            "    const label fIndex,\n"
            "    label overRideCase,\n"
            "    bool checkOnly,\n"
            "    bool forceOp\n"
            ")\n"
        )
            << " Invalid index: " << fIndex
            << abort(FatalError);
    }

    // Define the edges on the face to be collapsed
    FixedList<edge,4> checkEdge(edge(-1, -1));
    FixedList<label,4> checkEdgeIndex(-1);

    // Define checkEdges
    getCheckEdges(fIndex, (*this), map, checkEdge, checkEdgeIndex);

    // Determine the common vertices for the first and second edges
    label cv0 = checkEdge[1].commonVertex(checkEdge[0]);
    label cv1 = checkEdge[1].commonVertex(checkEdge[3]);
    label cv2 = checkEdge[2].commonVertex(checkEdge[0]);
    label cv3 = checkEdge[2].commonVertex(checkEdge[3]);

    // If coupled modification is set, and this is a
    // master face, collapse its slaves first.
    bool localCouple = false, procCouple = false;

    if (coupledModification_)
    {
        const face& fCheck = faces_[fIndex];

        const label faceEnum = coupleMap::FACE;
        const label pointEnum = coupleMap::POINT;

        // Is this a locally coupled edge (either master or slave)?
        if (locallyCoupledEntity(fIndex, true))
        {
            localCouple = true;
            procCouple = false;
        }
        else
        if (processorCoupledEntity(fIndex))
        {
            procCouple = true;
            localCouple = false;
        }

        if (localCouple && !procCouple)
        {
            // Determine the slave index.
            label sIndex = -1, pIndex = -1;

            forAll(patchCoupling_, patchI)
            {
                if (patchCoupling_(patchI))
                {
                    const coupleMap& cMap = patchCoupling_[patchI].patchMap();

                    if ((sIndex = cMap.findSlave(faceEnum, fIndex)) > -1)
                    {
                        pIndex = patchI;

                        break;
                    }

                    // The following bit happens only during the sliver
                    // exudation process, since slave edges are
                    // usually not added to the coupled edge-stack.
                    if ((sIndex = cMap.findMaster(faceEnum, fIndex)) > -1)
                    {
                        pIndex = patchI;

                        // Notice that we are collapsing a slave edge first.
                        collapsingSlave = true;

                        break;
                    }
                }
            }

            if (sIndex == -1)
            {
                FatalErrorIn
                (
                    "\n"
                    "const changeMap "
                    "dynamicTopoFvMesh::collapseQuadFace\n"
                    "(\n"
                    "    const label fIndex,\n"
                    "    label overRideCase,\n"
                    "    bool checkOnly,\n"
                    "    bool forceOp\n"
                    ")\n"
                )
                    << "Coupled maps were improperly specified." << nl
                    << " Slave index not found for: " << nl
                    << " Face: " << fIndex << ": " << fCheck << nl
                    << abort(FatalError);
            }
            else
            {
                // If we've found the slave, size up the list
                meshOps::sizeUpList
                (
                    changeMap(),
                    slaveMaps
                );

                // Save index and patch for posterity
                slaveMaps[0].index() = sIndex;
                slaveMaps[0].patchIndex() = pIndex;
            }

            if (debug > 1)
            {
                Pout<< nl << " >> Collapsing slave face: " << sIndex
                    << " for master face: " << fIndex << endl;
            }
        }
        else
        if (procCouple && !localCouple)
        {
            // Check slaves
            forAll(procIndices_, pI)
            {
                // Fetch reference to subMeshes
                const coupledPatchInfo& sendMesh = sendPatchMeshes_[pI];
                const coupledPatchInfo& recvMesh = recvPatchMeshes_[pI];

                const coupleMap& scMap = sendMesh.patchMap();
                const coupleMap& rcMap = recvMesh.patchMap();

                // If this face was sent to a lower-ranked
                // processor, skip it.
                if (procIndices_[pI] < Pstream::myProcNo())
                {
                    if (scMap.reverseEntityMap(faceEnum).found(fIndex))
                    {
                        if (debug > 3)
                        {
                            Pout<< "Face: " << fIndex
                                << "::" << fCheck
                                << " was sent to proc: "
                                << procIndices_[pI]
                                << ", so bailing out."
                                << endl;
                        }

                        return map;
                    }
                }

                label sIndex = -1;

                if ((sIndex = rcMap.findSlave(faceEnum, fIndex)) > -1)
                {
                    // Check if a lower-ranked processor is
                    // handling this face
                    if (procIndices_[pI] < Pstream::myProcNo())
                    {
                        if (debug > 3)
                        {
                            Pout<< "Face: " << fIndex
                                << "::" << fCheck
                                << " is handled by proc: "
                                << procIndices_[pI]
                                << ", so bailing out."
                                << endl;
                        }

                        return map;
                    }

                    label curIndex = slaveMaps.size();

                    // Size up the list
                    meshOps::sizeUpList
                    (
                        changeMap(),
                        slaveMaps
                    );

                    // Save index and patch for posterity
                    slaveMaps[curIndex].index() = sIndex;
                    slaveMaps[curIndex].patchIndex() = pI;

                    // Only one slave coupling is possible, so bail out
                    break;
                }
                else
                if
                (
                    (
                        rcMap.findSlave(pointEnum, cv0) > -1 &&
                        rcMap.findSlave(pointEnum, cv1) > -1
                    )
                 || (
                        rcMap.findSlave(pointEnum, cv2) > -1 &&
                        rcMap.findSlave(pointEnum, cv3) > -1
                    )
                )
                {
                    // An edge-only coupling exists.

                    // Check if a lower-ranked processor is
                    // handling this face
                    if (procIndices_[pI] < Pstream::myProcNo())
                    {
                         if (debug > 3)
                         {
                             Pout<< "Face edge on: " << fIndex
                                 << "::" << fCheck
                                 << " is handled by proc: "
                                 << procIndices_[pI]
                                 << ", so bailing out."
                                 << endl;
                         }

                         return map;
                    }

                    label p0 = rcMap.findSlave(pointEnum, cv0);
                    label p1 = rcMap.findSlave(pointEnum, cv1);

                    label p2 = rcMap.findSlave(pointEnum, cv2);
                    label p3 = rcMap.findSlave(pointEnum, cv3);

                    edge cEdge(-1, -1);
                    label edgeCouple = -1;

                    if ((p0 > -1 && p1 > -1) && (p2 == -1 && p3 == -1))
                    {
                        cEdge[0] = p0;
                        cEdge[1] = p1;

                        edgeCouple = 1;
                        sIndex = map.firstEdge();
                    }
                    else
                    if ((p0 == -1 && p1 == -1) && (p2 > -1 && p3 > -1))
                    {
                        cEdge[0] = p2;
                        cEdge[1] = p3;

                        edgeCouple = 2;
                        sIndex = map.secondEdge();
                    }

                    label curIndex = slaveMaps.size();

                    // Size up the list
                    meshOps::sizeUpList
                    (
                        changeMap(),
                        slaveMaps
                    );

                    // Unfortunately, since no edge maps are
                    // available in 2D, loop through boundary
                    // faces on subMesh and get the right edge
                    label seIndex = -1;

                    const dynamicTopoFvMesh& sMesh = recvMesh.subMesh();

                    for
                    (
                        label faceI = sMesh.nOldInternalFaces_;
                        faceI < sMesh.faceEdges_.size();
                        faceI++
                    )
                    {
                        const labelList& fEdges = sMesh.faceEdges_[faceI];

                        forAll(fEdges, edgeI)
                        {
                            if (sMesh.edges_[fEdges[edgeI]] == cEdge)
                            {
                                seIndex = fEdges[edgeI];
                                break;
                            }
                        }

                        if (seIndex > -1)
                        {
                            break;
                        }
                    }

                    if (seIndex == -1)
                    {
                        Pout<< "Face edge on: " << fIndex
                            << "::" << fCheck
                            << " sIndex: " << sIndex
                            << " could not be found."
                            << abort(FatalError);
                    }

                    // Save index and patch for posterity
                    //  - Negate the index to signify edge coupling
                    slaveMaps[curIndex].index() = -seIndex;
                    slaveMaps[curIndex].patchIndex() = pI;
                }
            }
        }
        else
        {
            // Something's wrong with coupling maps
            FatalErrorIn
            (
                "\n"
                "const changeMap "
                "dynamicTopoFvMesh::collapseQuadFace\n"
                "(\n"
                "    const label fIndex,\n"
                "    label overRideCase,\n"
                "    bool checkOnly,\n"
                "    bool forceOp\n"
                ")\n"
            )
                << "Coupled maps were improperly specified." << nl
                << " localCouple: " << localCouple << nl
                << " procCouple: " << procCouple << nl
                << " Face: " << fIndex << ": " << fCheck << nl
                << abort(FatalError);
        }

        // Temporarily turn off coupledModification
        unsetCoupledModification();

        // Test the master face for collapse, and decide on a case
        changeMap masterMap = collapseQuadFace(fIndex, -1, true, forceOp);

        // Turn it back on.
        setCoupledModification();

        // Master couldn't perform collapse
        if (masterMap.type() <= 0)
        {
            return masterMap;
        }

        // For edge-only coupling, define the points for checking
        List<FixedList<point, 2> > slaveMoveNewPoint(slaveMaps.size());
        List<FixedList<point, 2> > slaveMoveOldPoint(slaveMaps.size());

        // Now check each of the slaves for collapse feasibility
        forAll(slaveMaps, slaveI)
        {
            // Alias for convenience...
            changeMap& slaveMap = slaveMaps[slaveI];

            label slaveOverRide = -1;
            label sIndex = slaveMap.index();
            label pI = slaveMap.patchIndex();
            const coupleMap* cMapPtr = NULL;

            if (localCouple)
            {
                cMapPtr = &(patchCoupling_[pI].patchMap());
            }
            else
            if (procCouple)
            {
                const dynamicTopoFvMesh& sMesh =
                (
                    recvPatchMeshes_[pI].subMesh()
                );

                cMapPtr = &(recvPatchMeshes_[pI].patchMap());

                if (debug > 3)
                {
                    if (sIndex < 0)
                    {
                        Pout<< "Checking slave edge: " << mag(sIndex)
                            << "::" << sMesh.edges_[mag(sIndex)]
                            << " on proc: " << procIndices_[pI]
                            << " for master face: " << fIndex
                            << " using collapseCase: " << masterMap.type()
                            << endl;
                    }
                    else
                    {
                        Pout<< "Checking slave face: " << sIndex
                            << "::" << sMesh.faces_[sIndex]
                            << " on proc: " << procIndices_[pI]
                            << " for master face: " << fIndex
                            << " using collapseCase: " << masterMap.type()
                            << endl;
                    }
                }
            }

            // Define an overRide case for the slave
            FixedList<edge, 2> mEdge(edge(-1, -1)), sEdge(edge(-1, -1));

            if (collapsingSlave)
            {
                const Map<label>& rPointMap =
                (
                    cMapPtr->reverseEntityMap(pointEnum)
                );

                if (sIndex < 0)
                {

                }
                else
                {
                    mEdge[0][0] = rPointMap[cv0];
                    mEdge[0][1] = rPointMap[cv1];

                    mEdge[1][0] = rPointMap[cv2];
                    mEdge[1][1] = rPointMap[cv3];
                }
            }
            else
            {
                const Map<label>& pointMap =
                (
                    cMapPtr->entityMap(pointEnum)
                );

                if (sIndex < 0)
                {

                }
                else
                {
                    mEdge[0][0] = pointMap[cv0];
                    mEdge[0][1] = pointMap[cv1];

                    mEdge[1][0] = pointMap[cv2];
                    mEdge[1][1] = pointMap[cv3];
                }
            }

            // Determine checkEdges for the slave
            FixedList<edge,4> slaveCheckEdge(edge(-1, -1));
            FixedList<label,4> slaveCheckEdgeIndex(-1);

            if (localCouple)
            {
                getCheckEdges
                (
                    sIndex,
                    (*this),
                    slaveMap,
                    slaveCheckEdge,
                    slaveCheckEdgeIndex
                );

                sEdge[0] = edges_[slaveMap.firstEdge()];
                sEdge[1] = edges_[slaveMap.secondEdge()];
            }
            else
            if (procCouple)
            {
                const dynamicTopoFvMesh& sMesh =
                (
                    recvPatchMeshes_[pI].subMesh()
                );

                if (sIndex < 0)
                {
                    sEdge[0] = sMesh.edges_[mag(sIndex)];
                }
                else
                {
                    getCheckEdges
                    (
                        sIndex,
                        sMesh,
                        slaveMap,
                        slaveCheckEdge,
                        slaveCheckEdgeIndex
                    );

                    sEdge[0] = sMesh.edges_[slaveMap.firstEdge()];
                    sEdge[1] = sMesh.edges_[slaveMap.secondEdge()];
                }
            }

            // Compare edge orientations for edge-only coupling
            label compVal = -2;

            // Perform a topological comparison.
            switch (masterMap.type())
            {
                case 1:
                {
                    if (sIndex < 0)
                    {
                        slaveMoveNewPoint[slaveI][0] = points_[cv0];
                        slaveMoveNewPoint[slaveI][1] = points_[cv1];

                        slaveMoveOldPoint[slaveI][0] = oldPoints_[cv0];
                        slaveMoveOldPoint[slaveI][1] = oldPoints_[cv1];

                        // Check edge orientation
                        compVal = edge::compare(mEdge[0], sEdge[0]);
                    }
                    else
                    if (mEdge[0] == sEdge[0])
                    {
                        slaveOverRide = 1;
                    }
                    else
                    if (mEdge[1] == sEdge[0])
                    {
                        slaveOverRide = 2;
                    }
                    else
                    {
                        // Write out for for post-processing
                        writeVTK("mFace_" + Foam::name(fIndex), fIndex, 2);

                        FatalErrorIn
                        (
                            "\n"
                            "const changeMap "
                            "dynamicTopoFvMesh::collapseQuadFace\n"
                            "(\n"
                            "    const label fIndex,\n"
                            "    label overRideCase,\n"
                            "    bool checkOnly,\n"
                            "    bool forceOp\n"
                            ")\n"
                        )
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

                case 2:
                {
                    if (sIndex < 0)
                    {
                        slaveMoveNewPoint[slaveI][0] = points_[cv2];
                        slaveMoveNewPoint[slaveI][1] = points_[cv3];

                        slaveMoveOldPoint[slaveI][0] = oldPoints_[cv2];
                        slaveMoveOldPoint[slaveI][1] = oldPoints_[cv3];

                        // Check edge orientation
                        compVal = edge::compare(mEdge[0], sEdge[0]);
                    }
                    else
                    if (mEdge[1] == sEdge[1])
                    {
                        slaveOverRide = 2;
                    }
                    else
                    if (mEdge[0] == sEdge[1])
                    {
                        slaveOverRide = 1;
                    }
                    else
                    {
                        // Write out for for post-processing
                        writeVTK("mFace_" + Foam::name(fIndex), fIndex, 2);

                        FatalErrorIn
                        (
                            "\n"
                            "const changeMap "
                            "dynamicTopoFvMesh::collapseQuadFace\n"
                            "(\n"
                            "    const label fIndex,\n"
                            "    label overRideCase,\n"
                            "    bool checkOnly,\n"
                            "    bool forceOp\n"
                            ")\n"
                        )
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

                case 3:
                {
                    if (sIndex < 0)
                    {
                        slaveMoveNewPoint[slaveI][0] =
                        (
                            0.5 * (points_[cv0] + points_[cv2])
                        );

                        slaveMoveNewPoint[slaveI][1] =
                        (
                            0.5 * (points_[cv1] + points_[cv3])
                        );

                        slaveMoveNewPoint[slaveI][0] =
                        (
                            0.5 * (oldPoints_[cv0] + oldPoints_[cv2])
                        );

                        slaveMoveNewPoint[slaveI][1] =
                        (
                            0.5 * (oldPoints_[cv1] + oldPoints_[cv3])
                        );

                        // Check edge orientation
                        compVal = edge::compare(mEdge[0], sEdge[0]);
                    }
                    else
                    {
                        overRideCase = 3;
                    }

                    break;
                }
            }

            if (sIndex < 0)
            {
                // Swap components if necessary
                if (compVal == -1)
                {
                    Foam::Swap
                    (
                        slaveMoveNewPoint[slaveI][0],
                        slaveMoveNewPoint[slaveI][1]
                    );

                    Foam::Swap
                    (
                        slaveMoveOldPoint[slaveI][0],
                        slaveMoveOldPoint[slaveI][1]
                    );
                }
                else
                if (compVal != 1)
                {
                    FatalErrorIn
                    (
                        "\n"
                        "const changeMap "
                        "dynamicTopoFvMesh::collapseQuadFace\n"
                        "(\n"
                        "    const label fIndex,\n"
                        "    label overRideCase,\n"
                        "    bool checkOnly,\n"
                        "    bool forceOp\n"
                        ")\n"
                    )
                        << "Coupled topo-change for slave failed." << nl
                        << " Dissimilar edges: " << nl
                        << " mEdge: " << mEdge << nl
                        << " sEdge: " << sEdge << nl
                        << " masterMap type: " << masterMap.type() << nl
                        << abort(FatalError);
                }
            }

            // Temporarily turn off coupledModification
            unsetCoupledModification();

            // Test the slave face
            if (localCouple)
            {
                slaveMap =
                (
                    collapseQuadFace(sIndex, slaveOverRide, true, forceOp)
                );
            }
            else
            if (procCouple)
            {
                dynamicTopoFvMesh& sMesh = recvPatchMeshes_[pI].subMesh();

                if (sIndex < 0)
                {
                    // Edge-based coupling

                    // Build a hull of cells and tri-faces
                    // that are connected to the slave edge
                    labelList slaveHullCells;
                    labelList slaveHullTriFaces;

                    meshOps::constructPrismHull
                    (
                        mag(sIndex),
                        sMesh.faces_,
                        sMesh.cells_,
                        sMesh.owner_,
                        sMesh.neighbour_,
                        sMesh.edgeFaces_,
                        slaveHullTriFaces,
                        slaveHullCells
                    );

                    bool infeasible = false;

                    // Keep track of resulting cell quality,
                    // if collapse is indeed feasible
                    scalar slaveCollapseQuality(GREAT);

                    // Check whether the collapse is possible.
                    if
                    (
                        checkCollapse
                        (
                            sMesh,
                            slaveHullTriFaces,
                            FixedList<label, 2>(-1),
                            FixedList<label, 2>(-1),
                            sEdge[0],
                            slaveMoveNewPoint[slaveI],
                            slaveMoveOldPoint[slaveI],
                            slaveCollapseQuality,
                            false,
                            forceOp
                        )
                    )
                    {
                        infeasible = true;
                    }

                    if (infeasible)
                    {
                        slaveMap.type() = 0;
                    }
                    else
                    {
                        slaveMap.type() = 1;
                    }
                }
                else
                {
                    // Edge-based coupling
                    slaveMap =
                    (
                        sMesh.collapseQuadFace
                        (
                            sIndex,
                            slaveOverRide,
                            true,
                            forceOp
                        )
                    );
                }
            }

            // Turn it back on.
            setCoupledModification();

            if (slaveMap.type() <= 0)
            {
                // Slave couldn't perform collapse.
                map.type() = -2;

                return map;
            }

            // Save index and patch for posterity
            slaveMap.index() = sIndex;
            slaveMap.patchIndex() = pI;
        }

        // Next collapse each slave face (for real this time...)
        forAll(slaveMaps, slaveI)
        {
            // Alias for convenience...
            changeMap& slaveMap = slaveMaps[slaveI];

            label sIndex = slaveMap.index();
            label pI = slaveMap.patchIndex();
            label slaveOverRide = slaveMap.type();

            // Temporarily turn off coupledModification
            unsetCoupledModification();

            // Collapse the slave
            if (localCouple)
            {
                slaveMap =
                (
                    collapseQuadFace(sIndex, slaveOverRide, false, forceOp)
                );
            }
            else
            {
                dynamicTopoFvMesh& sMesh = recvPatchMeshes_[pI].subMesh();

                if (sIndex < 0)
                {
                    // Edge-based coupling

                    // Force operation to succeed
                    slaveMap.type() = 1;
                }
                else
                {
                    // Face-based coupling
                    slaveMap =
                    (
                        sMesh.collapseQuadFace
                        (
                            sIndex,
                            slaveOverRide,
                            false,
                            forceOp
                        )
                    );
                }
            }

            // Turn it back on.
            setCoupledModification();

            // The final operation has to succeed.
            if (slaveMap.type() <= 0)
            {
                FatalErrorIn
                (
                    "\n"
                    "const changeMap "
                    "dynamicTopoFvMesh::collapseQuadFace\n"
                    "(\n"
                    "    const label fIndex,\n"
                    "    label overRideCase,\n"
                    "    bool checkOnly,\n"
                    "    bool forceOp\n"
                    ")\n"
                )
                    << "Coupled topo-change for slave failed." << nl
                    << " Face: " << fIndex << ": " << fCheck << nl
                    << " Slave index: " << sIndex << nl
                    << " Patch index: " << pI << nl
                    << " Type: " << slaveMap.type() << nl
                    << abort(FatalError);
            }

            // Save index and patch for posterity
            slaveMap.index() = sIndex;
            slaveMap.patchIndex() = pI;
        }
    }

    // Build a hull of cells and tri-faces that are connected to each edge
    FixedList<labelList, 2> hullCells;
    FixedList<labelList, 2> hullTriFaces;

    meshOps::constructPrismHull
    (
        checkEdgeIndex[1],
        faces_,
        cells_,
        owner_,
        neighbour_,
        edgeFaces_,
        hullTriFaces[0],
        hullCells[0]
    );

    meshOps::constructPrismHull
    (
        checkEdgeIndex[2],
        faces_,
        cells_,
        owner_,
        neighbour_,
        edgeFaces_,
        hullTriFaces[1],
        hullCells[1]
    );

    // Determine the neighbouring cells
    label c0 = owner_[fIndex], c1 = neighbour_[fIndex];

    // Define variables for the prism-face calculation
    FixedList<label,2> c0BdyIndex(-1), c0IntIndex(-1);
    FixedList<label,2> c1BdyIndex(-1), c1IntIndex(-1);
    FixedList<face,2> c0BdyFace(face(3)), c0IntFace(face(4));
    FixedList<face,2> c1BdyFace(face(3)), c1IntFace(face(4));

    // Find the prism-faces
    meshOps::findPrismFaces
    (
        fIndex,
        c0,
        faces_,
        cells_,
        neighbour_,
        c0BdyFace,
        c0BdyIndex,
        c0IntFace,
        c0IntIndex
    );

    if (c1 != -1)
    {
        meshOps::findPrismFaces
        (
            fIndex,
            c1,
            faces_,
            cells_,
            neighbour_,
            c1BdyFace,
            c1BdyIndex,
            c1IntFace,
            c1IntIndex
        );
    }

    // Determine if either edge belongs to a boundary
    FixedList<bool, 2> nBoundCurves(false), edgeBoundary(false);

    edgeBoundary[0] = (whichEdgePatch(checkEdgeIndex[1]) > -1);
    edgeBoundary[1] = (whichEdgePatch(checkEdgeIndex[2]) > -1);

    // Decide on collapseCase
    label collapseCase = -1;

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
                    "\n"
                    "const changeMap "
                    "dynamicTopoFvMesh::collapseQuadFace\n"
                    "(\n"
                    "    const label fIndex,\n"
                    "    label overRideCase,\n"
                    "    bool checkOnly,\n"
                    "    bool forceOp\n"
                    ")\n"
                )   << "Collapsing an internal face that "
                    << "lies on two boundary patches. "
                    << "Face: " << fIndex << ": " << faces_[fIndex]
                    << endl;
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
        if (!nBoundCurves[0] && !nBoundCurves[1])
        {
            // No bounding curves. Collapse to mid-point.
            collapseCase = 3;
        }
        else
        {
            // Two bounding curves? This might cause pinching.
            // Bail out for now.
            return map;
        }
    }
    else
    {
        // Looks like this is an interior face.
        // Collapse case [3] by default
        collapseCase = 3;
    }

    // Perform an override if requested.
    if (overRideCase != -1)
    {
        collapseCase = overRideCase;
    }

    // Configure the new point-positions
    FixedList<bool, 2> check(false);
    FixedList<FixedList<label, 2>, 2> checkPoints(FixedList<label, 2>(-1));
    FixedList<point, 2> newPoint(vector::zero);
    FixedList<point, 2> oldPoint(vector::zero);

    // Replacement check points
    FixedList<label,2> original(-1), replacement(-1);

    switch (collapseCase)
    {
        case 1:
        {
            // Collapse to the first node
            original[0] = cv2; original[1] = cv3;
            replacement[0] = cv0; replacement[1] = cv1;

            newPoint[0] = points_[replacement[0]];
            newPoint[1] = points_[replacement[1]];

            oldPoint[0] = oldPoints_[replacement[0]];
            oldPoint[1] = oldPoints_[replacement[1]];

            // Define check-points
            check[1] = true;
            checkPoints[1][0] = original[0];
            checkPoints[1][1] = original[1];

            break;
        }

        case 2:
        {
            // Collapse to the second node
            original[0] = cv0; original[1] = cv1;
            replacement[0] = cv2; replacement[1] = cv3;

            newPoint[0] = points_[replacement[0]];
            newPoint[1] = points_[replacement[1]];

            oldPoint[0] = oldPoints_[replacement[0]];
            oldPoint[1] = oldPoints_[replacement[1]];

            // Define check-points
            check[0] = true;
            checkPoints[0][0] = original[0];
            checkPoints[0][1] = original[1];

            break;
        }

        case 3:
        {
            // Collapse to the mid-point
            original[0] = cv0; original[1] = cv1;
            replacement[0] = cv2; replacement[1] = cv3;

            // Define new point-positions
            newPoint[0] =
            (
                0.5 *
                (
                    points_[original[0]]
                  + points_[replacement[0]]
                )
            );

            newPoint[1] =
            (
                0.5 *
                (
                    points_[original[1]]
                  + points_[replacement[1]]
                )
            );

            // Define old point-positions
            oldPoint[0] =
            (
                0.5 *
                (
                    oldPoints_[original[0]]
                  + oldPoints_[replacement[0]]
                )
            );

            oldPoint[1] =
            (
                0.5 *
                (
                    oldPoints_[original[1]]
                  + oldPoints_[replacement[1]]
                )
            );

            // Define check-points
            check[0] = true;
            checkPoints[0][0] = original[0];
            checkPoints[0][1] = original[1];

            check[1] = true;
            checkPoints[1][0] = replacement[0];
            checkPoints[1][1] = replacement[1];

            break;
        }

        default:
        {
            // Don't think this will ever happen.
            FatalErrorIn
            (
                "\n"
                "const changeMap "
                "dynamicTopoFvMesh::collapseQuadFace\n"
                "(\n"
                "    const label fIndex,\n"
                "    label overRideCase,\n"
                "    bool checkOnly,\n"
                "    bool forceOp\n"
                ")\n"
            )
                << "Edge: " << fIndex << ": " << faces_[fIndex]
                << ". Couldn't decide on collapseCase."
                << abort(FatalError);

            break;
        }
    }

    // Keep track of resulting cell quality,
    // if collapse is indeed feasible
    scalar collapseQuality(GREAT);

    // Check collapsibility of cells around edges
    // with the re-configured point
    forAll(check, indexI)
    {
        if (!check[indexI])
        {
            continue;
        }

        // Check whether the collapse is possible.
        if
        (
            checkCollapse
            (
                (*this),
                hullTriFaces[indexI],
                c0BdyIndex,
                c1BdyIndex,
                checkPoints[indexI],
                newPoint,
                oldPoint,
                collapseQuality,
                (c1 != -1),
                forceOp
            )
        )
        {
            map.type() = 0;
            return map;
        }
    }

    // Are we only performing checks?
    if (checkOnly)
    {
        map.type() = collapseCase;

        if (debug > 2)
        {
            Pout<< "Face: " << fIndex
                << ":: " << faces_[fIndex] << nl
                << " collapseCase determined to be: " << collapseCase << nl
                << " Resulting quality: " << collapseQuality << nl
                << " edgeBoundary: " << edgeBoundary << nl
                << " nBoundCurves: " << nBoundCurves
                << endl;
        }

        return map;
    }

    if (debug > 1)
    {
        const labelList& fE = faceEdges_[fIndex];

        Pout<< nl << nl
            << "Face: " << fIndex << ": " << faces_[fIndex] << nl
            << "faceEdges: " << fE
            << " is to be collapsed. "
            << nl;

        label epIndex = whichPatch(fIndex);

        Pout<< "Patch: ";

        if (epIndex == -1)
        {
            Pout<< "Internal" << nl;
        }
        else
        {
            Pout<< boundaryMesh()[epIndex].name() << nl;
        }

        if (debug > 2)
        {
            Pout<< nl
                << "~~~~~~~~~~~~~~~~~~~~~~~~~" << nl
                << "Hulls before modification" << nl
                << "~~~~~~~~~~~~~~~~~~~~~~~~~" << nl;

            if (debug > 3)
            {
                Pout<< "Cell [0] Boundary faces: " << nl
                    << " Face: "<< c0BdyIndex[0]
                    << "::" << c0BdyFace[0] << nl
                    << " Face: "<< c0BdyIndex[1]
                    << "::" << c0BdyFace[1] << nl;

                Pout<< "Cell [0] Interior faces: " << nl
                    << " Face: "<< c0IntIndex[0]
                    << "::" << c0IntFace[0] << nl
                    << " Face: "<< c0IntIndex[1]
                    << "::" << c0IntFace[1] << nl;

                if (c1 != -1)
                {
                    Pout<< "Cell [1] Boundary faces: " << nl
                        << " Face: "<< c1BdyIndex[0]
                        << "::" << c1BdyFace[0] << nl
                        << " Face: "<< c1BdyIndex[1]
                        << "::" << c1BdyFace[1] << nl;

                    Pout<< "Cell [1] Interior faces: " << nl
                        << " Face: "<< c1IntIndex[0]
                        << "::" << c1IntFace[0] << nl
                        << " Face: "<< c1IntIndex[1]
                        << "::" << c1IntFace[1] << nl;
                }
            }

            Pout<< nl << "Cells belonging to first Edge Hull: "
                << hullCells[0] << nl;

            forAll(hullCells[0],cellI)
            {
                const cell& firstCurCell = cells_[hullCells[0][cellI]];

                Pout<< "Cell: " << hullCells[0][cellI]
                    << ": " << firstCurCell << nl;

                forAll(firstCurCell,faceI)
                {
                    Pout<< " Face: " << firstCurCell[faceI]
                        << " : " << faces_[firstCurCell[faceI]]
                        << " fE: " << faceEdges_[firstCurCell[faceI]]
                        << nl;

                    if (debug > 3)
                    {
                        const labelList& fE = faceEdges_[firstCurCell[faceI]];

                        forAll(fE, edgeI)
                        {
                            Pout<< "  Edge: " << fE[edgeI]
                                << " : " << edges_[fE[edgeI]]
                                << nl;
                        }
                    }
                }
            }

            const labelList& firstEdgeFaces = edgeFaces_[checkEdgeIndex[1]];

            Pout<< nl << "First Edge Face Hull: "
                << firstEdgeFaces << nl;

            forAll(firstEdgeFaces,indexI)
            {
                Pout<< " Face: " << firstEdgeFaces[indexI]
                    << " : " << faces_[firstEdgeFaces[indexI]]
                    << " fE: " << faceEdges_[firstEdgeFaces[indexI]]
                    << nl;

                if (debug > 3)
                {
                    const labelList& fE = faceEdges_[firstEdgeFaces[indexI]];

                    forAll(fE, edgeI)
                    {
                        Pout<< "  Edge: " << fE[edgeI]
                            << " : " << edges_[fE[edgeI]]
                            << nl;
                    }
                }
            }

            Pout<< nl << "Cells belonging to second Edge Hull: "
                << hullCells[1] << nl;

            forAll(hullCells[1], cellI)
            {
                const cell& secondCurCell = cells_[hullCells[1][cellI]];

                Pout<< "Cell: " << hullCells[1][cellI]
                    << ": " << secondCurCell << nl;

                forAll(secondCurCell, faceI)
                {
                    Pout<< " Face: " << secondCurCell[faceI]
                        << " : " << faces_[secondCurCell[faceI]]
                        << " fE: " << faceEdges_[secondCurCell[faceI]]
                        << nl;

                    if (debug > 3)
                    {
                        const labelList& fE = faceEdges_[secondCurCell[faceI]];

                        forAll(fE, edgeI)
                        {
                            Pout<< "  Edge: " << fE[edgeI]
                                << " : " << edges_[fE[edgeI]]
                                << nl;
                        }
                    }
                }
            }

            const labelList& secondEdgeFaces = edgeFaces_[checkEdgeIndex[2]];

            Pout<< nl << "Second Edge Face Hull: "
                << secondEdgeFaces << nl;

            forAll(secondEdgeFaces, indexI)
            {
                Pout<< " Face: " << secondEdgeFaces[indexI]
                    << " : " << faces_[secondEdgeFaces[indexI]]
                    << " fE: " << faceEdges_[secondEdgeFaces[indexI]]
                    << nl;

                if (debug > 3)
                {
                    const labelList& fE = faceEdges_[secondEdgeFaces[indexI]];

                    forAll(fE, edgeI)
                    {
                        Pout<< "  Edge: " << fE[edgeI]
                            << " : " << edges_[fE[edgeI]]
                            << nl;
                    }
                }
            }

            // Write out VTK files prior to change
            if (debug > 3)
            {
                labelHashSet vtkCells;

                forAll(hullCells[0], cellI)
                {
                    if (!vtkCells.found(hullCells[0][cellI]))
                    {
                        vtkCells.insert(hullCells[0][cellI]);
                    }
                }

                forAll(hullCells[1], cellI)
                {
                    if (!vtkCells.found(hullCells[1][cellI]))
                    {
                        vtkCells.insert(hullCells[1][cellI]);
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
    }

    // Edges / Quad-faces to throw or keep during collapse
    FixedList<label,2> ends(-1);
    FixedList<label,2> faceToKeep(-1), faceToThrow(-1);
    FixedList<label,4> edgeToKeep(-1), edgeToThrow(-1);

    // Maintain a list of modified faces for mapping
    DynamicList<label> modifiedFaces(10);

    // Case 2 & 3 use identical connectivity,
    // but different point locations
    if (collapseCase == 2 || collapseCase == 3)
    {
        const labelList& firstEdgeFaces = edgeFaces_[checkEdgeIndex[1]];

        // Collapse to the second node...
        forAll(firstEdgeFaces,faceI)
        {
            // Replace point indices on faces.
            meshOps::replaceLabel
            (
                cv0,
                cv2,
                faces_[firstEdgeFaces[faceI]]
            );

            meshOps::replaceLabel
            (
                cv1,
                cv3,
                faces_[firstEdgeFaces[faceI]]
            );

            // Add an entry for mapping
            if (findIndex(modifiedFaces, firstEdgeFaces[faceI]) == -1)
            {
                modifiedFaces.append(firstEdgeFaces[faceI]);
            }

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
        meshOps::findCommonEdge
        (
            c0BdyIndex[0],
            faceToKeep[0],
            faceEdges_,
            edgeToKeep[0]
        );

        meshOps::findCommonEdge
        (
            c0BdyIndex[1],
            faceToKeep[0],
            faceEdges_,
            edgeToKeep[1]
        );

        meshOps::findCommonEdge
        (
            c0BdyIndex[0],
            faceToThrow[0],
            faceEdges_,
            edgeToThrow[0]
        );

        meshOps::findCommonEdge
        (
            c0BdyIndex[1],
            faceToThrow[0],
            faceEdges_,
            edgeToThrow[1]
        );

        // Size down edgeFaces for the ends.
        meshOps::findCommonEdge
        (
            faceToThrow[0],
            faceToKeep[0],
            faceEdges_,
            ends[0]
        );

        meshOps::sizeDownList
        (
            faceToThrow[0],
            edgeFaces_[ends[0]]
        );

        if (c1 != -1)
        {
            meshOps::findCommonEdge
            (
                c1BdyIndex[0],
                faceToKeep[1],
                faceEdges_,
                edgeToKeep[2]
            );

            meshOps::findCommonEdge
            (
                c1BdyIndex[1],
                faceToKeep[1],
                faceEdges_,
                edgeToKeep[3]
            );

            meshOps::findCommonEdge
            (
                c1BdyIndex[0],
                faceToThrow[1],
                faceEdges_,
                edgeToThrow[2]
            );

            meshOps::findCommonEdge
            (
                c1BdyIndex[1],
                faceToThrow[1],
                faceEdges_,
                edgeToThrow[3]
            );

            // Size down edgeFaces for the ends.
            meshOps::findCommonEdge
            (
                faceToThrow[1],
                faceToKeep[1],
                faceEdges_,
                ends[1]
            );

            meshOps::sizeDownList
            (
                faceToThrow[1],
                edgeFaces_[ends[1]]
            );
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
                meshOps::sizeDownList
                (
                    origTriFace,
                    edgeFaces_[edgeToKeep[indexI]]
                );
            }
            else
            {
                meshOps::replaceLabel
                (
                    origTriFace,
                    retTriFace,
                    edgeFaces_[edgeToKeep[indexI]]
                );

                meshOps::replaceLabel
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
            meshOps::replaceLabel
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
                meshOps::sizeUpList
                (
                    firstEdgeFaces[faceI],
                    edgeFaces_[checkEdgeIndex[2]]
                );
            }
        }

        // Remove the current face from the replacement edge
        meshOps::sizeDownList
        (
            fIndex,
            edgeFaces_[checkEdgeIndex[2]]
        );

        // Replace point labels on all triangular boundary faces.
        forAll(hullCells[0],cellI)
        {
            const cell& cellToCheck = cells_[hullCells[0][cellI]];

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

        // Update map
        map.removeEdge(checkEdgeIndex[0]);
        map.removeEdge(checkEdgeIndex[1]);
        map.removeEdge(checkEdgeIndex[2]);

        forAll(edgeToThrow, indexI)
        {
            if (edgeToThrow[indexI] == -1)
            {
                continue;
            }

            removeEdge(edgeToThrow[indexI]);

            // Update map
            map.removeEdge(edgeToThrow[indexI]);
        }

        // Delete the two points...
        removePoint(cv0);
        removePoint(cv1);

        // Update map
        map.removePoint(cv0);
        map.removePoint(cv1);

        if (coupledModification_)
        {
            // Alias for convenience...
            const changeMap& slaveMap = slaveMaps[0];

            const coupleMap* cMapPtr = NULL;
            label pIndex = slaveMap.patchIndex();

            if (localCouple && !procCouple)
            {
                cMapPtr = &(patchCoupling_[pIndex].patchMap());
            }
            else
            if (procCouple && !localCouple)
            {
                coupledPatchInfo& recvMesh = recvPatchMeshes_[pIndex];
                cMapPtr = &(recvMesh.patchMap());
            }

            // Alias for convenience
            const coupleMap& cMap = *cMapPtr;

            // Remove the point entries.
            const label pointEnum = coupleMap::POINT;

            // Obtain references
            Map<label>& pointMap = cMap.entityMap(pointEnum);
            Map<label>& rPointMap = cMap.reverseEntityMap(pointEnum);

            if (collapsingSlave)
            {
                pointMap.erase(rPointMap[cv0]);
                rPointMap.erase(cv0);
            }
            else
            {
                rPointMap.erase(pointMap[cv0]);
                pointMap.erase(cv0);
            }

            if (collapsingSlave)
            {
                pointMap.erase(rPointMap[cv1]);
                rPointMap.erase(cv1);
            }
            else
            {
                rPointMap.erase(pointMap[cv1]);
                pointMap.erase(cv1);
            }

            // Remove the face entries
            const label faceEnum = coupleMap::FACE;

            // Obtain references
            Map<label>& faceMap = cMap.entityMap(faceEnum);
            Map<label>& rFaceMap = cMap.reverseEntityMap(faceEnum);

            if (collapsingSlave)
            {
                faceMap.erase(faceMap[fIndex]);
                rFaceMap.erase(fIndex);
            }
            else
            {
                rFaceMap.erase(faceMap[fIndex]);
                faceMap.erase(fIndex);
            }

            // Push operation into coupleMap
            switch (slaveMap.type())
            {
                case 1:
                {
                    cMap.pushOperation
                    (
                        slaveMap.index(),
                        coupleMap::COLLAPSE_FIRST
                    );

                    break;
                }

                case 2:
                {
                    cMap.pushOperation
                    (
                        slaveMap.index(),
                        coupleMap::COLLAPSE_SECOND
                    );

                    break;
                }

                case 3:
                {
                    cMap.pushOperation
                    (
                        slaveMap.index(),
                        coupleMap::COLLAPSE_MIDPOINT
                    );

                    break;
                }
            }
        }
    }
    else
    {
        const labelList& secondEdgeFaces = edgeFaces_[checkEdgeIndex[2]];

        // Collapse to the first node
        forAll(secondEdgeFaces,faceI)
        {
            // Replace point indices on faces.
            meshOps::replaceLabel
            (
                cv2,
                cv0,
                faces_[secondEdgeFaces[faceI]]
            );

            meshOps::replaceLabel
            (
                cv3,
                cv1,
                faces_[secondEdgeFaces[faceI]]
            );

            // Add an entry for mapping
            if (findIndex(modifiedFaces, secondEdgeFaces[faceI]) == -1)
            {
                modifiedFaces.append(secondEdgeFaces[faceI]);
            }

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
        meshOps::findCommonEdge
        (
            c0BdyIndex[0],
            faceToKeep[0],
            faceEdges_,
            edgeToKeep[0]
        );

        meshOps::findCommonEdge
        (
            c0BdyIndex[1],
            faceToKeep[0],
            faceEdges_,
            edgeToKeep[1]
        );

        meshOps::findCommonEdge
        (
            c0BdyIndex[0],
            faceToThrow[0],
            faceEdges_,
            edgeToThrow[0]
        );

        meshOps::findCommonEdge
        (
            c0BdyIndex[1],
            faceToThrow[0],
            faceEdges_,
            edgeToThrow[1]
        );

        // Size down edgeFaces for the ends.
        meshOps::findCommonEdge
        (
            faceToThrow[0],
            faceToKeep[0],
            faceEdges_,
            ends[0]
        );

        meshOps::sizeDownList
        (
            faceToThrow[0],
            edgeFaces_[ends[0]]
        );

        if (c1 != -1)
        {
            meshOps::findCommonEdge
            (
                c1BdyIndex[0],
                faceToKeep[1],
                faceEdges_,
                edgeToKeep[2]
            );

            meshOps::findCommonEdge
            (
                c1BdyIndex[1],
                faceToKeep[1],
                faceEdges_,
                edgeToKeep[3]
            );

            meshOps::findCommonEdge
            (
                c1BdyIndex[0],
                faceToThrow[1],
                faceEdges_,
                edgeToThrow[2]
            );

            meshOps::findCommonEdge
            (
                c1BdyIndex[1],
                faceToThrow[1],
                faceEdges_,
                edgeToThrow[3]
            );

            // Size down edgeFaces for the ends.
            meshOps::findCommonEdge
            (
                faceToThrow[1],
                faceToKeep[1],
                faceEdges_,
                ends[1]
            );

            meshOps::sizeDownList
            (
                faceToThrow[1],
                edgeFaces_[ends[1]]
            );
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
                meshOps::sizeDownList
                (
                    origTriFace,
                    edgeFaces_[edgeToKeep[indexI]]
                );
            }
            else
            {
                meshOps::replaceLabel
                (
                    origTriFace,
                    retTriFace,
                    edgeFaces_[edgeToKeep[indexI]]
                );

                meshOps::replaceLabel
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
            meshOps::replaceLabel
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
                meshOps::sizeUpList
                (
                    secondEdgeFaces[faceI],
                    edgeFaces_[checkEdgeIndex[1]]
                );
            }
        }

        // Remove the current face from the replacement edge
        meshOps::sizeDownList
        (
            fIndex,
            edgeFaces_[checkEdgeIndex[1]]
        );

        // Replace point labels on all triangular boundary faces.
        forAll(hullCells[1], cellI)
        {
            const cell& cellToCheck = cells_[hullCells[1][cellI]];

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

        // Update map
        map.removeEdge(checkEdgeIndex[0]);
        map.removeEdge(checkEdgeIndex[2]);
        map.removeEdge(checkEdgeIndex[3]);

        forAll(edgeToThrow, indexI)
        {
            if (edgeToThrow[indexI] == -1)
            {
                continue;
            }

            removeEdge(edgeToThrow[indexI]);

            // Update map
            map.removeEdge(edgeToThrow[indexI]);
        }

        // Delete the two points...
        removePoint(cv2);
        removePoint(cv3);

        // Update map
        map.removePoint(cv2);
        map.removePoint(cv3);

        if (coupledModification_)
        {
            // Alias for convenience...
            const changeMap& slaveMap = slaveMaps[0];

            const coupleMap* cMapPtr = NULL;
            label pIndex = slaveMap.patchIndex();

            if (localCouple && !procCouple)
            {
                cMapPtr = &(patchCoupling_[pIndex].patchMap());
            }
            else
            if (procCouple && !localCouple)
            {
                coupledPatchInfo& recvMesh = recvPatchMeshes_[pIndex];
                cMapPtr = &(recvMesh.patchMap());
            }

            // Alias for convenience
            const coupleMap& cMap = *cMapPtr;

            // Remove the point entries.
            const label pointEnum = coupleMap::POINT;

            // Obtain references
            Map<label>& pointMap = cMap.entityMap(pointEnum);
            Map<label>& rPointMap = cMap.reverseEntityMap(pointEnum);

            if (collapsingSlave)
            {
                pointMap.erase(rPointMap[cv2]);
                rPointMap.erase(cv2);
            }
            else
            {
                rPointMap.erase(pointMap[cv2]);
                pointMap.erase(cv2);
            }

            if (collapsingSlave)
            {
                pointMap.erase(rPointMap[cv3]);
                rPointMap.erase(cv3);
            }
            else
            {
                rPointMap.erase(pointMap[cv3]);
                pointMap.erase(cv3);
            }

            // Remove the face entries
            const label faceEnum = coupleMap::FACE;

            // Obtain references
            Map<label>& faceMap = cMap.entityMap(faceEnum);
            Map<label>& rFaceMap = cMap.reverseEntityMap(faceEnum);

            if (collapsingSlave)
            {
                faceMap.erase(faceMap[fIndex]);
                rFaceMap.erase(fIndex);
            }
            else
            {
                rFaceMap.erase(faceMap[fIndex]);
                faceMap.erase(fIndex);
            }

            // Push operation into coupleMap
            switch (slaveMap.type())
            {
                case 1:
                {
                    cMap.pushOperation
                    (
                        slaveMap.index(),
                        coupleMap::COLLAPSE_FIRST
                    );

                    break;
                }

                case 2:
                {
                    cMap.pushOperation
                    (
                        slaveMap.index(),
                        coupleMap::COLLAPSE_SECOND
                    );

                    break;
                }

                case 3:
                {
                    cMap.pushOperation
                    (
                        slaveMap.index(),
                        coupleMap::COLLAPSE_MIDPOINT
                    );

                    break;
                }
            }
        }
    }

    if (debug > 2)
    {
        Pout<< nl
            << "~~~~~~~~~~~~~~~~~~~~~~~~" << nl
            << "Hulls after modification" << nl
            << "~~~~~~~~~~~~~~~~~~~~~~~~" << nl;

        Pout<< nl << "Cells belonging to first Edge Hull: "
            << hullCells[0] << nl;

        forAll(hullCells[0], cellI)
        {
            const cell& firstCurCell = cells_[hullCells[0][cellI]];

            Pout<< "Cell: " << hullCells[0][cellI]
                << ": " << firstCurCell << nl;

            forAll(firstCurCell, faceI)
            {
                Pout<< " Face: " << firstCurCell[faceI]
                    << " : " << faces_[firstCurCell[faceI]]
                    << " fE: " << faceEdges_[firstCurCell[faceI]]
                    << nl;
            }
        }

        const labelList& firstEdgeFaces = edgeFaces_[checkEdgeIndex[1]];

        Pout<< nl << "First Edge Face Hull: " << firstEdgeFaces << nl;

        forAll(firstEdgeFaces, indexI)
        {
            Pout<< " Face: " << firstEdgeFaces[indexI]
                << " : " << faces_[firstEdgeFaces[indexI]]
                << " fE: " << faceEdges_[firstEdgeFaces[indexI]]
                << nl;
        }

        Pout<< nl << "Cells belonging to second Edge Hull: "
            << hullCells[1] << nl;

        forAll(hullCells[1], cellI)
        {
            const cell& secondCurCell = cells_[hullCells[1][cellI]];

            Pout<< "Cell: " << hullCells[1][cellI]
                << ": " << secondCurCell << nl;

            forAll(secondCurCell, faceI)
            {
                Pout<< " Face: " << secondCurCell[faceI]
                    << " : " << faces_[secondCurCell[faceI]]
                    << " fE: " << faceEdges_[secondCurCell[faceI]]
                    << nl;
            }
        }

        const labelList& secondEdgeFaces = edgeFaces_[checkEdgeIndex[2]];

        Pout<< nl << "Second Edge Face Hull: " << secondEdgeFaces << nl;

        forAll(secondEdgeFaces, indexI)
        {
            Pout<< " Face : " << secondEdgeFaces[indexI]
                << " : " << faces_[secondEdgeFaces[indexI]]
                << " fE: " << faceEdges_[secondEdgeFaces[indexI]]
                << nl;
        }

        Pout<< "Retained face: "
            << faceToKeep[0] << ": "
            << " owner: " << owner_[faceToKeep[0]]
            << " neighbour: " << neighbour_[faceToKeep[0]]
            << nl;

        Pout<< "Discarded face: "
            << faceToThrow[0] << ": "
            << " owner: " << owner_[faceToThrow[0]]
            << " neighbour: " << neighbour_[faceToThrow[0]]
            << nl;

        if (c1 != -1)
        {
            Pout<< "Retained face: "
                << faceToKeep[1] << ": "
                << " owner: " << owner_[faceToKeep[1]]
                << " neighbour: " << neighbour_[faceToKeep[1]]
                << nl;

            Pout<< "Discarded face: "
                << faceToThrow[1] << ": "
                << " owner: " << owner_[faceToThrow[1]]
                << " neighbour: " << neighbour_[faceToThrow[1]]
                << nl;
        }

        Pout<< endl;
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

                setFlip(faceToKeep[0]);
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

                    setFlip(faceToKeep[0]);
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

                setFlip(faceToKeep[0]);
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

                    setFlip(faceToKeep[1]);
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

                        setFlip(faceToKeep[1]);
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

                    setFlip(faceToKeep[1]);
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

        // Update map
        map.removeFace(faceToKeep[0]);
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
        label newFaceIndex =
        (
            insertFace
            (
                whichPatch(faceToThrow[0]),
                newFace,
                newOwn,
                -1
            )
        );

        // Add an entry for mapping
        modifiedFaces.append(newFaceIndex);

        // Update map
        map.addFace(newFaceIndex, labelList(1, faceToThrow[0]));

        // Add a faceEdges entry as well.
        // Edges don't have to change, since they're
        // all on the boundary anyway.
        faceEdges_.append(newFaceEdges);

        meshOps::replaceLabel
        (
            faceToKeep[0],
            newFaceIndex,
            cells_[newOwn]
        );

        // Correct edgeFaces with the new face label.
        forAll(newFaceEdges, edgeI)
        {
            meshOps::replaceLabel
            (
                faceToKeep[0],
                newFaceIndex,
                edgeFaces_[newFaceEdges[edgeI]]
            );
        }

        // Renumber the neighbour so that this face is removed correctly.
        neighbour_[faceToKeep[0]] = 0;
        removeFace(faceToKeep[0]);

        // Update map
        map.removeFace(faceToKeep[0]);
    }

    // Remove the unwanted faces in the cell(s) adjacent to this face,
    // and correct the cells that contain discarded faces
    const cell& cell_0 = cells_[c0];

    forAll(cell_0,faceI)
    {
        if (cell_0[faceI] != fIndex && cell_0[faceI] != faceToKeep[0])
        {
            removeFace(cell_0[faceI]);

            // Update map
            map.removeFace(cell_0[faceI]);
        }
    }

    if (cellCheck[0] != -1)
    {
        meshOps::replaceLabel
        (
            faceToThrow[0],
            faceToKeep[0],
            cells_[cellCheck[0]]
        );
    }

    // Remove the cell
    removeCell(c0);

    // Update map
    map.removeCell(c0);

    if (c1 == -1)
    {
        // Increment the surface-collapse counter
        statistics_[6]++;
    }
    else
    {
        // Remove orphaned faces
        if (owner_[faceToKeep[1]] == -1)
        {
            removeFace(faceToKeep[1]);

            // Update map
            map.removeFace(faceToKeep[1]);
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
            label newFaceIndex =
            (
                insertFace
                (
                    whichPatch(faceToThrow[1]),
                    newFace,
                    newOwn,
                    -1
                )
            );

            // Add an entry for mapping
            modifiedFaces.append(newFaceIndex);

            // Update map
            map.addFace(newFaceIndex, labelList(1, faceToThrow[1]));

            // Add a faceEdges entry as well.
            // Edges don't have to change, since they're
            // all on the boundary anyway.
            faceEdges_.append(newFaceEdges);

            meshOps::replaceLabel
            (
                faceToKeep[1],
                newFaceIndex,
                cells_[newOwn]
            );

            // Correct edgeFaces with the new face label.
            forAll(newFaceEdges, edgeI)
            {
                meshOps::replaceLabel
                (
                    faceToKeep[1],
                    newFaceIndex,
                    edgeFaces_[newFaceEdges[edgeI]]
                );
            }

            // Renumber the neighbour so that this face is removed correctly.
            neighbour_[faceToKeep[1]] = 0;
            removeFace(faceToKeep[1]);

            // Update map
            map.removeFace(faceToKeep[1]);
        }

        const cell& cell_1 = cells_[c1];

        forAll(cell_1, faceI)
        {
            if (cell_1[faceI] != fIndex && cell_1[faceI] != faceToKeep[1])
            {
                removeFace(cell_1[faceI]);

                // Update map
                map.removeFace(cell_1[faceI]);
            }
        }

        if (cellCheck[1] != -1)
        {
            meshOps::replaceLabel
            (
                faceToThrow[1],
                faceToKeep[1],
                cells_[cellCheck[1]]
            );
        }

        // Remove the cell
        removeCell(c1);

        // Update map
        map.removeCell(c1);
    }

    // Move old / new points
    oldPoints_[replacement[0]] = oldPoint[0];
    oldPoints_[replacement[1]] = oldPoint[1];

    points_[replacement[0]] = newPoint[0];
    points_[replacement[1]] = newPoint[1];

    // Finally remove the face
    removeFace(fIndex);

    // Update map
    map.removeFace(fIndex);

    // Write out VTK files after change
    if (debug > 3)
    {
        labelHashSet vtkCells;

        forAll(hullCells[0], cellI)
        {
            if (hullCells[0][cellI] == c0 || hullCells[0][cellI] == c1)
            {
                continue;
            }

            if (!vtkCells.found(hullCells[0][cellI]))
            {
                vtkCells.insert(hullCells[0][cellI]);
            }
        }

        forAll(hullCells[1], cellI)
        {
            if (hullCells[1][cellI] == c0 || hullCells[1][cellI] == c1)
            {
                continue;
            }

            if (!vtkCells.found(hullCells[1][cellI]))
            {
                vtkCells.insert(hullCells[1][cellI]);
            }
        }

        writeVTK
        (
            Foam::name(fIndex)
          + "_Collapse_1",
            vtkCells.toc()
        );
    }

    // Fill-in candidate mapping information
    labelList mC(2, -1);
    mC[0] = c0, mC[1] = c1;

    // Now that all old / new cells possess correct connectivity,
    // compute mapping information.
    forAll(hullCells, indexI)
    {
        forAll(hullCells[indexI], cellI)
        {
            label mcIndex = hullCells[indexI][cellI];

            // Skip collapsed cells
            if (mcIndex == c0 || mcIndex == c1)
            {
                continue;
            }

            // Set the mapping for this cell
            setCellMapping(mcIndex, mC);
        }
    }

    // Set face mapping information for modified faces
    forAll(modifiedFaces, faceI)
    {
        const label mfIndex = modifiedFaces[faceI];

        // Exclude deleted faces
        if (faces_[mfIndex].empty())
        {
            continue;
        }

        // Decide between default / weighted mapping
        // based on boundary information
        label fPatch = whichPatch(mfIndex);

        if (fPatch == -1)
        {
            setFaceMapping(mfIndex);
        }
        else
        {
            // Fill-in candidate mapping information
            labelList faceCandidates;

            const labelList& fEdges = faceEdges_[mfIndex];

            forAll(fEdges, edgeI)
            {
                if (whichEdgePatch(fEdges[edgeI]) == fPatch)
                {
                    // Loop through associated edgeFaces
                    const labelList& eFaces = edgeFaces_[fEdges[edgeI]];

                    forAll(eFaces, faceI)
                    {
                        if
                        (
                            (eFaces[faceI] != mfIndex) &&
                            (whichPatch(eFaces[faceI]) == fPatch)
                        )
                        {
                            faceCandidates.setSize
                            (
                                faceCandidates.size() + 1,
                                eFaces[faceI]
                            );
                        }
                    }
                }
            }

            // Set the mapping for this face
            setFaceMapping(mfIndex, faceCandidates);
        }
    }

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    statistics_[4]++;

    // Increment the number of modifications
    statistics_[0]++;

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
//     3: Collapsed to mid-point (default).
// - overRideCase is used to force a certain collapse configuration.
//    -1: Use this value to let collapseEdge decide a case.
//     1: Force collapse to first node.
//     2: Force collapse to second node.
//     3: Force collapse to mid-point.
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

    // For 2D meshes, perform face-collapse
    if (twoDMesh_)
    {
        return collapseQuadFace(eIndex, overRideCase, checkOnly);
    }

    // Figure out which thread this is...
    label tIndex = self();

    // Prepare the changeMaps
    changeMap map;
    List<changeMap> slaveMaps;
    bool collapsingSlave = false;

    if
    (
        (statistics_[0] > maxModifications_)
     && (maxModifications_ > -1)
    )
    {
        // Reached the max allowable topo-changes.
        stack(tIndex).clear();

        return map;
    }

    // Check if edgeRefinements are to be avoided on patch.
    if (!isSubMesh_)
    {
        if (lengthEstimator().checkRefinementPatch(whichEdgePatch(eIndex)))
        {
            return map;
        }
    }

    // Sanity check: Is the index legitimate?
    if (eIndex < 0)
    {
        FatalErrorIn
        (
            "\n"
            "const changeMap dynamicTopoFvMesh::collapseEdge\n"
            "(\n"
            "    const label eIndex,\n"
            "    label overRideCase,\n"
            "    bool checkOnly,\n"
            "    bool forceOp\n"
            ")\n"
        )
            << " Invalid index: " << eIndex
            << abort(FatalError);
    }

    // If coupled modification is set, and this is a
    // master edge, collapse its slaves as well.
    bool localCouple = false, procCouple = false;

    if (coupledModification_)
    {
        const edge& eCheck = edges_[eIndex];

        const label edgeEnum = coupleMap::EDGE;
        const label pointEnum = coupleMap::POINT;

        // Is this a locally coupled edge (either master or slave)?
        if (locallyCoupledEntity(eIndex, true))
        {
            localCouple = true;
            procCouple = false;
        }
        else
        if (processorCoupledEntity(eIndex))
        {
            procCouple = true;
            localCouple = false;
        }

        if (localCouple && !procCouple)
        {
            // Determine the slave index.
            label sIndex = -1, pIndex = -1;

            forAll(patchCoupling_, patchI)
            {
                if (patchCoupling_(patchI))
                {
                    const coupleMap& cMap = patchCoupling_[patchI].patchMap();

                    if ((sIndex = cMap.findSlave(edgeEnum, eIndex)) > -1)
                    {
                        pIndex = patchI;

                        break;
                    }

                    // The following bit happens only during the sliver
                    // exudation process, since slave edges are
                    // usually not added to the coupled edge-stack.
                    if ((sIndex = cMap.findMaster(edgeEnum, eIndex)) > -1)
                    {
                        pIndex = patchI;

                        // Notice that we are collapsing a slave edge first.
                        collapsingSlave = true;

                        break;
                    }
                }
            }

            if (sIndex == -1)
            {
                FatalErrorIn
                (
                    "\n"
                    "const changeMap dynamicTopoFvMesh::collapseEdge\n"
                    "(\n"
                    "    const label eIndex,\n"
                    "    label overRideCase,\n"
                    "    bool checkOnly,\n"
                    "    bool forceOp\n"
                    ")\n"
                )
                    << "Coupled maps were improperly specified." << nl
                    << " Slave index not found for: " << nl
                    << " Edge: " << eIndex << ": " << eCheck << nl
                    << abort(FatalError);
            }
            else
            {
                // If we've found the slave, size up the list
                meshOps::sizeUpList
                (
                    changeMap(),
                    slaveMaps
                );

                // Save index and patch for posterity
                slaveMaps[0].index() = sIndex;
                slaveMaps[0].patchIndex() = pIndex;
            }

            if (debug > 1)
            {
                Pout<< nl << " >> Collapsing slave edge: " << sIndex
                    << " for master edge: " << eIndex << endl;
            }
        }
        else
        if (procCouple && !localCouple)
        {
            // Check slaves
            forAll(procIndices_, pI)
            {
                // Fetch reference to subMeshes
                const coupledPatchInfo& sendMesh = sendPatchMeshes_[pI];
                const coupledPatchInfo& recvMesh = recvPatchMeshes_[pI];

                const coupleMap& scMap = sendMesh.patchMap();
                const coupleMap& rcMap = recvMesh.patchMap();

                // If this edge was sent to a lower-ranked
                // processor, skip it.
                if (procIndices_[pI] < Pstream::myProcNo())
                {
                    if (scMap.reverseEntityMap(edgeEnum).found(eIndex))
                    {
                        if (debug > 3)
                        {
                            Pout<< "Edge: " << eIndex
                                << "::" << eCheck
                                << " was sent to proc: "
                                << procIndices_[pI]
                                << ", so bailing out."
                                << endl;
                        }

                        return map;
                    }
                }

                label sIndex = -1;

                if ((sIndex = rcMap.findSlave(edgeEnum, eIndex)) > -1)
                {
                    // Check if a lower-ranked processor is
                    // handling this edge
                    if (procIndices_[pI] < Pstream::myProcNo())
                    {
                        if (debug > 3)
                        {
                            Pout<< "Edge: " << eIndex
                                << "::" << eCheck
                                << " is handled by proc: "
                                << procIndices_[pI]
                                << ", so bailing out."
                                << endl;
                        }

                        return map;
                    }

                    label curIndex = slaveMaps.size();

                    // Size up the list
                    meshOps::sizeUpList
                    (
                        changeMap(),
                        slaveMaps
                    );

                    // Save index and patch for posterity
                    slaveMaps[curIndex].index() = sIndex;
                    slaveMaps[curIndex].patchIndex() = pI;
                }
                else
                if
                (
                    (rcMap.findSlave(pointEnum, eCheck[0]) > -1) ||
                    (rcMap.findSlave(pointEnum, eCheck[1]) > -1)
                )
                {
                    // A point-only coupling exists.

                    // Check if a lower-ranked processor is
                    // handling this edge
                    if (procIndices_[pI] < Pstream::myProcNo())
                    {
                         if (debug > 3)
                         {
                             Pout<< "Edge point on: " << eIndex
                                 << "::" << eCheck
                                 << " is handled by proc: "
                                 << procIndices_[pI]
                                 << ", so bailing out."
                                 << endl;
                         }

                         return map;
                    }

                    label p0Index = rcMap.findSlave(pointEnum, eCheck[0]);
                    label p1Index = rcMap.findSlave(pointEnum, eCheck[1]);

                    if (p0Index > -1 && p1Index == -1)
                    {
                        sIndex = p0Index;
                    }
                    else
                    if (p0Index == -1 && p1Index > -1)
                    {
                        sIndex = p1Index;
                    }

                    label curIndex = slaveMaps.size();

                    // Size up the list
                    meshOps::sizeUpList
                    (
                        changeMap(),
                        slaveMaps
                    );

                    // Save index and patch for posterity
                    //  - Negate the index to signify point coupling
                    slaveMaps[curIndex].index() = -sIndex;
                    slaveMaps[curIndex].patchIndex() = pI;
                }
            }
        }
        else
        {
            // Something's wrong with coupling maps
            FatalErrorIn
            (
                "\n"
                "const changeMap dynamicTopoFvMesh::collapseEdge\n"
                "(\n"
                "    const label eIndex,\n"
                "    label overRideCase,\n"
                "    bool checkOnly,\n"
                "    bool forceOp\n"
                ")\n"
            )
                << "Coupled maps were improperly specified." << nl
                << " localCouple: " << localCouple << nl
                << " procCouple: " << procCouple << nl
                << " Edge: " << eIndex << ": " << eCheck << nl
                << abort(FatalError);
        }

        // Temporarily turn off coupledModification
        unsetCoupledModification();

        // Test the master edge for collapse, and decide on a case
        changeMap masterMap = collapseEdge(eIndex, -1, true, forceOp);

        // Turn it back on.
        setCoupledModification();

        // Master couldn't perform collapse
        if (masterMap.type() <= 0)
        {
            return masterMap;
        }

        // For point-only coupling, define the points for checking
        pointField slaveMoveNewPoint(slaveMaps.size(), vector::zero);
        pointField slaveMoveOldPoint(slaveMaps.size(), vector::zero);

        // Now check each of the slaves for collapse feasibility
        forAll(slaveMaps, slaveI)
        {
            // Alias for convenience...
            changeMap& slaveMap = slaveMaps[slaveI];

            label slaveOverRide = -1;
            label sIndex = slaveMap.index();
            label pI = slaveMap.patchIndex();
            const coupleMap* cMapPtr = NULL;

            edge mEdge(eCheck), sEdge(-1, -1);

            if (localCouple)
            {
                sEdge = edges_[sIndex];

                cMapPtr = &(patchCoupling_[pI].patchMap());
            }
            else
            if (procCouple)
            {
                const dynamicTopoFvMesh& sMesh =
                (
                    recvPatchMeshes_[pI].subMesh()
                );

                cMapPtr = &(recvPatchMeshes_[pI].patchMap());

                if (sIndex < 0)
                {
                    if (debug > 3)
                    {
                        Pout<< "Checking slave point: " << mag(sIndex)
                            << "::" << sMesh.points_[mag(sIndex)]
                            << " on proc: " << procIndices_[pI]
                            << " for master edge: " << eIndex
                            << " using collapseCase: " << masterMap.type()
                            << endl;
                    }
                }
                else
                {
                    sEdge = sMesh.edges_[sIndex];

                    if (debug > 3)
                    {
                        Pout<< "Checking slave edge: " << sIndex
                            << "::" << sMesh.edges_[sIndex]
                            << " on proc: " << procIndices_[pI]
                            << " for master edge: " << eIndex
                            << " using collapseCase: " << masterMap.type()
                            << endl;
                    }
                }
            }

            const Map<label>& pointMap = cMapPtr->entityMap(pointEnum);

            // Perform a topological comparison.
            switch (masterMap.type())
            {
                case 1:
                {
                    if (sIndex < 0)
                    {
                        slaveMoveNewPoint[slaveI] = points_[mEdge[0]];
                        slaveMoveOldPoint[slaveI] = oldPoints_[mEdge[0]];
                    }
                    else
                    if (pointMap[mEdge[0]] == sEdge[0])
                    {
                        slaveOverRide = 1;
                    }
                    else
                    if (pointMap[mEdge[1]] == sEdge[0])
                    {
                        slaveOverRide = 2;
                    }
                    else
                    {
                        // Write out for for post-processing
                        writeVTK("mEdge_" + Foam::name(eIndex), eIndex, 1);

                        FatalErrorIn
                        (
                            "\n"
                            "const changeMap dynamicTopoFvMesh"
                            "::collapseEdge\n"
                            "(\n"
                            "    const label eIndex,\n"
                            "    label overRideCase,\n"
                            "    bool checkOnly,\n"
                            "    bool forceOp\n"
                            ")\n"
                        )
                            << "Coupled collapse failed." << nl
                            << "Master: " << eIndex << " : " << mEdge << nl
                            << "Slave: " << sIndex << " : " << sEdge << nl
                            << abort(FatalError);
                    }

                    break;
                }

                case 2:
                {
                    if (sIndex < 0)
                    {
                        slaveMoveNewPoint[slaveI] = points_[mEdge[1]];
                        slaveMoveOldPoint[slaveI] = oldPoints_[mEdge[1]];
                    }
                    else
                    if (pointMap[mEdge[1]] == sEdge[1])
                    {
                        slaveOverRide = 2;
                    }
                    else
                    if (pointMap[mEdge[0]] == sEdge[1])
                    {
                        slaveOverRide = 1;
                    }
                    else
                    {
                        // Write out for for post-processing
                        writeVTK("mEdge_" + Foam::name(eIndex), eIndex, 1);

                        FatalErrorIn
                        (
                            "\n"
                            "const changeMap dynamicTopoFvMesh"
                            "::collapseEdge\n"
                            "(\n"
                            "    const label eIndex,\n"
                            "    label overRideCase,\n"
                            "    bool checkOnly,\n"
                            "    bool forceOp\n"
                            ")\n"
                        )
                            << "Coupled collapse failed." << nl
                            << "Master: " << eIndex << " : " << mEdge << nl
                            << "Slave: " << sIndex << " : " << sEdge << nl
                            << abort(FatalError);
                    }

                    break;
                }

                case 3:
                {
                    if (sIndex < 0)
                    {
                        slaveMoveNewPoint[slaveI] =
                        (
                            0.5 * (points_[mEdge[0]] + points_[mEdge[1]])
                        );

                        slaveMoveOldPoint[slaveI] =
                        (
                            0.5 * (oldPoints_[mEdge[0]] + oldPoints_[mEdge[1]])
                        );
                    }
                    else
                    {
                        slaveOverRide = 3;
                    }

                    break;
                }
            }

            // Temporarily turn off coupledModification
            unsetCoupledModification();

            // Test the slave edge
            if (localCouple)
            {
                slaveMap = collapseEdge(sIndex, slaveOverRide, true, forceOp);
            }
            else
            if (procCouple)
            {
                dynamicTopoFvMesh& sMesh = recvPatchMeshes_[pI].subMesh();

                if (sIndex < 0)
                {
                    // Point-based coupling
                    scalar slaveCollapseQuality(GREAT);
                    DynamicList<label> cellsChecked(10);

                    // Check cells connected to coupled point
                    const labelList& pEdges = sMesh.pointEdges_[mag(sIndex)];

                    bool infeasible = false;

                    forAll(pEdges, edgeI)
                    {
                        const labelList& eFaces =
                        (
                            sMesh.edgeFaces_[pEdges[edgeI]]
                        );

                        // Build a list of cells to check
                        forAll(eFaces, faceI)
                        {
                            label own = sMesh.owner_[eFaces[faceI]];
                            label nei = sMesh.neighbour_[eFaces[faceI]];

                            // Check owner cell
                            if (findIndex(cellsChecked, own) == -1)
                            {
                                // Check if point movement is feasible
                                if
                                (
                                    sMesh.checkCollapse
                                    (
                                        slaveMoveNewPoint[slaveI],
                                        slaveMoveOldPoint[slaveI],
                                        mag(sIndex),
                                        own,
                                        cellsChecked,
                                        slaveCollapseQuality,
                                        forceOp
                                    )
                                )
                                {
                                    infeasible = true;
                                    break;
                                }
                            }

                            if (nei == -1)
                            {
                                continue;
                            }

                            // Check neighbour cell
                            if (findIndex(cellsChecked, nei) == -1)
                            {
                                // Check if point movement is feasible
                                if
                                (
                                    sMesh.checkCollapse
                                    (
                                        slaveMoveNewPoint[slaveI],
                                        slaveMoveOldPoint[slaveI],
                                        mag(sIndex),
                                        nei,
                                        cellsChecked,
                                        slaveCollapseQuality,
                                        forceOp
                                    )
                                )
                                {
                                    infeasible = true;
                                    break;
                                }
                            }
                        }

                        if (infeasible)
                        {
                            break;
                        }
                    }

                    if (infeasible)
                    {
                        slaveMap.type() = 0;
                    }
                    else
                    {
                        slaveMap.type() = 1;
                    }
                }
                else
                {
                    // Edge-based coupling
                    slaveMap =
                    (
                        sMesh.collapseEdge
                        (
                            sIndex,
                            slaveOverRide,
                            true,
                            forceOp
                        )
                    );
                }
            }

            // Turn it back on.
            setCoupledModification();

            if (slaveMap.type() <= 0)
            {
                // Slave couldn't perform collapse.
                map.type() = -2;

                return map;
            }

            // Save index and patch for posterity
            slaveMap.index() = sIndex;
            slaveMap.patchIndex() = pI;
        }

        // Next collapse each slave edge (for real this time...)
        forAll(slaveMaps, slaveI)
        {
            // Alias for convenience...
            changeMap& slaveMap = slaveMaps[slaveI];

            label sIndex = slaveMap.index();
            label pI = slaveMap.patchIndex();
            label slaveOverRide = slaveMap.type();

            // Temporarily turn off coupledModification
            unsetCoupledModification();

            // Collapse the slave
            if (localCouple)
            {
                slaveMap = collapseEdge(sIndex, slaveOverRide, false, forceOp);
            }
            else
            {
                dynamicTopoFvMesh& sMesh = recvPatchMeshes_[pI].subMesh();

                if (sIndex < 0)
                {
                    // Point based coupling
                    const coupleMap& cMap = recvPatchMeshes_[pI].patchMap();

                    // Move points to new location,
                    // and update operation into coupleMap
                    sMesh.points_[mag(sIndex)] = slaveMoveNewPoint[slaveI];
                    sMesh.oldPoints_[mag(sIndex)] = slaveMoveOldPoint[slaveI];

                    cMap.pushOperation
                    (
                        mag(sIndex),
                        coupleMap::MOVE_POINT,
                        slaveMoveNewPoint[slaveI],
                        slaveMoveOldPoint[slaveI]
                    );

                    // Force operation to succeed
                    slaveMap.type() = 1;
                }
                else
                {
                    // Edge-based coupling
                    slaveMap =
                    (
                        sMesh.collapseEdge
                        (
                            sIndex,
                            slaveOverRide,
                            false,
                            forceOp
                        )
                    );
                }
            }

            // Turn it back on.
            setCoupledModification();

            // The final operation has to succeed.
            if (slaveMap.type() <= 0)
            {
                FatalErrorIn
                (
                    "\n"
                    "const changeMap dynamicTopoFvMesh::collapseEdge\n"
                    "(\n"
                    "    const label eIndex,\n"
                    "    label overRideCase,\n"
                    "    bool checkOnly,\n"
                    "    bool forceOp\n"
                    ")\n"
                )
                    << "Coupled topo-change for slave failed." << nl
                    << " Edge: " << eIndex << ": " << eCheck << nl
                    << " Slave index: " << sIndex << nl
                    << " Patch index: " << pI << nl
                    << " Type: " << slaveMap.type() << nl
                    << abort(FatalError);
            }

            // Save index and patch for posterity
            slaveMap.index() = sIndex;
            slaveMap.patchIndex() = pI;
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

    // Construct a hull around this edge
    meshOps::constructHull
    (
        eIndex,
        faces_,
        edges_,
        cells_,
        owner_,
        neighbour_,
        faceEdges_,
        edgeFaces_,
        edgePoints_,
        edgeHull,
        faceHull,
        cellHull,
        ringEntities
    );

    // Check whether points of the edge lies on a boundary
    FixedList<label, 2> nBoundCurves(0), checkPoints(-1);
    const FixedList<bool,2> edgeBoundary = checkEdgeBoundary(eIndex);

    // Decide on collapseCase
    label collapseCase = -1;

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

    // Configure the new point-position
    point newPoint = vector::zero;
    point oldPoint = vector::zero;

    label collapsePoint = -1, replacePoint = -1;

    switch (collapseCase)
    {
        case 1:
        {
            // Collapse to the first node
            replacePoint = edges_[eIndex][0];
            collapsePoint = edges_[eIndex][1];

            newPoint = points_[replacePoint];
            oldPoint = oldPoints_[replacePoint];

            checkPoints[0] = collapsePoint;

            break;
        }

        case 2:
        {
            // Collapse to the second node
            replacePoint = edges_[eIndex][1];
            collapsePoint = edges_[eIndex][0];

            newPoint = points_[replacePoint];
            oldPoint = oldPoints_[replacePoint];

            checkPoints[0] = collapsePoint;

            break;
        }

        case 3:
        {
            // Collapse to the mid-point
            replacePoint = edges_[eIndex][1];
            collapsePoint = edges_[eIndex][0];

            newPoint =
            (
                0.5 *
                (
                    points_[replacePoint]
                  + points_[collapsePoint]
                )
            );

            oldPoint =
            (
                0.5 *
                (
                    oldPoints_[replacePoint]
                  + oldPoints_[collapsePoint]
                )
            );

            checkPoints[0] = replacePoint;
            checkPoints[1] = collapsePoint;

            break;
        }

        default:
        {
            // Don't think this will ever happen.
            FatalErrorIn
            (
                "\n"
                "const changeMap dynamicTopoFvMesh::collapseEdge\n"
                "(\n"
                "    const label eIndex,\n"
                "    label overRideCase,\n"
                "    bool checkOnly,\n"
                "    bool forceOp\n"
                ")\n"
            )
                << "Edge: " << eIndex << ": " << edges_[eIndex]
                << ". Couldn't decide on collapseCase."
                << abort(FatalError);

            break;
        }
    }

    // Loop through edges and check for feasibility of collapse
    // Also, keep track of resulting cell quality,
    // if collapse is indeed feasible
    scalar collapseQuality(GREAT);
    DynamicList<label> cellsChecked(10);

    // Add all hull cells as 'checked',
    // and therefore, feasible
    forAll(cellHull, cellI)
    {
        if (cellHull[cellI] == -1)
        {
            continue;
        }

        cellsChecked.append(cellHull[cellI]);
    }

    // Check collapsibility of cells around edges
    // with the re-configured point
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
                if (findIndex(cellsChecked, own) == -1)
                {
                    // Check if a collapse is feasible
                    if
                    (
                        checkCollapse
                        (
                            newPoint,
                            oldPoint,
                            checkPoints[pointI],
                            own,
                            cellsChecked,
                            collapseQuality,
                            forceOp
                        )
                    )
                    {
                        map.type() = 0;
                        return map;
                    }
                }

                if (nei == -1)
                {
                    continue;
                }

                // Check neighbour cell
                if (findIndex(cellsChecked, nei) == -1)
                {
                    // Check if a collapse is feasible
                    if
                    (
                        checkCollapse
                        (
                            newPoint,
                            oldPoint,
                            checkPoints[pointI],
                            nei,
                            cellsChecked,
                            collapseQuality,
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

    // Add a map entry of the existing edge as 'added',
    // since it will be deleted at the end of this operation
    map.addEdge(eIndex, labelList(edges_[eIndex]));

    // Are we only performing checks?
    if (checkOnly)
    {
        // Return a succesful collapse possibility
        map.type() = collapseCase;

        // Make note of the removed point
        map.removePoint(collapsePoint);

        if (debug > 2)
        {
            Pout<< nl << "Edge: " << eIndex
                << ":: " << edges_[eIndex] << nl
                << " collapseCase determined to be: " << collapseCase << nl
                << " Resulting quality: " << collapseQuality << nl
                << " collapsePoint: " << collapsePoint << nl
                << " nBoundCurves: " << nBoundCurves << nl
                << endl;
        }

        return map;
    }

    // Define indices on the hull for removal / replacement
    label removeEdgeIndex = -1, replaceEdgeIndex = -1;
    label removeFaceIndex = -1, replaceFaceIndex = -1;

    if (replacePoint == edges_[eIndex][0])
    {
        replaceEdgeIndex = 0;
        replaceFaceIndex = 1;
        removeEdgeIndex = 2;
        removeFaceIndex = 3;
    }
    else
    if (replacePoint == edges_[eIndex][1])
    {
        removeEdgeIndex = 0;
        removeFaceIndex = 1;
        replaceEdgeIndex = 2;
        replaceFaceIndex = 3;
    }
    else
    {
        // Don't think this will ever happen.
        FatalErrorIn
        (
            "\n"
            "const changeMap dynamicTopoFvMesh::collapseEdge\n"
            "(\n"
            "    const label eIndex,\n"
            "    label overRideCase,\n"
            "    bool checkOnly,\n"
            "    bool forceOp\n"
            ")\n"
        )
            << "Edge: " << eIndex << ": " << edges_[eIndex]
            << ". Couldn't decide on removal / replacement indices."
            << abort(FatalError);
    }

    if (debug > 1)
    {
        Pout<< nl << nl
            << "Edge: " << eIndex
            << ": " << edges_[eIndex]
            << " is to be collapsed. " << nl;

        label epIndex = whichEdgePatch(eIndex);

        Pout<< "Patch: ";

        if (epIndex == -1)
        {
            Pout<< "Internal" << nl;
        }
        else
        {
            Pout<< boundaryMesh()[epIndex].name() << nl;
        }

        Pout<< " nBoundCurves: " << nBoundCurves << nl
            << " collapseCase: " << collapseCase << nl
            << " Resulting quality: " << collapseQuality << endl;

        if (debug > 2)
        {
            Pout<< " Vertices: " << edgePoints_[eIndex] << nl
                << " Edges: " << edgeHull << nl
                << " Faces: " << faceHull << nl
                << " Cells: " << cellHull << nl
                << " replacePoint: " << replacePoint << nl
                << " collapsePoint: " << collapsePoint << nl
                << " checkPoints: " << checkPoints << nl
                << " ringEntities (removed faces): " << endl;

            forAll(ringEntities[removeFaceIndex], faceI)
            {
                label fIndex = ringEntities[removeFaceIndex][faceI];

                if (fIndex == -1)
                {
                    continue;
                }

                Pout<< fIndex << ": " << faces_[fIndex] << nl;
            }

            Pout<< " ringEntities (removed edges): " << nl;

            forAll(ringEntities[removeEdgeIndex], edgeI)
            {
                label ieIndex = ringEntities[removeEdgeIndex][edgeI];

                if (ieIndex == -1)
                {
                    continue;
                }

                Pout<< ieIndex << ": " << edges_[ieIndex] << nl;
            }

            Pout<< " ringEntities (replacement faces): " << nl;

            forAll(ringEntities[replaceFaceIndex], faceI)
            {
                label fIndex = ringEntities[replaceFaceIndex][faceI];

                if (fIndex == -1)
                {
                    continue;
                }

                Pout<< fIndex << ": " << faces_[fIndex] << nl;
            }

            Pout<< " ringEntities (replacement edges): " << nl;

            forAll(ringEntities[replaceEdgeIndex], edgeI)
            {
                label ieIndex = ringEntities[replaceEdgeIndex][edgeI];

                if (ieIndex == -1)
                {
                    continue;
                }

                Pout<< ieIndex << ": " << edges_[ieIndex] << nl;
            }

            labelList& collapsePointEdges = pointEdges_[collapsePoint];

            Pout<< " pointEdges (collapsePoint): ";

            forAll(collapsePointEdges, edgeI)
            {
                Pout<< collapsePointEdges[edgeI] << " ";
            }

            Pout<< nl;

            // Write out VTK files prior to change
            if (debug > 3)
            {
                writeVTK
                (
                    Foam::name(eIndex)
                  + '(' + Foam::name(collapsePoint)
                  + ',' + Foam::name(replacePoint) + ')'
                  + "_Collapse_0",
                    cellsChecked
                );
            }
        }
    }

    // Prepare entities for coupledModification mapping
    FixedList<bool, 2> mappedFace(false);
    FixedList<bool, 2> mappedRmEdge(false);
    FixedList<bool, 2> mappedRpEdge(false);
    FixedList<label, 2> bFaces(-1), brmEdges(-1), brpEdges(-1);

    if (whichEdgePatch(eIndex) > -1)
    {
        // Update number of surface collapses, if necessary.
        statistics_[6]++;

        // Prepare coupled entities for mapping.
        if (coupledModification_)
        {
            label nBF = 0;

            // Identify boundary entities that need to be updated
            forAll(faceHull, indexI)
            {
                if (whichPatch(faceHull[indexI]) > -1)
                {
                    bFaces[nBF] = faceHull[indexI];
                    brmEdges[nBF] = ringEntities[removeEdgeIndex][indexI];
                    brpEdges[nBF] = ringEntities[replaceEdgeIndex][indexI];

                    // If this face isn't coupled, skip updates
                    if
                    (
                        (!locallyCoupledEntity(bFaces[nBF], true, false, true))
                     && (!processorCoupledEntity(bFaces[nBF], true))
                    )
                    {
                        mappedFace[nBF] = true;
                    }

                    // If edges aren't coupled, skip updates
                    if
                    (
                        (!locallyCoupledEntity(brmEdges[nBF], true, false))
                     && (!processorCoupledEntity(brmEdges[nBF]))
                    )
                    {
                        mappedRmEdge[nBF] = true;
                    }

                    if
                    (
                        (!locallyCoupledEntity(brpEdges[nBF], true, false))
                     && (!processorCoupledEntity(brpEdges[nBF]))
                    )
                    {
                        mappedRpEdge[nBF] = true;
                    }

                    nBF++;
                }
            }

            if (debug > 2)
            {
                Pout<< nl << "Coupled collapse map candidates: " << nl
                    << " bFaces: " << bFaces
                    << " mappedFace: " << mappedFace << nl
                    << " brmEdges: " << brmEdges
                    << " mappedRmEdge: " << mappedRmEdge << nl
                    << " brpEdges: " << brpEdges
                    << " mappedRpEdge: " << mappedRpEdge << nl
                    << endl;
            }
        }
    }

    // Maintain a list of modified faces for mapping
    DynamicList<label> modifiedFaces(10);

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
                meshOps::replaceLabel
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
            meshOps::replaceLabel
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
                // Renumber the face...
                faces_[rmvEdgeFaces[faceI]][replaceIndex] = replacePoint;

                // Add an entry for mapping
                if (findIndex(modifiedFaces, rmvEdgeFaces[faceI]) == -1)
                {
                    modifiedFaces.append(rmvEdgeFaces[faceI]);
                }
            }

            // Hull faces should be removed for the replacement edge
            if (rmvEdgeFaces[faceI] == faceHull[indexI])
            {
                meshOps::sizeDownList
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
                meshOps::sizeUpList
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
        meshOps::sizeDownList
        (
            faceToRemove,
            edgeFaces_[edgeHull[indexI]]
        );

        // Size down edgePoints for the ring edges
        meshOps::sizeDownList
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

                    setFlip(replaceFace);
                }
                else
                if
                (
                    (neighbour_[faceToRemove] == -1) &&
                    (neighbour_[replaceFace] != -1)
                )
                {
                    // This interior face would need to be converted
                    // to a boundary one, and flipped as well.
                    face newFace = faces_[replaceFace].reverseFace();
                    label newOwner = neighbour_[replaceFace];
                    label newNeighbour = neighbour_[faceToRemove];
                    labelList newFE = faceEdges_[replaceFace];

                    label newFaceIndex =
                    (
                        insertFace
                        (
                            whichPatch(faceToRemove),
                            newFace,
                            newOwner,
                            newNeighbour
                        )
                    );

                    // Set this face aside for mapping
                    modifiedFaces.append(newFaceIndex);

                    // Update map.
                    map.addFace(newFaceIndex, labelList(1, faceToRemove));

                    // Ensure that all edges of this face are
                    // on the boundary.
                    forAll(newFE, edgeI)
                    {
                        if (whichEdgePatch(newFE[edgeI]) == -1)
                        {
                            edge newEdge = edges_[newFE[edgeI]];
                            labelList newEF = edgeFaces_[newFE[edgeI]];
                            labelList newEP = edgePoints_[newFE[edgeI]];

                            // Need patch information for the new edge.
                            // Find the corresponding edge in ringEntities.
                            // Note that hullEdges doesn't need to be checked,
                            // since they are common to both faces.
                            label i =
                            (
                                findIndex
                                (
                                    ringEntities[replaceEdgeIndex],
                                    newFE[edgeI]
                                )
                            );

                            label repIndex =
                            (
                                whichEdgePatch
                                (
                                    ringEntities[removeEdgeIndex][i]
                                )
                            );

                            // Insert the new edge
                            label newEdgeIndex =
                            (
                                insertEdge
                                (
                                    repIndex,
                                    newEdge,
                                    newEF,
                                    newEP
                                )
                            );

                            // Replace faceEdges for all
                            // connected faces.
                            forAll(newEF, faceI)
                            {
                                meshOps::replaceLabel
                                (
                                    newFE[edgeI],
                                    newEdgeIndex,
                                    faceEdges_[newEF[faceI]]
                                );
                            }

                            // Remove the edge
                            removeEdge(newFE[edgeI]);

                            // Update map
                            map.removeEdge(newFE[edgeI]);

                            map.addEdge
                            (
                                newEdgeIndex,
                                labelList(1, newFE[edgeI])
                            );

                            // Replace faceEdges with new edge index
                            newFE[edgeI] = newEdgeIndex;

                            // Modify ringEntities
                            ringEntities[replaceEdgeIndex][i] = newEdgeIndex;
                        }
                    }

                    // Add the new faceEdges
                    faceEdges_.append(newFE);

                    // Replace edgeFaces with the new face index
                    const labelList& newFEdges = faceEdges_[newFaceIndex];

                    forAll(newFEdges, edgeI)
                    {
                        meshOps::replaceLabel
                        (
                            replaceFace,
                            newFaceIndex,
                            edgeFaces_[newFEdges[edgeI]]
                        );
                    }

                    // Remove the face
                    removeFace(replaceFace);

                    // Update map
                    map.removeFace(replaceFace);

                    // Replace label for the new owner
                    meshOps::replaceLabel
                    (
                        replaceFace,
                        newFaceIndex,
                        cells_[newOwner]
                    );

                    // Modify ringEntities and replaceFace
                    replaceFace = newFaceIndex;
                    ringEntities[replaceFaceIndex][indexI] = newFaceIndex;
                }
                else
                if
                (
                    (neighbour_[faceToRemove] == -1) &&
                    (neighbour_[replaceFace] == -1)
                )
                {
                    // Wierd overhanging cell. Since replaceFace
                    // would be an orphan if this continued, remove
                    // the face entirely.
                    labelList rmFE = faceEdges_[replaceFace];

                    forAll(rmFE, edgeI)
                    {
                        if
                        (
                            (edgeFaces_[rmFE[edgeI]].size() == 1) &&
                            (edgeFaces_[rmFE[edgeI]][0] == replaceFace)
                        )
                        {
                            // This edge has to be removed entirely.
                            removeEdge(rmFE[edgeI]);

                            // Update map
                            map.removeEdge(rmFE[edgeI]);

                            label i =
                            (
                                findIndex
                                (
                                    ringEntities[replaceEdgeIndex],
                                    rmFE[edgeI]
                                )
                            );

                            // Modify ringEntities
                            ringEntities[replaceEdgeIndex][i] = -1;
                        }
                        else
                        {
                            // Size-down edgeFaces
                            meshOps::sizeDownList
                            (
                                replaceFace,
                                edgeFaces_[rmFE[edgeI]]
                            );
                        }
                    }

                    // Remove the face
                    removeFace(replaceFace);

                    // Update map
                    map.removeFace(replaceFace);

                    // Modify ringEntities and replaceFace
                    replaceFace = -1;
                    ringEntities[replaceFaceIndex][indexI] = -1;
                }
                else
                {
                    // Keep orientation intact, and update the owner
                    owner_[replaceFace] = neighbour_[faceToRemove];
                }
            }
            else
            if (neighbour_[faceToRemove] == -1)
            {
                // This interior face would need to be converted
                // to a boundary one, but with orientation intact.
                face newFace = faces_[replaceFace];
                label newOwner = owner_[replaceFace];
                label newNeighbour = neighbour_[faceToRemove];
                labelList newFE = faceEdges_[replaceFace];

                label newFaceIndex =
                (
                    insertFace
                    (
                        whichPatch(faceToRemove),
                        newFace,
                        newOwner,
                        newNeighbour
                    )
                );

                // Set this face aside for mapping
                modifiedFaces.append(newFaceIndex);

                // Update map
                map.addFace(newFaceIndex, labelList(1, faceToRemove));

                // Ensure that all edges of this face are
                // on the boundary.
                forAll(newFE, edgeI)
                {
                    if (whichEdgePatch(newFE[edgeI]) == -1)
                    {
                        edge newEdge = edges_[newFE[edgeI]];
                        labelList newEF = edgeFaces_[newFE[edgeI]];
                        labelList newEP = edgePoints_[newFE[edgeI]];

                        // Need patch information for the new edge.
                        // Find the corresponding edge in ringEntities.
                        // Note that hullEdges doesn't need to be checked,
                        // since they are common to both faces.
                        label i =
                        (
                            findIndex
                            (
                                ringEntities[replaceEdgeIndex],
                                newFE[edgeI]
                            )
                        );

                        label repIndex =
                        (
                            whichEdgePatch
                            (
                                ringEntities[removeEdgeIndex][i]
                            )
                        );

                        // Insert the new edge
                        label newEdgeIndex =
                        (
                            insertEdge
                            (
                                repIndex,
                                newEdge,
                                newEF,
                                newEP
                            )
                        );

                        // Replace faceEdges for all
                        // connected faces.
                        forAll(newEF, faceI)
                        {
                            meshOps::replaceLabel
                            (
                                newFE[edgeI],
                                newEdgeIndex,
                                faceEdges_[newEF[faceI]]
                            );
                        }

                        // Remove the edge
                        removeEdge(newFE[edgeI]);

                        // Update map
                        map.removeEdge(newFE[edgeI]);

                        map.addEdge
                        (
                            newEdgeIndex,
                            labelList(1, newFE[edgeI])
                        );

                        // Replace faceEdges with new edge index
                        newFE[edgeI] = newEdgeIndex;

                        // Modify ringEntities
                        ringEntities[replaceEdgeIndex][i] = newEdgeIndex;
                    }
                }

                // Add the new faceEdges
                faceEdges_.append(newFE);

                // Replace edgeFaces with the new face index
                const labelList& newFEdges = faceEdges_[newFaceIndex];

                forAll(newFEdges, edgeI)
                {
                    meshOps::replaceLabel
                    (
                        replaceFace,
                        newFaceIndex,
                        edgeFaces_[newFEdges[edgeI]]
                    );
                }

                // Remove the face
                removeFace(replaceFace);

                // Update map
                map.removeFace(replaceFace);

                // Replace label for the new owner
                meshOps::replaceLabel
                (
                    replaceFace,
                    newFaceIndex,
                    cells_[newOwner]
                );

                // Modify ringEntities and replaceFace
                replaceFace = newFaceIndex;
                ringEntities[replaceFaceIndex][indexI] = newFaceIndex;
            }
            else
            {
                // Keep orientation intact, and update the neighbour
                neighbour_[replaceFace] = neighbour_[faceToRemove];
            }

            // Update the cell
            if (neighbour_[faceToRemove] != -1)
            {
                meshOps::replaceLabel
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

                    setFlip(replaceFace);
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
            meshOps::replaceLabel
            (
                faceToRemove,
                replaceFace,
                cells_[owner_[faceToRemove]]
            );
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

            // Update map
            map.removeFace(faceToRemove);

            // Remove the hull cell
            removeCell(cellToRemove);

            // Update map
            map.removeCell(cellToRemove);
        }

        // Remove the hull edge and associated edgeFaces
        removeEdge(edgeToRemove);

        // Update map
        map.removeEdge(edgeToRemove);

        // Remove the hull face
        removeFace(faceHull[indexI]);

        // Update map
        map.removeFace(faceHull[indexI]);
    }

    // Loop through pointEdges for the collapsePoint,
    // and replace all occurrences with replacePoint.
    // Size-up pointEdges for the replacePoint as well.
    const labelList& pEdges = pointEdges_[collapsePoint];

    forAll(pEdges, edgeI)
    {
        // Renumber edges
        label edgeIndex = pEdges[edgeI];

        if (edgeIndex != eIndex)
        {
            if (edges_[edgeIndex][0] == collapsePoint)
            {
                edges_[edgeIndex][0] = replacePoint;

                meshOps::sizeUpList
                (
                    edgeIndex,
                    pointEdges_[replacePoint]
                );
            }
            else
            if (edges_[edgeIndex][1] == collapsePoint)
            {
                edges_[edgeIndex][1] = replacePoint;

                meshOps::sizeUpList
                (
                    edgeIndex,
                    pointEdges_[replacePoint]
                );
            }
            else
            {
                // Looks like pointEdges is inconsistent
                FatalErrorIn
                (
                    "\n"
                    "const changeMap dynamicTopoFvMesh::collapseEdge\n"
                    "(\n"
                    "    const label eIndex,\n"
                    "    label overRideCase,\n"
                    "    bool checkOnly,\n"
                    "    bool forceOp\n"
                    ")\n"
                )
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
                const face& faceToCheck = faces_[eFaces[faceI]];

                if ((replaceIndex = faceToCheck.which(collapsePoint)) > -1)
                {
                    // Renumber the face...
                    faces_[eFaces[faceI]][replaceIndex] = replacePoint;

                    // Set this face aside for mapping
                    if (findIndex(modifiedFaces, eFaces[faceI]) == -1)
                    {
                        modifiedFaces.append(eFaces[faceI]);
                    }

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
                    meshOps::replaceLabel
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
        // If the ring edge was removed, don't bother.
        if (ringEntities[replaceEdgeIndex][edgeI] == -1)
        {
            continue;
        }

        buildEdgePoints(ringEntities[replaceEdgeIndex][edgeI]);
    }

    // Move old / new points
    oldPoints_[replacePoint] = oldPoint;
    points_[replacePoint] = newPoint;

    // Remove the collapse point
    removePoint(collapsePoint);

    // Update map
    map.removePoint(collapsePoint);

    // Remove the edge
    removeEdge(eIndex);

    // Update map
    map.removeEdge(eIndex);

    // For cell-mapping, exclude all hull-cells
    forAll(cellsChecked, indexI)
    {
        if (cells_[cellsChecked[indexI]].empty())
        {
            cellsChecked[indexI] = -1;
        }
        else
        if (findIndex(cellHull, cellsChecked[indexI]) > 0)
        {
            cellsChecked[indexI] = -1;
        }
    }

    // Write out VTK files after change
    if (debug > 3)
    {
        writeVTK
        (
            Foam::name(eIndex)
          + '(' + Foam::name(collapsePoint)
          + ',' + Foam::name(replacePoint) + ')'
          + "_Collapse_1",
            cellsChecked
        );
    }

    // Now that all old / new cells possess correct connectivity,
    // compute mapping information.
    forAll(cellsChecked, cellI)
    {
        if (cellsChecked[cellI] < 0)
        {
            continue;
        }

        // Fill-in candidate mapping information
        labelList mC(1, cellsChecked[cellI]);

        // Set the mapping for this cell
        setCellMapping(cellsChecked[cellI], mC);
    }

    // Set face mapping information for modified faces
    forAll(modifiedFaces, faceI)
    {
        const label mfIndex = modifiedFaces[faceI];

        // Exclude deleted faces
        if (faces_[mfIndex].empty())
        {
            continue;
        }

        // Decide between default / weighted mapping
        // based on boundary information
        label fPatch = whichPatch(mfIndex);

        if (fPatch == -1)
        {
            setFaceMapping(mfIndex);
        }
        else
        {
            // Fill-in candidate mapping information
            labelList faceCandidates;

            const labelList& fEdges = faceEdges_[mfIndex];

            forAll(fEdges, edgeI)
            {
                if (whichEdgePatch(fEdges[edgeI]) == fPatch)
                {
                    // Loop through associated edgeFaces
                    const labelList& eFaces = edgeFaces_[fEdges[edgeI]];

                    forAll(eFaces, faceI)
                    {
                        if
                        (
                            (eFaces[faceI] != mfIndex) &&
                            (whichPatch(eFaces[faceI]) == fPatch)
                        )
                        {
                            faceCandidates.setSize
                            (
                                faceCandidates.size() + 1,
                                eFaces[faceI]
                            );
                        }
                    }
                }
            }

            // Set the mapping for this face
            setFaceMapping(mfIndex, faceCandidates);
        }
    }

    // If modification is coupled, update mapping info.
    if (coupledModification_)
    {
        // Check if the collapse point is present
        // on a processor not involved in the current
        // operation, and update if necessary.
        if (procCouple && !localCouple)
        {
            forAll(procIndices_, pI)
            {
                bool involved = false;

                forAll(slaveMaps, slaveI)
                {
                    // Alias for convenience...
                    const changeMap& slaveMap = slaveMaps[slaveI];

                    if (slaveMap.patchIndex() == pI && slaveMap.index() >= 0)
                    {
                        // Involved in this operation. Break out.
                        involved = true;
                        break;
                    }
                }

                if (involved)
                {
                    continue;
                }

                // Check coupleMaps for point coupling
                const label pointEnum = coupleMap::POINT;

                const coupledPatchInfo& recvMesh = recvPatchMeshes_[pI];
                const coupleMap& cMap = recvMesh.patchMap();

                // Obtain non-const references
                Map<label>& pointMap = cMap.entityMap(pointEnum);
                Map<label>& rPointMap = cMap.reverseEntityMap(pointEnum);

                label sI = -1;

                if (collapsingSlave)
                {
                    if ((sI = cMap.findMaster(pointEnum, collapsePoint)) > -1)
                    {
                        if (rPointMap.found(replacePoint))
                        {
                            rPointMap[replacePoint] = sI;
                        }
                        else
                        {
                            rPointMap.insert(replacePoint, sI);
                        }

                        pointMap[sI] = replacePoint;
                    }
                }
                else
                {
                    if ((sI = cMap.findSlave(pointEnum, collapsePoint)) > -1)
                    {
                        if (pointMap.found(replacePoint))
                        {
                            pointMap[replacePoint] = sI;
                        }
                        else
                        {
                            pointMap.insert(replacePoint, sI);
                        }

                        rPointMap[sI] = replacePoint;
                    }
                }

                if (sI > -1 && debug > 2)
                {
                    Pout<< " Found " << collapsePoint
                        << " on proc: " << procIndices_[pI]
                        << endl;
                }
            }
        }

        forAll(slaveMaps, slaveI)
        {
            // Alias for convenience...
            const changeMap& slaveMap = slaveMaps[slaveI];

            // Skip updates for point-based coupling
            if (slaveMap.index() < 0)
            {
                continue;
            }

            label pI = slaveMap.patchIndex();

            // Fetch the appropriate coupleMap
            const coupleMap* cMapPtr = NULL;

            // Configure the slave edge.
            //  - collapseEdge store this as the first edge
            //    on the list of added edges.
            //  - Since slave edges have already been deleted
            //    at this point, use changeMap information instead.
            edge sEdge
            (
                slaveMap.addedEdgeList()[0].masterObjects()[0],
                slaveMap.addedEdgeList()[0].masterObjects()[1]
            );

            if (localCouple && !procCouple)
            {
                cMapPtr = &(patchCoupling_[pI].patchMap());
            }
            else
            if (procCouple && !localCouple)
            {
                cMapPtr = &(recvPatchMeshes_[pI].patchMap());
            }

            // Configure the slave replacement point
            label scPoint = slaveMap.removedPointList()[0];
            label srPoint = (sEdge[0] == scPoint) ? sEdge[1] : sEdge[0];

            // Alias for convenience...
            const coupleMap& cMap = *cMapPtr;

            label pointEnum = coupleMap::POINT;

            // Obtain references
            Map<label>& pointMap = cMap.entityMap(pointEnum);
            Map<label>& rPointMap = cMap.reverseEntityMap(pointEnum);

            if (collapsingSlave)
            {
                if (rPointMap[replacePoint] == scPoint)
                {
                    pointMap[srPoint] = replacePoint;
                    rPointMap[replacePoint] = srPoint;
                }

                pointMap.erase(rPointMap[collapsePoint]);
                rPointMap.erase(collapsePoint);
            }
            else
            {
                if (pointMap[replacePoint] == scPoint)
                {
                    rPointMap[srPoint] = replacePoint;
                    pointMap[replacePoint] = srPoint;
                }

                rPointMap.erase(pointMap[collapsePoint]);
                pointMap.erase(collapsePoint);
            }

            const label edgeEnum = coupleMap::EDGE;

            // Obtain references
            Map<label>& edgeMap = cMap.entityMap(edgeEnum);
            Map<label>& rEdgeMap = cMap.reverseEntityMap(edgeEnum);

            // Remove for main edge
            if (collapsingSlave)
            {
                edgeMap.erase(rEdgeMap[eIndex]);
                rEdgeMap.erase(eIndex);
            }
            else
            {
                rEdgeMap.erase(edgeMap[eIndex]);
                edgeMap.erase(eIndex);
            }

            // Update replacement coupled edges.
            // This needs to be done before the next step,
            // where edge map entries are removed entirely.
            forAll(brpEdges, edgeI)
            {
                // Skip if we're already done
                if (mappedRpEdge[edgeI])
                {
                    continue;
                }

                label rpIndex = brpEdges[edgeI];
                label rmIndex = brmEdges[edgeI];

                // Check if the coupled edge is in
                // the removed edges list. If so, the
                // edge needs to be updated
                const labelList& srmE = slaveMap.removedEdgeList();

                if (collapsingSlave)
                {
                    if (rEdgeMap.found(rpIndex))
                    {
                        // Map reassigment is necessary only if
                        // the slave edge was removed
                        if (findIndex(srmE, rEdgeMap[rpIndex]) > -1)
                        {
                            // Since rpIndex was found, rmIndex must
                            // also exist, since all faceEdges are mapped
                            edgeMap[rEdgeMap[rmIndex]] = rpIndex;
                            rEdgeMap[rpIndex] = rEdgeMap[rmIndex];
                        }

                        mappedRpEdge[edgeI] = true;
                    }
                }
                else
                {
                    if (edgeMap.found(rpIndex))
                    {
                        // Map reassigment is necessary only if
                        // the slave edge was removed
                        if (findIndex(srmE, edgeMap[rpIndex]) > -1)
                        {
                            // Since rpIndex was found, rmIndex must
                            // also exist, since all faceEdges are mapped
                            rEdgeMap[edgeMap[rmIndex]] = rpIndex;
                            edgeMap[rpIndex] = edgeMap[rmIndex];
                        }

                        mappedRpEdge[edgeI] = true;
                    }
                }
            }

            // Check and remove other edges
            forAll(brmEdges, edgeI)
            {
                // Skip if we're already done
                if (mappedRmEdge[edgeI])
                {
                    continue;
                }

                label rmIndex = brmEdges[edgeI];

                if (collapsingSlave)
                {
                    if (rEdgeMap.found(rmIndex))
                    {
                        edgeMap.erase(rEdgeMap[rmIndex]);
                        rEdgeMap.erase(rmIndex);

                        mappedRmEdge[edgeI] = true;
                    }
                }
                else
                {
                    if (edgeMap.found(rmIndex))
                    {
                        rEdgeMap.erase(edgeMap[rmIndex]);
                        edgeMap.erase(rmIndex);

                        mappedRmEdge[edgeI] = true;
                    }
                }
            }

            const label faceEnum = coupleMap::FACE;

            // Obtain references
            Map<label>& faceMap = cMap.entityMap(faceEnum);
            Map<label>& rFaceMap = cMap.reverseEntityMap(faceEnum);

            // Check and remove other faces
            forAll(bFaces, faceI)
            {
                // Skip if we're already done
                if (mappedFace[faceI])
                {
                    continue;
                }

                label fIndex = bFaces[faceI];

                if (collapsingSlave)
                {
                    if (rFaceMap.found(fIndex))
                    {
                        faceMap.erase(rFaceMap[fIndex]);
                        rFaceMap.erase(fIndex);

                        mappedFace[faceI] = true;
                    }
                }
                else
                {
                    if (faceMap.found(fIndex))
                    {
                        rFaceMap.erase(faceMap[fIndex]);
                        faceMap.erase(fIndex);

                        mappedFace[faceI] = true;
                    }
                }
            }

            // If any interior faces in the master map were
            // converted to boundaries, account for it
            const List<objectMap>& madF = map.addedFaceList();

            forAll(madF, faceI)
            {
                label fIndex = madF[faceI].index();
                const face& mFace = faces_[fIndex];

                triFace cFace
                (
                    pointMap.found(mFace[0]) ? pointMap[mFace[0]] : -1,
                    pointMap.found(mFace[1]) ? pointMap[mFace[1]] : -1,
                    pointMap.found(mFace[2]) ? pointMap[mFace[2]] : -1
                );

                // Skip mapping if all points were not found
                if (cFace[0] == -1 || cFace[1] == -1 || cFace[2] == -1)
                {
                    continue;
                }

                bool matchedFace = false;

                // Select appropriate mesh
                const dynamicTopoFvMesh* meshPtr = NULL;

                if (localCouple)
                {
                    meshPtr = this;
                }
                else
                if (procCouple)
                {
                    meshPtr = &(recvPatchMeshes_[pI].subMesh());
                }

                const dynamicTopoFvMesh& sMesh = *meshPtr;

                // Check all edgeFaces connected to slave points
                forAll(cFace, pointI)
                {
                    label cpIndex = cFace[pointI];

                    const labelList& pEdges = sMesh.pointEdges_[cpIndex];

                    forAll(pEdges, edgeI)
                    {
                        label ceIndex = pEdges[edgeI];

                        const labelList& eFaces = sMesh.edgeFaces_[ceIndex];

                        forAll(eFaces, faceJ)
                        {
                            label checkFace = eFaces[faceJ];

                            if (sMesh.whichPatch(checkFace) == -1)
                            {
                                continue;
                            }

                            const face& sFace = sMesh.faces_[checkFace];

                            if
                            (
                                triFace::compare
                                (
                                    triFace(sFace[0], sFace[1], sFace[2]),
                                    cFace
                                )
                            )
                            {
                                if (debug > 2)
                                {
                                    Pout<< " Found face: " << checkFace
                                        << "::" << cFace
                                        << " with fIndex: " << fIndex
                                        << "::" << mFace
                                        << endl;
                                }

                                // Update faceMaps
                                if (collapsingSlave)
                                {
                                    if (faceMap.found(checkFace))
                                    {
                                        faceMap[checkFace] = fIndex;
                                    }
                                    else
                                    {
                                        faceMap.insert(checkFace, fIndex);
                                    }

                                    if (rFaceMap.found(fIndex))
                                    {
                                        rFaceMap[fIndex] = checkFace;
                                    }
                                    else
                                    {
                                        rFaceMap.insert(fIndex, checkFace);
                                    }

                                }
                                else
                                {
                                    if (rFaceMap.found(checkFace))
                                    {
                                        rFaceMap[checkFace] = fIndex;
                                    }
                                    else
                                    {
                                        rFaceMap.insert(checkFace, fIndex);
                                    }

                                    if (faceMap.found(fIndex))
                                    {
                                        faceMap[fIndex] = checkFace;
                                    }
                                    else
                                    {
                                        faceMap.insert(fIndex, checkFace);
                                    }
                                }

                                matchedFace = true;

                                break;
                            }
                        }

                        if (matchedFace)
                        {
                            break;
                        }
                    }

                    if (matchedFace)
                    {
                        break;
                    }
                }
            }

            // Check and match interior master edges converted to boundaries
            const List<objectMap>& madE = map.addedEdgeList();

            for (label i = 1; i < madE.size(); i++)
            {
                label meIndex = madE[i].index();

                const edge& mE = edges_[meIndex];

                edge cE
                (
                    pointMap.found(mE[0]) ? pointMap[mE[0]] : -1,
                    pointMap.found(mE[1]) ? pointMap[mE[1]] : -1
                );

                // Skip mapping if all points were not found
                if (cE[0] == -1 || cE[1] == -1)
                {
                    continue;
                }

                bool matchedEdge = false;

                // Select appropriate mesh
                const dynamicTopoFvMesh* meshPtr = NULL;

                if (localCouple)
                {
                    meshPtr = this;
                }
                else
                if (procCouple)
                {
                    meshPtr = &(recvPatchMeshes_[pI].subMesh());
                }

                const dynamicTopoFvMesh& sMesh = *meshPtr;

                // Look at edges connected to points
                forAll(cE, pointI)
                {
                    const labelList& pEdges = sMesh.pointEdges_[cE[pointI]];

                    forAll(pEdges, edgeI)
                    {
                        label checkEdge = pEdges[edgeI];

                        const edge& sE = sMesh.edges_[checkEdge];

                        if (sE == cE)
                        {
                            if (debug > 2)
                            {
                                Pout<< " Found edge: " << checkEdge
                                    << "::" << cE
                                    << " with meIndex: " << meIndex
                                    << "::" << mE
                                    << endl;
                            }

                            // Update edgeMaps
                            if (collapsingSlave)
                            {
                                if (edgeMap.found(checkEdge))
                                {
                                    edgeMap[checkEdge] = meIndex;
                                }
                                else
                                {
                                    edgeMap.insert(checkEdge, meIndex);
                                }

                                if (rEdgeMap.found(meIndex))
                                {
                                    rEdgeMap[meIndex] = checkEdge;
                                }
                                else
                                {
                                    rEdgeMap.insert(meIndex, checkEdge);
                                }
                            }
                            else
                            {
                                if (rEdgeMap.found(checkEdge))
                                {
                                    rEdgeMap[checkEdge] = meIndex;
                                }
                                else
                                {
                                    rEdgeMap.insert(checkEdge, meIndex);
                                }

                                if (edgeMap.found(meIndex))
                                {
                                    edgeMap[meIndex] = checkEdge;
                                }
                                else
                                {
                                    edgeMap.insert(meIndex, checkEdge);
                                }
                            }

                            matchedEdge = true;

                            break;
                        }
                    }

                    if (matchedEdge)
                    {
                        break;
                    }
                }
            }

            // If any interior faces in the slave map were
            // converted to boundaries, account for it
            const List<objectMap>& sadF = slaveMap.addedFaceList();

            forAll(sadF, faceI)
            {
                const face* facePtr = NULL;
                label fIndex = sadF[faceI].index();

                if (localCouple && !procCouple)
                {
                    facePtr = &(faces_[fIndex]);
                }
                else
                if (procCouple && !localCouple)
                {
                    facePtr =
                    (
                        &(recvPatchMeshes_[pI].subMesh().faces_[fIndex])
                    );
                }

                const face& sFace = *facePtr;

                triFace cFace
                (
                    rPointMap.found(sFace[0]) ? rPointMap[sFace[0]] : -1,
                    rPointMap.found(sFace[1]) ? rPointMap[sFace[1]] : -1,
                    rPointMap.found(sFace[2]) ? rPointMap[sFace[2]] : -1
                );

                // Skip mapping if all points were not found
                if (cFace[0] == -1 || cFace[1] == -1 || cFace[2] == -1)
                {
                    continue;
                }

                bool matchedFace = false;

                forAll(cFace, pointI)
                {
                    label cpIndex = cFace[pointI];

                    const labelList& pEdges = pointEdges_[cpIndex];

                    forAll(pEdges, edgeI)
                    {
                        label ceIndex = pEdges[edgeI];

                        const labelList& eFaces = edgeFaces_[ceIndex];

                        forAll(eFaces, faceJ)
                        {
                            label checkFace = eFaces[faceJ];

                            if (whichPatch(checkFace) == -1)
                            {
                                continue;
                            }

                            const face& mFace = faces_[checkFace];

                            if
                            (
                                triFace::compare
                                (
                                    triFace(mFace[0], mFace[1], mFace[2]),
                                    cFace
                                )
                            )
                            {
                                if (debug > 2)
                                {
                                    Pout<< " Found face: " << checkFace
                                        << "::" << cFace
                                        << " with fIndex: " << fIndex
                                        << "::" << sFace
                                        << endl;
                                }

                                // Update faceMaps
                                if (collapsingSlave)
                                {
                                    if (rFaceMap.found(checkFace))
                                    {
                                        rFaceMap[checkFace] = fIndex;
                                    }
                                    else
                                    {
                                        rFaceMap.insert(checkFace, fIndex);
                                    }

                                    if (faceMap.found(fIndex))
                                    {
                                        faceMap[fIndex] = checkFace;
                                    }
                                    else
                                    {
                                        faceMap.insert(fIndex, checkFace);
                                    }
                                }
                                else
                                {
                                    if (faceMap.found(checkFace))
                                    {
                                        faceMap[checkFace] = fIndex;
                                    }
                                    else
                                    {
                                        faceMap.insert(checkFace, fIndex);
                                    }

                                    if (rFaceMap.found(fIndex))
                                    {
                                        rFaceMap[fIndex] = checkFace;
                                    }
                                    else
                                    {
                                        rFaceMap.insert(fIndex, checkFace);
                                    }
                                }

                                matchedFace = true;

                                break;
                            }
                        }

                        if (matchedFace)
                        {
                            break;
                        }
                    }

                    if (matchedFace)
                    {
                        break;
                    }
                }
            }

            // Check and match interior slave edges converted to boundaries
            const List<objectMap>& sadE = slaveMap.addedEdgeList();

            for (label i = 1; i < sadE.size(); i++)
            {
                label seIndex = sadE[i].index();

                edge sE(-1, -1);

                if (localCouple && !procCouple)
                {
                    sE = edges_[seIndex];
                }
                else
                if (procCouple && !localCouple)
                {
                    sE = recvPatchMeshes_[pI].subMesh().edges_[seIndex];
                }

                edge cE
                (
                    rPointMap.found(sE[0]) ? rPointMap[sE[0]] : -1,
                    rPointMap.found(sE[1]) ? rPointMap[sE[1]] : -1
                );

                // Skip mapping if all points were not found
                if (cE[0] == -1 || cE[1] == -1)
                {
                    continue;
                }

                bool matchedEdge = false;

                // Look at edges connected to points
                forAll(cE, pointI)
                {
                    const labelList& pEdges = pointEdges_[cE[pointI]];

                    forAll(pEdges, edgeI)
                    {
                        label checkEdge = pEdges[edgeI];

                        const edge& mE = edges_[checkEdge];

                        if (mE == cE)
                        {
                            if (debug > 2)
                            {
                                Pout<< " Found edge: " << checkEdge
                                    << "::" << cE
                                    << " with seIndex: " << seIndex
                                    << "::" << sE
                                    << endl;
                            }

                            // Update edgeMaps
                            if (collapsingSlave)
                            {
                                if (rEdgeMap.found(checkEdge))
                                {
                                    rEdgeMap[checkEdge] = seIndex;
                                }
                                else
                                {
                                    rEdgeMap.insert(checkEdge, seIndex);
                                }

                                if (edgeMap.found(seIndex))
                                {
                                    edgeMap[seIndex] = checkEdge;
                                }
                                else
                                {
                                    edgeMap.insert(seIndex, checkEdge);
                                }
                            }
                            else
                            {
                                if (edgeMap.found(checkEdge))
                                {
                                    edgeMap[checkEdge] = seIndex;
                                }
                                else
                                {
                                    edgeMap.insert(checkEdge, seIndex);
                                }

                                if (rEdgeMap.found(seIndex))
                                {
                                    rEdgeMap[seIndex] = checkEdge;
                                }
                                else
                                {
                                    rEdgeMap.insert(seIndex, checkEdge);
                                }
                            }

                            matchedEdge = true;

                            break;
                        }
                    }

                    if (matchedEdge)
                    {
                        break;
                    }
                }
            }

            // Push operation into coupleMap
            switch (slaveMap.type())
            {
                case 1:
                {
                    cMap.pushOperation
                    (
                        slaveMap.index(),
                        coupleMap::COLLAPSE_FIRST
                    );

                    break;
                }

                case 2:
                {
                    cMap.pushOperation
                    (
                        slaveMap.index(),
                        coupleMap::COLLAPSE_SECOND
                    );

                    break;
                }

                case 3:
                {
                    cMap.pushOperation
                    (
                        slaveMap.index(),
                        coupleMap::COLLAPSE_MIDPOINT
                    );

                    break;
                }
            }
        }

        // Should have found maps for edges by now.
        // (Both local and processor)
        if
        (
            !mappedRmEdge[0] ||
            !mappedRmEdge[1] ||
            !mappedRpEdge[0] ||
            !mappedRpEdge[1]
        )
        {
            FatalErrorIn
            (
                "\n"
                "const changeMap dynamicTopoFvMesh::collapseEdge\n"
                "(\n"
                "    const label eIndex,\n"
                "    label overRideCase,\n"
                "    bool checkOnly,\n"
                "    bool forceOp\n"
                ")\n"
            )
                << "Edge: " << eIndex
                << ": " << edges_[eIndex] << nl
                << " Failed to update edge maps." << nl
                << " Edges : " << nl
                << " brmEdges: " << brmEdges << nl
                << " brpEdges: " << brpEdges << nl
                << " mappedRmEdge: " << mappedRmEdge << nl
                << " mappedRpEdge: " << mappedRpEdge << nl
                << abort(FatalError);
        }

        // Should have found maps for faces by now.
        // (Both local and processor)
        if (!mappedFace[0] || !mappedFace[1])
        {
            FatalErrorIn
            (
                "\n"
                "const changeMap dynamicTopoFvMesh::collapseEdge\n"
                "(\n"
                "    const label eIndex,\n"
                "    label overRideCase,\n"
                "    bool checkOnly,\n"
                "    bool forceOp\n"
                ")\n"
            )
                << "Edge: " << eIndex
                << ": " << edges_[eIndex] << nl
                << " Failed to update face maps." << nl
                << " Faces : " << nl
                << bFaces[0] << " :: " << faces_[bFaces[0]] << nl
                << bFaces[1] << " :: " << faces_[bFaces[1]] << nl
                << " mappedFace: " << mappedFace << nl
                << abort(FatalError);
        }
    }

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    statistics_[4]++;

    // Increment the number of modifications
    statistics_[0]++;

    // Return a succesful collapse
    map.type() = collapseCase;

    return map;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
