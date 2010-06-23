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
    dynamicTopoFvMesh

Description
    Functions specific to coupled connectivity

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*----------------------------------------------------------------------------*/

#include "dynamicTopoFvMesh.H"

#include "triFace.H"
#include "globalMeshData.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Identify coupled patches.
//  - Also builds global shared point information.
//  - Returns true if no coupled patches were found.
bool dynamicTopoFvMesh::identifyCoupledPatches()
{
    bool coupledPatchesAbsent = true;

    // Check if patches are explicitly coupled
    if (patchCoupling_.size())
    {
        coupledPatchesAbsent = false;
    }

    // Maintain a separate list of processor IDs in procIndices.
    // This is done because this sub-domain may talk to processors
    // that share only edges/points.
    if (Pstream::parRun())
    {
        const polyBoundaryMesh& boundary = boundaryMesh();

        forAll(boundary, patchI)
        {
            if (isA<processorPolyPatch>(boundary[patchI]))
            {
                coupledPatchesAbsent = false;

                break;
            }
        }

        // Prepare a list of points for sub-mesh creation.
        //  - Obtain global shared-points information, if necessary.
        if (procIndices_.empty())
        {
            // Fetch the list of global points from polyMesh.
            const globalMeshData& gData = polyMesh::globalData();

            const labelList& spAddr = gData.sharedPointAddr();
            const labelList& spLabels = gData.sharedPointLabels();

            labelList nRecvPoints(Pstream::nProcs(), 0);
            labelListList spBuffer(Pstream::nProcs(), labelList(0));

            if (gData.nGlobalPoints())
            {
                if (debug)
                {
                    Info << " Found "
                         << gData.nGlobalPoints()
                         << " global points."
                         << endl;
                }

                // Send others my addressing.
                for (label proc = 0; proc < Pstream::nProcs(); proc++)
                {
                    if (proc != Pstream::myProcNo())
                    {
                        // Send number of entities first.
                        pWrite(proc, spAddr.size());

                        // Send the buffer.
                        if (spAddr.size())
                        {
                            pWrite(proc, spAddr);
                        }
                    }
                }

                // Receive addressing from others
                for (label proc = 0; proc < Pstream::nProcs(); proc++)
                {
                    if (proc != Pstream::myProcNo())
                    {
                        label procInfoSize = -1;

                        // How many entities am I going to be receiving?
                        pRead(proc, procInfoSize);

                        if (procInfoSize)
                        {
                            // Size the receive buffer.
                            spBuffer[proc].setSize(procInfoSize, -1);

                            // Schedule for receipt.
                            pRead(proc, spBuffer[proc]);
                        }
                    }
                }
            }
            else
            if (debug)
            {
                Info << "Did not find any global points." << endl;
            }

            labelHashSet immNeighbours;
            labelListList procPatchPoints(Pstream::nProcs());

            // Insert my immediate neighbours into the list.
            forAll(boundary, pI)
            {
                if (isA<processorPolyPatch>(boundary[pI]))
                {
                    const processorPolyPatch& pp =
                    (
                        refCast<const processorPolyPatch>(boundary[pI])
                    );

                    label neiProcNo = pp.neighbProcNo();

                    // Insert all boundary points.
                    procPatchPoints[neiProcNo] = pp.meshPoints();

                    // Keep track of immediate neighbours.
                    immNeighbours.insert(neiProcNo);
                }
            }

            if (gData.nGlobalPoints())
            {
                // Wait for all transfers to complete.
                waitForBuffers();

                // Now loop through all processor addressing, and check if
                // any labels coincide with my global shared points.
                // If this is true, we need to be talking to that neighbour
                // as well (if not already).
                for (label proc = 0; proc < Pstream::nProcs(); proc++)
                {
                    if
                    (
                        (proc != Pstream::myProcNo()) &&
                        (!immNeighbours.found(proc))
                    )
                    {
                        bool foundGlobalMatch = false;

                        labelHashSet neiSet;

                        forAll(spBuffer[proc], pointI)
                        {
                            forAll(spAddr, pointJ)
                            {
                                if (spAddr[pointJ] == spBuffer[proc][pointI])
                                {
                                    // Make an entry
                                    neiSet.insert(spLabels[pointJ]);

                                    foundGlobalMatch = true;

                                    break;
                                }
                            }
                        }

                        if (foundGlobalMatch && debug)
                        {
                            Pout << "Additionally talking to processor: "
                                 << proc << endl;
                        }

                        procPatchPoints[proc] = neiSet.toc();
                    }
                }
            }

            // Estimate an initial size
            procIndices_.setSize(Pstream::nProcs());

            // Patch sub meshes need to be prepared in ascending
            // order of neighbouring processors.
            label nTotalProcs = 0;

            forAll(procPatchPoints, procI)
            {
                if (procPatchPoints[procI].size())
                {
                    procIndices_[nTotalProcs++] = procI;
                }
            }

            // Shrink to actual size
            procIndices_.setSize(nTotalProcs);

            // Size the PtrLists.
            sendPatchMeshes_.setSize(nTotalProcs);
            recvPatchMeshes_.setSize(nTotalProcs);

            // Create send/recv patch meshes, and copy
            // the list of points for each processor.
            forAll(procIndices_, pI)
            {
                label proc = procIndices_[pI];

                if (proc < Pstream::myProcNo())
                {
                    sendPatchMeshes_.set
                    (
                        pI,
                        new coupledPatchInfo
                        (
                            *this,               // Reference to this mesh
                            false,               // Not local
                            true,                // Sent to neighbour
                            proc,                // Master index
                            Pstream::myProcNo()  // Slave index
                        )
                    );

                    sendPatchMeshes_[pI].patchMap().subMeshPoints() =
                    (
                        procPatchPoints[procIndices_[pI]]
                    );

                    recvPatchMeshes_.set
                    (
                        pI,
                        new coupledPatchInfo
                        (
                            *this,               // Reference to this mesh
                            false,               // Not local
                            false,               // Not sent to neighbour
                            proc,                // Master index
                            Pstream::myProcNo()  // Slave index
                        )
                    );

                    recvPatchMeshes_[pI].patchMap().subMeshPoints() =
                    (
                        procPatchPoints[procIndices_[pI]]
                    );
                }
                else
                {
                    sendPatchMeshes_.set
                    (
                        pI,
                        new coupledPatchInfo
                        (
                            *this,               // Reference to this mesh
                            false,               // Not local
                            true,                // Sent to neighbour
                            Pstream::myProcNo(), // Master index
                            proc                 // Slave index
                        )
                    );

                    sendPatchMeshes_[pI].patchMap().subMeshPoints() =
                    (
                        procPatchPoints[procIndices_[pI]]
                    );

                    recvPatchMeshes_.set
                    (
                        pI,
                        new coupledPatchInfo
                        (
                            *this,               // Reference to this mesh
                            false,               // Not local
                            false,               // Not sent to neighbour
                            Pstream::myProcNo(), // Master index
                            proc                 // Slave index
                        )
                    );

                    recvPatchMeshes_[pI].patchMap().subMeshPoints() =
                    (
                        procPatchPoints[procIndices_[pI]]
                    );
                }
            }

            if (debug > 3)
            {
                Pout << "Talking to processors: " << procIndices_ << endl;

                forAll(procIndices_, pI)
                {
                    label proc = procIndices_[pI];

                    // Write out points as a VTK
                    if (proc < Pstream::myProcNo())
                    {
                        writeVTK
                        (
                            "subMeshPoints_" +
                            Foam::name(Pstream::myProcNo()) + "to" +
                            Foam::name(proc),
                            sendPatchMeshes_[pI].patchMap().subMeshPoints(),
                            0
                        );
                    }
                    else
                    {
                        writeVTK
                        (
                            "subMeshPoints_" +
                            Foam::name(Pstream::myProcNo()) + "to" +
                            Foam::name(proc),
                            recvPatchMeshes_[pI].patchMap().subMeshPoints(),
                            0
                        );
                    }
                }
            }
        }
    }

    return coupledPatchesAbsent;
}


// Read coupled patch information from dictionary.
void dynamicTopoFvMesh::readCoupledPatches()
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    // Size the initial list
    patchCoupling_.clear();
    patchCoupling_.setSize(boundary.size());

    if (dict_.found("coupledPatches") || mandatory_)
    {
        const dictionary& coupledPatches = dict_.subDict("coupledPatches");

        // Determine master and slave patches
        forAllConstIter(dictionary, coupledPatches, dIter)
        {
            const dictionary& dictI = dIter().dict();

            // Lookup the master / slave patches
            word masterPatch = dictI.lookup("master");
            word slavePatch  = dictI.lookup("slave");

            // Determine patch indices
            label mPatch = boundary.findPatchID(masterPatch);
            label sPatch = boundary.findPatchID(slavePatch);

            if (mPatch == -1 && sPatch == -1)
            {
                // This pair doesn't exist. This might be
                // true for some sub-domains.
                continue;
            }

            if (debug)
            {
                Info << " Adding master: " << masterPatch
                     << " : " << mPatch
                     << " with slave: " << slavePatch
                     << " : " << sPatch << endl;
            }

            // Add to the list if entries are legitimate
            if
            (
                mPatch != sPatch &&
                boundary[mPatch].size() == boundary[sPatch].size()
            )
            {
                // Check whether patches are associated with zones.
                Switch specifyZones
                (
                    dictI.lookup("specifyZones")
                );

                label mZone = -1, sZone = -1;

                if (specifyZones)
                {
                    const faceZoneMesh& faceZones = polyMesh::faceZones();

                    mZone = faceZones.findZoneID
                    (
                        dictI.lookup("masterZone")
                    );

                    sZone = faceZones.findZoneID
                    (
                        dictI.lookup("slaveZone")
                    );
                }

                // Configure a regIOobject for check-in.
                coupleMap cMap
                (
                    IOobject
                    (
                        "coupleMap_"
                      + Foam::name(mPatch)
                      + "_To_"
                      + Foam::name(sPatch)
                      + "_Local",
                        time().timeName(),
                        *this,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE,
                        true
                    ),
                    true,
                    false,
                    mPatch,
                    sPatch
                );

                patchCoupling_.set
                (
                    mPatch,
                    new coupledPatchInfo
                    (
                        *this,
                        cMap,
                        mZone,
                        sZone
                    )
                );
            }
            else
            {
                FatalErrorIn("dynamicTopoFvMesh::readCoupledPatches()")
                    << " Coupled patches are either wrongly specified,"
                    << " or the sizes don't match." << nl
                    << " Master: " << mPatch << ":" << masterPatch
                    << " Size: " << boundary[mPatch].size() << nl
                    << " Slave: " << sPatch << ":" << slavePatch
                    << " Size: " << boundary[sPatch].size() << nl
                    << abort(FatalError);
            }
        }
    }

    // Initialize entitiesToAvoid to some arbitrary size
    entitiesToAvoid_.setSize(50, -1);
}


// Initialize coupled patch connectivity for topology modifications.
//  - Send and receive sub meshes for processor patches.
//  - Made static because this function may be run in a separate thread.
void dynamicTopoFvMesh::initCoupledConnectivity
(
    void *argument
)
{
    // Recast the argument
    dynamicTopoFvMesh *mesh = reinterpret_cast<dynamicTopoFvMesh*>(argument);

    // Identify coupled patches.
    if (mesh->identifyCoupledPatches())
    {
        return;
    }

    // Build and send patch sub-meshes.
    mesh->buildProcessorPatchMeshes();

    // Build inital maps for locally coupled patches.
    mesh->buildLocalCoupledMaps();

    // Build maps for coupled processor patches.
    mesh->buildProcessorCoupledMaps();

    // Synchronize before continuing
    mesh->synchronize();
}


// Handle topology changes for coupled patches
void dynamicTopoFvMesh::handleCoupledPatches()
{
    if (!patchCoupling_.size() && procIndices_.empty())
    {
        return;
    }

    if (debug)
    {
        // Check coupled-patch sizes first.
        forAll(patchCoupling_, patchI)
        {
            if (patchCoupling_(patchI))
            {
                const coupleMap& cMap = patchCoupling_[patchI].patchMap();

                label mSize = patchSizes_[cMap.masterIndex()];
                label sSize = patchSizes_[cMap.slaveIndex()];

                if (mSize != sSize)
                {
                    Pout << "Coupled patch-count is inconsistent." << nl
                         << " Master Patch: " << cMap.masterIndex()
                         << " Count: " << mSize << nl
                         << " Slave Patch: " << cMap.slaveIndex()
                         << " Count: " << sSize
                         << endl;

                    FatalErrorIn("dynamicTopoFvMesh::handleCoupledPatches()")
                        << " Failures were found in connectivity"
                        << " prior to coupled topo-changes."
                        << abort(FatalError);
                }
            }
        }

        Info << "Handling coupled patches...";
    }

    // Set coupled modifications.
    setCoupledModification();

    // Loop through the coupled stack and perform changes.
    if (twoDMesh_)
    {
        if (edgeRefinement_)
        {
            // Initialize the face stack
            initCoupledFaceStack();

            edgeBisectCollapse2D(&(handlerPtr_[0]));
        }

        // Cannot swap on surfaces in 2D, so don't bother.
    }
    else
    {
        if (edgeRefinement_)
        {
            // Initialize the edge stack
            initCoupledEdgeStack();

            edgeBisectCollapse3D(&(handlerPtr_[0]));
        }

        // Re-Initialize the edge stack
        initCoupledEdgeStack();

        swap3DEdges(&(handlerPtr_[0]));
    }

    // Build a list of entities that need to be avoided
    // by regular topo-changes.
    buildEntitiesToAvoid();

    // Reset coupled modifications.
    unsetCoupledModification();

    if (debug)
    {
        Info << "Done." << endl;

        // Check coupled-patch sizes after changes.
        forAll(patchCoupling_, patchI)
        {
            if (patchCoupling_(patchI))
            {
                const coupleMap& cMap = patchCoupling_[patchI].patchMap();

                label mSize = patchSizes_[cMap.masterIndex()];
                label sSize = patchSizes_[cMap.slaveIndex()];

                if (mSize != sSize)
                {
                    Pout << "Coupled patch-count is inconsistent." << nl
                         << " Master Patch: " << cMap.masterIndex()
                         << " Count: " << mSize << nl
                         << " Slave Patch: " << cMap.slaveIndex()
                         << " Count: " << sSize
                         << endl;

                    FatalErrorIn("dynamicTopoFvMesh::handleCoupledPatches()")
                        << " Failures were found in connectivity"
                        << " after coupled topo-changes."
                        << abort(FatalError);
                }
            }
        }
    }
}


// Build patch sub-meshes for processor patches
void dynamicTopoFvMesh::buildProcessorPatchMeshes()
{
    if (procIndices_.empty())
    {
        return;
    }

    // Maintain a list of cells common to multiple processors.
    labelHashSet commonCells;

    forAll(procIndices_, pI)
    {
        label proc = procIndices_[pI];

        coupledPatchInfo& sendMesh = sendPatchMeshes_[pI];

        // Build the subMesh.
        buildProcessorPatchMesh
        (
            sendMesh,
            commonCells
        );

        const coupleMap& scMap = sendMesh.patchMap();

        // Send my sub-mesh to the neighbour.
        pWrite(proc, scMap.nEntities());

        if (debug > 3)
        {
            Pout << "Sending to [" << proc << "]:: nEntities: "
                 << scMap.nEntities()
                 << endl;
        }

        // Send the pointBuffer
        pWrite(proc, scMap.pointBuffer());

        // Send connectivity (points, edges, faces, cells, etc)
        forAll(scMap.entityBuffer(), bufferI)
        {
            if (scMap.entityBuffer(bufferI).size())
            {
                pWrite(proc, scMap.entityBuffer(bufferI));
            }
        }

        // Obtain references
        coupledPatchInfo& recvMesh = recvPatchMeshes_[pI];

        const coupleMap& rcMap = recvMesh.patchMap();

        // First read entity sizes.
        pRead(proc, rcMap.nEntities());

        if (debug > 3)
        {
            Pout << "Receiving from [" << proc << "]:: nEntities: "
                 << rcMap.nEntities()
                 << endl;
        }

        // Size the buffers.
        rcMap.allocateBuffers();

        // Receive the pointBuffer
        pRead(proc, rcMap.pointBuffer());

        // Receive connectivity (points, edges, faces, cells, etc)
        forAll(rcMap.entityBuffer(), bufferI)
        {
            if (rcMap.entityBuffer(bufferI).size())
            {
                pRead(proc, rcMap.entityBuffer(bufferI));
            }
        }
    }

    // We won't wait for all transfers to complete for the moment.
    // Meanwhile, do some other useful work, if possible.
}


// Build patch sub-mesh for a specified processor
// - At this point, procIndices is available as a sorted list
//   of neighbouring processors.
void dynamicTopoFvMesh::buildProcessorPatchMesh
(
    coupledPatchInfo& subMesh,
    labelHashSet& commonCells
)
{
    label nP = 0, nE = 0, nF = 0, nC = 0;

    // Obtain references
    const coupleMap& cMap = subMesh.patchMap();
    const labelList& subMeshPoints = cMap.subMeshPoints();

    Map<label>& rPointMap = cMap.reverseEntityMap(coupleMap::POINT);
    Map<label>& rEdgeMap  = cMap.reverseEntityMap(coupleMap::EDGE);
    Map<label>& rFaceMap  = cMap.reverseEntityMap(coupleMap::FACE);
    Map<label>& rCellMap  = cMap.reverseEntityMap(coupleMap::CELL);

    Map<label>& pointMap = cMap.entityMap(coupleMap::POINT);
    Map<label>& edgeMap  = cMap.entityMap(coupleMap::EDGE);
    Map<label>& faceMap  = cMap.entityMap(coupleMap::FACE);
    Map<label>& cellMap  = cMap.entityMap(coupleMap::CELL);

    // Add all cells connected to points on the subMeshPoints list
    label proc = -1;

    if (cMap.masterIndex() == Pstream::myProcNo())
    {
        proc = cMap.slaveIndex();
    }
    else
    {
        proc = cMap.masterIndex();
    }

    // Check if this is a direct neighbour
    const polyBoundaryMesh& boundary = boundaryMesh();

    label patchIndex = -1;
    bool nonDirectNeighbour = true;

    forAll(boundary, patchI)
    {
        if (isA<processorPolyPatch>(boundary[patchI]))
        {
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[patchI])
            );

            if (pp.neighbProcNo() == proc)
            {
                patchIndex = patchI;
                nonDirectNeighbour = false;

                break;
            }
        }
    }

    // Add sub-mesh points first.
    // Additional halo points will be added later.
    forAll(subMeshPoints, pointI)
    {
        pointMap.insert(nP, subMeshPoints[pointI]);
        rPointMap.insert(subMeshPoints[pointI], nP);
        nP++;
    }

    // Size up the point-buffer with global index information.
    // Direct neighbours do not require any addressing.
    labelList& spBuffer = cMap.entityBuffer(coupleMap::POINT);
    spBuffer = identity(subMeshPoints.size());

    if (nonDirectNeighbour)
    {
        // Fetch the list of global points from polyMesh.
        const globalMeshData& gData = polyMesh::globalData();

        const labelList& spAddr = gData.sharedPointAddr();
        const labelList& spLabels = gData.sharedPointLabels();

        forAll(subMeshPoints, pointI)
        {
            label addrIndex = findIndex(spLabels, subMeshPoints[pointI]);

            // Store global shared-point index
            spBuffer[pointI] = spAddr[addrIndex];
        }
    }

    labelHashSet localCommonCells;

    // Detect all cells surrounding shared points.
    forAll(subMeshPoints, pointI)
    {
        // Loop through pointEdges for this point.
        const labelList& pEdges =
        (
            pointEdges_[subMeshPoints[pointI]]
        );

        forAll(pEdges, edgeI)
        {
            const labelList& eFaces = edgeFaces_[pEdges[edgeI]];

            forAll(eFaces, faceI)
            {
                label own = owner_[eFaces[faceI]];
                label nei = neighbour_[eFaces[faceI]];

                // Check owner cell
                if (!rCellMap.found(own))
                {
                    if (commonCells.found(own))
                    {
                        // Add locally common cells at the end.
                        localCommonCells.set(own, empty());
                    }
                    else
                    {
                        cellMap.insert(nC, own);
                        rCellMap.insert(own, nC);
                        nC++;

                        commonCells.insert(own);
                    }
                }

                // Check neighbour cell
                if (!rCellMap.found(nei) && nei != -1)
                {
                    if (commonCells.found(nei))
                    {
                        // Add locally common cells at the end.
                        localCommonCells.set(nei, empty());
                    }
                    else
                    {
                        cellMap.insert(nC, nei);
                        rCellMap.insert(nei, nC);
                        nC++;

                        commonCells.insert(nei);
                    }
                }
            }
        }
    }

    // Set the number of unique cells at this point.
    cMap.nEntities(coupleMap::UNIQUE_CELL) = nC;

    // Now add locally common cells.
    forAllConstIter(labelHashSet, localCommonCells, cIter)
    {
        cellMap.insert(nC, cIter.key());
        rCellMap.insert(cIter.key(), nC);
        nC++;
    }

    // Allocate the faceMap
    forAllIter(Map<label>, rCellMap, cIter)
    {
        const cell& thisCell = cells_[cIter.key()];

        forAll(thisCell, faceI)
        {
            if (!rFaceMap.found(thisCell[faceI]))
            {
                faceMap.insert(nF, thisCell[faceI]);
                rFaceMap.insert(thisCell[faceI], nF);
                nF++;
            }
        }
    }

    // Allocate the edgeMap. Interior edges need to be detected
    // first and added before boundary ones. Do this in two stages.
    // - The definition of an 'interior' edge here is one which
    //   has all its connected faces in the faceMap. In some sense,
    //   this is a reverse of the traditional definition.
    for (label stage = 0; stage < 2; stage++)
    {
        forAllIter(Map<label>::iterator, rFaceMap, fIter)
        {
            const labelList& fEdges = faceEdges_[fIter.key()];

            forAll(fEdges, edgeI)
            {
                if (!rEdgeMap.found(fEdges[edgeI]))
                {
                    if (stage == 0)
                    {
                        // Check if any of edgeFaces do not belong
                        // to the faceMap.
                        const labelList& eFaces = edgeFaces_[fEdges[edgeI]];

                        bool boundaryEdge = false;

                        forAll(eFaces, faceI)
                        {
                            if (!rFaceMap.found(eFaces[faceI]))
                            {
                                // This face was not found. Designate this
                                // edge as a 'boundary'.
                                boundaryEdge = true;

                                break;
                            }
                        }

                        if (!boundaryEdge)
                        {
                            edgeMap.insert(nE, fEdges[edgeI]);
                            rEdgeMap.insert(fEdges[edgeI], nE);
                            nE++;
                        }
                    }
                    else
                    {
                        // Adding only boundary edges at this stage.
                        edgeMap.insert(nE, fEdges[edgeI]);
                        rEdgeMap.insert(fEdges[edgeI], nE);
                        nE++;
                    }
                }
            }
        }

        // Set the number of internal edges at this point
        if (stage == 0)
        {
            cMap.nEntities(coupleMap::INTERNAL_EDGE) = nE;
        }
    }

    // Set additional points in the pointMap
    forAllIter(Map<label>::iterator, rEdgeMap, eIter)
    {
        const edge& thisEdge = edges_[eIter.key()];

        if (!rPointMap.found(thisEdge[0]))
        {
            pointMap.insert(nP, thisEdge[0]);
            rPointMap.insert(thisEdge[0], nP);
            nP++;
        }

        if (!rPointMap.found(thisEdge[1]))
        {
            pointMap.insert(nP, thisEdge[1]);
            rPointMap.insert(thisEdge[1], nP);
            nP++;
        }
    }

    // Assign sizes to the mesh
    cMap.nEntities(coupleMap::POINT) = nP;
    cMap.nEntities(coupleMap::EDGE)  = nE;
    cMap.nEntities(coupleMap::FACE)  = nF;
    cMap.nEntities(coupleMap::CELL)  = nC;
    cMap.nEntities(coupleMap::SHARED_POINT) = subMeshPoints.size();

    // Size up buffers and fill them
    cMap.allocateBuffers();

    pointField& pBuffer = cMap.pointBuffer();

    forAllIter(Map<label>::iterator, pointMap, pIter)
    {
        pBuffer[pIter.key()] = points_[pIter()];
    }

    // Edge buffer size: 2 points for every edge
    labelList& eBuffer = cMap.entityBuffer(coupleMap::EDGE);

    label index = 0;

    forAllIter(Map<label>::iterator, edgeMap, eIter)
    {
        edge& edgeToCheck = edges_[eIter()];

        eBuffer[index++] = rPointMap[edgeToCheck[0]];
        eBuffer[index++] = rPointMap[edgeToCheck[1]];
    }

    labelList& fpBuffer = cMap.entityBuffer(coupleMap::PATCH_ID);

    label faceIndex = 0;

    forAllIter(Map<label>::iterator, faceMap, fIter)
    {
        // Find the actual patchID for this face.
        label pIndex = whichPatch(fIter());

        // Fill it in, provided it isn't a processor boundary face.
        if (pIndex == -1)
        {
            fpBuffer[faceIndex++] = pIndex;
        }
        else
        if (isA<processorPolyPatch>(boundary[pIndex]))
        {
            fpBuffer[faceIndex++] = -2;
        }
        else
        {
            fpBuffer[faceIndex++] = pIndex;
        }
    }

    index = 0;
    face thisFace(3);
    labelList& fBuffer  = cMap.entityBuffer(coupleMap::FACE);
    labelList& feBuffer = cMap.entityBuffer(coupleMap::FACE_EDGE);

    forAllIter(Map<label>::iterator, faceMap, fIter)
    {
        if (rCellMap.found(owner_[fIter()]))
        {
            thisFace = faces_[fIter()];
        }
        else
        {
            // This face is pointed the wrong way.
            thisFace = faces_[fIter()].reverseFace();
        }

        const labelList& fEdges = faceEdges_[fIter()];

        forAll(fEdges, indexI)
        {
            fBuffer[index] = rPointMap[thisFace[indexI]];
            feBuffer[index++] = rEdgeMap[fEdges[indexI]];
        }
    }

    index = 0;
    labelList& cBuffer  = cMap.entityBuffer(coupleMap::CELL);

    forAllIter(Map<label>::iterator, cellMap, cIter)
    {
        const cell& cellToCheck = cells_[cIter()];

        forAll(cellToCheck, faceI)
        {
            cBuffer[index++] = rFaceMap[cellToCheck[faceI]];
        }
    }

    // Set maps as built.
    subMesh.setBuiltMaps();

    // For debugging purposes...
    if (debug > 3)
    {
        Pout << "Writing out coupledPatchInfo for processor: "
             << proc << endl;

        writeVTK
        (
            "psMesh_" + Foam::name(Pstream::myProcNo())
          + "to" + Foam::name(proc),
            rCellMap.toc()
        );
    }
}


// Build coupled maps for locally coupled patches.
//   - Performs a geometric match initially, since the mesh provides
//     no explicit information for topological coupling.
//   - For each subsequent topology change, maps are updated during
//     the re-ordering stage.
void dynamicTopoFvMesh::buildLocalCoupledMaps()
{
    if (!patchCoupling_.size())
    {
        return;
    }

    if (debug)
    {
        Info << "Building local coupled maps...";
    }

    // Check if a geometric tolerance has been specified.
    scalar gTol = debug::tolerances("processorMatchTol");

    if (dict_.found("gTol") || mandatory_)
    {
        gTol = readScalar(dict_.lookup("gTol"));
    }

    const polyBoundaryMesh& boundary = boundaryMesh();

    forAll(patchCoupling_, patchI)
    {
        if (!patchCoupling_(patchI))
        {
            continue;
        }

        // Build maps only for the first time.
        // coupledPatchInfo is subsequently mapped for each topo-change.
        if (patchCoupling_[patchI].builtMaps())
        {
            continue;
        }

        const coupleMap& cMap = patchCoupling_[patchI].patchMap();

        // Map points on each patch.
        const labelList& mP = boundary[cMap.masterIndex()].meshPoints();
        const labelList& sP = boundary[cMap.slaveIndex()].meshPoints();

        // Sanity check: Do patches have equal number of entities?
        if (mP.size() != sP.size())
        {
            FatalErrorIn("dynamicTopoFvMesh::buildLocalCoupledMaps()")
                << "Patch point sizes are not consistent."
                << abort(FatalError);
        }

        // Check if maps were read from disk.
        // If so, geometric checking is unnecessary.
        if (cMap.entityMap(coupleMap::POINT).size() == 0)
        {
            forAll(mP, indexI)
            {
                bool matched = false;
                scalar minDistance = GREAT;

                forAll(sP, indexJ)
                {
                    scalar distance =
                    (
                        mag(points_[mP[indexI]] - points_[sP[indexJ]])
                    );

                    minDistance = minDistance < distance
                                ? minDistance : distance;

                    if (distance < gTol)
                    {
                        // Add a map entry
                        cMap.mapSlave
                        (
                            coupleMap::POINT,
                            mP[indexI],
                            sP[indexJ]
                        );

                        cMap.mapMaster
                        (
                            coupleMap::POINT,
                            sP[indexJ],
                            mP[indexI]
                        );

                        matched = true;

                        break;
                    }
                }

                if (!matched)
                {
                    FatalErrorIn("dynamicTopoFvMesh::buildLocalCoupledMaps()")
                        << " Failed to match point: " << mP[indexI]
                        << ": " << points_[mP[indexI]]
                        << " within a tolerance of: " << gTol << nl
                        << " Missed by: " << minDistance
                        << abort(FatalError);
                }
            }
        }

        const labelListList& mpF = boundary[cMap.masterIndex()].pointFaces();
        const labelListList& spF = boundary[cMap.slaveIndex()].pointFaces();

        // Abbreviate for convenience
        const Map<label>& pMap = cMap.entityMap(coupleMap::POINT);
        const Map<label>& eMap = cMap.entityMap(coupleMap::EDGE);
        const Map<label>& fMap = cMap.entityMap(coupleMap::FACE);

        // Set up a comparison face.
        face cFace(4);

        // Now that all points are matched,
        // perform topological matching for higher entities.
        forAll(mP, indexI)
        {
            // Fetch global point indices
            label mp = mP[indexI];
            label sp = pMap[mp];

            // Fetch local point indices
            label lmp = boundary[cMap.masterIndex()].whichPoint(mp);
            label lsp = boundary[cMap.slaveIndex()].whichPoint(sp);

            // Fetch patch starts
            label mStart = boundary[cMap.masterIndex()].start();
            label sStart = boundary[cMap.slaveIndex()].start();

            // Match faces for both 2D and 3D.
            const labelList& mpFaces = mpF[lmp];
            const labelList& spFaces = spF[lsp];

            forAll(mpFaces, faceI)
            {
                // Fetch the global face index
                label mfIndex = (mStart + mpFaces[faceI]);

                if (fMap.found(mfIndex))
                {
                    continue;
                }

                const face& mFace = faces_[mfIndex];

                // Configure the face for comparison.
                forAll(mFace, pointI)
                {
                    cFace[pointI] = pMap[mFace[pointI]];
                }

                if (twoDMesh_)
                {
                    bool matched = false;

                    forAll(spFaces, faceJ)
                    {
                        // Fetch the global face index
                        label sfIndex = (sStart + spFaces[faceJ]);

                        const face& sFace = faces_[sfIndex];

                        if (face::compare(cFace, sFace))
                        {
                            // Found the slave. Add a map entry
                            cMap.mapSlave
                            (
                                coupleMap::FACE,
                                mfIndex,
                                sfIndex
                            );

                            cMap.mapMaster
                            (
                                coupleMap::FACE,
                                sfIndex,
                                mfIndex
                            );

                            matched = true;

                            break;
                        }
                    }

                    if (!matched)
                    {
                        FatalErrorIn("buildLocalCoupledMaps()")
                            << " Failed to match face: "
                            << mfIndex << ": " << mFace
                            << abort(FatalError);
                    }
                }
                else
                {
                    label slaveFaceIndex = -1;

                    forAll(spFaces, faceJ)
                    {
                        // Fetch the global face index
                        label sfIndex = (sStart + spFaces[faceJ]);

                        const face& sFace = faces_[sfIndex];

                        if (triFace::compare(triFace(cFace), triFace(sFace)))
                        {
                            // Found the slave. Add a map entry
                            cMap.mapSlave
                            (
                                coupleMap::FACE,
                                mfIndex,
                                sfIndex
                            );

                            cMap.mapMaster
                            (
                                coupleMap::FACE,
                                sfIndex,
                                mfIndex
                            );

                            slaveFaceIndex = sfIndex;

                            break;
                        }
                    }

                    if (slaveFaceIndex == -1)
                    {
                        FatalErrorIn("buildLocalCoupledMaps()")
                            << " Failed to match face: "
                            << mfIndex << ": " << mFace
                            << abort(FatalError);
                    }

                    // Match all edges on this face as well.
                    const labelList& mfEdges = faceEdges_[mfIndex];
                    const labelList& sfEdges = faceEdges_[slaveFaceIndex];

                    forAll(mfEdges, edgeI)
                    {
                        if (eMap.found(mfEdges[edgeI]))
                        {
                            continue;
                        }

                        const edge& mEdge = edges_[mfEdges[edgeI]];

                        // Configure a comparison edge.
                        edge cEdge(pMap[mEdge[0]], pMap[mEdge[1]]);

                        bool matchedEdge = false;

                        forAll(sfEdges, edgeJ)
                        {
                            const edge& sEdge = edges_[sfEdges[edgeJ]];

                            if (cEdge == sEdge)
                            {
                                // Found the slave. Add a map entry
                                cMap.mapSlave
                                (
                                    coupleMap::EDGE,
                                    mfEdges[edgeI],
                                    sfEdges[edgeJ]
                                );

                                cMap.mapMaster
                                (
                                    coupleMap::EDGE,
                                    sfEdges[edgeJ],
                                    mfEdges[edgeI]
                                );

                                matchedEdge = true;

                                break;
                            }
                        }

                        if (!matchedEdge)
                        {
                            FatalErrorIn("buildProcessorCoupledMaps()")
                                << " Failed to match edge: "
                                << mfEdges[edgeI] << ": "
                                << mEdge << nl
                                << " cEdge: " << cEdge
                                << abort(FatalError);
                        }
                    }
                }
            }
        }

        // Set maps as built
        patchCoupling_[patchI].setBuiltMaps();
    }

    if (debug)
    {
        Info << "Done." << endl;
    }
}


// Build coupled maps for coupled processor patches
void dynamicTopoFvMesh::buildProcessorCoupledMaps()
{
    if (procIndices_.empty())
    {
        return;
    }

    Map<label> nPrc;

    const polyBoundaryMesh& boundary = boundaryMesh();

    // Build a list of nearest neighbours.
    forAll(boundary, patchI)
    {
        if (isA<processorPolyPatch>(boundary[patchI]))
        {
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[patchI])
            );

            nPrc.insert(pp.neighbProcNo(), patchI);
        }
    }

    // Wait for all transfers to complete.
    waitForBuffers();

    // Put un-matched faces in a list.
    labelHashSet unMatchedFaces;

    forAll(procIndices_, pI)
    {
        label proc = procIndices_[pI];

        coupledPatchInfo& recvMesh = recvPatchMeshes_[pI];
        const coupleMap& cMap = recvMesh.patchMap();

        // Set the autoPtr.
        recvMesh.setMesh
        (
            proc,
            new dynamicTopoFvMesh
            (
                (*this),
                IOobject
                (
                    "subMesh",
                    time().timeName(),
                    time()
                ),
                cMap.pointBuffer(),
                cMap.nEntities(coupleMap::INTERNAL_EDGE),
                cMap.edges(),
                cMap.faces(),
                cMap.faceEdges(),
                cMap.owner(),
                cMap.neighbour(),
                cMap.cells()
            )
        );

        // Sanity check: Do sub-mesh point sizes match?
        if
        (
            cMap.subMeshPoints().size()
         != cMap.nEntities(coupleMap::SHARED_POINT)
        )
        {
            FatalErrorIn("dynamicTopoFvMesh::buildProcessorCoupledMaps()")
                << " Sub-mesh point sizes don't match." << nl
                << " My procID: " << Pstream::myProcNo() << nl
                << " Slave processor: " << proc << nl
                << abort(FatalError);
        }

        // First, topologically map points based on subMeshPoints.
        const labelList& mP = cMap.subMeshPoints();

        if (nPrc.found(proc))
        {
            // This is an immediate neighbour.
            // All patch points must be matched.

            // Fetch addressing for this patch.
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>
                (
                    boundary[nPrc[proc]]
                )
            );

            const labelList& neiPoints = pp.neighbPoints();

            if (debug)
            {
                if (findIndex(neiPoints, -1) > -1)
                {
                    FatalErrorIn("buildProcessorCoupledMaps()")
                        << " Multiply connected point." << nl
                        << " My procID: " << Pstream::myProcNo() << nl
                        << " Slave processor: " << proc << nl
                        << abort(FatalError);
                }
            }

            forAll(mP, pointI)
            {
                // Add a map entry
                cMap.mapSlave
                (
                    coupleMap::POINT,
                    mP[pointI],
                    neiPoints[pointI]
                );

                cMap.mapMaster
                (
                    coupleMap::POINT,
                    neiPoints[pointI],
                    mP[pointI]
                );
            }
        }
        else
        {
            // Disconnected neighbour.
            // Look at globalPoint addressing for its information.
            const globalMeshData& gD = polyMesh::globalData();
            const labelList& spA = gD.sharedPointAddr();
            const labelList& spL = gD.sharedPointLabels();

            // Fetch global-point addressing from buffer
            const labelList& spB = cMap.entityBuffer(coupleMap::POINT);

            // Match all points first.
            forAll(mP, pointI)
            {
                label maIndex = findIndex(spL, mP[pointI]);
                label saIndex = findIndex(spB, spA[maIndex]);

                // Add a map entry
                cMap.mapSlave
                (
                    coupleMap::POINT,
                    mP[pointI],
                    spB[saIndex]
                );

                cMap.mapMaster
                (
                    coupleMap::POINT,
                    spB[saIndex],
                    mP[pointI]
                );
            }
        }

        if (debug > 1)
        {
            const Map<label>& pMap = cMap.entityMap(coupleMap::POINT);
            const pointField& sPoints = recvMesh.subMesh().points();

            forAllConstIter(Map<label>, pMap, pIter)
            {
                scalar dist = mag(points_[pIter.key()] - sPoints[pIter()]);

                if (dist > debug::tolerances("processorMatchTol"))
                {
                    FatalErrorIn
                    (
                        "dynamicTopoFvMesh::buildProcessorCoupledMaps()"
                    )
                        << " Failed to match point: " << pIter.key()
                        << ": " << points_[pIter.key()]
                        << " with point: " << pIter()
                        << ": " << sPoints[pIter()]
                        << " Missed by: " << dist
                        << abort(FatalError);
                }
            }
        }

        // Set up a comparison face.
        face cFace(4);

        // Now match all faces connected to master points.
        if (nPrc.found(proc))
        {
            // Fetch a global faceList for the slave subMesh
            const faceList& slaveFaces = recvMesh.subMesh().faces();

            // This is an immediate neighbour.
            label mStart = boundary[nPrc[proc]].start();
            label mSize  = boundary[nPrc[proc]].size();

            // Abbreviate for convenience
            const dynamicTopoFvMesh& sMesh = recvMesh.subMesh();
            const Map<label>& pMap = cMap.entityMap(coupleMap::POINT);
            const Map<label>& eMap = cMap.entityMap(coupleMap::EDGE);
            const Map<label>& fMap = cMap.entityMap(coupleMap::FACE);

            // Fetch global pointFaces for the slave.
            const labelListList& spF = recvMesh.subMesh().pointFaces();

            // Match patch faces for both 2D and 3D.
            for(label i = 0; i < mSize; i++)
            {
                // Fetch the global face index
                label mfIndex = (mStart + i);

                if (fMap.found(mfIndex))
                {
                    continue;
                }

                const face& mFace = faces_[mfIndex];

                // Configure the face for comparison.
                forAll(mFace, pointI)
                {
                    cFace[pointI] = pMap[mFace[pointI]];
                }

                // Fetch pointFaces for the zeroth point.
                const labelList& spFaces = spF[cFace[0]];

                if (twoDMesh_)
                {
                    bool matched = false;

                    forAll(spFaces, faceJ)
                    {
                        // Fetch the global face index
                        label sfIndex = spFaces[faceJ];

                        const face& sFace = slaveFaces[sfIndex];

                        if (face::compare(cFace, sFace))
                        {
                            // Found the slave. Add a map entry
                            cMap.mapSlave
                            (
                                coupleMap::FACE,
                                mfIndex,
                                sfIndex
                            );

                            cMap.mapMaster
                            (
                                coupleMap::FACE,
                                sfIndex,
                                mfIndex
                            );

                            matched = true;

                            break;
                        }
                    }

                    if (!matched)
                    {
                        unMatchedFaces.insert(mfIndex);

                        continue;
                    }
                }
                else
                {
                    label slaveFaceIndex = -1;

                    forAll(spFaces, faceJ)
                    {
                        // Fetch the global face index
                        label sfIndex = spFaces[faceJ];

                        const face& sFace = slaveFaces[sfIndex];

                        if (triFace::compare(triFace(cFace), triFace(sFace)))
                        {
                            // Found the slave. Add a map entry
                            cMap.mapSlave
                            (
                                coupleMap::FACE,
                                mfIndex,
                                sfIndex
                            );

                            cMap.mapMaster
                            (
                                coupleMap::FACE,
                                sfIndex,
                                mfIndex
                            );

                            slaveFaceIndex = sfIndex;

                            break;
                        }
                    }

                    if (slaveFaceIndex == -1)
                    {
                        unMatchedFaces.insert(mfIndex);

                        continue;
                    }

                    // Match all edges on this face as well.
                    const labelList& mfEdges = faceEdges_[mfIndex];
                    const labelList& sfEdges =
                    (
                        sMesh.faceEdges_[slaveFaceIndex]
                    );

                    forAll(mfEdges, edgeI)
                    {
                        if (eMap.found(mfEdges[edgeI]))
                        {
                            continue;
                        }

                        const edge& mEdge = edges_[mfEdges[edgeI]];

                        // Configure a comparison edge.
                        edge cEdge(pMap[mEdge[0]], pMap[mEdge[1]]);

                        bool matchedEdge = false;

                        forAll(sfEdges, edgeJ)
                        {
                            const edge& sEdge =
                            (
                                sMesh.edges_[sfEdges[edgeJ]]
                            );

                            if (cEdge == sEdge)
                            {
                                // Found the slave. Add a map entry
                                cMap.mapSlave
                                (
                                    coupleMap::EDGE,
                                    mfEdges[edgeI],
                                    sfEdges[edgeJ]
                                );

                                cMap.mapMaster
                                (
                                    coupleMap::EDGE,
                                    sfEdges[edgeJ],
                                    mfEdges[edgeI]
                                );

                                matchedEdge = true;

                                break;
                            }
                        }

                        if (!matchedEdge)
                        {
                            // Write out the edge
                            writeVTK("mEdge", mfEdges[edgeI], 1);

                            FatalErrorIn("buildProcessorCoupledMaps()")
                                << " Failed to match edge: "
                                << mfEdges[edgeI] << ": "
                                << mEdge << nl
                                << " cEdge: " << cEdge
                                << " for processor: " << proc
                                << abort(FatalError);
                        }
                    }
                }
            }
        }
        else
        {
            // Not a nearest neighbour. Attempt to match
            // edges, provided any common ones exist.
        }

        if (unMatchedFaces.size())
        {
            // Write out the face
            writeVTK("mFaces", unMatchedFaces.toc(), 2);

            FatalErrorIn("dynamicTopoFvMesh::buildProcessorCoupledMaps()")
                << " Unmatched faces were found for processor: " << proc
                << abort(FatalError);
        }
    }
}


// Fill buffers with length-scale info
// and exchange across processors.
void dynamicTopoFvMesh::exchangeLengthBuffers()
{
    forAll(procIndices_, pI)
    {
        coupledPatchInfo& sendMesh = sendPatchMeshes_[pI];
        coupledPatchInfo& recvMesh = recvPatchMeshes_[pI];

        const coupleMap& scMap = sendMesh.patchMap();
        const coupleMap& rcMap = recvMesh.patchMap();

        if (scMap.slaveIndex() == Pstream::myProcNo())
        {
            // Fill in buffers to send.
            scalarList sendLength
            (
                scMap.nEntities(coupleMap::CELL),
                0.0
            );

            Map<label>& cellMap = scMap.entityMap(coupleMap::CELL);

            forAllIter(Map<label>, cellMap, cIter)
            {
                sendLength[cIter.key()] = lengthScale_[cIter()];
            }

            pWrite(scMap.masterIndex(), sendLength);

            if (debug > 4)
            {
                Pout << "Sending to: "
                     << scMap.masterIndex()
                     << " nCells: "
                     << scMap.nEntities(coupleMap::CELL)
                     << endl;
            }
        }

        if (rcMap.masterIndex() == Pstream::myProcNo())
        {
            // Schedule receipt from neighbour
            resizableList<scalar>& recvLength =
            (
                recvMesh.subMesh().lengthScale_
            );

            recvLength.setSize
            (
                rcMap.nEntities(coupleMap::CELL),
                0.0
            );

            // Schedule for receipt
            pRead(rcMap.slaveIndex(), recvLength);

            if (debug > 4)
            {
                Pout << "Receiving from: "
                     << rcMap.slaveIndex()
                     << " nCells: "
                     << rcMap.nEntities(coupleMap::CELL)
                     << endl;
            }
        }
    }
}


// Send length-scale info across processors
void dynamicTopoFvMesh::writeLengthScaleInfo
(
    const labelList& cellLevels,
    const resizableList<scalar>& lengthScale
)
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    // Clear existing buffers
    sendLblBuffer_.clear();
    recvLblBuffer_.clear();
    sendSclBuffer_.clear();
    recvSclBuffer_.clear();

    // Number of sent faces across processors
    labelList nSendFaces(boundary.size(), 0);
    labelList nRecvFaces(boundary.size(), 0);

    // Corresponding face labels
    sendLblBuffer_.setSize(boundary.size());
    recvLblBuffer_.setSize(boundary.size());

    // Length-scales corresponding to face-labels
    sendSclBuffer_.setSize(boundary.size());
    recvSclBuffer_.setSize(boundary.size());

    // Fill send buffers with cell-level and length-scale info.
    forAll(boundary, pI)
    {
        if (isA<processorPolyPatch>(boundary[pI]))
        {
            const labelList& fCells = boundary[pI].faceCells();

            // Set the initial buffer size
            sendLblBuffer_[pI].setSize(boundary[pI].size(), 0);
            sendSclBuffer_[pI].setSize(boundary[pI].size(), 0.0);

            forAll(fCells, faceI)
            {
                label cI = fCells[faceI];

                // Does the adjacent cell have a non-zero level?
                if (cellLevels[cI] > 0)
                {
                    // Fill the send buffer.
                    sendLblBuffer_[pI][nSendFaces[pI]] = faceI;
                    sendSclBuffer_[pI][nSendFaces[pI]] = lengthScale[cI];

                    nSendFaces[pI]++;
                }
            }

            // Resize to actual value
            sendLblBuffer_[pI].setSize(nSendFaces[pI]);
            sendSclBuffer_[pI].setSize(nSendFaces[pI]);

            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[pI])
            );

            label neiProcNo = pp.neighbProcNo();

            // First perform a blocking send/receive of the number of faces.
            pWrite(neiProcNo, nSendFaces[pI]);
            pRead(neiProcNo, nRecvFaces[pI]);
        }
    }

    // Send info to neighbouring processors.
    forAll(boundary, patchI)
    {
        if (isA<processorPolyPatch>(boundary[patchI]))
        {
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[patchI])
            );

            label neiProcNo = pp.neighbProcNo();

            if (debug > 4)
            {
                Pout << " Processor patch " << patchI << ' ' << pp.name()
                     << " communicating with " << neiProcNo
                     << "  Sending: " << nSendFaces[patchI]
                     << endl;
            }

            // Next, perform a non-blocking send/receive of buffers.
            // But only if the buffer size is non-zero.
            if (nSendFaces[patchI] != 0)
            {
                pWrite(neiProcNo, sendLblBuffer_[patchI]);
                pWrite(neiProcNo, sendSclBuffer_[patchI]);
            }

            if (nRecvFaces[patchI] != 0)
            {
                // Size the receive buffers
                recvLblBuffer_[patchI].setSize(nRecvFaces[patchI], 0);
                recvSclBuffer_[patchI].setSize(nRecvFaces[patchI], 0.0);

                pRead(neiProcNo, recvLblBuffer_[patchI]);
                pRead(neiProcNo, recvSclBuffer_[patchI]);
            }
        }
    }
}


// Receive length-scale info across processors
void dynamicTopoFvMesh::readLengthScaleInfo
(
    const label level,
    label& visitedCells,
    labelList& cellLevels,
    resizableList<scalar>& lengthScale,
    labelHashSet& levelCells
)
{
    const polyBoundaryMesh& boundary = boundaryMesh();
    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();

    // Wait for all transfers to complete.
    OPstream::waitRequests();
    IPstream::waitRequests();

    // Now re-visit cells and update length-scales.
    forAll(boundary, patchI)
    {
        if (recvLblBuffer_[patchI].size())
        {
            const labelList& fCells = boundary[patchI].faceCells();

            forAll(recvLblBuffer_[patchI], i)
            {
                label nLabel = recvLblBuffer_[patchI][i];

                label cI = fCells[nLabel];

                label& ngbLevel = cellLevels[cI];

                // For an unvisited cell, update the level
                if (ngbLevel == 0)
                {
                    ngbLevel = level + 1;
                    levelCells.insert(cI);
                    visitedCells++;
                }

                // Visit all neighbours and re-calculate length-scale
                const cell& cellCheck = cells()[cI];
                scalar sumLength = 0.0;
                label nTouchedNgb = 0;

                forAll(cellCheck, faceI)
                {
                    label sLevel = -1, sLCell = -1;
                    label pF = boundary.whichPatch(cellCheck[faceI]);

                    if (pF == -1)
                    {
                        // Internal face. Determine neighbour.
                        if (own[cellCheck[faceI]] == cI)
                        {
                            sLCell = nei[cellCheck[faceI]];
                        }
                        else
                        {
                            sLCell = own[cellCheck[faceI]];
                        }

                        sLevel = cellLevels[sLCell];

                        if ((sLevel < ngbLevel) && (sLevel > 0))
                        {
                            sumLength += lengthScale[sLCell];

                            nTouchedNgb++;
                        }
                    }
                    else
                    if (isA<processorPolyPatch>(boundary[pF]))
                    {
                        // Determine the local index.
                        label local = boundary[pF].whichFace(cellCheck[faceI]);

                        // Is this label present in the list?
                        forAll(recvLblBuffer_[pF], j)
                        {
                            if (recvLblBuffer_[pF][j] == local)
                            {
                                sumLength += recvSclBuffer_[pF][j];

                                nTouchedNgb++;

                                break;
                            }
                        }
                    }
                    else
                    if (fixedPatches_.found(boundary[pF].name()))
                    {
                        sumLength +=
                        (
                            fixedPatches_[boundary[pF].name()][0].scalarToken()
                        );

                        nTouchedNgb++;
                    }
                }

                sumLength /= nTouchedNgb;

                // Scale the length and assign to this cell
                scalar sLength = sumLength*growthFactor_;

                lengthScale[cI] = sLength;
            }
        }
    }
}


// Wait for buffer transfer completion.
void dynamicTopoFvMesh::waitForBuffers() const
{
    if (Pstream::parRun())
    {
        OPstream::waitRequests();
        IPstream::waitRequests();
    }
}


// Parallel blocking send
void dynamicTopoFvMesh::pWrite
(
    const label toID,
    const label& data
) const
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
void dynamicTopoFvMesh::pRead
(
    const label fromID,
    label& data
) const
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
void dynamicTopoFvMesh::pWrite
(
    const label toID,
    const FixedList<Type, Size>& data
) const
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
void dynamicTopoFvMesh::pRead
(
    const label fromID,
    FixedList<Type, Size>& data
) const
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
void dynamicTopoFvMesh::pWrite
(
    const label toID,
    const List<Type>& data
) const
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
void dynamicTopoFvMesh::pRead
(
    const label fromID,
    List<Type>& data
) const
{
    IPstream::read
    (
        Pstream::nonBlocking,
        fromID,
        reinterpret_cast<char*>(&data[0]),
        data.size()*sizeof(Type)
    );
}


// Parallel non-blocking send for resizableLists
template <class Type>
void dynamicTopoFvMesh::pWrite
(
    const label toID,
    const resizableList<Type>& data
) const
{
    OPstream::write
    (
        Pstream::nonBlocking,
        toID,
        reinterpret_cast<const char*>(&data[0]),
        data.size()*sizeof(Type)
    );
}


// Parallel non-blocking receive for resizableLists
template <class Type>
void dynamicTopoFvMesh::pRead
(
    const label fromID,
    resizableList<Type>& data
) const
{
    IPstream::read
    (
        Pstream::nonBlocking,
        fromID,
        reinterpret_cast<char*>(&data[0]),
        data.size()*sizeof(Type)
    );
}


} // End namespace Foam

// ************************************************************************* //
