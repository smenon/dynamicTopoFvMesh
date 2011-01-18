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

#include "Time.H"
#include "triFace.H"
#include "changeMap.H"
#include "matchPoints.H"
#include "globalMeshData.H"
#include "coupledPatchInfo.H"
#include "dynamicTopoFvMesh.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Set coupled modification
void dynamicTopoFvMesh::setCoupledModification() const
{
    coupledModification_ = true;
}


// Unset coupled modification
void dynamicTopoFvMesh::unsetCoupledModification() const
{
    coupledModification_ = false;
}


// Initialize the coupled stack
void dynamicTopoFvMesh::initCoupledStack
(
    const labelHashSet& entities,
    bool useEntities
)
{
    // Clear existing lists/stacks.
    stack(0).clear();

    if (useEntities)
    {
        bool emptyEntity;

        // Initialize the stack with entries
        // in the labelHashSet and return
        forAllConstIter(labelHashSet, entities, eIter)
        {
            // Add only valid entities
            emptyEntity =
            (
                twoDMesh_ ?
                faces_[eIter.key()].empty() :
                edgeFaces_[eIter.key()].empty()
            );

            if (emptyEntity)
            {
                continue;
            }

            stack(0).insert(eIter.key());
        }

        if (debug > 3 && Pstream::parRun())
        {
            Pout<< nl << "Entity stack size: " << stack(0).size() << endl;

            if (debug > 4)
            {
                // Write out stack entities
                labelList stackElements(stack(0).size(), -1);

                forAll(stackElements, elemI)
                {
                    stackElements[elemI] = stack(0)[elemI];
                }

                label elemType = twoDMesh_ ? 2 : 1;

                writeVTK
                (
                    "entityStack_"
                  + Foam::name(Pstream::myProcNo()),
                    stackElements,
                    elemType
                );
            }
        }

        return;
    }

    const polyBoundaryMesh& boundary = boundaryMesh();

    // Loop though boundary faces and check whether
    // they belong to master/slave coupled patches.
    for (label faceI = nOldInternalFaces_; faceI < faces_.size(); faceI++)
    {
        // Add only valid faces
        if (faces_[faceI].empty())
        {
            continue;
        }

        label pIndex = whichPatch(faceI);

        if (pIndex == -1)
        {
            continue;
        }

        // Check if this is a locally coupled master face.
        if (patchCoupling_.size())
        {
            if (patchCoupling_(pIndex))
            {
                // Add this to the coupled modification stack.
                if (twoDMesh_)
                {
                    stack(0).push(faceI);
                }
                else
                {
                    const labelList& mfEdges = faceEdges_[faceI];

                    forAll(mfEdges, edgeI)
                    {
                        // Add this to the coupled modification stack.
                        stack(0).push(mfEdges[edgeI]);
                    }
                }
            }
        }

        // Check if this is a processor patch.
        if (isA<processorPolyPatch>(boundary[pIndex]))
        {
            // Check if this is a master processor patch.
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[pIndex])
            );

            label neiProcID = pp.neighbProcNo();

            if (neiProcID > Pstream::myProcNo())
            {
                // Add this to the coupled modification stack.
                if (twoDMesh_)
                {
                    stack(0).push(faceI);
                }
                else
                {
                    const labelList& mfEdges = faceEdges_[faceI];

                    forAll(mfEdges, edgeI)
                    {
                        stack(0).push(mfEdges[edgeI]);
                    }
                }
            }
        }
    }

    if (debug > 3 && Pstream::parRun())
    {
        Pout<< nl << "Coupled stack size: " << stack(0).size() << endl;

        if (debug > 4)
        {
            // Write out stack entities
            labelList stackElements(stack(0).size(), -1);

            forAll(stackElements, elemI)
            {
                stackElements[elemI] = stack(0)[elemI];
            }

            label elemType = twoDMesh_ ? 2 : 1;

            writeVTK
            (
                "coupledStack_"
              + Foam::name(Pstream::myProcNo()),
                stackElements,
                elemType
            );
        }
    }
}


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
                    Info<< " Found " << gData.nGlobalPoints()
                        << " global points." << endl;
                }

                // Send others my addressing.
                for (label proc = 0; proc < Pstream::nProcs(); proc++)
                {
                    if (proc != Pstream::myProcNo())
                    {
                        // Send number of entities first.
                        meshOps::pWrite(proc, spAddr.size());

                        // Send the buffer.
                        if (spAddr.size())
                        {
                            meshOps::pWrite(proc, spAddr);
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
                        meshOps::pRead(proc, procInfoSize);

                        if (procInfoSize)
                        {
                            // Size the receive buffer.
                            spBuffer[proc].setSize(procInfoSize, -1);

                            // Schedule for receipt.
                            meshOps::pRead(proc, spBuffer[proc]);
                        }
                    }
                }
            }
            else
            if (debug)
            {
                Info<< "Did not find any global points." << endl;
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
                meshOps::waitForBuffers();

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
                            Pout<< "Additionally talking to processor: "
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
                            twoDMesh_,           // 2D or 3D
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
                            twoDMesh_,           // 2D or 3D
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
                            twoDMesh_,           // 2D or 3D
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
                            twoDMesh_,           // 2D or 3D
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
                Pout<< "Talking to processors: " << procIndices_ << endl;

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

    // Clear list
    patchCoupling_.clear();

    if (dict_.found("coupledPatches") || mandatory_)
    {
        // Size the initial list
        patchCoupling_.setSize(boundary.size());

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
                Info<< " Adding master: " << masterPatch << " : " << mPatch
                    << " with slave: " << slavePatch << " : " << sPatch
                    << endl;
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
                    twoDMesh_,
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
                FatalErrorIn("void dynamicTopoFvMesh::readCoupledPatches()")
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
    dynamicTopoFvMesh *mesh = static_cast<dynamicTopoFvMesh*>(argument);

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
}


// Move coupled patch subMesh points
void dynamicTopoFvMesh::moveCoupledSubMeshes()
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (debug)
    {
        Info<< " void dynamicTopoFvMesh::moveCoupledSubMeshes() :"
            << " Moving points for coupled subMeshes."
            << endl;
    }

    forAll(procIndices_, pI)
    {
        label proc = procIndices_[pI];

        const coupledPatchInfo& sPM = sendPatchMeshes_[pI];
        const coupledPatchInfo& rPM = recvPatchMeshes_[pI];

        // Fetch the coupleMap
        const coupleMap& scMap = sPM.patchMap();
        const coupleMap& rcMap = rPM.patchMap();

        Map<label>& pointMap = scMap.entityMap(coupleMap::POINT);

        // Fill point buffers
        pointField& pBuffer = scMap.pointBuffer();
        pointField& opBuffer = scMap.oldPointBuffer();

        forAllConstIter(Map<label>, pointMap, pIter)
        {
            pBuffer[pIter.key()] = points_[pIter()];
            opBuffer[pIter.key()] = oldPoints_[pIter()];
        }

        // Buffers have already been allocated
        // to the right size, so just transfer points

        // Send point buffers to neighbour
        meshOps::pWrite(proc, scMap.pointBuffer());
        meshOps::pWrite(proc, scMap.oldPointBuffer());

        // Receive point buffers from neighbour
        meshOps::pRead(proc, rcMap.pointBuffer());
        meshOps::pRead(proc, rcMap.oldPointBuffer());
    }

    // Wait for transfers to complete before moving on
    meshOps::waitForBuffers();
}


// Insert the cells around the coupled master entity to the mesh
// - Returns a changeMap with a type specifying:
//     1: Insertion was successful
//    -1: Insertion failed
//    -2: Failed because entity was being handled elsewhere
// - The changeMap index specifies the converted mIndex.
const changeMap dynamicTopoFvMesh::insertCells(const label mIndex)
{
    // Prepare the changeMaps
    changeMap map;

    // Maintain face counts for each inserted cell
    Map<label> nCellFaces;

    // Track inserted entities
    List<Map<label> > pointsToInsert(procIndices_.size());
    List<Map<label> > edgesToInsert(procIndices_.size());
    List<Map<label> > facesToInsert(procIndices_.size());
    List<Map<label> > cellsToInsert(procIndices_.size());

    // Track converted entities
    List<Map<label> > edgesToConvert(procIndices_.size());
    List<Map<label> > facesToConvert(procIndices_.size());

    labelList masterConvertPatch(procIndices_.size(), -1);
    labelList slaveConvertPatch(procIndices_.size(), -1);

    // First check to ensure that this case can be handled
    forAll(procIndices_, pI)
    {
        const label cplEnum = twoDMesh_ ? coupleMap::FACE : coupleMap::EDGE;

        // Fetch reference to maps
        const coupledPatchInfo& rPM = recvPatchMeshes_[pI];
        const coupleMap& cMap = rPM.patchMap();

        // Does a coupling exist?
        if (cMap.findSlave(cplEnum, mIndex) == -1)
        {
            continue;
        }

        // If this entity is being handled elsewhere, bail out
        if (procIndices_[pI] < Pstream::myProcNo())
        {
            map.type() = -2;

            return map;
        }
    }

    // Perform a few debug calls before modifications
    if (debug > 1)
    {
        Pout<< nl << nl
            << "Inserting cell(s) around coupled "
            << (twoDMesh_ ? "face: " : "edge: ") << mIndex
            << endl;
    }

    // Agglomerate cells from surrounding subMeshes
    // and add them to this processor.
    forAll(procIndices_, pI)
    {
        const label cplEnum = twoDMesh_ ? coupleMap::FACE : coupleMap::EDGE;

        // Fetch reference to maps
        const coupledPatchInfo& rPM = recvPatchMeshes_[pI];
        const coupleMap& cMap = rPM.patchMap();

        label sIndex = -1;

        if ((sIndex = cMap.findSlave(cplEnum, mIndex)) == -1)
        {
            continue;
        }

        // Fetch references
        const dynamicTopoFvMesh& mesh = rPM.subMesh();
        Map<label>& procCellMap = cellsToInsert[pI];

        if (twoDMesh_)
        {
            // Insert the owner cell
            procCellMap.insert(mesh.owner_[sIndex], -1);
        }
        else
        {
            // Insert all cells connected to this edge
            const labelList& eFaces = mesh.edgeFaces_[sIndex];

            forAll(eFaces, faceI)
            {
                const label own = mesh.owner_[eFaces[faceI]];
                const label nei = mesh.neighbour_[eFaces[faceI]];

                if (!procCellMap.found(own))
                {
                    procCellMap.insert(own, -1);
                }

                if (!procCellMap.found(nei) && (nei != -1))
                {
                    procCellMap.insert(nei, -1);
                }
            }
        }

        // Loop through inserted cells and
        // create an equivalent on this mesh
        forAllIter(Map<label>, procCellMap, cIter)
        {
            const cell& checkCell = mesh.cells_[cIter.key()];

            scalar sLengthScale = -1.0;

            if (edgeRefinement_)
            {
                sLengthScale = mesh.lengthScale_[cIter.key()];
            }

            // Add an empty cell for now, and update
            // with face information at a later stage.
            cIter() = insertCell(cell(checkCell.size()), sLengthScale);

            // Initialize a face counter
            nCellFaces.insert(cIter(), 0);

            if (debug > 3)
            {
                Pout<< " Map cell: " << cIter()
                    << " for cell: " << cIter.key()
                    << " on proc: " << procIndices_[pI]
                    << endl;
            }

            // Add this cell to the map.
            map.addCell(cIter());

            forAll(checkCell, fI)
            {
                label mfIndex = -1, sfIndex = checkCell[fI];

                // Check whether a face mapping exists for this face
                if ((mfIndex = cMap.findMaster(coupleMap::FACE, sfIndex)) > -1)
                {
                    // Mark as a candidate for conversion
                    facesToInsert[pI].insert(sfIndex, mfIndex);
                    facesToConvert[pI].insert(sfIndex, mfIndex);

                    // Add a patch entry as well
                    masterConvertPatch[pI] = whichPatch(mfIndex);
                    slaveConvertPatch[pI] = mesh.whichPatch(sfIndex);
                }
                else
                {
                    // Mark for insertion. Some faces may be duplicated
                    // across processors, but this will be dealt
                    // with at a later stage.
                    facesToInsert[pI].insert(sfIndex, -1);
                }

                // Loop through edges and check whether edge-mapping exists
                const labelList& fEdges = mesh.faceEdges_[sfIndex];

                forAll(fEdges, edgeI)
                {
                    const label eIndex = fEdges[edgeI];
                    const edge& sEdge = mesh.edges_[eIndex];

                    // Meshes in 2D don't have edge-mapping, so check
                    // point maps instead. If either point doesn't exist,
                    // this is an edge that needs to be inserted.
                    label cMs = cMap.findMaster(coupleMap::POINT, sEdge[0]);
                    label cMe = cMap.findMaster(coupleMap::POINT, sEdge[1]);

                    if (!pointsToInsert[pI].found(sEdge[0]))
                    {
                        pointsToInsert[pI].insert(sEdge[0], cMs);
                    }

                    if (!pointsToInsert[pI].found(sEdge[1]))
                    {
                        pointsToInsert[pI].insert(sEdge[1], cMe);
                    }

                    if (!edgesToInsert[pI].found(eIndex))
                    {
                        edgesToInsert[pI].insert(eIndex, -1);
                    }

                    // If both points have maps, this is a conversion edge
                    if (cMs > -1 && cMe > -1)
                    {
                        if (!edgesToConvert[pI].found(eIndex))
                        {
                            edgesToConvert[pI].insert(eIndex, -1);
                        }
                    }
                }
            }
        }

        // Push operation for the slave into coupleMap
        cMap.pushOperation
        (
            sIndex,
            coupleMap::REMOVE_CELL
        );
    }

    // Build a list of edges that need to be converted to interior.
    // - Do this by looking at edges of master face conversion candidates.
    // - Some edges may not need conversion, but deal with this later.
    forAll(facesToConvert, procI)
    {
        // Fetch reference to subMesh
        const coupledPatchInfo& rPM = recvPatchMeshes_[procI];

        const coupleMap& cMap = rPM.patchMap();
        const dynamicTopoFvMesh& mesh = rPM.subMesh();
        const Map<label>& procFaceMap = facesToConvert[procI];

        forAllConstIter(Map<label>, procFaceMap, fIter)
        {
            const labelList& mfEdges = faceEdges_[fIter()];
            const labelList& sfEdges = mesh.faceEdges_[fIter.key()];

            forAll(sfEdges, edgeI)
            {
                if (edgesToInsert[procI][sfEdges[edgeI]] != -1)
                {
                    // Already mapped this edge. Move on.
                    continue;
                }

                // Configure the comparison edge
                const edge& sEdge = mesh.edges_[sfEdges[edgeI]];

                label cMs = cMap.findMaster(coupleMap::POINT, sEdge[0]);
                label cMe = cMap.findMaster(coupleMap::POINT, sEdge[1]);

                edge cEdge(cMs, cMe);

                forAll(mfEdges, edgeJ)
                {
                    const edge& mEdge = edges_[mfEdges[edgeJ]];

                    if (mEdge == cEdge)
                    {
                        edgesToInsert[procI][sfEdges[edgeI]] = mfEdges[edgeJ];
                        edgesToConvert[procI][sfEdges[edgeI]] = mfEdges[edgeJ];
                        break;
                    }
                }
            }
        }
    }

    // Write out prior to modifications
    if (debug > 4)
    {
        forAll(cellsToInsert, procI)
        {
            label procIndex = procIndices_[procI];

            // Fetch reference to subMesh
            const coupledPatchInfo& rPM = recvPatchMeshes_[procI];
            const dynamicTopoFvMesh& mesh = rPM.subMesh();

            if (pointsToInsert[procI].size())
            {
                mesh.writeVTK
                (
                    "pointsToInsert_"
                  + Foam::name(mIndex) + '_'
                  + Foam::name(procIndex),
                    pointsToInsert[procI].toc(),
                    0
                );
            }

            if (edgesToInsert[procI].size())
            {
                mesh.writeVTK
                (
                    "edgesToInsert_"
                  + Foam::name(mIndex) + '_'
                  + Foam::name(procIndex),
                    edgesToInsert[procI].toc(),
                    1
                );
            }

            if (facesToInsert[procI].size())
            {
                mesh.writeVTK
                (
                    "facesToInsert_"
                  + Foam::name(mIndex) + '_'
                  + Foam::name(procIndex),
                    facesToInsert[procI].toc(),
                    2
                );
            }

            if (cellsToInsert[procI].size())
            {
                mesh.writeVTK
                (
                    "cellsToInsert_"
                  + Foam::name(mIndex) + '_'
                  + Foam::name(procIndex),
                    cellsToInsert[procI].toc(),
                    3
                );
            }

            if (edgesToConvert[procI].size())
            {
                mesh.writeVTK
                (
                    "edgesToConvert_"
                  + Foam::name(mIndex) + '_'
                  + Foam::name(procIndex),
                    edgesToConvert[procI].toc(),
                    1
                );
            }

            if (facesToConvert[procI].size())
            {
                mesh.writeVTK
                (
                    "facesToConvert_"
                  + Foam::name(mIndex) + '_'
                  + Foam::name(procIndex),
                    facesToConvert[procI].toc(),
                    2
                );
            }
        }
    }

    // Loop through all insertion points, and merge if necessary
    scalar mergeTol =
    (
        twoDMesh_ ?	0.0 : 1e-4 * mag(edges_[mIndex].vec(points_))
    );

    forAll(pointsToInsert, procI)
    {
        Map<label>& procIPointMap = pointsToInsert[procI];

        // Fetch reference to subMesh
        const coupledPatchInfo& rPMI = recvPatchMeshes_[procI];

        const coupleMap& cMapI = rPMI.patchMap();
        const dynamicTopoFvMesh& meshI = rPMI.subMesh();

        forAllIter(Map<label>, procIPointMap, pItI)
        {
            if (pItI() != -1)
            {
                continue;
            }

            const point& pointI = meshI.points_[pItI.key()];

            bool foundMerge = false;

            forAll(pointsToInsert, procJ)
            {
                if (procJ == procI)
                {
                    continue;
                }

                Map<label>& procJPointMap = pointsToInsert[procJ];

                // Fetch reference to subMesh
                const coupledPatchInfo& rPMJ = recvPatchMeshes_[procJ];

                const coupleMap& cMapJ = rPMJ.patchMap();
                const dynamicTopoFvMesh& meshJ = rPMJ.subMesh();

                // Compare points with this processor
                forAllIter(Map<label>, procJPointMap, pItJ)
                {
                    if (pItJ() != -1)
                    {
                        continue;
                    }

                    const point& pointJ = meshJ.points_[pItJ.key()];

                    if (mag(pointI - pointJ) < mergeTol)
                    {
                        // Make a merge entry
                        label mergePointIndex =
                        (
                            insertPoint
                            (
                                pointJ,
                                meshJ.oldPoints_[pItJ.key()],
                                labelList(1, -1)
                            )
                        );

                        pItI() = mergePointIndex;
                        pItJ() = mergePointIndex;

                        // Update maps for the new point
                        cMapI.mapSlave(coupleMap::POINT, pItI(), pItI.key());
                        cMapI.mapMaster(coupleMap::POINT, pItI.key(), pItI());

                        cMapJ.mapSlave(coupleMap::POINT, pItJ(), pItJ.key());
                        cMapJ.mapMaster(coupleMap::POINT, pItJ.key(), pItJ());

                        if (debug > 3)
                        {
                            Pout<< " Map point: " << mergePointIndex
                                << " pointI: " << pItI.key()
                                << " procI: " << procIndices_[procI]
                                << " pointJ: " << pItJ.key()
                                << " procJ: " << procIndices_[procJ]
                                << endl;
                        }

                        // Add this point to the map.
                        map.addPoint(mergePointIndex);

                        foundMerge = true;
                        break;
                    }
                }

                if (foundMerge)
                {
                    break;
                }
            }

            if (!foundMerge)
            {
                // Found a unique point
                label uniquePointIndex =
                (
                    insertPoint
                    (
                        pointI,
                        meshI.oldPoints_[pItI.key()],
                        labelList(1, -1)
                    )
                );

                pItI() = uniquePointIndex;

                // Update maps for the new point
                cMapI.mapSlave(coupleMap::POINT, pItI(), pItI.key());
                cMapI.mapMaster(coupleMap::POINT, pItI.key(), pItI());

                if (debug > 3)
                {
                    Pout<< " Map point: " << uniquePointIndex
                        << " for point: " << pItI.key()
                        << " on proc: " << procIndices_[procI]
                        << endl;
                }

                // Add this point to the map.
                map.addPoint(uniquePointIndex);
            }
        }
    }

    // Add edges to the mesh, noting that all
    // required points have already been added.
    forAll(edgesToInsert, procI)
    {
        // Fetch reference to subMesh
        const coupledPatchInfo& rPM = recvPatchMeshes_[procI];

        const coupleMap& cMap = rPM.patchMap();
        const dynamicTopoFvMesh& mesh = rPM.subMesh();

        Map<label>& procEdgeMap = edgesToInsert[procI];

        forAllIter(Map<label>, procEdgeMap, eIter)
        {
            // Skip if the edge is a conversion candidate
            if (eIter() != -1)
            {
                continue;
            }

            const edge& sEdge = mesh.edges_[eIter.key()];

            edge newEdge
            (
                cMap.findMaster(coupleMap::POINT, sEdge[0]),
                cMap.findMaster(coupleMap::POINT, sEdge[1])
            );

            // Ensure that this edge hasn't been added before
            label neIndex = -1;
            bool foundDuplicate = false;

            const List<objectMap>& addedEdges = map.addedEdgeList();

            forAll(addedEdges, edgeI)
            {
                neIndex = addedEdges[edgeI].index();

                if (edges_[neIndex] == newEdge)
                {
                    foundDuplicate = true;
                    break;
                }
            }

            if (foundDuplicate)
            {
                // Note the duplicate index for later
                eIter() = neIndex;

                continue;
            }

            // Insert edge with null edgeFaces / edgePoints for now.
            // This can be corrected later.
            eIter() =
            (
                insertEdge
                (
                    masterConvertPatch[procI],
                    newEdge,
                    labelList(0),
                    labelList(0)
                )
            );

            if (debug > 3)
            {
                Pout<< " Map edge: " << eIter() << "::" << newEdge
                    << " for edge: " << eIter.key() << "::" << sEdge
                    << endl;
            }

            // Add this edge to the map.
            map.addEdge(eIter());
        }
    }

    // Add faces to the mesh, noting that all required
    // points and edges have already been added.
    forAll(facesToInsert, procI)
    {
        // Fetch reference to subMesh
        const coupledPatchInfo& rPM = recvPatchMeshes_[procI];

        const coupleMap& cMap = rPM.patchMap();
        const dynamicTopoFvMesh& mesh = rPM.subMesh();

        Map<label>& procFaceMap = facesToInsert[procI];

        forAllIter(Map<label>, procFaceMap, fIter)
        {
            // Skip if the face is a conversion candidate
            if (fIter() != -1)
            {
                continue;
            }

            const face& sFace = mesh.faces_[fIter.key()];
            const labelList& sfEdges = mesh.faceEdges_[fIter.key()];

            face nFace(sFace.size()), cFace(sFace.size());
            labelList nFaceEdges(sfEdges.size());

            // Configure points from map
            forAll(nFace, pI)
            {
                nFace[pI] = cMap.findMaster(coupleMap::POINT, sFace[pI]);
            }

            // This may be a face shared by two other processors.
            // Ensure that duplicates are added only once.
            label dupeOwnCell = -1;
            bool foundDuplicate = false, foundInsertedDuplicate = false;

            forAll(facesToInsert, procJ)
            {
                if (procJ == procI)
                {
                    continue;
                }

                // Fetch reference to subMesh
                const coupledPatchInfo& rmJ = recvPatchMeshes_[procJ];

                const coupleMap& cMJ = rmJ.patchMap();
                const dynamicTopoFvMesh& mJ = rmJ.subMesh();
                const Map<label>& procFaceMapJ = facesToInsert[procJ];

                forAllConstIter(Map<label>, procFaceMapJ, fJ)
                {
                    const face& sFJ = mJ.faces_[fJ.key()];

                    // Discard dissimilar face sizes
                    if (sFJ.size() != cFace.size())
                    {
                        continue;
                    }

                    // Configure points from map
                    forAll(cFace, pI)
                    {
                        cFace[pI] = cMJ.findMaster(coupleMap::POINT, sFJ[pI]);
                    }

                    if (cFace.size() == 3)
                    {
                        // Optimized triangular face comparison
                        if
                        (
                            triFace::compare
                            (
                                triFace(nFace[0], nFace[1], nFace[2]),
                                triFace(cFace[0], cFace[1], cFace[2])
                            )
                        )
                        {
                            foundDuplicate = true;
                        }
                    }
                    else
                    {
                        // Regular face compare
                        if (face::compare(nFace, cFace))
                        {
                            foundDuplicate = true;
                        }
                    }

                    if (foundDuplicate)
                    {
                        // Record the owner for posterity
                        dupeOwnCell =
                        (
                            cellsToInsert[procJ][mJ.owner_[fJ.key()]]
                        );

                        // Was the duplicate face inserted before this?
                        if (fJ() != -1)
                        {
                            // Note the duplicate index for later
                            fIter() = fJ();

                            foundInsertedDuplicate = true;
                        }

                        break;
                    }
                }
            }

            // Face was inserted before, so don't insert again
            if (foundInsertedDuplicate)
            {
                continue;
            }

            // Configure edges from edgesToInsert
            forAll(sfEdges, eI)
            {
                label meIndex = -1;

                // Configure with the appropriate edge
                if (edgesToInsert[procI].found(sfEdges[eI]))
                {
                    meIndex = edgesToInsert[procI][sfEdges[eI]];
                }

                if (meIndex == -1)
                {
                    // Something is wrong here.
                    Pout<< "  Could not find correspondence for slave edge: "
                        << sfEdges[eI] << ":: " << mesh.edges_[sfEdges[eI]]
                        << nl << " mIndex: " << mIndex
                        << abort(FatalError);
                }

                nFaceEdges[eI] = meIndex;
            }

            // Determine patch, owner and neighbour for this face
            label nPatch = -1, nOwner = -1, nNeighbour = -1;

            const polyBoundaryMesh& slaveBoundary = mesh.boundaryMesh();

            label sfPatch = mesh.whichPatch(fIter.key());
            label sFaceOwn = mesh.owner_[fIter.key()];
            label sFaceNei = mesh.neighbour_[fIter.key()];

            label mFaceOwn =
            (
                cellsToInsert[procI].found(sFaceOwn) ?
                cellsToInsert[procI][sFaceOwn] : -1
            );

            label mFaceNei =
            (
                cellsToInsert[procI].found(sFaceNei) ?
                cellsToInsert[procI][sFaceNei] : -1
            );

            // If a duplicate face was found, over-ride neighbour.
            // Face-flipping will be taken care of automatically.
            if (foundDuplicate)
            {
                mFaceNei = dupeOwnCell;
            }

            if (mFaceOwn != -1 && mFaceNei == -1)
            {
                // Boundary face already has correct orientation
                nOwner = mFaceOwn;
                nNeighbour = -1;

                // Determine patch
                if (sfPatch == -1)
                {
                    // Slave face was an interior one
                    nPatch = masterConvertPatch[procI];
                }
                else
                if
                (
                    isA<processorPolyPatch>(slaveBoundary[sfPatch]) ||
                    (sfPatch == (slaveBoundary.size() - 1))
                )
                {
                    // Processor, or 'defaultPatch'
                    nPatch = masterConvertPatch[procI];
                }
                else
                {
                    // Physical type
                    nPatch = sfPatch;
                }
            }
            else
            if (mFaceOwn == -1 && mFaceNei != -1)
            {
                // Boundary face is inverted. Flip it
                nFace = nFace.reverseFace();
                nOwner = mFaceNei;
                nNeighbour = -1;

                // Determine patch
                if (sfPatch == -1)
                {
                    // Slave face was an interior one
                    nPatch = masterConvertPatch[procI];
                }
                else
                if
                (
                    isA<processorPolyPatch>(slaveBoundary[sfPatch]) ||
                    (sfPatch == (slaveBoundary.size() - 1))
                )
                {
                    // Processor, or 'defaultPatch'
                    nPatch = masterConvertPatch[procI];
                }
                else
                {
                    // Physical type
                    nPatch = sfPatch;
                }
            }
            else
            if (mFaceOwn != -1 && mFaceNei != -1)
            {
                // Interior face. Check if a flip is necessary.
                if (mFaceNei < mFaceOwn)
                {
                    nFace = nFace.reverseFace();
                    nOwner = mFaceNei;
                    nNeighbour = mFaceOwn;
                }
                else
                {
                    nOwner = mFaceOwn;
                    nNeighbour = mFaceNei;
                }

                nPatch = -1;
            }
            else
            if (mFaceOwn == -1 && mFaceNei == -1)
            {
                // Something is wrong here.
                Pout<< "Could not find correct owner / neighbour info: " << nl
                    << " Face: " << nFace << nl
                    << " Owner: " << mFaceOwn << nl
                    << " Neighbour: " << mFaceNei << nl
                    << " - Slave Face: " << sFace << nl
                    << " - Slave Patch: " << slaveBoundary[sfPatch].name() << nl
                    << " - Slave Owner: " << sFaceOwn << nl
                    << " - Slave Neighbour: " << sFaceNei << nl
                    << abort(FatalError);
            }

            // Insert the new face
            fIter() =
            (
                insertFace
                (
                    nPatch,
                    nFace,
                    nOwner,
                    nNeighbour
                )
            );

            // Add the faceEdges entry as well
            faceEdges_.append(nFaceEdges);

            // Size up edgeFaces for each edge
            forAll(nFaceEdges, edgeI)
            {
                meshOps::sizeUpList
                (
                    fIter(),
                    edgeFaces_[nFaceEdges[edgeI]]
                );
            }

            if (debug > 3)
            {
                Pout<< " Map face: " << fIter() << "::" << nFace
                    << " Own: " << nOwner << " Nei: " << nNeighbour
                    << " fE: " << nFaceEdges << nl
                    << " for face: " << fIter.key() << "::" << sFace
                    << " Own: " << sFaceOwn << " Nei: " << sFaceNei
                    << " fE: " << sfEdges
                    << endl;
            }

            // Add this face to the map.
            map.addFace(fIter());

            // Update cells
            cells_[nOwner][nCellFaces[nOwner]++] = fIter();

            if (nNeighbour > -1)
            {
                cells_[nNeighbour][nCellFaces[nNeighbour]++] = fIter();
            }
        }
    }

    // Loop through conversion faces
    forAll(facesToConvert, procI)
    {
        // Fetch reference to subMesh
        const coupledPatchInfo& rPM = recvPatchMeshes_[procI];

        const coupleMap& cMap = rPM.patchMap();
        const dynamicTopoFvMesh& mesh = rPM.subMesh();
        const Map<label>& procFaceMap = facesToConvert[procI];

        forAllConstIter(Map<label>, procFaceMap, fIter)
        {
            label fIndex = fIter();

            // Fetch the owner information
            label mFaceOwn = owner_[fIndex];
            label sFaceOwn = mesh.owner_[fIter.key()];
            label newNeighbour = cellsToInsert[procI][sFaceOwn];

            // Insert a new internal face.
            // Orientation should be correct, because this is a boundary
            // face converted to an interior, and adjacent to an added cell.
            label newFaceIndex =
            (
                insertFace
                (
                    -1,
                    faces_[fIndex],
                    mFaceOwn,
                    newNeighbour
                )
            );

            // Add the new faceEdges from the existing face. This may contain
            // edges that need to be converted, but that will be done later.
            faceEdges_.append(faceEdges_[fIndex]);

            // Update map
            map.addFace(newFaceIndex, labelList(1, fIndex));

            // Update the owner cell
            meshOps::replaceLabel
            (
                fIndex,
                newFaceIndex,
                cells_[mFaceOwn]
            );

            // Update the neighbour cell
            cells_[newNeighbour][nCellFaces[newNeighbour]++] = newFaceIndex;

            // Replace edgeFaces with the new face index
            const labelList& fEdges = faceEdges_[newFaceIndex];

            forAll(fEdges, edgeI)
            {
                meshOps::replaceLabel
                (
                    fIndex,
                    newFaceIndex,
                    edgeFaces_[fEdges[edgeI]]
                );
            }

            // Remove the old boundary face
            removeFace(fIndex);

            // Update map
            map.removeFace(fIndex);

            // Obtain references
            Map<label>& faceMap = cMap.entityMap(coupleMap::FACE);
            Map<label>& rFaceMap = cMap.reverseEntityMap(coupleMap::FACE);

            // Erase entries
            rFaceMap.erase(faceMap[fIndex]);
            faceMap.erase(fIndex);

            // For 2D meshes, the boundary face gets converted
            // to an interior one. Note the index for further operations.
            if ((mIndex == fIndex) && twoDMesh_)
            {
                map.index() = newFaceIndex;
            }
        }
    }

    // Loop through conversion edges
    forAll(edgesToConvert, procI)
    {
        // Fetch reference to subMesh
        const coupledPatchInfo& rPM = recvPatchMeshes_[procI];

        const coupleMap& cMap = rPM.patchMap();
        const Map<label>& procEdgeMap = edgesToConvert[procI];

        // Loop through conversion edges
        forAllConstIter(Map<label>, procEdgeMap, eIter)
        {
            label eIndex = eIter();

            bool allInterior = true;

            if (eIndex < 0)
            {
                // Not a conversion edge. Search added edge list for index
                eIndex = edgesToInsert[procI][eIter.key()];

                if (eIndex < 0)
                {
                    // Something's wrong
                    Pout<< " Edge: " << eIter.key()
                        << " :: " << rPM.subMesh().edges_[eIter.key()]
                        << " was never inserted."
                        << abort(FatalError);
                }
            }

            const labelList& eFaces = edgeFaces_[eIndex];

            if (eFaces.size())
            {
                forAll(eFaces, faceI)
                {
                    if (whichPatch(eFaces[faceI]) > -1)
                    {
                        allInterior = false;
                        break;
                    }
                }
            }
            else
            {
                // Edge was deleted before, so skip it
                allInterior = false;
            }

            if (allInterior)
            {
                // This edge needs to be converted to an interior one
                label newEdgeIndex =
                (
                    insertEdge
                    (
                        -1,
                        edges_[eIndex],
                        edgeFaces_[eIndex],
                        labelList(0)
                    )
                );

                // Update map
                map.addEdge(newEdgeIndex, labelList(1, eIndex));

                // Update faceEdges information for all connected faces
                forAll(eFaces, faceI)
                {
                    meshOps::replaceLabel
                    (
                        eIndex,
                        newEdgeIndex,
                        faceEdges_[eFaces[faceI]]
                    );
                }

                // Remove the old boundary edge
                removeEdge(eIndex);

                // Update map
                map.removeEdge(eIndex);

                // For 3D meshes, the boundary edge gets converted
                // to an interior one. Note the index for further operations.
                if ((mIndex == eIndex) && !twoDMesh_)
                {
                    map.index() = newEdgeIndex;
                }

                // Obtain references. Updates are actually only in 3D.
                Map<label>& edgeMap = cMap.entityMap(coupleMap::EDGE);
                Map<label>& rEdgeMap = cMap.reverseEntityMap(coupleMap::EDGE);

                // Erase entries
                rEdgeMap.erase(edgeMap[eIndex]);
                edgeMap.erase(eIndex);
            }
        }
    }

    if (!twoDMesh_)
    {
        // Fix edgePoints for all new / converted edges
        DynamicList<label> fixEdges(10);

        forAll(edgesToInsert, procI)
        {
            const Map<label>& procEdgeMap = edgesToInsert[procI];

            forAllConstIter(Map<label>, procEdgeMap, eIter)
            {
                if ((findIndex(fixEdges, eIter()) < 0) && (eIter() != mIndex))
                {
                    fixEdges.append(eIter());
                }
            }
        }

        // Fix the conversion edge as well
        fixEdges.append(map.index());

        // Sequentially fix all edges
        forAll(fixEdges, indexI)
        {
            buildEdgePoints(fixEdges[indexI]);
        }
    }

    List<changeMap> slaveMaps(procIndices_.size());

    // Loop through all processors, and remove cells
    forAll(cellsToInsert, procI)
    {
        const Map<label>& procCellMap = cellsToInsert[procI];

        if (procCellMap.empty())
        {
            // Set type to something recognizable
            slaveMaps[procI].type() = -7;

            continue;
        }

        // Prepare a list of cells
        const labelList cellsToRemove = procCellMap.toc();

        // Fetch non-const reference to subMesh
        coupledPatchInfo& rPM = recvPatchMeshes_[procI];
        dynamicTopoFvMesh& mesh = rPM.subMesh();

        // Now remove cells for this processor
        slaveMaps[procI] =
        (
            mesh.removeCells
            (
                cellsToRemove,
                slaveConvertPatch[procI]
            )
        );
    }

    // Now map entities from the removeCells operation
    forAll(slaveMaps, procI)
    {
        const changeMap& slaveMap = slaveMaps[procI];

        // Skip empty entities
        if (slaveMap.type() == -7)
        {
            continue;
        }

        const coupledPatchInfo& rPM = recvPatchMeshes_[procI];

        const coupleMap& cMap = rPM.patchMap();
        const dynamicTopoFvMesh& mesh = rPM.subMesh();

        // Map faces from patches.
        const Map<label>& procFaceMap = facesToInsert[procI];
        const List<objectMap>& asfList = slaveMap.addedFaceList();

        forAll(asfList, indexI)
        {
            label sfIndex = asfList[indexI].index(), mfIndex = -1;

            if (mesh.whichPatch(sfIndex) == slaveConvertPatch[procI])
            {
                // Configure a comparison face
                face cFace(mesh.faces_[sfIndex]);

                forAll(cFace, pI)
                {
                    cFace[pI] = cMap.findMaster(coupleMap::POINT, cFace[pI]);
                }

                bool foundMatch = false;

                forAllConstIter(Map<label>, procFaceMap, fIter)
                {
                    const face& mFace = faces_[fIter()];

                    // Discard dissimilar face sizes
                    if (mFace.size() != cFace.size())
                    {
                        continue;
                    }

                    if (cFace.size() == 3)
                    {
                        // Optimized triangular face comparison
                        if
                        (
                            triFace::compare
                            (
                                triFace(mFace[0], mFace[1], mFace[2]),
                                triFace(cFace[0], cFace[1], cFace[2])
                            )
                        )
                        {
                            mfIndex = fIter();
                            foundMatch = true;
                            break;
                        }
                    }
                    else
                    {
                        // Regular face compare
                        if (face::compare(mFace, cFace))
                        {
                            mfIndex = fIter();
                            foundMatch = true;
                            break;
                        }
                    }
                }

                if (foundMatch)
                {
                    // Update maps for the new face
                    cMap.mapSlave(coupleMap::FACE, mfIndex, sfIndex);
                    cMap.mapMaster(coupleMap::FACE, sfIndex, mfIndex);
                }
                else
                {
                    // Something is wrong here.
                    Pout<< " Could not find master face for: " << nl
                        << "  Slave face: " << sfIndex
                        << "  :: " << mesh.faces_[sfIndex] << nl
                        << "  cFace: " << cFace << nl
                        << abort(FatalError);
                }
            }
        }

        // Map edges for 3D meshes
        if (!twoDMesh_)
        {
            // Map edges from patches.
            const Map<label>& procEdgeMap = edgesToInsert[procI];
            const List<objectMap>& aseList = slaveMap.addedEdgeList();

            forAll(aseList, indexI)
            {
                label seIndex = aseList[indexI].index(), meIndex = -1;

                if (mesh.whichEdgePatch(seIndex) == slaveConvertPatch[procI])
                {
                    // Configure a comparison edge
                    edge cEdge(mesh.edges_[seIndex]);

                    cEdge[0] = cMap.findMaster(coupleMap::POINT, cEdge[0]);
                    cEdge[1] = cMap.findMaster(coupleMap::POINT, cEdge[1]);

                    bool foundMatch = false;

                    forAllConstIter(Map<label>, procEdgeMap, eIter)
                    {
                        const edge& mEdge = edges_[eIter()];

                        if (mEdge == cEdge)
                        {
                            meIndex = eIter();
                            foundMatch = true;
                            break;
                        }
                    }

                    if (foundMatch)
                    {
                        // Update maps for the new edge
                        cMap.mapSlave(coupleMap::EDGE, meIndex, seIndex);
                        cMap.mapMaster(coupleMap::EDGE, seIndex, meIndex);
                    }
                    else
                    {
                        // Something is wrong here.
                        Pout<< " Could not find master edge for: " << nl
                            << "  Slave edge: " << seIndex
                            << "  :: " << mesh.edges_[seIndex] << nl
                            << "  cEdge: " << cEdge << nl
                            << abort(FatalError);
                    }
                }
            }
        }
    }

    // Write out cells for debug, if necessary
    if (debug > 4)
    {
        const List<objectMap>& acList = map.addedCellList();

        labelList addedCells(acList.size(), -1);

        forAll(acList, indexI)
        {
            addedCells[indexI] = acList[indexI].index();
        }

        // Write it out
        writeVTK("insertCells(" + Foam::name(mIndex) + ')',	addedCells);
    }

    // Specify that the operation was successful
    map.type() = 1;

    // Return the changeMap
    return map;
}


// Remove the specified cells from the mesh,
// and add internal faces/edges to the specified patch
// - Returns a changeMap with a type specifying:
//     1: Operation was successful
//    -1: Operation failed
const changeMap dynamicTopoFvMesh::removeCells
(
    const labelList& cList,
    const label patch
)
{
    if (cList.empty() || patch < 0)
    {
        FatalErrorIn
        (
            "const changeMap dynamicTopoFvMesh::removeCells"
            "(const labelList& cList, const label patch)"
        )
            << " Wrong arguments. " << nl
            << " cList: " << cList << nl
            << " patch: " << patch << nl
            << abort(FatalError);
    }

    // Perform a few debug calls before modifications
    if (debug > 1)
    {
        Pout<< nl << nl
            << "Removing cell(s): " << cList
            << " and adding internal faces / edges to patch: "
            << boundaryMesh()[patch].name()
            << endl;
    }

    changeMap map;

    labelHashSet pointsToRemove, edgesToRemove, facesToRemove;
    Map<label> facesToConvert, edgesToConvert;

    // First loop through all cells and accumulate
    // a set of faces to be removed/converted.
    forAll(cList, cellI)
    {
        const cell& cellToCheck = cells_[cList[cellI]];

        forAll(cellToCheck, faceI)
        {
            label own = owner_[cellToCheck[faceI]];
            label nei = neighbour_[cellToCheck[faceI]];

            if (nei == -1)
            {
                if (!facesToRemove.found(cellToCheck[faceI]))
                {
                    facesToRemove.insert(cellToCheck[faceI]);
                }
            }
            else
            if
            (
                (findIndex(cList, own) != -1) &&
                (findIndex(cList, nei) != -1)
            )
            {
                if (!facesToRemove.found(cellToCheck[faceI]))
                {
                    facesToRemove.insert(cellToCheck[faceI]);
                }
            }
            else
            {
                facesToConvert.set(cellToCheck[faceI], -1);
            }
        }
    }

    // Add all edges as candidates for conversion.
    // Some of these will be removed altogether.
    forAllConstIter(labelHashSet, facesToRemove, fIter)
    {
        const labelList& fEdges = faceEdges_[fIter.key()];

        forAll(fEdges, edgeI)
        {
            if (whichEdgePatch(fEdges[edgeI]) == patch)
            {
                // Make an identical map
                edgesToConvert.set(fEdges[edgeI], fEdges[edgeI]);
            }
            else
            {
                edgesToConvert.set(fEdges[edgeI], -1);
            }
        }
    }

    forAllConstIter(Map<label>, facesToConvert, fIter)
    {
        const labelList& fEdges = faceEdges_[fIter.key()];

        forAll(fEdges, edgeI)
        {
            if (whichEdgePatch(fEdges[edgeI]) == patch)
            {
                // Make an identical map
                edgesToConvert.set(fEdges[edgeI], fEdges[edgeI]);
            }
            else
            {
                edgesToConvert.set(fEdges[edgeI], -1);
            }
        }
    }

    // Build a list of edges to be removed.
    forAllConstIter(Map<label>, edgesToConvert, eIter)
    {
        const labelList& eFaces = edgeFaces_[eIter.key()];

        bool allRemove = true;

        forAll(eFaces, faceI)
        {
            if (facesToConvert.found(eFaces[faceI]))
            {
                allRemove = false;
                break;
            }
        }

        if (allRemove)
        {
            if (!edgesToRemove.found(eIter.key()))
            {
                edgesToRemove.insert(eIter.key());
            }
        }
    }

    // Weed-out the conversion list.
    forAllConstIter(labelHashSet, edgesToRemove, eIter)
    {
        edgesToConvert.erase(eIter.key());
    }

    // Build a set of points to be removed.
    if (!twoDMesh_)
    {
        forAllConstIter(labelHashSet, edgesToRemove, eIter)
        {
            const edge& edgeToCheck = edges_[eIter.key()];

            forAll(edgeToCheck, pointI)
            {
                const labelList& pEdges = pointEdges_[edgeToCheck[pointI]];

                bool allRemove = true;

                forAll(pEdges, edgeI)
                {
                    if (!edgesToRemove.found(pEdges[edgeI]))
                    {
                        allRemove = false;
                        break;
                    }
                }

                if (allRemove)
                {
                    if (!pointsToRemove.found(edgeToCheck[pointI]))
                    {
                        pointsToRemove.insert(edgeToCheck[pointI]);
                    }
                }
            }
        }
    }

    forAllIter(Map<label>, edgesToConvert, eIter)
    {
        const labelList& eFaces = edgeFaces_[eIter.key()];

        label nConvFaces = 0;

        forAll(eFaces, faceI)
        {
            if (facesToConvert.found(eFaces[faceI]))
            {
                nConvFaces++;
            }
        }

        if (nConvFaces > 2)
        {
            Pout<< "Invalid conversion. Bailing out." << endl;
            return map;
        }
    }

    // Write out candidates for post-processing
    if (debug > 2)
    {
        writeVTK("pointsToRemove", pointsToRemove.toc(), 0);
        writeVTK("edgesToRemove", edgesToRemove.toc(), 1);
        writeVTK("facesToRemove", facesToRemove.toc(), 2);
        writeVTK("cellsToRemove", cList, 3);
        writeVTK("edgesToConvert", edgesToConvert.toc(), 1);
        writeVTK("facesToConvert", facesToConvert.toc(), 2);
    }

    // Loop through all faces for conversion, check orientation
    // and create new faces in their place.
    forAllIter(Map<label>, facesToConvert, fIter)
    {
        // Check if this internal face is oriented properly.
        face newFace;
        label newOwner = -1;
        labelList fEdges = faceEdges_[fIter.key()];

        if (findIndex(cList, neighbour_[fIter.key()]) != -1)
        {
            // Orientation is correct
            newFace = faces_[fIter.key()];
            newOwner = owner_[fIter.key()];
        }
        else
        if (findIndex(cList, owner_[fIter.key()]) != -1)
        {
            // Face is to be reversed.
            newFace = faces_[fIter.key()].reverseFace();
            newOwner = neighbour_[fIter.key()];

            setFlip(fIter.key());
        }
        else
        {
            // Something's terribly wrong
            FatalErrorIn
            (
                "\n"
                "const changeMap dynamicTopoFvMesh::removeCells\n"
                "(\n"
                "    const labelList& cList,\n"
                "    const label patch\n"
                ")\n"
            )
                << nl << " Invalid mesh. "
                << abort(FatalError);
        }

        // Insert the reconfigured face at the boundary.
        fIter() =
        (
            insertFace
            (
                patch,
                newFace,
                newOwner,
                -1
            )
        );

        // Add the faceEdges entry.
        // Edges will be corrected later.
        faceEdges_.append(fEdges);

        // Add this face to the map.
        map.addFace(fIter());

        // Replace cell with the new face label
        meshOps::replaceLabel
        (
            fIter.key(),
            fIter(),
            cells_[newOwner]
        );

        // Remove the internal face.
        removeFace(fIter.key());
    }

    // Create a new edge for each converted edge
    forAllIter(Map<label>, edgesToConvert, eIter)
    {
        if (eIter() == -1)
        {
            // Create copies before appending.
            edge newEdge = edges_[eIter.key()];
            labelList eFaces = edgeFaces_[eIter.key()];
            labelList ePoints;

            if (!twoDMesh_)
            {
                ePoints = edgePoints_[eIter.key()];
            }

            eIter() =
            (
                insertEdge
                (
                    patch,
                    newEdge,
                    eFaces,
                    ePoints
                )
            );

            // Add this edge to the map.
            map.addEdge(eIter());

            // Remove the edge
            removeEdge(eIter.key());
        }
    }

    // Loop through all faces for conversion, and replace edgeFaces.
    forAllConstIter(Map<label>, facesToConvert, fIter)
    {
        // Make a copy, because this list is going to
        // be modified within this loop.
        labelList fEdges = faceEdges_[fIter()];

        forAll(fEdges, edgeI)
        {
            if (edgesToConvert.found(fEdges[edgeI]))
            {
                meshOps::replaceLabel
                (
                    fIter.key(),
                    fIter(),
                    edgeFaces_[edgesToConvert[fEdges[edgeI]]]
                );

                meshOps::replaceLabel
                (
                    fEdges[edgeI],
                    edgesToConvert[fEdges[edgeI]],
                    faceEdges_[fIter()]
                );
            }
        }
    }

    // Loop through all edges for conversion, and size-down edgeFaces.
    forAllConstIter(Map<label>, edgesToConvert, eIter)
    {
        // Make a copy, because this list is going to
        // be modified within this loop.
        labelList eFaces = edgeFaces_[eIter()];

        forAll(eFaces, faceI)
        {
            if (facesToRemove.found(eFaces[faceI]))
            {
                meshOps::sizeDownList
                (
                    eFaces[faceI],
                    edgeFaces_[eIter()]
                );
            }

            // Replace old edges with new ones.
            labelList& fEdges = faceEdges_[eFaces[faceI]];

            forAll(fEdges, edgeI)
            {
                if (edgesToConvert.found(fEdges[edgeI]))
                {
                    fEdges[edgeI] = edgesToConvert[fEdges[edgeI]];
                }
            }
        }
    }

    // At this point, edgeFaces is consistent.
    // Correct edge-points for all converted edges
    if (!twoDMesh_)
    {
        forAllConstIter(Map<label>, edgesToConvert, eIter)
        {
            buildEdgePoints(eIter());
        }
    }

    // Remove unwanted faces
    forAllConstIter(labelHashSet, facesToRemove, fIter)
    {
        removeFace(fIter.key());
    }

    // Remove unwanted edges
    forAllConstIter(labelHashSet, edgesToRemove, eIter)
    {
        removeEdge(eIter.key());
    }

    // Remove unwanted points
    forAllConstIter(labelHashSet, pointsToRemove, pIter)
    {
        removePoint(pIter.key());
    }

    // Remove all cells
    forAll(cList, cellI)
    {
        removeCell(cList[cellI]);
    }

    // Set the flag
    topoChangeFlag_ = true;

    // Specify that the operation was successful
    map.type() = 1;

    return map;
}


// Handle topology changes for coupled patches
void dynamicTopoFvMesh::handleCoupledPatches
(
    labelHashSet& entities
)
{
    if (patchCoupling_.empty() && procIndices_.empty())
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
                    Pout<< "Coupled patch-count is inconsistent." << nl
                        << " Master Patch: " << cMap.masterIndex()
                        << " Count: " << mSize << nl
                        << " Slave Patch: " << cMap.slaveIndex()
                        << " Count: " << sSize
                        << endl;

                    FatalErrorIn
                    (
                        "void dynamicTopoFvMesh::handleCoupledPatches()"
                    )
                        << " Failures were found in connectivity"
                        << " prior to coupled topo-changes."
                        << abort(FatalError);
                }
            }
        }

        Info<< "Handling coupled patches...";
    }

    // Set coupled modifications.
    setCoupledModification();

    // Loop through the coupled stack and perform changes.
    if (edgeRefinement_)
    {
        // Initialize the coupled stack
        initCoupledStack(entities, false);

        edgeRefinementEngine(&(handlerPtr_[0]));
    }

    // Re-Initialize the stack
    initCoupledStack(entities, false);

    if (twoDMesh_)
    {
        swap2DEdges(&(handlerPtr_[0]));
    }
    else
    {
        swap3DEdges(&(handlerPtr_[0]));
    }

    // Build a list of entities that need to be avoided
    // by regular topo-changes.
    buildEntitiesToAvoid(entities, true);

    // Reset coupled modifications.
    unsetCoupledModification();

    if (debug)
    {
        Info<< "Done." << endl;

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
                    Pout<< "Coupled patch-count is inconsistent." << nl
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

    // Schedule transfer of topology operations across processors
    forAll(procIndices_, pI)
    {
        label proc = procIndices_[pI];

        coupledPatchInfo& sPM = sendPatchMeshes_[pI];
        coupledPatchInfo& rPM = recvPatchMeshes_[pI];

        if (proc < Pstream::myProcNo())
        {
            const coupleMap& cMap = sPM.patchMap();

            // How many entities am I receiving..
            label nEntities = -1;

            meshOps::pRead(proc, nEntities);

            if (debug > 3)
            {
                Pout<< " Op tranfer:"
                    << " Receiving from [" << proc << "]:: nEntities: "
                    << nEntities << endl;
            }

            if (nEntities)
            {
                // Size up the receipt buffers
                labelList& indices = cMap.entityIndices();
                labelList& operations = cMap.entityOperations();

                indices.setSize(nEntities);
                operations.setSize(nEntities);

                // Schedule indices and operations for receipt
                meshOps::pRead(proc, indices);
                meshOps::pRead(proc, operations);
            }
        }
        else
        {
            const coupleMap& cMap = rPM.patchMap();

            label nEntities = cMap.entityIndices().size();

            if (debug > 3)
            {
                Pout<< " Op tranfer:"
                    << " Sending to [" << proc << "]:: "
                    << " nEntities: " << nEntities << nl
                    << "  entityIndices: " << cMap.entityIndices() << nl
                    << "  entityOperations: " << cMap.entityOperations() << nl
                    << endl;
            }

            meshOps::pWrite(proc, nEntities);

            if (nEntities)
            {
                // Schedule transfer to processor
                meshOps::pWrite(proc, cMap.entityIndices());
                meshOps::pWrite(proc, cMap.entityOperations());
            }
        }
    }

    // We won't wait for transfers to complete for the moment,
    // and will deal with operations once the internal mesh
    // has been dealt with.
}


// Synchronize topology operations across processors
void dynamicTopoFvMesh::syncCoupledPatches(labelHashSet& entities)
{
    if (!Pstream::parRun())
    {
        return;
    }

    // Wait for all transfers to complete.
    meshOps::waitForBuffers();

    // Buffer for cell-removal
    DynamicList<label> rCellList(10);

    forAll(procIndices_, pI)
    {
        label proc = procIndices_[pI];

        const coupledPatchInfo& sPM = sendPatchMeshes_[pI];

        if (proc < Pstream::myProcNo())
        {
            const coupleMap& cMap = sPM.patchMap();
            const labelList& indices = cMap.entityIndices();
            const labelList& operations = cMap.entityOperations();

            // Fetch the appropriate map
            const Map<label>* entityMapPtr =
            (
                twoDMesh_ ?
                &(cMap.entityMap(coupleMap::FACE)) :
                &(cMap.entityMap(coupleMap::EDGE))
            );

            const Map<label>& entityMap = *entityMapPtr;

            // Find the appropriate processor patch
            // for cell-removal operations
            label procPatch = -1;

            const polyBoundaryMesh& boundary = boundaryMesh();

            forAll(boundary, pI)
            {
                if (!isA<processorPolyPatch>(boundary[pI]))
                {
                    continue;
                }

                const processorPolyPatch& pp =
                (
                    refCast<const processorPolyPatch>(boundary[pI])
                );

                if (pp.neighbProcNo() == proc)
                {
                    procPatch = pI;
                    break;
                }
            }

            if (procPatch == -1)
            {
                Pout<< " * * * Sync Operations * * * " << nl
                    << " Could not find patch for proc: " << proc << nl
                    << abort(FatalError);
            }

            // Keep track of added entities from initial set
            label nEntities =
            (
                twoDMesh_ ?
                cMap.nEntities(coupleMap::FACE) :
                cMap.nEntities(coupleMap::EDGE)
            );

            // Specify a mapping for added indices
            Map<label> addedEntityMap;

            // Sequentially execute operations
            forAll(indices, indexI)
            {
                label index = indices[indexI], op = operations[indexI];

                if (debug > 3)
                {
                    Pout<< " Recv Op index: " << index << endl;
                }

                // Determine the appropriate local index
                label localIndex =
                (
                    entityMap.found(index) ?
                    entityMap[index] :
                    addedEntityMap[index]
                );

                changeMap opMap;

                switch (op)
                {
                    case coupleMap::BISECTION:
                    {
                        opMap = bisectEdge(localIndex);

                        label nIndex =
                        (
                            twoDMesh_ ?
                            opMap.addedFaceList()[2].index() :
                            opMap.addedEdgeList()[0].index()
                        );

                        // Insert the added index into the map
                        addedEntityMap.insert(nEntities++, nIndex);

                        if (debug > 3)
                        {
                            Pout<< " Adding Op index: " << nIndex
                                << " for index: " << localIndex
                                << endl;
                        }

                        break;
                    }

                    case coupleMap::COLLAPSE_FIRST:
                    {
                        opMap = collapseEdge(localIndex, 1);
                        break;
                    }

                    case coupleMap::COLLAPSE_SECOND:
                    {
                        opMap = collapseEdge(localIndex, 2);
                        break;
                    }

                    case coupleMap::COLLAPSE_MIDPOINT:
                    {
                        opMap = collapseEdge(localIndex, 3);
                        break;
                    }

                    case coupleMap::REMOVE_CELL:
                    {
                        // Clear existing list
                        rCellList.clear();

                        if (twoDMesh_)
                        {
                            // Insert the owner cell
                            rCellList.append(owner_[localIndex]);
                        }
                        else
                        {
                            // Insert all cells connected to this edge
                            const labelList& eFaces = edgeFaces_[localIndex];

                            forAll(eFaces, faceI)
                            {
                                const label own = owner_[eFaces[faceI]];
                                const label nei = neighbour_[eFaces[faceI]];

                                if (findIndex(rCellList, own) == -1)
                                {
                                    rCellList.append(own);
                                }

                                if (nei != -1)
                                {
                                    if (findIndex(rCellList, nei) == -1)
                                    {
                                        rCellList.append(nei);
                                    }
                                }
                            }
                        }

                        opMap = removeCells(rCellList, procPatch);

                        break;
                    }
                }

                if (opMap.type() <= 0)
                {
                    Pout<< " * * * Sync Operations * * * " << nl
                        << " Operation failed." << nl
                        << " Index: " << index << nl
                        << " localIndex: " << localIndex << nl
                        << " operation: " << op << nl
                        << " opMap.type: " << opMap.type() << nl
                        << endl;
                }
            }
        }
    }

    // Re-Initialize the stack with avoided entities from subMeshes
    // and leave those on processor patches untouched
    labelHashSet procEntities;

    buildEntitiesToAvoid(procEntities, false);

    // First remove processor entries
    forAllConstIter(labelHashSet, procEntities, pIter)
    {
        if (entities.found(pIter.key()))
        {
            entities.erase(pIter.key());
        }
    }

    // Initialize the coupled stack, using supplied entities
    initCoupledStack(entities, true);

    // Loop through the coupled stack and perform changes.
    if (edgeRefinement_)
    {
        edgeRefinementEngine(&(handlerPtr_[0]));
    }

    // Re-Initialize the stack, using supplied entities
    initCoupledStack(entities, true);

    if (twoDMesh_)
    {
        swap2DEdges(&(handlerPtr_[0]));
    }
    else
    {
        swap3DEdges(&(handlerPtr_[0]));
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
    Map<labelList> commonCells;

    forAll(procIndices_, pI)
    {
        label proc = procIndices_[pI];

        coupledPatchInfo& sPM = sendPatchMeshes_[pI];

        // Build the subMesh.
        buildProcessorPatchMesh
        (
            sPM,
            commonCells
        );

        const coupleMap& scMap = sPM.patchMap();

        // Send my sub-mesh to the neighbour.
        meshOps::pWrite(proc, scMap.nEntities());

        if (debug > 3)
        {
            Pout<< "Sending to [" << proc << "]:: nEntities: "
                << scMap.nEntities() << endl;
        }

        // Send the pointBuffers
        meshOps::pWrite(proc, scMap.pointBuffer());
        meshOps::pWrite(proc, scMap.oldPointBuffer());

        // Send connectivity (points, edges, faces, cells, etc)
        forAll(scMap.entityBuffer(), bufferI)
        {
            if (scMap.entityBuffer(bufferI).size())
            {
                meshOps::pWrite(proc, scMap.entityBuffer(bufferI));
            }
        }

        // Obtain references
        const coupledPatchInfo& rPM = recvPatchMeshes_[pI];
        const coupleMap& rcMap = rPM.patchMap();

        // First read entity sizes.
        meshOps::pRead(proc, rcMap.nEntities());

        if (debug > 3)
        {
            Pout<< "Receiving from [" << proc << "]:: nEntities: "
                << rcMap.nEntities() << endl;
        }

        // Size the buffers.
        rcMap.allocateBuffers();

        // Receive the pointBuffers
        meshOps::pRead(proc, rcMap.pointBuffer());
        meshOps::pRead(proc, rcMap.oldPointBuffer());

        // Receive connectivity (points, edges, faces, cells, etc)
        forAll(rcMap.entityBuffer(), bufferI)
        {
            if (rcMap.entityBuffer(bufferI).size())
            {
                meshOps::pRead(proc, rcMap.entityBuffer(bufferI));
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
    Map<labelList>& commonCells
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
    if (twoDMesh_)
    {
        // No pointEdges structure, so loop through all cells
        // to check for those connected to subMeshPoints
        forAll(cells_, cellI)
        {
            const cell& cellToCheck = cells_[cellI];
            const labelList cellPoints = cellToCheck.labels(faces_);

            forAll(cellPoints, pointI)
            {
                if (rPointMap.found(cellPoints[pointI]))
                {
                    if (commonCells.found(cellI))
                    {
                        // Add locally common cells at the end.
                        if (!localCommonCells.found(cellI))
                        {
                            localCommonCells.insert(cellI);
                        }

                        // Check if the processor exists on the list
                        // and if not, add it.
                        labelList& procList = commonCells[cellI];

                        if (findIndex(procList, proc) == -1)
                        {
                            meshOps::sizeUpList(proc, procList);
                        }
                    }
                    else
                    {
                        cellMap.insert(nC, cellI);
                        rCellMap.insert(cellI, nC);
                        nC++;

                        commonCells.insert(cellI, labelList(1, proc));
                    }

                    break;
                }
            }
        }
    }
    else
    {
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
                            if (!localCommonCells.found(own))
                            {
                                localCommonCells.insert(own);
                            }

                            // Check if the processor exists on the list
                            // and if not, add it.
                            labelList& procList = commonCells[own];

                            if (findIndex(procList, proc) == -1)
                            {
                                meshOps::sizeUpList(proc, procList);
                            }
                        }
                        else
                        {
                            cellMap.insert(nC, own);
                            rCellMap.insert(own, nC);
                            nC++;

                            commonCells.insert(own, labelList(1, proc));
                        }
                    }

                    // Check neighbour cell
                    if (!rCellMap.found(nei) && nei != -1)
                    {
                        if (commonCells.found(nei))
                        {
                            // Add locally common cells at the end.
                            if (!localCommonCells.found(nei))
                            {
                                localCommonCells.insert(nei);
                            }

                            // Check if the processor exists on the list
                            // and if not, add it.
                            labelList& procList = commonCells[nei];

                            if (findIndex(procList, proc) == -1)
                            {
                                meshOps::sizeUpList(proc, procList);
                            }
                        }
                        else
                        {
                            cellMap.insert(nC, nei);
                            rCellMap.insert(nei, nC);
                            nC++;

                            commonCells.insert(nei, labelList(1, proc));
                        }
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

    // Keep track of inserted boundary face indices
    // - Add exposed internal faces to a 'default' patch
    //   at the end of the list.
    labelList bdyFaceSizes(boundary.size() + 1, 0);
    labelList bdyFaceStarts(boundary.size() + 1, 0);
    List<Map<label> > bdyFaceIndices(boundary.size() + 1);

    // Allocate the faceMap. Since the number of internal faces is unknown,
    // detect internal ones first and update boundaries later.
    label sumNFE = 0;

    forAllConstIter(Map<label>, rCellMap, cIter)
    {
        const cell& thisCell = cells_[cIter.key()];

        forAll(thisCell, faceI)
        {
            label fIndex = thisCell[faceI];

            if (!rFaceMap.found(fIndex))
            {
                // Determine the patch index
                label patchID = whichPatch(fIndex);

                if (patchID == -1)
                {
                    // Internal face. Check if this needs to
                    // be added to the 'default' patch.
                    label own = owner_[fIndex];
                    label nei = neighbour_[fIndex];

                    if (rCellMap.found(own) && rCellMap.found(nei))
                    {
                        faceMap.insert(nF, fIndex);
                        rFaceMap.insert(fIndex, nF);
                        nF++;
                    }
                    else
                    {
                        // Update boundary maps.
                        label bfI = bdyFaceSizes[boundary.size()]++;

                        // Skip faceMap and update the reverseMap for now.
                        // faceMap will be updated once all
                        // internal faces have been detected.
                        rFaceMap.insert(fIndex, bfI);
                        bdyFaceIndices[boundary.size()].insert(bfI, fIndex);
                    }
                }
                else
                {
                    // Update boundary maps.
                    label bfI = bdyFaceSizes[patchID]++;

                    // Skip faceMap and update the reverseMap for now.
                    // faceMap will be updated once all
                    // internal faces have been detected.
                    rFaceMap.insert(fIndex, bfI);
                    bdyFaceIndices[patchID].insert(bfI, fIndex);
                }

                // Accumulate face sizes
                sumNFE += faces_[fIndex].size();
            }
        }
    }

    // Set the number of internal faces at this point
    cMap.nEntities(coupleMap::INTERNAL_FACE) = nF;

    // Set patch starts
    bdyFaceStarts[0] = nF;

    for (label i = 1; i < bdyFaceStarts.size(); i++)
    {
        bdyFaceStarts[i] = bdyFaceStarts[i-1] + bdyFaceSizes[i-1];
    }

    // Update faceMap and reverseFaceMap for boundaries
    forAll(bdyFaceIndices, patchI)
    {
        label pStart = bdyFaceStarts[patchI];

        forAllConstIter(Map<label>, bdyFaceIndices[patchI], fIter)
        {
            faceMap.insert(fIter.key() + pStart, fIter());
            rFaceMap[fIter()] = fIter.key() + pStart;
        }

        // Update face-count
        nF += bdyFaceSizes[patchI];
    }

    // Keep track of inserted boundary edge indices
    // - Add exposed internal edges to a 'default' patch
    //   at the end of the list.
    labelList bdyEdgeSizes(boundary.size() + 1, 0);
    labelList bdyEdgeStarts(boundary.size() + 1, 0);
    List<Map<label> > bdyEdgeIndices(boundary.size() + 1);

    // Allocate the edgeMap. Since the number of internal edges is unknown,
    // detect internal ones first and update boundaries later.
    forAllConstIter(Map<label>, rFaceMap, fIter)
    {
        const labelList& fEdges = faceEdges_[fIter.key()];

        forAll(fEdges, edgeI)
        {
            label eIndex = fEdges[edgeI];

            if (!rEdgeMap.found(eIndex))
            {
                // Determine the patch index
                label patchID = whichEdgePatch(eIndex);

                if (patchID == -1)
                {
                    bool boundaryEdge = false;

                    // Check if any cells touching edgeFaces
                    // do not belong to the cellMap.
                    const labelList& eFaces = edgeFaces_[eIndex];

                    forAll(eFaces, faceI)
                    {
                        label fIndex = eFaces[faceI];

                        label own = owner_[fIndex];
                        label nei = neighbour_[fIndex];

                        if (!rCellMap.found(own) || !rCellMap.found(nei))
                        {
                            boundaryEdge = true;
                            break;
                        }
                    }

                    if (boundaryEdge)
                    {
                        // Update boundary maps.
                        label beI = bdyEdgeSizes[boundary.size()]++;

                        // Skip edgeMap and update the reverseMap for now.
                        // edgeMap will be updated once all
                        // internal edges have been detected.
                        rEdgeMap.insert(eIndex, beI);
                        bdyEdgeIndices[boundary.size()].insert(beI, eIndex);
                    }
                    else
                    {
                        // Internal edge
                        edgeMap.insert(nE, eIndex);
                        rEdgeMap.insert(eIndex, nE);
                        nE++;
                    }
                }
                else
                {
                    // Update boundary maps.
                    label beI = bdyEdgeSizes[patchID]++;

                    // Skip edgeMap and update the reverseMap for now.
                    // edgeMap will be updated once all
                    // internal edges have been detected.
                    rEdgeMap.insert(eIndex, beI);
                    bdyEdgeIndices[patchID].insert(beI, eIndex);
                }
            }
        }
    }

    // Set the number of internal edges at this point
    cMap.nEntities(coupleMap::INTERNAL_EDGE) = nE;

    // Set patch starts
    bdyEdgeStarts[0] = nE;

    for (label i = 1; i < bdyEdgeStarts.size(); i++)
    {
        bdyEdgeStarts[i] = bdyEdgeStarts[i-1] + bdyEdgeSizes[i-1];
    }

    // Update edgeMap and reverseEdgeMap for boundaries
    forAll(bdyEdgeIndices, patchI)
    {
        label pStart = bdyEdgeStarts[patchI];

        forAllConstIter(Map<label>, bdyEdgeIndices[patchI], eIter)
        {
            edgeMap.insert(eIter.key() + pStart, eIter());
            rEdgeMap[eIter()] = eIter.key() + pStart;
        }

        // Update edge-count
        nE += bdyEdgeSizes[patchI];
    }

    // Set additional points in the pointMap
    forAllConstIter(Map<label>, rEdgeMap, eIter)
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
    cMap.nEntities(coupleMap::EDGE) = nE;
    cMap.nEntities(coupleMap::FACE) = nF;
    cMap.nEntities(coupleMap::CELL) = nC;
    cMap.nEntities(coupleMap::SHARED_POINT) = subMeshPoints.size();
    cMap.nEntities(coupleMap::NFE_SIZE) = sumNFE;
    cMap.nEntities(coupleMap::NBDY) = boundary.size() + 1;

    // Size up buffers and fill them
    cMap.allocateBuffers();

    pointField& pBuffer = cMap.pointBuffer();
    pointField& opBuffer = cMap.oldPointBuffer();

    forAllConstIter(Map<label>, pointMap, pIter)
    {
        pBuffer[pIter.key()] = points_[pIter()];
        opBuffer[pIter.key()] = oldPoints_[pIter()];
    }

    label index = 0;

    // Edge buffer size: 2 points for every edge
    labelList& eBuffer = cMap.entityBuffer(coupleMap::EDGE);

    for (label i = 0; i < nE; i++)
    {
        label eIndex = edgeMap[i];
        const edge& edgeToCheck = edges_[eIndex];

        eBuffer[index++] = rPointMap[edgeToCheck[0]];
        eBuffer[index++] = rPointMap[edgeToCheck[1]];
    }

    index = 0;
    face thisFace;
    labelList& fBuffer = cMap.entityBuffer(coupleMap::FACE);
    labelList& feBuffer = cMap.entityBuffer(coupleMap::FACE_EDGE);

    for (label i = 0; i < nF; i++)
    {
        label fIndex = faceMap[i];
        label own = owner_[fIndex];
        label nei = neighbour_[fIndex];

        if (rCellMap.found(own))
        {
            // Check if this face is pointed the right way
            if (rCellMap.found(nei) && (rCellMap[nei] < rCellMap[own]))
            {
                thisFace = faces_[fIndex].reverseFace();
            }
            else
            {
                thisFace = faces_[fIndex];
            }
        }
        else
        {
            // This face is pointed the wrong way.
            thisFace = faces_[fIndex].reverseFace();
        }

        const labelList& fEdges = faceEdges_[fIndex];

        forAll(fEdges, indexI)
        {
            fBuffer[index] = rPointMap[thisFace[indexI]];
            feBuffer[index] = rEdgeMap[fEdges[indexI]];

            index++;
        }
    }

    // Fill variable size face-list sizes for 2D
    if (twoDMesh_)
    {
        labelList& nfeBuffer = cMap.entityBuffer(coupleMap::NFE_BUFFER);

        forAllConstIter(Map<label>, faceMap, fIter)
        {
            nfeBuffer[fIter.key()] = faces_[fIter()].size();
        }
    }

    index = 0;
    labelList& cBuffer = cMap.entityBuffer(coupleMap::CELL);

    for (label i = 0; i < nC; i++)
    {
        label cIndex = cellMap[i];
        const cell& cellToCheck = cells_[cIndex];

        forAll(cellToCheck, faceI)
        {
            cBuffer[index++] = rFaceMap[cellToCheck[faceI]];
        }
    }

    // Fill in boundary information
    cMap.entityBuffer(coupleMap::FACE_STARTS) = bdyFaceStarts;
    cMap.entityBuffer(coupleMap::FACE_SIZES) = bdyFaceSizes;
    cMap.entityBuffer(coupleMap::EDGE_STARTS) = bdyEdgeStarts;
    cMap.entityBuffer(coupleMap::EDGE_SIZES) = bdyEdgeSizes;

    labelList& ptBuffer = cMap.entityBuffer(coupleMap::PATCH_ID);

    // Fill types for all but the last one (which is default).
    forAll(boundary, patchI)
    {
        if (isA<processorPolyPatch>(boundary[patchI]))
        {
            // Fetch the neighbouring processor ID
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[patchI])
            );

            // Set processor patches to a special type
            ptBuffer[patchI] = -2 - pp.neighbProcNo();
        }
        else
        {
            // Conventional physical patch. Make an identical
            // map, since physical boundaries are present on
            // all processors.
            ptBuffer[patchI] = patchI;
        }
    }

    // Fill the default patch as 'internal'
    ptBuffer[boundary.size()] = -1;

    // Set maps as built.
    subMesh.setBuiltMaps();

    // For debugging purposes...
    if (debug > 3)
    {
        Pout<< "Writing out coupledPatchInfo for processor: "
            << proc << endl;

        writeVTK
        (
            "psMesh_" + Foam::name(Pstream::myProcNo())
          + "to" + Foam::name(proc),
            rCellMap.toc()
        );

        // Write out patch information
        Pout<< " faceStarts: " << bdyFaceStarts << nl
            << " faceSizes: " << bdyFaceSizes << nl
            << " edgeStarts: " << bdyEdgeStarts << nl
            << " edgeSizes: " << bdyEdgeSizes << nl
            << " patchTypes: " << ptBuffer << endl;
    }
}


// Build coupled maps for locally coupled patches.
//   - Performs a geometric match initially, since the mesh provides
//     no explicit information for topological coupling.
//   - For each subsequent topology change, maps are updated during
//     the re-ordering stage.
void dynamicTopoFvMesh::buildLocalCoupledMaps()
{
    if (patchCoupling_.empty())
    {
        return;
    }

    if (debug)
    {
        Info<< "Building local coupled maps...";
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
            FatalErrorIn("void dynamicTopoFvMesh::buildLocalCoupledMaps()")
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
                    FatalErrorIn
                    (
                        "void dynamicTopoFvMesh::buildLocalCoupledMaps()"
                    )
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

                if (twoDMesh_)
                {
                    // Set up a comparison face.
                    face cFace(4);

                    // Configure the face for comparison.
                    forAll(mFace, pointI)
                    {
                        cFace[pointI] = pMap[mFace[pointI]];
                    }

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
                        FatalErrorIn
                        (
                            "void dynamicTopoFvMesh::buildLocalCoupledMaps()"
                        )
                            << " Failed to match face: "
                            << mfIndex << ": " << mFace << nl
                            << " Comparison face: " << cFace
                            << abort(FatalError);
                    }
                }
                else
                {
                    label slaveFaceIndex = -1;

                    // Set up a comparison face.
                    triFace cFace
                    (
                        pMap[mFace[0]],
                        pMap[mFace[1]],
                        pMap[mFace[2]]
                    );

                    forAll(spFaces, faceJ)
                    {
                        // Fetch the global face index
                        label sfIndex = (sStart + spFaces[faceJ]);

                        const face& sFace = faces_[sfIndex];

                        if (triFace::compare(cFace, triFace(sFace)))
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
                        FatalErrorIn
                        (
                            "void dynamicTopoFvMesh::buildLocalCoupledMaps()"
                        )
                            << " Failed to match face: "
                            << mfIndex << ": " << mFace << nl
                            << " Comparison face: " << cFace
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
                            FatalErrorIn
                            (
                                "void dynamicTopoFvMesh::"
                                "buildLocalCoupledMaps()"
                            )
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
        Info<< "Done." << endl;
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
    meshOps::waitForBuffers();

    // Put un-matched faces in a list.
    labelHashSet unMatchedFaces;

    forAll(procIndices_, pI)
    {
        label proc = procIndices_[pI];

        coupledPatchInfo& rPM = recvPatchMeshes_[pI];
        const coupleMap& cMap = rPM.patchMap();
        const labelList& ptBuffer = cMap.entityBuffer(coupleMap::PATCH_ID);

        // Specify the list of patch names and types
        wordList patchNames(ptBuffer.size());
        wordList patchTypes(ptBuffer.size());

        // Specify additional info for processor types
        labelList sProcNo(ptBuffer.size(), -1);
        labelList nProcNo(ptBuffer.size(), -1);

        forAll(patchTypes, pI)
        {
            if (pI == patchTypes.size() - 1)
            {
                patchNames[pI] = "defaultPatch";
                patchTypes[pI] = "patch";
            }
            else
            if (ptBuffer[pI] <= -2)
            {
                // Back out the neighbouring processor ID
                label neiProcNo = Foam::mag(ptBuffer[pI] + 2);

                patchNames[pI] =
                (
                    "procBoundary"
                  + Foam::name(proc) + "to"
                  + Foam::name(neiProcNo)
                );

                // Fill in additional patch information
                sProcNo[pI] = proc;
                nProcNo[pI] = neiProcNo;

                patchTypes[pI] = "processor";
            }
            else
            {
                patchNames[pI] = boundary[ptBuffer[pI]].name();
                patchTypes[pI] = boundary[ptBuffer[pI]].type();
            }
        }

        // Set the autoPtr.
        rPM.setMesh
        (
            proc,
            new dynamicTopoFvMesh
            (
                (*this),
                IOobject
                (
                    fvMesh::defaultRegion,
                    time().timeName(),
                    time()
                ),
                xferCopy(cMap.pointBuffer()),
                xferCopy(cMap.oldPointBuffer()),
                xferCopy(cMap.edges()),
                xferCopy(cMap.faces()),
                xferCopy(cMap.faceEdges()),
                xferCopy(cMap.owner()),
                xferCopy(cMap.neighbour()),
                xferCopy(cMap.cells()),
                cMap.entityBuffer(coupleMap::FACE_STARTS),
                cMap.entityBuffer(coupleMap::FACE_SIZES),
                cMap.entityBuffer(coupleMap::EDGE_STARTS),
                cMap.entityBuffer(coupleMap::EDGE_SIZES),
                patchNames,
                patchTypes,
                sProcNo,
                nProcNo
            )
        );

        // Sanity check: Do sub-mesh point sizes match?
        if
        (
            cMap.subMeshPoints().size()
         != cMap.nEntities(coupleMap::SHARED_POINT)
        )
        {
            FatalErrorIn("void dynamicTopoFvMesh::buildProcessorCoupledMaps()")
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
                refCast<const processorPolyPatch>(boundary[nPrc[proc]])
            );

            const labelList& neiPoints = pp.neighbPoints();

            if (debug)
            {
                if (findIndex(neiPoints, -1) > -1)
                {
                    FatalErrorIn
                    (
                        "void dynamicTopoFvMesh::buildProcessorCoupledMaps()"
                    )
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
            const pointField& sPoints = rPM.subMesh().points();

            forAllConstIter(Map<label>, pMap, pIter)
            {
                scalar dist = mag(points_[pIter.key()] - sPoints[pIter()]);

                if (dist > debug::tolerances("processorMatchTol"))
                {
                    FatalErrorIn
                    (
                        "void dynamicTopoFvMesh::buildProcessorCoupledMaps()"
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
            const faceList& slaveFaces = rPM.subMesh().faces();

            // This is an immediate neighbour.
            label mStart = boundary[nPrc[proc]].start();
            label mSize  = boundary[nPrc[proc]].size();

            // Abbreviate for convenience
            const dynamicTopoFvMesh& sMesh = rPM.subMesh();
            const Map<label>& pMap = cMap.entityMap(coupleMap::POINT);
            const Map<label>& eMap = cMap.entityMap(coupleMap::EDGE);
            const Map<label>& fMap = cMap.entityMap(coupleMap::FACE);

            // Fetch global pointFaces for the slave.
            const labelListList& spF = rPM.subMesh().pointFaces();

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
                    label sFaceIndex = -1;

                    forAll(spFaces, faceJ)
                    {
                        // Fetch the global face index
                        label sfIndex = spFaces[faceJ];

                        const face& sFace = slaveFaces[sfIndex];

                        if
                        (
                            triFace::compare
                            (
                                triFace(cFace[0], cFace[1], cFace[2]),
                                triFace(sFace[0], sFace[1], sFace[2])
                            )
                        )
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

                            sFaceIndex = sfIndex;

                            break;
                        }
                    }

                    if (sFaceIndex == -1)
                    {
                        unMatchedFaces.insert(mfIndex);

                        continue;
                    }

                    // Match all edges on this face as well.
                    const labelList& mfEdges = faceEdges_[mfIndex];
                    const labelList& sfEdges = sMesh.faceEdges_[sFaceIndex];

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

                            FatalErrorIn
                            (
                                "void dynamicTopoFvMesh::"
                                "buildProcessorCoupledMaps()"
                            )
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
            notImplemented
            (
                "void dynamicTopoFvMesh::buildProcessorCoupledMaps()"
            );
        }

        if (unMatchedFaces.size())
        {
            // Write out the face
            writeVTK("mFaces", unMatchedFaces.toc(), 2);

            FatalErrorIn
            (
                "void dynamicTopoFvMesh::buildProcessorCoupledMaps()"
            )
                << " Unmatched faces were found for processor: " << proc
                << abort(FatalError);
        }

        // If the received mesh is of higher rank,
        // build associated edgePoints.
        if (proc > Pstream::myProcNo() && !twoDMesh_)
        {
            dynamicTopoFvMesh& mesh = rPM.subMesh();

            forAll(mesh.edges_, edgeI)
            {
                mesh.buildEdgePoints(edgeI);
            }
        }
    }
}


// Initialize coupled boundary ordering
// - Assumes that faces_ and points_ are consistent
// - Assumes that patchStarts_ and patchSizes_ are consistent
void dynamicTopoFvMesh::initCoupledBoundaryOrdering
(
    List<pointField>& centres,
    List<pointField>& anchors
) const
{
    if (!Pstream::parRun())
    {
        return;
    }

    const polyBoundaryMesh& boundary = boundaryMesh();

    forAll(boundary, pI)
    {
        if (isA<processorPolyPatch>(boundary[pI]))
        {
            // Check if this is a master processor patch.
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[pI])
            );

            label start = patchStarts_[pI];
            label size = patchSizes_[pI];

            // Prepare centres and anchors
            centres[pI].setSize(size, vector::zero);
            anchors[pI].setSize(size, vector::zero);

            if (Pstream::myProcNo() < pp.neighbProcNo())
            {
                forAll(centres[pI], fI)
                {
                    centres[pI][fI] = faces_[fI + start].centre(points_);
                    anchors[pI][fI] = points_[faces_[fI + start][0]];
                }

                if (debug)
                {
                    if (debug > 3)
                    {
                        Pout<< "Sending to: " << pp.neighbProcNo()
                            << " nCentres: " << size << endl;

                        // Write out my centres to disk
                        meshOps::writeVTK
                        (
                            (*this),
                            "centres_" + pp.name(),
                            size, size, size,
                            centres[pI]
                        );
                    }

                    // Ensure that we're sending the right size
                    meshOps::pWrite(pp.neighbProcNo(), size);
                }

                // Send information to neighbour
                meshOps::pWrite(pp.neighbProcNo(), centres[pI]);
                meshOps::pWrite(pp.neighbProcNo(), anchors[pI]);
            }
            else
            {
                if (debug)
                {
                    label nEntities = -1;

                    // Ensure that we're receiving the right size
                    meshOps::pRead(pp.neighbProcNo(), nEntities);

                    if (debug > 3)
                    {
                        Pout<< " Recving from: " << pp.neighbProcNo()
                            << " nCentres: " << size << endl;
                    }

                    if (nEntities != size)
                    {
                        FatalErrorIn
                        (
                            "void dynamicTopoFvMesh::"
                            "initCoupledBoundaryOrdering() const"
                        )
                            << "Incorrect send / recv sizes: " << nl
                            << " nEntities: " << nEntities << nl
                            << " size: " << size << nl
                            << abort(FatalError);
                    }
                }

                // Schedule receive from neighbour
                meshOps::pRead(pp.neighbProcNo(), centres[pI]);
                meshOps::pRead(pp.neighbProcNo(), anchors[pI]);
            }
        }
    }
}


// Synchronize coupled boundary ordering
bool dynamicTopoFvMesh::syncCoupledBoundaryOrdering
(
    List<pointField>& centres,
    List<pointField>& anchors,
    labelListList& patchMaps,
    labelListList& rotations
) const
{
    if (!Pstream::parRun())
    {
        return false;
    }

    bool anyChange = false;

    const polyBoundaryMesh& boundary = boundaryMesh();

    // Calculate centres and tolerances for any slave patches
    List<scalarField> slaveTols(boundary.size());
    List<pointField> slaveCentres(boundary.size());

    scalar matchTol = Foam::debug::tolerances("meshOpsMatchTol", 1e-4);

    forAll(boundary, pI)
    {
        if (isA<processorPolyPatch>(boundary[pI]))
        {
            // Check if this is a slave processor patch.
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[pI])
            );

            label start = patchStarts_[pI];
            label size = patchSizes_[pI];

            if (Pstream::myProcNo() > pp.neighbProcNo())
            {
                slaveTols[pI].setSize(size, 0.0);
                slaveCentres[pI].setSize(size, vector::zero);

                forAll(slaveCentres[pI], fI)
                {
                    point& fc = slaveCentres[pI][fI];

                    const face& checkFace = faces_[fI + start];

                    // Calculate centre
                    fc = checkFace.centre(points_);

                    scalar maxLen = -GREAT;

                    forAll(checkFace, fpI)
                    {
                        maxLen = max(maxLen, mag(points_[checkFace[fpI]] - fc));
                    }

                    slaveTols[pI][fI] = matchTol*maxLen;
                }

                // Write out my centres to disk
                if (debug > 3)
                {
                    meshOps::writeVTK
                    (
                        (*this),
                        "slaveCentres_" + pp.name(),
                        size, size, size,
                        slaveCentres[pI]
                    );
                }
            }
        }
    }

    // Wait for transfers before continuing.
    meshOps::waitForBuffers();

    forAll(boundary, pI)
    {
        if (isA<processorPolyPatch>(boundary[pI]))
        {
            // Check if this is a master processor patch.
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[pI])
            );

            labelList& patchMap = patchMaps[pI];
            labelList& rotation = rotations[pI];

            // Initialize map and rotation
            patchMap.setSize(patchSizes_[pI], -1);
            rotation.setSize(patchSizes_[pI], 0);

            if (Pstream::myProcNo() < pp.neighbProcNo())
            {
                // Do nothing (i.e. identical mapping, zero rotation).
                forAll(patchMap, pfI)
                {
                    patchMap[pfI] = pfI;
                }
            }
            else
            {
                // Try zero separation automatic matching
                bool matchedAll =
                (
                    matchPoints
                    (
                        slaveCentres[pI],
                        centres[pI],
                        slaveTols[pI],
                        true,
                        patchMap
                    )
                );

                // Write out master centres to disk
                if (debug > 3 || !matchedAll)
                {
                    label mSize = centres[pI].size();
                    label sSize = slaveCentres[pI].size();

                    meshOps::writeVTK
                    (
                        (*this),
                        "masterCentres_" + pp.name(),
                        mSize, mSize, mSize,
                        centres[pI]
                    );

                    meshOps::writeVTK
                    (
                        (*this),
                        "slaveCentres_" + pp.name(),
                        sSize, sSize, sSize,
                        slaveCentres[pI]
                    );

                    if (!matchedAll)
                    {
                        Pout<< " Patch: " << pp.name()
                            << " mSize: " << mSize << " sSize: " << sSize
                            << " failed on match for face centres."
                            << endl;
                    }
                }

                label start = patchStarts_[pI];

                // Set rotation.
                forAll(patchMap, oldFaceI)
                {
                    label newFaceI = patchMap[oldFaceI];

                    const point& anchor = anchors[pI][newFaceI];
                    const scalar& faceTol = slaveTols[pI][oldFaceI];
                    const face& checkFace = faces_[start + oldFaceI];

                    label anchorFp = -1;
                    scalar minDSqr = GREAT;

                    forAll(checkFace, fpI)
                    {
                        scalar dSqr = magSqr(anchor - points_[checkFace[fpI]]);

                        if (dSqr < minDSqr)
                        {
                            minDSqr = dSqr;
                            anchorFp = fpI;
                        }
                    }

                    if (anchorFp == -1 || mag(minDSqr) > faceTol)
                    {
                        FatalErrorIn
                        (
                            "void dynamicTopoFvMesh::"
                            "syncCoupledBoundaryOrdering() const"
                        )
                            << "Cannot find anchor: " << anchor << nl
                            << " Face: " << checkFace << nl
                            << " Vertices: "
                            << UIndirectList<point>(points_, checkFace) << nl
                            << " on patch: " << pp.name()
                            << abort(FatalError);
                    }
                    else
                    {
                        // Positive rotation
                        rotation[newFaceI] =
                        (
                            (checkFace.size() - anchorFp) % checkFace.size()
                        );
                    }
                }

                // Set the flag
                anyChange = true;
            }
        }
    }

    return anyChange;
}


// Fill buffers with length-scale info
// and exchange across processors.
void dynamicTopoFvMesh::exchangeLengthBuffers()
{
    if (!Pstream::parRun())
    {
        return;
    }

    forAll(procIndices_, pI)
    {
        coupledPatchInfo& sPM = sendPatchMeshes_[pI];
        coupledPatchInfo& rPM = recvPatchMeshes_[pI];

        const coupleMap& scMap = sPM.patchMap();
        const coupleMap& rcMap = rPM.patchMap();

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

            meshOps::pWrite(scMap.masterIndex(), sendLength);

            if (debug > 4)
            {
                Pout<< "Sending to: "
                    << scMap.masterIndex()
                    << " nCells: "
                    << scMap.nEntities(coupleMap::CELL)
                    << endl;
            }
        }

        if (rcMap.masterIndex() == Pstream::myProcNo())
        {
            // Schedule receipt from neighbour
            rPM.subMesh().lengthScale_.setSize
            (
                rcMap.nEntities(coupleMap::CELL),
                0.0
            );

            // Schedule for receipt
            meshOps::pRead
            (
                rcMap.slaveIndex(),
                rPM.subMesh().lengthScale_
            );

            if (debug > 4)
            {
                Pout<< "Receiving from: "
                    << rcMap.slaveIndex()
                    << " nCells: "
                    << rcMap.nEntities(coupleMap::CELL)
                    << endl;
            }
        }
    }

    // Wait for transfers before continuing.
    meshOps::waitForBuffers();

    if (debug > 4)
    {
        forAll(procIndices_, pI)
        {
            const coupledPatchInfo& rPM = recvPatchMeshes_[pI];
            const coupleMap& rcMap = rPM.patchMap();

            if (rcMap.masterIndex() == Pstream::myProcNo())
            {
                rPM.subMesh().writeVTK
                (
                    "lengthScale_" + Foam::name(rcMap.slaveIndex()),
                    identity(rPM.subMesh().nCells()),
                    3,
                    false,
                    false,
                    rPM.subMesh().lengthScale_
                );
            }
        }
    }
}


scalar dynamicTopoFvMesh::processorLengthScale(const label index) const
{
    scalar procScale = 0.0;

    if (twoDMesh_)
    {
        // First check the master processor
        procScale += lengthScale_[owner_[index]];

        // Next, check the slave processor
        bool foundSlave = false;

        forAll(procIndices_, pI)
        {
            // Fetch non-const reference to subMeshes
            const label faceEnum = coupleMap::FACE;
            const coupledPatchInfo& rPM = recvPatchMeshes_[pI];
            const coupleMap& cMap = rPM.patchMap();

            label sIndex = -1;

            if ((sIndex = cMap.findSlave(faceEnum, index)) > -1)
            {
                procScale +=
                (
                    rPM.subMesh().lengthScale_[rPM.subMesh().owner_[sIndex]]
                );

                foundSlave = true;
                break;
            }
        }

        // Should have found at least one slave
        if (!foundSlave)
        {
            FatalErrorIn
            (
                "scalar dynamicTopoFvMesh::processorLengthScale"
                "(const label index) const"
            )
                << "Processor lengthScale lookup failed: " << nl
                << " Master face: " << index
                << " :: " << faces_[index] << nl
                << abort(FatalError);
        }

        // Average the scale
        procScale *= 0.5;
    }
    else
    {
        const label edgeEnum = coupleMap::EDGE;
        const labelList& eFaces = edgeFaces_[index];

        // First check the master processor
        label nC = 0;

        forAll(eFaces, faceI)
        {
            label own = owner_[eFaces[faceI]];
            label nei = neighbour_[eFaces[faceI]];

            procScale += lengthScale_[own];
            nC++;

            if (nei > -1)
            {
                procScale += lengthScale_[nei];
                nC++;
            }
        }

        // Next check slaves
        bool foundSlave = false;

        forAll(procIndices_, pI)
        {
            const coupledPatchInfo& rPM = recvPatchMeshes_[pI];
            const coupleMap& cMap = rPM.patchMap();

            label sIndex = -1;

            if ((sIndex = cMap.findSlave(edgeEnum, index)) > -1)
            {
                // Fetch connectivity from patchSubMesh
                const labelList& peFaces =
                (
                    rPM.subMesh().edgeFaces_[sIndex]
                );

                foundSlave = true;

                forAll(peFaces, faceI)
                {
                    label own = rPM.subMesh().owner_[peFaces[faceI]];
                    label nei = rPM.subMesh().neighbour_[peFaces[faceI]];

                    procScale += lengthScale_[own];
                    nC++;

                    if (nei > -1)
                    {
                        procScale += lengthScale_[nei];
                        nC++;
                    }
                }
            }
        }

        // Should have found at least one slave
        if (!foundSlave)
        {
            FatalErrorIn
            (
                "scalar dynamicTopoFvMesh::processorLengthScale"
                "(const label index) const"
            )
                << "Processor lengthScale lookup failed: " << nl
                << " Master edge: " << index
                << " :: " << edges_[index] << nl
                << abort(FatalError);
        }

        // Average the final scale
        procScale /= nC;
    }

    return procScale;
}


// Method to determine whether the master face is locally coupled
bool dynamicTopoFvMesh::locallyCoupledEntity
(
    const label index,
    bool checkSlaves,
    bool checkProcs,
    bool checkFace
) const
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    // Bail out if no patchCoupling is present
    if (patchCoupling_.empty())
    {
        return false;
    }

    if (twoDMesh_ || checkFace)
    {
        label patch = whichPatch(index);

        if (patch == -1)
        {
            return false;
        }

        // Processor checks receive priority.
        if (isA<processorPolyPatch>(boundary[patch]) && checkProcs)
        {
            return false;
        }

        // Check coupled master patches.
        if (patchCoupling_(patch))
        {
            return true;
        }
        else
        if (checkSlaves)
        {
            // Check on slave patches as well.
            forAll(patchCoupling_, pI)
            {
                if (patchCoupling_(pI))
                {
                    const coupleMap& cMap = patchCoupling_[pI].patchMap();

                    if (cMap.slaveIndex() == patch)
                    {
                        return true;
                    }
                }
            }
        }
    }
    else
    {
        const labelList& eFaces = edgeFaces_[index];

        // Search for boundary faces, and determine boundary type.
        forAll(eFaces, faceI)
        {
            if (neighbour_[eFaces[faceI]] == -1)
            {
                label patch = whichPatch(eFaces[faceI]);

                // Processor checks receive priority.
                if (isA<processorPolyPatch>(boundary[patch]) && checkProcs)
                {
                    return false;
                }

                // Check coupled master patches.
                if (patchCoupling_(patch))
                {
                    return true;
                }

                if (checkSlaves)
                {
                    // Check on slave patches as well.
                    forAll(patchCoupling_, pI)
                    {
                        if (patchCoupling_(pI))
                        {
                            const coupleMap& cMap =
                            (
                                patchCoupling_[pI].patchMap()
                            );

                            if (cMap.slaveIndex() == patch)
                            {
                                return true;
                            }
                        }
                    }
                }
            }
        }
    }

    // Could not find any faces on locally coupled patches.
    return false;
}


// Method to determine the locally coupled patch index
label dynamicTopoFvMesh::locallyCoupledEdgePatch(const label eIndex) const
{
    const labelList& eFaces = edgeFaces_[eIndex];

    // Search for boundary faces, and determine boundary type.
    forAll(eFaces, faceI)
    {
        if (neighbour_[eFaces[faceI]] == -1)
        {
            label patch = whichPatch(eFaces[faceI]);

            // Check coupled master patches.
            if (patchCoupling_(patch))
            {
                return patch;
            }

            // Check on slave patches as well.
            forAll(patchCoupling_, pI)
            {
                if (patchCoupling_(pI))
                {
                    const coupleMap& cMap = patchCoupling_[pI].patchMap();

                    if (cMap.slaveIndex() == patch)
                    {
                        return patch;
                    }
                }
            }
        }
    }

    // Could not find any faces on locally coupled patches.
    FatalErrorIn
    (
        "label dynamicTopoFvMesh::locallyCoupledEdgePatch"
        "(const label cIndex) const"
    )
        << "Edge: " << eIndex << ":: " << edges_[eIndex]
        << " does not lie on any coupled patches."
        << abort(FatalError);

    return -1;
}


// Method to determine if the entity is on a processor boundary
//  - Also provides an additional check for 'pure' processor edges
//    i.e., edges that do not abut a physical patch. This is necessary
//    while deciding on collapse cases towards bounding curves.
bool dynamicTopoFvMesh::processorCoupledEntity
(
    const label index,
    bool checkFace,
    bool checkEdge,
    bool checkPure
) const
{
    // Skip check for serial runs
    if (!Pstream::parRun())
    {
        return false;
    }

    const polyBoundaryMesh& boundary = boundaryMesh();

    label patch = -2;

    if ((twoDMesh_ || checkFace) && !checkEdge)
    {
        patch = whichPatch(index);

        if (patch == -1)
        {
            return false;
        }

        if (isA<processorPolyPatch>(boundary[patch]))
        {
            return true;
        }
    }
    else
    {
        const labelList& eFaces = edgeFaces_[index];

        label nPhysical = 0, nProcessor = 0;

        // Search for boundary faces, and determine boundary type.
        forAll(eFaces, faceI)
        {
            if (neighbour_[eFaces[faceI]] == -1)
            {
                label patch = whichPatch(eFaces[faceI]);

                if (isA<processorPolyPatch>(boundary[patch]))
                {
                    // Increment the processor patch count
                    nProcessor++;

                    if (!checkPure)
                    {
                        // We don't have to validate if this
                        // is a 'pure' processor edge, so bail out.
                        return true;
                    }
                }
                else
                {
                    // Not a processor patch.
                    nPhysical++;
                }
            }
        }

        // Purity check
        if (checkPure)
        {
            if (nProcessor && !nPhysical)
            {
                return true;
            }
        }
    }

    // Could not find any faces on processor patches.
    return false;
}


// Build a list of entities that need to be avoided
// by regular topo-changes.
void dynamicTopoFvMesh::buildEntitiesToAvoid
(
    labelHashSet& entities,
    bool checkSubMesh
)
{
    entities.clear();

    // Build a set of entities to avoid during regular modifications,
    // and build a master stack for coupled modifications.
    const polyBoundaryMesh& boundary = boundaryMesh();

    // Determine locally coupled slave patches.
    labelHashSet localMasterPatches, localSlavePatches;

    forAll(patchCoupling_, patchI)
    {
        if (patchCoupling_(patchI))
        {
            const coupleMap& cMap = patchCoupling_[patchI].patchMap();

            localMasterPatches.insert(cMap.masterIndex());
            localSlavePatches.insert(cMap.slaveIndex());
        }
    }

    // Loop through boundary faces and check whether
    // they belong to master/slave coupled patches.
    for (label faceI = nOldInternalFaces_; faceI < faces_.size(); faceI++)
    {
        // Add only valid faces
        if (faces_[faceI].empty())
        {
            continue;
        }

        label pIndex = whichPatch(faceI);

        if (pIndex == -1)
        {
            continue;
        }

        // Check if this is a coupled face.
        if
        (
            localMasterPatches.found(pIndex) ||
            localSlavePatches.found(pIndex) ||
            isA<processorPolyPatch>(boundary[pIndex])
        )
        {
            if (twoDMesh_)
            {
                // Avoid this face during regular modification.
                if (!entities.found(faceI))
                {
                    entities.insert(faceI);
                }
            }
            else
            {
                const labelList& fEdges = faceEdges_[faceI];

                forAll(fEdges, edgeI)
                {
                    // Avoid this edge during regular modification.
                    if (!entities.found(fEdges[edgeI]))
                    {
                        entities.insert(fEdges[edgeI]);
                    }
                }
            }
        }
    }

    // Loop through entities contained in patchSubMeshes, if requested
    if (checkSubMesh)
    {
        forAll(procIndices_, pI)
        {
            const coupleMap& cMap = sendPatchMeshes_[pI].patchMap();
            const Map<label> rEdgeMap = cMap.reverseEntityMap(coupleMap::EDGE);

            if (cMap.slaveIndex() == Pstream::myProcNo())
            {
                forAllConstIter(Map<label>, rEdgeMap, eIter)
                {
                    if (twoDMesh_)
                    {
                        const labelList& eFaces = edgeFaces_[eIter.key()];

                        forAll(eFaces, faceI)
                        {
                            if (faces_[eFaces[faceI]].size() == 4)
                            {
                                if (!entities.found(eFaces[faceI]))
                                {
                                    entities.insert(eFaces[faceI]);
                                }
                            }
                        }
                    }
                    else
                    {
                        const edge& check = edges_[eIter.key()];
                        const labelList& eFaces = edgeFaces_[eIter.key()];

                        // Skip deleted edges
                        if (eFaces.size())
                        {
                            forAll(check, pI)
                            {
                                const labelList& pE = pointEdges_[check[pI]];

                                forAll(pE, edgeI)
                                {
                                    if (!entities.found(pE[edgeI]))
                                    {
                                        entities.insert(pE[edgeI]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (debug > 3)
    {
        Pout<< nl << "nEntitiesToAvoid: " << entities.size() << endl;

        if (debug > 4)
        {
            // Write out entities
            label elemType = twoDMesh_ ? 2 : 1;

            writeVTK
            (
                "entitiesToAvoid_"
              + Foam::name(Pstream::myProcNo()),
                entities.toc(),
                elemType
            );
        }
    }
}


// Check whether the specified edge is a coupled master edge.
bool dynamicTopoFvMesh::isCoupledMaster
(
    const label eIndex
) const
{
    if (!coupledModification_)
    {
        return true;
    }

    return locallyCoupledEntity(eIndex);
}


} // End namespace Foam

// ************************************************************************* //
