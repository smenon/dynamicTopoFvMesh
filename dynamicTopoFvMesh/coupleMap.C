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
    coupleMap

Description
    Implementation of the coupleMap class

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "coupleMap.H"
#include "boolList.H"
#include "demandDrivenData.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline coupleMap::coupleMap
(
    const IOobject& io,
    const bool twoDMesh,
    const bool isLocal,
    const bool isSend,
    const label patchIndex,
    const label masterIndex,
    const label slaveIndex
)
:
    regIOobject(io),
    twoDMesh_(twoDMesh),
    isLocal_(isLocal),
    isSend_(isSend),
    patchIndex_(patchIndex),
    masterIndex_(masterIndex),
    slaveIndex_(slaveIndex),
    nEntities_(-1),
    nInternalFaces_(-1),
    ownerPtr_(NULL),
    neighbourPtr_(NULL),
    edgesPtr_(NULL),
    facesPtr_(NULL),
    cellsPtr_(NULL),
    faceEdgesPtr_(NULL)
{
    if
    (
        (io.readOpt() == IOobject::MUST_READ)
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        // Construct an Istream and read from disk.
        readData(readStream(typeName));
        close();
    }
}


// Construct as copy
inline coupleMap::coupleMap(const coupleMap& cm)
:
    regIOobject(cm, true),
    twoDMesh_(cm.twoDMesh_),
    isLocal_(cm.isLocal_),
    isSend_(cm.isSend_),
    patchIndex_(cm.patchIndex_),
    masterIndex_(cm.masterIndex_),
    slaveIndex_(cm.slaveIndex_),
    nEntities_(cm.nEntities_),
    nInternalFaces_(-1),
    ownerPtr_(NULL),
    neighbourPtr_(NULL),
    edgesPtr_(NULL),
    facesPtr_(NULL),
    cellsPtr_(NULL),
    faceEdgesPtr_(NULL)
{
    if
    (
        (cm.readOpt() == IOobject::MUST_READ)
     || (cm.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        // Construct an Istream and read from disk.
        readData(readStream(typeName));
        close();
    }
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

inline coupleMap::~coupleMap()
{
    clearMaps();
    clearBuffers();
    clearAddressing();
}


// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

inline void coupleMap::clearAddressing() const
{
    nInternalFaces_ = -1;
    deleteDemandDrivenData(ownerPtr_);
    deleteDemandDrivenData(neighbourPtr_);
    deleteDemandDrivenData(edgesPtr_);
    deleteDemandDrivenData(facesPtr_);
    deleteDemandDrivenData(cellsPtr_);
    deleteDemandDrivenData(faceEdgesPtr_);
}


inline void coupleMap::makeAddressing() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ownerPtr_ || neighbourPtr_ || nInternalFaces_ > -1)
    {
        FatalErrorIn("inline void coupleMap::makeAddressing() const")
            << "Addressing has already been calculated."
            << abort(FatalError);
    }

    label nFaces = nEntities(coupleMap::FACE);

    if (nFaces < 0)
    {
        FatalErrorIn("inline void coupleMap::makeAddressing() const")
            << "Invalid buffers. Cannot continue."
            << abort(FatalError);
    }

    // Set sizes.
    ownerPtr_ = new labelList(nFaces, -1);
    neighbourPtr_ = new labelList(nFaces, -1);

    labelList& own = *ownerPtr_;
    labelList& nei = *neighbourPtr_;

    boolList markedFaces(nFaces, false);

    label nInternalFaces_ = 0;

    // Fetch the demand-driven cell list.
    const cellList& cList = cells();

    forAll(cList, cellI)
    {
        const cell& cellToCheck = cList[cellI];

        forAll(cellToCheck, faceI)
        {
            if (!markedFaces[cellToCheck[faceI]])
            {
                // First visit: owner
                own[cellToCheck[faceI]] = cellI;
                markedFaces[cellToCheck[faceI]] = true;
            }
            else
            {
                // Second visit: neighbour
                nei[cellToCheck[faceI]] = cellI;
                nInternalFaces_++;
            }
        }
    }

    nei.setSize(nInternalFaces_);
}


inline void coupleMap::makeEdges() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (edgesPtr_)
    {
        FatalErrorIn("inline void coupleMap::makeEdges() const")
            << "Edges have already been calculated."
            << abort(FatalError);
    }

    label nEdges = nEntities(coupleMap::EDGE);
    const labelList& eBuffer = entityBuffer(coupleMap::EDGE);

    edgesPtr_ = new edgeList(nEdges);

    edgeList& edges = *edgesPtr_;

    forAll(edges, edgeI)
    {
        edges[edgeI][0] = eBuffer[(2*edgeI)+0];
        edges[edgeI][1] = eBuffer[(2*edgeI)+1];
    }
}


inline void coupleMap::makeFaces() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (facesPtr_ || faceEdgesPtr_)
    {
        FatalErrorIn("inline void coupleMap::makeFaces() const")
            << "Faces have already been calculated."
            << abort(FatalError);
    }

    label nFaces = nEntities(coupleMap::FACE);

    const labelList& fBuffer = entityBuffer(coupleMap::FACE);
    const labelList& feBuffer = entityBuffer(coupleMap::FACE_EDGE);
    const labelList& nfeBuffer = entityBuffer(coupleMap::NFE_BUFFER);

    facesPtr_ = new faceList(nFaces);
    faceEdgesPtr_ = new labelListList(nFaces);

    faceList& faces = *facesPtr_;
    labelListList& faceEdges = *faceEdgesPtr_;

    label sumNFE = 0;

    forAll(faces, faceI)
    {
        face& f = faces[faceI];
        labelList& fe = faceEdges[faceI];

        // Fetch the buffer value for 2D meshes
        label nfe = twoDMesh_ ? nfeBuffer[faceI] : 3;

        // Size up the lists
        f.setSize(nfe, -1);
        fe.setSize(nfe, -1);

        for (label p = 0; p < nfe; p++)
        {
            f[p] = fBuffer[sumNFE + p];
            fe[p] = feBuffer[sumNFE + p];
        }

        sumNFE += nfe;
    }

    if (sumNFE != nEntities(coupleMap::NFE_SIZE))
    {
        FatalErrorIn("inline void coupleMap::makeFaces() const")
            << " Mismatched buffer." << nl
            << " sumNFE: " << sumNFE << nl
            << " NFE_SIZE: " << nEntities(coupleMap::NFE_SIZE) << nl
            << abort(FatalError);
    }
}


inline void coupleMap::makeCells() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (cellsPtr_)
    {
        FatalErrorIn("inline void coupleMap::makeCells() const")
            << "Cells have already been calculated."
            << abort(FatalError);
    }

    label nCells = nEntities(coupleMap::CELL);
    const labelList& cBuffer = entityBuffer(coupleMap::CELL);

    // Set number of cell faces
    label ncf = twoDMesh_ ? 5 : 4;

    cellsPtr_ = new cellList(nCells, cell(ncf));

    cellList& cells = *cellsPtr_;

    forAll(cells, cellI)
    {
        for (label f = 0; f < ncf; f++)
        {
            cells[cellI][f] = cBuffer[(ncf*cellI) + f];
        }
    }
}


inline void coupleMap::makeFaceMap() const
{
    // It is an error to attempt to recalculate
    // if the map is already calculated
    if (faceMap_.size())
    {
        FatalErrorIn("inline void coupleMap::makeFaceMap() const")
            << "faceMap has already been calculated."
            << abort(FatalError);
    }

    const Map<label>& mapFaces = entityMap(coupleMap::FACE);

    // Size the list
    faceMap_.setSize(mapFaces.size(), -1);

    // Fill-in entries
    forAllConstIter(Map<label>, mapFaces, fIter)
    {
        faceMap_[fIter.key()] = fIter();
    }
}


inline void coupleMap::makeCellMap() const
{
    // It is an error to attempt to recalculate
    // if the map is already calculated
    if (cellMap_.size())
    {
        FatalErrorIn("inline void coupleMap::makeCellMap() const")
            << "cellMap has already been calculated."
            << abort(FatalError);
    }

    const Map<label>& mapCells = entityMap(coupleMap::CELL);

    // Size the list
    cellMap_.setSize(mapCells.size(), -1);

    // Fill-in entries
    forAllConstIter(Map<label>, mapCells, cIter)
    {
        cellMap_[cIter.key()] = cIter();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline pointField& coupleMap::pointBuffer() const
{
    return pointBuffer_;
}


inline pointField& coupleMap::oldPointBuffer() const
{
    return oldPointBuffer_;
}


inline labelList& coupleMap::subMeshPoints() const
{
    return subMeshPoints_;
}


inline List<labelPair>& coupleMap::globalProcPoints() const
{
    return globalProcPoints_;
}


inline void coupleMap::allocateBuffers() const
{
    forAll(nEntities_, entityI)
    {
        if (nEntities_[entityI] < 0)
        {
            FatalErrorIn("inline void coupleMap::allocateBuffers() const")
                << " Entity sizes are not valid." << nl
                << " nEntities: " << nEntities_
                << abort(FatalError);
        }
    }

    // Size up point buffers
    pointBuffer().setSize(nEntities(coupleMap::POINT));
    oldPointBuffer().setSize(nEntities(coupleMap::POINT));

    // Size up connectivity buffers
    entityBuffer(coupleMap::POINT).setSize
    (
        nEntities(coupleMap::GLOBAL_POINT)
      - nEntities(coupleMap::SHARED_POINT)
    );

    entityBuffer(coupleMap::EDGE).setSize(2 * nEntities(coupleMap::EDGE));

    // Size up boundary buffers
    entityBuffer(coupleMap::FACE_STARTS).setSize(nEntities(coupleMap::NBDY));
    entityBuffer(coupleMap::FACE_SIZES).setSize(nEntities(coupleMap::NBDY));
    entityBuffer(coupleMap::EDGE_STARTS).setSize(nEntities(coupleMap::NBDY));
    entityBuffer(coupleMap::EDGE_SIZES).setSize(nEntities(coupleMap::NBDY));
    entityBuffer(coupleMap::PATCH_ID).setSize(nEntities(coupleMap::NBDY));

    // nFaceEdges buffer is required only for 2D,
    // due to a mix of triangle / quad faces
    if (twoDMesh_)
    {
        entityBuffer(coupleMap::NFE_BUFFER).setSize(nEntities(coupleMap::FACE));
    }

    // Set number of cell faces
    label ncf = twoDMesh_ ? 5 : 4;
    entityBuffer(coupleMap::CELL).setSize(ncf*nEntities(coupleMap::CELL));

    // Allocate for variable size face-lists
    // For 3D, this is always 3 times the number of faces
    entityBuffer(coupleMap::FACE).setSize(nEntities(coupleMap::NFE_SIZE));
    entityBuffer(coupleMap::FACE_EDGE).setSize(nEntities(coupleMap::NFE_SIZE));
}


inline label coupleMap::findSlave
(
    const label eType,
    const label Index
) const
{
    Map<label>::const_iterator it = entityMap_[eType].find(Index);

    if (it == entityMap_[eType].end())
    {
        return -1;
    }
    else
    {
        return it();
    }
}


inline label coupleMap::findMaster
(
    const label eType,
    const label Index
) const
{
    Map<label>::const_iterator it = reverseEntityMap_[eType].find(Index);

    if (it == reverseEntityMap_[eType].end())
    {
        return -1;
    }
    else
    {
        return it();
    }
}


inline void coupleMap::removeSlave
(
    const label eType,
    const label Index
) const
{
    Map<label>::iterator it = reverseEntityMap_[eType].find(Index);

    if (it != reverseEntityMap_[eType].end())
    {
        reverseEntityMap_[eType].erase(it);
    }
}


inline void coupleMap::removeMaster
(
    const label eType,
    const label Index
) const
{
    Map<label>::iterator it = entityMap_[eType].find(Index);

    if (it != entityMap_[eType].end())
    {
        entityMap_[eType].erase(it);
    }
}


inline void coupleMap::mapSlave
(
    const label eType,
    const label master,
    const label slave
) const
{
    entityMap_[eType].set(master, slave);
}


inline void coupleMap::mapMaster
(
    const label eType,
    const label slave,
    const label master
) const
{
    reverseEntityMap_[eType].set(slave, master);
}


inline void coupleMap::pushOperation
(
    const label index,
    const label opType,
    const point& newPoint,
    const point& oldPoint
) const
{
    entityIndices_.setSize(entityIndices_.size() + 1, index);
    entityOperations_.setSize(entityOperations_.size() + 1, opType);

    if
    (
        opType == coupleMap::MOVE_POINT ||
        opType == coupleMap::CONVERT_PATCH
    )
    {
        moveNewPoints_.setSize(moveNewPoints_.size() + 1, newPoint);
        moveOldPoints_.setSize(moveOldPoints_.size() + 1, oldPoint);
    }
}


inline void coupleMap::transferMaps
(
    const label eType,
    Map<label>& newEntityMap,
    Map<label>& newReverseEntityMap
) const
{
    entityMap_[eType].transfer(newEntityMap);
    reverseEntityMap_[eType].transfer(newReverseEntityMap);
}


inline void coupleMap::clearMaps() const
{
    faceMap_.clear();
    cellMap_.clear();

    forAll(entityMap_, mapI)
    {
        entityMap_[mapI].clear();
        reverseEntityMap_[mapI].clear();
    }
}


inline void coupleMap::clearBuffers() const
{
    pointBuffer_.clear();
    oldPointBuffer_.clear();

    subMeshPoints_.clear();
    globalProcPoints_.clear();

    forAll(entityBuffer_, bufferI)
    {
        entityBuffer_[bufferI].clear();
    }

    entityIndices_.clear();
    entityOperations_.clear();

    moveNewPoints_.clear();
    moveOldPoints_.clear();
}


inline label coupleMap::nInternalFaces() const
{
    if (nInternalFaces_ == -1)
    {
        makeAddressing();
    }

    return nInternalFaces_;
}


inline const labelList& coupleMap::owner() const
{
    if (!ownerPtr_)
    {
        makeAddressing();
    }

    return *ownerPtr_;
}


inline const labelList& coupleMap::neighbour() const
{
    if (!neighbourPtr_)
    {
        makeAddressing();
    }

    return *neighbourPtr_;
}


inline const edgeList& coupleMap::edges() const
{
    if (!edgesPtr_)
    {
        makeEdges();
    }

    return *edgesPtr_;
}


inline const faceList& coupleMap::faces() const
{
    if (!facesPtr_)
    {
        makeFaces();
    }

    return *facesPtr_;
}


inline const cellList& coupleMap::cells() const
{
    if (!cellsPtr_)
    {
        makeCells();
    }

    return *cellsPtr_;
}


inline const labelListList& coupleMap::faceEdges() const
{
    if (!faceEdgesPtr_)
    {
        makeFaces();
    }

    return *faceEdgesPtr_;
}


inline const labelList& coupleMap::faceMap() const
{
    if (faceMap_.empty())
    {
        makeFaceMap();
    }

    return faceMap_;
}


inline const labelList& coupleMap::cellMap() const
{
    if (cellMap_.empty())
    {
        makeCellMap();
    }

    return cellMap_;
}


inline bool coupleMap::readData(Istream& is)
{
    Map<label> tmpMap(is);

    entityMap(coupleMap::POINT).transfer(tmpMap);

    // Prepare the reversePointMap as well.
    const Map<label>& pMap = entityMap(coupleMap::POINT);
    Map<label>& rpMap = reverseEntityMap(coupleMap::POINT);

    forAllConstIter(Map<label>, pMap, pIter)
    {
        rpMap.set(pIter(), pIter.key());
    }

    return !is.bad();
}


inline bool coupleMap::writeData(Ostream& os) const
{
    // Only write-out point-map information
    // to avoid geometric checking.
    // The rest can be constructed topologically.
    return (os << entityMap(coupleMap::POINT)).good();;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void coupleMap::operator=(const coupleMap& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("coupleMap::operator=(const coupleMap&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}


} // End namespace Foam

// ************************************************************************* //
