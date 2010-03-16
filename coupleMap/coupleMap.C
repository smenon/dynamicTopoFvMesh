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

\*----------------------------------------------------------------------------*/

#include "coupleMap.H"
#include "boolList.H"
#include "demandDrivenData.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coupleMap, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupleMap::coupleMap
(
    const IOobject& io,
    const bool isLocal,
    const label masterIndex,
    const label slaveIndex
)
:
    regIOobject(io),
    isLocal_(isLocal),
    masterIndex_(masterIndex),
    slaveIndex_(slaveIndex),
    nEntities_(-1),
    nInternalFaces_(-1),
    ownerPtr_(NULL),
    neighbourPtr_(NULL)
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
coupleMap::coupleMap(const coupleMap& cm)
:
    regIOobject(cm, true),
    isLocal_(cm.isLocal_),
    masterIndex_(cm.masterIndex_),
    slaveIndex_(cm.slaveIndex_),
    nEntities_(cm.nEntities_),
    nInternalFaces_(-1),
    ownerPtr_(NULL),
    neighbourPtr_(NULL)
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

coupleMap::~coupleMap()
{
    clearMaps();

    nInternalFaces_ = -1;
    deleteDemandDrivenData(ownerPtr_);
    deleteDemandDrivenData(neighbourPtr_);
}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void coupleMap::makeAddressing() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ownerPtr_ || neighbourPtr_ || nInternalFaces_ > -1)
    {
        FatalErrorIn("coupleMap::makeAddressing()")
            << "Addressing has already been calculated."
            << abort(FatalError);
    }

    label nFaces = nEntities(coupleMap::FACE);
    label nCells = nEntities(coupleMap::CELL);

    const labelList& cBuffer = entityBuffer(coupleMap::CELL);

    if (nCells < 0 || nFaces < 0)
    {
        FatalErrorIn("coupleMap::makeAddressing()")
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

    for (label cellI = 0; cellI < nCells; cellI++)
    {
        for (label f = 0; f < 4; f++)
        {
            label faceI = cBuffer[(4*cellI)+f];

            if (!markedFaces[faceI])
            {
                // First visit: owner
                own[faceI] = cellI;
                markedFaces[faceI] = true;
            }
            else
            {
                // Second visit: neighbour
                nei[faceI] = cellI;
                nInternalFaces_++;
            }
        }
    }

    nei.setSize(nInternalFaces_);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label coupleMap::masterIndex() const
{
    return masterIndex_;
}

label coupleMap::slaveIndex() const
{
    return slaveIndex_;
}

bool coupleMap::isLocal() const
{
    return isLocal_;
}

pointField& coupleMap::pointBuffer() const
{
    return pointBuffer_;
}

labelList& coupleMap::subMeshPoints() const
{
    return subMeshPoints_;
}

label coupleMap::findSlaveIndex
(
    const label eType,
    const label Index
) const
{
    if (entityMap_[eType].found(Index))
    {
        return entityMap_[eType][Index];
    }
    else
    {
        return -1;
    }
}

label coupleMap::findMasterIndex
(
    const label eType,
    const label Index
) const
{
    if (reverseEntityMap_[eType].found(Index))
    {
        return reverseEntityMap_[eType][Index];
    }
    else
    {
        return -1;
    }
}

void coupleMap::removeSlaveIndex
(
    const label eType,
    const label Index
) const
{
    if (reverseEntityMap_[eType].found(Index))
    {
        reverseEntityMap_[eType].erase(Index);
    }
}

void coupleMap::removeMasterIndex
(
    const label eType,
    const label Index
) const
{
    if (entityMap_[eType].found(Index))
    {
        entityMap_[eType].erase(Index);
    }
}

void coupleMap::mapSlave
(
    const label eType,
    const label master,
    const label slave
) const
{
    entityMap_[eType].set(master, slave);
}

void coupleMap::mapMaster
(
    const label eType,
    const label slave,
    const label master
) const
{
    reverseEntityMap_[eType].set(slave, master);
}

void coupleMap::transferMaps
(
    const label eType,
    Map<label>& newEntityMap,
    Map<label>& newReverseEntityMap
) const
{
    entityMap_[eType].transfer(newEntityMap);
    reverseEntityMap_[eType].transfer(newReverseEntityMap);
}

void coupleMap::clearMaps() const
{
    forAll(entityMap_, mapI)
    {
        entityMap_[mapI].clear();
        reverseEntityMap_[mapI].clear();
    }
}

FixedList<label,6>& coupleMap::nEntities() const
{
    return nEntities_;
}

label& coupleMap::nEntities(const label eType) const
{
    return nEntities_[eType];
}

Map<label>& coupleMap::entityMap(const label eType) const
{
    return entityMap_[eType];
}

Map<label>& coupleMap::reverseEntityMap(const label eType) const
{
    return reverseEntityMap_[eType];
}

FixedList<labelList,5>& coupleMap::entityBuffer() const
{
    return entityBuffer_;
}

labelList& coupleMap::entityBuffer(const label eType) const
{
    return entityBuffer_[eType];
}

label coupleMap::nInternalFaces() const
{
    if (nInternalFaces_ == -1)
    {
        makeAddressing();
    }

    return nInternalFaces_;
}

const labelList& coupleMap::owner() const
{
    if (!ownerPtr_)
    {
        makeAddressing();
    }

    return *ownerPtr_;
}

const labelList& coupleMap::neighbour() const
{
    if (!neighbourPtr_)
    {
        makeAddressing();
    }

    return *neighbourPtr_;
}

bool coupleMap::readData(Istream& is)
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

bool coupleMap::writeData(Ostream& os) const
{
    // Only write-out point-map information
    // to avoid geometric checking.
    // The rest can be constructed topologically.
    return (os << entityMap(coupleMap::POINT)).good();;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void coupleMap::operator=(const coupleMap& rhs)
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
