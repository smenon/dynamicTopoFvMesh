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

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coupleMap, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupleMap::coupleMap
(
    const IOobject& io,
    const label masterIndex,
    const label slaveIndex
)
:
    regIOobject(io),
    masterIndex_(masterIndex),
    slaveIndex_(slaveIndex),
    nEntities_(-1)
{}

// Construct as copy
coupleMap::coupleMap(const coupleMap& cm)
:
    regIOobject(cm, true),
    masterIndex_(cm.masterIndex_),
    slaveIndex_(cm.slaveIndex_),
    nEntities_(cm.nEntities_)
{}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

coupleMap::~coupleMap()
{
    clearMaps();
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

pointField& coupleMap::pointBuffer() const
{
    return pointBuffer_;
}

const Map<label>& coupleMap::masterToSlaveMap() const
{
    return masterToSlave_;
}

const Map<label>& coupleMap::slaveToMasterMap() const
{
    return slaveToMaster_;
}

label coupleMap::findSlaveIndex(const label Index) const
{
    if (masterToSlave_.found(Index))
    {
        return masterToSlave_[Index];
    }
    else
    {
        return -1;
    }
}

label coupleMap::findMasterIndex(const label Index) const
{
    if (slaveToMaster_.found(Index))
    {
        return slaveToMaster_[Index];
    }
    else
    {
        return -1;
    }
}

void coupleMap::removeSlaveIndex(const label Index) const
{
    if (slaveToMaster_.found(Index))
    {
        slaveToMaster_.erase(Index);
    }
}

void coupleMap::removeMasterIndex(const label Index) const
{
    if (masterToSlave_.found(Index))
    {
        masterToSlave_.erase(Index);
    }
}

void coupleMap::mapSlave
(
    const label master,
    const label slave
) const
{
    masterToSlave_.insert(master, slave);
}

void coupleMap::mapMaster
(
    const label slave,
    const label master
) const
{
    slaveToMaster_.insert(slave, master);
}

void coupleMap::transferMaps
(
    Map<label>& newMasterToSlave,
    Map<label>& newSlaveToMaster
) const
{
    masterToSlave_.transfer(newMasterToSlave);
    slaveToMaster_.transfer(newSlaveToMaster);
}

void coupleMap::clearMaps() const
{
    masterToSlave_.clear();
    slaveToMaster_.clear();
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
