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

#include "dynamicTopoFvMesh.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Constructor for coupledPatchInfo
dynamicTopoFvMesh::coupledPatchInfo::coupledPatchInfo
(
    const dynamicTopoFvMesh& mesh,
    const label masterIndex,
    const label slaveIndex,
    const label mfzIndex,
    const label sfzIndex
)
:
    mesh_(mesh),
    builtMaps_(false),
    masterIndex_(masterIndex),
    slaveIndex_(slaveIndex),
    masterFaceZone_(mfzIndex),
    slaveFaceZone_(sfzIndex),
    nEntities_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dynamicTopoFvMesh::coupledPatchInfo::setMesh
(
    label index,
    dynamicTopoFvMesh* mesh
)
{
    subMesh_.set(mesh);

    // Modify the subMesh directory as well.
    subMesh_().meshSubDir = "proc_" + Foam::name(index);
}

bool dynamicTopoFvMesh::coupledPatchInfo::builtMaps() const
{
    return builtMaps_;
}

void dynamicTopoFvMesh::coupledPatchInfo::setBuiltMaps()
{
    builtMaps_ = true;
}

label dynamicTopoFvMesh::coupledPatchInfo::slaveIndex() const
{
    return slaveIndex_;
}

label dynamicTopoFvMesh::coupledPatchInfo::masterFaceZone() const
{
    return masterFaceZone_;
}

label dynamicTopoFvMesh::coupledPatchInfo::slaveFaceZone() const
{
    return slaveFaceZone_;
}

FixedList<label,6>&
dynamicTopoFvMesh::coupledPatchInfo::nEntities()
{
    return nEntities_;
}

label& dynamicTopoFvMesh::coupledPatchInfo::nEntities
(
    const label eType
)
{
    return nEntities_[eType];
}

Map<label>&
dynamicTopoFvMesh::coupledPatchInfo::entityMap
(
    const label eType
)
{
    return entityMap_[eType];
}

Map<label>&
dynamicTopoFvMesh::coupledPatchInfo::reverseEntityMap
(
    const label eType
)
{
    return reverseEntityMap_[eType];
}

const Map<label>&
dynamicTopoFvMesh::coupledPatchInfo::masterToSlaveMap()
{
    return masterToSlave_;
}

const Map<label>&
dynamicTopoFvMesh::coupledPatchInfo::slaveToMasterMap()
{
    return slaveToMaster_;
}

label dynamicTopoFvMesh::coupledPatchInfo::findSlaveIndex
(
    const label Index
) const
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

label dynamicTopoFvMesh::coupledPatchInfo::findMasterIndex
(
    const label Index
) const
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

// Remove a master/slave entry, if found.
void dynamicTopoFvMesh::coupledPatchInfo::removeSlaveIndex
(
    const label Index
)
{
    if (slaveToMaster_.found(Index))
    {
        slaveToMaster_.erase(Index);
    }
}

void dynamicTopoFvMesh::coupledPatchInfo::removeMasterIndex
(
    const label Index
)
{
    if (masterToSlave_.found(Index))
    {
        masterToSlave_.erase(Index);
    }
}

void dynamicTopoFvMesh::coupledPatchInfo::mapSlave
(
    const label master,
    const label slave
)
{
    masterToSlave_.insert(master, slave);
}

void dynamicTopoFvMesh::coupledPatchInfo::mapMaster
(
    const label slave,
    const label master
)
{
    slaveToMaster_.insert(slave, master);
}

void dynamicTopoFvMesh::coupledPatchInfo::transferMaps
(
    Map<label>& newMasterToSlave,
    Map<label>& newSlaveToMaster
)
{
    masterToSlave_.transfer(newMasterToSlave);
    slaveToMaster_.transfer(newSlaveToMaster);
}

void dynamicTopoFvMesh::coupledPatchInfo::clearMaps()
{
    masterToSlave_.clear();
    slaveToMaster_.clear();
}

FixedList<labelList,5>&
dynamicTopoFvMesh::coupledPatchInfo::entityBuffer()
{
    return entityBuffer_;
}

labelList& dynamicTopoFvMesh::coupledPatchInfo::entityBuffer
(
    const label eType
)
{
    return entityBuffer_[eType];
}

pointField& dynamicTopoFvMesh::coupledPatchInfo::pointBuffer()
{
    return pointBuffer_;
}

scalarList& dynamicTopoFvMesh::coupledPatchInfo::lengthBuffer()
{
    return lengthBuffer_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void dynamicTopoFvMesh::coupledPatchInfo::operator=
(
    const coupledPatchInfo& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "coupledPatchInfo::operator=(const Foam::coupledPatchInfo&)"
        )
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
