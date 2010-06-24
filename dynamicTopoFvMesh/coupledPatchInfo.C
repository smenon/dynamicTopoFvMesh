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
    const coupleMap& cMap,
    const label mfzIndex,
    const label sfzIndex
)
:
    mesh_(mesh),
    builtMaps_(false),
    patchMap_(cMap),
    masterFaceZone_(mfzIndex),
    slaveFaceZone_(sfzIndex)
{}

dynamicTopoFvMesh::coupledPatchInfo::coupledPatchInfo
(
    const dynamicTopoFvMesh& mesh,
    const bool isLocal,
    const bool isSend,
    const label mPatch,
    const label sPatch,
    const label mfzIndex,
    const label sfzIndex
)
:
    mesh_(mesh),
    builtMaps_(false),
    patchMap_
    (
        IOobject
        (
            "coupleMap_"
          + Foam::name(mPatch)
          + "_To_"
          + Foam::name(sPatch)
          + word(isLocal ? "_Local" : "_Proc")
          + word(isSend ? "_Send" : "_Recv"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        isLocal,
        isSend,
        mPatch,
        sPatch
    ),
    masterFaceZone_(mfzIndex),
    slaveFaceZone_(sfzIndex)
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

dynamicTopoFvMesh& dynamicTopoFvMesh::coupledPatchInfo::subMesh()
{
    if (!subMesh_.valid())
    {
        FatalErrorIn
        (
            "dynamicTopoFvMesh& "
            "dynamicTopoFvMesh::coupledPatchInfo::subMesh()"
        )
            << " Sub-mesh pointer has not been set."
            << abort(FatalError);
    }

    return subMesh_();
}

bool dynamicTopoFvMesh::coupledPatchInfo::builtMaps() const
{
    return builtMaps_;
}

void dynamicTopoFvMesh::coupledPatchInfo::setBuiltMaps()
{
    builtMaps_ = true;
}

const coupleMap& dynamicTopoFvMesh::coupledPatchInfo::patchMap() const
{
    return patchMap_;
}

label dynamicTopoFvMesh::coupledPatchInfo::masterFaceZone() const
{
    return masterFaceZone_;
}

label dynamicTopoFvMesh::coupledPatchInfo::slaveFaceZone() const
{
    return slaveFaceZone_;
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
            "void dynamicTopoFvMesh::"
            "coupledPatchInfo::operator="
            "(const Foam::coupledPatchInfo&)"
        )
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
