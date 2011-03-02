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

#include "Time.H"
#include "coupledInfo.H"
#include "dynamicTopoFvMesh.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Constructor for coupledInfo
coupledInfo::coupledInfo
(
    const dynamicTopoFvMesh& mesh,
    const coupleMap& cMap,
    const label mfzIndex,
    const label sfzIndex
)
:
    mesh_(mesh),
    builtMaps_(false),
    map_(cMap),
    masterFaceZone_(mfzIndex),
    slaveFaceZone_(sfzIndex)
{}


coupledInfo::coupledInfo
(
    const dynamicTopoFvMesh& mesh,
    const bool isTwoDMesh,
    const bool isLocal,
    const bool isSend,
    const label patchIndex,
    const label mPatch,
    const label sPatch,
    const label mfzIndex,
    const label sfzIndex
)
:
    mesh_(mesh),
    builtMaps_(false),
    map_
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
        isTwoDMesh,
        isLocal,
        isSend,
        patchIndex,
        mPatch,
        sPatch
    ),
    masterFaceZone_(mfzIndex),
    slaveFaceZone_(sfzIndex)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void coupledInfo::setMesh
(
    label index,
    dynamicTopoFvMesh* mesh
)
{
    subMesh_.set(mesh);
}


dynamicTopoFvMesh& coupledInfo::subMesh()
{
    if (!subMesh_.valid())
    {
        FatalErrorIn("dynamicTopoFvMesh& coupledInfo::subMesh()")
            << " Sub-mesh pointer has not been set."
            << abort(FatalError);
    }

    return subMesh_();
}


const dynamicTopoFvMesh& coupledInfo::subMesh() const
{
    if (!subMesh_.valid())
    {
        FatalErrorIn("const dynamicTopoFvMesh& coupledInfo::subMesh() const")
            << " Sub-mesh pointer has not been set."
            << abort(FatalError);
    }

    return subMesh_();
}


bool coupledInfo::builtMaps() const
{
    return builtMaps_;
}


void coupledInfo::setBuiltMaps()
{
    builtMaps_ = true;
}


coupleMap& coupledInfo::map()
{
    return map_;
}


const coupleMap& coupledInfo::map() const
{
    return map_;
}


label coupledInfo::masterFaceZone() const
{
    return masterFaceZone_;
}


label coupledInfo::slaveFaceZone() const
{
    return slaveFaceZone_;
}


template <class Type>
void coupledInfo::mapVolField
(
    const wordList& fieldNames,
    const word& fieldType,
    OSstream& strStream
) const
{
    strStream
        << fieldType << token::NL
        << token::BEGIN_BLOCK << token::NL;

    forAll(fieldNames, i)
    {
        const GeometricField<Type, fvPatchField, volMesh>& fld =
        (
            mesh_.lookupObject
            <
                GeometricField<Type, fvPatchField, volMesh>
            >(fieldNames[i])
        );

        // Create and map the internal-field values
        Field<Type> internalField
        (
            fld.internalField(),
            map().cellMap()
        );

        // Create and map the patch field values
//        PtrList<fvPatchField<Type> > patchFields(pm.size());
    }

    strStream
        << token::END_BLOCK << token::NL;
}


template <class Type>
void coupledInfo::mapSurfaceField
(
    const wordList& fieldNames,
    const word& fieldType,
    OSstream& strStream
) const
{
    strStream
        << fieldType << token::NL
        << token::BEGIN_BLOCK << token::NL;

    forAll(fieldNames, i)
    {
        const GeometricField<Type, fvsPatchField, surfaceMesh>& fld =
        (
            mesh_.lookupObject
            <
                GeometricField<Type, fvsPatchField, surfaceMesh>
            >(fieldNames[i])
        );

        // Create and map the internal-field values
        Field<Type> internalField
        (
            fld.internalField(),
            SubList<label>
            (
                map().faceMap(),
                subMesh().nInternalFaces()
            )
        );
    }

    strStream
        << token::END_BLOCK << token::NL;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void coupledInfo::operator=(const coupledInfo& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "void coupledInfo::operator=(const Foam::coupledInfo&)"
        )
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
