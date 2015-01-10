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

#include "subMeshProcessorPointPatchField.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

// Construct from patch and internal field
template<class Type>
Foam::subMeshProcessorPointPatchField<Type>::subMeshProcessorPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(p, iF),
    procPatch_(refCast<const subMeshProcessorPointPatch>(p))
{}


// Construct from patch and internal field and patch field
template<class Type>
Foam::subMeshProcessorPointPatchField<Type>::subMeshProcessorPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const Field<Type>& f
)
:
    coupledPointPatchField<Type>(p, iF, f),
    procPatch_(refCast<const subMeshProcessorPointPatch>(p))
{}


// Construct from patch, internal field and dictionary
template<class Type>
Foam::subMeshProcessorPointPatchField<Type>::subMeshProcessorPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    coupledPointPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const subMeshProcessorPointPatch>(p))
{}


// Construct by mapping given PointPatchField onto a new patch
template<class Type>
Foam::subMeshProcessorPointPatchField<Type>::subMeshProcessorPointPatchField
(
    const subMeshProcessorPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    coupledPointPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const subMeshProcessorPointPatch>(p))
{}


// Construct as copy
template<class Type>
Foam::subMeshProcessorPointPatchField<Type>::subMeshProcessorPointPatchField
(
    const subMeshProcessorPointPatchField<Type>& ptf
)
:
    coupledPointPatchField<Type>(ptf),
    procPatch_(refCast<const subMeshProcessorPointPatch>(ptf.patch()))
{}


// Construct as copy setting internal field reference
template<class Type>
Foam::subMeshProcessorPointPatchField<Type>::subMeshProcessorPointPatchField
(
    const subMeshProcessorPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(ptf, iF),
    procPatch_(refCast<const subMeshProcessorPointPatch>(ptf.patch()))
{}

// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
Foam::subMeshProcessorPointPatchField<Type>::~subMeshProcessorPointPatchField()
{}


// ************************************************************************* //
