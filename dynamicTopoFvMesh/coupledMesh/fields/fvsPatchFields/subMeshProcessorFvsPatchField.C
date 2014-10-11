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

#include "subMeshProcessorFvsPatchField.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

// Construct from patch and internal field
template<class Type>
Foam::subMeshProcessorFvsPatchField<Type>::subMeshProcessorFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(p, iF),
    procPatch_(refCast<const subMeshProcessorFvPatch>(p))
{}


// Construct from patch and internal field and patch field
template<class Type>
Foam::subMeshProcessorFvsPatchField<Type>::subMeshProcessorFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const Field<Type>& f
)
:
    coupledFvsPatchField<Type>(p, iF, f),
    procPatch_(refCast<const subMeshProcessorFvPatch>(p))
{}


// Construct from patch, internal field and dictionary
template<class Type>
Foam::subMeshProcessorFvsPatchField<Type>::subMeshProcessorFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    coupledFvsPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const subMeshProcessorFvPatch>(p))
{}


// Construct by mapping given fvsPatchField onto a new patch
template<class Type>
Foam::subMeshProcessorFvsPatchField<Type>::subMeshProcessorFvsPatchField
(
    const subMeshProcessorFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvsPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const subMeshProcessorFvPatch>(p))
{}


// Construct as copy
template<class Type>
Foam::subMeshProcessorFvsPatchField<Type>::subMeshProcessorFvsPatchField
(
    const subMeshProcessorFvsPatchField<Type>& ptf
)
:
    coupledFvsPatchField<Type>(ptf),
    procPatch_(refCast<const subMeshProcessorFvPatch>(ptf.patch()))
{}


// Construct as copy setting internal field reference
template<class Type>
Foam::subMeshProcessorFvsPatchField<Type>::subMeshProcessorFvsPatchField
(
    const subMeshProcessorFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(ptf, iF),
    procPatch_(refCast<const subMeshProcessorFvPatch>(ptf.patch()))
{}

// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
Foam::subMeshProcessorFvsPatchField<Type>::~subMeshProcessorFvsPatchField()
{}


// ************************************************************************* //
