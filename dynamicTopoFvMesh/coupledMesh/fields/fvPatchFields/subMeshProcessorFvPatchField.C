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

#include "subMeshProcessorFvPatchField.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

// Construct from patch and internal field
template<class Type>
Foam::subMeshProcessorFvPatchField<Type>::subMeshProcessorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    procPatch_(refCast<const subMeshProcessorFvPatch>(p))
{}


// Construct from patch and internal field and patch field
template<class Type>
Foam::subMeshProcessorFvPatchField<Type>::subMeshProcessorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    coupledFvPatchField<Type>(p, iF, f),
    procPatch_(refCast<const subMeshProcessorFvPatch>(p))
{}


// Construct from patch, internal field and dictionary
template<class Type>
Foam::subMeshProcessorFvPatchField<Type>::subMeshProcessorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const subMeshProcessorFvPatch>(p))
{}


// Construct by mapping given fvPatchField onto a new patch
template<class Type>
Foam::subMeshProcessorFvPatchField<Type>::subMeshProcessorFvPatchField
(
    const subMeshProcessorFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const subMeshProcessorFvPatch>(p))
{}


// Construct as copy
template<class Type>
Foam::subMeshProcessorFvPatchField<Type>::subMeshProcessorFvPatchField
(
    const subMeshProcessorFvPatchField<Type>& ptf
)
:
    coupledFvPatchField<Type>(ptf),
    procPatch_(refCast<const subMeshProcessorFvPatch>(ptf.patch()))
{}


// Construct as copy setting internal field reference
template<class Type>
Foam::subMeshProcessorFvPatchField<Type>::subMeshProcessorFvPatchField
(
    const subMeshProcessorFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    procPatch_(refCast<const subMeshProcessorFvPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
Foam::subMeshProcessorFvPatchField<Type>::~subMeshProcessorFvPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return neighbour field given internal field
template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::subMeshProcessorFvPatchField<Type>::patchNeighbourField() const
{
    return *this;
}


// Initialise the evaluation of the patch field
template<class Type>
void Foam::subMeshProcessorFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{}


// Evaluate the patch field
template<class Type>
void Foam::subMeshProcessorFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{}


// Return patch-normal gradient
template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::subMeshProcessorFvPatchField<Type>::snGrad
(
    const scalarField& deltaCoeffs
) const
{
    return deltaCoeffs * (*this - this->patchInternalField());
}


// Is all data available
template<class Type>
bool Foam::subMeshProcessorFvPatchField<Type>::ready() const
{
    return true;
}


// Initialise neighbour matrix update
template<class Type>
void Foam::subMeshProcessorFvPatchField<Type>::initInterfaceMatrixUpdate
(
    scalarField&,
    const scalarField& psiInternal,
    const scalarField&,
    const direction,
    const Pstream::commsTypes commsType
) const
{}


// Update result field based on interface functionality
template<class Type>
void Foam::subMeshProcessorFvPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
    const scalarField&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{}


// Initialise neighbour matrix update
template<class Type>
void Foam::subMeshProcessorFvPatchField<Type>::initInterfaceMatrixUpdate
(
    Field<Type>&,
    const Field<Type>& psiInternal,
    const scalarField&,
    const Pstream::commsTypes commsType
) const
{}


// Update result field based on interface functionality
template<class Type>
void Foam::subMeshProcessorFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const Field<Type>&,
    const scalarField& coeffs,
    const Pstream::commsTypes commsType
) const
{}


// ************************************************************************* //
