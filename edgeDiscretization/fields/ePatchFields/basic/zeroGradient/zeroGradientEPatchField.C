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

#include "zeroGradientEPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
zeroGradientEPatchField<Type>::zeroGradientEPatchField
(
    const ePatch& p,
    const DimensionedField<Type, EdgeMesh>& iF
)
:
    ePatchField<Type>(p, iF)
{}

/*
template<class Type>
zeroGradientEPatchField<Type>::zeroGradientEPatchField
(
    const zeroGradientEPatchField<Type>& ptf,
    const ePatch& p,
    const DimensionedField<Type, EdgeMesh>& iF,
    const ePatchFieldMapper& mapper
)
:
    ePatchField<Type>(ptf, p, iF, mapper)
{}
*/


template<class Type>
zeroGradientEPatchField<Type>::zeroGradientEPatchField
(
    const ePatch& p,
    const DimensionedField<Type, EdgeMesh>& iF,
    const dictionary&
)
:
    ePatchField<Type>(p, iF)
{
    // ePatchField<Type>::operator=(this->patchInternalField());
}


template<class Type>
zeroGradientEPatchField<Type>::zeroGradientEPatchField
(
    const zeroGradientEPatchField& zgpf
)
:
    ePatchField<Type>(zgpf)
{}


template<class Type>
zeroGradientEPatchField<Type>::zeroGradientEPatchField
(
    const zeroGradientEPatchField& zgpf,
    const DimensionedField<Type, EdgeMesh>& iF
)
:
    ePatchField<Type>(zgpf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void zeroGradientEPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // this->operator==(this->patchInternalField());
    ePatchField<Type>::evaluate();
}


template<class Type>
tmp<Field<Type> > zeroGradientEPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::one)
    );
}


template<class Type>
tmp<Field<Type> > zeroGradientEPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}


template<class Type>
tmp<Field<Type> > zeroGradientEPatchField<Type>::gradientInternalCoeffs() const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}


template<class Type>
tmp<Field<Type> > zeroGradientEPatchField<Type>::gradientBoundaryCoeffs() const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
