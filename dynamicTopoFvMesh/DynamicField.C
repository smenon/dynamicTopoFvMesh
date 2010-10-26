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

#include "DynamicField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
const char* const DynamicField<T, SizeInc, SizeMult, SizeDiv>::typeName
(
    "DynamicField"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
DynamicField<T, SizeInc, SizeMult, SizeDiv>::DynamicField(Istream& is)
:
    Field<T>(is),
    capacity_(Field<T>::size())
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
tmp<DynamicField<T, SizeInc, SizeMult, SizeDiv> >
DynamicField<T, SizeInc, SizeMult, SizeDiv>::clone() const
{
    return tmp<DynamicField<T, SizeInc, SizeMult, SizeDiv> >
    (
        new DynamicField<T, SizeInc, SizeMult, SizeDiv>(*this)
    );
}


// * * * * * * * * * * * * * * * IOstream Operator * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Ostream& operator<<
(
    Ostream& os,
    const DynamicField<T, SizeInc, SizeMult, SizeDiv>& f
)
{
    os << static_cast<const Field<T>&>(f);
    return os;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Ostream& operator<<
(
    Ostream& os,
    const tmp<DynamicField<T, SizeInc, SizeMult, SizeDiv> >& tf
)
{
    os << tf();
    tf.clear();
    return os;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Istream& operator>>
(
    Istream& is,
    DynamicField<T, SizeInc, SizeMult, SizeDiv>& lst
)
{
    is >> static_cast<Field<T>&>(lst);
    lst.capacity_ = lst.Field<T>::size();

    return is;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
