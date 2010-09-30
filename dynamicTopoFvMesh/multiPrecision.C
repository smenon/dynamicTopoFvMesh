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
    multiPrecision

Description
    Various utility functions that perform multi-precision operations,
    using the Multi-Precision Floating-Point Reliable (MPFR) library.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "multiPrecision.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const pTraits<mpScalar>::typeName = "mpScalar";

template<>
const mpScalar pTraits<mpScalar>::zero = mpScalar("0.0");

template<>
const mpScalar pTraits<mpScalar>::one = mpScalar("1.0");

template<>
const char* pTraits<mpScalar>::componentNames[] = { "x" };

template<>
pTraits<mpScalar>::pTraits(Istream& is)
{
    is >> p_;
}

template<>
const char* const mpVector::typeName = "mpVector";

template<>
const char* mpVector::componentNames[] = {"x", "y", "z"};

template<>
const mpVector mpVector::zero =
(
   mpVector
   (
       pTraits<mpScalar>::zero,
       pTraits<mpScalar>::zero,
       pTraits<mpScalar>::zero
   )
);

template<>
const mpVector mpVector::one =
(
   mpVector
   (
       pTraits<mpScalar>::one,
       pTraits<mpScalar>::one,
       pTraits<mpScalar>::one
   )
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
