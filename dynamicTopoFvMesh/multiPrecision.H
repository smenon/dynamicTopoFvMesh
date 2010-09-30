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

SourceFiles
    multiPrecisionI.H
    multiPrecision.C

\*---------------------------------------------------------------------------*/

#ifndef multiPrecision_H
#define multiPrecision_H

#include <mpfr.h>
#include "word.H"
#include "fieldTypes.H"

namespace Foam
{

template<unsigned P, mpfr_rnd_t R>
class multiPrecisionScalar;

template<unsigned P, mpfr_rnd_t R>
Ostream& operator<<(Ostream&, const multiPrecisionScalar<P,R>&);

template<unsigned P, mpfr_rnd_t R>
Istream& operator>>(Istream&, multiPrecisionScalar<P,R>&);

/*---------------------------------------------------------------------------*\
                      Class multiPrecisionScalar Declaration
\*---------------------------------------------------------------------------*/

template<unsigned P, mpfr_rnd_t R>
class multiPrecisionScalar
{
    //- Multi-precision float type
    mpfr_t v_;

#   ifdef FULLDEBUG
    //- Convert to word, useful for debugging
    word dv_;
#   endif

public:

    //- Constructors

        // Construct null
        inline multiPrecisionScalar();

        // Construct as copy
        inline multiPrecisionScalar(const multiPrecisionScalar<P,R>& v);

        // Construct from scalar
        inline multiPrecisionScalar(const scalar& v);

        // Construct from string constant
        inline multiPrecisionScalar(const char* s, const int base=10);

        // Construct from input stream
        inline multiPrecisionScalar(Istream& is);

    //- Destructor

        inline ~multiPrecisionScalar();

    //- Member functions

        inline const multiPrecisionScalar<P,R> mag() const;
        inline const multiPrecisionScalar<P,R> sqr() const;
        inline const multiPrecisionScalar<P,R> sqrt() const;

    //- Operator overloads

        // Assignment
        inline multiPrecisionScalar& operator=
        (
            const scalar& v
        );

        inline multiPrecisionScalar& operator=
        (
            const multiPrecisionScalar<P,R>& v
        );

        // Arithmetic
        inline multiPrecisionScalar& operator+=
        (
            const multiPrecisionScalar<P,R>& v
        );

        inline multiPrecisionScalar& operator-=
        (
            const multiPrecisionScalar<P,R>& v
        );

        inline multiPrecisionScalar& operator*=
        (
            const multiPrecisionScalar<P,R>& v
        );

        inline multiPrecisionScalar& operator/=
        (
            const multiPrecisionScalar<P,R>& v
        );

        // Comparison
        inline bool operator>(const multiPrecisionScalar<P,R>& v);
        inline bool operator<(const multiPrecisionScalar<P,R>& v);
        inline bool operator>(const scalar& v);
        inline bool operator<(const scalar& v);

        // Convert to string
        inline const word convertToWord
        (
            const size_t n,
            const int base = 10
        ) const;

        // Type-casts
        inline operator scalar() const;

        // IO
        friend Ostream& operator<<
        #ifndef __CINT__
        <P,R>
        #endif
        (Ostream&, const multiPrecisionScalar<P,R>&);

        friend Istream& operator>>
        #ifndef __CINT__
        <P,R>
        #endif
        (Istream&, multiPrecisionScalar<P,R>&);
};


// Define traits for multiPrecisionScalar
template<unsigned P, mpfr_rnd_t R>
class pTraits<multiPrecisionScalar<P,R> >
{
    multiPrecisionScalar<P,R> p_;

public:

    //- Component type
    typedef multiPrecisionScalar<P,R> cmptType;

    // Member constants

        enum
        {
            dim = 3,         // Dimensionality of space
            rank = 0,        // Rank of multiPrecisionScalar is 0
            nComponents = 1  // Number of components is 1
        };

    // Static data members

        static const char* const typeName;
        static const char* componentNames[];
        static const multiPrecisionScalar<P,R> zero;
        static const multiPrecisionScalar<P,R> one;

    // Constructors

        //- Construct from Istream
        pTraits(Istream& is);

    // Member Functions

        operator multiPrecisionScalar<P,R>() const
        {
            return p_;
        }
};


// Typedefs
typedef multiPrecisionScalar<64, MPFR_RNDN> mpScalar;

typedef Vector<mpScalar> mpVector;
typedef Tensor<mpScalar> mpTensor;
typedef SymmTensor<mpScalar> mpSymmTensor;
typedef SphericalTensor<mpScalar> mpSphericalTensor;

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiPrecisionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
