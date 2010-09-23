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

template<unsigned P, mpfr_rnd_t R>
const char*
const pTraits<multiPrecisionScalar<P,R> >::typeName = "multiPrecisionScalar";

template<unsigned P, mpfr_rnd_t R>
const multiPrecisionScalar<P,R>
pTraits<multiPrecisionScalar<P,R> >::zero = multiPrecisionScalar<P,R>("0.0");

template<unsigned P, mpfr_rnd_t R>
const multiPrecisionScalar<P,R>
pTraits<multiPrecisionScalar<P,R> >::one = multiPrecisionScalar<P,R>("1.0");

template<unsigned P, mpfr_rnd_t R>
const char* pTraits<multiPrecisionScalar<P,R> >::componentNames[] = { "x" };

template<unsigned P, mpfr_rnd_t R>
pTraits<multiPrecisionScalar<P,R> >::pTraits(Istream& is)
{
    is >> p_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>::multiPrecisionScalar()
{
    // Initialize with precision
    mpfr_init2(v_, P);

    // Set to zero and rounding mode
    mpfr_set_d(v_, 0.0, R);
}


// Construct as copy
template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>::multiPrecisionScalar
(
    const multiPrecisionScalar<P,R>& v
)
{
    // Initialize with precision
    mpfr_init2(v_, P);

    // Set to value and rounding mode
    mpfr_set(v_, v.v_, R);
}


// Construct from scalar
template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>::multiPrecisionScalar
(
    const scalar& v
)
{
    // Initialize with precision
    mpfr_init2(v_, P);

    // Set to value and rounding mode
    mpfr_set_d(v_, v, R);
}


// Construct from string constant
template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>::multiPrecisionScalar
(
    const char* s,
    const int base
)
{
    // Initialize with precision
    mpfr_init2(v_, P);

    // Set to value and rounding mode
    mpfr_set_str(v_, s, base, R);
}


// Construct from input stream
template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>::multiPrecisionScalar(Istream& is)
:
    multiPrecisionScalar<P,R>()
{
    is >> *this;
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>::~multiPrecisionScalar()
{
    mpfr_clear(v_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> multiPrecisionScalar<P,R>::mag() const
{
    multiPrecisionScalar<P,R> result(*this);

    mpfr_abs(result.v_, result.v_, R);

    return result;
}


template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> multiPrecisionScalar<P,R>::sqr() const
{
    multiPrecisionScalar<P,R> result(*this);

    mpfr_sqr(result.v_, result.v_, R);

    return result;
}


template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> multiPrecisionScalar<P,R>::sqrt() const
{
    multiPrecisionScalar<P,R> result(*this);

    mpfr_sqrt(result.v_, result.v_, R);

    return result;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>& multiPrecisionScalar<P,R>::operator=
(
    const scalar& v
)
{
    mpfr_set_d(v_, v, R);

    return *this;
}


template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>& multiPrecisionScalar<P,R>::operator=
(
    const multiPrecisionScalar<P,R>& v
)
{
    if (this == &v)
    {
        FatalErrorIn
        (
            "multiPrecisionScalar::operator="
            "(const multiPrecisionScalar&)"
        )
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    mpfr_set(v_, v.v_, R);

    return *this;
}


template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>& multiPrecisionScalar<P,R>::operator+=
(
    const multiPrecisionScalar<P,R>& v
)
{
    mpfr_add(v_, v_, v.v_, R);

    return *this;
}


template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>& multiPrecisionScalar<P,R>::operator-=
(
    const multiPrecisionScalar<P,R>& v
)
{
    mpfr_sub(v_, v_, v.v_, R);

    return *this;
}


template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>& multiPrecisionScalar<P,R>::operator*=
(
    const multiPrecisionScalar<P,R>& v
)
{
    mpfr_mul(v_, v_, v.v_, R);

    return *this;
}


template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>& multiPrecisionScalar<P,R>::operator/=
(
    const multiPrecisionScalar<P,R>& v
)
{
    mpfr_div(v_, v_, v.v_, R);

    return *this;
}


template<unsigned P, mpfr_rnd_t R>
inline bool multiPrecisionScalar<P,R>::operator>
(
    const multiPrecisionScalar<P,R>& v
)
{
    return ( mpfr_greater_p(v_, v.v_) != 0 );
}


template<unsigned P, mpfr_rnd_t R>
inline bool multiPrecisionScalar<P,R>::operator<
(
    const multiPrecisionScalar<P,R>& v
)
{
    return ( mpfr_less_p(v_, v.v_) != 0 );
}


template<unsigned P, mpfr_rnd_t R>
inline bool multiPrecisionScalar<P,R>::operator>
(
    const scalar& v
)
{
    return ( (*this) > multiPrecisionScalar<P,R>(v) );
}


template<unsigned P, mpfr_rnd_t R>
inline bool multiPrecisionScalar<P,R>::operator<
(
    const scalar& v
)
{
    return ( (*this) < multiPrecisionScalar<P,R>(v) );
}


template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>::operator scalar() const
{
    return scalar(mpfr_get_d(v_, R));
}


template<unsigned P, mpfr_rnd_t R>
Ostream& operator<<
(
    Ostream& os,
    const multiPrecisionScalar<P,R>& v
)
{
    os << v.operator scalar();

    os.check
    (
        "Ostream& operator<<(Ostream& os, const multiPrecisionScalar& s)"
    );

    return os;
}


template<unsigned P, mpfr_rnd_t R>
Istream& operator>>
(
    Istream& is,
    multiPrecisionScalar<P,R>& v
)
{
    token t(is);

    if (t.isNumber())
    {
        v = t.number();
    }
    else
    if (t.isString())
    {
        // Convert string to number
        v = multiPrecisionScalar<P,R>(t.stringToken().c_str());

        // Flag empty strings as an error
        if (!t.stringToken().size())
        {
            is.setBad();

            FatalIOErrorIn("operator>>(Istream&, multiPrecisionScalar&)", is)
                << "wrong token type - "
                << "expected string found non-string characters "
                << t.info()
                << exit(FatalIOError);

            return is;
        }
    }
    else
    {
        is.setBad();

        FatalIOErrorIn("operator>>(Istream&, multiPrecisionScalar&)", is)
            << "wrong token type - expected multiPrecisionScalar found "
            << t.info()
            << exit(FatalIOError);

        return is;
    }

    is.check("Istream& operator>>(Istream& is, multiPrecisionScalar& s)");

    return is;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> operator+
(
    const multiPrecisionScalar<P,R>& a,
    const multiPrecisionScalar<P,R>& b
)
{
    return multiPrecisionScalar<P,R>(a) += b;
}


template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> operator+
(
    const scalar& a,
    const multiPrecisionScalar<P,R>& b
)
{
    return multiPrecisionScalar<P,R>(a) += b;
}


template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> operator+
(
    const multiPrecisionScalar<P,R>& a,
    const scalar& b
)
{
    return multiPrecisionScalar<P,R>(a) += multiPrecisionScalar<P,R>(b);
}


template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> operator-
(
    const multiPrecisionScalar<P,R>& a,
    const multiPrecisionScalar<P,R>& b
)
{
    return multiPrecisionScalar<P,R>(a) -= b;
}


template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> operator-
(
    const scalar& a,
    const multiPrecisionScalar<P,R>& b
)
{
    return multiPrecisionScalar<P,R>(a) -= b;
}


template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> operator-
(
    const multiPrecisionScalar<P,R>& a,
    const scalar& b
)
{
    return multiPrecisionScalar<P,R>(a) -= multiPrecisionScalar<P,R>(b);
}


template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> operator*
(
    const multiPrecisionScalar<P,R>& a,
    const multiPrecisionScalar<P,R>& b
)
{
    return multiPrecisionScalar<P,R>(a) *= b;
}


template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> operator*
(
    const scalar& a,
    const multiPrecisionScalar<P,R>& b
)
{
    return multiPrecisionScalar<P,R>(a) *= b;
}


template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> operator*
(
    const multiPrecisionScalar<P,R>& a,
    const scalar& b
)
{
    return multiPrecisionScalar<P,R>(a) *= multiPrecisionScalar<P,R>(b);
}


template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> operator/
(
    const multiPrecisionScalar<P,R>& a,
    const multiPrecisionScalar<P,R>& b
)
{
    return multiPrecisionScalar<P,R>(a) /= b;
}


template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> operator/
(
    const scalar& a,
    const multiPrecisionScalar<P,R>& b
)
{
    return multiPrecisionScalar<P,R>(a) /= b;
}


template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> operator/
(
    const multiPrecisionScalar<P,R>& a,
    const scalar& b
)
{
    return multiPrecisionScalar<P,R>(a) /= multiPrecisionScalar<P,R>(b);
}


template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> mag
(
    const multiPrecisionScalar<P,R>& v
)
{
    return v.mag();
}


template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> sqr
(
    const multiPrecisionScalar<P,R>& v
)
{
    return v.sqr();
}


template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> sqrt
(
    const multiPrecisionScalar<P,R>& v
)
{
    return v.sqrt();
}


} // End namespace Foam

// ************************************************************************* //
