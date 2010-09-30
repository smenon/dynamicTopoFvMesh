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
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>::multiPrecisionScalar()
{
    // Initialize with precision
    mpfr_init2(v_, P);

    // Set to zero and rounding mode
    mpfr_set_zero(v_, 0);

#   ifdef FULLDEBUG
    dv_ = this->convertToWord(16);
#   endif
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

#   ifdef FULLDEBUG
    dv_ = this->convertToWord(16);
#   endif
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

    // Convert to a string intermediate
    std::ostringstream oss;
    oss.precision(16);
    oss << v;

    // Set string representation
    mpfr_set_str(v_, oss.str().c_str(), 10, R);

#   ifdef FULLDEBUG
    dv_ = this->convertToWord(16);
#   endif
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

#   ifdef FULLDEBUG
    dv_ = this->convertToWord(16);
#   endif
}


// Construct from input stream
template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>::multiPrecisionScalar(Istream& is)
:
    multiPrecisionScalar<P,R>()
{
    is >> *this;

#   ifdef FULLDEBUG
    dv_ = this->convertToWord(16);
#   endif
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>::~multiPrecisionScalar()
{
    mpfr_clear(v_);

#   ifdef FULLDEBUG
    dv_ = word("0.0");
#   endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> multiPrecisionScalar<P,R>::mag() const
{
    multiPrecisionScalar<P,R> result(*this);

    mpfr_abs(result.v_, result.v_, R);

#   ifdef FULLDEBUG
    result.dv_ = result.convertToWord(16);
#   endif

    return result;
}


template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> multiPrecisionScalar<P,R>::sqr() const
{
    multiPrecisionScalar<P,R> result(*this);

    mpfr_sqr(result.v_, result.v_, R);

#   ifdef FULLDEBUG
    result.dv_ = result.convertToWord(16);
#   endif

    return result;
}


template<unsigned P, mpfr_rnd_t R>
inline const multiPrecisionScalar<P,R> multiPrecisionScalar<P,R>::sqrt() const
{
    multiPrecisionScalar<P,R> result(*this);

    mpfr_sqrt(result.v_, result.v_, R);

#   ifdef FULLDEBUG
    result.dv_ = result.convertToWord(16);
#   endif

    return result;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>& multiPrecisionScalar<P,R>::operator=
(
    const scalar& v
)
{
    // Convert to a string intermediate
    std::ostringstream oss;
    oss.precision(16);
    oss << v;

    // Set string representation
    mpfr_set_str(v_, oss.str().c_str(), 10, R);

#   ifdef FULLDEBUG
    dv_ = this->convertToWord(16);
#   endif

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

#   ifdef FULLDEBUG
    dv_ = this->convertToWord(16);
#   endif

    return *this;
}


template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>& multiPrecisionScalar<P,R>::operator+=
(
    const multiPrecisionScalar<P,R>& v
)
{
    mpfr_add(v_, v_, v.v_, R);

#   ifdef FULLDEBUG
    dv_ = this->convertToWord(16);
#   endif

    return *this;
}


template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>& multiPrecisionScalar<P,R>::operator-=
(
    const multiPrecisionScalar<P,R>& v
)
{
    mpfr_sub(v_, v_, v.v_, R);

#   ifdef FULLDEBUG
    dv_ = this->convertToWord(16);
#   endif

    return *this;
}


template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>& multiPrecisionScalar<P,R>::operator*=
(
    const multiPrecisionScalar<P,R>& v
)
{
    mpfr_mul(v_, v_, v.v_, R);

#   ifdef FULLDEBUG
    dv_ = this->convertToWord(16);
#   endif

    return *this;
}


template<unsigned P, mpfr_rnd_t R>
inline multiPrecisionScalar<P,R>& multiPrecisionScalar<P,R>::operator/=
(
    const multiPrecisionScalar<P,R>& v
)
{
    mpfr_div(v_, v_, v.v_, R);

#   ifdef FULLDEBUG
    dv_ = this->convertToWord(16);
#   endif

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
inline const word multiPrecisionScalar<P,R>::convertToWord
(
    const size_t n,
    const int base
) const
{
    string out;
    size_t slen, nslen;
    char *s = NULL, *ns = NULL;

//    // Check for infinity
//    if (mpfr_inf_p(v_) != 0)
//    {
//        if (mpfr_sgn(v_) > 0)
//        {
//            return "+@Inf@";
//        }
//        else
//        {
//            return "-@Inf@";
//        }
//    }
//
//    // Check for zero
//    if (mpfr_zero_p(v_) != 0)
//    {
//        return "0";
//    }
//
//    // Check for NaN
//    if (mpfr_nan_p(v_) != 0)
//    {
//        return "@NaN@";
//    }

    // Fetch the string
    mp_exp_t exp;

    s = mpfr_get_str(NULL, &exp, base, 0, v_, R);
    ns = mpfr_get_str(NULL, &exp, base, n, v_, R);

    if (s != NULL)
    {
        slen = strlen(s);
        nslen = strlen(ns);

        if (nslen <= slen)
        {
            mpfr_free_str(s);

            s = ns;
            slen = nslen;
        }
        else
        {
            mpfr_free_str(ns);
        }

        if (exp > 0 && static_cast<size_t>(exp) < slen)
        {
            if (s[0] == '-')
            {
                // Remove zeros starting from right end
                char* ptr = s + slen - 1;

                while (*ptr == '0' && ptr > (s + exp))
                {
                    ptr--;
                }

                if (ptr == (s + exp))
                {
                    out = string(s, exp + 1);
                }
                else
                {
                    out =
                    (
                        string(s, exp + 1)
                      + '.'
                      + string(s + exp + 1, ptr - (s + exp + 1) + 1)
                    );
                }
            }
            else
            {
                // Remove zeros starting from right end
                char* ptr = s + slen - 1;

                while (*ptr == '0' && ptr > (s + exp - 1))
                {
                    ptr--;
                }

                if (ptr == (s + exp - 1))
                {
                    out = string(s, exp);
                }
                else
                {
                    out =
                    (
                        string(s, exp)
                      + '.'
                      + string(s + exp, ptr - (s + exp) + 1)
                    );
                }
            }
        }
        else
        {
            if (s[0] == '-')
            {
                // Remove zeros starting from right end
                char* ptr = s + slen - 1;

                while (*ptr == '0' && ptr > (s + 1))
                {
                    ptr--;
                }

                if (ptr == (s + 1))
                {
                    out = string(s, 2);
                }
                else
                {
                    out =
                    (
                        string(s, 2)
                      + '.'
                      + string(s + 2, ptr - (s + 2) + 1)
                    );
                }
            }
            else
            {
                // Remove zeros starting from right end
                char* ptr = s + slen - 1;

                while (*ptr=='0' && ptr>s)
                {
                    ptr--;
                }

                if (ptr == s)
                {
                    out = string(s,1);
                }
                else
                {
                    out =
                    (
                        string(s, 1)
                      + '.'
                      + string(s + 1, ptr - (s + 1) + 1)
                    );
                }
            }

            // Make final string
            if (--exp)
            {
                std::ostringstream oss;
                oss << exp << std::dec;

                if (exp > 0)
                {
                    out += "e+" + oss.str();
                }
                else
                {
                    out += "e" + oss.str();
                }
            }
        }

        mpfr_free_str(s);

        return word(out);
    }
    else
    {
        return "Conversion error!";
    }
}


template<unsigned P, mpfr_rnd_t R>
Ostream& operator<<
(
    Ostream& os,
    const multiPrecisionScalar<P,R>& v
)
{
    os << v.convertToWord(os.precision());

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

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class T1, class T2>
inline const Vector<T1> convert(const Vector<T2>& v)
{
    Vector<T1> tv;

    tv.x() = v.x();
    tv.y() = v.y();
    tv.z() = v.z();

    return tv;
}


template<class T1, class T2>
inline const Tensor<T1> convert(const Tensor<T2>& v)
{
    Tensor<T1> tv;

    tv.xx() = v.xx(); tv.xy() = v.xy(); tv.xz() = v.xz();
    tv.yx() = v.yx(); tv.yy() = v.yy(); tv.yz() = v.yz();
    tv.zx() = v.zx(); tv.zy() = v.zy(); tv.zz() = v.zz();

    return tv;
}


template<class T1, class T2>
inline const SymmTensor<T1> convert(const SymmTensor<T2>& v)
{
    SymmTensor<T1> tv;

    tv.xx() = v.xx(); tv.xy() = v.xy(); tv.xz() = v.xz();
                      tv.yy() = v.yy(); tv.yz() = v.yz();
                                        tv.zz() = v.zz();

    return tv;
}


template<class T1, class T2>
inline const SphericalTensor<T1> convert(const SphericalTensor<T2>& v)
{
    SphericalTensor<T1> tv;

    tv.ii() = v.ii();

    return tv;
}


template<class T1, class T2>
inline const Field<Vector<T1> > convert(const Field<Vector<T2> >& f)
{
    Field<Vector<T1> > tf(f.size());

    forAll(tf, indexI)
    {
        tf[indexI] = convert<T1>(f[indexI]);
    }

    return tf;
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
