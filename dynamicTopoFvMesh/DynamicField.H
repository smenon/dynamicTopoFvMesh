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
    Foam::DynamicField

Description
    Dynamically sized Field.

SourceFiles
    DynamicFieldI.H
    DynamicField.C

\*---------------------------------------------------------------------------*/

#ifndef DynamicField_H
#define DynamicField_H

#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
class DynamicField;

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Ostream& operator<<
(
    Ostream&,
    const DynamicField<T, SizeInc, SizeMult, SizeDiv>&
);

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Ostream& operator<<
(
    Ostream&,
    const tmp<DynamicField<T, SizeInc, SizeMult, SizeDiv> >&
);

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Istream& operator>>
(
    Istream&,
    DynamicField<T, SizeInc, SizeMult, SizeDiv>&
);


/*---------------------------------------------------------------------------*\
                           Class DynamicField Declaration
\*---------------------------------------------------------------------------*/

template<class T, unsigned SizeInc=0, unsigned SizeMult=2, unsigned SizeDiv=1>
class DynamicField
:
    public Field<T>
{

    // Private data

        //- The capacity (allocated size) of the underlying field.
        label capacity_;

    // Private Member Functions

        //- Construct as copy of tmp<DynamicField>
#       ifdef ConstructFromTmp
        DynamicField
        (
            const tmp<DynamicField<T, SizeInc, SizeMult, SizeDiv> >&
        );
#       endif

        //- Construct from a dictionary entry
        DynamicField(const word&, const dictionary&, const label);

public:

    // Static data members

        static const char* const typeName;


    // Static Member Functions

        //- Return a null field
        inline static const DynamicField<T, SizeInc, SizeMult, SizeDiv>& null();


    // Constructors

        //- Construct null
        //  Used for temporary fields which are initialised after construction
        DynamicField();

        //- Construct given size
        //  Used for temporary fields which are initialised after construction
        explicit inline DynamicField(const label);

        //- Construct given size and initial value
        DynamicField(const label, const T&);

        //- Construct as copy of a UList\<T\>
        explicit inline DynamicField(const UList<T>&);

        //- Construct by 1 to 1 mapping from the given field
        inline DynamicField
        (
            const UList<T>& mapF,
            const labelList& mapAddressing
        );

        //- Construct by interpolative mapping from the given field
        inline DynamicField
        (
            const UList<T>& mapF,
            const labelListList& mapAddressing,
            const scalarListList& weights
        );

        //- Construct by mapping from the given field
        inline DynamicField
        (
            const UList<T>& mapF,
            const FieldMapper& map
        );

        //- Construct as copy
        inline DynamicField
        (
            const DynamicField<T, SizeInc, SizeMult, SizeDiv>&
        );

        //- Construct as copy or re-use as specified.
        inline DynamicField
        (
            DynamicField<T, SizeInc, SizeMult, SizeDiv>&,
            bool reUse
        );

        //- Construct from Istream
        inline DynamicField(Istream&);

        //- Clone
        tmp<DynamicField<T, SizeInc, SizeMult, SizeDiv> > clone() const;


    // Member Functions

        // Access

        //- Size of the underlying storage.
        inline label capacity() const;

        // Edit

        //- Alter the size of the underlying storage.
        //  The addressed size will be truncated if needed to fit, but will
        //  remain otherwise untouched.
        //  Use this or reserve() in combination with append().
        inline void setCapacity(const label);

        //- Alter the addressed list size.
        //  New space will be allocated if required.
        //  Use this to resize the list prior to using the operator[] for
        //  setting values (as per List usage).
        void setSize(const label);

        //- Alter the addressed list size and fill new space with a constant.
        inline void setSize(const label, const T&);

        //- Alter the addressed list size.
        //  New space will be allocated if required.
        //  Use this to resize the list prior to using the operator[] for
        //  setting values (as per List usage).
        inline void resize(const label);

        //- Alter the addressed list size and fill new space with a constant.
        inline void resize(const label, const T&);

        //- Reserve allocation space for at least this size.
        //  Never shrinks the allocated size, use setCapacity() for that.
        inline void reserve(const label);

        //- Clear the addressed list, i.e. set the size to zero.
        //  Allocated size does not change
        inline void clear();

        //- Clear the list and delete storage.
        inline void clearStorage();

        //- Shrink the allocated space to the number of elements used.
        //  Returns a reference to the DynamicList.
        inline DynamicField<T, SizeInc, SizeMult, SizeDiv>& shrink();

    // Member operators

        //- Append an element at the end of the list
        inline void append(const T&);

        //- Append a List at the end of this field
        inline void append(const UList<T>&);

        //- Remove and return the top element
        inline T remove();

        //- Return non-const access to an element, resizing list if necessary
        inline T& operator()(const label);

        //- Assignment of all addressed entries to the given value
        inline void operator=(const T&);

        //- Assignment from DynamicField
        inline void operator=
        (
            const DynamicField<T, SizeInc, SizeMult, SizeDiv>&
        );

        //- Assignment from UList
        inline void operator=(const UList<T>&);

        //- Assignment from tmp<DynamicField>
        inline void operator=
        (
            const tmp<DynamicField<T, SizeInc, SizeMult, SizeDiv> >&
        );

        //- Return element of Field.
        inline T& operator[](const label i);

        //- Return element of constant Field.
        inline const T& operator[](const label) const;

    // IOstream operators

        friend Ostream& operator<< <T, SizeInc, SizeMult, SizeDiv>
        (
            Ostream&,
            const DynamicField<T, SizeInc, SizeMult, SizeDiv>&
        );

        friend Ostream& operator<< <T, SizeInc, SizeMult, SizeDiv>
        (
            Ostream&,
            const tmp<DynamicField<T, SizeInc, SizeMult, SizeDiv> >&
        );

        friend Istream& operator>> <T, SizeInc, SizeMult, SizeDiv>
        (
            Istream&,
            DynamicField<T, SizeInc, SizeMult, SizeDiv>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DynamicFieldI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DynamicField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
