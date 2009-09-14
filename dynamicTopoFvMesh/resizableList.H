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
    resizableList

Description
    Similar in intent to DynamicList, but without inheritance.

Author
    Sandeep Menon

SourceFiles
    resizableListI.H

\*---------------------------------------------------------------------------*/

#ifndef resizableList_H
#define resizableList_H

#include "label.H"
#include "bool.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
template<class T> class List;

template<class T> class resizableList;

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class T> Ostream& operator<<
(
    Ostream&,
    resizableList<T>&
);


/*---------------------------------------------------------------------------*\
                           Class resizableList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class resizableList
{
    // Private data

        //- Vector of values of type T.
        T* __restrict__ v_;

        // Size of the list
        label size_;

        // Keeps track of the last index
        label lastIndex_;

    // Private Member Functions

public:

    // Constructors

        //- Construct null
        resizableList();

        //- Construct given a size and value for all elements
        resizableList(const label size, const T& data);

        //- Construct from an existing List<T> structure
        resizableList(const List<T>& data);

        //- Construct as copy
        resizableList(const resizableList<T>& data);


    // Destructor

        ~resizableList();


    // Member Functions

        // Memory allocation

            //- Allocation for a null-constructed object, given size
            void setSize(const label size);

            //- Allocation for a null-constructed object, given size/value
            void setSize(const label size, const T& data);

            //- Append an element to the end of the list
            //- Returns the index of the added element
            inline label append(const T& data);

            //- Clear contents of the list
            void clear();

            //- Deallocate all memory
            void clearOut();

        // Access

            //- Return the allocated size of the list.
            inline label allocSize() const;

            //- Return the number of elements in the list.
            inline label size() const;

        // Check

            //- Return true if the list is empty (i.e., if size == 0).
            inline bool empty() const;

            //- Check the validity of an index
            inline bool checkIndex(const label index) const;

    // Member Operators

        resizableList<T>& operator=(const resizableList& data);

        void operator=(const List<T>& data);

        //- Return subscript-checked data reference of HashList.
        inline T& operator[](const label index);

        //- Return subscript-checked data reference of constant HashList.
        inline const T& operator[](const label index) const;

        //- Assignment of all entries to the given value
        inline void operator=(const T& data);

    // IOstream Operators

        friend Ostream& operator<< <T>
        (
            Ostream&,
            resizableList<T>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "resizableListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
