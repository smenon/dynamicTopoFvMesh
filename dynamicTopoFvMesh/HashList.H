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
    HashList<T>

Description
    General templated class for dynamically-resizing one-dimensional lists
    of type T. Indexing is provided using a hash-function to achieve O(1)
    access time. List traversal is also achieved in linear time.
    The container is also thread-safe.

Author
    Sandeep Menon

SourceFiles
    HashListI.H

\*---------------------------------------------------------------------------*/

#ifndef HashList_H
#define HashList_H

#include "label.H"
#include "bool.H"
#include "List.H"

// Multi-threading support
#include "multiThreader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
template<class T> class List;

template<class T> class HashList;

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class T> Ostream& operator<<
(
    Ostream&,
    HashList<T>&
);


/*---------------------------------------------------------------------------*\
                           Class HashList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class HashList
{
    // Private data

        // Structure that holds one element of the linked-list
        struct Item
        {
            // Index of the item
            label   idx_;

            // Template data item
            T       data_;

            // Next and previous item (in hash-table link-list)
            Item    *prev_, *next_;

            // Next and previous item (by index)
            Item    *previdx_, *nextidx_;

            // Mutex for multi-threading
            Mutex   itemMutex_;

            // Default initialization
            Item()
            {
                idx_     = -1;
                prev_    = next_    = 0;
                previdx_ = nextidx_ = 0;
            }

            // Lock mutex
            inline void lock() { itemMutex_.lock(); }

            // Unlock mutex
            inline void unlock() { itemMutex_.unlock(); }
        };

        // Size of the list
        label size_;

        // Pointer to element[0] and element[N] of the list
        Item *startIdx_, *endIdx_;

        // Hash-table
        // Holds the pointer to the start of the linked-list
        Item **table_;

        // Mutex list for the HashTable
        Mutex *tableMutex_;

        // Global mutex for the entire list
        Mutex listMutex_;

        // Hash-table size. (Must be prime!)
        // Change if collision rates are unmanageable
        label hash_table_size_;

    // Private Member Functions

        // Hash to the table-entry
        inline label hash(label) const;

        //- Return an item of the HashList
        inline Item* getLinkItem(const label);

public:

    // Static data members

    // Constructors

        //- Construct null
        HashList(const label tableSize = 257);

        //- Construct given a size and value for all elements
        HashList(const label, const T&, const label tableSize = 257);

        //- Construct from an existing List<T> structure
        HashList(const List<T>&, const label tableSize = 257);

        //- Construct as copy
        HashList(const HashList<T>&);

    // Iterator class definition

        class iterator
        {
            friend class HashList;

            // Private data
            Item *curItem;

        public:

            // Construct null
            iterator() : curItem(0) {}

            // Constructor with input
            iterator(Item* p) : curItem(p) {}

            // Member functions

            //- Index access
            inline label index();

            // Member operators

            //- Assignment overload
            inline void operator=(const iterator&);

            //- Function overload
            inline T& operator()();

            //- Equality overloads
            inline bool operator!=(const iterator&) const;
            inline bool operator==(const iterator&) const;

            //- Postfix/prefix overload
            inline iterator& operator++(int);
            inline iterator& operator++();
        };

    // Destructor

        ~HashList();

    // Member Functions

        // Memory allocation

            //- Allocation for a null-constructed object, given size
            void setSize(const label, const label tableSize = 257);

            //- Allocation for a null-constructed object, given size/value
            void setSize(const label, const T&, const label tableSize = 257);

            //- Insertion of a new element at the given label
            inline iterator insert(const label, const T&);

            //- Deletion of an element at the given label
            inline void remove(const label);

            //- Append an element to the end of the list
            //- Returns the index of the added element
            inline label append(const T&);

            //- Clear contents of the list
            void clear();

            //- Sort the list in ascending order of indices
            void sort();

        // Access

            //- Return iterator to element[0]
            inline iterator begin();

            //- Return iterator to element[N]
            inline iterator end();

            //- Return nearest iterator to the specified index
            inline iterator getIterator(const label);

            //- Renumber the index of an iterator to a new value
            inline void reNumber(const label, iterator&);

            //- Return the index of the last element
            inline label lastIndex() const;

            //- Return the table-size of the HashList.
            inline label table_size() const;

            //- Return the number of elements in the HashList.
            inline label size() const;

            //- Lock the entire list [multi-threading]
            inline void lock();

            //- Unlock the entire list [multi-threading]
            inline void unlock();

        // Check

            //- Return true if the HashList is empty (i.e., if size() == 0).
            inline bool empty() const;

            //- Check the validity of an index
            inline bool checkIndex(const label) const;

            //- Check if an index exists
            inline bool found(const label) const;

        // Output

            //- Print out the link-list structure
            void toc() const;

            //- Print out the link-list at a particular location
            void toc(const label) const;

    // Member Operators

        HashList<T>& operator=(const HashList&);

        //- Return subscript-checked iterator of HashList.
        inline iterator operator()(const label);

        //- Return subscript-checked iterator of constant HashList.
        inline const iterator operator()(const label) const;

        //- Return subscript-checked data reference of HashList.
        inline T& operator[](const label);

        //- Return subscript-checked data reference of constant HashList.
        inline const T& operator[](const label) const;

        //- Assignment of all entries to the given value
        inline void operator=(const T&);

    // IOstream Operator

        friend Ostream& operator<< <T>
        (
            Ostream&,
            HashList<T>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "HashListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
