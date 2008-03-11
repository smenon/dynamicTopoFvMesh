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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
template<class T> class List;

template<class T> class HashList;

// Structure that holds one element of the linked-list
template<class T> struct Item {
    label   idx_;
    T       data_;
    // Next and previous item (in hash-table link-list)
    Item<T> *prev_, *next_;
    // Next and previous item (by index)
    Item<T> *previdx_, *nextidx_;
    // Default initialization
    Item<T>() {
        idx_     = -1;
        prev_    = next_    = 0;
        previdx_ = nextidx_ = 0;
    }
};

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

        // Size of the list
        label size_;

        // Pointer to element[0] and element[N] of the list
        Item<T> *startIdx_, *endIdx_;
        
	// Hash-table
        // Holds the pointer to the start of the linked-list
	Item<T> **table_;

        // Hash-table size. (Must be prime!)
        // Change only if collision rates are unmanageable
        static label hash_table_size_;
        
        // Rehash flag for automatic rehashing
        bool rehash_;
        
        // Index from which rehashing begins
        label rehashIndex_;

    // Private Member Functions

	// Hash to the table-entry
	inline label hash(label) const;
        
        //- Return an item of the HashList
        inline Item<T>* getLinkItem(const label);

public:

    // Static data members

    // Constructors

        //- Construct null
        HashList();

        //- Construct given a size and value for all elements
        HashList(const label, const T&);

        //- Construct from an existing List<T> structure
        HashList(const List<T>&);

        //- Construct as copy
        HashList(const HashList<T>&);
        
    // Iterator class definition
        
        class iterator {
            
            friend class HashList;
            
            // Private data
            Item<T> *curItem;
            
        public:
            
            // Construct null
            iterator() : curItem(0) {}
            
            // Constructor with input
            iterator(Item<T>* p) : curItem(p) {}
            
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
            
            //- Postfix overload
            inline iterator& operator++(int);
        };

    // Destructor

        ~HashList();

    // Member Functions
        
        // Memory allocation
        
            //- Allocation for a null-constructed object, given size/value
            void setSize(const label, const T&);
            
            //- Insertion of a new element at the given label 
            //- This will require a rehash
            inline iterator insert(const label, const T&);           
            
            //- Deletion of an element at the given label
            //- This will require a rehash
            inline void remove(const label);
            
            //- Append an element to the end of the list
            //- Returns the index of the added element
            inline label append(const T&);
            
            //- Clear contents of the list
            void clear();
            
            //- Rehash indices after insertion/deletion
            //inline void rehash();
            
            //- Sort the list in ascending order of indices
            void sort();
        
        // Access
                       
            //- Return const iterator to element[0]
            inline iterator begin();
            
            //- Return const iterator to element[N]
            inline iterator end(); 
            
            //- Rehash the index of an iterator to a new value
            inline void reNumber(const label, iterator&);

            //- Return the index at which rehashing begins
            inline label rehashIndex() const;
            
            //- Return the table-size of the HashList.
            inline label table_size() const;        

            //- Return the number of elements in the HashList.
            inline label size() const;

        // Check
            
            //- Return the rehash status of the list
            inline bool rehashRequired() const;
        
            //- Return true if the HashList is empty (i.e., if size() == 0).
            inline bool empty() const;  
            
            //- Check the validity of an index
            inline bool checkIndex(const label) const;

        // Output 
            
            //- Print out the link-list structure
            void toc() const;

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

// Initialize static member data

// Larger table-sizes probably improve access-rates on larger lists, but can be 
// wasteful in terms of memory. If collision rates are too high, either increase
// this value(preferred), or modify the hash-function.
template<class T> label HashList<T>::hash_table_size_ = 257;

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "HashListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
