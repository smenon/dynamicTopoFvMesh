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
    Template implementations of the resizableList class

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "resizableList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
template<class T>
resizableList<T>::resizableList()
:
    v_(0),
    size_(0),
    lastIndex_(-1)
{}


//- Construct given a size and value for all elements
template<class T>
resizableList<T>::resizableList
(
    const label size,
    const T& data
)
:
    v_(0),
    size_(0),
    lastIndex_(-1)
{
    // Allocate storage
    setSize(size, data);

    // Set UList properties and lastIndex
    UList<T>::reset(v_, size_);
    lastIndex_ = (size_ - 1);
}


//- Construct from an existing UList<T> structure
template<class T>
resizableList<T>::resizableList
(
    const UList<T>& data
)
:
    v_(0),
    size_(0),
    lastIndex_(-1)
{
    if (data.size() == 0)
    {
        return;
    }

    // Allocate storage
    setSize(data.size());

    // Assign
    this->operator=(data);

    // Set UList properties and lastIndex
    UList<T>::reset(v_, size_);
    lastIndex_ = (size_ - 1);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

// Deallocate memory
template<class T>
resizableList<T>::~resizableList()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Allocation for a null-constructed object, given size
template<class T>
void resizableList<T>::setSize(const label size)
{
    if (size < 0)
    {
        FatalErrorIn
        (
            "resizableList<T>::resizableList(const label size)"
        ) << nl << " Cannot allocate to size " << size
          << abort(FatalError);
    }

    if (size)
    {
        if (size > size_)
        {
            T* nv = new T[label(size)];

            // Copy existing contents to a new list, if necessary.
            if (this->size_)
            {
                register label i = min(this->size_, size);

#               ifdef USEMEMCPY
                if (contiguous<T>())
                {
                    memcpy(nv, this->v_, i*sizeof(T));
                }
                else
#               endif
                {
                    register T* vv = &this->v_[i];
                    register T* av = &nv[i];
                    while (i--) *--av = *--vv;
                }

                delete[] this->v_;
            }

            // Take over member data
            v_ = nv;
            size_ = size;
        }

        // Set UList properties and lastIndex
        UList<T>::reset(v_, size);
        lastIndex_ = (size - 1);
    }
    else
    {
        clear();
    }
}


//- Construct given a size and value for all elements
template<class T>
void resizableList<T>::setSize(const label size, const T& data)
{
    label oldSize = this->size_;

    // Allocate storage
    setSize(size);

    // Padding for extra elements.
    if (size > oldSize)
    {
        register label i = size - oldSize;
        register T* vv = &this->v_[size];
        while (i--) *--vv = data;
    }

    // Set UList properties and lastIndex
    UList<T>::reset(v_, size);
    lastIndex_ = (size - 1);
}


//- Append an element to the end of the list
//- Returns the index of the added element
template<class T>
inline label resizableList<T>::append(const T& data)
{
    label index = (lastIndex_ + 1);

    if (index >= size_)
    {
        setSize
        (
            max
            (
                index,
                label(1.1*size_ + 1)
            )
        );
    }

    // The setSize operation changes lastIndex_ to (size_ - 1)
    // We don't want that to happen in this case.
    lastIndex_ = index;
    UList<T>::reset(v_, (lastIndex_ + 1));

    this->operator[](lastIndex_) = data;

    return lastIndex_;
}


//- Return the last element in the list
template<class T>
inline T resizableList<T>::remove()
{
    if (lastIndex_ == -1)
    {
        FatalErrorIn
        (
            "inline T resizableList<T>::remove() "
        ) << "List is empty." << abort(FatalError);
    }

    // Fetch the reference
    T& tRef = UList<T>::operator[](lastIndex_);

    // Reset UList parameters
    UList<T>::reset(v_, lastIndex_--);

    return tRef;
}


//- Clear contents of the list (i.e., reset lastIndex)
template<class T>
void resizableList<T>::clear()
{
    // Reset the lastIndex
    lastIndex_ = -1;

    // Reset UList parameters
    UList<T>::reset(v_, 0);
}


//- Deallocate all memory
template<class T>
void resizableList<T>::clearOut()
{
    // Reset UList parameters
    UList<T>::reset(0, 0);

    if (size_)
    {
        delete[] this->v_;
    }

    size_ = 0;
    lastIndex_ = -1;
}


//- Return the number of elements in the list.
template<class T>
inline label resizableList<T>::allocSize() const
{
    return size_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
resizableList<T>& resizableList<T>::operator=
(
    const resizableList& data
)
{
    if (data.size() < 0)
    {
        FatalErrorIn
        (
            "resizableList<T>::operator=(const resizableList<T>&)"
        ) << nl << " Cannot allocate to size " << data.size()
          << abort(FatalError);
    }

    if (data.size() != this->size_)
    {
        setSize(data.size());
    }

    if (this->size_)
    {
        lastIndex_ = -1;

        forAll(data, indexI)
        {
            v_[++lastIndex_] = data[indexI];
        }
    }

    // Reset UList parameters
    UList<T>::reset(v_, (lastIndex_ + 1));

    return (*this);
}


// Assigment operator to UList<T>
template<class T>
inline void resizableList<T>::operator=
(
    const UList<T>& data
)
{
    if (data.size() < 0)
    {
        FatalErrorIn
        (
            "resizableList<T>::operator=(const UList<T>&)"
        ) << nl << " Cannot allocate to size " << data.size()
          << abort(FatalError);
    }

    if (data.size() != this->size_)
    {
        setSize(data.size());
    }

    if (this->size_)
    {
        lastIndex_ = -1;

        forAll(data, indexI)
        {
            v_[++lastIndex_] = data[indexI];
        }
    }

    // Reset UList parameters
    UList<T>::reset(v_, (lastIndex_ + 1));
}


//- Assignment of all entries to the given value
template<class T>
inline void resizableList<T>::operator=
(
    const T& data
)
{
    if (size_)
    {
        for (label i = 0; i < size_; i++)
        {
            v_[i] = data;
        }

        lastIndex_ = (size_ - 1);
    }

    // Reset UList parameters
    UList<T>::reset(v_, (lastIndex_ + 1));
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class T>
Ostream& operator<<
(
    Ostream& os,
    resizableList<T>& L
)
{
    if (L.size() < 11)
    {
        // Write size of the list
        os << L.size() << token::BEGIN_LIST;

        // Write contents
        forAll(L, indexI)
        {
            if (indexI > 0) os << token::SPACE;
            os << L[indexI];
        }

        // Write end of contents
        os << token::END_LIST;
    }
    else
    {
        // Write size of list
        os << nl << L.size();

        // Write beginning of contents
        os << nl << token::BEGIN_LIST << nl;

        // Write contents
        forAll(L, indexI)
        {
            os << indexI << token::SPACE << L[indexI] << nl;
        }

        // Write end of contents
        os << token::END_LIST;
    }

    // Check state of IOstream
    os.check("Ostream& operator<<(Ostream&, const resizableList&)");

    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
