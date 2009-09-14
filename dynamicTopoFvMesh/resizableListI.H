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
    // Allocate
    setSize(size, data);
}

//- Construct from an existing List<T> structure
template<class T>
resizableList<T>::resizableList
(
    const List<T>& data
)
:
    v_(0),
    size_(0),
    lastIndex_(-1)
{
    if (data.size() < 0)
    {
        FatalErrorIn("resizableList<T>::resizableList(const List<T>&)") << nl
            << " Cannot allocate to size " << data.size()
            << abort(FatalError);
    }

    if (data.size() == 0)
    {
        return;
    }

    // Assign
    this->operator=(data);
}

//- Construct as copy
template<class T>
resizableList<T>::resizableList
(
    const resizableList<T>& data
)
:
    v_(0),
    size_(0),
    lastIndex_(-1)
{
    if (data.size() < 0)
    {
        FatalErrorIn
        (
            "resizableList<T>::resizableList(const resizableList<T>&)"
        ) << nl << " Cannot allocate to size " << data.size()
          << abort(FatalError);
    }

    if (data.size_)
    {
        size_ = data.size_;

        this->v_ = new T[this->size_];

        forAll(data, indexI)
        {
            v_[++lastIndex_] = data[indexI];
        }
    }
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

        lastIndex_ = size - 1;
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
    this->setSize(size);

    // Padding for extra elements.
    if (size > oldSize)
    {
        register label i = size - oldSize;
        register T* vv = &this->v_[size];
        while (i--) *--vv = data;
    }
}

//- Append an element to the end of the list
//- Returns the index of the added element
template<class T>
inline label resizableList<T>::append(const T& data)
{
    label index = ++lastIndex_;

    if (index >= size_)
    {
        setSize
        (
            max
            (
                index,
                label(2*size_)
            )
        );
    }

    // The setSize operation changes lastIndex_ to (size_ - 1)
    // We don't want that to happen in this case.
    lastIndex_ = index;

    this->operator[](lastIndex_) = data;

    return index;
}

//- Clear contents of the list (i.e., reset lastIndex)
template<class T>
void resizableList<T>::clear()
{
    lastIndex_ = -1;
}

//- Deallocate all memory
template<class T>
void resizableList<T>::clearOut()
{
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

//- Return the number of elements in the list.
template<class T>
inline label resizableList<T>::size() const
{
    return (lastIndex_ + 1);
}

//- Return true if the list is empty (i.e., if size == 0).
template<class T>
inline bool resizableList<T>::empty() const
{
    return (size() == 0);
}

//- Check the validity of an index
template<class T>
inline bool resizableList<T>::checkIndex(const label index) const
{
    // Out-of-bounds check
    if (!size_)
    {
        FatalErrorIn("resizableList<T>::checkIndex(const label)")
            << "attempt to access element from zero sized list"
            << abort(FatalError);
    }
    else
    if (index < 0 || index >= size_)
    {
        FatalErrorIn("resizableList<T>::checkIndex(const label)")
            << "index " << index << " out of range 0 ... " << (size_-1)
            << abort(FatalError);
    }

    return true;
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

    return (*this);
}

// Assigment operator to List<T>
template<class T>
void resizableList<T>::operator=
(
    const List<T>& data
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
}

//- Return subscript-checked data reference of HashList.
template<class T>
inline T& resizableList<T>::operator[]
(
    const label index
)
{
#   ifdef FULLDEBUG
    checkIndex(index);
#   endif
    return v_[index];
}

//- Return subscript-checked data reference of constant HashList.
template<class T>
inline const T& resizableList<T>::operator[]
(
    const label index
) const
{
#   ifdef FULLDEBUG
    checkIndex(index);
#   endif
    return v_[index];
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
