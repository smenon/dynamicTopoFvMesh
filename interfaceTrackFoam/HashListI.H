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
    HashList

Description
    Template implementations of the HashList class

Author
    Sandeep Menon

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Hash to the table-entry (Thomas Wang's 7-shift integer hash)
template<class T>
inline label HashList<T>::hash(label key) const
{
    key -= (key<<6);
    key ^= (key>>17);
    key -= (key<<9);
    key ^= (key<<4);
    key -= (key<<3);
    key ^= (key<<10);
    key ^= (key>>15);
    return key;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
template<class T>
HashList<T>::HashList()
:
    size_(0),
    startIdx_(0),
    endIdx_(0),
    table_(0),
    rehash_(false),
    rehashIndex_(-1)
{}

//- Allocation for a null-constructed object, given size/value
template<class T>
void HashList<T>::setSize
(
    const label size,
    const T& data
)
{
    if (size<=0)
    {
        FatalErrorIn("HashList<T>::allocate(const size, const T&)") << nl
            << " Cannot allocate to size " << size
            << abort(FatalError);
    }

    // Initialize member data
    size_   = size;
    rehash_ = false;
    rehashIndex_ = size_;

    // Allocate the hash-table and initialize
    table_ = new Item<T>*[hash_table_size_];
    for(label i=0; i<hash_table_size_; i++)
    {
        table_[i] = 0;
    }

    // Allocate
    Item<T> *thisPtr=0, *prevPtr=0;
    for (label i=0; i<size_; i++)
    {
        thisPtr = new Item<T>;
        thisPtr->data_ = data;
        thisPtr->idx_ = i;

        // Update the hash-table
        label loc = hash(i)%hash_table_size_;
        if (table_[loc])
        {
            table_[loc]->prev_ = thisPtr;
            thisPtr->next_ = table_[loc];
        }
        table_[loc] = thisPtr;

        // Update the numerical indexing
        if (i==0)
        {
            startIdx_ = thisPtr;
        }
        else
        {
            prevPtr->nextidx_ = thisPtr;
            thisPtr->previdx_ = prevPtr;
        }
        prevPtr = thisPtr;
    }

    // Update the last element
    endIdx_=prevPtr;
}

//- Construct given a size and value for all elements
template<class T>
HashList<T>::HashList
(
    const label size,
    const T& data
)
{
    if (size<=0)
    {
        FatalErrorIn("HashList<T>::HashList(const size, const T&)") << nl
            << " Cannot allocate to size " << size
            << abort(FatalError);
    }

    // Initialize member data
    size_   = size;
    rehash_ = false;
    rehashIndex_ = size_;

    // Allocate the hash-table and initialize
    table_ = new Item<T>*[hash_table_size_];
    for(label i=0; i<hash_table_size_; i++)
    {
        table_[i] = 0;
    }

    // Allocate
    Item<T> *thisPtr=0, *prevPtr=0;
    for (label i=0; i<size_; i++)
    {
        thisPtr = new Item<T>;
        thisPtr->data_ = data;
        thisPtr->idx_  = i;

        // Update the hash-table
        label loc = hash(i)%hash_table_size_;
        if (table_[loc])
        {
            table_[loc]->prev_ = thisPtr;
            thisPtr->next_ = table_[loc];
        }
        table_[loc] = thisPtr;

        // Update the numerical indexing
        if (i==0)
        {
            startIdx_ = thisPtr;
        }
        else
        {
            prevPtr->nextidx_ = thisPtr;
            thisPtr->previdx_ = prevPtr;
        }
        prevPtr = thisPtr;
    }

    // Update the last element
    endIdx_=prevPtr;
}

//- Construct from an existing List<T> structure
template<class T>
HashList<T>::HashList
(
    const List<T>& data
)
{
    if (data.size()<=0)
    {
        FatalErrorIn("HashList<T>::HashList(const List<T>& data)") << nl
            << " Cannot allocate to size " << data.size()
            << abort(FatalError);
    }

    // Initialize member data
    size_   = data.size();
    rehash_ = false;
    rehashIndex_ = size_;

    // Allocate the hash-table and initialize
    table_ = new Item<T>*[hash_table_size_];
    for(label i=0; i<hash_table_size_; i++)
    {
        table_[i] = 0;
    }

    // Allocate
    Item<T> *thisPtr=0, *prevPtr=0;
    for (label i=0; i<size_; i++)
    {
        thisPtr = new Item<T>;
        thisPtr->data_ = data[i];
        thisPtr->idx_  = i;

        // Update the hash-table
        label loc = hash(i)%hash_table_size_;
        if (table_[loc])
        {
            table_[loc]->prev_ = thisPtr;
            thisPtr->next_ = table_[loc];
        }
        table_[loc] = thisPtr;

        // Update the numerical indexing
        if (i==0)
        {
            startIdx_ = thisPtr;
        }
        else
        {
            prevPtr->nextidx_ = thisPtr;
            thisPtr->previdx_ = prevPtr;
        }
        prevPtr = thisPtr;
    }

    // Update the last element
    endIdx_=prevPtr;
}

//- Construct as copy
template<class T>
HashList<T>::HashList
(
    const HashList<T>& data
)
{
    if (data.size()<=0)
    {
        FatalErrorIn("HashList<T>::HashList(const HashList<T>& data)") << nl
            << " Cannot allocate to size " << data.size()
            << abort(FatalError);
    }

    // Initialize member data
    size_   = data.size();
    rehash_ = data.rehashRequired();
    rehashIndex_ = data.rehashIndex();

    // Allocate the hash-table and initialize
    table_ = new Item<T>*[hash_table_size_];
    for(label i=0; i<hash_table_size_; i++)
    {
        table_[i] = 0;
    }

    // Allocate
    Item<T> *thisPtr=0, *prevPtr=0;
    for (label i=0; i<size_; i++)
    {
        thisPtr = new Item<T>;
        thisPtr->data_ = data[i];
        thisPtr->idx_  = i;

        // Update the hash-table
        label loc = hash(i)%hash_table_size_;
        if (table_[loc])
        {
            table_[loc]->prev_ = thisPtr;
            thisPtr->next_ = table_[loc];
        }
        table_[loc] = thisPtr;

        // Update the numerical indexing
        if (i==0)
        {
            startIdx_ = thisPtr;
        }
        else
        {
            prevPtr->nextidx_ = thisPtr;
            thisPtr->previdx_ = prevPtr;
        }
        prevPtr = thisPtr;
    }

    // Update the last element
    endIdx_=prevPtr;
}

//- Insertion of a new element at the given label
//- This will require a rehash() call
template<class T>
inline typename HashList<T>::iterator HashList<T>::insert
(
    const label index,
    const T& data
)
{
    // Hash the index
    label loc = hash(index)%hash_table_size_, rhIndex;
    rehash_   = true;

    // Loop through chained linked-list
    Item<T> *thisPtr=0, *newPtr=0, *searchPtr=0; 
    for (thisPtr=table_[loc]; thisPtr; thisPtr=thisPtr->next_)
    {
        if (index == thisPtr->idx_)
        {
            // Allocate 
            newPtr = new Item<T>;
            newPtr->data_ = data;

            // Update index pointers and set the rehash index
            if (thisPtr->previdx_)
            {
                thisPtr->previdx_->nextidx_ = newPtr;
                // Look for a previous index which is sane
                searchPtr=thisPtr;
                do
                {
                    searchPtr = searchPtr->previdx_;
                    rhIndex   = searchPtr->idx_;
                } while ((rhIndex == -1) && (searchPtr != startIdx_));
                rehashIndex_ = rhIndex < rehashIndex_ ? rhIndex : rehashIndex_;
            }
            newPtr->previdx_ = thisPtr->previdx_;
            newPtr->nextidx_ = thisPtr;
            thisPtr->previdx_ = newPtr;

            // Increment size and return pointer
            size_++;
            return iterator(newPtr);
        }
    }

    // The following line should never happen, just keeps gcc happy
#   ifdef FULLDEBUG
    FatalErrorIn("HashList<T>::insert(const label, const T&) ") << nl
            << " Index: " << index
            << " Attempt to access an element that doesn't exist. "
            << abort(FatalError);
#   endif
    return iterator(newPtr);
}

//- Deletion of an element at the given label
//- This will require a rehash() call
template<class T>
inline void HashList<T>::remove
(
    const label index
)
{
    label rhIndex;
    Item<T> *thisPtr=getLinkItem(index), *searchPtr=thisPtr;

    // Set the rehash flag
    rehash_ = true;

    // Update pointers, set the rehash index, and delete the item
    if (thisPtr->prev_)
    {
        thisPtr->prev_->next_ = thisPtr->next_;
    }
    else
    {
        // The table entry points to this item
        label loc = hash(index)%hash_table_size_;
        table_[loc] = thisPtr->next_;
    }
    if (thisPtr->next_)
    {
        thisPtr->next_->prev_ = thisPtr->prev_;
    }
    if (thisPtr->previdx_)
    {
        thisPtr->previdx_->nextidx_ = thisPtr->nextidx_;
        // Look for a previous index which is sane
        do
        {
            searchPtr = searchPtr->previdx_;
            rhIndex   = searchPtr->idx_;
        } while ((rhIndex == -1) && (searchPtr != startIdx_));
        rehashIndex_ = rhIndex < rehashIndex_ ? rhIndex : rehashIndex_;
    }
    else
    {
        // The first item is to be removed; update the first item
        rehashIndex_ = -1;
        startIdx_ = thisPtr->nextidx_;
    }
    if (thisPtr->nextidx_)
    {
        thisPtr->nextidx_->previdx_ = thisPtr->previdx_;
    }
    else
    {
        // The last item is being removed; update the last item
        endIdx_ = thisPtr->previdx_;
    }
    delete thisPtr;

    // Decrement size
    size_--;
}

//- Append an element to the end of the list
template<class T>
inline label HashList<T>::append
(
    const T& data
)
{
    // Index this item as an increment from the last index
    // This will allow random deletions from the list, while
    // keeping appends as valid
    label index=0;
    if (endIdx_) index = endIdx_->idx_+1;

    // Allocate and set the index
    Item<T> *newPtr = new Item<T>;
    newPtr->data_ = data;
    newPtr->idx_ = index;

    // Update the last element index
    if (endIdx_)
    {
        newPtr->previdx_ = endIdx_;
        endIdx_->nextidx_ = newPtr;
    }
    else
    {
        // A null endIdx_ occurs only in null-constructed HashLists
        startIdx_ = newPtr;
    }
    endIdx_ = newPtr;

    // Hash the new index, and update the table
    label loc = hash(index)%hash_table_size_;
    if (table_[loc])
    {
        table_[loc]->prev_ = newPtr;
        newPtr->next_ = table_[loc];
    }
    table_[loc] = newPtr;

    // Increment the size and return the index
    size_++;
    return index;
}

//- Clear contents of the list
template<class T>
void HashList<T>::clear()
{
    Item<T> *current, *node;

    // Delete all links sequentially
    current = startIdx_;
    while (current)
    {
        node = current->nextidx_;
        delete current;
        current = node;
    }

    // Clean out the hash table
    for(label i=0; i<hash_table_size_; i++)
    {
        table_[i] = 0;
    }
    size_ = 0; startIdx_ = 0; endIdx_ = 0;
    rehash_ = false; rehashIndex_ = -1;
}

//- Renumber (and rehash) indices after insertion/deletion
/*
template<class T>
inline void HashList<T>::rehash()
{
    if (this->rehash_)
    {
        // A rehashIndex less than zero implies that the entire list is to be rehashed
        Item<T> *thisPtr = (rehashIndex_ < 0) ? startIdx_ : getLinkItem(rehashIndex_);
        for(
             label idx = (thisPtr == startIdx_) ? 0 : rehashIndex_;
             thisPtr;
             thisPtr=thisPtr->nextidx_
           )
        {
            // Renumber the HashList in a sequential fashion
            thisPtr->idx_ = idx;

            // Remove the old links
            if (thisPtr->prev_) thisPtr->prev_->next_ = thisPtr->next_;
            if (thisPtr->next_) thisPtr->next_->prev_ = thisPtr->prev_;

            // Hash the new index, and update the table
            label loc = hash(idx)%hash_table_size_;
            if (table_[loc] == thisPtr)
            {
                table_[loc] = thisPtr->next_;
            }
            else
            if (table_[loc])
            {
                table_[loc]->prev_ = thisPtr;
                thisPtr->next_ = table_[loc];
            }
            thisPtr->prev_ = 0;
            table_[loc] = thisPtr;
            idx++;
        }

        // Reset the rehash flag and index
        rehash_      = false;
        rehashIndex_ = size_;
    }
}
*/

//- Sort the list in ascending order of indices
template<class T>
void HashList<T>::sort()
{
    // Set the starting index
    Item<T> *tmpItem = getLinkItem(0);
    tmpItem->previdx_ = 0;
    startIdx_ = tmpItem;

    // Loop through the rest of the elements and set accordingly
    for(label i=1; i<size_; i++)
    {
        Item<T> *item = getLinkItem(i);
        tmpItem->nextidx_ = item;
        item->previdx_    = tmpItem;
        tmpItem = item;
    }

    // Set the ending index
    tmpItem->nextidx_ = 0;
    endIdx_ = tmpItem;
}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

// Deallocate the hash-table structure and entity-list
template<class T>
HashList<T>::~HashList()
{
    if (table_)
    {
        clear();
        delete [] table_;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return an item of the HashList
template<class T>
inline Item<T>* HashList<T>::getLinkItem
(
    const label index
)
{
    // Hash the index
    label loc = hash(index)%hash_table_size_;

    // Loop through chained linked-list
    Item<T> *item=0;
    for (item=table_[loc]; item; item=item->next_)
    {
        if (index == item->idx_) return item;
    }

    // The following line should never happen, just keeps gcc happy
#   ifdef FULLDEBUG
    FatalErrorIn("HashList<T>::getLinkItem(const label index) ") << nl
            << " Index: " << index
            << " Attempt to access an element that doesn't exist. "
            << abort(FatalError);
#   endif
    return item;
}

//- Member function for the iterator
template<class T>
inline label HashList<T>::iterator::index()
{
    return curItem->idx_;
}

//- Rehash the index of an iterator to a new value
template<class T>
inline void HashList<T>::reNumber
(
    const label newIndex, iterator& iter
)
{
    // Hash the old index
    label oldloc = hash(iter.index())%hash_table_size_;

    // Reset the index for the item
    iter.curItem->idx_ = newIndex;
    label newloc = hash(newIndex)%hash_table_size_;

    // Change the hash-table pointers
    // Remove the old links
    if (iter.curItem->prev_)
        iter.curItem->prev_->next_ = iter.curItem->next_;
    if (iter.curItem->next_)
        iter.curItem->next_->prev_ = iter.curItem->prev_;

    // Update the table
    if (table_[oldloc] == iter.curItem)
        table_[oldloc] = iter.curItem->next_;
    if (table_[newloc])
    {
        table_[newloc]->prev_  = iter.curItem;
        iter.curItem->next_ = table_[newloc];
    }
    iter.curItem->prev_ = 0;
    table_[newloc] = iter.curItem;
}

//- Function overload - Return reference to data
template<class T>
inline T& HashList<T>::iterator::operator()()
{
    return curItem->data_;
}

// Iterator operator overload - Assignment 
template<class T>
inline void HashList<T>::iterator::operator=
(
    const iterator& rhs
)
{
    curItem = rhs.curItem;
}

// Iterator operator overload - Check for inequality
template<class T>
inline bool HashList<T>::iterator::operator!=
(
    const iterator& rhs
) const
{
    return !(curItem == rhs.curItem);
}

// Iterator operator overload - Check for equality
template<class T>
inline bool HashList<T>::iterator::operator==
(
    const iterator& rhs
) const
{
    return (curItem == rhs.curItem);
}

// Iterator operator overload - Increment to next item in the list
template<class T>
inline typename
HashList<T>::iterator&
HashList<T>::iterator::operator++
(
    int
)
{
    curItem = curItem->nextidx_;
    return (*this);
}

//- Return an iterator to the start of the list
template<class T>
inline typename HashList<T>::iterator HashList<T>::begin()
{
    return iterator(startIdx_);
}

//- Return an iterator to the end of the list
//- All lists end with NULL
template<class T>
inline typename HashList<T>::iterator HashList<T>::end()
{
    return iterator(0);
}

//- Return the index of the last element
template<class T>
inline label HashList<T>::lastIndex() const
{
    return endIdx_->idx_;
}

//- Return the index at which rehashing begins
template<class T>
inline label HashList<T>::rehashIndex() const
{
    return rehashIndex_;
}

//- Return the table-size of the HashList.
template<class T>
inline label HashList<T>::table_size() const
{
    return hash_table_size_;
}

//- Return the number of elements in the HashList.
template<class T>
inline label HashList<T>::size() const
{
    return size_;
}

//- Return the rehash status of the list
template<class T>
inline bool HashList<T>::rehashRequired() const
{
    return rehash_;
}

//- Return true if the HashList is empty (i.e., if size() == 0).
template<class T>
inline bool HashList<T>::empty() const
{
    return (size_ == 0);
}

//- Check the validity of an index
template<class T>
inline bool HashList<T>::checkIndex
(
    const label index
) const
{
    // Hash the index
    label loc = hash(index)%hash_table_size_;

    // Loop through chained linked-list
    for (Item<T> *item=table_[loc]; item; item=item->next_)
    {
        if (index == item->idx_) return true;
    }
    return false;
}

//- Print out the link-list structure
template<class T>
void HashList<T>::toc() const
{
    for(label i=0; i<hash_table_size_; i++)
    {
        Info << i << ": ";
        if (table_[i])
        {
            for(Item<T> *item = table_[i]; item; item=item->next_)
            {
                Info << item->idx_ << "->";
            }
        }
        Info << "NULL" << endl;
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
HashList<T>& HashList<T>::operator=
(
    const HashList& HLrhs
)
{
    Item<T> *lhs, *rhs;
    for(lhs=startIdx_, rhs=HLrhs.startIdx_;
        lhs;
        lhs=lhs->nextidx_, rhs=rhs->nextidx_)
    {
        lhs->data_ = rhs->data_;
    }
    return (*this);
}

//- Return subscript-checked iterator of HashList.
template<class T>
inline typename HashList<T>::iterator HashList<T>::operator()
(
    const label index
)
{
    // Hash the index
    label loc = hash(index)%hash_table_size_;

    // Loop through chained linked-list
    Item<T> *item=0;
    for (item=table_[loc]; item; item=item->next_)
    {
        if (index == item->idx_) return iterator(item);
    }

    // The following line should never happen, just keeps gcc happy
#   ifdef FULLDEBUG
    FatalErrorIn("HashList<T>::operator()(const label index) ") << nl
            << " Index: " << index
            << " Attempt to access an element that doesn't exist. "
            << abort(FatalError);
#   endif
    return iterator(item);
}

//- Return subscript-checked iterator of constant HashList.
template<class T>
inline const typename
HashList<T>::iterator
HashList<T>::operator()
(
    const label index
) const
{
    // Hash the index
    label loc = hash(index)%hash_table_size_;

    // Loop through chained linked-list
    Item<T> *item=0;
    for (item=table_[loc]; item; item=item->next_)
    {
        if (index == item->idx_) return iterator(item);
    }

    // The following line should never happen, just keeps gcc happy
#   ifdef FULLDEBUG
    FatalErrorIn("HashList<T>::operator()(const label index) const") << nl
            << " Index: " << index
            << " Attempt to access an element that doesn't exist. "
            << abort(FatalError);
#   endif
    return iterator(item);
}

//- Return element of HashList.
template<class T>
inline T& HashList<T>::operator[]
(
    const label index
)
{
    // Hash the index
    label loc = hash(index)%hash_table_size_;

    // Loop through chained linked-list
    Item<T> *item=0; 
    for (item=table_[loc]; item; item=item->next_)
    {
        if (index == item->idx_) return item->data_;
    }

    // The following line should never happen, just keeps gcc happy
#   ifdef FULLDEBUG
    FatalErrorIn("HashList<T>::operator[](const label index) ") << nl
        << " Index: " << index 
        << " Attempt to access an element that doesn't exist. "
        << abort(FatalError);
#   endif
    return item->data_;
}

//- Return element of constant HashList.
template<class T>
inline const T& HashList<T>::operator[]
(
    const label index
) const
{
    // Hash the index
    label loc = hash(index)%hash_table_size_;

    // Loop through chained linked-list
    Item<T> *item=0; 
    for (item=table_[loc]; item; item=item->next_)
    {
        if (index == item->idx_) return item->data_;
    }

    // The following line should never happen, just keeps gcc happy
#   ifdef FULLDEBUG
    FatalErrorIn("HashList<T>::operator[](const label index) const ") << nl
        << " Index: " << index 
        << " Attempt to access an element that doesn't exist. "
        << abort(FatalError);
#   endif
    return item->data_;
}

//- Assignment of all entries to the given value
template<class T>
inline void HashList<T>::operator=
(
    const T& rhs
)
{
    for(Item<T> *item=startIdx_; item; item=item->nextidx_) item->data_ = rhs;
    return (*this);
}

// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class T>
Ostream& operator<<(Ostream& os, HashList<T>& L)
{
    // Write size of HashList
    os << nl << L.size();

    // Write beginning of contents
    os << nl << token::BEGIN_LIST << nl;

    // Write HashList contents
    for(typename HashList<T>::iterator iter = L.begin(); iter != L.end(); iter++)
        os << iter.index() << token::SPACE << iter() << nl;

    // Write end of contents
    os << token::END_LIST;

    // Check state of IOstream
    os.check("Ostream& operator<<(Ostream&, const HashList&)");

    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
