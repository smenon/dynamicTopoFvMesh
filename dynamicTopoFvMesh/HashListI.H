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

#include "HashList.H"


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
HashList<T>::HashList(const label tableSize)
:
    size_(0),
    startIdx_(0),
    endIdx_(0),
    table_(0)
{
    initTable(tableSize);
}

//- Allocation for a null-constructed object, given size
//  Assumes that a null constructor exists for T
template<class T>
void HashList<T>::setSize
(
    const label size,
    const label tableSize
)
{
    if (size < 0)
    {
        FatalErrorIn("HashList<T>::setSize(size, tableSize)") << nl
            << " Cannot allocate to size " << size
            << abort(FatalError);
    }

    if (size == 0)
    {
        clear();
        return;
    }

    if (tableSize <= 0)
    {
        FatalErrorIn("HashList<T>::HashList(size, tableSize)") << nl
            << " Cannot allocate table to " << tableSize
            << abort(FatalError);
    }

    // Allocate the hash-table and initialize
    if (table_ && (tableSize != hash_table_size_))
    {
        delete [] table_;
        table_ = new Item*[tableSize];

        // Allocate the mutex-list as well
        delete [] tableMutex_;
        tableMutex_ = new Mutex[tableSize];

        hash_table_size_ = tableSize;
    }

    // Clear out existing contents
    if (size_)
    {
        clear();
    }

    // Initialize member data
    size_ = size;

    // Allocate
    Item *thisPtr=0, *prevPtr=0;
    for (label i=0; i<size_; i++)
    {
        thisPtr = new Item;
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

//- Allocation for a null-constructed object, given size/value
template<class T>
void HashList<T>::setSize
(
    const label size,
    const T& data,
    const label tableSize
)
{
    if (size < 0)
    {
        FatalErrorIn("HashList<T>::setSize(size, data, tableSize)") << nl
            << " Cannot allocate to size " << size
            << abort(FatalError);
    }

    if (size == 0)
    {
        clear();
        return;
    }

    if (tableSize <= 0)
    {
        FatalErrorIn("HashList<T>::HashList(size, data, tableSize)") << nl
            << " Cannot allocate table to " << tableSize
            << abort(FatalError);
    }

    // Allocate the hash-table and initialize
    if (table_ && (tableSize != hash_table_size_))
    {
        delete [] table_;
        table_ = new Item*[tableSize];

        // Allocate the mutex-list as well
        delete [] tableMutex_;
        tableMutex_ = new Mutex[tableSize];

        hash_table_size_ = tableSize;
    }

    // Clear out existing contents
    if (size_)
    {
        clear();
    }

    // Initialize member data
    size_ = size;

    // Allocate
    Item *thisPtr=0, *prevPtr=0;
    for (label i=0; i<size_; i++)
    {
        thisPtr = new Item;
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
    const T& data,
    const label tableSize
)
:
    size_(0),
    startIdx_(0),
    endIdx_(0),
    table_(0)
{
    // Initialize
    initTable(tableSize);

    // Allocate
    setSize(size, data, tableSize);
}

//- Construct from an existing List<T> structure
template<class T>
HashList<T>::HashList
(
    const List<T>& data,
    const label tableSize
)
:
    size_(0),
    startIdx_(0),
    endIdx_(0),
    table_(0)
{
    initTable(tableSize);

    if (data.size() < 0)
    {
        FatalErrorIn("HashList<T>::HashList(const List<T>&)") << nl
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
HashList<T>::HashList
(
    const HashList<T>& data
)
:
    size_(0),
    startIdx_(0),
    endIdx_(0),
    table_(0)
{
    initTable(data.table_size());

    if (data.size() < 0)
    {
        FatalErrorIn("HashList<T>::HashList(const HashList<T>&)") << nl
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

//- Insertion of a new element at the given label
template<class T>
inline typename HashList<T>::iterator HashList<T>::insert
(
    const label index,
    const T& data
)
{
    // Hash the index
    label loc = hash(index)%hash_table_size_;

    // Lock the list
    lock();

    // Lock the table at the index
    // tableMutex_[loc].lock();

    // Loop through chained linked-list
    Item *thisPtr=0, *newPtr=0;
    for (thisPtr=table_[loc]; thisPtr; thisPtr=thisPtr->next_)
    {
        if (index == thisPtr->idx_)
        {
            // Allocate
            newPtr = new Item;
            newPtr->data_ = data;

            // Update index pointers
            if (thisPtr->previdx_)
            {
                // Lock item
                thisPtr->previdx_->lock();

                thisPtr->previdx_->nextidx_ = newPtr;

                // Unlock item
                thisPtr->previdx_->unlock();
            }

            newPtr->previdx_ = thisPtr->previdx_;
            newPtr->nextidx_ = thisPtr;
            thisPtr->previdx_ = newPtr;

            // Increment size and return pointer
            size_++;

            // Unlock the table at the index
            // tableMutex_[loc].unlock();

            // Unlock the list
            unlock();

            return iterator(newPtr);
        }
    }

    // The following line should never happen, just keeps gcc happy
    FatalErrorIn("HashList<T>::insert(const label, const T&) ") << nl
            << " Index: " << index
            << " Attempt to access an element that doesn't exist. "
            << abort(FatalError);

    return iterator(newPtr);
}

//- Deletion of an element at the given label
template<class T>
inline void HashList<T>::remove
(
    const label index
)
{
    Item *thisPtr=getLinkItem(index);

    // Lock the list
    lock();

    // Update pointers and delete the item
    if (thisPtr->prev_)
    {
        // Lock items
        thisPtr->lock();
        thisPtr->prev_->lock();

        thisPtr->prev_->next_ = thisPtr->next_;

        // Unlock items
        thisPtr->prev_->unlock();
        thisPtr->unlock();
    }
    else
    {
        // The table entry points to this item
        label loc = hash(index)%hash_table_size_;

        // Lock table
        // tableMutex_[loc].lock();

        table_[loc] = thisPtr->next_;

        // Unlock table
        // tableMutex_[loc].unlock();
    }

    if (thisPtr->next_)
    {
        // Lock items
        thisPtr->lock();
        thisPtr->next_->lock();

        thisPtr->next_->prev_ = thisPtr->prev_;

        // Unlock items
        thisPtr->next_->unlock();
        thisPtr->unlock();
    }

    if (thisPtr->previdx_)
    {
        // Lock items
        thisPtr->lock();
        thisPtr->previdx_->lock();

        thisPtr->previdx_->nextidx_ = thisPtr->nextidx_;

        // Unlock items
        thisPtr->previdx_->unlock();
        thisPtr->unlock();
    }
    else
    {
        // The first item is to be removed; update the first item
        startIdx_ = thisPtr->nextidx_;
    }

    if (thisPtr->nextidx_)
    {
        // Lock items
        thisPtr->lock();
        thisPtr->nextidx_->lock();

        thisPtr->nextidx_->previdx_ = thisPtr->previdx_;

        // Unlock items
        thisPtr->nextidx_->unlock();
        thisPtr->unlock();
    }
    else
    {
        // The last item is being removed; update the last item
        endIdx_ = thisPtr->previdx_;
    }

    delete thisPtr;

    // Decrement size
    size_--;

    // Unlock the list
    unlock();
}

//- Deletion of an element, given an iterator
template<class T>
inline void HashList<T>::remove(iterator& iter)
{
    Item *thisPtr = iter.curItem;
    label index = thisPtr->idx_;

    // Lock the list
    lock();

    // Update pointers and delete the item
    if (thisPtr->prev_)
    {
        // Lock items
        thisPtr->lock();
        thisPtr->prev_->lock();

        thisPtr->prev_->next_ = thisPtr->next_;

        // Unlock items
        thisPtr->prev_->unlock();
        thisPtr->unlock();
    }
    else
    {
        // The table entry points to this item
        label loc = hash(index)%hash_table_size_;

        // Lock table
        // tableMutex_[loc].lock();

        table_[loc] = thisPtr->next_;

        // Unlock table
        // tableMutex_[loc].unlock();
    }

    if (thisPtr->next_)
    {
        // Lock items
        thisPtr->lock();
        thisPtr->next_->lock();

        thisPtr->next_->prev_ = thisPtr->prev_;

        // Unlock items
        thisPtr->next_->unlock();
        thisPtr->unlock();
    }

    if (thisPtr->previdx_)
    {
        // Lock items
        thisPtr->lock();
        thisPtr->previdx_->lock();

        thisPtr->previdx_->nextidx_ = thisPtr->nextidx_;

        // Unlock items
        thisPtr->previdx_->unlock();
        thisPtr->unlock();
    }
    else
    {
        // The first item is to be removed; update the first item
        startIdx_ = thisPtr->nextidx_;
    }

    if (thisPtr->nextidx_)
    {
        // Lock items
        thisPtr->lock();
        thisPtr->nextidx_->lock();

        thisPtr->nextidx_->previdx_ = thisPtr->previdx_;

        // Unlock items
        thisPtr->nextidx_->unlock();
        thisPtr->unlock();
    }
    else
    {
        // The last item is being removed; update the last item
        endIdx_ = thisPtr->previdx_;
    }

    delete thisPtr;

    // Decrement size
    size_--;

    // Unlock the list
    unlock();
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

    // Lock the list
    lock();

    label index = 0;

    // Allocate and set the index
    Item *newPtr = new Item;
    newPtr->data_ = data;

    // Update the last element index
    if (endIdx_)
    {
        // Lock end item
        endIdx_->lock();
        
        // Update item
        index = endIdx_->idx_ + 1;
        newPtr->idx_ = index;
        newPtr->previdx_ = endIdx_;
        endIdx_->nextidx_ = newPtr;

        // Unlock end item
        endIdx_->unlock();
    }
    else
    {
        // A null endIdx_ occurs only in null-constructed HashLists
        newPtr->idx_ = index;
        startIdx_ = newPtr;
    }

    // Update the last index item
    endIdx_ = newPtr;

    // Hash the new index, and update the table
    label loc = hash(index)%hash_table_size_;

    // Lock table index
    // tableMutex_[loc].lock();

    if (table_[loc])
    {
        // Lock item
        table_[loc]->lock();

        table_[loc]->prev_ = newPtr;
        newPtr->next_ = table_[loc];

        // Unlock item
        table_[loc]->unlock();
    }

    // Start table with the new item
    table_[loc] = newPtr;

    // Unlock table index
    // tableMutex_[loc].unlock();

    // Increment the size and return the index
    size_++;

    // Unlock the list
    unlock();

    return index;
}

//- Append a null element to the end of the list
//- Valid only for template arguments that implement a null constructor.
template<class T>
inline label HashList<T>::append()
{
    // Index this item as an increment from the last index

    // Lock the list
    lock();

    label index = 0;

    // Allocate and set the index
    Item *newPtr = new Item;

    // Template data T is constructed as null.

    // Update the last element index
    if (endIdx_)
    {
        // Lock end item
        endIdx_->lock();

        // Update item
        index = endIdx_->idx_ + 1;
        newPtr->idx_ = index;
        newPtr->previdx_ = endIdx_;
        endIdx_->nextidx_ = newPtr;

        // Unlock end item
        endIdx_->unlock();
    }
    else
    {
        // A null endIdx_ occurs only in null-constructed HashLists
        newPtr->idx_ = index;
        startIdx_ = newPtr;
    }

    // Update the last index item
    endIdx_ = newPtr;

    // Hash the new index, and update the table
    label loc = hash(index)%hash_table_size_;

    // Lock table index
    // tableMutex_[loc].lock();

    if (table_[loc])
    {
        // Lock item
        table_[loc]->lock();

        table_[loc]->prev_ = newPtr;
        newPtr->next_ = table_[loc];

        // Unlock item
        table_[loc]->unlock();
    }

    // Start table with the new item
    table_[loc] = newPtr;

    // Unlock table index
    // tableMutex_[loc].unlock();

    // Increment the size and return the index
    size_++;

    // Unlock the list
    unlock();

    return index;
}

//- Clear contents of the list
template<class T>
void HashList<T>::clear()
{
    Item *current, *node;

    // Delete all links sequentially
    current = startIdx_;
    while (current)
    {
        node = current->nextidx_;
        delete current;
        current = node;
    }

    // Clean the hash table
    cleanTable();

    size_ = 0; startIdx_ = 0; endIdx_ = 0;
}

//- Sort the list in ascending order of indices
template<class T>
void HashList<T>::sort()
{
    // Set the starting index
    Item *tmpItem = getLinkItem(0);
    tmpItem->previdx_ = 0;
    startIdx_ = tmpItem;

    // Loop through the rest of the elements and set accordingly
    for(label i=1; i<size_; i++)
    {
        Item *item = getLinkItem(i);
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
        delete [] tableMutex_;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return an item of the HashList
template<class T>
inline typename HashList<T>::Item* HashList<T>::getLinkItem
(
    const label index
)
{
    // Hash the index
    label loc = hash(index)%hash_table_size_;

    // Lock the table index
    // tableMutex_[loc].lock();

    // Loop through chained linked-list
    Item *item=0;
    for (item=table_[loc]; item; item=item->next_)
    {
        if (index == item->idx_)
        {
            // Unlock the table index
            // tableMutex_[loc].unlock();

            return item;
        }
    }

    // The following line should never happen, just keeps gcc happy
    FatalErrorIn("HashList<T>::getLinkItem(const label index) ") << nl
            << " Index: " << index
            << " Attempt to access an element that doesn't exist. "
            << abort(FatalError);

    return item;
}

//- Initialize the hash-table
template<class T>
inline void HashList<T>::initTable
(
    const label tableSize
)
{
    if (tableSize <= 0)
    {
        FatalErrorIn("HashList<T>::HashList(tableSize)") << nl
            << " Cannot allocate table to " << tableSize
            << abort(FatalError);
    }
    else
    {
        hash_table_size_ = tableSize;
    }

    // Allocate the hash-table
    table_ = new Item*[hash_table_size_];

    // Clean the hash table
    cleanTable();

    // Allocate the mutex-list as well
    tableMutex_ = new Mutex[hash_table_size_];
}

//- Clean out the hash-table
template<class T>
inline void HashList<T>::cleanTable()
{
    for(label i=0; i<hash_table_size_; i++)
    {
        table_[i] = 0;
    }
}

//- Member function for the iterator
template<class T>
inline label HashList<T>::iterator::index()
{
    return curItem->idx_;
}

//- Lock the iterator
template<class T>
inline void HashList<T>::iterator::lock()
{
    return curItem->lock();
}

//- Unlock the iterator
template<class T>
inline void HashList<T>::iterator::unlock()
{
    return curItem->unlock();
}

//- Renumber the index of an iterator to a new value
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
    {
        iter.curItem->prev_->next_ = iter.curItem->next_;
    }

    if (iter.curItem->next_)
    {
        iter.curItem->next_->prev_ = iter.curItem->prev_;
    }

    // Update the table
    if (table_[oldloc] == iter.curItem)
    {
        table_[oldloc] = iter.curItem->next_;
    }

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

// Iterator operator overload - Increment to next item in the list
template<class T>
inline typename
HashList<T>::iterator&
HashList<T>::iterator::operator++()
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

//- Return nearest iterator to the specified index
template<class T>
inline typename HashList<T>::iterator HashList<T>::getIterator
(
    const label index
)
{
    label nearestIndex = index;
    while (!checkIndex(nearestIndex))
    {
        nearestIndex += 1;
    }

    return iterator(getLinkItem(nearestIndex));
}

//- Return the index of the last element
template<class T>
inline label HashList<T>::lastIndex() const
{
    return endIdx_->idx_;
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

//- Lock the entire list [multi-threading]
template<class T>
inline void HashList<T>::lock()
{
    listMutex_.lock();
}

//- Unlock the entire list [multi-threading]
template<class T>
inline void HashList<T>::unlock()
{
    listMutex_.unlock();
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
    // Out-of-bounds check
    if (!size_)
    {
        FatalErrorIn("HashList<T>::checkIndex(const label)")
            << "attempt to access element from zero sized list"
            << abort(FatalError);
    }
    else
    if (index < 0 || index >= size_)
    {
        FatalErrorIn("HashList<T>::checkIndex(const label)")
            << "index " << index << " out of range 0 ... " << size_-1
            << abort(FatalError);
    }

    return found(index);
}

//- Check if an index exists
template<class T>
inline bool HashList<T>::found
(
    const label index
) const
{
    // Hash the index
    label loc = hash(index)%hash_table_size_;

    // Lock the table index
    // tableMutex_[loc].lock();

    // Loop through chained linked-list
    for (Item *item=table_[loc]; item; item=item->next_)
    {
        if (index == item->idx_)
        {
            // Unlock the table index
            // tableMutex_[loc].unlock();

            return true;
        }
    }

    // Unlock the table index
    // tableMutex_[loc].unlock();

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
            for(Item *item = table_[i]; item; item=item->next_)
            {
                Info << item->idx_ << "->";
            }
        }

        Info << "NULL" << endl;
    }
}

//- Print out the link-list at a particular location
template<class T>
void HashList<T>::toc(const label loc) const
{
    Info << loc << ": ";

    if (table_[loc])
    {
        for(Item *item = table_[loc]; item; item=item->next_)
        {
            Info << item->idx_ << "->";
        }
    }

    Info << "NULL" << endl;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
HashList<T>& HashList<T>::operator=
(
    const HashList& data
)
{
    // Allocate the hash-table and initialize
    if (table_ && (hash_table_size_ != data.table_size()))
    {
        delete [] table_;
        table_ = new Item*[data.table_size()];

        // Allocate the mutex-list as well
        delete [] tableMutex_;
        tableMutex_ = new Mutex[data.table_size()];

        hash_table_size_ = data.table_size();
    }

    // Clear existing contents
    if (size_)
    {
        clear();
    }
    else
    {
        // Clean the hash table
        cleanTable();
    }

    // Initialize member data
    size_ = data.size();

    Item *thisPtr=0, *prevPtr=0;
    for (label i=0; i<size_; i++)
    {
        thisPtr = new Item;
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

    return (*this);
}

// Assigment operator to List<T>
template<class T>
void HashList<T>::operator=
(
    const List<T>& data
)
{
    // Hash table should be allocated at this point
    // Throw an error if otherwise
    if (hash_table_size_)
    {
        // Clear existing contents
        if (size_)
        {
            clear();
        }
        else
        {
            // Clean the hash table
            cleanTable();
        }

        // Reinitialize the size
        size_ = data.size();
    }
    else
    {
        FatalErrorIn("HashList<T>::operator=(const List<T>&)") << nl
            << " Hash table has not been allocated. "
            << abort(FatalError);
    }

    Item *thisPtr=0, *prevPtr=0;
    for (label i=0; i<size_; i++)
    {
        thisPtr = new Item;
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

//- Return subscript-checked iterator of HashList.
template<class T>
inline typename HashList<T>::iterator HashList<T>::operator()
(
    const label index
)
{
    // Hash the index
    label loc = hash(index)%hash_table_size_;

    // Lock the table index
    // tableMutex_[loc].lock();

    // Loop through chained linked-list
    Item *item=0;
    for (item=table_[loc]; item; item=item->next_)
    {
        if (index == item->idx_)
        {
            // Unlock the table index
            // tableMutex_[loc].unlock();

            return iterator(item);
        }
    }

    // The following line should never happen, just keeps gcc happy
    FatalErrorIn("HashList<T>::operator()(const label index) ") << nl
            << " Index: " << index
            << " Attempt to access an element that doesn't exist. "
            << abort(FatalError);

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

    // Lock the table index
    // tableMutex_[loc].lock();

    // Loop through chained linked-list
    Item *item=0;
    for (item=table_[loc]; item; item=item->next_)
    {
        if (index == item->idx_)
        {
            // Unlock the table index
            // tableMutex_[loc].unlock();

            return iterator(item);
        }
    }

    // The following line should never happen, just keeps gcc happy
    FatalErrorIn("HashList<T>::operator()(const label index) const") << nl
            << " Index: " << index
            << " Attempt to access an element that doesn't exist. "
            << abort(FatalError);

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

    // Lock the table index
    // tableMutex_[loc].lock();

    // Loop through chained linked-list
    Item *item=0;
    for (item=table_[loc]; item; item=item->next_)
    {
        if (index == item->idx_)
        {
            // Unlock the table index
            // tableMutex_[loc].unlock();

            return item->data_;
        }
    }

    // The following line should never happen, just keeps gcc happy
    FatalErrorIn("HashList<T>::operator[](const label index) ") << nl
            << " Index: " << index
            << " Attempt to access an element that doesn't exist. "
            << abort(FatalError);

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

    // Lock the table index
    // tableMutex_[loc].lock();

    // Loop through chained linked-list
    Item *item=0;
    for (item=table_[loc]; item; item=item->next_)
    {
        if (index == item->idx_)
        {
            // Unlock the table index
            // tableMutex_[loc].unlock();

            return item->data_;
        }
    }

    // The following line should never happen, just keeps gcc happy
    FatalErrorIn("HashList<T>::operator[](const label index) const ") << nl
            << " Index: " << index
            << " Attempt to access an element that doesn't exist. "
            << abort(FatalError);

    return item->data_;
}

//- Assignment of all entries to the given value
template<class T>
inline void HashList<T>::operator=
(
    const T& rhs
)
{
    for
    (
        Item *item=startIdx_;
        item;
        item=item->nextidx_
    )
    {
        item->data_ = rhs;
    }

    return (*this);
}

// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class T>
Ostream& operator<<(Ostream& os, HashList<T>& L)
{
    if (L.size() < 11)
    {
        // Write size of HashList
        os << L.size() << token::BEGIN_LIST;

        // Write HashList contents
        for
        (
            typename HashList<T>::iterator iter = L.begin();
            iter != L.end();
            iter++
        )
        {
            if (iter.index() > 0) os << token::SPACE;
            os << iter();
        }

        // Write end of contents
        os << token::END_LIST;
    }
    else
    {
        // Write size of HashList
        os << nl << L.size();

        // Write beginning of contents
        os << nl << token::BEGIN_LIST << nl;

        // Write HashList contents
        for
        (
            typename HashList<T>::iterator iter = L.begin();
            iter != L.end();
            iter++
        )
        {
            os << iter.index() << token::SPACE << iter() << nl;
        }

        // Write end of contents
        os << token::END_LIST;
    }

    // Check state of IOstream
    os.check("Ostream& operator<<(Ostream&, const HashList&)");

    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
