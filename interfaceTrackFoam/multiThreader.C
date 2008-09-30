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
    multiThreader

Description
    Implementation of the multiThreader class

Author
    Sandeep Menon

\*----------------------------------------------------------------------------*/

#include "multiThreader.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiThreader::multiThreader(const dictionary& threadDict)
:
    numThreads_(-1),
    threadedMethod_(NULL)
{
    // Initialize members
    for (label index=0; index < MT_MAX_THREADS; index++)
    {
        this->dataList_[index] = NULL;
        this->infoList_[index].ID = index;
    }
        
    // Obtain numThreads from the specified dictionary.
    if (threadDict.found("threads"))
    {
        numThreads_ = readLabel(threadDict.lookup("threads"));
        Info << "Initializing threading environment with " 
             << numThreads_ << " threads." << endl;
    }
    else
    {
        // Default number of threads at one (single-threaded)
        numThreads_ = 1;
        Info << "Defaulting threading environment to one thread." << endl;        
    }
    
    // Limit the number of threads to the global maximum
    numThreads_ = (numThreads_ < MT_MAX_THREADS) ? numThreads_ : MT_MAX_THREADS;
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::multiThreader> 
Foam::multiThreader::New
(
    const dictionary& threadDict
)
{
    return autoPtr<multiThreader>(new multiThreader(threadDict));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiThreader::~multiThreader()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Set the method for multithreaded execution
void Foam::multiThreader::setMethod
(
    threadFunctionType function 
)
{
    this->threadedMethod_ = function;    
}

//- Set the method argument for a specific index
void Foam::multiThreader::setData
(
    label index, 
    void *argument
)
{
    if (index >= numThreads_)
    {
        FatalErrorIn("multiThreader::setData(label index, void *argument)")
            << "Index out of range: 0 to "
            << (numThreads_ - 1)
            << abort(FatalError);
        return;        
    }
    else
    {
        this->dataList_[index] = argument;
    }
}

//- Execute threads
void Foam::multiThreader::executeThreads()
{
    if (!this->threadedMethod_)
    {
        FatalErrorIn("multiThreader::executeThreads()")
            << "Method was not set."
            << abort(FatalError);
        return;
    }
    
    // Initialize the pthread attributes
    pthread_attr_t attr;
    pthread_attr_init(&attr);    
    pthread_t pid[MT_MAX_THREADS];
    
    for (label tIndex = 1; tIndex < numThreads_; tIndex++)
    {
        // Set the appropriate argument structures
        this->infoList_[tIndex].methodArg = this->dataList_[tIndex];
        this->infoList_[tIndex].numThreads = this->numThreads_;
        
        int status = pthread_create
                     ( 
                        &(pid[tIndex]), 
                        &attr,
                        reinterpret_cast<externThreadFunctionType>
                        (
                            this->threadedMethod_
                        ),
                        reinterpret_cast<void *>
                        (
                            &this->infoList_[tIndex]
                        ) 
                     );
        
        if (status != 0)
        {
            FatalErrorIn("multiThreader::executeThreads()")
                << "pthread_create could not initialize thread: "
                << tIndex
                << abort(FatalError);
        }        
    }
    
    // Primary thread also calls the method (default single-threaded behaviour)
    this->infoList_[0].methodArg = this->dataList_[0];
    this->infoList_[0].numThreads = this->numThreads_;  
    this->threadedMethod_(reinterpret_cast<void *>(&this->infoList_[0]));
    
    // Wait for all threads to exit
    for (label tIndex = 1; tIndex < numThreads_; tIndex++)
    {
        pthread_join(pid[tIndex], NULL);
    }    
}

//- Return the number of threads
label Foam::multiThreader::getNumThreads()
{
    return numThreads_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::multiThreader::operator=(const multiThreader& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("multiThreader::operator=(const multiThreader&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

// ************************************************************************* //
