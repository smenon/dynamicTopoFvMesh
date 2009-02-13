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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::multiThreader::debug = false;
bool Foam::Mutex::debug = false;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiThreader::multiThreader(int numThreads)
:
    maxQueueSize_(10),
    poolInfo_(NULL)
{       
    if (numThreads > 0)
    {
        numThreads_ = numThreads;
        Info << "Initializing threading environment with " 
             << numThreads_ << " threads." << endl;        
    }
    else
    {
        // Default number of threads at one (single-threaded)
        numThreads_ = 1;
        Info << "Defaulting threading environment to one thread." << endl;         
    }
    
    // Initialize the thread pool
    initializeThreadPool();
}

Foam::Mutex::Mutex()
{
    // Set attributes based on debug flag
    pthread_mutexattr_t attribute;
    pthread_mutexattr_init(&attribute);

    if (debug)
    {
        pthread_mutexattr_settype(&attribute, PTHREAD_MUTEX_ERRORCHECK);
    }
    else
    {
        pthread_mutexattr_settype(&attribute, PTHREAD_MUTEX_NORMAL);
    }

    if (pthread_mutex_init(&lock_, &attribute))
    {
        FatalErrorIn("multiThreader::Mutex::Mutex()")
            << "Unable to initialize mutex"
            << abort(FatalError);        
    }

    // Destroy the attribute
    pthread_mutexattr_destroy(&attribute);
}

Foam::Conditional::Conditional()
{
    if (pthread_cond_init(&condition_, NULL))
    {
        FatalErrorIn("multiThreader::Conditional::Conditional()")
            << "Unable to initialize condition"
            << abort(FatalError);        
    }    
}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::multiThreader::~multiThreader()
{
    destroyThreadPool();
}

Foam::Mutex::~Mutex()
{
    if (pthread_mutex_destroy(&lock_))
    {
        FatalErrorIn("multiThreader::Mutex::~Mutex()")
            << "Unable to destroy mutex"
            << abort(FatalError);        
    }    
}

Foam::Conditional::~Conditional()
{
    if (pthread_cond_destroy(&condition_))
    {
        FatalErrorIn("multiThreader::Conditional::Conditional()")
            << "Unable to destroy condition"
            << abort(FatalError);        
    }     
}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::multiThreader::initializeThreadPool()
{
    // Allocate the threadPool structure
    poolInfo_ = new threadPool;
    
    // Initialize fields
    poolInfo_->threader = this;
    poolInfo_->numThreads = numThreads_;
    poolInfo_->queueSize = 0;
    poolInfo_->busyThreads = 0;
    poolInfo_->threads = new pthread_t[numThreads_];
    poolInfo_->head = NULL;
    poolInfo_->tail = NULL;
    
    // Initialize flags
    poolInfo_->queueClosed = false;
    poolInfo_->shutDown = false;
    
    // Initialize thread attributes
    pthread_attr_init(&(poolInfo_->attr));
    pthread_attr_setdetachstate(&(poolInfo_->attr), PTHREAD_CREATE_JOINABLE);
    
    // Create worker threads and have them wait for jobs
    for (int tIndex = 0; tIndex < numThreads_; tIndex++)
    {        
        int status = pthread_create
                     ( 
                         &(poolInfo_->threads[tIndex]), 
                         &(poolInfo_->attr),
                         reinterpret_cast<externThreadFunctionType>
                         (
                             poolThread
                         ),
                         reinterpret_cast<void *>
                         (
                             poolInfo_
                         ) 
                     );
        
        if (status != 0)
        {
            FatalErrorIn("multiThreader::initThreadPool()")
                << "pthread_create could not initialize thread: "
                << tIndex
                << abort(FatalError);
        }        
    }   
}

threadReturnType Foam::multiThreader::poolThread(void *arg)
{
    // Typecast the argument into the required structure
    threadPool *poolInfo = reinterpret_cast<threadPool *>(arg);
    
    // Work queue loop
    while (true)
    {
        // Lock the work queue
        poolInfo->queueLock.lock();
        
        // Wait for work to arrive in the queue
        while ((poolInfo->queueSize == 0) && (!poolInfo->shutDown)) 
        {
#           ifdef FULLDEBUG            
            if (debug)
            {
                Info << "poolThread::Wait on queueNotEmpty." << endl;
            }
#           endif             
            poolInfo->threader->waitForCondition
                                (
                                    poolInfo->queueNotEmpty,
                                    poolInfo->queueLock
                                );
        }  
        
        // Check for shutdown
        if (poolInfo->shutDown) 
        {
            poolInfo->queueLock.unlock();
            pthread_exit(NULL);
        }  
        
        // Pick an item off the queue, and get to work
        workQueueItem *myWorkItem = poolInfo->head;
        poolInfo->queueSize--;
        if (poolInfo->queueSize == 0)
        {
            poolInfo->head = poolInfo->tail = NULL;
        }
        else
        {
            poolInfo->head = myWorkItem->next;
        }
        
        // Handle a waiting destructor
        if (poolInfo->queueSize == 0)
        {
#           ifdef FULLDEBUG            
            if (debug)
            {
                Info << "poolThread::Signaling: Empty queue." << endl;
            }
#           endif             
            poolInfo->threader->signal(poolInfo->queueEmpty);
        }
        
        // Increment the busy queue
        poolInfo->busyThreads++;
        
        // Unlock the work queue
        poolInfo->queueLock.unlock();
        
        // Perform the work
        myWorkItem->function(myWorkItem->arg);
        
        // Free up the work item
        delete myWorkItem;
        
        // Lock the work queue
        poolInfo->queueLock.lock();
        
        // Finished the allotted work, decrement the busy queue
        poolInfo->busyThreads--;
        
        // Signal any conditions waiting on the busy queue
        if ((poolInfo->busyThreads == 0) && (poolInfo->queueSize == 0))
        {
#           ifdef FULLDEBUG            
            if (debug)
            {
                Info << "Signaling: No busy threads." << endl;
            }
#           endif            
            poolInfo->threader->signal(poolInfo->noBusyThreads);
        }     
        
        // Unlock the work queue
        poolInfo->queueLock.unlock();
    }
    
    return threadReturnValue;
}

void Foam::multiThreader::addToWorkQueue
(
    void (*tFunction)(void*), 
    void *arg
)
{
    // Lock the work queue
    poolInfo_->queueLock.lock();
    
    // If occupied, wait for the queue to free-up
    while
    (
         (poolInfo_->queueSize == maxQueueSize_) 
      && (!(poolInfo_->shutDown || poolInfo_->queueClosed))  
    ) 
    {
#       ifdef FULLDEBUG            
        if (debug)
        {
            Info << "addToWorkQueue:: Wait on queueNotFull." << endl;
        }
#       endif          
        waitForCondition(poolInfo_->queueNotFull, poolInfo_->queueLock);
    } 
    
    // Is the pool in the process of being destroyed?
    // Unlock the mutex and return to caller.
    if (poolInfo_->shutDown || poolInfo_->queueClosed) 
    {
        poolInfo_->queueLock.unlock();
        return;
    }  
    
    // Allocate a new work structure
    workQueueItem *newWorkItem = new workQueueItem;
    newWorkItem->function = tFunction;
    newWorkItem->arg = arg;
    newWorkItem->next = NULL;    
    
    // Add new work structure to the queue
    if (poolInfo_->queueSize == 0) 
    {
        poolInfo_->tail = poolInfo_->head = newWorkItem;
        broadCast(poolInfo_->queueNotEmpty);
    } 
    else 
    {
        poolInfo_->tail->next = newWorkItem;
        poolInfo_->tail = newWorkItem;
    }

    poolInfo_->queueSize++;    
    
    // Unlock the work queue
    poolInfo_->queueLock.unlock();    
}

//- Wait for all worker threads to complete
void Foam::multiThreader::waitForCompletion()
{
    // Lock the work queue
    poolInfo_->queueLock.lock();    
    
    // Wait for all threads to finish work
    waitForCondition(poolInfo_->noBusyThreads, poolInfo_->queueLock);
    
    // Unlock the work queue
    poolInfo_->queueLock.unlock();    
}

void Foam::multiThreader::destroyThreadPool()
{
    // Lock the work queue
    poolInfo_->queueLock.lock();
    
    // Is a shutdown already in progress?
    if (poolInfo_->queueClosed || poolInfo_->shutDown) 
    {
        // Unlock the mutex and return
        poolInfo_->queueLock.unlock();
        return;
    }    
    
    poolInfo_->queueClosed = true;

    // Wait for workers to drain the queue
    while (poolInfo_->queueSize != 0) 
    {
        waitForCondition(poolInfo_->queueEmpty, poolInfo_->queueLock);
    }

    poolInfo_->shutDown = true;        
    
    // Unlock the work queue
    poolInfo_->queueLock.unlock();
    
    // Wake up workers so that they check the shutdown flag
    broadCast(poolInfo_->queueNotEmpty);
    broadCast(poolInfo_->queueNotFull);

    // Wait for all workers to exit
    for(int i=0; i < numThreads_; i++) 
    {
        if (pthread_join(poolInfo_->threads[i],NULL))
        {
            FatalErrorIn("multiThreader::destroyThreadPool()")
                << "pthread_join failed."
                << abort(FatalError);            
        }
    }    
    
    // Destroy the attribute
    pthread_attr_destroy(&(poolInfo_->attr));
    
    // Deallocate the work-queue and pool structure
    delete [] poolInfo_->threads;
    
    workQueueItem *currentNode;
    while(poolInfo_->head != NULL) 
    {
        currentNode = poolInfo_->head->next; 
        poolInfo_->head = poolInfo_->head->next;
        delete currentNode;
    }
    
    delete poolInfo_;     
}

void Foam::multiThreader::waitForCondition
(
    Conditional& condition, 
    Mutex& mutex
)
{
    if (pthread_cond_wait(condition(),mutex()))
    {
        FatalErrorIn("multiThreader::waitForCondition(..)")
            << "Conditional wait failed."
            << abort(FatalError);            
    }    
}

void Foam::multiThreader::broadCast(Conditional& condition)
{
    if (pthread_cond_broadcast(condition()))
    {
        FatalErrorIn("multiThreader::broadCast()")
            << "Unable to broadcast."
            << abort(FatalError);                
    }    
}

void Foam::multiThreader::signal(Conditional& condition)
{
    if (pthread_cond_signal(condition()))
    {
        FatalErrorIn("multiThreader::signal()")
            << "Unable to signal."
            << abort(FatalError);                
    }    
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return the number of threads
int Foam::multiThreader::getNumThreads()
{
    return numThreads_;
}

//- Obtain the thread ID for a given index
const pthread_t Foam::multiThreader::getID(int index)
{
    if (poolInfo_ && index > -1 && index < numThreads_)
    {
        return poolInfo_->threads[index];
    }
    else
    {
        FatalErrorIn("multiThreader::getID(int index)")
            << "Invalid request for ID."
            << abort(FatalError);
    }

    // This should never happen anyway.
    return poolInfo_->threads[index];
}

//- Return true if the number of threads is more than one.
bool Foam::multiThreader::multiThreaded() const
{
    return (numThreads_ > 1);
}

//- Return the maxQueueSize
int Foam::multiThreader::getMaxQueueSize()
{
    return maxQueueSize_;
}

//- Set the maxQueueSize
void Foam::multiThreader::setMaxQueueSize(int size)
{
    if (size > 0)
    {
        maxQueueSize_ = size;
    }
    else
    {
        FatalErrorIn("multiThreader::setMaxQueueSize(int size)")
            << "Improper value for MaxQueueSize."
            << abort(FatalError);        
    }
}

void Foam::Mutex::lock()
{
    if (pthread_mutex_lock(&lock_))
    {
        FatalErrorIn("multiThreader::destroyThreadPool()")
            << "Unable to lock the work queue."
            << abort(FatalError);
    }     
}

bool Foam::Mutex::tryLock()
{
    label retVal;

    if ((retVal = pthread_mutex_trylock(&lock_)) != 0)
    {
#       ifdef FULLDEBUG
        if (debug)
        {
            if (retVal == EINVAL)
            {
                FatalErrorIn("multiThreader::Mutex::trylock()")
                    << "Mutex returned EINVAL."
                    << abort(FatalError);
            }
            if (retVal == EFAULT)
            {
                FatalErrorIn("multiThreader::Mutex::trylock()")
                    << "Mutex returned EFAULT."
                    << abort(FatalError);
            }
        }
#       endif
    }

    return retVal;
}

void Foam::Mutex::unlock()
{
    if (pthread_mutex_unlock(&lock_))
    {
        FatalErrorIn("multiThreader::Mutex::unlock()")
            << "Unable to unlock the mutex."
            << abort(FatalError);
    }    
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
