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
    maxQueueSize_(-1),
    poolInfo_(NULL)
{        
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
    
    // Obtain maxQueueSize from the specified dictionary.
    if (threadDict.found("maxQueueSize"))
    {
        maxQueueSize_ = readLabel(threadDict.lookup("maxQueueSize"));
    }
    else
    {
        maxQueueSize_ = 10;
    }
    
    // Initialize the thread pool
    initializeThreadPool();
}

Foam::multiThreader::Mutex::Mutex()
{
    if (pthread_mutex_init(&lock_, NULL))
    {
        FatalErrorIn("multiThreader::Mutex::Mutex()")
            << "Unable to initialize mutex"
            << abort(FatalError);        
    }    
}

Foam::multiThreader::Conditional::Conditional()
{
    if (pthread_cond_init(&condition_, NULL))
    {
        FatalErrorIn("multiThreader::Conditional::Conditional()")
            << "Unable to initialize condition"
            << abort(FatalError);        
    }    
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

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::multiThreader::~multiThreader()
{
    destroyThreadPool();
}

Foam::multiThreader::Mutex::~Mutex()
{
    if (pthread_mutex_destroy(&lock_))
    {
        FatalErrorIn("multiThreader::Mutex::~Mutex()")
            << "Unable to destroy mutex"
            << abort(FatalError);        
    }    
}

Foam::multiThreader::Conditional::~Conditional()
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
    poolInfo_->threads = new pthread_t[numThreads_];
    poolInfo_->head = NULL;
    poolInfo_->tail = NULL;
    
    // Initialize flags
    poolInfo_->queueClosed = false;
    poolInfo_->shutDown = false;
    
    // Create worker threads and have them wait for jobs
    for (label tIndex = 0; tIndex < numThreads_; tIndex++)
    {        
        int status = pthread_create
                     ( 
                         &(poolInfo_->threads[tIndex]), 
                         NULL,
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
            poolInfo->threader->signal(poolInfo->queueEmpty);
        }
        
        // Unlock the work queue
        poolInfo->queueLock.unlock();
        
        // Perform the work
        myWorkItem->function(myWorkItem->arg);
        
        // Free up the work item
        delete myWorkItem;
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
    for(label i=0; i < numThreads_; i++) 
    {
        if (pthread_join(poolInfo_->threads[i],NULL))
        {
            FatalErrorIn("multiThreader::destroyThreadPool()")
                << "pthread_join failed."
                << abort(FatalError);            
        }
    }    
    
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
label Foam::multiThreader::getNumThreads()
{
    return numThreads_;
}

//- Return the maxQueueSize
label Foam::multiThreader::getMaxQueueSize()
{
    return maxQueueSize_;
}

void Foam::multiThreader::Mutex::lock()
{
    if (pthread_mutex_lock(&lock_))
    {
        FatalErrorIn("multiThreader::destroyThreadPool()")
            << "Unable to lock the work queue."
            << abort(FatalError);
    }     
}

void Foam::multiThreader::Mutex::unlock()
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
