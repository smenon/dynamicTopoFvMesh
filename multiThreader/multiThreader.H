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
    Class which provides support for generic multithreaded execution of class
    methods using POSIX threads.

    Thread-pooling is implemented as described in:
   "Using POSIX Threads: Programming with Pthreads"
    by Brad nichols, Dick Buttlar, Jackie Farrell
      O'Reilly & Associates, Inc.
 
Author
    Sandeep Menon

SourceFiles
    multiThreader.C

\*---------------------------------------------------------------------------*/

#ifndef multiThreader_H
#define multiThreader_H

#include "error.H"

#include <pthread.h>
#include <errno.h>
#include <sys/types.h>
#include <unistd.h>
extern "C" 
{ 
    typedef void *(*externThreadFunctionType)(void *); 
}
typedef void *(*threadFunctionType)(void *);

#define threadReturnType void *
#define threadReturnValue NULL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Mutex
{
    // Private data
        pthread_mutex_t lock_;

public:

    // Debug switch
    static bool debug;

    // Constructor
    Mutex();

    // Destructor
    ~Mutex();

    // Lock the mutex
    void lock();

    // Try to acquire the mutex.
    // Returns zero if successful
    bool tryLock();

    // Unlock the mutex
    void unlock();

    // Operator overload
    inline pthread_mutex_t* operator()()
    {
        return &lock_;
    }
};

class rwMutex
{
    // Private data
        pthread_rwlock_t lock_;

public:

    // Enumerants for lock-types
    enum lockType
    {
        READ_LOCK,
        WRITE_LOCK
    };

    // Debug switch
    static bool debug;

    // Constructor
    rwMutex();

    // Destructor
    ~rwMutex();

    // Lock the mutex
    void lock(const lockType lType);

    // Try to acquire a lock
    bool tryLock(const lockType lType);

    // Unlock the mutex
    void unlock();

    // Operator overload
    inline pthread_rwlock_t* operator()()
    {
        return &lock_;
    }
};

class Conditional
{
    // Private data
        pthread_cond_t condition_;

public:

    // Constructor
    Conditional();

    // Destructor
    ~Conditional();

    // Operator overload
    inline pthread_cond_t* operator()()
    {
        return &condition_;
    }        
};

/*---------------------------------------------------------------------------*\
                           Class multiThreader Declaration
\*---------------------------------------------------------------------------*/

class multiThreader
{            
    // Private data
    
        int numThreads_;
        
        int maxQueueSize_;

        // Work-queue item: Holds a pointer to the method and its argument
        struct workQueueItem
        {
            void (*function)(void*);
            void *arg;
            workQueueItem *next;
        };            
        
        // Common structure for all threads in the pool.
        struct threadPool
        {
            multiThreader   *threader;
            int             numThreads;
            int             queueSize;
            int             busyThreads;
            pthread_t       *threads;
            pthread_attr_t  attr;
            workQueueItem   *head, *tail;
            
            // Flags
            bool queueClosed;
            bool shutDown;
            
            // Synchronization
            Mutex       queueLock;
            Conditional queueNotEmpty;
            Conditional queueNotFull;
            Conditional queueEmpty;
            Conditional noBusyThreads;
        } *poolInfo_;           

    // Private Member Functions
        
        void initializeThreadPool();
        
        static threadReturnType poolThread(void *arg);
        
        void destroyThreadPool();
        
        void waitForCondition(Conditional&, Mutex&);
        
        void broadCast(Conditional&);
        
        void signal(Conditional&);

        //- Disallow default bitwise copy construct
        multiThreader(const multiThreader&);

        //- Disallow default bitwise assignment
        void operator=(const multiThreader&);            
    
public:      
    
    // Debug switch
    
        static bool debug;
    
    // Constructor

        //- Construct with specified number of threads
        multiThreader(int numThreads);


    // Destructor

        ~multiThreader();
        

    // Member Functions
        
        //- Return the number of threads
        int getNumThreads();

        //- Obtain the thread ID for a given index
        pthread_t getID(int index);

        //- Return true if the number of threads is more than one.
        bool multiThreaded() const;
        
        //- Return the maxQueueSize
        int getMaxQueueSize();
        
        //- Set the maxQueueSize
        void setMaxQueueSize(int size);
        
        //- Add a function to the work queue
        void addToWorkQueue(void (*tFunction)(void*), void *arg);
        
        //- Wait for all worker threads to complete
        void waitForCompletion();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
