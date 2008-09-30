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

Author
    Sandeep Menon

SourceFiles
    multiThreader.C

\*---------------------------------------------------------------------------*/

#ifndef multiThreader_H
#define multiThreader_H

#include "IOdictionary.H"

#include <pthread.h>
#include <sys/types.h>
#include <unistd.h>
extern "C" 
{ 
    typedef void *(*externThreadFunctionType)(void *); 
}
typedef void *(*threadFunctionType)(void *);

#define MT_MAX_THREADS 32
#define threadReturnType void *
#define threadReturnValue NULL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
/*---------------------------------------------------------------------------*\
                           Class multiThreader Declaration
\*---------------------------------------------------------------------------*/

class multiThreader
{
public:    
    
        // Common structure that holds function / data information
        // for all instantiated threads. The methodArg pointer contains
        // any arguments that the multithreaded class method may require.
        struct threadInfo
        {
            label   ID;
            label   numThreads;
            void    *methodArg;
        };     
    
    // Selectors

        //- Construct from specified dictionary
        static autoPtr<multiThreader> New(const dictionary&);


    // Destructor

        ~multiThreader();
        

    // Member Functions

        //- Set the method for multithreaded execution
        void setMethod(threadFunctionType);
        
        //- Set the method argument for a specific index
        void setData(label, void *);
        
        //- Execute threads
        void executeThreads(); 
        
        //- Return the number of threads
        label getNumThreads();
        
private:
    
    // Private data
        label                   numThreads_;
        threadFunctionType      threadedMethod_;
        threadInfo              infoList_[MT_MAX_THREADS];
        void                    *dataList_[MT_MAX_THREADS];

    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiThreader(const multiThreader&);

        //- Disallow default bitwise assignment
        void operator=(const multiThreader&);

    // Constructors

        multiThreader(const dictionary&);    
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "multiThreader.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
