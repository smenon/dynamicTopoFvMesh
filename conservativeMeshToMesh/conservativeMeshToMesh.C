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

\*---------------------------------------------------------------------------*/

#include "clockTime.H"
#include "multiThreader.H"
#include "threadHandler.H"
#include "conservativeMeshToMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(conservativeMeshToMesh, 0);

// Explicitly specify typeNames for IOFields
defineCompoundTypeName(Field<scalar>, scalarField);
addCompoundToRunTimeSelectionTable(Field<scalar>, scalarField);

defineCompoundTypeName(Field<vector>, vectorField);
addCompoundToRunTimeSelectionTable(Field<vector>, vectorField);

defineTemplateTypeNameAndDebugWithName
(
    IOList<scalarField>, "scalarFieldList", 0
);

defineTemplateTypeNameAndDebugWithName
(
    IOList<vectorField>, "vectorFieldList", 0
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

conservativeMeshToMesh::conservativeMeshToMesh
(
    const fvMesh& meshFrom,
    const fvMesh& meshTo,
    const label nThreads,
    const bool forceRecalculation,
    const bool writeAddressing
)
:
    meshToMesh(meshFrom, meshTo),
    addressing_
    (
        IOobject
        (
            "addressing",
            meshTo.time().timeName(),
            meshTo,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        meshTo.nCells()
    ),
    weights_
    (
        IOobject
        (
            "weights",
            meshTo.time().timeName(),
            meshTo,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        meshTo.nCells()
    ),
    centres_
    (
        IOobject
        (
            "centres",
            meshTo.time().timeName(),
            meshTo,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        meshTo.nCells()
    )
{
    if (addressing_.headerOk() && weights_.headerOk() && centres_.headerOk())
    {
        // Check if sizes match. Otherwise, re-calculate.
        if
        (
            addressing_.size() == meshTo.nCells() &&
            weights_.size() == meshTo.nCells() &&
            centres_.size() == meshTo.nCells() &&
           !forceRecalculation
        )
        {
            Info << " Reading addressing from file." << endl;

            return;
        }

        Info << " Recalculating addressing." << endl;
    }
    else
    {
        Info << " Calculating addressing." << endl;
    }

    // Track calculation time
    clockTime calcTimer;

    if (nThreads == 1)
    {
        calcAddressingAndWeights(0, meshTo.nCells());
    }
    else
    {
        // Prior to multi-threaded operation,
        // force calculation of demand-driven data
        fromMesh().edges();
        fromMesh().cells();
        fromMesh().cellCentres();
        fromMesh().cellEdges();
        fromMesh().cellPoints();
        fromMesh().cellCells();

        toMesh().edges();
        toMesh().cells();
        toMesh().cellEdges();
        toMesh().cellPoints();

        multiThreader threader(nThreads);

        // Set one handler per thread
        PtrList<threadHandler<conservativeMeshToMesh> > hdl
        (
            threader.getNumThreads()
        );

        forAll(hdl, i)
        {
            hdl.set
            (
                i,
                new threadHandler<conservativeMeshToMesh>
                (
                    *this,
                    threader
                )
            );
        }

        // Simple, but inefficient load-balancing scheme
        labelList tStarts(threader.getNumThreads(), 0);
        labelList tSizes(threader.getNumThreads(), 0);

        label index = meshTo.nCells(), j = 0;

        while (index--)
        {
            tSizes[(j = tSizes.fcIndex(j))]++;
        }

        label total = 0;

        for (label i = 1; i < tStarts.size(); i++)
        {
            tStarts[i] = tSizes[i-1] + total;

            total += tSizes[i-1];
        }

        if (debug)
        {
            Info << " Load starts: " << tStarts << endl;
            Info << " Load sizes: " << tSizes << endl;
        }

        // Set the argument list for each thread
        forAll(hdl, i)
        {
            // Size up the argument list
            hdl[i].setSize(2);

            // Set the start/end cell indices
            hdl[i].set(0, reinterpret_cast<void *>(&tStarts[i]));
            hdl[i].set(1, reinterpret_cast<void *>(&tSizes[i]));

            // Lock the slave thread first
            hdl[i].lock(threadHandler<conservativeMeshToMesh>::START);
            hdl[i].unsetPredicate(threadHandler<conservativeMeshToMesh>::START);

            hdl[i].lock(threadHandler<conservativeMeshToMesh>::STOP);
            hdl[i].unsetPredicate(threadHandler<conservativeMeshToMesh>::STOP);
        }

        // Submit jobs to the work queue
        forAll(hdl, i)
        {
            threader.addToWorkQueue
            (
                &calcAddressingAndWeightsThreaded,
                reinterpret_cast<void *>(&(hdl[i]))
            );

            // Wait for a signal from this thread
            // before moving on.
            hdl[i].waitForSignal
            (
                threadHandler<conservativeMeshToMesh>::START
            );
        }

        // Synchronize all threads
        forAll(hdl, i)
        {
            hdl[i].waitForSignal(threadHandler<conservativeMeshToMesh>::STOP);
        }
    }

    Info << " Calculation time: " << calcTimer.elapsedTime() << endl;

    if (writeAddressing)
    {
        Info << " Writing addressing to disk." << endl;

        addressing_.write();
        weights_.write();
        centres_.write();
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

conservativeMeshToMesh::~conservativeMeshToMesh()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void conservativeMeshToMesh::calcAddressingAndWeightsThreaded
(
    void *argument
)
{
    // Recast the argument
    threadHandler<conservativeMeshToMesh> *thread =
    (
        reinterpret_cast<threadHandler<conservativeMeshToMesh>*>(argument)
    );

    if (thread->slave())
    {
        thread->sendSignal(threadHandler<conservativeMeshToMesh>::START);
    }

    conservativeMeshToMesh& interpolator = thread->reference();

    // Recast the pointers for the argument
    label& cellStart = *(reinterpret_cast<label*>(thread->operator()(0)));
    label& cellSize = *(reinterpret_cast<label*>(thread->operator()(1)));

    // Now calculate addressing
    interpolator.calcAddressingAndWeights(cellStart, cellSize);

    if (thread->slave())
    {
        thread->sendSignal(threadHandler<conservativeMeshToMesh>::STOP);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
