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
    ),
    counter_(0),
    matchTol_(1e-6),
    twoDMesh_(false)
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

    // Set mesh characteristics
    twoDMesh_ =
    (
        (meshFrom.nGeometricD() == 2 && meshTo.nGeometricD() == 2)
      ? true : false
    );

    // Check if the source mesh has a calculated addressing
    // If yes, try and invert that.
    IOobject srcHeader
    (
        "addressing",
        meshFrom.time().timeName(),
        meshFrom,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if (srcHeader.headerOk() && !addressing_.headerOk())
    {
        Info << " Found addressing in source directory."
             << " Checking for compatibility." << endl;

        if (invertAddressing())
        {
            Info << " Inversion successful. " << endl;
            return;
        }
    }

    // Track calculation time
    clockTime calcTimer;

    if (nThreads == 1)
    {
        calcAddressingAndWeights(0, meshTo.nCells(), true);
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
        fromMesh().faceEdges();
        fromMesh().cellCells();

        toMesh().edges();
        toMesh().cells();
        toMesh().cellEdges();
        toMesh().cellPoints();
        toMesh().faceEdges();

        multiThreader threader(nThreads);

        // Set one handler per thread
        PtrList<handler> hdl(threader.getNumThreads());

        forAll(hdl, i)
        {
            hdl.set
            (
                i,
                new handler(*this, threader)
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
            hdl[i].set(0, &tStarts[i]);
            hdl[i].set(1, &tSizes[i]);

            // Lock the slave thread first
            hdl[i].lock(handler::START);
            hdl[i].unsetPredicate(handler::START);

            hdl[i].lock(handler::STOP);
            hdl[i].unsetPredicate(handler::STOP);
        }

        // Submit jobs to the work queue
        forAll(hdl, i)
        {
            threader.addToWorkQueue
            (
                &calcAddressingAndWeightsThreaded,
                &(hdl[i])
            );

            // Wait for a signal from this thread
            // before moving on.
            hdl[i].waitForSignal(handler::START);
        }

        // Track progress of threads
        while (true)
        {
            sleep(3);

            ctrMutex_.lock();

            Info << "\r Progress: "
                 << 100.0 * (double(counter_) / toMesh().nCells())
                 << "% :  Cells processed: " << counter_
                 << " out of " << toMesh().nCells() << " total."
                 << "             "
                 << flush;

            ctrMutex_.unlock();

            if (counter_ == toMesh().nCells())
            {
                break;
            }
        }

        // Synchronize all threads
        forAll(hdl, i)
        {
            hdl[i].waitForSignal(handler::STOP);
        }
    }

    Info << nl << " Calculation time: " << calcTimer.elapsedTime() << endl;

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
    handler *thread = static_cast<handler*>(argument);

    if (thread->slave())
    {
        thread->sendSignal(handler::START);
    }

    conservativeMeshToMesh& interpolator = thread->reference();

    // Recast the pointers for the argument
    label& cellStart = *(static_cast<label*>(thread->operator()(0)));
    label& cellSize = *(static_cast<label*>(thread->operator()(1)));

    // Now calculate addressing
    interpolator.calcAddressingAndWeights(cellStart, cellSize);

    if (thread->slave())
    {
        thread->sendSignal(handler::STOP);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
