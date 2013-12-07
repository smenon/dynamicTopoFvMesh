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

#include "octree.H"
#include "clockTime.H"
#include "multiThreader.H"
#include "threadHandler.H"
#include "octreeDataFace.H"
#include "octreeDataCell.H"
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculate nearest cell addressing
void conservativeMeshToMesh::calcCellAddressing()
{
    if (debug)
    {
        Info<< "conservativeMeshToMesh::calcCellAddressing() : "
            << "calculating nearest cell addressing" << endl;
    }

    // Fetch references to source mesh data
    const cellList& srcCells = srcMesh().cells();
    const pointField& srcPoints = srcMesh().points();

    treeBoundBox meshBb(srcPoints);

    scalar typDim =
    (
        meshBb.avgDim()/(2.0*cbrt(scalar(srcCells.size())))
    );

    treeBoundBox shiftedBb
    (
        meshBb.min(),
        meshBb.max() + vector(typDim, typDim, typDim)
    );

    // Wrap data for octree into container
    octreeDataCell shapes(srcMesh());

    octree<octreeDataCell> oc
    (
        shiftedBb,  // overall bounding box
        shapes,     // all information needed to do checks on cells
        1,          // min. levels
        10.0,       // max. size of leaves
        10.0        // maximum ratio of cubes vs. cells
    );

    if (debug)
    {
        oc.printStats(Info);
    }

    // Loop through target mesh cells and
    // find the nearest on the source mesh
    const vectorField& tgtCentres = tgtMesh().cellCentres();

    scalar tightestDist;
    treeBoundBox tightest;

    forAll(cellAddressing_, cellI)
    {
        tightest = meshBb;                 // starting search bb
        tightestDist = Foam::GREAT;        // starting max distance

        cellAddressing_[cellI] =
        (
            oc.findNearest
            (
                tgtCentres[cellI],
                tightest,
                tightestDist
            )
        );
    }
}


// Output an entity as a VTK file
void conservativeMeshToMesh::writeVTK
(
    const word& name,
    const label entity,
    const label primitiveType,
    const bool useOldConnectivity
) const
{
    writeVTK
    (
        name,
        labelList(1, entity),
        primitiveType,
        useOldConnectivity
    );
}


// Output a list of primitives as a VTK file.
//  - primitiveType is:
//      0: List of points
//      1: List of edges
//      2: List of faces
//      3: List of cells
void conservativeMeshToMesh::writeVTK
(
    const word& name,
    const labelList& cList,
    const label primitiveType,
    const bool useOldConnectivity,
    const UList<scalar>& field
) const
{
    if (useOldConnectivity)
    {
        // Use points from polyMesh
        meshOps::writeVTK
        (
            tgtMesh(),
            name,
            cList,
            primitiveType,
            srcMesh().points(),
            srcMesh().edges(),
            srcMesh().faces(),
            srcMesh().cells(),
            srcMesh().faceOwner(),
            field
        );
    }
    else
    {
        meshOps::writeVTK
        (
            tgtMesh(),
            name,
            cList,
            primitiveType,
            tgtMesh().points(),
            tgtMesh().edges(),
            tgtMesh().faces(),
            tgtMesh().cells(),
            tgtMesh().faceOwner(),
            field
        );
    }
}


// Multi-threaded weighting factor computation
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

conservativeMeshToMesh::conservativeMeshToMesh
(
    const fvMesh& srcMesh,
    const fvMesh& tgtMesh,
    const label nThreads,
    const bool forceRecalculation,
    const bool writeAddressing
)
:
    meshSrc_(srcMesh),
    meshTgt_(tgtMesh),
    addressing_
    (
        IOobject
        (
            "addressing",
            tgtMesh.time().timeName(),
            tgtMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        tgtMesh.nCells()
    ),
    weights_
    (
        IOobject
        (
            "weights",
            tgtMesh.time().timeName(),
            tgtMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tgtMesh.nCells()
    ),
    centres_
    (
        IOobject
        (
            "centres",
            tgtMesh.time().timeName(),
            tgtMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        tgtMesh.nCells()
    ),
    cellAddressing_(tgtMesh.nCells()),
    counter_(0),
    boundaryAddressing_(tgtMesh.boundaryMesh().size())
{
    if (addressing_.headerOk() && weights_.headerOk() && centres_.headerOk())
    {
        // Check if sizes match. Otherwise, re-calculate.
        if
        (
            addressing_.size() == tgtMesh.nCells() &&
            weights_.size() == tgtMesh.nCells() &&
            centres_.size() == tgtMesh.nCells() &&
           !forceRecalculation
        )
        {
            Info<< " Reading addressing from file." << endl;

            return;
        }

        Info<< " Recalculating addressing." << endl;
    }
    else
    {
        Info<< " Calculating addressing." << endl;
    }

    // Check if the source mesh has a calculated addressing
    // If yes, try and invert that.
    IOobject srcHeader
    (
        "addressing",
        srcMesh.time().timeName(),
        srcMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if (srcHeader.headerOk() && !addressing_.headerOk())
    {
        Info<< " Found addressing in source directory."
            << " Checking for compatibility." << endl;

        if (invertAddressing())
        {
            Info<< " Inversion successful. " << endl;
            return;
        }
    }

    // Track calculation time
    clockTime calcTimer;

    // Compute nearest cell addressing
    calcCellAddressing();

    if (nThreads == 1)
    {
        calcAddressingAndWeights(0, tgtMesh.nCells(), true);
    }
    else
    {
        // Prior to multi-threaded operation,
        // force calculation of demand-driven data
        tgtMesh.cells();
        srcMesh.cells();
        srcMesh.cellCentres();
        srcMesh.cellCells();

        multiThreader threader(nThreads);

        // Set one handler per thread
        PtrList<handler> hdl(threader.getNumThreads());

        forAll(hdl, i)
        {
            hdl.set(i, new handler(*this, threader));
        }

        // Simple, but inefficient load-balancing scheme
        labelList tStarts(threader.getNumThreads(), 0);
        labelList tSizes(threader.getNumThreads(), 0);

        label index = tgtMesh.nCells(), j = 0;

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
            Info<< " Load starts: " << tStarts << endl;
            Info<< " Load sizes: " << tSizes << endl;
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
            sleep(1.0);

            ctrMutex_.lock();

            scalar percent = 100.0 * (double(counter_) / tgtMesh.nCells());

            Info<< "  Progress: " << percent << "% : "
                << "  Cells processed: " << counter_
                << "  out of " << tgtMesh.nCells() << " total."
                << "             \r"
                << flush;

            ctrMutex_.unlock();

            if (counter_ == tgtMesh.nCells())
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

    Info<< nl << " Calculation time: " << calcTimer.elapsedTime() << endl;

    if (writeAddressing)
    {
        Info<< " Writing addressing to disk." << endl;

        addressing_.write();
        weights_.write();
        centres_.write();
    }

    forAll (tgtMesh.boundaryMesh(), patchI)
    {
        const polyPatch& tgtPatch = tgtMesh.boundaryMesh()[patchI];

        label patchID = srcMesh.boundaryMesh().findPatchID(tgtPatch.name());

        if (patchID == -1)
        {
            FatalErrorIn
            (
                "\n\n"
                "void conservativeMeshToMesh::conservativeMeshToMesh()\n"
            )   << " Could not find " << tgtPatch.name()
                << " in the source mesh."
                << exit(FatalError);
        }

        const polyPatch& srcPatch = srcMesh.boundaryMesh()[patchID];

        if (srcPatch.size() == 0)
        {
            WarningIn("meshToMesh::calcAddressing()")
                << "Source patch " << srcPatch.name()
                << " has no faces. Not performing mapping for it."
                << endl;
            boundaryAddressing_[patchI] = -1;
        }
        else
        {
            treeBoundBox wallBb(srcPatch.localPoints());

            scalar typDim =
            (
                wallBb.avgDim()/(2.0*sqrt(scalar(srcPatch.size())))
            );

            treeBoundBox shiftedBb
            (
                wallBb.min(),
                wallBb.max() + vector(typDim, typDim, typDim)
            );

            // Wrap data for octree into container
            octreeDataFace shapes(srcPatch);

            octree<octreeDataFace> oc
            (
                shiftedBb,  // overall search domain
                shapes,     // all information needed to do checks on cells
                1,          // min levels
                20.0,       // maximum ratio of cubes v.s. cells
                2.0
            );

            const vectorField::subField centresToBoundary =
            (
                tgtPatch.faceCentres()
            );

            boundaryAddressing_[patchI].setSize(tgtPatch.size());

            scalar tightestDist;
            treeBoundBox tightest;

            forAll(tgtPatch, tgtI)
            {
                tightest = wallBb;                 // starting search bb
                tightestDist = Foam::GREAT;        // starting max distance

                boundaryAddressing_[patchI][tgtI] =
                (
                    oc.findNearest
                    (
                        centresToBoundary[tgtI],
                        tightest,
                        tightestDist
                    )
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

conservativeMeshToMesh::~conservativeMeshToMesh()
{}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

// Return source mesh
const fvMesh& conservativeMeshToMesh::srcMesh() const
{
    return meshSrc_;
}


// Return target mesh
const fvMesh& conservativeMeshToMesh::tgtMesh() const
{
    return meshTgt_;
}


// Fetch cell addressing
const labelList& conservativeMeshToMesh::cellAddressing() const
{
    return cellAddressing_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
