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

#if USE_MPFR
template<>
const char* const mpVector::typeName = "mpVector";

template<>
const char* mpVector::componentNames[] = {"x", "y", "z"};

template<>
const mpVector mpVector::zero =
(
   mpVector
   (
       pTraits<mpScalar>::zero,
       pTraits<mpScalar>::zero,
       pTraits<mpScalar>::zero
   )
);

template<>
const mpVector mpVector::one =
(
   mpVector
   (
       pTraits<mpScalar>::one,
       pTraits<mpScalar>::one,
       pTraits<mpScalar>::one
   )
);
#endif

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
    twoDMesh_(false),
    boundaryAddressing_(meshTo.boundaryMesh().size())
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
        toMesh().cells();
        fromMesh().cells();
        fromMesh().cellCentres();
        fromMesh().cellCells();

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
            sleep(0.5);

            ctrMutex_.lock();

            scalar percent = 100.0 * (double(counter_) / toMesh().nCells());

            Info << "  Progress: " << percent << "% : "
                 << "  Cells processed: " << counter_
                 << "  out of " << toMesh().nCells() << " total."
                 << "             \r"
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

    forAll (meshTo.boundaryMesh(), patchi)
    {
        const polyPatch& toPatch = meshTo.boundaryMesh()[patchi];

        label patchID = meshFrom.boundaryMesh().findPatchID(toPatch.name());

        if (patchID == -1)
        {
            FatalErrorIn
            (
                "\n\n"
                "void conservativeMeshToMesh::conservativeMeshToMesh()\n"
            )   << " Could not find " << toPatch.name()
                << " in the source mesh."
                << exit(FatalError);
        }

        const polyPatch& fromPatch = meshFrom.boundaryMesh()[patchID];

        if (fromPatch.size() == 0)
        {
            WarningIn("meshToMesh::calcAddressing()")
                << "Source patch " << fromPatch.name()
                << " has no faces. Not performing mapping for it."
                << endl;
            boundaryAddressing_[patchi] = -1;
        }
        else
        {
            treeBoundBox wallBb(fromPatch.localPoints());

            scalar typDim =
            (
                wallBb.avgDim()/(2.0*sqrt(scalar(fromPatch.size())))
            );

            treeBoundBox shiftedBb
            (
                wallBb.min(),
                wallBb.max() + vector(typDim, typDim, typDim)
            );

            // Wrap data for octree into container
            octreeDataFace shapes(fromPatch);

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
                toPatch.faceCentres()
            );

            boundaryAddressing_[patchi].setSize(toPatch.size());

            scalar tightestDist;
            treeBoundBox tightest;

            forAll(toPatch, toi)
            {
                tightest = wallBb;                 // starting search bb
                tightestDist = Foam::GREAT;        // starting max distance

                boundaryAddressing_[patchi][toi] = oc.findNearest
                (
                    centresToBoundary[toi],
                    tightest,
                    tightestDist
                );
            }
        }
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
            fromMesh(),
            name,
            cList,
            primitiveType,
            fromMesh().points(),
            fromMesh().edges(),
            fromMesh().faces(),
            fromMesh().cells(),
            fromMesh().faceOwner(),
            field
        );
    }
    else
    {
        meshOps::writeVTK
        (
            toMesh(),
            name,
            cList,
            primitiveType,
            toMesh().points(),
            toMesh().edges(),
            toMesh().faces(),
            toMesh().cells(),
            toMesh().faceOwner(),
            field
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
