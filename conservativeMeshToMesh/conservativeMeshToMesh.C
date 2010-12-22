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
#include "Random.H"
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
    const bool writeAddressing,
    const bool decompSource,
    const bool decompTarget
)
:
    meshFrom_(meshFrom),
    meshTo_(meshTo),
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
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        meshTo.nCells()
    ),
    volumes_
    (
        IOobject
        (
            "volumes",
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
    meshToMeshPtr_.set
    (
        new meshToMesh
        (
            decomposeMesh
            (
                meshFrom,
                decompSource,
                srcTetMesh_,
                srcTetFvMesh_,
                srcTetStarts_,
                srcTetSizes_
            ),
            decomposeMesh
            (
                meshTo,
                decompTarget,
                tgtTetMesh_,
                tgtTetFvMesh_,
                tgtTetStarts_,
                tgtTetSizes_
            )
        )
    );

    if (addressing_.headerOk() && volumes_.headerOk() && centres_.headerOk())
    {
        // Check if sizes match. Otherwise, re-calculate.
        if
        (
            addressing_.size() == meshTo.nCells() &&
            volumes_.size() == meshTo.nCells() &&
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

    // If the target mesh is decomposed, sizes would be wrong,
    // so fix it temporarily for the addressing calculation.
    // This will be readjusted later to actual addressing
    if (decompTarget)
    {
        addressing_.setSize(tgtTetFvMesh_().nCells());
        weights_.setSize(tgtTetFvMesh_().nCells());
        volumes_.setSize(tgtTetFvMesh_().nCells());
        centres_.setSize(tgtTetFvMesh_().nCells());
    }

    labelList srcPolyCell;

    // Build a map for each target tet to target polyhedral cell
    if (decompSource)
    {
        srcPolyCell.setSize(srcTetFvMesh_().nCells(), -1);

        forAll(srcTetStarts_, cellI)
        {
            label start = srcTetStarts_[cellI];
            label size = srcTetSizes_[cellI];

            // Select a random integer for this cell
            label randInt = Random(cellI).integer(0, 1000);

            for (label i = 0; i < size; i++)
            {
                srcPolyCell[start + i] = cellI;
                srcCellIndex_.insert(start + i, randInt);
            }
        }
    }
    else
    {
        // Identity map
        srcPolyCell = identity(meshFrom.nCells());

        // Assign a random index for each cell on the source
        forAll(srcPolyCell, cellI)
        {
            srcCellIndex_.insert(cellI, Random(cellI).integer(0, 1000));
        }
    }

    // Track calculation time
    clockTime calcTimer;

    if (nThreads == 1)
    {
        calcAddressingAndWeights(0, tgtMesh().nCells(), true);
    }
    else
    {
        // Prior to multi-threaded operation,
        // force calculation of demand-driven data
        tgtMesh().cells();
        srcMesh().cells();
        srcMesh().cellCentres();
        srcMesh().cellCells();

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

        label index = tgtMesh().nCells(), j = 0;

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
            sleep(1.0);

            ctrMutex_.lock();

            scalar percent = 100.0 * (double(counter_) / toMesh().nCells());

            Info << "  Progress: " << percent << "% : "
                 << "  Cells processed: " << counter_
                 << "  out of " << tgtMesh().nCells() << " total."
                 << "             \r"
                 << flush;

            ctrMutex_.unlock();

            if (counter_ == tgtMesh().nCells())
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

    if (decompTarget)
    {
        // Agglomerate intersections from tets for each polyhedral cell
        List<labelList> newAddressing(meshTo.nCells());
        List<scalarField> newVolumes(meshTo.nCells());
        List<vectorField> newCentres(meshTo.nCells());

        Map<scalar> polyVol;
        Map<vector> polyCtr;

        forAll(newAddressing, cellI)
        {
            // First count the number of intersections per tet
            labelList& cellAddr = newAddressing[cellI];
            scalarField& cellVols = newVolumes[cellI];
            vectorField& cellCtrs = newCentres[cellI];

            label start = tgtTetStarts_[cellI];
            label size = tgtTetSizes_[cellI];

            // Clear the map
            polyVol.clear();
            polyCtr.clear();

            for (label i = 0; i < size; i++)
            {
                const labelList& tetAddr = addressing_[start + i];
                const scalarField& tetVol = volumes_[start + i];
                const vectorField& tetCtr = centres_[start + i];

                forAll(tetAddr, j)
                {
                    // Add a new entry if necessary
                    if (!polyVol.found(srcPolyCell[tetAddr[j]]))
                    {
                        polyVol.insert(srcPolyCell[tetAddr[j]], 0.0);
                        polyCtr.insert(srcPolyCell[tetAddr[j]], vector::zero);
                    }

                    // Accumulate volumes and weighted-centroids
                    polyVol[srcPolyCell[tetAddr[j]]] += tetVol[j];
                    polyCtr[srcPolyCell[tetAddr[j]]] += tetCtr[j] * tetVol[j];
                }
            }

            // Set the addressing
            cellAddr = polyVol.toc();

            cellVols.setSize(cellAddr.size());
            cellCtrs.setSize(cellAddr.size());

            forAll(cellAddr, polyI)
            {
                cellVols[polyI] = polyVol[cellAddr[polyI]];
                cellCtrs[polyI] = polyCtr[cellAddr[polyI]] / cellVols[polyI];
            }
        }

        // Transfer new lists
        addressing_.transfer(newAddressing);
        volumes_.transfer(newVolumes);
        centres_.transfer(newCentres);
    }
    else
    if (decompSource)
    {
        // Source is decomposed. Renumber to correct addressing.
        // Volumes and centres are already okay.
        forAll(addressing_, cellI)
        {
            labelList& addr = addressing_[cellI];

            forAll(addr, i)
            {
                addr[i] = srcPolyCell[addr[i]];
            }
        }
    }

    if (writeAddressing)
    {
        Info << " Writing addressing to disk." << endl;

        addressing_.write();
        volumes_.write();
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

// Return basic source mesh (undecomposed)
const fvMesh& conservativeMeshToMesh::origSrcMesh() const
{
    return meshFrom_;
}


// Return basic target mesh (undecomposed)
const fvMesh& conservativeMeshToMesh::origTgtMesh() const
{
    return meshTo_;
}


// Return (decomposed) source mesh
const fvMesh& conservativeMeshToMesh::fromMesh() const
{
    return meshToMeshPtr_().fromMesh();
}


// Return (decomposed) target mesh
const fvMesh& conservativeMeshToMesh::toMesh() const
{
    return meshToMeshPtr_().toMesh();
}


// Return (decomposed) source mesh
const fvMesh& conservativeMeshToMesh::srcMesh() const
{
    if (srcTetFvMesh_.valid())
    {
        return srcTetFvMesh_();
    }

    return meshFrom_;
}


// Return (decomposed) target mesh
const fvMesh& conservativeMeshToMesh::tgtMesh() const
{
    if (tgtTetFvMesh_.valid())
    {
        return tgtTetFvMesh_();
    }

    return meshTo_;
}


// Decompose a given mesh, using tetDecomposition
const fvMesh& conservativeMeshToMesh::decomposeMesh
(
    const fvMesh& initialMesh,
    const bool decomp,
    autoPtr<tetPolyMesh>& tetDecomp,
    autoPtr<fvMesh>& tetDecompFvMesh,
    labelList& tetStarts,
    labelList& tetSizes
)
{
    if (!decomp)
    {
        // Do nothing. Return original mesh.
        return initialMesh;
    }

    // Reset pointer to null
    tetDecomp.clear();

    // Construct the tetPolyMesh
    tetDecomp.set(new tetPolyMesh(initialMesh));

    const tetPolyMesh& initTets = tetDecomp();

    label initialCells = initialMesh.nCells();

    Info<< " Creating decomposed tet mesh."
        << " Initial cells: " << initialCells
        << " Number of tets: " << initTets.nTets()
        << endl;

    tetStarts.setSize(initialCells, 0);
    tetSizes.setSize(initialCells, 0);

    // Set addressing from polyhedra to tets
    label nTotCells = 0;

    forAll(tetStarts, i)
    {
        label nTetsForCell = initTets.nTetsForCell(i);

        tetStarts[i] = nTotCells;
        tetSizes[i] = nTetsForCell;

        nTotCells += nTetsForCell;
    }

    if (nTotCells != initTets.nTets())
    {
        FatalErrorIn
        (
            "\n\n"
            "const fvMesh& conservativeMeshToMesh::decomposeMesh()\n"
        )   << " Incorrect addressing. " << nl
            << " nTets: " << initTets.nTets() << nl
            << " nTotCells: " << nTotCells << nl
            << exit(FatalError);
    }

    // Construct fvMesh from polyMesh, since a shapes
    // constructor doesn't exist
    polyMesh newMesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            initialMesh.time().timeName(),
            initialMesh.time()
        ),
        xferCopy(initTets.points()()),
        initTets.tetCells(),
        initTets.boundary().boundaryTriFaces(),
        initialMesh.boundaryMesh().names(),
        initialMesh.boundaryMesh().types(),
        polyPatch::typeName,
        "defaultPatch",
        initialMesh.boundaryMesh().physicalTypes()
    );

    // Now set the fvMesh from the decomposed polyMesh
    tetDecompFvMesh.set
    (
        new fvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                initialMesh.time().timeName(),
                initialMesh.time()
            ),
            xferCopy(newMesh.points()),
            xferCopy(newMesh.faces()),
            xferCopy(newMesh.faceOwner()),
            xferCopy(newMesh.faceNeighbour())
        )
    );

    // Add patches to newly created mesh
    const polyBoundaryMesh& newBoundary = newMesh.boundaryMesh();

    List<polyPatch*> patchList(newBoundary.size());

    forAll(newBoundary, patchI)
    {
        patchList[patchI] =
        (
            polyPatch::New
            (
                newBoundary[patchI].type(),
                newBoundary[patchI].name(),
                newBoundary[patchI].size(),
                newBoundary[patchI].start(),
                patchI,
                tetDecompFvMesh().boundaryMesh()
            ).ptr()
        );
    }

    // Add patches
    tetDecompFvMesh().addFvPatches(patchList);

    return tetDecompFvMesh();
}


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
            fromMesh(),
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
