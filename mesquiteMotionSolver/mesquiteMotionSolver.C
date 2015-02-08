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

#include "mesquiteMotionSolver.H"
#include "Random.H"
#include "IOmanip.H"
#include "SortableList.H"
#include "globalMeshData.H"
#include "cyclicPolyPatch.H"
#include "processorPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"

#include "tetrahedron.H"
#include "matchPoints.H"
#include "DimensionedField.H"
#include "pointPatchField.H"
#include "pointMesh.H"
#include "mapPolyMesh.H"

#include "hexMatcher.H"
#include "tetMatcher.H"
#include "pyrMatcher.H"
#include "prismMatcher.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mesquiteMotionSolver, 0);
    addToRunTimeSelectionTable
    (
        motionSolver,
        mesquiteMotionSolver,
        dictionary
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mesquiteMotionSolver::mesquiteMotionSolver
(
    const polyMesh& mesh
)
:
    motionSolver(mesh),
    MeshObject
    <
        polyMesh,
        Foam::UpdateableMeshObject,
        mesquiteMotionSolver
    >(mesh),
    Mesh_(mesh),
    twoDMesh_(mesh.nGeometricD() == 2 ? true : false),
    arraysInitialized_(false),
    nPoints_(mesh.nPoints()),
    nCells_(0),
    decompType_(-1),
    nDecompPoints_(0),
    surfaceSmoothing_(false),
    volumeCorrection_(false),
    volCorrTolerance_(1e-20),
    volCorrMaxIter_(100),
    tolerance_(1e-4),
    nSweeps_(1),
    surfInterval_(1),
    relax_(1.0),
    refPoints_
    (
        IOobject
        (
            "refPoints",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh.points()
    ),
    basePoints_(0),
    boundaryConditions_(0),
    oldVolume_(0.0)
{
    // Read options from the dictionary
    readOptions();
}


mesquiteMotionSolver::mesquiteMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    motionSolver(mesh),
    MeshObject
    <
        polyMesh,
        Foam::UpdateableMeshObject,
        mesquiteMotionSolver
    >(mesh),
    Mesh_(mesh),
    twoDMesh_(mesh.nGeometricD() == 2 ? true : false),
    arraysInitialized_(false),
    nPoints_(mesh.nPoints()),
    nCells_(0),
    decompType_(-1),
    nDecompPoints_(0),
    surfaceSmoothing_(false),
    volumeCorrection_(false),
    volCorrTolerance_(1e-20),
    volCorrMaxIter_(100),
    tolerance_(1e-4),
    nSweeps_(1),
    surfInterval_(1),
    relax_(1.0),
    refPoints_
    (
        IOobject
        (
            "refPoints",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh.points()
    ),
    basePoints_(0),
    boundaryConditions_(0),
    oldVolume_(0.0)
{
    // Read options from the dictionary
    readOptions();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mesquiteMotionSolver::~mesquiteMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Parallel blocking send
void mesquiteMotionSolver::parWrite
(
    const label toID,
    const label& data
)
{
    OPstream::write
    (
        Pstream::blocking,
        toID,
        reinterpret_cast<const char*>(&data),
        sizeof(label)
    );
}


// Parallel blocking receive
void mesquiteMotionSolver::parRead
(
    const label fromID,
    label& data
)
{
    IPstream::read
    (
        Pstream::blocking,
        fromID,
        reinterpret_cast<char*>(&data),
        sizeof(label)
    );
}


// Parallel non-blocking send for lists
template <class Type>
void mesquiteMotionSolver::parWrite
(
    const label toID,
    const UList<Type>& data
)
{
    OPstream::write
    (
        Pstream::nonBlocking,
        toID,
        reinterpret_cast<const char*>(&data[0]),
        data.size()*sizeof(Type)
    );
}


// Parallel non-blocking receive for lists
template <class Type>
void mesquiteMotionSolver::parRead
(
    const label fromID,
    UList<Type>& data
)
{
    IPstream::read
    (
        Pstream::nonBlocking,
        fromID,
        reinterpret_cast<char*>(&data[0]),
        data.size()*sizeof(Type)
    );
}


// Read options from the dictionary
void mesquiteMotionSolver::readOptions()
{
    if (debug)
    {
        Info<< "Reading options for mesquiteMotionSolver" << endl;
    }

    // Fetch the sub-dictionary
    const dictionary& optionsDict = subDict("mesquiteOptions");

    // Check if the 'pointDisplacement' file is to be used
    // for moving boundary conditions
    if (optionsDict.found("usePointDisplacement"))
    {
        const pointMesh& pMesh = pointMesh::New(Mesh_);

        // Check for existence of base points
        IOobject baseIO("basePoints", Mesh_.time().timeName(), Mesh_);

        if (baseIO.headerOk())
        {
            basePoints_.reset
            (
                new pointVectorField
                (
                    IOobject
                    (
                        "basePoints",
                        Mesh_.time().timeName(),
                        Mesh_,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    pMesh
                )
            );
        }
        else
        {
            basePoints_.reset
            (
                new pointVectorField
                (
                    IOobject
                    (
                        "basePoints",
                        Mesh_.time().timeName(),
                        Mesh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    pMesh,
                    dimensioned<vector>("0", dimLength, vector::zero)
                )
            );

            // Initialize with mesh points
            basePoints_().internalField() = Mesh_.points();
        }

        boundaryConditions_.reset
        (
            new pointVectorField
            (
                IOobject
                (
                    "pointDisplacement",
                    Mesh_.time().timeName(),
                    Mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                pMesh
            )
        );
    }

    // Check if any slip patches are specified
    if (optionsDict.found("slipPatches") || twoDMesh_)
    {
        labelHashSet slipPatchIDs;

        const polyBoundaryMesh& boundary = mesh().boundaryMesh();

        // For 2D meshes, add wedge / empty patches
        if (twoDMesh_)
        {
            forAll(boundary, patchI)
            {
                if
                (
                    boundary[patchI].type() == "wedge" ||
                    boundary[patchI].type() == "empty"
                )
                {
                    slipPatchIDs.insert(patchI);
                }
            }

            surfaceSmoothing_ = true;
        }
        else
        {
            wordList slipPatches = optionsDict.subDict("slipPatches").toc();

            forAll(slipPatches, wordI)
            {
                word& patchName = slipPatches[wordI];
                slipPatchIDs.insert(boundary.findPatchID(patchName));

                surfaceSmoothing_ = true;
            }
        }

        // Check if a tolerance has been specified
        if (optionsDict.found("tolerance"))
        {
            tolerance_ = readScalar(optionsDict.lookup("tolerance"));
        }

        // Check if volume correction is enabled
        if (optionsDict.found("volumeCorrection"))
        {
            volumeCorrection_ =
            (
                readBool(optionsDict.lookup("volumeCorrection"))
            );
        }

        // Check if volume correction tolerance is specified
        if (optionsDict.found("volCorrTolerance"))
        {
            volCorrTolerance_ =
            (
                readScalar(optionsDict.lookup("volCorrTolerance"))
            );
        }

        // Check if volume correction maxIter is specified
        if (optionsDict.found("volCorrMaxIter"))
        {
            volCorrMaxIter_ = readLabel(optionsDict.lookup("volCorrMaxIter"));
        }

        // Check if multiple sweeps have been requested
        if (optionsDict.found("nSweeps"))
        {
            nSweeps_ = readLabel(optionsDict.lookup("nSweeps"));
        }

        // Check if a surface smoothing interval has been specified
        if (optionsDict.found("surfInterval"))
        {
            surfInterval_ = readLabel(optionsDict.lookup("surfInterval"));
        }

        // Check if a relaxation factor is specified
        if (optionsDict.found("relaxationFactor"))
        {
            relax_ = readScalar(optionsDict.lookup("relaxationFactor"));
        }

        // Check if coupled patches exist.
        if (optionsDict.found("coupledPatches"))
        {
            const dictionary& coupledPatches =
            (
                optionsDict.subDict("coupledPatches")
            );

            // Determine master and slave patches
            forAllConstIter(dictionary, coupledPatches, dIter)
            {
                const dictionary& dictI = dIter().dict();

                // Lookup the master / slave patches
                word masterPatch = dictI.lookup("master");
                word slavePatch  = dictI.lookup("slave");

                // Determine patch indices
                label mPatch = boundary.findPatchID(masterPatch);
                label sPatch = boundary.findPatchID(slavePatch);

                if (mPatch == -1 && sPatch == -1)
                {
                    continue;
                }

                // Add to the patch-list, if the entry hasn't been found.
                if (!slipPatchIDs.found(sPatch) && (sPatch != -1))
                {
                    slipPatchIDs.insert(sPatch);
                }

                // Add to the list if entries are legitimate
                if
                (
                    mPatch != sPatch &&
                    boundary[mPatch].size() == boundary[sPatch].size()
                )
                {
                    patchCoupling_.insert(mPatch, sPatch);
                }
                else
                {
                    FatalErrorIn("void mesquiteMotionSolver::readOptions()")
                        << " Coupled patches are wrongly specified." << nl
                        << " Master: " << mPatch << ":" << masterPatch << nl
                        << " Slave: " << sPatch << ":" << slavePatch << nl
                        << abort(FatalError);
                }
            }
        }

        // Extract the final slip-patch list
        pIDs_ = slipPatchIDs.toc();
    }

    if (twoDMesh_)
    {
        return;
    }

    // The following applies only to 3D meshes
    Mesquite::MsqError err;

    // Add existing metrics to the hash-table
    qMetricTable_.insert
    (
        "AspectRatioGamma",
        autoPtr<Mesquite::QualityMetric>
        (
            new Mesquite::AspectRatioGammaQualityMetric
        )
    );

    qMetricTable_.insert
    (
        "ConditionNumber",
        autoPtr<Mesquite::QualityMetric>
        (
            new Mesquite::ConditionNumberQualityMetric
        )
    );

    qMetricTable_.insert
    (
        "InverseMeanRatio",
        autoPtr<Mesquite::QualityMetric>
        (
            new Mesquite::IdealWeightInverseMeanRatio
        )
    );

    qMetricTable_.insert
    (
        "MeanRatio",
        autoPtr<Mesquite::QualityMetric>
        (
            new Mesquite::IdealWeightMeanRatio
        )
    );

    qMetricTable_.insert
    (
        "VertexConditionNumber",
        autoPtr<Mesquite::QualityMetric>
        (
            new Mesquite::VertexConditionNumberQualityMetric
        )
    );

    qMetricTable_.insert
    (
        "EdgeLengthQuality",
        autoPtr<Mesquite::QualityMetric>
        (
            new Mesquite::EdgeLengthQualityMetric
        )
    );

    qMetricTable_.insert
    (
        "EdgeLength",
        autoPtr<Mesquite::QualityMetric>
        (
            new Mesquite::EdgeLengthMetric
        )
    );

    qMetricTable_.insert
    (
        "UntangleBeta",
        autoPtr<Mesquite::QualityMetric>
        (
            new Mesquite::UntangleBetaQualityMetric
        )
    );

    qMetricTable_.insert
    (
        "LocalSize",
        autoPtr<Mesquite::QualityMetric>
        (
            new Mesquite::LocalSizeQualityMetric
        )
    );

    // Size up metric list to contain at least one entry
    qMetrics_.setSize(1);

    // Read the optimization metric
    if (optionsDict.isDict("optMetric"))
    {
        // Metrics are in dictionary form
        const dictionary& metricDict = optionsDict.subDict("optMetric");

        // Read first metric
        qMetrics_[0] = word(metricDict.lookup("firstMetric"));

        // Optionally read the second metric
        if (metricDict.found("secondMetric"))
        {
            qMetrics_.setSize(2);
            qMetrics_[1] = word(metricDict.lookup("secondMetric"));
        }

        forAll(qMetrics_, wordI)
        {
            if (!qMetricTable_.found(qMetrics_[wordI]))
            {
                FatalErrorIn("void mesquiteMotionSolver::readOptions()")
                    << "Unrecognized quality metric: "
                    << qMetrics_[wordI] << nl
                    << "Available metrics are: "
                    << nl << qMetricTable_.toc()
                    << abort(FatalError);
            }
            else
            {
                Info<< "Selecting quality metric: "
                    << qMetrics_[wordI] << endl;
            }
        }
    }
    else
    {
        // Conventional keyword entry

        // Alias for convenience...
        word& qMetric = qMetrics_[0];

        qMetric = word(optionsDict.lookup("optMetric"));

        if (!qMetricTable_.found(qMetric))
        {
            FatalErrorIn("void mesquiteMotionSolver::readOptions()")
                << "Unrecognized quality metric: " << qMetric << nl
                << "Available metrics are: " << nl << qMetricTable_.toc()
                << abort(FatalError);
        }
        else
        {
            Info<< "Selecting quality metric: " << qMetric << endl;
        }
    }

    // Define the objective function table,
    // and add existing possibilities
    HashTable<label> ofTable;
    ofTable.insert("CompositeOFAdd", 0);
    ofTable.insert("CompositeOFMultiply", 1);
    ofTable.insert("CompositeOFScalarAdd", 2);
    ofTable.insert("CompositeOFScalarMultiply", 3);
    ofTable.insert("LInf", 4);
    ofTable.insert("LPtoP", 5);
    ofTable.insert("Max", 6);
    ofTable.insert("PMeanP", 7);
    ofTable.insert("StdDev", 8);
    ofTable.insert("Variance", 9);
    ofTable.insert("PatchPowerMeanP", 10);

    // Read the objective function
    word ofType(optionsDict.lookup("objFunction"));

    if (!ofTable.found(ofType))
    {
        FatalErrorIn("void mesquiteMotionSolver::readOptions()")
            << "Unrecognized objective function: " << ofType << nl
            << "Available types are: " << nl << ofTable.toc()
            << abort(FatalError);
    }
    else
    {
        Info<< "Selecting objective function: " << ofType << endl;
    }

    // Instantiate the appropriate objective function
    scalar scale = 0.0;
    label numTries = 0;
    label ofSelection = ofTable[ofType];
    FixedList<label, 2> types(-1);

    // Check if a composite function is requested
    if (ofSelection == 0 || ofSelection == 1)
    {
        numTries = 2;

        // Lookup both objective functions
        types[0] = ofTable[word(optionsDict.lookup("firstFunction"))];
        types[1] = ofTable[word(optionsDict.lookup("secondFunction"))];

        // Ensure that we're not making a composite of composites
        if (types[0] < 4 || types[1] < 4)
        {
            FatalErrorIn("void mesquiteMotionSolver::readOptions()")
                << "Cannot make a composite of composite functions."
                << abort(FatalError);
        }

        // Make sure we have two quality metrics
        if (qMetrics_.size() < 2)
        {
            FatalErrorIn("void mesquiteMotionSolver::readOptions()")
                << "Two quality metrics are needed for composite functions."
                << abort(FatalError);
        }
    }
    else
    // Check if a scaled function is requested
    if (ofSelection == 2 || ofSelection == 3)
    {
        numTries = 1;

        // Lookup the objective function to scale
        types[0] = ofTable[word(optionsDict.lookup("scaleFunction"))];

        // Lookup the scale value
        scale = readScalar(optionsDict.lookup("scale"));

        // Ensure that we're not making a composite of composites
        if (types[0] < 4)
        {
            FatalErrorIn("void mesquiteMotionSolver::readOptions()")
                << "Cannot make a composite of composite functions."
                << abort(FatalError);
        }
    }
    else
    {
        // Simple function
        types[0] = ofSelection;

        numTries = 1;
    }

    // Hold two additional pointers for composite functions
    List<autoPtr<Mesquite::ObjectiveFunction> > composite(numTries);

    for (label i = 0; i < numTries; i++)
    {
        switch (types[i])
        {
            case 4:
            {
                composite[i].set
                (
                    new Mesquite::LInfTemplate
                    (
                        &qMetricTable_[qMetrics_[i]]()
                    )
                );

                break;
            }

            case 5:
            {
                // Lookup the P value
                label pValue = readLabel(optionsDict.lookup("pValue"));

                composite[i].set
                (
                    new Mesquite::LPtoPTemplate
                    (
                        &qMetricTable_[qMetrics_[i]](),
                        pValue,
                        err
                    )
                );

                break;
            }

            case 6:
            {
                // Same as LInfTemplate, but absolute values are not used
                composite[i].set
                (
                    new Mesquite::MaxTemplate
                    (
                        &qMetricTable_[qMetrics_[i]]()
                    )
                );

                break;
            }

            case 7:
            {
                // Lookup the Power value
                scalar power = readScalar(optionsDict.lookup("power"));

                // Sum of each quality metric value raised to a power,
                // divided by the total number of quality metric values.
                composite[i].set
                (
                    new Mesquite::PMeanPTemplate
                    (
                        power,
                        &qMetricTable_[qMetrics_[i]]()
                    )
                );

                break;
            }

            case 8:
            {
                composite[i].set
                (
                    new Mesquite::StdDevTemplate
                    (
                        &qMetricTable_[qMetrics_[i]]()
                    )
                );

                break;
            }

            case 9:
            {
                composite[i].set
                (
                    new Mesquite::VarianceTemplate
                    (
                        &qMetricTable_[qMetrics_[i]]()
                    )
                );

                break;
            }

            case 10:
            {
                // Lookup the Power value
                scalar power = readScalar(optionsDict.lookup("power"));

                // p-mean^p of p-mean^p of metric values
                composite[i].set
                (
                    new Mesquite::PatchPowerMeanP
                    (
                        power,
                        &qMetricTable_[qMetrics_[i]]()
                    )
                );

                break;
            }

            default:
            {
                FatalErrorIn("void mesquiteMotionSolver::readOptions()")
                    << "Illegal function requested."
                    << abort(FatalError);

                break;
            }
        }
    }

    // Compile the final objective function,
    // autoPtr clears its contents for composite functions
    switch (ofSelection)
    {
        case 0:
        {
            // Composite of add
            objFunction_.set
            (
                new Mesquite::CompositeOFAdd
                (
                    composite[0].ptr(),
                    composite[1].ptr()
                )
            );

            break;
        }

        case 1:
        {
            // Composite of multiply
            objFunction_.set
            (
                new Mesquite::CompositeOFMultiply
                (
                    composite[0].ptr(),
                    composite[1].ptr()
                )
            );

            break;
        }

        case 2:
        {
            // Composite of scalar add
            objFunction_.set
            (
                new Mesquite::CompositeOFScalarAdd
                (
                    scale,
                    composite[0].ptr()
                )
            );

            break;
        }

        case 3:
        {
            // Composite of scalar multiply
            objFunction_.set
            (
                new Mesquite::CompositeOFScalarMultiply
                (
                    scale,
                    composite[0].ptr()
                )
            );

            break;
        }

        default:
        {
            // Simply copy the pointer
            objFunction_ = composite[0];

            break;
        }
    }

    // Define the optimization algorithm table,
    // and add existing possibilities
    HashTable<label> oaTable;
    oaTable.insert("ConjugateGradient", 0);
    oaTable.insert("FeasibleNewton", 1);
    oaTable.insert("NonSmoothDescent", 2);
    oaTable.insert("QuasiNewton", 3);
    oaTable.insert("Randomize", 4);
    oaTable.insert("Laplacian", 5);
    oaTable.insert("SmartLaplacian", 6);
    oaTable.insert("SteepestDescent", 7);
    oaTable.insert("TrustRegion", 8);

    // Read the optimization algorithm
    word oaType(optionsDict.lookup("optAlgorithm"));

    if (!oaTable.found(oaType))
    {
        FatalErrorIn("void mesquiteMotionSolver::readOptions()")
            << "Unrecognized optimization algorithm: " << oaType << nl
            << "Available types are: " << nl << oaTable.toc()
            << abort(FatalError);
    }
    else
    {
        Info<< "Selecting optimization algorithm: " << oaType << endl;
    }

    // Instantiate the appropriate objective function
    label oaSelection = oaTable[oaType];

    switch(oaSelection)
    {
        case 0:
        {
            optAlgorithm_.set(new Mesquite::ConjugateGradient(&objFunction_()));

            break;
        }

        case 1:
        {
            optAlgorithm_.set(new Mesquite::FeasibleNewton(&objFunction_()));

            break;
        }

        case 2:
        {
            // NonSmoothDescent optimizes on a quality metric only
            optAlgorithm_.set
            (
                new Mesquite::NonSmoothDescent
                (
                    dynamic_cast<Mesquite::ElementQM*>
                    (
                        &qMetricTable_[qMetrics_[0]]()
                    )
                )
            );

            break;
        }

        case 3:
        {
            optAlgorithm_.set(new Mesquite::QuasiNewton(&objFunction_()));

            break;
        }

        case 4:
        {
            // Not really a smoother; just randomizes positions
            optAlgorithm_.set(new Mesquite::Randomize());

            break;
        }

        case 5:
        {
            optAlgorithm_.set(new Mesquite::LaplacianSmoother(&objFunction_()));

            break;
        }

        case 6:
        {
            optAlgorithm_.set
            (
                new Mesquite::SmartLaplacianSmoother(&objFunction_())
            );

            break;
        }

        case 7:
        {
            optAlgorithm_.set(new Mesquite::SteepestDescent(&objFunction_()));

            break;
        }

        case 8:
        {
            optAlgorithm_.set(new Mesquite::TrustRegion(&objFunction_()));

            break;
        }
    }

    // Read termination criteria, if it exists.
    if (optionsDict.found("tcInner"))
    {
        const dictionary& innerDict = optionsDict.subDict("tcInner");

        label nSetOptions = 0;
        HashSet<word> options;

        options.insert("absGradL2");

        if (innerDict.found("absGradL2"))
        {
            tcInner_.add_absolute_gradient_L2_norm
            (
                readScalar(innerDict.lookup("absGradL2"))
            );

            nSetOptions++;
        }

        options.insert("absGradInf");

        if (innerDict.found("absGradInf"))
        {
            tcInner_.add_absolute_gradient_inf_norm
            (
                readScalar(innerDict.lookup("absGradInf"))
            );

            nSetOptions++;
        }

        options.insert("relGradL2");

        if (innerDict.found("relGradL2"))
        {
            tcInner_.add_relative_gradient_L2_norm
            (
                readScalar(innerDict.lookup("relGradL2"))
            );

            nSetOptions++;
        }

        options.insert("relGradInf");

        if (innerDict.found("relGradInf"))
        {
            tcInner_.add_relative_gradient_inf_norm
            (
                readScalar(innerDict.lookup("relGradInf"))
            );

            nSetOptions++;
        }

        options.insert("absQualImprovement");

        if (innerDict.found("absQualImprovement"))
        {
            tcInner_.add_absolute_quality_improvement
            (
                readScalar(innerDict.lookup("absQualImprovement"))
            );

            nSetOptions++;
        }

        options.insert("relQualImprovement");

        if (innerDict.found("relQualImprovement"))
        {
            tcInner_.add_relative_quality_improvement
            (
                readScalar(innerDict.lookup("relQualImprovement"))
            );

            nSetOptions++;
        }

        options.insert("absVertexMovement");

        if (innerDict.found("absVertexMovement"))
        {
            tcInner_.add_absolute_vertex_movement
            (
                readScalar(innerDict.lookup("absVertexMovement"))
            );

            nSetOptions++;
        }

        options.insert("relVertexMovement");

        if (innerDict.found("relVertexMovement"))
        {
            tcInner_.add_relative_vertex_movement
            (
                readScalar(innerDict.lookup("relVertexMovement"))
            );

            nSetOptions++;
        }

        options.insert("absSuccessiveImprovement");

        if (innerDict.found("absSuccessiveImprovement"))
        {
            tcInner_.add_absolute_successive_improvement
            (
                readScalar(innerDict.lookup("absSuccessiveImprovement"))
            );

            nSetOptions++;
        }

        options.insert("relSuccessiveImprovement");

        if (innerDict.found("relSuccessiveImprovement"))
        {
            tcInner_.add_relative_successive_improvement
            (
                readScalar(innerDict.lookup("relSuccessiveImprovement"))
            );

            nSetOptions++;
        }

        options.insert("iterationLimit");

        if (innerDict.found("iterationLimit"))
        {
            tcInner_.add_iteration_limit
            (
                readLabel(innerDict.lookup("iterationLimit"))
            );

            nSetOptions++;
        }

        options.insert("cpuTime");

        if (innerDict.found("cpuTime"))
        {
            tcInner_.add_cpu_time
            (
                readScalar(innerDict.lookup("cpuTime"))
            );

            nSetOptions++;
        }

        if (innerDict.size() == 0 || nSetOptions == 0)
        {
            FatalErrorIn("void mesquiteMotionSolver::readOptions()")
                << "Empty tcInner dictionary: " << nl
                << "Available types are: " << nl
                << options.toc()
                << abort(FatalError);
        }
    }
    else
    {
        Info<< "Inner termination criterion (tcInner) "
            << "was not found. Using default values."
            << endl;

        tcInner_.add_absolute_gradient_inf_norm(1e-4);
    }

    if (optionsDict.found("tcOuter"))
    {
        const dictionary& outerDict = optionsDict.subDict("tcOuter");

        label nSetOptions = 0;
        HashSet<word> options;

        options.insert("iterationLimit");

        if (outerDict.found("iterationLimit"))
        {
            tcOuter_.add_iteration_limit
            (
                readLabel(outerDict.lookup("iterationLimit"))
            );

            nSetOptions++;
        }

        options.insert("cpuTime");

        if (outerDict.found("cpuTime"))
        {
            tcOuter_.add_cpu_time
            (
                readScalar(outerDict.lookup("cpuTime"))
            );

            nSetOptions++;
        }

        if (outerDict.size() == 0 || nSetOptions == 0)
        {
            FatalErrorIn("void mesquiteMotionSolver::readOptions()")
                << "Empty tcOuter dictionary: " << nl
                << "Available types are: " << nl
                << options.toc()
                << abort(FatalError);
        }
    }
    else
    {
        Info<< "Outer termination criterion (tcOuter) "
            << "was not found. Using default values."
            << endl;

        tcOuter_.add_iteration_limit(1);
    }

    // Set decomposition type, if available
    if (optionsDict.found("decompType"))
    {
        word dType(optionsDict.lookup("decompType"));

        if (dType == "cell")
        {
            decompType_ = 1;
        }
        else
        if (dType == "face")
        {
            decompType_ = 2;
        }
        else
        {
            FatalErrorIn("void mesquiteMotionSolver::readOptions()")
                << "Unknown decomposition type: " << dType
                << "Available types are: cell or face"
                << abort(FatalError);
        }
    }
    else
    {
        // Set default
        decompType_ = 2;
    }
}


// Initialize connectivity arrays for Mesquite
void mesquiteMotionSolver::initArrays()
{
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    if (surfaceSmoothing_)
    {
        offsets_.setSize(pIDs_.size() + 1, 0);
        pNormals_.setSize(pIDs_.size());
        gradEdge_.setSize(pIDs_.size());
        localPts_.setSize(pIDs_.size());
        edgeMarker_.setSize(pIDs_.size());
        edgeConstant_.setSize(pIDs_.size());

        label totalSize = 0;

        forAll(pIDs_, patchI)
        {
            label nPts = boundary[pIDs_[patchI]].nPoints();
            label nEdg = boundary[pIDs_[patchI]].nEdges();

            pNormals_[patchI].setSize(nPts, vector::zero);
            localPts_[patchI].setSize(nPts, vector::zero);
            gradEdge_[patchI].setSize(nEdg, vector::zero);
            edgeMarker_[patchI].setSize(nEdg, 1.0);
            edgeConstant_[patchI].setSize(nEdg, 1.0);

            // Accumulate the total size
            totalSize += nPts;

            // Set offsets
            offsets_[patchI + 1] = totalSize;
        }

        // Initialize CG variables
        bV_.setSize(totalSize, vector::zero);
        xV_.setSize(totalSize, vector::zero);
        pV_.setSize(totalSize, vector::zero);
        rV_.setSize(totalSize, vector::zero);
        wV_.setSize(totalSize, vector::zero);
        bdy_.setSize(totalSize, vector::one);
        pointMarker_.setSize(totalSize, 1.0);

        // Prepare the boundary condition vectorField
        forAll(pIDs_, patchI)
        {
            const label pOffset = offsets_[patchI];
            const edgeList& edges = boundary[pIDs_[patchI]].edges();

            for
            (
                label i = boundary[pIDs_[patchI]].nInternalEdges();
                i < edges.size();
                i++
            )
            {
                bdy_[edges[i][0] + pOffset] = vector::zero;
                bdy_[edges[i][1] + pOffset] = vector::zero;
            }
        }

        origPoints_.setSize(refPoints_.size(), vector::zero);
    }

    if (twoDMesh_)
    {
        // Initialise parallel connectivity, if necessary
        initParallelConnectivity();

        // Optionally fix additional zone points
        initFixedZones();

        // Optionally fix additional boundary layer zone points
        initBoundaryLayerZones();

        // Set the flag
        arraysInitialized_ = true;

        return;
    }

    const dictionary& optionsDict = subDict("mesquiteOptions");

    label nPolyhedra = 0;
    label nCellPoints = 0, nFacePoints = 0;

    // Construct shape recognizers
    FixedList<label, 4> nTypes(0);
    FixedList<autoPtr<cellMatcher>, 4> matcher;
    FixedList<Mesquite::EntityTopology, 4> cellType;

    // Set types
    bool standardTypes = false;

    if (optionsDict.found("standardCellTypes"))
    {
        standardTypes = readBool(optionsDict.lookup("standardCellTypes"));
    }

    if (standardTypes)
    {
        matcher[0].set(new tetMatcher());
        cellType[0] = Mesquite::TETRAHEDRON;

        matcher[1].set(new hexMatcher());
        cellType[1] = Mesquite::HEXAHEDRON;

        matcher[2].set(new pyrMatcher());
        cellType[2] = Mesquite::PYRAMID;

        matcher[3].set(new prismMatcher());
        cellType[3] = Mesquite::PRISM;
    }
    else
    {
        // No standard types.
        // Default to tetrahedra, and decompose the rest.
        matcher[0].set(new tetMatcher());
        cellType[0] = Mesquite::TETRAHEDRON;

        matcher[1].set(NULL);
        matcher[2].set(NULL);
        matcher[3].set(NULL);
    }

    const faceList& meshFaces = mesh().faces();
    const cellList& meshCells = mesh().cells();
    const labelList& faceOwner = mesh().faceOwner();

    // Face marker list
    List<labelPair> faceMarker;

    if (decompType_ == 2)
    {
        faceMarker.setSize(meshFaces.size(), labelPair(-1, -1));
    }

    forAll(meshCells, cellI)
    {
        bool foundMatch = false;

        forAll(matcher, indexI)
        {
            if
            (
                matcher[indexI].valid() &&
                matcher[indexI]->isA(mesh(), cellI)
            )
            {
                nCells_++;
                nTypes[indexI]++;

                foundMatch = true;
                break;
            }
        }

        if (foundMatch)
        {
            continue;
        }

        // Assume polyhedron
        nPolyhedra++;

        // Decompose into tets.
        const cell& checkCell = meshCells[cellI];

        // Add a central point for the cell
        nCellPoints++;
        nDecompPoints_++;

        forAll(checkCell, faceI)
        {
            const label fIndex = checkCell[faceI];
            const face& checkFace = meshFaces[fIndex];

            switch (decompType_)
            {
                case 1:
                {
                    nCells_ += (checkFace.size() - 2);
                    nTypes[0] += (checkFace.size() - 2);

                    break;
                }

                case 2:
                {
                    if (checkFace.size() == 3)
                    {
                        nCells_++;
                        nTypes[0]++;
                    }
                    else
                    {
                        // Add a central point for the face
                        labelPair& fPair = faceMarker[fIndex];

                        if (fPair.first() == -1)
                        {
                            fPair.first() = 1;

                            nFacePoints++;
                            nDecompPoints_++;
                        }

                        nCells_ += checkFace.size();
                        nTypes[0] += checkFace.size();
                    }

                    break;
                }

                default:
                {
                    FatalErrorIn("void mesquiteMotionSolver::initArrays()")
                        << "Unknown polyhedron decomposition type."
                        << abort(FatalError);

                    break;
                }
            }
        }
    }

    if (debug)
    {
        Pout<< "nTet       = " << nTypes[0] << nl
            << "nHex       = " << nTypes[1] << nl
            << "nPyr       = " << nTypes[2] << nl
            << "nPrism     = " << nTypes[3] << nl
            << "nPolyhedra = " << nPolyhedra << endl;
    }

    label nCellToNode =
    (
        (4 * nTypes[0])
      + (8 * nTypes[1])
      + (5 * nTypes[2])
      + (6 * nTypes[3])
    );

    // Prepare arrays for mesquite
    vtxCoords_.setSize(3 * (nPoints_ + nDecompPoints_));
    cellToNode_.setSize(nCellToNode);
    fixFlags_.setSize(nPoints_ + nDecompPoints_);
    mixedTypes_.setSize(nCells_);
    decompCellCentres_.setSize(nCellPoints);
    decompFaceCentres_.setSize(nFacePoints);

    // Reset faceMarker
    if (nFacePoints)
    {
        forAll(faceMarker, faceI)
        {
            faceMarker[faceI].first() = -1;
        }
    }

    // Reset counters
    nDecompPoints_ = 0;
    nCellPoints = 0, nFacePoints = 0;

    // Set connectivity information
    label cIndex = 0, cellIndex = 0;

    forAll(meshCells, cellI)
    {
        bool foundMatch = false;

        forAll(matcher, indexI)
        {
            if
            (
                matcher[indexI].valid() &&
                matcher[indexI]->isA(mesh(), cellI)
            )
            {
                // Match shape
                matcher[indexI]->matchShape
                (
                    false,
                    meshFaces,
                    faceOwner,
                    cellI,
                    meshCells[cellI]
                );

                // Set cellToNode
                const labelList& cellToNode = matcher[indexI]->vertLabels();

                forAll(cellToNode, nodeI)
                {
                    cellToNode_[cIndex++] = cellToNode[nodeI];
                }

                // Set element type
                mixedTypes_[cellIndex++] = cellType[indexI];

                foundMatch = true;
                break;
            }
        }

        if (foundMatch)
        {
            continue;
        }

        // Decompose into tets.
        const cell& checkCell = meshCells[cellI];

        // Add a central point for the cell
        label cellPoint = nPoints_ + nDecompPoints_;

        decompCellCentres_[nCellPoints].first() = cellI;
        decompCellCentres_[nCellPoints].second() = cellPoint;

        nCellPoints++;
        nDecompPoints_++;

        switch (decompType_)
        {
            case 1:
            {
                forAll(checkCell, faceI)
                {
                    const label fIndex = checkCell[faceI];
                    const face& checkFace = meshFaces[fIndex];
                    const label& checkOwner = faceOwner[fIndex];

                    const label nP = checkFace.size();

                    if (checkOwner == cellI)
                    {
                        for (label nI = (nP - 1); nI > 1; nI--)
                        {
                            cellToNode_[cIndex++] = checkFace[0];
                            cellToNode_[cIndex++] = checkFace[nI];
                            cellToNode_[cIndex++] = checkFace[nI - 1];
                            cellToNode_[cIndex++] = cellPoint;

                            mixedTypes_[cellIndex++] = cellType[0];
                        }
                    }
                    else
                    {
                        for (label nI = 1; nI < (nP - 1); nI++)
                        {
                            cellToNode_[cIndex++] = checkFace[0];
                            cellToNode_[cIndex++] = checkFace[nI];
                            cellToNode_[cIndex++] = checkFace[nI + 1];
                            cellToNode_[cIndex++] = cellPoint;

                            mixedTypes_[cellIndex++] = cellType[0];
                        }
                    }
                }

                break;
            }

            case 2:
            {
                forAll(checkCell, faceI)
                {
                    const label fIndex = checkCell[faceI];
                    const face& checkFace = meshFaces[fIndex];
                    const label& checkOwner = faceOwner[fIndex];

                    const label nP = checkFace.size();

                    if (nP == 3)
                    {
                        if (checkOwner == cellI)
                        {
                            cellToNode_[cIndex++] = checkFace[2];
                            cellToNode_[cIndex++] = checkFace[1];
                            cellToNode_[cIndex++] = checkFace[0];
                            cellToNode_[cIndex++] = cellPoint;

                            mixedTypes_[cellIndex++] = cellType[0];
                        }
                        else
                        {
                            cellToNode_[cIndex++] = checkFace[0];
                            cellToNode_[cIndex++] = checkFace[1];
                            cellToNode_[cIndex++] = checkFace[2];
                            cellToNode_[cIndex++] = cellPoint;

                            mixedTypes_[cellIndex++] = cellType[0];
                        }
                    }
                    else
                    {
                        // Add a central point for the face
                        labelPair& fPair = faceMarker[fIndex];

                        if (fPair.first() == -1)
                        {
                            fPair.first() = fIndex;
                            fPair.second() = nPoints_ + nDecompPoints_;

                            nFacePoints++;
                            nDecompPoints_++;
                        }

                        if (checkOwner == cellI)
                        {
                            for (label pI = nP - 1, pJ = 0; pJ < nP; pI = pJ++)
                            {
                                cellToNode_[cIndex++] = checkFace[pJ];
                                cellToNode_[cIndex++] = checkFace[pI];
                                cellToNode_[cIndex++] = fPair.second();
                                cellToNode_[cIndex++] = cellPoint;

                                mixedTypes_[cellIndex++] = cellType[0];
                            }
                        }
                        else
                        {
                            for (label pI = nP - 1, pJ = 0; pJ < nP; pI = pJ++)
                            {
                                cellToNode_[cIndex++] = checkFace[pI];
                                cellToNode_[cIndex++] = checkFace[pJ];
                                cellToNode_[cIndex++] = fPair.second();
                                cellToNode_[cIndex++] = cellPoint;

                                mixedTypes_[cellIndex++] = cellType[0];
                            }
                        }
                    }
                }

                break;
            }
        }
    }

    // Condense the faceMarker
    if (nFacePoints)
    {
        nFacePoints = 0;

        forAll(faceMarker, faceI)
        {
            const labelPair& fPair = faceMarker[faceI];

            if (fPair.first() > -1)
            {
                decompFaceCentres_[nFacePoints++] = fPair;
            }
        }
    }

    faceMarker.clear();

    // Fix patch information, but blank out first
    forAll(fixFlags_, pointI)
    {
        fixFlags_[pointI] = 0;
    }

    forAll(boundary, patchI)
    {
        // Leave processor boundaries out.
        if (isA<processorPolyPatch>(boundary[patchI]))
        {
            continue;
        }

        const labelList& meshPointLabels = boundary[patchI].meshPoints();

        forAll(meshPointLabels, pointI)
        {
            fixFlags_[meshPointLabels[pointI]] = 1;
        }
    }

    // Initialise parallel connectivity, if necessary
    initParallelConnectivity();

    // Optionally fix additional zone points
    initFixedZones();

    // Optionally fix additional boundary layer zone points
    initBoundaryLayerZones();

    // Fix decomposed centroids
    fixDecomposedCentroids();

    // Clear demand-driven addressing
    clearDemandDrivenAddressing();

    // Set the flag
    arraysInitialized_ = true;
}


// Identify coupled patches
void mesquiteMotionSolver::identifyCoupledPatches()
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (debug)
    {
        Info<< "Identifying coupled patches" << endl;
    }

    // Maintain a separate list of processor IDs in procIndices.
    // This is done because this sub-domain may talk to processors
    // that share only edges/points.
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    // Fetch the list of global points from polyMesh.
    const globalMeshData& gData = mesh().globalData();

    const labelList& spA = gData.sharedPointAddr();
    const labelList& spL = gData.sharedPointLabels();

    labelListList spB(Pstream::nProcs(), labelList(0));

    if (gData.nGlobalPoints())
    {
        if (debug)
        {
            Info<< " mesquiteMotionSolver::identifyCoupledPatches :"
                << " Found " << gData.nGlobalPoints()
                << " global points." << endl;
        }

        // Send others my addressing.
        for (label proc = 0; proc < Pstream::nProcs(); proc++)
        {
            if (proc != Pstream::myProcNo())
            {
                // Send number of entities first.
                parWrite(proc, spA.size());

                // Send the buffer.
                if (spA.size())
                {
                    parWrite(proc, spA);
                }
            }
        }

        // Receive addressing from others
        for (label proc = 0; proc < Pstream::nProcs(); proc++)
        {
            if (proc != Pstream::myProcNo())
            {
                label procInfoSize = -1;

                // How many entities am I going to be receiving?
                parRead(proc, procInfoSize);

                if (procInfoSize)
                {
                    // Size the receive buffer.
                    spB[proc].setSize(procInfoSize, -1);

                    // Schedule for receipt.
                    parRead(proc, spB[proc]);
                }
            }
        }
    }
    else
    if (debug)
    {
        Info<< " mesquiteMotionSolver::identifyCoupledPatches :"
            << " Did not find any global points." << endl;
    }

    labelHashSet immNeighbours;
    labelListList procPoints(Pstream::nProcs());

    // Track globally shared points
    List<DynamicList<labelPair> > globalProcPoints
    (
        Pstream::nProcs(),
        DynamicList<labelPair>(5)
    );

    // Insert my immediate neighbours into the list.
    forAll(boundary, pI)
    {
        if (isA<processorPolyPatch>(boundary[pI]))
        {
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[pI])
            );

            label neiProcNo = pp.neighbProcNo();

            // Insert all boundary points.
            procPoints[neiProcNo] = pp.meshPoints();

            // Keep track of immediate neighbours.
            immNeighbours.insert(neiProcNo);
        }
    }

    if (gData.nGlobalPoints())
    {
        // Wait for all transfers to complete.
        OPstream::waitRequests();
        IPstream::waitRequests();

        // Now loop through all processor addressing, and check if
        // any labels coincide with my global shared points.
        // If this is true, we need to be talking to that neighbour
        // as well (if not already).
        for (label proc = 0; proc < Pstream::nProcs(); proc++)
        {
            if (proc == Pstream::myProcNo())
            {
                continue;
            }

            bool foundGlobalMatch = false;

            // Fetch reference to buffer
            const labelList& procBuffer = spB[proc];
            DynamicList<labelPair>& gPP = globalProcPoints[proc];

            forAll(procBuffer, pointI)
            {
                forAll(spA, pointJ)
                {
                    if (spA[pointJ] == procBuffer[pointI])
                    {
                        // Make an entry, if one wasn't made already
                        if (findIndex(procPoints[proc], spL[pointJ]) == -1)
                        {
                            gPP.append(labelPair(spL[pointJ], spA[pointJ]));
                        }

                        foundGlobalMatch = true;

                        break;
                    }
                }
            }

            // Sort in ascending order of global point labels
            if (gPP.size())
            {
                SortableList<label> spALabels(gPP.size(), -1);

                forAll(gPP, pointI)
                {
                    spALabels[pointI] = gPP[pointI].second();
                }

                spALabels.sort();

                // Reorder list and transfer
                List<labelPair> sortedGPP(gPP.size());

                const labelList& indices = spALabels.indices();

                forAll(sortedGPP, pointI)
                {
                    sortedGPP[pointI] = gPP[indices[pointI]];
                }

                gPP.transfer(sortedGPP);
            }

            if (!immNeighbours.found(proc))
            {
                if (foundGlobalMatch && debug)
                {
                    Pout<< " mesquiteMotionSolver::identifyCoupledPatches :"
                        << " Additionally talking to processor: "
                        << proc << endl;
                }
            }
        }
    }

    // Estimate an initial size
    procIndices_.setSize(Pstream::nProcs());

    // Patch sub meshes need to be prepared in ascending
    // order of neighbouring processors.
    label nTotalProcs = 0;

    forAll(procPoints, procI)
    {
        if (procPoints[procI].size())
        {
            procIndices_[nTotalProcs++] = procI;
        }

        // Check for point / edge coupling
        if (globalProcPoints[procI].size() && !procPoints[procI].size())
        {
            procIndices_[nTotalProcs++] = procI;
        }
    }

    // Shrink to actual size
    procIndices_.setSize(nTotalProcs);

    // Transfer to shortened lists
    procPointLabels_.setSize(nTotalProcs, labelList(0));
    globalProcPoints_.setSize(nTotalProcs, List<labelPair>(0));

    forAll(procIndices_, pI)
    {
        procPointLabels_[pI].transfer(procPoints[procIndices_[pI]]);
        globalProcPoints_[pI].transfer(globalProcPoints[procIndices_[pI]]);
    }
}


// Prepare for parallel surface smoothing
void mesquiteMotionSolver::initParallelSurfaceSmoothing()
{
    if (!surfaceSmoothing_)
    {
        return;
    }

    if (debug)
    {
        Info<< "Initializing parallel smoothing connectivity" << endl;
    }

    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    // Maintain a list of points for geometric matching
    List<DynamicList<point> > auxSurfPoints
    (
        procIndices_.size(),
        DynamicList<point>(20)
    );

    // Prepare a new Map for all processors
    sendSurfFields_.setSize(procIndices_.size(), vectorField(0));
    recvSurfFields_.setSize(procIndices_.size(), vectorField(0));
    sendSurfPointMap_.setSize(procIndices_.size(), Map<label>());
    recvSurfPointMap_.setSize(procIndices_.size(), Map<label>());
    unmatchedRecvPoints_.setSize(procIndices_.size(), Map<label>());

    // Build a reverse processor map for convenience
    labelList procMap(Pstream::nProcs(), -1);

    forAll(procIndices_, pI)
    {
        procMap[procIndices_[pI]] = pI;
    }

    forAll(pIDs_, patchI)
    {
        const label pOffset = offsets_[patchI];
        const pointField& points = boundary[pIDs_[patchI]].localPoints();

        forAll(boundary, patchJ)
        {
            // Determine processor id
            label pI = -1;

            if (isA<processorPolyPatch>(boundary[patchJ]))
            {
                const processorPolyPatch& pp =
                (
                    refCast<const processorPolyPatch>(boundary[patchJ])
                );

                // Fetch index from map
                pI = procMap[pp.neighbProcNo()];
            }
            else
            {
                // Disregard non-processor patches
                continue;
            }

            // Fetch references
            Map<label>& sspMap = sendSurfPointMap_[pI];
            Map<label>& rspMap = recvSurfPointMap_[pI];
            DynamicList<point>& asPoints = auxSurfPoints[pI];

            // Fetch mesh points for this patch
            const labelList& mp = boundary[patchJ].meshPoints();

            forAll(mp, pointI)
            {
                // Find local index in patch
                label local = boundary[pIDs_[patchI]].whichPoint(mp[pointI]);

                if (local == -1)
                {
                    continue;
                }

                // Locate index in field
                label fieldIndex = local + pOffset;

                // Fix boundary condition for this point.
                bdy_[fieldIndex] = vector::one;

                // Modify point marker
                if (procIndices_[pI] < Pstream::myProcNo())
                {
                    pointMarker_[fieldIndex] = 0.0;
                }

                // Add an entry if it wasn't done before
                if (!sspMap.found(fieldIndex))
                {
                    sspMap.insert(fieldIndex, asPoints.size());
                    rspMap.insert(fieldIndex, asPoints.size());

                    // Insert corresponding point
                    asPoints.append(points[local]);
                }
            }
        }

        // Loop through global-points information,
        // and check if additional points need to be added
        forAll(globalProcPoints_, pI)
        {
            // Fetch references
            Map<label>& sspMap = sendSurfPointMap_[pI];
            Map<label>& rspMap = recvSurfPointMap_[pI];
            DynamicList<point>& asPoints = auxSurfPoints[pI];

            const List<labelPair>& pLabels = globalProcPoints_[pI];

            forAll(pLabels, pointI)
            {
                // Find local index in patch
                label local =
                (
                    boundary[pIDs_[patchI]].whichPoint
                    (
                        pLabels[pointI].first()
                    )
                );

                if (local == -1)
                {
                    continue;
                }

                // Locate index in field
                label fieldIndex = local + pOffset;

                // Fix boundary condition for this point.
                bdy_[fieldIndex] = vector::one;

                // Modify point marker
                if (procIndices_[pI] < Pstream::myProcNo())
                {
                    pointMarker_[fieldIndex] = 0.0;
                }

                // Add an entry if it wasn't done before
                if (!sspMap.found(fieldIndex))
                {
                    if (debug)
                    {
                        Pout<< " Found local: " << local
                            << " for global: " << pLabels[pointI]
                            << endl;
                    }

                    sspMap.insert(fieldIndex, asPoints.size());
                    rspMap.insert(fieldIndex, asPoints.size());

                    // Insert corresponding point
                    asPoints.append(points[local]);
                }
            }
        }
    }

    // Check to ensure that field sizes match
    labelList nProcSize(procIndices_.size(), -1);

    forAll(procIndices_, pI)
    {
        // Send size to neighbour
        parWrite(procIndices_[pI], sendSurfPointMap_[pI].size());

        // Receive size from neighbour
        parRead(procIndices_[pI], nProcSize[pI]);

        // Require a sync for size mismatch
        bool requiresSync = false;

        if (nProcSize[pI] != sendSurfPointMap_[pI].size())
        {
            requiresSync = true;
        }

        if (debug)
        {
            Pout<< " neiProcNo: " << procIndices_[pI]
                << " mySize: " << sendSurfPointMap_[pI].size()
                << " neiSize: " << nProcSize[pI]
                << " requireSync: " << (Switch(requiresSync)).asText()
                << endl;
        }
    }

    // Prepare boundary markers for correction
    const labelListList& eFaces = mesh().edgeFaces();
    vectorField bdyMarker(pointMarker_.size(), vector::zero);

    // Prepare edgeMarkers
    forAll(pIDs_, patchI)
    {
        const label pOffset = offsets_[patchI];
        const edgeList& edges = boundary[pIDs_[patchI]].edges();
        const labelList& meshEdges = boundary[pIDs_[patchI]].meshEdges();

        for
        (
            label i = boundary[pIDs_[patchI]].nInternalEdges();
            i < edges.size();
            i++
        )
        {
            // If both points on this edge are marked,
            // this edge needs to be left out.
            label i0 = edges[i][0] + pOffset;
            label i1 = edges[i][1] + pOffset;

            bool p0 = (pointMarker_[i0] < 0.5);
            bool p1 = (pointMarker_[i1] < 0.5);

            if (p0 && p1)
            {
                edgeMarker_[patchI][i] = 0.0;
            }

            bool noProc = true;

            const labelList& eF = eFaces[meshEdges[i]];

            forAll(eF, fI)
            {
                label curPatchID = boundary.whichPatch(eF[fI]);

                // Disregard internal faces
                if (curPatchID == -1)
                {
                    continue;
                }

                if (isA<processorPolyPatch>(boundary[curPatchID]))
                {
                    // Disregard processor patches
                    noProc = false;
                    break;
                }
            }

            // Mark boundary condition vector for this edge
            if (noProc)
            {
                bdyMarker[i0] = vector::one;
                bdyMarker[i1] = vector::one;
            }
        }
    }

    // Configure a magnitude for distance matching
    const boundBox& box = mesh().bounds();

    const vector outsideBoxMax = box.max() + vector::one;

    // Send and receive points for geometric match
    forAll(procIndices_, pI)
    {
        label proc = procIndices_[pI];

        // Set size to max, in case of mismatch
        label maxSize = Foam::max(sendSurfPointMap_[pI].size(), nProcSize[pI]);

        // Fetch references
        const DynamicList<point>& bufPoints = auxSurfPoints[pI];

        vectorField& sendField = sendSurfFields_[pI];
        vectorField& recvField = recvSurfFields_[pI];

        // Size-up send / recv fields
        sendField.setSize(maxSize, vector::zero);
        recvField.setSize(maxSize, vector::zero);

        // Send and receive points
        if (recvField.size())
        {
            // Get field from neighbour
            parRead(proc, recvField);
        }

        if (sendField.size())
        {
            // Copy auxiliary points.
            // If there is a mismatch, send invalid values
            // so that geometric matches don't result in
            // false positives. This can occur in situations
            // where points on a neighbouring processor lie
            // on a slip patch, but this processor touches
            // the patch only via points or edges.
            forAll(bufPoints, pointI)
            {
                sendField[pointI] = bufPoints[pointI];
            }

            for (label i = bufPoints.size(); i < sendField.size(); ++i)
            {
                sendField[i] = outsideBoxMax;
            }

            // Send points to neighbour
            parWrite(proc, sendField);
        }
    }

    // Fetch relative tolerance
    scalar relTol = 1e-4;

    // Compute tolerance
    scalar tol = relTol * box.mag();

    // Wait for all transfers to complete.
    OPstream::waitRequests();
    IPstream::waitRequests();

    // Set up a reversePointMap for mismatches
    labelListList reversePointMap(procIndices_.size());

    forAll(procIndices_, pI)
    {
        // Fetch references
        vectorField& recvField = recvSurfFields_[pI];
        DynamicList<point>& bufPoints = auxSurfPoints[pI];

        // Match points geometrically
        labelList pointMap;

        matchPoints
        (
            bufPoints,
            recvField,
            List<scalar>(bufPoints.size(), tol),
            false,
            pointMap
        );

        // Alias for convenience
        labelList& reverseMap = reversePointMap[pI];

        // Size-up the reverse-map
        reverseMap.setSize(nProcSize[pI], -1);

        // Fill reversePointMap
        forAll(pointMap, pointI)
        {
            label mapIndex = pointMap[pointI];

            if (mapIndex < 0)
            {
                if (debug)
                {
                    Pout<< " Could not match sent point: " << pointI
                        << " :: " << bufPoints[pointI]
                        << " for proc: " << procIndices_[pI]
                        << endl;
                }

                continue;
            }

            reverseMap[mapIndex] = pointI;
        }

        // Check reversePointMap for unfilled points
        forAll(reverseMap, pointI)
        {
            // Compute tolSqr for matching
            scalar tolSqr = sqr(tol);

            if (reverseMap[pointI] == -1)
            {
                // Search my points for the unmatched point
                const vector& uPoint = recvField[pointI];
                const pointField& points = mesh().points();

                bool foundMatch = false;

                forAll(points, pointJ)
                {
                    scalar distSqr = magSqr(uPoint - points[pointJ]);

                    if (distSqr < tolSqr)
                    {
                        unmatchedRecvPoints_[pI].insert(pointI, pointJ);

                        foundMatch = true;
                        break;
                    }
                }

                if (!foundMatch)
                {
                    Pout<< " Could not match recv point: " << pointI
                        << " :: " << recvField[pointI]
                        << " for proc: " << procIndices_[pI]
                        << abort(FatalError);
                }
            }
        }

        // Renumber to shuffled indices
        Map<label>& rspMap = recvSurfPointMap_[pI];

        forAllIter(Map<label>, rspMap, pIter)
        {
            // Skip mismatches
            if (pointMap[pIter()] < 0)
            {
                rspMap.erase(pIter);
                continue;
            }

            pIter() = pointMap[pIter()];
        }
    }

    // Reset buffers
    forAll(procIndices_, pI)
    {
        sendSurfFields_[pI] = vector::zero;
        recvSurfFields_[pI] = vector::zero;
    }

    // Fix boundary conditions
    transferBuffers(bdyMarker);

    forAll(bdyMarker, indexI)
    {
        if (bdyMarker[indexI] > vector(0.5, 0.5, 0.5))
        {
            bdy_[indexI] = vector::zero;
        }
    }
}


// Prepare mesquite connectivity for parallel runs
void mesquiteMotionSolver::initMesquiteParallelArrays()
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (debug)
    {
        Info<< "Initializing parallel connectivity for Mesquite" << endl;
    }

    typedef DynamicList<label> DynamicLabelList;

    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    Map<label> nPrc;
    List<DynamicLabelList> procDecompFaces(boundary.size());

    // Build a list of nearest neighbours.
    forAll(boundary, patchI)
    {
        if (isA<processorPolyPatch>(boundary[patchI]))
        {
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[patchI])
            );

            nPrc.insert(pp.neighbProcNo(), patchI);
            procDecompFaces[patchI].reserve(pp.size());
        }
    }

    // Distribute face centres to respective processor patches
    forAll(decompFaceCentres_, indexI)
    {
        const label fIndex = decompFaceCentres_[indexI].first();
        const label patchI = boundary.whichPatch(fIndex);

        if (patchI < 0 || !isA<processorPolyPatch>(boundary[patchI]))
        {
            continue;
        }

        procDecompFaces[patchI].append(indexI);
    }

    // Auxiliary cellToNode addressing sizes
    labelList nAuxCellToNode(procIndices_.size(), -1);

    // Size up maps and buffers
    nAuxPoints_.setSize(procIndices_.size(), -1);
    nAuxCells_.setSize(procIndices_.size(), -1);

    sendPointMap_.setSize(procIndices_.size());
    recvPointMap_.setSize(procIndices_.size());

    sendCellToNode_.setSize(procIndices_.size());
    recvCellToNode_.setSize(procIndices_.size());
    sendPointBuffer_.setSize(procIndices_.size());
    recvPointBuffer_.setSize(procIndices_.size());

    // Prepare cellToNode / cellTypes for transfer
    labelList revCellMap;
    Map<label> fwdCellMap;
    labelList nSharedPoints(procIndices_.size(), 0);
    labelListList sendFixFlags(procIndices_.size());
    labelListList recvFixFlags(procIndices_.size());
    labelListList recvTypes(procIndices_.size());
    labelListList sendTypes(procIndices_.size());

    // Fetch demand-driven data
    const labelList& typeMap = cellTypeMap();
    const labelList& cellOffsets = cellOffset();
    const labelListList& pointCellList = pointCells();

    forAll(procIndices_, pI)
    {
        label proc = procIndices_[pI];

        const label invalidProcPatch = -2;
        const labelList& procPoints = procPointLabels_[pI];

        // Add all points / cells connected to procPoints
        label nProcPoints = 0, nProcCells = 0, procPatch = invalidProcPatch;

        // Determine addressing for all shared points.
        if (nPrc.found(proc))
        {
            // Note index for later
            procPatch = nPrc[proc];

            // Fetch addressing for this patch.
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[procPatch])
            );

            const labelList& neiPoints = pp.neighbPoints();

            forAll(procPoints, pointI)
            {
                const label pIndex = procPoints[pointI];
                const label nIndex = neiPoints[pointI];

                recvPointMap_[pI].insert(nIndex, pIndex);
                sendPointMap_[pI].insert(pIndex, nProcPoints);

                nProcPoints++;
            }
        }

        // Add global points as well
        //  - pLabels is already sorted in ascending order
        //    of global point indices, so this should match
        //    across processors
        const List<labelPair>& pLabels = globalProcPoints_[pI];

        forAll(pLabels, pointI)
        {
            const label pIndex = pLabels[pointI].first();

            recvPointMap_[pI].insert(nProcPoints, pIndex);
            sendPointMap_[pI].insert(pIndex, nProcPoints);

            nProcPoints++;
        }

        // Add decomposition face points
        if (procPatch != invalidProcPatch)
        {
            DynamicLabelList& procFaces = procDecompFaces[procPatch];

            forAll(procFaces, faceI)
            {
                const label indexI = procFaces[faceI];
                const label pIndex = decompFaceCentres_[indexI].second();

                recvPointMap_[pI].insert(nProcPoints, pIndex);
                sendPointMap_[pI].insert(pIndex, nProcPoints);

                nProcPoints++;
            }

            // Clear after use
            procFaces.clearStorage();
        }

        // Make a note of shared point size
        nSharedPoints[pI] = nProcPoints;

        // Fetch the list of shared processor points
        const labelList sharedPoints = sendPointMap_[pI].toc();

        // Estimate size for cellToNode and cell types
        label nCellToNode = 0;

        forAll(sharedPoints, pointI)
        {
            const labelList& pCells = pointCellList[sharedPoints[pointI]];

            forAll(pCells, cellI)
            {
                const label cIndex = pCells[cellI];

                if (fwdCellMap.found(cIndex))
                {
                    continue;
                }

                // Fetch cellToNode for this cell
                const label offset = cellOffsets[cIndex];
                const label cellType = mixedTypes_[cIndex];
                const label nCellNodes = typeMap[cellType];

                // Update counter
                nCellToNode += nCellNodes;

                // Update point map
                for (label nodeI = 0; nodeI < nCellNodes; nodeI++)
                {
                    const label pIndex = cellToNode_[offset + nodeI];

                    // Add discovered points to map
                    if (sendPointMap_[pI].found(pIndex))
                    {
                        continue;
                    }

                    sendPointMap_[pI].insert(pIndex, nProcPoints++);
                }

                // Add cell to list
                fwdCellMap.insert(cIndex, nProcCells++);
            }
        }

        // Invert the forward cell map
        revCellMap.setSize(nProcCells);

        forAllConstIter(Map<label>, fwdCellMap, cIter)
        {
            revCellMap[cIter()] = cIter.key();
        }

        // Populate cellToNode and cell types
        labelList& myTypes = sendTypes[pI];
        labelList& myCellToNode = sendCellToNode_[pI];

        myTypes.setSize(nProcCells);
        myCellToNode.setSize(nCellToNode);

        nCellToNode = 0;

        for (label cellI = 0; cellI < nProcCells; cellI++)
        {
            const label cIndex = revCellMap[cellI];

            // Fetch cellToNode for this cell
            const label offset = cellOffsets[cIndex];
            const label cellType = mixedTypes_[cIndex];
            const label nCellNodes = typeMap[cellType];

            for (label nodeI = 0; nodeI < nCellNodes; nodeI++)
            {
                const label pIndex = cellToNode_[offset + nodeI];

                myCellToNode[nCellToNode++] = sendPointMap_[pI][pIndex];
            }

            myTypes[cellI] = cellType;
        }

        // Transfer entity sizes
        parWrite(proc, nCellToNode);
        parWrite(proc, nProcPoints);
        parWrite(proc, nProcCells);

        // Schedule transfer of entities
        if (nProcPoints)
        {
            // Size up send buffer
            sendPointBuffer_[pI].setSize(nProcPoints, vector::zero);

            // Assign fixFlags and send to neighbour
            labelList& myFixFlags = sendFixFlags[pI];

            myFixFlags.setSize(nProcPoints, -1);

            forAllConstIter(Map<label>, sendPointMap_[pI], pIter)
            {
                myFixFlags[pIter()] = fixFlags_[pIter.key()];
            }

            parWrite(proc, myFixFlags);
        }

        if (nProcCells)
        {
            // Send cell types to neighbour
            parWrite(proc, myTypes);
            parWrite(proc, myCellToNode);
        }

        parRead(proc, nAuxCellToNode[pI]);
        parRead(proc, nAuxPoints_[pI]);
        parRead(proc, nAuxCells_[pI]);

        // Schedule transfer of entities
        if (nAuxPoints_[pI])
        {
            // Size up recv buffer
            recvPointBuffer_[pI].setSize(nAuxPoints_[pI], vector::zero);

            // Receive fix-flags from neighbour
            recvFixFlags[pI].setSize(nAuxPoints_[pI], -1);

            parRead(proc, recvFixFlags[pI]);
        }

        if (nAuxCells_[pI])
        {
            // Receive cell-types from neighbour
            recvTypes[pI].setSize(nAuxCells_[pI]);

            parRead(proc, recvTypes[pI]);

            // Receive cellToNode from neighbour
            recvCellToNode_[pI].setSize(nAuxCellToNode[pI]);

            parRead(proc, recvCellToNode_[pI]);
        }

        // Clear cellMaps for next processor
        fwdCellMap.clear();
        revCellMap.clear();
    }

    // Wait for transfers to complete
    OPstream::waitRequests();
    IPstream::waitRequests();

    // Prepare new arrays for mesquite
    const label nCellToNode = cellToNode_.size();

    // Initialize counters
    label cellIndex = nCells_;
    label cIndex = nCellToNode;
    label pIndex = nPoints_ + nDecompPoints_;

    // Count the number of extra entities
    label nExtraPoints = 0, nExtraCells = 0, nExtraCellToNode = 0;

    forAll(procIndices_, pI)
    {
        // Number of extra points excludes shared points
        nExtraPoints += (nAuxPoints_[pI] - nSharedPoints[pI]);

        nExtraCells += nAuxCells_[pI];

        nExtraCellToNode += nAuxCellToNode[pI];
    }

    List<double> vtxCoordsNew(3 * (nPoints_ + nDecompPoints_ + nExtraPoints));
    List<unsigned long> cellToNodeNew(nCellToNode + nExtraCellToNode);
    List<int> fixFlagsNew(nPoints_ + nDecompPoints_ + nExtraPoints);
    List<Mesquite::EntityTopology> mixedTypesNew(nCells_ + nExtraCells);

    // Copy existing arrays
    forAll(cellToNode_, i)
    {
        cellToNodeNew[i] = cellToNode_[i];
    }

    forAll(fixFlags_, i)
    {
        fixFlagsNew[i] = fixFlags_[i];
    }

    forAll(mixedTypes_, i)
    {
        mixedTypesNew[i] = mixedTypes_[i];
    }

    // Transfer pointers
    vtxCoords_.transfer(vtxCoordsNew);
    cellToNode_.transfer(cellToNodeNew);
    fixFlags_.transfer(fixFlagsNew);
    mixedTypes_.transfer(mixedTypesNew);

    // Now prepare maps for extra points in buffer,
    // and set the fix flag for all of them.
    forAll(procIndices_, pI)
    {
        const labelList& ctn = recvCellToNode_[pI];

        // Configure points
        label ptIndex = 0;

        forAll(ctn, pointI)
        {
            // Add discovered points to map
            const label recvIndex = ctn[pointI];

            Map<label>::const_iterator it = recvPointMap_[pI].find(recvIndex);

            if (it == recvPointMap_[pI].end())
            {
                recvPointMap_[pI].insert(recvIndex, pIndex);

                // Hold all new points as fixed.
                fixFlags_[pIndex] = 1;

                // Incremement point count
                ptIndex = pIndex++;
            }
            else
            {
                // Pick mapped value
                ptIndex = it();
            }

            // Assign cellToNode
            cellToNode_[cIndex++] = ptIndex;
        }

        // Configure neighbour fixFlags
        const labelList& neiFlags = recvFixFlags[pI];

        forAllConstIter(Map<label>, recvPointMap_[pI], pIter)
        {
            const label ptIndex = pIter();

            fixFlags_[ptIndex] = (fixFlags_[ptIndex] || neiFlags[pIter.key()]);
        }

        // Assign cell types from neighbour
        const labelList& cellTypes = recvTypes[pI];

        forAll(cellTypes, cellI)
        {
            mixedTypes_[cellIndex++] =
            (
                Mesquite::EntityTopology(cellTypes[cellI])
            );
        }
    }

    label expectedCells = (nCellToNode + nExtraCellToNode);
    label expectedPoints = (nPoints_ + nDecompPoints_ + nExtraPoints);

    if ( (cIndex != expectedCells) || (pIndex != expectedPoints) )
    {
        FatalErrorIn
        (
            "void mesquiteMotionSolver::initMesquiteParallelArrays()"
        )
            << " Wrong sizes: " << nl
            << " cIndex: " << cIndex
            << " Expected: " << expectedCells << nl
            << " pIndex: " << pIndex
            << " Expected: " << expectedPoints << nl
            << abort(FatalError);
    }

    // Override existing point / cell sizes
    nPoints_ = (nPoints_ + nExtraPoints);
    nCells_ = (nCells_ + nExtraCells);

    // Set up pointFractions
    pointFractions_.setSize(nPoints_, 1.0);

    forAll(procIndices_, pI)
    {
        const labelList& procPoints = procPointLabels_[pI];

        forAll(procPoints, pointI)
        {
            pointFractions_[procPoints[pointI]] += 1.0;
        }
    }

    // Loop through global-points information,
    // and check if additional points need to be added
    forAll(globalProcPoints_, pI)
    {
        const List<labelPair>& pLabels = globalProcPoints_[pI];

        forAll(pLabels, pointI)
        {
            pointFractions_[pLabels[pointI].first()] += 1.0;
        }
    }

    // Invert pointFractions
    pointFractions_ = (1.0 / pointFractions_);
}


// Prepare flags for boundary layer zones
void mesquiteMotionSolver::initBoundaryLayerZones()
{
    const dictionary& optionsDict = subDict("mesquiteOptions");

    if (!optionsDict.found("boundaryLayerZones"))
    {
        return;
    }

    const pointZoneMesh& pZones = mesh().pointZones();
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();
    const dictionary& bDict = optionsDict.subDict("boundaryLayerZones");

    wordList blPatch = bDict.toc();

    forAll(blPatch, wordI)
    {
        const word& patchName = blPatch[wordI];
        const label patchIndex = boundary.findPatchID(patchName);

        if (patchIndex == -1)
        {
            FatalErrorIn("void mesquiteMotionSolver::initBoundaryLayerZones()")
                << "Cannot find patch: " << patchName
                << abort(FatalError);
        }

        const wordList blZones(bDict.lookup(patchName));

        forAll(blZones, zoneI)
        {
            const label zoneID = pZones.findZoneID(blZones[zoneI]);

            if (zoneID == -1)
            {
                continue;
            }

            const pointZone& pLabels = pZones[zoneID];

            if (surfaceSmoothing_)
            {
                forAll(pIDs_, patchI)
                {
                    const label pOffset = offsets_[patchI];
                    const polyPatch& bPatch = boundary[pIDs_[patchI]];

                    forAll(pLabels, pointI)
                    {
                        // Find local index in patch
                        const label pLabel = pLabels[pointI];
                        const label local = bPatch.whichPoint(pLabel);

                        if (local == -1)
                        {
                            continue;
                        }

                        const label fieldIndex = local + pOffset;

                        // Modify boundary condition
                        bdy_[fieldIndex] = vector::zero;
                    }
                }
            }

            if (twoDMesh_)
            {
                continue;
            }

            forAll(pLabels, pointI)
            {
                fixFlags_[pLabels[pointI]] = 1;
            }
        }
    }
}


// Prepare flags for fixed zones
void mesquiteMotionSolver::initFixedZones()
{
    const dictionary& optionsDict = subDict("mesquiteOptions");

    if (!optionsDict.found("fixPointsZone"))
    {
        return;
    }

    const pointZoneMesh& pZones = mesh().pointZones();
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();
    const dictionary& bDict = optionsDict.subDict("fixPointsZone");

    wordList fixZones = bDict.toc();

    forAll(fixZones, zoneI)
    {
        label zoneID = pZones.findZoneID(fixZones[zoneI]);

        if (zoneID == -1)
        {
            continue;
        }

        const pointZone& pLabels = pZones[zoneID];

        if (surfaceSmoothing_)
        {
            forAll(pIDs_, patchI)
            {
                const label pOffset = offsets_[patchI];
                const polyPatch& bPatch = boundary[pIDs_[patchI]];

                forAll(pLabels, pointI)
                {
                    // Find local index in patch
                    const label pLabel = pLabels[pointI];
                    const label local = bPatch.whichPoint(pLabel);

                    if (local == -1)
                    {
                        continue;
                    }

                    const label fieldIndex = local + pOffset;

                    // Modify boundary condition
                    bdy_[fieldIndex] = vector::zero;
                }
            }
        }

        if (twoDMesh_)
        {
            continue;
        }

        forAll(pLabels, pointI)
        {
            fixFlags_[pLabels[pointI]] = 1;
        }
    }
}


// Fix decomposed centroids
void mesquiteMotionSolver::fixDecomposedCentroids()
{
    if (nDecompPoints_ == 0)
    {
        return;
    }

    label nFixedFaces = 0;

    const faceList& meshFaces = mesh().faces();

    // Fix points on decomposed boundary faces
    forAll(decompFaceCentres_, indexI)
    {
        const label fIndex = decompFaceCentres_[indexI].first();
        const label pIndex = decompFaceCentres_[indexI].second();

        const face& checkFace = meshFaces[fIndex];

        bool allFixed = true;

        forAll(checkFace, nI)
        {
            const label checkIndex = checkFace[nI];

            if (fixFlags_[checkIndex] == 0)
            {
                allFixed = false;
                break;
            }
        }

        if (allFixed)
        {
            nFixedFaces++;
            fixFlags_[pIndex] = 1;
        }
    }

    label nFixedCells = 0;

    // Fetch demand-driven data
    const labelList& typeMap = cellTypeMap();
    const labelList& cellOffsets = cellOffset();
    const labelListList& pointCellList = pointCells();

    // Loop through all cells connected to the decomposed cell centroid,
    // and check if all points on the cell (other than the centroid) are fixed.
    // If all of them are fixed, then fix the decomposed centroid as well.
    forAll(decompCellCentres_, indexI)
    {
        const label pIndex = decompCellCentres_[indexI].second();
        const labelList& pCells = pointCellList[pIndex];

        bool allFixed = true;

        forAll(pCells, cellI)
        {
            const label cIndex = pCells[cellI];

            // Fetch cellToNode for this cell
            const label offset = cellOffsets[cIndex];
            const label cellType = mixedTypes_[cIndex];
            const label nCellNodes = typeMap[cellType];

            for (label nodeI = 0; nodeI < nCellNodes; nodeI++)
            {
                const label checkIndex = cellToNode_[offset + nodeI];

                if (checkIndex == pIndex)
                {
                    continue;
                }

                if (fixFlags_[checkIndex] == 0)
                {
                    allFixed = false;
                    break;
                }
            }

            if (!allFixed)
            {
                break;
            }
        }

        if (allFixed)
        {
            nFixedCells++;
            fixFlags_[pIndex] = 1;
        }
    }

    if (debug)
    {
        Pout<< "void mesquiteMotionSolver::fixDecomposedCentroids() :"
            << " Fixed decomposed centroids on "
            << nFixedCells << " cells and "
            << nFixedFaces << " faces"
            << endl;
    }
}


// Compute and return cell typeMap
labelList mesquiteMotionSolver::cellTypeMap() const
{
    // Define a map for supported types
    const label types[][2] =
    {
        { Mesquite::TETRAHEDRON, 4 },
        { Mesquite::HEXAHEDRON, 8 },
        { Mesquite::PYRAMID, 5 },
        { Mesquite::PRISM, 6 }
    };

    // Allocate the type map
    label maxSize = 0;
    label nTypes(sizeof(types) / sizeof(types[0]));

    for (label i = 0; i < nTypes; i++)
    {
        maxSize = Foam::max(types[i][0], maxSize);
    }

    // Prepare the type map
    labelList typeMap(maxSize + 1);

    for (label i = 0; i < nTypes; i++)
    {
        typeMap[types[i][0]] = types[i][1];
    }

    return typeMap;
}


// Initialize pointCells addressing
void mesquiteMotionSolver::initPointCells() const
{
    if (debug)
    {
        Info<< "void mesquiteMotionSolver::initPointCells() const : "
            << "Calculating pointCells connectivity" << endl;
    }

    if (pointCells_.size() || cellOffset_.size())
    {
        FatalErrorIn
        (
            "void mesquiteMotionSolver::initPointCells() const"
        )   << "pointCells connectivity already allocated."
            << abort(FatalError);
    }

    const label nCells = nCells_;
    const label nPoints = nPoints_ + nDecompPoints_;

    // Count pointCells
    label counter = 0;
    labelList nPointCells(nPoints, 0);

    // Allocate the lists
    cellOffset_.setSize(nCells);
    pointCells_.setSize(nPoints);

    // Count cells per point
    const labelList& typeMap = cellTypeMap();

    for (label cellI = 0; cellI < nCells; ++cellI)
    {
        const label cellType = mixedTypes_[cellI];
        const label nCellNodes = typeMap[cellType];

        for (label nodeI = 0; nodeI < nCellNodes; nodeI++)
        {
            const unsigned long pIndex = cellToNode_[counter + nodeI];

            nPointCells[pIndex]++;
        }

        cellOffset_[cellI] = counter;
        counter += nCellNodes;
    }

    // Allocate for inidividual points
    forAll(pointCells_, pointI)
    {
        pointCells_[pointI].setSize(nPointCells[pointI]);
        nPointCells[pointI] = 0;
    }

    // Assign pointCells list
    counter = 0;

    for (label cellI = 0; cellI < nCells; ++cellI)
    {
        const label cellType = mixedTypes_[cellI];
        const label nCellNodes = typeMap[cellType];

        for (label nodeI = 0; nodeI < nCellNodes; nodeI++)
        {
            const unsigned long pIndex = cellToNode_[counter + nodeI];

            pointCells_[pIndex][nPointCells[pIndex]++] = cellI;
        }

        counter += nCellNodes;
    }
}


// Compute and return pointCells addressing
const labelListList& mesquiteMotionSolver::pointCells() const
{
    if (pointCells_.empty())
    {
        initPointCells();
    }

    return pointCells_;
}


// Compute and return cellOffset addressing
const labelList& mesquiteMotionSolver::cellOffset() const
{
    if (cellOffset_.empty())
    {
        initPointCells();
    }

    return cellOffset_;
}


// Clear demand-driven addressing
void mesquiteMotionSolver::clearDemandDrivenAddressing()
{
    pointCells_.clear();
    cellOffset_.clear();
}


// Dump current mesh for debugging
void mesquiteMotionSolver::writeDecomposition() const
{
    labelListList cpList(nCells_);
    vectorField mPoints(nPoints_ + nDecompPoints_);

    forAll(mPoints, pointI)
    {
        mPoints[pointI][0] = vtxCoords_[(3 * pointI) + 0];
        mPoints[pointI][1] = vtxCoords_[(3 * pointI) + 1];
        mPoints[pointI][2] = vtxCoords_[(3 * pointI) + 2];
    }

    label cIndex = 0;

    forAll(mixedTypes_, cellI)
    {
        label nnodes = -1;

        switch (mixedTypes_[cellI])
        {
            case Mesquite::TETRAHEDRON:
            {
                nnodes = 4;
                break;
            }

            default:
            {
                FatalErrorIn("void mesquiteMotionSolver::writeDecomposition()")
                    << "Unknown type."
                    << abort(FatalError);
            }
        }

        cpList[cellI].setSize(nnodes);

        forAll(cpList[cellI], nI)
        {
            cpList[cellI][nI] = cellToNode_[cIndex++];
        }
    }

    // Write out the mesh
    //    meshOps::writeVTK
    //    (
    //        mesh(),
    //        "tetDecompMesh",
    //        nPoints_ + nDecompPoints_,
    //        nCells_,
    //        cellToNode_.size(),
    //        mPoints,
    //        cpList,
    //        3,
    //        Map<label>(),
    //        Map<label>(),
    //        UList<scalar>(),
    //        fixFlags_
    //    );
}


// Private member function to construct parallel connectivity data
void mesquiteMotionSolver::initParallelConnectivity()
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (debug)
    {
        Info<< "Initializing parallel connectivity" << endl;
    }

    // First identify all processors that I talk to...
    identifyCoupledPatches();

    // Prepare for parallel surface smoothing
    initParallelSurfaceSmoothing();

    // If this is a 2D mesh, we're done
    if (twoDMesh_)
    {
        arraysInitialized_ = true;

        return;
    }

    // Initialize connectivity for mesquite
    initMesquiteParallelArrays();

    // Set the flag and return
    arraysInitialized_ = true;
}


// Compute centroids based on up-to-date point positions
void mesquiteMotionSolver::computeCentroids
(
    const pointField& p,
    pointField& fCtrs,
    pointField& cCtrs
)
{
    label nInvTets = 0;
    scalarField cVols(mesh().nCells(), 0.0);

    // Size the lists
    fCtrs.setSize(mesh().nFaces(), vector::zero);
    cCtrs.setSize(mesh().nCells(), vector::zero);

    const scalar oneThird = (1.0 / 3.0);
    const faceList& fs = mesh().faces();

    forAll(fs, faceI)
    {
        const face& f = fs[faceI];
        const label nPoints = f.size();

        if (nPoints == 3)
        {
            fCtrs[faceI] = oneThird * (p[f[0]] + p[f[1]] + p[f[2]]);
        }
        else
        {
            scalar sumA = 0.0;
            vector sumAc = vector::zero;

            point fCentre = p[f[0]];

            for (label pI = 1; pI < nPoints; pI++)
            {
                fCentre += p[f[pI]];
            }

            fCentre /= nPoints;

            for (label pI = nPoints - 1, pJ = 0; pJ < nPoints; pI = pJ++)
            {
                const point& currPoint = p[f[pI]];
                const point& nextPoint = p[f[pJ]];

                vector c = currPoint + nextPoint + fCentre;
                vector n = (nextPoint - currPoint) ^ (fCentre - currPoint);
                scalar a = mag(n);

                sumA += a;
                sumAc += a * c;
            }

            fCtrs[faceI] = oneThird * sumAc / (sumA + VSMALL);
        }
    }

    const labelList& own = mesh().faceOwner();
    const labelList& nei = mesh().faceNeighbour();

    labelField nCellFaces(mesh().nCells(), 0);
    vectorField cEst(mesh().nCells(), vector::zero);

    forAll(own, faceI)
    {
        cEst[own[faceI]] += fCtrs[faceI];
        nCellFaces[own[faceI]] += 1;
    }

    forAll(nei, faceI)
    {
        cEst[nei[faceI]] += fCtrs[faceI];
        nCellFaces[nei[faceI]] += 1;
    }

    forAll(cCtrs, cellI)
    {
        cEst[cellI] /= nCellFaces[cellI];
    }

    scalar tetVol;

    forAll(own, faceI)
    {
        const face& f = fs[faceI];
        const label nPoints = f.size();

        if (nPoints == 3)
        {
            tetPointRef tpr(p[f[2]], p[f[1]], p[f[0]], cEst[own[faceI]]);

            tetVol = tpr.mag();

            if (tetVol < 0.0)
            {
                nInvTets++;

                if (debug)
                {
                    Info<< "Warning: Inverted tetrahedron:" << nl
                        << " tetVol: " << tetVol << nl
                        << " tet: " << tpr << nl
                        << " Face: " << f << endl;
                }
            }

            // Accumulate volume-weighted tet centre
            cCtrs[own[faceI]] += tetVol * tpr.centre();

            // Accumulate tet volume
            cVols[own[faceI]] += tetVol;
        }
        else
        {
            for (label pI = nPoints - 1, pJ = 0; pJ < nPoints; pI = pJ++)
            {
                const point& currPoint = p[f[pI]];
                const point& nextPoint = p[f[pJ]];

                tetPointRef tpr
                (
                    nextPoint,
                    currPoint,
                    fCtrs[faceI],
                    cEst[own[faceI]]
                );

                tetVol = tpr.mag();

                if (tetVol < 0.0)
                {
                    nInvTets++;

                    if (debug)
                    {
                        Info<< "Warning: Inverted tetrahedron:" << nl
                            << " tetVol: " << tetVol << nl
                            << " tet: " << tpr << nl
                            << " Face: " << f << endl;
                    }
                }

                // Accumulate volume-weighted tet centre
                cCtrs[own[faceI]] += tetVol * tpr.centre();

                // Accumulate tet volume
                cVols[own[faceI]] += tetVol;
            }
        }
    }

    forAll(nei, faceI)
    {
        const face& f = fs[faceI];
        const label nPoints = f.size();

        if (nPoints == 3)
        {
            tetPointRef tpr(p[f[0]], p[f[1]], p[f[2]], cEst[nei[faceI]]);

            tetVol = tpr.mag();

            if (tetVol < 0.0)
            {
                nInvTets++;

                if (debug)
                {
                    Info<< "Warning: Inverted tetrahedron:" << nl
                        << " tetVol: " << tetVol << nl
                        << " tet: " << tpr << nl
                        << " Face: " << f << endl;
                }
            }

            // Accumulate volume-weighted tet centre
            cCtrs[nei[faceI]] += tetVol * tpr.centre();

            // Accumulate tet volume
            cVols[nei[faceI]] += tetVol;
        }
        else
        {
            for (label pI = nPoints - 1, pJ = 0; pJ < nPoints; pI = pJ++)
            {
                const point& currPoint = p[f[pI]];
                const point& nextPoint = p[f[pJ]];

                tetPointRef tpr
                (
                    currPoint,
                    nextPoint,
                    fCtrs[faceI],
                    cEst[nei[faceI]]
                );

                tetVol = tpr.mag();

                if (tetVol < 0.0)
                {
                    nInvTets++;

                    if (debug)
                    {
                        Info<< "Warning: Inverted tetrahedron:" << nl
                            << " tetVol: " << tetVol << nl
                            << " tet: " << tpr << nl
                            << " Face: " << f << endl;
                    }
                }

                // Accumulate volume-weighted tet centre
                cCtrs[nei[faceI]] += tetVol * tpr.centre();

                // Accumulate tet volume
                cVols[nei[faceI]] += tetVol;
            }
        }
    }

    if (nInvTets)
    {
        Info<< "Warning: Inverted tetrahedralisation was found:" << nl
            << " nInvTets: " << nInvTets << nl
            << " Enable debug flag for more information." << endl;
    }

    cCtrs /= cVols + VSMALL;
}


// Apply centroid positions
void mesquiteMotionSolver::applyCentroids()
{
    // Compute up-to-date centroids
    vectorField xF, xC;

    computeCentroids(refPoints_, xF, xC);

    // Copy points from decomposition
    forAll(decompCellCentres_, indexI)
    {
        label cIndex = decompCellCentres_[indexI].first();
        label pIndex = decompCellCentres_[indexI].second();

        const vector& v = xC[cIndex];

        vtxCoords_[(3 * pIndex) + 0] = v.x();
        vtxCoords_[(3 * pIndex) + 1] = v.y();
        vtxCoords_[(3 * pIndex) + 2] = v.z();
    }

    forAll(decompFaceCentres_, indexI)
    {
        label fIndex = decompFaceCentres_[indexI].first();
        label pIndex = decompFaceCentres_[indexI].second();

        const vector& v = xF[fIndex];

        vtxCoords_[(3 * pIndex) + 0] = v.x();
        vtxCoords_[(3 * pIndex) + 1] = v.y();
        vtxCoords_[(3 * pIndex) + 2] = v.z();
    }
}


// Copy auxiliary points to/from buffers
void mesquiteMotionSolver::copyAuxiliaryPoints(bool firstCopy)
{
    if (firstCopy && nDecompPoints_)
    {
        applyCentroids();
    }

    if (!Pstream::parRun())
    {
        return;
    }

    if (firstCopy)
    {
        forAll(procIndices_, pI)
        {
            label proc = procIndices_[pI];

            const Map<label>& pointMap = sendPointMap_[pI];

            // Fetch reference to send / recv point buffers
            vectorField& psField = sendPointBuffer_[pI];
            vectorField& prField = recvPointBuffer_[pI];

            label pIndex = -1;
            vector v = vector::zero;

            // Fill the send buffer
            forAllConstIter(Map<label>, pointMap, pIter)
            {
                pIndex = pIter.key();

                v.x() = vtxCoords_[(3 * pIndex) + 0];
                v.y() = vtxCoords_[(3 * pIndex) + 1];
                v.z() = vtxCoords_[(3 * pIndex) + 2];

                psField[pIter()] = v;
            }

            if (debug)
            {
                Pout<< " copyAuxiliaryPoints: Proc: " << proc
                    << " Sending: " << " : " << psField.size()
                    << " Recving: " << " : " << prField.size()
                    << endl;
            }

            // Send points to neighbour
            parWrite(proc, psField);

            // Receive points from neighbour
            parRead(proc, prField);
        }

        // Wait for all transfers to complete.
        OPstream::waitRequests();
        IPstream::waitRequests();

        forAll(procIndices_, pI)
        {
            const Map<label>& pointMap = recvPointMap_[pI];

            // Fetch reference to recv point buffer
            const vectorField& prField = recvPointBuffer_[pI];

            label pIndex = -1;
            vector v = vector::zero;

            // Copy from recv buffer
            forAllConstIter(Map<label>, pointMap, pIter)
            {
                pIndex = pIter();

                v = prField[pIter.key()];

                vtxCoords_[(3 * pIndex) + 0] = v.x();
                vtxCoords_[(3 * pIndex) + 1] = v.y();
                vtxCoords_[(3 * pIndex) + 2] = v.z();
            }
        }
    }
    else
    {
        // Perform a reverse copy by swapping send / recv buffers
        forAll(procIndices_, pI)
        {
            label proc = procIndices_[pI];

            const Map<label>& pointMap = recvPointMap_[pI];

            // Fetch reference to send / recv point buffers
            vectorField& psField = sendPointBuffer_[pI];
            vectorField& prField = recvPointBuffer_[pI];

            label pIndex = -1;

            // Copy to recv buffer
            forAllConstIter(Map<label>, pointMap, pIter)
            {
                pIndex = pIter();

                vector& v = prField[pIter.key()];

                v.x() = vtxCoords_[(3 * pIndex) + 0];
                v.y() = vtxCoords_[(3 * pIndex) + 1];
                v.z() = vtxCoords_[(3 * pIndex) + 2];
            }

            if (debug)
            {
                Pout<< " copyAuxiliaryPoints: Proc: " << proc
                    << " Sending: " << " : " << prField.size()
                    << " Recving: " << " : " << psField.size()
                    << endl;
            }

            // Send points to neighbour
            parWrite(proc, prField);

            // Receive points from neighbour
            parRead(proc, psField);
        }

        // Wait for all transfers to complete.
        OPstream::waitRequests();
        IPstream::waitRequests();

        forAll(procIndices_, pI)
        {
            const Map<label>& pointMap = sendPointMap_[pI];

            // Fetch reference to send buffer
            const vectorField& psField = sendPointBuffer_[pI];

            label pIndex = -1;
            vector v = vector::zero;

            // Count the number of shared points
            label nShared =
            (
                procPointLabels_[pI].size()
              + globalProcPoints_[pI].size()
            );

            forAllConstIter(Map<label>, pointMap, pIter)
            {
                pIndex = pIter.key();

                // Only update shared points
                if (pIter() < nShared)
                {
                    v = psField[pIter()];

                    // Accumulate positions
                    vtxCoords_[(3 * pIndex) + 0] += v.x();
                    vtxCoords_[(3 * pIndex) + 1] += v.y();
                    vtxCoords_[(3 * pIndex) + 2] += v.z();
                }
            }
        }

        // Now loop through all points
        // and average point positions
        forAll(pointFractions_, pointI)
        {
            scalar fract = pointFractions_[pointI];

            vtxCoords_[(3 * pointI) + 0] *= fract;
            vtxCoords_[(3 * pointI) + 1] *= fract;
            vtxCoords_[(3 * pointI) + 2] *= fract;
        }
    }
}


// Sparse matrix-vector multiply [3D]
void mesquiteMotionSolver::A
(
    const vectorField& p,
    vectorField& w
)
{
    w = vector::zero;

    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    // Gradient (n2e)
    forAll(pIDs_, patchI)
    {
        const label pOffset = offsets_[patchI];
        const edgeList& edges = boundary[pIDs_[patchI]].edges();

        forAll(edges, edgeI)
        {
            gradEdge_[patchI][edgeI] =
            (
                p[edges[edgeI][1] + pOffset]
              - p[edges[edgeI][0] + pOffset]
            );
        }
    }

    // Divergence (e2n)
    forAll(pIDs_, patchI)
    {
        const label pOffset = offsets_[patchI];
        const edgeList& edges = boundary[pIDs_[patchI]].edges();

        forAll(edges, edgeI)
        {
            gradEdge_[patchI][edgeI] *= edgeMarker_[patchI][edgeI];

            // Account for edge constant
            gradEdge_[patchI][edgeI] *= edgeConstant_[patchI][edgeI];

            w[edges[edgeI][0] + pOffset] += gradEdge_[patchI][edgeI];
            w[edges[edgeI][1] + pOffset] -= gradEdge_[patchI][edgeI];
        }
    }

    // Transfer buffers after divergence compute.
    transferBuffers(w);

    // Apply boundary conditions
    applyBCs(w);
}


// Transfer buffers for surface point fields
void mesquiteMotionSolver::transferBuffers
(
    vectorField& field,
    bool fix
)
{
    if (!Pstream::parRun())
    {
        return;
    }

    forAll(procIndices_, pI)
    {
        label neiProcNo = procIndices_[pI];
        const Map<label>& pointMap = recvSurfPointMap_[pI];

        // Fetch references
        vectorField& recvField = recvSurfFields_[pI];
        vectorField& sendField = sendSurfFields_[pI];

        if (recvField.size())
        {
            // Schedule receipt from neighbour
            parRead(neiProcNo, recvField);
        }

        if (sendField.size())
        {
            // Prepare the send buffer
            forAllConstIter(Map<label>, pointMap, pIter)
            {
                sendField[pIter()] = field[pIter.key()];
            }

            parWrite(neiProcNo, sendField);
        }
    }

    // Wait for all transfers to complete
    OPstream::waitRequests();
    IPstream::waitRequests();

    forAll(procIndices_, pI)
    {
        // Fetch references
        const Map<label>& pointMap = sendSurfPointMap_[pI];
        const vectorField& recvField = recvSurfFields_[pI];

        if (recvField.size())
        {
            if (fix)
            {
                // Fix field values
                forAllConstIter(Map<label>, pointMap, pIter)
                {
                    field[pIter.key()] = recvField[pIter()];
                }
            }
            else
            {
                // Correct field values
                forAllConstIter(Map<label>, pointMap, pIter)
                {
                    field[pIter.key()] += recvField[pIter()];
                }
            }
        }
    }
}


// Apply boundary conditions
void mesquiteMotionSolver::applyBCs
(
    vectorField& field
)
{
    forAll(pIDs_, patchI)
    {
        const label pOffset = offsets_[patchI];

        // Apply slip conditions for internal nodes
        forAll(pNormals_[patchI], pointI)
        {
            const vector& n = pNormals_[patchI][pointI];

            field[pointI + pOffset] -=
            (
                (field[pointI + pOffset] & n)*n
            );
        }
    }

    // Component-wise multiply the field with BCs.
    field = cmptMultiply(field, bdy_);

    if (twoDMesh_)
    {
        return;
    }

    // If no boundaries were fixed, fix a few points at random
    if (bdy_.size())
    {
        if (Foam::min(bdy_) > vector(0.5, 0.5, 0.5))
        {
            Random randomizer(std::time(NULL));

            label nFix = 0; //(field.size()*5)/100;

            for (label i = 0; i < nFix; i++)
            {
                field[randomizer.integer(0, field.size() - 1)] = vector::zero;
            }
        }
    }
}


// Vector dot-product
scalar mesquiteMotionSolver::dot
(
    const vectorField& f1,
    const vectorField& f2
)
{
    scalar s = 0.0;

    forAll(f1, indexI)
    {
        s += pointMarker_[indexI]*(f1[indexI] & f2[indexI]);
    }

    // Reduce across processors
    reduce(s, sumOp<scalar>());

    return s;
}


scalar mesquiteMotionSolver::normFactor
(
    const vectorField& x,
    const vectorField& b,
    const vectorField& w,
    vectorField& tmpField
)
{
    // Obtain an average for the input vector
    scalar sumPM = 0.0;
    vector xRef = vector::zero;

    forAll(x, pointI)
    {
        sumPM += pointMarker_[pointI];
        xRef += pointMarker_[pointI] * x[pointI];
    }

    reduce(xRef, sumOp<vector>());
    reduce(sumPM, sumOp<scalar>());

    xRef /= sumPM;

    A(vectorField(x.size(), xRef), tmpField);

    tmp<vectorField> nFw = (w - tmpField);
    tmp<vectorField> nFb = (b - tmpField);

    if (debug)
    {
        Info<< " xRef: " << xRef << nl
            << " cmptSumMag(nFw): " << cmptSumMag(nFw) << nl
            << " cmptSumMag(nFb): " << cmptSumMag(nFb) << endl;
    }

    return cmptSumMag(nFw) + cmptSumMag(nFb) + 1.0e-20;
}


// Component-wise sumMag
scalar mesquiteMotionSolver::cmptSumMag
(
    const vectorField& field
)
{
    scalar cSum = 0.0, m = 0.0;

    forAll(field, i)
    {
        m = pointMarker_[i];
        cSum += m*(mag(field[i].x()) + mag(field[i].y()) + mag(field[i].z()));
    }

    // Reduce across processors
    reduce(cSum, sumOp<scalar>());

    return cSum;
}


// CG solver
label mesquiteMotionSolver::CG
(
    const vectorField& b,
    vectorField& p,
    vectorField& r,
    vectorField& w,
    vectorField& x
)
{
    // Local variables
    scalar alpha, beta, rho, rhoOld, residual, wApA;
    label maxIter = x.size(), iter = 0;

    reduce(maxIter, sumOp<label>());

    // Compute initial residual
    A(x, w);

    // Compute the normFactor, using 'r' as scratch-space
    scalar norm = normFactor(x, b, w, r);

    if (debug)
    {
        Info<< "normFactor: " << norm << endl;
    }

    r = b - w;
    p = r;
    rho = dot(r, p);

    // Obtain the normalized residual
    residual = cmptSumMag(r)/norm;

    Info<< " Initial residual: " << residual;

    while ( (iter < maxIter) && (residual > tolerance_) )
    {
        A(p, w);

        wApA = dot(p, w);

        if ((mag(wApA)/norm) < VSMALL)
        {
            Info<< " Solution singularity.";

            return iter;
        }

        alpha = rho / wApA;

        forAll(x, i)
        {
            x[i] += (alpha*p[i]);
            r[i] -= (alpha*w[i]);
        }

        rhoOld = rho;

        rho = dot(r, r);

        beta = rho / rhoOld;

        forAll(p, i)
        {
            p[i] = r[i] + (beta*p[i]);
        }

        // Update the normalized residual
        residual = cmptSumMag(r)/norm;
        iter++;
    }

    Info<< " Final residual: " << residual;

    return iter;
}


// Apply fixed-value boundary conditions, if any.
void mesquiteMotionSolver::applyFixedValuePatches()
{
    if (debug)
    {
        Info<< "Applying fixed-value patches, if any" << endl;
    }

    // Fetch the sub-dictionary
    const dictionary& optionsDict = subDict("mesquiteOptions");

    // Fetch reference to boundary
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    // Construct a pointMesh.
    const pointMesh& pMesh = pointMesh::New(Mesh_);

    // Create a temporary dimensioned field
    DimensionedField<point, pointMesh> dPointField
    (
        IOobject
        (
            "dPointField",
            Mesh_.time().timeName(),
            Mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh,
        dimLength,
        pointField(refPoints_.size(), vector::zero)
    );

    if (boundaryConditions_.valid())
    {
        // Use the pointVectorField form of boundary conditions
        boundaryConditions_().correctBoundaryConditions();

        const pointField& bP = basePoints_().internalField();
        const pointField& iF = boundaryConditions_().internalField();

        // Mark boundary mesh points
        forAll(boundary, patchI)
        {
            const polyPatch& patch = boundary[patchI];

            // Skip processor patches
            if (isA<processorPolyPatch>(patch))
            {
                continue;
            }

            const labelList& mPoints = patch.meshPoints();

            forAll(mPoints, pointI)
            {
                const label pIndex = mPoints[pointI];
                dPointField[pIndex] = vector::one;
            }
        }

        // Displace boundary mesh points
        forAll(refPoints_, pI)
        {
            point& disp = dPointField[pI];
            point bpVal = (bP[pI] + iF[pI] - refPoints_[pI]);

            disp = cmptMultiply(disp, bpVal);
        }
    }
    else
    if (optionsDict.found("fixedValuePatches"))
    {
        // Check the dictionary for entries corresponding to constant
        // fixed-displacement BCs. This is done because a 'motionU'
        // field is not used to specify such BC types.
        const dictionary& fvpDict = optionsDict.subDict("fixedValuePatches");

        // Extract a list of patch names.
        wordList fixPatches = fvpDict.toc();

        // Accumulate a set of points, so that common-points
        // are not moved twice. If an overlap exists, the
        // last entry is used.
        forAll(fixPatches, wordI)
        {
            label patchI = boundary.findPatchID(fixPatches[wordI]);

            if (patchI == -1)
            {
                FatalErrorIn
                (
                    "void mesquiteMotionSolver::applyFixedValuePatches()"
                )
                    << "Cannot find patch: " << fixPatches[wordI]
                    << abort(FatalError);
            }

            // Create a patchField and evaluate.
            autoPtr<pointPatchField<point> > pField
            (
                pointPatchField<point>::New
                (
                    pMesh.boundary()[patchI],
                    dPointField,
                    fvpDict.subDict(fixPatches[wordI])
                )
            );

            pField().updateCoeffs();
        }
    }

    // Sync displacements in parallel
    if (Pstream::parRun())
    {
        vector smallVec(VSMALL, VSMALL, VSMALL);

        if (twoDMesh_)
        {
            // Size up the buffer
            vectorField parDisp(pointMarker_.size(), vector::zero);

            forAll(pIDs_, patchI)
            {
                const label pID = pIDs_[patchI];
                const label pOffset = offsets_[patchI];
                const labelList& meshPts = boundary[pID].meshPoints();

                forAll(meshPts, pointI)
                {
                    const label gIndex = meshPts[pointI];
                    parDisp[pOffset + pointI] = dPointField[gIndex];
                }
            }

            // Transfer buffers across processors
            transferBuffers(parDisp, true);

            forAll(pIDs_, patchI)
            {
                const label pID = pIDs_[patchI];
                const label pOffset = offsets_[patchI];
                const labelList& meshPts = boundary[pID].meshPoints();

                forAll(meshPts, pointI)
                {
                    const label gIndex = meshPts[pointI];
                    const vector& disp = parDisp[pOffset + pointI];

                    // Only update for non-zero displacement
                    if (cmptMag(dPointField[gIndex]) < smallVec)
                    {
                        dPointField[gIndex] = disp;
                    }
                }
            }
        }
        else
        {
            label nDomainPoints = refPoints_.size();

            forAll(procIndices_, pI)
            {
                label proc = procIndices_[pI];

                const Map<label>& pointMap = sendPointMap_[pI];

                // Fetch reference to send / recv point buffers
                vectorField& psField = sendPointBuffer_[pI];
                vectorField& prField = recvPointBuffer_[pI];

                // Fill the send buffer
                forAllConstIter(Map<label>, pointMap, pIter)
                {
                    if (pIter.key() < nDomainPoints)
                    {
                        psField[pIter()] = dPointField[pIter.key()];
                    }
                }

                // Send to neighbour
                parWrite(proc, psField);

                // Receive from neighbour
                parRead(proc, prField);
            }

            // Wait for all transfers to complete.
            OPstream::waitRequests();
            IPstream::waitRequests();

            forAll(procIndices_, pI)
            {
                const Map<label>& pointMap = recvPointMap_[pI];

                // Fetch reference to recv buffer
                const vectorField& prField = recvPointBuffer_[pI];

                forAllConstIter(Map<label>, pointMap, pIter)
                {
                    // Only update points in this domain
                    if (pIter() < nDomainPoints)
                    {
                        // Only update for non-zero displacement
                        if (cmptMag(dPointField[pIter()]) < smallVec)
                        {
                            dPointField[pIter()] = prField[pIter.key()];
                        }
                    }
                }
            }
        }
    }

    // Apply boundary layer displacement
    applyBoundaryLayerPatches(dPointField);

    // Now update refPoints with displacement values
    refPoints_ += dPointField;
}


// Apply boundary layer displacement
void mesquiteMotionSolver::applyBoundaryLayerPatches(pointField& displacement)
{
    const dictionary& optionsDict = subDict("mesquiteOptions");

    if (!optionsDict.found("boundaryLayerZones"))
    {
        return;
    }

    const pointZoneMesh& pZones = mesh().pointZones();
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();
    const dictionary& bDict = optionsDict.subDict("boundaryLayerZones");

    wordList blPatch = bDict.toc();

    forAll(blPatch, wordI)
    {
        const word& patchName = blPatch[wordI];
        const label patchIndex = boundary.findPatchID(patchName);

        // Fetch mesh points for this patch
        const labelList& mp = boundary[patchIndex].meshPoints();

        // Compute the average displacement for this patch
        vector avgDisp = vector::zero;

        forAll(mp, pointI)
        {
            avgDisp += displacement[mp[pointI]];
        }

        label patchSize = mp.size();
        reduce(patchSize, sumOp<label>());

        if (patchSize > 0)
        {
            avgDisp /= scalar(patchSize);
        }

        // Now apply average displacement to all zone points
        const wordList blZones(bDict.lookup(patchName));

        forAll(blZones, zoneI)
        {
            const label zoneID = pZones.findZoneID(blZones[zoneI]);

            if (zoneID == -1)
            {
                continue;
            }

            const pointZone& pLabels = pZones[zoneID];

            forAll(pLabels, pointI)
            {
                displacement[pLabels[pointI]] = avgDisp;
            }
        }
    }
}


// Private member function to perform Laplacian surface smoothing
void mesquiteMotionSolver::smoothSurfaces()
{
    if (debug)
    {
        Info<< "Smoothing surfaces" << endl;
    }

    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    // Copy refPoints prior to all surface-smoothing sweeps.
    origPoints_ = refPoints_;

    // Value to be used later in volume error correction
    if (volumeCorrection_)
    {
        oldVolume_ = sum(mesh().cellVolumes());
    }

    for (label i = 0; i < nSweeps_; i++)
    {
        // Prepare point-normals with updated point positions
        preparePointNormals();

        // Copy existing point-positions
        forAll(pIDs_, patchI)
        {
            const label pOffset = offsets_[patchI];
            const labelList& meshPts = boundary[pIDs_[patchI]].meshPoints();

            forAll(meshPts, pointI)
            {
                xV_[pointI + pOffset] = refPoints_[meshPts[pointI]];
            }
        }

        // Prepare edge constants
        prepareEdgeConstants(xV_);

        Info<< "Solving for point motion: ";

        label iters = CG(bV_, pV_, rV_, wV_, xV_);

        Info<< " No Iterations: " << iters << endl;

        // Update refPoints (with relaxation if necessary)
        forAll(pIDs_, patchI)
        {
            const label pOffset = offsets_[patchI];
            const labelList& meshPts = boundary[pIDs_[patchI]].meshPoints();

            forAll(meshPts, pointI)
            {
                refPoints_[meshPts[pointI]] =
                (
                    (relax_ * xV_[pointI + pOffset])
                  + ((1.0 - relax_) * origPoints_[meshPts[pointI]])
                );
            }
        }
    }

    // Check if a parallel sync is necessary
    bool requireParSync = false;

    forAll(unmatchedRecvPoints_, pI)
    {
        if (unmatchedRecvPoints_[pI].size())
        {
            requireParSync = true;
            break;
        }
    }

    // Reduce across processors.
    reduce(requireParSync, orOp<bool>());

    if (requireParSync)
    {
        // Fill buffers with current points, and transfer
        forAll(pIDs_, patchI)
        {
            const label pOffset = offsets_[patchI];
            const labelList& meshPts = boundary[pIDs_[patchI]].meshPoints();

            forAll(meshPts, pointI)
            {
                xV_[pointI + pOffset] = refPoints_[meshPts[pointI]];
            }
        }

        forAll(procIndices_, pI)
        {
            label neiProcNo = procIndices_[pI];
            const Map<label>& pointMap = sendSurfPointMap_[pI];

            // Fetch references
            vectorField& recvField = recvSurfFields_[pI];
            vectorField& sendField = sendSurfFields_[pI];

            if (recvField.size())
            {
                // Schedule receipt from neighbour
                parRead(neiProcNo, recvField);
            }

            if (sendField.size())
            {
                // Prepare the send buffer
                forAllConstIter(Map<label>, pointMap, pIter)
                {
                    sendField[pIter()] = xV_[pIter.key()];
                }

                parWrite(neiProcNo, sendField);
            }
        }

        // Wait for all transfers to complete
        OPstream::waitRequests();
        IPstream::waitRequests();

        // Correct unmatched points
        forAll(procIndices_, pI)
        {
            const vectorField& recvField = recvSurfFields_[pI];
            const Map<label>& pointMap = unmatchedRecvPoints_[pI];

            forAllConstIter(Map<label>, pointMap, pIter)
            {
                refPoints_[pIter()] = recvField[pIter.key()];
            }
        }
    }

    // Transform and update any cyclics
    forAll(boundary, patchI)
    {
        const polyPatch& patch = boundary[patchI];

        if (!isA<cyclicPolyPatch>(patch))
        {
            continue;
        }

        // Cast to cyclic
        const cyclicPolyPatch& half0 = refCast<const cyclicPolyPatch>(patch);

        // Skip the neighbor patch of the cyclic
        if (half0.neighbour())
        {
            continue;
        }

        // Fetch the other half
        const cyclicPolyPatch& half1 = half0.neighbPatch();

        bool rotate = (half0.transform() == cyclicPolyPatch::ROTATIONAL);
        bool translate = (half0.transform() == cyclicPolyPatch::TRANSLATIONAL);

        const label half0Start = half0.start();
        const label half1Start = half1.start();
        const faceList& meshFaces = mesh().faces();

        forAll(half0, faceI)
        {
            const label half0Index = (half0Start + faceI);
            const label half1Index = (half1Start + faceI);
            const face& half0Face = meshFaces[half0Index];
            const face& half1Face = meshFaces[half1Index];

            label fS = half0Face.size();

            forAll(half0Face, pointI)
            {
                const label p0Index = half0Face[pointI];
                const label p1Index = half1Face[(fS - pointI) % fS];

                if (translate)
                {
                    const vector& sepVector = half0.separationVector();

                    refPoints_[p1Index] = (refPoints_[p0Index] + sepVector);
                }
                else
                if (rotate)
                {
                    refPoints_[p1Index] = refPoints_[p0Index];

                    half0.transformPosition(refPoints_[p1Index], faceI);
                }
            }
        }
    }
}


// Find the quality of a tetrahedron.
// The function assumes points (a-b-c)
// are in counter-clockwise fashion when viewed from d.
inline scalar mesquiteMotionSolver::tetQuality
(
    const label cIndex,
    const pointField& pField,
    bool returnVolume
)
{
    const cell& cellToCheck = mesh().cells()[cIndex];

    // Disregard non-tetrahedral cells
    if (cellToCheck.size() != 4)
    {
        return 1.0;
    }

    const face& currFace = mesh().faces()[cellToCheck[0]];
    const face& nextFace = mesh().faces()[cellToCheck[1]];

    // Get the fourth point and compute cell quality
    forAll(nextFace, pointI)
    {
        if
        (
            nextFace[pointI] != currFace[0] &&
            nextFace[pointI] != currFace[1] &&
            nextFace[pointI] != currFace[2]
        )
        {
            // Compute cell-volume
            if (mesh().faceOwner()[cellToCheck[0]] == cIndex)
            {
                const point& a = pField[currFace[2]];
                const point& b = pField[currFace[1]];
                const point& c = pField[currFace[0]];
                const point& d = pField[nextFace[pointI]];

                // Obtain the magSqr edge-lengths
                scalar Le = ((b-a) & (b-a))
                          + ((c-a) & (c-a))
                          + ((d-a) & (d-a))
                          + ((c-b) & (c-b))
                          + ((d-b) & (d-b))
                          + ((d-c) & (d-c));

                scalar V = ((1.0/6.0)*(((b - a) ^ (c - a)) & (d - a)));

                if (returnVolume)
                {
                    return V;
                }
                else
                {
                    return sign(V)*((24.96100588*::cbrt(V*V))/Le);
                }
            }
            else
            {
                const point& a = pField[currFace[0]];
                const point& b = pField[currFace[1]];
                const point& c = pField[currFace[2]];
                const point& d = pField[nextFace[pointI]];

                // Obtain the magSqr edge-lengths
                scalar Le = ((b-a) & (b-a))
                          + ((c-a) & (c-a))
                          + ((d-a) & (d-a))
                          + ((c-b) & (c-b))
                          + ((d-b) & (d-b))
                          + ((d-c) & (d-c));

                scalar V = ((1.0/6.0)*(((b - a) ^ (c - a)) & (d - a)));

                if (returnVolume)
                {
                    return V;
                }
                else
                {
                    return sign(V)*((24.96100588*::cbrt(V*V))/Le);
                }
            }
        }
    }

    // Something's wrong with connectivity.
    FatalErrorIn("inline scalar mesquiteMotionSolver::tetQuality()")
        << "Cell: " << cIndex
        << " has inconsistent connectivity."
        << abort(FatalError);

    return 0.0;
}


// Private member function to check for invalid
// cells and correct if necessary.
void mesquiteMotionSolver::correctInvalidCells()
{
    // Fetch the sub-dictionary
    const dictionary& optionsDict = subDict("mesquiteOptions");

    bool checkTetValidity = true;

    if (optionsDict.found("checkTetValidity"))
    {
        checkTetValidity = readBool(optionsDict.lookup("checkTetValidity"));
    }

    // If not asked for explicitly, skip the check
    if (!checkTetValidity)
    {
        return;
    }

    // Loop through pointCells for all boundary points
    // and compute cell volume.
    const labelListList& pointCells = mesh().pointCells();
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    // Check if a minimum quality was specified.
    scalar thresh = 0.15;

    if (optionsDict.found("sliverThreshold"))
    {
        thresh = readScalar(optionsDict.lookup("sliverThreshold"));
    }

    // Obtain point-positions after smoothing
    pointField newField = refPoints_;

    DynamicList<label> invCells(50);

    forAll(pIDs_, patchI)
    {
        const labelList& meshPts = boundary[pIDs_[patchI]].meshPoints();

        forAll(meshPts, pointI)
        {
            const labelList& pCells = pointCells[meshPts[pointI]];

            forAll(pCells, cellI)
            {
                if (tetQuality(pCells[cellI], refPoints_) < thresh)
                {
                    // Add this cell to the list
                    if (findIndex(invCells, pCells[cellI]) == -1)
                    {
                        invCells.append(pCells[cellI]);
                    }
                }
            }
        }
    }

    label nInvCells = invCells.size();

    // Reduce across processors
    reduce(nInvCells, sumOp<label>());

    if (nInvCells == 0)
    {
        return;
    }

    bool valid = false;
    scalar lambda = 2.0, valFraction = 0.75;
    label nAttempts = 0;

    while (!valid)
    {
        // Assume as valid to begin with.
        valid = true;

        // Bisect the relaxation factor.
        lambda *= 0.5;

        // Update refPoints for the test points
        refPoints_ =
        (
            (lambda * newField)
          + ((1.0 - lambda) * origPoints_)
        );

        forAll(invCells, cellI)
        {
            // Compute the original value
            scalar origVal = tetQuality(invCells[cellI], origPoints_);

            // Compute the new value
            scalar newVal = tetQuality(invCells[cellI], refPoints_);

            if (newVal < (valFraction*origVal))
            {
                valid = false;

                break;
            }
        }

        // Reduce across processors
        reduce(valid, andOp<bool>());

        if (valid)
        {
            break;
        }

        nAttempts++;

        if (nAttempts > 50)
        {
            Info<< nl;

            WarningIn("void mesquiteMotionSolver::correctInvalidCells()")
                << " Failed to obtain a valid mesh." << nl
                << " Reverting to original point positions."
                << endl;

            refPoints_ = origPoints_;

            break;
        }
    }
}


//  Member function to adjust domain volume back to pre-smoothing value
//  +/- some tolerance. Uses the bisection method to identify an
//  approriate magnitude to displace all surface nodes (along point normals)
//  by some value to correct for volume loss/gain.
void mesquiteMotionSolver::correctGlobalVolume()
{
    if (!volumeCorrection_)
    {
        return;
    }

    const cellList& allCells = mesh().cells();

    // Obtain point-positions after smoothing
    pointField oldField = refPoints_;

    scalar lengthScale = mag(mesh().bounds().max());
    label iterations = 0;
    scalar magPlus, magMinus;

    // One of the initial bracket guesses will always be zero
    scalar magVal = 0;
    scalar domainVolume = 0;

    forAll(allCells, cellI)
    {
        domainVolume += tetQuality(cellI, refPoints_, true);
    }

    scalar error = (domainVolume - oldVolume_);

    // Check which way the error is going and
    // adjust left and right brackets accordingly
    if (error > 0.0)
    {
        magMinus = -1.0*lengthScale;
        magPlus = 0;
    }
    else
    {
        magPlus = lengthScale;
        magMinus = 0;
    }

    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    while (iterations < volCorrMaxIter_)
    {
        // Reset previously calculated volume
        domainVolume = 0;

        forAll(allCells, cellI)
        {
            domainVolume += tetQuality(cellI, refPoints_, true);
        }

        error = (domainVolume - oldVolume_);

        // Move one bracket side based on sign
        // of error of last iteration
        if (error > 0.0)
        {
            magPlus = magVal;
        }
        else
        {
            magMinus = magVal;
        }

        // Bring points back to old location just after smoothing
        refPoints_ = oldField;

        if (mag(error) < volCorrTolerance_)
        {
            // Final update of refPoints with optimal magnitude
            forAll(pIDs_, patchI)
            {
                const labelList& meshPts =
                (
                    boundary[pIDs_[patchI]].meshPoints()
                );

                forAll(meshPts, pointI)
                {
                    refPoints_[meshPts[pointI]] +=
                    (
                        magVal*pNormals_[patchI][pointI]
                    );
                }
            }

            domainVolume = 0;

            forAll(allCells, cellI)
            {
                domainVolume += tetQuality(cellI, refPoints_, true);
            }

            if (debug)
            {
                Info << nl
                     << "    Volume Correction Iterations: "
                     << iterations << endl;
                Info << "    Final Volume Error: "
                     << error << endl;
                Info << "    Magnitude of correction vector: "
                     << magVal << endl;
            }

            return;
        }

        // Bisection of updated guesses
        magVal = 0.5*(magPlus + magMinus);

        // Update refPoints
        forAll(pIDs_, patchI)
        {
            const labelList& meshPts = boundary[pIDs_[patchI]].meshPoints();

            forAll(meshPts, pointI)
            {
                // Move all surface points in the
                // point normal direction by magVal
                refPoints_[meshPts[pointI]] +=
                (
                    magVal*pNormals_[patchI][pointI]
                );
            }
        }

        iterations++;
    }

    // Output if loop doesn't exit by meeting error tolerance.
    WarningIn
    (
        "mesquiteSmoother::correctGlobalVolume()"
    )   << "Maximum volume correction iterations reached. "
        << endl;
}


// Enforce cylindrical constraints for slip-patches
void mesquiteMotionSolver::enforceCylindricalConstraints()
{
    if (!surfaceSmoothing_)
    {
        return;
    }

    if (debug)
    {
        Info<< "Enforcing cylindrical constraints, if any" << endl;
    }

    // Fetch the sub-dictionary
    const dictionary& optionsDict = subDict("mesquiteOptions");

    // Check for sub-dictionary entry
    if (optionsDict.found("cylindricalConstraints"))
    {
        const dictionary& constraintDict =
        (
            optionsDict.subDict("cylindricalConstraints")
        );

        const polyBoundaryMesh& boundary = mesh().boundaryMesh();

        // Read patch-information one-by one.
        wordList cstrPatches = constraintDict.toc();

        forAll(cstrPatches, wordI)
        {
            label pID = boundary.findPatchID(cstrPatches[wordI]);

            if (pID == -1 || findIndex(pIDs_, pID) == -1)
            {
                FatalErrorIn
                (
                    "void mesquiteMotionSolver::enforceCylindricalConstraints()"
                )
                    << " Cannot find patch: " << cstrPatches[wordI]
                    << " in the slipPatches sub-dictionary."
                    << abort(FatalError);
            }

            const dictionary& pD = constraintDict.subDict(cstrPatches[wordI]);

            // Read info.
            vector axisPoint(pD.lookup("axisPoint"));
            vector axisVector(pD.lookup("axisVector"));
            scalar radius = readScalar(pD.lookup("radius"));

            const labelList& meshPts = boundary[pID].meshPoints();

            axisVector /= mag(axisVector) + VSMALL;

            scalar maxViol = 0.0;
            label nViolations = 0;

            forAll(meshPts, pointI)
            {
                point& x = refPoints_[meshPts[pointI]];

                vector rx = (x - axisPoint);
                vector ra = (rx & axisVector)*axisVector;
                vector r  = (rx - ra);
                vector rn = r / mag(r);

                // Correct point position
                x -= (mag(r) - radius)*rn;

                if (debug)
                {
                    scalar viol = mag(mag(r) - radius);

                    if (viol > SMALL)
                    {
                        ++nViolations;

                        maxViol = (viol > maxViol) ? viol : maxViol;
                    }
                }
            }

            if (debug)
            {
                Pout<< " Patch: " << cstrPatches[wordI]
                    << " No. of violations: " << nViolations
                    << " Max constraint violation: " << maxViol
                    << endl;
            }
        }
    }
}


// Utility method to check validity of cells connected to a point.
bool mesquiteMotionSolver::checkValidity
(
    const vector& x,
    const labelList& jList,
    scalar& beta
)
{
    beta = 0.0;

    bool foundInvalid = false;

    for (label i = 0; i < jList.size(); i += 3)
    {
        const point& xm0 = refPoints_[jList[i+0]];
        const point& xm1 = refPoints_[jList[i+1]];
        const point& xm2 = refPoints_[jList[i+2]];

        // Prepare the Jacobian.
        tensor J
        (
            xm0.x() - x.x(), xm1.x() - x.x(), xm2.x() - x.x(),
            xm0.y() - x.y(), xm1.y() - x.y(), xm2.y() - x.y(),
            xm0.z() - x.z(), xm1.z() - x.z(), xm2.z() - x.z()
        );

        scalar alpha = det(J);

        beta += alpha;

        if (alpha < 0.0)
        {
            foundInvalid = true;
        }
    }

    // Prepare beta.
    beta = mag(beta) / (10.0 * jList.size());

    return foundInvalid;
}


// Prepare point-normals with updated point positions
void mesquiteMotionSolver::preparePointNormals()
{
    if (debug)
    {
        Info<< "Preparing point normals for surface smoothing" << endl;
    }

    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    forAll(pIDs_, patchI)
    {
        // First update localPoints with latest point positions
        const labelList& meshPts = boundary[pIDs_[patchI]].meshPoints();

        forAll(meshPts, pointI)
        {
            localPts_[patchI][pointI] = refPoints_[meshPts[pointI]];
        }

        // Now compute point normals from updated local points
        const faceList& faces = boundary[pIDs_[patchI]].localFaces();

        pNormals_[patchI] = vector::zero;

        forAll(faces, faceI)
        {
            const face& thisFace = faces[faceI];

            vector n = thisFace.normal(localPts_[patchI]);

            forAll(thisFace, pointI)
            {
                pNormals_[patchI][thisFace[pointI]] += n;
            }
        }
    }

    if (Pstream::parRun())
    {
        // Size up the buffer
        vectorField parNormals(pointMarker_.size(), vector::zero);

        // Copy point-normals to buffer
        forAll(pIDs_, patchI)
        {
            label pOffset = offsets_[patchI];
            vectorField& pN = pNormals_[patchI];

            forAll(pN, pointI)
            {
                parNormals[pOffset + pointI] = pN[pointI];
            }
        }

        // Transfer buffers across processors
        transferBuffers(parNormals);

        // Set updated normals
        forAll(pIDs_, patchI)
        {
            label pOffset = offsets_[patchI];
            vectorField& pN = pNormals_[patchI];

            forAll(pN, pointI)
            {
                pN[pointI] = parNormals[pOffset + pointI];
            }
        }
    }

    // Normalize all point-normals
    forAll(pIDs_, patchI)
    {
        pNormals_[patchI] /= mag(pNormals_[patchI]) + VSMALL;
    }
}


// Prepare non-uniform edge constants with updated point positions
void mesquiteMotionSolver::prepareEdgeConstants
(
    const vectorField& p
)
{
    // Fetch the sub-dictionary
    const dictionary& optionsDict = subDict("mesquiteOptions");

    // Check for sub-dictionary entry
    if (!optionsDict.found("nonUniformEdgeConstant"))
    {
        return;
    }

    if (debug)
    {
        Info<< "Preparing non-uniform edge constants" << endl;
    }

    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    forAll(pIDs_, patchI)
    {
        const label pOffset = offsets_[patchI];
        const edgeList& edges = boundary[pIDs_[patchI]].edges();

        scalarField& k = edgeConstant_[patchI];

        forAll(edges, edgeI)
        {
            k[edgeI] =
            (
                Foam::mag
                (
                    p[edges[edgeI][1] + pOffset]
                  - p[edges[edgeI][0] + pOffset]
                )
            );
        }
    }
}


//- Return point location obtained from the current motion field
tmp<pointField> mesquiteMotionSolver::curPoints() const
{
    tmp<pointField> tcurPoints(refPoints_);

    motionSolver::twoDCorrectPoints(tcurPoints());

    return tcurPoints;
}


//- Solve for new mesh points
void mesquiteMotionSolver::solve()
{
    // Initialize connectivity arrays for Mesquite
    if (!arraysInitialized_)
    {
        initArrays();
    }

    // Apply fixed-value motion BC's, if any.
    applyFixedValuePatches();

    // Perform surface smoothing first
    if
    (
        surfaceSmoothing_
     && (Mesh_.time().timeIndex() % surfInterval_ == 0)
     && (Mesh_.time().timeIndex() != 0)
    )
    {
        smoothSurfaces();

        if (twoDMesh_)
        {
            return;
        }

        // Check for invalid cells and correct if necessary.
        correctInvalidCells();

        // Correct global volume, if necessary
        correctGlobalVolume();

        // Enforce constraints, if necessary
        enforceCylindricalConstraints();
    }

    if (twoDMesh_)
    {
        return;
    }

    // Copy refPoints prior to internal smoothing
    origPoints_ = refPoints_;

    // Copy most recent point positions
    forAll(refPoints_, pointI)
    {
        vtxCoords_[(3 * pointI) + 0] = refPoints_[pointI][0];
        vtxCoords_[(3 * pointI) + 1] = refPoints_[pointI][1];
        vtxCoords_[(3 * pointI) + 2] = refPoints_[pointI][2];
    }

    // Copy auxiliary points
    copyAuxiliaryPoints(true);

    Mesquite::MsqError err;

    //- ArrayMesh object defined by Mesquite
    Mesquite::ArrayMesh
    msqMesh
    (
        3,                         // Number of coords per vertex
        nPoints_ + nDecompPoints_, // Number of vertices
        vtxCoords_.begin(),        // The vertex coordinates
        fixFlags_.begin(),         // Fixed vertex flags
        nCells_,                   // Number of elements
        mixedTypes_.begin(),       // Array with all Element types
        cellToNode_.begin(),       // Connectivity
        NULL,                      // Element offset connectivity
        false                      // Fortran-style array indexing
    );

    // Create an instruction queue
    Mesquite::InstructionQueue queue;

    // Apply termination criteria to the optimization algorithm
    optAlgorithm_->set_outer_termination_criterion(&tcOuter_);
    optAlgorithm_->set_inner_termination_criterion(&tcInner_);

    // Set up the quality assessor
    PtrList<Mesquite::QualityAssessor> qA(qMetrics_.size());

    forAll(qA, metricI)
    {
        qA.set
        (
            metricI,
            new Mesquite::QualityAssessor
            (
                &qMetricTable_[qMetrics_[metricI]]()
            )
        );

        // Assess the quality of the initial mesh before smoothing
        queue.add_quality_assessor(&qA[metricI], err);
    }

    // Set the master quality improver
    queue.set_master_quality_improver
    (
        &optAlgorithm_(),
        err
    );

    // Assess the quality of the final mesh after smoothing
    forAll(qA, metricI)
    {
        queue.add_quality_assessor(&qA[metricI], err);

        // Disable Mesquite output for parallel runs.
        if (Pstream::parRun())
        {
            qA[metricI].disable_printing_results();
        }
    }

    // Launches optimization on the mesh
    queue.run_instructions(&msqMesh, err);

    // Copy auxiliary points
    copyAuxiliaryPoints(false);

    // Copy updated positions
    forAll(refPoints_, pointI)
    {
        refPoints_[pointI][0] = vtxCoords_[(3 * pointI) + 0];
        refPoints_[pointI][1] = vtxCoords_[(3 * pointI) + 1];
        refPoints_[pointI][2] = vtxCoords_[(3 * pointI) + 2];

        // Relax points
        refPoints_[pointI] =
        (
            (relax_ * refPoints_[pointI])
          + ((1.0 - relax_) * origPoints_[pointI])
        );
    }
}


//- Move points
bool mesquiteMotionSolver::movePoints()
{
    return true;
}


//- Move points
void mesquiteMotionSolver::movePoints(const pointField&)
{}


//- Update topology
void mesquiteMotionSolver::updateMesh(const mapPolyMesh& mpm)
{
    update(mpm);
}


//- Update on topology change
void mesquiteMotionSolver::update(const mapPolyMesh& mpm)
{
    if (debug)
    {
        Info<< "Clearing out mesquiteMotionSolver for topo-changes" << endl;
    }

    motionSolver::updateMesh(mpm);

    if (surfaceSmoothing_)
    {
        // Clear out CG variables
        bV_.clear();
        xV_.clear();
        pV_.clear();
        rV_.clear();
        wV_.clear();
        bdy_.clear();
        pointMarker_.clear();

        localPts_.clear();
        gradEdge_.clear();
        pNormals_.clear();
        offsets_.clear();
        edgeMarker_.clear();
        edgeConstant_.clear();
    }

    nPoints_ = Mesh_.nPoints();
    nCells_ = 0;
    nDecompPoints_ = 0;

    decompCellCentres_.clear();
    decompFaceCentres_.clear();

    nAuxPoints_.clear();
    nAuxCells_.clear();

    // Reset refPoints
    refPoints_.clear();
    refPoints_ = Mesh_.points();

    // Boundary conditions will be mapped automatically
    // via the mesh registry

    // Clear the auxiliary maps / buffers
    procIndices_.clear();
    pointFractions_.clear();
    procPointLabels_.clear();
    globalProcPoints_.clear();

    unmatchedRecvPoints_.clear();
    sendSurfFields_.clear();
    recvSurfFields_.clear();
    sendSurfPointMap_.clear();
    recvSurfPointMap_.clear();

    sendPointMap_.clear();
    recvPointMap_.clear();
    sendCellToNode_.clear();
    recvCellToNode_.clear();
    sendPointBuffer_.clear();
    recvPointBuffer_.clear();

    // Clear Mesquite arrays
    vtxCoords_.clear();
    cellToNode_.clear();
    fixFlags_.clear();
    mixedTypes_.clear();
    pointCells_.clear();
    cellOffset_.clear();

    // Reset flag
    arraysInitialized_ = false;
}

} // End namespace Foam

// ************************************************************************* //
