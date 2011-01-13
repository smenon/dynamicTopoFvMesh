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
#include "globalMeshData.H"
#include "processorPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "coupleMap.H"

#include "matchPoints.H"
#include "DimensionedField.H"
#include "pointPatchField.H"
#include "pointMesh.H"

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
    MeshObject<polyMesh, mesquiteMotionSolver>(mesh),
    Mesh_(mesh),
    twoDMesh_(mesh.nGeometricD() == 2 ? true : false),
    arraysInitialized_(false),
    nPoints_(mesh.nPoints()),
    nCells_(mesh.nCells()),
    nAuxPoints_(0),
    nAuxCells_(0),
    surfaceSmoothing_(false),
    volumeCorrection_(false),
    volCorrTolerance_(1e-20),
    volCorrMaxIter_(100),
    tolerance_(1e-4),
    nSweeps_(1),
    surfInterval_(1),
    relax_(1.0),
    vtxCoords_(NULL),
    cellToNode_(NULL),
    fixFlags_(NULL),
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
    oldVolume_(0.0)
{
    // Read options from the dictionary
    readOptions();
}


mesquiteMotionSolver::mesquiteMotionSolver
(
    const polyMesh& mesh,
    Istream& msData
)
:
    motionSolver(mesh),
    MeshObject<polyMesh, mesquiteMotionSolver>(mesh),
    Mesh_(mesh),
    twoDMesh_(mesh.nGeometricD() == 2 ? true : false),
    arraysInitialized_(false),
    nPoints_(mesh.nPoints()),
    nCells_(mesh.nCells()),
    nAuxPoints_(0),
    nAuxCells_(0),
    surfaceSmoothing_(false),
    volumeCorrection_(false),
    volCorrTolerance_(1e-20),
    volCorrMaxIter_(100),
    tolerance_(1e-4),
    nSweeps_(1),
    surfInterval_(1),
    relax_(1.0),
    vtxCoords_(NULL),
    cellToNode_(NULL),
    fixFlags_(NULL),
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
    oldVolume_(0.0)
{
    // Read options from the dictionary
    readOptions();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mesquiteMotionSolver::~mesquiteMotionSolver()
{
    clearOut();
}


// Clear out addressing
void mesquiteMotionSolver::clearOut()
{
    if (debug)
    {
        Info<< "Clearing out mesquite arrays" << endl;
    }

    // Delete memory pointers
    delete [] vtxCoords_;
    delete [] cellToNode_;
    delete [] fixFlags_;

    // Reset to NULL
    vtxCoords_ = NULL;
    cellToNode_ = NULL;
    fixFlags_ = NULL;

    // Reset array flag
    arraysInitialized_ = false;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Read options from the dictionary
void mesquiteMotionSolver::readOptions()
{
    if (debug)
    {
        Info<< "Reading options for mesquiteMotionSolver" << endl;
    }

    // Fetch the sub-dictionary
    const dictionary& optionsDict = subDict("mesquiteOptions");

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

    // Read the optimization metric
    qMetric_ = word(optionsDict.lookup("optMetric"));

    if (!qMetricTable_.found(qMetric_))
    {
        FatalErrorIn("void mesquiteMotionSolver::readOptions()")
            << "Unrecognized quality metric: " << qMetric_ << nl
            << "Available metrics are: " << nl << qMetricTable_.toc()
            << abort(FatalError);
    }
    else
    {
        Info<< "Selecting quality metric: " << qMetric_ << endl;
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
                        &qMetricTable_[qMetric_]()
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
                        &qMetricTable_[qMetric_](),
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
                        &qMetricTable_[qMetric_]()
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
                        &qMetricTable_[qMetric_]()
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
                        &qMetricTable_[qMetric_]()
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
                        &qMetricTable_[qMetric_]()
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
                        &qMetricTable_[qMetric_]()
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
                        &qMetricTable_[qMetric_]()
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

        label totalSize = 0;

        forAll(pIDs_, patchI)
        {
            label nPts = boundary[pIDs_[patchI]].nPoints();
            label nEdg = boundary[pIDs_[patchI]].nEdges();

            pNormals_[patchI].setSize(nPts, vector::zero);
            localPts_[patchI].setSize(nPts, vector::zero);
            gradEdge_[patchI].setSize(nEdg, vector::zero);
            edgeMarker_[patchI].setSize(nEdg, 1.0);

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
            const edgeList& edges = boundary[pIDs_[patchI]].edges();

            for
            (
                label i = boundary[pIDs_[patchI]].nInternalEdges();
                i < edges.size();
                i++
            )
            {
                bdy_[edges[i][0] + offsets_[patchI]] = vector::zero;
                bdy_[edges[i][1] + offsets_[patchI]] = vector::zero;
            }
        }

        origPoints_.setSize(refPoints_.size(), vector::zero);
    }

    if (twoDMesh_)
    {
        // Initialise parallel connectivity, if necessary
        initParallelConnectivity();

        // Set the flag
        arraysInitialized_ = true;

        return;
    }

    // Prepare arrays for mesquite
    vtxCoords_ = new double[3 * nPoints_];
    cellToNode_ = new unsigned long[4 * nCells_];
    fixFlags_ = new int[nPoints_];

    // Set connectivity information
    label cIndex = 0;

    const faceList& meshFaces = mesh().faces();
    const cellList& meshCells = mesh().cells();
    const labelList& owner = mesh().faceOwner();

    forAll(meshCells, cellI)
    {
        const cell& curCell  = meshCells[cellI];
        const face& currFace = meshFaces[curCell[0]];
        const face& nextFace = meshFaces[curCell[1]];

        // Get the fourth point
        forAll(nextFace, pointI)
        {
            if
            (
                nextFace[pointI] != currFace[0]
             && nextFace[pointI] != currFace[1]
             && nextFace[pointI] != currFace[2]
            )
            {
                // Fill in cellPoints in order
                if (owner[curCell[0]] == cellI)
                {
                    cellToNode_[cIndex++] = currFace[2];
                    cellToNode_[cIndex++] = currFace[1];
                    cellToNode_[cIndex++] = currFace[0];
                    cellToNode_[cIndex++] = nextFace[pointI];
                }
                else
                {
                    cellToNode_[cIndex++] = currFace[0];
                    cellToNode_[cIndex++] = currFace[1];
                    cellToNode_[cIndex++] = currFace[2];
                    cellToNode_[cIndex++] = nextFace[pointI];
                }

                break;
            }
        }
    }

    // Fix patch information, but blank out first
    forAll(refPoints_, pointI)
    {
        fixFlags_[pointI] = 0;
    }

    forAll(boundary, patchI)
    {
        const labelList& meshPointLabels = boundary[patchI].meshPoints();

        // Leave processor boundaries out.
        if (boundary[patchI].type() == "processor")
        {
            continue;
        }

        forAll(meshPointLabels, pointI)
        {
            fixFlags_[meshPointLabels[pointI]] = 1;
        }
    }

    // Initialise parallel connectivity, if necessary
    initParallelConnectivity();

    // Set the flag
    arraysInitialized_ = true;
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

    Map<label> nPrc;
    label pIndex = nPoints_, cIndex = (4*nCells_);

    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

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
        }
    }

    // Search the registry for all mapping objects.
    HashTable<const coupleMap*> coupleMaps = Mesh_.lookupClass<coupleMap>();

    if (debug)
    {
        forAllIter(HashTable<const coupleMap*>, coupleMaps, cmIter)
        {
            const coupleMap& cMap = *(cmIter());

            Pout<< " Name: " << cMap.name()
                << " Send: " << cMap.isSend()
                << " Recv: " << cMap.isRecv()
                << endl;
        }
    }

    if (!twoDMesh_)
    {
        // Prepare arrays for mesquite
        forAllIter(HashTable<const coupleMap*>, coupleMaps, cmIter)
        {
            const coupleMap& cMap = *(cmIter());

            // Pick maps for which this sub-domain is master
            if
            (
                (cMap.masterIndex() == Pstream::myProcNo()) &&
                (cMap.isRecv()) && (!cMap.isLocal())
            )
            {
                nAuxPoints_ +=
                (
                    cMap.nEntities(coupleMap::POINT)
                  - cMap.nEntities(coupleMap::SHARED_POINT)
                );

                nAuxCells_ += cMap.nEntities(coupleMap::CELL);
            }
        }

        double *vtxCoordsNew = NULL;
        unsigned long *cellToNodeNew = NULL;
        int *fixFlagsNew = NULL;

        vtxCoordsNew = new double[3 * (nPoints_ + nAuxPoints_)];
        cellToNodeNew = new unsigned long[4 * (nCells_ + nAuxCells_)];
        fixFlagsNew = new int[(nPoints_ + nAuxPoints_)];

        // Copy existing arrays
        for (label i = 0; i < cIndex; i++)
        {
            cellToNodeNew[i] = cellToNode_[i];
        }

        for (label i = 0; i < pIndex; i++)
        {
            fixFlagsNew[i] = fixFlags_[i];
        }

        // Delete demand-driven data
        clearOut();

        // Transfer pointers
        vtxCoords_ = vtxCoordsNew;
        cellToNode_ = cellToNodeNew;
        fixFlags_ = fixFlagsNew;

        vtxCoordsNew = NULL;
        cellToNodeNew = NULL;
        fixFlagsNew = NULL;
    }

    // Search the registry for all mapping objects.
    forAllIter(HashTable<const coupleMap*>, coupleMaps, cmIter)
    {
        const coupleMap& cMap = *(cmIter());

        if (cMap.isLocal() || cMap.isSend())
        {
            continue;
        }

        const labelList& pBuffer = cMap.entityBuffer(coupleMap::POINT);

        label neiProcNo = -1;

        if (cMap.masterIndex() == Pstream::myProcNo())
        {
            neiProcNo = cMap.slaveIndex();
        }
        else
        if (cMap.slaveIndex() == Pstream::myProcNo())
        {
            neiProcNo = cMap.masterIndex();
        }
        else
        {
            FatalErrorIn
            (
                "void mesquiteMotionSolver::initParallelConnectivity()"
            )
                << " Wrongly configured coupleMap: " << nl
                << " My proc: " << Pstream::myProcNo() << nl
                << " masterIndex: " << cMap.masterIndex() << nl
                << " slaveIndex: " << cMap.slaveIndex() << nl
                << abort(FatalError);
        }

        // Prepare a new Map for this processor
        auxPointMap_.insert(neiProcNo, Map<label>());

        // Determine addressing for all shared points.
        if (nPrc.found(neiProcNo))
        {
            // Fetch addressing for this patch.
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[nPrc[neiProcNo]])
            );

            const labelList& neiPoints = pp.neighbPoints();
            const labelList& mPoints = boundary[nPrc[neiProcNo]].meshPoints();

            forAll(mPoints, pointI)
            {
                auxPointMap_[neiProcNo].insert
                (
                    neiPoints[pointI],
                    mPoints[pointI]
                );
            }
        }
        else
        {
            // Disconnected neighbour.
            // Look at globalPoint addressing
            // for its information.
            const globalMeshData& gD = mesh().globalData();
            const labelList& spA = gD.sharedPointAddr();
            const labelList& spL = gD.sharedPointLabels();

            forAll(pBuffer, pointI)
            {
                label addrIndex = findIndex(spA, pBuffer[pointI]);

                // Find the global index.
                auxPointMap_[neiProcNo].insert(pointI, spL[addrIndex]);
            }
        }

        if (twoDMesh_)
        {
            continue;
        }

        if (cMap.masterIndex() == Pstream::myProcNo())
        {
            // Fetch demand-driven addressing from coupleMap
            const labelList& own = cMap.owner();
            const cellList& cList = cMap.cells();
            const faceList& fList = cMap.faces();

            forAll(cList, cellI)
            {
                const cell& cellToCheck = cList[cellI];

                // Assume all tetrahedral cells
                FixedList<label, 4> p(-1);

                // Fetch the first two faces for this cell.
                const face& tFI = fList[cellToCheck[0]];
                const face& tFJ = fList[cellToCheck[1]];

                // Determine a unique point on the second face.
                forAll(tFJ, pJ)
                {
                    if
                    (
                        (tFJ[pJ] != tFI[0]) &&
                        (tFJ[pJ] != tFI[1]) &&
                        (tFJ[pJ] != tFI[2])
                    )
                    {
                        // Found all necessary points.
                        // Now check orientation.
                        if (own[cellToCheck[0]] == cellI)
                        {
                            p[0] = tFI[2];
                            p[1] = tFI[1];
                            p[2] = tFI[0];
                            p[3] = tFJ[pJ];
                        }
                        else
                        {
                            p[0] = tFI[0];
                            p[1] = tFI[1];
                            p[2] = tFI[2];
                            p[3] = tFJ[pJ];
                        }

                        break;
                    }
                }

                // Prepare maps for all new points and fill in cellToNode.
                forAll(p, i)
                {
                    if (!auxPointMap_[neiProcNo].found(p[i]))
                    {
                        auxPointMap_[neiProcNo].insert
                        (
                            p[i],
                            pIndex++
                        );

                        // Hold all new points as fixed.
                        fixFlags_[p[i]] = 1;
                    }

                    // Assign cellToNode
                    cellToNode_[cIndex++] = auxPointMap_[neiProcNo][p[i]];
                }
            }
        }
        else
        {
            // Fetch the neighbour's procIndex
            const labelList& smPoints = cMap.subMeshPoints();

            // Hold all subMeshPoints as fixed.
            forAll(smPoints, pointI)
            {
                fixFlags_[smPoints[pointI]] = 1;
            }
        }
    }

    // Initialize connectivity for surface smoothing.
    if (!surfaceSmoothing_)
    {
        return;
    }

    Map<DynamicList<point> > auxSurfPoints;
    boolList bdyMarker(pointMarker_.size(), false);

    forAllIter(HashTable<const coupleMap*>, coupleMaps, cmIter)
    {
        const coupleMap& cMap = *(cmIter());

        if (cMap.isLocal())
        {
            continue;
        }

        label neiProcNo = -1;
        bool slaveProc = false;

        // Fetch the neighbour's procIndex
        if (cMap.masterIndex() == Pstream::myProcNo())
        {
            neiProcNo = cMap.slaveIndex();
        }
        else
        if (cMap.slaveIndex() == Pstream::myProcNo())
        {
            neiProcNo = cMap.masterIndex();
            slaveProc = true;
        }

        // Prepare a new Map for this processor
        auxSurfPointMap_.insert(neiProcNo, Map<label>());
        auxNeiBufferMap_.insert(neiProcNo, Map<label>());
        auxSurfPoints.insert(neiProcNo, DynamicList<point>(20));

        // Fetch reference
        Map<label>& aspMap = auxSurfPointMap_[neiProcNo];
        Map<label>& anbMap = auxNeiBufferMap_[neiProcNo];
        DynamicList<point>& asPoints = auxSurfPoints[neiProcNo];

        // Fetch buffers from coupleMap
        const pointField& mPts = cMap.pointBuffer();
        const labelList& starts = cMap.entityBuffer(coupleMap::FACE_STARTS);
        const labelList& sizes = cMap.entityBuffer(coupleMap::FACE_SIZES);
        const label nShared = cMap.nEntities(coupleMap::SHARED_POINT);
        const faceList& fList = cMap.faces();

        forAll(pIDs_, patchI)
        {
            // Physical patches are ordered the same way
            // in coupleMaps, so pIDs_ can be used directly.
            label fStart = starts[pIDs_[patchI]];
            label fSize = sizes[pIDs_[patchI]];

            for (label faceI = fStart; faceI < fStart + fSize; faceI++)
            {
                const face& thisFace = fList[faceI];

                forAll(thisFace, pI)
                {
                    if (thisFace[pI] >= nShared)
                    {
                        continue;
                    }

                    // Fetch the global index.
                    label gIndex = auxPointMap_[neiProcNo][thisFace[pI]];
                    label lIndex = boundary[pIDs_[patchI]].whichPoint(gIndex);
                    label fieldIndex = lIndex + offsets_[patchI];

                    if (lIndex > -1)
                    {
                        // Fix boundary condition for this point.
                        bdy_[fieldIndex] = vector::one;

                        // Mark for posterity.
                        bdyMarker[fieldIndex] = true;

                        // Modify point marker
                        if (slaveProc)
                        {
                            pointMarker_[fieldIndex] = 0.0;
                        }

                        // Add a mapping entry.
                        if (!aspMap.found(fieldIndex))
                        {
                            aspMap.insert(fieldIndex, asPoints.size());
                            anbMap.insert(fieldIndex, asPoints.size());

                            // Insert corresponding point
                            asPoints.append(mPts[thisFace[pI]]);
                        }
                    }
                    else
                    {
                        // Something is wrong here.
                        Pout<< " Unable to find local index. " << nl
                            << " gIndex: " << gIndex << nl
                            << " lIndex: " << lIndex << nl
                            << abort(FatalError);
                    }
                }
            }
        }
    }

    if (debug)
    {
        Map<label> nProcSize;

        // Check to ensure that field sizes match

        // Transfer to neighbour
        forAllConstIter(Map<Map<label> >, auxSurfPointMap_, mIter)
        {
            OPstream toProc(Pstream::blocking, mIter.key(), sizeof(label));

            toProc << mIter().size();

            // Set an entry for this processor
            nProcSize.insert(mIter.key(), -1);
        }

        // Receive from neighbour
        forAllConstIter(Map<Map<label> >, auxSurfPointMap_, mIter)
        {
            IPstream fmProc(Pstream::blocking, mIter.key(), sizeof(label));

            fmProc >> nProcSize[mIter.key()];

            Pout<< " neiProcNo: " << mIter.key()
                << " mySize: " << mIter().size()
                << " neiSize: " << nProcSize[mIter.key()]
                << endl;

            if (nProcSize[mIter.key()] != mIter().size())
            {
                FatalErrorIn
                (
                    "void mesquiteMotionSolver::initParallelConnectivity()"
                )
                    << " Buffer size mismatch: " << nl
                    << " My proc: " << Pstream::myProcNo() << nl
                    << " neiProcNo: " << mIter.key()
                    << " mySize: " << mIter().size()
                    << " neiSize: " << nProcSize[mIter.key()]
                    << abort(FatalError);
            }
        }
    }

    // Match points geometrically
    forAllConstIter(Map<Map<label> >, auxSurfPointMap_, mIter)
    {
        label proc = mIter.key();
        label bSize = mIter().size();

        // Size-up send / recv fields
        sendFields_.insert(proc, vectorField(bSize, vector::zero));
        recvFields_.insert(proc, vectorField(bSize, vector::zero));

        // Fetch references
        vectorField& recvField = recvFields_[proc];

        // Send and receive points
        IPstream::read
        (
            Pstream::nonBlocking,
            proc,
            reinterpret_cast<char*>(recvField.begin()),
            recvField.byteSize()
        );

        DynamicList<point>& bufPoints = auxSurfPoints[proc];

        OPstream::write
        (
            Pstream::nonBlocking,
            proc,
            reinterpret_cast<const char*>(&bufPoints[0]),
            bufPoints.size()*sizeof(vector)
        );
    }

    // Wait for all transfers to complete.
    OPstream::waitRequests();
    IPstream::waitRequests();

    forAllIter(Map<Map<label> >, auxSurfPointMap_, mIter)
    {
        label proc = mIter.key();

        // Fetch references
        vectorField& recvField = recvFields_[proc];
        DynamicList<point>& bufPoints = auxSurfPoints[proc];

        labelList pMap(recvField.size(), -1);

        bool matchedAll =
        (
            matchPoints
            (
                bufPoints,
                recvField,
                List<scalar>(recvField.size(), SMALL),
                true,
                pMap
            )
        );

        if (matchedAll)
        {
            // Renumber to shuffled indices
            Map<label>& anbMap = auxNeiBufferMap_[proc];

            forAllIter(Map<label>, anbMap, pIter)
            {
                pIter() = pMap[pIter()];
            }
        }
        else
        {
            FatalErrorIn
            (
                "void mesquiteMotionSolver::initParallelConnectivity()"
            )
                << " Failed point match: " << nl
                << " My proc: " << Pstream::myProcNo() << nl
                << " neiProcNo: " << proc << nl
                << abort(FatalError);
        }
    }

    // Prepare edgeMarkers
    forAll(pIDs_, patchI)
    {
        const edgeList& edges = boundary[pIDs_[patchI]].edges();

        for
        (
            label i = boundary[pIDs_[patchI]].nInternalEdges();
            i < edges.size();
            i++
        )
        {
            // If both points on this edge are marked,
            // this edge needs to be left out.
            label i0 = edges[i][0] + offsets_[patchI];
            label i1 = edges[i][1] + offsets_[patchI];

            bool p0 = (pointMarker_[i0] < 0.5);
            bool p1 = (pointMarker_[i1] < 0.5);

            if (p0 && p1)
            {
                edgeMarker_[patchI][i] = 0.0;
            }

            // Check for corner points.
            // If only one vertex is marked,
            // fix the boundary condition vector
            if (bdyMarker[i0] && !bdyMarker[i1])
            {
                bdy_[i0] = vector::zero;
            }
            else
            if (!bdyMarker[i0] && bdyMarker[i1])
            {
                bdy_[i1] = vector::zero;
            }
        }
    }
}


// Copy auxiliary points to/from buffers
void mesquiteMotionSolver::copyAuxiliaryPoints(bool copyBack)
{
    if (!Pstream::parRun())
    {
        return;
    }

    // Search the registry for all mapping objects.
    HashTable<const coupleMap*> coupleMaps = Mesh_.lookupClass<coupleMap>();

    forAllIter(HashTable<const coupleMap*>, coupleMaps, cmIter)
    {
        coupleMap& cMap = const_cast<coupleMap&>(*cmIter());

        if (cMap.isLocal())
        {
            continue;
        }

        // Obtain a reference to the point buffer
        pointField& pField = cMap.pointBuffer();

        // Pick maps for which this sub-domain is master
        if (cMap.masterIndex() == Pstream::myProcNo())
        {
            // Fetch the neighbour's procIndex
            label neiProcNo = cMap.slaveIndex();

            const Map<label>& pointMap = auxPointMap_[neiProcNo];

            label pIndex = -1;
            vector v = vector::zero;

            if (copyBack)
            {
                if (cMap.isRecv())
                {
                    forAll(pField, pointI)
                    {
                        pIndex = pointMap[pointI];

                        v.x() = vtxCoords_[(3*pIndex)+0];
                        v.y() = vtxCoords_[(3*pIndex)+1];
                        v.z() = vtxCoords_[(3*pIndex)+2];

                        pField[pointI] = v;
                    }

                    if (debug)
                    {
                        Pout << "Sending to proc: "
                             << neiProcNo
                             << endl;
                    }

                    // Send points to the slave.
                    OPstream::write
                    (
                        Pstream::nonBlocking,
                        neiProcNo,
                        reinterpret_cast<const char*>(&pField[0]),
                        pField.size()*sizeof(vector)
                    );
                }
            }
            else
            {
                if (cMap.isRecv())
                {
                    forAll(pField, pointI)
                    {
                        pIndex = pointMap[pointI];

                        v = pField[pointI];

                        vtxCoords_[(3*pIndex)+0] = v.x();
                        vtxCoords_[(3*pIndex)+1] = v.y();
                        vtxCoords_[(3*pIndex)+2] = v.z();
                    }
                }
            }
        }
        else
        if (copyBack && cMap.isSend())
        {
            if (debug)
            {
                Pout << "Receiving from proc: "
                     << cMap.masterIndex()
                     << endl;
            }

            // This sub-domain is a slave.
            // Prepare for buffer receipt.
            IPstream::read
            (
                Pstream::nonBlocking,
                cMap.masterIndex(),
                reinterpret_cast<char*>(pField.begin()),
                pField.byteSize()
            );
        }
    }

    if (copyBack)
    {
        // Wait for all transfers to complete.
        OPstream::waitRequests();
        IPstream::waitRequests();

        DynamicList<label> procRanks(10);

        // Update points in descending order of processors,
        // so that the lowest rank is updated last.
        forAllIter(HashTable<const coupleMap*>, coupleMaps, cmIter)
        {
            coupleMap& cMap = const_cast<coupleMap&>(*cmIter());

            if (cMap.isLocal() || cMap.isRecv())
            {
                continue;
            }

            if (cMap.slaveIndex() == Pstream::myProcNo())
            {
                procRanks.append(cMap.masterIndex());
            }
        }

        procRanks.shrink();
        sort(procRanks);

        forAllReverse(procRanks, indexI)
        {
            // Now copy updated locations to refPoints.
            forAllIter(HashTable<const coupleMap*>, coupleMaps, cmIter)
            {
                coupleMap& cMap = const_cast<coupleMap&>(*cmIter());

                if (cMap.isLocal() || cMap.isRecv())
                {
                    continue;
                }

                // Obtain a reference to the point buffer
                const pointField& pField = cMap.pointBuffer();

                // Pick maps for which this sub-domain is a slave
                if
                (
                    cMap.slaveIndex() == Pstream::myProcNo() &&
                    cMap.masterIndex() == procRanks[indexI]
                )
                {
                    const labelList& smPoints = cMap.subMeshPoints();

                    // Only update shared-point positions.
                    forAll(smPoints, pointI)
                    {
                        refPoints_[smPoints[pointI]] = pField[pointI];
                    }
                }
            }
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
        const edgeList& edges = boundary[pIDs_[patchI]].edges();

        forAll(edges, edgeI)
        {
            gradEdge_[patchI][edgeI] =
            (
                p[edges[edgeI][1] + offsets_[patchI]]
              - p[edges[edgeI][0] + offsets_[patchI]]
            );
        }
    }

    // Divergence (e2n)
    forAll(pIDs_, patchI)
    {
        const edgeList& edges = boundary[pIDs_[patchI]].edges();

        forAll(edges, edgeI)
        {
            gradEdge_[patchI][edgeI] *= edgeMarker_[patchI][edgeI];

            w[edges[edgeI][0] + offsets_[patchI]] += gradEdge_[patchI][edgeI];
            w[edges[edgeI][1] + offsets_[patchI]] -= gradEdge_[patchI][edgeI];
        }
    }

    // Transfer buffers after divergence compute.
    transferBuffers(w);

    // Apply boundary conditions
    applyBCs(w);
}


// Transfer buffers after divergence compute.
void mesquiteMotionSolver::transferBuffers
(
    vectorField& field
)
{
    if (!Pstream::parRun())
    {
        return;
    }

    forAllConstIter(Map<Map<label> >, auxNeiBufferMap_, mIter)
    {
        label neiProcNo = mIter.key();
        const Map<label>& pointMap = mIter();

        // Fetch references
        vectorField& recvField = recvFields_[neiProcNo];
        vectorField& sendField = sendFields_[neiProcNo];

        // Schedule receipt from neighbour
        IPstream::read
        (
            Pstream::nonBlocking,
            neiProcNo,
            reinterpret_cast<char*>(recvField.begin()),
            recvField.byteSize()
        );

        // Prepare a send buffer.
        forAllConstIter(Map<label>, pointMap, pIter)
        {
            sendField[pIter()] = field[pIter.key()];
        }

        OPstream::write
        (
            Pstream::nonBlocking,
            neiProcNo,
            reinterpret_cast<const char*>(&sendField[0]),
            sendField.size()*sizeof(vector)
        );
    }

    // Wait for all transfers to complete.
    OPstream::waitRequests();
    IPstream::waitRequests();

    forAllConstIter(Map<Map<label> >, auxSurfPointMap_, mIter)
    {
        label neiProcNo = mIter.key();
        const Map<label>& pointMap = mIter();

        // Fetch references
        const vectorField& recvField = recvFields_[neiProcNo];

        // Prepare a send buffer.
        forAllConstIter(Map<label>, pointMap, pIter)
        {
            field[pIter.key()] += recvField[pIter()];
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
        // Apply slip conditions for internal nodes
        forAll(pNormals_[patchI], pointI)
        {
            const vector& n = pNormals_[patchI][pointI];

            field[pointI + offsets_[patchI]] -=
            (
                (field[pointI + offsets_[patchI]] & n)*n
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
    if (min(bdy_) > vector(0.5,0.5,0.5))
    {
        Random randomizer(1);
        label nFix = (field.size()*5)/100;

        for(label i = 0; i < nFix; i++)
        {
            field[randomizer.integer(0, field.size()-1)] = vector::zero;
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

    vectorField nFw = (w - tmpField);
    vectorField nFb = (b - tmpField);

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

    forAll(field,i)
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
    scalar alpha, beta, rho, rhoOld, residual;
    label maxIter = x.size(), iter = 0;

    // Compute initial residual
    A(x,w);

    // Compute the normFactor, using 'r' as scratch-space
    scalar norm = normFactor(x,b,w,r);

    if (debug)
    {
        Info<< "normFactor: " << norm << endl;
    }

    r = b - w;
    p = r;
    rho = dot(r,p);

    // Obtain the normalized residual
    residual = cmptSumMag(r)/norm;

    Info<< " Initial residual: " << residual;

    while ( (iter < maxIter) && (residual > tolerance_) )
    {
        A(p,w);

        alpha = rho / dot(p,w);

        forAll (x, i)
        {
            x[i] += (alpha*p[i]);
            r[i] -= (alpha*w[i]);
        }

        rhoOld = rho;

        rho = dot(r,r);

        beta = rho / rhoOld;

        forAll (p, i)
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

    // Check the dictionary for entries corresponding to constant
    // fixed-displacement BCs. This is done because a 'motionU'
    // field is not used to specify such BC types.
    if (optionsDict.found("fixedValuePatches"))
    {
        const polyBoundaryMesh& boundary = mesh().boundaryMesh();
        const dictionary& fvpDict = optionsDict.subDict("fixedValuePatches");

        // Extract a list of patch names.
        wordList fixPatches = fvpDict.toc();

        // Construct a pointMesh.
        pointMesh pMesh(Mesh_);

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

        // Now update refPoints with patch values
        refPoints_ += dPointField;

        // Apply values to points common with coupleMaps

        // Search the registry for all mapping objects.
        if (Pstream::parRun() && !twoDMesh_)
        {
            HashTable<const coupleMap*> cMaps = Mesh_.lookupClass<coupleMap>();

            forAllConstIter(HashTable<const coupleMap*>, cMaps, cmIter)
            {
                const coupleMap& cMap = *cmIter();

                if (cMap.isLocal() || cMap.isSend())
                {
                    continue;
                }

                // Obtain a reference to the point buffer
                pointField& pField = cMap.pointBuffer();

                // Pick maps for which this sub-domain is master
                if (cMap.masterIndex() == Pstream::myProcNo())
                {
                    // Fetch the neighbour's procIndex
                    label neiProcNo = cMap.slaveIndex();

                    const Map<label>& pointMap = auxPointMap_[neiProcNo];

                    label pIndex = -1;

                    forAll(pField, pointI)
                    {
                        pIndex = pointMap[pointI];

                        if (pIndex < label(nPoints_))
                        {
                            pField[pointI] = refPoints_[pIndex];
                        }
                    }
                }
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
            const labelList& meshPts = boundary[pIDs_[patchI]].meshPoints();

            forAll(meshPts,pointI)
            {
                xV_[pointI + offsets_[patchI]] = refPoints_[meshPts[pointI]];
            }
        }

        Info<< "Solving for point motion: ";

        label iters = CG(bV_, pV_, rV_, wV_, xV_);

        Info<< " No Iterations: " << iters << endl;

        // Update refPoints (with relaxation if necessary)
        forAll(pIDs_, patchI)
        {
            const labelList& meshPts = boundary[pIDs_[patchI]].meshPoints();

            forAll(meshPts,pointI)
            {
                refPoints_[meshPts[pointI]] =
                (
                    (relax_ * xV_[pointI + offsets_[patchI]])
                  + ((1.0 - relax_) * origPoints_[meshPts[pointI]])
                );
            }
        }

        // Update locally coupled patches
        if (patchCoupling_.size())
        {
            // Search the registry for all mapping objects.
            HashTable<const coupleMap*> coupleMaps =
            (
                Mesh_.lookupClass<coupleMap>()
            );

            forAllIter(Map<label>, patchCoupling_, patchI)
            {
                const labelList& meshPts = boundary[patchI.key()].meshPoints();

                bool foundMap = false;

                forAllIter(HashTable<const coupleMap*>, coupleMaps, cmIter)
                {
                    const coupleMap& cMap = *(cmIter());

                    // Ensure that both master and slave patches match.
                    if
                    (
                        (cMap.masterIndex() == patchI.key()) &&
                        (cMap.slaveIndex() == patchI()) &&
                        (cMap.isLocal())
                    )
                    {
                        const Map<label>& mtsMap =
                        (
                            cMap.entityMap(coupleMap::POINT)
                        );

                        forAll(meshPts, pointI)
                        {
                            refPoints_[mtsMap[meshPts[pointI]]] =
                            (
                                refPoints_[meshPts[pointI]]
                            );
                        }

                        foundMap = true;

                        break;
                    }
                }

                if (!foundMap)
                {
                    FatalErrorIn("void mesquiteMotionSolver::smoothSurfaces()")
                        << "Could not find coupling map: " << nl
                        << "Master: " << patchI.key() << nl
                        << "Slave: " << patchI()
                        << abort(FatalError);
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
    const pointField& pField
)
{
    const cell& cellToCheck = mesh().cells()[cIndex];

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

                return sign(V)*((24.96100588*::cbrt(V*V))/Le);
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

                return sign(V)*((24.96100588*::cbrt(V*V))/Le);
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

    // Loop through pointCells for all boundary points
    // and compute cell volume.
    const labelListList& pointCells = mesh().pointCells();
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    // Check if a minimum quality was specified.
    scalar thresh = 0.45;

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

    if (invCells.size() == 0)
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
                    << abort(FatalError);
            }

            const dictionary& pD = constraintDict.subDict(cstrPatches[wordI]);

            // Read info.
            vector axisPoint(pD.lookup("axisPoint"));
            vector axisVector(pD.lookup("axisVector"));
            scalar radius = readScalar(pD.lookup("radius"));

            const labelList& meshPts = boundary[pID].meshPoints();

            axisVector /= mag(axisVector) + VSMALL;

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
                        WarningIn
                        (
                            "void mesquiteMotionSolver::"
                            "enforceCylindricalConstraints()"
                        )
                            << " Constraint violation: " << viol << endl;
                    }
                }
            }

            // If this is a coupled patch, correct the slave as well.
            if (patchCoupling_.found(pID))
            {
                // Search the registry for all mapping objects.
                HashTable<const coupleMap*> coupleMaps =
                (
                    Mesh_.lookupClass<coupleMap>()
                );

                const labelList& mMeshPts = boundary[pID].meshPoints();

                bool foundMap = false;

                forAllIter
                (
                    HashTable<const coupleMap*>,
                    coupleMaps,
                    cmIter
                )
                {
                    const coupleMap& cMap = *(cmIter());

                    // Ensure that both master and slave patches match.
                    if
                    (
                        (cMap.masterIndex() == pID) &&
                        (cMap.slaveIndex() == patchCoupling_[pID]) &&
                        (cMap.isLocal())
                    )
                    {
                        const Map<label>& mtsMap =
                        (
                            cMap.entityMap(coupleMap::POINT)
                        );

                        forAll(mMeshPts, pointI)
                        {
                            refPoints_[mtsMap[mMeshPts[pointI]]] =
                            (
                                refPoints_[mMeshPts[pointI]]
                            );
                        }

                        foundMap = true;

                        break;
                    }
                }

                if (!foundMap)
                {
                    FatalErrorIn
                    (
                        "void mesquiteMotionSolver::"
                        "enforceCylindricalConstraints()"
                    )
                        << "Could not find coupling map: " << nl
                        << "Master: " << pID << nl
                        << "Slave: " << patchCoupling_[pID]
                        << abort(FatalError);
                }
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

    // Search the registry for all mapping objects.
    HashTable<const coupleMap*> coupleMaps = Mesh_.lookupClass<coupleMap>();

    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    forAll(pIDs_, patchI)
    {
        // First update localPoints with latest point positions
        const labelList& meshPts = boundary[pIDs_[patchI]].meshPoints();

        forAll(meshPts,pointI)
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

        // Add contributions from neighbouring processors
        forAllIter(HashTable<const coupleMap*>, coupleMaps, cmIter)
        {
            const coupleMap& cMap = *(cmIter());

            if (cMap.isLocal() || cMap.isSend())
            {
                continue;
            }

            label neiProcNo = -1;
            bool slaveProc = false;

            // Fetch the neighbour's procIndex
            if (cMap.masterIndex() == Pstream::myProcNo())
            {
                neiProcNo = cMap.slaveIndex();
            }
            else
            if (cMap.slaveIndex() == Pstream::myProcNo())
            {
                neiProcNo = cMap.masterIndex();
                slaveProc = true;
            }

            // Fetch reference
            const Map<label>& apMap = auxPointMap_[neiProcNo];

            // Fetch buffers from coupleMap
            const pointField& mPts = cMap.pointBuffer();
            const labelList& starts = cMap.entityBuffer(coupleMap::FACE_STARTS);
            const labelList& sizes = cMap.entityBuffer(coupleMap::FACE_SIZES);
            const label nShared = cMap.nEntities(coupleMap::SHARED_POINT);
            const faceList& fList = cMap.faces();

            // Physical patches are ordered the same way
            // in coupleMaps, so pIDs_ can be used directly.
            label fStart = starts[pIDs_[patchI]];
            label fSize = sizes[pIDs_[patchI]];

            for (label faceI = fStart; faceI < fStart + fSize; faceI++)
            {
                const face& thisFace = fList[faceI];

                vector n = thisFace.normal(mPts);

                forAll(thisFace, pI)
                {
                    if (thisFace[pI] >= nShared)
                    {
                        continue;
                    }

                    // Fetch the global index.
                    label gIndex = apMap[thisFace[pI]];
                    label lIndex = boundary[pIDs_[patchI]].whichPoint(gIndex);

                    if (lIndex > -1)
                    {
                        pNormals_[patchI][lIndex] += n;
                    }
                    else
                    {
                        // Something is wrong here.
                        Pout<< " Unable to find local index. " << nl
                            << " gIndex: " << gIndex << nl
                            << " lIndex: " << lIndex << nl
                            << abort(FatalError);
                    }
                }
            }
        }
    }

    // Normalize all point-normals
    forAll(pIDs_, patchI)
    {
        pNormals_[patchI] /= mag(pNormals_[patchI]) + VSMALL;
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

        // Enforce constraints, if necessary
        enforceCylindricalConstraints();
    }

    // Copy refPoints prior to internal smoothing
    origPoints_ = refPoints_;

    // Copy most recent point positions
    forAll(refPoints_, pointI)
    {
        vtxCoords_[(3*pointI)+0] = refPoints_[pointI][0];
        vtxCoords_[(3*pointI)+1] = refPoints_[pointI][1];
        vtxCoords_[(3*pointI)+2] = refPoints_[pointI][2];
    }

    // Copy auxiliary points from slaves
    copyAuxiliaryPoints();

    Mesquite::MsqError err;

    //- ArrayMesh object defined by Mesquite
    Mesquite::ArrayMesh
    msqMesh
    (
        3,                         // Number of coords per vertex
        nPoints_,                  // Number of vertices
        vtxCoords_,                // The vertex coordinates
        fixFlags_,                 // Fixed vertex flags
        nCells_,                   // Number of elements
        Mesquite::TETRAHEDRON,     // Element type
        cellToNode_,               // Connectivity
        false,                     // Fortran-style array indexing
        4                          // Number of nodes per element
    );

    // Create an instruction queue
    Mesquite::InstructionQueue queue;

    // Apply termination criteria to the optimization algorithm
    optAlgorithm_->set_outer_termination_criterion(&tcOuter_);
    optAlgorithm_->set_inner_termination_criterion(&tcInner_);

    // Set up the quality assessor
    Mesquite::QualityAssessor qA(&qMetricTable_[qMetric_]());

    // Assess the quality of the initial mesh before smoothing
    queue.add_quality_assessor(&qA, err);

    // Set the master quality improver
    queue.set_master_quality_improver
    (
        &optAlgorithm_(),
        err
    );

    // Assess the quality of the final mesh after smoothing
    queue.add_quality_assessor(&qA, err);

    // Disable slave output for parallel runs.
    if (Pstream::parRun() && !Pstream::master())
    {
        qA.disable_printing_results();
    }

    // Launches optimization on the mesh
    queue.run_instructions(&msqMesh, err);

    // Copy updated positions
    forAll(refPoints_, pointI)
    {
        refPoints_[pointI][0] = vtxCoords_[(3*pointI)+0];
        refPoints_[pointI][1] = vtxCoords_[(3*pointI)+1];
        refPoints_[pointI][2] = vtxCoords_[(3*pointI)+2];

        // Relax points
        refPoints_[pointI] =
        (
            (relax_ * refPoints_[pointI])
          + ((1.0 - relax_) * origPoints_[pointI])
        );
    }

    // Copy auxiliary points back to slaves
    copyAuxiliaryPoints(true);
}


//- Move points
bool mesquiteMotionSolver::movePoints() const
{
    return true;
}


//- Update topology
void mesquiteMotionSolver::updateMesh(const mapPolyMesh& mpm)
{
    update(mpm);
}


//- Update topology (using meshObjectBase)
bool mesquiteMotionSolver::updateMesh(const mapPolyMesh& mpm) const
{
    const_cast<mesquiteMotionSolver*>(this)->update(mpm);

    return true;
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
    }

    nPoints_ = Mesh_.nPoints();
    nCells_  = Mesh_.nCells();

    nAuxPoints_ = 0;
    nAuxCells_  = 0;

    // Reset refPoints
    refPoints_.clear();
    refPoints_ = Mesh_.points();

    // Clear the auxiliary point map
    auxPointMap_.clear();
    auxSurfPointMap_.clear();
    auxNeiBufferMap_.clear();
    sendFields_.clear();
    recvFields_.clear();

    // Clear Mesquite arrays
    clearOut();

    // Reset flag
    arraysInitialized_ = false;
}

} // End namespace Foam

// ************************************************************************* //
