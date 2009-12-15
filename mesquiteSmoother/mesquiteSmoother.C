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

#include "mesquiteSmoother.H"
#include "IOmanip.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "tensor2D.H"

// Forward declarations for LAPACK routines
extern "C"
void dgeqrf_
(
    const int* m,
    const int* n,
    double a[],
    const int* lda,
    double tau[],
    double work[],
    int* lwork,
    int* info
);

extern "C"
void dorgqr_
(
    const int* m,
    const int* n,
    const int* k,
    double a[],
    const int* lda,
    double tau[],
    double work[],
    int* lwork,
    int* info
);

extern "C"
void dpotrf_
(
    const char *uplo,
    const int* n,
    double a[],
    const int* lda,
    int* info
);

extern "C"
void dgesvd_
(
    const char* jobu,
    const char* jobvt,
    const int* m,
    const int* n,
    double a[],
    const int* lda,
    double s[],
    double u[],
    const int* ldu,
    double vt[],
    const int* ldvt,
    double work[],
    const int* lwork,
    int* info
);

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mesquiteSmoother, 0);
    addToRunTimeSelectionTable
    (
        motionSolver,
        mesquiteSmoother,
        dictionary
    );
}

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mesquiteSmoother::mesquiteSmoother
(
    const polyMesh& mesh
)
:
    motionSolver(mesh),
    Mesh_(mesh),
    nPoints_(mesh.nPoints()),
    nCells_(mesh.nCells()),
    surfaceSmoothing_(false),
    prioritizeVertices_(false),
    tolerance_(1e-4),
    nSweeps_(1),
    surfInterval_(1),
    vtxCoords_(NULL),
    cellToNode_(NULL),
    fixFlags_(NULL),
    refPoints_
    (
        IOobject
        (
            "refPoints",
            mesh.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh.points()
    )
{
    // Read options from the dictionary
    readOptions();

    // Initialize connectivity arrays for Mesquite
    initArrays();
}

mesquiteSmoother::mesquiteSmoother
(
    const polyMesh& mesh,
    Istream& msData
)
:
    motionSolver(mesh),
    Mesh_(mesh),
    nPoints_(mesh.nPoints()),
    nCells_(mesh.nCells()),
    surfaceSmoothing_(false),
    prioritizeVertices_(false),
    tolerance_(1e-4),
    nSweeps_(1),
    surfInterval_(1),
    vtxCoords_(NULL),
    cellToNode_(NULL),
    fixFlags_(NULL),
    refPoints_
    (
        IOobject
        (
            "refPoints",
            mesh.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh.points()
    )
{
    // Read options from the dictionary
    readOptions();

    // Initialize connectivity arrays for Mesquite
    initArrays();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mesquiteSmoother::~mesquiteSmoother()
{
    delete [] vtxCoords_;
    delete [] cellToNode_;
    delete [] fixFlags_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Read options from the dictionary
void mesquiteSmoother::readOptions()
{
    // Check if any slip patches are specified
    if (found("slipPatches"))
    {
        wordList slipPatches = subDict("slipPatches").toc();

        forAll(slipPatches, wordI)
        {
            word& patchName = slipPatches[wordI];

            forAll(mesh().boundaryMesh(), patchI)
            {
                if (mesh().boundaryMesh()[patchI].name() == patchName)
                {
                    slipPatchIDs_.insert(patchI);
                }
            }
        }

        // Toggle surface smoothing if slip patches are present
        if (!slipPatches.empty())
        {
            surfaceSmoothing_ = true;
        }

        // Check if a tolerance has been specified
        if (found("tolerance"))
        {
            tolerance_ = readScalar(lookup("tolerance"));
        }

        // Check if multiple sweeps have been requested
        if (found("nSweeps"))
        {
            nSweeps_ = readLabel(lookup("nSweeps"));
        }

        // Check if a surface smoothing interval has been specified
        if (found("surfInterval"))
        {
            surfInterval_ = readLabel(lookup("surfInterval"));
        }

        // Check if coupled patches exist.
        if (found("coupledPatches"))
        {
            dictionary coupledPatches = subDict("coupledPatches");

            const polyBoundaryMesh& boundary = mesh().boundaryMesh();

            // Determine master and slave patches
            wordList masterPatches = coupledPatches.toc();

            forAll(masterPatches, wordI)
            {
                // Lookup the slave patch
                word masterPatch = masterPatches[wordI];
                word slavePatch  = coupledPatches.lookup(masterPatch);

                // Determine patch indices
                label mPatch = -1, sPatch = -1;

                forAll(boundary,patchI)
                {
                    if (boundary[patchI].name() == masterPatch)
                    {
                        mPatch = patchI;
                    }

                    if (boundary[patchI].name() == slavePatch)
                    {
                        sPatch = patchI;
                    }
                }

                // It is considered an error to have slave-patches
                // on the slip-patches list.
                if (slipPatchIDs_.found(sPatch))
                {
                    FatalErrorIn
                    (
                        "mesquiteSmoother::readOptions()"
                    )
                        << " Slave Patch: " << slavePatch << nl
                        << " is specified in the slip-patches list." << nl
                        << " Please remove the entry."
                        << abort(FatalError);
                }

                if (mPatch == -1 && sPatch == -1)
                {
                    continue;
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
                    FatalErrorIn
                    (
                        "mesquiteSmoother::readOptions()"
                    )
                        << " Coupled patches are wrongly specified." << nl
                        << " Master: " << mPatch << ":" << masterPatch << nl
                        << " Slave: " << sPatch << ":" << slavePatch << nl
                        << abort(FatalError);
                }
            }
        }
    }

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
    qMetric_ = word(lookup("optMetric"));

    if (!qMetricTable_.found(qMetric_))
    {
        FatalErrorIn
        (
            "mesquiteSmoother::readOptions()"
        )
            << "Unrecognized quality metric: "
            << qMetric_ << nl
            << "Available metrics are: " << nl
            << qMetricTable_.toc()
            << abort(FatalError);
    }
    else
    {
        Info << "Selecting quality metric: " << qMetric_ << endl;
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
    word ofType(lookup("objFunction"));

    if (!ofTable.found(ofType))
    {
        FatalErrorIn
        (
            "mesquiteSmoother::readOptions()"
        )
            << "Unrecognized objective function: "
            << ofType << nl
            << "Available types are: " << nl
            << ofTable.toc()
            << abort(FatalError);
    }
    else
    {
        Info << "Selecting objective function: " << ofType << endl;
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
        types[0] = ofTable[word(lookup("firstFunction"))];
        types[1] = ofTable[word(lookup("secondFunction"))];

        // Ensure that we're not making a composite of composites
        if (types[0] < 4 || types[1] < 4)
        {
            FatalErrorIn
            (
                "mesquiteSmoother::readOptions()"
            )
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
        types[0] = ofTable[word(lookup("scaleFunction"))];

        // Lookup the scale value
        scale = readScalar(lookup("scale"));

        // Ensure that we're not making a composite of composites
        if (types[0] < 4)
        {
            FatalErrorIn
            (
                "mesquiteSmoother::readOptions()"
            )
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

    for(label i = 0; i < numTries; i++)
    {
        switch(types[i])
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
                label pValue = readLabel(lookup("pValue"));

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
                scalar power = readScalar(lookup("power"));

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
                scalar power = readScalar(lookup("power"));

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
                FatalErrorIn
                (
                    "mesquiteSmoother::readOptions()"
                )
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
    word oaType(lookup("optAlgorithm"));

    if (!oaTable.found(oaType))
    {
        FatalErrorIn
        (
            "mesquiteSmoother::readOptions()"
        )
            << "Unrecognized optimization algorithm: "
            << oaType << nl
            << "Available types are: " << nl
            << oaTable.toc()
            << abort(FatalError);
    }
    else
    {
        Info << "Selecting optimization algorithm: " << oaType << endl;
    }

    // Instantiate the appropriate objective function
    label oaSelection = oaTable[oaType];

    switch(oaSelection)
    {
        case 0:
            optAlgorithm_.set
            (
                new Mesquite::ConjugateGradient(&objFunction_())
            );

            break;

        case 1:
            optAlgorithm_.set
            (
                new Mesquite::FeasibleNewton(&objFunction_())
            );

            break;

        case 2:
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

        case 3:
            optAlgorithm_.set
            (
                new Mesquite::QuasiNewton(&objFunction_())
            );

            break;

        case 4:
            // Not really a smoother; just randomizes positions
            optAlgorithm_.set
            (
                new Mesquite::Randomize()
            );

            break;

        case 5:
            optAlgorithm_.set
            (
                new Mesquite::LaplacianSmoother(&objFunction_())
            );

            break;

        case 6:
            optAlgorithm_.set
            (
                new Mesquite::SmartLaplacianSmoother(&objFunction_())
            );

            break;

        case 7:
            optAlgorithm_.set
            (
                new Mesquite::SteepestDescent(&objFunction_())
            );

            break;

        case 8:
            optAlgorithm_.set
            (
                new Mesquite::TrustRegion(&objFunction_())
            );

            break;
    }

    // Read termination criteria, if it exists.
    if (found("tcInner"))
    {
        const dictionary& innerDict = subDict("tcInner");

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
            FatalErrorIn
            (
                "mesquiteSmoother::readOptions()"
            )
                << "Empty tcInner dictionary: " << nl
                << "Available types are: " << nl
                << options.toc()
                << abort(FatalError);
        }
    }
    else
    {
        Info << "Inner termination criterion (tcInner) "
             << "was not found. Using default values."
             << endl;

        tcInner_.add_absolute_gradient_inf_norm(1e-4);
    }

    if (found("tcOuter"))
    {
        const dictionary& outerDict = subDict("tcOuter");

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
            FatalErrorIn
            (
                "mesquiteSmoother::readOptions()"
            )
                << "Empty tcOuter dictionary: " << nl
                << "Available types are: " << nl
                << options.toc()
                << abort(FatalError);
        }
    }
    else
    {
        Info << "Outer termination criterion (tcOuter) "
             << "was not found. Using default values."
             << endl;

        tcOuter_.add_iteration_limit(1);
    }
}

// Initialize connectivity arrays
void mesquiteSmoother::initArrays()
{
    vtxCoords_ = new double[nPoints_*3];
    cellToNode_ = new unsigned long[nCells_*4];
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

    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    forAll(boundary, patchI)
    {
        const labelList& meshPointLabels = boundary[patchI].meshPoints();

        forAll(meshPointLabels, pointI)
        {
            fixFlags_[meshPointLabels[pointI]] = 1;
        }
    }

    if (surfaceSmoothing_)
    {
        // Extract the patch list
        pIDs_ = slipPatchIDs_.toc();

        pNormals_.setSize(pIDs_.size());
        localPts_.setSize(pIDs_.size());
        jPoints_.setSize(pIDs_.size());
        T_.setSize(pIDs_.size());
        a_.setSize(pIDs_.size());
        pFlags_.setSize(pIDs_.size());

        // Allocate for function values if we're prioritizing.
        if (prioritizeVertices_)
        {
            meanFnValues_.setSize(pIDs_.size());
            pFnValues_.setSize(pIDs_.size());
        }

        // Obtain pointCells connectivity
        const labelListList& pointCells = mesh().pointCells();

        forAll(pIDs_, patchI)
        {
            label nPts = boundary[pIDs_[patchI]].nPoints();
            label nIntEdges = boundary[pIDs_[patchI]].nInternalEdges();

            pNormals_[patchI].setSize(nPts, vector::zero);
            localPts_[patchI].setSize(nPts, vector::zero);

            // Initialize the transform tensor field
            T_[patchI].setSize(nPts, tensor::zero);

            // Allocate coefficients.
            // Six parametric coefficients per point.
            a_[patchI].setSize(nPts, scalarField(6, 0.0));

            // Flag to denote whether the point
            // lies on a patch boundary (true) or the interior (false).
            pFlags_[patchI].setSize(nPts, false);

            forAll(pFlags_[patchI], pointI)
            {
                const labelList& pEdges =
                (
                    boundary[pIDs_[patchI]].pointEdges()[pointI]
                );

                forAll(pEdges, edgeI)
                {
                    if (pEdges[edgeI] >= nIntEdges)
                    {
                        // Assign the flag
                        pFlags_[patchI][pointI] = true;
                        break;
                    }
                }
            }

            if (prioritizeVertices_)
            {
                pFnValues_[patchI].setSize(nPts, 0.0);
            }

            // Connectivity information for point Jacobians
            jPoints_[patchI].setSize(nPts, labelList(0));

            const labelList& meshPts = boundary[pIDs_[patchI]].meshPoints();

            forAll(meshPts, pointI)
            {
                label pIndex = meshPts[pointI];

                const labelList& pCells = pointCells[pIndex];

                // First size the list. Assume tet-mesh.
                jPoints_[patchI][pointI].setSize(3*pCells.size());

                label i = 0;

                forAll(pCells, cellI)
                {
                    const cell& curCell = meshCells[pCells[cellI]];

                    // Find a face which doesn't contain pIndex.
                    forAll(curCell, faceI)
                    {
                        const face& thisFace = meshFaces[curCell[faceI]];

                        if (thisFace.which(pIndex) == -1)
                        {
                            if (owner[curCell[faceI]] == pCells[cellI])
                            {
                                jPoints_[patchI][pointI][i++] = thisFace[2];
                                jPoints_[patchI][pointI][i++] = thisFace[1];
                                jPoints_[patchI][pointI][i++] = thisFace[0];
                            }
                            else
                            {
                                jPoints_[patchI][pointI][i++] = thisFace[0];
                                jPoints_[patchI][pointI][i++] = thisFace[1];
                                jPoints_[patchI][pointI][i++] = thisFace[2];
                            }

                            break;
                        }
                    }
                }
            }
        }

        // Build coupling maps.
        if (patchCoupling_.size())
        {
            masterToSlave_.setSize(boundary.size());

            const pointField& points = mesh().points();

            forAllIter(Map<label>, patchCoupling_, patchI)
            {
                const labelList& mLabels = boundary[patchI.key()].meshPoints();
                const labelList& sLabels = boundary[patchI()].meshPoints();

                label nMatchedPoints = 0;

                forAll(mLabels, pointI)
                {
                    forAll(sLabels, pointJ)
                    {
                        if
                        (
                            mag
                            (
                                points[mLabels[pointI]]
                              - points[sLabels[pointJ]]
                            ) < 1e-20
                        )
                        {
                            // Add a map entry
                            masterToSlave_[patchI.key()].insert
                            (
                                mLabels[pointI],
                                sLabels[pointJ]
                            );

                            nMatchedPoints++;

                            break;
                        }
                    }
                }

                // Make sure we were successful.
                if (nMatchedPoints != mLabels.size())
                {
                    FatalErrorIn("mesquiteSmoother::initArrays()")
                        << " Failed to match all points." << nl
                        << " Number of points required for match: "
                        << mLabels.size() << nl
                        << " Number of matched edges: " << nMatchedPoints
                        << abort(FatalError);
                }
            }
        }
    }
}

// Private member function to perform Laplacian surface smoothing
void mesquiteSmoother::smoothSurfaces()
{
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    Info << "Solving for surface point motion...";

    if (prioritizeVertices_)
    {
        preparePriorityList();
    }

    for (label i = 0; i < nSweeps_; i++)
    {
        // Prepare point-normals with updated point positions
        preparePointNormals();

        // Solve for surface point positions
        solveForSurfacePoints();
    }

    Info << "Done.";

    // Update coupled patches
    if (patchCoupling_.size())
    {
        forAllIter(Map<label>, patchCoupling_, patchI)
        {
            const labelList& meshPts = boundary[patchI.key()].meshPoints();

            forAll(meshPts, pointI)
            {
                refPoints_[masterToSlave_[patchI.key()][meshPts[pointI]]] =
                (
                    refPoints_[meshPts[pointI]]
                );
            }
        }
    }
}

// Solve for surface point positions
void mesquiteSmoother::solveForSurfacePoints()
{
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    forAll(pIDs_, patchI)
    {
        const labelList& meshPts = boundary[pIDs_[patchI]].meshPoints();

        // Set the current patch
        patchIndex_ = patchI;

        forAll(meshPts, pointI)
        {
            // Set the current point.
            pointIndex_ = pointI;

            // Optimize the point-position in local coordinates.
            optimizePoint();
        }
    }
}

// Optimize the point-position in local coordinates.
// Solve using the Sequential Quadratic Programming (SQP) approach.
void mesquiteSmoother::optimizePoint()
{
    vector x = vector::zero, xNew = vector::zero;

    // Evaluate the constraint on this point as well.
    scalar ceq = 0.0, f = 0.0;
    vector gnc = vector::zero, gf = vector::zero;

    // First evaluate the function and constraint
    // (including gradients) at the start point.
    evaluateFunction(2, x, f, gf);
    evaluateConstraint(2, x, ceq, gnc);

    // Start with the test point
    xNew = x;

    // Start with an identity matrix as the Hessian estimate.
    tensor H = I;

    // Set a few default options
    scalar eps = 1e-16, tolX = 1e-6, tolFn = 1e-6, tolCon = 1e-6;
    label iter = 0, nFnEvals = 1, nGradEvals = 1, verbosity = 1;
    label maxFnEvals = 300, maxIter = 400;
    vector searchDir = vector::zero, sDiff = vector::zero;
    scalar stepLength = 1.0;
    bool done = false, noFurtherProgress = false;

    // Obtain the initial constraint violation
    scalar mg = ceq, nc = ceq, ga;

    // Initialize old variables to current estimates.
    scalar ncOld = nc, fOld = f, GT = 0.0;
    scalar optErr = 0.0, feasErr = 0.0, optScal, feasScal;
    vector xOld = xNew, gfOld = vector::zero;
    vector an = vector::zero, anOld = vector::zero, matX = vector::zero;
    vector gOld = vector::zero, gNew = vector::zero, yL = vector::zero;

    // Initialize Lagrange-multipliers
    scalar lambda = 0.0, oldLambda = 0.0, Lambda = 0.0;
    scalar lambdaNLP = 0.0, newLambda = 0.0, lOld = 0.0;

    // Print out iteration info if necessary...
    if (verbosity > 0)
    {
        Info << endl;
        Info << "\tIter: "
             << "\tFn: "
             << "\tf(x): "
             << "\t\tMax Const.: "
             << "\tstepLength: "
             << "\tDir. deriv.: "
             << "\tOptimality: "
             << endl;
    }

    // Main loop
    while (!done)
    {
        // Add in the constraints...
        an = gnc;

        if (iter > 0)
        {
            optErr = cmptMax(cmptMag(gf + (lambdaNLP*an)));
            feasErr = mg;
            optScal = 1.0; feasScal = 1.0;

            // Print out iteration info if necessary...
            if (verbosity > 0)
            {
                Info << '\t' << iter << '\t' << nFnEvals
                     << '\t' << setw(6) << f << '\t'
                     << '\t' << setw(8) << mg
                     << '\t' << setw(8) << stepLength
                     << '\t' << setw(8) << (gf & searchDir)
                     << '\t' << setw(8) << optErr
                     << endl;
            }

            // Test for convergence...
            noFurtherProgress =
            (
                (
                    (cmptMax(cmptMag(searchDir)) < 2.0*tolX) ||
                    (mag(searchDir & gf) < 2.0*tolFn)
                )
                && (mg < tolCon)
            );

            if ((optErr < tolFn*optScal) && (feasErr < tolCon*feasScal))
            {
                done = true;
            }
            else
            if (noFurtherProgress)
            {
                if (mg < tolCon)
                {
                    tensor Q = tensor::zero, Z = tensor::zero;
                    vector R = vector::zero;

                    qrDecompose(an, Q, R, Z);

                    vector qTgf = (Q.T() & gf);

                    lambdaNLP = -(pInv(R) & qTgf);

                    optErr = cmptMax(cmptMag(gf + (lambdaNLP*an)));

                    optScal = 1.0;

                    if (optErr < (tolFn*optScal))
                    {
                        if (verbosity > 0)
                        {
                            Info << "Terminated because:" << endl;
                            Info << "\t First-order optimality < "
                                 << tolFn << nl
                                 << "\t and max. contraint violation < "
                                 << tolCon
                                 << endl;
                        }
                    }
                    else
                    if (cmptMax(cmptMag(searchDir)) < 2.0*tolX)
                    {
                        if (verbosity > 0)
                        {
                            Info << "Terminated because:" << endl;
                            Info << "\t | searchDir | < "
                                 << 2.0*tolX << nl
                                 << "\t and max. contraint violation < "
                                 << tolCon
                                 << endl;
                        }
                    }
                    else
                    {
                        if (verbosity > 0)
                        {
                            Info << "Terminated because:" << endl;
                            Info << "\t | grad(f) | < "
                                 << 2.0*tolFn << nl
                                 << "\t and max. contraint violation < "
                                 << tolCon
                                 << endl;
                        }
                    }
                }
                else
                {
                    if (cmptMax(cmptMag(searchDir)) < 2.0*tolX)
                    {
                        if (verbosity > 0)
                        {
                            Info << "Terminated because:" << endl;
                            Info << "No feasible solution." << endl;
                            Info << "\t | searchDir | < " << 2.0*tolX << nl
                                 << "\t but contraint was not satisfied."
                                 << endl;
                        }
                    }
                    else
                    {
                        if (verbosity > 0)
                        {
                            Info << "Terminated because:" << endl;
                            Info << "No feasible solution." << endl;
                            Info << "\t | grad(f) | < " << 2.0*tolFn << nl
                                 << "\t but contraint was not satisfied."
                                 << endl;
                        }
                    }
                }

                done = true;
            }
            else
            {
                // Reached maximum Function evaluations
                if (nFnEvals > maxFnEvals)
                {
                    xNew = matX;
                    f = fOld;
                    gf = gfOld;
                    done = true;
                }

                // Reached maximum SQP iterations
                if (iter >= maxIter)
                {
                    xNew = matX;
                    f = fOld;
                    gf = gfOld;
                    done = true;
                }
            }
        }

        // If we're not done yet, or this is the zeroth iteration:
        if (!done)
        {
            iter++;

            // Figure out the search direction.
            scalar schg = (an & gf);

            // Gradient should be opposite to the function gradient.
            if (schg > 0.0)
            {
                an = -an;
                nc = -nc;
            }

            // Check for the first evaluation.
            if (nGradEvals > 1)
            {
                newLambda = Lambda;
                gNew = gf + (newLambda*an);
                gOld = gfOld + (Lambda*anOld);
                yL = (gNew - gOld);
                sDiff = (xNew - xOld);

                // Make sure the Hessian is positive definite during update.
                if ((yL & sDiff) < (stepLength*stepLength*1e-3))
                {
                    while ((yL & sDiff) < -1e-5)
                    {
                        vector yMul = cmptMultiply(yL, sDiff);

                        if ((yMul.x() < yMul.y()) && (yMul.x() < yMul.z()))
                        {
                            yL.x() /= 2.0;
                        }
                        else
                        if ((yMul.y() < yMul.x()) && (yMul.y() < yMul.z()))
                        {
                            yL.y() /= 2.0;
                        }
                        else
                        if ((yMul.z() < yMul.x()) && (yMul.z() < yMul.y()))
                        {
                            yL.z() /= 2.0;
                        }
                    }

                    if ((yL & sDiff) < eps*sqrt(tr(H.T() & H)))
                    {
                        WarningIn
                        (
                            "mesquiteSmoother::optimizePoint()"
                        )
                            << " Algorithm requires a second Hessian update,"
                            << " but this has not been implemented."
                            << endl;
                    }
                }

                // Perform the BFGS update of the Hessian.
                if ((yL & sDiff) > eps)
                {
                    H +=
                    (
                        ((yL*yL)/(yL & sDiff))
                      - (((H & sDiff)*(H & sDiff)) / (sDiff & (H & sDiff)))
                    );
                }
            }
            else
            {
                oldLambda = (eps + magSqr(gf)) / ((an & an) + eps);
            }

            // Increment the gradient evaluation count.
            nGradEvals++;

            // Store old values
            lOld = Lambda; anOld = an; xOld = xNew;
            gfOld = gf; ncOld = nc; fOld = f;

            // Prepare inputs for the QP solver
            // for an updated search direction.
            GT = nc;
            H = 0.5*(H + H.T());

            // Reset the search-direction first.
            searchDir = vector::zero;

            // Call the QP solver.
            QPsolver(H, gf, an, -GT, searchDir, lambda);

            // Update the Lagrange Multipliers
            lambdaNLP = lambda;
            lambda = mag(lambda);
            ga = mag(nc);
            mg = ga;

            Lambda = lambda;
            oldLambda = Foam::max(Lambda, 0.5*(Lambda+oldLambda));

            // Start the line-search
            matX = xNew;
            scalar matL = f + (oldLambda*ga) + 1e-30, matL2;

            if (mg > 0)
            {
                matL2 = mg;
            }
            else
            if (f >= 0)
            {
                matL2 = -1.0 / (1.0 + f);
            }
            else
            {
                matL2 = 0.0;
            }

            if (f < 0.0)
            {
                matL2 += (f - 1.0);
            }

            scalar merit = (matL + 1), merit2 = (matL2 + 1);

            stepLength = 2.0;

            while
            (
                (merit2 > matL2) &&
                (merit > matL) &&
                (nFnEvals < maxFnEvals)
            )
            {
                stepLength /= 2.0;

                if (stepLength < 1e-4)
                {
                    stepLength = -stepLength;
                }

                xNew = matX + (stepLength*searchDir);

                // Evaluate the function value at this point.
                evaluateFunction(0, xNew, f, gf);

                nFnEvals++;

                // Evaluate the constraint at this point.
                evaluateConstraint(0, xNew, nc, gnc);
                ga = mag(nc);
                mg = ga;

                merit = f + (oldLambda*ga);

                if (mg > 0)
                {
                    merit2 = mg;
                }
                else
                if (f >= 0)
                {
                    merit2 = -1.0 / (1.0 + f);
                }
                else
                {
                    merit2 = 0.0;
                }

                if (f < 0.0)
                {
                    merit2 += (f - 1.0);
                }
            }

            scalar mf = mag(stepLength);
            Lambda = (mf*Lambda) + ((1 - mf)*lOld);

            // Now evaluate the function and constraint
            // (including gradients) at the updated point.
            evaluateFunction(2, xNew, f, gf);
            evaluateConstraint(2, xNew, ceq, gnc);

            // Update statistics
            nFnEvals++; nGradEvals++;
        }
    }

    // Now update the point-position...
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();
    const labelList& meshPts = boundary[pIDs_[patchIndex_]].meshPoints();

    // Obtain the transform tensor for this point.
    const tensor& T = T_[patchIndex_][pointIndex_].T();

    // Transform the point to global coordinates.
    refPoints_[meshPts[pointIndex_]] =
    (
        (T & xNew) + localPts_[patchIndex_][pointIndex_]
    );
}

//  QPsolver(H,f,A,b) solves the quadratic programming problem:
//
//  min 0.5*x'Hx + f'x   subject to:  Ax <= b
//   x
void mesquiteSmoother::QPsolver
(
    const tensor& H,
    const vector& f,
    const vector& Ac,
    const scalar& bc,
    vector& x,
    scalar& lambda,
    const bool normalize
)
{
    label iter = 0, maxIter = 30, type = -1;
    bool del_cstr = false;

    tensor Q = tensor::zero, Z = tensor::zero;
    vector A = Ac, actSet = vector::zero, R = vector::zero;
    scalar normf = 1, normA = 1;
    scalar n = 0.0, b = bc, rlambda = 0.0, actlambda = 0.0;

    if (normalize)
    {
        // Normalize the constraints.
        n = mag(A);
        if (n > SMALL)
        {
            A /= n;
            b /= n;
            normA = n;
        }
    }
    else
    {
        n = 1.0;
    }

    lambda = 0.0;

    // Set the active-set.
    actSet = A;

    if (normalize)
    {
        qrDecompose(actSet, Q, R, Z);

        scalar res = b - (A & x);

        vector minNormStep = (vector(Q.xx(),Q.yx(),Q.zx()) * (res/R.x()));

        x += minNormStep;
    }

    // Find the initial feasible solution.
    vector gf = (H & x) + f, sd = vector::zero;

    // Obtain the search direction.
    sd = computeSearchDirection(Z, H, gf, type);

    // Main loop
    while (iter < maxIter)
    {
        iter++;

        // Find the maximum possible step that can be taken
        // without violating constraints.
        scalar stepMin = 1e16;

        if (type == 1)
        {
            // Taking a Newton step.
            stepMin = 1.0;
            x += sd;
            del_cstr = true;
        }

        // Calculate gradient at this point.
        gf = (H & x) + f;

        if (del_cstr)
        {
            vector qTgf = (Q.T() & gf);

            // Solve for rlambda using SVD.
            rlambda = -(pInv(R) & qTgf);

            actlambda = mag(rlambda);

            if (actlambda > 0.0)
            {
                lambda = normf * (rlambda/normA);

                return;
            }
        }
    }
}

// QR Decomposition using LAPACK
void mesquiteSmoother::qrDecompose
(
    const vector& a,
    tensor& Q,
    vector& R,
    tensor& Z
)
{
    // Perform the QR decomposition of the Active Set transpose.
    int m = 3, n = 1, k = 3, lda = 3, lwork = 3, info;

    // Remember that Fortran transposes indices...
    Q = tensor
    (
        a.x(), a.y(), a.z(),
        0, 0, 0,
        0, 0, 0
    );

    R = vector::zero;
    tensor tau = tensor::zero, work = tensor::zero;

    // Call the LAPACK QR Decomposition routines
    ::dgeqrf_(&m, &n, &Q.xx(), &lda, &tau.xx(), &work.xx(), &lwork, &info);

    // Set R before Q gets overwritten.
    n = 3;
    R.x() = Q.xx(); R.y() = Q.yx(); R.z() = Q.zx();

    ::dorgqr_(&m, &n, &k, &Q.xx(), &lda, &tau.xx(), &work.xx(), &lwork, &info);

    // Store Z
    Z = tensor::zero;

    Z.xx() = Q.xy(); Z.xy() = Q.xz();
    Z.yx() = Q.yy(); Z.yy() = Q.yz();
    Z.zx() = Q.zy(); Z.zy() = Q.zz();
}

// Pseudo-inverse using LAPACK (for vectors)
vector mesquiteSmoother::pInv
(
    const vector& a
)
{
    const char jobu = 'S', jobvt = 'S';
    int m = 3, n = 1, lda = 3, ldu = 3, ldvt = 1, lwork = 6, info;

    double work[6];
    double s[1] = {0.0}, vt[1] = {0.0};
    vector u = vector::zero, cpy = a;

    dgesvd_
    (
        &jobu,
        &jobvt,
        &m, &n,
        &cpy.x(),
        &lda,
        s,
        &u.x(), &ldu,
        vt, &ldvt,
        work, &lwork,
        &info
    );

    vector x = vector::zero;

    if (s[0] > SMALL)
    {
        x = vt[0]*(1.0/s[0])*u;
    }

    return x;
}

// Pseudo-inverse using LAPACK (for 2D tensors)
tensor2D mesquiteSmoother::pInv
(
    const tensor2D& a
)
{
    const char jobu = 'A', jobvt = 'A';
    int m = 2, n = 2, lda = 2, ldu = 2, ldvt = 2, lwork = 10, info;

    double work[10];
    tensor2D sf = tensor2D::zero, vt = tensor2D::zero;
    tensor2D u = tensor2D::zero, cpy = a;

    dgesvd_
    (
        &jobu,
        &jobvt,
        &m, &n,
        &cpy.xx(), &lda,
        &sf.xx(),
        &u.xx(), &ldu,
        &vt.xx(), &ldvt,
        work, &lwork,
        &info
    );

    tensor2D vtsu = tensor2D::zero;

    if (sf.xx() > SMALL && sf.xy() > SMALL)
    {
        // Make sf diagonal.
        tensor2D s(1.0/sf.xx(), 0, 0, 1.0/sf.xy());

        // Compute the inverse from the decomposition.
        vtsu = vt & (s & u);
    }
    else
    {
        WarningIn
        (
            "mesquiteSmoother::pInv(tensor2D& a)"
        )
            << " Returning zero inverse."
            << endl;
    }

    return vtsu;
}

// Compute a search direction vector
vector mesquiteSmoother::computeSearchDirection
(
    const tensor& Z,
    const tensor& H,
    const vector& gf,
    label& type
)
{
    vector sd = vector::zero;
    tensor projH = Z.T() & H & Z;

    // Copy element values for factorization
    tensor2D a(projH.xx(), projH.yx(), projH.xy(), projH.yy());

    char factType = 'L';
    int n = 2, lda = 2, info = -1;

    // Compute the Cholesky factorization.
    dpotrf_(&factType, &n, &a.xx(), &lda, &info);

    if (!info)
    {
        // Factorization is positive definite. Take a Newton direction.
        vector zTgf = (Z.T() & gf);

        // Copy values for the matrix solve.
        vector2D b(zTgf.x(), zTgf.y());

        // Explicitly zero-out the off-diagonal
        a.yx() = 0;

        // c = inv(a) * (inv(a') * b)
        // The transpose is for output from LAPACK
        vector2D c = pInv(a).T() & (pInv(a.T()).T() & b);

        // Specify the type as Newton.
        type = 1;

        sd = -1.0*(Z & vector(c.x(), c.y(), 0));
    }
    else
    {
        WarningIn
        (
            "mesquiteSmoother::computeSearchDirection()"
        )
            << " cholTrap has not been implemented."
            << endl;
    }

    // Make sure that this is a descent direction
    if ((gf & sd) > 0.0)
    {
        sd *= -1.0;
    }

    return sd;
}

// Compute the accumulated cell metrics at a point.
void mesquiteSmoother::evaluateFunction
(
    const label type,
    const vector& x,
    scalar& fn,
    vector& fnGrad
)
{
    scalar f = 0.0;
    vector fGrad = vector::zero;

    // Obtain the set of ordered cell-points for the Jacobian.
    const labelList& jList = jPoints_[patchIndex_][pointIndex_];

    // Obtain the transform tensor for this point.
    const tensor& T = T_[patchIndex_][pointIndex_];

    // Fetch the current position for this point.
    const point& xi = localPts_[patchIndex_][pointIndex_];

    // Create a temporary tensor whose columns are the test point.
    tensor xT
    (
        x.x(), x.x(), x.x(),
        x.y(), x.y(), x.y(),
        x.z(), x.z(), x.z()
    );

    for(label i = 0; i < jList.size(); i += 3)
    {
        const point& xm0 = refPoints_[jList[i+0]];
        const point& xm1 = refPoints_[jList[i+1]];
        const point& xm2 = refPoints_[jList[i+2]];

        // Prepare a tensor in the local coordinate system,
        // whose origin lies at point 'xi'.
        tensor jL
        (
            xm0.x() - xi.x(), xm1.x() - xi.x(), xm2.x() - xi.x(),
            xm0.y() - xi.y(), xm1.y() - xi.y(), xm2.y() - xi.y(),
            xm0.z() - xi.z(), xm1.z() - xi.z(), xm2.z() - xi.z()
        );

        // Prepare the Jacobian.
        tensor J = (T & jL) - xT;

        // Obtain the determinant...
        scalar alpha = det(J);

        // Obtain the Frobenius norm
        scalar fro = tr(J.T() & J);

        // Obtain the denominator of the metric
        scalar den = cbrt(alpha * alpha);

        // Compute the element metric and accumulate.
        f = f + (fro/den);

        // Check if a gradient if being requested
        if (type > 1)
        {
            // Compute the analytical gradient
            vector gfMat =
            (
                (2.0/(3.0*den))*((3.0*J) - (fro*inv(J).T()))
              & vector::one
            );

            fGrad += gfMat;
        }
    }

    // Set return value based on type
    if (type == 0)
    {
        fn = f;
    }
    else
    if (type == 1)
    {
        fnGrad = (-1.0*fGrad);
    }
    else
    if (type == 2)
    {
        fn = f;
        fnGrad = (-1.0*fGrad);
    }
}

// Evaluate the point constraint for a given point.
void mesquiteSmoother::evaluateConstraint
(
    const label type,
    const vector& x,
    scalar& cstr,
    vector& grad
)
{
    // Fetch the coefficients for this point.
    const scalarField& a = a_[patchIndex_][pointIndex_];

    if (type == 0 || type == 2)
    {
        cstr =
        (
            a[0]
          + a[1]*x.x()
          + a[2]*x.y()
          + a[3]*x.x()*x.y()
          + a[4]*x.x()*x.x()
          + a[5]*x.y()*x.y()
          - x.z()
        );
    }

    if (type == 1 || type == 2)
    {
        grad.x() = a[1] + a[3]*x.y() + 2.0*a[4]*x.x();
        grad.y() = a[2] + a[3]*x.x() + 2.0*a[5]*x.y();
        grad.z() = -1.0;
    }
}

// Prepare a priority list of surface vertices for smoothing.
void mesquiteSmoother::preparePriorityList()
{
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    // Dummy gradient
    vector gf = vector::zero;

    forAll(pIDs_, patchI)
    {
        const labelList& meshPts = boundary[pIDs_[patchI]].meshPoints();

        patchIndex_ = patchI;

        meanFnValues_[patchI] = 0.0;

        forAll(meshPts, pointI)
        {
            pointIndex_ = pointI;

            // Evaluate the objective function for this point
            evaluateFunction
            (
                0,
                refPoints_[meshPts[pointI]],
                pFnValues_[patchI][pointI],
                gf
            );

            meanFnValues_[patchI] += pFnValues_[patchI][pointI];
        }

        // Evaluate the mean...
        meanFnValues_[patchI] /= meshPts.size();
    }
}

// Prepare point-normals with updated point positions.
//  - Use this new information to set-up a local coordinate system,
//    such that the z-axis corresponds to the point normal.
//  - Compute the coefficients for a parametric parabolic surface / curve
//    which passes through the point and its nearest neighbours.
void mesquiteSmoother::preparePointNormals()
{
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    // Define the global coordinate system vectors for convenience.
    vector x1(1,0,0);
    vector y1(0,1,0);
    vector z1(0,0,1);

    // Allocate a coefficient matrix for inversion.
    scalarMatrix A(6);

    forAll(pIDs_, patchI)
    {
        // First update localPoints with latest point positions
        const labelList& meshPts = boundary[pIDs_[patchI]].meshPoints();

        forAll(meshPts,pointI)
        {
            localPts_[patchI][pointI] = refPoints_[meshPts[pointI]];
        }

        // Now compute point normals from updated local points
        const labelListList& pFaces = boundary[pIDs_[patchI]].pointFaces();
        const labelListList& pEdges = boundary[pIDs_[patchI]].pointEdges();

        const edgeList& edges = boundary[pIDs_[patchI]].edges();
        const faceList& faces = boundary[pIDs_[patchI]].localFaces();

        pNormals_[patchI] = vector::zero;

        // Test patch planarity on-the-fly
        bool planePatch = true;
        vector cumVector = vector::zero;

        forAll(pNormals_[patchI], pointI)
        {
            vector& n = pNormals_[patchI][pointI];

            forAll(pFaces[pointI], faceI)
            {
                n += faces[pFaces[pointI][faceI]].normal(localPts_[patchI]);
            }

            // Normalize the vector
            n /= mag(n) + VSMALL;

            // Accumulate the vector and normalize it.
            if (planePatch)
            {
                cumVector += n;

                cumVector /= mag(cumVector) + VSMALL;

                if ((cumVector & n) < 0.99)
                {
                    // Fails the planarity test.
                    planePatch = false;
                }
            }
        }

        // If the patch is not planar, prepare surface coefficients
        // for each point on the patch.
        if (planePatch)
        {
            // Compute and assign a transform tensor for
            // all points on this patch.

            // Pick the normal at the first point.
            vector z2 = pNormals_[patchI][0];

            // Choose a vector with arbitrary orientation
            vector v(0.5, 0.5, 0.5);

            // Prepare the x-axis
            vector x2 = v - (v & z2)*z2;
            x2 /= mag(x2) + VSMALL;

            // Prepare the y-axis
            vector y2 = (x2 ^ z2);
            y2 /= mag(y2) + VSMALL;

            // Prepare the transform tensor for this point.
            tensor T = tensor::zero;

            T.xx() = (x1 & x2); T.xy() = (x1 & y2); T.xz() = (x1 & z2);
            T.yx() = (y1 & x2); T.yy() = (y1 & y2); T.yz() = (y1 & z2);
            T.zx() = (z1 & x2); T.zy() = (z1 & y2); T.zz() = (z1 & z2);

            // Assign the transform tensor / coefficients for all points.
            forAll(meshPts, pointI)
            {
                T_[patchI][pointI] = T;

                // The initArrays() function already assigns null-fields
                // for these points, so they won't be initialized here.
            }
        }
        else
        {
            label targetCount = 8;

            labelList coeffList(targetCount);
            pointField coeffPoints(targetCount, vector::zero);

            // Least-squares variables
            scalarListList lsqCoeffs(targetCount, scalarList(6));
            scalarField lsqSource(targetCount, 0.0);
            scalarField lsqWeights(targetCount, 0.0);

            forAll(meshPts, pointI)
            {
                // Avoid points on bounding curves for now.
                if (pFlags_[patchI][pointI])
                {
                    continue;
                }

                // Construct a local coordinate system at this point,
                // with the z-axis as the normal. The x-axis and y-axis
                // can be arbitrary.
                vector z2 = pNormals_[patchI][pointI];

                // Choose a vector with arbitrary orientation
                vector v(0.5, 0.5, 0.5);

                // Prepare the x-axis
                vector x2 = v - (v & z2)*z2;
                x2 /= mag(x2) + VSMALL;

                // Prepare the y-axis
                vector y2 = (x2 ^ z2);
                y2 /= mag(y2) + VSMALL;

                // Prepare the transform tensor for this point.
                // This tensor is attractive because it is orthogonal.
                // i.e., its inverse is also its transpose.
                tensor& T = T_[patchI][pointI];

                T.xx() = (x1 & x2); T.xy() = (x1 & y2); T.xz() = (x1 & z2);
                T.yx() = (y1 & x2); T.yy() = (y1 & y2); T.yz() = (y1 & z2);
                T.zx() = (z1 & x2); T.zy() = (z1 & y2); T.zz() = (z1 & z2);

                // Loop through edges connected to this point,
                // Attempt to obtain 6 points for a parabolic surface.
                const point& origin = localPts_[patchI][pointI];
                const labelList& peList = pEdges[pointI];

                // Reset the point list
                coeffList = -1;

                // The first point is always the origin.
                label nPoints = 0;

                lsqWeights[nPoints] = GREAT;
                coeffList[nPoints++] = pointI;

                forAll(peList, edgeI)
                {
                    // Get the other point on this edge.
                    label other = edges[peList[edgeI]].otherVertex(pointI);

                    const point& oPoint = refPoints_[meshPts[other]];

                    // Translate / rotate to the local origin.
                    coeffList[nPoints] = other;
                    coeffPoints[nPoints] = T & (oPoint - origin);

                    // Assign weights
                    lsqWeights[nPoints] = GREAT;
                    /*
                    (
                        1.0
                      / magSqr(oPoint - origin)
                    );
                    */

                    nPoints++;

                    // Break out if we're done.
                    if (nPoints == targetCount)
                    {
                        break;
                    }
                }

                while (nPoints < targetCount)
                {
                    // Need to find a few more points.
                    // Look at existing points and their neighbours.
                    for (label i = 1; i < nPoints; i++)
                    {
                        // Is this a point on a bounding curve?
                        // Don't use it.
                        if (pFlags_[patchI][coeffList[i]])
                        {
                            continue;
                        }

                        const labelList& opeList = pEdges[coeffList[i]];

                        forAll(opeList, edgeI)
                        {
                            // Get the other point on this edge.
                            label other =
                            (
                                edges[opeList[edgeI]].otherVertex(coeffList[i])
                            );

                            const point& oPoint = refPoints_[meshPts[other]];

                            if (findIndex(coeffList, other) == -1)
                            {
                                // Translate / rotate to the local origin.
                                coeffList[nPoints] = other;
                                coeffPoints[nPoints] = T & (oPoint - origin);

                                // Assign weights
                                lsqWeights[nPoints] = GREAT / 5.0;
                                /*
                                (
                                    1.0
                                  / magSqr(oPoint - origin)
                                );
                                */

                                nPoints++;

                                // Break out if we're done.
                                if (nPoints == targetCount)
                                {
                                    break;
                                }
                            }
                        }

                        // Have we achieved the requirement?
                        if (nPoints == targetCount)
                        {
                            break;
                        }
                    }
                }

                // Fetch the coefficient list.
                scalarField& b = a_[patchI][pointI];

                // Prepare the least-squares coefficients
                for (label i = 0; i < nPoints; i++)
                {
                    scalar x = coeffPoints[i][0];
                    scalar y = coeffPoints[i][1];
                    scalar z = coeffPoints[i][2];

                    lsqCoeffs[i][0] = 1.0;
                    lsqCoeffs[i][1] = x;
                    lsqCoeffs[i][2] = y;
                    lsqCoeffs[i][3] = x * y;
                    lsqCoeffs[i][4] = x * x;
                    lsqCoeffs[i][5] = y * y;
                    lsqSource[i] = z;
                }

                // Prepare the matrix
                prepareLeastSquaresMatrix
                (
                    nPoints,
                    lsqCoeffs,
                    lsqWeights,
                    lsqSource,
                    A,
                    b
                );

                /*
                // Prepare matrix coefficients and source.
                //  - Direct solution using only 6 points.
                for (label i = 0; i < 6; i++)
                {
                    scalar x = coeffPoints[i][0];
                    scalar y = coeffPoints[i][1];
                    scalar z = coeffPoints[i][2];

                    A[i][0] = 1.0;
                    A[i][1] = x;
                    A[i][2] = y;
                    A[i][3] = x * y;
                    A[i][4] = x * x;
                    A[i][5] = y * y;
                    b[i] = coeffPoints[i][2];
                }
                */

                // Solve for coefficients. Solution is stored in 'b'.
                scalarMatrix::solve(A,b);
            }
        }
    }
}

// Utility method to prepare the least-squares matrix
// for either a parabolic curve / surface fit.
void mesquiteSmoother::prepareLeastSquaresMatrix
(
    const label nPoints,
    const scalarListList& lsqCoeffs,
    const scalarField& lsqWeights,
    const scalarField& lsqSource,
    scalarMatrix& lsqMatrix,
    scalarField& source
)
{
    // Obtain number of columns (coefficients)
    label nCols = lsqCoeffs[0].size();

    // Since matrix is symmetric, only perform upper-diagonal
    for (label i = -1; i < nCols; i++)
    {
        for (label j = i+1; j < nCols; j++)
        {
            lsqMatrix[i+1][j] = 0.0;

            for (label k = 0; k < nPoints; k++)
            {
                lsqMatrix[i+1][j] +=
                (
                    lsqWeights[k]*lsqCoeffs[k][i+1]*lsqCoeffs[k][j]
                );
            }

            // Copy to lower diagonal.
            lsqMatrix[j][i+1] = lsqMatrix[i+1][j];
        }
    }

    // Perform the matrix-vector multiply
    for (label i = 0; i < nCols; i++)
    {
        source[i] = 0.0;

        for (label k = 0; k < nPoints; k++)
        {
            source[i] += (lsqWeights[k]*lsqCoeffs[k][i]*lsqSource[k]);
        }
    }
}

tmp<pointField> mesquiteSmoother::newPoints()
{
    solve();

    return curPoints();
}

//- Return point location obtained from the current motion field
tmp<pointField> mesquiteSmoother::curPoints() const
{
    tmp<pointField> tcurPoints(refPoints_);

    return tcurPoints;
}

void mesquiteSmoother::solve()
{
    // Perform surface smoothing first
    if
    (
        surfaceSmoothing_
     && (Mesh_.time().timeIndex() % surfInterval_ == 0)
     && (Mesh_.time().timeIndex() != 0)
    )
    {
        smoothSurfaces();
    }

    // Copy most recent point positions
    forAll(refPoints_, pointI)
    {
        vtxCoords_[(3*pointI)+0] = refPoints_[pointI][0];
        vtxCoords_[(3*pointI)+1] = refPoints_[pointI][1];
        vtxCoords_[(3*pointI)+2] = refPoints_[pointI][2];
    }

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
    Mesquite::QualityAssessor qA
    (
        &qMetricTable_[qMetric_](),
        err
    );

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

    // Launches optimization on the mesh
    queue.run_instructions(&msqMesh, err);

    // Copy updated positions back
    forAll(refPoints_, pointI)
    {
        refPoints_[pointI][0] = vtxCoords_[(3*pointI)+0];
        refPoints_[pointI][1] = vtxCoords_[(3*pointI)+1];
        refPoints_[pointI][2] = vtxCoords_[(3*pointI)+2];
    }
}

void mesquiteSmoother::updateMesh(const mapPolyMesh& mpm)
{
    motionSolver::updateMesh(mpm);

    // Clear Mesquite arrays
    delete [] vtxCoords_;
    delete [] cellToNode_;
    delete [] fixFlags_;

    nPoints_ = Mesh_.nPoints();
    nCells_  = Mesh_.nCells();

    // Reset refPoints
    refPoints_.clear();
    refPoints_ = Mesh_.points();

    if (surfaceSmoothing_)
    {
        masterToSlave_.clear();

        pNormals_.clear();
        localPts_.clear();
        jPoints_.clear();
        T_.clear();
        a_.clear();
        pFlags_.clear();

        if (prioritizeVertices_)
        {
            meanFnValues_.clear();
            pFnValues_.clear();
        }
    }

    // Initialize data structures
    initArrays();
}

} // End namespace Foam

// ************************************************************************* //
