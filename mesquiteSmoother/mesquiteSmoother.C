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
#include "Random.H"
#include "IOmanip.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"

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
    volumeCorrection_(false),
    volCorrTolerance_(1e-20),
    volCorrMaxIter_(100),
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
    ),
    oldVolume_(0.0)
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
    volumeCorrection_(false),
    volCorrTolerance_(1e-20),
    volCorrMaxIter_(100),
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
    ),
    oldVolume_(0.0)
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

        labelHashSet slipPatchIDs;

        forAll(slipPatches, wordI)
        {
            word& patchName = slipPatches[wordI];

            slipPatchIDs.insert(mesh().boundaryMesh().findPatchID(patchName));
        }

        // Extract the patch list
        pIDs_ = slipPatchIDs.toc();

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

        // Check if volume correction is enabled
        if (found("volumeCorrection"))
        {
            volumeCorrection_ = lookup("volumeCorrection");
        }

        // Check if volume correction tolerance is specified
        if (found("volCorrTolerance"))
        {
            volCorrTolerance_ = readScalar(lookup("volCorrTolerance"));
        }

        // Check if volume correction maxIter is specified
        if (found("volCorrMaxIter"))
        {
            volCorrMaxIter_ = readLabel(lookup("volCorrMaxIter"));
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
                if (findIndex(pIDs_, sPatch) > -1)
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

// Initialize connectivity arrays for Mesquite
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
        offsets_.setSize(pIDs_.size() + 1, 0);
        pNormals_.setSize(pIDs_.size());
        gradEdgeV_.setSize(pIDs_.size());
        localPts_.setSize(pIDs_.size());

        label totalSize = 0;

        forAll(pIDs_, patchI)
        {
            label nPts = boundary[pIDs_[patchI]].nPoints();
            label nEdg = boundary[pIDs_[patchI]].nEdges();

            pNormals_[patchI].setSize(nPts, vector::zero);
            gradEdgeV_[patchI].setSize(nEdg, vector::zero);
            localPts_[patchI].setSize(nEdg, vector::zero);

            // Accumulate the total size
            totalSize += nPts;

            // Set offsets
            offsets_[patchI + 1] = totalSize;
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

        // Initialize CG variables
        bV_.setSize(totalSize, vector::zero);
        xV_.setSize(totalSize, vector::zero);
        pV_.setSize(totalSize, vector::zero);
        rV_.setSize(totalSize, vector::zero);
        wV_.setSize(totalSize, vector::zero);

        origPoints_.setSize(refPoints_.size(), vector::zero);
    }
}

// Sparse matrix-vector multiply [3D]
void mesquiteSmoother::A
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
            gradEdgeV_[patchI][edgeI] =
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
            w[edges[edgeI][0] + offsets_[patchI]] += gradEdgeV_[patchI][edgeI];
            w[edges[edgeI][1] + offsets_[patchI]] -= gradEdgeV_[patchI][edgeI];
        }
    }

    // Apply boundary conditions
    applyBCs(w);
}

// Apply boundary conditions
void mesquiteSmoother::applyBCs
(
    vectorField& field
)
{
    // Blank out residuals at boundary nodes
    label nFixedBC = 0;
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    forAll(pIDs_, patchI)
    {
        const edgeList& edges = boundary[pIDs_[patchI]].edges();

        // Apply slip conditions for internal nodes
        forAll(pNormals_[patchI], pointI)
        {
            const vector& n = pNormals_[patchI][pointI];

            field[pointI + offsets_[patchI]] -=
            (
                (field[pointI + offsets_[patchI]] & n)*n
            );
        }

        // Blank out residuals for bounding curves
        for
        (
            label i = boundary[pIDs_[patchI]].nInternalEdges();
            i < edges.size();
            i++
        )
        {
            field[edges[i][0] + offsets_[patchI]] = vector::zero;
            field[edges[i][1] + offsets_[patchI]] = vector::zero;

            nFixedBC++;
        }
    }

    // If no boundaries were fixed, fix a few points at random
    if (nFixedBC == 0)
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
scalar mesquiteSmoother::dot
(
    const vectorField& f1,
    const vectorField& f2
)
{
    scalar s = 0.0;

    forAll(f1, indexI)
    {
        s += (f1[indexI] & f2[indexI]);
    }

    return s;
}

scalar mesquiteSmoother::normFactor
(
    const vectorField& x,
    const vectorField& b,
    const vectorField& w,
    vectorField& tmpField
)
{
    vector xRef = average(x);

    A(vectorField(x.size(), xRef),tmpField);

    vectorField nFw = (w - tmpField);
    vectorField nFb = (b - tmpField);

    return cmptSumMag(nFw) + cmptSumMag(nFb) + 1.0e-20;
}

// Component-wise sumMag
scalar mesquiteSmoother::cmptSumMag
(
    const vectorField& field
)
{
    scalar cSum = 0.0;

    forAll(field,i)
    {
        cSum += mag(field[i].x()) + mag(field[i].y()) + mag(field[i].z());
    }

    return cSum;
}

// CG solver
label mesquiteSmoother::CG
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

    r = b - w;
    p = r;
    rho = dot(r,p);

    // Obtain the normalized residual
    residual = cmptSumMag(r)/norm;

    Info << " Initial residual: " << residual;

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

    Info << " Final residual: " << residual;

    return iter;
}

// Private member function to perform Laplacian surface smoothing
void mesquiteSmoother::smoothSurfaces()
{
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    // Copy refPoints prior to all surface-smoothing sweeps.
    origPoints_ = refPoints_;

    // Value to be used later in volume error correction
    if (volumeCorrection_)
    {
        oldVolume_ = gSum(mesh().cellVolumes());
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

        Info << "Solving for point motion: ";

        label iters = CG(bV_, pV_, rV_, wV_, xV_);

        Info << " No Iterations: " << iters << endl;

        // Update refPoints
        forAll(pIDs_, patchI)
        {
            const labelList& meshPts = boundary[pIDs_[patchI]].meshPoints();

            forAll(meshPts,pointI)
            {
                refPoints_[meshPts[pointI]] = xV_[pointI + offsets_[patchI]];
            }
        }

        // Update coupled patches
        if (patchCoupling_.size())
        {
            forAllIter(Map<label>, patchCoupling_, patchI)
            {
                const labelList& meshPts = boundary[patchI.key()].meshPoints();

                forAll(meshPts,pointI)
                {
                    refPoints_[masterToSlave_[patchI.key()][meshPts[pointI]]] =
                    (
                        refPoints_[meshPts[pointI]]
                    );
                }
            }
        }
    }
}

// Find the volume of a tetrahedron.
// The function assumes points (a-b-c)
// are in counter-clockwise fashion when viewed from d.
inline scalar mesquiteSmoother::tetVolume
(
    const label cIndex
)
{
    const cell& cellToCheck = mesh().cells()[cIndex];

    const face& currFace = mesh().faces()[cellToCheck[0]];
    const face& nextFace = mesh().faces()[cellToCheck[1]];

    // Get the fourth point and compute cell volume
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
                const point& a = refPoints_[currFace[2]];
                const point& b = refPoints_[currFace[1]];
                const point& c = refPoints_[currFace[0]];
                const point& d = refPoints_[nextFace[pointI]];

                return ((1.0/6.0)*(((b - a) ^ (c - a)) & (d - a)));
            }
            else
            {
                const point& a = refPoints_[currFace[0]];
                const point& b = refPoints_[currFace[1]];
                const point& c = refPoints_[currFace[2]];
                const point& d = refPoints_[nextFace[pointI]];

                return ((1.0/6.0)*(((b - a) ^ (c - a)) & (d - a)));
            }
        }
    }

    // Something's wrong with connectivity.
    FatalErrorIn("mesquiteSmoother::tetVolume()")
        << "Cell: " << cIndex
        << " has inconsistent connectivity."
        << abort(FatalError);

    return 0.0;
}

// Private member function to check for invalid
// cells and correct if necessary.
void mesquiteSmoother::correctInvalidCells()
{
    // Loop through pointCells for all boundary points
    // and compute cell volume.
    const labelListList& pointCells = mesh().pointCells();
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    // Obtain point-positions after smoothing
    pointField oldField = refPoints_;

    DynamicList<label> invCells(50);

    forAll(pIDs_, patchI)
    {
        const labelList& meshPts = boundary[pIDs_[patchI]].meshPoints();

        forAll(meshPts, pointI)
        {
            const labelList& pCells = pointCells[meshPts[pointI]];

            forAll(pCells, cellI)
            {
                if (tetVolume(pCells[cellI]) < 0.0)
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

    InfoIn("mesquiteSmoother::correctInvalidCells()")
        << "Found " << invCells.size() << " invalid cells. "
        << "Attempting to correct..." << flush;

    // Looks like some inverted cells are present.
    // Check for free-vertices connected to invalid cells.
    const labelListList& cellPoints = mesh().cellPoints();
    const labelListList& pointFaces = mesh().pointFaces();

    labelHashSet freePoints;

    forAll(invCells, cellI)
    {
        const labelList& cPoints = cellPoints[invCells[cellI]];

        forAll(cPoints, pointI)
        {
            const labelList& pFaces = pointFaces[cPoints[pointI]];

            bool boundaryPoint = false;

            forAll(pFaces, faceI)
            {
                if (!mesh().isInternalFace(pFaces[faceI]))
                {
                    boundaryPoint = true;
                    break;
                }
            }

            if (!boundaryPoint && !freePoints.found(cPoints[pointI]))
            {
                freePoints.insert(cPoints[pointI]);
            }
        }
    }

    // Obtain mesh connectivity
    const faceList& meshFaces = mesh().faces();
    const cellList& meshCells = mesh().cells();
    const labelList& owner = mesh().faceOwner();

    // For each point, build a local hull for optimization
    labelList freePointList = freePoints.toc();
    labelListList jPoints(freePointList.size(), labelList(0));

    forAll(freePointList, pointI)
    {
        const labelList& pCells = pointCells[freePointList[pointI]];

        // First size the list. Assume tet-mesh.
        jPoints[pointI].setSize(3*pCells.size());

        label i = 0;

        forAll(pCells, cellI)
        {
            const cell& curCell = meshCells[pCells[cellI]];

            // Find a face which doesn't contain pIndex.
            forAll(curCell, faceI)
            {
                const face& thisFace = meshFaces[curCell[faceI]];

                if (thisFace.which(freePointList[pointI]) == -1)
                {
                    if (owner[curCell[faceI]] == pCells[cellI])
                    {
                        jPoints[pointI][i++] = thisFace[0];
                        jPoints[pointI][i++] = thisFace[1];
                        jPoints[pointI][i++] = thisFace[2];
                    }
                    else
                    {
                        jPoints[pointI][i++] = thisFace[2];
                        jPoints[pointI][i++] = thisFace[1];
                        jPoints[pointI][i++] = thisFace[0];
                    }

                    break;
                }
            }
        }
    }

    if (debug)
    {
        forAll(freePointList, pointI)
        {
            writePoint(freePointList[pointI], jPoints[pointI]);
        }
    }

    // Loop through various points, and check for a negative
    // Jacobian value. If an inverted cell exists, proceed
    // with the LBFGS algorithm.
    int N = 3;
    scalar beta = 0.0;
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *m_x = lbfgs_malloc(N);

    forAll(freePointList, pointI)
    {
        const vector& x = refPoints_[freePointList[pointI]];

        bool foundInvalid = checkValidity(x, jPoints[pointI], beta);

        if (foundInvalid)
        {
            // Prepare the node-position
            m_x[0] = x.x();
            m_x[1] = x.y();
            m_x[2] = x.z();

            if (debug)
            {
                Info << "x0: "
                     << m_x[0] << " "
                     << m_x[1] << " "
                     << m_x[2] << " "
                     << endl;
            }

            // Prepare data
            optInfo data
            (
                (*this),
                beta,
                jPoints[pointI],
                refPoints_
            );

            // Call the lbfgs algorithm
            lbfgs
            (
                N,
                m_x,
                &fx,
                _evaluate,
                NULL,
                reinterpret_cast<void *>(&data),
                NULL
            );

            if (debug)
            {
                Info << "xNew: "
                     << m_x[0] << " "
                     << m_x[1] << " "
                     << m_x[2] << " "
                     << endl;
            }

            foundInvalid = checkValidity
            (
                vector(m_x[0],m_x[1],m_x[2]),
                jPoints[pointI],
                beta
            );

            // Update position only if untangling was successful
            if (!foundInvalid)
            {
                refPoints_[freePointList[pointI]] =
                (
                    vector(m_x[0],m_x[1],m_x[2])
                );
            }
        }
    }

    // Free allocated memory
    lbfgs_free(m_x);

    // Perform a final check to ensure everything went okay.
    bool foundInvalid = false;

    forAll(freePointList, pointI)
    {
        const vector& x = refPoints_[freePointList[pointI]];

        foundInvalid = checkValidity(x, jPoints[pointI], beta);

        if (foundInvalid)
        {
            break;
        }
    }

    if (foundInvalid)
    {
        WarningIn
        (
            "mesquiteSmoother::correctInvalidCells()"
        )
            << " Failed to untangle mesh."
            << endl;

        Info << "Relaxing points to untangle mesh...";
    }
    else
    {
        Info << "Success." << endl;

        return;
    }

    // Clear out the invalid cell list, and obtain a new list
    invCells.clear();

    forAll(pIDs_, patchI)
    {
        const labelList& meshPts = boundary[pIDs_[patchI]].meshPoints();

        forAll(meshPts, pointI)
        {
            const labelList& pCells = pointCells[meshPts[pointI]];

            forAll(pCells, cellI)
            {
                if (tetVolume(pCells[cellI]) < 0.0)
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

    bool valid = false;
    scalar lambda = 2.0;
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
            (lambda * oldField)
          + ((1.0 - lambda) * origPoints_)
        );

        forAll(invCells, cellI)
        {
            if (tetVolume(invCells[cellI]) < 0.0)
            {
                valid = false;

                break;
            }
        }

        // Stop further testing if invalid.
        if (!valid)
        {
            break;
        }

        nAttempts++;

        if (nAttempts > 50)
        {
            WarningIn("mesquiteSmoother::correctInvalidCells()")
                    << " Failed to obtain an untangled mesh." << nl
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
void mesquiteSmoother::correctGlobalVolume()
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
        domainVolume += tetVolume(cellI);
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
            domainVolume += tetVolume(cellI);
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
                const labelList& meshPts = boundary[pIDs_[patchI]].meshPoints();

                forAll(meshPts,pointI)
                {
                    refPoints_[meshPts[pointI]] +=
                    (
                        magVal*pNormals_[patchI][pointI]
                    );
                }
            }

            domainVolume = 0;

            forAll(allCells,cellI)
            {
                domainVolume += tetVolume(cellI);
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

            forAll(meshPts,pointI)
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
void mesquiteSmoother::enforceCylindricalConstraints()
{
    if (!surfaceSmoothing_)
    {
        return;
    }

    // Check for sub-dictionary entry
    if (found("cylindricalConstraints"))
    {
        const dictionary& constraintDict = subDict("cylindricalConstraints");

        // Read patch-information one-by one.
        wordList cstrPatches = constraintDict.toc();

        forAll(cstrPatches, wordI)
        {
            label pID = mesh().boundaryMesh().findPatchID(cstrPatches[wordI]);

            if (pID == -1 || findIndex(pIDs_, pID) == -1)
            {
                FatalErrorIn
                (
                    "mesquiteSmoother::enforceCylindricalConstraints()"
                )
                    << " Cannot find patch: " << cstrPatches[wordI]
                    << abort(FatalError);
            }

            const dictionary& pD = constraintDict.subDict(cstrPatches[wordI]);

            // Read info.
            vector axisPoint(pD.lookup("axisPoint"));
            vector axisVector(pD.lookup("axisVector"));
            scalar radius = readScalar(pD.lookup("radius"));

            const labelList meshPts = mesh().boundaryMesh()[pID].meshPoints();

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
                            "mesquiteSmoother::enforceCylindricalConstraints()"
                        )
                            << " Constraint violation: " << viol << endl;
                    }
                }
            }
        }
    }
}

// Utility method to check validity of cells connected to a point.
bool mesquiteSmoother::checkValidity
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

// Static function for callback
lbfgsfloatval_t mesquiteSmoother::_evaluate
(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
)
{
    // Recast the argument
    optInfo *data =
    (
        reinterpret_cast<optInfo*>(instance)
    );

    // Evaluate the function first
    lbfgsfloatval_t fnVal =
    (
        data->motionSolver().evaluate
        (
            data->orderedPoints(),
            data->referencePoints(),
            data->beta(),
            x,
            g
        )
    );

    return fnVal;
}

// Actual function evaulation routine
lbfgsfloatval_t mesquiteSmoother::evaluate
(
    const labelList& jList,
    const pointField& refPoints,
    const scalar beta,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g
) const
{
    // Initialize the function and gradient
    lbfgsfloatval_t fx = 0.0;
    vector fxGrad = vector::zero;

    for (label i = 0; i < jList.size(); i += 3)
    {
        const point& xm0 = refPoints_[jList[i+0]];
        const point& xm1 = refPoints_[jList[i+1]];
        const point& xm2 = refPoints_[jList[i+2]];

        // Prepare the Jacobian.
        tensor J
        (
            xm0.x() - x[0], xm1.x() - x[0], xm2.x() - x[0],
            xm0.y() - x[1], xm1.y() - x[1], xm2.y() - x[1],
            xm0.z() - x[2], xm1.z() - x[2], xm2.z() - x[2]
        );

        // Obtain the determinant...
        scalar alpha = det(J);
        scalar amb = (alpha - beta);
        scalar u = (mag(amb) - amb);

        // Compute the element metric and accumulate.
        fx = fx + (u*u);

        // Accumulate the gradient as well.
        fxGrad += ((((-2*u*u*alpha) / mag(amb))*inv(J).T()) & vector::one);
    }

    // Negate the gradient
    g[0] = -fxGrad.x(); g[1] = -fxGrad.y(); g[2] = -fxGrad.z();

    return fx;
}

// Write a particular point out for post-processing
void mesquiteSmoother::writePoint
(
    const label pointIndex,
    const labelList& jList
)
{
    // Locally order cell-points and write them out.
    Map<label> pointMap;

    label nPoints = 0, nCells = 0;
    labelListList cpList(jList.size() / 3, labelList(4, -1));
    pointField points(jList.size(), vector::zero);

    // Zeroth point is the point under consideration
    points[nPoints++] = refPoints_[pointIndex];

    for (label i = 0; i < jList.size(); i += 3)
    {
        label pointJ = 0;

        if (!pointMap.found(jList[i+0]))
        {
            points[nPoints] = refPoints_[jList[i+0]];
            pointMap.insert(jList[i+0], nPoints);
            nPoints++;
        }

        cpList[nCells][pointJ++] = pointMap[jList[i+0]];

        if (!pointMap.found(jList[i+1]))
        {
            points[nPoints] = refPoints_[jList[i+1]];
            pointMap.insert(jList[i+1], nPoints);
            nPoints++;
        }

        cpList[nCells][pointJ++] = pointMap[jList[i+1]];

        if (!pointMap.found(jList[i+2]))
        {
            points[nPoints] = refPoints_[jList[i+2]];
            pointMap.insert(jList[i+2], nPoints);
            nPoints++;
        }

        cpList[nCells][pointJ++] = pointMap[jList[i+2]];

        // The last point is the point under consideration
        cpList[nCells][pointJ++] = 0;

        nCells++;
    }

    // Correct the number of points
    points.setSize(nPoints);

    writeVTK
    (
        "pointCells_"
      + Foam::name(mesh().time().timeIndex()) + '_'
      + Foam::name(pointIndex),
        mesh(),
        points,
        cpList
    );
}

// Output a list of points / cells as a VTK file.
//  - CellPoints are required to be ordered.
void mesquiteSmoother::writeVTK
(
    const word& name,
    const polyMesh& mesh,
    const pointField& points,
    const labelListList& cellPoints
)
{
    // Make the directory
    // fileName dirName(mesh.time().path()/"VTK"/mesh.time().timeName());
    fileName dirName(mesh.time().path()/"VTK");

    mkDir(dirName);

    // Open stream for output
    OFstream file(dirName/name+".vtk");

    // Write out the header
    file << "# vtk DataFile Version 2.0" << nl
         << name << ".vtk" << nl
         << "ASCII" << nl
         << "DATASET UNSTRUCTURED_GRID" << nl
         << "POINTS " << points.size() << " double" << nl;

    forAll(points, i)
    {
        file << setprecision(10)
             << points[i].x() << ' '
             << points[i].y() << ' '
             << points[i].z() << ' '
             << nl;
    }

    if (cellPoints.size())
    {
        file << "CELLS " << cellPoints.size()
             << " " << (4*cellPoints.size()) + cellPoints.size()
             << endl;

        forAll(cellPoints, cellI)
        {
            file << 4 << ' '
                 << cellPoints[cellI][0] << ' '
                 << cellPoints[cellI][1] << ' '
                 << cellPoints[cellI][2] << ' '
                 << cellPoints[cellI][3] << ' '
                 << nl;
        }

        file << "CELL_TYPES " << cellPoints.size() << endl;

        forAll(cellPoints, cellI)
        {
            // Tetrahedron
            file << "10" << nl;
        }
    }
    else
    {
        file << "CELLS " << points.size() << " " << 2*points.size() << endl;

        forAll(points, i)
        {
            file << 1 << ' ' << i << ' ' << nl;
        }

        file << "CELL_TYPES " << points.size() << endl;

        forAll(points, i)
        {
            // Vertex
            file << "1" << nl;
        }
    }
}

// Prepare point-normals with updated point positions
void mesquiteSmoother::preparePointNormals()
{
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
        const labelListList& pFaces = boundary[pIDs_[patchI]].pointFaces();
        const faceList& faces = boundary[pIDs_[patchI]].localFaces();

        pNormals_[patchI] = vector::zero;

        forAll(pNormals_[patchI], pointI)
        {
            vector& n = pNormals_[patchI][pointI];

            forAll(pFaces[pointI], faceI)
            {
                n += faces[pFaces[pointI][faceI]].normal(localPts_[patchI]);
            }

            // Normalize the vector
            n /= mag(n) + VSMALL;
        }
    }
}

tmp<pointField> mesquiteSmoother::newPoints()
{
    solve();

    return curPoints();
}

//- Return point location obtained from the current motion field
tmp<pointField>
mesquiteSmoother::curPoints() const
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

        // Check for invalid cells and correct if necessary.
        correctInvalidCells();

        // Correct for volume change due to smoothing
        correctGlobalVolume();

        // Enforce constraints, if necessary
        enforceCylindricalConstraints();
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
        // Clear out CG variables
        bV_.clear();
        xV_.clear();
        pV_.clear();
        rV_.clear();
        wV_.clear();

        masterToSlave_.clear();

        localPts_.clear();
        gradEdgeV_.clear();
        pNormals_.clear();
        offsets_.clear();
    }

    // Initialize data structures
    initArrays();
}

} // End namespace Foam

// ************************************************************************* //
