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
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mesquiteSmoother::mesquiteSmoother
(
    const polyMesh& mesh
)
:
    motionSolver(mesh),
    Mesh_(mesh),
    nPoints_(mesh.nPoints()),
    nCells_(mesh.nCells()),
    surfaceSmoothing_(false),
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

Foam::mesquiteSmoother::mesquiteSmoother
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

Foam::mesquiteSmoother::~mesquiteSmoother()
{
    delete [] vtxCoords_;
    delete [] cellToNode_;
    delete [] fixFlags_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Read options from the dictionary
void Foam::mesquiteSmoother::readOptions()
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

    // Set termination criteria for untangling
    //untInner_.add_absolute_quality_improvement(0.0);
    //untInner_.add_absolute_successive_improvement(1e-4);
    //untOuter_.add_iteration_limit(1);

    //untangleMetric_.set(new Mesquite::UntangleBetaQualityMetric(1e-8));
    //untangleFunc_.set(new Mesquite::LPtoPTemplate(&untangleMetric_(), 2, err));
    //untangleGlobal_.set(new Mesquite::ConjugateGradient(&untangleFunc_(),err));
    //untangleGlobal_->use_global_patch();
}

// Initialize connectivity arrays for Mesquite
void Foam::mesquiteSmoother::initArrays()
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
    }
}

// Sparse matrix-vector multiply [3D]
void Foam::mesquiteSmoother::A
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
void Foam::mesquiteSmoother::applyBCs
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
        label nFix = (field.size()*10)/100;

        for(label i = 0; i < nFix; i++)
        {
            field[randomizer.integer(0, field.size()-1)] = vector::zero;
        }
    }
}

// Vector dot-product
Foam::scalar Foam::mesquiteSmoother::dot
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

Foam::scalar Foam::mesquiteSmoother::normFactor
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
Foam::scalar Foam::mesquiteSmoother::cmptSumMag
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
Foam::label Foam::mesquiteSmoother::CG
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
void Foam::mesquiteSmoother::smoothSurfaces()
{
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

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

// Prepare point-normals with updated point positions
void Foam::mesquiteSmoother::preparePointNormals()
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

Foam::tmp<Foam::pointField> Foam::mesquiteSmoother::newPoints()
{
    solve();

    return curPoints();
}

//- Return point location obtained from the current motion field
Foam::tmp<Foam::pointField>
Foam::mesquiteSmoother::curPoints() const
{
    tmp<pointField> tcurPoints(refPoints_);

    return tcurPoints;
}

void Foam::mesquiteSmoother::solve()
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

    // Apply termination criteria for untangling
    //untangleGlobal_->set_inner_termination_criterion(&untInner_);
    //untangleGlobal_->set_outer_termination_criterion(&untOuter_);

    //untangleGlobal_->loop_over_mesh(&msqMesh, 0, 0, err);

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

void Foam::mesquiteSmoother::updateMesh(const mapPolyMesh& mpm)
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

// ************************************************************************* //
