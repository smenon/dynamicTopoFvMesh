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

Class
    dynamicTopoFvMesh

Description
    An implementation of dynamic changes to mesh-topology

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*----------------------------------------------------------------------------*/

#include "HashList.H"
#include "clockTime.H"
#include "dynamicTopoFvMesh.H"
#include "dynamicTopoFvMeshMapper.H"
#include "multiThreader.H"
#include "tetDecompositionMotionSolver.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"
#include "motionSolver.H"
#include "mapPolyMesh.H"

#include "fvPatchFields.H"
#include "fvsPatchFields.H"
#include "MapFvFields.H"
#include "MeshObject.H"

#include <dlfcn.h>

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynamicTopoFvMesh,0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dynamicTopoFvMesh::dynamicTopoFvMesh(const IOobject& io)
:
    fvMesh(io),
    numPatches_(polyMesh::boundaryMesh().size()),
    topoChangeFlag_(false),
    dict_
    (
        IOobject
        (
            "dynamicMeshDict",
            polyMesh::time().constant(),
            (*this),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    twoDMesh_(polyMesh::nGeometricD() == 2 ? true : false),
    edgeModification_
    (
        dict_.subDict("dynamicTopoFvMesh").lookup("edgeModification")
    ),
    solveForMotion_
    (
        dict_.subDict("dynamicTopoFvMesh").lookup("solveForMotion")
    ),
    mapper_(NULL),
    meshPoints_(polyMesh::points()),
    faces_(polyMesh::faces()),
    owner_(polyMesh::faceOwner()),
    neighbour_(polyMesh::faceNeighbour()),
    cells_(primitiveMesh::cells()),
    IOedges_
    (
        IOobject
        (
            "edges",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    IOpointEdges_
    (
        IOobject
        (
            "pointEdges",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    IOedgeFaces_
    (
        IOobject
        (
            "edgeFaces",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    IOfaceEdges_
    (
        IOobject
        (
            "faceEdges",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    IOedgePatchStarts_
    (
        IOobject
        (
            "edgePatchStarts",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    IOedgePatchSizes_
    (
        IOobject
        (
            "edgePatchSizes",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    oldPatchSizes_(numPatches_,0),
    patchSizes_(numPatches_,0),
    oldPatchStarts_(numPatches_,-1),
    patchStarts_(numPatches_,-1),
    oldEdgePatchSizes_(numPatches_,0),
    edgePatchSizes_(numPatches_,0),
    oldEdgePatchStarts_(numPatches_,-1),
    edgePatchStarts_(numPatches_,-1),
    oldPatchNMeshPoints_(numPatches_,-1),
    patchNMeshPoints_(numPatches_,-1),
    nOldPoints_(primitiveMesh::nPoints()),
    nPoints_(primitiveMesh::nPoints()),
    nOldEdges_(0),
    nEdges_(0),
    nOldFaces_(primitiveMesh::nFaces()),
    nFaces_(primitiveMesh::nFaces()),
    nOldCells_(primitiveMesh::nCells()),
    nCells_(primitiveMesh::nCells()),
    nOldInternalFaces_(primitiveMesh::nInternalFaces()),
    nInternalFaces_(primitiveMesh::nInternalFaces()),
    nOldInternalEdges_(0),
    nInternalEdges_(0),
    ratioMin_(0.0),
    ratioMax_(0.0),
    growthFactor_(1.0),
    maxLengthScale_(GREAT),
    sliverThreshold_(0.05),
    bisectInteriorFace_(-1),
    maxTetsPerEdge_(-1),
    allowTableResize_(false)
{
    // For backward compatibility, check the size of owner/neighbour
    if (owner_.size() != neighbour_.size())
    {
        // Append -1 for neighbours
        for(label i = nInternalFaces_; i < nFaces_; i++)
        {
            neighbour_.append(-1);
        }
    }

    // Initialize the motion-solver, if it was requested
    if (solveForMotion_)
    {
        motionPtr_.set(motionSolver::New(*this).ptr());
    }

    // Initialize the multiThreading environment
    if (dict_.subDict("dynamicTopoFvMesh").found("threads"))
    {
        threader_.set
        (
            new multiThreader
            (
                readLabel
                (
                    dict_.subDict("dynamicTopoFvMesh").lookup("threads")
                )
            )
        );
    }
    else
    {
        threader_.set(new multiThreader(1));
    }

    // Get the number of threads and allocate topoMeshStructures
    label nThreads = threader_->getNumThreads();
    structPtr_ = new topoMeshStruct[nThreads];
    for (label i=0; i<nThreads; i++)
    {
        structPtr_[i].mesh_ = this;
        structPtr_[i].nThreads_ = nThreads;
        structPtr_[i].pthreadID_ = threader_->getID(i);
    }

    // Size the stacks
    if (twoDMesh_)
    {
        faceStack_.setSize(nThreads);
    }
    else
    {
        edgeStack_.setSize(nThreads);
    }


    // Initialize patch-size information
    const polyBoundaryMesh& boundary = polyMesh::boundaryMesh();
    for(label i=0; i<numPatches_; i++)
    {
        patchNMeshPoints_[i] = boundary[i].meshPoints().size();
        oldPatchSizes_[i]  = patchSizes_[i]  = boundary[i].size();
        oldPatchStarts_[i] = patchStarts_[i] = boundary[i].start();
        oldPatchNMeshPoints_[i] = patchNMeshPoints_[i];
    }

    // Set sizes for the reverse maps
    reversePointMap_.setSize(nPoints_);
    reverseFaceMap_.setSize(nFaces_);
    reverseCellMap_.setSize(nCells_);

    // For tetrahedral meshes...
    if (!twoDMesh_)
    {
        // Open the tetMetric dynamic-link library
        void * metricLibPtr = NULL;
        char * error;

        if
        (
            dict_.subDict("dynamicTopoFvMesh").found("tetMetricLib")
        )
        {
            metricLibPtr =
                dlopen
                (
                    word
                    (
                        dict_.subDict("dynamicTopoFvMesh").lookup("tetMetricLib")
                    ).c_str(),
                    RTLD_LAZY|RTLD_GLOBAL
                );
        }
        else
        {
            metricLibPtr = dlopen("libtetMetrics.so", RTLD_LAZY|RTLD_GLOBAL);
        }

        if (!metricLibPtr)
        {
            FatalErrorIn
            (
                "dynamicTopoFvMesh::dynamicTopoFvMesh(const IOobject& io) "
            ) << nl
            << " Could not open the tetMetric library. "
            << abort(FatalError);
        }

        // Obtain the tetrahedral metric to be used.
        word tetMetric
        (
            dict_.subDict("dynamicTopoFvMesh").lookup("tetMetric")
        );

        // Obtain function address from the dll.
        tetMetric_ =
            reinterpret_cast<tetMetricReturnType>
            (
                dlsym
                (
                    metricLibPtr,
                    tetMetric.c_str()
                )
            );

        if ((error = dlerror()) != NULL)
        {
            FatalErrorIn
            (
                "dynamicTopoFvMesh::dynamicTopoFvMesh(const IOobject& io) "
            ) << nl
            << " Unrecognized tet-quality metric: " << tetMetric
            << " Reported dlsym() error: " << error
            << abort(FatalError);
        }

        // Check if swapping is to be avoided on any patches
        if (dict_.subDict("dynamicTopoFvMesh").found("noSwapPatches"))
        {
            wordList noSwapPatches =
                dict_.subDict
                (
                    "dynamicTopoFvMesh"
                ).subDict("noSwapPatches").toc();

            forAll(noSwapPatches, wordI)
            {
                word& patchName = noSwapPatches[wordI];

                forAll(boundaryMesh(), patchI)
                {
                    if (boundaryMesh()[patchI].name() == patchName)
                    {
                        noSwapPatchIDs_.insert(patchI);
                    }
                }
            }
        }

        // Check if a limit has been imposed on maxTetsPerEdge
        if (dict_.subDict("dynamicTopoFvMesh").found("maxTetsPerEdge"))
        {
            maxTetsPerEdge_ =
                readLabel
                (
                    dict_.subDict("dynamicTopoFvMesh").lookup("maxTetsPerEdge")
                );
        }
        else
        {
            maxTetsPerEdge_ = 7;
        }

        // Check if programming tables can be resized at runtime
        if (dict_.subDict("dynamicTopoFvMesh").found("allowTableResize"))
        {
            allowTableResize_ =
                readBool
                (
                    dict_.subDict("dynamicTopoFvMesh").lookup("allowTableResize")
                );
        }

        // Initialize edge-related connectivity structures
        initEdges();
    }

    // Define edgeModification options
    if (edgeModification_)
    {
        const dictionary& edgeOptionDict =
            dict_.subDict("dynamicTopoFvMesh").subDict("edgeOptions");
        ratioMax_ = readScalar(edgeOptionDict.lookup("bisectionRatio"));
        ratioMin_ = readScalar(edgeOptionDict.lookup("collapseRatio"));
        growthFactor_ = readScalar(edgeOptionDict.lookup("growthFactor"));

        if (edgeOptionDict.found("maxLengthScale"))
        {
            maxLengthScale_ =
                readScalar(edgeOptionDict.lookup("maxLengthScale"));
        }

        if (edgeOptionDict.found("fixedLengthScalePatches"))
        {
            fixedLengthScalePatches_ =
                edgeOptionDict.subDict("fixedLengthScalePatches");
        }

        if (edgeOptionDict.found("sliverThreshold"))
        {
            sliverThreshold_ =
                readScalar(edgeOptionDict.lookup("sliverThreshold"));
        }

        initLengthScale();
    }

    // Initialize mutex lists for multi-threading, if necessary
    initMutexLists();
}

// Default topoMeshStruct constructor
dynamicTopoFvMesh::topoMeshStruct::topoMeshStruct()
{
    mesh_ = 0;
    nThreads_ = pthreadID_ = -1;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dynamicTopoFvMesh::~dynamicTopoFvMesh()
{
    delete [] structPtr_;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Access a particular boundary displacement patch
void dynamicTopoFvMesh::setMotionBC
(
    const label& index,
    const vectorField& dispField
)
{
    if (solveForMotion_)
    {
        // Determine the kind of motion solver in use
        word solverType(dict_.lookup("solver"));

        //- Cell decomposition FEM motion solver
        if
        (
            (solverType == "laplaceCellDecomposition")
         || (solverType == "pseudoSolidCellDecomposition")
        )
        {
            // Boundary motion specified for the tetDecompositionMotionSolver
            tetPointVectorField& motionU = const_cast<tetPointVectorField&>
                (objectRegistry::lookupObject<tetPointVectorField>("motionU"));

            // Assign boundary conditions to the motion solver
            motionU.boundaryField()[index] == dispField/time().deltaT().value();
        }

        //- Face decomposition FEM motion solver
        if
        (
            (solverType == "laplaceFaceDecomposition")
         || (solverType == "pseudoSolidFaceDecomposition")
        )
        {
            // Boundary motion specified for the tetDecompositionMotionSolver
            tetPointVectorField& motionU = const_cast<tetPointVectorField&>
                (objectRegistry::lookupObject<tetPointVectorField>("motionU"));

            // Assign boundary conditions to the motion solver

            // The face-decomposition solver includes points at face-centres,
            // thus point motion has to be interpolated to these points
            tetPolyPatchInterpolation interpolator
            (
                refCast<const faceTetPolyPatch>
                (
                    motionU.boundaryField()[index].patch()
                )
            );

            motionU.boundaryField()[index] ==
                interpolator.pointToPointInterpolate
                (
                    dispField/time().deltaT().value()
                );
        }

        //- Spring-based Laplacian motion solver
        if
        (
            (solverType == "springMotionSolver")
        )
        {
            // Boundary motion specified for the springMotionSolver
            pointField& refPoints = const_cast<pointField&>
                (objectRegistry::lookupObject<pointField>("refPoints"));

            // Assign boundary conditions to the motion solver
            const labelList& meshPts = boundaryMesh()[index].meshPoints();
            forAll(meshPts,pointI)
            {
                refPoints[meshPts[pointI]] += dispField[pointI];
            }
        }
    }
}

// Return the mesh-mapper
const mapPolyMesh& dynamicTopoFvMesh::meshMap()
{
    if (mapper_.valid())
    {
        return mapper_();
    }
    else
    {
        FatalErrorIn("dynamicTopoFvMesh::meshMap() ") << nl
                << " Illegal request for the mesh mapper. " << nl
                << abort(FatalError);
    }

    return mapper_();
}

// Return old cell-centre information (prior to a topology change)
const vectorField& dynamicTopoFvMesh::oldCellCentres() const
{
    if (cellCentresPtr_.valid())
    {
        return cellCentresPtr_();
    }
    else
    {
        FatalErrorIn("dynamicTopoFvMesh::oldCellCentres() ") << nl
                << " Illegal request for old cell centres. " << nl
                << abort(FatalError);
    }

    return vectorField::null();
}

// Return mesh length-scale values
tmp<scalarField> dynamicTopoFvMesh::lengthScale()
{
    tmp<scalarField> tlengthScale
    (
        new scalarField(nCells(), 0.0)
    );

    scalarField& internalField = tlengthScale();

    // Obtain length-scale values from the mesh
    forAllIter(HashList<scalar>, lengthScale_, lIter)
    {
        internalField[lIter.index()] = lIter();
    }

    return tlengthScale;
}

// Return mesh cell-quality values
// Valid for 3D tetrahedral meshes only...
tmp<scalarField> dynamicTopoFvMesh::meshQuality()
{
    tmp<scalarField> tQuality
    (
        new scalarField(nCells(), 0.0)
    );

    // Valid for 3D tetrahedral meshes only...
    if (!twoDMesh_)
    {
        scalarField& iF = tQuality();

        // Compute statistics on the fly
        scalar maxQuality = -GREAT;
        scalar minQuality =  GREAT;
        scalar meanQuality = 0.0;

        const pointField& meshPoints = points();
        const faceList& meshFaces = faces();
        const cellList& meshCells = cells();

        const labelList& owner = faceOwner();

        // Loop through all cells in the mesh and compute cell quality
        forAll(meshCells, cellI)
        {
            const cell& curCell = meshCells[cellI];
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
                    // Compute cell-quality and write-out
                    if (owner[curCell[0]] == cellI)
                    {
                        iF[cellI] = (*tetMetric_)
                        (
                            meshPoints[currFace[2]],
                            meshPoints[currFace[1]],
                            meshPoints[currFace[0]],
                            meshPoints[nextFace[pointI]]
                        );
                    }
                    else
                    {
                        iF[cellI] = (*tetMetric_)
                        (
                            meshPoints[currFace[0]],
                            meshPoints[currFace[1]],
                            meshPoints[currFace[2]],
                            meshPoints[nextFace[pointI]]
                        );
                    }

                    // Update statistics
                    maxQuality = iF[cellI] > maxQuality ? iF[cellI] : maxQuality;
                    minQuality = iF[cellI] < minQuality ? iF[cellI] : minQuality;
                    meanQuality += iF[cellI];

                    break;
                }
            }
        }

        // Output statistics:
        if (debug)
        {
            Info << " ~~~ Mesh Quality Statistics ~~~ " << endl;
            Info << " Min: " << minQuality << endl;
            Info << " Max: " << maxQuality << endl;
            Info << " Mean: " << meanQuality/iF.size() << endl;
            Info << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << endl;
        }
    }

    return tQuality;
}

// Find the circumcenter, given three points
inline vector dynamicTopoFvMesh::circumCenter
(
    const point& a,
    const point& b,
    const point& c,
    const label& one,
    const label& two,
    const label& three
)
{
    scalar d1 =  (c - a)&(b - a);
    scalar d2 = -(c - b)&(b - a);
    scalar d3 =  (c - a)&(c - b);

    scalar c1 = d2*d3;
    scalar c2 = d3*d1;
    scalar c3 = d1*d2;

    scalar cd = c1 + c2 + c3;

#   ifdef FULLDEBUG
    if (cd < VSMALL)
    {
        FatalErrorIn
        (
            "dynamicTopoFvMesh::circumCenter(point& a, point& b, point& c) "
        ) << nl << " Encountered a co-linear set of points: " << nl
                << " Point a :: " << one << ": " << a << nl
                << " Point b :: " << two << ": " << b << nl
                << " Point c :: " << three << ": " << c << nl
                << abort(FatalError);
    }
#   endif

    return ((c2 + c3)*a + (c3 + c1)*b + (c1 + c2)*c)/(2*cd);
}

// Find the area of a triangle face.
// This function also assumes face right-handedness
inline scalar dynamicTopoFvMesh::triFaceArea
(
    const face& triFace
)
{
    return Foam::mag(triFaceNormal(triFace));
}

// Find the normal of a triangle face.
// This function also assumes face right-handedness
inline vector dynamicTopoFvMesh::triFaceNormal
(
    const face& triFace
)
{
    vector v = meshPoints_[triFace[1]] - meshPoints_[triFace[0]];
    vector w = meshPoints_[triFace[2]] - meshPoints_[triFace[0]];

    return 0.5 * (v ^ w);
}

// Find the volume of a tetrahedron. This function assumes proper orientation
// of the vertices, and will give negative values otherwise.
inline scalar dynamicTopoFvMesh::tetVolume
(
    const point& a,
    const point& b,
    const point& c,
    const point& d
)
{
    return (1.0/6.0)*(((b - a) ^ (c - a)) & (d - a));
}

// Method to determine the old boundary patch index for a given face
// Similar to the polyBoundaryMesh routine, but works on local information
inline label dynamicTopoFvMesh::whichPatch
(
    const label& index
) const
{
    if (index < nOldInternalFaces_) return -1;

    for(label i=0; i<numPatches_; i++)
    {
        if
        (
            index >= oldPatchStarts_[i]
         && index < oldPatchStarts_[i] + oldPatchSizes_[i]
        )
        {
            return i;
        }
    }

    // If not in any of the above, it's possible that the face was added
    // at the end of the list. Check addedFacePatches_ for the patch info
    if (addedFacePatches_.found(index))
    {
        return addedFacePatches_[index];
    }
    else
    {
        FatalErrorIn
        (
            "label dynamicTopoFvMesh::whichPatch(const label& index) const"
        )
        << "Cannot find patch information for face index " << index
        << " It appears that face ordering is inconsistent with patch information."
        << abort(FatalError);
    }

    return -2;
}

// Method to determine the old boundary patch index for a given edge
inline label dynamicTopoFvMesh::whichEdgePatch
(
    const label& index
) const
{
    if (index < nOldInternalEdges_) return -1;

    for(label i=0; i<numPatches_; i++)
    {
        if
        (
            index >= oldEdgePatchStarts_[i]
         && index < oldEdgePatchStarts_[i] + oldEdgePatchSizes_[i]
        )
        {
            return i;
        }
    }

    // If not in any of the above, it's possible that the edge was added
    // at the end of the list. Check addedEdgePatches_ for the patch info
    if (addedEdgePatches_.found(index))
    {
        return addedEdgePatches_[index];
    }
    else
    {
        FatalErrorIn
        (
            "label dynamicTopoFvMesh::whichEdgePatch(const label& index) const"
        )
        << "Cannot find patch information for edge index " << index
        << " It appears that edge ordering is inconsistent with patch information."
        << abort(FatalError);
    }

    return -2;
}

// Utility method to find the interior/boundary faces
// for an input quad-face and adjacent triangle-prism cell.
inline void dynamicTopoFvMesh::findPrismFaces
(
    const label findex,
    const label cindex,
    FixedList<face,2>& bdyf,
    FixedList<label,2>& bidx,
    FixedList<face,2>& intf,
    FixedList<label,2>& iidx
)
{
    label indexO=0, indexI=0;

    cell& c = cells_[cindex];

    forAll(c, i)
    {
        label faceIndex = c[i];
        // Don't count the face under consideration
        if (faceIndex != findex)
        {
            face& fi=faces_[faceIndex];
            if (neighbour_[faceIndex] == -1)
            {
                if (fi.size() == 3)
                {
                    // Triangular face on the boundary
                    bidx[indexO] = faceIndex;
                    bdyf[indexO++] = fi;
                }
                else
                {
                    // This seems to be a non-triangular face on the boundary
                    // Consider this as "interior" and move on
                    iidx[indexI] = faceIndex;
                    intf[indexI++] = fi;
                }
            }
            else
            {
                // Face on the interior
                iidx[indexI] = faceIndex;
                intf[indexI++] = fi;
            }
        }
    }
}

// Utility method to find the common edge between two faces.
// If an edge is found, returns the common edge on the first face in the argument
bool dynamicTopoFvMesh::findCommonEdge
(
    const face& first,
    const face& second,
    edge& common
)
{
    bool found=false;
    edgeList efi = first.edges();
    edgeList efj = second.edges();
    forAll(efi, edgeI)
    {
        forAll(efj, edgeJ)
        {
            if (efi[edgeI] == efj[edgeJ])
            {
                common = efi[edgeI];
                found = true; break;
            }
        }
        if (found) break;
    }
    return found;
}

// Utility method to find the isolated point on a triangular face
// that doesn't lie on the specified edge. Also returns the point next to it.
inline void dynamicTopoFvMesh::findIsolatedPoint
(
    const face& f,
    const edge& e,
    label& ptIndex,
    label& nextPtIndex
)
{
    bool found = false;
    forAll(f, pointI)
    {
        if ( f[pointI] != e[0] && f[pointI] != e[1] )
        {
            ptIndex = f[pointI];
            nextPtIndex = f[(pointI+1)%3];
            found = true;
            break;
        }
    }
    if (!found)
    {
        FatalErrorIn
        (
            "label dynamicTopoFvMesh::findIsolatedPoint()"
        )
        << "Cannot find isolated point in face " << f << endl
        << " Using edge: " << e
        << abort(FatalError);
    }
}

// Utility method to replace a label in a given list
inline void dynamicTopoFvMesh::replaceLabel
(
     const label& original,
     const label& replacement,
     labelList& list
)
{
    bool found = false;
    forAll(list, indexI)
    {
        if (list[indexI] == original)
        {
            list[indexI] = replacement;
            found = true;
            break;
        }
    }
    if (!found)
    {
        FatalErrorIn
        (
            "label dynamicTopoFvMesh::replaceLabel()"
        )   << "Cannot find " << original << " in list: " << list << endl
            << " Label: " << replacement << " was not used in replacement."
            << abort(FatalError);
    }
}

// Utility method for face-insertion
label dynamicTopoFvMesh::insertFace
(
    const label patch,
    const face& newFace,
    const label newOwner,
    const label newNeighbour,
    const edge& edgeToWatch = edge(0,0)
)
{
    // Append the specified face to each face-related list.
    // This will avoid rehashing of existing structures, but ordering is not
    // maintained. Reordering is performed after all pending changes
    // (flips, bisections, contractions, etc) have been made to the mesh
    label newFaceIndex = faces_.append(newFace);
    owner_.append(newOwner);
    neighbour_.append(newNeighbour);

    // Add a new unlocked face mutex
    if (threader_->multiThreaded())
    {
        faceMutex_.append();
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        Info << "Inserting face: "
             << newFaceIndex << ": "
             << newFace << endl;
    }
#   endif

    if (twoDMesh_)
    {
        // Push this onto the stack as well
        faceStack_[getThreadID(pthread_self())].push(newFaceIndex);

        if (edgeModification_)
        {
            edgeToWatch_.append(edgeToWatch);
        }
    }

    // Keep track of added boundary faces in a separate hash-table
    // This information will be required at the reordering stage
    addedFacePatches_.insert(newFaceIndex,patch);

    if (newNeighbour == -1)
    {
        // Modify patch information for this boundary face
        patchSizes_[patch]++;
        for(label i=patch+1; i<numPatches_; i++)
        {
            patchStarts_[i]++;
        }
    }
    else
    {
        // Increment the number of internal faces, and subsequent patch-starts
        nInternalFaces_++;
        for(label i=0; i<numPatches_; i++)
        {
            patchStarts_[i]++;
        }
    }

    // Increment the total face count
    nFaces_++;

    return newFaceIndex;
}

// Remove the specified face from the mesh
void dynamicTopoFvMesh::removeFace
(
    const label index
)
{
#   ifdef FULLDEBUG
    if (debug)
    {
        Info << "Removed face: "
             << index << ": "
             << faces_[index] << endl;
    }
#   endif

    faces_.remove(index);
    owner_.remove(index);

    // Remove the face mutex
    if (threader_->multiThreaded())
    {
        // Unlock it first
        faceMutex_[index].unlock();
        faceMutex_.remove(index);
    }

    if (neighbour_[index] == -1)
    {
        // Modify patch information for this boundary face
        label rmFacePatch = whichPatch(index);
        patchSizes_[rmFacePatch]--;
        for(label i=rmFacePatch+1; i<numPatches_; i++)
        {
            patchStarts_[i]--;
        }
    }
    else
    {
        // Decrement the internal face count, and subsequent patch-starts
        nInternalFaces_--;
        for(label i=0; i<numPatches_; i++)
        {
            patchStarts_[i]--;
        }
    }

    neighbour_.remove(index);

    if (twoDMesh_)
    {
        // Remove from the stack as well
        forAll(faceStack_, stackI)
        {
            faceStack_[stackI].remove(index);
        }

        if (edgeModification_)
        {
            edgeToWatch_.remove(index);
        }
    }

    // Update the reverse face-map, but only if this is a face that existed
    // at time [n]. Added faces which are deleted during the topology change
    // needn't be updated.
    if (index < nOldFaces_)
    {
        reverseFaceMap_[index] = -1;
    }

    // Decrement the total face-count
    nFaces_--;
}

// Insert the specified edge to the mesh
label dynamicTopoFvMesh::insertEdge
(
    const label patch,
    const edge& newEdge,
    const labelList& edgeFaces
)
{
    label newEdgeIndex = edges_.append(newEdge);
    edgeFaces_.append(edgeFaces);

    // Add a new unlocked edge mutex
    if (threader_->multiThreaded())
    {
        edgeMutex_.append();
    }

    // Add to the stack as well
    edgeStack_[getThreadID(pthread_self())].push(newEdgeIndex);

#   ifdef FULLDEBUG
    if (debug)
    {
        Info << "Inserting edge: "
             << newEdgeIndex << ": "
             << newEdge << endl;
    }
#   endif

    // Keep track of added edges in a separate hash-table
    // This information will be required at the reordering stage
    addedEdgePatches_.insert(newEdgeIndex,patch);

    if (patch >= 0)
    {
        // Modify patch information for this boundary edge
        edgePatchSizes_[patch]++;
        for(label i=patch+1; i<numPatches_; i++)
        {
            edgePatchStarts_[i]++;
        }
    }
    else
    {
        // Increment the number of internal edges, and subsequent patch-starts
        nInternalEdges_++;
        for(label i=0; i<numPatches_; i++)
        {
            edgePatchStarts_[i]++;
        }
    }

    // Size-up the pointEdges list
    sizeUpList(newEdgeIndex, pointEdges_[newEdge[0]]);
    sizeUpList(newEdgeIndex, pointEdges_[newEdge[1]]);

    // Increment the total edge count
    nEdges_++;

    return newEdgeIndex;
}

// Remove the specified edge from the mesh
void dynamicTopoFvMesh::removeEdge
(
    const label index
)
{
    edge& thisEdge = edges_[index];

#   ifdef FULLDEBUG
    if (debug)
    {
        Info << "Removing edge: "
             << index << ": "
             << thisEdge << endl;
    }
#   endif

    // Size-down the pointEdges list
    if (pointEdges_[thisEdge[0]].size())
    {
        sizeDownList(index, pointEdges_[thisEdge[0]]);
    }
    else
    {
        pointEdges_.remove(thisEdge[0]);
    }

    if (pointEdges_[thisEdge[1]].size())
    {
        sizeDownList(index, pointEdges_[thisEdge[1]]);
    }
    else
    {
        pointEdges_.remove(thisEdge[1]);
    }

    edges_.remove(index);
    edgeFaces_.remove(index);

    // Remove the edge mutex
    if (threader_->multiThreaded())
    {
        // Unlock it first
        edgeMutex_[index].unlock();
        edgeMutex_.remove(index);
    }

    // Remove from the stack as well
    forAll(edgeStack_, stackI)
    {
        edgeStack_[stackI].remove(index);
    }

    // Identify the patch for this edge
    label patch = whichEdgePatch(index);

    if (patch >= 0)
    {
        // Modify patch information for this boundary edge
        edgePatchSizes_[patch]--;
        for(label i=patch+1; i<numPatches_; i++)
        {
            edgePatchStarts_[i]--;
        }
    }
    else
    {
        // Decrement the internal edge count, and subsequent patch-starts
        nInternalEdges_--;
        for(label i=0; i<numPatches_; i++)
        {
            edgePatchStarts_[i]--;
        }
    }

    // Update reverse edge-map, but only if this is an edge that existed
    // at time [n]. Added edges which are deleted during the topology change
    // needn't be updated.
    if (index < nOldEdges_)
    {
        reverseEdgeMap_[index] = -1;
    }

    // Decrement the total edge-count
    nEdges_--;
}

// Utility method to size-up the list to include an item
inline void dynamicTopoFvMesh::sizeUpList
(
    const label item,
    labelList& list
)
{
    // Create a new list
    labelList newList(list.size()+1);

    // Copy individual items
    forAll(list, itemI)
    {
        newList[itemI] = list[itemI];
    }

    // Set the last element and overwrite
    newList[list.size()] = item;
    list.transfer(newList);
}

// Utility method to size-down the list to remove an item
inline void dynamicTopoFvMesh::sizeDownList
(
    const label item,
    labelList& list
)
{
    // Create a new list
    labelList newList(list.size()-1);

    // Copy individual items
    label n = 0;
    forAll(list, itemI)
    {
        if (list[itemI] != item)
        {
            newList[n++] = list[itemI];
        }
    }

    // Overwrite
    list.transfer(newList);
}

// Utility method to build a hull of faces/cells that are connected to the edge
// This will also determine whether the edge lies on a boundary
bool dynamicTopoFvMesh::constructPrismHull
(
    const edge& edgeToCheck,
    const label startFaceIndex,
    DynamicList<label>& hullFaces,
    DynamicList<label>& hullCells,
    DynamicList<label>& hullTriFaces,
    DynamicList<label>& edgePatches,
    bool requiresTriFaces
)
{
    // Get the two cells on either side...
    label c0 = owner_[startFaceIndex], c1 = neighbour_[startFaceIndex];

    bool isBoundary=false, foundQuadFace, foundTriFace;
    label faceToExclude, cellIndex;

    // Start a search from cell[0] and add to the list as we go along
    faceToExclude=startFaceIndex, cellIndex=c0, hullCells.append(c0);
    do
    {
        cell& cellToCheck = cells_[cellIndex];
        foundQuadFace = false; foundTriFace = false;
        forAll(cellToCheck,faceI)
        {
            if (cellToCheck[faceI] != faceToExclude)
            {
                face& faceToCheck = faces_[cellToCheck[faceI]];
                if (faceToCheck.nEdges() == 4 && !foundQuadFace)
                {
                    edgeList indexEdges = faceToCheck.edges();
                    forAll(indexEdges,edgeI)
                    {
                        if (indexEdges[edgeI] == edgeToCheck)
                        {
                            // Bingo... We have a match. Add to the dynamic list
                            hullFaces.append(cellToCheck[faceI]);
                            faceToExclude = cellToCheck[faceI];
                            foundQuadFace=true; break;
                        }
                    }
                }

                if
                (
                    requiresTriFaces && faceToCheck.nEdges() == 3
                 && !foundTriFace
                )
                {
                    hullTriFaces.append(cellToCheck[faceI]);
                    foundTriFace=true;
                }
            }

            // Found the faces we were looking for, break-out
            if (requiresTriFaces)
            {
                if (foundQuadFace && foundTriFace) break;
            }
            else
            {
                if (foundQuadFace) break;
            }

        }
#       ifdef FULLDEBUG
        if (requiresTriFaces)
        {
            if (!foundQuadFace || !foundTriFace)
            {
                FatalErrorIn
                (
                    "dynamicTopoFvMesh::constructPrismHull(...)"
                )
                << nl << " Failed to find a suitable quad/tri face. "
                << " Possibly not a prismatic mesh. " << nl
                << abort(FatalError);
            }
        }
        else
        {
            if (!foundQuadFace)
            {
                FatalErrorIn
                (
                    "dynamicTopoFvMesh::constructPrismHull(...)"
                )
                << nl << " Failed to find a suitable quad face. "
                << " Possibly not a prismatic mesh. " << nl
                << abort(FatalError);
            }
        }
#       endif

        // Decide which cell to check next
        if (owner_[faceToExclude] == cellIndex)
        {
            cellIndex = neighbour_[faceToExclude];
        }
        else
        {
            cellIndex = owner_[faceToExclude];
        }

        if (cellIndex == -1)
        {
            isBoundary=true;
            // Add the patchID to the list
            edgePatches.append(whichPatch(faceToExclude));
            break;
        }
        else
        {
            if (cellIndex != c0)
            {
                hullCells.append(cellIndex);
            }
        }

    } while ( faceToExclude != startFaceIndex );

    if (c1 == -1)
    {
        isBoundary = true;
        // Add the patchID to the list
        edgePatches.append(whichPatch(startFaceIndex));
    }
    else
    {
        // Check if the previous search hit a boundary.
        // If yes, start another search in the reverse direction.
        if (isBoundary)
        {
            // Start a search from cell[1] and add to the list as we go along
            faceToExclude=startFaceIndex, cellIndex=c1, hullCells.append(c1);
            do
            {
                cell& cellToCheck = cells_[cellIndex];
                foundQuadFace = false; foundTriFace = false;
                forAll(cellToCheck, faceI)
                {
                    if (cellToCheck[faceI] != faceToExclude)
                    {
                        face& faceToCheck = faces_[cellToCheck[faceI]];

                        if (faceToCheck.nEdges() == 4 && !foundQuadFace)
                        {
                            edgeList indexEdges = faceToCheck.edges();
                            forAll(indexEdges, edgeI)
                            {
                                if (indexEdges[edgeI] == edgeToCheck)
                                {
                                    // Bingo... We have a match. Add to the list
                                    hullFaces.append(cellToCheck[faceI]);
                                    faceToExclude = cellToCheck[faceI];
                                    foundQuadFace=true; break;
                                }
                            }
                        }

                        if
                        (
                            requiresTriFaces
                         && faceToCheck.nEdges() == 3
                         && !foundTriFace
                        )
                        {
                            hullTriFaces.append(cellToCheck[faceI]);
                            foundTriFace=true;
                        }
                    }

                    // Found the faces we were looking for, break-out
                    if (requiresTriFaces)
                    {
                        if (foundQuadFace && foundTriFace) break;
                    }
                    else
                    {
                        if (foundQuadFace) break;
                    }

                }
#               ifdef FULLDEBUG
                if (requiresTriFaces)
                {
                    if (!foundQuadFace || !foundTriFace)
                    {
                        FatalErrorIn
                        (
                            "dynamicTopoFvMesh::constructPrismHull(...)"
                        )
                        << nl << " Failed to find a suitable quad/tri face. "
                        << " Possibly not a prismatic mesh. " << nl
                        << abort(FatalError);
                    }
                }
                else
                {
                    if (!foundQuadFace)
                    {
                        FatalErrorIn
                        (
                            "dynamicTopoFvMesh::constructPrismHull(...)"
                        )
                        << nl << " Failed to find a suitable quad face. "
                        << " Possibly not a prismatic mesh. " << nl
                        << abort(FatalError);
                    }
                }
#               endif

                // Decide which cell to check next
                if (owner_[faceToExclude] == cellIndex)
                {
                    cellIndex = neighbour_[faceToExclude];
                }
                else
                {
                    cellIndex = owner_[faceToExclude];
                }

                if (cellIndex == -1)
                {
                    break;
                }
                else
                {
                    if (cellIndex != c0)
                    {
                        hullCells.append(cellIndex);
                    }
                }

            } while ( faceToExclude != startFaceIndex );

            // Add the starting Face index to the list as well
            hullFaces.append(startFaceIndex);
        }
    }

    return isBoundary;
}

// Utility method to build a counter-clockwise ring of vertices
// around the edge a-b (when viewed from vertex 'a')
inline bool dynamicTopoFvMesh::constructVertexRing
(
    const label eIndex,
    DynamicList<label>& hullCells,
    DynamicList<label>& hullFaces,
    DynamicList<label>& hullVertices,
    scalar& minQuality,
    labelHashSet& pLocks,
    labelHashSet& eLocks,
    labelHashSet& fLocks,
    labelHashSet& cLocks,
    bool requiresQuality = true
)
{
    bool found;
    minQuality = GREAT;
    label otherPoint = -1, nextPoint = -1, cellIndex = -1;
    label faceToExclude = -1, numPoints = 0, numFaces = 0;
    scalar cQuality = 0.0;

    // Try to lock this edge and its two points.
    if (threader_->multiThreaded())
    {
        if (edgeMutex_[eIndex].tryReadLock())
        {
            // Failed to acquire this edge.
            return true;
        }
        else
        {
            eLocks.insert(eIndex);

            edge& lockEdge = edges_[eIndex];

            // Check first point
            if (pointMutex_[lockEdge[0]].tryReadLock())
            {
                // Can't lock this point. Get out.
                unlockMutexLists(pLocks, eLocks, fLocks, cLocks);
                return true;
            }
            else
            {
                pLocks.insert(lockEdge[0]);
            }

            // Check second point
            if (pointMutex_[lockEdge[1]].tryReadLock())
            {
                // Can't lock this point. Get out.
                unlockMutexLists(pLocks, eLocks, fLocks, cLocks);
                return true;
            }
            else
            {
                pLocks.insert(lockEdge[1]);
            }
        }
    }

    // Obtain a reference to this edge
    edge& edgeToCheck = edges_[eIndex];

    // Decide which face to start with...
    labelList& eFaces = edgeFaces_[eIndex];
    label startFaceIndex = -1;

    // Check for an interior/boundary edge
    if (whichEdgePatch(eIndex) < 0)
    {
        // No orientation check, start with first entry
        startFaceIndex = eFaces[0];

        // Determine the orientation of the start-face
        findIsolatedPoint
        (
            faces_[startFaceIndex],
            edgeToCheck,
            otherPoint,
            nextPoint
        );
    }
    else
    {
        // Need to find a properly oriented start-face
        forAll(eFaces, faceI)
        {
            if (neighbour_[eFaces[faceI]] == -1)
            {
                findIsolatedPoint
                (
                    faces_[eFaces[faceI]],
                    edgeToCheck,
                    otherPoint,
                    nextPoint
                );

                if (nextPoint == edgeToCheck[0])
                {
                    startFaceIndex = eFaces[faceI];
                    break;
                }
            }
        }
    }

    // Try to lock the start-face and other-point
    if (threader_->multiThreaded())
    {
        if (faceMutex_[startFaceIndex].tryReadLock())
        {
            // Failed to acquire this face.
            unlockMutexLists(pLocks, eLocks, fLocks, cLocks);
            return true;
        }
        else
        {
            fLocks.insert(startFaceIndex);
        }
    }

    // Figure out the next cell to check
    if (nextPoint == edgeToCheck[0])
    {
        cellIndex = owner_[startFaceIndex];
    }
    else
    if (nextPoint == edgeToCheck[1])
    {
        cellIndex = neighbour_[startFaceIndex];
    }
    else
    {
        // Something's terribly wrong
        FatalErrorIn
        (
            "dynamicTopoFvMesh::constructVertexRing(...)"
        )
        << nl << " Failed to determine a vertex ring. "
        << " Possibly not a tetrahedral mesh. " << nl
        << abort(FatalError);
    }

    // Start a search and add to the list as we go along
    faceToExclude = startFaceIndex;

    do
    {
        if (threader_->multiThreaded())
        {
            // Try to acquire this cell
            if (cellMutex_[cellIndex].tryReadLock())
            {
                // Failed to acquire this cell.
                unlockMutexLists(pLocks, eLocks, fLocks, cLocks);
                return true;
            }
            else
            {
                cLocks.insert(cellIndex);
            }
        }

        cell& cellToCheck = cells_[cellIndex];

        // Add this point to the hull
        hullVertices(numPoints++) = otherPoint;

        // Add this face to the hull
        hullFaces(numFaces++) = faceToExclude;

        // Add this cell to the hull
        hullCells.append(cellIndex);

        found = false;

        // Loop through faces of this cell
        forAll(cellToCheck, faceI)
        {
            // Loop through edgeFaces and get the next face
            forAll(eFaces, faceII)
            {
                if
                (
                    eFaces[faceII] != faceToExclude
                 && eFaces[faceII] == cellToCheck[faceI]
                )
                {
                    face& faceToCheck = faces_[cellToCheck[faceI]];

                    if (requiresQuality)
                    {
                        // Check face-orientation and compute cell-quality
                        if (owner_[cellToCheck[faceI]] == cellIndex)
                        {
                            cQuality = (*tetMetric_)
                            (
                                meshPoints_[faceToCheck[2]],
                                meshPoints_[faceToCheck[1]],
                                meshPoints_[faceToCheck[0]],
                                meshPoints_[otherPoint]
                            );
                        }
                        else
                        {
                            cQuality = (*tetMetric_)
                            (
                                meshPoints_[faceToCheck[0]],
                                meshPoints_[faceToCheck[1]],
                                meshPoints_[faceToCheck[2]],
                                meshPoints_[otherPoint]
                            );
                        }

                        // Check if the quality is worse
                        minQuality = cQuality < minQuality
                                   ? cQuality : minQuality;

#                       ifdef FULLDEBUG
                        if (minQuality < 0)
                        {
                            Info << nl << nl
                                 << "*** Detected negative cell-quality! ***"
                                 << nl << "Quality: " << minQuality << nl
                                 << "Cell: " << cellIndex << ": " << cellToCheck
                                 << nl << "when using face: " << faceToCheck
                                 << " and point: " << otherPoint
                                 << nl << "Owner: "
                                 << owner_[cellToCheck[faceI]]
                                 << nl << "Neighbour: "
                                 << neighbour_[cellToCheck[faceI]]
                                 << endl;

                            Info << "Points: " << endl;
                            forAll(faceToCheck, pointI)
                            {
                                Info << faceToCheck[pointI]
                                     << ": " << meshPoints_[faceToCheck[pointI]]
                                     << endl;
                            }
                            Info << otherPoint
                                 << ": " << meshPoints_[otherPoint]
                                 << endl;

                            Info << "Faces: " << endl;
                            forAll(cellToCheck, faceI)
                            {
                                Info << cellToCheck[faceI] << ": "
                                     << faces_[cellToCheck[faceI]]
                                     << endl;
                            }

                            // Something's terribly wrong
                            FatalErrorIn
                            (
                                "dynamicTopoFvMesh::constructVertexRing(...)"
                            )
                            << " Not a valid tetrahedral mesh. "
                            << abort(FatalError);
                        }
#                       endif
                    }

                    // Find the isolated point
                    findIsolatedPoint
                    (
                        faceToCheck,
                        edgeToCheck,
                        otherPoint,
                        nextPoint
                    );

                    // Update faceToExclude
                    faceToExclude = cellToCheck[faceI];

                    // Try to lock the face
                    if (threader_->multiThreaded())
                    {
                        if (faceMutex_[faceToExclude].tryReadLock())
                        {
                            // Failed to acquire this face.
                            unlockMutexLists(pLocks, eLocks, fLocks, cLocks);
                            return true;
                        }
                        else
                        {
                            fLocks.insert(faceToExclude);
                        }
                    }

                    found = true; break;
                }
            }

            if (found) break;
        }

        if (!found)
        {
            // Something's terribly wrong
            FatalErrorIn
            (
                "dynamicTopoFvMesh::constructVertexRing(...)"
            )
            << " Failed to determine a vertex ring. " << nl
            << " edgeFaces connectivity is inconsistent. " << nl
            << "Edge: " << eIndex << ":: " << edgeToCheck << nl
            << "edgeFaces: " << eFaces
            << abort(FatalError);
        }

        // Decide which cell to check next
        if (nextPoint == edgeToCheck[0])
        {
            cellIndex = owner_[faceToExclude];
        }
        else if (nextPoint == edgeToCheck[1])
        {
            cellIndex = neighbour_[faceToExclude];
        }
        else
        {
            // Something's terribly wrong
            FatalErrorIn
            (
                "dynamicTopoFvMesh::constructVertexRing(...)"
            )
            << " Failed to determine a vertex ring. " << nl
            << " Possibly not a valid tetrahedral mesh. "
            << abort(FatalError);
        }

        // Check if this is a boundary
        if (cellIndex == -1)
        {
            // Add this point to the hull
            hullVertices(numPoints++) = otherPoint;

            // Add this face to the hull
            hullFaces(numFaces++) = faceToExclude;

            // Add this cell to the hull
            hullCells.append(cellIndex);

            break;
        }

    } while ( faceToExclude != startFaceIndex );

#   ifdef FULLDEBUG
    if (debug)
    {
        // Print out the ring
        // Info << endl;
        // Info << "Edge: " << edgeToCheck << endl;
        // Info << "Points: " << hullVertices << endl;
        // Info << "Faces: " << hullFaces << endl;
        // Info << "Cells: " << hullCells << endl;
    }
#   endif

    // Shrink all dynamic lists
    hullVertices.shrink();
    hullFaces.shrink();
    hullCells.shrink();

    // Return a successful lock
    return false;
}

// Utility to obtain write permissions for the cell hull.
// Returns zero on success
inline bool dynamicTopoFvMesh::obtainWritePriority
(
    const DynamicList<label>& hullCells,
    labelHashSet& pLocks,
    labelHashSet& eLocks,
    labelHashSet& fLocks,
    labelHashSet& cLocks
)
{
    bool haveContention = false;

    forAll(hullCells, cellI)
    {
        label cIndex = hullCells[cellI];

        if (cIndex != -1)
        {
            // Try to acquire this cell
            if (cellMutex_[cIndex].tryWriteLock())
            {
                // Failed to obtain write permission.
                // Check if contentionMutex is available.
                if (!haveContention)
                {
                    if (contentionMutex_.tryLock())
                    {
                        // Another thread has contention. 
                        unlockMutexLists
                        (
                            pLocks,
                            eLocks,
                            fLocks,
                            cLocks
                        );

                        return true;
                    }
                    else
                    {
                        // We have contention.
                        haveContention = true;
                    }
                }

                // Contention is available. Wait and lock.
                cellMutex_[cIndex].writeLock();
            }

            cLocks.insert(cIndex);

            // Try to acquire all faces of this cell
            cell& cellToCheck = cells_[cIndex];

            forAll(cellToCheck, faceI)
            {
                label fIndex = cellToCheck[faceI];

                if (!fLocks.found(fIndex))
                {
                    if (faceMutex_[fIndex].tryWriteLock())
                    {
                        // Failed to obtain write permission.
                        // Check if contentionMutex is available.
                        if (!haveContention)
                        {
                            if (contentionMutex_.tryLock())
                            {
                                // Another thread has contention.
                                unlockMutexLists
                                (
                                    pLocks,
                                    eLocks,
                                    fLocks,
                                    cLocks
                                );

                                return true;
                            }
                            else
                            {
                                // We have contention.
                                haveContention = true;
                            }
                        }

                        // Contention is available. Wait and lock.
                        faceMutex_[fIndex].writeLock();
                    }

                    fLocks.insert(fIndex);

                    // Try to acquire all edges of this face
                    labelList& fEdges = faceEdges_[fIndex];
                    forAll (fEdges, edgeI)
                    {
                        label eIndex = fEdges[edgeI];

                        if (!eLocks.found(eIndex))
                        {
                            if (edgeMutex_[eIndex].tryWriteLock())
                            {
                                // Failed to obtain write permission.
                                // Check if contentionMutex is available.
                                if (!haveContention)
                                {
                                    if (contentionMutex_.tryLock())
                                    {
                                        // Another thread has contention.
                                        unlockMutexLists
                                        (
                                            pLocks,
                                            eLocks,
                                            fLocks,
                                            cLocks
                                        );

                                        return true;
                                    }
                                    else
                                    {
                                        // We have contention.
                                        haveContention = true;
                                    }
                                }

                                // Contention is available. Wait and lock.
                                edgeMutex_[eIndex].writeLock();
                            }

                            eLocks.insert(eIndex);
                        }
                    }

                    // Try to acquire all points of this face
                    face& faceToCheck = faces_[fIndex];
                    forAll (faceToCheck, pointI)
                    {
                        label pIndex = faceToCheck[pointI];

                        if (!pLocks.found(pIndex))
                        {
                            if (pointMutex_[pIndex].tryWriteLock())
                            {
                                // Failed to obtain write permission.
                                // Check if contentionMutex is available.
                                if (!haveContention)
                                {
                                    if (contentionMutex_.tryLock())
                                    {
                                        // Another thread has contention.
                                        unlockMutexLists
                                        (
                                            pLocks,
                                            eLocks,
                                            fLocks,
                                            cLocks
                                        );

                                        return true;
                                    }
                                    else
                                    {
                                        // We have contention.
                                        haveContention = true;
                                    }
                                }

                                // Contention is available. Wait and lock.
                                pointMutex_[pIndex].writeLock();
                            }

                            pLocks.insert(pIndex);
                        }
                    }
                }
            }
        }
    }

    if (haveContention)
    {
        contentionMutex_.unlock();
    }

    // Return successful lock
    return false;
}

// Check whether the given edge is on a bounding curve
inline bool dynamicTopoFvMesh::checkBoundingCurve(label eIndex)
{
    // Internal edges don't count
    label edgePatch = -1;
    if ((edgePatch = whichEdgePatch(eIndex)) < 0)
    {
        return false;
    }
    else
    {
        // Check whether this edge shouldn't be swapped
        if (noSwapPatchIDs_.found(edgePatch))
        {
            return true;
        }
    }

    // Check if two boundary faces lie on different face-patches
    label fPatch, firstPatch = -1, secondPatch = -1;
    labelList& edgeFaces = edgeFaces_[eIndex];
    forAll(edgeFaces, faceI)
    {
        if
        (
            (fPatch = whichPatch(edgeFaces[faceI])) > -1
        )
        {
            if (firstPatch == -1)
            {
                firstPatch = fPatch;
            }
            else
            {
                secondPatch = fPatch;
                break;
            }
        }
    }

    if (firstPatch != secondPatch)
    {
        return true;
    }

    // Not on a bounding curve
    return false;
}

// Allocate dynamic programming tables
inline void dynamicTopoFvMesh::initTables
(
    const label mMax,
    scalarListList& Q,
    labelListList& K,
    labelListList& triangulations
)
{
    Q.setSize((mMax-2),scalarList(mMax,-1.0));
    K.setSize((mMax-2),labelList(mMax,-1));
    triangulations.setSize(3,labelList((mMax-2),-1));
}

// Utility method to fill the dynamic programming tables
// Returns the number of triangulations
inline label dynamicTopoFvMesh::fillTables
(
    const label eIndex,
    const DynamicList<label>& hullVertices,
    const scalar minQuality,
    scalarListList& Q,
    labelListList& K
)
{
    label m = hullVertices.size();
    edge& edgeToCheck = edges_[eIndex];

    for (label i = m-3; i >= 0; i--)
    {
        for (label j = i+2; j < m; j++)
        {
            for (label k = i+1; k < j; k++)
            {
                scalar qA = (*tetMetric_)
                (
                    meshPoints_[hullVertices[i]],
                    meshPoints_[hullVertices[k]],
                    meshPoints_[hullVertices[j]],
                    meshPoints_[edgeToCheck[0]]
                );

                /*
                if (qA < minQuality)
                {
                    Q[i][j] = qA;
                    break;
                }
                */

                scalar qB = (*tetMetric_)
                (
                    meshPoints_[hullVertices[j]],
                    meshPoints_[hullVertices[k]],
                    meshPoints_[hullVertices[i]],
                    meshPoints_[edgeToCheck[1]]
                );

                scalar q = Foam::min(qA,qB);

                if (k < j-1)
                {
                    q = Foam::min(q,Q[k][j]);
                }

                if (k > i+1)
                {
                    q = Foam::min(q,Q[i][k]);
                }

                if ((k == i+1) || (q > Q[i][j]))
                {
                    Q[i][j] = q;
                    K[i][j] = k;
                }
            }
        }
    }

    return m;
}

// Remove the edge according to the swap sequence
void dynamicTopoFvMesh::removeEdgeFlips
(
    const label m,
    const label eIndex,
    const labelListList& K,
    const DynamicList<label>& hullCells,
    const DynamicList<label>& hullFaces,
    const DynamicList<label>& hullVertices,
    labelListList& triangulations
)
{
    label numTriangulations = 0, isolatedVertex = -1;
    edge& edgeToCheck = edges_[eIndex];

    // Extract the appropriate triangulations
    extractTriangulation(0, (m-1), K, numTriangulations, triangulations);

    // Determine the 3-2 swap triangulation
    label t32 = identify32Swap(m, edgeToCheck, hullVertices, triangulations);

    // Perform a series of 2-3 swaps
    label numSwaps = 0;
    while (numSwaps < (m-3))
    {
        for (label i = 0; i < (m-2); i++)
        {
            if ( (i != t32) && (triangulations[0][i] != -1) )
            {
                // Check if triangulation is on the boundary
                if
                (
                    boundaryTriangulation
                    (
                        i,
                        isolatedVertex,
                        triangulations
                    )
                )
                {
                    // Perform 2-3 swap
                    swap23
                    (
                        isolatedVertex,
                        eIndex,
                        edgeToCheck,
                        i,
                        triangulations,
                        hullCells,
                        hullFaces,
                        hullVertices
                    );

                    // Done with this face, so reset it
                    triangulations[0][i] = -1;
                    triangulations[1][i] = -1;
                    triangulations[2][i] = -1;

                    numSwaps++;
                }
            }
        }

        if (numSwaps == 0)
        {
            Info << "Triangulations: " << endl;
            forAll(triangulations, row)
            {
                Info << triangulations[row] << endl;
            }

            // Should have performed at least one swap
            FatalErrorIn("dynamicTopoFvMesh::removeEdgeFlips()") << nl
                << "Did not perform any 2-3 swaps" << nl
                << abort(FatalError);
        }
    }

    // Perform the final 3-2 swap
    swap32
    (
        eIndex,
        edgeToCheck,
        t32,
        triangulations,
        hullCells,
        hullFaces,
        hullVertices
    );

    // Done with this face, so reset it
    triangulations[0][t32] = -1;
    triangulations[1][t32] = -1;
    triangulations[2][t32] = -1;

    // Lock the edge mutex
    eMutex_.lock();

    // Finally remove the edge
    removeEdge(eIndex);

    // Unlock the edge mutex
    eMutex_.unlock();
}

// Extract triangulations from the programming table
void dynamicTopoFvMesh::extractTriangulation
(
    const label i,
    const label j,
    const labelListList& K,
    label& numTriangulations,
    labelListList& triangulations
)
{
    if ( j >= (i+2) )
    {
        label k = K[i][j];

        // Fill in the triangulation list
        triangulations[0][numTriangulations] = i;
        triangulations[1][numTriangulations] = k;
        triangulations[2][numTriangulations] = j;

        // Increment triangulation count
        numTriangulations++;

        // Recursively call the function for the two sub-triangulations
        extractTriangulation(i,k,K,numTriangulations,triangulations);
        extractTriangulation(k,j,K,numTriangulations,triangulations);
    }
}

// Identify the 3-2 swap from the triangulation sequence
label dynamicTopoFvMesh::identify32Swap
(
    const label m,
    const edge& edgeToCheck,
    const DynamicList<label>& hullVertices,
    const labelListList& triangulations
)
{
    for (label i = 0; i < (m-2); i++)
    {
        pointHit curHit = triPointRef
        (
            meshPoints_[hullVertices[triangulations[0][i]]],
            meshPoints_[hullVertices[triangulations[1][i]]],
            meshPoints_[hullVertices[triangulations[2][i]]]
        ).ray
        (
            meshPoints_[edgeToCheck[0]],
            (meshPoints_[edgeToCheck[1]]-meshPoints_[edgeToCheck[0]])
        );

        if (curHit.hit())
        {
            return i;
        }
        else
        if (curHit.eligibleMiss())
        {
            return i;
        }
    }

    // This bit should never occur.
    Info << "Hull Vertices: " << endl;

    forAll(hullVertices, vertexI)
    {
        Info << hullVertices[vertexI] << ": "
             << meshPoints_[hullVertices[vertexI]]
             << endl;
    }

    FatalErrorIn("dynamicTopoFvMesh::identify32Swap()") << nl
        << "Could not determine 3-2 swap triangulation." << nl
        << "Edge: " << edgeToCheck << nl
        << "Edge Points: "
        << meshPoints_[edgeToCheck[0]] << ","
        << meshPoints_[edgeToCheck[1]] << nl
        << abort(FatalError);

    return -1;
}

// Routine to check whether the triangulation at the
// index lies on the boundary of the vertex ring.
bool dynamicTopoFvMesh::boundaryTriangulation
(
    const label index,
    label& isolatedVertex,
    labelListList& triangulations
)
{
    label first = 0, second = 0, third = 0;

    // Count for occurrences
    forAll(triangulations, row)
    {
        forAll(triangulations[row], col)
        {
            if (triangulations[row][col] == triangulations[0][index])
            {
                first++;
            }

            if (triangulations[row][col] == triangulations[1][index])
            {
                second++;
            }

            if (triangulations[row][col] == triangulations[2][index])
            {
                third++;
            }
        }
    }

    if (first == 1)
    {
        isolatedVertex = triangulations[0][index];
        return true;
    }

    if (second == 1)
    {
        isolatedVertex = triangulations[1][index];
        return true;
    }

    if (third == 1)
    {
        isolatedVertex = triangulations[2][index];
        return true;
    }

    // This isn't a boundary triangulation
    return false;
}

// Routine to perform 2-3 swaps
void dynamicTopoFvMesh::swap23
(
    const label isolatedVertex,
    const label edgeToCheckIndex,
    const edge& edgeToCheck,
    const label triangulationIndex,
    const labelListList& triangulations,
    const DynamicList<label>& hullCells,
    const DynamicList<label>& hullFaces,
    const DynamicList<label>& hullVertices
)
{
    // A 2-3 swap performs the following operations:
    //      [1] Remove face: [ edge[0] edge[1] isolatedVertex ]
    //      [2] Remove two cells on either side of removed face
    //      [3] Add one edge
    //      [4] Add three new faces
    //      [5] Add three new cells
    //      Update faceEdges and edgeFaces information

#   ifdef FULLDEBUG
    if (debug)
    {
        // Print out arguments
        Info << endl;
        Info << "== Swapping 2-3 ==" << endl;
        Info << "Edge: " << edgeToCheckIndex << ": " << edgeToCheck << endl;
        Info << "Ring: " << hullVertices << endl;
        Info << "Faces: " << hullFaces << endl;
        Info << "Cells: " << hullCells << endl;
        Info << "Triangulation: "
             << triangulations[0][triangulationIndex] << " "
             << triangulations[1][triangulationIndex] << " "
             << triangulations[2][triangulationIndex] << " "
             << endl;
        Info << "Isolated vertex: " << isolatedVertex << endl;
    }
#   endif

    label faceForRemoval = hullFaces[isolatedVertex];
    label vertexForRemoval = hullVertices[isolatedVertex];

    // Determine the two cells to be removed
    FixedList<label,2> cellsForRemoval;
    cellsForRemoval[0] = owner_[faceForRemoval];
    cellsForRemoval[1] = neighbour_[faceForRemoval];

    // Lock the cell mutex
    cMutex_.lock();

    // Add three new cells to the end of the cell list
    FixedList<label,3> newCellIndex(-1);
    newCellIndex[0] = cells_.append(cell(4));
    newCellIndex[1] = cells_.append(cell(4));
    newCellIndex[2] = cells_.append(cell(4));

    cell &newTetCell0 = cells_[newCellIndex[0]];
    cell &newTetCell1 = cells_[newCellIndex[1]];
    cell &newTetCell2 = cells_[newCellIndex[2]];

    // Update length-scale info
    if (edgeModification_)
    {
        scalar avgScale =
        (
             lengthScale_[cellsForRemoval[0]]
           + lengthScale_[cellsForRemoval[1]]
        )/2.0;

        for (label i = 0; i < 3; i++)
        {
            lengthScale_.append(avgScale);
        }
    }

    // Add three unlocked cell mutexes
    if (threader_->multiThreaded())
    {
        for (label i = 0; i < 3; i++)
        {
            cellMutex_.append();
        }
    }

    // Unlock the cell mutex
    cMutex_.unlock();

    // Obtain point-ordering for the other vertices
    // otherVertices[0] is the point before isolatedVertex
    // otherVertices[1] is the point after isolatedVertex
    FixedList<label,2> otherVertices;

    if (triangulations[0][triangulationIndex] == isolatedVertex)
    {
        otherVertices[0] = hullVertices[triangulations[2][triangulationIndex]];
        otherVertices[1] = hullVertices[triangulations[1][triangulationIndex]];
    }
    else
    if (triangulations[1][triangulationIndex] == isolatedVertex)
    {
        otherVertices[0] = hullVertices[triangulations[0][triangulationIndex]];
        otherVertices[1] = hullVertices[triangulations[2][triangulationIndex]];
    }
    else
    {
        otherVertices[0] = hullVertices[triangulations[1][triangulationIndex]];
        otherVertices[1] = hullVertices[triangulations[0][triangulationIndex]];
    }

    // Lock the face mutex
    fMutex_.lock();

    // Insert three new internal faces
    FixedList<label,3> newFaceIndex;
    face tmpTriFace(3);

    // First face: The actual triangulation
    tmpTriFace[0] = otherVertices[0];
    tmpTriFace[1] = vertexForRemoval;
    tmpTriFace[2] = otherVertices[1];

    newFaceIndex[0] = insertFace
                      (
                          -1,
                          tmpTriFace,
                          newCellIndex[0],
                          newCellIndex[1]
                      );

    // Second face: Triangle involving edgeToCheck[0]
    tmpTriFace[0] = otherVertices[0];
    tmpTriFace[1] = edgeToCheck[0];
    tmpTriFace[2] = otherVertices[1];

    newFaceIndex[1] = insertFace
                      (
                          -1,
                          tmpTriFace,
                          newCellIndex[1],
                          newCellIndex[2]
                      );

    // Third face: Triangle involving edgeToCheck[1]
    tmpTriFace[0] = otherVertices[1];
    tmpTriFace[1] = edgeToCheck[1];
    tmpTriFace[2] = otherVertices[0];

    newFaceIndex[2] = insertFace
                      (
                          -1,
                          tmpTriFace,
                          newCellIndex[0],
                          newCellIndex[2]
                      );

    // Append three dummy faceEdges entries.
    for (label i = 0; i < 3; i++)
    {
        faceEdges_.append();
    }

    // Unlock the face mutex
    fMutex_.unlock();

    // Add an entry to edgeFaces
    labelList newEdgeFaces(3);
    newEdgeFaces[0] = newFaceIndex[0];
    newEdgeFaces[1] = newFaceIndex[1];
    newEdgeFaces[2] = newFaceIndex[2];

    // Lock the edge mutex
    eMutex_.lock();

    // Add a new internal edge to the mesh
    label newEdgeIndex = insertEdge
                         (
                             -1,
                             edge
                             (
                                 otherVertices[0],
                                 otherVertices[1]
                             ),
                             newEdgeFaces
                         );

    // Unlock the edge mutex
    eMutex_.unlock();

    // Define the six edges to check while building faceEdges:
    FixedList<edge,6> check;

    check[0][0] = vertexForRemoval; check[0][1] = otherVertices[0];
    check[1][0] = vertexForRemoval; check[1][1] = otherVertices[1];
    check[2][0] = edgeToCheck[0];   check[2][1] = otherVertices[0];
    check[3][0] = edgeToCheck[1];   check[3][1] = otherVertices[0];
    check[4][0] = edgeToCheck[0];   check[4][1] = otherVertices[1];
    check[5][0] = edgeToCheck[1];   check[5][1] = otherVertices[1];

    // Add three new entries to faceEdges
    label nE0 = 0, nE1 = 0, nE2 = 0;
    labelListList newFaceEdges(3,labelList(3));

    newFaceEdges[0][nE0++] = newEdgeIndex;
    newFaceEdges[1][nE1++] = newEdgeIndex;
    newFaceEdges[2][nE2++] = newEdgeIndex;

    // Fill-in information for the three new cells,
    // and correct info on existing neighbouring cells
    label nF0 = 0, nF1 = 0, nF2 = 0;
    FixedList<bool,2> foundEdge;

    // Add the newly created faces to cells
    newTetCell0[nF0++] = newFaceIndex[0];
    newTetCell0[nF0++] = newFaceIndex[2];
    newTetCell1[nF1++] = newFaceIndex[0];
    newTetCell1[nF1++] = newFaceIndex[1];
    newTetCell2[nF2++] = newFaceIndex[1];
    newTetCell2[nF2++] = newFaceIndex[2];

    forAll(cellsForRemoval, cellI)
    {
        label cellIndex = cellsForRemoval[cellI];
        cell& cellToCheck = cells_[cellIndex];

        forAll(cellToCheck, faceI)
        {
            label faceIndex = cellToCheck[faceI];
            face& faceToCheck = faces_[faceIndex];

            foundEdge[0] = false; foundEdge[1] = false;

            // Find the face that contains only
            // edgeToCheck[0] or edgeToCheck[1]
            forAll(faceToCheck, pointI)
            {
                if (faceToCheck[pointI] == edgeToCheck[0])
                {
                    foundEdge[0] = true;
                }

                if (faceToCheck[pointI] == edgeToCheck[1])
                {
                    foundEdge[1] = true;
                }
            }

            // Face is connected to edgeToCheck[0]
            if (foundEdge[0] && !foundEdge[1])
            {
                // Check if a face-flip is necessary
                if(owner_[faceIndex] == cellIndex)
                {
                    if (neighbour_[faceIndex] == -1)
                    {
                        // Change the owner
                        owner_[faceIndex] = newCellIndex[1];
                    }
                    else
                    {
                        // Flip this face
                        faces_[faceIndex] = faceToCheck.reverseFace();
                        owner_[faceIndex] = neighbour_[faceIndex];
                        neighbour_[faceIndex] = newCellIndex[1];
                    }
                }
                else
                {
                    // Flip is unnecessary. Just update neighbour
                    neighbour_[faceIndex] = newCellIndex[1];
                }

                // Add this face to the cell
                newTetCell1[nF1++] = faceIndex;

                // Update faceEdges and edgeFaces, and add them to the stack
                const labelList& fEdges = faceEdges_[faceIndex];
                forAll(fEdges, edgeI)
                {
                    if (edges_[fEdges[edgeI]] == check[0])
                    {
                        newFaceEdges[0][nE0++] = fEdges[edgeI];
                        sizeUpList(newFaceIndex[0], edgeFaces_[fEdges[edgeI]]);

                        edgeStack(getThreadID(pthread_self())).push(fEdges[edgeI]);
                    }

                    if (edges_[fEdges[edgeI]] == check[1])
                    {
                        newFaceEdges[0][nE0++] = fEdges[edgeI];
                        sizeUpList(newFaceIndex[0], edgeFaces_[fEdges[edgeI]]);

                        edgeStack(getThreadID(pthread_self())).push(fEdges[edgeI]);
                    }

                    if (edges_[fEdges[edgeI]] == check[2])
                    {
                        newFaceEdges[1][nE1++] = fEdges[edgeI];
                        sizeUpList(newFaceIndex[1], edgeFaces_[fEdges[edgeI]]);

                        edgeStack(getThreadID(pthread_self())).push(fEdges[edgeI]);
                    }

                    if (edges_[fEdges[edgeI]] == check[4])
                    {
                        newFaceEdges[1][nE1++] = fEdges[edgeI];
                        sizeUpList(newFaceIndex[1], edgeFaces_[fEdges[edgeI]]);

                        edgeStack(getThreadID(pthread_self())).push(fEdges[edgeI]);
                    }
                }
            }

            // Face is connected to edgeToCheck[1]
            if (!foundEdge[0] && foundEdge[1])
            {
                // Check if a face-flip is necessary
                if(owner_[faceIndex] == cellIndex)
                {
                    if (neighbour_[faceIndex] == -1)
                    {
                        // Change the owner
                        owner_[faceIndex] = newCellIndex[0];
                    }
                    else
                    {
                        // Flip this face
                        faces_[faceIndex] = faceToCheck.reverseFace();
                        owner_[faceIndex] = neighbour_[faceIndex];
                        neighbour_[faceIndex] = newCellIndex[0];
                    }
                }
                else
                {
                    // Flip is unnecessary. Just update neighbour
                    neighbour_[faceIndex] = newCellIndex[0];
                }

                // Add this face to the cell
                newTetCell0[nF0++] = faceIndex;

                // Update faceEdges and edgeFaces
                const labelList& fEdges = faceEdges_[faceIndex];
                forAll(fEdges, edgeI)
                {
                    if (edges_[fEdges[edgeI]] == check[3])
                    {
                        newFaceEdges[2][nE2++] = fEdges[edgeI];
                        sizeUpList(newFaceIndex[2], edgeFaces_[fEdges[edgeI]]);

                        edgeStack(getThreadID(pthread_self())).push(fEdges[edgeI]);
                    }

                    if (edges_[fEdges[edgeI]] == check[5])
                    {
                        newFaceEdges[2][nE2++] = fEdges[edgeI];
                        sizeUpList(newFaceIndex[2], edgeFaces_[fEdges[edgeI]]);

                        edgeStack(getThreadID(pthread_self())).push(fEdges[edgeI]);
                    }
                }
            }

            // Face is connected to both edgeToCheck [0] and [1]
            if
            (
                (foundEdge[0] && foundEdge[1])
             && (faceIndex != faceForRemoval)
            )
            {
                // Check if a face-flip is necessary
                if(owner_[faceIndex] == cellIndex)
                {
                    if (neighbour_[faceIndex] == -1)
                    {
                        // Change the owner
                        owner_[faceIndex] = newCellIndex[2];
                    }
                    else
                    {
                        // Flip this face
                        faces_[faceIndex] = faceToCheck.reverseFace();
                        owner_[faceIndex] = neighbour_[faceIndex];
                        neighbour_[faceIndex] = newCellIndex[2];
                    }
                }
                else
                {
                    // Flip is unnecessary. Just update neighbour
                    neighbour_[faceIndex] = newCellIndex[2];
                }

                // Add this face to the cell
                newTetCell2[nF2++] = faceIndex;
            }
        }
    }

    // Now update faceEdges for the three new faces
    forAll(newFaceEdges, faceI)
    {
        faceEdges_[newFaceIndex[faceI]] = newFaceEdges[faceI];
    }

    // Generate mapping information for the three new cells
    // Prepare a list of master-objects to map from
    const labelListList& cc = cellCells();
    labelHashSet masterObjects;
    FixedList<label,2> parents(-1);

    forAll (cellsForRemoval, indexI)
    {
        // Determine an appropriate parent cell
        if (cellsForRemoval[indexI] < nOldCells_)
        {
            parents[indexI] = cellsForRemoval[indexI];
        }
        else
        {
            parents[indexI] = cellParents_[cellsForRemoval[indexI]];
        }

        // Find the cell's neighbours in the old mesh
        masterObjects.insert(parents[indexI]);
        forAll(cc[parents[indexI]],cellI)
        {
            if (!masterObjects.found(cc[parents[indexI]][cellI]))
            {
                masterObjects.insert(cc[parents[indexI]][cellI]);
            }
        }
    }

    // Lock the cell mutex
    cMutex_.lock();

    forAll(newCellIndex, cellI)
    {
        // Insert the parent cell [from first by default]
        cellParents_.insert(newCellIndex[cellI], parents[0]);

        // Insert mapping info into the HashTable
        cellsFromCells_.insert
        (
            newCellIndex[cellI],
            objectMap
            (
                newCellIndex[cellI],
                masterObjects.toc()
            )
        );
    }

    // Unlock the cell mutex
    cMutex_.unlock();

    // Lock the face mutex
    fMutex_.lock();

    // Remove the face
    removeFace(faceForRemoval);

    // Update edgeFaces for edges of the removed face
    labelList& fEdges = faceEdges_[faceForRemoval];
    forAll(fEdges, edgeI)
    {
        sizeDownList(faceForRemoval, edgeFaces_[fEdges[edgeI]]);

        edgeStack(getThreadID(pthread_self())).push(fEdges[edgeI]);
    }

    // Now remove the faceEdges entry
    faceEdges_.remove(faceForRemoval);

    // Unlock the face mutex
    fMutex_.unlock();

    // Lock the cell mutex
    cMutex_.lock();

    // Update the number of cells, and the reverse cell map
    nCells_++;

    forAll(cellsForRemoval, cellI)
    {
        label cIndex = cellsForRemoval[cellI];

#       ifdef FULLDEBUG
        if (debug)
        {
            Info << "Removing cell: "
                 << cIndex << ": "
                 << cells_[cIndex]
                 << endl;
        }
#       endif

        cells_.remove(cIndex);

        // Remove the cell mutex
        if (threader_->multiThreaded())
        {
            // Unlock it first
            cellMutex_[cIndex].unlock();
            cellMutex_.remove(cIndex);
        }

        if (edgeModification_)
        {
            lengthScale_.remove(cIndex);
        }

        if (cIndex < nOldCells_)
        {
            reverseCellMap_[cIndex] = -1;
        }

        // Check if the cell was added in the current morph, and delete
        if (cellsFromCells_.found(cIndex))
        {
            cellsFromCells_.erase(cIndex);
        }
    }

    // Unlock the cell mutex
    cMutex_.unlock();

#   ifdef FULLDEBUG
    if (debug)
    {
        Info << "Added edge: " << endl;
        Info << newEdgeIndex << ":: "
             << edges_[newEdgeIndex]
             << " edgeFaces: "
             << edgeFaces_[newEdgeIndex]
             << endl;
        Info << "Added faces: " << endl;
        forAll(newFaceIndex, faceI)
        {
            Info << newFaceIndex[faceI] << ":: "
                 << faces_[newFaceIndex[faceI]]
                 << " faceEdges: "
                 << faceEdges_[newFaceIndex[faceI]]
                 << endl;
        }
        Info << "Added cells: " << endl;
        forAll(newCellIndex, cellI)
        {
            Info << newCellIndex[cellI] << ":: "
                 << cells_[newCellIndex[cellI]]
                 << endl;
        }
    }
#   endif
}

// Routine to perform 3-2 or 2-2 swaps
void dynamicTopoFvMesh::swap32
(
    const label edgeToCheckIndex,
    const edge& edgeToCheck,
    const label triangulationIndex,
    const labelListList& triangulations,
    const DynamicList<label>& hullCells,
    const DynamicList<label>& hullFaces,
    const DynamicList<label>& hullVertices
)
{
    // A 2-2 / 3-2 swap performs the following operations:
    //      [1] Remove three faces surrounding edgeToCheck
    //      [2] Remove two (2-2 swap) or three(3-2 swap)
    //          cells surrounding edgeToCheck
    //      [3] Add one internal face
    //      [4] Add two new cells
    //      [5] If edgeToCheck is on a boundary,
    //          add two boundary faces and a boundary edge (2-2 swap)
    //      edgeToCheck is removed later by swap3DEdges
    //      Update faceEdges and edgeFaces information

    // Determine the patch this edge belongs to
    label edgePatch = whichEdgePatch(edgeToCheckIndex);

#   ifdef FULLDEBUG
    if (debug)
    {
        // Print out arguments
        Info << endl;
        if (edgePatch < 0)
        {
            Info << "== Swapping 3-2 ==" << endl;
        }
        else
        {
            Info << "== Swapping 2-2 ==" << endl;
        }
        Info << "Edge: " << edgeToCheckIndex << ": " << edgeToCheck << endl;
        Info << "Ring: " << hullVertices << endl;
        Info << "Faces: " << hullFaces << endl;
        Info << "Cells: " << hullCells << endl;
        Info << "Triangulation: "
             << triangulations[0][triangulationIndex] << " "
             << triangulations[1][triangulationIndex] << " "
             << triangulations[2][triangulationIndex] << " "
             << endl;
    }
#   endif

    // Determine the three faces to be removed
    FixedList<label,3> facesForRemoval;
    labelHashSet cellsForRemoval(3);

    forAll(facesForRemoval, faceI)
    {
        facesForRemoval[faceI]
            = hullFaces[triangulations[faceI][triangulationIndex]];

        label own = owner_[facesForRemoval[faceI]];
        label nei = neighbour_[facesForRemoval[faceI]];

        // Check and add cells as well
        if (!cellsForRemoval.found(own))
        {
            cellsForRemoval.insert(own);
        }

        if (!cellsForRemoval.found(nei) && nei != -1)
        {
            cellsForRemoval.insert(nei);
        }
    }

    labelList cellRemovalList = cellsForRemoval.toc();

    // Lock the cell mutex
    cMutex_.lock();

    // Add two new cells to the end of the cell list
    FixedList<label,2> newCellIndex(-1);
    newCellIndex[0] = cells_.append(cell(4));
    newCellIndex[1] = cells_.append(cell(4));

    cell &newTetCell0 = cells_[newCellIndex[0]];
    cell &newTetCell1 = cells_[newCellIndex[1]];

    // Update length-scale info
    if (edgeModification_)
    {
        scalar avgScale = 0.0;

        forAll(cellRemovalList, cellI)
        {
            avgScale += lengthScale_[cellRemovalList[cellI]];
        }

        avgScale /= cellRemovalList.size();

        for (label i = 0; i < 2; i++)
        {
            lengthScale_.append(avgScale);
        }
    }

    // Add three unlocked cell mutexes
    if (threader_->multiThreaded())
    {
        for (label i = 0; i < 2; i++)
        {
            cellMutex_.append();
        }
    }

    // Unlock the cell mutex
    cMutex_.unlock();

    // Lock the face mutex
    fMutex_.lock();

    // Insert a new internal face
    face newTriFace(3);

    newTriFace[0] = hullVertices[triangulations[0][triangulationIndex]];
    newTriFace[1] = hullVertices[triangulations[1][triangulationIndex]];
    newTriFace[2] = hullVertices[triangulations[2][triangulationIndex]];

    label newFaceIndex = insertFace
                         (
                             -1,
                             newTriFace,
                             newCellIndex[0],
                             newCellIndex[1]
                         );

    // Add faceEdges for the new face as well.
    faceEdges_.append(labelList(3));

    // Unlock the face mutex
    fMutex_.unlock();

    // Define the three edges to check while building faceEdges:
    FixedList<edge,3> check;

    check[0][0] = newTriFace[0]; check[0][1] = newTriFace[1];
    check[1][0] = newTriFace[1]; check[1][1] = newTriFace[2];
    check[2][0] = newTriFace[2]; check[2][1] = newTriFace[0];

    // New faceEdge entry for the interior face
    label nE = 0;
    labelList& newFaceEdges = faceEdges_[newFaceIndex];

    // For 2-2 swaps, two faces are introduced
    FixedList<label,2> nBE(0);
    labelListList bdyFaceEdges(2, labelList(3, -1));

    // Fill-in information for the two new cells,
    // and correct info on existing neighbouring cells
    label nF0 = 0, nF1 = 0;
    FixedList<bool,2> foundEdge;

    // For a 2-2 swap on a boundary edge,
    // add two boundary faces and an edge
    FixedList<label,2> newBdyFaceIndex(-1);
    label newEdgeIndex = -1;

    if (edgePatch > -1)
    {
        // Temporary local variables
        label nBE0 = 0, nBE1 = 0;
        label otherPoint = -1, nextPoint = -1, facePatch = -1;
        FixedList<label,2> bdyEdges0(-1), bdyEdges1(-1);
        FixedList<face,2> newBdyTriFace(face(3));
        edge newEdge(-1, -1);

        // Get a cue for face orientation from existing faces
        forAll(facesForRemoval, faceI)
        {
            if (neighbour_[facesForRemoval[faceI]] == -1)
            {
                facePatch = whichPatch(facesForRemoval[faceI]);
                face& faceToCheck = faces_[facesForRemoval[faceI]];

                findIsolatedPoint
                (
                    faceToCheck,
                    edgeToCheck,
                    otherPoint,
                    nextPoint
                );

                if (nextPoint == edgeToCheck[0])
                {
                    newEdge[1] = otherPoint;
                    newBdyTriFace[0][0] = otherPoint;
                    newBdyTriFace[0][1] = edgeToCheck[0];
                    newBdyTriFace[1][2] = otherPoint;
                }
                else
                {
                    newEdge[0] = otherPoint;
                    newBdyTriFace[1][0] = otherPoint;
                    newBdyTriFace[1][1] = edgeToCheck[1];
                    newBdyTriFace[0][2] = otherPoint;
                }

                // Also update faceEdges for the new boundary faces
                labelList& fEdges = faceEdges_[facesForRemoval[faceI]];

                forAll(fEdges, edgeI)
                {
                    if
                    (
                        edges_[fEdges[edgeI]]
                     == edge(edgeToCheck[0], otherPoint)
                    )
                    {
                        bdyFaceEdges[0][nBE[0]++] = fEdges[edgeI];
                        bdyEdges0[nBE0++] = fEdges[edgeI];

                        edgeStack(getThreadID(pthread_self())).push(fEdges[edgeI]);
                    }

                    if
                    (
                        edges_[fEdges[edgeI]]
                     == edge(edgeToCheck[1], otherPoint)
                    )
                    {
                        bdyFaceEdges[1][nBE[1]++] = fEdges[edgeI];
                        bdyEdges1[nBE1++] = fEdges[edgeI];

                        edgeStack(getThreadID(pthread_self())).push(fEdges[edgeI]);
                    }
                }
            }
        }

        // Lock the face mutex
        fMutex_.lock();

        // Insert the two new faces
        newBdyFaceIndex[0] = insertFace
                             (
                                 facePatch,
                                 newBdyTriFace[0],
                                 newCellIndex[1],
                                 -1
                             );

        newBdyFaceIndex[1] = insertFace
                             (
                                 facePatch,
                                 newBdyTriFace[1],
                                 newCellIndex[0],
                                 -1
                             );

        // Update the new cells
        newTetCell0[nF0++] = newBdyFaceIndex[1];
        newTetCell1[nF1++] = newBdyFaceIndex[0];

        // Insert the new edge
        labelList newBdyEdgeFaces(3, -1);
        newBdyEdgeFaces[0] = newBdyFaceIndex[0];
        newBdyEdgeFaces[1] = newFaceIndex;
        newBdyEdgeFaces[2] = newBdyFaceIndex[1];

        // Lock the edge mutex
        eMutex_.lock();

        // Insert the edge
        newEdgeIndex = insertEdge(edgePatch, newEdge, newBdyEdgeFaces);

        // Unlock the edge mutex
        eMutex_.unlock();

        // Update faceEdges with the new edge
        newFaceEdges[nE++] = newEdgeIndex;
        bdyFaceEdges[0][nBE[0]++] = newEdgeIndex;
        bdyFaceEdges[1][nBE[1]++] = newEdgeIndex;

        // Update edgeFaces with the two new faces
        forAll(bdyEdges0, edgeI)
        {
            sizeUpList(newBdyFaceIndex[0], edgeFaces_[bdyEdges0[edgeI]]);
            sizeUpList(newBdyFaceIndex[1], edgeFaces_[bdyEdges1[edgeI]]);
        }

        // Add faceEdges for the two new boundary faces
        faceEdges_.append(bdyFaceEdges[0]);
        faceEdges_.append(bdyFaceEdges[1]);

        // Unlock the face mutex
        fMutex_.unlock();
    }

    newTetCell0[nF0++] = newFaceIndex;
    newTetCell1[nF1++] = newFaceIndex;

    forAll(cellRemovalList, cellI)
    {
        label cellIndex = cellRemovalList[cellI];
        cell& cellToCheck = cells_[cellIndex];

        forAll(cellToCheck, faceI)
        {
            label faceIndex = cellToCheck[faceI];
            face& faceToCheck = faces_[faceIndex];

            foundEdge[0] = false; foundEdge[1] = false;

            // Find the face that contains only
            // edgeToCheck[0] or edgeToCheck[1]
            forAll(faceToCheck, pointI)
            {
                if (faceToCheck[pointI] == edgeToCheck[0])
                {
                    foundEdge[0] = true;
                }

                if (faceToCheck[pointI] == edgeToCheck[1])
                {
                    foundEdge[1] = true;
                }
            }

            // Face is connected to edgeToCheck[0]
            if (foundEdge[0] && !foundEdge[1])
            {
                // Check if a face-flip is necessary
                if(owner_[faceIndex] == cellIndex)
                {
                    if (neighbour_[faceIndex] == -1)
                    {
                        // Change the owner
                        owner_[faceIndex] = newCellIndex[1];
                    }
                    else
                    {
                        // Flip this face
                        faces_[faceIndex] = faceToCheck.reverseFace();
                        owner_[faceIndex] = neighbour_[faceIndex];
                        neighbour_[faceIndex] = newCellIndex[1];
                    }
                }
                else
                {
                    // Flip is unnecessary. Just update neighbour
                    neighbour_[faceIndex] = newCellIndex[1];
                }

                // Add this face to the cell
                newTetCell1[nF1++] = faceIndex;

                // Update faceEdges and edgeFaces
                const labelList& fEdges = faceEdges_[faceIndex];
                forAll(fEdges, edgeI)
                {
                    if
                    (
                        (edges_[fEdges[edgeI]] == check[0])
                     || (edges_[fEdges[edgeI]] == check[1])
                     || (edges_[fEdges[edgeI]] == check[2])
                    )
                    {
                        newFaceEdges[nE++] = fEdges[edgeI];
                        sizeUpList(newFaceIndex, edgeFaces_[fEdges[edgeI]]);

                        edgeStack(getThreadID(pthread_self())).push(fEdges[edgeI]);

                        break;
                    }
                }
            }

            // Face is connected to edgeToCheck[1]
            if (!foundEdge[0] && foundEdge[1])
            {
                // Check if a face-flip is necessary
                if(owner_[faceIndex] == cellIndex)
                {
                    if (neighbour_[faceIndex] == -1)
                    {
                        // Change the owner
                        owner_[faceIndex] = newCellIndex[0];
                    }
                    else
                    {
                        // Flip this face
                        faces_[faceIndex] = faceToCheck.reverseFace();
                        owner_[faceIndex] = neighbour_[faceIndex];
                        neighbour_[faceIndex] = newCellIndex[0];
                    }
                }
                else
                {
                    // Flip is unnecessary. Just update neighbour
                    neighbour_[faceIndex] = newCellIndex[0];
                }

                // Add this face to the cell
                newTetCell0[nF0++] = faceIndex;
            }
        }
    }

    // Generate mapping information for the two new cells
    const labelListList& cc = cellCells();
    labelHashSet masterObjects;
    FixedList<label,3> parents(-1);

    forAll (cellRemovalList, indexI)
    {
        // Determine an appropriate parent cell
        if (cellRemovalList[indexI] < nOldCells_)
        {
            parents[indexI] = cellRemovalList[indexI];
        }
        else
        {
            parents[indexI] = cellParents_[cellRemovalList[indexI]];
        }

        // Find the cell's neighbours in the old mesh
        masterObjects.insert(parents[indexI]);
        forAll(cc[parents[indexI]],cellI)
        {
            if (!masterObjects.found(cc[parents[indexI]][cellI]))
            {
                masterObjects.insert(cc[parents[indexI]][cellI]);
            }
        }
    }

    // Lock the cell mutex
    cMutex_.lock();

    forAll(newCellIndex, cellI)
    {
        // Insert the parent cell [from first by default]
        cellParents_.insert(newCellIndex[cellI], parents[0]);

        // Insert mapping info into the HashTable
        cellsFromCells_.insert
        (
            newCellIndex[cellI],
            objectMap
            (
                newCellIndex[cellI],
                masterObjects.toc()
            )
        );
    }

    // Unlock the cell mutex
    cMutex_.unlock();

    // Lock the face mutex
    fMutex_.lock();

    // Remove the faces and update associated edges
    forAll(facesForRemoval, faceI)
    {
        removeFace(facesForRemoval[faceI]);

        // Update edgeFaces
        labelList& fEdges = faceEdges_[facesForRemoval[faceI]];
        forAll(fEdges, edgeI)
        {
            label edgeIndex = fEdges[edgeI];

            if (edgeIndex != edgeToCheckIndex)
            {
                sizeDownList(facesForRemoval[faceI],edgeFaces_[edgeIndex]);

                edgeStack(getThreadID(pthread_self())).push(edgeIndex);
            }
        }

        // Now remove the faceEdges entry
        faceEdges_.remove(facesForRemoval[faceI]);
    }

    // Unlock the face mutex
    fMutex_.unlock();

    // Lock the cell mutex
    cMutex_.lock();

    if (edgePatch < 0)
    {
        // Update the number of cells only for 3-2 swaps
        nCells_--;
    }

    forAll(cellRemovalList, cellI)
    {
        label cIndex = cellRemovalList[cellI];

#       ifdef FULLDEBUG
        if (debug)
        {
            Info << "Removing cell: "
                 << cIndex << ": "
                 << cells_[cIndex]
                 << endl;
        }
#       endif

        cells_.remove(cIndex);

        // Remove the cell mutex
        if (threader_->multiThreaded())
        {
            // Unlock it first
            cellMutex_[cIndex].unlock();
            cellMutex_.remove(cIndex);
        }

        if (edgeModification_)
        {
            lengthScale_.remove(cIndex);
        }

        if (cIndex < nOldCells_)
        {
            reverseCellMap_[cIndex] = -1;
        }

        // Check if the cell was added in the current morph, and delete
        if (cellsFromCells_.found(cIndex))
        {
            cellsFromCells_.erase(cIndex);
        }
    }

    // Unlock the cell mutex
    cMutex_.unlock();

#   ifdef FULLDEBUG
    if (debug)
    {
        if (edgePatch > -1)
        {
            Info << "Added edge: " << endl;
            Info << newEdgeIndex << ":: "
                 << edges_[newEdgeIndex]
                 << " edgeFaces: "
                 << edgeFaces_[newEdgeIndex]
                 << endl;
        }
        Info << "Added face(s): " << endl;
        Info << newFaceIndex << ":: "
             << faces_[newFaceIndex];
        Info << " faceEdges: "
             << faceEdges_[newFaceIndex]
             << endl;
        if (edgePatch > -1)
        {
            forAll(newBdyFaceIndex, faceI)
            {
                Info << newBdyFaceIndex[faceI] << ":: "
                     << faces_[newBdyFaceIndex[faceI]]
                     << " faceEdges: "
                     << faceEdges_[newBdyFaceIndex[faceI]]
                     << endl;
            }
        }
        Info << "Added cells: " << endl;
        forAll(newCellIndex, cellI)
        {
            Info << newCellIndex[cellI] << ":: "
                 << cells_[newCellIndex[cellI]]
                 << endl;
        }
    }
#   endif
}

// Reorder points after a topology change
void dynamicTopoFvMesh::reOrderPoints
(
    pointField& points
)
{
    // *** Point renumbering *** //
    // If points were deleted during topology change, the numerical order ceases
    // to be continuous. Loop through all points and renumber sequentially.
    // Possible scope for bandwidth-reduction on the motion-solver.

    // Allocate for the mapping information
    pointMap_.setSize(nPoints_, -1);

    label pointRenum = 0;

    addedPointRenumbering_.clear();

    HashList<point>::iterator ptIter = meshPoints_.begin();
    HashList<labelList>::iterator peIter;

    if (!twoDMesh_)
    {
        peIter = pointEdges_.begin();
    }

    while(ptIter != meshPoints_.end())
    {
        // Obtain the index for this point
        label pIndex = ptIter.index();

        // Update the point info
        points[pointRenum] = ptIter();

        // Renumber the point index
        meshPoints_.reNumber(pointRenum, ptIter);

        if (!twoDMesh_)
        {
            pointEdges_.reNumber(pointRenum, peIter);
            peIter++;
        }

        // Added points are always numbered after nOldPoints_
        // (by virtue of the HashList append method)
        if (pIndex < nOldPoints_)
        {
            pointMap_[pointRenum]    = pIndex;
            reversePointMap_[pIndex] = pointRenum;
        }
        else
        {
            addedPointRenumbering_.insert(pIndex,pointRenum);
        }

        // Update the counter
        pointRenum++;

        // Update the iterators
        ptIter++;
    }
}

// Reorder edges after a topology change
void dynamicTopoFvMesh::reOrderEdges()
{
    // *** Edge renumbering *** //
    // If edges were deleted during topology change, the numerical order ceases
    // to be continuous. Edges are added to respective internal/boundary patches

    // Allocate for mapping information
    edgeMap_.setSize(nEdges_, -1);

    label edgeInOrder = 0, allEdges = edges_.lastIndex() + 1;
    edgeList oldEdges(allEdges);
    labelListList oldEdgeFaces(allEdges);

    addedEdgeRenumbering_.clear();
    Map<label> addedEdgeReverseRenumbering;

    // Make a copy of the old edge-based HashLists, and clear them
    HashList<edge>::iterator eIter = edges_.begin();
    HashList<labelList>::iterator efIter = edgeFaces_.begin();

    while (eIter != edges_.end())
    {
        oldEdges[eIter.index()] = eIter();
        oldEdgeFaces[efIter.index()] = efIter();
        eIter++; efIter++;
    }

    edges_.clear(); edgeFaces_.clear();

    // Keep track of inserted boundary edge indices
    labelList boundaryPatchIndices(edgePatchStarts_);

    // Loop through all edges and add internal ones first
    forAll(oldEdges, edgeI)
    {
        // Ensure that we're adding valid edges
        if (oldEdgeFaces[edgeI].size() > 0)
        {
            // Determine which patch this edge belongs to
            label patch = whichEdgePatch(edgeI);

            // Obtain references
            edge& thisEdge = oldEdges[edgeI];
            labelList& thisEF = oldEdgeFaces[edgeI];

            // Renumber edges
            if (thisEdge[0] < nOldPoints_)
            {
                thisEdge[0] = reversePointMap_[thisEdge[0]];
            }
            else
            {
                thisEdge[0] = addedPointRenumbering_[thisEdge[0]];
            }

            if (thisEdge[1] < nOldPoints_)
            {
                thisEdge[1] = reversePointMap_[thisEdge[1]];
            }
            else
            {
                thisEdge[1] = addedPointRenumbering_[thisEdge[1]];
            }

            // Renumber edgeFaces
            forAll(thisEF,faceI)
            {
                if (thisEF[faceI] < nOldFaces_)
                {
                    thisEF[faceI] = reverseFaceMap_[thisEF[faceI]];
                }
                else
                {
                    thisEF[faceI] = addedFaceRenumbering_[thisEF[faceI]];
                }
            }

            // Update maps for boundary edges. Edge insertion for
            // boundaries will be done after internal edges.
            if (patch >= 0)
            {
                label bEdgeIndex = boundaryPatchIndices[patch]++;

                // Update the maps
                if (edgeI < nOldEdges_)
                {
                    edgeMap_[bEdgeIndex] = edgeI;
                    reverseEdgeMap_[edgeI] = bEdgeIndex;
                }
                else
                {
                    addedEdgeRenumbering_.insert(edgeI,bEdgeIndex);
                    addedEdgeReverseRenumbering.insert(bEdgeIndex,edgeI);
                }
            }
            else
            {
                // Renumber internal edges and add normally.
                if (edgeI < nOldEdges_)
                {
                    edgeMap_[edgeInOrder] = edgeI;
                    reverseEdgeMap_[edgeI] = edgeInOrder;
                }
                else
                {
                    addedEdgeRenumbering_.insert(edgeI,edgeInOrder);
                }

                // Insert entities into HashLists...
                edges_.append(thisEdge);
                edgeFaces_.append(thisEF);

                edgeInOrder++;
            }
        }
    }

    // All internal edges have been inserted. Now insert boundary edges.
    label oldIndex;
    for(label i=nInternalEdges_; i<nEdges_; i++)
    {
        if (edgeMap_[i] == -1)
        {
            // This boundary edge was added during the topology change
            oldIndex = addedEdgeReverseRenumbering[i];
        }
        else
        {
            oldIndex = edgeMap_[i];
        }

        // Insert entities into HashLists...
        edges_.append(oldEdges[oldIndex]);
        edgeFaces_.append(oldEdgeFaces[oldIndex]);

        edgeInOrder++;
    }

    // Renumber faceEdges
    forAllIter(HashList<labelList>::iterator, faceEdges_, feIter)
    {
        labelList& faceEdges = feIter();
        forAll(faceEdges,edgeI)
        {
            if (faceEdges[edgeI] < nOldEdges_)
            {
                faceEdges[edgeI] = reverseEdgeMap_[faceEdges[edgeI]];
            }
            else
            {
                faceEdges[edgeI] = addedEdgeRenumbering_[faceEdges[edgeI]];
            }
        }
    }

    // Renumber pointEdges
    forAllIter(HashList<labelList>::iterator, pointEdges_, peIter)
    {
        labelList& pointEdges = peIter();
        forAll(pointEdges,edgeI)
        {
            if (pointEdges[edgeI] < nOldEdges_)
            {
                pointEdges[edgeI] = reverseEdgeMap_[pointEdges[edgeI]];
            }
            else
            {
                pointEdges[edgeI] = addedEdgeRenumbering_[pointEdges[edgeI]];
            }
        }
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        checkEdgeConnectivity();
    }
#   endif

    // Set edge connectivity in IOLists
    setEdgeConnectivity();
}

// Reorder faces in upper-triangular order after a topology change
void dynamicTopoFvMesh::reOrderFaces
(
    faceList& faces,
    labelList& owner,
    labelList& neighbour
)
{
    // *** Face renumbering *** //
    // Faces have to be renumbered if any were added/deleted/modified
    // Boundary faces are added to respective patches.
    // Internal faces, however, have to be added in upper-triangular ordering;
    // i.e., in the increasing order of neighbours

    // Allocate for mapping information
    faceMap_.setSize(nFaces_, -1);

    label faceInOrder = 0, allFaces = faces_.lastIndex() + 1;
    faceList oldFaces(allFaces);
    labelList oldOwner(allFaces), oldNeighbour(allFaces), visited(allFaces,0);
    edgeList oldEdgeToWatch(0);
    labelListList oldFaceEdges(0);

    addedFaceRenumbering_.clear();
    Map<label> addedFaceReverseRenumbering;

    // Make a copy of the old face-based HashLists, and clear them
    HashList<face>::iterator fIter = faces_.begin();
    HashList<label>::iterator oIter = owner_.begin();
    HashList<label>::iterator nIter = neighbour_.begin();

    while(fIter != faces_.end())
    {
        oldFaces[fIter.index()] = fIter();
        oldOwner[oIter.index()] = oIter();
        oldNeighbour[nIter.index()] = nIter();
        fIter++; oIter++; nIter++;
    }

    faces_.clear(); owner_.clear(); neighbour_.clear();

    if (edgeModification_ && twoDMesh_)
    {
        oldEdgeToWatch.setSize(allFaces);

        forAllIter(HashList<edge>::iterator, edgeToWatch_, eIter)
        {
            oldEdgeToWatch[eIter.index()] = eIter();
        }

        edgeToWatch_.clear();
    }

    if (!twoDMesh_)
    {
        oldFaceEdges.setSize(allFaces);

        forAllIter(HashList<labelList>::iterator, faceEdges_, feIter)
        {
            oldFaceEdges[feIter.index()] = feIter();
        }

        faceEdges_.clear();
    }

    // Mark the internal faces with -2 so that they are inserted first
    forAllIter(HashList<cell>::iterator, cells_, cIter)
    {
        const cell& curFaces = cIter();
        forAll(curFaces, faceI)
        {
            visited[curFaces[faceI]]--;
        }
    }

    // Upper-triangular ordering of faces:

    // Keep track of inserted boundary face indices
    labelList boundaryPatchIndices(patchStarts_);

    // Insertion cannot be done in one go as the faces need to be
    // added into the list in the increasing order of neighbour
    // cells.  Therefore, all neighbours will be detected first
    // and then added in the correct order.
    forAllIter(HashList<cell>::iterator, cells_, cIter)
    {
        // Record the neighbour cell
        label cellI = cIter.index();
        const cell& curFaces = cIter();
        labelList neiCells(curFaces.size(), -1);

        label nNeighbours = 0;

        forAll (curFaces, faceI)
        {
            if (visited[curFaces[faceI]] == -2)
            {
                // Face is internal and gets reordered
                label own =   oldOwner[curFaces[faceI]] < nOldCells_
                            ? reverseCellMap_[oldOwner[curFaces[faceI]]]
                            : addedCellRenumbering_[oldOwner[curFaces[faceI]]];
                label nei =   oldNeighbour[curFaces[faceI]] < nOldCells_
                            ? reverseCellMap_[oldNeighbour[curFaces[faceI]]]
                            : addedCellRenumbering_[oldNeighbour[curFaces[faceI]]];

                label smallerIndex = own < nei ? own : nei;
                label largerIndex  = own > nei ? own : nei;

                if (cellI == smallerIndex)
                {
                    neiCells[faceI] = largerIndex;
                    nNeighbours++;
                }
            }

            // Boundary faces are inserted normally. Update maps for now.
            // Face insertion for boundaries will be done after internal faces.
            if (visited[curFaces[faceI]] == -1)
            {
                label patchID = whichPatch(curFaces[faceI]);
                label bFaceIndex = boundaryPatchIndices[patchID]++;

                // Renumber the point-labels for this boundary-face
                face& faceRenumber = oldFaces[curFaces[faceI]];
                forAll(faceRenumber,pointI)
                {
                    if (faceRenumber[pointI] < nOldPoints_)
                    {
                        faceRenumber[pointI] =
                            reversePointMap_[faceRenumber[pointI]];
                    }
                    else
                    {
                        faceRenumber[pointI] =
                            addedPointRenumbering_[faceRenumber[pointI]];
                    }
                }

                // Renumber the edges in edgeToWatch
                if (edgeModification_ && twoDMesh_)
                {
                    edge& edgeRenumber = oldEdgeToWatch[curFaces[faceI]];

                    if (edgeRenumber[0] < nOldPoints_)
                    {
                        edgeRenumber[0] = reversePointMap_[edgeRenumber[0]];
                    }
                    else
                    {
                        edgeRenumber[0] = addedPointRenumbering_[edgeRenumber[0]];
                    }

                    if (edgeRenumber[1] < nOldPoints_)
                    {
                        edgeRenumber[1] = reversePointMap_[edgeRenumber[1]];
                    }
                    else
                    {
                        edgeRenumber[1] = addedPointRenumbering_[edgeRenumber[1]];
                    }
                }

                // Update the maps
                if (curFaces[faceI] < nOldFaces_)
                {
                    faceMap_[bFaceIndex] = curFaces[faceI];
                    reverseFaceMap_[curFaces[faceI]] = bFaceIndex;
                }
                else
                {
                    addedFaceRenumbering_.insert(curFaces[faceI],bFaceIndex);
                    addedFaceReverseRenumbering.insert(bFaceIndex,curFaces[faceI]);
                }

                // Mark this face as visited
                visited[curFaces[faceI]] = 0;
            }
        }

        // Add internal faces in the increasing order of neighbours
        for (label neiSearch = 0; neiSearch < nNeighbours; neiSearch++)
        {
            // Find the lowest neighbour which is still valid
            label nextNei = -1;
            label minNei = nCells_;

            forAll (neiCells, ncI)
            {
                if (neiCells[ncI] > -1 && neiCells[ncI] < minNei)
                {
                    nextNei = ncI;
                    minNei = neiCells[ncI];
                }
            }

            if (nextNei > -1)
            {
                // Face is internal and gets reordered
                if (curFaces[nextNei] < nOldFaces_)
                {
                    faceMap_[faceInOrder] = curFaces[nextNei];
                    reverseFaceMap_[curFaces[nextNei]] = faceInOrder;
                }
                else
                {
                    addedFaceRenumbering_.insert(curFaces[nextNei],faceInOrder);
                }

                // Renumber the point labels in this face
                face& faceRenumber = oldFaces[curFaces[nextNei]];
                forAll(faceRenumber, pointI)
                {
                    if (faceRenumber[pointI] < nOldPoints_)
                    {
                        faceRenumber[pointI] =
                            reversePointMap_[faceRenumber[pointI]];
                    }
                    else
                    {
                        faceRenumber[pointI] =
                            addedPointRenumbering_[faceRenumber[pointI]];
                    }
                }

                // Renumber owner/neighbour
                label oldOwn = oldOwner[curFaces[nextNei]];
                label oldNei = oldNeighbour[curFaces[nextNei]];

                label ownerRenumber =
                    oldOwner[curFaces[nextNei]] < nOldCells_
                  ? reverseCellMap_[oldOwn] : addedCellRenumbering_[oldOwn];

                label neighbourRenumber =
                    oldNeighbour[curFaces[nextNei]] < nOldCells_
                  ? reverseCellMap_[oldNei] : addedCellRenumbering_[oldNei];

                // Cell-reordering may have cause flipped faces.. Correct them.
                if (neighbourRenumber < ownerRenumber)
                {
                    faceRenumber = faceRenumber.reverseFace();
                }

                // Insert entities into HashLists...
                faces_.append(faceRenumber);
                owner_.append(cellI);
                neighbour_.append(minNei);

                if (edgeModification_ && twoDMesh_)
                {
                    edge& edgeRenumber = oldEdgeToWatch[curFaces[nextNei]];

                    if (edgeRenumber[0] < nOldPoints_)
                    {
                        edgeRenumber[0] = reversePointMap_[edgeRenumber[0]];
                    }
                    else
                    {
                        edgeRenumber[0] = addedPointRenumbering_[edgeRenumber[0]];
                    }

                    if (edgeRenumber[1] < nOldPoints_)
                    {
                        edgeRenumber[1] = reversePointMap_[edgeRenumber[1]];
                    }
                    else
                    {
                        edgeRenumber[1] = addedPointRenumbering_[edgeRenumber[1]];
                    }

                    edgeToWatch_.append(edgeRenumber);
                }

                if (!twoDMesh_)
                {
                    faceEdges_.append(oldFaceEdges[curFaces[nextNei]]);
                }

                // Insert entities into mesh-reset lists
                faces[faceInOrder] = faceRenumber;
                owner[faceInOrder] = cellI;
                neighbour[faceInOrder] = minNei;

                // Stop the neighbour from being used again
                neiCells[nextNei] = -1;

                // Mark this face as visited
                visited[curFaces[nextNei]] = 0;

                faceInOrder++;
            }
            else
            {
                FatalErrorIn("dynamicTopoFvMesh::reOrderFaces()") << nl
                    << "Error in internal face insertion" << nl
                    << abort(FatalError);
            }
        }
    }

    // All internal faces have been inserted. Now insert boundary faces.
    label oldIndex;
    for(label i=nInternalFaces_; i<nFaces_; i++)
    {
        if (faceMap_[i] == -1)
        {
            // This boundary face was added during the topology change
            oldIndex = addedFaceReverseRenumbering[i];
        }
        else
        {
            oldIndex = faceMap_[i];
        }

        // Renumber owner/neighbour
        label ownerRenumber =   oldOwner[oldIndex] < nOldCells_
                              ? reverseCellMap_[oldOwner[oldIndex]]
                              : addedCellRenumbering_[oldOwner[oldIndex]];

        // Insert entities into HashLists...
        faces_.append(oldFaces[oldIndex]);
        owner_.append(ownerRenumber);
        neighbour_.append(-1);

        if (edgeModification_ && twoDMesh_)
        {
            edgeToWatch_.append(oldEdgeToWatch[oldIndex]);
        }

        if (!twoDMesh_)
        {
            faceEdges_.append(oldFaceEdges[oldIndex]);
        }

        // Insert entities into mesh-reset lists
        faces[faceInOrder] = oldFaces[oldIndex];
        owner[faceInOrder] = ownerRenumber;
        neighbour[faceInOrder] = -1;

        faceInOrder++;
    }

    // Renumber all cells
    forAllIter(HashList<cell>::iterator, cells_, cIter)
    {
        cell& cellFaces = cIter();
        forAll(cellFaces,faceI)
        {
            if (cellFaces[faceI] < nOldFaces_)
            {
                cellFaces[faceI] = reverseFaceMap_[cellFaces[faceI]];
            }
            else
            {
                cellFaces[faceI] = addedFaceRenumbering_[cellFaces[faceI]];
            }
        }
    }

    // Final check to ensure everything went okay
    if (debug)
    {
        if (sum(visited) != 0)
        {
            FatalErrorIn("dynamicTopoFvMesh::reOrderFaces()") << nl
                    << " Algorithm did not visit every face in the mesh."
                    << " Something's messed up." << nl
                    << abort(FatalError);
        }
    }
}

// Reorder & renumber cells with bandwidth reduction after a topology change
void dynamicTopoFvMesh::reOrderCells()
{
    // *** Cell renumbering *** //
    // If cells were deleted during topology change, the numerical order ceases
    // to be continuous. Also, cells are always added at the end of the list by
    // virtue of the HashList append method. Thus, cells would now have to be
    // reordered so that bandwidth is reduced and renumbered to be sequential.

    // Allocate for mapping information
    cellMap_.setSize(nCells_, -1);

    label currentCell, cellInOrder = 0, allCells = cells_.lastIndex() + 1;
    SLList<label> nextCell;
    labelList ncc(allCells, 0);
    labelList visited(allCells, 0);
    labelListList cellCellAddr(allCells);
    cellList oldCells(allCells);
    scalarField oldLengthScale(0);

    addedCellRenumbering_.clear();

    // Make a copy of the old cell-based HashLists, and clear them
    forAllIter(HashList<cell>::iterator, cells_, cIter)
    {
        oldCells[cIter.index()] = cIter();
    }
    cells_.clear();

    if (edgeModification_)
    {
        oldLengthScale.setSize(allCells);
        forAllIter(HashList<scalar>::iterator, lengthScale_, cIter)
        {
            oldLengthScale[cIter.index()] = cIter();
        }
        lengthScale_.clear();
    }

    // Build a cell-cell addressing list
    HashList<label>::iterator ownIter = owner_.begin();
    HashList<label>::iterator neiIter = neighbour_.begin();

    while(ownIter != owner_.end())
    {
        if (neiIter() != -1)
        {
            ncc[ownIter()]++;
            ncc[neiIter()]++;
        }
        ownIter++; neiIter++;
    }

    forAll (cellCellAddr, cellI)
    {
        cellCellAddr[cellI].setSize(ncc[cellI]);

        // Mark off deleted cells as "visited"
        if (ncc[cellI] == 0)
        {
            visited[cellI] = 1;
        }
    }

    ncc = 0;
    ownIter = owner_.begin(); neiIter = neighbour_.begin();
    while(ownIter != owner_.end())
    {
        if (neiIter() != -1)
        {
            cellCellAddr[ownIter()][ncc[ownIter()]++] = neiIter();
            cellCellAddr[neiIter()][ncc[neiIter()]++] = ownIter();
        }
        ownIter++; neiIter++;
    }

    // Let's get to the "business bit" of the band-compression
    forAll (visited, cellI)
    {
        // Find the first cell that has not been visited yet
        if (visited[cellI] == 0)
        {
            // Use this cell as a start
            currentCell = cellI;
            nextCell.append(currentCell);

            // Loop through the nextCell list. Add the first cell into the
            // cell order if it has not already been visited and ask for its
            // neighbours. If the neighbour in question has not been visited,
            // add it to the end of the nextCell list
            while (nextCell.size() > 0)
            {
                currentCell = nextCell.removeHead();

                if (visited[currentCell] == 0)
                {
                    // Mark as visited and update cell mapping info
                    visited[currentCell] = 1;
                    if (currentCell < nOldCells_)
                    {
                        cellMap_[cellInOrder] = currentCell;
                        reverseCellMap_[currentCell] = cellInOrder;
                    }
                    else
                    {
                        addedCellRenumbering_.insert(currentCell,cellInOrder);
                    }

                    // Insert entities into HashLists...
                    cells_.append(oldCells[currentCell]);

                    if (edgeModification_)
                    {
                        lengthScale_.append(oldLengthScale[currentCell]);
                    }

                    cellInOrder++;

                    // Find if the neighbours have been visited
                    const labelList& neighbours = cellCellAddr[currentCell];

                    forAll (neighbours, nI)
                    {
                        if (visited[neighbours[nI]] == 0)
                        {
                            // Not visited, add to the list
                            nextCell.append(neighbours[nI]);
                        }
                    }
                }
            }
        }
    }

    // Loop through the cellsFromCells list, and renumber the map indices
    // HashTable keys, however, are not altered.
    forAllIter(Map<objectMap>, cellsFromCells_, cellI)
    {
        objectMap& thisMap = cellI();
        if (thisMap.index() < nOldCells_)
        {
            thisMap.index() = reverseCellMap_[thisMap.index()];
        }
        else
        {
            thisMap.index() = addedCellRenumbering_[thisMap.index()];
        }
    }

    if(debug)
    {
        if (sum(visited) != allCells)
        {
            FatalErrorIn("dynamicTopoFvMesh::reOrderCells()") << nl
                    << " Algorithm did not visit every cell in the mesh."
                    << " Something's messed up." << nl
                    << abort(FatalError);
        }
    }
}

// Reorder the faces in upper-triangular order, and generate mapping information
void dynamicTopoFvMesh::reOrderMesh
(
    pointField& points,
    faceList& faces,
    labelList& owner,
    labelList& neighbour
)
{
    if (debug)
    {
        Info << endl;
        Info << "=================" << endl;
        Info << " Mesh reOrdering " << endl;
        Info << "=================" << endl;
        Info << "Mesh Info [n]:" << endl;
        Info << "Points: " << nOldPoints_ << endl;
        Info << "Edges: " << nOldEdges_ << endl;
        Info << "Faces: " << nOldFaces_ << endl;
        Info << "Cells: " << nOldCells_ << endl;
        Info << "Internal Edges: " << nOldInternalEdges_ << endl;
        Info << "Internal Faces: " << nOldInternalFaces_ << endl;
        Info << "Patch Starts: " << oldPatchStarts_ << endl;
        Info << "Patch Sizes: " << oldPatchSizes_ << endl;
        Info << "=================" << endl;
        Info << "Mesh Info [n+1]:" << endl;
        Info << "Points: " << nPoints_ << endl;
        Info << "Edges: " << nEdges_ << endl;
        Info << "Faces: " << nFaces_ << endl;
        Info << "Cells: " << nCells_ << endl;
        Info << "Internal Edges: " << nInternalEdges_ << endl;
        Info << "Internal Faces: " << nInternalFaces_ << endl;
        Info << "Patch Starts: " << patchStarts_ << endl;
        Info << "Patch Sizes: " << patchSizes_ << endl;
        Info << "=================" << endl;
    }

    // Reorder the points
    if (debug) Info << "ReOrdering points..." << endl;
    reOrderPoints(points);

    // Reorder the cells
    if (debug) Info << "ReOrdering cells..." << endl;
    reOrderCells();

    // Reorder the faces
    if (debug) Info << "ReOrdering faces..." << endl;
    reOrderFaces(faces, owner, neighbour);

    // Reorder the edges
    if (!twoDMesh_)
    {
        if (debug) Info << "ReOrdering edges..." << endl;
        reOrderEdges();
    }
}

// Copy edge-based connectivity from HashLists
void dynamicTopoFvMesh::setEdgeConnectivity()
{
    IOedges_.setSize(nEdges_);
    forAllIter(HashList<edge>::iterator, edges_, eIter)
    {
        IOedges_[eIter.index()] = eIter();
    }

    IOpointEdges_.setSize(nPoints_);
    forAllIter(HashList<labelList>::iterator, pointEdges_, peIter)
    {
        IOpointEdges_[peIter.index()] = peIter();
    }

    IOedgeFaces_.setSize(nEdges_);
    forAllIter(HashList<labelList>::iterator, edgeFaces_, efIter)
    {
        IOedgeFaces_[efIter.index()] = efIter();
    }

    IOfaceEdges_.setSize(nFaces_);
    forAllIter(HashList<labelList>::iterator, faceEdges_, feIter)
    {
        IOfaceEdges_[feIter.index()] = feIter();
    }

    IOedgePatchSizes_ = edgePatchSizes_;
    IOedgePatchStarts_ = edgePatchStarts_;

    IOedges_.writeOpt() = IOobject::AUTO_WRITE;
    IOedges_.instance() = time().timeName();

    IOpointEdges_.writeOpt() = IOobject::AUTO_WRITE;
    IOpointEdges_.instance() = time().timeName();

    IOedgeFaces_.writeOpt() = IOobject::AUTO_WRITE;
    IOedgeFaces_.instance() = time().timeName();

    IOfaceEdges_.writeOpt() = IOobject::AUTO_WRITE;
    IOfaceEdges_.instance() = time().timeName();

    IOedgePatchSizes_.writeOpt() = IOobject::AUTO_WRITE;
    IOedgePatchSizes_.instance() = time().timeName();

    IOedgePatchStarts_.writeOpt() = IOobject::AUTO_WRITE;
    IOedgePatchStarts_.instance() = time().timeName();
}

// Check the state of edge-based connectivity HashLists
void dynamicTopoFvMesh::checkEdgeConnectivity()
{
    Info << "Checking edge-face connectivity...";
    labelList nEdgeFaces(nEdges_, 0);
    forAllIter(HashList<labelList>::iterator, faceEdges_, feIter)
    {
        // Check consistency of face-edge-points as well
        edgeList eList = faces_[feIter.index()].edges();

        labelList& faceEdges = feIter();

        forAll(faceEdges,edgeI)
        {
            nEdgeFaces[faceEdges[edgeI]]++;

            // Check if this edge actually belongs to this face
            bool found = false;
            edge& edgeToCheck = edges_[faceEdges[edgeI]];

            forAll(eList, edgeII)
            {
                if (edgeToCheck == eList[edgeII])
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                Info << "Edge: " << faceEdges[edgeI] << ": " << edgeToCheck
                     << " Was not found in face: " << feIter.index()
                     << ": " << faces_[feIter.index()] << endl;
                FatalErrorIn("dynamicTopoFvMesh::checkEdgeConnectivity()") << nl
                     << "Edge-Face connectivity is inconsistent." << nl
                     << abort(FatalError);
            }
        }
    }
    forAllIter(HashList<labelList>::iterator, edgeFaces_, efIter)
    {
        labelList& edgeFaces = efIter();
        if (edgeFaces.size() != nEdgeFaces[efIter.index()])
        {
            Info << "Edge: " << efIter.index()
                 << " Old index: " << edgeMap_[efIter.index()] << nl
                 << "edgeFaces: " << edgeFaces << endl;
            FatalErrorIn("dynamicTopoFvMesh::checkEdgeConnectivity()") << nl
                 << "Edge-Face connectivity is inconsistent." << nl
                 << abort(FatalError);
        }
    }
    Info << "Done." << endl;

    Info << "Checking point-edge connectivity...";
    labelList nPointEdges(nPoints_, 0);
    forAllIter(HashList<edge>::iterator, edges_, eIter)
    {
        nPointEdges[eIter()[0]]++;
        nPointEdges[eIter()[1]]++;
    }
    forAllIter(HashList<labelList>::iterator, pointEdges_, peIter)
    {
        labelList& pointEdges = peIter();
        if (pointEdges.size() != nPointEdges[peIter.index()])
        {
            Info << "Point: " << peIter.index()
                 << " Old index: " << pointMap_[peIter.index()] << nl
                 << "pointEdges: " << pointEdges << endl;
            FatalErrorIn("dynamicTopoFvMesh::checkEdgeConnectivity()") << nl
                 << "Point-Edge connectivity is inconsistent."
                 << abort(FatalError);
        }
    }
    Info << "Done." << endl;
}

// Check tet-specific connectivity structures
void dynamicTopoFvMesh::checkTetConnectivity()
{
    // Loop through all cells and construct cell-to-node
    label cIndex = 0;
    labelList cellIndex(nCells_);
    List<labelHashSet> cellToNode(nCells_);

    forAllIter(HashList<cell>, cells_, cIter)
    {
        label cellI = cIter.index();
        cellIndex[cIndex] = cellI;
        cell& thisCell = cells_[cellI];

        forAll(thisCell, faceI)
        {
            labelList& fEdges = faceEdges_[thisCell[faceI]];

            forAll(fEdges, edgeI)
            {
                edge& thisEdge = edges_[fEdges[edgeI]];

                if (!cellToNode[cIndex].found(thisEdge[0]))
                {
                    cellToNode[cIndex].insert(thisEdge[0]);
                }

                if (!cellToNode[cIndex].found(thisEdge[1]))
                {
                    cellToNode[cIndex].insert(thisEdge[1]);
                }
            }
        }

        cIndex++;
    }

    // Preliminary check for size
    label nFailedChecks = 0;
    forAll(cellToNode, cellI)
    {
        if (cellToNode[cellI].size() != 4)
        {
            Info << "Warning: Cell: "
                 << cellIndex[cellI] << " is inconsistent. " << endl;

            cell& failedCell = cells_[cellIndex[cellI]];
            Info << "Cell faces: " << failedCell << endl;
            forAll(failedCell, faceI)
            {
                Info << "\tFace: " << failedCell[faceI]
                     << " :: " << faces_[failedCell[faceI]] << endl;

                labelList& fEdges = faceEdges_[failedCell[faceI]];
                forAll(fEdges, edgeI)
                {
                    Info << "\t\tEdge: " << fEdges[edgeI]
                         << " :: " << edges_[fEdges[edgeI]] << endl;
                }
            }

            nFailedChecks++;
        }
    }

    if (nFailedChecks)
    {
        FatalErrorIn
        (
            "dynamicTopoFvMesh::checkTetConnectivity()"
        )
            << nFailedChecks << " cells failed connectivity checks."
            << abort(FatalError);
    }
}

// Calculate the edge length-scale for the mesh
void dynamicTopoFvMesh::calculateLengthScale()
{
    if (edgeModification_)
    {
        label level = 1, visitedCells = 0;
        labelList cellLevels(nCells(),0);
        scalarField lengthScale(nCells(),0.0);

        // Obtain the cellCells addressing list
        const labelListList& cc = cellCells();

        // Obtain the list of patches for which the length-scale is fixed
        wordList toc = fixedLengthScalePatches_.toc();

        // Loop through all boundaries and mark adjacent cells
        const polyBoundaryMesh& bdy = boundaryMesh();
        const labelList& own = faceOwner();
        const pointField& pList = points();
        forAll(bdy,patchI)
        {
            const polyPatch& bdyPatch = bdy[patchI];
            // Loop through all fixed length-scale patches first
            bool fixed = false;
            forAll(toc,wordI)
            {
                word& patchName = toc[wordI];
                if (bdy[patchI].name() == patchName)
                {
                    label pStart = bdyPatch.start();
                    forAll(bdyPatch,faceI)
                    {
                        label ownCell = own[pStart+faceI];
                        if (cellLevels[ownCell] == 0)
                        {
                            cellLevels[ownCell] = level;
                            lengthScale[ownCell] =
                                fixedLengthScalePatches_[patchName][0].scalarToken();
                            visitedCells++;
                        }
                    }
                    fixed = true; break;
                }
            }

            // Set boundary patch size if no fixed-length scale is specified.
            if
            (
                (toc.size() == 0)
                && (bdyPatch.type() != "wedge")
                && (bdyPatch.type() != "empty")
                && (bdyPatch.type() != "symmetryPlane")
            )
            {
                label pStart = bdyPatch.start();
                forAll(bdyPatch,faceI)
                {
                    label ownCell = own[pStart+faceI];

                    if (cellLevels[ownCell] == 0)
                    {
                        cellLevels[ownCell] = level;

                        if (twoDMesh_)
                        {
                            // Get length-scale from edgeToWatch
                            edge& e = edgeToWatch_[pStart+faceI];
                            lengthScale[ownCell] = mag(pList[e[0]] - pList[e[1]]);
                        }
                        else
                        {
                            // Average edge-lengths for this face
                            scalar edgeLength = 0.0;
                            labelList& fEdges = faceEdges_[pStart+faceI];
                            forAll(fEdges, edgeI)
                            {
                                edge& e = edges_[fEdges[edgeI]];
                                edgeLength += mag(pList[e[0]] - pList[e[1]]);
                            }
                            lengthScale[ownCell] = (edgeLength/fEdges.size());
                        }

                        visitedCells++;
                    }
                }
            }
        }

        // Perform multiple sweeps through the mesh...
        while (visitedCells < nCells())
        {
            // Loop through cells, and increment neighbour
            // cells of the current level
            forAll(cellLevels,cellI)
            {
                if (cellLevels[cellI] == level)
                {
                    // Obtain the cells neighbouring this one
                    const labelList& cList = cc[cellI];
                    forAll(cList, indexI)
                    {
                        label& ngbLevel = cellLevels[cList[indexI]];
                        if (ngbLevel == 0)
                        {
                            ngbLevel = level+1;

                            // Compute the mean of the existing
                            // neighbour length-scales
                            const labelList& ncList = cc[cList[indexI]];
                            scalar sumLength = 0.0;
                            label nTouchedNgb = 0;
                            forAll(ncList, indexJ)
                            {
                                label sLevel = cellLevels[ncList[indexJ]];
                                if ((sLevel < ngbLevel) && (sLevel > 0))
                                {
                                    sumLength += lengthScale[ncList[indexJ]];
                                    nTouchedNgb++;
                                }
                            }
                            sumLength /= nTouchedNgb;

                            // Scale the length and assign to this cell
                            scalar sLength = sumLength*growthFactor_;
                            sLength = (sLength < maxLengthScale_)
                                    ? sLength : maxLengthScale_;
                            lengthScale[cList[indexI]] = sLength;
                            visitedCells++;
                        }
                    }
                }
            }

            // Move on to the next level
            level++;
        }

        if (debug)
        {
            Info << "Max Length Scale: " << maxLengthScale_ << endl;
            Info << "Length Scale sweeps: " << level << endl;
        }

        // Check if everything went okay
        if (visitedCells != nCells())
        {
            FatalErrorIn("dynamicTopoFvMesh::calculateLengthScale()")
                    << " Algorithm did not visit every cell in the mesh."
                    << " Something's messed up." << nl
                    << " Visited cells: " << visitedCells
                    << " nCells: " << nCells()
                    << abort(FatalError);
        }

        // Copy the most recent length-scale values
        forAllIter(HashList<scalar>, lengthScale_, lIter)
        {
            lIter() = lengthScale[lIter.index()];
        }
    }
}

// Return the appropriate length-scale for boundary face
scalar dynamicTopoFvMesh::boundaryLengthScale
(
    const label faceIndex
)
{
    label bFacePatch = whichPatch(faceIndex);

    // Loop through all fixed length-scale patches, and return the fixed value
    wordList toc = fixedLengthScalePatches_.toc();
    forAll(toc,wordI)
    {
        word& patchName = toc[wordI];

        if (boundaryMesh()[bFacePatch].name() == patchName)
        {
            return (fixedLengthScalePatches_[patchName][0].scalarToken());
        }
    }

    return lengthScale_[owner_[faceIndex]];
}

// 2D Edge-swapping engine
void dynamicTopoFvMesh::swap2DEdges(void *argument)
{
    // Recast the argument
    topoMeshStruct *thread = reinterpret_cast<topoMeshStruct*>(argument);
    dynamicTopoFvMesh *mesh = thread->mesh_;

    bool found, foundinner;
    label otherPointIndex = -1, nextPoint = -1;
    FixedList<label,2> cellLabels(-1);
    FixedList<label,2> c0BdyIndex, c0IntIndex, c1BdyIndex, c1IntIndex;
    FixedList<face,2>  c0BdyFace,  c0IntFace,  c1BdyFace,  c1IntFace;
    FixedList<label,4> commonFaceIndex;
    FixedList<face,4>  commonFaces;
    FixedList<edge,2>  commonEdges;

    // Pick items off the stack
    while (!mesh->faceStack(mesh->getThreadID(pthread_self())).empty())
    {
        // Retrieve the index for this face
        label findex = mesh->faceStack(mesh->getThreadID(pthread_self())).pop();

        // Get the two cells on either side...
        mesh->faceCells(findex, cellLabels);
        label c0 = cellLabels[0];
        label c1 = cellLabels[1];

        // Consider only internal faces..
        if (c1 == -1) continue;

        // Find the interior/boundary faces.
        mesh->findPrismFaces
        (
            findex,
            c0,
            c0BdyFace,
            c0BdyIndex,
            c0IntFace,
            c0IntIndex
        );

        mesh->findPrismFaces
        (
            findex,
            c1,
            c1BdyFace,
            c1BdyIndex,
            c1IntFace,
            c1IntIndex
        );

        // Find the common faces / edges on the boundary
        // At the end of this loop, commonFaces [0] & [1] share commonEdge [0]
        // and commonFaces [2] & [3] share commonEdge [1]
        // Also, commonFaces[0] & [2] are connected to cell[0],
        // and commonFaces[1] & [3] are connected to cell[1]
        found = false; foundinner = false;
        edgeList e1 = c0BdyFace[0].edges(), e2 = c1BdyFace[0].edges();
        edgeList e3 = c0BdyFace[1].edges(), e4 = c1BdyFace[1].edges();
        forAll(e1,edgeI)
        {
            forAll(e2,edgeJ)
            {
                if (e1[edgeI] == e2[edgeJ])
                {
                    // These two faces share an edge, store for posterity
                    commonFaces[0] = c0BdyFace[0];
                    commonFaces[1] = c1BdyFace[0];
                    commonFaces[2] = c0BdyFace[1];
                    commonFaces[3] = c1BdyFace[1];
                    commonFaceIndex[0] = c0BdyIndex[0];
                    commonFaceIndex[1] = c1BdyIndex[0];
                    commonFaceIndex[2] = c0BdyIndex[1];
                    commonFaceIndex[3] = c1BdyIndex[1];
                    commonEdges[0] = e1[edgeI];
                    forAll(e3,edgeK)
                    {
                        forAll(e4,edgeL)
                        {
                            if (e3[edgeK] == e4[edgeL])
                            {
                                commonEdges[1] = e3[edgeK];
                                foundinner = true; break;
                            }
                        }
                        if (foundinner) break;
                    }
                    found = true; break;
                }
            }
            if (found) break;
        }
        if (!found)
        {
            // A match was obviously not found before,
            // but we know the common faces now.
            commonFaces[0] = c0BdyFace[0];
            commonFaces[1] = c1BdyFace[1];
            commonFaces[2] = c0BdyFace[1];
            commonFaces[3] = c1BdyFace[0];
            commonFaceIndex[0] = c0BdyIndex[0];
            commonFaceIndex[1] = c1BdyIndex[1];
            commonFaceIndex[2] = c0BdyIndex[1];
            commonFaceIndex[3] = c1BdyIndex[0];
            // Start a new search for common edges
            forAll(e1,edgeI)
            {
                forAll(e4,edgeJ)
                {
                    if (e1[edgeI] == e4[edgeJ])
                    {
                        commonEdges[0] = e1[edgeI];
                        forAll(e2,edgeK)
                        {
                            forAll(e3,edgeL)
                            {
                                if (e2[edgeK] == e3[edgeL])
                                {
                                    commonEdges[1] = e2[edgeK];
                                    foundinner = true; break;
                                }
                            }
                            if (foundinner) break;
                        }
                        found = true; break;
                    }
                }
                if (found) break;
            }
        }

        // Construct a circle passing through the three
        // points of the first face, define its radius
        label zeroIndex  = commonFaces[0][0];
        label oneIndex   = commonFaces[0][1];
        label twoIndex   = commonFaces[0][2];
        point& pointZero = mesh->meshPoints()[zeroIndex];
        point& pointOne  = mesh->meshPoints()[oneIndex];
        point& pointTwo  = mesh->meshPoints()[twoIndex];

        // Determine the circumcenter
        point center = mesh->circumCenter
                       (
                           pointZero,
                           pointOne,
                           pointTwo,
                           zeroIndex,
                           oneIndex,
                           twoIndex
                       );

        scalar radius = (pointZero - center)&(pointZero - center);

        // Find the isolated point on the other face, and the point next to it
        mesh->findIsolatedPoint
        (
            commonFaces[1],
            commonEdges[0],
            otherPointIndex,
            nextPoint
        );

        // ...and determine whether it lies in this circle
        point otherPoint = mesh->meshPoints()[otherPointIndex];

        if
        (
            ((otherPoint - center)&(otherPoint - center)) < radius
        )
        {
            mesh->swapQuadFace
            (
                findex,
                commonFaceIndex,
                commonFaces,
                commonEdges
            );
        }
    }
}

// Initialize the length-scale field
void dynamicTopoFvMesh::initLengthScale()
{
    lengthScale_.setSize(nCells_, 0.0);

    if (twoDMesh_)
    {
        // Allocate fields
        edgeToWatch_.setSize(nFaces_,edge(0,0));

        // Loop through all quad-faces and build initial edge-lengths
        bool found;
        for(label findex=0; findex<nFaces_; findex++)
        {
            face& fi = faces_[findex];
            if (fi.size() == 4)
            {
                label c0 = owner_[findex];
                cell& cell_0 = cells_[c0];
                edgeList efi = fi.edges();

                // Look for a triangular face on the boundary
                found=false;
                forAll(cell_0, i)
                {
                    label faceIndex = cell_0[i];
                    face& fj=faces_[faceIndex];
                    if (neighbour_[faceIndex] == -1 && fj.size() == 3)
                    {
                        // Match an edge, and compute its length
                        edgeList efj = fj.edges();
                        forAll(efi, indexI)
                        {
                            forAll(efj, indexJ)
                            {
                                if (efi[indexI] == efj[indexJ])
                                {
                                    edgeToWatch_[findex] = efi[indexI];
                                    found=true; break;
                                }
                            }
                            if (found) break;
                        }
                        break;
                    }
                }
            }
        }
    }
}

// Initialize edge related connectivity lists
void dynamicTopoFvMesh::initEdges()
{
    if (IOedges_.headerOk())
    {
        // Connectivity is already read in from disk, so copy to HashLists.
        forAll(IOedges_, edgeI)
        {
            edges_.append(IOedges_[edgeI]);
        }

        forAll(IOpointEdges_, pointI)
        {
            pointEdges_.append(IOpointEdges_[pointI]);
        }

        forAll(IOedgeFaces_, edgeI)
        {
            edgeFaces_.append(IOedgeFaces_[edgeI]);
        }

        forAll(IOfaceEdges_, faceI)
        {
            faceEdges_.append(IOfaceEdges_[faceI]);
        }

        edgePatchSizes_ = IOedgePatchSizes_;
        edgePatchStarts_ = IOedgePatchStarts_;

        // Obtain nEdges_ and nInternalEdges_ from patches and set sizes
        label lastPatch = numPatches_-1;
        nEdges_ = edgePatchStarts_[lastPatch] + edgePatchSizes_[lastPatch];
        nInternalEdges_ = edgePatchStarts_[0];
        reverseEdgeMap_.setSize(nEdges_);

        // Check to see that everything is okay
#       ifdef FULLDEBUG
        checkEdgeConnectivity();
#       endif
    }
    else
    {
        // Set sizes first.
        nEdges_ = primitiveMesh::nEdges();
        reverseEdgeMap_.setSize(nEdges_);

        // Obtain connectivity from primitive mesh
        const edgeList& edges = primitiveMesh::edges();
        const labelListList& pEdges = primitiveMesh::pointEdges();
        const labelListList& fEdges = primitiveMesh::faceEdges();
        const labelListList& eFaces = primitiveMesh::edgeFaces();

        // Allocate lists for re-ordering
        labelList edgePatch(nEdges_, -1);

        // Edge-patches are the same as faces
        for(label i = nInternalFaces_; i < nFaces_; i++)
        {
            const labelList& fEdge = fEdges[i];
            forAll(fEdge, edgeI)
            {
                edgePatch[fEdge[edgeI]] = whichPatch(i);
            }
        }

        // Loop through edgePatch and renumber internal edges
        forAll(edgePatch, edgeI)
        {
            if (edgePatch[edgeI] == -1)
            {
                reverseEdgeMap_[edgeI] = nInternalEdges_++;
            }
            else
            {
                edgePatchSizes_[edgePatch[edgeI]]++;
            }
        }

        // Calculate patch-starts
        label startCount = nInternalEdges_;
        forAll(edgePatchStarts_, patchI)
        {
            edgePatchStarts_[patchI] = startCount;
            startCount += edgePatchSizes_[patchI];
        }

        // Now renumber boundary edges
        labelList patchCount(edgePatchStarts_);
        forAll(edgePatch, edgeI)
        {
            if (edgePatch[edgeI] >= 0)
            {
                reverseEdgeMap_[edgeI] = patchCount[edgePatch[edgeI]]++;
            }
        }

        // Renumber and fill in faceEdges
        forAll(fEdges, faceI)
        {
            const labelList& fEdge = fEdges[faceI];
            labelList renumberFaceEdge(fEdge.size(),-1);

            forAll(fEdge, edgeI)
            {
                renumberFaceEdge[edgeI] = reverseEdgeMap_[fEdge[edgeI]];
            }

            faceEdges_.append(renumberFaceEdge);
        }

        // Renumber and fill in pointEdges
        forAll(pEdges, pointI)
        {
            const labelList& pEdge = pEdges[pointI];
            labelList renumberPointEdges(pEdge.size(), -1);

            forAll(pEdge, edgeI)
            {
                renumberPointEdges[edgeI] = reverseEdgeMap_[pEdge[edgeI]];
            }

            pointEdges_.append(renumberPointEdges);
        }

        // Renumber and fill in edges and edgeFaces
        edges_.setSize(nEdges_, edge(-1,-1));
        edgeFaces_.setSize(nEdges_, labelList(0));

        forAll(edges, edgeI)
        {
            edges_[reverseEdgeMap_[edgeI]] = edges[edgeI];
            edgeFaces_[reverseEdgeMap_[edgeI]] = eFaces[edgeI];
        }

        // Now that connectivity is constructed, copy and set the output-option
        setEdgeConnectivity();
    }
}

// Initialize mutex lists
void dynamicTopoFvMesh::initMutexLists()
{
    if (threader_->multiThreaded())
    {
        pointMutex_.clear();
        pointMutex_.setSize(nPoints_);

        edgeMutex_.clear();
        edgeMutex_.setSize(nEdges_);

        faceMutex_.clear();
        faceMutex_.setSize(nFaces_);

        cellMutex_.clear();
        cellMutex_.setSize(nCells_);
    }
}

// Unlock mutex lists with HashSets
inline void dynamicTopoFvMesh::unlockMutexLists
(
    labelHashSet& pLocks,
    labelHashSet& eLocks,
    labelHashSet& fLocks,
    labelHashSet& cLocks
)
{
    if (threader_->multiThreaded())
    {
        labelList lockedCells = cLocks.toc();
        forAll(lockedCells, cellI)
        {
            cellMutex_[lockedCells[cellI]].unlock();
        }
        cLocks.clear();

        labelList lockedFaces = fLocks.toc();
        forAll(lockedFaces, faceI)
        {
            faceMutex_[lockedFaces[faceI]].unlock();
        }
        fLocks.clear();

        labelList lockedEdges = eLocks.toc();
        forAll(lockedEdges, edgeI)
        {
            edgeMutex_[lockedEdges[edgeI]].unlock();
        }
        eLocks.clear();

        labelList lockedPoints = pLocks.toc();
        forAll(lockedPoints, pointI)
        {
            pointMutex_[lockedPoints[pointI]].unlock();
        }
        pLocks.clear();
    }
}

// Return a reference to the multiThreader
const multiThreader& dynamicTopoFvMesh::threader() const
{
    return threader_();
}

// Push items on to the stack
inline void dynamicTopoFvMesh::stack::push(const label index)
{
    stackMutex_.lock();

    if (!stack_.found(index))
    {
        stack_.insert(index);
    }

    stackMutex_.unlock();
}

//- Insert item onto stack (no checking)
inline void dynamicTopoFvMesh::stack::insert(const label index)
{
    stack_.insert(index);
}

// Pop an item off the stack
inline label dynamicTopoFvMesh::stack::pop()
{
    stackMutex_.lock();

    const label index = stack_.begin().key();

    stack_.erase(index);

    stackMutex_.unlock();

    return index;
}

// Remove a specific item off the stack
inline void dynamicTopoFvMesh::stack::remove(const label index)
{
    stackMutex_.lock();

    if (stack_.found(index))
    {
        stack_.erase(index);
    }

    stackMutex_.unlock();
}

// Return if the stack is empty or not
inline bool dynamicTopoFvMesh::stack::empty()
{
    return (stack_.size() == 0);
}

//- Return the size of the stack
inline label dynamicTopoFvMesh::stack::size()
{
    return stack_.size();
}

//- Print out the stack
inline void dynamicTopoFvMesh::stack::print()
{
    Info << stack_ << endl;
}

// Return length-scale at an edge-location in the mesh [3D]
inline scalar dynamicTopoFvMesh::meshEdgeLengthScale
(
    const label eIndex
)
{
    scalar scale = 0.0;
    labelList& eFaces = edgeFaces_[eIndex];

    // Determine whether the edge is internal
    if (whichEdgePatch(eIndex) < 0)
    {
        forAll(eFaces, faceI)
        {
            scale += lengthScale_[owner_[eFaces[faceI]]];
            scale += lengthScale_[neighbour_[eFaces[faceI]]];
        }

        scale /= (2.0*eFaces.size());
    }
    else
    {
        // Search for boundary faces, and average their scale
        forAll(eFaces, faceI)
        {
            if (neighbour_[eFaces[faceI]] == -1)
            {
                scale += boundaryLengthScale(eFaces[faceI]);
            }
        }

        scale /= 2.0;
    }

    return scale;
}

// Check if a given face is a quad
inline bool dynamicTopoFvMesh::checkQuadFace(const label fIndex)
{
    return (faces_[fIndex].size() == 4);
}

// Return the length of an edge
inline scalar dynamicTopoFvMesh::edgeLength(const label eIndex)
{
    edge& thisEdge = edges_[eIndex];

    return mag(meshPoints_[thisEdge[1]] - meshPoints_[thisEdge[0]]);
}

// Return cell-labels on either side of the face
inline void dynamicTopoFvMesh::faceCells
(
    const label fIndex,
    FixedList<label,2> cellLabels
)
{
    cellLabels[0] = owner_[fIndex];
    cellLabels[1] = neighbour_[fIndex];
}

// Return length-scale at an cell-location in the mesh
inline scalar dynamicTopoFvMesh::meshCellLengthScale
(
    const label cIndex
)
{
    return lengthScale_[cIndex];
}

// Does the mesh perform edge-modification?
bool dynamicTopoFvMesh::edgeModification()
{
    return edgeModification_;
}

// 2D Edge-bisection/collapse engine
void dynamicTopoFvMesh::edgeBisectCollapse2D
(
    void *argument
)
{
    // Loop through all quad-faces and bisect/collapse
    // edges (quad-faces) by the criterion:
    // Bisect when boundary edge-length > ratioMax_*originalLength
    // Collapse when boundary edge-length < ratioMin_*originalLength

    // Recast the argument
    topoMeshStruct *thread = reinterpret_cast<topoMeshStruct*>(argument);
    dynamicTopoFvMesh *mesh = thread->mesh_;

    // Cell labels for faces
    FixedList<label,2> cellLabels(-1);

    // Pick items off the stack
    while (!mesh->faceStack(mesh->getThreadID(pthread_self())).empty())
    {
        // Retrieve the index for this face
        label findex = mesh->faceStack(mesh->getThreadID(pthread_self())).pop();

        // Select only quad-faces
        if (mesh->checkQuadFace(findex))
        {
            // Measure the boundary edge-length of the face in question
            edge& checkEdge = mesh->edgeToWatch()[findex];
            point& a = mesh->meshPoints()[checkEdge[0]];
            point& b = mesh->meshPoints()[checkEdge[1]];
            scalar length = mag(b-a);

            // Get the two cells on either side...
            mesh->faceCells(findex, cellLabels);
            label c0 = cellLabels[0];
            label c1 = cellLabels[1];

            // Determine the length-scale at this face
            scalar scale=0;
            if (c1 == -1)
            {
                scale = mesh->boundaryLengthScale(findex);

                // Check if this boundary face is adjacent to a sliver-cell,
                // and remove it by a two-step bisection/collapse operation.
                bool sliverRemoved = mesh->remove2DSliver
                                     (
                                         findex
                                     );

                if (sliverRemoved)
                {
                    // Set the flag
                    mesh->topoChangeFlag() = true;

                    continue;
                }
            }
            else
            {
                scale = 0.5*
                        (
                            mesh->meshCellLengthScale(c0)
                          + mesh->meshCellLengthScale(c1)
                        );
            }

            //== Edge Bisection ==//
            if(length > mesh->ratioMax()*scale)
            {
                // Set the flag
                mesh->topoChangeFlag() = true;

                // Bisect this face
                mesh->bisectQuadFace(findex);
            }
            else
            //== Edge Collapse ==//
            if(length < mesh->ratioMin()*scale)
            {
                // Collapse this face
                bool success = mesh->collapseQuadFace(findex);

                if (success)
                {
                    // Set the flag
                    mesh->topoChangeFlag() = true;
                }
            }
        }
    }
}

// 3D Edge-swapping engine
void dynamicTopoFvMesh::swap3DEdges
(
    void *argument
)
{
    // Recast the argument
    topoMeshStruct *thread = reinterpret_cast<topoMeshStruct*>(argument);
    dynamicTopoFvMesh *mesh = thread->mesh_;

    // Obtain maxTetsPerEdge
    label mMax = mesh->maxTetsPerEdge();

    // Hull variables
    scalar minQuality;
    DynamicList<label> cellHull(mMax);
    DynamicList<label> faceHull(mMax);
    DynamicList<label> vertexHull(mMax);
    labelHashSet pLocks, eLocks, fLocks, cLocks;

    // Dynamic programming variables
    scalarListList Q;
    labelListList K, triangulations;

    // Allocate dynamic programming tables
    mesh->initTables(mMax, Q, K, triangulations);

    // Pick edges off the stack
    while (!mesh->edgeStack(mesh->getThreadID(pthread_self())).empty())
    {
        // Retrieve an edge from the stack
        label eIndex = mesh->edgeStack(mesh->getThreadID(pthread_self())).pop();

        // Check if this edge is on a bounding curve
        if (mesh->checkBoundingCurve(eIndex))
        {
            continue;
        }

        // Obtain a ring of vertices around this edge, and try to lock it
        while
        (
            mesh->constructVertexRing
            (
                eIndex,
                cellHull,
                faceHull,
                vertexHull,
                minQuality,
                pLocks,
                eLocks,
                fLocks,
                cLocks
            )
        )
        {
            if
            (
                !mesh->edgeStack
                (
                    mesh->getThreadID(pthread_self())
                ).empty()
            )
            {
                label oldIndex = eIndex;

                // Pop another edge
                eIndex = mesh->edgeStack
                         (
                             mesh->getThreadID(pthread_self())
                         ).pop();

                // Put the old edge back on the stack
                mesh->edgeStack
                (
                    mesh->getThreadID(pthread_self())
                ).push(oldIndex);
            }

            // Clear lists
            cellHull.clear(); faceHull.clear(); vertexHull.clear();
        }

        // Check if a table-resize is necessary
        if (vertexHull.size() > mMax)
        {
            if (mesh->allowTableResize())
            {
                // Resize the tables to account for
                // more tets per edge
                mMax = vertexHull.size();
                Q.clear(); K.clear(); triangulations.clear();
                mesh->initTables(mMax, Q, K, triangulations);
            }
            else
            {
                // Move on to the next edge
                mesh->unlockMutexLists(pLocks, eLocks, fLocks, cLocks);
                cellHull.clear(); faceHull.clear(); vertexHull.clear();
                continue;
            }
        }

        // Fill the dynamic programming tables
        label m = mesh->fillTables(eIndex, vertexHull, minQuality, Q, K);

        // Remove this edge if necessary...
        if (Q[0][m-1] > minQuality)
        {
#           ifdef FULLDEBUG
            if (debug)
            {
                Info << endl;
                Info << "Old Hull Quality: " << minQuality << endl;
                Info << "New Hull Quality: " << Q[0][m-1] << endl;
            }
#           endif

            // Try to acquire write priority for this hull
            if (mesh->threader().multiThreaded())
            {
                // Unlock all entities (from read-lock)
                mesh->unlockMutexLists(pLocks, eLocks, fLocks, cLocks);

                if
                (
                    mesh->obtainWritePriority
                    (
                        cellHull,
                        pLocks,
                        eLocks,
                        fLocks,
                        cLocks
                    )
                )
                {
                    // This thread failed in contention. Push the edge back.
                    mesh->edgeStack
                    (
                        mesh->getThreadID(pthread_self())
                    ).push(eIndex);

                    cellHull.clear(); faceHull.clear(); vertexHull.clear();

                    continue;
                }
            }

            // Remove this edge according to the swap sequence
            mesh->removeEdgeFlips
            (
                m,
                eIndex,
                K,
                cellHull,
                faceHull,
                vertexHull,
                triangulations
            );

            // Remove hull entities from the write-lock list.
            if (mesh->threader().multiThreaded())
            {
                forAll(cellHull, cellI)
                {
                    cLocks.erase(cellHull[cellI]);
                }

                forAll(faceHull, faceI)
                {
                    fLocks.erase(faceHull[faceI]);
                }

                eLocks.erase(eIndex);
            }

            // Set the flag
            mesh->topoChangeFlag() = true;
        }

        // Move on to the next edge.
        cellHull.clear(); faceHull.clear(); vertexHull.clear();

        // Unlock all entities
        mesh->unlockMutexLists(pLocks, eLocks, fLocks, cLocks);
    }
}

// 3D Edge-bisection/collapse engine
void dynamicTopoFvMesh::edgeBisectCollapse3D
(
    void *argument
)
{
    // Loop through all edges and bisect/collapse by the criterion:
    // Bisect when edge-length > ratioMax_*originalLength
    // Collapse when edge-length < ratioMin_*originalLength

    // Recast the argument
    topoMeshStruct *thread = reinterpret_cast<topoMeshStruct*>(argument);
    dynamicTopoFvMesh *mesh = thread->mesh_;

    while (!mesh->edgeStack(mesh->getThreadID(pthread_self())).empty())
    {
        // Retrieve an edge from the stack
        label eIndex = mesh->edgeStack(mesh->getThreadID(pthread_self())).pop();

        // Measure the edge-length
        scalar length = mesh->edgeLength(eIndex);

        // Determine the length-scale at this point in the mesh
        scalar scale = mesh->meshEdgeLengthScale(eIndex);

        //== Edge Bisection ==//
        if(length > mesh->ratioMax()*scale)
        {
            // Set the flag
            mesh->topoChangeFlag() = true;

            // Bisect this edge
            mesh->bisectEdge(eIndex);
        }
        else
        //== Edge Collapse ==//
        if(length < mesh->ratioMin()*scale)
        {
            // Collapse this edge
            bool success = mesh->collapseEdge(eIndex);

            // The edge can safely be deleted, since the iterator points
            // to the next valid edge on the list.
            if (success)
            {
                // Set the flag
                mesh->topoChangeFlag() = true;
            }
        }
    }
}

// Return the face-stack
inline dynamicTopoFvMesh::stack& dynamicTopoFvMesh::faceStack
(
    const label threadID
)
{
    return faceStack_[threadID];
}

// Return the edge-stack for a particular thread
inline dynamicTopoFvMesh::stack& dynamicTopoFvMesh::edgeStack
(
    const label threadID
)
{
    return edgeStack_[threadID];
}

// Return the integer threadID for a given pthreadID
inline label dynamicTopoFvMesh::getThreadID
(
    const pthread_t& pthreadID
)
{
    for (label i = 0; i < threader_->getNumThreads(); i++)
    {
        if (pthread_equal(structPtr_[i].pthreadID_, pthreadID))
        {
            return i;
        }
    }

    return -1;
}

// Initialize edge-stacks
inline void dynamicTopoFvMesh::initEdgeStacks()
{
    label tID = 0;

    for
    (
        HashList<edge>::iterator iter = edges_.begin();
        iter != edges_.end();
        iter++
    )
    {
        edgeStack_[tID].insert(iter.index());
        tID = edgeStack_.fcIndex(tID);
    }
}

// Initialize face-stacks
inline void dynamicTopoFvMesh::initFaceStacks()
{
    label tID = 0;

    for
    (
        HashList<face>::iterator iter = faces_.begin();
        iter != faces_.end();
        iter++
    )
    {
        faceStack_[tID].insert(iter.index());
        tID = faceStack_.fcIndex(tID);
    }
}

// Method for the swapping of a quad-face in 2D
void dynamicTopoFvMesh::swapQuadFace
(
    const label findex,
    const FixedList<label,4>& commonFaceIndex,
    const FixedList<face,4>&  commonFaces,
    const FixedList<edge,2>&  commonEdges
)
{
    face f;
    edge firstEdge(0,0);
    bool found = false;
    FixedList<label,4> otherPointIndex(-1), nextToOtherPoint(-1);
    FixedList<label,2> c0BdyIndex, c0IntIndex, c1BdyIndex, c1IntIndex;
    FixedList<face,2>  c0BdyFace,  c0IntFace,  c1BdyFace,  c1IntFace;
    FixedList<label,4> commonIntFaceIndex;
    FixedList<face,4>  commonIntFaces;

    // Obtain a reference for this face...
    face& thisFace = faces_[findex];

    // Get the two cells on either side...
    label c0 = owner_[findex];
    label c1 = neighbour_[findex];

    // Get cell references
    cell &cell_0 = cells_[c0];
    cell &cell_1 = cells_[c1];

    // This face needs to be flipped...
    if (debug)
    {
        Info << nl << nl << "Face: " << findex
             << " needs to be flipped. " << endl;
        Info << "Cell[0]: " << c0 << ": " << cell_0 << endl;
        Info << "Cell[1]: " << c1 << ": " << cell_1 << endl;
        Info << "Common Faces: Set 1: "
             << commonFaceIndex[0] << ": " << commonFaces[0] << ", "
             << commonFaceIndex[1] << ": " << commonFaces[1] << endl;
        Info << "Common Faces: Set 2: "
             << commonFaceIndex[2] << ": " << commonFaces[2] << ", "
             << commonFaceIndex[3] << ": " << commonFaces[3] << endl;
        Info << "Old face: " << faces_[findex] << endl;
    }

    // Find the interior/boundary faces.
    findPrismFaces
    (
        findex,
        c0,
        c0BdyFace,
        c0BdyIndex,
        c0IntFace,
        c0IntIndex
    );

    findPrismFaces
    (
        findex,
        c1,
        c1BdyFace,
        c1BdyIndex,
        c1IntFace,
        c1IntIndex
    );

    // Set the flag
    topoChangeFlag_ = true;

    // Find the points that don't lie on shared edges
    // and the points next to them (for orientation)
    findIsolatedPoint
    (
        commonFaces[1],
        commonEdges[0],
        otherPointIndex[1],
        nextToOtherPoint[1]
    );

    findIsolatedPoint
    (
        commonFaces[0],
        commonEdges[0],
        otherPointIndex[0],
        nextToOtherPoint[0]
    );

    findIsolatedPoint
    (
        commonFaces[2],
        commonEdges[1],
        otherPointIndex[2],
        nextToOtherPoint[2]
    );

    findIsolatedPoint
    (
        commonFaces[3],
        commonEdges[1],
        otherPointIndex[3],
        nextToOtherPoint[3]
    );

    // Find the other two edges on the face being flipped
    // First edge detected belongs to cell[0] by default
    edgeList eThis = thisFace.edges();
    forAll(eThis,edgeI)
    {
        if
        (
            eThis[edgeI] != commonEdges[0]
         && eThis[edgeI] != commonEdges[1]
        )
        {
            if
            (
                eThis[edgeI][0] == nextToOtherPoint[0]
             || eThis[edgeI][1] == nextToOtherPoint[0]
            )
            {
                firstEdge = eThis[edgeI];
            }
        }
    }

    // Find the interior faces that share the first edge
    // At the end of this loop, commonIntFaces [0] & [1] share firstEdge
    // and commonIntFaces [2] & [3] share the secondEdge,
    // where [0],[2] lie on cell[0] and [1],[3] lie on cell[1]
    found = false;
    edgeList e1 = c0IntFace[0].edges();
    forAll(e1,edgeI)
    {
        if ( e1[edgeI] == firstEdge )
        {
            commonIntFaces[0] = c0IntFace[0];
            commonIntFaces[2] = c0IntFace[1];
            commonIntFaceIndex[0] = c0IntIndex[0];
            commonIntFaceIndex[2] = c0IntIndex[1];
            found = true; break;
        }
    }
    if (!found)
    {
        // The edge was obviously not found before
        commonIntFaces[0] = c0IntFace[1];
        commonIntFaces[2] = c0IntFace[0];
        commonIntFaceIndex[0] = c0IntIndex[1];
        commonIntFaceIndex[2] = c0IntIndex[0];
    }

    found = false;
    edgeList e3 = c1IntFace[0].edges();
    forAll(e3,edgeI)
    {
        if ( e3[edgeI] == firstEdge )
        {
            commonIntFaces[1] = c1IntFace[0];
            commonIntFaces[3] = c1IntFace[1];
            commonIntFaceIndex[1] = c1IntIndex[0];
            commonIntFaceIndex[3] = c1IntIndex[1];
            found = true; break;
        }
    }
    if (!found)
    {
        // The edge was obviously not found before
        commonIntFaces[1] = c1IntFace[1];
        commonIntFaces[3] = c1IntFace[0];
        commonIntFaceIndex[1] = c1IntIndex[1];
        commonIntFaceIndex[3] = c1IntIndex[0];
    }

    // Modify the five faces belonging to this hull
    face& newFace = faces_[findex];
    face& newBdyFace0 = faces_[commonFaceIndex[0]];
    face& newBdyFace1 = faces_[commonFaceIndex[1]];
    face& newBdyFace2 = faces_[commonFaceIndex[2]];
    face& newBdyFace3 = faces_[commonFaceIndex[3]];
    label c0count=0, c1count=0;

    // Define parameters for the new flipped face
    newFace[0] = otherPointIndex[0];
    newFace[1] = otherPointIndex[1];
    newFace[2] = otherPointIndex[3];
    newFace[3] = otherPointIndex[2];
    cell_0[c0count++] = findex;
    cell_1[c1count++] = findex;
    owner_[findex] = c0;
    neighbour_[findex] = c1;

    // Modify the edge-to-watch
    if (edgeModification_)
    {
        edge& edgeToModify = edgeToWatch_[findex];
        edgeToModify[0] = otherPointIndex[0];
        edgeToModify[1] = otherPointIndex[1];
    }

    // Four modified boundary faces need to be constructed,
    // but right-handedness is also important.
    // Take a cue from the existing boundary-face orientation

    // Zeroth boundary face - Owner c[0], Neighbour -1
    newBdyFace0[0] = otherPointIndex[0];
    newBdyFace0[1] = nextToOtherPoint[0];
    newBdyFace0[2] = otherPointIndex[1];
    cell_0[c0count++] = commonFaceIndex[0];
    owner_[commonFaceIndex[0]] = c0;
    neighbour_[commonFaceIndex[0]] = -1;

    // First boundary face - Owner c[1], Neighbour -1
    newBdyFace1[0] = otherPointIndex[1];
    newBdyFace1[1] = nextToOtherPoint[1];
    newBdyFace1[2] = otherPointIndex[0];
    cell_1[c1count++] = commonFaceIndex[1];
    owner_[commonFaceIndex[1]] = c1;
    neighbour_[commonFaceIndex[1]] = -1;

    // Second boundary face - Owner c[0], Neighbour -1
    newBdyFace2[0] = otherPointIndex[3];
    newBdyFace2[1] = nextToOtherPoint[3];
    newBdyFace2[2] = otherPointIndex[2];
    cell_0[c0count++] = commonFaceIndex[2];
    owner_[commonFaceIndex[2]] = c0;
    neighbour_[commonFaceIndex[2]] = -1;

    // Third boundary face - Owner c[1], Neighbour -1
    newBdyFace3[0] = otherPointIndex[2];
    newBdyFace3[1] = nextToOtherPoint[2];
    newBdyFace3[2] = otherPointIndex[3];
    cell_1[c1count++] = commonFaceIndex[3];
    owner_[commonFaceIndex[3]] = c1;
    neighbour_[commonFaceIndex[3]] = -1;

    if (debug)
    {
        Info << "New flipped face: " << newFace << endl;
        Info << "New boundary face[0]" << commonFaceIndex[0]
             << ": " << newBdyFace0 << endl;
        Info << "New boundary face[1]" << commonFaceIndex[1]
             << ": " << newBdyFace1 << endl;
        Info << "New boundary face[2]" << commonFaceIndex[2]
             << ": " << newBdyFace2 << endl;
        Info << "New boundary face[3]" << commonFaceIndex[3]
             << ": " << newBdyFace3 << endl;
    }

    // Check the orientation of the two quad faces, and modify as necessary
    label newOwn=0, newNei=0;

    // The quad face belonging to cell[1] now becomes a part of cell[0]
    if ( neighbour_[commonIntFaceIndex[1]] == -1 )
    {
        // Boundary face
        // Face doesn't need to be flipped, just update the owner
        f          = commonIntFaces[1];
        newOwn     = c0;
        newNei     = -1;
    }
    else
    if ( owner_[commonIntFaceIndex[1]] == c1 )
    {
        // This face is on the interior, check for previous owner
        // Upper-triangular ordering has to be maintained, however...
        if ( c0 > neighbour_[commonIntFaceIndex[1]] )
        {
            // Flip is necessary
            f          = commonIntFaces[1].reverseFace();
            newOwn     = neighbour_[commonIntFaceIndex[1]];
            newNei     = c0;
        }
        else
        {
            // Flip isn't necessary, just change the owner
            f          = commonIntFaces[1];
            newOwn     = c0;
            newNei     = neighbour_[commonIntFaceIndex[1]];
        }
    }
    else
    if ( neighbour_[commonIntFaceIndex[1]] == c1 )
    {
        // This face is on the interior, check for previous neighbour
        // Upper-triangular ordering has to be maintained, however...
        if ( c0 < owner_[commonIntFaceIndex[1]] )
        {
            // Flip is necessary
            f          = commonIntFaces[1].reverseFace();
            newOwn     = c0;
            newNei     = owner_[commonIntFaceIndex[1]];
        }
        else
        {
            // Flip isn't necessary, just change the neighbour
            f          = commonIntFaces[1];
            newOwn     = owner_[commonIntFaceIndex[1]];
            newNei     = c0;
        }
    }

    faces_[commonIntFaceIndex[1]] = f;
    cell_0[c0count++] = commonIntFaceIndex[0];
    cell_0[c0count++] = commonIntFaceIndex[1];
    owner_[commonIntFaceIndex[1]] = newOwn;
    neighbour_[commonIntFaceIndex[1]] = newNei;

    // The quad face belonging to cell[0] now becomes a part of cell[1]
    if ( neighbour_[commonIntFaceIndex[2]] == -1 )
    {
        // Boundary face
        // Face doesn't need to be flipped, just update the owner
        f          = commonIntFaces[2];
        newOwn     = c1;
        newNei     = -1;
    }
    else
    if ( owner_[commonIntFaceIndex[2]] == c0 )
    {
        // This face is on the interior, check for previous owner
        // Upper-triangular ordering has to be maintained, however...
        if ( c1 > neighbour_[commonIntFaceIndex[2]] )
        {
            // Flip is necessary
            f          = commonIntFaces[2].reverseFace();
            newOwn     = neighbour_[commonIntFaceIndex[2]];
            newNei     = c1;
        }
        else
        {
            // Flip isn't necessary, just change the owner
            f          = commonIntFaces[2];
            newOwn     = c1;
            newNei     = neighbour_[commonIntFaceIndex[2]];
        }
    }
    else
    if ( neighbour_[commonIntFaceIndex[2]] == c0 )
    {
        // This face is on the interior, check for previous neighbour
        // Upper-triangular ordering has to be maintained, however...
        if ( c1 < owner_[commonIntFaceIndex[2]] )
        {
            // Flip is necessary
            f          = commonIntFaces[2].reverseFace();
            newOwn     = c1;
            newNei     = owner_[commonIntFaceIndex[2]];
        }
        else
        {
            // Flip isn't necessary, just change the neighbour
            f          = commonIntFaces[2];
            newOwn     = owner_[commonIntFaceIndex[2]];
            newNei     = c1;
        }
    }

    faces_[commonIntFaceIndex[2]] = f;
    cell_1[c1count++] = commonIntFaceIndex[2];
    cell_1[c1count++] = commonIntFaceIndex[3];
    owner_[commonIntFaceIndex[2]] = newOwn;
    neighbour_[commonIntFaceIndex[2]] = newNei;

    // Generate mapping information for both cells
    label firstParent, secondParent;
    const labelListList& cc = cellCells();
    labelHashSet c0MasterObjects(6);
    labelHashSet c1MasterObjects(6);

    if (c0 < nOldCells_)
    {
        firstParent = c0;
    }
    else
    {
        firstParent = cellParents_[c0];
    }

    if (c1 < nOldCells_)
    {
        secondParent = c1;
    }
    else
    {
        secondParent = cellParents_[c1];
    }

    // Find the cell's neighbours in the old mesh
    c0MasterObjects.insert(firstParent);
    c1MasterObjects.insert(firstParent);

    forAll(cc[firstParent],cellI)
    {
        if (!c0MasterObjects.found(cc[firstParent][cellI]))
        {
            c0MasterObjects.insert(cc[firstParent][cellI]);
        }
        if (!c1MasterObjects.found(cc[firstParent][cellI]))
        {
            c1MasterObjects.insert(cc[firstParent][cellI]);
        }
    }

    c0MasterObjects.insert(secondParent);
    c1MasterObjects.insert(secondParent);

    forAll(cc[secondParent],cellI)
    {
        if (!c0MasterObjects.found(cc[secondParent][cellI]))
        {
            c0MasterObjects.insert(cc[secondParent][cellI]);
        }
        if (!c1MasterObjects.found(cc[secondParent][cellI]))
        {
            c1MasterObjects.insert(cc[secondParent][cellI]);
        }
    }

    // Insert mapping info into the HashTable
    cellsFromCells_.insert(c0,objectMap(c0,c0MasterObjects.toc()));
    cellsFromCells_.insert(c1,objectMap(c1,c1MasterObjects.toc()));
}

// Method for the bisection of a quad-face in 2D
void dynamicTopoFvMesh::bisectQuadFace
(
    const label findex
)
{
    // Local variables
    bool found;
    label replaceFace;
    FixedList<edge,2> commonEdges;
    FixedList<label,4> otherPointIndex(-1), nextToOtherPoint(-1);
    FixedList<label,2> c0BdyIndex, c0IntIndex, c1BdyIndex, c1IntIndex;
    FixedList<face,2>  c0BdyFace,  c0IntFace,  c1BdyFace,  c1IntFace;
    edge tmpEdge(0,0), firstEdge(0,0), secondEdge(0,0);

    // Obtain a reference for this face...
    face& thisFace = faces_[findex];

    // Get the two cells on either side...
    label c0 = owner_[findex], c1 = neighbour_[findex];

    // Find the prism faces for cell[0].
    cell &cell_0 = cells_[c0];
    findPrismFaces
    (
        findex,
        c0,
        c0BdyFace,
        c0BdyIndex,
        c0IntFace,
        c0IntIndex
    );

    if (debug)
    {
        Info << nl << nl << "Face: " << findex
             << ": " << thisFace << " is to be bisected. " << endl;
        Info << "Cell[0]: " << c0 << ": " << cell_0 << endl;
        forAll(cell_0, faceI)
        {
            Info << cell_0[faceI] << ": " << faces_[cell_0[faceI]] << endl;
        }
    }

    // Find the common-edge between the triangular boundary faces
    // and the face under consideration.
    findCommonEdge(c0BdyFace[0],thisFace,commonEdges[0]);
    findCommonEdge(c0BdyFace[1],thisFace,commonEdges[1]);

    // Find the isolated point on both boundary faces of cell[0]
    findIsolatedPoint
    (
        c0BdyFace[0],
        commonEdges[0],
        otherPointIndex[0],
        nextToOtherPoint[0]
    );

    findIsolatedPoint
    (
        c0BdyFace[1],
        commonEdges[1],
        otherPointIndex[1],
        nextToOtherPoint[1]
    );

    // Add two new points to the end of the list
    label newPtIndex0 = meshPoints_.append
                        (
                            0.5*
                            (
                                meshPoints_[commonEdges[0][0]]
                              + meshPoints_[commonEdges[0][1]]
                            )
                        );

    label newPtIndex1 = meshPoints_.append
                        (
                            0.5*
                            (
                                meshPoints_[commonEdges[1][0]]
                              + meshPoints_[commonEdges[1][1]]
                            )
                        );
    nPoints_ += 2;

    // Add a new prism cell to the end of the list
    label newCellIndex0 = cells_.append(cell(5));
    cell &newCell0 = cells_[newCellIndex0];

    // Generate mapping information for this new cell
    label firstParent;
    const labelListList& cc = cellCells();
    labelHashSet c0MasterObjects;

    if (c0 < nOldCells_)
    {
        firstParent = c0;
    }
    else
    {
        firstParent = cellParents_[c0];
    }

    // Insert the parent cell
    cellParents_.insert(newCellIndex0,firstParent);

    // Find the cell's neighbours in the old mesh
    c0MasterObjects.insert(firstParent);
    forAll(cc[firstParent],cellI)
    {
        if (!c0MasterObjects.found(cc[firstParent][cellI]))
        {
            c0MasterObjects.insert(cc[firstParent][cellI]);
        }
    }

    // Insert mapping info into the HashTable
    cellsFromCells_.insert
    (
        newCellIndex0,
        objectMap
        (
            newCellIndex0,
            c0MasterObjects.toc()
        )
    );

    // Add a new element to the lengthScale field
    // (Currently the same as cell[0])
    lengthScale_.append(lengthScale_[c0]);

    // Modify the two existing triangle boundary faces

    // Zeroth boundary face - Owner = c[0] & Neighbour [-1] (unchanged)
    face& BdyFace0 = faces_[c0BdyIndex[0]];
    BdyFace0[0] = otherPointIndex[0];
    BdyFace0[1] = nextToOtherPoint[0];
    BdyFace0[2] = newPtIndex0;

    // First boundary face - Owner = newCell[0], Neighbour = -1
    face& BdyFace1 = faces_[c0BdyIndex[1]];
    BdyFace1[0] = otherPointIndex[1];
    BdyFace1[1] = nextToOtherPoint[1];
    BdyFace1[2] = newPtIndex1;
    owner_[c0BdyIndex[1]] = newCellIndex0;
    replaceLabel(c0BdyIndex[1],-1,cell_0);

    // Detect edges other than commonEdges
    edgeList eThis = thisFace.edges();
    forAll(eThis, edgeI)
    {
        if (eThis[edgeI] != commonEdges[0] && eThis[edgeI] != commonEdges[1])
        {
            if
            (
                eThis[edgeI][0] == nextToOtherPoint[0]
             || eThis[edgeI][1] == nextToOtherPoint[0]
            )
            {
                firstEdge = eThis[edgeI];
            }

            if
            (
                eThis[edgeI][0] == nextToOtherPoint[1]
             || eThis[edgeI][1] == nextToOtherPoint[1]
            )
            {
                secondEdge = eThis[edgeI];
            }
        }
    }

    // Modify point-labels on the quad face under consideration
    replaceLabel
    (
        commonEdges[0].otherVertex(nextToOtherPoint[0]),
        newPtIndex0,
        thisFace
    );

    replaceLabel
    (
        nextToOtherPoint[1],
        newPtIndex1,
        thisFace
    );

    // Change the edge-length criteria for this face
    tmpEdge[0] = nextToOtherPoint[0];
    tmpEdge[1] = newPtIndex0;
    edgeToWatch_[findex] = tmpEdge;

    if (debug)
    {
        Info << "Modified thisFace: " << findex << ": " << thisFace << endl;
    }

    // Find the interior face that contains secondEdge
    found = false;
    edgeList e1 = c0IntFace[0].edges();
    forAll(e1, edgeI)
    {
        if ( e1[edgeI] == secondEdge )
        {
            replaceLabel(c0IntIndex[0],-1,cell_0);
            replaceFace = c0IntIndex[0];
            found = true; break;
        }
    }
    if (!found)
    {
        // The edge was obviously not found before
        replaceLabel(c0IntIndex[1],-1,cell_0);
        replaceFace = c0IntIndex[1];
    }

    // Check if face reversal is necessary for the replacement
    if (owner_[replaceFace] == c0)
    {
        if (neighbour_[replaceFace] == -1)
        {
            // Change the owner
            owner_[replaceFace] = newCellIndex0;
        }
        else
        {
            // This face has to be reversed
            faces_[replaceFace] = faces_[replaceFace].reverseFace();
            owner_[replaceFace] = neighbour_[replaceFace];
            neighbour_[replaceFace] = newCellIndex0;
        }
    }
    else
    {
        // Keep owner, but change neighbour
        neighbour_[replaceFace] = newCellIndex0;
    }

    // Define the faces for the new cell
    newCell0[0] = c0BdyIndex[1];
    newCell0[1] = replaceFace;

    // Define the set of new faces and insert...

    // Temporary data
    label newFaceIndex;
    face tmpQuadFace(4), tmpTriFace(3);
    edge edgeToWatch;

    // New interior face; Owner = cell[0] & Neighbour = newCell[0]
    tmpQuadFace[0] = otherPointIndex[0];
    tmpQuadFace[1] = newPtIndex0;
    tmpQuadFace[2] = newPtIndex1;
    tmpQuadFace[3] = otherPointIndex[1];
    edgeToWatch[0] = otherPointIndex[0];
    edgeToWatch[1] = newPtIndex0;
    newFaceIndex = insertFace
                   (
                       -1,
                       tmpQuadFace,
                       c0,
                       newCellIndex0,
                       edgeToWatch
                   );
    replaceLabel(-1, newFaceIndex, newCell0);
    replaceLabel(-1, newFaceIndex, cell_0);

    // remove2DSliver requires this face index for removal
    bisectInteriorFace_ = newFaceIndex;

    // Second boundary face; Owner = newCell[0] & Neighbour = [-1]
    tmpTriFace[0] = otherPointIndex[0];
    tmpTriFace[1] = newPtIndex0;
    tmpTriFace[2] = commonEdges[0].otherVertex(nextToOtherPoint[0]);
    edgeToWatch[0] = edgeToWatch[1] = 0;
    newFaceIndex = insertFace
                   (
                       whichPatch(c0BdyIndex[0]),
                       tmpTriFace,
                       newCellIndex0,
                       -1,
                       edgeToWatch
                   );
    replaceLabel(-1, newFaceIndex, newCell0);

    // Third boundary face; Owner = c[0] & Neighbour = [-1]
    tmpTriFace[0] = otherPointIndex[1];
    tmpTriFace[1] = newPtIndex1;
    tmpTriFace[2] = commonEdges[1].otherVertex(nextToOtherPoint[1]);
    edgeToWatch[0] = edgeToWatch[1] = 0;
    newFaceIndex = insertFace
                   (
                       whichPatch(c0BdyIndex[1]),
                       tmpTriFace,
                       c0,
                       -1,
                       edgeToWatch
                   );
    replaceLabel(-1, newFaceIndex, cell_0);

    if (c1 == -1)
    {
        // The quad boundary face resulting from bisection;
        // Owner = newCell[0] & Neighbour = [-1]
        tmpQuadFace[0] = newPtIndex1;
        tmpQuadFace[1] = nextToOtherPoint[1];
        tmpQuadFace[2] = commonEdges[0].otherVertex(nextToOtherPoint[0]);
        tmpQuadFace[3] = newPtIndex0;
        edgeToWatch[0] = newPtIndex1;
        edgeToWatch[1] = nextToOtherPoint[1];
        newFaceIndex = insertFace
                       (
                           whichPatch(findex),
                           tmpQuadFace,
                           newCellIndex0,
                           -1,
                           edgeToWatch
                       );
        replaceLabel(-1, newFaceIndex, newCell0);

        if (debug)
        {
            Info << "Modified Cell[0]: " << c0 << ": " << cell_0 << endl;
            forAll(cell_0, faceI)
            {
                Info << cell_0[faceI]
                     << ": " << faces_[cell_0[faceI]] << endl;
            }
            Info << "New Cell[0]: " << newCellIndex0 << ": " << newCell0 << endl;
            forAll(newCell0, faceI)
            {
                Info << newCell0[faceI]
                     << ": " << faces_[newCell0[faceI]] << endl;
            }
        }
    }
    else
    {
        cell &cell_1 = cells_[c1];

        // Find the prism faces for cell[1].
        findPrismFaces
        (
            findex,
            c1,
            c1BdyFace,
            c1BdyIndex,
            c1IntFace,
            c1IntIndex
        );

        // Add a new prism cell to the end of the list
        label newCellIndex1 = cells_.append(cell(5));
        cell &newCell1 = cells_[newCellIndex1];

        // Generate mapping information for this new cell
        label secondParent;
        labelHashSet c1MasterObjects;

        if (c1 < nOldCells_)
        {
            secondParent = c1;
        }
        else
        {
            secondParent = cellParents_[c1];
        }

        // Insert the parent cell
        cellParents_.insert(newCellIndex1,secondParent);

        // Find the cell's neighbours in the old mesh
        c1MasterObjects.insert(secondParent);
        forAll(cc[secondParent],cellI)
        {
            if (!c1MasterObjects.found(cc[secondParent][cellI]))
            {
                c1MasterObjects.insert(cc[secondParent][cellI]);
            }
        }

        // Insert mapping info into the HashTable
        cellsFromCells_.insert
        (
            newCellIndex1,
            objectMap
            (
                newCellIndex1,
                c1MasterObjects.toc()
            )
        );

        // Add a new element to the lengthScale field
        // (Currently the same as cell[1])
        lengthScale_.append(lengthScale_[c1]);

        if (debug)
        {
            Info << "Cell[1]: " << c1 << ": " << cell_1 << endl;
            forAll(cell_1, faceI)
                Info << cell_1[faceI] << ": " << faces_[cell_1[faceI]] << endl;
        }

        // Find the interior face that contains secondEdge
        found = false;
        edgeList e2 = c1IntFace[0].edges();
        forAll(e2, edgeI)
        {
            if ( e2[edgeI] == secondEdge )
            {
                replaceLabel(c1IntIndex[0], -1, cell_1);
                replaceFace = c1IntIndex[0];
                found = true; break;
            }
        }
        if (!found)
        {
            // The edge was obviously not found before
            replaceLabel(c1IntIndex[1], -1, cell_1);
            replaceFace = c1IntIndex[1];
        }

        // Check if face reversal is necessary for the replacement
        if (owner_[replaceFace] == c1)
        {
            if (neighbour_[replaceFace] == -1)
            {
                // Change the owner
                owner_[replaceFace] = newCellIndex1;
            }
            else
            {
                // This face has to be reversed
                faces_[replaceFace] = faces_[replaceFace].reverseFace();
                owner_[replaceFace] = neighbour_[replaceFace];
                neighbour_[replaceFace] = newCellIndex1;
            }
        }
        else
        {
            // Keep owner, but change neighbour
            neighbour_[replaceFace] = newCellIndex1;
        }

        // Define attributes for the new prism cell
        newCell1[0] = replaceFace;

        // The interior face resulting from bisection;
        // Owner = newCell[0] & Neighbour = newCell[1]
        tmpQuadFace[0] = newPtIndex1;
        tmpQuadFace[1] = nextToOtherPoint[1];
        tmpQuadFace[2] = commonEdges[0].otherVertex(nextToOtherPoint[0]);
        tmpQuadFace[3] = newPtIndex0;
        edgeToWatch[0] = newPtIndex1;
        edgeToWatch[1] = nextToOtherPoint[1];
        newFaceIndex = insertFace
                       (
                           -1,
                           tmpQuadFace,
                           newCellIndex0,
                           newCellIndex1,
                           edgeToWatch
                       );
        replaceLabel(-1, newFaceIndex, newCell0);
        replaceLabel(-1, newFaceIndex, newCell1);
        newCell1[1] = newFaceIndex;

        // Check for common edges among the two boundary faces
        // Find the isolated point on both boundary faces of cell[1]
        label patch_0 = -2, patch_1 = -2;
        if(findCommonEdge(c1BdyFace[0],c0BdyFace[0],commonEdges[0]))
        {
            findCommonEdge(c1BdyFace[1],c0BdyFace[1],commonEdges[1]);

            findIsolatedPoint
            (
                c1BdyFace[0],
                commonEdges[0],
                otherPointIndex[2],
                nextToOtherPoint[2]
            );

            findIsolatedPoint
            (
                c1BdyFace[1],
                commonEdges[1],
                otherPointIndex[3],
                nextToOtherPoint[3]
            );

            // Modify the two existing triangle boundary faces
            // Zeroth boundary face - Owner = newCell[1], Neighbour = -1
            face& BdyFace2 = faces_[c1BdyIndex[0]];
            BdyFace2[0] = otherPointIndex[2];
            BdyFace2[1] = nextToOtherPoint[2];
            BdyFace2[2] = newPtIndex0;
            owner_[c1BdyIndex[0]] = newCellIndex1;
            replaceLabel(c1BdyIndex[0], -1, cell_1);
            newCell1[2] = c1BdyIndex[0];

            // First boundary face - Owner = c[1] & Neighbour [-1] (unchanged)
            face& BdyFace3 = faces_[c1BdyIndex[1]];
            BdyFace3[0] = otherPointIndex[3];
            BdyFace3[1] = nextToOtherPoint[3];
            BdyFace3[2] = newPtIndex1;

            // Obtain patch info
            patch_0 = whichPatch(c1BdyIndex[0]);
            patch_1 = whichPatch(c1BdyIndex[1]);
        }
        else
        {
            findCommonEdge(c1BdyFace[0],c0BdyFace[1],commonEdges[1]);
            findCommonEdge(c1BdyFace[1],c0BdyFace[0],commonEdges[0]);

            findIsolatedPoint
            (
                c1BdyFace[0],
                commonEdges[1],
                otherPointIndex[3],
                nextToOtherPoint[3]
            );

            findIsolatedPoint
            (
                c1BdyFace[1],
                commonEdges[0],
                otherPointIndex[2],
                nextToOtherPoint[2]
            );

            // Modify the two existing triangle boundary faces
            // Zeroth boundary face - Owner = newCell[1], Neighbour = -1
            face& BdyFace2 = faces_[c1BdyIndex[1]];
            BdyFace2[0] = otherPointIndex[2];
            BdyFace2[1] = nextToOtherPoint[2];
            BdyFace2[2] = newPtIndex0;
            owner_[c1BdyIndex[1]] = newCellIndex1;
            replaceLabel(c1BdyIndex[1], -1, cell_1);
            newCell1[2] = c1BdyIndex[1];

            // First boundary face - Owner = c[1] & Neighbour [-1] (unchanged)
            face& BdyFace3 = faces_[c1BdyIndex[0]];
            BdyFace3[0] = otherPointIndex[3];
            BdyFace3[1] = nextToOtherPoint[3];
            BdyFace3[2] = newPtIndex1;

            // Obtain patch info
            patch_0 = whichPatch(c1BdyIndex[1]);
            patch_1 = whichPatch(c1BdyIndex[0]);
        }

        // New interior face; Owner = cell[1] & Neighbour = newCell[1]
        tmpQuadFace[0] = newPtIndex0;
        tmpQuadFace[1] = otherPointIndex[2];
        tmpQuadFace[2] = otherPointIndex[3];
        tmpQuadFace[3] = newPtIndex1;
        edgeToWatch[0] = newPtIndex1;
        edgeToWatch[1] = otherPointIndex[3];
        newFaceIndex = insertFace
                       (
                           -1,
                           tmpQuadFace,
                           c1,
                           newCellIndex1,
                           edgeToWatch
                       );
        replaceLabel(-1, newFaceIndex, newCell1);
        replaceLabel(-1, newFaceIndex, cell_1);

        // Second boundary face; Owner = cell[1] & Neighbour [-1]
        tmpTriFace[0] = otherPointIndex[2];
        tmpTriFace[1] = newPtIndex0;
        tmpTriFace[2] = commonEdges[0].otherVertex(nextToOtherPoint[2]);
        edgeToWatch[0] = edgeToWatch[1] = 0;
        newFaceIndex = insertFace
                       (
                           patch_0,
                           tmpTriFace,
                           c1,
                           -1,
                           edgeToWatch
                       );
        replaceLabel(-1, newFaceIndex, cell_1);

        // Third boundary face; Owner = newCell[1] & Neighbour [-1]
        tmpTriFace[0] = otherPointIndex[3];
        tmpTriFace[1] = newPtIndex1;
        tmpTriFace[2] = commonEdges[1].otherVertex(nextToOtherPoint[3]);
        edgeToWatch[0] = edgeToWatch[1] = 0;
        newFaceIndex = insertFace
                       (
                           patch_1,
                           tmpTriFace,
                           newCellIndex1,
                           -1,
                           edgeToWatch
                       );
        replaceLabel(-1, newFaceIndex, newCell1);

        if (debug)
        {
            Info << nl << "Modified Cell[0]: " << c0 << ": " << cell_0 << endl;
            forAll(cell_0, faceI)
            {
                Info << cell_0[faceI]
                     << ": " << faces_[cell_0[faceI]] << endl;
            }

            Info << "New Cell[0]: " << newCellIndex0 << ": " << newCell0 << endl;
            forAll(newCell0, faceI)
            {
                Info << newCell0[faceI] << ": "
                     << faces_[newCell0[faceI]] << endl;
            }

            Info << nl << "Modified Cell[1]: " << c1 << ": " << cell_1 << endl;
            forAll(cell_1, faceI)
            {
                Info << cell_1[faceI] << ": "
                     << faces_[cell_1[faceI]] << endl;
            }

            Info << "New Cell[1]: " << newCellIndex1 << ": " << newCell1 << endl;
            forAll(newCell1, faceI)
            {
                Info << newCell1[faceI] << ": "
                     << faces_[newCell1[faceI]] << endl;
            }
        }
    }

    // Update the number of cells
    nCells_++;
    if (c1 != -1) nCells_++;
}

// Method for the collapse of a quad-face in 2D
// Returns a boolean value indicating whether the collapse was valid
bool dynamicTopoFvMesh::collapseQuadFace
(
    const label findex
)
{
    // Obtain a reference for this face...
    face& thisFace = faces_[findex];

    // This face is to be collapsed...
    if (debug)
    {
        Info << nl << nl
             << "Face: " << findex << ": " << thisFace
             << " is to be collapsed. " << endl;
    }

    // Local variables
    FixedList<label,2> c0BdyIndex, c0IntIndex, c1BdyIndex, c1IntIndex;
    FixedList<face,2> c0BdyFace, c0IntFace, c1BdyFace, c1IntFace;
    edge tmpEdge(0, 0), firstEdge, secondEdge;
    face tmpTriFace(3);

    // Find the two edges from checkEdge
    edge& checkEdge = edgeToWatch_[findex];
    edgeList thisEdges = thisFace.edges();
    forAll(thisEdges,edgeI)
    {
        if
        (
            (
                checkEdge[0] == thisEdges[edgeI][0]
             || checkEdge[0] == thisEdges[edgeI][1]
            )
         && !(checkEdge == thisEdges[edgeI])
        )
        {
            firstEdge = thisEdges[edgeI];
        }
        else
        if
        (
            (
                checkEdge[1] == thisEdges[edgeI][0]
             || checkEdge[1] == thisEdges[edgeI][1]
            )
         && !(checkEdge == thisEdges[edgeI])
        )
        {
            secondEdge = thisEdges[edgeI];
        }
        else
        if (!(checkEdge == thisEdges[edgeI]))
        {
            // This is the fourth edge
            tmpEdge = thisEdges[edgeI];
        }
    }

    // Build a hull of faces that are connected to each edge
    // This will also determine whether the edge lies on a boundary
    DynamicList<label> firstEdgeHull(10), secondEdgeHull(10);
    DynamicList<label> firstCells(10),    secondCells(10);
    DynamicList<label> firstTriFaces(10), secondTriFaces(10);
    DynamicList<label> firstEdgePatches(2), secondEdgePatches(2);

    bool firstEdgeBoundary  = constructPrismHull
                              (
                                  firstEdge,
                                  findex,
                                  firstEdgeHull,
                                  firstCells,
                                  firstTriFaces,
                                  firstEdgePatches,
                                  true
                              );

    bool secondEdgeBoundary = constructPrismHull
                              (
                                  secondEdge,
                                  findex,
                                  secondEdgeHull,
                                  secondCells,
                                  secondTriFaces,
                                  secondEdgePatches,
                                  true
                              );

    if (debug)
    {
        Info << endl;
        Info << "-------------------------" << endl;
        Info << "Hulls before modification" << endl;
        Info << "-------------------------" << endl;

        Info << nl << "Cells belonging to first Edge Hull: "
             << firstCells << endl;

        forAll(firstCells,cellI)
        {
            cell &firstCurCell = cells_[firstCells[cellI]];
            Info << "Cell: " << firstCells[cellI]
                 << ": " << firstCurCell << endl;
            forAll(firstCurCell,faceI)
            {
                Info << firstCurCell[faceI]
                     << ": " << faces_[firstCurCell[faceI]] << endl;
            }
        }

        Info << nl << "First Edge Hull: " << firstEdgeHull << endl;
        forAll(firstEdgeHull,indexI)
        {
            Info << firstEdgeHull[indexI]
                 << ": " << faces_[firstEdgeHull[indexI]] << endl;
        }

        Info << nl << "Cells belonging to second Edge Hull: "
             << secondCells << endl;

        forAll(secondCells, cellI)
        {
            cell &secondCurCell = cells_[secondCells[cellI]];
            Info << "Cell: " << secondCells[cellI]
                 << ": " << secondCurCell << endl;
            forAll(secondCurCell, faceI)
            {
                Info << secondCurCell[faceI]
                     << ": " << faces_[secondCurCell[faceI]] << endl;
            }
        }

        Info << nl << "Second Edge Hull: " << secondEdgeHull << endl;
        forAll(secondEdgeHull, indexI)
        {
            Info << secondEdgeHull[indexI]
                 << ": " << faces_[secondEdgeHull[indexI]] << endl;
        }
    }

    // Determine the common vertices for the first and second edges
    label cv0 = firstEdge.commonVertex(checkEdge);
    label cv1 = firstEdge.commonVertex(tmpEdge);
    label cv2 = secondEdge.commonVertex(checkEdge);
    label cv3 = secondEdge.commonVertex(tmpEdge);

    // Determine the neighbouring cells
    label c0 = owner_[findex], c1 = neighbour_[findex];

    // Find the prism-faces
    FixedList<label,2> faceToKeep(0), faceToThrow(0);

    findPrismFaces
    (
        findex,
        c0,
        c0BdyFace,
        c0BdyIndex,
        c0IntFace,
        c0IntIndex
    );

    if (c1 != -1)
    {
        findPrismFaces
        (
            findex,
            c1,
            c1BdyFace,
            c1BdyIndex,
            c1IntFace,
            c1IntIndex
        );
    }

    // Collapse preferentially towards a symmetryPlane.
    if (firstEdgeBoundary && secondEdgeBoundary)
    {
        WarningIn
        (
            "dynamicTopoFvMesh::collapseQuadFace"
            "(const label, face&)"
        )   << "Collapsing a face that lies on two boundary patches. "
            << "Algorithm will look for a symmetryPlane and collapse "
            << "the face preferentially towards it.\n"
            << "Face: " << findex << ": " << thisFace << endl;

        forAll(firstEdgePatches, patchI)
        {
            if
            (
                boundaryMesh()[firstEdgePatches[patchI]].type()
                == "symmetryPlane"
            )
            {
                secondEdgeBoundary = false;
            }
        }

        if (secondEdgeBoundary)
        {
            forAll(secondEdgePatches, patchI)
            {
                if
                (
                    boundaryMesh()[secondEdgePatches[patchI]].type()
                    == "symmetryPlane"
                )
                {
                    firstEdgeBoundary = false;
                }
            }
        }
    }

    if (!firstEdgeBoundary && secondEdgeBoundary)
    {
        // Check whether the collapse is possible.
        forAll(firstTriFaces, indexI)
        {
            if
            (
                (firstTriFaces[indexI] == c0BdyIndex[0])
             || (firstTriFaces[indexI] == c0BdyIndex[1])
            )
            {
                continue;
            }

            if (c1 != -1)
            {
                if
                (
                    (firstTriFaces[indexI] == c1BdyIndex[0])
                 || (firstTriFaces[indexI] == c1BdyIndex[1])
                )
                {
                    continue;
                }
            }

            face &triFace = faces_[firstTriFaces[indexI]];
            forAll(triFace, pointI)
            {
                tmpTriFace[pointI] = triFace[pointI];
                if (triFace[pointI] == cv0)
                {
                    tmpTriFace[pointI] = cv2;
                }
                if (triFace[pointI] == cv1)
                {
                    tmpTriFace[pointI] = cv3;
                }
            }

            // Compute the area and check if it's zero/negative
            scalar origArea = triFaceArea(triFace);
            scalar newArea  = triFaceArea(tmpTriFace);
            if
            (
                (Foam::sign(origArea) != Foam::sign(newArea))
             || (mag(newArea) < VSMALL)
            )
            {
                return false;
            }
        }

        // Collapse to the second node...
        forAll(firstEdgeHull,faceI)
        {
            face& replacementFace = faces_[firstEdgeHull[faceI]];
            replaceLabel(cv0,cv2,replacementFace);
            replaceLabel(cv1,cv3,replacementFace);

            // Modify edgeToWatch as well
            edge& replacementCheckEdge = edgeToWatch_[firstEdgeHull[faceI]];
            if (replacementCheckEdge[0] == cv0) replacementCheckEdge[0] = cv2;
            if (replacementCheckEdge[1] == cv0) replacementCheckEdge[1] = cv2;
            if (replacementCheckEdge[0] == cv1) replacementCheckEdge[0] = cv3;
            if (replacementCheckEdge[1] == cv1) replacementCheckEdge[1] = cv3;

            // Determine the quad-face in cell[0] & cell[1]
            // that belongs to firstEdgeHull
            if (firstEdgeHull[faceI] == c0IntIndex[0])
            {
                faceToKeep[0]  = c0IntIndex[1];
                faceToThrow[0] = c0IntIndex[0];
            }

            if (firstEdgeHull[faceI] == c0IntIndex[1])
            {
                faceToKeep[0]  = c0IntIndex[0];
                faceToThrow[0] = c0IntIndex[1];
            }

            if (c1 != -1)
            {
                if (firstEdgeHull[faceI] == c1IntIndex[0])
                {
                    faceToKeep[1]  = c1IntIndex[1];
                    faceToThrow[1] = c1IntIndex[0];
                }
                if (firstEdgeHull[faceI] == c1IntIndex[1])
                {
                    faceToKeep[1]  = c1IntIndex[0];
                    faceToThrow[1] = c1IntIndex[1];
                }
            }
        }

        // All triangular boundary faces also need to have point labels replaced
        forAll(firstCells,cellI)
        {
            cell& cellToCheck = cells_[firstCells[cellI]];
            forAll(cellToCheck,faceI)
            {
                face& faceToCheck = faces_[cellToCheck[faceI]];
                if (faceToCheck.size() == 3)
                {
                    forAll(faceToCheck,pointI)
                    {
                        if (faceToCheck[pointI] == cv0)
                        {
                            faceToCheck[pointI] = cv2;
                        }
                        if (faceToCheck[pointI] == cv1)
                        {
                            faceToCheck[pointI] = cv3;
                        }
                    }
                }
            }
        }

        // Delete the two points...
        meshPoints_.remove(cv0);
        meshPoints_.remove(cv1);
        nPoints_ -= 2;

        // Update the reverse point map
        if (cv0 < nOldPoints_) reversePointMap_[cv0] = -1;
        if (cv1 < nOldPoints_) reversePointMap_[cv1] = -1;
    }
    else
    {
        // Check whether the collapse is possible.
        forAll(secondTriFaces, indexI)
        {
            if
            (
                (secondTriFaces[indexI] == c0BdyIndex[0])
             || (secondTriFaces[indexI] == c0BdyIndex[1])
            )
            {
                continue;
            }

            if (c1 != -1)
            {
                if
                (
                    (secondTriFaces[indexI] == c1BdyIndex[0])
                 || (secondTriFaces[indexI] == c1BdyIndex[1])
                )
                {
                    continue;
                }
            }

            face &triFace = faces_[secondTriFaces[indexI]];
            forAll(triFace, pointI)
            {
                tmpTriFace[pointI] = triFace[pointI];
                if (triFace[pointI] == cv2)
                {
                    tmpTriFace[pointI] = cv0;
                }
                if (triFace[pointI] == cv3)
                {
                    tmpTriFace[pointI] = cv1;
                }
            }

            // Compute the area and check if it's zero/negative
            scalar origArea = triFaceArea(triFace);
            scalar newArea  = triFaceArea(tmpTriFace);
            if
            (
                (Foam::sign(origArea) != Foam::sign(newArea))
             || (mag(newArea) < VSMALL)
            )
            {
                return false;
            }
        }

        // Collapse to the first node by default...
        forAll(secondEdgeHull,faceI)
        {
            face& replacementFace = faces_[secondEdgeHull[faceI]];
            replaceLabel(cv2, cv0, replacementFace);
            replaceLabel(cv3, cv1, replacementFace);

            // Modify edgeToWatch as well
            edge& replacementCheckEdge = edgeToWatch_[secondEdgeHull[faceI]];
            if (replacementCheckEdge[0] == cv2) replacementCheckEdge[0] = cv0;
            if (replacementCheckEdge[1] == cv2) replacementCheckEdge[1] = cv0;
            if (replacementCheckEdge[0] == cv3) replacementCheckEdge[0] = cv1;
            if (replacementCheckEdge[1] == cv3) replacementCheckEdge[1] = cv1;

            // Determine the quad-face(s) in cell[0] & cell[1]
            // that belongs to secondEdgeHull
            if (secondEdgeHull[faceI] == c0IntIndex[0])
            {
                faceToKeep[0]  = c0IntIndex[1];
                faceToThrow[0] = c0IntIndex[0];
            }

            if (secondEdgeHull[faceI] == c0IntIndex[1])
            {
                faceToKeep[0]  = c0IntIndex[0];
                faceToThrow[0] = c0IntIndex[1];
            }

            if (c1 != -1)
            {
                if (secondEdgeHull[faceI] == c1IntIndex[0])
                {
                    faceToKeep[1]  = c1IntIndex[1];
                    faceToThrow[1] = c1IntIndex[0];
                }
                if (secondEdgeHull[faceI] == c1IntIndex[1])
                {
                    faceToKeep[1]  = c1IntIndex[0];
                    faceToThrow[1] = c1IntIndex[1];
                }
            }
        }
        // All triangular boundary faces also need to have point labels replaced
        forAll(secondCells, cellI)
        {
            cell& cellToCheck = cells_[secondCells[cellI]];
            forAll(cellToCheck, faceI)
            {
                face& faceToCheck = faces_[cellToCheck[faceI]];
                if (faceToCheck.size() == 3)
                {
                    forAll(faceToCheck, pointI)
                    {
                        if (faceToCheck[pointI] == cv2)
                        {
                            faceToCheck[pointI] = cv0;
                        }
                        if (faceToCheck[pointI] == cv3)
                        {
                            faceToCheck[pointI] = cv1;
                        }
                    }
                }
            }
        }

        // Delete the two points...
        meshPoints_.remove(cv2);
        meshPoints_.remove(cv3);
        nPoints_ -= 2;

        // Update the reverse point map
        if (cv2 < nOldPoints_) reversePointMap_[cv2] = -1;
        if (cv3 < nOldPoints_) reversePointMap_[cv3] = -1;
    }

    if (debug)
    {
        Info << endl;
        Info << "------------------------" << endl;
        Info << "Hulls after modification" << endl;
        Info << "------------------------" << endl;

        Info << nl << "Cells belonging to first Edge Hull: "
             << firstCells << endl;

        forAll(firstCells, cellI)
        {
            cell &firstCurCell = cells_[firstCells[cellI]];
            Info << "Cell: " << firstCells[cellI]
                 << ": " << firstCurCell << endl;
            forAll(firstCurCell, faceI)
            {
                Info << firstCurCell[faceI]
                     << ": " << faces_[firstCurCell[faceI]] << endl;
            }
        }

        Info << nl << "First Edge Hull: " << firstEdgeHull << endl;
        forAll(firstEdgeHull, indexI)
        {
            Info << firstEdgeHull[indexI]
                 << ": " << faces_[firstEdgeHull[indexI]] << endl;
        }

        Info << nl << "Cells belonging to second Edge Hull: "
             << secondCells << endl;

        forAll(secondCells, cellI)
        {
            cell &secondCurCell = cells_[secondCells[cellI]];
            Info << "Cell: " << secondCells[cellI]
                 << ": " << secondCurCell << endl;
            forAll(secondCurCell, faceI)
            {
                Info << secondCurCell[faceI]
                     << ": " << faces_[secondCurCell[faceI]] << endl;
            }
        }

        Info << nl << "Second Edge Hull: " << secondEdgeHull << endl;
        forAll(secondEdgeHull, indexI)
        {
            Info << secondEdgeHull[indexI]
                 << ": " << faces_[secondEdgeHull[indexI]] << endl;
        }

        Info << endl;
        Info << "Retained face: "
             << faceToKeep[0] << ": "
             << " owner: " << owner_[faceToKeep[0]]
             << " neighbour: " << neighbour_[faceToKeep[0]]
             << endl;
        Info << "Discarded face: "
             << faceToThrow[0] << ": "
             << " owner: " << owner_[faceToThrow[0]]
             << " neighbour: " << neighbour_[faceToThrow[0]]
             << endl;
        if (c1 != -1)
        {
            Info << "Retained face: "
                 << faceToKeep[1] << ": "
                 << " owner: " << owner_[faceToKeep[1]]
                 << " neighbour: " << neighbour_[faceToKeep[1]]
                 << endl;
            Info << "Discarded face: "
                 << faceToThrow[1] << ": "
                 << " owner: " << owner_[faceToThrow[1]]
                 << " neighbour: " << neighbour_[faceToThrow[1]]
                 << endl;
        }
    }

    // Ensure proper orientation for the two retained faces
    FixedList<label,2> cellCheck(0);
    if (owner_[faceToThrow[0]] == c0)
    {
        cellCheck[0] = neighbour_[faceToThrow[0]];
        if (owner_[faceToKeep[0]] == c0)
        {
            if
            (
                (neighbour_[faceToThrow[0]] > neighbour_[faceToKeep[0]])
             && (neighbour_[faceToKeep[0]] != -1)
            )
            {
                // This face is to be flipped
                faces_[faceToKeep[0]] = faces_[faceToKeep[0]].reverseFace();
                owner_[faceToKeep[0]] = neighbour_[faceToKeep[0]];
                neighbour_[faceToKeep[0]] = neighbour_[faceToThrow[0]];
            }
            else
            {
                if (neighbour_[faceToThrow[0]] != -1)
                {
                    // Keep orientation intact, and update the owner
                    owner_[faceToKeep[0]] = neighbour_[faceToThrow[0]];
                }
                else
                {
                    // This face will need to be flipped and converted
                    // to a boundary face. Flip it now, so that conversion
                    // happens later.
                    faces_[faceToKeep[0]] = faces_[faceToKeep[0]].reverseFace();
                    owner_[faceToKeep[0]] = neighbour_[faceToKeep[0]];
                    neighbour_[faceToKeep[0]] = -1;
                }
            }
        }
        else
        {
            // Keep orientation intact, and update the neighbour
            neighbour_[faceToKeep[0]] = neighbour_[faceToThrow[0]];
        }
    }
    else
    {
        cellCheck[0] = owner_[faceToThrow[0]];
        if (neighbour_[faceToKeep[0]] == c0)
        {
            if (owner_[faceToThrow[0]] < owner_[faceToKeep[0]])
            {
                // This face is to be flipped
                faces_[faceToKeep[0]] = faces_[faceToKeep[0]].reverseFace();
                neighbour_[faceToKeep[0]] = owner_[faceToKeep[0]];
                owner_[faceToKeep[0]] = owner_[faceToThrow[0]];
            }
            else
            {
                // Keep orientation intact, and update the neighbour
                neighbour_[faceToKeep[0]] = owner_[faceToThrow[0]];
            }
        }
        else
        {
            // Keep orientation intact, and update the owner
            owner_[faceToKeep[0]] = owner_[faceToThrow[0]];
        }
    }

    if (c1 != -1)
    {
        if (owner_[faceToThrow[1]] == c1)
        {
            cellCheck[1] = neighbour_[faceToThrow[1]];
            if (owner_[faceToKeep[1]] == c1)
            {
                if
                (
                    (neighbour_[faceToThrow[1]] > neighbour_[faceToKeep[1]])
                 && (neighbour_[faceToKeep[1]] != -1)
                )
                {
                    // This face is to be flipped
                    faces_[faceToKeep[1]] = faces_[faceToKeep[1]].reverseFace();
                    owner_[faceToKeep[1]] = neighbour_[faceToKeep[1]];
                    neighbour_[faceToKeep[1]] = neighbour_[faceToThrow[1]];
                }
                else
                {
                    if (neighbour_[faceToThrow[1]] != -1)
                    {
                        // Keep orientation intact, and update the owner
                        owner_[faceToKeep[1]] = neighbour_[faceToThrow[1]];
                    }
                    else
                    {
                        // This face will need to be flipped and converted
                        // to a boundary face. Flip it now, so that conversion
                        // happens later.
                        faces_[faceToKeep[1]] =
                            faces_[faceToKeep[1]].reverseFace();
                        owner_[faceToKeep[1]] = neighbour_[faceToKeep[1]];
                        neighbour_[faceToKeep[1]] = -1;
                    }
                }
            }
            else
            {
                // Keep orientation intact, and update the neighbour
                neighbour_[faceToKeep[1]] = neighbour_[faceToThrow[1]];
            }
        }
        else
        {
            cellCheck[1] = owner_[faceToThrow[1]];
            if (neighbour_[faceToKeep[1]] == c1)
            {
                if (owner_[faceToThrow[1]] < owner_[faceToKeep[1]])
                {
                    // This face is to be flipped
                    faces_[faceToKeep[1]] = faces_[faceToKeep[1]].reverseFace();
                    neighbour_[faceToKeep[1]] = owner_[faceToKeep[1]];
                    owner_[faceToKeep[1]] = owner_[faceToThrow[1]];
                }
                else
                {
                    // Keep orientation intact, and update the neighbour
                    neighbour_[faceToKeep[1]] = owner_[faceToThrow[1]];
                }
            }
            else
            {
                // Keep orientation intact, and update the owner
                owner_[faceToKeep[1]] = owner_[faceToThrow[1]];
            }
        }
    }

    // Remove orphaned faces
    if (owner_[faceToKeep[0]] == -1)
    {
        removeFace(faceToKeep[0]);
    }
    else
    if
    (
        (neighbour_[faceToKeep[0]] == -1)
     && (whichPatch(faceToKeep[0]) < 0)
    )
    {
        // This face is being converted from interior to boundary. Remove
        // from the interior list and add as a boundary face to the end.
        label newFaceIndex0 = insertFace
                              (
                                  whichPatch(faceToThrow[0]),
                                  faces_[faceToKeep[0]],
                                  owner_[faceToKeep[0]],
                                  -1,
                                  edgeToWatch_[faceToKeep[0]]
                              );

        replaceLabel
        (
            faceToKeep[0],
            newFaceIndex0,
            cells_[owner_[faceToKeep[0]]]
        );

        // Renumber the neighbour so that this face is removed correctly.
        neighbour_[faceToKeep[0]] = 0;
        removeFace(faceToKeep[0]);
    }

    // Remove the unwanted faces in the cell(s) adjacent to this face,
    // and correct the cells that contain discarded faces
    cell &cell_0 = cells_[c0];
    forAll(cell_0,faceI)
    {
        if (cell_0[faceI] != findex && cell_0[faceI] != faceToKeep[0])
        {
           removeFace(cell_0[faceI]);
        }
    }
    cells_.remove(c0);
    lengthScale_.remove(c0);
    if (cellCheck[0] != -1)
    {
        replaceLabel(faceToThrow[0], faceToKeep[0], cells_[cellCheck[0]]);
    }

    // Update the number of cells, and the reverse map
    nCells_--;
    if (c0 < nOldCells_)
    {
        reverseCellMap_[c0] = -1;
    }

    // Check if the cell was added in the current morph, and delete
    if (cellsFromCells_.found(c0))
    {
        cellsFromCells_.erase(c0);
    }

    if (c1 != -1)
    {
        // Remove orphaned faces
        if (owner_[faceToKeep[1]] == -1)
        {
            removeFace(faceToKeep[1]);
        }
        else
        if
        (
            (neighbour_[faceToKeep[1]] == -1)
         && (whichPatch(faceToKeep[1]) < 0)
        )
        {
            // This face is being converted from interior to boundary. Remove
            // from the interior list and add as a boundary face to the end.
            label newFaceIndex1 = insertFace
                                  (
                                      whichPatch(faceToThrow[1]),
                                      faces_[faceToKeep[1]],
                                      owner_[faceToKeep[1]],
                                      -1,
                                      edgeToWatch_[faceToKeep[1]]
                                  );

            replaceLabel
            (
                faceToKeep[1],
                newFaceIndex1,
                cells_[owner_[faceToKeep[1]]]
            );

            // Renumber the neighbour so that this face is removed correctly.
            neighbour_[faceToKeep[1]] = 0;
            removeFace(faceToKeep[1]);
        }

        cell &cell_1 = cells_[c1];
        forAll(cell_1, faceI)
        {
            if (cell_1[faceI] != findex && cell_1[faceI] != faceToKeep[1])
            {
               removeFace(cell_1[faceI]);
            }
        }
        cells_.remove(c1);
        lengthScale_.remove(c1);
        if (cellCheck[1] != -1)
        {
            replaceLabel(faceToThrow[1], faceToKeep[1], cells_[cellCheck[1]]);
        }

        // Update the number of cells, and the reverse map
        nCells_--;
        if (c1 < nOldCells_)
        {
            reverseCellMap_[c1] = -1;
        }

        // Check if the cell was added in the current morph, and delete
        if (cellsFromCells_.found(c1))
        {
            cellsFromCells_.erase(c1);
        }
    }

    // Finally remove the face
    removeFace(findex);

    // Return a successful collapse
    return true;
}

// Method for the bisection of an edge in 3D
void dynamicTopoFvMesh::bisectEdge
(
    const label eIndex
)
{
    // Edge bisection performs the following operations:
    //      [1] Add a point at middle of the edge
    //      [2] Bisect all faces surrounding this edge
    //      [3] Bisect all cells surrounding this edge
    //      [4] Create internal/external edges for each bisected face
    //      [5] Create internal faces for each bisected cell
    //      Update faceEdges and edgeFaces information

    face tmpTriFace(3);
    labelList tmpEdgeFaces(3,-1);
    labelList tmpIntEdgeFaces(4,-1);
    labelList tmpFaceEdges(3,-1);
    edge& thisEdge = edges_[eIndex];

#   ifdef FULLDEBUG
    if (debug)
    {
        Info << nl << nl << "Edge: " << eIndex
             << ": " << thisEdge << " is to be bisected. " << endl;
    }
#   endif

    // Obtain maxTetsPerEdge
    label mMax = maxTetsPerEdge();

    // Hull variables
    scalar minQuality;
    DynamicList<label> cellHull(mMax);
    DynamicList<label> faceHull(mMax);
    DynamicList<label> vertexHull(mMax);

    // Lock hull-entities to prevent them from being modified by other threads.
    labelHashSet pLocks, eLocks, fLocks, cLocks;

    // Obtain a ring of vertices around this edge
    if
    (
        constructVertexRing
        (
            eIndex,
            cellHull,
            faceHull,
            vertexHull,
            minQuality,
            pLocks,
            eLocks,
            fLocks,
            cLocks,
            false
        )
    )
    {
        // Put this edge back on the stack and bail out
        edgeStack_[getThreadID(pthread_self())].push(eIndex);
        return;
    }

    // Lock the point mutex
    pMutex_.lock();

    // Add a new point to the end of the list
    label newPointIndex =
        meshPoints_.append
        (
            0.5*
            (
                meshPoints_[thisEdge[0]]
              + meshPoints_[thisEdge[1]]
            )
        );

    // Add an entry to pointEdges as well
    pointEdges_.append(labelList(0));

    // Add an unlocked point mutex
    if (threader_->multiThreaded())
    {
        pointMutex_.append();
    }

    nPoints_++;

    // Unlock the point mutex
    pMutex_.unlock();

    // Lock the edge mutex
    eMutex_.lock();

    // Add a new edge to the end of the list
    label newEdgeIndex =
        insertEdge
        (
            whichEdgePatch(eIndex),
            edge(newPointIndex,thisEdge[1]),
            labelList(faceHull.size(),-1)
        );

    // Remove the existing edge from the pointEdges list
    // of the modified point, and add it to the new point
    sizeDownList(eIndex, pointEdges_[thisEdge[1]]);
    sizeUpList(eIndex, pointEdges_[newPointIndex]);

    // Modify the existing edge
    thisEdge[1] = newPointIndex;

    // Obtain new references
    edge& newEdge = edges_[newEdgeIndex];
    labelList& newEdgeFaces = edgeFaces_[newEdgeIndex];

    // Unlock the edge mutex
    eMutex_.unlock();

    // Keep track of added entities
    labelList addedCellIndices(cellHull.size(),-1);
    labelList addedFaceIndices(faceHull.size(),-1);
    labelList addedEdgeIndices(faceHull.size(),-1);
    labelList addedIntFaceIndices(faceHull.size(),-1);

    // Obtain cellCells for mapping information
    const labelListList& cc = cellCells();

    // Now loop through the hull and bisect individual entities
    forAll(vertexHull, indexI)
    {
        // Fetch the existing face
        face& currFace = faces_[faceHull[indexI]];

        // Modify the existing face
        replaceLabel
        (
            newEdge[1],
            newPointIndex,
            currFace
        );

        // Obtain circular indices
        label nextI = vertexHull.fcIndex(indexI);
        label prevI = vertexHull.rcIndex(indexI);

        // Check if this is an interior/boundary face
        if (cellHull[indexI] != -1)
        {
            cell& currCell = cells_[cellHull[indexI]];

            // Lock the cell mutex
            cMutex_.lock();

            // Create a new cell
            addedCellIndices[indexI] = cells_.append(cell(4));
            cell& newCell = cells_[addedCellIndices[indexI]];
            nCells_++;

            // Add an unlocked cell mutex
            if (threader_->multiThreaded())
            {
                cellMutex_.append();
            }

            // Generate mapping information for this new cell
            label parent;
            labelHashSet masterObjects;

            if (cellHull[indexI] < nOldCells_)
            {
                parent = cellHull[indexI];
            }
            else
            {
                parent = cellParents_[cellHull[indexI]];
            }

            // Insert the parent cell
            cellParents_.insert(addedCellIndices[indexI], parent);

            // Unlock the cell mutex
            cMutex_.unlock();

            // Find the cell's neighbours in the old mesh
            masterObjects.insert(parent);
            forAll(cc[parent], cellI)
            {
                if (!masterObjects.found(cc[parent][cellI]))
                {
                    masterObjects.insert(cc[parent][cellI]);
                }
            }

            // Lock the cell mutex
            cMutex_.lock();

            // Insert mapping info into the HashTable
            cellsFromCells_.insert
            (
                addedCellIndices[indexI],
                objectMap
                (
                    addedCellIndices[indexI],
                    masterObjects.toc()
                )
            );

            // Add a new element to the lengthScale field
            lengthScale_.append(lengthScale_[cellHull[indexI]]);

            // Unlock the cell mutex
            cMutex_.unlock();

            // Configure the interior face
            tmpTriFace[0] = vertexHull[nextI];
            tmpTriFace[1] = vertexHull[indexI];
            tmpTriFace[2] = newPointIndex;

            // Lock the face mutex
            fMutex_.lock();

            // Insert the face
            addedIntFaceIndices[indexI] =
                insertFace
                (
                    -1,
                    tmpTriFace,
                    cellHull[indexI],
                    addedCellIndices[indexI]
                );

            // Add a faceEdges entry as well
            faceEdges_.append(tmpFaceEdges);

            // Unlock the face mutex
            fMutex_.unlock();

            // Add to the new cell
            newCell[0] = addedIntFaceIndices[indexI];

            // Modify the existing face
            forAll(currCell, faceI)
            {
                // Check if this face contains newEdge[1]
                if
                (
                    (currCell[faceI] != faceHull[indexI])
                 && (currCell[faceI] != faceHull[nextI])
                 && (faces_[currCell[faceI]].which(newEdge[1]) > -1)
                )
                {
                    label replaceFace = currCell[faceI];

                    // Check if face reversal is necessary
                    if (owner_[replaceFace] == cellHull[indexI])
                    {
                        if (neighbour_[replaceFace] == -1)
                        {
                            // Change the owner
                            owner_[replaceFace] = addedCellIndices[indexI];
                        }
                        else
                        {
                            // This face has to be reversed
                            faces_[replaceFace] = faces_[replaceFace].reverseFace();
                            owner_[replaceFace] = neighbour_[replaceFace];
                            neighbour_[replaceFace] = addedCellIndices[indexI];
                        }
                    }
                    else
                    {
                        // Keep owner, but change neighbour
                        neighbour_[replaceFace] = addedCellIndices[indexI];
                    }

                    // Look for the edge on the ring
                    labelList& rFaceEdges = faceEdges_[replaceFace];

                    forAll(rFaceEdges, edgeI)
                    {
                        if
                        (
                            edges_[rFaceEdges[edgeI]]
                         == edge(vertexHull[indexI],vertexHull[nextI])
                        )
                        {
                            // Modify edgeFaces to add the
                            // new interior face
                            sizeUpList
                            (
                                addedIntFaceIndices[indexI],
                                edgeFaces_[rFaceEdges[edgeI]]
                            );

                            // Add this edge to faceEdges
                            // for the new interior face
                            faceEdges_[addedIntFaceIndices[indexI]][0] =
                                rFaceEdges[edgeI];
                        }
                    }

                    // Replace face labels
                    replaceLabel
                    (
                        replaceFace,
                        addedIntFaceIndices[indexI],
                        currCell
                    );

                    // Add to the new cell
                    newCell[1] = replaceFace;

                    break;
                }
            }

            // Check if this is a boundary face
            if (cellHull[prevI] == -1)
            {
                // Configure the boundary face
                tmpTriFace[0] = newPointIndex;
                tmpTriFace[1] = newEdge[1];
                tmpTriFace[2] = vertexHull[indexI];

                // Lock the face mutex
                fMutex_.lock();

                // Insert the face
                addedFaceIndices[indexI] =
                    insertFace
                    (
                        whichPatch(faceHull[indexI]),
                        tmpTriFace,
                        addedCellIndices[indexI],
                        -1
                    );

                // Configure edgeFaces
                tmpEdgeFaces[0] = faceHull[indexI];
                tmpEdgeFaces[1] = addedIntFaceIndices[indexI];
                tmpEdgeFaces[2] = addedFaceIndices[indexI];

                // Lock the edge mutex
                eMutex_.lock();

                // Add an edge
                addedEdgeIndices[indexI] =
                    insertEdge
                    (
                        whichPatch(faceHull[indexI]),
                        edge(newPointIndex,vertexHull[indexI]),
                        tmpEdgeFaces
                    );

                // Unlock the edge mutex
                eMutex_.unlock();

                // Add this edge to the interior-face faceEdges entry
                faceEdges_[addedIntFaceIndices[indexI]][1] =
                    addedEdgeIndices[indexI];

                // Configure faceEdges for this boundary face
                tmpFaceEdges[0] = addedEdgeIndices[indexI];
                tmpFaceEdges[1] = newEdgeIndex;

                // Find the third edge
                labelList& rFaceEdges = faceEdges_[faceHull[indexI]];

                forAll(rFaceEdges, edgeI)
                {
                    if
                    (
                        edges_[rFaceEdges[edgeI]]
                     == edge(newEdge[1],vertexHull[indexI])
                    )
                    {
                        // Add this edge to faceEdges for this face
                        tmpFaceEdges[2] = rFaceEdges[edgeI];

                        // Modify faceEdges
                        replaceLabel
                        (
                            rFaceEdges[edgeI],
                            addedEdgeIndices[indexI],
                            rFaceEdges
                        );

                        // Modify edgeFaces
                        replaceLabel
                        (
                            faceHull[indexI],
                            addedFaceIndices[indexI],
                            edgeFaces_[tmpFaceEdges[2]]
                        );

                        break;
                    }
                }

                // Add the faceEdges entry
                faceEdges_.append(tmpFaceEdges);

                // Unlock the face mutex
                fMutex_.unlock();

                // Add an entry to newEdgeFaces
                newEdgeFaces[indexI] = addedFaceIndices[indexI];

                // Add an entry for this cell
                newCell[2] = addedFaceIndices[indexI];
            }
            else
            // Check if a cell was added before this
            if (addedCellIndices[prevI] != -1)
            {
                // Configure the interior face
                tmpTriFace[0] = vertexHull[indexI];
                tmpTriFace[1] = newEdge[1];
                tmpTriFace[2] = newPointIndex;

                // Lock the face mutex
                fMutex_.lock();

                // Insert the face
                addedFaceIndices[indexI] =
                    insertFace
                    (
                        -1,
                        tmpTriFace,
                        addedCellIndices[prevI],
                        addedCellIndices[indexI]
                    );

                // Configure edgeFaces
                tmpIntEdgeFaces[0] = faceHull[indexI];
                tmpIntEdgeFaces[1] = addedIntFaceIndices[indexI];
                tmpIntEdgeFaces[2] = addedFaceIndices[indexI];
                tmpIntEdgeFaces[3] = addedIntFaceIndices[prevI];

                // Lock the edge mutex
                eMutex_.lock();

                // Add an internal edge
                addedEdgeIndices[indexI] =
                    insertEdge
                    (
                        -1,
                        edge(newPointIndex,vertexHull[indexI]),
                        tmpIntEdgeFaces
                    );

                // Unlock the edge mutex
                eMutex_.unlock();

                // Add this edge to the interior-face faceEdges entry..
                faceEdges_[addedIntFaceIndices[indexI]][1] =
                    addedEdgeIndices[indexI];

                // ... and to the previous interior face as well
                faceEdges_[addedIntFaceIndices[prevI]][2] =
                    addedEdgeIndices[indexI];

                // Configure faceEdges for this split interior face
                tmpFaceEdges[0] = addedEdgeIndices[indexI];
                tmpFaceEdges[1] = newEdgeIndex;

                // Find the third edge
                labelList& rFaceEdges = faceEdges_[faceHull[indexI]];

                forAll(rFaceEdges, edgeI)
                {
                    if
                    (
                        edges_[rFaceEdges[edgeI]]
                     == edge(newEdge[1],vertexHull[indexI])
                    )
                    {
                        // Add this edge to faceEdges for this face
                        tmpFaceEdges[2] = rFaceEdges[edgeI];

                        // Modify faceEdges
                        replaceLabel
                        (
                            rFaceEdges[edgeI],
                            addedEdgeIndices[indexI],
                            rFaceEdges
                        );

                        // Modify edgeFaces
                        replaceLabel
                        (
                            faceHull[indexI],
                            addedFaceIndices[indexI],
                            edgeFaces_[tmpFaceEdges[2]]
                        );

                        break;
                    }
                }

                // Add the faceEdges entry
                faceEdges_.append(tmpFaceEdges);

                // Unlock the face mutex
                fMutex_.unlock();

                // Add an entry to newEdgeFaces
                newEdgeFaces[indexI] = addedFaceIndices[indexI];

                // Add an entry for this cell
                newCell[2] = addedFaceIndices[indexI];

                // Make the final entry for the previous cell
                cells_[addedCellIndices[prevI]][3] = addedFaceIndices[indexI];
            }

            // Do the first interior face at the end
            if (indexI == vertexHull.size() - 1)
            {
                // Configure the interior face
                tmpTriFace[0] = newPointIndex;
                tmpTriFace[1] = newEdge[1];
                tmpTriFace[2] = vertexHull[0];

                // Lock the face mutex
                fMutex_.lock();

                // Insert the face
                addedFaceIndices[0] =
                    insertFace
                    (
                        -1,
                        tmpTriFace,
                        addedCellIndices[0],
                        addedCellIndices[indexI]
                    );

                // Configure edgeFaces
                tmpIntEdgeFaces[0] = faceHull[0];
                tmpIntEdgeFaces[1] = addedIntFaceIndices[0];
                tmpIntEdgeFaces[2] = addedFaceIndices[0];
                tmpIntEdgeFaces[3] = addedIntFaceIndices[indexI];

                // Lock the edge mutex
                eMutex_.lock();

                // Add an internal edge
                addedEdgeIndices[0] =
                    insertEdge
                    (
                        -1,
                        edge(newPointIndex,vertexHull[0]),
                        tmpIntEdgeFaces
                    );

                // Unlock the edge mutex
                eMutex_.unlock();

                // Add this edge to the interior-face faceEdges entry..
                faceEdges_[addedIntFaceIndices[0]][1] =
                    addedEdgeIndices[0];

                // ... and to the previous interior face as well
                faceEdges_[addedIntFaceIndices[indexI]][2] =
                    addedEdgeIndices[0];

                // Configure faceEdges for the first split face
                tmpFaceEdges[0] = addedEdgeIndices[0];
                tmpFaceEdges[1] = newEdgeIndex;

                // Find the third edge
                labelList& rFaceEdges = faceEdges_[faceHull[0]];

                forAll(rFaceEdges, edgeI)
                {
                    if
                    (
                        edges_[rFaceEdges[edgeI]]
                     == edge(newEdge[1],vertexHull[0])
                    )
                    {
                        // Add this edge to faceEdges for this face
                        tmpFaceEdges[2] = rFaceEdges[edgeI];

                        // Modify faceEdges
                        replaceLabel
                        (
                            rFaceEdges[edgeI],
                            addedEdgeIndices[0],
                            rFaceEdges
                        );

                        // Modify edgeFaces
                        replaceLabel
                        (
                            faceHull[0],
                            addedFaceIndices[0],
                            edgeFaces_[tmpFaceEdges[2]]
                        );

                        break;
                    }
                }

                // Add the faceEdges entry
                faceEdges_.append(tmpFaceEdges);

                // Unlock the face mutex
                fMutex_.unlock();

                // Add an entry to newEdgeFaces
                newEdgeFaces[0] = addedFaceIndices[0];

                // Add an entry for this cell
                newCell[3] = addedFaceIndices[0];

                // Make the final entry for the first cell
                cells_[addedCellIndices[0]][2] = addedFaceIndices[0];
            }
        }
        else
        {
            // Configure the final boundary face
            tmpTriFace[0] = vertexHull[indexI];
            tmpTriFace[1] = newEdge[1];
            tmpTriFace[2] = newPointIndex;

            // Lock the face mutex
            fMutex_.lock();

            // Insert the face
            addedFaceIndices[indexI] =
                insertFace
                (
                    whichPatch(faceHull[indexI]),
                    tmpTriFace,
                    addedCellIndices[prevI],
                    -1
                );

            // Configure edgeFaces
            tmpEdgeFaces[0] = addedFaceIndices[indexI];
            tmpEdgeFaces[1] = addedIntFaceIndices[prevI];
            tmpEdgeFaces[2] = faceHull[indexI];

            // Lock the edge mutex
            eMutex_.lock();

            // Add an edge
            addedEdgeIndices[indexI] =
                insertEdge
                (
                    whichPatch(faceHull[indexI]),
                    edge(newPointIndex,vertexHull[indexI]),
                    tmpEdgeFaces
                );

            // Unlock the edge mutex
            eMutex_.unlock();

            // Add a faceEdges entry to the previous interior face
            faceEdges_[addedIntFaceIndices[prevI]][2] =
                addedEdgeIndices[indexI];

            // Configure faceEdges for the final boundary face
            tmpFaceEdges[0] = addedEdgeIndices[indexI];
            tmpFaceEdges[1] = newEdgeIndex;

            // Find the third edge
            labelList& rFaceEdges = faceEdges_[faceHull[indexI]];

            forAll(rFaceEdges, edgeI)
            {
                if
                (
                    edges_[rFaceEdges[edgeI]]
                 == edge(newEdge[1],vertexHull[indexI])
                )
                {
                    // Add this edge to faceEdges for this face
                    tmpFaceEdges[2] = rFaceEdges[edgeI];

                    // Modify faceEdges
                    replaceLabel
                    (
                        rFaceEdges[edgeI],
                        addedEdgeIndices[indexI],
                        rFaceEdges
                    );

                    // Modify edgeFaces
                    replaceLabel
                    (
                        faceHull[indexI],
                        addedFaceIndices[indexI],
                        edgeFaces_[tmpFaceEdges[2]]
                    );

                    break;
                }
            }

            // Add the faceEdges entry
            faceEdges_.append(tmpFaceEdges);

            // Unlock the face mutex
            fMutex_.unlock();

            // Add an entry to newEdgeFaces
            newEdgeFaces[indexI] = addedFaceIndices[indexI];

            // Make the final entry for the previous cell
            cells_[addedCellIndices[prevI]][3] = addedFaceIndices[indexI];
        }
    }

    // Unlock all entities
    unlockMutexLists(pLocks, eLocks, fLocks, cLocks);

#   ifdef FULLDEBUG
    if (debug)
    {
        Info << "Modified cells: " << endl;
        forAll(cellHull, cellI)
        {
            if (cellHull[cellI] != -1)
            {
                Info << cellHull[cellI] << ":: "
                     << cells_[cellHull[cellI]]
                     << endl;
            }
        }

        Info << "Added cells: " << endl;
        forAll(addedCellIndices, cellI)
        {
            if (addedCellIndices[cellI] != -1)
            {
                Info << addedCellIndices[cellI] << ":: "
                     << cells_[addedCellIndices[cellI]]
                     << endl;
            }
        }

        Info << "Modified faces: " << endl;
        forAll(faceHull, faceI)
        {
            Info << faceHull[faceI] << ":: "
                 << faces_[faceHull[faceI]] << ": "
                 << owner_[faceHull[faceI]] << ": "
                 << neighbour_[faceHull[faceI]] << " "
                 << "faceEdges:: " << faceEdges_[faceHull[faceI]]
                 << endl;
        }

        Info << "Added faces: " << endl;
        forAll(addedFaceIndices, faceI)
        {
            Info << addedFaceIndices[faceI] << ":: "
                 << faces_[addedFaceIndices[faceI]] << ": "
                 << owner_[addedFaceIndices[faceI]] << ": "
                 << neighbour_[addedFaceIndices[faceI]] << " "
                 << "faceEdges:: " << faceEdges_[addedFaceIndices[faceI]]
                 << endl;
        }
        forAll(addedIntFaceIndices, faceI)
        {
            if (addedIntFaceIndices[faceI] != -1)
            {
                Info << addedIntFaceIndices[faceI] << ":: "
                     << faces_[addedIntFaceIndices[faceI]] << ": "
                     << owner_[addedIntFaceIndices[faceI]] << ": "
                     << neighbour_[addedIntFaceIndices[faceI]] << " "
                     << "faceEdges:: "
                     << faceEdges_[addedIntFaceIndices[faceI]]
                     << endl;
            }
        }

        Info << "New edge:: " << newEdgeIndex
             << ": " << edges_[newEdgeIndex]
             << " edgeFaces:: " << newEdgeFaces << endl;

        Info << "Added edges: " << endl;
        forAll(addedEdgeIndices, edgeI)
        {
            Info << addedEdgeIndices[edgeI]
                 << ":: " << edges_[addedEdgeIndices[edgeI]]
                 << " edgeFaces:: " << edgeFaces_[addedEdgeIndices[edgeI]]
                 << endl;
        }

        Info << "New Point:: " << newPointIndex << endl;
        Info << "pointEdges:: " << pointEdges_[newPointIndex] << endl;
    }
#   endif
}

// Method for the collapse of an edge in 3D
// Returns a boolean value indicating whether the collapse was valid
bool dynamicTopoFvMesh::collapseEdge
(
    const label eIndex
)
{
    // Edge collapse performs the following operations:
    //      [1] Checks if either vertex of the edge is on a boundary
    //      [2] Checks whether cells attached to deleted vertices will be valid
    //          after the edge-collapse operation
    //      [3] Deletes all cells surrounding this edge
    //      [4] Deletes all faces surrounding this edge
    //      [5] Deletes all faces surrounding the deleted vertex attached
    //          to the cells in [3]
    //      [6] Checks the orientation of faces connected to the retained
    //          vertices
    //      [7] Remove one of the vertices of the edge
    //      Update faceEdges and edgeFaces information

    bool found = false;
    label replaceIndex = -1;
    edge& thisEdge = edges_[eIndex];
    FixedList<bool,2> edgeBoundary(false);

#   ifdef FULLDEBUG
    if (debug)
    {
        Info << nl << nl << "Edge: " << eIndex
             << ": " << thisEdge << " is to be collapsed. " << endl;
    }
#   endif

    // Obtain maxTetsPerEdge
    label mMax = maxTetsPerEdge();

    // Hull variables
    scalar minQuality;
    DynamicList<label> cellHull(mMax);
    DynamicList<label> faceHull(mMax);
    DynamicList<label> vertexHull(mMax);

    // Lock hull-entities to prevent them from being modified by other threads.
    labelHashSet pLocks, eLocks, fLocks, cLocks;

    // Obtain a ring of vertices around this edge
    if
    (
        constructVertexRing
        (
            eIndex,
            cellHull,
            faceHull,
            vertexHull,
            minQuality,
            pLocks,
            eLocks,
            fLocks,
            cLocks,
            false
        )
    )
    {
        // Put this edge back on the stack and bail out
        edgeStack_[getThreadID(pthread_self())].push(eIndex);
        return false;
    }

    // Determine ring edges and the hull faces connected to them
    labelList edgeHull(vertexHull.size(), -1);
    labelListList hullEdgesAndFaces(4, labelList(faceHull.size(), -1));

    constructEdgeRing
    (
        eIndex,
        thisEdge,
        vertexHull,
        faceHull,
        cellHull,
        edgeHull,
        hullEdgesAndFaces,
        edgeBoundary
    );

    // Configure the new point-position
    point newPoint = vector::zero;

    // Decide which point to remove
    FixedList<label,2> checkPoints(-1);
    label collapsePoint = -1, replacePoint = -1;
    label removeEdgeIndex = -1, removeFaceIndex = -1;
    label replaceEdgeIndex = -1, replaceFaceIndex = -1;

    if (edgeBoundary[0] && !edgeBoundary[1])
    {
        // Collapse to the first node
        replacePoint = thisEdge[0];
        collapsePoint = thisEdge[1];
        replaceEdgeIndex = 0;
        replaceFaceIndex = 1;
        removeEdgeIndex = 2;
        removeFaceIndex = 3;
        checkPoints[0] = collapsePoint;
        newPoint = meshPoints_[thisEdge[0]];
    }
    else
    if (!edgeBoundary[0] && edgeBoundary[1])
    {
        // Collapse to the second node
        replacePoint = thisEdge[1];
        collapsePoint = thisEdge[0];
        removeEdgeIndex = 0;
        removeFaceIndex = 1;
        replaceEdgeIndex = 2;
        replaceFaceIndex = 3;
        checkPoints[0] = collapsePoint;
        newPoint = meshPoints_[thisEdge[1]];
    }
    else
    {
        // Collapse to the second node by default
        replacePoint = thisEdge[1];
        collapsePoint = thisEdge[0];
        removeEdgeIndex = 0;
        removeFaceIndex = 1;
        replaceEdgeIndex = 2;
        replaceFaceIndex = 3;
        checkPoints[0] = collapsePoint;
        checkPoints[1] = replacePoint;
        // Position collapse to the mid-point
        newPoint = 0.5*(meshPoints_[thisEdge[0]] + meshPoints_[thisEdge[1]]);
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        Info << "Vertices: " << vertexHull << endl;
        Info << "Edges: " << edgeHull << endl;
        Info << "Faces: " << faceHull << endl;
        Info << "Cells: " << cellHull << endl;
        Info << "replacePoint: " << replacePoint << endl;
        Info << "collapsePoint: " << collapsePoint << endl;
        Info << "hullEdgesAndFaces (removed faces): " << endl;

        forAll(hullEdgesAndFaces[removeFaceIndex], faceI)
        {
            label fIndex = hullEdgesAndFaces[removeFaceIndex][faceI];

            if (fIndex != -1)
            {
                Info << fIndex << ": " << faces_[fIndex] << endl;
            }
        }

        Info << "hullEdgesAndFaces (removed edges): " << endl;
        forAll(hullEdgesAndFaces[removeEdgeIndex], edgeI)
        {
            label ieIndex = hullEdgesAndFaces[removeEdgeIndex][edgeI];

            if (ieIndex != -1)
            {
                Info << ieIndex << ": " << edges_[ieIndex] << endl;
            }
        }

        Info << "hullEdgesAndFaces (replaced faces): " << endl;
        forAll(hullEdgesAndFaces[replaceFaceIndex], faceI)
        {
            label fIndex = hullEdgesAndFaces[replaceFaceIndex][faceI];

            if (fIndex != -1)
            {
                Info << fIndex << ": " << faces_[fIndex] << endl;
            }
        }

        Info << "hullEdgesAndFaces (replaced edges): " << endl;
        forAll(hullEdgesAndFaces[replaceEdgeIndex], edgeI)
        {
            label ieIndex = hullEdgesAndFaces[replaceEdgeIndex][edgeI];

            if (ieIndex != -1)
            {
                Info << ieIndex << ": " << edges_[ieIndex] << endl;
            }
        }

        labelList& collapsePointEdges = pointEdges_[collapsePoint];
        Info << "pointEdges (collapsePoint): ";
        forAll(collapsePointEdges, edgeI)
        {
            Info << collapsePointEdges[edgeI] << " ";
        }
        Info << endl;
    }
#   endif

    // Loop through edges and check for feasibility of collapse
    labelHashSet cellsChecked;

    // Add all hull cells as 'checked'
    forAll(cellHull, cellI)
    {
        if (cellHull[cellI] != -1)
        {
            cellsChecked.insert(cellHull[cellI]);
        }
    }

    // Check collapsibility of cells around edges with the re-configured point
    forAll(checkPoints, pointI)
    {
        if (checkPoints[pointI] == -1)
        {
            continue;
        }

        labelList& checkPointEdges = pointEdges_[checkPoints[pointI]];

        forAll(checkPointEdges, edgeI)
        {
            labelList& eFaces = edgeFaces_[checkPointEdges[edgeI]];

            // Build a list of cells to check
            forAll(eFaces, faceI)
            {
                label own = owner_[eFaces[faceI]];
                label nei = neighbour_[eFaces[faceI]];

                // Check owner cell
                if (!cellsChecked.found(own))
                {
                    // Check if a collapse is feasible
                    if
                    (
                        checkCollapse
                        (
                            newPoint,
                            checkPoints[pointI],
                            own,
                            edgeBoundary,
                            cellsChecked
                        )
                    )
                    {
                        // Unlock all entities
                        unlockMutexLists(pLocks, eLocks, fLocks, cLocks);

                        return false;
                    }
                }

                // Check neighbour cell
                if (!cellsChecked.found(nei) && nei != -1)
                {
                    // Check if a collapse is feasible
                    if
                    (
                        checkCollapse
                        (
                            newPoint,
                            checkPoints[pointI],
                            nei,
                            edgeBoundary,
                            cellsChecked
                        )
                    )
                    {
                        // Unlock all entities
                        unlockMutexLists(pLocks, eLocks, fLocks, cLocks);

                        return false;
                    }
                }
            }
        }
    }

    // Renumber all hull faces and edges
    forAll(faceHull, indexI)
    {
        // Loop through all faces of the edge to be removed
        // and reassign them to the replacement edge
        label edgeToRemove = hullEdgesAndFaces[removeEdgeIndex][indexI];
        label faceToRemove = hullEdgesAndFaces[removeFaceIndex][indexI];
        label cellToRemove = cellHull[indexI];
        label replaceEdge = hullEdgesAndFaces[replaceEdgeIndex][indexI];
        label replaceFace = hullEdgesAndFaces[replaceFaceIndex][indexI];

        labelList& rmvEdgeFaces = edgeFaces_[edgeToRemove];
        labelList& rplEdgeFaces = edgeFaces_[replaceEdge];

        // Replace edge labels
        forAll(rmvEdgeFaces, faceI)
        {
            replaceLabel
            (
                edgeToRemove,
                replaceEdge,
                faceEdges_[rmvEdgeFaces[faceI]]
            );

            // Loop through faces associated with this edge,
            // and renumber them as well.
            face& faceToCheck = faces_[rmvEdgeFaces[faceI]];
            if ((replaceIndex = faceToCheck.which(collapsePoint)) > -1)
            {
#               ifdef FULLDEBUG
                if (debug)
                {
                    Info << "Renumbering face: "
                         << rmvEdgeFaces[faceI] << ": "
                         << faceToCheck << endl;
                }
#               endif

                faceToCheck[replaceIndex] = replacePoint;
            }
        }

        // Size up existing edgeFaces
        forAll(rmvEdgeFaces, faceI)
        {
            // Avoid hull faces
            if (rmvEdgeFaces[faceI] == faceHull[indexI])
            {
                sizeDownList
                (
                    faceHull[indexI],
                    rplEdgeFaces
                );

                continue;
            }

            found = false;

            forAll(hullEdgesAndFaces[removeFaceIndex], faceII)
            {
                if
                (
                    rmvEdgeFaces[faceI]
                 == hullEdgesAndFaces[removeFaceIndex][faceII]
                )
                {
                    found = true;
                    break;
                }
            }

            // Size-up the list if the face hasn't been found
            if (!found)
            {
                sizeUpList
                (
                    rmvEdgeFaces[faceI],
                    rplEdgeFaces
                );
            }
        }

        if (cellToRemove != -1)
        {
            // Size down edgeFaces for the ring edges
            sizeDownList
            (
                faceToRemove,
                edgeFaces_[edgeHull[indexI]]
            );

            // Ensure proper orientation of retained faces
            if (owner_[faceToRemove] == cellToRemove)
            {
                if (owner_[replaceFace] == cellToRemove)
                {
                    if
                    (
                        (neighbour_[faceToRemove] > neighbour_[replaceFace])
                     && (neighbour_[replaceFace] != -1)
                    )
                    {
                        // This face is to be flipped
                        faces_[replaceFace] = faces_[replaceFace].reverseFace();
                        owner_[replaceFace] = neighbour_[replaceFace];
                        neighbour_[replaceFace] = neighbour_[faceToRemove];
                    }
                    else
                    {
                        // Keep orientation intact, and update the owner
                        owner_[replaceFace] = neighbour_[faceToRemove];
                    }
                }
                else
                {
                    // Keep orientation intact, and update the neighbour
                    neighbour_[replaceFace] = neighbour_[faceToRemove];
                }

                // Update the cell
                if (neighbour_[faceToRemove] != -1)
                {
                    replaceLabel
                    (
                        faceToRemove,
                        replaceFace,
                        cells_[neighbour_[faceToRemove]]
                    );
                }
            }
            else
            {
                if (neighbour_[replaceFace] == cellToRemove)
                {
                    if (owner_[faceToRemove] < owner_[replaceFace])
                    {
                        // This face is to be flipped
                        faces_[replaceFace] = faces_[replaceFace].reverseFace();
                        neighbour_[replaceFace] = owner_[replaceFace];
                        owner_[replaceFace] = owner_[faceToRemove];
                    }
                    else
                    {
                        // Keep orientation intact, and update the neighbour
                        neighbour_[replaceFace] = owner_[faceToRemove];
                    }
                }
                else
                {
                    // Keep orientation intact, and update the owner
                    owner_[replaceFace] = owner_[faceToRemove];
                }

                // Update the cell
                if (owner_[faceToRemove] != -1)
                {
                    replaceLabel
                    (
                        faceToRemove,
                        replaceFace,
                        cells_[owner_[faceToRemove]]
                    );
                }
            }

            // Check orientation of faces
            if (owner_[replaceFace] == -1)
            {
                FatalErrorIn("dynamicTopoFvMesh::collapseEdge()")
                    << "Face: " << replaceFace << ": " << faces_[replaceFace]
                    << " is an orphan, i.e, no owner cell."
                    << abort(FatalError);
            }
            else
            if
            (
                (neighbour_[replaceFace] == -1)
             && (whichPatch(replaceFace) < 0)
            )
            {
                FatalErrorIn("dynamicTopoFvMesh::collapseEdge()")
                    << "Face: " << replaceFace << faces_[replaceFace]
                    << " is being converted to a boundary."
                    << abort(FatalError);
            }
        }
    }

    // Lock mutexes
    eMutex_.lock();
    fMutex_.lock();
    cMutex_.lock();

    // Remove all hull entities
    forAll(faceHull, indexI)
    {
        label edgeToRemove = hullEdgesAndFaces[removeEdgeIndex][indexI];
        label faceToRemove = hullEdgesAndFaces[removeFaceIndex][indexI];
        label cellToRemove = cellHull[indexI];

        if (cellToRemove != -1)
        {
            // Remove faceToRemove and associated faceEdges
            removeFace(faceToRemove);
            faceEdges_.remove(faceToRemove);

            // Remove from list of locked faces
            if (fLocks.found(faceToRemove))
            {
                fLocks.erase(faceToRemove);
            }

            // Remove the hull cell
            cells_.remove(cellToRemove);
            lengthScale_.remove(cellToRemove);

            // Remove from list of locked cells
            if (cLocks.found(cellToRemove))
            {
                cLocks.erase(cellToRemove);
            }

            // Remove the cell mutex
            if (threader_->multiThreaded())
            {
                // Unlock it first
                cellMutex_[cellToRemove].unlock();
                cellMutex_.remove(cellToRemove);
            }

            // Update the number of cells, and the reverse cell map
            nCells_--;

            if (cellToRemove < nOldCells_)
            {
                reverseCellMap_[cellToRemove] = -1;
            }

            // Check if the cell was added in the current morph, and delete
            if (cellsFromCells_.found(cellToRemove))
            {
                cellsFromCells_.erase(cellToRemove);
            }
        }

        // Remove the hull edge and associated edgeFaces
        removeEdge(edgeToRemove);

        // Remove from list of locked edges
        if (eLocks.found(edgeToRemove))
        {
            eLocks.erase(edgeToRemove);
        }

        // Remove the hull face and associated faceEdges
        removeFace(faceHull[indexI]);
        faceEdges_.remove(faceHull[indexI]);

        // Remove from list of locked faces
        if (fLocks.found(faceHull[indexI]))
        {
            fLocks.erase(faceHull[indexI]);
        }
    }

    // Unlock mutexes
    cMutex_.unlock();
    fMutex_.unlock();
    eMutex_.unlock();

    // Loop through pointEdges for the collapsePoint,
    // and replace all occurrences with replacePoint.
    // Size-up pointEdges for the replacePoint as well.
    labelList& pEdges = pointEdges_[collapsePoint];

    forAll(pEdges, edgeI)
    {
        // Renumber edges
        edge& edgeToCheck = edges_[pEdges[edgeI]];
        labelList& eFaces = edgeFaces_[pEdges[edgeI]];

        if (pEdges[edgeI] != eIndex)
        {
#           ifdef FULLDEBUG
            if (debug)
            {
                Info << "Renumbering [edge]: "
                     << pEdges[edgeI] << ": "
                     << edgeToCheck << endl;
            }
#           endif

            if (edgeToCheck[0] == collapsePoint)
            {
                edgeToCheck[0] = replacePoint;

                sizeUpList
                (
                    pEdges[edgeI],
                    pointEdges_[replacePoint]
                );
            }
            else
            if (edgeToCheck[1] == collapsePoint)
            {
                edgeToCheck[1] = replacePoint;

                sizeUpList
                (
                    pEdges[edgeI],
                    pointEdges_[replacePoint]
                );
            }
            else
            {
                // Looks like pointEdges is inconsistent
                FatalErrorIn("dynamicTopoFvMesh::collapseEdge()") << nl
                    << "pointEdges is inconsistent." << nl
                    << "Point: " << collapsePoint << nl
                    << "pointEdges: " << pEdges << nl
                    << abort(FatalError);
            }

            // Loop through faces associated with this edge,
            // and renumber them as well.
            forAll(eFaces, faceI)
            {
                face& faceToCheck = faces_[eFaces[faceI]];

                if ((replaceIndex = faceToCheck.which(collapsePoint)) > -1)
                {
#                   ifdef FULLDEBUG
                    if (debug)
                    {
                        Info << "Renumbering face: "
                             << eFaces[faceI] << ": "
                             << faceToCheck << endl;
                    }
#                   endif

                    faceToCheck[replaceIndex] = replacePoint;
                }
            }
        }
    }

    // Lock the point mutex
    pMutex_.lock();

    // Move to the new point
    meshPoints_[replacePoint] = newPoint;

    // Remove the collapse point
    meshPoints_.remove(collapsePoint);
    nPoints_--;

    // Remove from list of locked points
    if (pLocks.found(collapsePoint))
    {
        pLocks.erase(collapsePoint);
    }

    // Remove the point mutex
    if (threader_->multiThreaded())
    {
        // Unlock it first
        pointMutex_[collapsePoint].unlock();
        pointMutex_.remove(collapsePoint);
    }

    // Null pointEdges so that removeEdge deletes it.
    pointEdges_[collapsePoint] = labelList(0);

    // Unlock the point mutex
    pMutex_.unlock();

    // Update the reverse point map
    if (collapsePoint < nOldPoints_)
    {
        reversePointMap_[collapsePoint] = -1;
    }

    // Lock the edge mutex
    eMutex_.lock();

    // Remove the edge
    removeEdge(eIndex);

    // Remove from list of locked edges
    if (eLocks.found(eIndex))
    {
        eLocks.erase(eIndex);
    }

    // Unlock the edge mutex
    eMutex_.unlock();

    // Unlock all entities
    unlockMutexLists(pLocks, eLocks, fLocks, cLocks);

    // Return a successful collapse
    return true;
}

// Utility to determine ring edges and the hull edges/faces connected to them.
void dynamicTopoFvMesh::constructEdgeRing
(
    const label eIndex,
    const edge& edgeToCheck,
    const labelList& hullVertices,
    const labelList& hullFaces,
    const labelList& hullCells,
    labelList& ringEdges,
    labelListList& hullEdgesAndFaces,
    FixedList<bool,2>& edgeBoundary
)
{
    forAll(hullVertices, indexI)
    {
        // Obtain circular indices
        label nextI = hullVertices.fcIndex(indexI);

        // Obtain edges connected to top and bottom
        // vertices of edgeToCheck
        labelList& fEdges = faceEdges_[hullFaces[indexI]];
        forAll(fEdges, edgeI)
        {
            if
            (
                edges_[fEdges[edgeI]]
             == edge(edgeToCheck[0],hullVertices[indexI])
            )
            {
                hullEdgesAndFaces[0][indexI] = fEdges[edgeI];
            }

            if
            (
                edges_[fEdges[edgeI]]
             == edge(edgeToCheck[1],hullVertices[indexI])
            )
            {
                hullEdgesAndFaces[2][indexI] = fEdges[edgeI];
            }
        }

        if (hullCells[indexI] != -1)
        {
            cell& currCell = cells_[hullCells[indexI]];

            // Look for the ring-edge
            forAll(currCell, faceI)
            {
                // Check if this face contains edgeToCheck[0]
                if
                (
                    (currCell[faceI] != hullFaces[indexI])
                 && (currCell[faceI] != hullFaces[nextI])
                 && (faces_[currCell[faceI]].which(edgeToCheck[0]) > -1)
                )
                {
                    hullEdgesAndFaces[1][indexI] = currCell[faceI];

                    // Look for the edge on the ring
                    labelList& rFaceEdges = faceEdges_[currCell[faceI]];

                    bool found = false;

                    forAll(rFaceEdges, edgeI)
                    {
                        if
                        (
                            edges_[rFaceEdges[edgeI]]
                         == edge(hullVertices[indexI],hullVertices[nextI])
                        )
                        {
                            ringEdges[indexI] = rFaceEdges[edgeI];

                            found = true;
                            break;
                        }
                    }

                    if (!found)
                    {
                        // Looks like faceEdges is inconsistent
                        FatalErrorIn("dynamicTopoFvMesh::constructEdgeRing()")
                            << nl << "Unable to find ring edge: " << nl
                            << edge(hullVertices[indexI],hullVertices[nextI])
                            << " in face: " << currCell[faceI]
                            << ": " << faces_[currCell[faceI]] << nl
                            << " faceEdges: " << rFaceEdges << nl
                            << abort(FatalError);
                    }
                }

                // Check if this face contains edgeToCheck[1]
                if
                (
                    (currCell[faceI] != hullFaces[indexI])
                 && (currCell[faceI] != hullFaces[nextI])
                 && (faces_[currCell[faceI]].which(edgeToCheck[1]) > -1)
                )
                {
                    hullEdgesAndFaces[3][indexI] = currCell[faceI];
                }
            }
        }
    }

    // Loop through edges connected to both points and check if any of them
    // lie on boundaries
    forAll(edgeToCheck, pointI)
    {
        labelList& pEdges = pointEdges_[edgeToCheck[pointI]];

        forAll(pEdges, edgeI)
        {
            // Determine the patch this edge belongs to
            if (whichEdgePatch(pEdges[edgeI]) > -1)
            {
                edgeBoundary[pointI] = true;
                break;
            }
        }
    }

    // Check if either point lies on a bounding curve
    if (edgeBoundary[0] && edgeBoundary[1])
    {
        edgeBoundary = false;

        forAll(edgeToCheck, pointI)
        {
            labelList& pEdges = pointEdges_[edgeToCheck[pointI]];

            forAll(pEdges, edgeI)
            {
                if (checkBoundingCurve(pEdges[edgeI]))
                {
                    edgeBoundary[pointI] = true;
                    break;
                }
            }
        }
    }
}

// Utility method to check whether the cell given by 'cellIndex' will yield
// a valid cell when 'pointIndex' is moved to 'newPoint'. The routine performs
// volume-based checks. Returns 'true' if the collapse in NOT feasible, and
// makes entries in cellsChecked to avoid repetitive checks.
bool dynamicTopoFvMesh::checkCollapse
(
    const point& newPoint,
    const label pointIndex,
    const label cellIndex,
    const FixedList<bool,2>& edgeBoundary,
    labelHashSet& cellsChecked
)
{
    label faceIndex = -1;
    scalar cellVolume = 0.0;
    cell& cellToCheck = cells_[cellIndex];

    // Look for a face that doesn't contain 'pointIndex'
    forAll(cellToCheck, faceI)
    {
        face& currFace = faces_[cellToCheck[faceI]];

        if (currFace.which(pointIndex) < 0)
        {
            faceIndex = cellToCheck[faceI];
            break;
        }
    }

    // Compute cell-volume
    face& faceToCheck = faces_[faceIndex];

    if (owner_[faceIndex] == cellIndex)
    {
        cellVolume = tetVolume
        (
            meshPoints_[faceToCheck[2]],
            meshPoints_[faceToCheck[1]],
            meshPoints_[faceToCheck[0]],
            newPoint
        );
    }
    else
    {
        cellVolume = tetVolume
        (
            meshPoints_[faceToCheck[0]],
            meshPoints_[faceToCheck[1]],
            meshPoints_[faceToCheck[2]],
            newPoint
        );
    }

    // Final cell-volume check
    if (cellVolume < VSMALL)
    {
#       ifdef FULLDEBUG
        if (debug)
        {
            WarningIn
            (
                "dynamicTopoFvMesh::checkCollapse"
            )   << "\nCollapsing cell: " << cellIndex << " containing points:\n"
                << faceToCheck[0] << "," << faceToCheck[1] << ","
                << faceToCheck[2] << "," << pointIndex << nl
                << "will yield a negative volume: " << cellVolume
                << ", when " << pointIndex << " is moved to location: " << nl
                << newPoint << endl;
        }
#       endif
        return true;
    }

    // No problems, so a collapse is feasible
    cellsChecked.insert(cellIndex);
    return false;
}

// Check if the boundary face is adjacent to a sliver-cell,
// and remove it by a two-step bisection/collapse operation.
bool dynamicTopoFvMesh::remove2DSliver
(
    const label findex
)
{
    FixedList<label,2> c0BdyIndex, c0IntIndex;
    FixedList<face,2>  c0BdyFace,  c0IntFace;

    // Measure the boundary edge-length of the face in question
    edge& checkEdge = edgeToWatch_[findex];
    point& a = meshPoints_[checkEdge[0]];
    point& b = meshPoints_[checkEdge[1]];
    scalar length = mag(b-a);

    // Determine the neighbouring cell
    label c0 = owner_[findex];

    // Find the prism-faces
    findPrismFaces
    (
        findex,
        c0,
        c0BdyFace,
        c0BdyIndex,
        c0IntFace,
        c0IntIndex
    );

    // Determine the boundary triangular face area
    scalar area = triFaceArea(c0BdyFace[0]);

    // This cell has to be removed...
    if (mag(area) < (sliverThreshold_*length*length))
    {
        // Step 1: Bisect the boundary quad face
        bisectInteriorFace_ = -1;
        bisectQuadFace(findex);

        // Step 2: Collapse the newly created internal quad face
        bool success = collapseQuadFace(bisectInteriorFace_);

        if (!success)
        {
            WarningIn
            (
                "dynamicTopoFvMesh::remove2DSliver(const label, face&)"
            )
            << "Attempt to remove sliver cell: "
            << c0 << " failed. Simulation will continue."
            << endl;

            return false;
        }
        else
        {
            return true;
        }
    }

    return false;
}

// Update mesh corresponding to the given map
void dynamicTopoFvMesh::updateMesh(const mapPolyMesh& mpm)
{
    // Delete oldPoints in polyMesh
    polyMesh::resetMotion();

    // Update polyMesh.
    polyMesh::updateMesh(mpm);

    // Clear out surface-interpolation
    surfaceInterpolation::movePoints();

    // Clear-out fvMesh geometry and addressing
    fvMesh::clearOut();

    // Update topology for all registered classes
    meshObjectBase::allUpdateTopology<fvMesh>(*this);

    // Map all fields
    mapFields(mpm);
}

// Map all fields in time using the given map
void dynamicTopoFvMesh::mapFields(const mapPolyMesh& meshMap)
{
    if (debug)
    {
        Info << "void dynamicTopoFvMesh::mapFields(const mapPolyMesh& meshMap): "
             << "Mapping fvFields."
             << endl;
    }

    //- Field mapping class
    dynamicTopoFvMeshMapper fieldMapper(*this,meshMap);

    // Map all scalar volFields in the objectRegistry
    MapGeometricFields
    <
        scalar,
        fvPatchField,
        dynamicTopoFvMeshMapper,
        volMesh
    >
    (
        fieldMapper
    );

    // Map all vector volFields in the objectRegistry
    MapGeometricFields
    <
        vector,
        fvPatchField,
        dynamicTopoFvMeshMapper,
        volMesh
    >
    (
        fieldMapper
    );

    // Map all the surfaceFields in the objectRegistry
    MapGeometricFields
    <
        scalar,
        fvsPatchField,
        dynamicTopoFvMeshMapper,
        surfaceMesh
    >
    (
        fieldMapper
    );

    // Old volumes are not mapped since interpolation is
    // performed at the same time level.
}

// Update the mesh for motion
// This routine assumes that all boundary motions have been defined
// and incorporated into the mesh for the current time-step.
void dynamicTopoFvMesh::updateMotion()
{
    if (solveForMotion_)
    {
        // Solve for motion
        movePoints(motionPtr_->newPoints());
    }
}

// MultiThreaded topology modifier [2D]
void dynamicTopoFvMesh::threadedTopoModifier2D()
{
    if (edgeModification_)
    {
        // Initialize the face stacks
        initFaceStacks();

        // Submit jobs to the work queue
        for (label i = 0; i < threader_().getNumThreads(); i++)
        {
            threader_().addToWorkQueue
                        (
                            &edgeBisectCollapse2D,
                            reinterpret_cast<void *>(&(structPtr_[i]))
                        );
        }

        // Wait for all work to complete
        threader_().waitForCompletion();
    }

    if (debug) Info << nl << "2D Edge Bisection/Collapse complete." << endl;

    // Re-initialize the face stacks
    initFaceStacks();

    // Submit jobs to the work queue
    for (label i = 0; i < threader_().getNumThreads(); i++)
    {
        threader_().addToWorkQueue
                    (
                        &swap2DEdges,
                        reinterpret_cast<void *>(&(structPtr_[i]))
                    );
    }

    // Wait for all work to complete
    threader_().waitForCompletion();

    if (debug) Info << nl << "2D Edge Swapping complete." << endl;
}

// MultiThreaded topology modifier [3D]
void dynamicTopoFvMesh::threadedTopoModifier3D()
{
    if (edgeModification_)
    {
        // Initialize the edge stacks
        initEdgeStacks();

        // Submit jobs to the work queue
        for (label i = 0; i < threader_().getNumThreads(); i++)
        {
            threader_().addToWorkQueue
                        (
                            &edgeBisectCollapse3D,
                            reinterpret_cast<void *>(&(structPtr_[i]))
                        );
        }

        // Wait for all work to complete
        threader_().waitForCompletion();
    }

    if (debug) Info << nl << "3D Edge Bisection/Collapse complete." << endl;

    // Re-initialize the edge stacks
    initEdgeStacks();

    // Submit jobs to the work queue
    for (label i = 0; i < threader_().getNumThreads(); i++)
    {
        threader_().addToWorkQueue
                    (
                        &swap3DEdges,
                        reinterpret_cast<void *>(&(structPtr_[i]))
                    );
    }

    // Wait for all work to complete
    threader_().waitForCompletion();

    if (debug) Info << nl << "3D Edge Swapping complete." << endl;
}

// Return the number of edges in the mesh.
// Override of primitiveMesh member function
label dynamicTopoFvMesh::nEdges() const
{
    if (!nEdges_)
    {
        return primitiveMesh::nEdges();
    }

    return nEdges_;
}

// Return the number of internal edges in the mesh.
// Override of primitiveMesh member function
label dynamicTopoFvMesh::nInternalEdges() const
{
    if (!nInternalEdges_)
    {
        FatalErrorIn
        (
            "dynamicTopoFvMesh::nInternalEdges()"
        )
        << "Internal edges has not been allocated."
        << abort(FatalError);
    }

    return nInternalEdges_;
}

// Return the list of edges in the mesh.
// Override of primitiveMesh member function.
const edgeList& dynamicTopoFvMesh::edges() const
{
    if (!nEdges_)
    {
        return primitiveMesh::edges();
    }

    return IOedges_;
}

// Update the mesh for topology changes
// Return true if changes have occurred
bool dynamicTopoFvMesh::updateTopology()
{
    // Calculate the edge length-scale for the mesh
    calculateLengthScale();

    // Keep a copy of existing sizes
    nOldPoints_ = nPoints_;
    nOldEdges_  = nEdges_;
    nOldFaces_  = nFaces_;
    nOldCells_  = nCells_;
    nOldInternalFaces_ = nInternalFaces_;
    nOldInternalEdges_ = nInternalEdges_;
    for(label i=0; i<numPatches_; i++)
    {
        oldPatchSizes_[i] = patchSizes_[i];
        oldPatchStarts_[i] = patchStarts_[i];
        oldEdgePatchSizes_[i] = edgePatchSizes_[i];
        oldEdgePatchStarts_[i] = edgePatchStarts_[i];
        oldPatchNMeshPoints_[i] = patchNMeshPoints_[i];
    }

    // Print out the mesh bandwidth
    if (debug)
    {
        label band=0;
        const labelList& oldOwner = faceOwner();
        const labelList& oldNeighbour = faceNeighbour();

        for(label faceI = 0; faceI < nInternalFaces_; faceI++)
        {
            label diff = oldNeighbour[faceI] - oldOwner[faceI];
            if (diff > band) band = diff;
        }

        Info << "Mesh size: " << nCells()
             << "    Bandwidth before renumbering: " << band << endl;
    }

    // Obtain the most recent point-positions
    const pointField& currentPoints = points();
    HashList<point>::iterator pIter = meshPoints_.begin();
    while (pIter != meshPoints_.end())
    {
        pIter() = currentPoints[pIter.index()];
        pIter++;
    }

    // Track mesh topology modification time
    clockTime topologyTimer;

    //== Connectivity changes ==//

    // Reset the flag
    topoChangeFlag_ = false;

    // Invoke the threaded topoModifier
    if ( twoDMesh_ )
    {
        threadedTopoModifier2D();
    }
    else
    {
        threadedTopoModifier3D();
    }

    Info << "Topo modifier time: " << topologyTimer.elapsedTime() << endl;

    clockTime reOrderingTimer;

    // Apply all pending topology changes, if necessary
    if (topoChangeFlag_)
    {

        // Allocate temporary lists for mesh-reset
        pointField points(nPoints_);
        faceList faces(nFaces_);
        labelList owner(nFaces_);
        labelList neighbour(nFaces_);

        // Null temporaries
        List<objectMap> pointsFromPoints(0);
        List<objectMap> facesFromPoints(0);
        List<objectMap> facesFromEdges(0);
        List<objectMap> facesFromFaces(0);
        List<objectMap> cellsFromPoints(0);
        List<objectMap> cellsFromEdges(0);
        List<objectMap> cellsFromFaces(0);
        List<objectMap> cellsFromCells(cellsFromCells_.size());
        labelHashSet flipFaceFlux(0);
        labelListList pointZoneMap(0);
        labelListList faceZonePointMap(0);
        labelListList faceZoneFaceMap(0);
        labelListList cellZoneMap(0);
        pointField preMotionPoints(0);

        // Reorder the mesh and obtain current topological information
        reOrderMesh(points, faces, owner, neighbour);

        // Copy cell-mapping information
        label indexI=0;
        forAllIter(Map<objectMap>, cellsFromCells_, cellI)
        {
            cellsFromCells[indexI++] = cellI();
        }

        // Obtain the patch-point labels for mapping before resetting the mesh
        labelListList oldMeshPointLabels(numPatches_);
        for(label i=0; i<numPatches_; i++)
        {
            oldMeshPointLabels[i] = boundaryMesh()[i].meshPoints();
        }

        // Obtain geometry information for mapping as well
        cellCentresPtr_.set
        (
            new vectorField(polyMesh::cellCentres())
        );

        faceCentresPtr_.set
        (
            new vectorField(polyMesh::faceCentres())
        );

        pointPositionsPtr_.set
        (
            new vectorField(currentPoints)
        );

        // Reset the mesh
        polyMesh::resetPrimitives
        (
            nFaces_,
            points,
            faces,
            owner,
            neighbour,
            patchSizes_,
            patchStarts_
        );

        // Generate mapping for points on boundary patches
        labelListList patchPointMap(numPatches_);
        for(label i=0; i<numPatches_; i++)
        {
            const labelList& meshPointLabels = boundaryMesh()[i].meshPoints();
            patchNMeshPoints_[i] = meshPointLabels.size();
            patchPointMap[i].setSize(meshPointLabels.size(), -1);
            forAll(meshPointLabels, pointI)
            {
                // Check if the position has been maintained (This saves a search).
                // Otherwise, perform a search for the old position in the patch.
                if (pointI < oldPatchNMeshPoints_[i])
                {
                    if (meshPointLabels[pointI] == oldMeshPointLabels[i][pointI])
                    {
                        // Good. Position is maintained. Make an entry
                        patchPointMap[i][pointI] = pointI;
                    }
                    else
                    {
                        // Start a linear search for the old position
                        bool foundOldPos=false;
                        forAll(oldMeshPointLabels[i],pointJ)
                        {
                            if
                            (
                                oldMeshPointLabels[i][pointJ]
                             == meshPointLabels[pointI]
                            )
                            {
                                patchPointMap[i][pointI] = pointJ;
                                foundOldPos=true; break;
                            }
                        }
                        if (!foundOldPos)
                        {
                            // Couldn't find a match. Must be a new label.
                            patchPointMap[i][pointI] = -1;
                        }
                    }
                }
            }
        }

        // Clear the existing mapper
        if (mapper_.valid()) mapper_.clear();

        // Generate new mesh mapping information
        mapper_.set
        (
            new mapPolyMesh
            (
                (*this),
                nOldPoints_,
                nOldFaces_,
                nOldCells_,
                pointMap_,
                pointsFromPoints,
                faceMap_,
                facesFromPoints,
                facesFromEdges,
                facesFromFaces,
                cellMap_,
                cellsFromPoints,
                cellsFromEdges,
                cellsFromFaces,
                cellsFromCells,
                reversePointMap_,
                reverseFaceMap_,
                reverseCellMap_,
                flipFaceFlux,
                patchPointMap,
                pointZoneMap,
                faceZonePointMap,
                faceZoneFaceMap,
                cellZoneMap,
                preMotionPoints,
                oldPatchStarts_,
                oldPatchNMeshPoints_,
                true
            )
        );

        // Update the underlying mesh, and map all related fields
        updateMesh(mapper_);

        // Discard old information after mapping
        pointPositionsPtr_.clear();
        cellCentresPtr_.clear();
        faceCentresPtr_.clear();

        // Update the motion-solver, if necessary
        if (motionPtr_.valid())
        {
            motionPtr_().updateMesh(mapper_);
        }

        // Print out the mesh bandwidth
        if (debug)
        {
            label band=0;

            for(label faceI = 0; faceI < nInternalFaces_; faceI++)
            {
                label diff = neighbour[faceI] - owner[faceI];
                if (diff > band) band = diff;
            }

            Info << "Mesh size: " << nCells()
                 << "    Bandwidth after renumbering: " << band << endl;

            Info << "----------- " << endl;
            Info << "Patch Info: " << endl;
            Info << "----------- " << endl;
            Info << "Old Patch MeshPoints: " << oldPatchNMeshPoints_ << endl;
            Info << "New Patch MeshPoints: " << patchNMeshPoints_ << endl;
            Info << "Old Patch Starts: " << mapper_->oldPatchStarts() << endl;
            Info << "Old Patch Sizes: " << mapper_->oldPatchSizes() << endl;
        }

        // Clear the current and reverse maps
        pointMap_.clear();
        edgeMap_.clear();
        faceMap_.clear();
        cellMap_.clear();
        reversePointMap_.clear();
        reverseEdgeMap_.clear();
        reverseFaceMap_.clear();
        reverseCellMap_.clear();
        addedFacePatches_.clear();
        addedEdgePatches_.clear();
        cellsFromCells_.clear();
        cellParents_.clear();

        // Set new sizes for the reverse maps
        reversePointMap_.setSize(nPoints_);
        reverseEdgeMap_.setSize(nEdges_);
        reverseFaceMap_.setSize(nFaces_);
        reverseCellMap_.setSize(nCells_);

        // Re-initialize mutex lists for multi-threading, if necessary
        initMutexLists();

        // Basic checks for mesh-validity
        if (debug)
        {
            checkMesh(true);
        }
    }

    Info << "Reordering time: " << reOrderingTimer.elapsedTime() << endl;

    return topoChangeFlag_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void dynamicTopoFvMesh::operator=(const dynamicTopoFvMesh& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "dynamicTopoFvMesh::operator=(const dynamicTopoFvMesh&)"
        )
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
