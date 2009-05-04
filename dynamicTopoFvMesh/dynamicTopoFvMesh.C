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
#include "GeometricFields.H"
#include "dynamicTopoFvMesh.H"
#include "dynamicTopoFvMeshMapper.H"
#include "multiThreader.H"
#include "triPointRef.H"
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
    mapper_(NULL),
    meshPoints_(polyMesh::points()),
    faces_(polyMesh::faces()),
    owner_(polyMesh::faceOwner()),
    neighbour_(polyMesh::faceNeighbour()),
    cells_(primitiveMesh::cells()),
    eMeshPtr_(NULL),
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
    nModifications_(0),
    maxModifications_(-1),
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
    structPtr_.setSize(nThreads);

    for (label i=0; i<nThreads; i++)
    {
        structPtr_.set(i, new topoMeshStruct);
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

    // Initialize lock-lists
    pLocks_.setSize(nThreads);
    eLocks_.setSize(nThreads);
    fLocks_.setSize(nThreads);
    cLocks_.setSize(nThreads);

    // Initialize patch-size information
    const polyBoundaryMesh& boundary = polyMesh::boundaryMesh();
    for(label i=0; i<numPatches_; i++)
    {
        patchNMeshPoints_[i] = boundary[i].meshPoints().size();
        oldPatchSizes_[i]  = patchSizes_[i]  = boundary[i].size();
        oldPatchStarts_[i] = patchStarts_[i] = boundary[i].start();
        oldPatchNMeshPoints_[i] = patchNMeshPoints_[i];
    }

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
                        dict_.subDict
                        (
                            "dynamicTopoFvMesh"
                        ).lookup("tetMetricLib")
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
            ) << nl << " Could not open the tetMetric library. "
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
            ) << nl << " Unrecognized tet-quality metric: " << tetMetric
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
                    dict_.subDict
                    (
                        "dynamicTopoFvMesh"
                    ).lookup("maxTetsPerEdge")
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
                    dict_.subDict
                    (
                        "dynamicTopoFvMesh"
                    ).lookup("allowTableResize")
                );
        }
    }

    // Initialize edge-related connectivity structures
    initEdges();

    // Set sizes for the reverse maps
    reversePointMap_.setSize(nPoints_);
    reverseEdgeMap_.setSize(nEdges_);
    reverseFaceMap_.setSize(nFaces_);
    reverseCellMap_.setSize(nCells_);

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

        if (edgeOptionDict.found("maxModifications"))
        {
            maxModifications_ =
                readLabel(edgeOptionDict.lookup("maxModifications"));
        }

        // Initialize the lengthScale field
        lengthScale_.setSize(nCells_, 0.0);
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
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return the mesh-mapper
const mapPolyMesh& dynamicTopoFvMesh::meshMap()
{
    if (mapper_.valid())
    {
        return mapper_();
    }
    else
    {
        FatalErrorIn
        (
            "dynamicTopoFvMesh::meshMap()"
        )
            << nl << " Illegal request for the mesh mapper."
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
        FatalErrorIn
        (
            "dynamicTopoFvMesh::oldCellCentres() "
        )
            << nl << " Illegal request for old cell centres."
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

    if (edgeModification_)
    {
        scalarField& internalField = tlengthScale();

        // Obtain length-scale values from the mesh
        forAllIter(HashList<scalar>, lengthScale_, lIter)
        {
            internalField[lIter.index()] = lIter();
        }
    }

    return tlengthScale;
}

// Return mesh cell-quality values
// Valid for 3D tetrahedral meshes only...
tmp<scalarField> dynamicTopoFvMesh::meshQuality
(
    bool outputOption
)
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
                    maxQuality = iF[cellI] > maxQuality
                               ? iF[cellI] : maxQuality;
                    minQuality = iF[cellI] < minQuality
                               ? iF[cellI] : minQuality;
                    meanQuality += iF[cellI];

                    break;
                }
            }
        }

        // Output statistics:
        if (outputOption)
        {
            Info << " ~~~ Mesh Quality Statistics ~~~ " << endl;
            Info << " Min: " << minQuality << endl;
            Info << " Max: " << maxQuality << endl;
            Info << " Mean: " << meanQuality/iF.size() << endl;
            Info << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << endl;

            if (minQuality < 0.0)
            {
                WarningIn
                (
                    "dynamicTopoFvMesh::meshQuality()"
                )   << nl << "Minimum cell quality is: " << minQuality << endl;
            }
        }
    }

    return tQuality;
}

// Perform a Delaunay test on an internal face
inline void dynamicTopoFvMesh::testDelaunay
(
    const label fIndex,
    bool& failed
)
{
    failed = false;
    label eIndex = -1, pIndex = -1;
    FixedList<bool,2> foundTriFace(false);
    FixedList<FixedList<label,3>,2> triFaces;

    // Figure out which thread this is...
    label tIndex = self();

    // Boundary faces are discarded.
    if (whichPatch(fIndex) > -1)
    {
        return;
    }

    // Attempt to read-lock this face
    if (tryFaceLock(fIndex))
    {
        faceStack(tIndex).push(fIndex);
        return;
    }

    const labelList& fEdges = faceEdges_[fIndex];

    forAll(fEdges, edgeI)
    {
        // Break out if both triangular faces are found
        if (foundTriFace[0] && foundTriFace[1])
        {
            break;
        }

        // Obtain edgeFaces for this edge
        const labelList& eFaces = edgeFaces_[fEdges[edgeI]];

        forAll(eFaces, faceI)
        {
            face& thisFace = faces_[eFaces[faceI]];

            if (thisFace.size() == 3)
            {
                if (foundTriFace[0])
                {
                    // Update the second face.
                    triFaces[1][0] = thisFace[0];
                    triFaces[1][1] = thisFace[1];
                    triFaces[1][2] = thisFace[2];

                    foundTriFace[1] = true;

                    // Take this edge
                    eIndex = fEdges[edgeI];
                }
                else
                {
                    // Update the first face.
                    triFaces[0][0] = thisFace[0];
                    triFaces[0][1] = thisFace[1];
                    triFaces[0][2] = thisFace[2];

                    foundTriFace[0] = true;
                }
            }
        }
    }

    // Obtain point references for the first face
    point& a = meshPoints_[triFaces[0][0]];
    point& b = meshPoints_[triFaces[0][1]];
    point& c = meshPoints_[triFaces[0][2]];

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
            "dynamicTopoFvMesh::testDelaunay(const label fIndex) "
        ) << nl << " Encountered a co-linear set of points: " << nl
                << " Point a :: " << triFaces[0][0] << ": " << a << nl
                << " Point b :: " << triFaces[0][1] << ": " << b << nl
                << " Point c :: " << triFaces[0][2] << ": " << c << nl
                << abort(FatalError);
    }
#   endif

    point cCenter = ((c2 + c3)*a + (c3 + c1)*b + (c1 + c2)*c)/(2*cd);
    scalar rSquared = (a - cCenter)&(a - cCenter);

    // Find the isolated point on the second face
    edge& e = edges_[eIndex];

    // Check the first point
    if (triFaces[1][0] != e.start() && triFaces[1][0] != e.end())
    {
        pIndex = triFaces[1][0];
    }

    // Check the second point
    if (triFaces[1][1] != e.start() && triFaces[1][1] != e.end())
    {
        pIndex = triFaces[1][1];
    }

    // Check the third point
    if (triFaces[1][2] != e.start() && triFaces[1][2] != e.end())
    {
        pIndex = triFaces[1][2];
    }

    // ...and determine whether it lies in this circle
    point& otherPoint = meshPoints_[pIndex];

    if
    (
        ((otherPoint - cCenter)&(otherPoint - cCenter)) < rSquared
    )
    {
        // Failed the test.
        failed = true;
    }
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

// Compare two triangular faces.
// Identical to triFace definition.
inline label dynamicTopoFvMesh::compare
(
    const face& a,
    const face& b
)
{
    if
    (
        (a[0] == b[0] && a[1] == b[1] && a[2] == b[2])
     || (a[0] == b[1] && a[1] == b[2] && a[2] == b[0])
     || (a[0] == b[2] && a[1] == b[0] && a[2] == b[1])
    )
    {
        // Identical
        return 1;
    }
    else if
    (
        (a[0] == b[2] && a[1] == b[1] && a[2] == b[0])
     || (a[0] == b[1] && a[1] == b[0] && a[2] == b[2])
     || (a[0] == b[0] && a[1] == b[2] && a[2] == b[1])
    )
    {
        // Same face, but reversed orientation
        return -1;
    }
    else
    {
        // Faces don't match.
        return 0;
    }
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
)
{
    if (index < nOldInternalFaces_)
    {
        return -1;
    }

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
            << "Cannot find patch information for face index: " << index << nl
            << " It appears that face ordering is"
            << " inconsistent with patch information."
            << abort(FatalError);
    }

    return -2;
}

// Method to determine the old boundary patch index for a given edge
inline label dynamicTopoFvMesh::whichEdgePatch
(
    const label& index
)
{
    if (index < nOldInternalEdges_)
    {
        return -1;
    }

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
            << "Cannot find patch information for edge index: " << index << nl
            << " It appears that edge ordering is"
            << " inconsistent with patch information."
            << abort(FatalError);
    }

    return -2;
}

// Utility method to find the interior/boundary faces
// for an input quad-face and adjacent triangle-prism cell.
void dynamicTopoFvMesh::findPrismFaces
(
    const label fIndex,
    const label cIndex,
    FixedList<face,2>& bdyf,
    FixedList<label,2>& bidx,
    FixedList<face,2>& intf,
    FixedList<label,2>& iidx
)
{
    label indexO = 0, indexI = 0;

    cell& c = cells_[cIndex];

    forAll(c, i)
    {
        label faceIndex = c[i];

        // Don't count the face under consideration
        if (faceIndex != fIndex)
        {
            face& fi = faces_[faceIndex];

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
bool dynamicTopoFvMesh::findCommonEdge
(
    const label first,
    const label second,
    label& common
)
{
    bool found = false;

    labelList& fEi = faceEdges_[first];
    labelList& fEj = faceEdges_[second];

    forAll(fEi, edgeI)
    {
        forAll(fEj, edgeJ)
        {
            if (fEi[edgeI] == fEj[edgeJ])
            {
                common = fEi[edgeI];

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
    // Check the first point
    if ( f[0] != e.start() && f[0] != e.end() )
    {
        ptIndex = f[0];
        nextPtIndex = f[1];
        return;
    }

    // Check the second point
    if ( f[1] != e.start() && f[1] != e.end() )
    {
        ptIndex = f[1];
        nextPtIndex = f[2];
        return;
    }

    // Check the third point
    if ( f[2] != e.start() && f[2] != e.end() )
    {
        ptIndex = f[2];
        nextPtIndex = f[0];
        return;
    }

    // This bit should never happen.
    FatalErrorIn
    (
        "label dynamicTopoFvMesh::findIsolatedPoint()"
    )
        << "Cannot find isolated point in face " << f << endl
        << " Using edge: " << e
        << abort(FatalError);
}

// Utility method to replace a label in a given list
inline void dynamicTopoFvMesh::replaceLabel
(
     const label original,
     const label replacement,
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

// Method to insert a label between two labels in a list
// Assumes that all labels are unique.
inline void dynamicTopoFvMesh::insertLabel
(
    const label newLabel,
    const label labelA,
    const label labelB,
    labelList& list
)
{
    // Create a new list
    bool found = false;
    label origSize = list.size();
    labelList newList(origSize + 1);

    label index = 0, nextI = -1;

    // Start a linear search
    forAll(list, itemI)
    {
        newList[index++] = list[itemI];

        nextI = list.fcIndex(itemI);

        if
        (
           ((list[itemI] == labelA && list[nextI] == labelB)
         || (list[itemI] == labelB && list[nextI] == labelA))
         && !found
        )
        {
            found = true;
            newList[index++] = newLabel;
        }
    }

    if (!found)
    {
        FatalErrorIn
        (
            "label dynamicTopoFvMesh::insertLabel()"
        )   << "Cannot insert " << newLabel << " in list: " << list << endl
            << " Labels: "
            << labelA << " and " << labelB << " were not found in sequence."
            << abort(FatalError);
    }

    // Transfer the list
    list.transfer(newList);
}

// Utility method for face-insertion
label dynamicTopoFvMesh::insertFace
(
    const label patch,
    const face& newFace,
    const label newOwner,
    const label newNeighbour
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
        faceStack_[self()].push(newFaceIndex);
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
    faceEdges_.remove(index);

    // Remove the face mutex
    if (threader_->multiThreaded())
    {
        // Unlock it first
        faceMutex_[index].unlock();
        faceMutex_.remove(index);
    }

    // Identify the patch for this face
    label patch = whichPatch(index);

    if (patch >= 0)
    {
        // Modify patch information for this boundary face
        patchSizes_[patch]--;

        for(label i = (patch + 1); i < numPatches_; i++)
        {
            patchStarts_[i]--;
        }
    }
    else
    {
        // Decrement the internal face count, and subsequent patch-starts
        nInternalFaces_--;

        forAll(patchStarts_, patchI)
        {
            patchStarts_[patchI]--;
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
    const labelList& edgeFaces,
    const labelList& edgePoints
)
{
    label newEdgeIndex = edges_.append(newEdge);
    edgeFaces_.append(edgeFaces);

    if (!twoDMesh_)
    {
        edgePoints_.append(edgePoints);
    }

    // Add a new unlocked edge mutex
    if (threader_->multiThreaded())
    {
        edgeMutex_.append();
    }

    // Add to the stack as well
    edgeStack(self()).push(newEdgeIndex);

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

    if (!twoDMesh_)
    {
        // Remove the edgePoints entry
        edgePoints_.remove(index);

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
        edgeStack(stackI).remove(index);
    }

    // Identify the patch for this edge
    label patch = whichEdgePatch(index);

    if (patch >= 0)
    {
        // Modify patch information for this boundary edge
        edgePatchSizes_[patch]--;

        for(label i = (patch + 1); i < numPatches_; i++)
        {
            edgePatchStarts_[i]--;
        }
    }
    else
    {
        // Decrement the internal edge count, and subsequent patch-starts
        nInternalEdges_--;

        forAll(edgePatchStarts_, patchI)
        {
            edgePatchStarts_[patchI]--;
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

// Check for the occurrence of a label in the list
inline bool dynamicTopoFvMesh::foundInList
(
    const label item,
    const labelList& list
)
{
    forAll(list, itemI)
    {
        if (list[itemI] == item)
        {
            return true;
        }
    }

    return false;
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

// Utility method to build a hull of cells connected to the edge [2D]
void dynamicTopoFvMesh::constructPrismHull
(
    const label eIndex,
    labelHashSet& hullTriFaces,
    labelHashSet& hullCells
)
{
    // Obtain references
    labelList& eFaces = edgeFaces_[eIndex];

    // Loop through edgeFaces and add cells
    forAll(eFaces, faceI)
    {
        label c0 = owner_[eFaces[faceI]];
        label c1 = neighbour_[eFaces[faceI]];

        if
        (
            !hullCells.found(c0)
        )
        {
            // Add this cell
            hullCells.insert(c0);

            // Find associated triFaces and add them too
            cell& cellToCheck = cells_[c0];

            forAll(cellToCheck, faceJ)
            {
                face& faceToCheck = faces_[cellToCheck[faceJ]];

                if
                (
                    (faceToCheck.size() == 3)
                 && !(hullTriFaces.found(cellToCheck[faceJ]))
                )
                {
                    hullTriFaces.insert(cellToCheck[faceJ]);
                }
            }
        }

        if
        (
            !hullCells.found(c1)
         && (c1 != -1)
        )
        {
            // Add this cell
            hullCells.insert(c1);

            // Find associated triFaces and add them too
            cell& cellToCheck = cells_[c1];

            forAll(cellToCheck, faceJ)
            {
                face& faceToCheck = faces_[cellToCheck[faceJ]];

                if
                (
                    (faceToCheck.size() == 3)
                 && !(hullTriFaces.found(cellToCheck[faceJ]))
                )
                {
                    hullTriFaces.insert(cellToCheck[faceJ]);
                }
            }
        }
    }
}

// Utility method to build a hull of cells (and faces) around an edge.
// Assumes that eIndex is already locked.
bool dynamicTopoFvMesh::constructHull
(
    const label eIndex,
    labelList& hullEdges,
    labelList& hullFaces,
    labelList& hullCells,
    labelListList& ringEntities,
    const rwMutex::lockType lType
)
{
    // [1] hullEdges is an ordered list of edge-labels around eIndex,
    //     but not connected to it.
    //      - Ordering is in the same manner as edgePoints.
    // [2] hullFaces is an ordered list of face-labels connected to eIndex.
    //      - Ordering is in the same manner as edgePoints.
    // [3] hullCells is an ordered list of cell-labels connected to eIndex.
    //      - For boundary hulls, the last cell label is -1
    // [4] ringEntities are edges and faces connected to eIndex[0] and eIndex[1]
    //      - ringEntities[0]: edges connected to eIndex[0]
    //      - ringEntities[1]: faces connected to eIndex[0]
    //      - ringEntities[2]: edges connected to eIndex[1]
    //      - ringEntities[3]: faces connected to eIndex[1]

    bool found;
    label otherPoint = -1, nextPoint = -1;

    // Figure out which thread this is...
    label tIndex = self();

    // Obtain a reference to this edge, and its edgeFaces
    edge& edgeToCheck = edges_[eIndex];
    labelList& eFaces = edgeFaces_[eIndex];
    labelList& hullVertices = edgePoints_[eIndex];

    // Try to lock the two points of this edge
    if (tryEdgePointLock(eIndex, lType))
    {
        edgeStack(tIndex).push(eIndex);
        return true;
    }

    // Try to lock hull-points around this edge
    if (tryEdgeHullLock(eIndex, lType))
    {
        edgeStack(tIndex).push(eIndex);
        return true;
    }

    // Temporary tri-face for comparison
    face oFace(3);

    // Loop through all faces of this edge and add them to hullFaces
    forAll(eFaces, faceI)
    {
        face& faceToCheck = faces_[eFaces[faceI]];

        // Find the isolated point on this face,
        // and compare it with hullVertices
        findIsolatedPoint
        (
            faceToCheck,
            edgeToCheck,
            otherPoint,
            nextPoint
        );

        found = false;

        forAll(hullVertices, indexI)
        {
            if (hullVertices[indexI] == otherPoint)
            {
                // Fill in the position of this face on the hull
                hullFaces[indexI] = eFaces[faceI];

                // Try to lock edges of this face.
                // Also, avoid locking eIndex multiple times.
                if (tryFaceEdgeLock(hullFaces[indexI], lType, eIndex))
                {
                    edgeStack(tIndex).push(eIndex);
                    return true;
                }

                // Obtain edges connected to top and bottom
                // vertices of edgeToCheck
                labelList& fEdges = faceEdges_[hullFaces[indexI]];

                forAll(fEdges, edgeI)
                {
                    if
                    (
                        edges_[fEdges[edgeI]]
                     == edge(edgeToCheck[0], otherPoint)
                    )
                    {
                        ringEntities[0][indexI] = fEdges[edgeI];
                    }

                    if
                    (
                        edges_[fEdges[edgeI]]
                     == edge(edgeToCheck[1], otherPoint)
                    )
                    {
                        ringEntities[2][indexI] = fEdges[edgeI];
                    }
                }

                // Depending on the orientation of this face,
                // fill in hull cell indices as well
                if (nextPoint == edgeToCheck[0])
                {
                    hullCells[indexI] = owner_[eFaces[faceI]];
                }
                else
                if (nextPoint == edgeToCheck[1])
                {
                    hullCells[indexI] = neighbour_[eFaces[faceI]];
                }
                else
                {
                    // Something's terribly wrong
                    FatalErrorIn
                    (
                        "dynamicTopoFvMesh::constructHull()"
                    )
                        << nl << " Failed to construct hull. "
                        << nl << " Possibly not a tetrahedral mesh. "
                        << abort(FatalError);
                }

                if (hullCells[indexI] != -1)
                {
                    label nextI = hullVertices.fcIndex(indexI);
                    label nextHullPoint = hullVertices[nextI];
                    cell& currCell = cells_[hullCells[indexI]];

                    // Look for the ring-faces
                    forAll(currCell, faceI)
                    {
                        face& cFace = faces_[currCell[faceI]];

                        // Build a comparison face
                        oFace[0] = edgeToCheck[0];
                        oFace[1] = otherPoint;
                        oFace[2] = nextHullPoint;

                        // Check if this face contains edgeToCheck[0]
                        if (compare(cFace, oFace))
                        {
                            ringEntities[1][indexI] = currCell[faceI];
                        }

                        // Build a comparison face
                        oFace[0] = edgeToCheck[1];
                        oFace[1] = nextHullPoint;
                        oFace[2] = otherPoint;

                        // Check if this face contains edgeToCheck[1]
                        if (compare(cFace, oFace))
                        {
                            ringEntities[3][indexI] = currCell[faceI];
                        }
                    }

                    // Scan one the faces for the ring-edge
                    labelList& rFaceEdges = faceEdges_[ringEntities[1][indexI]];

                    forAll(rFaceEdges, edgeI)
                    {
                        if
                        (
                            edges_[rFaceEdges[edgeI]]
                         == edge(otherPoint,nextHullPoint)
                        )
                        {
                            hullEdges[indexI] = rFaceEdges[edgeI];
                            break;
                        }
                    }

                    // Try to lock the ring-edge.
                    if (tryEdgeLock(hullEdges[indexI],lType))
                    {
                        edgeStack(tIndex).push(eIndex);
                        return true;
                    }
                }

                // Done with this index. Break out.
                found = true; break;
            }
        }

        // Throw an error if the point wasn't found
        if (!found)
        {
            // Something's terribly wrong
            FatalErrorIn
            (
                "dynamicTopoFvMesh::constructHull()"
            )
                << " Failed to construct hull. " << nl
                << " edgeFaces connectivity is inconsistent. " << nl
                << " Edge: " << eIndex << ":: " << edgeToCheck << nl
                << " edgeFaces: " << eFaces << nl
                << " edgePoints: " << hullVertices
                << abort(FatalError);
        }
    }

    // Return a successful lock
    return false;
}

// Utility method to build edgePoints for an edge [3D].
// Assumes that edgeFaces information is consistent.
void dynamicTopoFvMesh::buildEdgePoints
(
    const label eIndex
)
{
    bool found = false;
    label faceIndex = -1, cellIndex = -1;
    label otherPoint = -1, nextPoint = -1;

    // Obtain references
    edge& edgeToCheck = edges_[eIndex];
    labelList& eFaces = edgeFaces_[eIndex];
    labelList& ePoints = edgePoints_[eIndex];

    // Re-size the list first
    ePoints.setSize(eFaces.size(), -1);

    if (whichEdgePatch(eIndex) == -1)
    {
        // Internal edge.
        // Pick the first face and start with that
        faceIndex = eFaces[0];
    }
    else
    {
        // Need to find a properly oriented start-face
        forAll(eFaces, faceI)
        {
            if (whichPatch(eFaces[faceI]) > -1)
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
                    faceIndex = eFaces[faceI];
                    break;
                }
            }
        }
    }

    // Shuffle vertices to appear in CCW order
    forAll(ePoints, indexI)
    {

        findIsolatedPoint
        (
            faces_[faceIndex],
            edgeToCheck,
            otherPoint,
            nextPoint
        );

        // Add the isolated point
        ePoints[indexI] = otherPoint;

        // Figure out how this edge is oriented.
        if (nextPoint == edgeToCheck[0])
        {
            // Counter-clockwise. Pick the owner.
            cellIndex = owner_[faceIndex];
        }
        else
        if (whichPatch(faceIndex) == -1)
        {
            // Clockwise. Pick the neighbour.
            cellIndex = neighbour_[faceIndex];
        }
        else
        {
            // Looks like we've hit a boundary face. Break out.
            break;
        }

        const cell& cellToCheck = cells_[cellIndex];

        found = false;

        // Assuming tet-cells,
        // Loop through edgeFaces and get the next face
        forAll(eFaces, faceI)
        {
            if
            (
                eFaces[faceI] != faceIndex
             && eFaces[faceI] == cellToCheck[0]
            )
            {
                faceIndex = cellToCheck[0];
                found = true; break;
            }

            if
            (
                eFaces[faceI] != faceIndex
             && eFaces[faceI] == cellToCheck[1]
            )
            {
                faceIndex = cellToCheck[1];
                found = true; break;
            }

            if
            (
                eFaces[faceI] != faceIndex
             && eFaces[faceI] == cellToCheck[2]
            )
            {
                faceIndex = cellToCheck[2];
                found = true; break;
            }

            if
            (
                eFaces[faceI] != faceIndex
             && eFaces[faceI] == cellToCheck[3]
            )
            {
                faceIndex = cellToCheck[3];
                found = true; break;
            }
        }

#       ifdef FULLDEBUG
        if (!found)
        {
            // Something's terribly wrong
            FatalErrorIn
            (
                "void dynamicTopoFvMesh::buildEdgePoints"
            )
                << " Failed to determine a vertex ring. " << nl
                << " edgeFaces connectivity is inconsistent. " << nl
                << "Edge: " << eIndex << ":: " << edgeToCheck << nl
                << "edgeFaces: " << eFaces
                << abort(FatalError);
        }
#       endif
    }
}

// Utility to check whether points of an edge lie on a boundary.
// Assumes that eIndex is already locked.
void dynamicTopoFvMesh::checkEdgeBoundary
(
    const label eIndex,
    FixedList<bool,2>& edgeBoundary
)
{
    edge& edgeToCheck = edges_[eIndex];

    // Loop through edges connected to both points,
    // and check if any of them lie on boundaries.
    // Used to ensure that collapses happen towards boundaries.
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
    // Used to ensure that collapses happen towards bounding curves.
    // Note that if the edge itself is on a bounding curve, collapse is valid.
    if (edgeBoundary[0] && edgeBoundary[1])
    {
        FixedList<label, 2> nBoundCurves(0);

        forAll(edgeToCheck, pointI)
        {
            labelList& pEdges = pointEdges_[edgeToCheck[pointI]];

            forAll(pEdges, edgeI)
            {
                if (checkBoundingCurve(pEdges[edgeI]))
                {
                    nBoundCurves[pointI]++;
                }
            }
        }

        // Pick the point which is connected to more bounding curves
        if (nBoundCurves[0] > nBoundCurves[1])
        {
            edgeBoundary[0] = true;
            edgeBoundary[1] = false;
        }
        else
        if (nBoundCurves[1] > nBoundCurves[0])
        {
            edgeBoundary[1] = true;
            edgeBoundary[0] = false;
        }
        else
        {
            // Bounding edge: collapseEdge can collapse this to the mid-point
            edgeBoundary = true;
        }
    }
}

// Utility method to compute the minimum quality of a vertex hull
inline bool dynamicTopoFvMesh::computeMinQuality
(
    const label eIndex,
    scalar& minQuality
)
{
    minQuality = GREAT;
    scalar cQuality = 0.0;

    // Figure out which thread this is...
    label tIndex = self();

    // Try to read-lock this edge.
    if (tryEdgeLock(eIndex))
    {
        edgeStack(tIndex).push(eIndex);
        return true;
    }

    // Obtain a reference to this edge and corresponding edgePoints
    edge& edgeToCheck = edges_[eIndex];
    labelList& hullVertices = edgePoints_[eIndex];

    // Try to read-lock the two points of this edge.
    if (tryEdgePointLock(eIndex))
    {
        edgeStack(tIndex).push(eIndex);
        return true;
    }

    // Try to read-lock hull-points around this edge
    if (tryEdgeHullLock(eIndex))
    {
        edgeStack(tIndex).push(eIndex);
        return true;
    }

    // Assume that hullVertices is ordered
    // in a CCW order around edgeToCheck[0]
    forAll(hullVertices, indexI)
    {
        // Compute quality as well
        if (indexI != 0)
        {
            cQuality = (*tetMetric_)
            (
                meshPoints_[edgeToCheck[0]],
                meshPoints_[hullVertices[indexI-1]],
                meshPoints_[edgeToCheck[1]],
                meshPoints_[hullVertices[indexI]]
            );

            // Check if the quality is worse
            minQuality = cQuality < minQuality ? cQuality : minQuality;
        }
    }

    // Ensure that the mesh is valid
    if (minQuality < 0.0)
    {
        FatalErrorIn("dynamicTopoFvMesh::computeMinQuality()")
            << "Encountered negative cell-quality!" << nl
            << "Edge: " << eIndex << ": " << edgeToCheck << nl
            << "EdgePoints: " << hullVertices << nl
            << "Minimum Quality: " << minQuality
            << abort(FatalError);
    }

    // Undo all read-locks
    unlockMutexLists(tIndex);

    // Return a successful lock
    return false;
}

// Utility to obtain a read/write lock for an entity.
// Returns true on success.
inline bool dynamicTopoFvMesh::obtainLock
(
    const label tIndex,
    const rwMutex::lockType lType,
    rwMutex& entityMutex
)
{
    label numRetries = 0;

    while (entityMutex.tryLock(lType))
    {
        // Failed to acquire this mutex. Retry.
        if (lType == rwMutex::WRITE_LOCK && numRetries > 10)
        {
            // Unsuccessful after several retries. Bail out.
            return false;
        }
        else
        {
            numRetries++;
        }
    }

    // Successful lock.
    return true;
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

    // Figure out which thread this is...
    label tIndex = self();

    // Try to read-lock this edge.
    if (tryEdgeLock(eIndex))
    {
        // Couldn't lock this edge. Put it back on the stack and bail out.
        edgeStack(tIndex).push(eIndex);
        return true;
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
        // Undo all read-locks and return
        unlockMutexLists(tIndex);
        return true;
    }

    // Undo all read-locks
    unlockMutexLists(tIndex);

    // Not on a bounding curve
    return false;
}

// Allocate dynamic programming tables
void dynamicTopoFvMesh::initTables
(
    scalarListList& Q,
    labelListList& K,
    labelListList& triangulations
)
{
    label mMax = maxTetsPerEdge_;

    Q.setSize((mMax-2),scalarList(mMax,-1.0));
    K.setSize((mMax-2),labelList(mMax,-1));
    triangulations.setSize(3,labelList((mMax-2),-1));
}

// Utility method to fill the dynamic programming tables
//  - Returns true if the operation completed successfully.
//  - Returns false if read-locks were not obtained, or
//    tables could not be resized.
bool dynamicTopoFvMesh::fillTables
(
    const label eIndex,
    const scalar minQuality,
    label& m,
    scalarListList& Q,
    labelListList& K,
    labelListList& triangulations
)
{
    // Figure out which thread this is...
    label tIndex = self();

    // Try to read-lock this edge.
    if (tryEdgeLock(eIndex))
    {
        edgeStack(tIndex).push(eIndex);
        return false;
    }

    // Try to read-lock the two points of this edge.
    if (tryEdgePointLock(eIndex))
    {
        edgeStack(tIndex).push(eIndex);
        return true;
    }

    // Try to read-lock hull-points around this edge
    if (tryEdgeHullLock(eIndex))
    {
        edgeStack(tIndex).push(eIndex);
        return true;
    }

    edge& edgeToCheck = edges_[eIndex];
    labelList& hullVertices = edgePoints_[eIndex];

    // Fill in the size
    m = hullVertices.size();

    // Check if a table-resize is necessary
    if (m > maxTetsPerEdge_)
    {
        if (allowTableResize_)
        {
            // Resize the tables to account for
            // more tets per edge
            maxTetsPerEdge_ = m;
            Q.clear(); K.clear(); triangulations.clear();
            initTables(Q, K, triangulations);
        }
        else
        {
            // Undo all read-locks
            unlockMutexLists(tIndex);

            // Can't resize. Bail out.
            return false;
        }
    }

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

    // Undo all read-locks
    unlockMutexLists(tIndex);

    return true;
}

// Remove the edge according to the swap sequence
void dynamicTopoFvMesh::removeEdgeFlips
(
    const label eIndex,
    const labelListList& K,
    labelListList& triangulations
)
{
    // Figure out which thread this is...
    label tIndex = self();

    // Try to write-lock this edge.
    if (tryEdgeLock(eIndex, rwMutex::WRITE_LOCK))
    {
        edgeStack(tIndex).push(eIndex);
        return;
    }

    // Make a copy of edgePoints, since it will be
    // modified during swaps
    labelList hullVertices(edgePoints_[eIndex]);

    label m = hullVertices.size();
    labelList hullFaces(m, -1);
    labelList hullCells(m, -1);
    labelList hullEdges(m, -1);
    labelListList ringEntities(4, labelList(m, -1));

    // Construct the hull, and try to write-lock entities.
    if
    (
        constructHull
        (
            eIndex,
            hullEdges,
            hullFaces,
            hullCells,
            ringEntities,
            rwMutex::WRITE_LOCK
        )
    )
    {
        // Unable to write-lock.
        edgeStack(tIndex).push(eIndex);
        return;
    }

    label numTriangulations = 0, isolatedVertex = -1;

    // Extract the appropriate triangulations
    extractTriangulation(0, (m-1), K, numTriangulations, triangulations);

    // Determine the 3-2 swap triangulation
    label t32 = identify32Swap(eIndex, hullVertices, triangulations);

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
                        i,
                        triangulations,
                        hullVertices,
                        hullFaces,
                        hullCells
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

    // Perform the final 3-2 / 2-2 swap
    swap32
    (
        eIndex,
        t32,
        triangulations,
        hullVertices,
        hullFaces,
        hullCells
    );

    // Done with this face, so reset it
    triangulations[0][t32] = -1;
    triangulations[1][t32] = -1;
    triangulations[2][t32] = -1;

    // Write lock the edge mutex
    eMutex_.lock(rwMutex::WRITE_LOCK);

    // Finally remove the edge
    removeEdge(eIndex);

    // Unlock the edge mutex
    eMutex_.unlock();

    // Set the flag
    topoChangeFlag_ = true;
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
    const label eIndex,
    const labelList& hullVertices,
    const labelListList& triangulations
)
{
    label m = hullVertices.size();
    edge& edgeToCheck = edges_[eIndex];

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
    const label eIndex,
    const label triangulationIndex,
    const labelListList& triangulations,
    const labelList& hullVertices,
    const labelList& hullFaces,
    const labelList& hullCells
)
{
    // A 2-3 swap performs the following operations:
    //      [1] Remove face: [ edge[0] edge[1] isolatedVertex ]
    //      [2] Remove two cells on either side of removed face
    //      [3] Add one edge
    //      [4] Add three new faces
    //      [5] Add three new cells
    //      Update faceEdges, edgeFaces and edgePoints information

    // Figure out which edge this is...
    label tIndex = self();

    // Obtain edge reference
    // Assume that write-locks have been obtained
    edge& edgeToCheck = edges_[eIndex];

#   ifdef FULLDEBUG
    if (debug)
    {
        // Print out arguments
        Info << endl;
        Info << "== Swapping 2-3 ==" << endl;
        Info << "Edge: " << eIndex << ": " << edgeToCheck << endl;
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

    // Write lock the cell mutex
    cMutex_.lock(rwMutex::WRITE_LOCK);

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

    // Unlock the cell mutex from write-lock
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

    // Write lock the face mutex
    fMutex_.lock(rwMutex::WRITE_LOCK);

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

    // Unlock the face mutex from write lock
    fMutex_.unlock();

    // Add an entry to edgeFaces
    labelList newEdgeFaces(3);
    newEdgeFaces[0] = newFaceIndex[0];
    newEdgeFaces[1] = newFaceIndex[1];
    newEdgeFaces[2] = newFaceIndex[2];

    // Add an entry for edgePoints as well
    labelList newEdgePoints(3);
    newEdgePoints[0] = vertexForRemoval;
    newEdgePoints[1] = edgeToCheck[0];
    newEdgePoints[2] = edgeToCheck[1];

    // Write lock the edge mutex
    eMutex_.lock(rwMutex::WRITE_LOCK);

    // Add a new internal edge to the mesh
    label newEdgeIndex = insertEdge
                         (
                             -1,
                             edge
                             (
                                 otherVertices[0],
                                 otherVertices[1]
                             ),
                             newEdgeFaces,
                             newEdgePoints
                         );

    // Unlock the edge mutex from write lock
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

            // Check if face contains edgeToCheck[0]
            if
            (
                (faceToCheck[0] == edgeToCheck[0])
             || (faceToCheck[1] == edgeToCheck[0])
             || (faceToCheck[2] == edgeToCheck[0])
            )
            {
                foundEdge[0] = true;
            }

            // Check if face contains edgeToCheck[1]
            if
            (
                (faceToCheck[0] == edgeToCheck[1])
             || (faceToCheck[1] == edgeToCheck[1])
             || (faceToCheck[2] == edgeToCheck[1])
            )
            {
                foundEdge[1] = true;
            }

            // Face is connected to edgeToCheck[0]
            if (foundEdge[0] && !foundEdge[1])
            {
                // Check if a face-flip is necessary
                if (owner_[faceIndex] == cellIndex)
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

                // Update faceEdges, edgeFaces, and edgePoints.
                // Add them to the stack as well.
                const labelList& fEdges = faceEdges_[faceIndex];

                forAll(fEdges, edgeI)
                {
                    if (edges_[fEdges[edgeI]] == check[0])
                    {
                        newFaceEdges[0][nE0++] = fEdges[edgeI];

                        sizeUpList
                        (
                            newFaceIndex[0],
                            edgeFaces_[fEdges[edgeI]]
                        );

                        insertLabel
                        (
                            otherVertices[1],
                            edgeToCheck[0],
                            edgeToCheck[1],
                            edgePoints_[fEdges[edgeI]]
                        );

                        edgeStack(tIndex).push(fEdges[edgeI]);
                    }

                    if (edges_[fEdges[edgeI]] == check[1])
                    {
                        newFaceEdges[0][nE0++] = fEdges[edgeI];

                        sizeUpList
                        (
                            newFaceIndex[0],
                            edgeFaces_[fEdges[edgeI]]
                        );

                        insertLabel
                        (
                            otherVertices[0],
                            edgeToCheck[0],
                            edgeToCheck[1],
                            edgePoints_[fEdges[edgeI]]
                        );

                        edgeStack(tIndex).push(fEdges[edgeI]);
                    }

                    if (edges_[fEdges[edgeI]] == check[2])
                    {
                        newFaceEdges[1][nE1++] = fEdges[edgeI];

                        sizeUpList
                        (
                            newFaceIndex[1],
                            edgeFaces_[fEdges[edgeI]]
                        );

                        insertLabel
                        (
                            otherVertices[1],
                            vertexForRemoval,
                            edgeToCheck[1],
                            edgePoints_[fEdges[edgeI]]
                        );

                        edgeStack(tIndex).push(fEdges[edgeI]);
                    }

                    if (edges_[fEdges[edgeI]] == check[4])
                    {
                        newFaceEdges[1][nE1++] = fEdges[edgeI];

                        sizeUpList
                        (
                            newFaceIndex[1],
                            edgeFaces_[fEdges[edgeI]]
                        );

                        insertLabel
                        (
                            otherVertices[0],
                            vertexForRemoval,
                            edgeToCheck[1],
                            edgePoints_[fEdges[edgeI]]
                        );

                        edgeStack(tIndex).push(fEdges[edgeI]);
                    }
                }
            }

            // Face is connected to edgeToCheck[1]
            if (!foundEdge[0] && foundEdge[1])
            {
                // Check if a face-flip is necessary
                if (owner_[faceIndex] == cellIndex)
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

                // Update faceEdges, edgeFaces, and edgePoints.
                // Add them to the stack as well.
                const labelList& fEdges = faceEdges_[faceIndex];

                forAll(fEdges, edgeI)
                {
                    if (edges_[fEdges[edgeI]] == check[3])
                    {
                        newFaceEdges[2][nE2++] = fEdges[edgeI];

                        sizeUpList
                        (
                            newFaceIndex[2],
                            edgeFaces_[fEdges[edgeI]]
                        );

                        insertLabel
                        (
                            otherVertices[1],
                            vertexForRemoval,
                            edgeToCheck[0],
                            edgePoints_[fEdges[edgeI]]
                        );

                        edgeStack(tIndex).push(fEdges[edgeI]);
                    }

                    if (edges_[fEdges[edgeI]] == check[5])
                    {
                        newFaceEdges[2][nE2++] = fEdges[edgeI];

                        sizeUpList
                        (
                            newFaceIndex[2],
                            edgeFaces_[fEdges[edgeI]]
                        );

                        insertLabel
                        (
                            otherVertices[0],
                            vertexForRemoval,
                            edgeToCheck[0],
                            edgePoints_[fEdges[edgeI]]
                        );

                        edgeStack(tIndex).push(fEdges[edgeI]);
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
                if (owner_[faceIndex] == cellIndex)
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

    forAll(cellsForRemoval, indexI)
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

    // Write lock the cell mutex
    cMutex_.lock(rwMutex::WRITE_LOCK);

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

    // Unlock the cell mutex from write lock
    cMutex_.unlock();

    // Write lock the face mutex
    fMutex_.lock(rwMutex::WRITE_LOCK);

    // Update edgeFaces and edgePoints for edges of the removed face
    label otherPoint = -1, nextPoint = -1;
    face& checkFace = faces_[faceForRemoval];
    labelList& fEdges = faceEdges_[faceForRemoval];

    forAll(fEdges, edgeI)
    {
        sizeDownList
        (
            faceForRemoval,
            edgeFaces_[fEdges[edgeI]]
        );

        // Find the isolated point and remove it
        findIsolatedPoint
        (
            checkFace,
            edges_[fEdges[edgeI]],
            otherPoint,
            nextPoint
        );

        sizeDownList
        (
            otherPoint,
            edgePoints_[fEdges[edgeI]]
        );

        edgeStack(tIndex).push(fEdges[edgeI]);
    }

    // Remove the face
    removeFace(faceForRemoval);

    // Unlock the face mutex from write lock
    fMutex_.unlock();

    // Write lock the cell mutex
    cMutex_.lock(rwMutex::WRITE_LOCK);

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

    // Unlock the cell mutex from write lock
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
    const label eIndex,
    const label triangulationIndex,
    const labelListList& triangulations,
    const labelList& hullVertices,
    const labelList& hullFaces,
    const labelList& hullCells
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
    //      Update faceEdges, edgeFaces and edgePoints information

    // Figure out which edge this is...
    label tIndex = self();

    // Obtain edge reference
    // Assume that write-locks have been obtained
    edge& edgeToCheck = edges_[eIndex];

    // Determine the patch this edge belongs to
    label edgePatch = whichEdgePatch(eIndex);

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

        Info << "Edge: " << eIndex << ": " << edgeToCheck << endl;
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

    // Write lock the cell mutex
    cMutex_.lock(rwMutex::WRITE_LOCK);

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

    // Add two unlocked cell mutexes
    if (threader_->multiThreaded())
    {
        for (label i = 0; i < 2; i++)
        {
            cellMutex_.append();
        }
    }

    // Unlock the cell mutex from write lock
    cMutex_.unlock();

    // Write lock the face mutex
    fMutex_.lock(rwMutex::WRITE_LOCK);

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

    // Unlock the face mutex from write lock
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
    label otherPoint = -1, nextPoint = -1;
    FixedList<bool,2> foundEdge;

    // For a 2-2 swap on a boundary edge,
    // add two boundary faces and an edge
    FixedList<label,2> newBdyFaceIndex(-1);
    label newEdgeIndex = -1;

    if (edgePatch > -1)
    {
        // Temporary local variables
        label facePatch = -1;
        edge newEdge(-1, -1);
        FixedList<label,2> nBEdge(0);
        FixedList<FixedList<label,2>,2> bdyEdges;
        FixedList<face,2> newBdyTriFace(face(3));

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
                        bdyEdges[0][nBEdge[0]++] = fEdges[edgeI];

                        edgeStack(tIndex).push(fEdges[edgeI]);
                    }

                    if
                    (
                        edges_[fEdges[edgeI]]
                     == edge(edgeToCheck[1], otherPoint)
                    )
                    {
                        bdyFaceEdges[1][nBE[1]++] = fEdges[edgeI];
                        bdyEdges[1][nBEdge[1]++] = fEdges[edgeI];

                        edgeStack(tIndex).push(fEdges[edgeI]);
                    }
                }
            }
        }

        // Write lock the face mutex
        fMutex_.lock(rwMutex::WRITE_LOCK);

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

        // Add an edgeFaces entry
        labelList newBdyEdgeFaces(3, -1);
        newBdyEdgeFaces[0] = newBdyFaceIndex[0];
        newBdyEdgeFaces[1] = newFaceIndex;
        newBdyEdgeFaces[2] = newBdyFaceIndex[1];

        // Find the point other than the new edge
        // on the new triangular face
        findIsolatedPoint
        (
            newTriFace,
            newEdge,
            otherPoint,
            nextPoint
        );

        // Add an edgePoints entry
        labelList newBdyEdgePoints(3, -1);
        newBdyEdgePoints[0] = edgeToCheck[0];
        newBdyEdgePoints[1] = otherPoint;
        newBdyEdgePoints[2] = edgeToCheck[1];

        // Write lock the edge mutex
        eMutex_.lock(rwMutex::WRITE_LOCK);

        // Insert the edge
        newEdgeIndex = insertEdge
                       (
                           edgePatch,
                           newEdge,
                           newBdyEdgeFaces,
                           newBdyEdgePoints
                       );

        // Unlock the edge mutex from write lock
        eMutex_.unlock();

        // Update faceEdges with the new edge
        newFaceEdges[nE++] = newEdgeIndex;
        bdyFaceEdges[0][nBE[0]++] = newEdgeIndex;
        bdyFaceEdges[1][nBE[1]++] = newEdgeIndex;

        // Update edgeFaces and edgePoints with the two new faces
        forAll(bdyEdges[0], edgeI)
        {
            sizeUpList(newBdyFaceIndex[0], edgeFaces_[bdyEdges[0][edgeI]]);
            sizeUpList(newBdyFaceIndex[1], edgeFaces_[bdyEdges[1][edgeI]]);

            // Replace the edgePoints label, and preserve position on the list
            findIsolatedPoint
            (
                newBdyTriFace[0],
                edges_[bdyEdges[0][edgeI]],
                otherPoint,
                nextPoint
            );

            replaceLabel
            (
                edgeToCheck[1],
                otherPoint,
                edgePoints_[bdyEdges[0][edgeI]]
            );

            // Size up edgePoints again, so that it is sized down later
            sizeUpList(edgeToCheck[1], edgePoints_[bdyEdges[0][edgeI]]);

            // Replace the edgePoints label, and preserve position on the list
            findIsolatedPoint
            (
                newBdyTriFace[1],
                edges_[bdyEdges[1][edgeI]],
                otherPoint,
                nextPoint
            );

            replaceLabel
            (
                edgeToCheck[0],
                otherPoint,
                edgePoints_[bdyEdges[1][edgeI]]
            );

            // Size up edgePoints again, so that it is sized down later
            sizeUpList(edgeToCheck[0], edgePoints_[bdyEdges[1][edgeI]]);
        }

        // Add faceEdges for the two new boundary faces
        faceEdges_.append(bdyFaceEdges[0]);
        faceEdges_.append(bdyFaceEdges[1]);

        // Unlock the face mutex from write lock
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
                if (owner_[faceIndex] == cellIndex)
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

                // Update faceEdges, edgeFaces and edgePoints
                labelList& fEdges = faceEdges_[faceIndex];

                forAll(fEdges, edgeI)
                {
                    edge& checkEdge = edges_[fEdges[edgeI]];

                    if
                    (
                        (checkEdge == check[0])
                     || (checkEdge == check[1])
                     || (checkEdge == check[2])
                    )
                    {
                        newFaceEdges[nE++] = fEdges[edgeI];

                        sizeUpList
                        (
                            newFaceIndex,
                            edgeFaces_[fEdges[edgeI]]
                        );

                        // Find the isolated point and insert it
                        findIsolatedPoint
                        (
                            newTriFace,
                            checkEdge,
                            otherPoint,
                            nextPoint
                        );

                        insertLabel
                        (
                            otherPoint,
                            edgeToCheck[0],
                            edgeToCheck[1],
                            edgePoints_[fEdges[edgeI]]
                        );

                        edgeStack(tIndex).push(fEdges[edgeI]);

                        break;
                    }
                }
            }

            // Face is connected to edgeToCheck[1]
            if (!foundEdge[0] && foundEdge[1])
            {
                // Check if a face-flip is necessary
                if (owner_[faceIndex] == cellIndex)
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

    forAll(cellRemovalList, indexI)
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

    // Write lock the cell mutex
    cMutex_.lock(rwMutex::WRITE_LOCK);

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

    // Unlock the cell mutex from write lock
    cMutex_.unlock();

    // Write lock the face mutex
    fMutex_.lock(rwMutex::WRITE_LOCK);

    // Remove the faces and update associated edges
    forAll(facesForRemoval, faceI)
    {
        // Update edgeFaces and edgePoints
        face& checkFace = faces_[facesForRemoval[faceI]];
        labelList& fEdges = faceEdges_[facesForRemoval[faceI]];

        forAll(fEdges, edgeI)
        {
            label edgeIndex = fEdges[edgeI];

            if (edgeIndex != eIndex)
            {
                sizeDownList
                (
                    facesForRemoval[faceI],
                    edgeFaces_[edgeIndex]
                );

                // Find the isolated point and remove it
                findIsolatedPoint
                (
                    checkFace,
                    edges_[edgeIndex],
                    otherPoint,
                    nextPoint
                );

                sizeDownList
                (
                    otherPoint,
                    edgePoints_[edgeIndex]
                );

                edgeStack(tIndex).push(edgeIndex);
            }
        }

        // Now remove the face...
        removeFace(facesForRemoval[faceI]);
    }

    // Unlock the face mutex from write lock
    fMutex_.unlock();

    // Write lock the cell mutex
    cMutex_.lock(rwMutex::WRITE_LOCK);

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

    // Unlock the cell mutex from write lock
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

    // Allocate for the mapping information
    pointMap_.setSize(nPoints_, -1);

    label pointRenum = 0;

    addedPointRenumbering_.clear();

    HashList<point>::iterator ptIter = meshPoints_.begin();

    while (ptIter != meshPoints_.end())
    {
        // Obtain the index for this point
        label pIndex = ptIter.index();

        // Update the point info
        points[pointRenum] = ptIter();

        // Renumber the point index
        meshPoints_.reNumber(pointRenum, ptIter);

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
void dynamicTopoFvMesh::reOrderEdges
(
    edgeList& edges,
    labelListList& edgeFaces
)
{
    // *** Edge renumbering *** //
    // If edges were deleted during topology change, the numerical order ceases
    // to be continuous. Edges are added to respective internal/boundary patches

    // Allocate for mapping information
    edgeMap_.setSize(nEdges_, -1);

    label edgeInOrder = 0, allEdges = edges_.lastIndex() + 1;
    edgeList oldEdges(allEdges);
    labelListList oldEdgeFaces(allEdges);
    labelListList oldEdgePoints(allEdges);

    addedEdgeRenumbering_.clear();
    Map<label> addedEdgeReverseRenumbering;

    // Transfer old edge-based HashLists, and clear them
    HashList<edge>::iterator eIter = edges_.begin();
    HashList<labelList>::iterator efIter = edgeFaces_.begin();

    while (eIter != edges_.end())
    {
        oldEdges[eIter.index()] = eIter();
        oldEdgeFaces[efIter.index()].transfer(efIter());
        eIter++; efIter++;
    }

    edges_.clear(); edgeFaces_.clear();

    if (!twoDMesh_)
    {
        HashList<labelList>::iterator epIter = edgePoints_.begin();

        while (epIter != edgePoints_.end())
        {
            oldEdgePoints[epIter.index()].transfer(epIter());
            epIter++;
        }

        edgePoints_.clear();
    }

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

            // Renumber edgePoints
            if (!twoDMesh_)
            {
                labelList& thisEP = oldEdgePoints[edgeI];

                forAll(thisEP,pointI)
                {
                    if (thisEP[pointI] < nOldPoints_)
                    {
                        thisEP[pointI] = reversePointMap_[thisEP[pointI]];
                    }
                    else
                    {
                        thisEP[pointI] = addedPointRenumbering_[thisEP[pointI]];
                    }
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

                // Insert entities into mesh-reset lists
                edges[edgeInOrder] = thisEdge;
                edgeFaces[edgeInOrder].transfer(thisEF);

                if (!twoDMesh_)
                {
                    edgePoints_.append(oldEdgePoints[edgeI]);
                }

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

        // Insert entities into mesh-reset lists
        edges[edgeInOrder] = oldEdges[oldIndex];
        edgeFaces[edgeInOrder].transfer(oldEdgeFaces[oldIndex]);

        if (!twoDMesh_)
        {
            edgePoints_.append(oldEdgePoints[oldIndex]);
        }

        edgeInOrder++;
    }
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

    addedFaceRenumbering_.clear();
    Map<label> addedFaceReverseRenumbering;

    // Make a copy of the old face-based HashLists, and clear them
    HashList<face>::iterator fIter = faces_.begin();
    HashList<label>::iterator oIter = owner_.begin();
    HashList<label>::iterator nIter = neighbour_.begin();

    while(fIter != faces_.end())
    {
        oldFaces[fIter.index()].transfer(fIter());
        oldOwner[oIter.index()] = oIter();
        oldNeighbour[nIter.index()] = nIter();
        fIter++; oIter++; nIter++;
    }

    faces_.clear(); owner_.clear(); neighbour_.clear();

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

        forAll(curFaces, faceI)
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

            forAll(neiCells, ncI)
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

                // Cell-reordering may cause flipped faces.. Correct them.
                if (neighbourRenumber < ownerRenumber)
                {
                    faceRenumber = faceRenumber.reverseFace();
                }

                // Insert entities into HashLists...
                faces_.append(faceRenumber);
                owner_.append(cellI);
                neighbour_.append(minNei);

                // Insert entities into mesh-reset lists
                faces[faceInOrder].transfer(faceRenumber);
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

        // Insert entities into mesh-reset lists
        // NOTE: From OF-1.5 onwards, neighbour array
        //       does not store -1 for boundary faces
        faces[faceInOrder].transfer(oldFaces[oldIndex]);
        owner[faceInOrder] = ownerRenumber;

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
        oldCells[cIter.index()].transfer(cIter());
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

    forAll(cellCellAddr, cellI)
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
    forAll(visited, cellI)
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

                    forAll(neighbours, nI)
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

    if (debug)
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
    edgeList& edges,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    labelListList& edgeFaces
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
        Info << "Patch Starts [Face]: " << oldPatchStarts_ << endl;
        Info << "Patch Sizes  [Face]: " << oldPatchSizes_ << endl;
        Info << "Patch Starts [Edge]: " << oldEdgePatchStarts_ << endl;
        Info << "Patch Sizes  [Edge]: " << oldEdgePatchSizes_ << endl;
        Info << "=================" << endl;
        Info << "Mesh Info [n+1]:" << endl;
        Info << "Points: " << nPoints_ << endl;
        Info << "Edges: " << nEdges_ << endl;
        Info << "Faces: " << nFaces_ << endl;
        Info << "Cells: " << nCells_ << endl;
        Info << "Internal Edges: " << nInternalEdges_ << endl;
        Info << "Internal Faces: " << nInternalFaces_ << endl;
        Info << "Patch Starts [Face]: " << patchStarts_ << endl;
        Info << "Patch Sizes: [Face]: " << patchSizes_ << endl;
        Info << "Patch Starts [Edge]: " << edgePatchStarts_ << endl;
        Info << "Patch Sizes: [Edge]: " << edgePatchSizes_ << endl;
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
    if (debug) Info << "ReOrdering edges..." << endl;
    reOrderEdges(edges, edgeFaces);
}

// Copy edge-based connectivity from eMesh to HashLists
void dynamicTopoFvMesh::setEdgeConnectivity()
{
    if (!twoDMesh_)
    {
        pointEdges_ = eMeshPtr_->pointEdges();
    }

    faceEdges_ = eMeshPtr_->faceEdges();
}

// Check the state of connectivity HashLists
void dynamicTopoFvMesh::checkConnectivity()
{
    Info << "Checking edge-face connectivity...";

    label nFailedChecks = 0;
    label allEdges = edges_.lastIndex() + 1;
    labelList nEdgeFaces(allEdges, 0);

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
                Info << nl << nl << "Edge: " << faceEdges[edgeI]
                     << ": " << edgeToCheck << nl
                     << "was not found in face: " << feIter.index()
                     << ": " << faces_[feIter.index()] << nl
                     << "faceEdges: " << faceEdges
                     << endl;

                nFailedChecks++;

                WarningIn
                (
                    "dynamicTopoFvMesh::checkConnectivity()"
                )
                     << nl << "Edge-Face connectivity is inconsistent."
                     << endl;
            }
        }
    }

    forAllIter(HashList<labelList>::iterator, edgeFaces_, efIter)
    {
        labelList& edgeFaces = efIter();

        if (edgeFaces.size() != nEdgeFaces[efIter.index()])
        {
            Info << nl << nl << "Edge: " << efIter.index()
                 << "edgeFaces: " << edgeFaces << endl;

            nFailedChecks++;

            WarningIn
            (
                "dynamicTopoFvMesh::checkConnectivity()"
            )
                << nl << "Edge-Face connectivity is inconsistent."
                << endl;
        }

        // Check if this edge belongs to faceEdges for each face
        forAll(edgeFaces, faceI)
        {
            if
            (
                !foundInList
                (
                    efIter.index(), faceEdges_[edgeFaces[faceI]]
                )
            )
            {
                Info << nl << nl << "Edge: " << efIter.index()
                     << ", edgeFaces: " << edgeFaces << nl
                     << "was not found in faceEdges of face: "
                     << edgeFaces[faceI] << nl
                     << "faceEdges: " << faceEdges_[edgeFaces[faceI]]
                     << endl;

                nFailedChecks++;

                WarningIn
                (
                    "dynamicTopoFvMesh::checkConnectivity()"
                )
                    << nl << "Edge-Face connectivity is inconsistent."
                    << endl;
            }
        }
    }

    Info << "Done." << endl;

    if (!twoDMesh_)
    {
        Info << "Checking point-edge connectivity...";

        label allPoints = meshPoints_.lastIndex() + 1;
        List<labelHashSet> hlPointEdges(allPoints);

        forAllIter(HashList<edge>::iterator, edges_, eIter)
        {
            hlPointEdges[eIter()[0]].insert(eIter.index());
            hlPointEdges[eIter()[1]].insert(eIter.index());
        }

        forAllIter(HashList<labelList>::iterator, pointEdges_, peIter)
        {
            labelList& pointEdges = peIter();

            forAll(pointEdges, edgeI)
            {
                if (!hlPointEdges[peIter.index()].found(pointEdges[edgeI]))
                {
                    Info << nl << nl << "Point: " << peIter.index()
                         << "pointEdges: " << pointEdges << nl
                         << "hlPointEdges: " << hlPointEdges[peIter.index()]
                         << endl;

                    nFailedChecks++;

                    WarningIn
                    (
                        "dynamicTopoFvMesh::checkConnectivity()"
                    )
                        << nl << "Point-Edge connectivity is inconsistent."
                        << endl;
                }
            }

            // Do a size check as well
            if (hlPointEdges[peIter.index()].size() != pointEdges.size())
            {
                Info << nl << nl << "Point: " << peIter.index()
                     << "pointEdges: " << pointEdges << nl
                     << "hlPointEdges: " << hlPointEdges[peIter.index()]
                     << endl;

                nFailedChecks++;

                WarningIn
                (
                    "dynamicTopoFvMesh::checkConnectivity()"
                )
                    << "Size inconsistency."
                    << nl << "Point-Edge connectivity is inconsistent."
                    << endl;
            }
        }

        Info << "Done." << endl;

        Info << "Checking edge-points connectivity...";

        label otherPoint = -1, nextPoint = -1;

        forAllIter(HashList<labelList>::iterator, edgePoints_, epIter)
        {
            // Do a preliminary size check
            labelList& edgePoints = epIter();
            labelList& edgeFaces = edgeFaces_[epIter.index()];

            if (edgePoints.size() != edgeFaces.size())
            {
                Info << nl << nl
                     << "Edge: " << epIter.index()
                     << " " << edges_[epIter.index()] << endl;

                Info << "edgeFaces: " << edgeFaces << endl;
                forAll(edgeFaces, faceI)
                {
                    Info << edgeFaces[faceI] << ": "
                         << faces_[edgeFaces[faceI]]
                         << endl;
                }

                Info << "edgePoints: " << edgePoints << endl;

                nFailedChecks++;

                WarningIn
                (
                    "dynamicTopoFvMesh::checkConnectivity()"
                )
                    << nl << "Edge-Points connectivity is inconsistent."
                    << endl;
            }

            // Now check to see that both lists are consistent.
            edge& edgeToCheck = edges_[epIter.index()];

            forAll(edgeFaces, faceI)
            {
                face& faceToCheck = faces_[edgeFaces[faceI]];

                findIsolatedPoint
                (
                    faceToCheck,
                    edgeToCheck,
                    otherPoint,
                    nextPoint
                );

                if (!foundInList(otherPoint, edgePoints))
                {
                    Info << nl << nl
                         << "Edge: " << epIter.index()
                         << " " << edges_[epIter.index()] << endl;

                    Info << "edgeFaces: " << edgeFaces << endl;
                    forAll(edgeFaces, faceI)
                    {
                        Info << edgeFaces[faceI] << ": "
                             << faces_[edgeFaces[faceI]]
                             << endl;
                    }

                    Info << "edgePoints: " << edgePoints << endl;

                    nFailedChecks++;

                    WarningIn
                    (
                        "dynamicTopoFvMesh::checkConnectivity()"
                    )
                        << nl << "Edge-Points connectivity is inconsistent."
                        << endl;
                }
            }
        }

        Info << "Done." << endl;

        Info << "Checking cell-point connectivity...";

        // Loop through all cells and construct cell-to-node
        label cIndex = 0;
        label allCells = cells_.lastIndex() + 1;
        labelList cellIndex(allCells);
        List<labelHashSet> cellToNode(allCells);

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

        // Resize the lists
        cellIndex.setSize(cIndex);
        cellToNode.setSize(cIndex);

        // Preliminary check for size
        forAll(cellToNode, cellI)
        {
            if (cellToNode[cellI].size() != 4)
            {
                Info << nl << "Warning: Cell: "
                     << cellIndex[cellI] << " is inconsistent. "
                     << endl;

                cell& failedCell = cells_[cellIndex[cellI]];

                Info << "Cell faces: " << failedCell << endl;

                forAll(failedCell, faceI)
                {
                    Info << "\tFace: " << failedCell[faceI]
                         << " :: " << faces_[failedCell[faceI]]
                         << endl;

                    labelList& fEdges = faceEdges_[failedCell[faceI]];

                    forAll(fEdges, edgeI)
                    {
                        Info << "\t\tEdge: " << fEdges[edgeI]
                             << " :: " << edges_[fEdges[edgeI]]
                             << endl;
                    }
                }

                nFailedChecks++;
            }
        }

        Info << "Done." << endl;
    }

    if (nFailedChecks)
    {
        FatalErrorIn
        (
            "dynamicTopoFvMesh::checkConnectivity()"
        )
            << nFailedChecks << " failures were found in connectivity."
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

        // HashSet to keep track of cells in each level
        labelHashSet levelCells;

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
            forAll(toc,wordI)
            {
                word& pName = toc[wordI];

                if (bdy[patchI].name() == pName)
                {
                    label pStart = bdyPatch.start();

                    forAll(bdyPatch,faceI)
                    {
                        label ownCell = own[pStart+faceI];

                        if (cellLevels[ownCell] == 0)
                        {
                            cellLevels[ownCell] = level;
                            lengthScale[ownCell] =
                                fixedLengthScalePatches_[pName][0].scalarToken();

                            levelCells.insert(ownCell);

                            visitedCells++;
                        }
                    }

                    break;
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
                            label eIndex = getTriBoundaryEdge(pStart+faceI);
                            edge& e = edges_[eIndex];

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

                        levelCells.insert(ownCell);

                        visitedCells++;
                    }
                }
            }
        }

        // Perform multiple sweeps through the mesh...
        while (visitedCells < nCells())
        {
            // Loop through cells of the current level
            labelList currentLevelCells = levelCells.toc();
            levelCells.clear();

            // Loop through cells, and increment neighbour
            // cells of the current level
            forAll(currentLevelCells,cellI)
            {
                // Obtain the cells neighbouring this one
                const labelList& cList = cc[currentLevelCells[cellI]];

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

                        levelCells.insert(cList[indexI]);

                        visitedCells++;
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

// Given a boundary quad face, return a boundary triangular face.
// For 2D simplical meshes only.
label dynamicTopoFvMesh::getTriBoundaryFace
(
    const label faceIndex
)
{
    const labelList& fEdges = faceEdges_[faceIndex];

    forAll(fEdges, edgeI)
    {
        // Obtain edgeFaces for this edge
        const labelList& eFaces = edgeFaces_[fEdges[edgeI]];

        forAll(eFaces, faceI)
        {
            if (faces_[eFaces[faceI]].size() == 3)
            {
                // Found a triangular face. Return this face.
                return eFaces[faceI];
            }
        }
    }

    // This bit should never happen.
    FatalErrorIn
    (
        "label dynamicTopoFvMesh::getTriBoundaryFace()"
    )
        << "Cannot find a triangular face bordering face: "
        << faceIndex << " :: " << faces_[faceIndex]
        << abort(FatalError);

    return -1;
}

// Given a boundary quad face, pick out a boundary edge that
// contains a triangular face. For 2D simplical meshes only.
label dynamicTopoFvMesh::getTriBoundaryEdge
(
    const label faceIndex
)
{
    const labelList& fEdges = faceEdges_[faceIndex];

    forAll(fEdges, edgeI)
    {
        // Obtain edgeFaces for this edge
        const labelList& eFaces = edgeFaces_[fEdges[edgeI]];

        forAll(eFaces, faceI)
        {
            if (faces_[eFaces[faceI]].size() == 3)
            {
                // Found a triangular face. Return this edge.
                return fEdges[edgeI];
            }
        }
    }

    // This bit should never happen.
    FatalErrorIn
    (
        "label dynamicTopoFvMesh::getTriBoundaryEdge()"
    )
        << "Cannot find a triangular face bordering face: "
        << faceIndex << " :: " << faces_[faceIndex]
        << abort(FatalError);

    return -1;
}

// 2D Edge-swapping engine
void dynamicTopoFvMesh::swap2DEdges(void *argument)
{
    // Recast the argument
    topoMeshStruct *thread = reinterpret_cast<topoMeshStruct*>(argument);
    dynamicTopoFvMesh *mesh = thread->mesh_;

    // Figure out which thread this is...
    label tIndex = mesh->self();

    bool failed = false;

    // Pick items off the stack
    while (!mesh->faceStack(tIndex).empty())
    {
        // Retrieve the index for this face
        label fIndex = mesh->faceStack(tIndex).pop();

        // Perform a Delaunay test and check if a flip is necesary.
        mesh->testDelaunay
        (
            fIndex,
            failed
        );

        if (failed)
        {
            // Swap this face.
            mesh->swapQuadFace(fIndex);
        }
    }
}

// Initialize edge related connectivity lists
void dynamicTopoFvMesh::initEdges()
{
    // Initialize eMesh, and copy to HashLists
    eMeshPtr_.set(new eMesh(*this));

    // Obtain information
    nEdges_ = eMeshPtr_->nEdges();
    nInternalEdges_ = eMeshPtr_->nInternalEdges();
    edgePatchSizes_ = eMeshPtr_->boundary().patchSizes();
    edgePatchStarts_ = eMeshPtr_->boundary().patchStarts();

    // Set HashLists with edge connectivity information
    edges_ = eMeshPtr_->edges();
    edgeFaces_ = eMeshPtr_->edgeFaces();
    faceEdges_ = eMeshPtr_->faceEdges();

    if (!twoDMesh_)
    {
        pointEdges_ = eMeshPtr_->pointEdges();
        edgePoints_ = eMeshPtr_->edgePoints();
    }
}

// Initialize mutex lists
void dynamicTopoFvMesh::initMutexLists()
{
    if (threader_->multiThreaded())
    {
        pointMutex_.setSize(nPoints_);
        edgeMutex_.setSize(nEdges_);
        faceMutex_.setSize(nFaces_);
        cellMutex_.setSize(nCells_);
    }
}

// Unlock mutex lists
inline void dynamicTopoFvMesh::unlockMutexLists(const label tIndex)
{
    if (threader_->multiThreaded())
    {
        forAll(cLocks_[tIndex], cellI)
        {
            cellMutex_[cLocks_[tIndex][cellI]].unlock();
        }
        cLocks_[tIndex].clear();

        forAll(fLocks_[tIndex], faceI)
        {
            faceMutex_[fLocks_[tIndex][faceI]].unlock();
        }
        fLocks_[tIndex].clear();

        forAll(eLocks_[tIndex], edgeI)
        {
            edgeMutex_[eLocks_[tIndex][edgeI]].unlock();
        }
        eLocks_[tIndex].clear();

        forAll(pLocks_[tIndex], pointI)
        {
            pointMutex_[pLocks_[tIndex][pointI]].unlock();
        }
        pLocks_[tIndex].clear();
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

//- Clear out the stack
inline void dynamicTopoFvMesh::stack::clear()
{
    stackMutex_.lock();

    stack_.clear();

    stackMutex_.unlock();
}

//- Print out the stack
inline void dynamicTopoFvMesh::stack::print()
{
    Info << stack_ << endl;
}

// Return length-scale at an face-location in the mesh [2D]
inline scalar dynamicTopoFvMesh::meshFaceLengthScale
(
    const label fIndex
)
{
    scalar scale = 0.0;

    // Determine whether the face is internal
    if (whichPatch(fIndex) < 0)
    {
        scale =
            0.5*
            (
                lengthScale_[owner_[fIndex]]
              + lengthScale_[neighbour_[fIndex]]
            );
    }
    else
    {
        scale = boundaryLengthScale(fIndex);
    }

    return scale;
}

// Compute length-scale at an edge-location in the mesh [3D]
inline bool dynamicTopoFvMesh::meshEdgeLengthScale
(
    const label eIndex,
    scalar& scale
)
{
    // Figure out which thread this is...
    label tIndex = self();

    // Try to read-lock this edge.
    if (tryEdgeLock(eIndex))
    {
        edgeStack(tIndex).push(eIndex);
        return true;
    }

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

    // Undo all read-locks
    unlockMutexLists(tIndex);

    // Return a successful lock
    return false;
}

// Check if a given face is a quad
inline bool dynamicTopoFvMesh::checkQuadFace(const label fIndex)
{
    return (faces_[fIndex].size() == 4);
}

// Compute the length of an edge
inline bool dynamicTopoFvMesh::edgeLength
(
    const label eIndex,
    scalar& length
)
{
    // Figure out which thread this is...
    label tIndex = self();

    // Try to read-lock this edge.
    // Calling function is expected to push it back on stack if it fails.
    if (tryEdgeLock(eIndex))
    {
        return true;
    }

    // Try to read-lock the two points of this edge.
    if (tryEdgePointLock(eIndex))
    {
        return true;
    }

    edge& thisEdge = edges_[eIndex];

    length = mag(meshPoints_[thisEdge[1]] - meshPoints_[thisEdge[0]]);

    // Undo all read-locks
    unlockMutexLists(tIndex);

    // Return a successful lock
    return false;
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
    // Bisect when boundary edge-length > ratioMax_*lengthScale
    // Collapse when boundary edge-length < ratioMin_*lengthScale

    // Recast the argument
    topoMeshStruct *thread = reinterpret_cast<topoMeshStruct*>(argument);
    dynamicTopoFvMesh *mesh = thread->mesh_;

    // Figure out which thread this is...
    label tIndex = mesh->self();
    scalar length = 0.0, scale = 0.0;

    // Pick items off the stack
    while (!mesh->faceStack(tIndex).empty())
    {
        // Retrieve the index for this face
        label fIndex = mesh->faceStack(tIndex).pop();

        // Select only quad-faces
        if (mesh->checkQuadFace(fIndex))
        {
            // Measure the boundary edge-length of the face in question
            if
            (
                mesh->edgeLength
                (
                    mesh->getTriBoundaryEdge(fIndex),
                    length
                )
            )
            {
                continue;
            }

            // Determine the length-scale at this face
            mesh->meshFaceLengthScale(fIndex);

            // Check if this boundary face is adjacent to a sliver-cell,
            // and remove it by a two-step bisection/collapse operation.
            if (mesh->whichPatch(fIndex) != -1)
            {
                mesh->remove2DSliver(fIndex);
            }

            if (length > mesh->ratioMax()*scale)
            {
                // Bisect this face
                mesh->bisectQuadFace(fIndex);
            }
            else
            if (length < mesh->ratioMin()*scale)
            {
                // Collapse this face
                mesh->collapseQuadFace(fIndex);
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

    // Figure out which thread this is...
    label tIndex = mesh->self();

    // Hull variables
    label m = -1;
    scalar minQuality;

    // Dynamic programming variables
    scalarListList Q;
    labelListList K, triangulations;

    // Allocate dynamic programming tables
    mesh->initTables(Q, K, triangulations);

    // Pick edges off the stack
    while (!mesh->edgeStack(tIndex).empty())
    {
        // Retrieve an edge from the stack
        label eIndex = mesh->edgeStack(tIndex).pop();

        // Check if this edge is on a bounding curve
        if (mesh->checkBoundingCurve(eIndex))
        {
            continue;
        }

        // Compute the minimum quality of cells around this edge
        if
        (
            mesh->computeMinQuality
            (
                eIndex,
                minQuality
            )
        )
        {
            continue;
        }

        // Fill the dynamic programming tables
        if (mesh->fillTables(eIndex, minQuality, m, Q, K, triangulations))
        {
            // Check if edge-swapping is required.
            if (Q[0][m-1] > minQuality)
            {
                // Remove this edge according to the swap sequence
                mesh->removeEdgeFlips(eIndex, K, triangulations);
            }
        }
    }
}

// 3D Edge-bisection/collapse engine
void dynamicTopoFvMesh::edgeBisectCollapse3D
(
    void *argument
)
{
    // Loop through all edges and bisect/collapse by the criterion:
    // Bisect when edge-length > ratioMax_*lengthScale
    // Collapse when edge-length < ratioMin_*lengthScale

    // Recast the argument
    topoMeshStruct *thread = reinterpret_cast<topoMeshStruct*>(argument);
    dynamicTopoFvMesh *mesh = thread->mesh_;

    // Figure out which thread this is...
    label tIndex = mesh->self();
    scalar length = 0.0, scale = 0.0;

    while (!mesh->edgeStack(tIndex).empty())
    {
        // Retrieve an edge from the stack
        label eIndex = mesh->edgeStack(tIndex).pop();

        // Measure the edge-length
        if (mesh->edgeLength(eIndex, length))
        {
            mesh->edgeStack(tIndex).push(eIndex);
            continue;
        }

        // Determine the length-scale at this point in the mesh
        if (mesh->meshEdgeLengthScale(eIndex, scale))
        {
            continue;
        }

        if (length > mesh->ratioMax()*scale)
        {
            // Bisect this edge
            mesh->bisectEdge(eIndex);
        }
        else
        if (length < mesh->ratioMin()*scale)
        {
            // Collapse this edge
            mesh->collapseEdge(eIndex);
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

// Try to lock a point.
bool dynamicTopoFvMesh::tryPointLock
(
    const label pIndex,
    const rwMutex::lockType lType
)
{
    if (threader_->multiThreaded())
    {
        // Figure out which thread this is...
        label tIndex = self();

        if (obtainLock(tIndex, lType, pointMutex_[pIndex]))
        {
            pLocks_[tIndex].append(pIndex);
        }
        else
        {
            unlockMutexLists(tIndex);
            return true;
        }
    }

    // Return a successful lock
    return false;
}

// Try to lock a face.
// Returns false on success.
bool dynamicTopoFvMesh::tryFaceLock
(
    const label fIndex,
    const rwMutex::lockType lType
)
{
    if (threader_->multiThreaded())
    {
        // Figure out which thread this is...
        label tIndex = self();

        if (obtainLock(tIndex, lType, faceMutex_[fIndex]))
        {
            fLocks_[tIndex].append(fIndex);
        }
        else
        {
            unlockMutexLists(tIndex);
            return true;
        }
    }

    // Return a successful lock
    return false;
}

// Try to lock all points of a face.
// Assumes face is already locked.
// Returns false on success.
bool dynamicTopoFvMesh::tryFacePointLock
(
    const label fIndex,
    const rwMutex::lockType lType
)
{
    if (threader_->multiThreaded())
    {
        // Figure out which thread this is...
        label tIndex = self();

        // Try to lock all points of this face
        face& lockFace = faces_[fIndex];

        forAll(lockFace, pointI)
        {
            if (obtainLock(tIndex, lType, pointMutex_[lockFace[pointI]]))
            {
                pLocks_[tIndex].append(lockFace[pointI]);
            }
            else
            {
                unlockMutexLists(tIndex);
                return true;
            }
        }
    }

    // Return a successful lock
    return false;
}

// Try to lock all edges of a face.
// Assumes face is already locked.
// Returns false on success.
bool dynamicTopoFvMesh::tryFaceEdgeLock
(
    const label fIndex,
    const rwMutex::lockType lType,
    const label edgeToAvoid
)
{
    if (threader_->multiThreaded())
    {
        // Figure out which thread this is...
        label tIndex = self();

        // Try to lock all edges of this face
        labelList& fEdges = faceEdges_[fIndex];

        forAll(fEdges, edgeI)
        {
            if (edgeToAvoid != -1)
            {
                if (obtainLock(tIndex, lType, edgeMutex_[fEdges[edgeI]]))
                {
                    eLocks_[tIndex].append(fEdges[edgeI]);
                }
                else
                {
                    unlockMutexLists(tIndex);
                    return true;
                }
            }
        }
    }

    // Return a successful lock
    return false;
}

// Try to lock an edge.
// Returns false on success.
bool dynamicTopoFvMesh::tryEdgeLock
(
    const label eIndex,
    const rwMutex::lockType lType
)
{
    if (threader_->multiThreaded())
    {
        // Figure out which thread this is...
        label tIndex = self();

        if (obtainLock(tIndex, lType, edgeMutex_[eIndex]))
        {
            eLocks_[tIndex].append(eIndex);
        }
        else
        {
            unlockMutexLists(tIndex);
            return true;
        }
    }

    // Return a successful lock
    return false;
}

// Try to lock the two points of an edge.
// Assumes edge is already locked.
// Returns false on success.
bool dynamicTopoFvMesh::tryEdgePointLock
(
    const label eIndex,
    const rwMutex::lockType lType
)
{
    if (threader_->multiThreaded())
    {
        // Figure out which thread this is...
        label tIndex = self();

        // Obtain a reference to this edge...
        edge& lockEdge = edges_[eIndex];

        if (obtainLock(tIndex, lType, pointMutex_[lockEdge[0]]))
        {
            pLocks_[tIndex].append(lockEdge[0]);
        }
        else
        {
            unlockMutexLists(tIndex);
            return true;
        }

        if (obtainLock(tIndex, lType, pointMutex_[lockEdge[1]]))
        {
            pLocks_[tIndex].append(lockEdge[1]);
        }
        else
        {
            unlockMutexLists(tIndex);
            return true;
        }
    }

    // Return a successful lock
    return false;
}

// Try to lock hull-points around an edge
// Assumes edge is already locked.
// Returns false on success.
bool dynamicTopoFvMesh::tryEdgeHullLock
(
    const label eIndex,
    const rwMutex::lockType lType
)
{
    if (threader_->multiThreaded())
    {
        // Figure out which thread this is...
        label tIndex = self();

        // Try to lock all hull points of this edge
        labelList& ePoints = edgePoints_[eIndex];

        forAll(ePoints, pointI)
        {
            if (obtainLock(tIndex, lType, pointMutex_[ePoints[pointI]]))
            {
                pLocks_[tIndex].append(ePoints[pointI]);
            }
            else
            {
                unlockMutexLists(tIndex);
                return true;
            }
        }
    }

    // Return a successful lock
    return false;
}

// Try to lock the faces around an edge.
// Assumes edge is already locked.
// Returns false on success.
bool dynamicTopoFvMesh::tryEdgeFaceLock
(
    const label eIndex,
    const rwMutex::lockType lType
)
{
    if (threader_->multiThreaded())
    {
        // Figure out which thread this is...
        label tIndex = self();

        // Try to lock all faces of this edge
        labelList& eFaces = edgeFaces_[eIndex];

        forAll(eFaces, faceI)
        {
            if (obtainLock(tIndex, lType, faceMutex_[eFaces[faceI]]))
            {
                fLocks_[tIndex].append(eFaces[faceI]);
            }
            else
            {
                unlockMutexLists(tIndex);
                return true;
            }
        }
    }

    // Return a successful lock
    return false;
}

// Remove an index from the point lock list
void dynamicTopoFvMesh::removePointLock
(
    const label pIndex
)
{
    if (threader_->multiThreaded())
    {
        // Figure out which thread this is...
        label tIndex = self();

        if (foundInList(pIndex, pLocks_[tIndex]))
        {
            sizeDownList(pIndex, pLocks_[tIndex]);
        }
    }
}

// Remove an index from the edge lock list
void dynamicTopoFvMesh::removeEdgeLock
(
    const label eIndex
)
{
    if (threader_->multiThreaded())
    {
        // Figure out which thread this is...
        label tIndex = self();

        if (foundInList(eIndex, eLocks_[tIndex]))
        {
            sizeDownList(eIndex, eLocks_[tIndex]);
        }
    }
}

// Remove an index from the face lock list
void dynamicTopoFvMesh::removeFaceLock
(
    const label fIndex
)
{
    if (threader_->multiThreaded())
    {
        // Figure out which thread this is...
        label tIndex = self();

        if (foundInList(fIndex, fLocks_[tIndex]))
        {
            sizeDownList(fIndex, fLocks_[tIndex]);
        }
    }
}

// Remove an index from the cell lock list
void dynamicTopoFvMesh::removeCellLock
(
    const label cIndex
)
{
    if (threader_->multiThreaded())
    {
        // Figure out which thread this is...
        label tIndex = self();

        if (foundInList(cIndex, cLocks_[tIndex]))
        {
            sizeDownList(cIndex, cLocks_[tIndex]);
        }
    }
}

// Return the integer threadID for a given pthread
// Return zero for single-threaded operation
inline label dynamicTopoFvMesh::self()
{
    if (threader_->multiThreaded())
    {
        for (label i = 0; i < threader_->getNumThreads(); i++)
        {
            if (pthread_equal(structPtr_[i].pthreadID_, pthread_self()))
            {
                return i;
            }
        }

        // This bit should never happen.
        FatalErrorIn
        (
            "label dynamicTopoFvMesh::self()"
        )
            << "Cannot find a corresponding thread ID." << endl
            << abort(FatalError);
    }

    return 0;
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
    const label fIndex
)
{
    face f;
    bool found = false;
    label commonIndex = -1;
    FixedList<label,4> otherPointIndex(-1), nextToOtherPoint(-1);
    FixedList<label,2> c0BdyIndex, c0IntIndex, c1BdyIndex, c1IntIndex;
    FixedList<face,2>  c0BdyFace,  c0IntFace,  c1BdyFace,  c1IntFace;
    FixedList<label,2> commonEdgeIndex(-1);
    FixedList<edge,2>  commonEdges;
    FixedList<label,4> otherEdgeIndex(-1);
    FixedList<label,4> commonFaceIndex(-1), modifiedEdgeIndex(-1);
    FixedList<face,4>  commonFaces(face(3)), commonIntFaces(face(4));
    FixedList<label,4> commonIntFaceIndex(-1);
    FixedList<bool,2> foundTriFace0(false), foundTriFace1(false);
    FixedList<face,2> triFaces0(face(3)), triFaces1(face(3));

    // Get the two cells on either side...
    label c0 = owner_[fIndex];
    label c1 = neighbour_[fIndex];

    // Get cell references
    cell& cell_0 = cells_[c0];
    cell& cell_1 = cells_[c1];

    // Need to find common faces and edges...
    // At the end of this loop, commonFaces [0] & [1] share commonEdge[0]
    // and commonFaces [2] & [3] share commonEdge[1]
    // Also, commonFaces[0] & [2] are connected to cell[0],
    // and commonFaces[1] & [3] are connected to cell[1]
    labelList& fEdges = faceEdges_[fIndex];

    forAll(fEdges, edgeI)
    {
        // Break out if all triangular faces are found
        if
        (
            foundTriFace0[0] && foundTriFace0[1]
         && foundTriFace1[0] && foundTriFace1[1]
        )
        {
            break;
        }

        // Obtain edgeFaces for this edge
        labelList& eFaces = edgeFaces_[fEdges[edgeI]];

        forAll(eFaces, faceI)
        {
            face& eFace = faces_[eFaces[faceI]];

            if (eFace.size() == 3)
            {
                // Found a triangular face. Determine which cell it belongs to.
                if (owner_[eFaces[faceI]] == c0)
                {
                    if (foundTriFace0[0])
                    {
                        // Update the second face on cell[0].
                        commonIndex = 2;
                        foundTriFace0[1] = true;

                        if (foundTriFace1[1])
                        {
                            commonEdgeIndex[1] = fEdges[edgeI];
                            commonEdges[1] = edges_[fEdges[edgeI]];
                        }
                    }
                    else
                    {
                        // Update the first face on cell[0].
                        commonIndex = 0;
                        foundTriFace0[0] = true;

                        if (foundTriFace1[0])
                        {
                            commonEdgeIndex[0] = fEdges[edgeI];
                            commonEdges[0] = edges_[fEdges[edgeI]];
                        }
                    }
                }
                else
                {
                    if (foundTriFace1[0])
                    {
                        // Update the second face on cell[1].
                        commonIndex = 3;
                        foundTriFace1[1] = true;

                        if (foundTriFace0[1])
                        {
                            commonEdgeIndex[1] = fEdges[edgeI];
                            commonEdges[1] = edges_[fEdges[edgeI]];
                        }
                    }
                    else
                    {
                        // Update the first face on cell[1].
                        commonIndex = 1;
                        foundTriFace1[0] = true;

                        if (foundTriFace0[0])
                        {
                            commonEdgeIndex[0] = fEdges[edgeI];
                            commonEdges[0] = edges_[fEdges[edgeI]];
                        }
                    }
                }

                // Store the face and index
                commonFaces[commonIndex][0] = eFace[0];
                commonFaces[commonIndex][1] = eFace[1];
                commonFaces[commonIndex][2] = eFace[2];

                commonFaceIndex[commonIndex] = eFaces[faceI];
            }
        }
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        Info << nl << nl << "Face: " << fIndex
             << " needs to be flipped. " << endl;

        Info << "Cell[0]: " << c0 << ": " << cell_0 << endl;
        Info << "Cell[1]: " << c1 << ": " << cell_1 << endl;

        Info << "Common Faces: Set 1: "
             << commonFaceIndex[0] << ": " << commonFaces[0] << ", "
             << commonFaceIndex[1] << ": " << commonFaces[1] << endl;

        Info << "Common Faces: Set 2: "
             << commonFaceIndex[2] << ": " << commonFaces[2] << ", "
             << commonFaceIndex[3] << ": " << commonFaces[3] << endl;

        Info << "Old face: " << faces_[fIndex] << endl;
    }
#   endif

    // Find the interior/boundary faces.
    findPrismFaces
    (
        fIndex,
        c0,
        c0BdyFace,
        c0BdyIndex,
        c0IntFace,
        c0IntIndex
    );

    findPrismFaces
    (
        fIndex,
        c1,
        c1BdyFace,
        c1BdyIndex,
        c1IntFace,
        c1IntIndex
    );

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
    forAll(fEdges, edgeI)
    {
        if
        (
            fEdges[edgeI] != commonEdgeIndex[0]
         && fEdges[edgeI] != commonEdgeIndex[1]
        )
        {
            // Obtain a reference to this edge
            edge& eThis = edges_[fEdges[edgeI]];

            if
            (
                eThis[0] == nextToOtherPoint[0]
             || eThis[1] == nextToOtherPoint[0]
            )
            {
                otherEdgeIndex[0] = fEdges[edgeI];
            }
            else
            {
                otherEdgeIndex[1] = fEdges[edgeI];
            }
        }
    }

    // At the end of this loop, commonIntFaces [0] & [1] share otherEdges[0]
    // and commonIntFaces [2] & [3] share the otherEdges[1],
    // where [0],[2] lie on cell[0] and [1],[3] lie on cell[1]
    found = false;

    labelList& e1 = faceEdges_[c0IntIndex[0]];

    forAll(e1,edgeI)
    {
        if (e1[edgeI] == otherEdgeIndex[0])
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
        // The edge was not found before
        commonIntFaces[0] = c0IntFace[1];
        commonIntFaces[2] = c0IntFace[0];
        commonIntFaceIndex[0] = c0IntIndex[1];
        commonIntFaceIndex[2] = c0IntIndex[0];
    }

    found = false;

    labelList& e3 = faceEdges_[c1IntIndex[0]];

    forAll(e3,edgeI)
    {
        if (e3[edgeI] == otherEdgeIndex[0])
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
        // The edge was not found before
        commonIntFaces[1] = c1IntFace[1];
        commonIntFaces[3] = c1IntFace[0];
        commonIntFaceIndex[1] = c1IntIndex[1];
        commonIntFaceIndex[3] = c1IntIndex[0];
    }

    // Find two common edges between interior/interior faces
    findCommonEdge
    (
        c0IntIndex[0],
        c0IntIndex[1],
        otherEdgeIndex[2]
    );

    findCommonEdge
    (
        c1IntIndex[0],
        c1IntIndex[1],
        otherEdgeIndex[3]
    );

    // Find four common edges between interior/boundary faces
    findCommonEdge
    (
        commonFaceIndex[1],
        commonIntFaceIndex[1],
        modifiedEdgeIndex[0]
    );

    findCommonEdge
    (
        commonFaceIndex[3],
        commonIntFaceIndex[1],
        modifiedEdgeIndex[1]
    );

    findCommonEdge
    (
        commonFaceIndex[0],
        commonIntFaceIndex[2],
        modifiedEdgeIndex[2]
    );

    findCommonEdge
    (
        commonFaceIndex[2],
        commonIntFaceIndex[2],
        modifiedEdgeIndex[3]
    );

    // Modify the five faces belonging to this hull
    face& newFace = faces_[fIndex];
    face& newBdyFace0 = faces_[commonFaceIndex[0]];
    face& newBdyFace1 = faces_[commonFaceIndex[1]];
    face& newBdyFace2 = faces_[commonFaceIndex[2]];
    face& newBdyFace3 = faces_[commonFaceIndex[3]];
    labelList& newFEdges = faceEdges_[fIndex];
    label c0count=0, c1count=0;

    // Size down edgeFaces for the original face
    sizeDownList
    (
        fIndex,
        edgeFaces_[otherEdgeIndex[0]]
    );

    sizeDownList
    (
        fIndex,
        edgeFaces_[otherEdgeIndex[1]]
    );

    // Size up edgeFaces for the face after flipping
    sizeUpList
    (
        fIndex,
        edgeFaces_[otherEdgeIndex[2]]
    );

    sizeUpList
    (
        fIndex,
        edgeFaces_[otherEdgeIndex[3]]
    );

    // Replace edgeFaces and faceEdges
    replaceLabel
    (
        modifiedEdgeIndex[0],
        modifiedEdgeIndex[1],
        faceEdges_[commonFaceIndex[1]]
    );

    replaceLabel
    (
        commonFaceIndex[1],
        commonFaceIndex[0],
        edgeFaces_[modifiedEdgeIndex[0]]
    );

    replaceLabel
    (
        modifiedEdgeIndex[1],
        modifiedEdgeIndex[0],
        faceEdges_[commonFaceIndex[3]]
    );

    replaceLabel
    (
        commonFaceIndex[3],
        commonFaceIndex[2],
        edgeFaces_[modifiedEdgeIndex[1]]
    );

    replaceLabel
    (
        modifiedEdgeIndex[2],
        modifiedEdgeIndex[3],
        faceEdges_[commonFaceIndex[0]]
    );

    replaceLabel
    (
        commonFaceIndex[0],
        commonFaceIndex[1],
        edgeFaces_[modifiedEdgeIndex[2]]
    );

    replaceLabel
    (
        modifiedEdgeIndex[3],
        modifiedEdgeIndex[2],
        faceEdges_[commonFaceIndex[2]]
    );

    replaceLabel
    (
        commonFaceIndex[2],
        commonFaceIndex[3],
        edgeFaces_[modifiedEdgeIndex[3]]
    );

    // Define parameters for the new flipped face
    newFace[0] = otherPointIndex[0];
    newFace[1] = otherPointIndex[1];
    newFace[2] = otherPointIndex[3];
    newFace[3] = otherPointIndex[2];
    newFEdges[0] = otherEdgeIndex[2];
    newFEdges[1] = commonEdgeIndex[0];
    newFEdges[2] = otherEdgeIndex[3];
    newFEdges[3] = commonEdgeIndex[1];
    cell_0[c0count++] = fIndex;
    cell_1[c1count++] = fIndex;
    owner_[fIndex] = c0;
    neighbour_[fIndex] = c1;

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

#   ifdef FULLDEBUG
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
#   endif

    // Check the orientation of the two quad faces, and modify as necessary
    label newOwn=0, newNei=0;

    // The quad face belonging to cell[1] now becomes a part of cell[0]
    if (neighbour_[commonIntFaceIndex[1]] == -1)
    {
        // Boundary face
        // Face doesn't need to be flipped, just update the owner
        f          = commonIntFaces[1];
        newOwn     = c0;
        newNei     = -1;
    }
    else
    if (owner_[commonIntFaceIndex[1]] == c1)
    {
        // This face is on the interior, check for previous owner
        // Upper-triangular ordering has to be maintained, however...
        if (c0 > neighbour_[commonIntFaceIndex[1]])
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
    if (neighbour_[commonIntFaceIndex[1]] == c1)
    {
        // This face is on the interior, check for previous neighbour
        // Upper-triangular ordering has to be maintained, however...
        if (c0 < owner_[commonIntFaceIndex[1]])
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
    if (neighbour_[commonIntFaceIndex[2]] == -1)
    {
        // Boundary face
        // Face doesn't need to be flipped, just update the owner
        f          = commonIntFaces[2];
        newOwn     = c1;
        newNei     = -1;
    }
    else
    if (owner_[commonIntFaceIndex[2]] == c0)
    {
        // This face is on the interior, check for previous owner
        // Upper-triangular ordering has to be maintained, however...
        if (c1 > neighbour_[commonIntFaceIndex[2]])
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
    if (neighbour_[commonIntFaceIndex[2]] == c0)
    {
        // This face is on the interior, check for previous neighbour
        // Upper-triangular ordering has to be maintained, however...
        if (c1 < owner_[commonIntFaceIndex[2]])
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

    // Set the flag
    topoChangeFlag_ = true;
}

// Method for the bisection of a quad-face in 2D
void dynamicTopoFvMesh::bisectQuadFace
(
    const label fIndex
)
{
    // Quad-face bisection performs the following operations:
    //      [1] Add two points at the middle of the face
    //      [2] Create a new internal face for each bisected cell
    //      [3] Modify existing face and create a new half-face
    //      [4] Modify triangular boundary faces, and create new ones as well
    //      [5] Create edges for new faces
    //      Update faceEdges and edgeFaces information

    bool found;
    label replaceFace = -1;
    face tmpQuadFace(4), tmpTriFace(3);
    FixedList<label,7> newFaceIndex(-1), newEdgeIndex(-1);
    FixedList<edge,2> commonEdges;
    FixedList<label,4> modifiedEdgeIndex(-1);
    FixedList<label,2> commonEdgeIndex(-1), commonFaceIndex(-1);
    FixedList<label,2> newPtIndex(-1), newCellIndex(-1), otherEdgePoint(-1);
    FixedList<label,4> otherEdgeIndex(-1);
    FixedList<label,4> otherPointIndex(-1), nextToOtherPoint(-1);
    FixedList<label,2> c0BdyIndex, c0IntIndex, c1BdyIndex, c1IntIndex;
    FixedList<face,2>  c0BdyFace,  c0IntFace,  c1BdyFace,  c1IntFace;

    // Obtain a reference for this face...
    face& thisFace = faces_[fIndex];

    // Get the two cells on either side...
    label c0 = owner_[fIndex], c1 = neighbour_[fIndex];

    // Find the prism faces for cell[0].
    cell &cell_0 = cells_[c0];

    findPrismFaces
    (
        fIndex,
        c0,
        c0BdyFace,
        c0BdyIndex,
        c0IntFace,
        c0IntIndex
    );

#   ifdef FULLDEBUG
    if (debug)
    {
        Info << nl << nl << "Face: " << fIndex
             << ": " << thisFace << " is to be bisected. " << endl;

        Info << "Cell[0]: " << c0 << ": " << cell_0 << endl;

        forAll(cell_0, faceI)
        {
            Info << cell_0[faceI] << ": " << faces_[cell_0[faceI]] << endl;
        }
    }
#   endif

    // Find the common-edge between the triangular boundary faces
    // and the face under consideration.
    findCommonEdge(c0BdyIndex[0], fIndex, commonEdgeIndex[0]);
    findCommonEdge(c0BdyIndex[1], fIndex, commonEdgeIndex[1]);

    commonEdges[0] = edges_[commonEdgeIndex[0]];
    commonEdges[1] = edges_[commonEdgeIndex[1]];

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

    // For convenience...
    otherEdgePoint[0] = commonEdges[0].otherVertex(nextToOtherPoint[0]);
    otherEdgePoint[1] = commonEdges[1].otherVertex(nextToOtherPoint[1]);

    // Add two new points to the end of the list
    newPtIndex[0] = meshPoints_.append
                    (
                        0.5*
                        (
                            meshPoints_[commonEdges[0][0]]
                          + meshPoints_[commonEdges[0][1]]
                        )
                    );

    newPtIndex[1] = meshPoints_.append
                    (
                        0.5*
                        (
                            meshPoints_[commonEdges[1][0]]
                          + meshPoints_[commonEdges[1][1]]
                        )
                    );
    nPoints_ += 2;

    // Add a new prism cell to the end of the list
    newCellIndex[0] = cells_.append(cell(5));
    cell &newCell0 = cells_[newCellIndex[0]];

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
    cellParents_.insert(newCellIndex[0],firstParent);

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
        newCellIndex[0],
        objectMap
        (
            newCellIndex[0],
            c0MasterObjects.toc()
        )
    );

    // Add a new element to the lengthScale field
    // (Currently the same as cell[0])
    lengthScale_.append(lengthScale_[c0]);

    // Modify the two existing triangle boundary faces

    // Zeroth boundary face - Owner = c[0] & Neighbour [-1] (unchanged)
    replaceLabel
    (
        otherEdgePoint[0],
        newPtIndex[0],
        c0BdyFace[0]
    );

    // First boundary face - Owner = newCell[0], Neighbour = -1
    replaceLabel
    (
        otherEdgePoint[1],
        newPtIndex[1],
        c0BdyFace[1]
    );

    owner_[c0BdyIndex[1]] = newCellIndex[0];
    replaceLabel(c0BdyIndex[1],-1,cell_0);

    // Detect edges other than commonEdges
    labelList& fEdges = faceEdges_[fIndex];

    forAll(fEdges, edgeI)
    {
        if
        (
            fEdges[edgeI] != commonEdgeIndex[0]
         && fEdges[edgeI] != commonEdgeIndex[1]
        )
        {
            edge& eThis = edges_[fEdges[edgeI]];

            if
            (
                eThis[0] == nextToOtherPoint[0]
             || eThis[1] == nextToOtherPoint[0]
            )
            {
                otherEdgeIndex[0] = fEdges[edgeI];
            }
            else
            {
                otherEdgeIndex[1] = fEdges[edgeI];
            }
        }
    }

    // Modify point-labels on the quad face under consideration
    replaceLabel
    (
        otherEdgePoint[0],
        newPtIndex[0],
        thisFace
    );

    replaceLabel
    (
        nextToOtherPoint[1],
        newPtIndex[1],
        thisFace
    );

#   ifdef FULLDEBUG
    if (debug)
    {
        Info << "Modified thisFace: " << fIndex
             << ": " << thisFace << endl;
    }
#   endif

    // Find the interior face that contains otherEdgeIndex[1]
    found = false;

    labelList& e1 = faceEdges_[c0IntIndex[0]];

    forAll(e1, edgeI)
    {
        if (e1[edgeI] == otherEdgeIndex[1])
        {
            replaceLabel(c0IntIndex[0], -1, cell_0);
            replaceFace = c0IntIndex[0];
            found = true; break;
        }
    }

    if (!found)
    {
        // The edge was not found before
        replaceLabel(c0IntIndex[1], -1, cell_0);
        replaceFace = c0IntIndex[1];
    }

    // Check if face reversal is necessary for the replacement
    if (owner_[replaceFace] == c0)
    {
        if (neighbour_[replaceFace] == -1)
        {
            // Change the owner
            owner_[replaceFace] = newCellIndex[0];
        }
        else
        {
            // This face has to be reversed
            faces_[replaceFace] = faces_[replaceFace].reverseFace();
            owner_[replaceFace] = neighbour_[replaceFace];
            neighbour_[replaceFace] = newCellIndex[0];
        }
    }
    else
    {
        // Keep owner, but change neighbour
        neighbour_[replaceFace] = newCellIndex[0];
    }

    // Define the faces for the new cell
    newCell0[0] = c0BdyIndex[1];
    newCell0[1] = replaceFace;

    // Define the set of new faces and insert...

    // New interior face; Owner = cell[0] & Neighbour = newCell[0]
    tmpQuadFace[0] = otherPointIndex[0];
    tmpQuadFace[1] = newPtIndex[0];
    tmpQuadFace[2] = newPtIndex[1];
    tmpQuadFace[3] = otherPointIndex[1];

    newFaceIndex[0] = insertFace
                      (
                          -1,
                          tmpQuadFace,
                          c0,
                          newCellIndex[0]
                      );

    replaceLabel(-1, newFaceIndex[0], newCell0);
    replaceLabel(-1, newFaceIndex[0], cell_0);

    // remove2DSliver requires this face index for removal
    bisectInteriorFace_ = newFaceIndex[0];

    // Second boundary face; Owner = newCell[0] & Neighbour = [-1]
    tmpTriFace[0] = otherPointIndex[0];
    tmpTriFace[1] = newPtIndex[0];
    tmpTriFace[2] = otherEdgePoint[0];

    newFaceIndex[1] = insertFace
                      (
                          whichPatch(c0BdyIndex[0]),
                          tmpTriFace,
                          newCellIndex[0],
                          -1
                      );

    replaceLabel(-1, newFaceIndex[1], newCell0);

    // Third boundary face; Owner = c[0] & Neighbour = [-1]
    tmpTriFace[0] = otherPointIndex[1];
    tmpTriFace[1] = newPtIndex[1];
    tmpTriFace[2] = otherEdgePoint[1];

    newFaceIndex[2] = insertFace
                      (
                          whichPatch(c0BdyIndex[1]),
                          tmpTriFace,
                          c0,
                          -1
                      );

    replaceLabel(-1, newFaceIndex[2], cell_0);

    // Create / modify edges...
    labelList tmpTriEdgeFaces(3, -1);

    // The edge bisecting the zeroth boundary triangular face
    tmpTriEdgeFaces[0] = c0BdyIndex[0];
    tmpTriEdgeFaces[1] = newFaceIndex[0];
    tmpTriEdgeFaces[2] = newFaceIndex[1];

    newEdgeIndex[1] = insertEdge
                      (
                          whichEdgePatch(commonEdgeIndex[0]),
                          edge(newPtIndex[0], otherPointIndex[0]),
                          tmpTriEdgeFaces
                      );

    // The edge bisecting the first boundary triangular face
    tmpTriEdgeFaces[0] = c0BdyIndex[1];
    tmpTriEdgeFaces[1] = newFaceIndex[0];
    tmpTriEdgeFaces[2] = newFaceIndex[2];

    newEdgeIndex[2] = insertEdge
                      (
                          whichEdgePatch(commonEdgeIndex[1]),
                          edge(newPtIndex[1], otherPointIndex[1]),
                          tmpTriEdgeFaces
                      );

    if (c1 == -1)
    {
        // The quad boundary face resulting from bisection;
        // Owner = newCell[0] & Neighbour = [-1]
        tmpQuadFace[0] = newPtIndex[1];
        tmpQuadFace[1] = nextToOtherPoint[1];
        tmpQuadFace[2] = otherEdgePoint[0];
        tmpQuadFace[3] = newPtIndex[0];

        newFaceIndex[3] = insertFace
                          (
                              whichPatch(fIndex),
                              tmpQuadFace,
                              newCellIndex[0],
                              -1
                          );

        // Add a faceEdges entry as well
        faceEdges_.append(labelList(4, -1));

        replaceLabel(-1, newFaceIndex[3], newCell0);

        labelList tmpBiEdgeFaces(2, -1);

        // The edge bisecting the face
        tmpTriEdgeFaces[0] = fIndex;
        tmpTriEdgeFaces[1] = newFaceIndex[0];
        tmpTriEdgeFaces[2] = newFaceIndex[3];

        newEdgeIndex[0] = insertEdge
                          (
                              whichEdgePatch(otherEdgeIndex[0]),
                              edge(newPtIndex[0], newPtIndex[1]),
                              tmpTriEdgeFaces
                          );

        // Create / replace side edges created from face bisection
        tmpBiEdgeFaces[0] = newFaceIndex[1];
        tmpBiEdgeFaces[1] = newFaceIndex[3];

        newEdgeIndex[3] = insertEdge
                          (
                              whichEdgePatch(commonEdgeIndex[0]),
                              edge(newPtIndex[0], otherEdgePoint[0]),
                              tmpBiEdgeFaces
                          );

        tmpBiEdgeFaces[0] = c0BdyIndex[1];
        tmpBiEdgeFaces[1] = newFaceIndex[3];

        newEdgeIndex[4] = insertEdge
                          (
                              whichEdgePatch(commonEdgeIndex[1]),
                              edge(newPtIndex[1], nextToOtherPoint[1]),
                              tmpBiEdgeFaces
                          );

        // Replace an edge on the bisected face
        replaceLabel
        (
            otherEdgeIndex[1],
            newEdgeIndex[0],
            faceEdges_[fIndex]
        );

        // Now that edges are defined, configure faceEdges
        replaceLabel
        (
            commonEdgeIndex[1],
            newEdgeIndex[4],
            faceEdges_[c0BdyIndex[1]]
        );

#       ifdef FULLDEBUG
        if (debug)
        {
            Info << "Modified Cell[0]: "
                 << c0 << ": " << cell_0 << endl;

            forAll(cell_0, faceI)
            {
                Info << cell_0[faceI]
                     << ": " << faces_[cell_0[faceI]] << endl;
            }

            Info << "New Cell[0]: " << newCellIndex[0]
                 << ": " << newCell0 << endl;

            forAll(newCell0, faceI)
            {
                Info << newCell0[faceI]
                     << ": " << faces_[newCell0[faceI]] << endl;
            }
        }
#       endif
    }
    else
    {
        cell &cell_1 = cells_[c1];

        // Find the prism faces for cell[1].
        findPrismFaces
        (
            fIndex,
            c1,
            c1BdyFace,
            c1BdyIndex,
            c1IntFace,
            c1IntIndex
        );

        // Add a new prism cell to the end of the list
        newCellIndex[1] = cells_.append(cell(5));
        cell &newCell1 = cells_[newCellIndex[1]];

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
        cellParents_.insert(newCellIndex[1],secondParent);

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
            newCellIndex[1],
            objectMap
            (
                newCellIndex[1],
                c1MasterObjects.toc()
            )
        );

        // Add a new element to the lengthScale field
        // (Currently the same as cell[1])
        lengthScale_.append(lengthScale_[c1]);

#       ifdef FULLDEBUG
        if (debug)
        {
            Info << "Cell[1]: " << c1 << ": " << cell_1 << endl;

            forAll(cell_1, faceI)
            {
                Info << cell_1[faceI] << ": "
                     << faces_[cell_1[faceI]] << endl;
            }
        }
#       endif

        // Find the interior face that contains otherEdgeIndex[1]
        found = false;

        labelList& e2 = faceEdges_[c1IntIndex[0]];

        forAll(e2, edgeI)
        {
            if (e2[edgeI] == otherEdgeIndex[1])
            {
                replaceLabel(c1IntIndex[0], -1, cell_1);
                replaceFace = c1IntIndex[0];
                found = true; break;
            }
        }

        if (!found)
        {
            // The edge was not found before
            replaceLabel(c1IntIndex[1], -1, cell_1);
            replaceFace = c1IntIndex[1];
        }

        // Check if face reversal is necessary for the replacement
        if (owner_[replaceFace] == c1)
        {
            if (neighbour_[replaceFace] == -1)
            {
                // Change the owner
                owner_[replaceFace] = newCellIndex[1];
            }
            else
            {
                // This face has to be reversed
                faces_[replaceFace] = faces_[replaceFace].reverseFace();
                owner_[replaceFace] = neighbour_[replaceFace];
                neighbour_[replaceFace] = newCellIndex[1];
            }
        }
        else
        {
            // Keep owner, but change neighbour
            neighbour_[replaceFace] = newCellIndex[1];
        }

        // Define attributes for the new prism cell
        newCell1[0] = replaceFace;

        // The interior face resulting from bisection;
        // Owner = newCell[0] & Neighbour = newCell[1]
        tmpQuadFace[0] = newPtIndex[1];
        tmpQuadFace[1] = nextToOtherPoint[1];
        tmpQuadFace[2] = otherEdgePoint[0];
        tmpQuadFace[3] = newPtIndex[0];

        newFaceIndex[3] = insertFace
                          (
                              -1,
                              tmpQuadFace,
                              newCellIndex[0],
                              newCellIndex[1]
                          );

        replaceLabel(-1, newFaceIndex[3], newCell0);
        replaceLabel(-1, newFaceIndex[3], newCell1);
        newCell1[1] = newFaceIndex[3];

        // Check for common edges among the two boundary faces
        // Find the isolated point on both boundary faces of cell[1]
        if
        (
            findCommonEdge(c1BdyIndex[0], c0BdyIndex[0], commonEdgeIndex[0])
        )
        {
            findCommonEdge(c1BdyIndex[1], c0BdyIndex[1], commonEdgeIndex[1]);

            commonFaceIndex[0] = c1BdyIndex[0];
            commonFaceIndex[1] = c1BdyIndex[1];
        }
        else
        {
            findCommonEdge(c1BdyIndex[0], c0BdyIndex[1], commonEdgeIndex[1]);
            findCommonEdge(c1BdyIndex[1], c0BdyIndex[0], commonEdgeIndex[0]);

            commonFaceIndex[0] = c1BdyIndex[1];
            commonFaceIndex[1] = c1BdyIndex[0];
        }

        commonEdges[0] = edges_[commonEdgeIndex[0]];
        commonEdges[1] = edges_[commonEdgeIndex[1]];

        findIsolatedPoint
        (
            faces_[commonFaceIndex[0]],
            commonEdges[0],
            otherPointIndex[2],
            nextToOtherPoint[2]
        );

        findIsolatedPoint
        (
            faces_[commonFaceIndex[1]],
            commonEdges[1],
            otherPointIndex[3],
            nextToOtherPoint[3]
        );

        // For convenience...
        otherEdgePoint[0] = commonEdges[0].otherVertex(nextToOtherPoint[2]);
        otherEdgePoint[1] = commonEdges[1].otherVertex(nextToOtherPoint[3]);

        // Modify the two existing triangle boundary faces
        // Zeroth boundary face - Owner = newCell[1], Neighbour = -1
        replaceLabel
        (
            otherEdgePoint[0],
            newPtIndex[0],
            faces_[commonFaceIndex[0]]
        );

        owner_[commonFaceIndex[0]] = newCellIndex[1];
        replaceLabel(commonFaceIndex[0], -1, cell_1);
        newCell1[2] = commonFaceIndex[0];

        // First boundary face - Owner = c[1] & Neighbour [-1] (unchanged)
        replaceLabel
        (
            otherEdgePoint[1],
            newPtIndex[1],
            faces_[commonFaceIndex[1]]
        );

        // New interior face; Owner = cell[1] & Neighbour = newCell[1]
        tmpQuadFace[0] = newPtIndex[0];
        tmpQuadFace[1] = otherPointIndex[2];
        tmpQuadFace[2] = otherPointIndex[3];
        tmpQuadFace[3] = newPtIndex[1];

        newFaceIndex[4] = insertFace
                          (
                              -1,
                              tmpQuadFace,
                              c1,
                              newCellIndex[1]
                          );

        replaceLabel(-1, newFaceIndex[4], newCell1);
        replaceLabel(-1, newFaceIndex[4], cell_1);

        // Second boundary face; Owner = cell[1] & Neighbour [-1]
        tmpTriFace[0] = otherPointIndex[2];
        tmpTriFace[1] = newPtIndex[0];
        tmpTriFace[2] = otherEdgePoint[0];

        newFaceIndex[5] = insertFace
                          (
                              whichPatch(commonFaceIndex[0]),
                              tmpTriFace,
                              c1,
                              -1
                          );

        replaceLabel(-1, newFaceIndex[5], cell_1);

        // Third boundary face; Owner = newCell[1] & Neighbour [-1]
        tmpTriFace[0] = otherPointIndex[3];
        tmpTriFace[1] = newPtIndex[1];
        tmpTriFace[2] = otherEdgePoint[1];

        newFaceIndex[6] = insertFace
                          (
                              whichPatch(commonFaceIndex[1]),
                              tmpTriFace,
                              newCellIndex[1],
                              -1
                          );

        replaceLabel(-1, newFaceIndex[6], newCell1);

        // Create/modify edges...
        labelList tmpQuadEdgeFaces(4, -1);

        // The internal edge bisecting the face
        tmpQuadEdgeFaces[0] = fIndex;
        tmpQuadEdgeFaces[1] = newFaceIndex[0];
        tmpQuadEdgeFaces[2] = newFaceIndex[3];
        tmpQuadEdgeFaces[3] = newFaceIndex[4];

        newEdgeIndex[0] = insertEdge
                          (
                              -1,
                              edge(newPtIndex[0], newPtIndex[1]),
                              tmpQuadEdgeFaces
                          );

        // Create / replace side edges created from face bisection
        tmpTriEdgeFaces[0] = commonFaceIndex[0];
        tmpTriEdgeFaces[1] = newFaceIndex[3];
        tmpTriEdgeFaces[2] = newFaceIndex[1];

        newEdgeIndex[3] = insertEdge
                          (
                              whichEdgePatch(commonEdgeIndex[0]),
                              edge(newPtIndex[0], nextToOtherPoint[0]),
                              tmpTriEdgeFaces
                          );

        tmpTriEdgeFaces[0] = c0BdyIndex[1];
        tmpTriEdgeFaces[1] = newFaceIndex[3];
        tmpTriEdgeFaces[2] = newFaceIndex[6];

        newEdgeIndex[4] = insertEdge
                          (
                              whichEdgePatch(commonEdgeIndex[1]),
                              edge(newPtIndex[1], otherEdgePoint[1]),
                              tmpTriEdgeFaces
                          );

        // The edge bisecting the second boundary triangular face
        tmpTriEdgeFaces[0] = commonFaceIndex[0];
        tmpTriEdgeFaces[1] = newFaceIndex[4];
        tmpTriEdgeFaces[2] = newFaceIndex[5];

        newEdgeIndex[5] = insertEdge
                          (
                              whichEdgePatch(commonEdgeIndex[0]),
                              edge(newPtIndex[0], otherPointIndex[2]),
                              tmpTriEdgeFaces
                          );

        // The edge bisecting the third boundary triangular face
        tmpTriEdgeFaces[0] = commonFaceIndex[1];
        tmpTriEdgeFaces[1] = newFaceIndex[4];
        tmpTriEdgeFaces[2] = newFaceIndex[6];

        newEdgeIndex[6] = insertEdge
                          (
                              whichEdgePatch(commonEdgeIndex[1]),
                              edge(newPtIndex[1], otherPointIndex[3]),
                              tmpTriEdgeFaces
                          );

#       ifdef FULLDEBUG
        if (debug)
        {
            Info << nl << "Modified Cell[0]: "
                 << c0 << ": " << cell_0 << endl;

            forAll(cell_0, faceI)
            {
                Info << cell_0[faceI]
                     << ": " << faces_[cell_0[faceI]] << endl;
            }

            Info << "New Cell[0]: "
                 << newCellIndex[0] << ": " << newCell0 << endl;

            forAll(newCell0, faceI)
            {
                Info << newCell0[faceI] << ": "
                     << faces_[newCell0[faceI]] << endl;
            }

            Info << nl << "Modified Cell[1]: "
                 << c1 << ": " << cell_1 << endl;

            forAll(cell_1, faceI)
            {
                Info << cell_1[faceI] << ": "
                     << faces_[cell_1[faceI]] << endl;
            }

            Info << "New Cell[1]: "
                 << newCellIndex[1] << ": " << newCell1 << endl;

            forAll(newCell1, faceI)
            {
                Info << newCell1[faceI] << ": "
                     << faces_[newCell1[faceI]] << endl;
            }
        }
#       endif
    }

    // Set the flag
    topoChangeFlag_ = true;

    // Update the number of cells
    nCells_++;

    if (c1 != -1)
    {
        nCells_++;
    }
}

// Method for the collapse of a quad-face in 2D
void dynamicTopoFvMesh::collapseQuadFace
(
    const label fIndex
)
{
    // Obtain a reference for this face...
    face& thisFace = faces_[fIndex];

    // This face is to be collapsed...
    if (debug)
    {
        Info << nl << nl
             << "Face: " << fIndex << ": " << thisFace
             << " is to be collapsed. " << endl;
    }

    // Local variables
    FixedList<label,2> c0BdyIndex, c0IntIndex, c1BdyIndex, c1IntIndex;
    FixedList<face,2> c0BdyFace, c0IntFace, c1BdyFace, c1IntFace;
    FixedList<edge,4> checkEdge(edge(-1,-1));
    FixedList<label,4> checkEdgeIndex;
    face tmpTriFace(3);

    // Define checkEdges
    checkEdgeIndex[0] = getTriBoundaryEdge(fIndex);
    checkEdge[0] = edges_[checkEdgeIndex[0]];

    labelList& fEdges = faceEdges_[fIndex];

    forAll(fEdges, edgeI)
    {
        if (checkEdgeIndex[0] != fEdges[edgeI])
        {
            edge& thisEdge = edges_[fEdges[edgeI]];

            if
            (
                (
                    checkEdge[0].start() == thisEdge[0]
                 || checkEdge[0].start() == thisEdge[1]
                )
            )
            {
                checkEdgeIndex[1] = fEdges[edgeI];
                checkEdge[1] = thisEdge;
            }
            else
            if
            (
                (
                    checkEdge[0].end() == thisEdge[0]
                 || checkEdge[0].end() == thisEdge[1]
                )
            )
            {
                checkEdgeIndex[2] = fEdges[edgeI];
                checkEdge[2] = thisEdge;
            }
            else
            {
                checkEdgeIndex[3] = fEdges[edgeI];
                checkEdge[3] = thisEdge;
            }
        }
    }

    // Determine if either edge belongs to a boundary
    bool firstEdgeBoundary  = (whichEdgePatch(checkEdgeIndex[1]) > -1);
    bool secondEdgeBoundary = (whichEdgePatch(checkEdgeIndex[2]) > -1);

    // Build a hull of cells and tri-faces that are connected to each edge
    labelHashSet firstHullCells, secondHullCells;
    labelHashSet firstHullTriFaces, secondHullTriFaces;

    constructPrismHull
    (
        checkEdgeIndex[1],
        firstHullTriFaces,
        firstHullCells
    );

    constructPrismHull
    (
        checkEdgeIndex[2],
        secondHullTriFaces,
        secondHullCells
    );

    // Obtain lists from hashSets
    labelList firstCells = firstHullCells.toc();
    labelList secondCells = secondHullCells.toc();
    labelList firstTriFaces = firstHullTriFaces.toc();
    labelList secondTriFaces = secondHullTriFaces.toc();

    // Obtain references to edgeFaces
    labelList& firstEdgeFaces = edgeFaces_[checkEdgeIndex[1]];
    labelList& secondEdgeFaces = edgeFaces_[checkEdgeIndex[2]];

#   ifdef FULLDEBUG
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

        Info << nl << "First Edge Face Hull: " << firstEdgeFaces << endl;

        forAll(firstEdgeFaces,indexI)
        {
            Info << firstEdgeFaces[indexI]
                 << ": " << faces_[firstEdgeFaces[indexI]] << endl;
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

        Info << nl << "Second Edge Face Hull: " << secondEdgeFaces << endl;

        forAll(secondEdgeFaces, indexI)
        {
            Info << secondEdgeFaces[indexI]
                 << ": " << faces_[secondEdgeFaces[indexI]] << endl;
        }
    }
#   endif

    // Determine the common vertices for the first and second edges
    label cv0 = checkEdge[1].commonVertex(checkEdge[0]);
    label cv1 = checkEdge[1].commonVertex(checkEdge[3]);
    label cv2 = checkEdge[2].commonVertex(checkEdge[0]);
    label cv3 = checkEdge[2].commonVertex(checkEdge[3]);

    // Determine the neighbouring cells
    label c0 = owner_[fIndex], c1 = neighbour_[fIndex];

    // Find the prism-faces
    FixedList<label,2> faceToKeep(0), faceToThrow(0);

    findPrismFaces
    (
        fIndex,
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
            fIndex,
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
            << "Face: " << fIndex << ": " << thisFace << endl;

        if
        (
            boundaryMesh()[whichEdgePatch(checkEdgeIndex[1])].type()
            == "symmetryPlane"
        )
        {
            secondEdgeBoundary = false;
        }

        if (secondEdgeBoundary)
        {
            if
            (
                boundaryMesh()[whichEdgePatch(checkEdgeIndex[2])].type()
                == "symmetryPlane"
            )
            {
                firstEdgeBoundary = false;
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
                return;
            }
        }

        // Collapse to the second node...
        forAll(firstEdgeFaces,faceI)
        {
            face& replacementFace = faces_[firstEdgeFaces[faceI]];
            replaceLabel(cv0,cv2,replacementFace);
            replaceLabel(cv1,cv3,replacementFace);

            // Determine the quad-face in cell[0] & cell[1]
            // that belongs to firstEdgeFaces
            if (firstEdgeFaces[faceI] == c0IntIndex[0])
            {
                faceToKeep[0]  = c0IntIndex[1];
                faceToThrow[0] = c0IntIndex[0];
            }

            if (firstEdgeFaces[faceI] == c0IntIndex[1])
            {
                faceToKeep[0]  = c0IntIndex[0];
                faceToThrow[0] = c0IntIndex[1];
            }

            if (c1 != -1)
            {
                if (firstEdgeFaces[faceI] == c1IntIndex[0])
                {
                    faceToKeep[1]  = c1IntIndex[1];
                    faceToThrow[1] = c1IntIndex[0];
                }

                if (firstEdgeFaces[faceI] == c1IntIndex[1])
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
        if (cv0 < nOldPoints_)
        {
            reversePointMap_[cv0] = -1;
        }

        if (cv1 < nOldPoints_)
        {
            reversePointMap_[cv1] = -1;
        }
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
                return;
            }
        }

        // Collapse to the first node by default...
        forAll(secondEdgeFaces,faceI)
        {
            face& replacementFace = faces_[secondEdgeFaces[faceI]];
            replaceLabel(cv2, cv0, replacementFace);
            replaceLabel(cv3, cv1, replacementFace);

            // Determine the quad-face(s) in cell[0] & cell[1]
            // that belongs to secondEdgeFaces
            if (secondEdgeFaces[faceI] == c0IntIndex[0])
            {
                faceToKeep[0]  = c0IntIndex[1];
                faceToThrow[0] = c0IntIndex[0];
            }

            if (secondEdgeFaces[faceI] == c0IntIndex[1])
            {
                faceToKeep[0]  = c0IntIndex[0];
                faceToThrow[0] = c0IntIndex[1];
            }

            if (c1 != -1)
            {
                if (secondEdgeFaces[faceI] == c1IntIndex[0])
                {
                    faceToKeep[1]  = c1IntIndex[1];
                    faceToThrow[1] = c1IntIndex[0];
                }

                if (secondEdgeFaces[faceI] == c1IntIndex[1])
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

#   ifdef FULLDEBUG
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

        Info << nl << "First Edge Face Hull: " << firstEdgeFaces << endl;

        forAll(firstEdgeFaces, indexI)
        {
            Info << firstEdgeFaces[indexI]
                 << ": " << faces_[firstEdgeFaces[indexI]] << endl;
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

        Info << nl << "Second Edge Face Hull: " << secondEdgeFaces << endl;

        forAll(secondEdgeFaces, indexI)
        {
            Info << secondEdgeFaces[indexI]
                 << ": " << faces_[secondEdgeFaces[indexI]] << endl;
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
#   endif

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
        label newFaceIndex = insertFace
                             (
                                 whichPatch(faceToThrow[0]),
                                 faces_[faceToKeep[0]],
                                 owner_[faceToKeep[0]],
                                 -1
                             );

        replaceLabel
        (
            faceToKeep[0],
            newFaceIndex,
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
        if (cell_0[faceI] != fIndex && cell_0[faceI] != faceToKeep[0])
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
            label newFaceIndex = insertFace
                                 (
                                     whichPatch(faceToThrow[1]),
                                     faces_[faceToKeep[1]],
                                     owner_[faceToKeep[1]],
                                     -1
                                 );

            replaceLabel
            (
                faceToKeep[1],
                newFaceIndex,
                cells_[owner_[faceToKeep[1]]]
            );

            // Renumber the neighbour so that this face is removed correctly.
            neighbour_[faceToKeep[1]] = 0;
            removeFace(faceToKeep[1]);
        }

        cell &cell_1 = cells_[c1];

        forAll(cell_1, faceI)
        {
            if (cell_1[faceI] != fIndex && cell_1[faceI] != faceToKeep[1])
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
    removeFace(fIndex);

    // Set the flag
    topoChangeFlag_ = true;
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
    //      Update faceEdges, edgeFaces and edgePoints information

    // Figure out which thread this is...
    label tIndex = self();

    // Try to write-lock this edge.
    if (tryEdgeLock(eIndex, rwMutex::WRITE_LOCK))
    {
        edgeStack(tIndex).push(eIndex);
        return;
    }

    if
    (
        (nModifications_ > maxModifications_)
     && (maxModifications_ > -1)
    )
    {
        // Reached the max allowable topo-changes.
        edgeStack(tIndex).clear();
        return;
    }

    // Hull variables
    face tmpTriFace(3);
    labelList tmpEdgeFaces(3,-1);
    labelList tmpIntEdgeFaces(4,-1);
    labelList tmpEdgePoints(3,-1);
    labelList tmpIntEdgePoints(4,-1);
    labelList tmpFaceEdges(3,-1);
    edge& thisEdge = edges_[eIndex];
    labelList& vertexHull = edgePoints_[eIndex];
    label m = vertexHull.size();

#   ifdef FULLDEBUG
    if (debug)
    {
        Info << nl << nl << "Edge: " << eIndex
             << ": " << thisEdge << " is to be bisected. " << endl;
    }
#   endif

    // Size up the hull lists
    labelList cellHull(m, -1);
    labelList faceHull(m, -1);
    labelList edgeHull(m, -1);
    labelListList ringEntities(4, labelList(m, -1));

    // Construct a hull around this edge, and write-lock entities
    if
    (
        constructHull
        (
            eIndex,
            edgeHull,
            faceHull,
            cellHull,
            ringEntities,
            rwMutex::WRITE_LOCK
        )
    )
    {
        // Put this edge back on the stack and bail out
        edgeStack(tIndex).push(eIndex);
        return;
    }

    // Write lock the point mutex
    pMutex_.lock(rwMutex::WRITE_LOCK);

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

    // Unlock the point mutex from write lock
    pMutex_.unlock();

    // Write lock the edge mutex
    eMutex_.lock(rwMutex::WRITE_LOCK);

    // Add a new edge to the end of the list
    label newEdgeIndex =
        insertEdge
        (
            whichEdgePatch(eIndex),
            edge(newPointIndex,thisEdge[1]),
            labelList(faceHull.size(),-1),
            vertexHull
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

    // Unlock the edge mutex from write lock
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

        // Modify edgePoints for the edge connected to thisEdge[0]
        replaceLabel
        (
            newEdge[1],
            newPointIndex,
            edgePoints_[ringEntities[0][indexI]]
        );

        // Obtain circular indices
        label nextI = vertexHull.fcIndex(indexI);
        label prevI = vertexHull.rcIndex(indexI);

        // Check if this is an interior/boundary face
        if (cellHull[indexI] != -1)
        {
            cell& currCell = cells_[cellHull[indexI]];

            // Write lock the cell mutex
            cMutex_.lock(rwMutex::WRITE_LOCK);

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

            // Unlock the cell mutex from write lock
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

            // Write lock the cell mutex
            cMutex_.lock(rwMutex::WRITE_LOCK);

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

            // Unlock the cell mutex from write lock
            cMutex_.unlock();

            // Configure the interior face
            tmpTriFace[0] = vertexHull[nextI];
            tmpTriFace[1] = vertexHull[indexI];
            tmpTriFace[2] = newPointIndex;

            // Write lock the face mutex
            fMutex_.lock(rwMutex::WRITE_LOCK);

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

            // Unlock the face mutex from write lock
            fMutex_.unlock();

            // Add to the new cell
            newCell[0] = addedIntFaceIndices[indexI];

            // Modify the existing ring face connected to newEdge[1]
            label replaceFace = ringEntities[3][indexI];

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

            // Modify the edge on the ring.
            // Add the new interior face to edgeFaces.
            sizeUpList
            (
                addedIntFaceIndices[indexI],
                edgeFaces_[edgeHull[indexI]]
            );

            // Insert the new point to edgePoints for the ring edge
            insertLabel
            (
                newPointIndex,
                thisEdge[0],
                newEdge[1],
                edgePoints_[edgeHull[indexI]]
            );

            // Add this edge to faceEdges for the new interior face
            faceEdges_[addedIntFaceIndices[indexI]][0] = edgeHull[indexI];

            // Replace face labels
            replaceLabel
            (
                replaceFace,
                addedIntFaceIndices[indexI],
                currCell
            );

            // Add to the new cell
            newCell[1] = replaceFace;

            // Check if this is a boundary face
            if (cellHull[prevI] == -1)
            {
                // Configure the boundary face
                tmpTriFace[0] = newPointIndex;
                tmpTriFace[1] = newEdge[1];
                tmpTriFace[2] = vertexHull[indexI];

                // Write lock the face mutex
                fMutex_.lock(rwMutex::WRITE_LOCK);

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

                // Configure edgePoints
                tmpEdgePoints[0] = thisEdge[0];
                tmpEdgePoints[1] = vertexHull[nextI];
                tmpEdgePoints[2] = newEdge[1];

                // Write lock the edge mutex
                eMutex_.lock(rwMutex::WRITE_LOCK);

                // Add an edge
                addedEdgeIndices[indexI] =
                    insertEdge
                    (
                        whichPatch(faceHull[indexI]),
                        edge(newPointIndex,vertexHull[indexI]),
                        tmpEdgeFaces,
                        tmpEdgePoints
                    );

                // Unlock the edge mutex from write lock
                eMutex_.unlock();

                // Add this edge to the interior-face faceEdges entry
                faceEdges_[addedIntFaceIndices[indexI]][1] =
                    addedEdgeIndices[indexI];

                // Configure faceEdges for this boundary face
                tmpFaceEdges[0] = addedEdgeIndices[indexI];
                tmpFaceEdges[1] = newEdgeIndex;
                tmpFaceEdges[2] = ringEntities[2][indexI];

                // Modify faceEdges for the hull face
                replaceLabel
                (
                    ringEntities[2][indexI],
                    addedEdgeIndices[indexI],
                    faceEdges_[faceHull[indexI]]
                );

                // Modify edgeFaces for the edge connected to newEdge[1]
                replaceLabel
                (
                    faceHull[indexI],
                    addedFaceIndices[indexI],
                    edgeFaces_[ringEntities[2][indexI]]
                );

                // Modify edgePoints for the edge connected to newEdge[1]
                replaceLabel
                (
                    thisEdge[0],
                    newPointIndex,
                    edgePoints_[ringEntities[2][indexI]]
                );

                // Add the faceEdges entry
                faceEdges_.append(tmpFaceEdges);

                // Unlock the face mutex from write lock
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

                // Write lock the face mutex
                fMutex_.lock(rwMutex::WRITE_LOCK);

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

                // Configure edgePoints
                tmpIntEdgePoints[0] = thisEdge[0];
                tmpIntEdgePoints[1] = vertexHull[nextI];
                tmpIntEdgePoints[2] = newEdge[1];
                tmpIntEdgePoints[3] = vertexHull[prevI];

                // Write lock the edge mutex
                eMutex_.lock(rwMutex::WRITE_LOCK);

                // Add an internal edge
                addedEdgeIndices[indexI] =
                    insertEdge
                    (
                        -1,
                        edge(newPointIndex,vertexHull[indexI]),
                        tmpIntEdgeFaces,
                        tmpIntEdgePoints
                    );

                // Unlock the edge mutex from write lock
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
                tmpFaceEdges[2] = ringEntities[2][indexI];

                // Modify faceEdges for the hull face
                replaceLabel
                (
                    ringEntities[2][indexI],
                    addedEdgeIndices[indexI],
                    faceEdges_[faceHull[indexI]]
                );

                // Modify edgeFaces for the edge connected to newEdge[1]
                replaceLabel
                (
                    faceHull[indexI],
                    addedFaceIndices[indexI],
                    edgeFaces_[ringEntities[2][indexI]]
                );

                // Modify edgePoints for the edge connected to newEdge[1]
                replaceLabel
                (
                    thisEdge[0],
                    newPointIndex,
                    edgePoints_[ringEntities[2][indexI]]
                );

                // Add the faceEdges entry
                faceEdges_.append(tmpFaceEdges);

                // Unlock the face mutex from write lock
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

                // Write lock the face mutex
                fMutex_.lock(rwMutex::WRITE_LOCK);

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

                // Configure edgePoints
                tmpIntEdgePoints[0] = thisEdge[0];
                tmpIntEdgePoints[1] = vertexHull[1];
                tmpIntEdgePoints[2] = newEdge[1];
                tmpIntEdgePoints[3] = vertexHull[indexI];

                // Write lock the edge mutex
                eMutex_.lock(rwMutex::WRITE_LOCK);

                // Add an internal edge
                addedEdgeIndices[0] =
                    insertEdge
                    (
                        -1,
                        edge(newPointIndex,vertexHull[0]),
                        tmpIntEdgeFaces,
                        tmpIntEdgePoints
                    );

                // Unlock the edge mutex from write lock
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
                tmpFaceEdges[2] = ringEntities[2][0];

                // Modify faceEdges for the hull face
                replaceLabel
                (
                    ringEntities[2][0],
                    addedEdgeIndices[0],
                    faceEdges_[faceHull[0]]
                );

                // Modify edgeFaces for the edge connected to newEdge[1]
                replaceLabel
                (
                    faceHull[0],
                    addedFaceIndices[0],
                    edgeFaces_[ringEntities[2][0]]
                );

                // Modify edgePoints for the edge connected to newEdge[1]
                replaceLabel
                (
                    thisEdge[0],
                    newPointIndex,
                    edgePoints_[ringEntities[2][0]]
                );

                // Add the faceEdges entry
                faceEdges_.append(tmpFaceEdges);

                // Unlock the face mutex from write lock
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

            // Write lock the face mutex
            fMutex_.lock(rwMutex::WRITE_LOCK);

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

            // Configure edgePoints
            tmpEdgePoints[0] = newEdge[1];
            tmpEdgePoints[1] = vertexHull[prevI];
            tmpEdgePoints[2] = thisEdge[0];

            // Write lock the edge mutex
            eMutex_.lock(rwMutex::WRITE_LOCK);

            // Add an edge
            addedEdgeIndices[indexI] =
                insertEdge
                (
                    whichPatch(faceHull[indexI]),
                    edge(newPointIndex,vertexHull[indexI]),
                    tmpEdgeFaces,
                    tmpEdgePoints
                );

            // Unlock the edge mutex from write lock
            eMutex_.unlock();

            // Add a faceEdges entry to the previous interior face
            faceEdges_[addedIntFaceIndices[prevI]][2] =
                addedEdgeIndices[indexI];

            // Configure faceEdges for the final boundary face
            tmpFaceEdges[0] = addedEdgeIndices[indexI];
            tmpFaceEdges[1] = newEdgeIndex;
            tmpFaceEdges[2] = ringEntities[2][indexI];

            // Modify faceEdges for the hull face
            replaceLabel
            (
                ringEntities[2][indexI],
                addedEdgeIndices[indexI],
                faceEdges_[faceHull[indexI]]
            );

            // Modify edgeFaces for the edge connected to newEdge[1]
            replaceLabel
            (
                faceHull[indexI],
                addedFaceIndices[indexI],
                edgeFaces_[ringEntities[2][indexI]]
            );

            // Modify edgePoints for the edge connected to newEdge[1]
            replaceLabel
            (
                thisEdge[0],
                newPointIndex,
                edgePoints_[ringEntities[2][indexI]]
            );

            // Add the faceEdges entry
            faceEdges_.append(tmpFaceEdges);

            // Unlock the face mutex from write lock
            fMutex_.unlock();

            // Add an entry to newEdgeFaces
            newEdgeFaces[indexI] = addedFaceIndices[indexI];

            // Make the final entry for the previous cell
            cells_[addedCellIndices[prevI]][3] = addedFaceIndices[indexI];
        }
    }

    // Unlock all entities
    unlockMutexLists(tIndex);

#   ifdef FULLDEBUG
    if (debug)
    {
        Info << "Vertices: " << vertexHull << endl;
        Info << "Edges: " << edgeHull << endl;
        Info << "Faces: " << faceHull << endl;
        Info << "Cells: " << cellHull << endl;

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
             << ": " << edges_[newEdgeIndex] << nl
             << " edgeFaces:: " << edgeFaces_[newEdgeIndex] << nl
             << " edgePoints:: " << edgePoints_[newEdgeIndex]
             << endl;

        Info << "Added edges: " << endl;
        forAll(addedEdgeIndices, edgeI)
        {
            Info << addedEdgeIndices[edgeI]
                 << ":: " << edges_[addedEdgeIndices[edgeI]] << nl
                 << " edgeFaces:: " << edgeFaces_[addedEdgeIndices[edgeI]] << nl
                 << " edgePoints:: " << edgePoints_[addedEdgeIndices[edgeI]]
                 << endl;
        }

        Info << "New Point:: " << newPointIndex << endl;
        Info << "pointEdges:: " << pointEdges_[newPointIndex] << endl;
    }
#   endif

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the number of modifications
    nModifications_++;
}

// Method for the collapse of an edge in 3D
void dynamicTopoFvMesh::collapseEdge
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
    //      Update faceEdges, edgeFaces and edgePoints information

    // Figure out which thread this is...
    label tIndex = self();

    // Try to write-lock this edge.
    if (tryEdgeLock(eIndex, rwMutex::WRITE_LOCK))
    {
        edgeStack(tIndex).push(eIndex);
        return;
    }

    if
    (
        (nModifications_ > maxModifications_)
     && (maxModifications_ > -1)
    )
    {
        // Reached the max allowable topo-changes.
        edgeStack(tIndex).clear();
        return;
    }

    // Hull variables
    bool found = false;
    edge& thisEdge = edges_[eIndex];
    labelList& vertexHull = edgePoints_[eIndex];
    label replaceIndex = -1, m = vertexHull.size();
    FixedList<bool,2> edgeBoundary(false);

#   ifdef FULLDEBUG
    if (debug)
    {
        Info << nl << nl << "Edge: " << eIndex
             << ": " << thisEdge << " is to be collapsed. " << endl;
    }
#   endif

    // Size up the hull lists
    labelList cellHull(m, -1);
    labelList faceHull(m, -1);
    labelList edgeHull(m, -1);
    labelListList ringEntities(4, labelList(m, -1));

    // Construct a hull around this edge, and write-lock entities
    if
    (
        constructHull
        (
            eIndex,
            edgeHull,
            faceHull,
            cellHull,
            ringEntities,
            rwMutex::WRITE_LOCK
        )
    )
    {
        // Put this edge back on the stack and bail out
        edgeStack(tIndex).push(eIndex);
        return;
    }

    // Check whether points of the edge lies on a boundary
    checkEdgeBoundary(eIndex, edgeBoundary);

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
        Info << "ringEntities (removed faces): " << endl;

        forAll(ringEntities[removeFaceIndex], faceI)
        {
            label fIndex = ringEntities[removeFaceIndex][faceI];

            if (fIndex != -1)
            {
                Info << fIndex << ": " << faces_[fIndex] << endl;
            }
        }

        Info << "ringEntities (removed edges): " << endl;
        forAll(ringEntities[removeEdgeIndex], edgeI)
        {
            label ieIndex = ringEntities[removeEdgeIndex][edgeI];

            if (ieIndex != -1)
            {
                Info << ieIndex << ": " << edges_[ieIndex] << endl;
            }
        }

        Info << "ringEntities (replacement faces): " << endl;
        forAll(ringEntities[replaceFaceIndex], faceI)
        {
            label fIndex = ringEntities[replaceFaceIndex][faceI];

            if (fIndex != -1)
            {
                Info << fIndex << ": " << faces_[fIndex] << endl;
            }
        }

        Info << "ringEntities (replacement edges): " << endl;
        forAll(ringEntities[replaceEdgeIndex], edgeI)
        {
            label ieIndex = ringEntities[replaceEdgeIndex][edgeI];

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
                        unlockMutexLists(tIndex);

                        return;
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
                        unlockMutexLists(tIndex);

                        return;
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
        label edgeToRemove = ringEntities[removeEdgeIndex][indexI];
        label faceToRemove = ringEntities[removeFaceIndex][indexI];
        label cellToRemove = cellHull[indexI];
        label replaceEdge = ringEntities[replaceEdgeIndex][indexI];
        label replaceFace = ringEntities[replaceFaceIndex][indexI];

        labelList& rmvEdgeFaces = edgeFaces_[edgeToRemove];
        labelList& rplEdgeFaces = edgeFaces_[replaceEdge];

        // Replace edgePoints for all edges emanating from hullVertices
        // except ring-edges; those are sized-down later
        labelList& hullPointEdges = pointEdges_[vertexHull[indexI]];

        forAll(hullPointEdges, edgeI)
        {
            labelList& hullEdgePoints = edgePoints_[hullPointEdges[edgeI]];

            if
            (
                 foundInList(collapsePoint, hullEdgePoints)
             && !foundInList(replacePoint, hullEdgePoints)
            )
            {
                replaceLabel
                (
                    collapsePoint,
                    replacePoint,
                    hullEdgePoints
                );
            }
        }

        forAll(rmvEdgeFaces, faceI)
        {
            // Replace edge labels for faces
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

            // Hull faces should be removed for the replacement edge
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

            // Need to avoid ring faces as well.
            forAll(ringEntities[removeFaceIndex], faceII)
            {
                if
                (
                    rmvEdgeFaces[faceI]
                 == ringEntities[removeFaceIndex][faceII]
                )
                {
                    found = true;
                    break;
                }
            }

            // Size-up the replacement edge list if the face hasn't been found.
            // These faces are connected to the edge slated for
            // removal, but do not belong to the hull.
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

            // Size down edgePoints for the ring edges
            sizeDownList
            (
                collapsePoint,
                edgePoints_[edgeHull[indexI]]
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

    // Write lock mutexes
    eMutex_.lock(rwMutex::WRITE_LOCK);
    fMutex_.lock(rwMutex::WRITE_LOCK);
    cMutex_.lock(rwMutex::WRITE_LOCK);

    // Remove all hull entities
    forAll(faceHull, indexI)
    {
        label edgeToRemove = ringEntities[removeEdgeIndex][indexI];
        label faceToRemove = ringEntities[removeFaceIndex][indexI];
        label cellToRemove = cellHull[indexI];

        if (cellToRemove != -1)
        {
            // Remove faceToRemove and associated faceEdges
            removeFace(faceToRemove);

            // Remove from list of locked faces
            removeFaceLock(faceToRemove);

            // Remove the hull cell
            cells_.remove(cellToRemove);
            lengthScale_.remove(cellToRemove);

            // Remove from list of locked cells
            removeCellLock(cellToRemove);

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
        removeEdgeLock(edgeToRemove);

        // Remove the hull face
        removeFace(faceHull[indexI]);

        // Remove from list of locked faces
        removeFaceLock(faceHull[indexI]);
    }

    // Unlock mutexes from write lock
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

                    // Look for an edge on this face that doesn't
                    // contain collapsePoint or replacePoint.
                    label rplIndex = -1;
                    labelList& fEdges = faceEdges_[eFaces[faceI]];

                    forAll(fEdges, edgeI)
                    {
                        edge& eCheck = edges_[fEdges[edgeI]];

                        if
                        (
                            eCheck[0] != collapsePoint
                         && eCheck[1] != collapsePoint
                         && eCheck[0] != replacePoint
                         && eCheck[1] != replacePoint
                        )
                        {
                            rplIndex = fEdges[edgeI];
                            break;
                        }
                    }

                    // Modify edgePoints for this edge
                    replaceLabel
                    (
                        collapsePoint,
                        replacePoint,
                        edgePoints_[rplIndex]
                    );
                }
            }
        }
    }

    // At this point, edgePoints for the replacement edges are broken,
    // but edgeFaces are consistent. So use this information to re-build
    // edgePoints for all replacement edges.
    forAll(ringEntities[replaceEdgeIndex], edgeI)
    {
#       ifdef FULLDEBUG
        if (debug)
        {
            Info << "Building edgePoints for edge: "
                 << ringEntities[replaceEdgeIndex][edgeI] << ": "
                 << edges_[ringEntities[replaceEdgeIndex][edgeI]]
                 << endl;
        }
#       endif

        buildEdgePoints(ringEntities[replaceEdgeIndex][edgeI]);
    }

    // Write lock the point mutex
    pMutex_.lock(rwMutex::WRITE_LOCK);

    // Move to the new point
    meshPoints_[replacePoint] = newPoint;

    // Remove the collapse point
    meshPoints_.remove(collapsePoint);
    nPoints_--;

    // Remove from list of locked points
    removePointLock(collapsePoint);

    // Remove the point mutex
    if (threader_->multiThreaded())
    {
        // Unlock it first
        pointMutex_[collapsePoint].unlock();
        pointMutex_.remove(collapsePoint);
    }

    // Null pointEdges so that removeEdge deletes it.
    pointEdges_[collapsePoint] = labelList(0);

    // Unlock the point mutex from write lock
    pMutex_.unlock();

    // Update the reverse point map
    if (collapsePoint < nOldPoints_)
    {
        reversePointMap_[collapsePoint] = -1;
    }

    // Write lock the edge mutex
    eMutex_.lock(rwMutex::WRITE_LOCK);

    // Remove the edge
    removeEdge(eIndex);

    // Remove from list of locked edges
    removeEdgeLock(eIndex);

    // Unlock the edge mutex from write lock
    eMutex_.unlock();

    // Unlock all entities (from write lock)
    unlockMutexLists(tIndex);

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the number of modifications
    nModifications_++;
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
void dynamicTopoFvMesh::remove2DSliver
(
    const label fIndex
)
{
    // Measure the boundary edge-length of the face in question
    scalar length = 0.0;

    edgeLength
    (
        getTriBoundaryEdge(fIndex),
        length
    );

    // Determine the boundary triangular face area
    scalar area = triFaceArea(faces_[getTriBoundaryFace(fIndex)]);

    // This cell has to be removed...
    if (mag(area) < (sliverThreshold_*length*length))
    {
        // Step 1: Bisect the boundary quad face
        bisectInteriorFace_ = -1;
        bisectQuadFace(fIndex);

        // Step 2: Collapse the newly created internal quad face
        collapseQuadFace(bisectInteriorFace_);
    }
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
        Info << "void dynamicTopoFvMesh::mapFields(const mapPolyMesh& meshMap):"
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

// MultiThreaded topology modifier [2D]
void dynamicTopoFvMesh::threadedTopoModifier2D()
{
    if (edgeModification_)
    {
        // Initialize the face stacks
        initFaceStacks();

        if (threader_->multiThreaded())
        {
            // Submit jobs to the work queue
            for (label i = 0; i < threader_->getNumThreads(); i++)
            {
                threader_->addToWorkQueue
                           (
                               &edgeBisectCollapse2D,
                               reinterpret_cast<void *>(&(structPtr_[i]))
                           );
            }

            // Wait for all work to complete
            threader_->waitForCompletion();
        }
        else
        {
            edgeBisectCollapse2D(&(structPtr_[0]));
        }

        if (debug)
        {
            Info << nl << "2D Edge Bisection/Collapse complete." << endl;
        }
    }

    // Re-initialize the face stacks
    initFaceStacks();

    if (threader_->multiThreaded())
    {
        // Submit jobs to the work queue
        for (label i = 0; i < threader_->getNumThreads(); i++)
        {
            threader_->addToWorkQueue
                       (
                           &swap2DEdges,
                           reinterpret_cast<void *>(&(structPtr_[i]))
                       );
        }

        // Wait for all work to complete
        threader_->waitForCompletion();
    }
    else
    {
        swap2DEdges(reinterpret_cast<void *>(&(structPtr_[0])));
    }

    if (debug)
    {
        Info << nl << "2D Edge Swapping complete." << endl;
    }
}

// MultiThreaded topology modifier [3D]
void dynamicTopoFvMesh::threadedTopoModifier3D()
{
    if (edgeModification_)
    {
        // Initialize the edge stacks
        initEdgeStacks();

        if (threader_->multiThreaded())
        {
            // Submit jobs to the work queue
            for (label i = 0; i < threader_->getNumThreads(); i++)
            {
                threader_->addToWorkQueue
                           (
                               &edgeBisectCollapse3D,
                               reinterpret_cast<void *>(&(structPtr_[i]))
                           );
            }

            // Wait for all work to complete
            threader_->waitForCompletion();
        }
        else
        {
            edgeBisectCollapse3D(&(structPtr_[0]));
        }

        if (debug)
        {
            Info << nl << "3D Edge Bisection/Collapse complete." << endl;
        }
    }

    // Re-initialize the edge stacks
    initEdgeStacks();

    if (threader_->multiThreaded())
    {
        // Submit jobs to the work queue
        for (label i = 0; i < threader_->getNumThreads(); i++)
        {
            threader_->addToWorkQueue
                       (
                           &swap3DEdges,
                           reinterpret_cast<void *>(&(structPtr_[i]))
                       );
        }

        // Wait for all work to complete
        threader_->waitForCompletion();
    }
    else
    {
        swap3DEdges(reinterpret_cast<void *>(&(structPtr_[0])));
    }

    if (debug)
    {
        Info << nl << "3D Edge Swapping complete." << endl;
    }
}

// Return reference to the dictionary
const dictionary& dynamicTopoFvMesh::dynamicMeshDict() const
{
    return dict_;
}

// Return reference to the edge mesh
eMesh& dynamicTopoFvMesh::EdgeMesh()
{
    if (!eMeshPtr_.valid())
    {
        FatalErrorIn
        (
            "dynamicTopoFvMesh::edges()"
        )
            << "eMesh has not been allocated."
            << abort(FatalError);
    }

    return eMeshPtr_();
}

// Return the number of edges in the mesh.
// Override of primitiveMesh member function
label dynamicTopoFvMesh::nEdges() const
{
    if (!eMeshPtr_.valid())
    {
        FatalErrorIn
        (
            "dynamicTopoFvMesh::edges()"
        )
            << "eMesh has not been allocated."
            << abort(FatalError);
    }

    return eMeshPtr_->nEdges();
}

// Return the number of internal edges in the mesh.
// Override of primitiveMesh member function
label dynamicTopoFvMesh::nInternalEdges() const
{
    if (!eMeshPtr_.valid())
    {
        FatalErrorIn
        (
            "dynamicTopoFvMesh::edges()"
        )
            << "eMesh has not been allocated."
            << abort(FatalError);
    }

    return eMeshPtr_->nInternalEdges();
}

// Return the ordered list of edges in the mesh.
// Override of primitiveMesh member function.
const edgeList& dynamicTopoFvMesh::edges() const
{
    if (!eMeshPtr_.valid())
    {
        FatalErrorIn
        (
            "dynamicTopoFvMesh::edges()"
        )
            << "eMesh has not been allocated."
            << abort(FatalError);
    }

    return eMeshPtr_->edges();
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

    forAllIter(HashList<point>::iterator, meshPoints_, pIter)
    {
        pIter() = currentPoints[pIter.index()];
    }

    // Track mesh topology modification time
    clockTime topologyTimer;

    //== Connectivity changes ==//

    // Reset the flag
    topoChangeFlag_ = false;
    nModifications_ = 0;

    // Invoke the threaded topoModifier
    if (twoDMesh_)
    {
        threadedTopoModifier2D();
    }
    else
    {
        threadedTopoModifier3D();
    }

    Info << "Topo modifier time: " << topologyTimer.elapsedTime() << endl;

    if (debug)
    {
        // Check connectivity structures for consistency
        checkConnectivity();
    }

    clockTime reOrderingTimer;

    // Apply all pending topology changes, if necessary
    if (topoChangeFlag_)
    {

        // Allocate temporary lists for mesh-reset
        pointField points(nPoints_);
        edgeList edges(nEdges_);
        faceList faces(nFaces_);
        labelList owner(nFaces_);
        labelList neighbour(nInternalFaces_);
        labelListList edgeFaces(nEdges_);
        labelList oldPatchStarts(oldPatchStarts_);
        labelList oldPatchNMeshPoints(oldPatchNMeshPoints_);

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
        reOrderMesh
        (
            points,
            edges,
            faces,
            owner,
            neighbour,
            edgeFaces
        );

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

        // Reset the edge mesh and set pointEdges/faceEdges
        eMeshPtr_->resetPrimitives
        (
            edges,
            edgeFaces,
            edgePatchSizes_,
            edgePatchStarts_
        );

        // Copy edge-based connectivity from eMesh to HashLists
        setEdgeConnectivity();

        // Generate mapping for points on boundary patches
        labelListList patchPointMap(numPatches_);
        for(label i=0; i<numPatches_; i++)
        {
            const labelList& meshPointLabels = boundaryMesh()[i].meshPoints();
            patchNMeshPoints_[i] = meshPointLabels.size();
            patchPointMap[i].setSize(meshPointLabels.size(), -1);
            forAll(meshPointLabels, pointI)
            {
                // Check if the position has been maintained.
                // Otherwise, perform a search for the old position in the patch
                if (pointI < oldPatchNMeshPoints_[i])
                {
                    if
                    (
                        meshPointLabels[pointI]
                     == oldMeshPointLabels[i][pointI]
                    )
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
                oldPatchStarts,
                oldPatchNMeshPoints,
                true
            )
        );

        // Update the underlying mesh, and map all related fields
        updateMesh(mapper_);

        // Discard old information after mapping
        pointPositionsPtr_.clear();
        cellCentresPtr_.clear();
        faceCentresPtr_.clear();

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
        }

        // Clear the current and reverse edge maps
        // Other entities are taken over by mapPolyMesh
        edgeMap_.clear();
        reverseEdgeMap_.clear();
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
