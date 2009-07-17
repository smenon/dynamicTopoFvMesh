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
#include "mapPolyMesh.H"

#include "fvPatchFields.H"
#include "fvsPatchFields.H"
#include "MapFvFields.H"
#include "MeshObject.H"

#include <dlfcn.h>
#include <mpi.h>

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
    interval_(1),
    mapper_(NULL),
    points_(polyMesh::points()),
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
    curvatureRatio_(1.0),
    nModifications_(0),
    nBisections_(-1),
    nCollapses_(-1),
    nSwaps_(-1),
    maxModifications_(-1),
    bisectInterior_(-1),
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

    // Enable/disable run-time debug level
    if (dict_.found("debug"))
    {
        debug = readLabel(dict_.lookup("debug"));
    }

    if (dict_.subDict("dynamicTopoFvMesh").found("interval"))
    {
        interval_ = readLabel
                    (
                        dict_.subDict
                        ("dynamicTopoFvMesh").lookup("interval")
                    );
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

    if (nThreads == 1)
    {
        structPtr_.setSize(1);
        structPtr_.set(0, new topoMeshStruct(this, nThreads));
        structPtr_[0].setMaster();

        // Size the stacks
        if (twoDMesh_)
        {
            faceStack_.setSize(1);
        }
        else
        {
            edgeStack_.setSize(1);
        }
    }
    else
    {
        // Index '0' is master, rest are slaves
        structPtr_.setSize(nThreads + 1);

        // Size the stacks
        if (twoDMesh_)
        {
            faceStack_.setSize(nThreads + 1);
        }
        else
        {
            edgeStack_.setSize(nThreads + 1);
        }

        for (label i = 0; i <= nThreads; i++)
        {
            structPtr_.set(i, new topoMeshStruct(this, nThreads));

            if (i == 0)
            {
                structPtr_[0].ID() = -1;
                structPtr_[0].setMaster();
            }
            else
            {
                structPtr_[i].ID() = threader_->getID(i-1);
                structPtr_[i].setSlave();
            }
        }
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
            (
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
                )
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
        (
            reinterpret_cast<tetMetricReturnType>
            (
                dlsym
                (
                    metricLibPtr,
                    tetMetric.c_str()
                )
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
            (
                dict_.subDict
                (
                    "dynamicTopoFvMesh"
                ).subDict("noSwapPatches").toc()
            );

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
            (
                readLabel
                (
                    dict_.subDict
                    (
                        "dynamicTopoFvMesh"
                    ).lookup("maxTetsPerEdge")
                )
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
            (
                readBool
                (
                    dict_.subDict
                    (
                        "dynamicTopoFvMesh"
                    ).lookup("allowTableResize")
                )
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
        readEdgeOptions();

        // Set curvature patches
        curvatureFields_.setSize(numPatches_, scalarList(0));

        // Initialize the lengthScale field
        lengthScale_.setSize(nCells_, 0.0);
    }
}

// Constructor for topoMeshStruct
dynamicTopoFvMesh::topoMeshStruct::topoMeshStruct
(
    dynamicTopoFvMesh *mesh,
    label nThreads
)
:
    mesh_(mesh),
    nThreads_(nThreads),
    pthreadID_(-1),
    master_(false),
    isCoupled_(false)
{}

// Constructor for patchSubMesh
dynamicTopoFvMesh::patchSubMesh::patchSubMesh()
:
    nPoints_(-1),
    nEdges_(-1),
    nFaces_(-1),
    nCells_(-1)
{}

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

// Set curvature information for a particular patch
void dynamicTopoFvMesh::setCurvatureField
(
    const label pID,
    const scalarField& field
)
{
    curvatureFields_[pID] = 1.0/mag(field);
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

                    // Add to the list of slivers
                    if
                    (
                        (iF[cellI] < sliverThreshold_)
                     && (iF[cellI] > 0.0)
                    )
                    {
                        thresholdSlivers_.insert(cellI);
                    }

                    break;
                }
            }
        }

        // Output statistics:
        if (outputOption || (debug > 0))
        {
            if (minQuality < 0.0)
            {
                WarningIn
                (
                    "dynamicTopoFvMesh::meshQuality()"
                )   << nl << "Minimum cell quality is: " << minQuality << endl;
            }

            if (thresholdSlivers_.size())
            {
                WarningIn
                (
                    "dynamicTopoFvMesh::meshQuality()"
                )   << nl << thresholdSlivers_.size()
                    << " cells are below the specified quality threshold of: "
                    << sliverThreshold_ << endl;
            }

            // Reduce stats across procs.
            label nCells = iF.size();

            reduce(minQuality, minOp<scalar>());
            reduce(maxQuality, maxOp<scalar>());
            reduce(meanQuality, sumOp<scalar>());
            reduce(nCells, sumOp<label>());

            Info << " ~~~ Mesh Quality Statistics ~~~ " << endl;
            Info << " Min: " << minQuality << endl;
            Info << " Max: " << maxQuality << endl;
            Info << " Mean: " << meanQuality/nCells << endl;
            Info << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << endl;
        }
    }

    return tQuality;
}

// Perform a Delaunay test on an internal face
void dynamicTopoFvMesh::testDelaunay
(
    const label fIndex,
    bool& failed
)
{
    failed = false;
    label eIndex = -1, pIndex = -1, fLabel = -1;
    FixedList<bool,2> foundTriFace(false);
    FixedList<FixedList<label,3>,2> triFaces(FixedList<label,3>(-1));

    // Boundary faces are discarded.
    if (whichPatch(fIndex) > -1)
    {
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

                    fLabel = eFaces[faceI];
                }
            }
        }
    }

    // Obtain point references for the first face
    point& a = points_[triFaces[0][0]];

    point cCenter = circumCenter(fLabel);
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
    point& otherPoint = points_[pIndex];

    if
    (
        ((otherPoint - cCenter)&(otherPoint - cCenter)) < rSquared
    )
    {
        // Failed the test.
        failed = true;
    }
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

    if (debug > 2)
    {
        Info << "Inserting face: "
             << newFaceIndex << ": "
             << newFace << endl;
    }

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

// Remove the specified cell from the mesh,
// and add internal faces to the specified patch
void dynamicTopoFvMesh::removeCell
(
    const label index,
    const label patch
)
{
    if (debug > 2)
    {
        Info << "Removed cell: "
             << index << ": "
             << cells_[index] << endl;
    }

    label ptIndex = -1, nextPoint = -1;
    cell& cellToCheck = cells_[index];

    forAll(cellToCheck, faceI)
    {
        labelList& faceEdges = faceEdges_[cellToCheck[faceI]];

        // Delete this face if it's on a boundary
        if (neighbour_[cellToCheck[faceI]] == -1)
        {
            forAll(faceEdges, edgeI)
            {
                // Size down edgeFaces...
                sizeDownList
                (
                    cellToCheck[faceI],
                    edgeFaces_[faceEdges[edgeI]]
                );

                edge& edgeToCheck = edges_[faceEdges[edgeI]];

                // Size-down edgePoints as well
                if (!twoDMesh_)
                {
                    findIsolatedPoint
                    (
                        faces_[cellToCheck[faceI]],
                        edgeToCheck,
                        ptIndex,
                        nextPoint
                    );

                    sizeDownList
                    (
                        ptIndex,
                        edgePoints_[faceEdges[edgeI]]
                    );
                }

                if (edgeFaces_[faceEdges[edgeI]].size() == 0)
                {
                    // Hanging edge. Check its points and remove
                    // them if necessary.
                    forAll(edgeToCheck, pointI)
                    {
                        // Check for hanging nodes...
                        if (pointEdges_[edgeToCheck[pointI]].size() == 1)
                        {
                            // Null pointEdges so that removeEdge deletes it.
                            pointEdges_[edgeToCheck[pointI]] = labelList(0);

                            // Remove the point
                            points_.remove(edgeToCheck[pointI]);

                            nPoints_--;

                            // Update the reverse point map
                            if (edgeToCheck[pointI] < nOldPoints_)
                            {
                                reversePointMap_[edgeToCheck[pointI]] = -1;
                            }
                        }
                    }

                    removeEdge(faceEdges[edgeI]);
                }
            }

            removeFace(cellToCheck[faceI]);
        }
        else
        if (patch != -1)
        {
            // Check if this internal face is oriented properly.
            face newFace;
            label newOwner = -1;

            if (neighbour_[cellToCheck[faceI]] == index)
            {
                // Orientation is correct
                newFace = faces_[cellToCheck[faceI]];
                newOwner = owner_[cellToCheck[faceI]];
            }
            else
            if (owner_[cellToCheck[faceI]] == index)
            {
                newFace = faces_[cellToCheck[faceI]].reverseFace();
                newOwner = neighbour_[cellToCheck[faceI]];
            }
            else
            {
                // Something's terribly wrong
                FatalErrorIn
                (
                    "dynamicTopoFvMesh::removeCell()"
                )
                    << nl << " Invalid mesh. "
                    << abort(FatalError);
            }

            // Insert a new boundary face
            label newFaceIndex =
                insertFace
                (
                    patch,
                    newFace,
                    newOwner,
                    -1
                );

            // Add the faceEdges entry
            faceEdges_.append(faceEdges);

            // Replace edgeFaces with the new face label
            forAll(faceEdges, edgeI)
            {
                replaceLabel
                (
                    cellToCheck[faceI],
                    newFaceIndex,
                    edgeFaces_[faceEdges[edgeI]]
                );
            }

            // Replace cell with the new face label
            replaceLabel
            (
                cellToCheck[faceI],
                newFaceIndex,
                cells_[newOwner]
            );

            // Remove the internal face.
            removeFace(cellToCheck[faceI]);
        }
    }

    // Update cell info
    cells_.remove(index);
    nCells_--;

    if (index < nOldCells_)
    {
        reverseCellMap_[index] = -1;
    }

    // Check if the cell was added in the current morph, and delete
    if (cellsFromCells_.found(index))
    {
        cellsFromCells_.erase(index);
    }

    // Set the flag
    topoChangeFlag_ = true;
}

// Remove the specified face from the mesh
void dynamicTopoFvMesh::removeFace
(
    const label index
)
{
    if (debug > 2)
    {
        Info << "Removed face: "
             << index << ": "
             << faces_[index] << endl;
    }

    faces_.remove(index);
    owner_.remove(index);
    faceEdges_.remove(index);

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

    // Check and remove from the list of added face patches
    if (addedFacePatches_.found(index))
    {
        addedFacePatches_.erase(index);
    }

    // Decrement the total face-count
    nFaces_--;
}

// Split an internal face into two boundary faces
void dynamicTopoFvMesh::splitInternalFace
(
    const label internalFace
)
{

}

// Merge two boundary faces into one internal face
label dynamicTopoFvMesh::mergeBoundaryFaces
(
    const label firstFace,
    const label secondFace
)
{
    if (debug > 2)
    {
        Info << "Merging faces: "
             << firstFace << " and "
             << secondFace << endl;
    }

    // Sanity check: Are these actually boundary faces?
    if (neighbour_[firstFace] != -1 || neighbour_[secondFace] != -1)
    {
        FatalErrorIn
        (
            "dynamicTopoFvMesh::mergeBoundaryFaces()"
        )
            << nl << " Faces: "
            << firstFace << " and " << secondFace
            << " are not on boundaries. "
            << abort(FatalError);
    }

    // Perform distance-based checks to determine corresponding points
    Map<label> mapPoints;
    face& firstPolyFace = faces_[firstFace];
    face& secondPolyFace = faces_[secondFace];

    forAll(firstPolyFace, pointI)
    {
        forAll(secondPolyFace, pointJ)
        {
            if
            (
                magSqr
                (
                    points_[firstPolyFace[pointI]]
                  - points_[secondPolyFace[pointJ]]
                )
                < VSMALL
            )
            {
                // Found the point. Add it.
                mapPoints.insert(pointI, pointJ);
                break;
            }
        }
    }

    // Sanity check: Were all points matched up?
    if
    (
        (twoDMesh_ && mapPoints.size() != 4) ||
        (!twoDMesh_ && mapPoints.size() != 3)
    )
    {
        FatalErrorIn
        (
            "dynamicTopoFvMesh::mergeBoundaryFaces()"
        )
            << nl << " Faces: "
            << firstFace << " and " << secondFace
            << " do not match up. "
            << abort(FatalError);
    }

    // Obtain owners for both faces, and compare their labels
    face newFace;
    label removedFace = -1, retainedFace = -1;
    label newOwner = -1, newNeighbour = -1;

    if (owner_[firstFace] < owner_[secondFace])
    {
        // Retain the first face
        newFace = firstPolyFace;
        newOwner = owner_[firstFace];
        newNeighbour = owner_[secondFace];
    }
    else
    {
        // Retain the second face
        newFace = secondPolyFace;
        newOwner = owner_[secondFace];
        newNeighbour = owner_[firstFace];
    }

    labelList& faceEdges = faceEdges_[retainedFace];

    // Replace cell with the new face label
    replaceLabel
    (
        removedFace,
        retainedFace,
        cells_[newNeighbour]
    );

    // Remove the boundary face
    removeFace(removedFace);

    // Insert a new interior face
    label newFaceIndex = insertFace(-1, newFace, newOwner, newNeighbour);

    // Add the faceEdges entry
    faceEdges_.append(faceEdges);

    return newFaceIndex;
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

    // Add to the stack as well
    edgeStack(self()).push(newEdgeIndex);

    if (debug > 2)
    {
        Info << "Inserting edge: "
             << newEdgeIndex << ": "
             << newEdge << endl;
    }

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

    if (debug > 2)
    {
        Info << "Removing edge: "
             << index << ": "
             << thisEdge << endl;
    }

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

    // Check and remove from the list of added edge patches
    if (addedEdgePatches_.found(index))
    {
        addedEdgePatches_.erase(index);
    }

    // Decrement the total edge-count
    nEdges_--;
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
void dynamicTopoFvMesh::constructHull
(
    const label eIndex,
    labelList& hullEdges,
    labelList& hullFaces,
    labelList& hullCells,
    labelListList& ringEntities
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

    // Obtain a reference to this edge, and its edgeFaces
    edge& edgeToCheck = edges_[eIndex];
    labelList& eFaces = edgeFaces_[eIndex];
    labelList& hullVertices = edgePoints_[eIndex];

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
    }
}

// Utility to check whether points of an edge lie on a boundary.
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

// Check whether the given edge is on a bounding curve
bool dynamicTopoFvMesh::checkBoundingCurve(label eIndex)
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

// Check triangulation quality for an edge index
bool dynamicTopoFvMesh::checkQuality
(
    const label eIndex,
    const label m,
    const scalarListList& Q,
    const scalar minQuality
)
{
    if (Q[0][m-1] > minQuality)
    {
        return true;
    }

    return false;
}

// Utility method to fill the dynamic programming tables
//  - Returns true if the operation completed successfully.
//  - Returns false if tables could not be resized.
bool dynamicTopoFvMesh::fillTables
(
    const label eIndex,
    label& m,
    scalarListList& Q,
    labelListList& K,
    labelListList& triangulations
)
{
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
                    points_[hullVertices[i]],
                    points_[hullVertices[k]],
                    points_[hullVertices[j]],
                    points_[edgeToCheck[0]]
                );

                scalar qB = (*tetMetric_)
                (
                    points_[hullVertices[j]],
                    points_[hullVertices[k]],
                    points_[hullVertices[i]],
                    points_[edgeToCheck[1]]
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
    // Make a copy of edgePoints, since it will be
    // modified during swaps
    labelList hullVertices(edgePoints_[eIndex]);

    label m = hullVertices.size();
    labelList hullFaces(m, -1);
    labelList hullCells(m, -1);
    labelList hullEdges(m, -1);
    labelListList ringEntities(4, labelList(m, -1));

    // Construct the hull
    constructHull
    (
        eIndex,
        hullEdges,
        hullFaces,
        hullCells,
        ringEntities
    );

    label numTriangulations = 0, isolatedVertex = -1;

    // Extract the appropriate triangulations
    extractTriangulation(0, (m-1), K, numTriangulations, triangulations);

    // Determine the 3-2 swap triangulation
    label t32 = identify32Swap(eIndex, hullVertices, triangulations);

    // Check that the triangulation is valid
    if (t32 == -1)
    {
        // Reset all triangulations and bail out
        triangulations[0] = -1;
        triangulations[1] = -1;
        triangulations[2] = -1;

        return;
    }

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
                        numTriangulations,
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
        numTriangulations,
        triangulations,
        hullVertices,
        hullFaces,
        hullCells
    );

    // Done with this face, so reset it
    triangulations[0][t32] = -1;
    triangulations[1][t32] = -1;
    triangulations[2][t32] = -1;

    // Finally remove the edge
    removeEdge(eIndex);

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    nSwaps_++;
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
// Algorithm taken from:
//   ALGORITHMS TO TEST RAY-TRIANGLE INTERSECTION.
//   R. J. Segura and F. R. Feito,
//   Journal of WSCG, pp. 200-1, 2001.
label dynamicTopoFvMesh::identify32Swap
(
    const label eIndex,
    const labelList& hullVertices,
    const labelListList& triangulations
)
{
    label m = hullVertices.size();
    scalar tolerance = VSMALL;
    edge& edgeToCheck = edges_[eIndex];
    FixedList<label, 3> sign(-2);

    // Relax the tolerance for boundary edges
    if (whichEdgePatch(eIndex) > -1)
    {
        tolerance = 0.01*mag(tangentToEdge(eIndex));
    }

    for (label i = 0; i < (m-2); i++)
    {
        sign[0] =
        (
            tetVolumeSign
            (
                points_[hullVertices[triangulations[0][i]]],
                points_[hullVertices[triangulations[1][i]]],
                points_[edgeToCheck[1]],
                points_[edgeToCheck[0]],
                tolerance
            )
        );

        sign[1] =
        (
            tetVolumeSign
            (
                points_[hullVertices[triangulations[1][i]]],
                points_[hullVertices[triangulations[2][i]]],
                points_[edgeToCheck[1]],
                points_[edgeToCheck[0]],
                tolerance
            )
        );

        sign[2] =
        (
            tetVolumeSign
            (
                points_[hullVertices[triangulations[2][i]]],
                points_[hullVertices[triangulations[0][i]]],
                points_[edgeToCheck[1]],
                points_[edgeToCheck[0]],
                tolerance
            )
        );

        // Intersects at edge AC
        if ((sign[0]==0) && (sign[1]==sign[2]))
        {
            return i;
        }

        // Intersects at edge BC
        if ((sign[1]==0) && (sign[0]==sign[2]))
        {
            return i;
        }

        // Intersects at edge AB
        if ((sign[2]==0) && (sign[0]==sign[1]))
        {
            return i;
        }

        // Intersects inside
        if ((sign[0]==sign[1]) && (sign[1]==sign[2]))
        {
            return i;
        }
    }

    // Could not find an intersecting triangulation
    if (debug > 1)
    {
        Info << "Hull Vertices: " << endl;

        forAll(hullVertices, vertexI)
        {
            Info << hullVertices[vertexI] << ": "
                 << points_[hullVertices[vertexI]]
                 << endl;
        }

        InfoIn("dynamicTopoFvMesh::identify32Swap()") << nl
            << "Could not determine 3-2 swap triangulation." << nl
            << "Edge: " << edgeToCheck << nl
            << "Edge Points: "
            << points_[edgeToCheck[0]] << ","
            << points_[edgeToCheck[1]] << nl
            << abort(FatalError);
    }

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

// Output a list of cells as a VTK file.
// Uses the current state of connectivity.
void dynamicTopoFvMesh::writeVTK
(
    const word& name,
    const labelList& cList
)
{
    label nCells = 0;
    label nTotalCells = 0;
    label nPoints = 0;

    // Estimate a size for points and cellPoints
    List<vector> points(6*cList.size());
    labelListList cpList(cList.size());

    // Create a map for local points
    Map<label> pointMap;

    forAll(cList, cellI)
    {
        if (cList[cellI] == -1)
        {
            continue;
        }

        cell& thisCell = cells_[cList[cellI]];

        // Point-ordering for tetrahedra is different
        if (thisCell.size() == 4)
        {
            const face& currFace = faces_[thisCell[0]];
            const face& nextFace = faces_[thisCell[1]];

            // Size the list
            cpList[nCells].setSize(4);

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
                    // Write-out in order
                    if (owner_[thisCell[0]] == cList[cellI])
                    {
                        cpList[nCells][0] = currFace[2];
                        cpList[nCells][1] = currFace[1];
                        cpList[nCells][2] = currFace[0];
                        cpList[nCells][3] = nextFace[pointI];
                    }
                    else
                    {
                        cpList[nCells][0] = currFace[0];
                        cpList[nCells][1] = currFace[1];
                        cpList[nCells][2] = currFace[2];
                        cpList[nCells][3] = nextFace[pointI];
                    }

                    break;
                }
            }

            // Renumber to local ordering
            forAll(cpList[nCells], pointI)
            {
                // Check if this point was added to the map
                if (!pointMap.found(cpList[nCells][pointI]))
                {
                    // Point was not found, so add it
                    points[nPoints] = points_[cpList[nCells][pointI]];

                    // Update the map
                    pointMap.insert(cpList[nCells][pointI], nPoints);

                    // Increment the number of points
                    nPoints++;
                }

                // Renumber it.
                cpList[nCells][pointI] = pointMap[cpList[nCells][pointI]];
            }

            nTotalCells += 4;
            nCells++;
        }
    }

    // Make the directory
    fileName dirName(time().path()/time().timeName()/"VTK");

    mkDir(dirName);

    // Open stream for output
    OFstream file(dirName/name+".vtk");

    // Write out the header
    file << "# vtk DataFile Version 2.0" << nl
         << name << ".vtk" << nl
         << "ASCII" << nl
         << "DATASET UNSTRUCTURED_GRID" << nl
         << "POINTS " << nPoints << " float" << nl;

    for(label i = 0; i < nPoints; i++)
    {
        file << float(points[i].x()) << ' '
             << float(points[i].y()) << ' '
             << float(points[i].z()) << ' '
             << nl;
    }

    file << "CELLS " << nCells << " " << nTotalCells + nCells << endl;

    forAll(cpList, i)
    {
        if (cList[i] == -1)
        {
            continue;
        }

        file << cpList[i].size() << ' ';

        forAll(cpList[i], j)
        {
            file << cpList[i][j] << ' ';
        }

        file << nl;
    }

    file << "CELL_TYPES " << nCells << endl;
    forAll(cpList, i)
    {
        if (cpList[i].size() == 4)
        {
            // Tetrahedron
            file << "10" << nl;
        }

        if (cpList[i].size() == 6)
        {
            // Wedge
            file << "13" << nl;
        }
    }

    file << nl;
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
                Pout << nl << nl << "Edge: " << faceEdges[edgeI]
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

    label nInternalEdges = 0;
    labelList patchInfo(numPatches_, 0);

    forAllIter(HashList<labelList>::iterator, edgeFaces_, efIter)
    {
        labelList& edgeFaces = efIter();

        if (edgeFaces.size() != nEdgeFaces[efIter.index()])
        {
            Pout << nl << nl << "Edge: " << efIter.index()
                 << "edgeFaces: " << edgeFaces << endl;

            nFailedChecks++;

            WarningIn
            (
                "dynamicTopoFvMesh::checkConnectivity()"
            )
                << nl << "Edge-Face connectivity is inconsistent."
                << endl;
        }

        label nBF = 0;

        // Check if this edge belongs to faceEdges for each face
        forAll(edgeFaces, faceI)
        {
            label i = -1;

            if
            (
                !foundInList
                (
                    efIter.index(), faceEdges_[edgeFaces[faceI]], i
                )
            )
            {
                Pout << nl << nl << "Edge: " << efIter.index()
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

            if (neighbour_[edgeFaces[faceI]] == -1)
            {
                nBF++;
            }
        }

        if (nBF == 0)
        {
            nInternalEdges++;

            // Check if this edge is actually internal.
            if (whichEdgePatch(efIter.index()) >= 0)
            {
                Pout << "Edge: " << efIter.index()
                     << ": " << edges_[efIter.index()] << " is internal, "
                     << " but patch is specified as: "
                     << whichEdgePatch(efIter.index())
                     << endl;

                nFailedChecks++;
            }
        }
        else
        {
            label patchID = whichEdgePatch(efIter.index());

            // Check if this edge is actually on a boundary.
            if (patchID < 0)
            {
                Pout << "Edge: " << efIter.index()
                     << ": " << edges_[efIter.index()]
                     << " is on a boundary, but patch is specified as: "
                     << patchID << endl;

                nFailedChecks++;
            }
            else
            {
                patchInfo[patchID]++;
            }
        }
    }

    if (nInternalEdges != nInternalEdges_)
    {
        Pout << nl << "Internal edge-count is inconsistent." << nl
             << " Counted internal edges: " << nInternalEdges
             << " Actual count: " << nInternalEdges_ << endl;

        nFailedChecks++;
    }

    forAll(patchInfo, patchI)
    {
        if (patchInfo[patchI] != edgePatchSizes_[patchI])
        {
            Pout << "Patch-count is inconsistent." << nl
                 << " Patch: " << patchI
                 << " Counted edges: " << patchInfo[patchI]
                 << " Actual count: " << edgePatchSizes_[patchI] << endl;

            nFailedChecks++;
        }
    }

    // Check added edge patches to ensure that it is consistent
    forAllIter(Map<label>, addedEdgePatches_, aepIter)
    {
        label key = aepIter.key();
        label patch = aepIter();

        label nBF = 0;
        labelList& edgeFaces = edgeFaces_[key];

        // Check if any faces on boundaries
        forAll(edgeFaces, faceI)
        {
            if (neighbour_[edgeFaces[faceI]] == -1)
            {
                nBF++;
            }
        }

        if ((patch < 0) && (nBF > 0))
        {
            Pout << nl << nl << "Edge: " << key
                 << ", edgeFaces: " << edgeFaces
                 << " is internal, but contains boundary faces."
                 << endl;

            nFailedChecks++;
        }

        if ((patch >= 0) && (nBF != 2))
        {
            Pout << nl << nl << "Edge: " << key
                 << ", edgeFaces: " << edgeFaces
                 << " is on a boundary patch, but doesn't contain"
                 << " two boundary faces."
                 << endl;

            nFailedChecks++;
        }
    }

    Info << "Done." << endl;

    if (!twoDMesh_)
    {
        Info << "Checking point-edge connectivity...";

        label allPoints = points_.lastIndex() + 1;
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
                    Pout << nl << nl << "Point: " << peIter.index()
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
                Pout << nl << nl << "Point: " << peIter.index()
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
                Pout << nl << nl
                     << "Edge: " << epIter.index()
                     << " " << edges_[epIter.index()] << endl;

                Pout << "edgeFaces: " << edgeFaces << endl;
                forAll(edgeFaces, faceI)
                {
                    Info << edgeFaces[faceI] << ": "
                         << faces_[edgeFaces[faceI]]
                         << endl;
                }

                Pout << "edgePoints: " << edgePoints << endl;

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

                label i = -1;

                if (!foundInList(otherPoint, edgePoints, i))
                {
                    Pout << nl << nl
                         << "Edge: " << epIter.index()
                         << " " << edges_[epIter.index()] << endl;

                    Pout << "edgeFaces: " << edgeFaces << endl;
                    forAll(edgeFaces, faceI)
                    {
                        Info << edgeFaces[faceI] << ": "
                             << faces_[edgeFaces[faceI]]
                             << endl;
                    }

                    Pout << "edgePoints: " << edgePoints << endl;

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
                Pout << nl << "Warning: Cell: "
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

    reduce(nFailedChecks, orOp<bool>());

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
        wordList toc = fixedPatches_.toc();

        // Obtain the list of patches for which
        // curvature-based length-scale is specified
        wordList cToc = curvaturePatches_.toc();

        // Do a preliminary sanity check to avoid duplication
        forAll(toc, wordI)
        {
            word& pName = toc[wordI];

            forAll(cToc, wordI)
            {
                if (pName == cToc[wordI])
                {
                    FatalErrorIn
                    (
                        "dynamicTopoFvMesh::calculateLengthScale()"
                    )
                        << " Conflicting fixed length-scale patch: "
                        << pName << abort(FatalError);
                }
            }
        }

        // Loop through all boundaries and mark adjacent cells
        const polyBoundaryMesh& boundary = boundaryMesh();
        const labelList& own = faceOwner();
        const pointField& pList = points();

        forAll(boundary,patchI)
        {
            const polyPatch& bdyPatch = boundary[patchI];

            // Loop through all fixed curvature patches
            forAll(cToc, wordI)
            {
                word& pName = cToc[wordI];

                if (boundary[patchI].name() == pName)
                {
                    label pStart = bdyPatch.start();

                    if (curvatureFields_[patchI].empty())
                    {
                        FatalErrorIn
                        (
                            "dynamicTopoFvMesh::calculateLengthScale()"
                        )
                            << " Curvature field for patch: "
                            << pName << " is empty."
                            << abort(FatalError);
                    }

                    forAll(bdyPatch,faceI)
                    {
                        label ownCell = own[pStart+faceI];

                        if (cellLevels[ownCell] != 0)
                        {
                            continue;
                        }

                        cellLevels[ownCell] = level;

                        lengthScale[ownCell] =
                        (
                            curvatureRatio_
                           *curvatureFields_[patchI][faceI]
                           *growthFactor_
                        );

                        levelCells.insert(ownCell);

                        visitedCells++;
                    }

                    break;
                }
            }

            // Loop through all fixed length-scale patches
            forAll(toc,wordI)
            {
                word& pName = toc[wordI];

                if (boundary[patchI].name() == pName)
                {
                    label pStart = bdyPatch.start();

                    forAll(bdyPatch,faceI)
                    {
                        label ownCell = own[pStart+faceI];

                        if (cellLevels[ownCell] != 0)
                        {
                            continue;
                        }

                        cellLevels[ownCell] = level;

                        lengthScale[ownCell] =
                        (
                            fixedPatches_[pName][0].scalarToken()
                           *growthFactor_
                        );

                        levelCells.insert(ownCell);

                        visitedCells++;
                    }

                    break;
                }
            }

            // Set boundary patch size if no fixed-length scale is specified.
            if
            (
                (toc.size() == 0) &&
                (cToc.size() == 0) &&
                (bdyPatch.type() != "processor") &&
                (bdyPatch.type() != "wedge") &&
                (bdyPatch.type() != "empty") &&
                (bdyPatch.type() != "symmetryPlane")
            )
            {
                label pStart = bdyPatch.start();

                forAll(bdyPatch,faceI)
                {
                    label ownCell = own[pStart+faceI];

                    if (cellLevels[ownCell] != 0)
                    {
                        continue;
                    }

                    cellLevels[ownCell] = level;

                    if (twoDMesh_)
                    {
                        label eIndex = getTriBoundaryEdge(pStart+faceI);
                        edge& e = edges_[eIndex];

                        lengthScale[ownCell] =
                        (
                            mag(pList[e[0]] - pList[e[1]])*growthFactor_
                        );
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

                        lengthScale[ownCell] =
                        (
                            (edgeLength/fEdges.size())*growthFactor_
                        );
                    }

                    levelCells.insert(ownCell);

                    visitedCells++;
                }
            }
        }

        bool doneWithSweeps = false;

        // Perform multiple sweeps through the mesh...
        while (!doneWithSweeps)
        {
            if (Pstream::parRun())
            {
                writeLengthScaleInfo
                (
                    cellLevels,
                    lengthScale
                );
            }

            // Loop through cells of the current level
            labelList currLvlCells = levelCells.toc();
            levelCells.clear();

            // Loop through cells, and increment neighbour
            // cells of the current level
            forAll(currLvlCells,cellI)
            {
                // Obtain the cells neighbouring this one
                const labelList& cList = cc[currLvlCells[cellI]];

                forAll(cList, indexI)
                {
                    label& ngbLevel = cellLevels[cList[indexI]];

                    if (ngbLevel == 0)
                    {
                        ngbLevel = level + 1;

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

            if (Pstream::parRun())
            {
                readLengthScaleInfo
                (
                    level,
                    visitedCells,
                    cellLevels,
                    lengthScale,
                    levelCells
                );
            }

            if (debug > 2)
            {
                Pout << "Processed level: " << level << nl
                     << " Visited: " << visitedCells
                     << " out of " << nCells() << endl;
            }

            // Move on to the next level
            level++;

            if (visitedCells >= nCells())
            {
                doneWithSweeps = true;
            }

            // Wait for everyone to complete.
            reduce(doneWithSweeps, andOp<bool>());
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

// Send length-scale info across processors
void dynamicTopoFvMesh::writeLengthScaleInfo
(
    const labelList& cellLevels,
    const scalarList& lengthScale
)
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    // Clear existing buffers
    sendLblBuffer_.clear();
    recvLblBuffer_.clear();
    sendSclBuffer_.clear();
    recvSclBuffer_.clear();

    // Number of sent faces across processors
    labelList nSendFaces(boundary.size(), 0);
    labelList nRecvFaces(boundary.size(), 0);

    // Corresponding face labels
    sendLblBuffer_.setSize(boundary.size());
    recvLblBuffer_.setSize(boundary.size());

    // Length-scales corresponding to face-labels
    sendSclBuffer_.setSize(boundary.size());
    recvSclBuffer_.setSize(boundary.size());

    // Fill send buffers with cell-level and length-scale info.
    forAll(boundary, pI)
    {
        if (isA<processorPolyPatch>(boundary[pI]))
        {
            const labelList& fCells = boundary[pI].faceCells();

            // Set the initial buffer size
            sendLblBuffer_[pI].setSize(boundary[pI].size(), 0);
            sendSclBuffer_[pI].setSize(boundary[pI].size(), 0.0);

            forAll(fCells, faceI)
            {
                label cI = fCells[faceI];

                // Does the adjacent cell have a non-zero level?
                if (cellLevels[cI] > 0)
                {
                    // Fill the send buffer.
                    sendLblBuffer_[pI][nSendFaces[pI]] = faceI;
                    sendSclBuffer_[pI][nSendFaces[pI]] = lengthScale[cI];

                    nSendFaces[pI]++;
                }
            }

            // Resize to actual value
            sendLblBuffer_[pI].setSize(nSendFaces[pI]);
            sendSclBuffer_[pI].setSize(nSendFaces[pI]);

            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[pI])
            );

            label neiProcNo = pp.neighbProcNo();

            // First perform a blocking send/receive of the number of faces.
            pWrite(neiProcNo, nSendFaces[pI]);
            pRead(neiProcNo, nRecvFaces[pI]);
        }
    }

    // Send info to neighbouring processors.
    forAll(boundary, patchI)
    {
        if (isA<processorPolyPatch>(boundary[patchI]))
        {
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[patchI])
            );

            label neiProcNo = pp.neighbProcNo();

            if (debug > 3)
            {
                Pout << " Processor patch " << patchI << ' ' << pp.name()
                     << " communicating with " << neiProcNo
                     << "  Sending: " << nSendFaces[patchI]
                     << endl;
            }

            // Next, perform a non-blocking send/receive of buffers.
            // But only if the buffer size is non-zero.
            if (nSendFaces[patchI] != 0)
            {
                pWrite(neiProcNo, sendLblBuffer_[patchI]);
                pWrite(neiProcNo, sendSclBuffer_[patchI]);
            }

            if (nRecvFaces[patchI] != 0)
            {
                // Size the receive buffers
                recvLblBuffer_[patchI].setSize(nRecvFaces[patchI], 0);
                recvSclBuffer_[patchI].setSize(nRecvFaces[patchI], 0.0);

                pRead(neiProcNo, recvLblBuffer_[patchI]);
                pRead(neiProcNo, recvSclBuffer_[patchI]);
            }
        }
    }
}

// Receive length-scale info across processors
void dynamicTopoFvMesh::readLengthScaleInfo
(
    const label level,
    label& visitedCells,
    labelList& cellLevels,
    scalarList& lengthScale,
    labelHashSet& levelCells
)
{
    const polyBoundaryMesh& boundary = boundaryMesh();
    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();

    // Wait for all transfers to complete.
    OPstream::waitRequests();
    IPstream::waitRequests();

    // Now re-visit cells and update length-scales.
    forAll(boundary, patchI)
    {
        if (recvLblBuffer_[patchI].size())
        {
            const labelList& fCells = boundary[patchI].faceCells();

            forAll(recvLblBuffer_[patchI], i)
            {
                label nLabel = recvLblBuffer_[patchI][i];

                label cI = fCells[nLabel];

                label& ngbLevel = cellLevels[cI];

                // For an unvisited cell, update the level
                if (ngbLevel == 0)
                {
                    ngbLevel = level + 1;
                    levelCells.insert(cI);
                    visitedCells++;
                }

                // Visit all neighbours and re-calculate length-scale
                const cell& cellCheck = cells()[cI];
                scalar sumLength = 0.0;
                label nTouchedNgb = 0;

                forAll(cellCheck, faceI)
                {
                    label sLevel = -1, sLCell = -1;
                    label pF = boundary.whichPatch(cellCheck[faceI]);

                    if (pF == -1)
                    {
                        // Internal face. Determine neighbour.
                        if (own[cellCheck[faceI]] == cI)
                        {
                            sLCell = nei[cellCheck[faceI]];
                        }
                        else
                        {
                            sLCell = own[cellCheck[faceI]];
                        }

                        sLevel = cellLevels[sLCell];

                        if ((sLevel < ngbLevel) && (sLevel > 0))
                        {
                            sumLength += lengthScale[sLCell];

                            nTouchedNgb++;
                        }
                    }
                    else
                    if (isA<processorPolyPatch>(boundary[pF]))
                    {
                        // Determine the local index.
                        label local = boundary[pF].whichFace(cellCheck[faceI]);

                        // Is this label present in the list?
                        forAll(recvLblBuffer_[pF], j)
                        {
                            if (recvLblBuffer_[pF][j] == local)
                            {
                                sumLength += recvSclBuffer_[pF][j];

                                nTouchedNgb++;

                                break;
                            }
                        }
                    }
                    else
                    if (fixedPatches_.found(boundary[pF].name()))
                    {
                        sumLength +=
                        (
                            fixedPatches_[boundary[pF].name()][0].scalarToken()
                        );

                        nTouchedNgb++;
                    }
                }

                sumLength /= nTouchedNgb;

                // Scale the length and assign to this cell
                scalar sLength = sumLength*growthFactor_;

                sLength = (sLength < maxLengthScale_)
                        ? sLength : maxLengthScale_;

                lengthScale[cI] = sLength;
            }
        }
    }
}

// Parallel blocking send
void dynamicTopoFvMesh::pWrite
(
    const label toID,
    const label& data
)
{
    OPstream::write
    (
        Pstream::blocking,
        toID,
        reinterpret_cast<const char*>
        (
            &data
        ),
        sizeof(label)
    );
}

// Parallel blocking receive
void dynamicTopoFvMesh::pRead
(
    const label fromID,
    label& data
)
{
    IPstream::read
    (
        Pstream::blocking,
        fromID,
        reinterpret_cast<char*>
        (
            &data
        ),
        sizeof(label)
    );
}

// Parallel non-blocking send for lists
template <class Type>
void dynamicTopoFvMesh::pWrite
(
    const label toID,
    const List<Type>& data
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
void dynamicTopoFvMesh::pRead
(
    const label fromID,
    List<Type>& data
)
{
    IPstream::read
    (
        Pstream::nonBlocking,
        fromID,
        reinterpret_cast<char*>(data.begin()),
        data.byteSize()
    );
}

// Read edge-modification options from the dictionary
void dynamicTopoFvMesh::readEdgeOptions()
{
    const dictionary& edgeOptionDict =
    (
        dict_.subDict("dynamicTopoFvMesh").subDict("edgeOptions")
    );

    ratioMax_ = readScalar(edgeOptionDict.lookup("bisectionRatio"));
    ratioMin_ = readScalar(edgeOptionDict.lookup("collapseRatio"));
    growthFactor_ = readScalar(edgeOptionDict.lookup("growthFactor"));

    if (edgeOptionDict.found("maxLengthScale"))
    {
        maxLengthScale_ =
        (
            readScalar(edgeOptionDict.lookup("maxLengthScale"))
        );
    }

    if (edgeOptionDict.found("fixedLengthScalePatches"))
    {
        fixedPatches_ =
        (
            edgeOptionDict.subDict("fixedLengthScalePatches")
        );
    }

    if (edgeOptionDict.found("curvaturePatches"))
    {
        curvaturePatches_ =
        (
            edgeOptionDict.subDict("curvaturePatches")
        );

        curvatureRatio_ =
        (
            readScalar(edgeOptionDict.lookup("curvatureRatio"))
        );
    }

    if (edgeOptionDict.found("sliverThreshold"))
    {
        sliverThreshold_ =
        (
            readScalar(edgeOptionDict.lookup("sliverThreshold"))
        );
    }

    if (edgeOptionDict.found("maxModifications"))
    {
        maxModifications_ =
        (
            readLabel(edgeOptionDict.lookup("maxModifications"))
        );
    }
}

// Return the appropriate length-scale for boundary face
scalar dynamicTopoFvMesh::boundaryLengthScale
(
    const label faceIndex
)
{
    label bFacePatch = whichPatch(faceIndex);

    // Check curvature patches
    if (curvaturePatches_.found(boundaryMesh()[bFacePatch].name()))
    {
        // Find the local index in the patch
        label lIndex = -1;

        if (faceIndex < nOldFaces_)
        {
            lIndex = faceIndex - boundaryMesh()[bFacePatch].start();
        }
        else
        {
            lIndex = faceParents_[faceIndex];
        }

        return (curvatureRatio_*curvatureFields_[bFacePatch][lIndex]);
    }

    // Check fixed length-scale patches
    if (fixedPatches_.found(boundaryMesh()[bFacePatch].name()))
    {
        return
        (
            fixedPatches_[boundaryMesh()[bFacePatch].name()][0].scalarToken()
        );
    }

    return lengthScale_[owner_[faceIndex]];
}

// Compute the face-curvature for a given boundary face
// For 3D simplical meshes only
scalar dynamicTopoFvMesh::faceCurvature
(
    const label fIndex
)
{
    label nf = 0;
    scalar curvature = 0.0, kAvg = 0.0;
    vector te = vector::zero, m = vector::zero;
    FixedList<scalar, 2> k(0.0);
    FixedList<label, 2> ebFaces(-1);
    FixedList<point, 2> r, rcg, e;
    labelList& fEdges = faceEdges_[fIndex];

    // Loop through edges of this face and compute
    forAll(fEdges, edgeI)
    {
        label eIndex = fEdges[edgeI];
        labelList& eFaces = edgeFaces_[eIndex];

        nf = 0; ebFaces = -1;

        // Find the two boundary faces for this edge
        forAll(eFaces, faceI)
        {
            if (neighbour_[eFaces[faceI]] == -1)
            {
                ebFaces[nf++] = eFaces[faceI];
            }
        }

        // Compute the circumcenters
        r[0] = circumCenter(ebFaces[0]);
        r[1] = circumCenter(ebFaces[1]);

        // Compute face centers
        rcg[0] = triFaceCenter(ebFaces[0]);
        rcg[1] = triFaceCenter(ebFaces[1]);

        // Obtain tangent-to-edge
        te = tangentToEdge(eIndex);

        // Compute face tangential vectors
        e[0] = (rcg[0] - ((te&rcg[0])*te))/(te&te);
        e[1] = (rcg[1] - ((te&rcg[1])*te))/(te&te);

        // Normalize them...
        e[0] /= mag(e[0]);
        e[1] /= mag(e[1]);

        // Now compute the binormal vector
        m = ((r[1]&e[1])*e[0]) + ((r[0]&e[0])*e[1]);

        // Compute the correction factors for this edge
        k[0] =
            Foam::sqrt
            (
                1 + (1 - ((m&e[0])/(m&m)))
             * ((te&te)/(4*(r[0]&r[0])))
            );

        k[1] =
            Foam::sqrt
            (
                1 + (1 - ((m&e[1])/(m&m)))
             * ((te&te)/(4*(r[1]&r[1])))
            );

        kAvg = 0.5*(k[0] + k[1]);
    }

    return curvature;
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

    // If this is a slave, acquire the lock for this structure
    if (thread->slave())
    {
        thread->lock();
    }

    dynamicTopoFvMesh& mesh = thread->mesh();

    // Figure out which thread this is...
    label tIndex = mesh.self();

    bool failed = false;

    // Pick items off the stack
    while (!mesh.faceStack(tIndex).empty())
    {
        // Retrieve the index for this face
        label fIndex = mesh.faceStack(tIndex).pop();

        // Perform a Delaunay test and check if a flip is necesary.
        mesh.testDelaunay
        (
            fIndex,
            failed
        );

        if (failed)
        {
            if (thread->master())
            {
                // Swap this face.
                mesh.swapQuadFace(fIndex);
            }
            else
            {
                // Push this on to the master stack
                mesh.faceStack(0).push(fIndex);
            }
        }
    }

    // If this is a slave, relinquish the lock for this structure
    if (thread->slave())
    {
        thread->unlock();
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

    const polyBoundaryMesh& boundary = boundaryMesh();

    // Check if patches are explicitly coupled
    if (dict_.subDict("dynamicTopoFvMesh").found("coupledPatches"))
    {
        dictionary coupledPatches =
            dict_.subDict("dynamicTopoFvMesh").subDict("coupledPatches");

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

            if (mPatch == -1 && sPatch == -1)
            {
                // This pair doesn't exist. This might be
                // true for some sub-domains.
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
                FatalErrorIn("dynamicTopoFvMesh::initEdges()")
                        << " Coupled patches are wrongly specified." << nl
                        << " Master: " << mPatch << ":" << masterPatch << nl
                        << " Slave: " << sPatch << ":" << slavePatch << nl
                        << abort(FatalError);
            }
        }
    }

    // Add processor patches to the list
    forAll(boundary, patchI)
    {
        if (isA<processorPolyPatch>(boundary[patchI]))
        {
            patchCoupling_.insert(patchI, -1);
        }
    }
}

// Return a reference to the multiThreader
const multiThreader& dynamicTopoFvMesh::threader() const
{
    return threader_();
}

// Does the mesh perform edge-modification?
bool dynamicTopoFvMesh::edgeModification()
{
    return edgeModification_;
}

// Send and receive patchSubMeshes for coupled patches
void dynamicTopoFvMesh::sendAndRecvCoupledMeshes()
{
    if (!patchCoupling_.size())
    {
        return;
    }

    // Build patch sub-meshes (and clear existing ones).
    buildPatchSubMeshes();

    // Build maps for locally coupled patches.
    buildCoupledMaps(true);

    // Wait for all transfers to complete.
    waitForBuffers();

    // Build maps for coupled processor patches.
    buildCoupledMaps(false);

    // Prepare the master stack
    if (!twoDMesh_)
    {
        forAll(masterToSlave_, patchI)
        {
            forAllIter(Map<label>::iterator, masterToSlave_[patchI], indexI)
            {
                edgeStack(0).push(indexI.key());
            }
        }
    }
}

// Handle topology changes for coupled patches
void dynamicTopoFvMesh::handleCoupledPatches()
{
    if (!patchCoupling_.size())
    {
        return;
    }

    // Set coupled modifications.
    structPtr_[0].isCoupled() = true;

    // Loop through the coupled stack and perform changes.
    if (twoDMesh_)
    {
        if (edgeModification_)
        {
            edgeBisectCollapse2D(&(structPtr_[0]));
        }

        swap2DEdges(&(structPtr_[0]));
    }
    else
    {
        if (edgeModification_)
        {
            edgeBisectCollapse3D(&(structPtr_[0]));
        }

        swap3DEdges(&(structPtr_[0]));
    }

    // Reset coupled modifications.
    structPtr_[0].isCoupled() = false;
}

// Build patch sub-meshes for coupled patches
void dynamicTopoFvMesh::buildPatchSubMeshes()
{
    if (debug)
    {
        Info << "Building patch sub-meshes for coupled patches...";
    }

    const polyBoundaryMesh& boundary = boundaryMesh();

    // Size the PtrLists.
    sendPatchMeshes_.clear();
    recvPatchMeshes_.clear();

    // Clear the maps
    masterToSlave_.clear();
    slaveToMaster_.clear();

    sendPatchMeshes_.setSize(boundary.size());
    recvPatchMeshes_.setSize(boundary.size());

    masterToSlave_.setSize(boundary.size());
    slaveToMaster_.setSize(boundary.size());

    forAllIter(Map<label>, patchCoupling_, patchI)
    {
        label mPatch = patchI.key();
        label sPatch = patchI();

        if (sPatch == -1)
        {
            // This is a processor patch.

            // Find out which neighbour it talks to.
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[mPatch])
            );

            label neiProc = pp.neighbProcNo();

            if (neiProc < Pstream::myProcNo())
            {
                sendPatchMeshes_.set(mPatch, new patchSubMesh());

                buildPatchSubMesh(mPatch, sendPatchMeshes_[mPatch]);

                // Send my sub-mesh to the neighbour.
                pWrite(neiProc, sendPatchMeshes_[mPatch].nPoints());
                pWrite(neiProc, sendPatchMeshes_[mPatch].nEdges());
                pWrite(neiProc, sendPatchMeshes_[mPatch].nFaces());
                pWrite(neiProc, sendPatchMeshes_[mPatch].nCells());

                pWrite(neiProc, sendPatchMeshes_[mPatch].pointBuffer());
                pWrite(neiProc, sendPatchMeshes_[mPatch].edgeBuffer());
                pWrite(neiProc, sendPatchMeshes_[mPatch].faceBuffer());
                pWrite(neiProc, sendPatchMeshes_[mPatch].cellBuffer());

                if (edgeModification_)
                {
                    pWrite(neiProc, sendPatchMeshes_[mPatch].lengthBuffer());
                }
            }
            else
            {
                recvPatchMeshes_.set(mPatch, new patchSubMesh());

                // First read entity sizes.
                pRead(neiProc, recvPatchMeshes_[mPatch].nPoints());
                pRead(neiProc, recvPatchMeshes_[mPatch].nEdges());
                pRead(neiProc, recvPatchMeshes_[mPatch].nFaces());
                pRead(neiProc, recvPatchMeshes_[mPatch].nCells());

                // Obtain references.
                pointField& pBuffer = recvPatchMeshes_[mPatch].pointBuffer();
                labelList& eBuffer = recvPatchMeshes_[mPatch].edgeBuffer();
                labelList& fBuffer = recvPatchMeshes_[mPatch].faceBuffer();
                labelList& cBuffer = recvPatchMeshes_[mPatch].cellBuffer();

                if (!twoDMesh_)
                {
                    // Size the buffers.
                    pBuffer.setSize(recvPatchMeshes_[mPatch].nPoints());
                    eBuffer.setSize(2*recvPatchMeshes_[mPatch].nEdges());
                    fBuffer.setSize(3*recvPatchMeshes_[mPatch].nFaces());
                    cBuffer.setSize(4*recvPatchMeshes_[mPatch].nCells());

                    // Receive buffers
                    pRead(neiProc, pBuffer);
                    pRead(neiProc, eBuffer);
                    pRead(neiProc, fBuffer);
                    pRead(neiProc, cBuffer);
                }

                if (edgeModification_)
                {
                    scalarList& lB = recvPatchMeshes_[mPatch].lengthBuffer();
                    lB.setSize(recvPatchMeshes_[mPatch].nCells());
                    pRead(neiProc, lB);
                }
            }
        }
    }

    if (debug)
    {
        Info << "Done.";
    }
}

// Build patch sub-mesh for a specified coupled patch
void dynamicTopoFvMesh::buildPatchSubMesh
(
    const label patchID,
    patchSubMesh& subMesh
)
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    // Obtain mesh information
    label nP = 0, nE = 0, nF = 0, nC = 0;
    const labelListList& cCells = cellCells();
    const labelList& fCells = boundary[patchID].faceCells();

    // Obtain references
    Map<label>& rPointMap = subMesh.reversePointMap();
    Map<label>& rEdgeMap = subMesh.reverseEdgeMap();
    Map<label>& rFaceMap = subMesh.reverseFaceMap();
    Map<label>& rCellMap = subMesh.reverseCellMap();

    Map<label>& pointMap = subMesh.pointMap();
    Map<label>& edgeMap = subMesh.edgeMap();
    Map<label>& faceMap = subMesh.faceMap();
    Map<label>& cellMap = subMesh.cellMap();

    // Allocate the cellMap
    forAll(fCells, faceI)
    {
        cellMap.insert(nC, fCells[faceI]);
        rCellMap.insert(fCells[faceI], nC);
        nC++;
    }

    forAllIter(Map<label>::iterator, rCellMap, cIter)
    {
        const labelList& cList = cCells[cIter.key()];

        forAll(cList, cellI)
        {
            if (!rCellMap.found(cList[cellI]))
            {
                cellMap.insert(nC, cList[cellI]);
                rCellMap.insert(cList[cellI], nC);
                nC++;
            }
        }
    }

    // Allocate the faceMap
    forAllIter(Map<label>::iterator, rCellMap, cIter)
    {
        const cell& thisCell = cells_[cIter.key()];

        forAll(thisCell, faceI)
        {
            if (!rFaceMap.found(thisCell[faceI]))
            {
                faceMap.insert(nF, thisCell[faceI]);
                rFaceMap.insert(thisCell[faceI], nF);
                nF++;
            }
        }
    }

    // Allocate the edgeMap
    forAllIter(Map<label>::iterator, rFaceMap, fIter)
    {
        const labelList& fEdges = faceEdges_[fIter.key()];

        forAll(fEdges, edgeI)
        {
            if (!rEdgeMap.found(fEdges[edgeI]))
            {
                edgeMap.insert(nE, fEdges[edgeI]);
                rEdgeMap.insert(fEdges[edgeI], nE);
                nE++;
            }
        }
    }

    // Allocate the pointMap
    forAllIter(Map<label>::iterator, rEdgeMap, eIter)
    {
        const edge& thisEdge = edges_[eIter.key()];

        if (!rPointMap.found(thisEdge[0]))
        {
            pointMap.insert(nP, thisEdge[0]);
            rPointMap.insert(thisEdge[0], nP);
            nP++;
        }

        if (!rPointMap.found(thisEdge[1]))
        {
            pointMap.insert(nP, thisEdge[1]);
            rPointMap.insert(thisEdge[1], nP);
            nP++;
        }
    }

    // Assign sizes to the mesh
    subMesh.nPoints() = nP;
    subMesh.nEdges() = nE;
    subMesh.nFaces() = nF;
    subMesh.nCells() = nC;

    // Size up buffers and fill them
    pointField& pBuffer = subMesh.pointBuffer();
    pBuffer.setSize(subMesh.nPoints(), vector::zero);

    forAllIter(Map<label>::iterator, pointMap, pIter)
    {
        pBuffer[pIter.key()] = points_[pIter()];
    }

    // Edge buffer size: 2 points for every edge
    labelList& eBuffer = subMesh.edgeBuffer();
    eBuffer.setSize(2*subMesh.nEdges(), -1);

    label index = 0;

    forAllIter(Map<label>::iterator, edgeMap, eIter)
    {
        edge& edgeToCheck = edges_[eIter()];

        eBuffer[index++] = rPointMap[edgeToCheck[0]];
        eBuffer[index++] = rPointMap[edgeToCheck[1]];
    }

    labelList& fBuffer = subMesh.faceBuffer();
    labelList& cBuffer = subMesh.cellBuffer();

    if (!twoDMesh_)
    {
        // Face buffer size: 3 edges for every face in 3D
        fBuffer.setSize(3*subMesh.nFaces(), -1);

        index = 0;

        forAllIter(Map<label>::iterator, faceMap, fIter)
        {
            labelList& fEdges = faceEdges_[fIter()];

            fBuffer[index++] = rEdgeMap[fEdges[0]];
            fBuffer[index++] = rEdgeMap[fEdges[1]];
            fBuffer[index++] = rEdgeMap[fEdges[2]];
        }

        // Cell buffer size: 4 faces for every cell in 3D
        cBuffer.setSize(4*subMesh.nCells(), -1);

        index = 0;

        forAllIter(Map<label>::iterator, cellMap, cIter)
        {
            cell& cellToCheck = cells_[cIter()];

            cBuffer[index++] = rFaceMap[cellToCheck[0]];
            cBuffer[index++] = rFaceMap[cellToCheck[1]];
            cBuffer[index++] = rFaceMap[cellToCheck[2]];
            cBuffer[index++] = rFaceMap[cellToCheck[3]];
        }
    }

    if (edgeModification_)
    {
        scalarList& lBuffer = subMesh.lengthBuffer();

        lBuffer.setSize(subMesh.nCells(), 0.0);

        forAllIter(Map<label>::iterator, cellMap, cIter)
        {
            lBuffer[cIter.key()] = lengthScale_[cIter()];
        }
    }
}

// Build coupled maps
void dynamicTopoFvMesh::buildCoupledMaps
(
    bool localOnly
)
{
    if (!patchCoupling_.size())
    {
        return;
    }

    const polyBoundaryMesh& boundary = boundaryMesh();

    label neiProc = -1, start = -1;
    labelHashSet mList, sList;

    forAllIter(Map<label>, patchCoupling_, patchI)
    {
        if (patchI() == -1)
        {
            // Are we doing only local coupled patches?
            if (localOnly)
            {
                continue;
            }
            else
            {
                // Is this patch a processor slave?
                // If it is, skip it.
                const processorPolyPatch& pp =
                (
                    refCast<const processorPolyPatch>
                    (
                        boundary[patchI.key()]
                    )
                );

                neiProc = pp.neighbProcNo();

                if (neiProc < Pstream::myProcNo())
                {
                    continue;
                }
            }
        }

        // Build a list of master edges
        start = boundary[patchI.key()].start();

        for (label i = 0; i < boundary[patchI.key()].size(); i++)
        {
            labelList& mfEdges = faceEdges_[start+i];

            forAll(mfEdges, edgeI)
            {
                if (!mList.found(mfEdges[edgeI]))
                {
                    mList.insert(mfEdges[edgeI]);
                }
            }
        }

        // Build a list of edge-centres for geometric matching.
        labelList mEdges = mList.toc();
        pointField mCentres(mEdges.size(), vector::zero);

        forAll(mEdges, edgeI)
        {
            const edge& mE = edges_[mEdges[edgeI]];
            mCentres[edgeI] = 0.5*(points_[mE[0]] + points_[mE[1]]);
        }

        // Empty slave lists. Sizes are set depending on patch type.
        labelList sEdges(0);
        pointField sCentres(0);

        // Build a list of slave edges
        if (patchI() == -1)
        {
            // Processor patch.
            // Build a list from the patchSubMesh.

            // If any edges here coincide with locally coupled edges,
            // processor edges need to be given priority.
        }
        else
        {
            // Locally coupled patch.
            // Build a list from local edges.
            start = boundary[patchI()].start();

            for (label i = 0; i < boundary[patchI()].size(); i++)
            {
                labelList& mfEdges = faceEdges_[start+i];

                forAll(mfEdges, edgeI)
                {
                    if (!sList.found(mfEdges[edgeI]))
                    {
                        sList.insert(mfEdges[edgeI]);
                    }
                }
            }

            // Sanity check: Do patches have equal number of edges?
            if (mList.size() != sList.size())
            {
                FatalErrorIn("dynamicTopoFvMesh::buildCoupledMaps()")
                    << "Patch edge sizes are not consistent."
                    << abort(FatalError);
            }

            // Build a list of edge-centres for geometric matching.
            sEdges = sList.toc();
            sCentres.setSize(sEdges.size(), vector::zero);

            forAll(mEdges, edgeI)
            {
                const edge& sE = edges_[sEdges[edgeI]];
                sCentres[edgeI] = 0.5*(points_[sE[0]] + points_[sE[1]]);
            }
        }

        label nMatchedEdges = 0;

        forAll(mCentres, edgeI)
        {
            forAll(sCentres, edgeJ)
            {
                if (mag(mCentres[edgeI] - sCentres[edgeJ]) < SMALL)
                {
                    // Add a map entry
                    masterToSlave_[patchI.key()].insert
                    (
                        mEdges[edgeI],
                        sEdges[edgeJ]
                    );

                    slaveToMaster_[patchI.key()].insert
                    (
                        sEdges[edgeJ],
                        mEdges[edgeI]
                    );

                    nMatchedEdges++;

                    break;
                }
            }
        }

        // Make sure we were successful.
        if (nMatchedEdges != mCentres.size())
        {
            FatalErrorIn("dynamicTopoFvMesh::buildCoupledMaps()")
                << "Failed to match all patch edges."
                << abort(FatalError);
        }

        // Clear the HashSets
        mList.clear();
        sList.clear();
    }
}

// Wait for buffer transfer completion.
void dynamicTopoFvMesh::waitForBuffers()
{
    if (!patchCoupling_.size())
    {
        return;
    }

    OPstream::waitRequests();
    IPstream::waitRequests();
}

// Synchronize and exit for parallel runs.
void dynamicTopoFvMesh::synchronizeAndExit()
{
    if (Pstream::parRun())
    {
        Info << "Synchronizing all processes...";

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();

        Info << "Done." << endl;
        Info << "Terminating normally." << endl;
    }

    ::exit(0);
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

    // If this is a slave, acquire the lock for this structure
    if (thread->slave())
    {
        thread->lock();
    }

    dynamicTopoFvMesh& mesh = thread->mesh();

    // Figure out which thread this is...
    label tIndex = mesh.self();
    scalar length = 0.0, scale = 0.0;

    // Pick items off the stack
    while (!mesh.faceStack(tIndex).empty())
    {
        // Retrieve the index for this face
        label fIndex = mesh.faceStack(tIndex).pop();

        // Select only quad-faces
        if (mesh.checkQuadFace(fIndex))
        {
            // Measure the boundary edge-length of the face in question
            mesh.edgeLength
            (
                mesh.getTriBoundaryEdge(fIndex),
                length
            );

            // Determine the length-scale at this face
            mesh.meshFaceLengthScale(fIndex, scale);

            // Check if this boundary face is adjacent to a sliver-cell,
            // and remove it by a two-step bisection/collapse operation.
            if (mesh.whichPatch(fIndex) != -1)
            {
                mesh.remove2DSliver(fIndex);
            }

            if (length > mesh.ratioMax()*scale)
            {
                if (thread->master())
                {
                    // Bisect this face
                    mesh.bisectQuadFace(fIndex);
                }
                else
                {
                    // Push this on to the master stack
                    mesh.faceStack(0).push(fIndex);
                }
            }
            else
            if (length < mesh.ratioMin()*scale)
            {
                if (thread->master())
                {
                    // Collapse this face
                    mesh.collapseQuadFace(fIndex);
                }
                else
                {
                    // Push this on to the master stack
                    mesh.faceStack(0).push(fIndex);
                }
            }
        }
    }

    // If this is a slave, relinquish the lock for this structure
    if (thread->slave())
    {
        thread->unlock();
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

    // If this is a slave, acquire the lock for this structure
    if (thread->slave())
    {
        thread->lock();
    }

    dynamicTopoFvMesh& mesh = thread->mesh();

    // Figure out which thread this is...
    label tIndex = mesh.self();

    // Hull variables
    label m = -1;
    scalar minQuality;

    // Dynamic programming variables
    scalarListList Q;
    labelListList K, triangulations;

    // Allocate dynamic programming tables
    mesh.initTables(Q, K, triangulations);

    // Pick edges off the stack
    while (!mesh.edgeStack(tIndex).empty())
    {
        // Retrieve an edge from the stack
        label eIndex = mesh.edgeStack(tIndex).pop();

        // Compute the minimum quality of cells around this edge
        mesh.computeMinQuality
        (
            eIndex,
            minQuality
        );

        // Check if this edge is on a bounding curve
        if (mesh.checkBoundingCurve(eIndex))
        {
            continue;
        }

        // Fill the dynamic programming tables
        if (mesh.fillTables(eIndex, m, Q, K, triangulations))
        {
            // Check if edge-swapping is required.
            if (mesh.checkQuality(eIndex, m, Q, minQuality))
            {
                if (thread->master())
                {
                    // Remove this edge according to the swap sequence
                    mesh.removeEdgeFlips(eIndex, K, triangulations);
                }
                else
                {
                    // Push this on to the master stack
                    mesh.edgeStack(0).push(eIndex);
                }
            }
        }
    }

    // If this is a slave, relinquish the lock for this structure
    if (thread->slave())
    {
        thread->unlock();
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

    // If this is a slave, acquire the lock for this structure
    if (thread->slave())
    {
        thread->lock();
    }

    dynamicTopoFvMesh& mesh = thread->mesh();

    // Figure out which thread this is...
    label tIndex = mesh.self();
    scalar length = 0.0, scale = 0.0;

    while (!mesh.edgeStack(tIndex).empty())
    {
        // Retrieve an edge from the stack
        label eIndex = mesh.edgeStack(tIndex).pop();

        // Measure the edge-length
        mesh.edgeLength(eIndex, length);

        // Determine the length-scale at this point in the mesh
        mesh.meshEdgeLengthScale(eIndex, scale);

        if (length > mesh.ratioMax()*scale)
        {
            if (thread->master())
            {
                // Bisect this edge
                mesh.bisectEdge(eIndex);
            }
            else
            {
                // Push this on to the master stack
                mesh.edgeStack(0).push(eIndex);
            }
        }
        else
        if (length < mesh.ratioMin()*scale)
        {
            if (thread->master())
            {
                // Collapse this edge
                mesh.collapseEdge(eIndex);
            }
            else
            {
                // Push this on to the master stack
                mesh.edgeStack(0).push(eIndex);
            }
        }
    }

    // If this is a slave, relinquish the lock for this structure
    if (thread->slave())
    {
        thread->unlock();
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
            points_[faceToCheck[2]],
            points_[faceToCheck[1]],
            points_[faceToCheck[0]],
            newPoint
        );
    }
    else
    {
        cellVolume = tetVolume
        (
            points_[faceToCheck[0]],
            points_[faceToCheck[1]],
            points_[faceToCheck[2]],
            newPoint
        );
    }

    // Final cell-volume check
    if (cellVolume < VSMALL)
    {
        if (debug > 1)
        {
            InfoIn
            (
                "dynamicTopoFvMesh::checkCollapse"
            )   << "\nCollapsing cell: " << cellIndex
                << " containing points:\n"
                << faceToCheck[0] << "," << faceToCheck[1] << ","
                << faceToCheck[2] << "," << pointIndex << nl
                << "will yield a negative volume: " << cellVolume
                << ", when " << pointIndex
                << " is moved to location: " << nl
                << newPoint << endl;
        }

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
        if (self() == 0)
        {
            // Step 1: Bisect the boundary quad face
            bisectInterior_ = -1;
            bisectQuadFace(fIndex);

            // Step 2: Collapse the newly created internal quad face
            collapseQuadFace(bisectInterior_);
        }
        else
        {
            // Push this on to the master stack
            faceStack(0).push(fIndex);
        }
    }
}

// Remove sliver cells in 3D
void dynamicTopoFvMesh::removeSlivers()
{
    forAllIter(labelHashSet, thresholdSlivers_, iter)
    {
        label nIntFaces = 0, nBdyFaces = 0;
        label commonIntEdge = -1, commonBdyEdge = -1;
        FixedList<label, 4> intFaces(-1);
        FixedList<label, 4> bdyFaces(-1);

        cell& cellToCheck = cells_[iter.key()];

        // Determine the number of interior/boundary faces
        forAll(cellToCheck, faceI)
        {
            if (neighbour_[cellToCheck[faceI]] == -1)
            {
                bdyFaces[nBdyFaces++] = cellToCheck[faceI];
            }
            else
            {
                intFaces[nIntFaces++] = cellToCheck[faceI];
            }
        }

        // Check if this is a surface sliver
        if (nIntFaces == 2 && nBdyFaces == 2)
        {
            findCommonEdge
            (
                bdyFaces[0],
                bdyFaces[1],
                commonBdyEdge
            );

            findCommonEdge
            (
                intFaces[0],
                intFaces[1],
                commonIntEdge
            );

            WarningIn
            (
                "dynamicTopoFvMesh::removeSlivers()"
            )   << nl << "Removing Cell: " << iter.key()
                << nl << "Boundary Edge: " << commonBdyEdge
                << ": " << edges_[commonBdyEdge] << endl;

            // Step 1: Bisect the interior edge
            bisectEdge(commonIntEdge);

            // Reset the interior edge index
            bisectInterior_ = -1;

            // Step 2: Bisect the boundary edge
            bisectEdge(commonBdyEdge);

            // Step 3: Collapse the temporary interior edge
            collapseEdge(bisectInterior_);
        }
    }

    // Clear out the list
    thresholdSlivers_.clear();
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

    // Map all spherical tensor volFields in the objectRegistry
    MapGeometricFields
    <
        sphericalTensor,
        fvPatchField,
        dynamicTopoFvMeshMapper,
        volMesh
    >
    (
        fieldMapper
    );

    // Map all symmTensor volFields in the objectRegistry
    MapGeometricFields
    <
        symmTensor,
        fvPatchField,
        dynamicTopoFvMeshMapper,
        volMesh
    >
    (
        fieldMapper
    );

    // Map all tensor volFields in the objectRegistry
    MapGeometricFields
    <
        tensor,
        fvPatchField,
        dynamicTopoFvMeshMapper,
        volMesh
    >
    (
        fieldMapper
    );

    // Map all the scalar surfaceFields in the objectRegistry
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

    // Map all the vector surfaceFields in the objectRegistry
    MapGeometricFields
    <
        vector,
        fvsPatchField,
        dynamicTopoFvMeshMapper,
        surfaceMesh
    >
    (
        fieldMapper
    );

    // Map all the sphericalTensor surfaceFields in the objectRegistry
    MapGeometricFields
    <
        sphericalTensor,
        fvsPatchField,
        dynamicTopoFvMeshMapper,
        surfaceMesh
    >
    (
        fieldMapper
    );

    // Map all the symmTensor surfaceFields in the objectRegistry
    MapGeometricFields
    <
        symmTensor,
        fvsPatchField,
        dynamicTopoFvMeshMapper,
        surfaceMesh
    >
    (
        fieldMapper
    );

    // Map all the tensor surfaceFields in the objectRegistry
    MapGeometricFields
    <
        tensor,
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
            // Lock slaves
            lockSlaveThreads();

            // Submit jobs to the work queue
            for (label i = 1; i <= threader_->getNumThreads(); i++)
            {
                threader_->addToWorkQueue
                           (
                               &edgeBisectCollapse2D,
                               reinterpret_cast<void *>(&(structPtr_[i]))
                           );
            }

            // Unlock slaves
            unlockSlaveThreads();

            // Synchronize threads
            synchronizeThreads();
        }

        // Set the master thread to implement modifications
        edgeBisectCollapse2D(&(structPtr_[0]));

        if (debug)
        {
            Info << nl << "2D Edge Bisection/Collapse complete." << endl;
        }
    }

    // Re-initialize the face stacks
    initFaceStacks();

    if (threader_->multiThreaded())
    {
        // Lock slaves
        lockSlaveThreads();

        // Submit jobs to the work queue
        for (label i = 1; i <= threader_->getNumThreads(); i++)
        {
            threader_->addToWorkQueue
                       (
                           &swap2DEdges,
                           reinterpret_cast<void *>(&(structPtr_[i]))
                       );
        }

        // Unlock slaves
        unlockSlaveThreads();

        // Synchronize threads
        synchronizeThreads();
    }

    // Set the master thread to implement modifications
    swap2DEdges(reinterpret_cast<void *>(&(structPtr_[0])));

    if (debug)
    {
        Info << nl << "2D Edge Swapping complete." << endl;
    }
}

// MultiThreaded topology modifier [3D]
void dynamicTopoFvMesh::threadedTopoModifier3D()
{
    // Send and receive patch sub-meshes for coupled patches
    sendAndRecvCoupledMeshes();

    // Handle coupled patches first
    handleCoupledPatches();

    if (edgeModification_)
    {
        // Initialize the edge stacks
        initEdgeStacks();

        if (threader_->multiThreaded())
        {
            // Lock slaves
            lockSlaveThreads();

            // Submit jobs to the work queue
            for (label i = 1; i <= threader_->getNumThreads(); i++)
            {
                threader_->addToWorkQueue
                           (
                               &edgeBisectCollapse3D,
                               reinterpret_cast<void *>(&(structPtr_[i]))
                           );
            }

            // Unlock slaves
            unlockSlaveThreads();

            // Synchronize threads
            synchronizeThreads();
        }

        // Set the master thread to implement modifications
        edgeBisectCollapse3D(&(structPtr_[0]));

        if (debug)
        {
            Info << nl << "3D Edge Bisection/Collapse complete." << endl;
        }
    }

    // Re-initialize the edge stacks
    initEdgeStacks();

    if (threader_->multiThreaded())
    {
        // Lock slaves
        lockSlaveThreads();

        // Submit jobs to the work queue
        for (label i = 1; i <= threader_->getNumThreads(); i++)
        {
            threader_->addToWorkQueue
                       (
                           &swap3DEdges,
                           reinterpret_cast<void *>(&(structPtr_[i]))
                       );
        }

        // Unlock slaves
        unlockSlaveThreads();

        // Synchronize threads
        synchronizeThreads();
    }

    // Set the master thread to implement modifications
    swap3DEdges(reinterpret_cast<void *>(&(structPtr_[0])));

    if (debug)
    {
        Info << nl << "3D Edge Swapping complete." << endl;
    }
}

// Lock all slave threads
void dynamicTopoFvMesh::lockSlaveThreads()
{
    for (label i = 1; i <= threader_->getNumThreads(); i++)
    {
        structPtr_[i].lock();
    }
}

// Unlock all slave threads
void dynamicTopoFvMesh::unlockSlaveThreads()
{
    for (label i = 1; i <= threader_->getNumThreads(); i++)
    {
        structPtr_[i].unlock();
    }

    // Allow slaves to acquire the mutex.
    // Should be a better mechanism to do this.
    sleep(1);
}

// Synchronize all slave threads
void dynamicTopoFvMesh::synchronizeThreads()
{
    for (label i = 1; i <= threader_->getNumThreads(); i++)
    {
        structPtr_[i].lock();

        // Now that the lock is acquired, unlock it for later use
        structPtr_[i].unlock();
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
    // Re-read options if they have been modified at run-time
    if (dict_.readIfModified())
    {
        // Enable/disable run-time debug level
        if (dict_.found("debug"))
        {
            debug = readLabel(dict_.lookup("debug"));
        }

        if (dict_.subDict("dynamicTopoFvMesh").found("interval"))
        {
            interval_ = readLabel
                        (
                            dict_.subDict
                            ("dynamicTopoFvMesh").lookup("interval")
                        );
        }

        // Read edge options
        readEdgeOptions();
    }

    // Return if re-meshing is not at interval
    if (time().timeIndex() % interval_ != 0)
    {
        return false;
    }

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
    if (debug > 1)
    {
        label band=0;

        const labelList& oldOwner = faceOwner();
        const labelList& oldNeighbour = faceNeighbour();

        for(label faceI = 0; faceI < nInternalFaces_; faceI++)
        {
            label diff = oldNeighbour[faceI] - oldOwner[faceI];

            if (diff > band)
            {
                band = diff;
            }
        }

        reduce(band, maxOp<label>());

        Info << "Mesh size: " << nCells()
             << "    Bandwidth before renumbering: " << band << endl;
    }

    // Obtain the most recent point-positions
    const pointField& currentPoints = points();

    forAllIter(HashList<point>::iterator, points_, pIter)
    {
        pIter() = currentPoints[pIter.index()];
    }

    // Track mesh topology modification time
    clockTime topologyTimer;

    //== Connectivity changes ==//

    // Reset the flag
    thresholdSlivers_.clear();
    topoChangeFlag_ = false;
    nModifications_ = 0;
    nBisections_ = 0;
    nCollapses_ = 0;
    nSwaps_ = 0;

    // Obtain mesh stats before topo-changes
    meshQuality(true);

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
        if (debug > 1)
        {
            label band=0;

            for(label faceI = 0; faceI < nInternalFaces_; faceI++)
            {
                label diff = neighbour[faceI] - owner[faceI];

                if (diff > band)
                {
                    band = diff;
                }
            }

            reduce(band, maxOp<label>());

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
        faceParents_.clear();

        // Set new sizes for the reverse maps
        reversePointMap_.setSize(nPoints_);
        reverseEdgeMap_.setSize(nEdges_);
        reverseFaceMap_.setSize(nFaces_);
        reverseCellMap_.setSize(nCells_);

        // Basic checks for mesh-validity
        if (debug)
        {
            checkMesh(true);
        }
    }

    // Print out topo-stats
    Info << "Reordering time: " << reOrderingTimer.elapsedTime() << endl;
    Info << "nBisections: " << nBisections_ << endl;
    Info << "nCollapses: " << nCollapses_ << endl;
    Info << "nSwaps: " << nSwaps_ << endl;

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
