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

#include <iomanip>
#include "IOmanip.H"
#include <dlfcn.h>

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynamicTopoFvMesh,0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from IOobject
dynamicTopoFvMesh::dynamicTopoFvMesh(const IOobject& io)
:
    fvMesh(io),
    numPatches_(polyMesh::boundaryMesh().size()),
    topoChangeFlag_(false),
    isSubMesh_(false),
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
    mandatory_
    (
        dict_.subDict("dynamicTopoFvMesh").lookup("allOptionsMandatory")
    ),
    twoDMesh_(polyMesh::nGeometricD() == 2 ? true : false),
    edgeRefinement_
    (
        dict_.subDict("dynamicTopoFvMesh").lookup("edgeRefinement")
    ),
    coupledModification_(false),
    slaveModification_(false),
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
    curvatureDeviation_(0.0),
    minLengthScale_(VSMALL),
    maxLengthScale_(GREAT),
    sliverThreshold_(0.1),
    nModifications_(0),
    nBisections_(0),
    nCollapses_(0),
    nSwaps_(0),
    maxModifications_(-1),
    bisectInterior_(-1),
    proximityBins_(0),
    sliceThreshold_(VSMALL),
    sliceHoldOff_(0),
    sliceBoxes_(0),
    slicePairs_(0),
    maxTetsPerEdge_(-1),
    allowTableResize_(false),
    gTol_(1e-20)
{
    // For backward compatibility, check the size of owner/neighbour
    if (owner_.size() != neighbour_.size())
    {
        // Size up to number of faces
        neighbour_.setSize(nFaces_);

        // Padding with -1 for neighbours
        for(label i = nInternalFaces_; i < nFaces_; i++)
        {
            neighbour_[i] = -1;
        }
    }

    // Initialize the multiThreading environment
    initializeThreadingEnvironment();

    // Read optional parameters.
    readOptionalParameters();

    // Initialize patch-size information
    const polyBoundaryMesh& boundary = polyMesh::boundaryMesh();
    for(label i=0; i<numPatches_; i++)
    {
        patchNMeshPoints_[i] = boundary[i].meshPoints().size();
        oldPatchSizes_[i]  = patchSizes_[i]  = boundary[i].size();
        oldPatchStarts_[i] = patchStarts_[i] = boundary[i].start();
        oldPatchNMeshPoints_[i] = patchNMeshPoints_[i];
    }

    // Open the tetMetric dynamic-link library (for 3D only)
    loadMetricLibrary();

    // Initialize edge-related connectivity structures
    initEdges();

    // Set sizes for the reverse maps
    reversePointMap_.setSize(nPoints_, -7);
    reverseEdgeMap_.setSize(nEdges_, -7);
    reverseFaceMap_.setSize(nFaces_, -7);
    reverseCellMap_.setSize(nCells_, -7);

    // Define edgeRefinement options
    readRefinementOptions();
}

//- Construct from components. Used for subMeshes only.
dynamicTopoFvMesh::dynamicTopoFvMesh
(
    const dynamicTopoFvMesh& mesh,
    const IOobject& io,
    const pointField& points,
    const label nInternalEdges,
    const edgeList& edges,
    const faceList& faces,
    const labelListList& faceEdges,
    const cellList& cells
)
:
    fvMesh(io, points, faces, cells, false),
    numPatches_(1),
    topoChangeFlag_(false),
    isSubMesh_(true),
    dict_(mesh.dynamicMeshDict()),
    mandatory_(mesh.mandatory_),
    twoDMesh_(mesh.twoDMesh_),
    edgeRefinement_(mesh.edgeRefinement_),
    coupledModification_(false),
    slaveModification_(false),
    interval_(1),
    mapper_(NULL),
    points_(points),
    faces_(faces),
    cells_(cells),
    edges_(edges),
    faceEdges_(faceEdges),
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
    nOldPoints_(points.size()),
    nPoints_(points.size()),
    nOldEdges_(edges.size()),
    nEdges_(edges.size()),
    nOldFaces_(faces.size()),
    nFaces_(faces.size()),
    nOldCells_(cells.size()),
    nCells_(cells.size()),
    nOldInternalFaces_(0),
    nInternalFaces_(0),
    nOldInternalEdges_(nInternalEdges),
    nInternalEdges_(nInternalEdges),
    ratioMin_(mesh.ratioMin_),
    ratioMax_(mesh.ratioMax_),
    growthFactor_(mesh.growthFactor_),
    curvatureDeviation_(mesh.curvatureDeviation_),
    minLengthScale_(mesh.minLengthScale_),
    maxLengthScale_(mesh.maxLengthScale_),
    sliverThreshold_(mesh.sliverThreshold_),
    nModifications_(0),
    nBisections_(0),
    nCollapses_(0),
    nSwaps_(0),
    maxModifications_(mesh.maxModifications_),
    bisectInterior_(-1),
    proximityBins_(0),
    sliceThreshold_(VSMALL),
    sliceHoldOff_(0),
    sliceBoxes_(0),
    slicePairs_(0),
    maxTetsPerEdge_(mesh.maxTetsPerEdge_),
    allowTableResize_(mesh.allowTableResize_),
    gTol_(mesh.gTol_),
    tetMetric_(mesh.tetMetric_)
{
    // Initialize owner and neighbour
    owner_.setSize(faces.size(), -1);
    neighbour_.setSize(faces.size(), -1);

    boolList markedFaces(nFaces_, false);

    forAll(cells, cellI)
    {
        const labelList& thisCell = cells[cellI];

        forAll(thisCell, faceI)
        {
            if (!markedFaces[thisCell[faceI]])
            {
                // First visit: Owner
                owner_[thisCell[faceI]] = cellI;

                markedFaces[thisCell[faceI]] = true;
            }
            else
            {
                // Second visit: Neighbour
                neighbour_[thisCell[faceI]] = cellI;

                nInternalFaces_++;
            }
        }
    }

    nOldInternalFaces_ = nInternalFaces_;

    // Initialize the multiThreading environment.
    // Force to single-threaded.
    initializeThreadingEnvironment(1);

    const polyBoundaryMesh& boundary = polyMesh::boundaryMesh();

    // Add a default patch as boundary for polyMesh.
    polyMesh::addPatches
    (
        List<polyPatch*>
        (
            1,
            new polyPatch
            (
                "defaultPatch",
                (nFaces_ - nInternalFaces_),
                nInternalFaces_,
                0,
                boundary
            )
        )
    );

    // Initialize patch-size information
    oldPatchSizes_[0] = patchSizes_[0] = (nFaces_ - nInternalFaces_);
    oldPatchStarts_[0] = patchStarts_[0] = nInternalFaces_;
    oldEdgePatchSizes_[0] = edgePatchSizes_[0] = (nEdges_ - nInternalEdges_);
    oldEdgePatchStarts_[0] = edgePatchStarts_[0] = nInternalEdges_;

    // Set sizes for the reverse maps
    reversePointMap_.setSize(nPoints_, -7);
    reverseEdgeMap_.setSize(nEdges_, -7);
    reverseFaceMap_.setSize(nFaces_, -7);
    reverseCellMap_.setSize(nCells_, -7);

    // Now build edgeFaces and pointEdges information.
    invertConnectivity(nEdges_, faceEdges_, edgeFaces_);
    invertConnectivity(nPoints_, edges_, pointEdges_);

    // Size-up edgePoints for now, but explicitly construct
    // for each edge later, based on point coupling.
    if (!twoDMesh_)
    {
        edgePoints_.setSize(nEdges_, labelList(0));
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dynamicTopoFvMesh::~dynamicTopoFvMesh()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return the mesh-mapper
const mapPolyMesh& dynamicTopoFvMesh::meshMap() const
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

    if (edgeRefinement_)
    {
        scalarField& internalField = tlengthScale();

        // Re-calculate lengthScale
        calculateLengthScale();

        // Obtain length-scale values from the mesh
        forAll(lengthScale_, cellI)
        {
            internalField[cellI] = lengthScale_[cellI];
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
        new scalarField(cells_.size(), 0.0)
    );

    // Valid for 3D tetrahedral meshes only...
    if (!twoDMesh_)
    {
        scalarField& iF = tQuality();

        // Compute statistics on the fly
        label nCells = 0;
        scalar maxQuality = -GREAT;
        scalar minQuality =  GREAT;
        scalar meanQuality = 0.0;

        // Loop through all cells in the mesh and compute cell quality
        forAll(cells_, cellI)
        {
            if (cells_[cellI].empty())
            {
                continue;
            }

            // Compute cell quality
            iF[cellI] = tetQuality(cellI);

            // Update statistics
            maxQuality = Foam::max(iF[cellI], maxQuality);
            minQuality = Foam::min(iF[cellI], minQuality);
            meanQuality += iF[cellI];
            nCells++;

            // Add to the list of slivers
            if
            (
                (iF[cellI] < sliverThreshold_)
             && (iF[cellI] > 0.0)
            )
            {
                thresholdSlivers_.insert(cellI, iF[cellI]);
            }
        }

        // Output statistics:
        if (outputOption || (debug > 0))
        {
            // Reduce statistics across processors.
            reduce(minQuality, minOp<scalar>());
            reduce(maxQuality, maxOp<scalar>());
            reduce(meanQuality, sumOp<scalar>());
            reduce(nCells, sumOp<label>());

            if (minQuality < 0.0)
            {
                WarningIn
                (
                    "dynamicTopoFvMesh::meshQuality()"
                )   << nl << "Minimum cell quality is: "
                    << minQuality << endl;
            }

            Info << " ~~~ Mesh Quality Statistics ~~~ " << endl;
            Info << " Min: " << minQuality << endl;
            Info << " Max: " << maxQuality << endl;
            Info << " Mean: " << meanQuality/nCells << endl;
            Info << " Cells: " << nCells << endl;
            Info << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << endl;
        }
    }

    return tQuality;
}

// Perform a Delaunay test on an internal face
bool dynamicTopoFvMesh::testDelaunay
(
    const label fIndex
) const
{
    bool failed = false, procCouple = false;
    label eIndex = -1, pIndex = -1, fLabel = -1;
    FixedList<bool,2> foundTriFace(false);
    FixedList<FixedList<label,3>,2> triFaces(FixedList<label,3>(-1));

    // Boundary faces are discarded.
    if (whichPatch(fIndex) > -1)
    {
        procCouple = processorCoupledFace(fIndex);

        if (!procCouple)
        {
            return failed;
        }
    }

    if (procCouple)
    {
        // Detect faces across processor boundaries.

    }
    else
    {
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
                const face& thisFace = faces_[eFaces[faceI]];

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
    }

    // Obtain point references for the first face
    point a = points_[triFaces[0][0]];

    point cCenter = circumCenter(fLabel), otherPoint = vector::zero;
    scalar rSquared = (a - cCenter)&(a - cCenter);

    // Find the isolated point on the second face
    if (procCouple)
    {
        // Find the other point across the processor boundary.

    }
    else
    {
        const edge& e = edges_[eIndex];

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
        otherPoint = points_[pIndex];
    }

    if
    (
        ((otherPoint - cCenter)&(otherPoint - cCenter)) < rSquared
    )
    {
        // Failed the test.
        failed = true;
    }

    return failed;
}

// Utility method to find the interior (quad) / boundary (tri) faces
// for an input quad-face and adjacent triangle-prism cell.
void dynamicTopoFvMesh::findPrismFaces
(
    const label fIndex,
    const label cIndex,
    FixedList<face,2>& bdyf,
    FixedList<label,2>& bidx,
    FixedList<face,2>& intf,
    FixedList<label,2>& iidx
) const
{
    label indexO = 0, indexI = 0;

    const cell& c = cells_[cIndex];

    forAll(c, i)
    {
        label faceIndex = c[i];

        // Don't count the face under consideration
        if (faceIndex != fIndex)
        {
            const face& fi = faces_[faceIndex];

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
) const
{
    bool found = false;

    const labelList& fEi = faceEdges_[first];
    const labelList& fEj = faceEdges_[second];

    forAll(fEi, edgeI)
    {
        forAll(fEj, edgeJ)
        {
            if (fEi[edgeI] == fEj[edgeJ])
            {
                common = fEi[edgeI];

                found = true;
                break;
            }
        }

        if (found)
        {
            break;
        }
    }

    return found;
}

// Insert the specified cell to the mesh
label dynamicTopoFvMesh::insertCell
(
    const cell& newCell,
    const scalar lengthScale,
    const label mappingCell,
    const label zoneID
)
{
    label newCellIndex = cells_.size();

    if (debug > 2)
    {
        Info << "Inserting cell: "
             << newCellIndex << ": "
             << newCell << endl;
    }

    cells_.append(newCell);

    if (edgeRefinement_)
    {
        lengthScale_.append(lengthScale);
    }

    // Generate mapping information for this new cell
    const labelListList& cc = cellCells();

    label parent;
    labelHashSet masterObjects;

    if (mappingCell < nOldCells_)
    {
        parent = mappingCell;
    }
    else
    {
        parent = cellParents_[mappingCell];
    }

    // Insert the parent cell
    cellParents_.insert(newCellIndex, parent);

    // Find the cell's neighbours in the old mesh
    masterObjects.insert(parent);

    forAll(cc[parent], cellI)
    {
        if (!masterObjects.found(cc[parent][cellI]))
        {
            masterObjects.insert(cc[parent][cellI]);
        }
    }

    // Insert mapping info into the HashTable
    cellsFromCells_.insert
    (
        newCellIndex,
        objectMap
        (
            newCellIndex,
            masterObjects.toc()
        )
    );

    // Add to the zone if necessary
    if (zoneID >= 0)
    {
        addedCellZones_.insert(newCellIndex, zoneID);
    }

    nCells_++;

    return newCellIndex;
}

// Remove the specified cell from the mesh
void dynamicTopoFvMesh::removeCell
(
    const label cIndex
)
{
    if (debug > 2)
    {
        Info << "Removing cell: "
             << cIndex << ": "
             << cells_[cIndex]
             << endl;
    }

    cells_[cIndex].clear();

    if (edgeRefinement_)
    {
        lengthScale_[cIndex] = -1.0;
    }

    // Update the number of cells, and the reverse cell map
    nCells_--;

    if (cIndex < nOldCells_)
    {
        reverseCellMap_[cIndex] = -1;
    }
    else
    {
        // Store this information for the reOrdering stage
        deletedCells_.insert(cIndex);
    }

    // Check if the cell was added in the current morph, and delete
    if (cellsFromCells_.found(cIndex))
    {
        cellsFromCells_.erase(cIndex);
    }

    // Check if this cell was added to a zone
    if (addedCellZones_.found(cIndex))
    {
        addedCellZones_.erase(cIndex);
    }
}

// Utility method for face-insertion
label dynamicTopoFvMesh::insertFace
(
    const label patch,
    const face& newFace,
    const label newOwner,
    const label newNeighbour,
    const label zoneID
)
{
    // Append the specified face to each face-related list.
    // Reordering is performed after all pending changes
    // (flips, bisections, contractions, etc) have been made to the mesh
    label newFaceIndex = faces_.size();

    faces_.append(newFace);
    owner_.append(newOwner);
    neighbour_.append(newNeighbour);

    if (debug > 2)
    {
        Info << "Inserting face: "
             << newFaceIndex << ": "
             << newFace;

        Info << " Patch: ";

        if (patch == -1)
        {
            Info << "Internal" << endl;
        }
        else
        {
            Info << boundaryMesh()[patch].name() << endl;
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

    // Add to the zone if explicitly specified
    if (zoneID >= 0)
    {
        addedFaceZones_.insert(newFaceIndex, zoneID);
    }
    else
    {
        // No zone was specified. Check if this
        // face is added to a coupled patch associated
        // with a faceZone.
        forAllIter(Map<coupledPatchInfo>, patchCoupling_, pIter)
        {
            if (pIter.key() == patch)
            {
                if (pIter().masterFaceZone() > -1)
                {
                    addedFaceZones_.insert
                    (
                        newFaceIndex,
                        pIter().masterFaceZone()
                    );
                }

                break;
            }

            if (pIter().slaveIndex() == patch)
            {
                if (pIter().slaveFaceZone() > -1)
                {
                    addedFaceZones_.insert
                    (
                        newFaceIndex,
                        pIter().slaveFaceZone()
                    );
                }

                break;
            }
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
    const cell& cellToCheck = cells_[index];

    forAll(cellToCheck, faceI)
    {
        const labelList& faceEdges = faceEdges_[cellToCheck[faceI]];

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

                const edge& edgeToCheck = edges_[faceEdges[edgeI]];

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

                if (edgeFaces_[faceEdges[edgeI]].empty())
                {
                    // Hanging edge. Check its points and remove
                    // them if necessary.
                    forAll(edgeToCheck, pointI)
                    {
                        // Check for hanging nodes...
                        if (pointEdges_[edgeToCheck[pointI]].size() == 1)
                        {
                            removePoint(edgeToCheck[pointI]);
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
            (
                insertFace
                (
                    patch,
                    newFace,
                    newOwner,
                    -1
                )
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
    removeCell(index);

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

    // Clear entities.
    faces_[index].clear();
    owner_[index] = -1;
    neighbour_[index] = -1;
    faceEdges_[index].clear();

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
    else
    {
        // Store this information for the reOrdering stage
        deletedFaces_.insert(index);
    }

    // Check and remove from the list of added face patches
    if (addedFacePatches_.found(index))
    {
        addedFacePatches_.erase(index);
    }

    // Check if this face was added to a zone
    if (addedFaceZones_.found(index))
    {
        addedFaceZones_.erase(index);
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
    label newEdgeIndex = edges_.size();

    edges_.append(newEdge);
    edgeFaces_.append(edgeFaces);

    if (!twoDMesh_)
    {
        edgePoints_.append(edgePoints);
    }

    if (debug > 2)
    {
        Info << "Inserting edge: "
             << newEdgeIndex << ": "
             << newEdge;

        Info << " Patch: ";

        if (patch == -1)
        {
            Info << "Internal" << endl;
        }
        else
        {
            Info << boundaryMesh()[patch].name() << endl;
        }

        if (findIndex(edgePoints, -1) != -1)
        {
            FatalErrorIn("dynamicTopoFvMesh::insertEdge()")
                << " EdgePoints is incorrectly specified." << nl
                << " edgePoints: " << edgePoints << nl
                << abort(FatalError);
        }
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
    if (!twoDMesh_)
    {
        sizeUpList(newEdgeIndex, pointEdges_[newEdge[0]]);
        sizeUpList(newEdgeIndex, pointEdges_[newEdge[1]]);
    }

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
    if (debug > 2)
    {
        Info << "Removing edge: "
             << index << ": "
             << edges_[index] << endl;
    }

    if (!twoDMesh_)
    {
        // Remove the edgePoints entry
        edgePoints_[index].clear();

        // Size-down the pointEdges list
        if (pointEdges_[edges_[index][0]].size())
        {
            sizeDownList(index, pointEdges_[edges_[index][0]]);
        }

        if (pointEdges_[edges_[index][1]].size())
        {
            sizeDownList(index, pointEdges_[edges_[index][1]]);
        }

        // Remove from the stack as well
        forAll(edgeStack_, stackI)
        {
            edgeStack(stackI).remove(index);
        }
    }

    edges_[index] = edge(-1, -1);
    edgeFaces_[index].clear();

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
    else
    {
        // Store this information for the reOrdering stage
        deletedEdges_.insert(index);
    }

    // Check and remove from the list of added edge patches
    if (addedEdgePatches_.found(index))
    {
        addedEdgePatches_.erase(index);
    }

    // Decrement the total edge-count
    nEdges_--;
}

// Insert the specified point to the mesh
label dynamicTopoFvMesh::insertPoint
(
    const point& newPoint,
    const labelList& mappingPoints,
    const label zoneID
)
{
    // Add a new point to the end of the list
    label newPointIndex = points_.size();

    points_.append(newPoint);

    if (debug > 2)
    {
        Info << "Inserting point: "
             << newPointIndex << ": "
             << newPoint
             << "  Mapped from: "
             << mappingPoints << endl;
    }

    // Add an empty entry to pointEdges as well.
    // This entry can be sized-up appropriately at a later stage.
    if (!twoDMesh_)
    {
        pointEdges_.append(labelList(0));
    }

    labelHashSet masterObjects;

    forAll(mappingPoints, pointI)
    {
        label parent;

        if (mappingPoints[pointI] < nOldPoints_)
        {
            parent = mappingPoints[pointI];
        }
        else
        {
            parent = pointParents_[mappingPoints[pointI]];
        }

        if (!masterObjects.found(parent))
        {
            masterObjects.insert(parent);
        }
    }

    if (twoDMesh_)
    {
        pointParents_.insert(newPointIndex, masterObjects.begin().key());
    }
    else
    {
        // Insert the parent point.
        // This has to be done carefully: If a mapping point is on
        // the surface, then the new point must also preferentially
        // map from the surface point.
        label surfPointIndex = -1;

        forAllIter(labelHashSet, masterObjects, oIter)
        {
            const labelList& pEdges = pointEdges_[oIter.key()];

            bool foundBoundaryEdge = false;

            forAll(pEdges, edgeI)
            {
                if (whichEdgePatch(pEdges[edgeI]) > -1)
                {
                    surfPointIndex =
                    (
                        edges_[pEdges[edgeI]].otherVertex(oIter.key())
                    );

                    foundBoundaryEdge = true;
                    break;
                }
            }

            if (foundBoundaryEdge)
            {
                break;
            }
        }

        if (surfPointIndex == -1)
        {
            pointParents_.insert
            (
                newPointIndex, masterObjects.begin().key()
            );
        }
        else
        {
            pointParents_.insert(newPointIndex, surfPointIndex);
        }
    }

    // Insert mapping info into the HashTable
    pointsFromPoints_.insert
    (
        newPointIndex,
        objectMap
        (
            newPointIndex,
            masterObjects.toc()
        )
    );

    // Add to the zone if necessary
    if (zoneID >= 0)
    {
        addedPointZones_.insert(newPointIndex, zoneID);
    }

    nPoints_++;

    return newPointIndex;
}

// Remove the specified point from the mesh
void dynamicTopoFvMesh::removePoint
(
    const label index
)
{
    if (debug > 2)
    {
        Info << "Removing point: "
             << index << ": "
             << points_[index] << endl;
    }

    // Remove the point
    points_[index] = point();

    // Remove pointEdges as well
    if (!twoDMesh_)
    {
        pointEdges_[index].clear();
    }

    // Update the reverse point map
    if (index < nOldPoints_)
    {
        reversePointMap_[index] = -1;
    }
    else
    {
        deletedPoints_.insert(index);
    }

    // Check if the point was added in the current morph, and delete
    if (pointsFromPoints_.found(index))
    {
        pointsFromPoints_.erase(index);
    }

    // Check if this point was added to a zone
    if (addedPointZones_.found(index))
    {
        addedPointZones_.erase(index);
    }

    // Decrement the total point-count
    nPoints_--;
}

// Utility method to build a hull of cells connected to the edge [2D]
void dynamicTopoFvMesh::constructPrismHull
(
    const label eIndex,
    labelHashSet& hullTriFaces,
    labelHashSet& hullCells
) const
{
    // Obtain references
    const labelList& eFaces = edgeFaces_[eIndex];

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
            const cell& cellToCheck = cells_[c0];

            forAll(cellToCheck, faceJ)
            {
                const face& faceToCheck = faces_[cellToCheck[faceJ]];

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
            !hullCells.found(c1) &&
            (c1 != -1)
        )
        {
            // Add this cell
            hullCells.insert(c1);

            // Find associated triFaces and add them too
            const cell& cellToCheck = cells_[c1];

            forAll(cellToCheck, faceJ)
            {
                const face& faceToCheck = faces_[cellToCheck[faceJ]];

                if
                (
                    (faceToCheck.size() == 3) &&
                   !(hullTriFaces.found(cellToCheck[faceJ]))
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
) const
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
    const edge& edgeToCheck = edges_[eIndex];
    const labelList& eFaces = edgeFaces_[eIndex];
    const labelList& hullVertices = edgePoints_[eIndex];

    // Temporary tri-face for comparison
    face oFace(3);

    // Loop through all faces of this edge and add them to hullFaces
    forAll(eFaces, faceI)
    {
        const face& faceToCheck = faces_[eFaces[faceI]];

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
                const labelList& fEdges = faceEdges_[hullFaces[indexI]];

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
                    const cell& currCell = cells_[hullCells[indexI]];

                    // Look for the ring-faces
                    forAll(currCell, faceI)
                    {
                        const face& cFace = faces_[currCell[faceI]];

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
                    const labelList& rFaceEdges =
                    (
                        faceEdges_[ringEntities[1][indexI]]
                    );

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
    const label eIndex,
    const label checkIndex
)
{
    bool found = false;
    label faceIndex = -1, cellIndex = -1;
    label otherPoint = -1, nextPoint = -1;

    // Obtain references
    const edge& edgeToCheck = edges_[eIndex];
    const labelList& eFaces = edgeFaces_[eIndex];

    // Re-size the list first
    labelList& ePoints = edgePoints_[eIndex];
    ePoints.setSize(eFaces.size(), -1);

    if (whichEdgePatch(eIndex) == -1 && !isSubMesh_)
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

                if (nextPoint == edgeToCheck[checkIndex])
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
        if (nextPoint == edgeToCheck[checkIndex])
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
                << " Edge: " << eIndex << ":: " << edgeToCheck << nl
                << " edgeFaces: " << eFaces
                << abort(FatalError);
        }
    }
}

// Utility to invert a connectivity list (pointEdges from edges, etc)
// Re-write of invertManyToMany, to work on resizableLists
template<class InList, class OutList>
void dynamicTopoFvMesh::invertConnectivity
(
    const label nEntities,
    const resizableList<InList>& inEntities,
    resizableList<OutList>& outEntities
) const
{
    labelList nInPerOut(nEntities, 0);

    forAll(inEntities, indexI)
    {
        const InList& inEntity = inEntities[indexI];

        forAll(inEntity, indexJ)
        {
            nInPerOut[inEntity[indexJ]]++;
        }
    }

    // Size outEntities
    outEntities.setSize(nEntities);

    forAll(nInPerOut, indexI)
    {
        outEntities[indexI].setSize(nInPerOut[indexI]);
    }

    nInPerOut = 0;

    // Fill outEntities
    forAll(inEntities, indexI)
    {
        const InList& inEntity = inEntities[indexI];

        forAll(inEntity, indexJ)
        {
            label entityI = inEntity[indexJ];

            outEntities[entityI][nInPerOut[entityI]++] = indexI;
        }
    }
}

// Utility to check whether points of an edge lie on a boundary.
const FixedList<bool,2>
dynamicTopoFvMesh::checkEdgeBoundary
(
    const label eIndex
) const
{
    FixedList<bool,2> edgeBoundary(false);

    const edge& edgeToCheck = edges_[eIndex];

    // Loop through edges connected to both points,
    // and check if any of them lie on boundaries.
    // Used to ensure that collapses happen towards boundaries.
    forAll(edgeToCheck, pointI)
    {
        const labelList& pEdges = pointEdges_[edgeToCheck[pointI]];

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

    return edgeBoundary;
}

// Check whether the given edge is on a bounding curve
bool dynamicTopoFvMesh::checkBoundingCurve(const label eIndex) const
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
        if (findIndex(noSwapPatchIDs_, edgePatch) > -1)
        {
            return true;
        }
    }

    if (coupledModification_)
    {
        return false;
    }

    // Check if two boundary faces lie on different face-patches
    FixedList<vector, 2> fNorm;
    label fPatch, firstPatch = -1, secondPatch = -1, count = 0;
    const labelList& edgeFaces = edgeFaces_[eIndex];

    forAll(edgeFaces, faceI)
    {
        if ((fPatch = whichPatch(edgeFaces[faceI])) > -1)
        {
            // Obtain the normal.
            fNorm[count] = triFaceNormal(faces_[edgeFaces[faceI]]);

            // Normalize it.
            fNorm[count] /= mag(fNorm[count]);

            count++;

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

    scalar deviation = (fNorm[0] & fNorm[1]);

    // Check if the curvature is too high
    if (mag(deviation) < 0.85)
    {
        return true;
    }

    // Check if the edge borders two different patches
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
    labelList& m,
    PtrList<scalarListList>& Q,
    PtrList<labelListList>& K,
    PtrList<labelListList>& triangulations,
    const label checkIndex
) const
{
    label mMax = maxTetsPerEdge_;

    // Check if resizing is necessary only for a particular index.
    if (checkIndex != -1)
    {
        m[checkIndex] = -1;
        Q[checkIndex].setSize((mMax-2),scalarList(mMax,-1.0));
        K[checkIndex].setSize((mMax-2),labelList(mMax,-1));
        triangulations[checkIndex].setSize(3,labelList((mMax-2),-1));

        return;
    }

    // Size all elements by default.
    label numIndices = -1;

    if (coupledModification_)
    {
        numIndices = getMaxCouplingIndex() + 1;
    }
    else
    {
        numIndices = 1;
    }

    m.setSize(numIndices, -1);
    Q.setSize(numIndices);
    K.setSize(numIndices);
    triangulations.setSize(numIndices);

    forAll(Q, indexI)
    {
        Q.set
        (
            indexI,
            new scalarListList((mMax-2),scalarList(mMax,-1.0))
        );

        K.set
        (
            indexI,
            new labelListList((mMax-2),labelList(mMax,-1))
        );

        triangulations.set
        (
            indexI,
            new labelListList(3,labelList((mMax-2),-1))
        );
    }
}

// Check triangulation quality for an edge index
bool dynamicTopoFvMesh::checkQuality
(
    const label eIndex,
    const labelList& m,
    const PtrList<scalarListList>& Q,
    const scalar minQuality,
    const label checkIndex
) const
{
    bool myResult = false;

    // Non-coupled check
    if (Q[checkIndex][0][m[checkIndex]-1] > minQuality)
    {
        myResult = true;

        if (debug > 2)
        {
            Info << " eIndex: " << eIndex
                 << " minQuality: " << minQuality
                 << " newQuality: " << Q[checkIndex][0][m[checkIndex]-1]
                 << endl;
        }
    }

    if (coupledModification_)
    {
        if (locallyCoupledEdge(eIndex))
        {
            // Check the quality of the slave edge as well.
            label slaveIndex = -1;

            // Loop through masterToSlave and determine the slave index.
            forAllConstIter(Map<coupledPatchInfo>, patchCoupling_, patchI)
            {
                if ((slaveIndex = patchI().findSlaveIndex(eIndex)) > -1)
                {
                    break;
                }
            }

            // Turn off switch temporarily.
            unsetCoupledModification();

            // Recursively call for the slave edge.
            myResult =
            (
                myResult && checkQuality(slaveIndex, m, Q, minQuality, 1)
            );

            // Turn it back on.
            setCoupledModification();
        }
        else
        if (processorCoupledEdge(eIndex))
        {

        }
    }

    return myResult;
}

// Utility method to fill the dynamic programming tables
//  - Returns true if the operation completed successfully.
//  - Returns false if tables could not be resized.
bool dynamicTopoFvMesh::fillTables
(
    const label eIndex,
    const scalar minQuality,
    labelList& m,
    PtrList<scalarListList>& Q,
    PtrList<labelListList>& K,
    PtrList<labelListList>& triangulations,
    const label checkIndex
) const
{
    const edge& edgeToCheck = edges_[eIndex];
    const labelList& hullVertices = edgePoints_[eIndex];

    // Fill in the size
    m[checkIndex] = hullVertices.size();

    // Check if a table-resize is necessary
    if (m[checkIndex] > maxTetsPerEdge_)
    {
        if (allowTableResize_)
        {
            // Resize the tables to account for
            // more tets per edge
            maxTetsPerEdge_ = m[checkIndex];

            // Clear tables for this index.
            Q[checkIndex].clear();
            K[checkIndex].clear();
            triangulations[checkIndex].clear();

            // Resize for this index.
            initTables(m, Q, K, triangulations, checkIndex);
        }
        else
        {
            // Can't resize. Bail out.
            return false;
        }
    }

    for (label i = (m[checkIndex]-3); i >= 0; i--)
    {
        for (label j = i+2; j < m[checkIndex]; j++)
        {
            for (label k = i+1; k < j; k++)
            {
                scalar q = (*tetMetric_)
                (
                    points_[hullVertices[i]],
                    points_[hullVertices[k]],
                    points_[hullVertices[j]],
                    points_[edgeToCheck[0]]
                );

                // For efficiency, check the bottom triangulation
                // only when the top one if less than the hull quality.
                if (q > minQuality)
                {
                    q =
                    (
                        Foam::min
                        (
                            q,
                            (*tetMetric_)
                            (
                                points_[hullVertices[j]],
                                points_[hullVertices[k]],
                                points_[hullVertices[i]],
                                points_[edgeToCheck[1]]
                            )
                        )
                    );
                }

                if (k < j-1)
                {
                    q = Foam::min(q,Q[checkIndex][k][j]);
                }

                if (k > i+1)
                {
                    q = Foam::min(q,Q[checkIndex][i][k]);
                }

                if ((k == i+1) || (q > Q[checkIndex][i][j]))
                {
                    Q[checkIndex][i][j] = q;
                    K[checkIndex][i][j] = k;
                }
            }
        }
    }

    if (coupledModification_)
    {
        if (locallyCoupledEdge(eIndex))
        {
            // Fill tables for the slave edge as well.
            label slaveIndex = -1;

            // Determine the slave index.
            forAllConstIter(Map<coupledPatchInfo>, patchCoupling_, patchI)
            {
                if ((slaveIndex = patchI().findSlaveIndex(eIndex)) > -1)
                {
                    break;
                }
            }

            // Turn off switch temporarily.
            unsetCoupledModification();

            // Recursively call for the slave edge.
            bool success =
            (
                fillTables(slaveIndex, minQuality, m, Q, K, triangulations, 1)
            );

            // Turn it back on.
            setCoupledModification();

            return success;
        }
        else
        if (processorCoupledEdge(eIndex))
        {

        }
    }

    // Print out tables for debugging
    if (debug > 3)
    {
        printTables(m, Q, K, checkIndex);
    }

    return true;
}

// Print out tables for debugging
void dynamicTopoFvMesh::printTables
(
    const labelList& m,
    const PtrList<scalarListList>& Q,
    const PtrList<labelListList>& K,
    const label checkIndex
) const
{
    Info << "m: " << m[checkIndex] << endl;

    // Print out Q
    Info << "===" << endl;
    Info << " Q " << endl;
    Info << "===" << endl;

    Info << "   ";

    for (label j = 0; j < m[checkIndex]; j++)
    {
        std::cout << std::setfill('-')
                  << std::setw(12) << j;
    }

    Info << nl;

    for (label i = 0; i < (m[checkIndex]-2); i++)
    {
        Info << i << ": ";

        for (label j = 0; j < m[checkIndex]; j++)
        {
            std::cout << std::setfill(' ')
                      << std::setw(12) << Q[checkIndex][i][j];
        }

        Info << nl;
    }

    // Print out K
    Info << "===" << endl;
    Info << " K " << endl;
    Info << "===" << endl;

    Info << "   ";

    for (label j = 0; j < m[checkIndex]; j++)
    {
        std::cout << std::setfill('-')
                  << std::setw(12) << j;
    }

    Info << nl;

    for (label i = 0; i < (m[checkIndex]-2); i++)
    {
        Info << i << ": ";

        for (label j = 0; j < m[checkIndex]; j++)
        {
            std::cout << std::setfill(' ')
                      << std::setw(12) << K[checkIndex][i][j];
        }

        Info << nl;
    }

    Info << endl;
}

// Remove the edge according to the swap sequence.
// Returns true if the swap-sequence was performed successfully.
bool dynamicTopoFvMesh::removeEdgeFlips
(
    const label eIndex,
    const scalar minQuality,
    const PtrList<labelListList>& K,
    PtrList<labelListList>& triangulations,
    const label checkIndex
)
{
    changeMap map;
    scalar swapQuality = GREAT;

    if (debug > 2)
    {
        Info << " Removing edge : " << eIndex << " by flipping."
             << " Edge: " << edges_[eIndex]
             << " minQuality: " << minQuality << endl;
    }

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
    extractTriangulation
    (
        0,
        m-1,
        K[checkIndex],
        numTriangulations,
        triangulations[checkIndex]
    );

    scalar tolF = 0.1;

    // Determine the final swap triangulation
    label tF =
    (
        identify32Swap
        (
            eIndex,
            hullVertices,
            triangulations[checkIndex],
            tolF
        )
    );

    // Check that the triangulation is valid
    label pIndex = -1;

    if (tF == -1)
    {
        // Reset all triangulations and bail out
        triangulations[checkIndex][0] = -1;
        triangulations[checkIndex][1] = -1;
        triangulations[checkIndex][2] = -1;

        return false;
    }
    else
    if ((pIndex = whichEdgePatch(eIndex)) > -1)
    {
        // Boundary edges may encounter precision issues
        // when trying to identify a 2-2 swap. Ensure that the right
        // decision was made.
        const labelListList& bT = triangulations[checkIndex];

        // Check if at least two points of the triangulation
        // are on the mesh boundary
        label nAttempts = 0, nBdyPoints = 0;
        bool foundPoints = false;

        while (!foundPoints)
        {
            forAll(bT, indexI)
            {
                const labelList& pEdges =
                (
                    pointEdges_[hullVertices[bT[indexI][tF]]]
                );

                forAll(pEdges, edgeI)
                {
                    if (whichEdgePatch(pEdges[edgeI]) == pIndex)
                    {
                        nBdyPoints++;
                        break;
                    }
                }
            }

            if (nBdyPoints >= 2)
            {
                foundPoints = true;
            }

            if (!foundPoints)
            {
                // Try again with a reduced tolerance.
                tolF *= 0.5;

                tF = identify32Swap(eIndex, hullVertices, bT, tolF);

                nAttempts++;
            }

            if (nAttempts > 2 || tF == -1)
            {
                // Reset all triangulations and bail out
                triangulations[checkIndex][0] = -1;
                triangulations[checkIndex][1] = -1;
                triangulations[checkIndex][2] = -1;

                return false;
            }
        }
    }

    if (debug > 2)
    {
        Info << " Identified tF as: " << tF << endl;
        Info << " Triangulation: "
             << triangulations[checkIndex][0][tF] << " "
             << triangulations[checkIndex][1][tF] << " "
             << triangulations[checkIndex][2][tF] << " "
             << endl;
        Info << " All triangulations: " << nl
             << ' ' << triangulations[checkIndex][0] << nl
             << ' ' << triangulations[checkIndex][1] << nl
             << ' ' << triangulations[checkIndex][2] << nl
             << endl;
    }

    if (coupledModification_)
    {
        if (locallyCoupledEdge(eIndex))
        {
            // Flip the slave edge as well.
            label slaveIndex = -1;

            // Determine the slave index.
            forAllIter(Map<coupledPatchInfo>, patchCoupling_, patchI)
            {
                if ((slaveIndex = patchI().findSlaveIndex(eIndex)) > -1)
                {
                    break;
                }
            }

            if (debug > 2)
            {
                Info << nl << "Removing slave edge: " << slaveIndex
                     << " for master edge: " << eIndex << endl;
            }

            // Turn off switch temporarily.
            unsetCoupledModification();
            setSlaveModification();

            // Recursively call for the slave edge.
            bool success =
            (
                removeEdgeFlips(slaveIndex, minQuality, K, triangulations, 1)
            );

            // Turn it back on.
            setCoupledModification();
            unsetSlaveModification();

            // Bail out if the slave failed.
            if (!success)
            {
                // Reset all triangulations and bail out
                triangulations[checkIndex][0] = -1;
                triangulations[checkIndex][1] = -1;
                triangulations[checkIndex][2] = -1;

                return false;
            }
        }
    }

    // Perform a series of 2-3 swaps
    label numSwaps = 0;

    while (numSwaps < (m-3))
    {
        for (label i = 0; i < (m-2); i++)
        {
            if ( (i != tF) && (triangulations[checkIndex][0][i] != -1) )
            {
                // Check if triangulation is on the boundary
                if
                (
                    boundaryTriangulation
                    (
                        i,
                        isolatedVertex,
                        triangulations[checkIndex]
                    )
                )
                {
                    // Perform 2-3 swap
                    map =
                    (
                        swap23
                        (
                            isolatedVertex,
                            eIndex,
                            i,
                            numTriangulations,
                            triangulations[checkIndex],
                            hullVertices,
                            hullFaces,
                            hullCells
                        )
                    );

                    if (debug > 2)
                    {
                        scalar triQuality =
                        (
                            Foam::min
                            (
                                tetQuality(owner_[map.opposingFace()]),
                                tetQuality(neighbour_[map.opposingFace()])
                            )
                        );

                        swapQuality = Foam::min(triQuality, swapQuality);
                    }

                    // Done with this face, so reset it
                    triangulations[checkIndex][0][i] = -1;
                    triangulations[checkIndex][1][i] = -1;
                    triangulations[checkIndex][2][i] = -1;

                    numSwaps++;
                }
            }
        }

        if (numSwaps == 0)
        {
            Info << "Triangulations: " << endl;
            forAll(triangulations[checkIndex], row)
            {
                Info << triangulations[checkIndex][row] << endl;
            }

            // Should have performed at least one swap
            FatalErrorIn("dynamicTopoFvMesh::removeEdgeFlips()") << nl
                << "Did not perform any 2-3 swaps" << nl
                << abort(FatalError);
        }
    }

    // Perform the final 3-2 / 2-2 swap
    map =
    (
        swap32
        (
            eIndex,
            tF,
            numTriangulations,
            triangulations[checkIndex],
            hullVertices,
            hullFaces,
            hullCells
        )
    );

    if (debug > 2)
    {
        scalar triQuality =
        (
            Foam::min
            (
                tetQuality(owner_[map.opposingFace()]),
                tetQuality(neighbour_[map.opposingFace()])
            )
        );

        swapQuality = Foam::min(triQuality, swapQuality);

        if (swapQuality < minQuality)
        {
            WarningIn("dynamicTopoFvMesh::removeEdgeFlips()") << nl
                << " Swap failed to improve quality." << nl
                << " MinQuality: " << minQuality << nl
                << " SwapQuality: " << swapQuality
                << abort(FatalError);
        }
    }

    // Done with this face, so reset it
    triangulations[checkIndex][0][tF] = -1;
    triangulations[checkIndex][1][tF] = -1;
    triangulations[checkIndex][2][tF] = -1;

    // Finally remove the edge
    removeEdge(eIndex);

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    nSwaps_++;

    // Return a successful operation.
    return true;
}

// Extract triangulations from the programming table
void dynamicTopoFvMesh::extractTriangulation
(
    const label i,
    const label j,
    const labelListList& K,
    label& numTriangulations,
    labelListList& triangulations
) const
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
    const labelListList& triangulations,
    const scalar tolFraction
) const
{
    label m = hullVertices.size();
    scalar tolerance = VSMALL;

    const edge& edgeToCheck = edges_[eIndex];
    FixedList<label, 3> sign(-2);
    FixedList<scalar, 3> vol(0.0);

    // Relax the tolerance for boundary edges
    if (whichEdgePatch(eIndex) > -1)
    {
        tolerance = tolFraction*mag(tangentToEdge(eIndex));
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
                tolerance,
                vol[0]
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
                tolerance,
                vol[1]
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
                tolerance,
                vol[2]
            )
        );

        if (debug > 2)
        {
            Info << " tetVolumeSign: " << sign << endl;
            Info << " tetVolume: " << vol << endl;
        }

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
            << endl;
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
) const
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

// Output a list of primitives as a VTK file.
// Uses the current state of connectivity.
// primitiveType is:
//   0: List of points
//   1: List of edges
//   2: List of faces
//   3: List of cells
//   4: List of cells w/ associated connectivity as fields
void dynamicTopoFvMesh::writeVTK
(
    const word& name,
    const labelList& cList,
    const label primitiveType
) const
{
    label nTotalCells = 0;
    label nPoints = 0, nEdges = 0, nFaces = 0, nCells = 0;

    // Estimate a size for points and cellPoints
    List<vector> points(6*cList.size());

    // Connectivity lists
    labelListList cpList(cList.size());
    labelListList feList, epList, fpList;

    // Create a map for local points
    Map<label> pointMap, edgeMap, faceMap;

    // Track surface-points, if requested.
    labelHashSet surfPoints;

    forAll(cList, cellI)
    {
        if (cList[cellI] < 0)
        {
            continue;
        }

        // Are we looking at points?
        if (primitiveType == 0)
        {
            // Size the list
            cpList[nCells].setSize(1);

            cpList[nCells] = cList[cellI];

            nTotalCells++;
        }

        // Are we looking at edges?
        if (primitiveType == 1)
        {
            // Size the list
            cpList[nCells].setSize(2);

            const edge& thisEdge = edges_[cList[cellI]];

            cpList[nCells][0] = thisEdge[0];
            cpList[nCells][1] = thisEdge[1];

            nTotalCells += 2;
        }

        // Are we looking at faces?
        if (primitiveType == 2)
        {
            const face& thisFace = faces_[cList[cellI]];

            if (thisFace.size() == 3)
            {
                // Size the list
                cpList[nCells].setSize(3);

                // Write out in order
                cpList[nCells][0] = thisFace[0];
                cpList[nCells][1] = thisFace[1];
                cpList[nCells][2] = thisFace[2];

                nTotalCells += 3;
            }
            else
            if (thisFace.size() == 4)
            {
                // Size the list
                cpList[nCells].setSize(4);

                // Write out in order
                cpList[nCells][0] = thisFace[0];
                cpList[nCells][1] = thisFace[1];
                cpList[nCells][2] = thisFace[2];
                cpList[nCells][3] = thisFace[3];

                nTotalCells += 4;
            }
        }

        // Are we looking at cells?
        if (primitiveType == 3 || primitiveType == 4)
        {
            const cell& thisCell = cells_[cList[cellI]];

            if (thisCell.size() == 4)
            {
                // Point-ordering for tetrahedra
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

                nTotalCells += 4;
            }
            else
            if (thisCell.size() == 5)
            {
                // Point-ordering for wedge cells
                label firstTriFace = -1;

                // Size the list
                cpList[nCells].setSize(6);

                // Figure out one triangle face
                forAll(thisCell, faceI)
                {
                    const face& currFace = faces_[thisCell[faceI]];

                    if (currFace.size() == 3)
                    {
                        if (firstTriFace == -1)
                        {
                            firstTriFace = thisCell[faceI];

                            // Right-handedness is assumed here.
                            // Tri-faces are always on the boundary.
                            cpList[nCells][0] = currFace[0];
                            cpList[nCells][1] = currFace[1];
                            cpList[nCells][2] = currFace[2];
                        }
                        else
                        {
                            // Detect the three other points.
                            forAll(thisCell, faceJ)
                            {
                                const face& nextFace = faces_[thisCell[faceJ]];

                                if (nextFace.size() == 4)
                                {
                                    // Search for vertices on currFace
                                    // in this face.
                                    label i = -1, p = -1, n = -1;

                                    if ((i=nextFace.which(currFace[0])) != -1)
                                    {
                                        p = nextFace.prevLabel(i);
                                        n = nextFace.nextLabel(i);

                                        if
                                        (
                                            p != currFace[1] &&
                                            p != currFace[2]
                                        )
                                        {
                                            cpList[nCells][3] = p;
                                        }
                                        else
                                        if
                                        (
                                            n != currFace[1] &&
                                            n != currFace[2]
                                        )
                                        {
                                            cpList[nCells][3] = n;
                                        }
                                    }

                                    if ((i=nextFace.which(currFace[1])) != -1)
                                    {
                                        p = nextFace.prevLabel(i);
                                        n = nextFace.nextLabel(i);

                                        if
                                        (
                                            p != currFace[0] &&
                                            p != currFace[2]
                                        )
                                        {
                                            cpList[nCells][4] = p;
                                        }
                                        else
                                        if
                                        (
                                            n != currFace[0] &&
                                            n != currFace[2]
                                        )
                                        {
                                            cpList[nCells][4] = n;
                                        }
                                    }

                                    if ((i=nextFace.which(currFace[2])) != -1)
                                    {
                                        p = nextFace.prevLabel(i);
                                        n = nextFace.nextLabel(i);

                                        if
                                        (
                                            p != currFace[0] &&
                                            p != currFace[1]
                                        )
                                        {
                                            cpList[nCells][5] = p;
                                        }
                                        else
                                        if
                                        (
                                            n != currFace[0] &&
                                            n != currFace[1]
                                        )
                                        {
                                            cpList[nCells][5] = n;
                                        }
                                    }
                                }
                            }

                            break;
                        }
                    }
                }

                nTotalCells += 6;
            }

            // Add to the list of surface points
            if (primitiveType == 4)
            {
                forAll(thisCell, faceI)
                {
                    const face& thisFace = faces_[thisCell[faceI]];

                    if (whichPatch(thisCell[faceI]) > -1)
                    {
                        forAll(thisFace, pointI)
                        {
                            if (!surfPoints.found(thisFace[pointI]))
                            {
                                surfPoints.insert(thisFace[pointI]);
                            }
                        }
                    }

                    // Check if this face was added to the map
                    if (!faceMap.found(thisCell[faceI]))
                    {
                        faceMap.insert(thisCell[faceI], nFaces);

                        nFaces++;
                    }

                    // Check if edges of this face were added
                    const labelList& fEdges = faceEdges_[thisCell[faceI]];

                    forAll(fEdges, edgeI)
                    {
                        if (!edgeMap.found(fEdges[edgeI]))
                        {
                            edgeMap.insert(fEdges[edgeI], nEdges);

                            nEdges++;
                        }
                    }
                }
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

        nCells++;
    }

    // Make the directory
    fileName dirName(time().path()/"VTK"/time().timeName());

    mkDir(dirName);

    // Open stream for output
    OFstream file(dirName/name+".vtk");

    // Write out the header
    file << "# vtk DataFile Version 2.0" << nl
         << name << ".vtk" << nl
         << "ASCII" << nl
         << "DATASET UNSTRUCTURED_GRID" << nl
         << "POINTS " << nPoints << " double" << nl;

    for (label i = 0; i < nPoints; i++)
    {
        file << setprecision(10)
             << points[i].x() << ' '
             << points[i].y() << ' '
             << points[i].z() << ' '
             << nl;
    }

    file << "CELLS " << nCells << " " << nTotalCells + nCells << endl;

    forAll(cpList, i)
    {
        if (cpList[i].size())
        {
            file << cpList[i].size() << ' ';

            forAll(cpList[i], j)
            {
                file << cpList[i][j] << ' ';
            }

            file << nl;
        }
    }

    file << "CELL_TYPES " << nCells << endl;

    forAll(cpList, i)
    {
        if (cpList[i].size() == 1)
        {
            // Vertex
            file << "1" << nl;
        }

        if (cpList[i].size() == 2)
        {
            // Edge
            file << "3" << nl;
        }

        if (cpList[i].size() == 3)
        {
            // Triangle face
            file << "5" << nl;
        }

        if
        (
            (cpList[i].size() == 4) &&
            (primitiveType == 2)
        )
        {
            // Quad face
            file << "9" << nl;
        }

        if
        (
            (cpList[i].size() == 4) &&
            (primitiveType == 3 || primitiveType == 4)
        )
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

    if (primitiveType < 4)
    {
        return;
    }

    // Write out auxiliary connectivity fields, if necessary
    label nAuxFields = 0;

    if (surfPoints.size())
    {
        nAuxFields++;
    }

    if (nCells)
    {
        nAuxFields++;
    }

    if (nFaces)
    {
        nAuxFields += 2;

        // Size the feList and fpList.
        // Assume a uniform mesh
        fpList.setSize(nFaces, labelList(3, -1));
        feList.setSize(nFaces, labelList(3, -1));

        forAllIter(Map<label>, faceMap, fIter)
        {
            const face& thisFace = faces_[fIter.key()];

            forAll(thisFace, pointI)
            {
                fpList[fIter()][pointI] = pointMap[thisFace[pointI]];
            }

            const labelList& fEdges = faceEdges_[fIter.key()];

            forAll(fEdges, edgeI)
            {
                feList[fIter()][edgeI] = edgeMap[fEdges[edgeI]];
            }
        }
    }

    if (nEdges)
    {
        nAuxFields++;

        // Size the epList
        epList.setSize(nEdges, labelList(2, -1));

        forAllIter(Map<label>, edgeMap, eIter)
        {
            const edge& thisEdge = edges_[eIter.key()];

            forAll(thisEdge, pointI)
            {
                epList[eIter()][pointI] = pointMap[thisEdge[pointI]];
            }
        }
    }

    if (nAuxFields)
    {
        file << "FIELD ConnData " << nAuxFields << endl;

        if (nFaces)
        {
            // Write out cellFaces.
            file << "cellFaces " << 4 << ' ' << nCells << " int" << endl;

            forAll(cList, cellI)
            {
                if (cList[cellI] < 0)
                {
                    continue;
                }

                const cell& thisCell = cells_[cList[cellI]];

                forAll(thisCell, faceI)
                {
                    file << faceMap[thisCell[faceI]] << " ";
                }

                file << nl;
            }

            // Write out facePoints
            file << "facePoints " << 3 << ' ' << nFaces << " int" << endl;

            forAll(fpList, faceI)
            {
                const labelList& list = fpList[faceI];

                forAll(list, indexI)
                {
                    file << list[indexI] << " ";
                }

                file << nl;
            }

            // Write out faceEdges
            file << "faceEdges " << 3 << ' ' << nFaces << " int" << endl;

            forAll(feList, faceI)
            {
                const labelList& list = feList[faceI];

                forAll(list, indexI)
                {
                    file << list[indexI] << " ";
                }

                file << nl;
            }
        }

        // Write out edges
        if (nEdges)
        {
            file << "edgePoints " << 2 << ' ' << nEdges << " int" << endl;

            forAll(epList, edgeI)
            {
                const labelList& list = epList[edgeI];

                forAll(list, indexI)
                {
                    file << list[indexI] << " ";
                }

                file << nl;
            }
        }

        // Write out surface points
        if (surfPoints.size())
        {
            file << "surfacePoints 1 " << surfPoints.size() << " int" << endl;

            // Write out surface points
            forAllIter(labelHashSet, surfPoints, spIter)
            {
                file << pointMap[spIter.key()] << nl;
            }
        }
    }

    file << nl;
}

// Check the state of connectivity lists
void dynamicTopoFvMesh::checkConnectivity
(
    label maxErrors
) const
{
    label nFailedChecks = 0;

    messageStream ConnectivityWarning
    (
        "dynamicTopoFvMesh Connectivity Warning",
        messageStream::SERIOUS,
        maxErrors
    );

    // Check face-label ranges
    Info << "Checking index ranges...";

    forAll(edges_, edgeI)
    {
        const edge& curEdge = edges_[edgeI];

        if (curEdge == edge(-1, -1))
        {
            continue;
        }

        if
        (
            curEdge[0] < 0 || curEdge[0] > (points_.size()-1) ||
            curEdge[1] < 0 || curEdge[1] > (points_.size()-1)
        )
        {
            Pout << "Edge " << edgeI
                 << " contains vertex labels out of range: "
                 << curEdge
                 << " Max point index = " << (points_.size()-1) << endl;

            nFailedChecks++;

            ConnectivityWarning()
                << nl << "Edge-point connectivity is inconsistent."
                << endl;
        }
    }

    forAll(faces_, faceI)
    {
        const face& curFace = faces_[faceI];

        if (curFace.empty())
        {
            continue;
        }

        if (min(curFace) < 0 || max(curFace) > (points_.size()-1))
        {
            Pout << "Face " << faceI
                 << " contains vertex labels out of range: "
                 << curFace
                 << " Max point index = " << (points_.size()-1) << endl;

            nFailedChecks++;

            ConnectivityWarning()
                << nl << "Face-point connectivity is inconsistent."
                << endl;
        }
    }

    forAll(cells_, cellI)
    {
        const cell& curCell = cells_[cellI];

        if (curCell.empty())
        {
            continue;
        }

        if (min(curCell) < 0 || max(curCell) > (faces_.size()-1))
        {
            Pout << "Cell " << cellI
                 << " contains vertex labels out of range: "
                 << curCell
                 << " Max point index = " << (faces_.size()-1) << endl;

            nFailedChecks++;

            ConnectivityWarning()
                << nl << "Cell-Face connectivity is inconsistent."
                << endl;
        }
    }

    Info << "Done." << endl;

    Info << "Checking edge-face connectivity...";

    label allEdges = edges_.size();
    labelList nEdgeFaces(allEdges, 0);

    forAll(faceEdges_, faceI)
    {
        const labelList& faceEdges = faceEdges_[faceI];

        if (faceEdges.empty())
        {
            continue;
        }

        // Check consistency of face-edge-points as well
        edgeList eList = faces_[faceI].edges();

        forAll(faceEdges,edgeI)
        {
            nEdgeFaces[faceEdges[edgeI]]++;

            // Check if this edge actually belongs to this face
            bool found = false;
            const edge& edgeToCheck = edges_[faceEdges[edgeI]];

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
                     << "was not found in face: " << faceI
                     << ": " << faces_[faceI] << nl
                     << "faceEdges: " << faceEdges
                     << endl;

                nFailedChecks++;

                ConnectivityWarning()
                     << nl << "Edge-Face connectivity is inconsistent."
                     << endl;
            }
        }
    }

    label nInternalEdges = 0;
    labelList patchInfo(numPatches_, 0);

    forAll(edgeFaces_, edgeI)
    {
        const labelList& edgeFaces = edgeFaces_[edgeI];

        if (edgeFaces.empty())
        {
            continue;
        }

        if (edgeFaces.size() != nEdgeFaces[edgeI])
        {
            Pout << nl << nl << "Edge: " << edgeI
                 << ": edgeFaces: " << edgeFaces << endl;

            nFailedChecks++;

            ConnectivityWarning()
                << nl << "Edge-Face connectivity is inconsistent."
                << endl;
        }

        label nBF = 0;

        // Check if this edge belongs to faceEdges for each face
        forAll(edgeFaces, faceI)
        {
            if (findIndex(faceEdges_[edgeFaces[faceI]], edgeI) == -1)
            {
                Pout << nl << nl << "Edge: " << edgeI << ": " << edges_[edgeI]
                     << ", edgeFaces: " << edgeFaces << nl
                     << "was not found in faceEdges of face: "
                     << edgeFaces[faceI] << ": " << faces_[edgeFaces[faceI]]
                     << nl << "faceEdges: " << faceEdges_[edgeFaces[faceI]]
                     << endl;

                nFailedChecks++;

                ConnectivityWarning()
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
            if (whichEdgePatch(edgeI) >= 0)
            {
                Pout << "Edge: " << edgeI
                     << ": " << edges_[edgeI] << " is internal, "
                     << " but patch is specified as: "
                     << whichEdgePatch(edgeI)
                     << endl;

                nFailedChecks++;
            }
        }
        else
        {
            label patchID = whichEdgePatch(edgeI);

            // Check if this edge is actually on a boundary.
            if (patchID < 0)
            {
                Pout << "Edge: " << edgeI
                     << ": " << edges_[edgeI]
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
    forAllConstIter(Map<label>, addedEdgePatches_, aepIter)
    {
        label key = aepIter.key();
        label patch = aepIter();

        label nBF = 0;
        const labelList& edgeFaces = edgeFaces_[key];

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

        label allPoints = points_.size();
        List<labelHashSet> hlPointEdges(allPoints);

        forAll(edges_, edgeI)
        {
            if (edgeFaces_[edgeI].size())
            {
                hlPointEdges[edges_[edgeI][0]].insert(edgeI);
                hlPointEdges[edges_[edgeI][1]].insert(edgeI);
            }
        }

        forAll(pointEdges_, pointI)
        {
            const labelList& pointEdges = pointEdges_[pointI];

            if (pointEdges.empty())
            {
                continue;
            }

            forAll(pointEdges, edgeI)
            {
                if (!hlPointEdges[pointI].found(pointEdges[edgeI]))
                {
                    Pout << nl << nl << "Point: " << pointI << nl
                         << "pointEdges: " << pointEdges << nl
                         << "hlPointEdges: " << hlPointEdges[pointI]
                         << endl;

                    nFailedChecks++;

                    ConnectivityWarning()
                        << nl << "Point-Edge connectivity is inconsistent."
                        << endl;
                }
            }

            // Do a size check as well
            if
            (
                hlPointEdges[pointI].size() != pointEdges.size() ||
                pointEdges.size() == 1
            )
            {
                Pout << nl << nl << "Point: " << pointI << nl
                     << "pointEdges: " << pointEdges << nl
                     << "hlPointEdges: " << hlPointEdges[pointI]
                     << endl;

                nFailedChecks++;

                ConnectivityWarning()
                    << "Size inconsistency."
                    << nl << "Point-Edge connectivity is inconsistent."
                    << endl;
            }
        }

        Info << "Done." << endl;

        Info << "Checking edge-points connectivity...";

        label otherPoint = -1, nextPoint = -1;

        forAll(edgePoints_, edgeI)
        {
            // Do a preliminary size check
            const labelList& edgePoints = edgePoints_[edgeI];
            const labelList& edgeFaces = edgeFaces_[edgeI];

            if (edgeFaces.empty())
            {
                continue;
            }

            if (edgePoints.size() != edgeFaces.size())
            {
                Pout << nl << nl
                     << "Edge: " << edgeI
                     << " " << edges_[edgeI] << endl;

                Pout << "edgeFaces: " << edgeFaces << endl;
                forAll(edgeFaces, faceI)
                {
                    Info << edgeFaces[faceI] << ": "
                         << faces_[edgeFaces[faceI]]
                         << endl;
                }

                Pout << "edgePoints: " << edgePoints << endl;

                nFailedChecks++;

                ConnectivityWarning()
                    << nl << "Edge-Points connectivity is inconsistent."
                    << endl;
            }

            // Now check to see that both lists are consistent.
            const edge& edgeToCheck = edges_[edgeI];

            forAll(edgeFaces, faceI)
            {
                const face& faceToCheck = faces_[edgeFaces[faceI]];

                findIsolatedPoint
                (
                    faceToCheck,
                    edgeToCheck,
                    otherPoint,
                    nextPoint
                );

                if (findIndex(edgePoints, otherPoint) == -1)
                {
                    Pout << nl << nl
                         << "Edge: " << edgeI
                         << " " << edges_[edgeI] << endl;

                    Pout << "edgeFaces: " << edgeFaces << endl;
                    forAll(edgeFaces, faceI)
                    {
                        Info << edgeFaces[faceI] << ": "
                             << faces_[edgeFaces[faceI]]
                             << endl;
                    }

                    Pout << "edgePoints: " << edgePoints << endl;

                    nFailedChecks++;

                    ConnectivityWarning()
                        << nl << "Edge-Points connectivity is inconsistent."
                        << endl;
                }
            }
        }

        Info << "Done." << endl;
    }

    Info << "Checking cell-point connectivity...";

    // Loop through all cells and construct cell-to-node
    label cIndex = 0;
    label allCells = cells_.size();
    labelList cellIndex(allCells);
    List<labelHashSet> cellToNode(allCells);

    forAll(cells_, cellI)
    {
        const cell& thisCell = cells_[cellI];

        if (thisCell.empty())
        {
            continue;
        }

        cellIndex[cIndex] = cellI;

        forAll(thisCell, faceI)
        {
            const labelList& fEdges = faceEdges_[thisCell[faceI]];

            forAll(fEdges, edgeI)
            {
                const edge& thisEdge = edges_[fEdges[edgeI]];

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
        if
        (
            (cellToNode[cellI].size() != 6 && twoDMesh_) ||
            (cellToNode[cellI].size() != 4 && !twoDMesh_)
        )
        {
            Pout << nl << "Warning: Cell: "
                 << cellIndex[cellI] << " is inconsistent. "
                 << endl;

            const cell& failedCell = cells_[cellIndex[cellI]];

            Info << "Cell faces: " << failedCell << endl;

            forAll(failedCell, faceI)
            {
                Info << "\tFace: " << failedCell[faceI]
                     << " :: " << faces_[failedCell[faceI]]
                     << endl;

                const labelList& fEdges = faceEdges_[failedCell[faceI]];

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

// Perform spatial hashing on a set of points
void dynamicTopoFvMesh::spatialHash
(
    const pointField& pointLocations,
    const labelList& pointIndices,
    const boundBox& box,
    const label resolution,
    labelListList& bins,
    label removeIndex
) const
{
    label binSize = bins.size(), nD = resolution;

    const point& bMin = box.min();
    const point& bMax = box.max();

    // Extend bounding-box dimensions a bit to avoid edge-effects.
    scalar ext = 0.02*(mag(bMax - bMin));

    // Define an inverse grid-cell size.
    scalar xL = nD/(bMax.x() - bMin.x() + ext);
    scalar yL = nD/(bMax.y() - bMin.y() + ext);
    scalar zL = nD/(bMax.z() - bMin.z() + ext);

    // Loop through all points and bin them.
    forAll(pointLocations, pointI)
    {
        // Translate to boundBox minimum.
        point p = pointLocations[pointI] - bMin;

        // Hash the position.
        label i = label(mag(::floor(p.x()*xL)));
        label j = label(mag(::floor(p.y()*yL)));
        label k = label(mag(::floor(p.z()*zL)));

        label pos = ((k*nD*nD)+(j*nD)+i) % binSize;

        if (removeIndex)
        {
            // Remove the index.
            sizeDownList
            (
                pointIndices[pointI],
                bins[pos]
            );
        }
        else
        {
            // Store the index.
            sizeUpList
            (
                pointIndices[pointI],
                bins[pos]
            );
        }
    }
}

// Prepare for proximity-based refinement, if necessary
void dynamicTopoFvMesh::prepareProximityPatches()
{
    if (!proximityPatches_.size())
    {
        return;
    }

    if (debug)
    {
        Info << "Preparing patches for proximity-based refinement...";
    }

    // Clear out existing lists.
    proximityBins_.clear();

    proximityBins_.setSize(997, labelList(0));

    const polyBoundaryMesh& boundary = polyMesh::boundaryMesh();

    bool setSpatialRes = false;

    // Loop through all proximity patches and spatially hash patch faces.
    forAll(boundary, patchI)
    {
        if (proximityPatches_.found(boundary[patchI].name()))
        {
            const polyPatch& proxPatch = boundary[patchI];

            // Construct a bounding-box of face centres.
            // Do not synchronize in parallel, since the patch
            // may not be present on all sub-domains.
            boundBox box(proxPatch.faceCentres(), false);

            const point& bMin = box.min();
            const point& bMax = box.max();

            // Further hashing requires this information.
            proxBoundBox_.min() = bMin;
            proxBoundBox_.max() = bMax;

            // Build a list of face indices
            labelList faceIndices
            (
                identity(proxPatch.size()) + proxPatch.start()
            );

            // For spatial resolution, pick an edge on this patch.
            if (!setSpatialRes)
            {
                spatialRes_ =
                (
                    label
                    (
                        ::floor
                        (
                            mag(bMax - bMin)
                          / (3.0*edgeLength(faceEdges_[proxPatch.start()][0]))
                        )
                    )
                );

                setSpatialRes = true;
            }

            spatialHash
            (
                proxPatch.faceCentres(),
                faceIndices,
                box,
                spatialRes_,
                proximityBins_
            );
        }
    }

    if (debug)
    {
        Info << "Done." << endl;
    }
}

void dynamicTopoFvMesh::handleMeshSlicing()
{
    if (slicePairs_.empty())
    {
        return;
    }

    if (sliceHoldOff_)
    {
        // Hold-off mesh slicing for a few time-steps.
        sliceHoldOff_--;
        return;
    }

    Info << "Slicing Mesh...";

    // Loop through candidates and weed-out invalid points
    forAll(slicePairs_, pairI)
    {
        const labelPair& pairToCheck = slicePairs_[pairI];

        bool available = true;

        forAll(pairToCheck, indexI)
        {
            if (deletedPoints_.found(pairToCheck[indexI]))
            {
                available = false;
                break;
            }

            if (pairToCheck[indexI] < nOldPoints_)
            {
                if (reversePointMap_[pairToCheck[indexI]] == -1)
                {
                    available = false;
                    break;
                }
            }
        }

        if (available)
        {
            // Slice the mesh at this point.
            sliceMesh(pairToCheck);
        }
    }

    Info << "Done." << endl;

    checkConnectivity(10);

    // Clear out data.
    sliceBoxes_.clear();
    slicePairs_.clear();

    // Set the sliceHoldOff value
    sliceHoldOff_ = 50;
}

// Test an edge for proximity with other faces on proximity patches
// and return the scalar distance to an oppositely-oriented face.
scalar dynamicTopoFvMesh::testProximity
(
    const label eIndex,
    label& proximityFace
) const
{
    const edge& thisEdge = edges_[eIndex];
    const labelList& eFaces = edgeFaces_[eIndex];

    // Reset the proximity face
    proximityFace = -1;

    // Obtain the edge centre.
    point eCentre =
    (
        0.5 * (points_[thisEdge[0]] + points_[thisEdge[1]])
    );

    // Obtain the edge-normal
    vector eNorm = vector::zero;

    forAll(eFaces, faceI)
    {
        if (neighbour_[eFaces[faceI]] == -1)
        {
            // Obtain the normal.
            eNorm += triFaceNormal(faces_[eFaces[faceI]]);
        }
    }

    eNorm /= (mag(eNorm) + VSMALL);

    DynamicList<label> posIndices(20);
    scalar testStep = edgeLength(eIndex);
    scalar minDistance = GREAT, minDeviation = -0.9;
    label nD = spatialRes_, binSize = proximityBins_.size();

    const point& bMin = proxBoundBox_.min();
    const point& bMax = proxBoundBox_.max();

    // Extend bounding-box dimensions a bit to avoid edge-effects.
    scalar ext = 0.02*(mag(bMax - bMin));

    // Define an inverse grid-cell size.
    scalar xL = nD/(bMax.x() - bMin.x() + ext);
    scalar yL = nD/(bMax.y() - bMin.y() + ext);
    scalar zL = nD/(bMax.z() - bMin.z() + ext);

    // Now take multiple steps in both edge-normal directions,
    // and add to the list of boxes to be checked.
    for (scalar dir = -1.0; dir < 2.0; dir += 2.0)
    {
        for (scalar step = 0.0; step < 5.0*testStep; step += testStep)
        {
            // Hash the point-location
            point p = (eCentre + (dir*step*eNorm)) - bMin;

            label i = label(mag(::floor(p.x()*xL)));
            label j = label(mag(::floor(p.y()*yL)));
            label k = label(mag(::floor(p.z()*zL)));

            label pos = ((k*nD*nD)+(j*nD)+i) % binSize;

            if (findIndex(posIndices, pos) == -1)
            {
                posIndices.append(pos);
            }
        }
    }

    // Obtain old-mesh face geometry for reference.
    const vectorField& faceAreas = primitiveMesh::faceAreas();
    const vectorField& faceCentres = primitiveMesh::faceCentres();

    forAll(posIndices, indexI)
    {
        const labelList& posBin = proximityBins_[posIndices[indexI]];

        forAll(posBin, faceI)
        {
            // Step 1: Measure the distance to the face.
            vector rFace = (faceCentres[posBin[faceI]] - eCentre);

            scalar distance = mag(rFace);

            // Step 2: Check if this face is oriented away from edge.
            const vector& fNorm = faceAreas[posBin[faceI]];

            scalar deviation = (eNorm & (fNorm/mag(fNorm)));

            if
            (
                (deviation < minDeviation) &&
                (distance < minDistance)
            )
            {
                // Update statistics
                proximityFace = posBin[faceI];
                minDistance = distance;
                // minDeviation = deviation;

                // Define whether rFace went through
                // domain by comparing with the edge-normal.
                if ((eNorm & (rFace/distance)) > 0.0)
                {
                    // Outside the domain
                }
                else
                {
                    // Inside the domain
                }
            }
        }
    }

    if (proximityFace != -1)
    {
        // Check if we need to mark points for mesh-slicing.
        if (minDistance < sliceThreshold_)
        {
            // Check if any points on this face are still around.
            // If yes, mark one of them as the end point
            // for Dijkstra's algorithm. The start point will be a point
            // on this edge.
            labelPair proxPoints(thisEdge[0], -1);

            const face& proxFace = polyMesh::faces()[proximityFace];

            bool foundPoint = false;

            forAll(proxFace, pointI)
            {
                if (reversePointMap_[proxFace[pointI]] != -1)
                {
                    proxPoints.second() = proxFace[pointI];
                    foundPoint = true;
                    break;
                }
            }

            if (foundPoint)
            {
                // Lock the edge mutex
                entityMutex_[1].lock();

                // Add this entry as a candidate for mesh slicing.
                label curSize = slicePairs_.size();

                slicePairs_.setSize(curSize + 1);

                slicePairs_[curSize] = proxPoints;

                // Unlock the edge mutex
                entityMutex_[1].unlock();
            }
        }
    }

    return minDistance;
}

// Calculate the edge length-scale for the mesh
void dynamicTopoFvMesh::calculateLengthScale()
{
    if (!edgeRefinement_)
    {
        return;
    }

    label level = 1, visitedCells = 0;
    labelList cellLevels(nCells(), 0);

    // Size the local field
    lengthScale_.setSize(nCells(), 0.0);

    // HashSet to keep track of cells in each level
    labelHashSet levelCells;

    // Prepare for proximity-based refinement, if necessary
    prepareProximityPatches();

    // Obtain the cellCells addressing list
    const labelListList& cc = polyMesh::cellCells();
    const polyBoundaryMesh& boundary = polyMesh::boundaryMesh();
    const labelList& own = polyMesh::faceOwner();

    forAll(boundary,patchI)
    {
        const polyPatch& bdyPatch = boundary[patchI];

        if
        (
            (!freePatches_.found(bdyPatch.name())) &&
            (bdyPatch.type() != "processor") &&
            (bdyPatch.type() != "cyclic") &&
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

                lengthScale_[ownCell] =
                (
                    boundaryLengthScale(pStart+faceI)*growthFactor_
                );

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
                lengthScale_
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
                            sumLength += lengthScale_[ncList[indexJ]];

                            nTouchedNgb++;
                        }
                    }

                    sumLength /= nTouchedNgb;

                    // Scale the length and assign to this cell
                    scalar sLength = sumLength*growthFactor_;

                    lengthScale_[cList[indexI]] = sLength;

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
                lengthScale_,
                levelCells
            );
        }

        if (debug > 4)
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
}

// Calculate geometric tolerance for the mesh
void dynamicTopoFvMesh::calculateGeometricTolerance()
{

}

// Compute the growth factor of an existing mesh
scalar dynamicTopoFvMesh::computeGrowthFactor()
{
    if (!edgeRefinement_)
    {
        return -1.0;
    }

    scalar growthFactor = 0.0;

    label level = 1, visitedCells = 0;
    labelList cellLevels(nCells(),0);

    // Obtain the addressing lists
    const edgeList& edges = polyMesh::edges();
    const pointField& points = polyMesh::points();
    const labelList& owner = polyMesh::faceOwner();
    const labelListList& cc = polyMesh::cellCells();
    const labelListList& cEdges = polyMesh::cellEdges();
    const polyBoundaryMesh& boundary = polyMesh::boundaryMesh();

    // HashSet to keep track of cells in each level
    labelHashSet levelCells;

    // Obtain the list of patches for which the length-scale is fixed
    wordList toc = fixedPatches_.toc();

    forAll(boundary,patchI)
    {
        const polyPatch& bdyPatch = boundary[patchI];

        // Loop through all fixed length-scale patches
        forAll(toc,wordI)
        {
            const word& pName = toc[wordI];

            if (boundary[patchI].name() == pName)
            {
                label pStart = bdyPatch.start();

                forAll(bdyPatch,faceI)
                {
                    label ownCell = owner[pStart+faceI];

                    if (cellLevels[ownCell] != 0)
                    {
                        continue;
                    }

                    cellLevels[ownCell] = level;

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
            (bdyPatch.type() != "processor") &&
            (bdyPatch.type() != "cyclic") &&
            (bdyPatch.type() != "wedge") &&
            (bdyPatch.type() != "empty") &&
            (bdyPatch.type() != "symmetryPlane")
        )
        {
            label pStart = bdyPatch.start();

            forAll(bdyPatch,faceI)
            {
                label ownCell = owner[pStart+faceI];

                if (cellLevels[ownCell] != 0)
                {
                    continue;
                }

                cellLevels[ownCell] = level;

                levelCells.insert(ownCell);

                visitedCells++;
            }
        }
    }

    scalar prevAvg = 0.0;

    while (visitedCells < nCells())
    {
        // Loop through cells of the current level
        labelList currLvlCells = levelCells.toc();

        levelCells.clear();

        scalar avgEdgeLength = 0.0;

        // Loop through cells, and increment neighbour
        // cells of the current level
        forAll(currLvlCells,cellI)
        {
            const labelList& eList = cEdges[currLvlCells[cellI]];

            scalar avgCellLength = 0.0;

            forAll(eList, edgeI)
            {
                avgCellLength += edges[eList[edgeI]].mag(points);
            }

            avgCellLength /= eList.size();

            avgEdgeLength += avgCellLength;

            // Obtain the cells neighbouring this one,
            // and increment their level.
            const labelList& cList = cc[currLvlCells[cellI]];

            forAll(cList, indexI)
            {
                label& ngbLevel = cellLevels[cList[indexI]];

                if (ngbLevel == 0)
                {
                    ngbLevel = level + 1;

                    levelCells.insert(cList[indexI]);

                    visitedCells++;
                }
            }
        }

        avgEdgeLength /= currLvlCells.size();

        if (level == 1)
        {
            prevAvg = avgEdgeLength;
        }
        else
        {
            growthFactor += (avgEdgeLength/prevAvg);

            prevAvg = avgEdgeLength;
        }

        // Move on to the next level
        level++;
    }

    // Take the average growth factor
    growthFactor /= (level - 2);

    return growthFactor;
}

// Send length-scale info across processors
void dynamicTopoFvMesh::writeLengthScaleInfo
(
    const labelList& cellLevels,
    const resizableList<scalar>& lengthScale
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
    resizableList<scalar>& lengthScale,
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
) const
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
) const
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

// Parallel non-blocking send for fixed lists
template <class Type, label Size>
void dynamicTopoFvMesh::pWrite
(
    const label toID,
    const FixedList<Type, Size>& data
) const
{
    OPstream::write
    (
        Pstream::blocking,
        toID,
        reinterpret_cast<const char*>(&data[0]),
        Size*sizeof(Type)
    );
}

// Parallel non-blocking receive for fixed lists
template <class Type, label Size>
void dynamicTopoFvMesh::pRead
(
    const label fromID,
    FixedList<Type, Size>& data
) const
{
    IPstream::read
    (
        Pstream::blocking,
        fromID,
        reinterpret_cast<char*>(data.begin()),
        Size*sizeof(Type)
    );
}

// Parallel non-blocking send for lists
template <class Type>
void dynamicTopoFvMesh::pWrite
(
    const label toID,
    const List<Type>& data
) const
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
) const
{
    IPstream::read
    (
        Pstream::nonBlocking,
        fromID,
        reinterpret_cast<char*>(data.begin()),
        data.byteSize()
    );
}

// Read optional dictionary parameters
void dynamicTopoFvMesh::readOptionalParameters()
{
    // Enable/disable run-time debug level
    if (dict_.found("debug") || mandatory_)
    {
        debug = readLabel(dict_.lookup("debug"));
    }
    else
    {
        debug = 0;
    }

    if (dict_.subDict("dynamicTopoFvMesh").found("interval") || mandatory_)
    {
        interval_ =
        (
            readLabel
            (
                dict_.subDict
                ("dynamicTopoFvMesh").lookup("interval")
            )
        );
    }
    else
    {
        interval_ = 1;
    }

    // For tetrahedral meshes...
    if (!twoDMesh_)
    {
        // Check if swapping is to be avoided on any patches
        if
        (
            dict_.subDict("dynamicTopoFvMesh").found("noSwapPatches") ||
            mandatory_
        )
        {
            wordList noSwapPatches =
            (
                dict_.subDict
                (
                    "dynamicTopoFvMesh"
                ).subDict("noSwapPatches").toc()
            );

            // Ensure that patches are legitimate.
            checkPatches(noSwapPatches);

            noSwapPatchIDs_.setSize(noSwapPatches.size());

            label indexI = 0;

            forAll(noSwapPatches, wordI)
            {
                const word& patchName = noSwapPatches[wordI];

                noSwapPatchIDs_[indexI++] =
                (
                    boundaryMesh().findPatchID(patchName)
                );
            }
        }

        // Check if a limit has been imposed on maxTetsPerEdge
        if
        (
            dict_.subDict("dynamicTopoFvMesh").found("maxTetsPerEdge") ||
            mandatory_
        )
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
        if
        (
            dict_.subDict("dynamicTopoFvMesh").found("allowTableResize") ||
            mandatory_
        )
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
        else
        {
            allowTableResize_ = false;
        }
    }
}

// Read edge refinement options from the dictionary
void dynamicTopoFvMesh::readRefinementOptions
(
    bool reRead
)
{
    if (!edgeRefinement_)
    {
        return;
    }

    const dictionary& edgeOptionDict =
    (
        dict_.subDict("dynamicTopoFvMesh").subDict("refinementOptions")
    );

    scalar oldRatioMax = ratioMax_;
    scalar oldRatioMin = ratioMin_;
    scalar oldGrowthFactor = growthFactor_;

    ratioMax_ = readScalar(edgeOptionDict.lookup("bisectionRatio"));
    ratioMin_ = readScalar(edgeOptionDict.lookup("collapseRatio"));

    if (reRead)
    {
        // Check if values have changed, and report it.
        if (mag(oldRatioMax - ratioMax_) > SMALL)
        {
            Info << "\tOld ratioMax: " << oldRatioMax << nl
                 << "\tNew ratioMax: " << ratioMax_ << endl;
        }

        if (mag(oldRatioMin - ratioMin_) > SMALL)
        {
            Info << "\tOld ratioMin: " << oldRatioMin << nl
                 << "\tNew ratioMin: " << ratioMin_ << endl;
        }

        if (mag(oldGrowthFactor - growthFactor_) > SMALL)
        {
            Info << "\tOld growthFactor: " << oldGrowthFactor << nl
                 << "\tNew growthFactor: " << growthFactor_ << endl;
        }
    }

    if (edgeOptionDict.found("maxLengthScale") || mandatory_)
    {
        maxLengthScale_ =
        (
            readScalar(edgeOptionDict.lookup("maxLengthScale"))
        );
    }

    if (edgeOptionDict.found("minLengthScale") || mandatory_)
    {
        minLengthScale_ =
        (
            readScalar(edgeOptionDict.lookup("minLengthScale"))
        );
    }

    // Sanity check: Are length scales correctly specified?
    if (minLengthScale_ > maxLengthScale_)
    {
        FatalErrorIn("dynamicTopoFvMesh::readRefinementOptions()")
            << " Length-scales are incorrectly specified." << nl
            << " minLengthScale: " << minLengthScale_ << nl
            << " maxLengthScale: " << maxLengthScale_ << nl
            << abort(FatalError);
    }

    if (edgeOptionDict.found("fixedLengthScalePatches") || mandatory_)
    {
        fixedPatches_ =
        (
            edgeOptionDict.subDict("fixedLengthScalePatches")
        );

        // Ensure that patches are legitimate.
        checkPatches(fixedPatches_.toc());
    }

    // Check if swapping is to be avoided on any patches
    if (edgeOptionDict.found("noModificationPatches") || mandatory_)
    {
        wordList noModPatches =
        (
            edgeOptionDict.subDict("noModificationPatches").toc()
        );

        // Ensure that patches are legitimate.
        checkPatches(noModPatches);

        noModPatchIDs_.setSize(noModPatches.size());

        label indexI = 0;

        forAll(noModPatches, wordI)
        {
            const word& patchName = noModPatches[wordI];

            noModPatchIDs_[indexI++] =
            (
                boundaryMesh().findPatchID(patchName)
            );
        }
    }

    // Check local coupled patches for fixed length-scales
    if (dict_.found("coupledPatches") || mandatory_)
    {
        const dictionary coupledPatches = dict_.subDict("coupledPatches");

        // Determine master and slave patches
        label indexI = 0;
        wordList masterPatches(coupledPatches.size());
        wordList slavePatches(coupledPatches.size());

        forAllConstIter(dictionary, coupledPatches, dIter)
        {
            const dictionary& dictI = dIter().dict();

            masterPatches[indexI] = word(dictI.lookup("master"));
            slavePatches[indexI] = word(dictI.lookup("slave"));

            indexI++;
        }

        // Ensure that patches are legitimate.
        checkPatches(masterPatches);

        // Check whether coupled patches are fixedPatches as well.
        forAll(masterPatches, wordI)
        {
            word pName(masterPatches[wordI]);

            if (fixedPatches_.found(pName))
            {
                // Add the slave patch to the list as well.
                // If it already exists, over-ride the value.
                fixedPatches_.add
                (
                    slavePatches[wordI],
                    fixedPatches_[pName][0].scalarToken(),
                    true
                );
            }
        }
    }

    if (edgeOptionDict.found("freeLengthScalePatches") || mandatory_)
    {
        freePatches_ =
        (
            edgeOptionDict.subDict("freeLengthScalePatches")
        );

        // Ensure that patches are legitimate.
        checkPatches(freePatches_.toc());

        // Check if fixed and free patches are conflicting
        if (fixedPatches_.size() && freePatches_.size())
        {
            wordList fixedPatchList = fixedPatches_.toc();
            wordList freePatchList = freePatches_.toc();

            forAll(fixedPatchList, wordI)
            {
                forAll(freePatchList, wordJ)
                {
                    if (fixedPatchList[wordI] == freePatchList[wordJ])
                    {
                        FatalErrorIn
                        (
                            "dynamicTopoFvMesh::readRefinementOptions()"
                        )
                            << " Conflicting fixed/free patches." << nl
                            << " Fixed patch: " << fixedPatchList[wordI] << nl
                            << " Free patch: " << freePatchList[wordJ] << nl
                            << abort(FatalError);
                    }
                }
            }
        }
    }

    if (edgeOptionDict.found("computeGrowthFactor") || mandatory_)
    {
        if (!reRead)
        {
            // Compute the growth factor from the mesh for the first time.
            growthFactor_ = computeGrowthFactor();
        }
    }
    else
    {
        growthFactor_ = readScalar(edgeOptionDict.lookup("growthFactor"));
    }

    if (edgeOptionDict.found("curvaturePatches") || mandatory_)
    {
        curvaturePatches_ =
        (
            edgeOptionDict.subDict("curvaturePatches")
        );

        // Ensure that patches are legitimate.
        checkPatches(curvaturePatches_.toc());

        curvatureDeviation_ =
        (
            readScalar(edgeOptionDict.lookup("curvatureDeviation"))
        );

        if
        (
            (curvatureDeviation_ > 1.0 || curvatureDeviation_ < 0.0)
        )
        {
            FatalErrorIn("dynamicTopoFvMesh::readRefinementOptions()")
                << " Curvature deviation out of range [0..1]"
                << abort(FatalError);
        }
    }

    if (edgeOptionDict.found("proximityPatches") || mandatory_)
    {
        proximityPatches_ =
        (
            edgeOptionDict.subDict("proximityPatches")
        );

        // Ensure that patches are legitimate.
        checkPatches(proximityPatches_.toc());

        // Check if a threshold for slicing has been specified.
        if (edgeOptionDict.found("sliceThreshold") || mandatory_)
        {
            sliceThreshold_ =
            (
                readScalar(edgeOptionDict.lookup("sliceThreshold"))
            );

            // Cap the threshold value
            sliceThreshold_ = Foam::max(sliceThreshold_, minLengthScale_);
        }
    }

    if (edgeOptionDict.found("sliverThreshold") || mandatory_)
    {
        sliverThreshold_ =
        (
            readScalar(edgeOptionDict.lookup("sliverThreshold"))
        );

        if
        (
            (sliverThreshold_ > 1.0 || sliverThreshold_ < 0.0)
        )
        {
            FatalErrorIn("dynamicTopoFvMesh::readRefinementOptions()")
                << " Sliver threshold out of range [0..1]"
                << abort(FatalError);
        }
    }

    if (edgeOptionDict.found("maxModifications") || mandatory_)
    {
        maxModifications_ =
        (
            readLabel(edgeOptionDict.lookup("maxModifications"))
        );
    }
}

// Check for legitimacy of patches
void dynamicTopoFvMesh::checkPatches
(
    const wordList& patchList
) const
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    forAll(patchList, wordI)
    {
        bool foundPatch = false;

        forAll(boundary, patchI)
        {
            if (boundary[patchI].name() == patchList[wordI])
            {
                foundPatch = true;
                break;
            }
        }

        if (!foundPatch)
        {
            FatalErrorIn("dynamicTopoFvMesh::checkPatches()")
                << " Could not find patch: "
                << patchList[wordI] << nl
                << abort(FatalError);
        }
    }
}

// Given a boundary quad face, return a boundary triangular face.
// For 2D simplical meshes only.
label dynamicTopoFvMesh::getTriBoundaryFace
(
    const label fIndex
) const
{
    const labelList& fEdges = faceEdges_[fIndex];

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
        << fIndex << " :: " << faces_[fIndex]
        << abort(FatalError);

    return -1;
}

// Given a boundary quad face, pick out a boundary edge that
// contains a triangular face. For 2D simplical meshes only.
label dynamicTopoFvMesh::getTriBoundaryEdge
(
    const label fIndex
) const
{
    const labelList& fEdges = faceEdges_[fIndex];

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
        << fIndex << " :: " << faces_[fIndex]
        << abort(FatalError);

    return -1;
}

// 2D Edge-swapping engine
void dynamicTopoFvMesh::swap2DEdges(void *argument)
{
    // Recast the argument
    threadHandler<dynamicTopoFvMesh> *thread =
    (
        reinterpret_cast<threadHandler<dynamicTopoFvMesh>*>(argument)
    );

    if (thread->slave())
    {
        thread->sendSignal(threadHandler<dynamicTopoFvMesh>::START);
    }

    dynamicTopoFvMesh& mesh = thread->reference();

    // Figure out which thread this is...
    label tIndex = mesh.self();

    // Pick items off the stack
    while (!mesh.faceStack(tIndex).empty())
    {
        // Retrieve the index for this face
        label fIndex = mesh.faceStack(tIndex).pop();

        // Perform a Delaunay test and check if a flip is necesary.
        bool failed = mesh.testDelaunay(fIndex);

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

    if (thread->slave())
    {
        thread->sendSignal(threadHandler<dynamicTopoFvMesh>::STOP);
    }
}

// Initialize edge related connectivity lists
void dynamicTopoFvMesh::initEdges()
{
    // Initialize eMesh, and copy to local lists
    eMeshPtr_.set(new eMesh(*this));

    // Obtain information
    nEdges_ = eMeshPtr_->nEdges();
    nInternalEdges_ = eMeshPtr_->nInternalEdges();
    edgePatchSizes_ = eMeshPtr_->boundary().patchSizes();
    edgePatchStarts_ = eMeshPtr_->boundary().patchStarts();

    // Set old edge information
    nOldEdges_ = nEdges_;
    nOldInternalEdges_ = nInternalEdges_;
    oldEdgePatchSizes_ = edgePatchSizes_;
    oldEdgePatchStarts_ = edgePatchStarts_;

    // Set local lists with edge connectivity information
    edges_ = eMeshPtr_->edges();
    edgeFaces_ = eMeshPtr_->edgeFaces();
    faceEdges_ = eMeshPtr_->faceEdges();

    if (!twoDMesh_)
    {
        pointEdges_ = eMeshPtr_->pointEdges();
        edgePoints_ = eMeshPtr_->edgePoints();
    }

    // Read coupled patch information from dictionary.
    readCoupledPatches();
}

// Load the mesh-quality metric library
void dynamicTopoFvMesh::loadMetricLibrary()
{
    if (twoDMesh_)
    {
        return;
    }

    void * metricLibPtr = NULL;
    char * error;

    if
    (
        dict_.subDict("dynamicTopoFvMesh").found("tetMetricLib") ||
        mandatory_
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
            "dynamicTopoFvMesh::loadMetricLibrary() "
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
        typedef void (*returnType) ();

        // Load the list of symbols
        returnType availableMetrics =
        (
            reinterpret_cast<returnType>
            (
                dlsym(metricLibPtr,"reportMetrics")
            )
        );

        Info << " Available metrics: " << endl;

        // Execute the reportMetrics function.
        (*availableMetrics)();

        FatalErrorIn
        (
            "dynamicTopoFvMesh::loadMetricLibrary() "
        ) << nl << " Unrecognized tet-quality metric: " << tetMetric
          << abort(FatalError);
    }
}

// Initialize the threading environment.
//  - Provides an override option to avoid reading from the dictionary.
void dynamicTopoFvMesh::initializeThreadingEnvironment
(
    const label specThreads
)
{
    // Initialize an IOobject for the IOmultiThreader object
    IOobject io
    (
        "threader",
        this->time().timeName(),
        (*this),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    );

    if (specThreads > 0)
    {
        threader_.set(new IOmultiThreader(io, specThreads));
    }
    else
    {
        if (dict_.subDict("dynamicTopoFvMesh").found("threads") || mandatory_)
        {
            threader_.set
            (
                new IOmultiThreader
                (
                    io,
                    readLabel
                    (
                        dict_.subDict("dynamicTopoFvMesh").lookup("threads")
                    )
                )
            );
        }
        else
        {
            threader_.set(new IOmultiThreader(io, 1));
        }
    }

    // Get the number of threads and allocate threadHandlers
    label nThreads = threader_->getNumThreads();

    if (nThreads == 1)
    {
        handlerPtr_.setSize(1);

        handlerPtr_.set
        (
            0,
            new threadHandler<dynamicTopoFvMesh>
            (
                (*this),
                threader()
            )
        );

        handlerPtr_[0].setMaster();

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
        handlerPtr_.setSize(nThreads + 1);

        // Set the thread scheduling sequence
        topoSeq_.setSize(nThreads);

        // Linear sequence from 1 to nThreads
        forAll(topoSeq_, indexI)
        {
            topoSeq_[indexI] = indexI + 1;
        }

        // Size the stacks
        if (twoDMesh_)
        {
            faceStack_.setSize(nThreads + 1);
        }
        else
        {
            edgeStack_.setSize(nThreads + 1);
        }

        forAll(handlerPtr_, threadI)
        {
            handlerPtr_.set
            (
                threadI,
                new threadHandler<dynamicTopoFvMesh>
                (
                    (*this),
                    threader()
                )
            );

            if (threadI == 0)
            {
                handlerPtr_[0].setID(-1);
                handlerPtr_[0].setMaster();
            }
            else
            {
                handlerPtr_[threadI].setID(threader_->getID(threadI-1));
                handlerPtr_[threadI].setSlave();
            }
        }

        // For reOrdering, one handler for each reOrdering method
        reOrderPtr_.setSize(4);

        // Initialize reOrdering handlers
        forAll(reOrderPtr_, memberI)
        {
            reOrderPtr_.set
            (
                memberI,
                new threadHandler<dynamicTopoFvMesh>
                (
                    (*this),
                    threader()
                )
            );
        }

        // Set the thread scheduling sequence
        reOrderSeq_.setSize(4);

        // Points, cells, faces and edges (in that order)
        reOrderSeq_[0] = 0;
        reOrderSeq_[1] = 3;
        reOrderSeq_[2] = 2;
        reOrderSeq_[3] = 1;

        // Set argument sizes for individual members

        // Points take two arguments
        // (One pointField and one labelListList)
        reOrderPtr_[0].setSize(2);

        // Edges take three arguments
        // (One edgeList and two labelListLists)
        reOrderPtr_[1].setSize(3);

        // Faces take five arguments
        // (One faceList, two labelLists, and two labelListLists)
        reOrderPtr_[2].setSize(5);

        // Cells take one argument
        // (One labelListList)
        reOrderPtr_[3].setSize(1);
    }
}

// Return a reference to the multiThreader
const multiThreader& dynamicTopoFvMesh::threader() const
{
    return threader_();
}

// Does the mesh perform edge refinement?
bool dynamicTopoFvMesh::edgeRefinement() const
{
    return edgeRefinement_;
}

// Build the global shared-points list
void dynamicTopoFvMesh::buildGlobalPoints()
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    // Local variables
    label nSendPoints = 0;
    DynamicList<label> sharedPoints(boundary[0].size());
    labelList startBuffer(0), procBuffer(0);

    forAll(boundary, patchI)
    {
        if (isA<processorPolyPatch>(boundary[patchI]))
        {
            // Build a list of points and send them to the master.
            label start = patchStarts_[patchI];

            for (label i = 0; i < patchSizes_[patchI]; i++)
            {
                const face& faceToCheck = faces_[i + start];

                forAll(faceToCheck, pointI)
                {
                    if (findIndex(sharedPoints, faceToCheck[pointI]) == -1)
                    {
                        sharedPoints.append(faceToCheck[pointI]);
                    }
                }
            }
        }
    }

    // Prepare send sizes.
    nSendPoints = sharedPoints.size();

    if (Pstream::master())
    {
        labelList nRecvPoints(Pstream::nProcs(), 0);
        List<pointField> pBuffer(Pstream::nProcs());

        // Loop through all slave processors and receive info.
        for (label proc = 1; proc < Pstream::nProcs(); proc++)
        {
            // How many entities am I going to be receiving?
            pRead(proc, nRecvPoints[proc]);

            if (nRecvPoints[proc])
            {
                // Size the buffer
                pBuffer[proc].setSize(nRecvPoints[proc], vector::zero);

                // Schedule for receipt.
                pRead(proc, pBuffer[proc]);
            }
        }

        // As a master, copy my information to the receive buffer
        // for spatial comparisons.
        nRecvPoints[0] = nSendPoints;

        if (nRecvPoints[0])
        {
            // Size the buffer
            pBuffer[0].setSize(nRecvPoints[0], vector::zero);

            forAll(sharedPoints, pointI)
            {
                pBuffer[0][pointI] = points_[sharedPoints[pointI]];
            }
        }

        // Wait for all transfers to complete.
        waitForBuffers();

        // Global shared point detection algorithm
        //  - Performs a geometric comparison of positions.
        //  - Since this has the potential of being quite a large number,
        //    perform spatial-binning to minimize the search.
        //  - As a first-pass, determine the bounding-box.
        //  - Sub-divide the bounding box into regular intervals,
        //    and bin positions based on computed intervals.
        //  - Now perform a linear search between points in each bin.

        // Find a bound-box start point.
        bool foundStart = false;
        vector minLoc = vector::zero;
        vector maxLoc = vector::zero;

        for (label proc = 0; proc < Pstream::nProcs(); proc++)
        {
            if (nRecvPoints[proc])
            {
                // Assign the centre.
                minLoc = pBuffer[proc][0];
                maxLoc = pBuffer[proc][0];

                foundStart = true;
                break;
            }

            if (foundStart)
            {
                break;
            }
        }

        label index = 0;
        labelList pStarts(Pstream::nProcs(), 0);
        pointField points(sum(nRecvPoints), vector::zero);
        labelList pIndices(sum(nRecvPoints), 0);

        for (label proc = 0; proc < Pstream::nProcs(); proc++)
        {
            if (nRecvPoints[proc])
            {
                for(label pointI = 0; pointI < nRecvPoints[proc]; pointI++)
                {
                    minLoc = ::Foam::min(minLoc, pBuffer[proc][pointI]);
                    maxLoc = ::Foam::max(maxLoc, pBuffer[proc][pointI]);

                    points[index] = pBuffer[proc][pointI];
                    pIndices[index] = index;

                    index++;
                }
            }

            // Set point-starts to determine processor IDs
            if (proc > 0)
            {
                pStarts[proc] = pStarts[proc-1] + nRecvPoints[proc-1];
            }
        }

        // Write out received buffers for debugging
        if (debug > 2)
        {
            Info << "nRecvPoints: " << nRecvPoints << endl;
            Info << "Bounding box max: " << maxLoc << endl;
            Info << "Bounding box min: " << minLoc << endl;
        }

        labelListList pointBins(997, labelList(0));

        // Prepare a boundBox for spatial hashing
        boundBox box(minLoc, maxLoc);

        // Perform a spatial hash of all point locations
        spatialHash(points, pIndices, box, 10, pointBins);

        // Mapping between points and processors.
        Map<labelList> pMap;

        // Now that points have been grouped together, perform comparisons
        // between points in the same bin.
        forAll(pointBins, binI)
        {
            const labelList& bin = pointBins[binI];

            forAll(bin, pointI)
            {
                bool foundAMatch = false;

                forAll(bin, pointJ)
                {
                    if
                    (
                        mag(points[bin[pointI]] - points[bin[pointJ]]) < gTol_
                     && (bin[pointI] != bin[pointJ])
                    )
                    {
                        // Found a match. Which processor does each
                        // point belong to ?
                        FixedList<label, 2> proc(-1), loc(-1);

                        forAll(pStarts, i)
                        {
                            if
                            (
                                bin[pointI] >= pStarts[i]
                             && bin[pointI] < pStarts[i] + nRecvPoints[i]
                            )
                            {
                                proc[0] = i;
                                loc[0] = bin[pointI];
                            }

                            if
                            (
                                bin[pointJ] >= pStarts[i]
                             && bin[pointJ] < pStarts[i] + nRecvPoints[i]
                            )
                            {
                                proc[1] = i;
                                loc[1] = bin[pointJ];
                            }

                            if (proc[0] > -1 && proc[1] > -1)
                            {
                                break;
                            }
                        }

                        // Make an entry stating that the specified
                        // point talks to the corresponding processor.
                        forAll(loc, i)
                        {
                            if (pMap.found(loc[i]))
                            {
                                forAll(proc, j)
                                {
                                    if (findIndex(pMap[loc[i]], proc[j]) == -1)
                                    {
                                        // Size up the list.
                                        sizeUpList(proc[j], pMap[loc[i]]);
                                    }
                                }
                            }
                            else
                            {
                                // Copy proc contents.
                                pMap.insert(loc[i], proc);
                            }
                        }

                        foundAMatch = true;
                    }
                }

                if (!foundAMatch)
                {
                    // Should've found at least one match.
                    FatalErrorIn
                    (
                        "dynamicTopoFvMesh::identifyCoupledPatches()"
                    )
                        << "Could not match point: "
                        << points[bin[pointI]]
                        << abort(FatalError);
                }
            }
        }

        // Send assembled processor sharing info back to slaves.

        // Size-up buffers first.
        labelListList labelStarts(Pstream::nProcs());
        labelListList procInfo(Pstream::nProcs());

        for (label proc = 0; proc < Pstream::nProcs(); proc++)
        {
            if (nRecvPoints[proc])
            {
                labelStarts[proc].setSize(nRecvPoints[proc], 0);

                // Count the number of entities that we will be sending.
                label j = 0, totalSize = 0, start = pStarts[proc];

                for (label i = 0; i < nRecvPoints[proc]; i++)
                {
                    totalSize += pMap[i+start].size();
                }

                procInfo[proc].setSize(totalSize, 0);

                // Now fill the buffer.
                for (label i = 0; i < nRecvPoints[proc]; i++)
                {
                    const labelList& procList = pMap[i+start];

                    labelStarts[proc][i] = j;

                    forAll(procList, edgeI)
                    {
                        procInfo[proc][j++] = procList[edgeI];
                    }
                }

                if (proc == 0)
                {
                    // As a master, copy my own processor
                    // sharing info for later analysis.
                    startBuffer.setSize(nSendPoints, -1);
                    procBuffer.setSize(totalSize, -1);

                    startBuffer = labelStarts[0];
                    procBuffer = procInfo[0];
                }
                else
                {
                    // Send back to the slave.
                    pWrite(proc, totalSize);
                    pWrite(proc, labelStarts[proc]);
                    pWrite(proc, procInfo[proc]);
                }
            }
        }

        // Wait for all transfers to complete.
        waitForBuffers();
    }
    else
    {
        // Send number of entities first.
        pWrite(Pstream::masterNo(), nSendPoints);

        if (debug > 2)
        {
            Pout << "Sending " << nSendPoints << " points. " << endl;
        }

        if (nSendPoints)
        {
            // Size the buffer and copy info into it.
            pointField pBuffer(nSendPoints, vector::zero);

            forAll(sharedPoints, pointI)
            {
                pBuffer[pointI] = points_[sharedPoints[pointI]];
            }

            // Send buffers to the master.
            pWrite(Pstream::masterNo(), pBuffer);

            // Receive processor sharing info from the master.
            label procInfoSize = -1;

            // How many entities am I going to be receiving?
            pRead(Pstream::masterNo(), procInfoSize);

            // Size the receive buffer.
            startBuffer.setSize(nSendPoints, -1);
            procBuffer.setSize(procInfoSize, -1);

            // Schedule for receipt.
            pRead(Pstream::masterNo(), startBuffer);
            pRead(Pstream::masterNo(), procBuffer);

            // Wait for all transfers to complete.
            waitForBuffers();
        }
    }

    procIndices_.clear();
    subMeshPoints_.clear();

    Map<labelHashSet> procPatchPoints;

    forAll(startBuffer, pointI)
    {
        label start = startBuffer[pointI];

        label end =
        (
            (pointI == (startBuffer.size() - 1)) ?
            procBuffer.size() : startBuffer[pointI + 1]
        );

        for (label i = start; i < end; i++)
        {
            if (procBuffer[i] != Pstream::myProcNo())
            {
                // Was this processor detected before?
                if (!procPatchPoints.found(procBuffer[i]))
                {
                    // Add a new entry.
                    procPatchPoints.insert(procBuffer[i], labelHashSet(10));
                }

                // Add this point.
                procPatchPoints[procBuffer[i]].insert(sharedPoints[pointI]);
            }
        }
    }

    if (procPatchPoints.size())
    {
        // Copy the indices.
        procIndices_ = procPatchPoints.toc();

        // Patch sub meshes need to be prepared in ascending
        // order of neighbouring processors.
        sort(procIndices_);

        subMeshPoints_.setSize(procIndices_.size());

        // Copy the list of points for each processor.
        forAll(procIndices_, patchI)
        {
            subMeshPoints_[patchI].setSize
            (
                procPatchPoints[procIndices_[patchI]].size(), -1
            );

            subMeshPoints_[patchI] =
            (
                procPatchPoints[procIndices_[patchI]].toc()
            );
        }
    }

    if (debug > 3)
    {
        forAll(procIndices_, patchI)
        {
            // Write out points as a VTK
            writeVTK
            (
                "subMeshPoints_" +
                Foam::name(Pstream::myProcNo()) + "to" +
                Foam::name(procIndices_[patchI]),
                subMeshPoints_[patchI],
                0
            );
        }
    }
}

// Identify coupled patches.
//  - Also builds global shared point information.
//  - Returns true if no coupled patches were found.
bool dynamicTopoFvMesh::identifyCoupledPatches()
{
    bool coupledPatchesAbsent = true;

    procIndices_.clear();

    // Check if patches are explicitly coupled
    if (patchCoupling_.size())
    {
        coupledPatchesAbsent = false;
    }

    // Maintain a separate list of processor IDs.
    // This is done because this sub-domain may talk to processors
    // that share only edges/points.
    if (Pstream::parRun())
    {
        const polyBoundaryMesh& boundary = boundaryMesh();

        forAll(boundary, patchI)
        {
            if (isA<processorPolyPatch>(boundary[patchI]))
            {
                coupledPatchesAbsent = false;

                break;
            }
        }

        // Build the global shared-points information, if necessary.
        buildGlobalPoints();
    }

    return coupledPatchesAbsent;
}

// Read coupled patch information from dictionary.
void dynamicTopoFvMesh::readCoupledPatches()
{
    patchCoupling_.clear();

    if (dict_.found("coupledPatches") || mandatory_)
    {
        const dictionary& coupledPatches = dict_.subDict("coupledPatches");

        const polyBoundaryMesh& boundary = boundaryMesh();

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
                // Check whether patches are associated with zones.
                Switch specifyZones
                (
                    dictI.lookup("specifyZones")
                );

                label mZone = -1, sZone = -1;

                if (specifyZones)
                {
                    const faceZoneMesh& faceZones = polyMesh::faceZones();

                    mZone = faceZones.findZoneID
                    (
                        dictI.lookup("masterZone")
                    );

                    sZone = faceZones.findZoneID
                    (
                        dictI.lookup("slaveZone")
                    );
                }

                patchCoupling_.insert
                (
                    mPatch,
                    coupledPatchInfo(sPatch, mZone, sZone)
                );
            }
            else
            {
                FatalErrorIn("dynamicTopoFvMesh::readCoupledPatches()")
                        << " Coupled patches are either wrongly specified,"
                        << " or the sizes don't match." << nl
                        << " Master: " << mPatch << ":" << masterPatch
                        << " Size: " << boundary[mPatch].size() << nl
                        << " Slave: " << sPatch << ":" << slavePatch
                        << " Size: " << boundary[sPatch].size() << nl
                        << abort(FatalError);
            }
        }
    }

    // Check if a geometric tolerance has been specified.
    if (dict_.found("gTol") || mandatory_)
    {
        gTol_ = readScalar(dict_.lookup("gTol"));
    }
    else
    {
        // If not, attempt to calculate it from the mesh.
        calculateGeometricTolerance();
    }

    // Initialize entitiesToAvoid to some arbitrary size
    entitiesToAvoid_.setSize(50, -1);
}

// Initialize coupled patches for topology modifications.
//  - Send and receive sub meshes for processor patches
void dynamicTopoFvMesh::initCoupledPatches()
{
    // Identify coupled patches.
    if (identifyCoupledPatches())
    {
        return;
    }

    // Build and send patch sub-meshes (and clear existing ones).
    buildCoupledPatchMeshes();

    // Build maps for locally coupled patches.
    buildLocalCoupledMaps();

    // Build maps for coupled processor patches.
    buildProcCoupledMaps();
}

// Handle topology changes for coupled patches
void dynamicTopoFvMesh::handleCoupledPatches()
{
    if (!(patchCoupling_.size() || procIndices_.size()))
    {
        return;
    }

    if (debug)
    {
        Info << "Handling coupled patches...";
    }

    // Set coupled modifications.
    setCoupledModification();
    unsetSlaveModification();

    // Loop through the coupled stack and perform changes.
    if (twoDMesh_)
    {
        if (edgeRefinement_)
        {
            // Initialize the face stack
            initCoupledFaceStack();

            edgeBisectCollapse2D(&(handlerPtr_[0]));
        }

        // Re-Initialize the face stack
        initCoupledFaceStack();

        swap2DEdges(&(handlerPtr_[0]));
    }
    else
    {
        if (edgeRefinement_)
        {
            // Initialize the edge stack
            initCoupledEdgeStack();

            edgeBisectCollapse3D(&(handlerPtr_[0]));
        }

        // Re-Initialize the edge stack
        initCoupledEdgeStack();

        swap3DEdges(&(handlerPtr_[0]));
    }

    // Build a list of entities that need to be avoided
    // by regular topo-changes.
    buildEntitiesToAvoid();

    // Reset coupled modifications.
    unsetCoupledModification();
    unsetSlaveModification();

    if (debug)
    {
        Info << "Done." << endl;
    }
}

// Build patch sub-meshes for processor patches
void dynamicTopoFvMesh::buildCoupledPatchMeshes()
{
    if (!procIndices_.size())
    {
        return;
    }

    if (debug)
    {
        Info << "Building patch sub-meshes for processor patches...";
    }

    // Size the PtrLists.
    sendPatchMeshes_.clear();
    recvPatchMeshes_.clear();

    sendPatchMeshes_.setSize(Pstream::nProcs());
    recvPatchMeshes_.setSize(Pstream::nProcs());

    // Lists of entities that need to be avoided while building sub meshes.
    labelHashSet cellsToAvoid;

    forAll(procIndices_, procI)
    {
        label proc = procIndices_[procI];

        if (proc < Pstream::myProcNo())
        {
            sendPatchMeshes_.set(proc, new coupledPatchInfo());

            coupledPatchInfo& sendMesh = sendPatchMeshes_[proc];

            // Build the subMesh.
            buildCoupledPatchMesh
            (
                proc,
                sendMesh,
                cellsToAvoid
            );

            // Send my sub-mesh to the neighbour.
            pWrite(proc, sendMesh.nEntities());

            if (debug > 3)
            {
                Pout << "Sending:" << nl
                     << "\t nEntities: " << sendMesh.nEntities()
                     << endl;
            }

            // Send the pointBuffer
            pWrite(proc, sendMesh.pointBuffer());

            // Send connectivity (points, edges, faces, cells, etc)
            forAll(sendMesh.entityBuffer(), bufferI)
            {
                pWrite(proc, sendMesh.entityBuffer(bufferI));
            }

            if (edgeRefinement_)
            {
                pWrite(proc, sendMesh.lengthBuffer());
            }
        }
        else
        {
            recvPatchMeshes_.set(proc, new coupledPatchInfo());

            coupledPatchInfo& recvMesh = recvPatchMeshes_[proc];

            // First read entity sizes.
            pRead(proc, recvMesh.nEntities());

            if (debug > 3)
            {
                Pout << "Receiving:" << nl
                        << "\t nEntities: " << recvMesh.nEntities()
                        << endl;
            }

            // Size the buffers.
            recvMesh.pointBuffer().setSize
            (
                recvMesh.nEntities(coupledPatchInfo::POINT)
            );

            recvMesh.entityBuffer(coupledPatchInfo::POINT).setSize
            (
                recvMesh.nEntities(coupledPatchInfo::SHARED_POINT)
            );

            recvMesh.entityBuffer(coupledPatchInfo::EDGE).setSize
            (
                2*recvMesh.nEntities(coupledPatchInfo::EDGE)
            );

            recvMesh.entityBuffer(coupledPatchInfo::FACE).setSize
            (
                3*recvMesh.nEntities(coupledPatchInfo::FACE)
            );

            recvMesh.entityBuffer(coupledPatchInfo::CELL).setSize
            (
                4*recvMesh.nEntities(coupledPatchInfo::CELL)
            );

            recvMesh.entityBuffer(coupledPatchInfo::FACE_EDGE).setSize
            (
                3*recvMesh.nEntities(coupledPatchInfo::FACE)
            );

            // Receive the pointBuffer
            pRead(proc, recvMesh.pointBuffer());

            // Receive connectivity (points, edges, faces, cells, etc)
            forAll(recvMesh.entityBuffer(), bufferI)
            {
                pRead(proc, recvMesh.entityBuffer(bufferI));
            }

            if (edgeRefinement_)
            {
                recvMesh.lengthBuffer().setSize
                (
                    recvMesh.nEntities(coupledPatchInfo::CELL)
                );

                pRead(proc, recvMesh.lengthBuffer());
            }
        }
    }

    // We won't wait for all transfers to complete for the moment.
    // Meanwhile, build local coupling maps (or do some other useful work)

    if (debug)
    {
        Info << "Done." << endl;
    }
}

// Build patch sub-mesh for a specified processor
void dynamicTopoFvMesh::buildCoupledPatchMesh
(
    const label proc,
    coupledPatchInfo& subMesh,
    labelHashSet& cellsToAvoid
)
{
    label nP = 0, nE = 0, nF = 0, nC = 0, sP = 0;

    // Shorten enumerants for convenience
    const label pointEnum = coupledPatchInfo::POINT;
    const label edgeEnum  = coupledPatchInfo::EDGE;
    const label faceEnum  = coupledPatchInfo::FACE;
    const label cellEnum  = coupledPatchInfo::CELL;
    const label iEdgeEnum = coupledPatchInfo::INTERNAL_EDGE;
    const label fEdgeEnum = coupledPatchInfo::FACE_EDGE;
    const label shPointEnum = coupledPatchInfo::SHARED_POINT;

    // Obtain references
    Map<label>& rPointMap = subMesh.reverseEntityMap(pointEnum);
    Map<label>& rEdgeMap  = subMesh.reverseEntityMap(edgeEnum);
    Map<label>& rFaceMap  = subMesh.reverseEntityMap(faceEnum);
    Map<label>& rCellMap  = subMesh.reverseEntityMap(cellEnum);

    Map<label>& pointMap = subMesh.entityMap(pointEnum);
    Map<label>& edgeMap  = subMesh.entityMap(edgeEnum);
    Map<label>& faceMap  = subMesh.entityMap(faceEnum);
    Map<label>& cellMap  = subMesh.entityMap(cellEnum);

    // Add all cells connected to points on the subMeshPoints list
    label procIndex = -1;

    forAll(procIndices_, procI)
    {
        if (proc == procIndices_[procI])
        {
            procIndex = procI;

            // Loop through points detected by identifyCoupledPatches
            forAll(subMeshPoints_[procIndex], pointI)
            {
                // Loop through pointEdges for this point.
                const labelList& pEdges =
                (
                    pointEdges_[subMeshPoints_[procIndex][pointI]]
                );

                forAll(pEdges, edgeI)
                {
                    const labelList& eFaces = edgeFaces_[pEdges[edgeI]];

                    forAll(eFaces, faceI)
                    {
                        label own = owner_[eFaces[faceI]];
                        label nei = neighbour_[eFaces[faceI]];

                        // Check owner cell
                        if
                        (
                           !rCellMap.found(own) &&
                           !cellsToAvoid.found(own)
                        )
                        {
                            cellMap.insert(nC, own);
                            rCellMap.insert(own, nC);
                            nC++;

                            cellsToAvoid.insert(own);
                        }

                        // Check neighbour cell
                        if
                        (
                           !rCellMap.found(nei) &&
                           !cellsToAvoid.found(nei) &&
                            nei != -1
                        )
                        {
                            cellMap.insert(nC, nei);
                            rCellMap.insert(nei, nC);
                            nC++;

                            cellsToAvoid.insert(nei);
                        }
                    }
                }
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

    // Allocate the edgeMap. Interior edges need to be detected
    // first and added before boundary ones. Do this in two stages.
    for (label stage = 0; stage < 2; stage++)
    {
        forAllIter(Map<label>::iterator, rFaceMap, fIter)
        {
            const labelList& fEdges = faceEdges_[fIter.key()];

            forAll(fEdges, edgeI)
            {
                if (!rEdgeMap.found(fEdges[edgeI]))
                {
                    if (stage == 0)
                    {
                        // Check if any of edgeFaces do not belong
                        // to the faceMap.
                        const labelList& eFaces = edgeFaces_[fEdges[edgeI]];

                        bool boundaryEdge = false;

                        forAll(eFaces, faceI)
                        {
                            if (!rFaceMap.found(eFaces[faceI]))
                            {
                                // This face was not found. Designate this
                                // edge as a 'boundary'.
                                boundaryEdge = true;

                                break;
                            }
                        }

                        if (!boundaryEdge)
                        {
                            edgeMap.insert(nE, fEdges[edgeI]);
                            rEdgeMap.insert(fEdges[edgeI], nE);
                            nE++;
                        }
                    }
                    else
                    {
                        // Adding only boundary edges at this stage.
                        edgeMap.insert(nE, fEdges[edgeI]);
                        rEdgeMap.insert(fEdges[edgeI], nE);
                        nE++;
                    }
                }
            }
        }

        // Set the number of internal edges at this point
        if (stage == 0)
        {
            subMesh.nEntities(iEdgeEnum) = nE;
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

    // Loop through subMeshPoints for the processor
    // and fill a mapped buffer for them as well.
    // This allows the neighbour to match-up edges easily.
    labelList& spBuffer = subMesh.entityBuffer(pointEnum);

    spBuffer.setSize(subMeshPoints_[procIndex].size(), -1);

    forAll(subMeshPoints_[procIndex], pointI)
    {
        if (rPointMap.found(subMeshPoints_[procIndex][pointI]))
        {
            spBuffer[sP++] = rPointMap[subMeshPoints_[procIndex][pointI]];
        }
    }

    // Shorten the buffer to actual size.
    spBuffer.setSize(sP);

    // Assign sizes to the mesh
    subMesh.nEntities(pointEnum) = nP;
    subMesh.nEntities(edgeEnum)  = nE;
    subMesh.nEntities(faceEnum)  = nF;
    subMesh.nEntities(cellEnum)  = nC;
    subMesh.nEntities(shPointEnum) = sP;

    // Size up buffers and fill them
    pointField& pBuffer = subMesh.pointBuffer();

    pBuffer.setSize(nP, vector::zero);

    forAllIter(Map<label>::iterator, pointMap, pIter)
    {
        pBuffer[pIter.key()] = points_[pIter()];
    }

    // Edge buffer size: 2 points for every edge
    labelList& eBuffer = subMesh.entityBuffer(edgeEnum);

    eBuffer.setSize(2 * nE, -1);

    label index = 0;

    forAllIter(Map<label>::iterator, edgeMap, eIter)
    {
        edge& edgeToCheck = edges_[eIter()];

        eBuffer[index++] = rPointMap[edgeToCheck[0]];
        eBuffer[index++] = rPointMap[edgeToCheck[1]];
    }

    labelList& fBuffer  = subMesh.entityBuffer(faceEnum);
    labelList& cBuffer  = subMesh.entityBuffer(cellEnum);
    labelList& feBuffer = subMesh.entityBuffer(fEdgeEnum);

    if (!twoDMesh_)
    {
        // Face buffer size: 3 points/edges for every face in 3D
        fBuffer.setSize(3 * nF, -1);
        feBuffer.setSize(3 * nF, -1);

        index = 0;

        face thisFace(3);

        forAllIter(Map<label>::iterator, faceMap, fIter)
        {
            if (!rCellMap.found(owner_[fIter()]))
            {
                // This face is pointed the wrong way.
                thisFace = faces_[fIter()].reverseFace();
            }
            else
            {
                thisFace = faces_[fIter()];
            }

            const labelList& fEdges = faceEdges_[fIter()];

            fBuffer[index] = rPointMap[thisFace[0]];
            feBuffer[index++] = rEdgeMap[fEdges[0]];

            fBuffer[index] = rPointMap[thisFace[1]];
            feBuffer[index++] = rEdgeMap[fEdges[1]];

            fBuffer[index] = rPointMap[thisFace[2]];
            feBuffer[index++] = rEdgeMap[fEdges[2]];
        }

        // Cell buffer size: 4 faces for every cell in 3D
        cBuffer.setSize(4 * nC, -1);

        index = 0;

        forAllIter(Map<label>::iterator, cellMap, cIter)
        {
            const cell& cellToCheck = cells_[cIter()];

            cBuffer[index++] = rFaceMap[cellToCheck[0]];
            cBuffer[index++] = rFaceMap[cellToCheck[1]];
            cBuffer[index++] = rFaceMap[cellToCheck[2]];
            cBuffer[index++] = rFaceMap[cellToCheck[3]];
        }
    }

    if (edgeRefinement_)
    {
        scalarList& lBuffer = subMesh.lengthBuffer();

        lBuffer.setSize(nC, 0.0);

        forAllIter(Map<label>::iterator, cellMap, cIter)
        {
            lBuffer[cIter.key()] = lengthScale_[cIter()];
        }
    }

    // For debugging purposes...
    if (debug > 3)
    {
        Pout << "Writing out coupledPatchInfo for processor: "
             << proc << endl;

        writeVTK
        (
            "psMesh_" + Foam::name(Pstream::myProcNo())
          + "to" + Foam::name(proc),
            rCellMap.toc()
        );
    }
}

// Build coupled maps
void dynamicTopoFvMesh::buildLocalCoupledMaps()
{
    if (!patchCoupling_.size())
    {
        return;
    }

    if (debug)
    {
        Info << "Building local coupled maps...";
    }

    const polyBoundaryMesh& boundary = boundaryMesh();

    DynamicList<label> mList(50), sList(50);

    forAllIter(Map<coupledPatchInfo>, patchCoupling_, patchI)
    {
        // Clear existing maps.
        patchI().clearMaps();

        // Build a list of master entities [faces in 2D / edges in 3D]
        label start = -1;

        if (twoDMesh_)
        {
            // Build a list of master entities
            start = boundary[patchI.key()].start();

            for (label i = 0; i < boundary[patchI.key()].size(); i++)
            {
                mList.append(start+i);
            }

            // Build a list of slave entities
            start = boundary[patchI().slaveIndex()].start();

            for (label i = 0; i < boundary[patchI().slaveIndex()].size(); i++)
            {
                sList.append(start+i);
            }
        }
        else
        {
            // Build a list of master entities
            start = boundary[patchI.key()].start();

            for (label i = 0; i < boundary[patchI.key()].size(); i++)
            {
                const labelList& mfEdges = faceEdges_[start+i];

                forAll(mfEdges, edgeI)
                {
                    if (findIndex(mList, mfEdges[edgeI]) == -1)
                    {
                        mList.append(mfEdges[edgeI]);
                    }
                }
            }

            // Build a list of slave entities
            start = boundary[patchI().slaveIndex()].start();

            for (label i = 0; i < boundary[patchI().slaveIndex()].size(); i++)
            {
                const labelList& sfEdges = faceEdges_[start+i];

                forAll(sfEdges, edgeI)
                {
                    if (findIndex(sList, sfEdges[edgeI]) == -1)
                    {
                        sList.append(sfEdges[edgeI]);
                    }
                }
            }
        }

        // Sanity check: Do patches have equal number of entities?
        if (mList.size() != sList.size())
        {
            FatalErrorIn("dynamicTopoFvMesh::buildLocalCoupledMaps()")
                << "Patch sizes are not consistent."
                << abort(FatalError);
        }

        // Build a list of entity-centres for geometric matching.
        pointField mCentres(mList.size(), vector::zero);
        pointField sCentres(sList.size(), vector::zero);

        if (twoDMesh_)
        {
            forAll(mList, faceI)
            {
                const face& mfaceToCheck = faces_[mList[faceI]];
                const face& sfaceToCheck = faces_[sList[faceI]];

                // Assume quad-faces for 2D patches.
                mCentres[faceI] = 0.25 *
                (
                    points_[mfaceToCheck[0]]
                  + points_[mfaceToCheck[1]]
                  + points_[mfaceToCheck[2]]
                  + points_[mfaceToCheck[3]]
                );

                sCentres[faceI] = 0.25 *
                (
                    points_[sfaceToCheck[0]]
                  + points_[sfaceToCheck[1]]
                  + points_[sfaceToCheck[2]]
                  + points_[sfaceToCheck[3]]
                );
            }
        }
        else
        {
            forAll(mList, edgeI)
            {
                const edge& mE = edges_[mList[edgeI]];
                const edge& sE = edges_[sList[edgeI]];

                mCentres[edgeI] = 0.5*(points_[mE[0]] + points_[mE[1]]);
                sCentres[edgeI] = 0.5*(points_[sE[0]] + points_[sE[1]]);
            }
        }

        label nMatchedEdges = 0;

        forAll(mCentres, indexI)
        {
            bool matched = false;
            scalar minDistance = GREAT;

            forAll(sCentres, indexJ)
            {
                scalar distance = mag(mCentres[indexI] - sCentres[indexJ]);

                minDistance = minDistance < distance
                            ? minDistance : distance;

                if (distance < gTol_)
                {
                    // Add a map entry
                    patchI().mapSlave
                    (
                        mList[indexI],
                        sList[indexJ]
                    );

                    patchI().mapMaster
                    (
                        sList[indexJ],
                        mList[indexI]
                    );

                    matched = true;

                    nMatchedEdges++;

                    break;
                }
            }

            if (!matched)
            {
                FatalErrorIn("dynamicTopoFvMesh::buildCoupledMaps()")
                    << " Failed to match edge: " << mList[indexI]
                    << ": " << edges_[mList[indexI]]
                    << " within a tolerance of: "
                    << gTol_ << nl << " Missed by: " << minDistance
                    << abort(FatalError);
            }
        }

        // Make sure we were successful.
        if (nMatchedEdges != mCentres.size())
        {
            FatalErrorIn("dynamicTopoFvMesh::buildLocalCoupledMaps()")
                << " Failed to match all patch edges within a tolerance of: "
                << gTol_ << nl
                << " Number of edges required for match: "
                << mCentres.size() << nl
                << " Number of matched edges: " << nMatchedEdges
                << abort(FatalError);
        }

        // Clear the lists
        mList.clear();
        sList.clear();
    }

    if (debug)
    {
        Info << "Done." << endl;
    }
}

// Build coupled maps for coupled processor patches
void dynamicTopoFvMesh::buildProcCoupledMaps()
{
    if (!procIndices_.size())
    {
        return;
    }

    // Wait for all transfers to complete.
    waitForBuffers();

    // Shorten enumerants for convenience
    const label edgeEnum  = coupledPatchInfo::EDGE;
    const label faceEnum  = coupledPatchInfo::FACE;
    const label cellEnum  = coupledPatchInfo::CELL;
    const label iEdgeEnum = coupledPatchInfo::INTERNAL_EDGE;
    const label fEdgeEnum = coupledPatchInfo::FACE_EDGE;

    forAll(procIndices_, procI)
    {
        label proc = procIndices_[procI];

        if (proc > Pstream::myProcNo())
        {
            coupledPatchInfo& recvMesh = recvPatchMeshes_[proc];

            pointField smPoints(recvMesh.pointBuffer());
            edgeList smEdges(recvMesh.nEntities(edgeEnum));
            faceList smFaces;
            cellList smCells;
            labelListList smFEdges;

            // Set connectivity from buffers.
            labelList& eBuffer = recvMesh.entityBuffer(edgeEnum);

            forAll(smEdges, edgeI)
            {
                smEdges[edgeI][0] = eBuffer[(2*edgeI) + 0];
                smEdges[edgeI][1] = eBuffer[(2*edgeI) + 1];
            }

            if (!twoDMesh_)
            {
                // Set sizes.
                smFaces.setSize(recvMesh.nEntities(faceEnum), face(3));
                smCells.setSize(recvMesh.nEntities(cellEnum), cell(4));
                smFEdges.setSize
                (
                    recvMesh.nEntities(fEdgeEnum),
                    labelList(3, -1)
                );

                // Copy connectivity from buffers.
                labelList& fBuffer = recvMesh.entityBuffer(faceEnum);
                labelList& feBuffer = recvMesh.entityBuffer(fEdgeEnum);

                forAll(smFaces, faceI)
                {
                    smFaces[faceI][0] = fBuffer[(3*faceI) + 0];
                    smFaces[faceI][1] = fBuffer[(3*faceI) + 1];
                    smFaces[faceI][2] = fBuffer[(3*faceI) + 2];

                    smFEdges[faceI][0] = feBuffer[(3*faceI) + 0];
                    smFEdges[faceI][1] = feBuffer[(3*faceI) + 1];
                    smFEdges[faceI][2] = feBuffer[(3*faceI) + 2];
                }

                labelList& cBuffer = recvMesh.entityBuffer(cellEnum);

                forAll(smCells, cellI)
                {
                    smCells[cellI][0] = cBuffer[(4*cellI) + 0];
                    smCells[cellI][1] = cBuffer[(4*cellI) + 1];
                    smCells[cellI][2] = cBuffer[(4*cellI) + 2];
                    smCells[cellI][3] = cBuffer[(4*cellI) + 3];
                }
            }

            // Set the autoPtr.
            recvMesh.setMesh
            (
                proc,
                new dynamicTopoFvMesh
                (
                    (*this),
                    IOobject
                    (
                        "subMesh",
                        time().timeName(),
                        time()
                    ),
                    smPoints,
                    recvMesh.nEntities(iEdgeEnum),
                    smEdges,
                    smFaces,
                    smFEdges,
                    smCells
                )
            );

            // If any edges here coincide with locally coupled edges,
            // processor edges need to be given priority.
        }
    }
}

// Wait for buffer transfer completion.
void dynamicTopoFvMesh::waitForBuffers() const
{
    if (Pstream::parRun())
    {
        OPstream::waitRequests();
        IPstream::waitRequests();
    }
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
    threadHandler<dynamicTopoFvMesh> *thread =
    (
        reinterpret_cast<threadHandler<dynamicTopoFvMesh>*>(argument)
    );

    if (thread->slave())
    {
        thread->sendSignal(threadHandler<dynamicTopoFvMesh>::START);
    }

    dynamicTopoFvMesh& mesh = thread->reference();

    // Figure out which thread this is...
    label tIndex = mesh.self();

    // Pick items off the stack
    while (!mesh.faceStack(tIndex).empty())
    {
        // Retrieve the index for this face
        label fIndex = mesh.faceStack(tIndex).pop();

        // Check if edgeRefinements are to be avoided.
        if (mesh.checkFaceModification(fIndex))
        {
            continue;
        }

        // Check if this boundary face is adjacent to a sliver-cell,
        // and remove it by a two-step bisection/collapse operation.
        mesh.remove2DSliver(fIndex);

        if (mesh.checkFaceBisection(fIndex))
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
        if (mesh.checkFaceCollapse(fIndex))
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

    if (thread->slave())
    {
        thread->sendSignal(threadHandler<dynamicTopoFvMesh>::STOP);
    }
}

// 3D Edge-swapping engine
void dynamicTopoFvMesh::swap3DEdges
(
    void *argument
)
{
    // Recast the argument
    threadHandler<dynamicTopoFvMesh> *thread =
    (
        reinterpret_cast<threadHandler<dynamicTopoFvMesh>*>(argument)
    );

    if (thread->slave())
    {
        thread->sendSignal(threadHandler<dynamicTopoFvMesh>::START);
    }

    dynamicTopoFvMesh& mesh = thread->reference();

    // Figure out which thread this is...
    label tIndex = mesh.self();

    // Dynamic programming variables
    labelList m;
    PtrList<scalarListList> Q;
    PtrList<labelListList> K, triangulations;

    // Allocate dynamic programming tables
    mesh.initTables(m, Q, K, triangulations);

    // Pick edges off the stack
    while (!mesh.edgeStack(tIndex).empty())
    {
        // Retrieve an edge from the stack
        label eIndex = mesh.edgeStack(tIndex).pop();

        // Compute the minimum quality of cells around this edge
        scalar minQuality = mesh.computeMinQuality(eIndex);

        // Check if this edge is on a bounding curve
        if (mesh.checkBoundingCurve(eIndex))
        {
            continue;
        }

        // Fill the dynamic programming tables
        if (mesh.fillTables(eIndex, minQuality, m, Q, K, triangulations))
        {
            // Check if edge-swapping is required.
            if (mesh.checkQuality(eIndex, m, Q, minQuality))
            {
                if (thread->master())
                {
                    // Remove this edge according to the swap sequence
                    mesh.removeEdgeFlips(eIndex, minQuality, K, triangulations);
                }
                else
                {
                    // Push this on to the master stack
                    mesh.edgeStack(0).push(eIndex);
                }
            }
        }
    }

    if (thread->slave())
    {
        thread->sendSignal(threadHandler<dynamicTopoFvMesh>::STOP);
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
    threadHandler<dynamicTopoFvMesh> *thread =
    (
        reinterpret_cast<threadHandler<dynamicTopoFvMesh>*>(argument)
    );

    if (thread->slave())
    {
        thread->sendSignal(threadHandler<dynamicTopoFvMesh>::START);
    }

    dynamicTopoFvMesh& mesh = thread->reference();

    // Figure out which thread this is...
    label tIndex = mesh.self();

    while (!mesh.edgeStack(tIndex).empty())
    {
        // Retrieve an edge from the stack
        label eIndex = mesh.edgeStack(tIndex).pop();

        // Check if edgeRefinements are to be avoided.
        if (mesh.checkEdgeModification(eIndex))
        {
            continue;
        }

        if (mesh.checkEdgeBisection(eIndex))
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
        if (mesh.checkEdgeCollapse(eIndex))
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

    if (thread->slave())
    {
        thread->sendSignal(threadHandler<dynamicTopoFvMesh>::STOP);
    }
}

// Utility method to check for invalid face-collapse.
// Returns 'true' if the collapse in NOT feasible.
bool dynamicTopoFvMesh::checkCollapse
(
    const labelList& triFaces,
    const FixedList<label,2>& c0BdyIndex,
    const FixedList<label,2>& c1BdyIndex,
    const FixedList<label,2>& original,
    const FixedList<label,2>& replacement,
    const bool checkNeighbour
) const
{
    face tmpTriFace(3);

    forAll(triFaces, indexI)
    {
        if
        (
            (triFaces[indexI] == c0BdyIndex[0])
         || (triFaces[indexI] == c0BdyIndex[1])
        )
        {
            continue;
        }

        if (checkNeighbour)
        {
            if
            (
                (triFaces[indexI] == c1BdyIndex[0])
             || (triFaces[indexI] == c1BdyIndex[1])
            )
            {
                continue;
            }
        }

        const face &triFace = faces_[triFaces[indexI]];

        forAll(triFace, pointI)
        {
            tmpTriFace[pointI] = triFace[pointI];

            if (triFace[pointI] == original[0])
            {
                tmpTriFace[pointI] = replacement[0];
            }

            if (triFace[pointI] == original[1])
            {
                tmpTriFace[pointI] = replacement[1];
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
            // Inverted and/or degenerate.
            return true;
        }
    }

    // No problems, so a collapse is feasible.
    return false;
}

// Utility method to check whether the cell given by 'cellIndex' will yield
// a valid cell when 'pointIndex' is moved to 'newPoint'. The routine performs
// metric-based checks. Returns 'true' if the collapse in NOT feasible, and
// makes entries in cellsChecked to avoid repetitive checks.
bool dynamicTopoFvMesh::checkCollapse
(
    const point& newPoint,
    const label pointIndex,
    const label cellIndex,
    labelHashSet& cellsChecked,
    bool forceOp
) const
{
    label faceIndex = -1;
    scalar cQuality = 0.0;
    const cell& cellToCheck = cells_[cellIndex];

    // Look for a face that doesn't contain 'pointIndex'
    forAll(cellToCheck, faceI)
    {
        const face& currFace = faces_[cellToCheck[faceI]];

        if (currFace.which(pointIndex) < 0)
        {
            faceIndex = cellToCheck[faceI];
            break;
        }
    }

    // Compute cell-volume
    const face& faceToCheck = faces_[faceIndex];

    if (owner_[faceIndex] == cellIndex)
    {
        cQuality =
        (
            (*tetMetric_)
            (
                points_[faceToCheck[2]],
                points_[faceToCheck[1]],
                points_[faceToCheck[0]],
                newPoint
            )
        );
    }
    else
    {
        cQuality =
        (
            (*tetMetric_)
            (
                points_[faceToCheck[0]],
                points_[faceToCheck[1]],
                points_[faceToCheck[2]],
                newPoint
            )
        );
    }

    // Final quality check
    if (cQuality < sliverThreshold_ && !forceOp)
    {
        if (debug > 2)
        {
            InfoIn("dynamicTopoFvMesh::checkCollapse()")
                << "\nCollapsing cell: " << cellIndex
                << " containing points:\n"
                << faceToCheck[0] << "," << faceToCheck[1] << ","
                << faceToCheck[2] << "," << pointIndex << nl
                << "will yield a quality of: " << cQuality
                << ", when " << pointIndex
                << " is moved to location: " << nl
                << newPoint << endl;
        }

        return true;
    }

    // Quality below 0.0 is a no-no
    if (cQuality < 0.0)
    {
        return true;
    }

    // No problems, so a collapse is feasible
    cellsChecked.insert(cellIndex);

    return false;
}

// Utility method to compute the quality of a vertex hull
// around an edge after bisection.
scalar dynamicTopoFvMesh::computeBisectionQuality
(
    const label eIndex
) const
{
    scalar minQuality = GREAT;
    scalar cQuality = 0.0;

    // Obtain a reference to this edge and corresponding edgePoints
    const edge& edgeToCheck = edges_[eIndex];
    const labelList& hullVertices = edgePoints_[eIndex];

    // Obtain point references
    const point& a = points_[edgeToCheck[0]];
    const point& c = points_[edgeToCheck[1]];

    // Compute the mid-point of the edge
    point midPoint = 0.5*(a + c);

    if (whichEdgePatch(eIndex) < 0)
    {
        // Internal edge.
        forAll(hullVertices, indexI)
        {
            label prevIndex = hullVertices.rcIndex(indexI);

            // Pick vertices off the list
            const point& b = points_[hullVertices[prevIndex]];
            const point& d = points_[hullVertices[indexI]];

            // Compute the quality of the upper half.
            cQuality = (*tetMetric_)(a, b, midPoint, d);

            // Check if the quality is worse
            minQuality = Foam::min(cQuality, minQuality);

            // Compute the quality of the lower half.
            cQuality = (*tetMetric_)(midPoint, b, c, d);

            // Check if the quality is worse
            minQuality = Foam::min(cQuality, minQuality);
        }
    }
    else
    {
        // Boundary edge.
        for(label indexI = 1; indexI < hullVertices.size(); indexI++)
        {
            // Pick vertices off the list
            const point& b = points_[hullVertices[indexI-1]];
            const point& d = points_[hullVertices[indexI]];

            // Compute the quality of the upper half.
            cQuality = (*tetMetric_)(a, b, midPoint, d);

            // Check if the quality is worse
            minQuality = Foam::min(cQuality, minQuality);

            // Compute the quality of the lower half.
            cQuality = (*tetMetric_)(midPoint, b, c, d);

            // Check if the quality is worse
            minQuality = Foam::min(cQuality, minQuality);
        }
    }

    // Ensure that the mesh is valid
    if (minQuality < sliverThreshold_)
    {
        if (debug > 3 && minQuality < 0.0)
        {
            // Write out cells for post processing.
            labelHashSet iCells;

            const labelList& eFaces = edgeFaces_[eIndex];

            forAll(eFaces, faceI)
            {
                if (!iCells.found(owner_[eFaces[faceI]]))
                {
                    iCells.insert(owner_[eFaces[faceI]]);
                }

                if (!iCells.found(neighbour_[eFaces[faceI]]))
                {
                    iCells.insert(neighbour_[eFaces[faceI]]);
                }
            }

            writeVTK(Foam::name(eIndex) + "_iCells", iCells.toc());
        }

        if (debug > 2)
        {
            InfoIn("dynamicTopoFvMesh::computeBisectionQuality()")
                << "Bisecting edge will fall below the "
                << "sliver threshold of: " << sliverThreshold_ << nl
                << "Edge: " << eIndex << ": " << edgeToCheck << nl
                << "EdgePoints: " << hullVertices << nl
                << "Minimum Quality: " << minQuality << nl
                << "Mid point: " << midPoint
                << endl;
        }
    }

    return minQuality;
}

// Utility method to compute the quality of cells
// around a face after trisection.
scalar dynamicTopoFvMesh::computeTrisectionQuality
(
    const label fIndex
) const
{
    scalar minQuality = GREAT;
    scalar cQuality = 0.0;

    point midPoint = triFaceCenter(fIndex);

    FixedList<label,2> apexPoint(-1);

    // Find the apex point
    apexPoint[0] = tetApexPoint(owner_[fIndex], fIndex);

    const face& faceToCheck = faces_[fIndex];

    forAll(faceToCheck, pointI)
    {
        // Pick vertices off the list
        const point& b = points_[faceToCheck[pointI]];
        const point& c = points_[apexPoint[0]];
        const point& d = points_[faceToCheck[faceToCheck.fcIndex(pointI)]];

        // Compute the quality of the upper half.
        cQuality = (*tetMetric_)(midPoint, b, c, d);

        // Check if the quality is worse
        minQuality = Foam::min(cQuality, minQuality);
    }

    if (whichPatch(fIndex) == -1)
    {
        apexPoint[1] = tetApexPoint(neighbour_[fIndex], fIndex);

        forAll(faceToCheck, pointI)
        {
            // Pick vertices off the list
            const point& b = points_[faceToCheck[pointI]];
            const point& c = points_[apexPoint[1]];
            const point& d = points_[faceToCheck[faceToCheck.rcIndex(pointI)]];

            // Compute the quality of the upper half.
            cQuality = (*tetMetric_)(midPoint, b, c, d);

            // Check if the quality is worse
            minQuality = Foam::min(cQuality, minQuality);
        }
    }

    return minQuality;
}

// Check if the boundary face is adjacent to a sliver-cell,
// and remove it by a two-step bisection/collapse operation.
void dynamicTopoFvMesh::remove2DSliver
(
    const label fIndex
)
{
    if (!edgeRefinement_)
    {
        return;
    }

    // Only boundary faces are considered.
    if (whichPatch(fIndex) == -1)
    {
        return;
    }

    // Measure the boundary edge-length of the face in question
    scalar length = edgeLength(getTriBoundaryEdge(fIndex));

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

// Indentify the sliver type in 3D
const changeMap dynamicTopoFvMesh::identifySliverType
(
    const label cIndex
)
{
    changeMap map;

    // Ensure that this cell actually exists.
    if (cells_[cIndex].empty())
    {
        return map;
    }

    label fourthPoint = -1;
    scalar minDistance = GREAT;
    face triFace(3), testFace(3), faceToCheck(3);
    FixedList<edge, 2> edgeToCheck(edge(-1,-1));

    const cell& cellToCheck = cells_[cIndex];

    // Find the point-face pair with minimum perpendicular distance
    forAll(cellToCheck, faceI)
    {
        label isolatedPoint = -1;
        label nextFaceI = cellToCheck.fcIndex(faceI);

        // Pick two faces from this cell.
        const face& currFace = faces_[cellToCheck[faceI]];
        const face& nextFace = faces_[cellToCheck[nextFaceI]];

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
                isolatedPoint = nextFace[pointI];

                // Configure a triangular face with correct orientation.
                if (owner_[cellToCheck[faceI]] == cIndex)
                {
                    testFace[0] = currFace[2];
                    testFace[1] = currFace[1];
                    testFace[2] = currFace[0];
                }
                else
                {
                    testFace[0] = currFace[0];
                    testFace[1] = currFace[1];
                    testFace[2] = currFace[2];
                }

                break;
            }
        }

        // Obtain the unit normal.
        vector testNormal = triFaceNormal(testFace);

        testNormal /= (mag(testNormal) + VSMALL);

        // Project the isolated point onto the face.
        vector p = points_[isolatedPoint] - points_[testFace[0]];
        vector q = p - ((p & testNormal)*testNormal);

        // Compute the distance
        scalar distance = mag(p - q);

        // Is it the least so far?
        if (distance < minDistance)
        {
            // Use this point-face pair.
            fourthPoint = isolatedPoint;
            triFace = testFace;
            minDistance = distance;
        }
    }

    // Obtain the face-normal.
    vector refArea = triFaceNormal(triFace);

    // Normalize it.
    vector n = refArea/mag(refArea);

    // Define edge-vectors.
    vector r1 = points_[triFace[1]] - points_[triFace[0]];
    vector r2 = points_[triFace[2]] - points_[triFace[1]];
    vector r3 = points_[triFace[0]] - points_[triFace[2]];

    // Project the fourth point onto the face.
    vector r4 = points_[fourthPoint] - points_[triFace[0]];

    r4 = r4 - ((r4 & n)*n);

    // Define the two other vectors.
    vector r5 = r4 - r1;
    vector r6 = r5 - r2;

    // Calculate three signed triangle areas, using triFace[0] as the origin.
    scalar t1 = n & (0.5 * (r1 ^ r4));
    scalar t2 = n & (0.5 * (r2 ^ r5));
    scalar t3 = n & (0.5 * (r3 ^ r6));

    // Determine sliver types based on are magnitudes.
    if (t1 > 0 && t2 > 0 && t3 > 0)
    {
        // Region R0: Cap cell.
        map.type() = 2;
        map.apexPoint() = fourthPoint;

        faceToCheck[0] = triFace[0];
        faceToCheck[1] = triFace[1];
        faceToCheck[2] = triFace[2];
    }

    if (t1 < 0 && t2 > 0 && t3 > 0)
    {
        // Region R1: Sliver cell.
        map.type() = 1;

        edgeToCheck[0][0] = triFace[2];
        edgeToCheck[0][1] = fourthPoint;
        edgeToCheck[1][0] = triFace[0];
        edgeToCheck[1][1] = triFace[1];
    }

    if (t1 > 0 && t2 < 0 && t3 > 0)
    {
        // Region R2: Sliver cell.
        map.type() = 1;

        edgeToCheck[0][0] = triFace[0];
        edgeToCheck[0][1] = fourthPoint;
        edgeToCheck[1][0] = triFace[1];
        edgeToCheck[1][1] = triFace[2];
    }

    if (t1 > 0 && t2 > 0 && t3 < 0)
    {
        // Region R3: Sliver cell.
        map.type() = 1;

        edgeToCheck[0][0] = triFace[1];
        edgeToCheck[0][1] = fourthPoint;
        edgeToCheck[1][0] = triFace[2];
        edgeToCheck[1][1] = triFace[0];
    }

    if (t1 < 0 && t2 > 0 && t3 < 0)
    {
        // Region R4: Cap cell.
        map.type() = 2;
        map.apexPoint() = triFace[0];

        faceToCheck[0] = triFace[1];
        faceToCheck[1] = triFace[2];
        faceToCheck[2] = fourthPoint;
    }

    if (t1 < 0 && t2 < 0 && t3 > 0)
    {
        // Region R5: Cap cell.
        map.type() = 2;
        map.apexPoint() = triFace[1];

        faceToCheck[0] = triFace[2];
        faceToCheck[1] = triFace[0];
        faceToCheck[2] = fourthPoint;
    }

    if (t1 > 0 && t2 < 0 && t3 < 0)
    {
        // Region R6: Cap cell.
        map.type() = 2;
        map.apexPoint() = triFace[2];

        faceToCheck[0] = triFace[0];
        faceToCheck[1] = triFace[1];
        faceToCheck[2] = fourthPoint;
    }

    // See if an over-ride to wedge/spade is necessary.
    // Obtain a reference area magnitude.
    scalar refMag = 0.1*(refArea & n);

    if (mag(t1) < refMag)
    {
        if (mag(t3) < refMag)
        {
            // Wedge case: Too close to point [0]
            map.type() = 4;

            edgeToCheck[0][0] = triFace[0];
            edgeToCheck[0][1] = fourthPoint;
        }
        else
        if (mag(t2) < refMag)
        {
            // Wedge case: Too close to point [1]
            map.type() = 4;

            edgeToCheck[0][0] = triFace[1];
            edgeToCheck[0][1] = fourthPoint;
        }
        else
        if ((mag(t2) > refMag) && (mag(t3) > refMag))
        {
            // Spade case: Too close to edge vector r1
            map.type() = 3;
            map.apexPoint() = fourthPoint;

            edgeToCheck[0][0] = triFace[0];
            edgeToCheck[0][1] = triFace[1];
        }
    }

    if (mag(t2) < refMag)
    {
        if (mag(t3) < refMag)
        {
            // Wedge case: Too close to point [2]
            map.type() = 4;

            edgeToCheck[0][0] = triFace[2];
            edgeToCheck[0][1] = fourthPoint;
        }
        else
        if ((mag(t1) > refMag) && (mag(t3) > refMag))
        {
            // Spade case: Too close to edge vector r2
            map.type() = 3;
            map.apexPoint() = fourthPoint;

            edgeToCheck[0][0] = triFace[1];
            edgeToCheck[0][1] = triFace[2];
        }
    }

    if (mag(t3) < refMag)
    {
        if ((mag(t1) > refMag) && (mag(t2) > refMag))
        {
            // Spade case: Too close to edge vector r3
            map.type() = 3;
            map.apexPoint() = fourthPoint;

            edgeToCheck[0][0] = triFace[2];
            edgeToCheck[0][1] = triFace[0];
        }
    }

    // Determine appropriate information for sliver exudation.
    if (map.type() == 1)
    {
        FixedList<bool, 2> foundEdge(false);

        // Search the cell-faces for first and second edges.
        forAll(cellToCheck, faceI)
        {
            const labelList& fEdges = faceEdges_[cellToCheck[faceI]];

            forAll(fEdges, edgeI)
            {
                const edge& thisEdge = edges_[fEdges[edgeI]];

                if (thisEdge == edgeToCheck[0])
                {
                    map.firstEdge() = fEdges[edgeI];

                    foundEdge[0] = true;
                }

                if (thisEdge == edgeToCheck[1])
                {
                    map.secondEdge() = fEdges[edgeI];

                    foundEdge[1] = true;
                }
            }

            if (foundEdge[0] && foundEdge[1])
            {
                break;
            }
        }
    }
    else
    if (map.type() == 2)
    {
        // Search the cell-faces for opposing faces.
        forAll(cellToCheck, faceI)
        {
            const face& thisFace = faces_[cellToCheck[faceI]];

            if (compare(thisFace, faceToCheck) != 0)
            {
                map.opposingFace() = cellToCheck[faceI];

                break;
            }
        }
    }
    else
    if (map.type() == 3 || map.type() == 4)
    {
        bool foundEdge = false;

        // Search the cell-faces for first edge.
        forAll(cellToCheck, faceI)
        {
            const labelList& fEdges = faceEdges_[cellToCheck[faceI]];

            forAll(fEdges, edgeI)
            {
                const edge& thisEdge = edges_[fEdges[edgeI]];

                if (thisEdge == edgeToCheck[0])
                {
                    map.firstEdge() = fEdges[edgeI];

                    foundEdge = true;
                }
            }

            if (foundEdge)
            {
                break;
            }
        }
    }

    if (debug > 2)
    {
        Pout << "Cell: " << cIndex
             << " Identified sliver type as: "
             << map.type() << endl;
    }

    // Return the result.
    return map;
}

// Remove sliver cells in 3D
void dynamicTopoFvMesh::removeSlivers()
{
    if (!edgeRefinement_)
    {
        return;
    }

    // If coupled patches exist, set the flag
    if (patchCoupling_.size() || procIndices_.size())
    {
        coupledModification_ = true;
    }

    forAllIter(Map<scalar>, thresholdSlivers_, iter)
    {
        // First check if this sliver cell is handled elsewhere.
        if (procIndices_.size())
        {
            bool foundInSubMesh = false;

            forAll(procIndices_, procI)
            {
                label proc = procIndices_[procI];

                if (proc < Pstream::myProcNo())
                {
                    Map<label>& rCellMap =
                    (
                        sendPatchMeshes_[proc].reverseEntityMap
                        (
                            coupledPatchInfo::CELL
                        )
                    );

                    if (rCellMap.found(iter.key()))
                    {
                        // This cell was sent to another sub-domain.
                        foundInSubMesh = true;
                        break;
                    }
                }
            }

            if (foundInSubMesh)
            {
                continue;
            }
        }

        // Small polyhedron reconnection for sliver removal
        bool useSPR = false;

        if (useSPR)
        {
            DynamicList<face> polyCellFaces(25);
            DynamicList<label> agCells(25), polyCell(25);

            // Agglomerate cells around the sliver for SPR
            agglomerateTetCells
            (
                iter.key(),
                agCells,
                polyCell,
                polyCellFaces
            );

            // Recursively call the SPR algorithm for the
            // best possible tetrahedralization
            bool success =
            (
                smallPolyhedronReconnection
                (
                    iter(),
                    polyCellFaces
                )
            );

            if (success)
            {
                continue;
            }
        }

        // Identify the sliver type.
        changeMap map = identifySliverType(iter.key());

        if (debug)
        {
            WarningIn
            (
                "dynamicTopoFvMesh::removeSlivers()"
            )   << nl << "Removing Cell: " << iter.key()
                << " of sliver type: " << map.type()
                << " with quality: " << iter() << endl;
        }

        // Take action based on the type of sliver.
        if (map.type() == 1)
        {
            // Sliver cell.
            // Determine which edges need to be bisected.
            label firstEdge = map.firstEdge();
            label secondEdge = map.secondEdge();

            // Force bisection on both edges.
            changeMap firstMap  = bisectEdge(firstEdge, false, true);
            changeMap secondMap = bisectEdge(secondEdge, false, true);

            // Collapse the intermediate edge.
            // Since we don't know which edge it is, search
            // through recently added edges and compare.
            edge edgeToCheck(firstMap.addedPoint(), secondMap.addedPoint());

            bool foundCollapseEdge = false;
            const labelList firstMapEdges = firstMap.addedEdgeList();
            const labelList secondMapEdges = secondMap.addedEdgeList();

            // Loop through the first list.
            forAll(firstMapEdges, edgeI)
            {
                const edge& thisEdge = edges_[firstMapEdges[edgeI]];

                if (thisEdge == edgeToCheck)
                {
                    // Collapse this edge.
                    collapseEdge
                    (
                        firstMapEdges[edgeI],
                        -1,
                        false,
                        true
                    );

                    foundCollapseEdge = true;
                    break;
                }
            }

            // Loop through the second list.
            if (!foundCollapseEdge)
            {
                forAll(secondMapEdges, edgeI)
                {
                    const edge& thisEdge = edges_[secondMapEdges[edgeI]];

                    if (thisEdge == edgeToCheck)
                    {
                        // Collapse this edge.
                        collapseEdge
                        (
                            secondMapEdges[edgeI],
                            -1,
                            false,
                            true
                        );

                        break;
                    }
                }
            }
        }
        else
        if (map.type() == 2)
        {
            // Cap cell.
            label opposingFace = map.opposingFace();

            // Force trisection of the opposing face.
            changeMap faceMap = trisectFace(opposingFace, false, true);

            // Collapse the intermediate edge.
            // Since we don't know which edge it is, search
            // through recently added edges and compare.
            edge edgeToCheck(map.apexPoint(), faceMap.addedPoint());

            const labelList faceMapEdges = faceMap.addedEdgeList();

            forAll(faceMapEdges, edgeI)
            {
                const edge& thisEdge = edges_[faceMapEdges[edgeI]];

                if (thisEdge == edgeToCheck)
                {
                    // Collapse this edge.
                    collapseEdge(faceMapEdges[edgeI], -1, false, true);

                    break;
                }
            }
        }
        else
        if (map.type() == 3)
        {
            // Spade cell.

            // Force bisection on the first edge.
            changeMap firstMap = bisectEdge(map.firstEdge(), false, true);

            // Collapse the intermediate edge.
            // Since we don't know which edge it is, search
            // through recently added edges and compare.
            edge edgeToCheck(firstMap.addedPoint(), firstMap.apexPoint());

            const labelList firstMapEdges = firstMap.addedEdgeList();

            // Loop through the first list.
            forAll(firstMapEdges, edgeI)
            {
                const edge& thisEdge = edges_[firstMapEdges[edgeI]];

                if (thisEdge == edgeToCheck)
                {
                    // Collapse this edge.
                    collapseEdge(firstMapEdges[edgeI], -1, false, true);

                    break;
                }
            }
        }
        else
        if (map.type() == 4)
        {
            // Wedge cell.

            // Collapse the first edge.
            collapseEdge(map.firstEdge(), -1, false, true);
        }
    }

    // Clear out the list
    thresholdSlivers_.clear();

    // If coupled patches exist, reset the flag
    if (patchCoupling_.size() || procIndices_.size())
    {
        coupledModification_ = false;
    }
}

// Agglomerate tetrahedral cells around a certain cell index
void dynamicTopoFvMesh::agglomerateTetCells
(
    const label cIndex,
    DynamicList<label>& aggCells,
    DynamicList<label>& polyCell,
    DynamicList<face>& polyCellFaces
) const
{
    // First insert the cell itself
    aggCells.append(cIndex);

    label nAggCells = 1, maxCells = 15;

    while (nAggCells < maxCells)
    {
        bool newCellsNotFound = true;

        forAll(aggCells, cellI)
        {
            // First loop through all faces of this cell,
            // and add neighbouring cells
            const cell& thisCell = cells_[aggCells[cellI]];

            forAll(thisCell, faceI)
            {
                label own = owner_[thisCell[faceI]];
                label nei = neighbour_[thisCell[faceI]];

                if (findIndex(aggCells, own) == -1)
                {
                    // Add this cell
                    aggCells.append(own);

                    newCellsNotFound = false;
                    nAggCells++;
                }

                if ((findIndex(aggCells, nei) == -1) && (nei != -1))
                {
                    // Add this cell
                    aggCells.append(nei);

                    newCellsNotFound = false;
                    nAggCells++;
                }
            }

            // Have we achieved the target?
            if (nAggCells >= maxCells)
            {
                break;
            }

            if (newCellsNotFound)
            {
                break;
            }
        }
    }

    // Shrink list to actual size
    aggCells.shrink();

    if (debug > 2)
    {
        writeVTK("aggCells_" + Foam::name(cIndex), aggCells);
    }

    // Construct a polyhedral cell with the agglomeration.
    forAll(aggCells, cellI)
    {
        const cell& thisCell = cells_[aggCells[cellI]];

        // Check if connected cells are not on the list.
        // These must be the bounding faces of the polyCell.
        forAll(thisCell, faceI)
        {
            label own = owner_[thisCell[faceI]];
            label nei = neighbour_[thisCell[faceI]];

            if (findIndex(aggCells, own) == -1)
            {
                if (findIndex(polyCell, thisCell[faceI]) == -1)
                {
                    polyCell.append(thisCell[faceI]);
                    polyCellFaces.append(faces_[thisCell[faceI]]);
                }
            }

            if (findIndex(aggCells, nei) == -1)
            {
                if (findIndex(polyCell, thisCell[faceI]) == -1)
                {
                    polyCell.append(thisCell[faceI]);
                    polyCellFaces.append(faces_[thisCell[faceI]].reverseFace());
                }
            }
        }
    }

    // Shrink lists to actual size
    polyCell.shrink();
    polyCellFaces.shrink();

    if (debug > 2)
    {
        writeVTK("polyCell_" + Foam::name(cIndex), polyCell, 2);
    }
}

bool dynamicTopoFvMesh::smallPolyhedronReconnection
(
    const scalar q0,
    DynamicList<face>& polyCellFaces
)
{
    scalar qc = q0;

    DynamicList<label> polyPoints(3*polyCellFaces.size());

    forAll(polyCellFaces, faceI)
    {
        const face& thisFace = polyCellFaces[faceI];

        if (thisFace.size() == 0)
        {
            continue;
        }

        // Pick this face, and build a candidate
        // list of points to loop through for
        // this polyCell
        polyPoints.clear();

        forAll(polyCellFaces, faceJ)
        {
            const face& nextFace = polyCellFaces[faceJ];

            forAll(nextFace, pointJ)
            {
                label pIndex = nextFace[pointJ];

                if
                (
                    (thisFace.which(pIndex) == -1) &&
                    (findIndex(polyPoints, pIndex) == -1)
                )
                {
                    polyPoints.append(pIndex);
                }
            }
        }

        if (polyPoints.size() == 1)
        {
            // This is the final tet cell.
            // Compute the quality and bail out.
            scalar cQuality = (*tetMetric_)
            (
                points_[thisFace[0]],
                points_[thisFace[1]],
                points_[thisFace[2]],
                points_[polyPoints[0]]
            );

            if (cQuality > qc)
            {
                qc = cQuality;
            }

            break;
        }
        else
        {
            // Loop through all points, and compute quality
            forAll(polyPoints, pointI)
            {
                scalar cQuality = (*tetMetric_)
                (
                    points_[thisFace[0]],
                    points_[thisFace[1]],
                    points_[thisFace[2]],
                    points_[polyPoints[pointI]]
                );

                if (cQuality < qc)
                {
                    continue;
                }

                // Looks like a better quality is available
                // from this face. Check the sub-mesh as well.

            }
        }
    }

    // Look through all accumulated sub-quality values,
    // and determine whether a better quality was available.
    if (qc > q0)
    {
        return true;
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
    // Initialize coupled patches for topology modifications.
    initCoupledPatches();

    // Handle coupled patches.
    handleCoupledPatches();

    if (edgeRefinement_)
    {
        // Initialize the face stacks
        initFaceStacks();

        if (threader_->multiThreaded())
        {
            // Lock slaves
            lockSlaveThreads(topoSeq_, handlerPtr_);

            // Submit jobs to the work queue
            forAll(topoSeq_, i)
            {
                threader_->addToWorkQueue
                (
                    &edgeBisectCollapse2D,
                    reinterpret_cast<void *>(&(handlerPtr_[topoSeq_[i]]))
                );

                // Wait for a signal from this thread
                // before moving on.
                handlerPtr_[topoSeq_[i]].waitForSignal
                (
                    threadHandler<dynamicTopoFvMesh>::START
                );
            }

            // Synchronize threads
            synchronizeThreads(topoSeq_, handlerPtr_);
        }

        // Set the master thread to implement modifications
        edgeBisectCollapse2D(&(handlerPtr_[0]));

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
        lockSlaveThreads(topoSeq_, handlerPtr_);

        // Submit jobs to the work queue
        forAll(topoSeq_, i)
        {
            threader_->addToWorkQueue
            (
                &swap2DEdges,
                reinterpret_cast<void *>(&(handlerPtr_[topoSeq_[i]]))
            );

            // Wait for a signal from this thread
            // before moving on.
            handlerPtr_[topoSeq_[i]].waitForSignal
            (
                threadHandler<dynamicTopoFvMesh>::START
            );
        }

        // Synchronize threads
        synchronizeThreads(topoSeq_, handlerPtr_);
    }

    // Set the master thread to implement modifications
    swap2DEdges(reinterpret_cast<void *>(&(handlerPtr_[0])));

    if (debug)
    {
        Info << nl << "2D Edge Swapping complete." << endl;
    }
}

// MultiThreaded topology modifier [3D]
void dynamicTopoFvMesh::threadedTopoModifier3D()
{
    // Initialize coupled patches for topology modifications.
    initCoupledPatches();

    // Remove sliver cells first.
    // removeSlivers();

    // Handle coupled patches.
    handleCoupledPatches();

    if (edgeRefinement_)
    {
        // Initialize the cell stacks
        initEdgeStacks();

        if (threader_->multiThreaded())
        {
            // Lock slave threads
            lockSlaveThreads(topoSeq_, handlerPtr_);

            // Submit jobs to the work queue
            forAll(topoSeq_, i)
            {
                threader_->addToWorkQueue
                (
                    &edgeBisectCollapse3D,
                    reinterpret_cast<void *>(&(handlerPtr_[topoSeq_[i]]))
                );

                // Wait for a signal from this thread
                // before moving on.
                handlerPtr_[topoSeq_[i]].waitForSignal
                (
                    threadHandler<dynamicTopoFvMesh>::START
                );
            }

            // Synchronize threads
            synchronizeThreads(topoSeq_, handlerPtr_);
        }

        // Set the master thread to implement modifications
        edgeBisectCollapse3D(&(handlerPtr_[0]));

        // Handle mesh slicing events, if necessary
        handleMeshSlicing();

        if (debug)
        {
            Info << nl << "3D Edge Bisection/Collapse complete." << endl;
        }
    }

    // Initialize the edge stacks
    initEdgeStacks();

    if (threader_->multiThreaded())
    {
        // Lock slave threads
        lockSlaveThreads(topoSeq_, handlerPtr_);

        // Submit jobs to the work queue
        forAll(topoSeq_, i)
        {
            threader_->addToWorkQueue
            (
                &swap3DEdges,
                reinterpret_cast<void *>(&(handlerPtr_[topoSeq_[i]]))
            );

            // Wait for a signal from this thread
            // before moving on.
            handlerPtr_[topoSeq_[i]].waitForSignal
            (
                threadHandler<dynamicTopoFvMesh>::START
            );
        }

        // Synchronize threads
        synchronizeThreads(topoSeq_, handlerPtr_);
    }

    // Set the master thread to implement modifications
    swap3DEdges(reinterpret_cast<void *>(&(handlerPtr_[0])));

    if (debug)
    {
        Info << nl << "3D Edge Swapping complete." << endl;
    }
}

// Lock all slave threads
template <class Type>
void dynamicTopoFvMesh::lockSlaveThreads
(
    const labelList& sequence,
    PtrList<threadHandler<Type> >& handler
)
{
    forAll(sequence, i)
    {
        handler[sequence[i]].lock(threadHandler<Type>::START);
        handler[sequence[i]].lock(threadHandler<Type>::STOP);

        handler[sequence[i]].unsetPredicate(threadHandler<Type>::START);
        handler[sequence[i]].unsetPredicate(threadHandler<Type>::STOP);
    }
}

// Synchronize all slave threads
template <class Type>
void dynamicTopoFvMesh::synchronizeThreads
(
    const labelList& sequence,
    PtrList<threadHandler<Type> >& handler
)
{
    forAll(sequence, i)
    {
        // Wait for a signal from this thread
        // before moving on.
        handler[sequence[i]].waitForSignal(threadHandler<Type>::STOP);
    }
}

// Return reference to the dictionary
const IOdictionary& dynamicTopoFvMesh::dynamicMeshDict() const
{
    return dict_;
}

// Update the mesh for topology changes
// Return true if changes have occurred
bool dynamicTopoFvMesh::updateTopology()
{
    // Re-read options if they have been modified at run-time
    if (dict_.readIfModified())
    {
        // Re-read optional parameters
        readOptionalParameters();

        // Re-read edge options
        readRefinementOptions(true);
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

    forAll(currentPoints, pointI)
    {
        points_[pointI] = currentPoints[pointI];
    }

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

    bool sliversAbsent = true;

    if (thresholdSlivers_.size())
    {
        sliversAbsent = false;
    }

    // Reduce across processors.
    reduce(sliversAbsent, andOp<bool>());

    // Return if the interval is invalid.
    // Handy while using only mesh-motion.
    if (interval_ < 0)
    {
        return false;
    }

    // Return if re-meshing is not at interval,
    // or sliver cells are absent.
    if ((time().timeIndex() % interval_ != 0) && sliversAbsent)
    {
        return false;
    }

    // Calculate the edge length-scale for the mesh
    calculateLengthScale();

    // Track mesh topology modification time
    clockTime topologyTimer;

    // Invoke the threaded topoModifier
    if (twoDMesh_)
    {
        threadedTopoModifier2D();
    }
    else
    {
        threadedTopoModifier3D();
    }

    Info << " Topo modifier time: " << topologyTimer.elapsedTime() << endl;

    clockTime reOrderingTimer;

    // Reduce across processors.
    reduce(topoChangeFlag_, orOp<bool>());

    // Apply all pending topology changes, if necessary
    if (topoChangeFlag_)
    {
        // Obtain references to zones, if any
        pointZoneMesh& pointZones = polyMesh::pointZones();
        faceZoneMesh& faceZones = polyMesh::faceZones();
        cellZoneMesh& cellZones = polyMesh::cellZones();

        // Allocate temporary lists for mesh-reset
        pointField points(nPoints_);
        edgeList edges(nEdges_);
        faceList faces(nFaces_);
        labelList owner(nFaces_);
        labelList neighbour(nInternalFaces_);
        labelListList faceEdges(nFaces_);
        labelListList edgeFaces(nEdges_);
        labelList oldPatchStarts(oldPatchStarts_);
        labelList oldPatchNMeshPoints(oldPatchNMeshPoints_);
        labelListList pointZoneMap(pointZones.size());
        labelListList faceZonePointMap(faceZones.size());
        labelListList faceZoneFaceMap(faceZones.size());
        labelListList cellZoneMap(cellZones.size());

        // Null temporaries
        List<objectMap> pointsFromPoints(pointsFromPoints_.size());
        List<objectMap> facesFromPoints(0);
        List<objectMap> facesFromEdges(0);
        List<objectMap> facesFromFaces(0);
        List<objectMap> cellsFromPoints(0);
        List<objectMap> cellsFromEdges(0);
        List<objectMap> cellsFromFaces(0);
        List<objectMap> cellsFromCells(cellsFromCells_.size());
        labelHashSet flipFaceFlux(0);
        pointField preMotionPoints(0);

        // Obtain faceZone point maps before reordering
        List<Map<label> > oldFaceZonePointMaps(faceZones.size());

        forAll(faceZones, fzI)
        {
            oldFaceZonePointMaps[fzI] = faceZones[fzI]().meshPointMap();
        }

        // Reorder the mesh and obtain current topological information
        reOrderMesh
        (
            points,
            edges,
            faces,
            owner,
            neighbour,
            faceEdges,
            edgeFaces,
            pointZoneMap,
            faceZoneFaceMap,
            cellZoneMap
        );

        // Copy point-mapping information
        label indexI = 0;

        forAllIter(Map<objectMap>, pointsFromPoints_, pointI)
        {
            pointsFromPoints[indexI++] = pointI();
        }

        // Copy cell-mapping information
        label indexC = 0;

        forAllIter(Map<objectMap>, cellsFromCells_, cellI)
        {
            cellsFromCells[indexC++] = cellI();
        }

        // Obtain the patch-point maps before resetting the mesh
        List<Map<label> > oldPatchPointMaps(numPatches_);

        forAll(oldPatchPointMaps, patchI)
        {
            oldPatchPointMaps[patchI] = boundaryMesh()[patchI].meshPointMap();
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

        // Reset the edge mesh
        eMeshPtr_->resetPrimitives
        (
            edges,
            faceEdges,
            edgeFaces,
            edgePatchSizes_,
            edgePatchStarts_
        );

        // Generate mapping for points on boundary patches
        labelListList patchPointMap(numPatches_);

        for (label i = 0; i < numPatches_; i++)
        {
            // Obtain new patch mesh points after reset.
            const labelList& meshPointLabels = boundaryMesh()[i].meshPoints();

            patchNMeshPoints_[i] = meshPointLabels.size();

            patchPointMap[i].setSize(meshPointLabels.size(), -1);

            forAll(meshPointLabels, pointI)
            {
                label oldIndex = pointMap_[meshPointLabels[pointI]];

                // Check if the point existed before...
                if (oldIndex > -1)
                {
                    // Look for the old position on this patch.
                    Map<label>::const_iterator oIter =
                    (
                        oldPatchPointMaps[i].find(oldIndex)
                    );

                    // Add an entry if the point was found
                    if (oIter != oldPatchPointMaps[i].end())
                    {
                        patchPointMap[i][pointI] = oIter();
                    }
                }
            }
        }

        // Generate mapping for faceZone points
        forAll(faceZones, fzI)
        {
            // Obtain new face zone mesh points after reset.
            const labelList& meshPointLabels = faceZones[fzI]().meshPoints();

            faceZonePointMap[fzI].setSize(meshPointLabels.size(), -1);

            forAll(meshPointLabels, pointI)
            {
                label oldIndex = pointMap_[meshPointLabels[pointI]];

                // Check if the point existed before...
                if (oldIndex > -1)
                {
                    // Look for the old position on this patch.
                    Map<label>::const_iterator oIter =
                    (
                        oldFaceZonePointMaps[fzI].find(oldIndex)
                    );

                    // Add an entry if the point was found
                    if (oIter != oldFaceZonePointMaps[fzI].end())
                    {
                        faceZonePointMap[fzI][pointI] = oIter();
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

        // Clear unwanted member data
        addedFacePatches_.clear();
        addedEdgePatches_.clear();
        addedPointZones_.clear();
        addedFaceZones_.clear();
        addedCellZones_.clear();
        cellsFromCells_.clear();
        pointsFromPoints_.clear();
        cellParents_.clear();
        pointParents_.clear();

        // Set new sizes for the reverse maps
        reversePointMap_.setSize(nPoints_, -7);
        reverseEdgeMap_.setSize(nEdges_, -7);
        reverseFaceMap_.setSize(nFaces_, -7);
        reverseCellMap_.setSize(nCells_, -7);

        // Update "old" information
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

        // Basic checks for mesh-validity
        if (debug > 1)
        {
            checkMesh(true);
        }
    }

    // Print out topo-stats
    Info << " Reordering time: " << reOrderingTimer.elapsedTime() << endl;
    Info << " nBisections: " << nBisections_ << endl;
    Info << " nCollapses: " << nCollapses_ << endl;
    Info << " nSwaps: " << nSwaps_ << endl;

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
