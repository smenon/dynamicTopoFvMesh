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

#include "dynamicTopoFvMesh.H"
#include "addToRunTimeSelectionTable.H"

#include "clockTime.H"
#include "mapPolyMesh.H"
#include "volFields.H"
#include "motionSolver.H"
#include "MapTopoFvFields.H"
#include "fvPatchFields.H"
#include "fvsPatchFields.H"
#include "MeshObject.H"
#include "topoMapper.H"
#include "SortableList.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynamicTopoFvMesh,0);

addToRunTimeSelectionTable(dynamicFvMesh, dynamicTopoFvMesh, IOobject);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from IOobject
dynamicTopoFvMesh::dynamicTopoFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
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
    bandWidthReduction_(false),
    coupledModification_(false),
    interval_(1),
    mapper_(NULL),
    mPtr_(NULL),
    oldPoints_(polyMesh::points()),
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
    proximityBins_(0),
    sliceThreshold_(VSMALL),
    sliceHoldOff_(0),
    sliceBoxes_(0),
    slicePairs_(0),
    maxTetsPerEdge_(-1),
    allowTableResize_(false)
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

    const polyBoundaryMesh& boundary = polyMesh::boundaryMesh();

    // Initialize the multiThreading environment
    initializeThreadingEnvironment();

    // Read optional parameters.
    readOptionalParameters();

    // Initialize patch-size information
    for(label i=0; i<numPatches_; i++)
    {
        patchNMeshPoints_[i] = boundary[i].meshPoints().size();
        oldPatchSizes_[i] = patchSizes_[i]  = boundary[i].size();
        oldPatchStarts_[i] = patchStarts_[i] = boundary[i].start();
        oldPatchNMeshPoints_[i] = patchNMeshPoints_[i];
    }

    // Open the tetMetric dynamic-link library (for 3D only)
    loadMetric();

    // Initialize edge-related connectivity structures
    initEdges();

    // Initialize coupled patch connectivity for topology modifications.
    // This needs to be done before the motionSolver is initialized.
    initCoupledConnectivity(this);

    // Load the mesh-motion solver
    loadMotionSolver();

    // Load the field-mapper
    loadFieldMapper();

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
    const labelList& owner,
    const labelList& neighbour,
    const cellList& cells
)
:
    dynamicFvMesh(io, points, faces, owner, neighbour, false),
    numPatches_(1),
    topoChangeFlag_(false),
    isSubMesh_(true),
    dict_(mesh.dict_),
    mandatory_(mesh.mandatory_),
    twoDMesh_(mesh.twoDMesh_),
    edgeRefinement_(mesh.edgeRefinement_),
    bandWidthReduction_(mesh.bandWidthReduction_),
    coupledModification_(false),
    interval_(1),
    mapper_(NULL),
    mPtr_(NULL),
    oldPoints_(points),
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
    nOldInternalFaces_(primitiveMesh::nInternalFaces()),
    nInternalFaces_(primitiveMesh::nInternalFaces()),
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
    proximityBins_(0),
    sliceThreshold_(VSMALL),
    sliceHoldOff_(0),
    sliceBoxes_(0),
    slicePairs_(0),
    maxTetsPerEdge_(mesh.maxTetsPerEdge_),
    allowTableResize_(mesh.allowTableResize_),
    tetMetric_(mesh.tetMetric_)
{
    // Initialize owner and neighbour
    owner_.setSize(faces.size(), -1);
    neighbour_.setSize(faces.size(), -1);

    // Set owner and neighbour from polyMesh
    const labelList& own = polyMesh::faceOwner();
    const labelList& nei = polyMesh::faceNeighbour();

    forAll(own, faceI)
    {
        owner_[faceI] = own[faceI];
    }

    forAll(nei, faceI)
    {
        neighbour_[faceI] = nei[faceI];
    }

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

    // Size-up edgePoints for now, but explicitly construct
    // for each edge later, based on point coupling.
    if (!twoDMesh_)
    {
        invertConnectivity(nPoints_, edges_, pointEdges_);
        edgePoints_.setSize(nEdges_, labelList(0));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dynamicTopoFvMesh::~dynamicTopoFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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


// Obtain map weighting factors for a tetrahedral cell.
//  - Returns true when weights are consistent (i.e., sum to 1.0)
bool dynamicTopoFvMesh::computeTetWeights
(
    const label cIndex,
    const scalar cellVolume,
    const labelList& mappingCells,
    const scalar searchFactor,
    labelList& parents,
    scalarField& weights
) const
{
    // Obtain candidate parents for this cell
    labelList candidates = cellParents(cIndex, searchFactor, mappingCells);

    // Track actual intersections
    label nIntersects = 0;

    // Compute intersection volumes with candidates
    scalarField intVolumes(candidates.size(), 0.0);

    forAll(candidates, indexI)
    {
        intVolumes[indexI] =
        (
            tetIntersection
            (
                cIndex,
                candidates[indexI]
            )
        );

        if (intVolumes[indexI] > 0.0)
        {
            nIntersects++;
        }
    }

    // Now copy only valid intersections.
    parents.setSize(nIntersects, -1);
    weights.setSize(nIntersects, 0.0);

    // Reset counter
    nIntersects = 0;

    forAll(intVolumes, indexI)
    {
        if (intVolumes[indexI] > 0.0)
        {
            parents[nIntersects] = candidates[indexI];
            weights[nIntersects] = intVolumes[indexI];

            nIntersects++;
        }
    }

    // Test the weights for consistency
    scalarField testWeights = (weights/cellVolume);

    if (mag(1.0 - sum(testWeights)) > 1e-10)
    {
        // Inconsistent weights. Check whether any edges
        // lie on boundary patches. These cells can have
        // relaxed weights to account for mild convexity.
        bool foundBoundary = false;

        const cell& cellToCheck = cells_[cIndex];

        forAll(cellToCheck, faceI)
        {
            const labelList& fEdges = faceEdges_[cellToCheck[faceI]];

            forAll(fEdges, edgeI)
            {
                if (whichEdgePatch(fEdges[edgeI]) > -1)
                {
                    foundBoundary = true;
                    break;
                }
            }

            if (foundBoundary)
            {
                break;
            }
        }

        if (foundBoundary)
        {
            // Normalize by sum of intersections
            weights /= sum(weights);
        }
        else
        {
            // Weights are inconsistent. Notify caller.
            weights /= cellVolume;

            return false;
        }
    }
    else
    {
        // Normalize by volume
        weights /= cellVolume;
    }

    // Return consistent weights
    return true;
}


// Determine the intersection volume between two tetrahedra
scalar dynamicTopoFvMesh::tetIntersection
(
    const label newCellIndex,
    const label oldCellIndex
) const
{
    scalar intVol = 0.0, tolFactor = 1e-8;

    // For post-processing purposes, define a name
    word cvxName
    (
        "cvxSet_"
      + Foam::name(newCellIndex)
      + '_'
      + Foam::name(oldCellIndex)
    );

    bool intersects = false;

    // Check indices first
    if (newCellIndex >= cells_.size() || newCellIndex < 0)
    {
        FatalErrorIn("dynamicTopoFvMesh::tetIntersection")
            << " Wrong newCellIndex: " << newCellIndex << nl
            << " nCells: " << nCells_
            << abort(FatalError);
    }

    if (oldCellIndex >= nOldCells_ || oldCellIndex < 0)
    {
        FatalErrorIn("dynamicTopoFvMesh::tetIntersection")
            << " Wrong oldCellIndex: " << oldCellIndex << nl
            << " nOldCells: " << nOldCells_
            << abort(FatalError);
    }

    // Fetch connectivity.
    const cellList& oldCells = polyMesh::cells();
    const faceList& oldFaces = polyMesh::faces();

    const cell& newCell = cells_[newCellIndex];
    const cell& oldCell = oldCells[oldCellIndex];

    // Check topologically for common points / edges / faces.
    // These count as intersections.
    FixedList<label, 4> oldCellPoints(-1), newCellPoints(-1);
    FixedList<edge, 6> oldEdges(edge(-1,-1));
    FixedList<label, 6> newEdges(-1);

    const face& oldBaseFace = oldFaces[oldCell[0]];
    const face& newBaseFace = faces_[newCell[0]];

    label oldApex =
    (
        findIsolatedPoint
        (
            oldBaseFace,
            oldFaces[oldCell[1]]
        )
    );

    label newApex =
    (
        findIsolatedPoint
        (
            newBaseFace,
            faces_[newCell[1]]
        )
    );

    // Build the old / new points list
    oldCellPoints[0] = oldBaseFace[0];
    oldCellPoints[1] = oldBaseFace[1];
    oldCellPoints[2] = oldBaseFace[2];
    oldCellPoints[3] = oldApex;

    newCellPoints[0] = newBaseFace[0];
    newCellPoints[1] = newBaseFace[1];
    newCellPoints[2] = newBaseFace[2];
    newCellPoints[3] = newApex;

    // Build the old edge list
    oldEdges[0] = edge(oldBaseFace[0], oldBaseFace[1]);
    oldEdges[1] = edge(oldBaseFace[1], oldBaseFace[2]);
    oldEdges[2] = edge(oldBaseFace[2], oldBaseFace[0]);
    oldEdges[3] = edge(oldBaseFace[0], oldApex);
    oldEdges[4] = edge(oldBaseFace[1], oldApex);
    oldEdges[5] = edge(oldBaseFace[2], oldApex);

    // Get a list of edge indices for the new cell.
    label nEdg = 0;

    forAll(newCell, fI)
    {
        const labelList& fEdges = faceEdges_[newCell[fI]];

        forAll(fEdges, edgeI)
        {
            if (findIndex(newEdges, fEdges[edgeI]) == -1)
            {
                newEdges[nEdg++] = fEdges[edgeI];
            }
        }

        if (nEdg == 6)
        {
            break;
        }
    }

    // Track all possible intersections from here on.
    label nInts = 0;
    vectorField tP(0);
    vector intPoint = vector::zero;
    FixedList<vector,2> segment(vector::zero);

    // Check whether any old points are within
    // the new cell. Count these as 'intersections'.
    bool allCommon = true;

    forAll(oldCellPoints, pI)
    {
        bool foundUnique = true;

        forAll(newCellPoints, pJ)
        {
            if (oldCellPoints[pI] == newCellPoints[pJ])
            {
                foundUnique = false;
                break;
            }
        }

        if (foundUnique)
        {
            const point& checkPoint = oldPoints_[oldCellPoints[pI]];

            if (pointInTet(newCellIndex, checkPoint, false, true))
            {
                sizeUpList(checkPoint, tP);

                nInts++;
            }

            // Set the flag
            allCommon = false;
        }
    }

    if (allCommon)
    {
        // Looks like this cell is identical to the old cell.
        tP.setSize(4, vector::zero);

        forAll(newCellPoints, pI)
        {
            tP[pI] = oldPoints_[newCellPoints[pI]];
        }

        // Compute intersection volume.
        intVol =
        (
            convexSetVolume
            (
                cvxName,
                tolFactor,
                tP
            )
        );

        intersects = true;

        return intVol;
    }

    // Common flags
    bool foundCommonFace = false;
    bool foundCommonEdge = false;
    bool foundCommonPoint = false;

    // Check for a common face.
    // Note that two common faces cannot occur,
    // since we've already checked points.
    label nCo = 0, nCn = 0;
    FixedList<label, 3> commonOldEdgeIndices(-1);
    FixedList<label, 3> commonNewEdgeIndices(-1);

    forAll(oldCell, faceI)
    {
        const face& oldFace = oldFaces[oldCell[faceI]];

        forAll(newCell, faceJ)
        {
            const face& newFace = faces_[newCell[faceJ]];

            if (triFaceCompare(oldFace, newFace))
            {
                // Fill in first three points from the common face.
                const face& commonFace = faces_[newCell[faceJ]];

                // Add to the list.
                forAll(commonFace, pI)
                {
                    sizeUpList(oldPoints_[commonFace[pI]], tP);
                }

                nInts += 3;

                foundCommonFace = true;

                // Also identify common edges.
                const labelList& fEdges = faceEdges_[newCell[faceJ]];

                forAll(fEdges, edgeI)
                {
                    forAll(oldEdges, edgeJ)
                    {
                        if (oldEdges[edgeJ] == edges_[fEdges[edgeI]])
                        {
                            commonOldEdgeIndices[nCo++] = edgeJ;
                            commonNewEdgeIndices[nCn++] = fEdges[edgeI];

                            break;
                        }
                    }
                }

                break;
            }
        }

        if (foundCommonFace)
        {
            break;
        }
    }

    // Check for one common edge.
    // Note that two common edges cannot occur,
    // since we've already checked faces.
    if (!foundCommonFace)
    {
        forAll(newEdges, edgeI)
        {
            forAll(oldEdges, edgeJ)
            {
                if (oldEdges[edgeJ] == edges_[newEdges[edgeI]])
                {
                    commonOldEdgeIndices[nCo++] = edgeJ;
                    commonNewEdgeIndices[nCn++] = newEdges[edgeI];

                    // Fill in first two points from the common edge.
                    const edge& commonEdge = edges_[newEdges[edgeI]];

                    // Add to the list.
                    forAll(commonEdge, pI)
                    {
                        sizeUpList(oldPoints_[commonEdge[pI]], tP);
                    }

                    nInts += 2;

                    foundCommonEdge = true;

                    break;
                }
            }

            if (foundCommonEdge)
            {
                break;
            }
        }
    }

    // If a common edge wasn't found, look for a common point.
    // Obviously, two common points cannot occur,
    // since we've already checked edges.
    if (!foundCommonFace && !foundCommonEdge)
    {
        forAll(newCellPoints, pointI)
        {
            forAll(oldCellPoints, pointJ)
            {
                if (oldCellPoints[pointJ] == newCellPoints[pointI])
                {
                    // Add the common point to the list.
                    sizeUpList(oldPoints_[oldCellPoints[pointJ]], tP);

                    nInts++;

                    foundCommonPoint = true;

                    break;
                }
            }

            if (foundCommonPoint)
            {
                break;
            }
        }
    }

    // Loop through all new edges, and find possible intersections.
    forAll(newEdges, edgeI)
    {
        // Avoid common edges.
        if (findIndex(commonNewEdgeIndices, newEdges[edgeI]) > -1)
        {
            continue;
        }

        const edge& edgeToCheck = edges_[newEdges[edgeI]];

        // Deal with segment-point intersections first.
        bool foundIntersection = false;

        forAll(oldCellPoints, pI)
        {
            // Skip common points
            if
            (
                (oldCellPoints[pI] == edgeToCheck[0]) ||
                (oldCellPoints[pI] == edgeToCheck[1])
            )
            {
                continue;
            }

            segment[0] = oldPoints_[edgeToCheck[0]];
            segment[1] = oldPoints_[edgeToCheck[1]];

            foundIntersection =
            (
                segmentPointIntersection
                (
                    segment,
                    edgeToCheck,
                    oldCellPoints[pI],
                    tolFactor,
                    true
                )
            );

            if (foundIntersection)
            {
                // Add to the list.
                sizeUpList(oldPoints_[oldCellPoints[pI]], tP);

                nInts++;

                break;
            }
        }

        if (!foundIntersection)
        {
            forAll(oldCell, fI)
            {
                const face& oldCheckFace = oldFaces[oldCell[fI]];

                // Ensure that face doesn't contain edgeToCheck
                if
                (
                    (oldCheckFace.which(edgeToCheck[0]) > -1) &&
                    (oldCheckFace.which(edgeToCheck[1]) > -1)
                )
                {
                    continue;
                }

                segment[0] = oldPoints_[edgeToCheck[0]];
                segment[1] = oldPoints_[edgeToCheck[1]];

                // Reset flag
                foundIntersection = false;

                foundIntersection =
                (
                    segmentTriFaceIntersection
                    (
                        segment,
                        edgeToCheck,
                        oldCheckFace,
                        tolFactor,
                        intPoint,
                        true
                    )
                );

                if (foundIntersection)
                {
                    // Add to the list.
                    sizeUpList(intPoint, tP);

                    nInts++;
                }
            }
        }
    }

    // Loop through all old edges, and find possible
    // intersections with new cell faces.
    forAll(oldEdges, edgeI)
    {
        // Avoid common edges.
        if (findIndex(commonOldEdgeIndices, edgeI) > -1)
        {
            continue;
        }

        const edge& edgeToCheck = oldEdges[edgeI];

        // Deal with segment-point intersections first.
        bool foundIntersection = false;

        forAll(newCellPoints, pI)
        {
            // Skip common points
            if
            (
                (newCellPoints[pI] == edgeToCheck[0]) ||
                (newCellPoints[pI] == edgeToCheck[1])
            )
            {
                continue;
            }

            segment[0] = oldPoints_[edgeToCheck[0]];
            segment[1] = oldPoints_[edgeToCheck[1]];

            foundIntersection =
            (
                segmentPointIntersection
                (
                    segment,
                    edgeToCheck,
                    newCellPoints[pI],
                    tolFactor,
                    true
                )
            );

            if (foundIntersection)
            {
                // Add to the list.
                sizeUpList(oldPoints_[newCellPoints[pI]], tP);

                nInts++;

                break;
            }
        }

        if (!foundIntersection)
        {
            forAll(newCell, fI)
            {
                const face& newCheckFace = faces_[newCell[fI]];

                // Ensure that face doesn't contain edgeToCheck
                if
                (
                    (newCheckFace.which(edgeToCheck[0]) > -1) &&
                    (newCheckFace.which(edgeToCheck[1]) > -1)
                )
                {
                    continue;
                }

                segment[0] = oldPoints_[edgeToCheck[0]];
                segment[1] = oldPoints_[edgeToCheck[1]];

                // Reset flag
                foundIntersection = false;

                foundIntersection =
                (
                    segmentTriFaceIntersection
                    (
                        segment,
                        edgeToCheck,
                        newCheckFace,
                        tolFactor,
                        intPoint,
                        true
                    )
                );

                if (foundIntersection)
                {
                    // Add to the list.
                    sizeUpList(intPoint, tP);

                    nInts++;
                }
            }
        }
    }

    // Found a polyhedral intersecting volume.
    // Compute the volume from points and return.
    if (nInts >= 4)
    {
        intVol =
        (
            convexSetVolume
            (
                cvxName,
                tolFactor,
                tP
            )
        );

        intersects = true;
    }

    // Return intersection volume.
    return intVol;
}

// Compute the volume of a polyhedron
// formed by a convex set of points.
scalar dynamicTopoFvMesh::convexSetVolume
(
    const word& cvxSetName,
    const scalar tolFraction,
    const vectorField& cvxSet
) const
{
    scalar cVol = 0.0;
    face tmpFace(3);

    DynamicList<face> testFaces(10);

    // Loop through all points, and build faces with every
    // other point in the set
    forAll(cvxSet, i)
    {
        forAll(cvxSet, j)
        {
            // Skip duplicates.
            if (j == i)
            {
                continue;
            }

            forAll(cvxSet, k)
            {
                // Skip duplicates.
                if (k == i || k == j)
                {
                    continue;
                }

                // Configure the face.
                tmpFace[0] = i;
                tmpFace[1] = j;
                tmpFace[2] = k;

                // Compute the normal to this face
                vector n = tmpFace.normal(cvxSet);

                n /= mag(n) + VSMALL;

                // Include all other co-planar points
                forAll(cvxSet, l)
                {
                    // Skip duplicates.
                    if (findIndex(tmpFace, l) > -1)
                    {
                        continue;
                    }

                    vector rfVec = (cvxSet[l] - cvxSet[i]);
                    scalar dotProd = (rfVec/mag(rfVec)) & n;

                    if (mag(dotProd) < tolFraction*mag(rfVec))
                    {
                        // Need to configure a new face.
                        face newFace(3);
                        bool foundLocation = false;

                        forAll(tmpFace, pI)
                        {
                            label nI = tmpFace.fcIndex(pI);

                            newFace[0] = tmpFace[pI];
                            newFace[1] = l;
                            newFace[2] = tmpFace[nI];

                            // Compute the normal.
                            vector nNew = newFace.normal(cvxSet);

                            if ((n & nNew) > 0.0)
                            {
                                // Insert the point.
                                insertLabel
                                (
                                    l,
                                    tmpFace[pI],
                                    tmpFace[nI],
                                    tmpFace
                                );

                                foundLocation = true;

                                break;
                            }
                        }

                        if (!foundLocation)
                        {
                            FatalErrorIn
                            (
                                "label dynamicTopoFvMesh::convexSetVolume()"
                            )   << "Cannot find appropriate configuration."
                                << nl << "Face: " << tmpFace.points(cvxSet)
                                << nl << " with point: " << cvxSet[l]
                                << nl << " Set: " << cvxSet
                                << abort(FatalError);
                        }
                    }
                }

                label curFaceSign = 0;
                bool foundInternalFace = false;

                // Check all other points in the set,
                // and decide if all points lie on one side.
                forAll(cvxSet, l)
                {
                    // Skip duplicates.
                    if (findIndex(tmpFace, l) > -1)
                    {
                        continue;
                    }

                    vector rfVec = (cvxSet[l] - cvxSet[i]);
                    scalar dotProd = (rfVec/mag(rfVec)) & n;

                    // Obtain the sign of this point.
                    label fSign = Foam::sign(dotProd);

                    // Update the current sign if necessary.
                    if (curFaceSign == 0)
                    {
                        curFaceSign = fSign;
                    }
                    else
                    if (curFaceSign != fSign)
                    {
                        // Interior face. Bail out.
                        foundInternalFace = true;
                        break;
                    }
                }

                if (!foundInternalFace)
                {
                    // Looks like we found a face on the boundary.
                    // Check its sign to ensure that it points outward.
                    if (curFaceSign == 1)
                    {
                        tmpFace = tmpFace.reverseFace();
                    }

                    // Ensure that the face wasn't checked in.
                    bool alreadyCheckedIn = false;

                    forAll(testFaces, faceI)
                    {
                        const face& checkFace = testFaces[faceI];

                        if (face::compare(tmpFace, checkFace))
                        {
                            alreadyCheckedIn = true;
                            break;
                        }
                    }

                    // Add this face to the list of faces.
                    if (!alreadyCheckedIn)
                    {
                        testFaces.append(tmpFace);
                    }
                }

                // Reset the face size.
                tmpFace.setSize(3, -1);
            }
        }
    }

    // Weed-out faces that are sub-sets of larger faces
    label nFaces = 0;
    faceList cellFaces(testFaces.size());

    forAll(testFaces, faceI)
    {
        bool subset = false;

        const face& checkFace = testFaces[faceI];

        forAll(testFaces, faceJ)
        {
            const face& testFace = testFaces[faceJ];

            if (testFace.size() > checkFace.size())
            {
                bool foundUniquePoint = false;

                forAll(checkFace, pI)
                {
                    if (findIndex(testFace, checkFace[pI]) == -1)
                    {
                        foundUniquePoint = true;
                        break;
                    }
                }

                if (!foundUniquePoint)
                {
                    subset = true;
                    break;
                }
            }
        }

        if (!subset)
        {
            cellFaces[nFaces++] = testFaces[faceI];
        }
    }

    // Set to actual size
    cellFaces.setSize(nFaces);

    // Find cell-centroid
    vector xC = average(cvxSet);

    // Calculate volume from all accumulated faces.
    forAll(cellFaces, faceI)
    {
        const face& checkFace = cellFaces[faceI];

        vector xF = checkFace.centre(cvxSet);
        vector Sf = checkFace.normal(cvxSet);

        cVol += Sf & (xF - xC);
    }

    cVol *= (1.0/3.0);

    if (debug > 4)
    {
        // Write out points for post-processing
        labelListList cpList(cvxSet.size(), labelList(1));

        forAll(cpList, i)
        {
            cpList[i][0] = i;
        }

        writeVTK
        (
            cvxSetName,
            cvxSet.size(),
            cvxSet.size(),
            cvxSet.size(),
            cvxSet,
            cpList,
            0
        );

        Info << " Convex set: " << cvxSetName
             << " cellFaces: " << cellFaces
             << " Volume: " << cVol
             << endl;
    }

    // Return the computed volume.
    return cVol;
}


// Obtain a list of possible parent cells from the old mesh.
labelList dynamicTopoFvMesh::cellParents
(
    const label newCellIndex,
    const scalar searchFactor,
    const labelList& mappingCells
) const
{
    labelHashSet masterCells, finalCells;

    // Fetch connectivity from the old mesh.
    const cellList& cells = polyMesh::cells();
    const labelList& owner = polyMesh::faceOwner();
    const labelList& neighbour = polyMesh::faceNeighbour();

    forAll(mappingCells, cellI)
    {
        if (mappingCells[cellI] < 0)
        {
            continue;
        }

        if (mappingCells[cellI] < nOldCells_)
        {
            masterCells.set(mappingCells[cellI], empty());
        }
        else
        if (cellParents_.found(mappingCells[cellI]))
        {
            const labelList& nParents = cellParents_[mappingCells[cellI]];

            forAll(nParents, cI)
            {
                masterCells.set(nParents[cI], empty());
            }
        }
    }

    for (label attempt = 0; attempt < 5; attempt++)
    {
        // Fetch the initial set of candidates
        labelList initList = masterCells.toc();

        // Accumulate a larger stencil of cell neighbours
        forAll(initList, indexI)
        {
            label cellIndex = initList[indexI];

            const cell& cellToCheck = cells[cellIndex];

            forAll(cellToCheck, faceI)
            {
                // Add owner to the list.
                masterCells.set(owner[cellToCheck[faceI]], empty());

                if (cellToCheck[faceI] < neighbour.size())
                {
                    // Add the neighbour.
                    masterCells.set(neighbour[cellToCheck[faceI]], empty());
                }
            }
        }
    }

    // Fetch the new cell, and determine its bounds.
    const cell& newCell = cells_[newCellIndex];

    DynamicList<label> pointLabels(6);

    forAll(newCell, faceI)
    {
        const face& faceToCheck = faces_[newCell[faceI]];

        forAll(faceToCheck, pointI)
        {
            if (findIndex(pointLabels, faceToCheck[pointI]) == -1)
            {
                pointLabels.append(faceToCheck[pointI]);
            }
        }
    }

    pointField cellPoints(pointLabels.size());

    forAll(cellPoints, pointI)
    {
        cellPoints[pointI] = oldPoints_[pointLabels[pointI]];
    }

    // Prepare an axis-aligned bounding box around the cell,
    // and add cells according to cell-centre positions.
    boundBox cellBox(cellPoints, false);

    // Define a search radius
    vector bC = cellBox.midpoint();
    vector bMax = searchFactor * (cellBox.max() - bC);

    const vectorField& cellCentres = polyMesh::cellCentres();

    forAllIter(labelHashSet, masterCells, mIter)
    {
        vector xC = (cellCentres[mIter.key()] - bC);

        if ((xC & xC) < (bMax & bMax))
        {
            finalCells.insert(mIter.key());
        }
    }

    return finalCells.toc();
}


// Set fill-in mapping information for a particular cell
void dynamicTopoFvMesh::setCellMapping
(
    const label cIndex,
    const labelList& mapCells,
    const scalarField& mapWeights
)
{
    if (debug > 3)
    {
        Info << "Inserting mapping cell: " << cIndex << nl
             << " mapCells: " << mapCells << nl
             << " cellWeights: " << mapWeights
             << endl;
    }

    // Ensure compatible sizes.
    if (mapCells.size() != mapWeights.size())
    {
        FatalErrorIn("dynamicTopoFvMesh::setCellMapping()")
            << nl << " Incompatible mapping for cell: "
            << cIndex << ":: " << cells_[cIndex] << nl
            << " mapCells: " << mapCells
            << " cellWeights: " << mapWeights
            << abort(FatalError);
    }

    // Insert weights into the list, and overwrite if necessary
    bool index = -1;

    forAll(cellsFromCells_, indexI)
    {
        if (cellsFromCells_[indexI].index() == cIndex)
        {
            index = indexI;
            break;
        }
    }

    if (index == -1)
    {
        sizeUpList(objectMap(cIndex, mapCells), cellsFromCells_);
    }
    else
    {
        cellsFromCells_[index].masterObjects() = mapCells;
    }

    cellWeights_.set(cIndex, mapWeights);
}


// Set fill-in mapping information for a particular face
void dynamicTopoFvMesh::setFaceMapping
(
    const label fIndex,
    const labelList& mapFaces,
    const scalarField& mapWeights
)
{
    label patch = whichPatch(fIndex);

    if (debug > 3)
    {
        Info << "Inserting mapping face: " << fIndex << nl
             << " patch: " << patch << nl
             << " mapFaces: " << mapFaces << nl
             << " faceWeights: " << mapWeights
             << endl;
    }

    // Ensure compatible sizes.
    if (mapFaces.size() != mapWeights.size())
    {
        FatalErrorIn("dynamicTopoFvMesh::setFaceMapping()")
            << nl << " Incompatible mapping: " << nl
            << " mapFaces: " << mapFaces
            << " faceWeights: " << mapWeights
            << abort(FatalError);
    }

    // Insert weights into the list, and overwrite if necessary
    bool index = -1;

    forAll(facesFromFaces_, indexI)
    {
        if (facesFromFaces_[indexI].index() == fIndex)
        {
            index = indexI;
            break;
        }
    }

    if (index == -1)
    {
        sizeUpList(objectMap(fIndex, mapFaces), facesFromFaces_);
    }
    else
    {
        facesFromFaces_[index].masterObjects() = mapFaces;
    }

    faceWeights_.set(fIndex, mapWeights);
}


// Insert the specified cell to the mesh.
label dynamicTopoFvMesh::insertCell
(
    const cell& newCell,
    const scalar lengthScale,
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

    // Check if this cell was added to a zone
    if (addedCellZones_.found(cIndex))
    {
        addedCellZones_.erase(cIndex);
    }

    // Check if the cell was added in the current morph, and delete
    forAll(cellsFromPoints_, indexI)
    {
        if (cellsFromPoints_[indexI].index() == cIndex)
        {
            // Remove entry from the list
            removeIndex(indexI, cellsFromPoints_);
            break;
        }
    }

    forAll(cellsFromEdges_, indexI)
    {
        if (cellsFromEdges_[indexI].index() == cIndex)
        {
            // Remove entry from the list
            removeIndex(indexI, cellsFromEdges_);
            break;
        }
    }

    forAll(cellsFromFaces_, indexI)
    {
        if (cellsFromFaces_[indexI].index() == cIndex)
        {
            // Remove entry from the list
            removeIndex(indexI, cellsFromFaces_);
            break;
        }
    }

    forAll(cellsFromCells_, indexI)
    {
        if (cellsFromCells_[indexI].index() == cIndex)
        {
            // Remove entry from the list
            removeIndex(indexI, cellsFromCells_);
            break;
        }
    }

    // Check if any explicit cell weights were specified
    if (cellWeights_.found(cIndex))
    {
        cellWeights_.erase(cIndex);
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
             << newFace
             << " Owner: " << newOwner
             << " Neighbour: " << newNeighbour;

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
        // Increment the number of internal faces,
        // and subsequent patch-starts
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
        forAll(patchCoupling_, patchI)
        {
            if (patchCoupling_(patchI))
            {
                if (patchI == patch)
                {
                    if (patchCoupling_[patchI].masterFaceZone() > -1)
                    {
                        addedFaceZones_.insert
                        (
                            newFaceIndex,
                            patchCoupling_[patchI].masterFaceZone()
                        );
                    }

                    break;
                }

                if (patchCoupling_[patchI].patchMap().slaveIndex() == patch)
                {
                    if (patchCoupling_[patchI].slaveFaceZone() > -1)
                    {
                        addedFaceZones_.insert
                        (
                            newFaceIndex,
                            patchCoupling_[patchI].slaveFaceZone()
                        );
                    }

                    break;
                }
            }
        }
    }

    // Increment the total face count
    nFaces_++;

    return newFaceIndex;
}


// Remove the specified face from the mesh
void dynamicTopoFvMesh::removeFace
(
    const label fIndex
)
{
    if (debug > 2)
    {
        Info << "Removed face: "
             << fIndex << ": "
             << faces_[fIndex] << endl;
    }

    // Identify the patch for this face
    label patch = whichPatch(fIndex);

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
    faces_[fIndex].clear();
    owner_[fIndex] = -1;
    neighbour_[fIndex] = -1;
    faceEdges_[fIndex].clear();

    if (twoDMesh_)
    {
        // Remove from the stack as well
        forAll(faceStack_, stackI)
        {
            faceStack_[stackI].remove(fIndex);
        }
    }

    // Update coupled face maps, if necessary.
    forAll(patchCoupling_, patchI)
    {
        if (patchCoupling_(patchI))
        {
            const coupleMap& cMap = patchCoupling_[patchI].patchMap();

            cMap.removeMasterIndex(coupleMap::FACE, fIndex);
            cMap.removeSlaveIndex(coupleMap::FACE, fIndex);
        }
    }

    // Update the reverse face-map, but only if this is a face that existed
    // at time [n]. Added faces which are deleted during the topology change
    // needn't be updated.
    if (fIndex < nOldFaces_)
    {
        reverseFaceMap_[fIndex] = -1;
    }
    else
    {
        // Store this information for the reOrdering stage
        deletedFaces_.insert(fIndex);
    }

    // Check and remove from the list of added face patches
    if (addedFacePatches_.found(fIndex))
    {
        addedFacePatches_.erase(fIndex);
    }

    // Check if this face was added to a zone
    if (addedFaceZones_.found(fIndex))
    {
        addedFaceZones_.erase(fIndex);
    }

    // Check if the face was added in the current morph, and delete
    forAll(facesFromPoints_, indexI)
    {
        if (facesFromPoints_[indexI].index() == fIndex)
        {
            // Remove entry from the list
            removeIndex(indexI, facesFromPoints_);
            break;
        }
    }

    forAll(facesFromEdges_, indexI)
    {
        if (facesFromEdges_[indexI].index() == fIndex)
        {
            // Remove entry from the list
            removeIndex(indexI, facesFromEdges_);
            break;
        }
    }

    forAll(facesFromFaces_, indexI)
    {
        if (facesFromFaces_[indexI].index() == fIndex)
        {
            // Remove entry from the list
            removeIndex(indexI, facesFromFaces_);
            break;
        }
    }

    // Remove from the flipFaces list, if necessary
    if (flipFaces_.found(fIndex))
    {
        flipFaces_.erase(fIndex);
    }

    // Check if any explicit face weights were specified
    if (faceWeights_.found(fIndex))
    {
        faceWeights_.erase(fIndex);
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

        if (!twoDMesh_)
        {
            if (findIndex(edgePoints, -1) != -1)
            {
                FatalErrorIn("dynamicTopoFvMesh::insertEdge()")
                    << " EdgePoints is incorrectly specified." << nl
                    << " edgePoints: " << edgePoints << nl
                    << abort(FatalError);
            }
        }
    }

    // Keep track of added edges in a separate hash-table
    // This information will be required at the reordering stage
    addedEdgePatches_.insert(newEdgeIndex,patch);

    if (patch >= 0)
    {
        // Modify patch information for this boundary edge
        edgePatchSizes_[patch]++;

        for (label i = patch + 1; i < numPatches_; i++)
        {
            edgePatchStarts_[i]++;
        }
    }
    else
    {
        // Increment the number of internal edges, and subsequent patch-starts
        nInternalEdges_++;

        for (label i = 0; i < numPatches_; i++)
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
    const label eIndex
)
{
    if (debug > 2)
    {
        Info << "Removing edge: "
             << eIndex << ": "
             << edges_[eIndex] << endl;
    }

    if (!twoDMesh_)
    {
        // Remove the edgePoints entry
        edgePoints_[eIndex].clear();

        // Size-down the pointEdges list
        if (pointEdges_[edges_[eIndex][0]].size())
        {
            sizeDownList(eIndex, pointEdges_[edges_[eIndex][0]]);
        }

        if (pointEdges_[edges_[eIndex][1]].size())
        {
            sizeDownList(eIndex, pointEdges_[edges_[eIndex][1]]);
        }

        // Remove from the stack as well
        forAll(edgeStack_, stackI)
        {
            edgeStack(stackI).remove(eIndex);
        }

        // Update coupled face maps, if necessary.
        forAll(patchCoupling_, patchI)
        {
            if (patchCoupling_(patchI))
            {
                const coupleMap& cMap = patchCoupling_[patchI].patchMap();

                cMap.removeMasterIndex(coupleMap::EDGE, eIndex);
                cMap.removeSlaveIndex(coupleMap::EDGE, eIndex);
            }
        }
    }

    edges_[eIndex] = edge(-1, -1);
    edgeFaces_[eIndex].clear();

    // Identify the patch for this edge
    label patch = whichEdgePatch(eIndex);

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
    if (eIndex < nOldEdges_)
    {
        reverseEdgeMap_[eIndex] = -1;
    }
    else
    {
        // Store this information for the reOrdering stage
        deletedEdges_.insert(eIndex);
    }

    // Check and remove from the list of added edge patches
    if (addedEdgePatches_.found(eIndex))
    {
        addedEdgePatches_.erase(eIndex);
    }

    // Decrement the total edge-count
    nEdges_--;
}


// Insert the specified point to the mesh
label dynamicTopoFvMesh::insertPoint
(
    const point& newPoint,
    const point& oldPoint,
    const labelList& mappingPoints,
    const label zoneID
)
{
    // Add a new point to the end of the list
    label newPointIndex = points_.size();

    points_.append(newPoint);
    oldPoints_.append(oldPoint);

    if (debug > 2)
    {
        Info << "Inserting new point: "
             << newPointIndex << ": "
             << newPoint
             << " and old point: "
             << oldPoint
             << "  Mapped from: "
             << mappingPoints << endl;
    }

    // Add an empty entry to pointEdges as well.
    // This entry can be sized-up appropriately at a later stage.
    if (!twoDMesh_)
    {
        pointEdges_.append(labelList(0));
    }

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
    const label pIndex
)
{
    if (debug > 2)
    {
        Info << "Removing point: "
             << pIndex << ": "
             << points_[pIndex] << endl;
    }

    // Remove the point
    // (or just make sure that it's never used anywhere else)
    // points_[pIndex] = point();
    // oldPoints_[pIndex] = point();

    // Remove pointEdges as well
    if (!twoDMesh_)
    {
        pointEdges_[pIndex].clear();
    }

    // Update the reverse point map
    if (pIndex < nOldPoints_)
    {
        reversePointMap_[pIndex] = -1;
    }
    else
    {
        deletedPoints_.insert(pIndex);
    }

    // Check if this point was added to a zone
    if (addedPointZones_.found(pIndex))
    {
        addedPointZones_.erase(pIndex);
    }

    // Update coupled point maps, if necessary.
    forAll(patchCoupling_, patchI)
    {
        if (patchCoupling_(patchI))
        {
            // Obtain references
            const coupleMap & cMap = patchCoupling_[patchI].patchMap();

            Map<label>& pointMap = cMap.entityMap(coupleMap::POINT);
            Map<label>& rPointMap = cMap.reverseEntityMap(coupleMap::POINT);

            if (pointMap.found(pIndex))
            {
                // Erase the reverse map first
                rPointMap.erase(pointMap[pIndex]);

                // Update pointMap
                pointMap.erase(pIndex);
            }
        }
    }

    // Check if the point was added in the current morph, and delete
    forAll(pointsFromPoints_, indexI)
    {
        if (pointsFromPoints_[indexI].index() == pIndex)
        {
            // Remove entry from the list
            removeIndex(indexI, pointsFromPoints_);
            break;
        }
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
                        if (triFaceCompare(cFace, oFace))
                        {
                            ringEntities[1][indexI] = currCell[faceI];
                        }

                        // Build a comparison face
                        oFace[0] = edgeToCheck[1];
                        oFace[1] = nextHullPoint;
                        oFace[2] = otherPoint;

                        // Check if this face contains edgeToCheck[1]
                        if (triFaceCompare(cFace, oFace))
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
            Info << "edgeFaces: " << endl;

            forAll(eFaces, faceI)
            {
                Info << " Face: " << eFaces[faceI]
                     << ":: " << faces_[eFaces[faceI]]
                     << " Owner: " << owner_[eFaces[faceI]]
                     << " Neighbour: " << neighbour_[eFaces[faceI]]
                     << endl;
            }

            writeVTK("vRingEdgeFaces", eFaces, 2);
            writeVTK("vRingCellToCheck", cellIndex);

            // Something's terribly wrong
            FatalErrorIn
            (
                "void dynamicTopoFvMesh::buildEdgePoints"
            )
                << " Failed to determine a vertex ring. " << nl
                << " edgeFaces connectivity is inconsistent. " << nl
                << " Edge: " << eIndex << ":: " << edgeToCheck << nl
                << " edgeFaces: " << eFaces << nl
                << " Patch: " << whichEdgePatch(eIndex) << nl
                << " cellIndex: " << cellIndex
                << " :: " << cellToCheck << nl
                << " Current edgePoints: " << ePoints
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

        label pos = mag(((k*nD*nD)+(j*nD)+i) % binSize);

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
        if
        (
            (proximityPatches_.found(boundary[patchI].name())) ||
            (twoDMesh_ && boundary[patchI].type() == "symmetryPlane")
        )
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

        if (twoDMesh_)
        {
            forAll(pairToCheck, indexI)
            {
                if (deletedFaces_.found(pairToCheck[indexI]))
                {
                    available = false;
                    break;
                }

                if (pairToCheck[indexI] < nOldFaces_)
                {
                    if (reverseFaceMap_[pairToCheck[indexI]] == -1)
                    {
                        available = false;
                        break;
                    }
                }
            }
        }
        else
        {
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
        }

        if (available)
        {
            // Slice the mesh at this point.
            sliceMesh(pairToCheck);
        }
    }

    Info << "Done." << endl;

    // Clear out data.
    sliceBoxes_.clear();
    slicePairs_.clear();

    // Set the sliceHoldOff value
    sliceHoldOff_ = 50;
}


// Test an edge / face for proximity with other faces on proximity patches
// and return the scalar distance to an oppositely-oriented face.
scalar dynamicTopoFvMesh::testProximity
(
    const label index,
    label& proximityFace
) const
{
    DynamicList<label> posIndices(20);
    scalar minDistance = GREAT, minDeviation = -0.9, testStep = 0.0;
    label nD = spatialRes_, binSize = proximityBins_.size();
    vector gCentre = vector::zero, gNorm = vector::zero;

    // Obtain min/max extents
    const point& bMin = proxBoundBox_.min();
    const point& bMax = proxBoundBox_.max();

    // Extend bounding-box dimensions a bit to avoid edge-effects.
    scalar ext = 0.02*(mag(bMax - bMin));

    // Define an inverse grid-cell size.
    scalar xL = nD/(bMax.x() - bMin.x() + ext);
    scalar yL = nD/(bMax.y() - bMin.y() + ext);
    scalar zL = nD/(bMax.z() - bMin.z() + ext);

    // Reset the proximity face
    proximityFace = -1;

    if (twoDMesh_)
    {
        // Obtain the face-normal.
        gNorm = quadFaceNormal(faces_[index]);

        gNorm /= (mag(gNorm) + VSMALL);

        // Obtain the face centre.
        gCentre = quadFaceCentre(index);

        // Calculate a test step-size
        testStep = edgeLength(getTriBoundaryEdge(index));
    }
    else
    {
        const edge& thisEdge = edges_[index];
        const labelList& eFaces = edgeFaces_[index];

        // Obtain the edge centre.
        gCentre = (0.5 * (points_[thisEdge[0]] + points_[thisEdge[1]]));

        // Obtain the edge-normal
        forAll(eFaces, faceI)
        {
            if (neighbour_[eFaces[faceI]] == -1)
            {
                // Obtain the normal.
                gNorm += triFaceNormal(faces_[eFaces[faceI]]);
            }
        }

        gNorm /= (mag(gNorm) + VSMALL);

        // Calculate a test step-size
        testStep = edgeLength(index);
    }

    // Now take multiple steps in both normal directions,
    // and add to the list of boxes to be checked.
    for (scalar dir = -1.0; dir < 2.0; dir += 2.0)
    {
        for (scalar step = 0.0; step < 5.0*testStep; step += testStep)
        {
            // Hash the point-location
            point p = (gCentre + (dir*step*gNorm)) - bMin;

            label i = label(mag(::floor(p.x()*xL)));
            label j = label(mag(::floor(p.y()*yL)));
            label k = label(mag(::floor(p.z()*zL)));

            label pos = mag(((k*nD*nD)+(j*nD)+i) % binSize);

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
            vector rFace = (faceCentres[posBin[faceI]] - gCentre);

            scalar distance = mag(rFace);

            // Step 2: Check if this face is oriented away from face / edge.
            const vector& fNorm = faceAreas[posBin[faceI]];

            scalar deviation = (gNorm & (fNorm/(mag(fNorm) + VSMALL)));

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
                // domain by comparing with the normal.
                if ((gNorm & (rFace/distance)) > 0.0)
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
            labelPair proxPoints(-1, -1);
            bool foundPoint = false;

            if (twoDMesh_)
            {
                proxPoints.first() = index;
                proxPoints.second() = proximityFace;

                if
                (
                    (faces_[index].size() == 4) &&
                    (polyMesh::faces()[proximityFace].size() == 4)
                )
                {
                    foundPoint = true;
                }
            }
            else
            {
                const edge& thisEdge = edges_[index];
                const face& proxFace = polyMesh::faces()[proximityFace];

                // Check if any points on this face are still around.
                // If yes, mark one of them as the end point
                // for Dijkstra's algorithm. The start point will be a point
                // on this edge.
                proxPoints.first() = thisEdge[0];

                forAll(proxFace, pointI)
                {
                    if (reversePointMap_[proxFace[pointI]] != -1)
                    {
                        proxPoints.second() = proxFace[pointI];
                        foundPoint = true;
                        break;
                    }
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
void dynamicTopoFvMesh::calculateLengthScale(bool dump)
{
    if (!edgeRefinement_)
    {
        return;
    }

    Switch dumpLengthScale(false);

    if
    (
        dict_.subDict("dynamicTopoFvMesh").found("dumpLengthScale") ||
        mandatory_
    )
    {
        dumpLengthScale =
        (
            dict_.subDict("dynamicTopoFvMesh").lookup("dumpLengthScale")
        );
    }

    volScalarField *lsPtr(NULL);

    if (dumpLengthScale && time().outputTime() && dump)
    {
        lsPtr = new volScalarField
        (
            IOobject
            (
                "lengthScale",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            dimensionedScalar("scalar", dimLength, 0)
        );
    }

    // Bail-out if a dumping was not requested in dictionary.
    if (dump && !dumpLengthScale && !time().outputTime())
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

    // Exchange length-scale buffers across processors.
    if (Pstream::parRun())
    {
        exchangeLengthBuffers();
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

    // Check if length-scale is to be dumped to disk.
    if (dumpLengthScale && time().outputTime() && dump)
    {
        scalarField& lengthScale = lsPtr->internalField();

        // Obtain length-scale values from the mesh
        forAll(lengthScale_, cellI)
        {
            lengthScale[cellI] = lengthScale_[cellI];
        }

        lsPtr->write();
    }

    // Wait for transfers before continuing.
    waitForBuffers();
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

    if
    (
        dict_.subDict("dynamicTopoFvMesh").found("bandWidthReduction") ||
        mandatory_
    )
    {
        bandWidthReduction_ =
        (
            dict_.subDict
            ("dynamicTopoFvMesh").lookup("bandWidthReduction")
        );
    }
    else
    {
        bandWidthReduction_ = false;
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

    // Read growthFactor from dictionary
    growthFactor_ = readScalar(edgeOptionDict.lookup("growthFactor"));

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
            // sliceThreshold_ = Foam::max(sliceThreshold_, minLengthScale_);
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


// Load the mesh-quality metric from the library
void dynamicTopoFvMesh::loadMetric()
{
    if (twoDMesh_)
    {
        return;
    }

    // Specify the dictionary we would be looking in...
    const dictionary& meshDict = dict_.subDict("dynamicTopoFvMesh");

    // Select an appropriate metric
    tetMetric_ = tetMetric::New(meshDict, meshDict.lookup("tetMetric"));
}


// Load the mesh-motion solver
void dynamicTopoFvMesh::loadMotionSolver()
{
    if (mPtr_.valid())
    {
        FatalErrorIn
        (
            "dynamicTopoFvMesh::loadMotionSolver() "
        ) << nl << " Motion solver already loaded. "
          << abort(FatalError);
    }
    else
    if (dict_.found("solver"))
    {
        mPtr_ = motionSolver::New(*this);
    }
}


// Load the field mapper
void dynamicTopoFvMesh::loadFieldMapper()
{
    if (mapper_.valid())
    {
        FatalErrorIn
        (
            "dynamicTopoFvMesh::loadFieldMapper() "
        ) << nl << " Field mapper already loaded. "
          << abort(FatalError);
    }
    else
    {
        mapper_.set(new topoMapper(*this));
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
    }
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


// Method for the swapping of an edge in 3D
//  - To be used mainly for testing purposes only.
//  - Use swap3DEdges on the entire mesh for efficiency.
void dynamicTopoFvMesh::swapEdge
(
    const label eIndex,
    bool forceOp
)
{
    // Dynamic programming variables
    labelList m;
    PtrList<scalarListList> Q;
    PtrList<labelListList> K, triangulations;

    // Allocate dynamic programming tables
    initTables(m, Q, K, triangulations);

    // Compute the minimum quality of cells around this edge
    scalar minQuality = computeMinQuality(eIndex);

    // Check if this edge is on a bounding curve
    if (checkBoundingCurve(eIndex))
    {
        FatalErrorIn("dynamicTopoFvMesh::swapEdge(const label eIndex)")
            << nl << " Cannot swap edges on bounding curves. "
            << abort(FatalError);
    }

    // Fill the dynamic programming tables
    if (fillTables(eIndex, minQuality, m, Q, K, triangulations))
    {
        // Check if edge-swapping is required.
        scalar newQuality = Q[0][0][m[0]-1];

        if (newQuality > minQuality)
        {
            // Remove this edge according to the swap sequence
            removeEdgeFlips(eIndex, minQuality, K, triangulations);
        }
        else
        if (forceOp)
        {
            if (newQuality < 0.0)
            {
                FatalErrorIn("dynamicTopoFvMesh::swapEdge(const label eIndex)")
                    << " Forcing swap on edge: " << eIndex
                    << ":: " << edges_[eIndex]
                    << " will yield an invalid cell quality: "
                    << newQuality << " Old Quality: " << minQuality
                    << abort(FatalError);
            }
            else
            {
                removeEdgeFlips(eIndex, minQuality, K, triangulations);
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


// Utility method to compute the quality of a vertex hull
// around an edge after bisection.
scalar dynamicTopoFvMesh::computeBisectionQuality
(
    const label eIndex
) const
{
    scalar minQuality = GREAT, minVolume = GREAT;
    scalar cQuality = 0.0, oldVolume = 0.0;

    // Obtain a reference to this edge and corresponding edgePoints
    const edge& edgeToCheck = edges_[eIndex];
    const labelList& hullVertices = edgePoints_[eIndex];

    // Obtain point references
    const point& a = points_[edgeToCheck[0]];
    const point& c = points_[edgeToCheck[1]];

    const point& aOld = oldPoints_[edgeToCheck[0]];
    const point& cOld = oldPoints_[edgeToCheck[1]];

    // Compute the mid-point of the edge
    point midPoint = 0.5*(a + c);
    point oldPoint = 0.5*(aOld + cOld);

    if (whichEdgePatch(eIndex) < 0)
    {
        // Internal edge.
        forAll(hullVertices, indexI)
        {
            label prevIndex = hullVertices.rcIndex(indexI);

            // Pick vertices off the list
            const point& b = points_[hullVertices[prevIndex]];
            const point& d = points_[hullVertices[indexI]];

            const point& bOld = oldPoints_[hullVertices[prevIndex]];
            const point& dOld = oldPoints_[hullVertices[indexI]];

            // Compute the quality of the upper half.
            cQuality = tetMetric_(a, b, midPoint, d);

            // Compute old volume of the upper half.
            oldVolume = tetVolume(aOld, bOld, oldPoint, dOld);

            // Check if the volume / quality is worse
            minQuality = Foam::min(cQuality, minQuality);
            minVolume = Foam::min(oldVolume, minVolume);

            // Compute the quality of the lower half.
            cQuality = tetMetric_(midPoint, b, c, d);

            // Compute old volume of the lower half.
            oldVolume = tetVolume(oldPoint, bOld, cOld, dOld);

            // Check if the quality is worse
            minQuality = Foam::min(cQuality, minQuality);
            minVolume = Foam::min(oldVolume, minVolume);
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

            const point& bOld = oldPoints_[hullVertices[indexI-1]];
            const point& dOld = oldPoints_[hullVertices[indexI]];

            // Compute the quality of the upper half.
            cQuality = tetMetric_(a, b, midPoint, d);

            // Compute old volume of the upper half.
            oldVolume = tetVolume(aOld, bOld, oldPoint, dOld);

            // Check if the quality is worse
            minQuality = Foam::min(cQuality, minQuality);
            minVolume = Foam::min(oldVolume, minVolume);

            // Compute the quality of the lower half.
            cQuality = tetMetric_(midPoint, b, c, d);

            // Compute old volume of the lower half.
            oldVolume = tetVolume(oldPoint, bOld, cOld, dOld);

            // Check if the quality is worse
            minQuality = Foam::min(cQuality, minQuality);
            minVolume = Foam::min(oldVolume, minVolume);
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

    // If a negative old-volume was encountered,
    // return an invalid quality.
    if (minVolume < 0.0)
    {
        return minVolume;
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

    point midPoint = triFaceCentre(fIndex);

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
        cQuality = tetMetric_(midPoint, b, c, d);

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
            cQuality = tetMetric_(midPoint, b, c, d);

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

    label triEdge = getTriBoundaryEdge(fIndex);
    label triFace = getTriBoundaryFace(fIndex);

    // Measure the boundary edge-length of the face in question
    scalar length = edgeLength(triEdge);

    // Determine the boundary triangular face area
    scalar area = triFaceArea(faces_[triFace]);

    // This cell has to be removed...
    if (mag(area) < (0.2*length*length))
    {
        if (self() == 0)
        {
            if (debug > 1)
            {
                InfoIn("dynamicTopoFvMesh::remove2DSliver() ")
                    << nl
                    << " Considering face: " << fIndex
                    << ":: " << faces_[fIndex]
                    << " for sliver removal."
                    << endl;
            }

            // Find the isolated point.
            label ptIndex = -1, nextPtIndex = -1;

            const edge& edgeToCheck = edges_[triEdge];

            findIsolatedPoint
            (
                faces_[triFace],
                edgeToCheck,
                ptIndex,
                nextPtIndex
            );

            // Find the prism faces
            label fOwner = owner_[fIndex];
            FixedList<label,2> c0BdyIndex, c0IntIndex;
            FixedList<face,2>  c0BdyFace,  c0IntFace;

            findPrismFaces
            (
                fIndex,
                fOwner,
                c0BdyFace,
                c0BdyIndex,
                c0IntFace,
                c0IntIndex
            );

            // Determine the interior faces connected to each edge-point.
            label firstFace = -1, secondFace = -1;

            if (c0IntFace[0].which(edgeToCheck[0]) > -1)
            {
                firstFace  = c0IntIndex[0];
                secondFace = c0IntIndex[1];
            }
            else
            {
                firstFace  = c0IntIndex[1];
                secondFace = c0IntIndex[0];
            }

            point ec =
            (
                0.5 * (points_[edgeToCheck[0]] + points_[edgeToCheck[1]])
            );

            FixedList<vector, 2> p(vector::zero), q(vector::zero);
            FixedList<scalar, 2> proj(0.0);

            // Find the projection on the edge.
            forAll(edgeToCheck, pointI)
            {
                p[pointI] = (points_[ptIndex] - points_[edgeToCheck[pointI]]);
                q[pointI] = (ec - points_[edgeToCheck[pointI]]);

                q[pointI] /= (mag(q[pointI]) + VSMALL);

                proj[pointI] = (p[pointI] & q[pointI]);
            }

            // Take action based on the magnitude of the projection.
            if (mag(proj[0]) < VSMALL)
            {
                collapseQuadFace(firstFace);
                return;
            }

            if (mag(proj[1]) < VSMALL)
            {
                collapseQuadFace(secondFace);
                return;
            }

            if (proj[0] > 0.0 && proj[1] < 0.0)
            {
                changeMap map = bisectQuadFace(firstFace);

                // Loop through added faces, and collapse
                // the appropriate one
                const List<FixedList<label,2> >& aF = map.addedFaceList();

                forAll(aF, faceI)
                {
                    if
                    (
                        (owner_[aF[faceI][0]] == fOwner) &&
                        (aF[faceI][0] != firstFace)
                    )
                    {
                        collapseQuadFace(aF[faceI][0]);
                        break;
                    }
                }

                return;
            }

            if (proj[0] < 0.0 && proj[1] > 0.0)
            {
                changeMap map = bisectQuadFace(secondFace);

                // Loop through added faces, and collapse
                // the appropriate one
                const List<FixedList<label,2> >& aF = map.addedFaceList();

                forAll(aF, faceI)
                {
                    if
                    (
                        (owner_[aF[faceI][0]] == fOwner) &&
                        (aF[faceI][0] != secondFace)
                    )
                    {
                        collapseQuadFace(aF[faceI][0]);
                        break;
                    }
                }

                return;
            }

            if (proj[0] > 0.0 && proj[1] > 0.0)
            {
                changeMap map = bisectQuadFace(fIndex);

                // Loop through added faces, and collapse
                // the appropriate one
                const List<FixedList<label,2> >& aF = map.addedFaceList();

                forAll(aF, faceI)
                {
                    if
                    (
                        (owner_[aF[faceI][0]] == fOwner) &&
                        (aF[faceI][0] != fIndex)
                    )
                    {
                        collapseQuadFace(aF[faceI][0]);
                        break;
                    }
                }

                return;
            }
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

            if (triFaceCompare(thisFace, faceToCheck))
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

    // Check if a removeSlivers entry was found in the dictionary
    if (dict_.subDict("dynamicTopoFvMesh").found("removeSlivers"))
    {
        Switch rs =
        (
            dict_.subDict("dynamicTopoFvMesh").lookup("removeSlivers")
        );

        if (!rs)
        {
            return;
        }
    }

    // If coupled patches exist, set the flag
    if (patchCoupling_.size() || procIndices_.size())
    {
        setCoupledModification();
    }

    // Sort by sliver-quality.
    labelList cIndices(thresholdSlivers_.toc());
    SortableList<scalar> values(cIndices.size());

    // Fill-in values to sort by...
    forAll(cIndices, indexI)
    {
        values[indexI] = thresholdSlivers_[cIndices[indexI]];
    }

    // Explicitly sort by quality value.
    values.sort();

    const labelList& indices = values.indices();

    if (debug && thresholdSlivers_.size())
    {
        Info << "Sliver list: " << endl;

        forAll(indices, indexI)
        {
            label cIndex = cIndices[indices[indexI]];

            Info << " Cell: " << cIndex
                 << " Quality: " << thresholdSlivers_[cIndex]
                 << endl;
        }

        if (debug > 1)
        {
            writeVTK("sliverCells", cIndices, 3);
        }
    }

    forAll(indices, indexI)
    {
        // Fetch the cell index
        label cIndex = cIndices[indices[indexI]];

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
                        sendPatchMeshes_[proc].patchMap().reverseEntityMap
                        (
                            coupleMap::CELL
                        )
                    );

                    if (rCellMap.found(cIndex))
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

        // Identify the sliver type.
        changeMap map = identifySliverType(cIndex);

        if (debug)
        {
            WarningIn("dynamicTopoFvMesh::removeSlivers()")
                << nl << "Removing Cell: " << cIndex
                << " of sliver type: " << map.type()
                << " with quality: " << thresholdSlivers_[cIndex]
                << endl;
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
            edge edgeToCheck
            (
                firstMap.addedPointList()[0][0],
                secondMap.addedPointList()[0][0]
            );

            bool foundCollapseEdge = false;

            const List<FixedList<label,2> >& firstMapEdges =
            (
                firstMap.addedEdgeList()
            );

            const List<FixedList<label,2> >& secondMapEdges =
            (
                secondMap.addedEdgeList()
            );

            // Loop through the first list.
            forAll(firstMapEdges, edgeI)
            {
                const edge& thisEdge = edges_[firstMapEdges[edgeI][0]];

                if (thisEdge == edgeToCheck)
                {
                    // Collapse this edge.
                    collapseEdge
                    (
                        firstMapEdges[edgeI][0],
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
                    const edge& thisEdge = edges_[secondMapEdges[edgeI][0]];

                    if (thisEdge == edgeToCheck)
                    {
                        // Collapse this edge.
                        collapseEdge
                        (
                            secondMapEdges[edgeI][0],
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
            edge edgeToCheck
            (
                map.apexPoint(),
                faceMap.addedPointList()[0][0]
            );

            const List<FixedList<label,2> >& faceMapEdges =
            (
                faceMap.addedEdgeList()
            );

            forAll(faceMapEdges, edgeI)
            {
                const edge& thisEdge = edges_[faceMapEdges[edgeI][0]];

                if (thisEdge == edgeToCheck)
                {
                    // Collapse this edge.
                    collapseEdge(faceMapEdges[edgeI][0], -1, false, true);

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
            edge edgeToCheck
            (
                map.apexPoint(),
                firstMap.addedPointList()[0][0]
            );

            const List<FixedList<label,2> >& firstMapEdges =
            (
                firstMap.addedEdgeList()
            );

            // Loop through the first list.
            forAll(firstMapEdges, edgeI)
            {
                const edge& thisEdge = edges_[firstMapEdges[edgeI][0]];

                if (thisEdge == edgeToCheck)
                {
                    // Collapse this edge.
                    collapseEdge(firstMapEdges[edgeI][0], -1, false, true);

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
        unsetCoupledModification();
    }
}


// MultiThreaded topology modifier [2D]
void dynamicTopoFvMesh::threadedTopoModifier2D()
{
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

        // Handle mesh slicing events, if necessary
        handleMeshSlicing();

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
    // Remove sliver cells first.
    removeSlivers();

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


// Map all fields in time using a customized mapper
void dynamicTopoFvMesh::mapFields(const mapPolyMesh& meshMap)
{
    if (debug)
    {
        Info << "void dynamicTopoFvMesh::mapFields(const mapPolyMesh&): "
             << "Mapping fv fields."
             << endl;
    }

    // Set the mapPolyMesh object in the mapper
    mapper_().setMapper(meshMap);

    // Set weighting information.
    // This takes over the weight data.
    mapper_().setFaceWeights(faceWeights_);
    mapper_().setCellWeights(cellWeights_);

    // Map all the volFields in the objectRegistry
    MapGeometricFields<scalar,fvPatchField,topoMapper,volMesh>
        (mapper_());
    MapGeometricFields<vector,fvPatchField,topoMapper,volMesh>
        (mapper_());
    MapGeometricFields<sphericalTensor,fvPatchField,topoMapper,volMesh>
        (mapper_());
    MapGeometricFields<symmTensor,fvPatchField,topoMapper,volMesh>
        (mapper_());
    MapGeometricFields<tensor,fvPatchField,topoMapper,volMesh>
        (mapper_());

    // Map all the surfaceFields in the objectRegistry
    MapGeometricFields<scalar,fvsPatchField,topoMapper,surfaceMesh>
        (mapper_());
    MapGeometricFields<vector,fvsPatchField,topoMapper,surfaceMesh>
        (mapper_());
    MapGeometricFields<sphericalTensor,fvsPatchField,topoMapper,surfaceMesh>
        (mapper_());
    MapGeometricFields<symmTensor,fvsPatchField,topoMapper,surfaceMesh>
        (mapper_());
    MapGeometricFields<tensor,fvsPatchField,topoMapper,surfaceMesh>
        (mapper_());

    // Clear mapper
    mapper_().clear();
}


// Update the mesh for topology changes.
// Return true if changes have occurred
bool dynamicTopoFvMesh::update()
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

    // Set old-points before moving the mesh
    const pointField& oldPoints = points();

    forAll(oldPoints, pointI)
    {
        oldPoints_[pointI] = oldPoints[pointI];
    }

    // Set old cell-centre information for the mapping stage
    mapper_().setOldCellCentres(fvMesh::C());

    // Invoke mesh-motion solver and move points
    if (mPtr_.valid())
    {
        movePoints(mPtr_->newPoints());
    }

    // Obtain the most recent point-positions
    const pointField& newPoints = points();

    forAll(newPoints, pointI)
    {
        points_[pointI] = newPoints[pointI];
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

    // Return if the interval is invalid or first time-step (no V0).
    // Handy while using only mesh-motion.
    if (interval_ < 0 || time().timeIndex() < 1)
    {
        // Dump procIDs to disk, if requested.
        writeProcIDs();

        return false;
    }

    // Return if re-meshing is not at interval,
    // or sliver cells are absent.
    if ((time().timeIndex() % interval_ != 0) && sliversAbsent)
    {
        // Dump procIDs to disk, if requested.
        writeProcIDs();

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

    meshQuality(true);

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
        pointField preMotionPoints(nPoints_);
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
            preMotionPoints,
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

        // Obtain the patch-point maps before resetting the mesh
        List<Map<label> > oldPatchPointMaps(numPatches_);

        forAll(oldPatchPointMaps, patchI)
        {
            oldPatchPointMaps[patchI] = boundaryMesh()[patchI].meshPointMap();
        }

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

        // Generate new mesh mapping information
        mapPolyMesh mpm
        (
            (*this),
            nOldPoints_,
            nOldFaces_,
            nOldCells_,
            pointMap_,
            pointsFromPoints_,
            faceMap_,
            facesFromPoints_,
            facesFromEdges_,
            facesFromFaces_,
            cellMap_,
            cellsFromPoints_,
            cellsFromEdges_,
            cellsFromFaces_,
            cellsFromCells_,
            reversePointMap_,
            reverseFaceMap_,
            reverseCellMap_,
            flipFaces_,
            patchPointMap,
            pointZoneMap,
            faceZonePointMap,
            faceZoneFaceMap,
            cellZoneMap,
            preMotionPoints,
            oldPatchStarts,
            oldPatchNMeshPoints,
            true
        );

        // Move points to positions before mesh-motion
        movePoints(mpm.preMotionPoints());

        // Update the underlying mesh, and map all related fields
        updateMesh(mpm);

        // Reset old-volumes / mesh-fluxes.
        // This overrides mapped V0 values.
        resetMotion();
        setV0();

        // Correct volume fluxes on the old mesh
        mapper_().correctFluxes();

        // Update the mesh-motion solver
        if (mPtr_.valid())
        {
            mPtr_->updateMesh(mpm);
        }

        // Now move back to new points and
        // compute correct mesh-fluxes.
        movePoints(points);

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

        // Now that all connectivity changes are successful,
        // update coupled maps (in a separate thread, if available).
        if (Pstream::parRun())
        {
            if (threader_->multiThreaded())
            {
                threader_->addToWorkQueue
                (
                    &initCoupledConnectivity,
                    reinterpret_cast<void *>(this)
                );
            }
            else
            {
                initCoupledConnectivity(this);
            }
        }

        // Clear unwanted member data
        addedFacePatches_.clear();
        addedEdgePatches_.clear();
        addedPointZones_.clear();
        addedFaceZones_.clear();
        addedCellZones_.clear();
        faceParents_.clear();
        cellParents_.clear();

        // Clear the deleted entity map
        deletedPoints_.clear();
        deletedEdges_.clear();
        deletedFaces_.clear();
        deletedCells_.clear();

        // Clear flipFaces
        flipFaces_.clear();

        // Set new sizes for the reverse maps
        reversePointMap_.setSize(nPoints_, -7);
        reverseEdgeMap_.setSize(nEdges_, -7);
        reverseFaceMap_.setSize(nFaces_, -7);
        reverseCellMap_.setSize(nCells_, -7);

        // Update "old" information
        nOldPoints_ = nPoints_;
        nOldEdges_ = nEdges_;
        nOldFaces_ = nFaces_;
        nOldCells_ = nCells_;
        nOldInternalFaces_ = nInternalFaces_;
        nOldInternalEdges_ = nInternalEdges_;

        for (label i=0; i<numPatches_; i++)
        {
            oldPatchSizes_[i] = patchSizes_[i];
            oldPatchStarts_[i] = patchStarts_[i];
            oldEdgePatchSizes_[i] = edgePatchSizes_[i];
            oldEdgePatchStarts_[i] = edgePatchStarts_[i];
            oldPatchNMeshPoints_[i] = patchNMeshPoints_[i];
        }

        // Basic checks for mesh-validity
        if (debug > 2)
        {
            checkMesh(true);
        }
    }

    // Dump length-scale to disk, if requested.
    calculateLengthScale(true);

    // Dump procIDs to disk, if requested.
    writeProcIDs();

    // Print out topo-stats
    Info << " Reordering time: " << reOrderingTimer.elapsedTime() << endl;
    Info << " nBisections: " << nBisections_[0];
    Info << " nSurfBisections: " << nBisections_[1] << endl;
    Info << " nCollapses: " << nCollapses_[0];
    Info << " nSurfCollapses: " << nCollapses_[1] << endl;
    Info << " nSwaps: " << nSwaps_[0];
    Info << " nSurfSwaps: " << nSwaps_[1] << endl;

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
