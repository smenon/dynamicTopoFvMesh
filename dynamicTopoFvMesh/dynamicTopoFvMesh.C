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

#include "IOmanip.H"
#include "triFace.H"
#include "clockTime.H"
#include "mapPolyMesh.H"
#include "volFields.H"
#include "motionSolver.H"
#include "MapFvFields.H"
#include "fvPatchFields.H"
#include "fvsPatchFields.H"
#include "MeshObject.H"
#include "topoMapper.H"
#include "SortableList.H"
#include "StaticHashTable.H"
#include "lengthScaleEstimator.H"

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
    motionSolver_(NULL),
    lengthEstimator_(NULL),
    eMeshPtr_(NULL),
    oldPoints_(polyMesh::points()),
    points_(polyMesh::points()),
    faces_(polyMesh::faces()),
    owner_(polyMesh::faceOwner()),
    neighbour_(polyMesh::faceNeighbour()),
    cells_(primitiveMesh::cells()),
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
    nModifications_(0),
    maxModifications_(-1),
    nBisections_(0),
    nCollapses_(0),
    nSwaps_(0),
    sliverThreshold_(0.1),
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
    for (label i = 0; i < numPatches_; i++)
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

    // Load the length-scale estimator,
    // and read refinement options
    loadLengthScaleEstimator();
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
    motionSolver_(NULL),
    lengthEstimator_(NULL),
    eMeshPtr_(NULL),
    oldPoints_(points),
    points_(points),
    faces_(faces),
    cells_(cells),
    edges_(edges),
    faceEdges_(faceEdges),
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
    nModifications_(0),
    maxModifications_(mesh.maxModifications_),
    nBisections_(0),
    nCollapses_(0),
    nSwaps_(0),
    sliverThreshold_(mesh.sliverThreshold_),
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
    edgeFaces_ = invertManyToMany<labelList, labelList>(nEdges_, faceEdges_);

    // Size-up edgePoints for now, but explicitly construct
    // for each edge later, based on point coupling.
    if (!twoDMesh_)
    {
        pointEdges_ = invertManyToMany<edge, labelList>(nPoints_, edges_);
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


// Obtain map weighting factors for a face
void dynamicTopoFvMesh::computeFaceWeights
(
    const label fIndex,
    const labelList& mapCandidates,
    labelList& parents,
    scalarField& weights,
    vectorField& centres
) const
{
    scalar searchFactor = 1.0;

    label nOldIntersects = -1, nIntersects = 0, nAttempts = 0;

    // Maintain a list of candidates and intersection points
    labelList candidates;

    // Output option for the convex set algorithm
    bool output = false;

    while (nAttempts < 10)
    {
        // Obtain candidate parents for this face
        candidates =
        (
            faceParents
            (
                fIndex,
                searchFactor,
                mapCandidates
            )
        );

        // Set sizes
        boolList intersects(candidates.size(), false);

        // Test for intersections
        forAll(candidates, indexI)
        {
            vectorField tP;

            intersects[indexI] =
            (
                faceIntersection
                (
                    fIndex,
                    candidates[indexI],
                    tP
                )
            );

            if (intersects[indexI])
            {
                nIntersects++;
            }
        }

        if ((nIntersects == nOldIntersects) && (nIntersects != 0))
        {
            if (debug > 3)
            {
                Info << " Face: " << fIndex << nl
                     << " nCandidates: " << candidates.size() << nl
                     << " nIntersects: " << nIntersects
                     << endl;

                // Specify output option as well
                output = true;
            }

            // Set sizes
            parents.setSize(nIntersects, -1);
            weights.setSize(nIntersects, 0.0);
            centres.setSize(nIntersects, vector::zero);

            // Reset counter
            nIntersects = 0;

            forAll(intersects, indexI)
            {
                if (intersects[indexI])
                {
                    // Compute actual intersections
                    vectorField tP;

                    intersects[indexI] =
                    (
                        faceIntersection
                        (
                            fIndex,
                            candidates[indexI],
                            tP
                        )
                    );

                    // Skip false positives
                    if (intersects[indexI])
                    {
                        parents[nIntersects] = candidates[indexI];

                        // Compute weights
                        convexSetArea
                        (
                            fIndex,
                            parents[nIntersects],
                            tP,
                            weights[nIntersects],
                            centres[nIntersects],
                            output
                        );

                        nIntersects++;
                    }
                }
            }

            // Shorten to actual sizes
            parents.setSize(nIntersects);
            weights.setSize(nIntersects);
            centres.setSize(nIntersects);

            break;
        }
        else
        {
            nAttempts++;

            // Copy / reset parameters
            nOldIntersects = nIntersects;
            nIntersects = 0;

            // Expand the search radius and try again.
            searchFactor *= 1.4;
        }
    }

    // Fetch the face area
    scalar fArea = mag(faceNormal(faces_[fIndex], oldPoints_));
    vector fCentre = faceCentre(faces_[fIndex], oldPoints_);

    // Test weights for consistency
    if (mag(fArea - sum(weights)) > 1e-16)
    {
        // Write out for post-processing
        label uIdx = 0;
        labelList unMatched(candidates.size() - parents.size(), -1);

        forAll(candidates, cI)
        {
            if (findIndex(parents, candidates[cI]) == -1)
            {
                unMatched[uIdx++] = candidates[cI];
            }
        }

        writeVTK("nCell_" + Foam::name(fIndex), fIndex, 2, false, true);
        writeVTK("oCell_" + Foam::name(fIndex), candidates, 2, true, true);
        writeVTK("mCell_" + Foam::name(fIndex), parents, 2, true, true);
        writeVTK("uCell_" + Foam::name(fIndex), unMatched, 2, true, true);

        // Write out weights
        forAll(parents, indexI)
        {
            Info << parents[indexI] << ": "
                 << setprecision(16)
                 << weights[indexI]
                 << endl;
        }

        FatalErrorIn
        (
            "\n"
            "void dynamicTopoFvMesh::computeFaceWeights\n"
            "(\n"
            "    const label fIndex,\n"
            "    const labelList& mapCandidates,\n"
            "    labelList& parents,\n"
            "    scalarField& weights,\n"
            "    vectorField& centres\n"
            ") const"
        )
            << "Encountered non-conservative weighting factors." << nl
            << " Face: " << fIndex << nl
            << " mapCandidates: " << mapCandidates << nl
            << " nCandidates: " << candidates.size() << nl
            << " nParents: " << parents.size() << nl
            << " nAttempts: " << nAttempts << nl
            << setprecision(16)
            << " Face area: " << fArea << nl
            << " Sum(Weights): " << sum(weights) << nl
            << " Error: " << (fArea - sum(weights)) << nl
            << " Norm Sum(Weights): " << sum(weights/fArea) << nl
            << " Norm Error: " << mag(1.0 - sum(weights/fArea))
            << abort(FatalError);
    }

    // Return normalized weights
    weights /= fArea;
}


// Obtain map weighting factors for a cell
void dynamicTopoFvMesh::computeCellWeights
(
    const label cIndex,
    const labelList& mapCandidates,
    labelList& parents,
    scalarField& weights,
    vectorField& centres
) const
{
    scalar searchFactor = 1.0;

    label nOldIntersects = -1, nIntersects = 0, nAttempts = 0;

    // Maintain a list of candidates and intersection points
    labelList candidates;

    // Output option for the convex set algorithm
    bool output = false;

    while (nAttempts < 10)
    {
        // Obtain candidate parents for this cell
        candidates =
        (
            cellParents
            (
                cIndex,
                searchFactor,
                mapCandidates
            )
        );

        // Set sizes
        boolList intersects(candidates.size(), false);

        // Test for intersections
        forAll(candidates, indexI)
        {
            intersects[indexI] =
            (
                testCellIntersection
                (
                    cIndex,
                    candidates[indexI]
                )
            );

            if (intersects[indexI])
            {
                nIntersects++;
            }
        }

        if ((nIntersects == nOldIntersects) && (nIntersects != 0))
        {
            if (debug > 3)
            {
                Info << " Cell: " << cIndex << nl
                     << " nCandidates: " << candidates.size() << nl
                     << " nIntersects: " << nIntersects
                     << endl;

                // Specify output option as well
                output = true;
            }

            // Set sizes
            parents.setSize(nIntersects, -1);
            weights.setSize(nIntersects, 0.0);
            centres.setSize(nIntersects, vector::zero);

            // Reset counter
            nIntersects = 0;

            forAll(intersects, indexI)
            {
                if (intersects[indexI])
                {
                    // Compute actual intersections
                    vectorField tP(0);

                    intersects[indexI] =
                    (
                        cellIntersection
                        (
                            cIndex,
                            candidates[indexI],
                            tP
                        )
                    );

                    // Skip false positives
                    if (intersects[indexI])
                    {
                        parents[nIntersects] = candidates[indexI];

                        // Compute weights
                        convexSetVolume
                        (
                            cIndex,
                            parents[nIntersects],
                            tP,
                            weights[nIntersects],
                            centres[nIntersects],
                            output
                        );

                        nIntersects++;
                    }
                }
            }

            // Shorten to actual sizes
            parents.setSize(nIntersects);
            weights.setSize(nIntersects);
            centres.setSize(nIntersects);

            break;
        }
        else
        {
            nAttempts++;

            // Copy / reset parameters
            nOldIntersects = nIntersects;
            nIntersects = 0;

            // Expand the search radius and try again.
            searchFactor *= 1.4;
        }
    }

    // Fetch the volume of the cell
    scalar cellVolume = 0.0;
    vector cellCentre = vector::zero;

    cellCentreAndVolume
    (
        cIndex,
        oldPoints_,
        faces_,
        cells_,
        owner_,
        cellCentre,
        cellVolume
    );

    // Test weights for consistency
    if (mag(cellVolume - sum(weights)) > 1e-16)
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
                // Normalize by sum of weights
                weights /= sum(weights);

                return;
            }
        }

        // Write out for post-processing
        label uIdx = 0;
        labelList unMatched(candidates.size() - parents.size(), -1);

        forAll(candidates, cI)
        {
            if (findIndex(parents, candidates[cI]) == -1)
            {
                unMatched[uIdx++] = candidates[cI];
            }
        }

        writeVTK("nCell_" + Foam::name(cIndex), cIndex, 3, false, true);
        writeVTK("oCell_" + Foam::name(cIndex), candidates, 3, true, true);
        writeVTK("mCell_" + Foam::name(cIndex), parents, 3, true, true);
        writeVTK("uCell_" + Foam::name(cIndex), unMatched, 3, true, true);

        // Write out weights
        forAll(parents, indexI)
        {
            Info << parents[indexI] << ": "
                 << setprecision(16)
                 << weights[indexI]
                 << endl;
        }

        FatalErrorIn
        (
            "\n"
            "void dynamicTopoFvMesh::computeCellWeights\n"
            "(\n"
            "    const label cIndex,\n"
            "    const labelList& mapCandidates,\n"
            "    labelList& parents,\n"
            "    scalarField& weights,\n"
            "    vectorField& centres\n"
            ") const"
        )
            << "Encountered non-conservative weighting factors." << nl
            << " Cell: " << cIndex << nl
            << " mapCandidates: " << mapCandidates << nl
            << " nCandidates: " << candidates.size() << nl
            << " nParents: " << parents.size() << nl
            << " nAttempts: " << nAttempts << nl
            << setprecision(16)
            << " Cell volume: " << cellVolume << nl
            << " Sum(Weights): " << sum(weights) << nl
            << " Error: " << (cellVolume - sum(weights)) << nl
            << " Norm Sum(Weights): " << sum(weights/cellVolume) << nl
            << " Norm Error: " << mag(1.0 - sum(weights/cellVolume))
            << abort(FatalError);
    }

    // Return normalized weights
    weights /= cellVolume;
}


// Test for intersection between cells in old/new meshes
//   - Uses the static separating axis test for polyhedra,
//     outlined in work by David Eberly
//     'Intersection of Convex Objects: The Method of Separating Axes'
//     http://www.geometrictools.com/
bool dynamicTopoFvMesh::testCellIntersection
(
    const label newCellIndex,
    const label oldCellIndex
) const
{
    // Direction vector
    vector dir(vector::zero);

    // Check indices first
    if (newCellIndex >= cells_.size() || newCellIndex < 0)
    {
        FatalErrorIn
        (
            "\n"
            "bool dynamicTopoFvMesh::testCellIntersection\n"
            "(\n"
            "    const label newCellIndex,\n"
            "    const label oldCellIndex\n"
            ") const"
        )
            << " Wrong newCellIndex: " << newCellIndex << nl
            << " nCells: " << nCells_
            << abort(FatalError);
    }

    if (oldCellIndex >= nOldCells_ || oldCellIndex < 0)
    {
        FatalErrorIn
        (
            "\n"
            "bool dynamicTopoFvMesh::testCellIntersection\n"
            "(\n"
            "    const label newCellIndex,\n"
            "    const label oldCellIndex\n"
            ") const"
        )
            << " Wrong oldCellIndex: " << oldCellIndex << nl
            << " nOldCells: " << nOldCells_
            << abort(FatalError);
    }

    // Fetch references for each mesh
    const faceList& fromFaces = polyMesh::faces();
    const labelList& fromOwner = polyMesh::faceOwner();
    const cell& fromCell = polyMesh::cells()[oldCellIndex];
    const edgeList fromCellEdges = fromCell.edges(polyMesh::faces());
    const labelList fromCellPoints = fromCell.labels(polyMesh::faces());

    const cell& toCell = cells_[newCellIndex];
    const edgeList toCellEdges = toCell.edges(faces_);
    const labelList toCellPoints = toCell.labels(faces_);

    // Test faces of oldCell for separation
    forAll(fromCell, faceI)
    {
        const label fIndex = fromCell[faceI];

        dir = faceNormal(fromFaces[fIndex], oldPoints_);

        // Reverse normal if necessary
        if (fromOwner[fIndex] != oldCellIndex)
        {
            dir *= -1.0;
        }

        if
        (
            whichSide
            (
                toCellPoints,
                oldPoints_,
                dir,
                oldPoints_[fromFaces[fIndex][0]]
            ) > 0
        )
        {
            return false;
        }
    }

    // Test faces of newCell for separation
    forAll(toCell, faceI)
    {
        const label fIndex = toCell[faceI];

        dir = faceNormal(faces_[fIndex], oldPoints_);

        // Reverse normal if necessary
        if (owner_[fIndex] != newCellIndex)
        {
            dir *= -1.0;
        }

        if
        (
            whichSide
            (
                fromCellPoints,
                oldPoints_,
                dir,
                oldPoints_[faces_[fIndex][0]]
            ) > 0
        )
        {
            return false;
        }
    }

    // Test edges of both cells for separation
    forAll(fromCellEdges, edgeI)
    {
        const edge& fromEdge = fromCellEdges[edgeI];

        vector fromVec =
        (
            oldPoints_[fromEdge[1]] - oldPoints_[fromEdge[0]]
        );

        forAll(toCellEdges, edgeJ)
        {
            const edge& toEdge = toCellEdges[edgeJ];

            vector toVec =
            (
                oldPoints_[toEdge[1]] - oldPoints_[toEdge[0]]
            );

            dir = (fromVec ^ toVec);

            label firstSide =
            (
                whichSide
                (
                    fromCellPoints,
                    oldPoints_,
                    dir,
                    oldPoints_[fromEdge[0]]
                )
            );

            if (firstSide == 0)
            {
                continue;
            }

            label secondSide =
            (
                whichSide
                (
                    toCellPoints,
                    oldPoints_,
                    dir,
                    oldPoints_[fromEdge[0]]
                )
            );

            if (secondSide == 0)
            {
                continue;
            }

            if ((firstSide * secondSide) < 0)
            {
                return false;
            }
        }
    }

    return true;
}


// Return the intersection points between faces in old/new meshes
bool dynamicTopoFvMesh::faceIntersection
(
    const label newFaceIndex,
    const label oldFaceIndex,
    vectorField& tP
) const
{
    // Reset inputs
    tP.clear();

    // Fetch references for each mesh
    const face& fromFace = polyMesh::faces()[oldFaceIndex];
    const face& toFace = faces_[newFaceIndex];

    // Obtain face centre and projection normal
    vector xf = faceCentre(toFace, oldPoints_);
    vector nf = faceNormal(toFace, oldPoints_);

    nf /= mag(nf) + VSMALL;

    // Track all possible intersections from here on.
    label nInts = 0;
    Map<vector> intersections;
    vector intPoint = vector::zero;

    // Topologically check for common points,
    // and project uniques ones on to the face plane
    Map<label> projPoints;
    Map<labelList> commonPoints;
    vectorField projections(fromFace.size(), vector::zero);

    forAll(fromFace, pointI)
    {
        label pIndex = findIndex(toFace, fromFace[pointI]);

        vector r = oldPoints_[fromFace[pointI]];

        if (pIndex == -1)
        {
            // Project this point on to the toFace plane.
            projections[pointI] = xf + ((r - xf) - ((r - xf) & nf)*nf);

            projPoints.insert(fromFace[pointI], pointI);
        }
        else
        {
            commonPoints.insert(toFace[pIndex], labelList(0));

            projections[pointI] = r;

            intersections.set(++nInts, r);
        }
    }

    // Add all new points as well, if they resulted
    // from bisections of old face edges.
    forAll(toFace, pointI)
    {
        label pIndex = toFace[pointI];

        if (pIndex >= nOldPoints_)
        {
            // Check pointsFromPoints info
            label index = -1;

            forAll(pointsFromPoints_, indexI)
            {
                if (pointsFromPoints_[indexI].index() == pIndex)
                {
                    index = indexI;
                    break;
                }
            }

            const labelList& mObj = pointsFromPoints_[index].masterObjects();

            // Check if the old face contains all master points
            bool allMaster = true;

            forAll(mObj, pointJ)
            {
                if (findIndex(fromFace, mObj[pointJ]) == -1)
                {
                    allMaster = false;
                    break;
                }
            }

            if (allMaster)
            {
                commonPoints.insert(toFace[pointI], mObj);

                intersections.set(++nInts, oldPoints_[toFace[pointI]]);
            }
        }
    }

    // If all points are common, this is identical to the old face.
    if (nInts == fromFace.size())
    {
        // Copy intersections
        tP.setSize(nInts, vector::zero);

        nInts = 0;

        forAllConstIter(Map<vector>, intersections, pI)
        {
            tP[nInts++] = pI();
        }

        if (debug)
        {
            if (checkPointNearness(tP, 1e-20))
            {
                writeVTK(Foam::name(newFaceIndex),newFaceIndex,2,false,true);
                writeVTK(Foam::name(oldFaceIndex),oldFaceIndex,2,true,true);
                writeVTK
                (
                    "ccSet_"
                  + Foam::name(newFaceIndex)
                  + '<' + Foam::name(oldFaceIndex) + '>',
                    tP.size(),
                    tP.size(),
                    tP.size(),
                    tP
                );
            }
        }

        return true;
    }

    // Check whether any old projections are within
    // the new face. Count these as 'intersections'.
    forAll(fromFace, pointI)
    {
        if (commonPoints.found(fromFace[pointI]))
        {
            continue;
        }

        const point& checkPoint = projections[pointI];

        if (pointInFace(toFace, oldPoints_, checkPoint))
        {
            intersections.set(++nInts, checkPoint);
        }
    }

    // Check whether and new points are within
    // projected old faces. Count these as 'intersections'.
    face ifFace(identity(fromFace.size()));

    forAll(toFace, pointI)
    {
        if (commonPoints.found(toFace[pointI]))
        {
            continue;
        }

        const point& checkPoint = oldPoints_[toFace[pointI]];

        if (pointInFace(ifFace, projections, checkPoint))
        {
            intersections.set(++nInts, checkPoint);
        }
    }

    // Loop through all new edges, and find possible intersections
    // with (projections of) old face edges,
    forAll(toFace, pointI)
    {
        edge toEdge = toFace.faceEdge(pointI);
        label nextLabel = toFace.nextLabel(pointI);

        forAll(fromFace, pointJ)
        {
            label nextJ = fromFace.fcIndex(pointJ);
            edge fromEdge = fromFace.faceEdge(pointJ);

            // Avoid common points
            if (toEdge.commonVertex(fromEdge) > -1)
            {
                continue;
            }

            // Also check for bisection points

            point p1 = projections[pointJ];
            point p2 = projections[nextJ];
            point p3 = oldPoints_[toFace[pointI]];
            point p4 = oldPoints_[nextLabel];

            // Compute edge normal and tangent-to-edge
            vector te = (p4 - p3);
            vector n = (nf ^ te);
            n /= mag(n) + VSMALL;

            // Compute uValues
            scalar numOld = n & (p3 - p1);
            scalar denOld = n & (p2 - p1);

            // Check if the edges are parallel
            scalar tolerance = (1e-4 * mag(p2 - p1));

            if (mag(denOld) < tolerance)
            {
                continue;
            }

            scalar u = (numOld / denOld);
            vector checkPoint = p1 + u*(p2 - p1);
            scalar v = (te & (checkPoint - p3)) / (te & te);

            // Check for intersection along lines.
            if
            (
                ((u > tolerance) && (u < (1.0 - tolerance))) &&
                ((v > tolerance) && (v < (1.0 - tolerance)))
            )
            {
                intersections.set(++nInts, checkPoint);
            }
        }
    }

    // Copy intersections
    tP.setSize(nInts, vector::zero);

    nInts = 0;

    forAllConstIter(Map<vector>, intersections, pI)
    {
        tP[nInts++] = pI();
    }

    // Check for concurrent points.
    if (debug)
    {
        if (checkPointNearness(tP, 1e-20))
        {
            writeVTK(Foam::name(newFaceIndex),newFaceIndex,2,false,true);
            writeVTK(Foam::name(oldFaceIndex),oldFaceIndex,2,true,true);
            writeVTK
            (
                "ccSet_"
              + Foam::name(newFaceIndex)
              + '<' + Foam::name(oldFaceIndex) + '>',
                tP.size(),
                tP.size(),
                tP.size(),
                tP
            );
        }
    }

    // Found a convex set of points.
    if (nInts >= 3)
    {
        return true;
    }

    // Does not intersect
    return false;
}


// Return the intersection points between cells in old/new meshes
bool dynamicTopoFvMesh::cellIntersection
(
    const label newCellIndex,
    const label oldCellIndex,
    vectorField& tP
) const
{
    // Reset inputs
    tP.clear();

    // Fetch references for each mesh
    const cell& fromCell = polyMesh::cells()[oldCellIndex];
    const edgeList fromCellEdges = fromCell.edges(polyMesh::faces());
    const labelList fromCellPoints = fromCell.labels(polyMesh::faces());

    const cell& toCell = cells_[newCellIndex];
    const edgeList toCellEdges = toCell.edges(faces_);
    const labelList toCellPoints = toCell.labels(faces_);

    // Track all possible intersections from here on.
    label nInts = 0;
    Map<vector> intersections;
    vector intPoint = vector::zero;

    // Topologically check for common points
    Map<labelList> commonPoints;

    forAll(fromCellPoints, pointI)
    {
        label pIndex = findIndex(toCellPoints, fromCellPoints[pointI]);

        if (pIndex > -1)
        {
            commonPoints.insert(toCellPoints[pIndex], labelList(0));

            intersections.set(++nInts, oldPoints_[toCellPoints[pIndex]]);
        }
    }

    // Add all new points as well, if they resulted
    // from bisections of old cell edges.
    forAll(toCellPoints, pointI)
    {
        label pIndex = toCellPoints[pointI];

        if (pIndex >= nOldPoints_)
        {
            // Check pointsFromPoints info
            label index = -1;

            forAll(pointsFromPoints_, indexI)
            {
                if (pointsFromPoints_[indexI].index() == pIndex)
                {
                    index = indexI;
                    break;
                }
            }

            const labelList& mObj = pointsFromPoints_[index].masterObjects();

            // Check if the old cell contains all master points
            bool allMaster = true;

            forAll(mObj, pointJ)
            {
                if (findIndex(fromCellPoints, mObj[pointJ]) == -1)
                {
                    allMaster = false;
                    break;
                }
            }

            if (allMaster)
            {
                commonPoints.insert(toCellPoints[pointI], mObj);

                intersections.set(++nInts, oldPoints_[toCellPoints[pointI]]);
            }
        }
    }

    // If all points are common, this is identical to the old cell.
    if (nInts == fromCellPoints.size())
    {
        // Copy intersections
        tP.setSize(nInts, vector::zero);

        nInts = 0;

        forAllConstIter(Map<vector>, intersections, pI)
        {
            tP[nInts++] = pI();
        }

        if (debug)
        {
            if (checkPointNearness(tP, 1e-20))
            {
                writeVTK(Foam::name(newCellIndex),newCellIndex,3,false,true);
                writeVTK(Foam::name(oldCellIndex),oldCellIndex,3,true,true);
                writeVTK
                (
                    "ccSet_"
                  + Foam::name(newCellIndex)
                  + '<' + Foam::name(oldCellIndex) + '>',
                    tP.size(),
                    tP.size(),
                    tP.size(),
                    tP
                );
            }
        }

        return true;
    }

    // Check whether any old points are within
    // the new cell. Count these as 'intersections'.
    forAll(fromCellPoints, pointI)
    {
        if (commonPoints.found(fromCellPoints[pointI]))
        {
            continue;
        }

        const point& checkPoint = oldPoints_[fromCellPoints[pointI]];

        if
        (
            pointInCell
            (
                newCellIndex,
                toCell,
                faces_,
                owner_,
                oldPoints_,
                checkPoint
            )
        )
        {
            intersections.set(++nInts, checkPoint);
        }
    }

    // Check whether any new points are within
    // the old cell. Count these as 'intersections'.
    forAll(toCellPoints, pointI)
    {
        if (commonPoints.found(toCellPoints[pointI]))
        {
            continue;
        }

        const point& checkPoint = oldPoints_[toCellPoints[pointI]];

        if
        (
            pointInCell
            (
                oldCellIndex,
                fromCell,
                polyMesh::faces(),
                polyMesh::faceOwner(),
                oldPoints_,
                checkPoint
            )
        )
        {
            intersections.set(++nInts, checkPoint);
        }
    }

    bool foundIntersection = false, edgeIntersections = false;

    // Loop through edges from each cell, and check whether they intersect.
    List<Pair<edge> > FeToTe, TeToFe;

    forAll(fromCellEdges, edgeI)
    {
        const edge& fromEdge = fromCellEdges[edgeI];

        forAll(toCellEdges, edgeJ)
        {
            const edge& toEdge = toCellEdges[edgeJ];

            foundIntersection = false;

            foundIntersection =
            (
                segmentSegmentIntersection
                (
                    fromEdge,
                    toEdge,
                    oldPoints_,
                    intPoint
                )
            );

            if (foundIntersection)
            {
                FeToTe.setSize
                (
                    FeToTe.size() + 1,
                    Pair<edge>(fromEdge, toEdge)
                );

                TeToFe.setSize
                (
                    TeToFe.size() + 1,
                    Pair<edge>(toEdge, fromEdge)
                );

                intersections.set(++nInts, intPoint);

                // Note for later that edge-intersections exist.
                edgeIntersections = true;
            }
        }
    }

    // Loop through all old edges, and find possible
    // intersections with faces of the new cell.
    forAll(fromCellEdges, edgeI)
    {
        const edge& edgeToCheck = fromCellEdges[edgeI];

        forAll(toCell, faceI)
        {
            const face& faceToCheck = faces_[toCell[faceI]];

            // Avoid edge-edge intersections, if any.
            if (edgeIntersections)
            {
                // Is edgeToCheck in the list?
                bool foundEdge = false;
                const edgeList fEdges = faceToCheck.edges();

                forAll(FeToTe, indexI)
                {
                    if (FeToTe[indexI].first() == edgeToCheck)
                    {
                        // Check whether the intersecting edge
                        // exists on this face.
                        forAll(fEdges, edgeJ)
                        {
                            if (fEdges[edgeJ] == FeToTe[indexI].second())
                            {
                                foundEdge = true;
                                break;
                            }
                        }

                        if (foundEdge)
                        {
                            break;
                        }
                    }
                }

                if (foundEdge)
                {
                    continue;
                }
            }

            bool foundCommon = false;

            forAllConstIter(Map<labelList>, commonPoints, pIter)
            {
                if (faceToCheck.which(pIter.key()) > -1)
                {
                    // Avoid common points, since this implies that
                    // the edge intersects at a face point
                    if
                    (
                        (edgeToCheck[0] == pIter.key()) ||
                        (edgeToCheck[1] == pIter.key())
                    )
                    {
                        foundCommon = true;
                        break;
                    }

                    // Also check for bisection points
                    if
                    (
                        (findIndex(pIter(), edgeToCheck[0]) > -1) &&
                        (findIndex(pIter(), edgeToCheck[1]) > -1)
                    )
                    {
                        foundCommon = true;
                        break;
                    }
                }
            }

            if (foundCommon)
            {
                continue;
            }

            foundIntersection = false;

            foundIntersection =
            (
                segmentFaceIntersection
                (
                    edgeToCheck,
                    faceToCheck,
                    oldPoints_,
                    intPoint
                )
            );

            if (foundIntersection)
            {
                // Add to the list.
                intersections.set(++nInts, intPoint);
            }
        }
    }

    // Loop through all new edges, and find possible
    // intersections with faces of the old cell.
    forAll(toCellEdges, edgeI)
    {
        const edge& edgeToCheck = toCellEdges[edgeI];

        forAll(fromCell, faceI)
        {
            const face& faceToCheck = polyMesh::faces()[fromCell[faceI]];

            // Avoid edge-edge intersections, if any.
            if (edgeIntersections)
            {
                // Is edgeToCheck in the list?
                bool foundEdge = false;
                const edgeList fEdges = faceToCheck.edges();

                forAll(TeToFe, indexI)
                {
                    if (TeToFe[indexI].first() == edgeToCheck)
                    {
                        // Check whether the intersecting edge
                        // exists on this face.
                        forAll(fEdges, edgeJ)
                        {
                            if (fEdges[edgeJ] == TeToFe[indexI].second())
                            {
                                foundEdge = true;
                                break;
                            }
                        }

                        if (foundEdge)
                        {
                            break;
                        }
                    }
                }

                if (foundEdge)
                {
                    continue;
                }
            }

            bool foundCommon = false;

            forAllConstIter(Map<labelList>, commonPoints, pIter)
            {
                // Avoid common points, since this implies that
                // the edge intersects at a face point
                if (faceToCheck.which(pIter.key()) > -1)
                {
                    if
                    (
                        (edgeToCheck[0] == pIter.key()) ||
                        (edgeToCheck[1] == pIter.key())
                    )
                    {
                        foundCommon = true;
                        break;
                    }

                    // Also check for bisection points
                    if
                    (
                        (findIndex(pIter(), edgeToCheck[0]) > -1) &&
                        (findIndex(pIter(), edgeToCheck[1]) > -1)
                    )
                    {
                        foundCommon = true;
                        break;
                    }
                }
            }

            if (foundCommon)
            {
                continue;
            }

            foundIntersection = false;

            foundIntersection =
            (
                segmentFaceIntersection
                (
                    edgeToCheck,
                    faceToCheck,
                    oldPoints_,
                    intPoint
                )
            );

            if (foundIntersection)
            {
                // Add to the list.
                intersections.set(++nInts, intPoint);
            }
        }
    }

    // Copy intersections
    tP.setSize(nInts, vector::zero);

    nInts = 0;

    forAllConstIter(Map<vector>, intersections, pI)
    {
        tP[nInts++] = pI();
    }

    // Check for concurrent points.
    if (debug)
    {
        if (checkPointNearness(tP, 1e-20))
        {
            writeVTK(Foam::name(newCellIndex),newCellIndex,3,false,true);
            writeVTK(Foam::name(oldCellIndex),oldCellIndex,3,true,true);
            writeVTK
            (
                "ccSet_"
              + Foam::name(newCellIndex)
              + '<' + Foam::name(oldCellIndex) + '>',
                tP.size(),
                tP.size(),
                tP.size(),
                tP
            );
        }
    }

    // Found a convex set of points.
    if (nInts >= 4)
    {
        return true;
    }

    // Does not intersect
    return false;
}


// Compute the area / centre of a polygon
// formed by a convex set of points.
void dynamicTopoFvMesh::convexSetArea
(
    const label newFaceIndex,
    const label oldFaceIndex,
    const vectorField& cvxSet,
    scalar& fArea,
    vector& fCentre,
    bool output
) const
{
    // Reset inputs
    fArea = 0.0;
    fCentre = vector::zero;

    // Try the trivial case for a triangle.
    if (cvxSet.size() == 3)
    {
        triPointRef tpr(cvxSet[0], cvxSet[1], cvxSet[2]);

        fArea = tpr.mag();
        fCentre = tpr.centre();

        return;
    }

    // We need a reference normal. Use the new face.
    vector refNorm = faceNormal(faces_[newFaceIndex], oldPoints_);
    refNorm /= mag(refNorm) + VSMALL;

    // Track edges
    label nEdges = 0;
    edgeList testEdges(0);

    // Loop through all points, and build edges with every
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

            // Define the edge
            edge tmpEdge(i, j);

            // If this is an existing edge, skip it.
            bool foundExisting = false;

            forAll(testEdges, edgeI)
            {
                if (testEdges[edgeI] == tmpEdge)
                {
                    foundExisting = true;
                    break;
                }
            }

            if (foundExisting)
            {
                continue;
            }

            // Specify a tolerance for collinearity
            scalar tolerance = 1e-14;

            // Compute the normal to this edge
            vector n = (tmpEdge.vec(cvxSet) ^ refNorm);

            n /= mag(n) + VSMALL;

            label curEdgeSign = 0;
            bool foundInternalEdge = false;

            // Quick-reject test:
            //   Check all other points in the set,
            //   and decide if all points lie on one side.
            forAll(cvxSet, k)
            {
                // Skip duplicates.
                if (tmpEdge[0] == k || tmpEdge[1] == k)
                {
                    continue;
                }

                vector rfVec = (cvxSet[k] - cvxSet[i]);
                scalar dotProd = (rfVec/(mag(rfVec) + VSMALL)) & n;

                // Skip nearly collinear points.
                if (mag(dotProd) < tolerance)
                {
                    continue;
                }

                // Obtain the sign of this point.
                label eSign = Foam::sign(dotProd);

                // Update the current sign if necessary.
                if (curEdgeSign == 0)
                {
                    curEdgeSign = eSign;
                }
                else
                if (curEdgeSign != eSign)
                {
                    // Interior edge. Bail out.
                    foundInternalEdge = true;
                    break;
                }
            }

            if (foundInternalEdge)
            {
                continue;
            }

            // Looks like we found an edge on the boundary.
            // Check its sign to ensure that it points outward.
            if (curEdgeSign == 1)
            {
                n *= -1.0;
                tmpEdge = tmpEdge.reverseEdge();
            }

            // Add to the list of edges.
            testEdges.setSize(++nEdges, tmpEdge);
        }
    }

    // Sanity check - do points match edges?
    if (testEdges.size() != cvxSet.size())
    {
        FatalErrorIn
        (
            "\n"
            "void dynamicTopoFvMesh::convexSetArea\n"
            "(\n"
            "    const label newFaceIndex,\n"
            "    const label oldFaceIndex,\n"
            "    const vectorField& cvxSet,\n"
            "    scalar& fArea,\n"
            "    vector& fCentre,\n"
            "    bool output\n"
            ") const"
        )
            << " Points do not match edges. " << nl
            << " nPoints: " << cvxSet.size() << nl
            << " nEdges: " << testEdges.size() << nl
            << " Edge list: " << testEdges
            << abort(FatalError);
    }

    // Find an approximate face-centroid
    scalar sumA = 0.0;
    vector sumAc = vector::zero;
    vector xC = average(cvxSet);

    forAll(testEdges, edgeI)
    {
        const edge& e = testEdges[edgeI];

        vector c = cvxSet[e[0]] + cvxSet[e[1]] + xC;
        scalar a = mag(e.vec(cvxSet) ^ (xC - cvxSet[e[0]]));

        sumA += a;
        sumAc += a*c;
    }

    fCentre = (1.0/3.0)*sumAc/(sumA + VSMALL);
    fArea = 0.5*sumA;

    if (output)
    {
        Info << " newFaceIndex: " << newFaceIndex
             << " oldFaceIndex: " << oldFaceIndex << nl
             << " Edges: " << testEdges << nl
             << " Area: " << fArea << nl
             << " Centre: " << fCentre << nl
             << endl;
    }
}


// Compute the volume / centre of a polyhedron
// formed by a convex set of points.
void dynamicTopoFvMesh::convexSetVolume
(
    const label newCellIndex,
    const label oldCellIndex,
    const vectorField& cvxSet,
    scalar& cVolume,
    vector& cCentre,
    bool output
) const
{
    // Reset inputs
    cVolume = 0.0;
    cCentre = vector::zero;

    // Try the trivial case for a tetrahedron.
    // No checking for orientation here.
    if (cvxSet.size() == 4)
    {
        cCentre = average(cvxSet);

        cVolume =
        (
            mag
            (
                tetPointRef
                (
                    cvxSet[0],
                    cvxSet[1],
                    cvxSet[2],
                    cvxSet[3]
                ).mag()
            )
        );

        return;
    }

    // Track faces
    face tmpFace(3);
    label nFaces = 0;
    faceList testFaces(0);
    labelHashSet uniquePts;

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

                // Quick-reject test:
                //   If this is a subset of an existing face, skip it.
                bool foundSubSet = false;

                forAll(testFaces, faceI)
                {
                    const face& checkFace = testFaces[faceI];

                    if (checkFace.size() >= tmpFace.size())
                    {
                        bool foundUniquePoint = false;

                        forAll(tmpFace, pI)
                        {
                            if (findIndex(checkFace, tmpFace[pI]) == -1)
                            {
                                foundUniquePoint = true;
                                break;
                            }
                        }

                        if (!foundUniquePoint)
                        {
                            foundSubSet = true;
                            break;
                        }
                    }
                }

                if (foundSubSet)
                {
                    continue;
                }

                // Specify a tolerance for planarity
                scalar tolerance = 1e-14;

                // Compute the normal to this face
                vector n = tmpFace.normal(cvxSet);

                n /= mag(n) + VSMALL;

                label curFaceSign = 0;
                bool foundInternalFace = false;

                // Quick-reject test:
                //   Check all other points in the set,
                //   and decide if all points lie on one side.
                forAll(cvxSet, l)
                {
                    // Skip duplicates.
                    if (tmpFace[0] == l || tmpFace[1] == l || tmpFace[2] == l)
                    {
                        continue;
                    }

                    vector rfVec = (cvxSet[l] - cvxSet[i]);
                    scalar dotProd = (rfVec/(mag(rfVec) + VSMALL)) & n;

                    // Skip nearly co-planar points.
                    if (mag(dotProd) < tolerance)
                    {
                        continue;
                    }

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

                if (foundInternalFace)
                {
                    continue;
                }

                // Looks like we found a face on the boundary.
                // Check its sign to ensure that it points outward.
                if (curFaceSign == 1)
                {
                    n *= -1.0;
                    tmpFace = tmpFace.reverseFace();
                }

                // Ensure that the face wasn't checked in.
                bool alreadyCheckedIn = false;

                forAll(testFaces, faceI)
                {
                    // Fetch a non-const reference, since this face
                    // might be modified in this loop.
                    face& checkFace = testFaces[faceI];

                    label nCommon = 0;

                    uniquePts.clear();

                    forAll(tmpFace, pI)
                    {
                        if (findIndex(checkFace, tmpFace[pI]) > -1)
                        {
                            nCommon++;
                        }
                        else
                        {
                            uniquePts.insert(tmpFace[pI]);
                        }
                    }

                    if (nCommon >= 2)
                    {
                        if (checkFace.size() >= tmpFace.size())
                        {
                            // Check for unique points
                            if (uniquePts.size() > 0)
                            {
                                // Compute the existing normal
                                vector eNorm = checkFace.normal(cvxSet);

                                scalar dotProd =
                                (
                                    n & (eNorm/(mag(eNorm) + VSMALL))
                                );

                                if
                                (
                                    (mag(1.0 - dotProd) < tolerance) &&
                                    (dotProd > 0.0)
                                )
                                {
                                    // Add all unique points to checkFace
                                    insertPointLabels
                                    (
                                        n,
                                        cvxSet,
                                        uniquePts,
                                        checkFace
                                    );

                                    alreadyCheckedIn = true;
                                    break;
                                }
                            }
                            else
                            {
                                // Subset face
                                alreadyCheckedIn = true;
                                break;
                            }
                        }
                        else
                        {
                            // checkFace is a subset. Replace it.
                            checkFace = tmpFace;

                            alreadyCheckedIn = true;
                            break;
                        }
                    }
                }

                // Add this face to the list of faces.
                if (!alreadyCheckedIn)
                {
                    testFaces.setSize(++nFaces, tmpFace);
                }

                // Reset the face size.
                tmpFace.setSize(3, -1);
            }
        }
    }

    // Account for planarity test failure.
    //  - Check for subsets.
    forAll(testFaces, faceI)
    {
        // Fetch a non-const reference, since this face
        // might be modified in this loop.
        face& checkFace = testFaces[faceI];

        // Account for deleted testFaces
        if (checkFace.empty())
        {
            continue;
        }

        // Compute the normal to this face
        vector n = checkFace.normal(cvxSet);

        forAll(testFaces, faceJ)
        {
            if (faceI == faceJ)
            {
                continue;
            }

            // Fetch a non-const reference, since this face
            // might be modified in this loop.
            face& testFace = testFaces[faceJ];

            if (checkFace.size() >= testFace.size())
            {
                label nCommon = 0;

                uniquePts.clear();

                forAll(testFace, pI)
                {
                    if (findIndex(checkFace, testFace[pI]) > -1)
                    {
                        nCommon++;
                    }
                    else
                    {
                        uniquePts.insert(testFace[pI]);
                    }
                }

                if (nCommon >= 3)
                {
                    // Delete the test face
                    testFace.clear();

                    // Add all unique points to checkFace
                    // Failed the tolerance test before,
                    // so don't check for it now
                    if (uniquePts.size())
                    {
                        insertPointLabels
                        (
                            n,
                            cvxSet,
                            uniquePts,
                            checkFace
                        );
                    }
                }
            }
        }
    }

    // Find an approximate cell-centroid
    vector xC = average(cvxSet);

    // Calculate volume from all accumulated faces.
    forAll(testFaces, faceI)
    {
        const face& checkFace = testFaces[faceI];

        if (checkFace.empty())
        {
            continue;
        }

        vector xF = checkFace.centre(cvxSet);
        vector Sf = checkFace.normal(cvxSet);

        // Calculate 3*face-pyramid volume
        scalar pyr3Vol = Sf & (xF - xC);

        // Calculate face-pyramid centre
        vector pc = (3.0/4.0)*xF + (1.0/4.0)*xC;

        // Accumulate volume-weighted face-pyramid centre
        cCentre += pyr3Vol*pc;

        // Accumulate face-pyramid volume
        cVolume += pyr3Vol;
    }

    cCentre /= cVolume + VSMALL;
    cVolume *= (1.0/3.0);

    if (output)
    {
        // Write out faces as a standalone patch
        fileName dirName(time().path()/"VTK"/time().timeName());
        mkDir(dirName);

        PrimitivePatch<face, List, pointField>::writeVTK
        (
            dirName/fileName
            (
                "int_"
              + Foam::name(newCellIndex) + '_'
              + Foam::name(oldCellIndex)
            ),
            testFaces,
            cvxSet
        );

        Info << " newCellIndex: " << newCellIndex
             << " oldCellIndex: " << oldCellIndex << nl
             << " Faces: " << testFaces << nl
             << " Volume: " << cVolume << nl
             << " Centre: " << cCentre << nl
             << endl;
    }
}


// Obtain a list of possible parent faces from the old mesh
labelList dynamicTopoFvMesh::faceParents
(
    const label fIndex,
    const scalar searchFactor,
    const labelList& oldCandidates
) const
{
    labelStaticHashSet masterFaces;

    const face& faceToCheck = faces_[fIndex];
    const polyBoundaryMesh& boundary = boundaryMesh();

    // Determine the patch that this face belongs to..
    label patchIndex = whichPatch(fIndex);

    if (patchIndex == -1 || faceToCheck.empty())
    {
        FatalErrorIn
        (
            "\n"
            "labelList dynamicTopoFvMesh::faceParents\n"
            "(\n"
            "    const label fIndex,\n"
            "    const scalar searchFactor,\n"
            "    const labelList& oldCandidates\n"
            ") const"
        )
            << nl << " Illegal request for face: "
            << fIndex << ":: " << faces_[fIndex] << nl
            << " oldCandidates: " << oldCandidates
            << abort(FatalError);
    }

    // Fetch old patch start
    label patchStart = boundary[patchIndex].start();

    // Insert the old candidates first
    // Assume that candidates / parents are global face indices,
    // and all addressing into the same patch.
    forAll(oldCandidates, faceI)
    {
        if (oldCandidates[faceI] < 0)
        {
            continue;
        }

        if (oldCandidates[faceI] < nOldFaces_)
        {
            if (whichPatch(oldCandidates[faceI]) == patchIndex)
            {
                masterFaces.insert
                (
                    (oldCandidates[faceI] - patchStart),
                    empty()
                );
            }
        }
        else
        if (faceParents_.found(oldCandidates[faceI]))
        {
            const labelList& nParents = faceParents_[oldCandidates[faceI]];

            forAll(nParents, fI)
            {
                if (whichPatch(nParents[fI]) == patchIndex)
                {
                    masterFaces.insert
                    (
                        (nParents[fI] - patchStart),
                        empty()
                    );
                }
            }
        }
    }

    // Fetch connectivity from the old mesh.
    const labelListList& oldFaceFaces = boundary[patchIndex].faceFaces();

    for (label attempt = 0; attempt < 5; attempt++)
    {
        // Fetch the initial set of candidates
        labelList initList = masterFaces.toc();

        // Accumulate a larger stencil of face neighbours
        forAll(initList, indexI)
        {
            const labelList& ff = oldFaceFaces[initList[indexI]];

            forAll(ff, faceI)
            {
                masterFaces.insert(ff[faceI], empty());
            }
        }
    }

    // Fetch the list of face points
    pointField facePoints(faceToCheck.size());

    forAll(facePoints, pointI)
    {
        facePoints[pointI] = oldPoints_[faceToCheck[pointI]];
    }

    // Prepare an axis-aligned bounding box around the face,
    // and add faces according to face-centre positions.
    boundBox faceBox(facePoints, false);

    // Define a search radius
    vector bC = faceBox.midpoint();
    vector bMax = searchFactor * (faceBox.max() - bC);

    const vectorField& faceCentres = boundary[patchIndex].faceCentres();

    label nEntries = 0;
    labelList finalFaces(masterFaces.size(), -1);

    forAllConstIter(labelStaticHashSet, masterFaces, fIter)
    {
        vector xC = (faceCentres[fIter.key()] - bC);

        if ((xC & xC) < (bMax & bMax))
        {
            // Store global indices in finalFaces
            finalFaces[nEntries++] = (fIter.key() + patchStart);
        }
    }

    // Shrink to actual size
    finalFaces.setSize(nEntries);

    if (debug > 3)
    {
        Info << " Face: " << fIndex
             << " No. of parent candidates: "
             << nEntries
             << " searchFactor: "
             << searchFactor
             << endl;
    }

    return finalFaces;
}


// Obtain a list of possible parent cells from the old mesh.
labelList dynamicTopoFvMesh::cellParents
(
    const label cIndex,
    const scalar searchFactor,
    const labelList& oldCandidates
) const
{
    labelStaticHashSet masterCells;

    // Insert the old candidates first
    forAll(oldCandidates, cellI)
    {
        if (oldCandidates[cellI] < 0)
        {
            continue;
        }

        if (oldCandidates[cellI] < nOldCells_)
        {
            masterCells.insert(oldCandidates[cellI], empty());
        }
        else
        if (cellParents_.found(oldCandidates[cellI]))
        {
            const labelList& nParents = cellParents_[oldCandidates[cellI]];

            forAll(nParents, cI)
            {
                masterCells.insert(nParents[cI], empty());
            }
        }
    }

    // Fetch connectivity from the old mesh.
    const labelListList& oldCellCells = polyMesh::cellCells();

    for (label attempt = 0; attempt < 5; attempt++)
    {
        // Fetch the initial set of candidates
        labelList initList = masterCells.toc();

        // Accumulate a larger stencil of cell neighbours
        forAll(initList, indexI)
        {
            const labelList& cc = oldCellCells[initList[indexI]];

            forAll(cc, cellI)
            {
                masterCells.insert(cc[cellI], empty());
            }
        }
    }

    // Fetch the new cell, and determine its bounds.
    const cell& newCell = cells_[cIndex];

    labelList pLabels = newCell.labels(faces_);

    pointField cellPoints(pLabels.size());

    forAll(cellPoints, pointI)
    {
        cellPoints[pointI] = oldPoints_[pLabels[pointI]];
    }

    // Prepare an axis-aligned bounding box around the cell,
    // and add cells according to cell-centre positions.
    boundBox cellBox(cellPoints, false);

    // Define a search radius
    vector bC = cellBox.midpoint();
    vector bMax = searchFactor * (cellBox.max() - bC);

    const vectorField& cellCentres = polyMesh::cellCentres();

    label nEntries = 0;
    labelList finalCells(masterCells.size(), -1);

    forAllConstIter(labelStaticHashSet, masterCells, cIter)
    {
        vector xC = (cellCentres[cIter.key()] - bC);

        if ((xC & xC) < (bMax & bMax))
        {
            finalCells[nEntries++] = cIter.key();
        }
    }

    // Shrink to actual size
    finalCells.setSize(nEntries);

    if (debug > 3)
    {
        Info << " Cell: " << cIndex
             << " No. of parent candidates: "
             << nEntries
             << " searchFactor: "
             << searchFactor
             << endl;
    }

    return finalCells;
}


// Set fill-in mapping information for a particular cell
void dynamicTopoFvMesh::setCellMapping
(
    const label cIndex,
    const labelList& mapCells,
    const scalarField& mapWeights,
    const vectorField& mapCentres
)
{
    if (debug > 3)
    {
        Info << "Inserting mapping cell: " << cIndex << nl
             << " mapCells: " << mapCells << nl
             << " mapWeights: " << mapWeights << nl
             << " mapCentres: " << mapCentres
             << endl;
    }

    // Ensure compatible sizes.
    if (mapCells.size() != mapWeights.size() || mapCells.empty())
    {
        FatalErrorIn
        (
            "\n"
            "void dynamicTopoFvMesh::setCellMapping\n"
            "(\n"
            "    const label cIndex,\n"
            "    const labelList& mapCells,\n"
            "    const scalarField& mapWeights,\n"
            "    const vectorField& mapCentres\n"
            ")"
        )
            << nl << " Incompatible mapping for cell: "
            << cIndex << ":: " << cells_[cIndex] << nl
            << " mapCells: " << mapCells << nl
            << " mapWeights: " << mapWeights << nl
            << " mapCentres: " << mapCentres
            << abort(FatalError);
    }

    // Insert weights into the list, and overwrite if necessary
    label index = -1;

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

    // Set weights and centres
    cellWeights_.set(cIndex, mapWeights);
    cellCentres_.set(cIndex, mapCentres);
}


// Set fill-in mapping information for a particular face
void dynamicTopoFvMesh::setFaceMapping
(
    const label fIndex,
    const labelList& mapFaces,
    const scalarField& mapWeights,
    const vectorField& mapCentres
)
{
    label patch = whichPatch(fIndex);

    if (debug > 3)
    {
        Info << "Inserting mapping face: " << fIndex << nl
             << " patch: " << patch << nl
             << " mapFaces: " << mapFaces << nl
             << " mapWeights: " << mapWeights << nl
             << " mapCentres: " << mapCentres
             << endl;
    }

    bool foundError = false;

    // Check to ensure that internal faces are not mapped
    // from any faces in the mesh
    if (patch == -1 && mapFaces.size())
    {
        foundError = true;
    }

    // Check to ensure that boundary faces map
    // only from other faces on the same patch
    if (patch > -1 && mapFaces.empty())
    {
        foundError = true;
    }

    // Ensure compatible sizes.
    if
    (
        (mapFaces.size() != mapWeights.size()) ||
        (mapFaces.size() != mapCentres.size())
    )
    {
        foundError = true;
    }

    if (foundError)
    {
        FatalErrorIn
        (
            "\n"
            "void dynamicTopoFvMesh::setFaceMapping\n"
            "(\n"
            "    const label fIndex,\n"
            "    const labelList& mapFaces,\n"
            "    const scalarField& mapWeights,\n"
            "    const vectorField& mapCentres\n"
            ")"
        )
            << nl << " Incompatible mapping: " << nl
            << " Face: " << fIndex << nl
            << " mapFaces: " << mapFaces << nl
            << " mapWeights: " << mapWeights << nl
            << " mapCentres: " << mapCentres
            << abort(FatalError);
    }

    // Insert addressing into the list, and overwrite if necessary
    label index = -1;

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

    if (mapFaces.size())
    {
        // Set weights and centres
        faceWeights_.set(fIndex, mapWeights);
        faceCentres_.set(fIndex, mapCentres);
    }
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
        cellCentres_.erase(cIndex);
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
        forAll(entityStack_, stackI)
        {
            Stack(stackI).remove(fIndex);
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
        faceCentres_.erase(fIndex);
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
                FatalErrorIn
                (
                    "\n"
                    "label dynamicTopoFvMesh::insertEdge\n"
                    "(\n"
                    "    const label patch,\n"
                    "    const edge& newEdge,\n"
                    "    const labelList& edgeFaces,\n"
                    "    const labelList& edgePoints\n"
                    ")"
                )
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
        forAll(entityStack_, stackI)
        {
            Stack(stackI).remove(eIndex);
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
    const labelList& mapPoints,
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
             << mapPoints << endl;
    }

    // Make a pointsFromPoints entry
    sizeUpList
    (
        objectMap(newPointIndex, mapPoints),
        pointsFromPoints_
    );

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
                        "\n"
                        "void dynamicTopoFvMesh::constructHull\n"
                        "(\n"
                        "    const label eIndex,\n"
                        "    labelList& hullEdges,\n"
                        "    labelList& hullFaces,\n"
                        "    labelList& hullCells,\n"
                        "    labelListList& ringEntities\n"
                        ") const"
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
                        if (triFace::compare(triFace(cFace), triFace(oFace)))
                        {
                            ringEntities[1][indexI] = currCell[faceI];
                        }

                        // Build a comparison face
                        oFace[0] = edgeToCheck[1];
                        oFace[1] = nextHullPoint;
                        oFace[2] = otherPoint;

                        // Check if this face contains edgeToCheck[1]
                        if (triFace::compare(triFace(cFace), triFace(oFace)))
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
                "\n"
                "void dynamicTopoFvMesh::constructHull\n"
                "(\n"
                "    const label eIndex,\n"
                "    labelList& hullEdges,\n"
                "    labelList& hullFaces,\n"
                "    labelList& hullCells,\n"
                "    labelListList& ringEntities\n"
                ") const"
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
                "\n"
                "void dynamicTopoFvMesh::buildEdgePoints\n"
                "(\n"
                "    const label eIndex,\n"
                "    const label checkIndex\n"
                ")"
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


void dynamicTopoFvMesh::handleMeshSlicing()
{
    if (slicePairs_.empty())
    {
        return;
    }

    if (lengthEstimator().holdOff())
    {
        // Clear out data.
        slicePairs_.clear();
        lengthEstimator().clearBoxes();

        // Hold-off mesh slicing for a few time-steps.
        lengthEstimator().decrementHoldOff();

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
    slicePairs_.clear();
    lengthEstimator().clearBoxes();

    // Set the sliceHoldOff value
    lengthEstimator().setHoldOff(50);
}


// Test an edge / face for proximity with other faces on proximity patches
// and return the scalar distance to an oppositely-oriented face.
scalar dynamicTopoFvMesh::testProximity
(
    const label index,
    label& proximityFace
) const
{
    scalar proxDistance = GREAT, testStep = 0.0;
    vector gCentre = vector::zero, gNormal = vector::zero;

    if (twoDMesh_)
    {
        // Obtain the face-normal.
        gNormal = faceNormal(faces_[index], points_);

        // Obtain the face centre.
        gCentre = faceCentre(faces_[index], points_);

        // Calculate a test step-size
        testStep = edgeLength(edges_[getTriBoundaryEdge(index)], points_);
    }
    else
    {
        const edge& thisEdge = edges_[index];
        const labelList& eFaces = edgeFaces_[index];

        // Obtain the edge centre.
        gCentre = 0.5 * (points_[thisEdge[0]] + points_[thisEdge[1]]);

        // Obtain the edge-normal
        forAll(eFaces, faceI)
        {
            if (neighbour_[eFaces[faceI]] == -1)
            {
                // Obtain the normal.
                gNormal += faceNormal(faces_[eFaces[faceI]], points_);
            }
        }

        // Calculate a test step-size
        testStep = edgeLength(edges_[index], points_);
    }

    // Normalize
    gNormal /= (mag(gNormal) + VSMALL);

    // Test for proximity, and mark slice-pairs
    // is the distance is below threshold.
    if
    (
        lengthEstimator().testProximity
        (
            gCentre,
            gNormal,
            testStep,
            proximityFace,
            proxDistance
        )
    )
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

            label newSize = slicePairs_.size() + 1;

            // Const-cast slicePairs for resize
            List<labelPair>& sP = const_cast<List<labelPair>&>(slicePairs_);

            // Add this entry as a candidate for mesh slicing.
            sP.setSize(newSize, proxPoints);

            // Unlock the edge mutex
            entityMutex_[1].unlock();
        }
    }

    return proxDistance;
}


// Calculate the edge length-scale for the mesh
void dynamicTopoFvMesh::calculateLengthScale(bool dump)
{
    if (!edgeRefinement_)
    {
        return;
    }

    Switch dumpLengthScale(false);

    const dictionary& meshDict = dict_.subDict("dynamicTopoFvMesh");

    if (meshDict.found("dumpLengthScale") || mandatory_)
    {
        dumpLengthScale = meshDict.lookup("dumpLengthScale");
    }

    autoPtr<volScalarField> lsfPtr(NULL);

    if (dumpLengthScale && time().outputTime() && dump)
    {
        lsfPtr.set
        (
            new volScalarField
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
            )
        );
    }

    // Bail-out if a dumping was not requested in dictionary.
    if (dump && !dumpLengthScale && !time().outputTime())
    {
        return;
    }

    // Size the field and calculate length-scale
    lengthScale_.setSize(nCells_, 0.0);

    lengthEstimator().calculateLengthScale(lengthScale_);

    // Check if length-scale is to be dumped to disk.
    if (dumpLengthScale && time().outputTime() && dump)
    {
        // Obtain length-scale values from the mesh
        lsfPtr->internalField() = lengthScale_;

        lsfPtr->write();
    }

    if (Pstream::parRun())
    {
        // Exchange length-scale buffers across processors.
        exchangeLengthBuffers();

        // Wait for transfers before continuing.
        waitForBuffers();
    }
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

    // Set debug option for underlying classes as well.
    lengthScaleEstimator::debug = debug;

    if (debug > 2)
    {
        fvMesh::debug = true;
        polyMesh::debug = true;
    }

    const dictionary& meshSubDict = dict_.subDict("dynamicTopoFvMesh");

    if (meshSubDict.found("interval") || mandatory_)
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

    if (meshSubDict.found("bandWidthReduction") || mandatory_)
    {
        bandWidthReduction_ = meshSubDict.lookup("bandWidthReduction");
    }
    else
    {
        bandWidthReduction_ = false;
    }

    if (meshSubDict.found("sliverThreshold") || mandatory_)
    {
        sliverThreshold_ = readScalar(meshSubDict.lookup("sliverThreshold"));

        if (sliverThreshold_ > 1.0 || sliverThreshold_ < 0.0)
        {
            FatalErrorIn("void dynamicTopoFvMesh::readOptionalParameters()")
                << " Sliver threshold out of range [0..1]"
                << abort(FatalError);
        }
    }

    if (meshSubDict.found("maxModifications") || mandatory_)
    {
        maxModifications_ = readLabel(meshSubDict.lookup("maxModifications"));
    }

    // For tetrahedral meshes...
    if (!twoDMesh_)
    {
        // Check if swapping is to be avoided on any patches
        if (meshSubDict.found("noSwapPatches") || mandatory_)
        {
            wordList noSwapPatches = meshSubDict.subDict("noSwapPatches").toc();

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
        if (meshSubDict.found("maxTetsPerEdge") || mandatory_)
        {
            maxTetsPerEdge_ = readLabel(meshSubDict.lookup("maxTetsPerEdge"));
        }
        else
        {
            maxTetsPerEdge_ = 7;
        }

        // Check if programming tables can be resized at runtime
        if (meshSubDict.found("allowTableResize") || mandatory_)
        {
            allowTableResize_ =
            (
                readBool(meshSubDict.lookup("allowTableResize"))
            );
        }
        else
        {
            allowTableResize_ = false;
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
        "label dynamicTopoFvMesh::getTriBoundaryFace"
        "(const label fIndex) const"
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
        "label dynamicTopoFvMesh::getTriBoundaryEdge"
        "(const label fIndex) const"
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
    if (motionSolver_.valid())
    {
        FatalErrorIn
        (
            "void dynamicTopoFvMesh::loadMotionSolver() "
        ) << nl << " Motion solver already loaded. "
          << abort(FatalError);
    }
    else
    if (dict_.found("solver"))
    {
        motionSolver_ = motionSolver::New(*this);
    }
}


// Load the field mapper
void dynamicTopoFvMesh::loadFieldMapper()
{
    if (mapper_.valid())
    {
        FatalErrorIn
        (
            "void dynamicTopoFvMesh::loadFieldMapper() "
        ) << nl << " Field mapper already loaded. "
          << abort(FatalError);
    }
    else
    {
        mapper_.set(new topoMapper(*this));
    }
}


// Load the length scale estimator
void dynamicTopoFvMesh::loadLengthScaleEstimator()
{
    if (edgeRefinement_)
    {
        if (lengthEstimator_.valid())
        {
            FatalErrorIn
            (
                "void dynamicTopoFvMesh::loadLengthScaleEstimator() "
            ) << nl << " Length scale estimator already loaded. "
              << abort(FatalError);
        }
        else
        {
            const dictionary& meshDict = dict_.subDict("dynamicTopoFvMesh");

            lengthEstimator_.set
            (
                new lengthScaleEstimator
                (
                    *this,
                    meshDict.subDict("refinementOptions")
                )
            );
        }

        // Read options
        lengthEstimator().readRefinementOptions(false, mandatory_);

        // Set coupled patch options, if available
        if (dict_.found("coupledPatches") || mandatory_)
        {
            lengthEstimator().setCoupledPatches
            (
                dict_.subDict("coupledPatches")
            );
        }
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
            new meshHandler(*this, threader())
        );

        handlerPtr_[0].setMaster();

        // Size the stacks
        entityStack_.setSize(1);
    }
    else
    {
        // Index '0' is master, rest are slaves
        handlerPtr_.setSize(nThreads + 1);

        // Size the stacks
        entityStack_.setSize(nThreads + 1);

        forAll(handlerPtr_, threadI)
        {
            handlerPtr_.set
            (
                threadI,
                new meshHandler(*this, threader())
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
    meshHandler *thread = static_cast<meshHandler*>(argument);

    if (thread->slave())
    {
        thread->sendSignal(meshHandler::START);
    }

    dynamicTopoFvMesh& mesh = thread->reference();

    // Figure out which thread this is...
    label tIndex = mesh.self();

    // Pick items off the stack
    while (!mesh.Stack(tIndex).empty())
    {
        // Retrieve the index for this face
        label fIndex = mesh.Stack(tIndex).pop();

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
                mesh.Stack(0).push(fIndex);
            }
        }
    }

    if (thread->slave())
    {
        thread->sendSignal(meshHandler::STOP);
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
        FatalErrorIn
        (
            "void dynamicTopoFvMesh::swapEdge"
            "(const label eIndex, bool forceOp)"
        )
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
                FatalErrorIn
                (
                    "void dynamicTopoFvMesh::swapEdge"
                    "(const label eIndex, bool forceOp)"
                )
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
    meshHandler *thread = static_cast<meshHandler*>(argument);

    if (thread->slave())
    {
        thread->sendSignal(meshHandler::START);
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
    while (!mesh.Stack(tIndex).empty())
    {
        // Retrieve an edge from the stack
        label eIndex = mesh.Stack(tIndex).pop();

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
                    mesh.Stack(0).push(eIndex);
                }
            }
        }
    }

    if (thread->slave())
    {
        thread->sendSignal(meshHandler::STOP);
    }
}


// Edge refinement engine
void dynamicTopoFvMesh::edgeRefinementEngine
(
    void *argument
)
{
    // Loop through all edges and bisect/collapse by the criterion:
    // Bisect when edge-length > ratioMax_*lengthScale
    // Collapse when edge-length < ratioMin_*lengthScale

    // Recast the argument
    meshHandler *thread = static_cast<meshHandler*>(argument);

    if (thread->slave())
    {
        thread->sendSignal(meshHandler::START);
    }

    dynamicTopoFvMesh& mesh = thread->reference();

    // Figure out which thread this is...
    label tIndex = mesh.self();

    while (!mesh.Stack(tIndex).empty())
    {
        // Retrieve an entity from the stack
        label eIndex = mesh.Stack(tIndex).pop();

        if (mesh.checkBisection(eIndex))
        {
            if (thread->master())
            {
                // Bisect this edge
                mesh.bisectEdge(eIndex);
            }
            else
            {
                // Push this on to the master stack
                mesh.Stack(0).push(eIndex);
            }
        }
        else
        if (mesh.checkCollapse(eIndex))
        {
            if (thread->master())
            {
                // Collapse this edge
                mesh.collapseEdge(eIndex);
            }
            else
            {
                // Push this on to the master stack
                mesh.Stack(0).push(eIndex);
            }
        }
    }

    if (thread->slave())
    {
        thread->sendSignal(meshHandler::STOP);
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
            oldVolume = tetPointRef(aOld, bOld, oldPoint, dOld).mag();

            // Check if the volume / quality is worse
            minQuality = Foam::min(cQuality, minQuality);
            minVolume = Foam::min(oldVolume, minVolume);

            // Compute the quality of the lower half.
            cQuality = tetMetric_(midPoint, b, c, d);

            // Compute old volume of the lower half.
            oldVolume = tetPointRef(oldPoint, bOld, cOld, dOld).mag();

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
            oldVolume = tetPointRef(aOld, bOld, oldPoint, dOld).mag();

            // Check if the quality is worse
            minQuality = Foam::min(cQuality, minQuality);
            minVolume = Foam::min(oldVolume, minVolume);

            // Compute the quality of the lower half.
            cQuality = tetMetric_(midPoint, b, c, d);

            // Compute old volume of the lower half.
            oldVolume = tetPointRef(oldPoint, bOld, cOld, dOld).mag();

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
            InfoIn
            (
                "scalar dynamicTopoFvMesh::computeBisectionQuality"
                "(const label eIndex) const"
            )
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

    point midPoint = faceCentre(faces_[fIndex], points_);

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
    scalar length = edgeLength(edges_[triEdge], points_);

    // Determine the boundary triangular face area
    scalar area = mag(faceNormal(faces_[triFace], points_));

    // This cell has to be removed...
    if (mag(area) < (0.2*length*length))
    {
        if (self() == 0)
        {
            if (debug > 1)
            {
                InfoIn
                (
                    "void dynamicTopoFvMesh::remove2DSliver"
                    "(const label fIndex)"
                )
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
            Stack(0).push(fIndex);
        }
    }
}


// Indentify the sliver type in 3D
const changeMap dynamicTopoFvMesh::identifySliverType
(
    const label cIndex
) const
{
    changeMap map;

    // Ensure that this cell actually exists.
    if (cells_[cIndex].empty())
    {
        return map;
    }

    label fourthPoint = -1;
    scalar minDistance = GREAT;
    face tFace(3), testFace(3), faceToCheck(3);
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
        vector testNormal = faceNormal(testFace, points_);

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
            tFace = testFace;
            minDistance = distance;
        }
    }

    // Obtain the face-normal.
    vector refArea = faceNormal(tFace, points_);

    // Normalize it.
    vector n = refArea/mag(refArea);

    // Define edge-vectors.
    vector r1 = points_[tFace[1]] - points_[tFace[0]];
    vector r2 = points_[tFace[2]] - points_[tFace[1]];
    vector r3 = points_[tFace[0]] - points_[tFace[2]];

    // Project the fourth point onto the face.
    vector r4 = points_[fourthPoint] - points_[tFace[0]];

    r4 = r4 - ((r4 & n)*n);

    // Define the two other vectors.
    vector r5 = r4 - r1;
    vector r6 = r5 - r2;

    // Calculate three signed triangle areas, using tFace[0] as the origin.
    scalar t1 = n & (0.5 * (r1 ^ r4));
    scalar t2 = n & (0.5 * (r2 ^ r5));
    scalar t3 = n & (0.5 * (r3 ^ r6));

    // Determine sliver types based on are magnitudes.
    if (t1 > 0 && t2 > 0 && t3 > 0)
    {
        // Region R0: Cap cell.
        map.type() = 2;
        map.apexPoint() = fourthPoint;

        faceToCheck[0] = tFace[0];
        faceToCheck[1] = tFace[1];
        faceToCheck[2] = tFace[2];
    }

    if (t1 < 0 && t2 > 0 && t3 > 0)
    {
        // Region R1: Sliver cell.
        map.type() = 1;

        edgeToCheck[0][0] = tFace[2];
        edgeToCheck[0][1] = fourthPoint;
        edgeToCheck[1][0] = tFace[0];
        edgeToCheck[1][1] = tFace[1];
    }

    if (t1 > 0 && t2 < 0 && t3 > 0)
    {
        // Region R2: Sliver cell.
        map.type() = 1;

        edgeToCheck[0][0] = tFace[0];
        edgeToCheck[0][1] = fourthPoint;
        edgeToCheck[1][0] = tFace[1];
        edgeToCheck[1][1] = tFace[2];
    }

    if (t1 > 0 && t2 > 0 && t3 < 0)
    {
        // Region R3: Sliver cell.
        map.type() = 1;

        edgeToCheck[0][0] = tFace[1];
        edgeToCheck[0][1] = fourthPoint;
        edgeToCheck[1][0] = tFace[2];
        edgeToCheck[1][1] = tFace[0];
    }

    if (t1 < 0 && t2 > 0 && t3 < 0)
    {
        // Region R4: Cap cell.
        map.type() = 2;
        map.apexPoint() = tFace[0];

        faceToCheck[0] = tFace[1];
        faceToCheck[1] = tFace[2];
        faceToCheck[2] = fourthPoint;
    }

    if (t1 < 0 && t2 < 0 && t3 > 0)
    {
        // Region R5: Cap cell.
        map.type() = 2;
        map.apexPoint() = tFace[1];

        faceToCheck[0] = tFace[2];
        faceToCheck[1] = tFace[0];
        faceToCheck[2] = fourthPoint;
    }

    if (t1 > 0 && t2 < 0 && t3 < 0)
    {
        // Region R6: Cap cell.
        map.type() = 2;
        map.apexPoint() = tFace[2];

        faceToCheck[0] = tFace[0];
        faceToCheck[1] = tFace[1];
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

            edgeToCheck[0][0] = tFace[0];
            edgeToCheck[0][1] = fourthPoint;
        }
        else
        if (mag(t2) < refMag)
        {
            // Wedge case: Too close to point [1]
            map.type() = 4;

            edgeToCheck[0][0] = tFace[1];
            edgeToCheck[0][1] = fourthPoint;
        }
        else
        if ((mag(t2) > refMag) && (mag(t3) > refMag))
        {
            // Spade case: Too close to edge vector r1
            map.type() = 3;
            map.apexPoint() = fourthPoint;

            edgeToCheck[0][0] = tFace[0];
            edgeToCheck[0][1] = tFace[1];
        }
    }

    if (mag(t2) < refMag)
    {
        if (mag(t3) < refMag)
        {
            // Wedge case: Too close to point [2]
            map.type() = 4;

            edgeToCheck[0][0] = tFace[2];
            edgeToCheck[0][1] = fourthPoint;
        }
        else
        if ((mag(t1) > refMag) && (mag(t3) > refMag))
        {
            // Spade case: Too close to edge vector r2
            map.type() = 3;
            map.apexPoint() = fourthPoint;

            edgeToCheck[0][0] = tFace[1];
            edgeToCheck[0][1] = tFace[2];
        }
    }

    if (mag(t3) < refMag)
    {
        if ((mag(t1) > refMag) && (mag(t2) > refMag))
        {
            // Spade case: Too close to edge vector r3
            map.type() = 3;
            map.apexPoint() = fourthPoint;

            edgeToCheck[0][0] = tFace[2];
            edgeToCheck[0][1] = tFace[0];
        }
    }

    // Determine appropriate information for sliver exudation.
    switch (map.type())
    {
        case 1:
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

            break;
        }

        case 2:
        {
            // Search the cell-faces for opposing faces.
            forAll(cellToCheck, faceI)
            {
                const face& thisFace = faces_[cellToCheck[faceI]];

                if (triFace::compare(triFace(thisFace), triFace(faceToCheck)))
                {
                    map.opposingFace() = cellToCheck[faceI];

                    break;
                }
            }

            break;
        }

        case 3:
        case 4:
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

            break;
        }

        default:
        {
            WarningIn
            (
                "void dynamicTopoFvMesh::identifySliverType"
                "(const label cIndex) const"
            )
                << nl << "Could not identify sliver type for cell: "
                << cIndex
                << endl;
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


// Remove sliver cells
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
            InfoIn("void dynamicTopoFvMesh::removeSlivers()")
                << nl << "Removing Cell: " << cIndex
                << " of sliver type: " << map.type()
                << " with quality: " << thresholdSlivers_[cIndex]
                << endl;
        }

        // Take action based on the type of sliver.
        switch (map.type())
        {
            case 1:
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
                    const edge& thisEdge =
                    (
                        edges_[firstMapEdges[edgeI][0]]
                    );

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
                        const edge& thisEdge =
                        (
                            edges_[secondMapEdges[edgeI][0]]
                        );

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

                break;
            }

            case 2:
            {
                // Cap cell.
                label opposingFace = map.opposingFace();

                // Force trisection of the opposing face.
                changeMap faceMap =
                (
                    trisectFace(opposingFace, false, true)
                );

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
                        collapseEdge
                        (
                            faceMapEdges[edgeI][0],
                            -1,
                            false,
                            true
                        );

                        break;
                    }
                }

                break;
            }

            case 3:
            {
                // Spade cell.

                // Force bisection on the first edge.
                changeMap firstMap =
                (
                    bisectEdge(map.firstEdge(), false, true)
                );

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
                        collapseEdge
                        (
                            firstMapEdges[edgeI][0],
                            -1,
                            false,
                            true
                        );

                        break;
                    }
                }

                break;
            }

            case 4:
            {
                // Wedge cell.

                // Collapse the first edge.
                collapseEdge
                (
                    map.firstEdge(),
                    -1,
                    false,
                    true
                );

                break;
            }

            default:
            {
                WarningIn("void dynamicTopoFvMesh::removeSlivers()")
                    << nl << "Could not identify sliver type for cell: "
                    << cIndex
                    << endl;
            }
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


// MultiThreaded topology modifier
void dynamicTopoFvMesh::threadedTopoModifier()
{
    // Remove sliver cells first.
    removeSlivers();

    // Handle coupled patches.
    handleCoupledPatches();

    // Set the thread scheduling sequence
    labelList topoSequence(threader_->getNumThreads());

    // Linear sequence from 1 to nThreads
    forAll(topoSequence, indexI)
    {
        topoSequence[indexI] = indexI + 1;
    }

    if (edgeRefinement_)
    {
        // Initialize stacks
        initStacks();

        if (threader_->multiThreaded())
        {
            // Lock slave threads
            lockThreads(topoSequence, handlerPtr_);

            // Submit jobs to the work queue
            forAll(topoSequence, i)
            {
                threader_->addToWorkQueue
                (
                    &edgeRefinementEngine,
                    &(handlerPtr_[topoSequence[i]])
                );

                // Wait for a signal from this thread
                // before moving on.
                handlerPtr_[topoSequence[i]].waitForSignal
                (
                    meshHandler::START
                );
            }

            // Synchronize threads
            synchronizeThreads(topoSequence, handlerPtr_);
        }

        // Set the master thread to implement modifications
        edgeRefinementEngine(&(handlerPtr_[0]));

        // Handle mesh slicing events, if necessary
        handleMeshSlicing();

        if (debug)
        {
            Info << nl << "Edge Bisection/Collapse complete." << endl;
        }
    }

    // Re-Initialize stacks
    initStacks();

    if (threader_->multiThreaded())
    {
        // Lock slave threads
        lockThreads(topoSequence, handlerPtr_);

        // Submit jobs to the work queue
        forAll(topoSequence, i)
        {
            if (twoDMesh_)
            {
                threader_->addToWorkQueue
                (
                    &swap2DEdges,
                    &(handlerPtr_[topoSequence[i]])
                );
            }
            else
            {
                threader_->addToWorkQueue
                (
                    &swap3DEdges,
                    &(handlerPtr_[topoSequence[i]])
                );
            }

            // Wait for a signal from this thread
            // before moving on.
            handlerPtr_[topoSequence[i]].waitForSignal
            (
                meshHandler::START
            );
        }

        // Synchronize threads
        synchronizeThreads(topoSequence, handlerPtr_);
    }

    // Set the master thread to implement modifications
    if (twoDMesh_)
    {
        swap2DEdges(&(handlerPtr_[0]));
    }
    else
    {
        swap3DEdges(&(handlerPtr_[0]));
    }

    if (debug)
    {
        Info << nl << "Edge Swapping complete." << endl;
    }
}


// Reset the mesh and generate mapping information
void dynamicTopoFvMesh::resetMesh()
{
    // Reduce across processors.
    reduce(topoChangeFlag_, orOp<bool>());

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

        clockTime reOrderingTimer;

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

        // Print out topo-stats
        Info << " Reordering time: " << reOrderingTimer.elapsedTime() << endl;

        // Obtain the patch-point maps before resetting the mesh
        List<Map<label> > oldPatchPointMaps(numPatches_);

        forAll(oldPatchPointMaps, patchI)
        {
            oldPatchPointMaps[patchI] = boundaryMesh()[patchI].meshPointMap();
        }

        topoMapper& fieldMapper = mapper_();

        // Set face/cell centres and gradient information
        // for the mapping stage, prior to mesh reset
        fieldMapper.storeCentres();
        fieldMapper.storeGradients();

        // Reset the mesh with pre-motion points
        polyMesh::resetPrimitives
        (
            nFaces_,
            preMotionPoints,
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

        // Update the underlying mesh, and map all related fields
        updateMesh(mpm);

        // If this is the first time-step,
        // perform a dummy movePoints to force V0 creation
        if (time().timeIndex() == 1)
        {
            if (debug > 2)
            {
                InfoIn("void dynamicTopoFvMesh::resetMesh()")
                    << " Setting preMotionPoints for first move."
                    << nl << endl;
            }

            movePoints(mpm.preMotionPoints());
        }

        // Reset old-volumes
        resetMotion();
        setV0();

        // Correct volume fluxes on the old mesh
        fieldMapper.correctFluxes();

        // Now move mesh to new points and
        // compute correct mesh-fluxes.
        movePoints(points);

        // Update the mesh-motion solver
        if (motionSolver_.valid())
        {
            motionSolver_->updateMesh(mpm);
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
                    this
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

        for (label i = 0; i < numPatches_; i++)
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
    else
    {
        // No topology changes were made.
        // Only execute mesh-motion.
        if (motionSolver_.valid())
        {
            movePoints(motionSolver_->curPoints());
        }
    }

    // Dump length-scale to disk, if requested.
    calculateLengthScale(true);

    // Dump procIDs to disk, if requested.
    writeProcIDs();
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

    topoMapper& fieldMapper = mapper_();

    // Set the mapPolyMesh object in the mapper
    fieldMapper.setMapper(meshMap);

    // Set weighting information.
    // This takes over the weight data.
    fieldMapper.setFaceWeights(faceWeights_, faceCentres_);
    fieldMapper.setCellWeights(cellWeights_, cellCentres_);

    // Conservatively map scalar/vector volFields
    fieldMapper.conservativeMapFields<scalar>();
    fieldMapper.conservativeMapFields<vector>();

    // Map all the volFields in the objectRegistry
    MapGeometricFields<sphericalTensor,fvPatchField,topoMapper,volMesh>
        (fieldMapper);
    MapGeometricFields<symmTensor,fvPatchField,topoMapper,volMesh>
        (fieldMapper);
    MapGeometricFields<tensor,fvPatchField,topoMapper,volMesh>
        (fieldMapper);

    // Map all the surfaceFields in the objectRegistry
    MapGeometricFields<scalar,fvsPatchField,topoMapper,surfaceMesh>
        (fieldMapper);
    MapGeometricFields<vector,fvsPatchField,topoMapper,surfaceMesh>
        (fieldMapper);
    MapGeometricFields<sphericalTensor,fvsPatchField,topoMapper,surfaceMesh>
        (fieldMapper);
    MapGeometricFields<symmTensor,fvsPatchField,topoMapper,surfaceMesh>
        (fieldMapper);
    MapGeometricFields<tensor,fvsPatchField,topoMapper,surfaceMesh>
        (fieldMapper);

    // Clear mapper
    fieldMapper.clear();
}


// Update the mesh for motion / topology changes.
//  - Return true if topology changes have occurred
bool dynamicTopoFvMesh::update()
{
    // Re-read options if they have been modified at run-time
    if (dict_.readIfModified())
    {
        // Re-read optional parameters
        readOptionalParameters();

        // Read edge refinement options
        if (edgeRefinement_)
        {
            lengthEstimator().readRefinementOptions(true, mandatory_);
        }
    }

    // Set old point positions
    oldPoints_ = polyMesh::points();

    // Invoke mesh-motion solver and store new points
    if (motionSolver_.valid())
    {
        points_ = motionSolver_->newPoints()();
    }

    // Reset statistics
    topoChangeFlag_ = false;
    nModifications_ = 0;
    nBisections_ = 0;
    nCollapses_ = 0;
    nSwaps_ = 0;

    // Obtain mesh stats before topo-changes
    bool noSlivers = meshQuality(true);

    // Return if the interval is invalid,
    // not at re-mesh interval, or slivers are absent.
    // Handy while using only mesh-motion.
    if (interval_ < 0 || ((time().timeIndex() % interval_ != 0) && noSlivers))
    {
        // Move mesh to new positions
        if (motionSolver_.valid())
        {
            movePoints(motionSolver_->curPoints());
        }

        // Dump procIDs to disk, if requested.
        writeProcIDs();

        // Motion only.
        return false;
    }

    // Calculate the edge length-scale for the mesh
    calculateLengthScale();

    // Track mesh topology modification time
    clockTime topologyTimer;

    // Invoke the threaded topoModifier
    threadedTopoModifier();

    Info << " Topo modifier time: " << topologyTimer.elapsedTime() << endl;

    // Write out statistics
    Info << " Bisections :: Interior: " << nBisections_[0]
         << ", Surface: " << nBisections_[1] << endl;
    Info << " Collapses  :: Interior: " << nCollapses_[0]
         << ", Surface: " << nCollapses_[1] << endl;
    Info << " Swaps      :: Interior: " << nSwaps_[0]
         << ", Surface: " << nSwaps_[1] << endl;

    // Obtain mesh stats after topo-changes
    meshQuality(true);

    // Apply all topology changes (if any) and reset mesh.
    resetMesh();

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
