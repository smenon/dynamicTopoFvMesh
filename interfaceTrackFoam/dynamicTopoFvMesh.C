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
    An implementation of dynamic changes to mesh-topology (w/ smoothing)
 
Author
    Sandeep Menon
\*----------------------------------------------------------------------------*/

#include "dynamicTopoFvMesh.H"
#include "dynamicTopoFvMeshMapper.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynamicTopoFvMesh,0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::dynamicTopoFvMesh::dynamicTopoFvMesh(const IOobject& io)
:
    fvMesh(io),
    numPatches_(this->boundaryMesh().size()),
    topoChangeFlag_(false),
    dict_
    (
        IOobject
        (
        "dynamicMeshDict",
        this->time().constant(),
        (*this),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        )
    ),
    twoDMotion_(this->nGeometricD() == 2 ? true : false),
    edgeModification_(dict_.subDict("dynamicTopoFvMesh").lookup("edgeModification")),
    solveForMotion_(dict_.subDict("dynamicTopoFvMesh").lookup("solveForMotion")),
    mapper_(NULL),
    meshPoints_(this->points()),
    faces_(this->faces()),
    owner_(this->allOwner()),
    neighbour_(this->allNeighbour()),
    cells_(this->cells()),
    oldPatchSizes_(numPatches_,0),
    patchSizes_(numPatches_,0),
    oldPatchStarts_(numPatches_,-1),
    patchStarts_(numPatches_,-1),
    oldPatchNMeshPoints_(numPatches_,-1),
    patchNMeshPoints_(numPatches_,-1),
    nOldPoints_(this->nPoints()),
    nPoints_(this->nPoints()),
    nOldFaces_(this->nFaces()),
    nFaces_(this->nFaces()),
    nOldCells_(this->nCells()),
    nCells_(this->nCells()),
    nOldInternalFaces_(this->nInternalFaces()),
    nInternalFaces_(this->nInternalFaces()), 
    ratioMin_(0.0),
    ratioMax_(0.0),
    growthFactor_(1.0),
    maxLengthScale_(GREAT),
    bisectInteriorFace_(-1)
{
    // Initialize the motion-solver, if it was requested
    if (solveForMotion_)
    {
        motionPtr_.set(motionSolver::New(*this).ptr());
    }
    
    // Initialize the multiThreading environment
    threader_.set
    (
        multiThreader::New
        (
            dict_.subDict("dynamicTopoFvMesh")
        ).ptr()
    );
    
    // Get the number of threads and allocate topoMeshStructures
    label nThreads = threader_->getNumThreads();
    structPtr_ = new topoMeshStruct[nThreads];
    for (label i=0; i<nThreads; i++)
    {
        structPtr_[i].mesh = this;
        structPtr_[i].nThreads = nThreads;
        structPtr_[i].threadID = i;
    }

    // For tetrahedral meshes...
    if (!twoDMotion_)
    {
        // Obtain the tetrahedral metric to be used.
        word tetMetric(dict_.subDict("dynamicTopoFvMesh").lookup("tetMetric"));

        if (tetMetric == "Knupp")
        {
            tetMetric_.set(new Knupp);
        }
        else
        if (tetMetric == "Dihedral")
        {
            tetMetric_.set(new Dihedral);
        }
        else
        {
            FatalErrorIn("dynamicTopoFvMesh::dynamicTopoFvMesh(const IOobject& io) ") << nl
                    << " Unrecognized tet-quality metric: " << tetMetric
                    << abort(FatalError);
        }
    }
    
    // Initialize edge-related connectivity structures
    if (debug) Info << "Building edges..." << endl;
    //initEdges();

    // Define edgeModification options
    if (edgeModification_)
    {
        const dictionary& edgeOptionDict = dict_.subDict("dynamicTopoFvMesh").subDict("edgeOptions");
        ratioMax_ = readScalar(edgeOptionDict.lookup("bisectionRatio"));
        ratioMin_ = readScalar(edgeOptionDict.lookup("collapseRatio"));
        growthFactor_ = readScalar(edgeOptionDict.lookup("growthFactor"));
        
        if (edgeOptionDict.found("maxLengthScale"))
        {
            maxLengthScale_ = readScalar(edgeOptionDict.lookup("maxLengthScale"));
        }

        if (edgeOptionDict.found("fixedLengthScalePatches"))
        {
            fixedLengthScalePatches_ = edgeOptionDict.subDict("fixedLengthScalePatches");
        }

        initEdgeLengths();
    }

    // Set sizes for the reverse maps
    reversePointMap_.setSize(nPoints_);
    reverseFaceMap_.setSize(nFaces_);
    reverseCellMap_.setSize(nCells_);

    // Initialize patch-size information
    const polyBoundaryMesh& boundary = this->boundaryMesh();
    for(label i=0; i<numPatches_; i++)
    {
        oldPatchSizes_[i]  = patchSizes_[i]  = boundary[i].size();
        oldPatchStarts_[i] = patchStarts_[i] = boundary[i].start();
        oldPatchNMeshPoints_[i] = patchNMeshPoints_[i] = boundary[i].meshPoints().size();
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicTopoFvMesh::~dynamicTopoFvMesh()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Access a particular boundary displacement patch
void Foam::dynamicTopoFvMesh::setMotionBC
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
const Foam::autoPtr<mapPolyMesh> Foam::dynamicTopoFvMesh::meshMap()
{
    return mapper_;
}

// Return old cell-centre information (prior to a topology change)
const Foam::vectorField& Foam::dynamicTopoFvMesh::oldCellCentres() const
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
    
    return Foam::vectorField::null();
}

// Return mesh length-scale values
Foam::tmp<volScalarField> Foam::dynamicTopoFvMesh::lengthScale()
{
    tmp<volScalarField> tlengthScale
    (
        new volScalarField
        (
            IOobject
            (
                "lengthScale",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            *this,
            dimensionedScalar("lScale", dimLength, 0),
            zeroGradientFvPatchScalarField::typeName            
        )
    );

    scalarField& internalField = tlengthScale().internalField();
    
    // Obtain length-scale values from the mesh
    forAllIter(HashList<scalar>, lengthScale_, lIter)
    {
        internalField[lIter.index()] = lIter();
    }    
    
    return tlengthScale;
}

// Find the circumcenter, given three points
inline Foam::vector Foam::dynamicTopoFvMesh::circumCenter
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
    FatalErrorIn("dynamicTopoFvMesh::circumCenter(point& a, point& b, point& c) ") << nl
            << " Encountered a co-linear set of points: " << nl
            << " Point a :: " << one << ": " << a << nl
            << " Point b :: " << two << ": " << b << nl
            << " Point c :: " << three << ": " << c << nl
            << abort(FatalError);
#   endif

    return ((c2 + c3)*a + (c3 + c1)*b + (c1 + c2)*c)/(2*cd);
}

// Find the area of a triangle face. This function also assumes face right-handedness
inline Foam::scalar Foam::dynamicTopoFvMesh::triFaceArea
(
    const face& triFace
)
{
    vector v = meshPoints_[triFace[1]] - meshPoints_[triFace[0]];
    vector w = meshPoints_[triFace[2]] - meshPoints_[triFace[0]];

    vector area = (v^w);

    // Return the cross-product of the vectors
    return 0.5 * (area.z());
}

// Method to determine the old boundary patch index for a given face
// Similar to the polyBoundaryMesh routine, but works on local information
inline Foam::label Foam::dynamicTopoFvMesh::whichPatch
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
    // at the end of the list. Check boundaryPatches_ for the patch info
    if (boundaryPatches_.found(index))
        return boundaryPatches_[index];
    else
        FatalErrorIn
            (
            "label dynamicTopoFvMesh::whichPatch(const label& index) const"
            ) << "Cannot find patch information for face index " << index
            << nl
            << " It appears that face ordering is inconsistent with patch information."
            << abort(FatalError);

    return -2;
}

// Utility method to find the interior/boundary faces
// for an input quad-face and adjacent triangle-prism cell.
inline void Foam::dynamicTopoFvMesh::findPrismFaces
(
    const label& findex,
    const cell& c,
    faceList& bdyf,
    labelList& bidx,
    faceList& intf,
    labelList& iidx
)
{
    label indexO=0, indexI=0;

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

// Method to find the interior/boundary faces
// for an input tri-face and adjacent tet/pyramid cell.
inline void Foam::dynamicTopoFvMesh::findTetPyramidFaces
(
    const label& findex,
    const cell& c,
    faceList& bdyf,
    labelList& bidx,
    faceList& intf,
    labelList& iidx
)
{
    forAll(c, i)
    {
        
    }    
}

// Utility method to find the common edge between two faces.
// If an edge is found, returns the common edge on the first face in the argument
bool Foam::dynamicTopoFvMesh::findCommonEdge
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
void Foam::dynamicTopoFvMesh::findIsolatedPoint
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
            "label dynamicTopoFvMesh::findIsolatedPoint(const face&,const edge&,label&,label&)"
            ) << "Cannot find isolated point in face " << f << endl
            << " Using edge: " << e
            << abort(FatalError);
    }
}

// Utility method to replace a face-label in a given cell
inline void Foam::dynamicTopoFvMesh::replaceFaceLabel
(
     const label& original,
     const label& replacement,
     cell& c
)
{
    bool found = false;
    forAll(c, faceI)
    {
        if (c[faceI] == original)
        {
            c[faceI] = replacement;
            found = true;
            break; 
        }
    }
    if (!found)
    {
        FatalErrorIn
            (
            "label dynamicTopoFvMesh::replaceFaceLabel(const label&,const label&,cell&)"
            )   << "Cannot find face " << original << " in cell: " << c << endl
            << " Face: " << replacement << " was not used in replacement."
            << abort(FatalError);
    }
}

// Utility method to replace a point-label in a given face
inline void Foam::dynamicTopoFvMesh::replacePointLabel
(
     const label& original,
     const label& replacement,
     face& f
)
{
    bool found = false;
    forAll(f, pointI)
    {
        if (f[pointI] == original)
        {
            f[pointI] = replacement;
            found = true;
            break;
        }
    }
    if (!found)
    {
        FatalErrorIn
            (
            "label dynamicTopoFvMesh::replacePointLabel(const label&,const label&,face&)"
            )   << "Cannot find point " << original << " in face: " << f << endl
            << " Point: " << replacement << " was not used in replacement."
            << abort(FatalError);
    }
}

// Utility method for face-insertion
Foam::label Foam::dynamicTopoFvMesh::insertFace
(
    const label patch,
    const face& newFace,
    const label newOwner,
    const label newNeighbour,
    const edge& edgeToWatch
)
{
    label newFaceIndex;

    // Append the specified face to each face-related list.
    // This will avoid rehashing of existing structures, but ordering is not maintained
    // Reordering is performed after all pending changes 
    // (flips, bisections, contractions, etc) have been made to the mesh
    newFaceIndex = faces_.append(newFace);
    owner_.append(newOwner);
    neighbour_.append(newNeighbour);
    if (edgeModification_ && twoDMotion_)
    {
        edgeToWatch_.append(edgeToWatch);
    }

    // Keep track of added boundary faces in a separate hash-table
    // This information will be required at the reordering stage
    if (newNeighbour == -1)
    {
        boundaryPatches_.insert(newFaceIndex,patch);
        // Modify patch information for this boundary face
        patchSizes_[patch]++;
        for(label i=patch+1; i<numPatches_; i++)
            patchStarts_[i]++;
    }
    else
    {
        // Increment the number of internal faces, and subsequent patch-starts
        nInternalFaces_++;
        for(label i=0; i<numPatches_; i++)
            patchStarts_[i]++;
    }

    // Increment the total face count
    nFaces_++;

    return newFaceIndex;
}

// Remove the specified face from the mesh
void Foam::dynamicTopoFvMesh::removeFace
(
    const label index
)
{

    if (debug)
    {
        Info << "Removed face: " << index << endl;
        Info << faces_[index] << endl;
    }

    faces_.remove(index);
    owner_.remove(index);
    if (neighbour_[index] == -1)
    {
        // Modify patch information for this boundary face
        label rmFacePatch = whichPatch(index);
        patchSizes_[rmFacePatch]--;
        for(label i=rmFacePatch+1; i<numPatches_; i++)
            patchStarts_[i]--;
    }
    else
    {
        // Decrement the internal face count, and subsequent patch-starts
        nInternalFaces_--;
        for(label i=0; i<numPatches_; i++)
            patchStarts_[i]--;
    }
    neighbour_.remove(index);
    if (edgeModification_ && twoDMotion_)
    {
        edgeToWatch_.remove(index);
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

// Utility method to build a hull of faces/cells that are connected to the edge
// This will also determine whether the edge lies on a boundary
bool Foam::dynamicTopoFvMesh::constructPrismHull
(
    const edge& edgeToCheck,
    const label startFaceIndex,
    DynamicList<label>& hullFaces,
    DynamicList<label>& hullCells,
    DynamicList<label>& hullTriFaces,
    bool requiresTriFaces,
    label& patchID        
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
                if (requiresTriFaces && faceToCheck.nEdges() == 3 && !foundTriFace)
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
                FatalErrorIn
                (
                    "dynamicTopoFvMesh::constructPrismHull(...)"
                ) 
                << nl << " Failed to find a suitable quad/tri face. "
                << " Possibly not a prismatic mesh. " << nl
                << abort(FatalError);
        }
        else
        {
            if (!foundQuadFace)
                FatalErrorIn
                (
                    "dynamicTopoFvMesh::constructPrismHull(...)"
                ) 
                << nl << " Failed to find a suitable quad face. "
                << " Possibly not a prismatic mesh. " << nl
                << abort(FatalError);            
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
            // Determine the patchID
            patchID = whichPatch(faceToExclude);
            break; 
        }
        else
        {
            if (cellIndex != c0)
                hullCells.append(cellIndex);
        }
    } while ( faceToExclude != startFaceIndex );

    if (c1 == -1)
    {
        isBoundary = true;
        // Determine the patchID
        patchID = whichPatch(startFaceIndex);        
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
                        FatalErrorIn
                        (
                            "dynamicTopoFvMesh::constructPrismHull(...)"
                        ) 
                        << nl << " Failed to find a suitable quad/tri face. "
                        << " Possibly not a prismatic mesh. " << nl
                        << abort(FatalError);
                }
                else
                {
                    if (!foundQuadFace)
                        FatalErrorIn
                        (
                            "dynamicTopoFvMesh::constructPrismHull(...)"
                        ) 
                        << nl << " Failed to find a suitable quad face. "
                        << " Possibly not a prismatic mesh. " << nl
                        << abort(FatalError);
                }
#               endif
                // Decide which cell to check next
                if (owner_[faceToExclude] == cellIndex)
                    cellIndex = neighbour_[faceToExclude];
                else
                    cellIndex = owner_[faceToExclude];
                if (cellIndex == -1)
                {
                    break;
                }
                else
                {
                    if (cellIndex != c0) hullCells.append(cellIndex);
                }
            } while ( faceToExclude != startFaceIndex );
            // Add the starting Face index to the list as well
            hullFaces.append(startFaceIndex);
        }
    }

    return isBoundary;
}

// Reorder points after a topology change
void Foam::dynamicTopoFvMesh::reOrderPoints
(
    pointField& points
)
{
    // *** Point renumbering *** //
    // If points were deleted during topology change, the numerical order ceases 
    // to be continuous. Loop through all points and renumber sequentially. 
    // Possible scope for bandwidth-reduction on the motion-solver.

    label pointRenum = 0;

    addedPointRenumbering_.clear();

    HashList<point>::iterator ptIter = meshPoints_.begin();
    while(ptIter != meshPoints_.end())
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

// Reorder faces in upper-triangular order after a topology change
void Foam::dynamicTopoFvMesh::reOrderFaces
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

    label faceInOrder = 0, allFaces = faces_.lastIndex() + 1;
    faceList oldFaces(allFaces);
    labelList oldOwner(allFaces), oldNeighbour(allFaces), visited(allFaces,0);
    edgeList oldEdgeToWatch(0);
    
    HashTable<label,label> addedFaceRenumbering;
    HashTable<label,label> addedFaceReverseRenumbering;

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
    if (edgeModification_ && twoDMotion_)
    {
        oldEdgeToWatch.setSize(allFaces);
        forAllIter(HashList<edge>::iterator, edgeToWatch_, eIter)
            oldEdgeToWatch[eIter.index()] = eIter();
        edgeToWatch_.clear();
    }

    // Mark the internal faces with -2 so that they are inserted first
    forAllIter(HashList<cell>::iterator, cells_, cIter)
    {
        const cell& curFaces = cIter();
        forAll(curFaces, faceI) 
            visited[curFaces[faceI]]--;
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
                        faceRenumber[pointI] = reversePointMap_[faceRenumber[pointI]];
                    else
                        faceRenumber[pointI] = addedPointRenumbering_[faceRenumber[pointI]];
                }

                // Renumber the edges in edgeToWatch
                if (edgeModification_ && twoDMotion_)
                {
                    edge& edgeRenumber = oldEdgeToWatch[curFaces[faceI]];
                    if (edgeRenumber[0] < nOldPoints_)
                        edgeRenumber[0] = reversePointMap_[edgeRenumber[0]];
                    else
                        edgeRenumber[0] = addedPointRenumbering_[edgeRenumber[0]];
                    if (edgeRenumber[1] < nOldPoints_)
                        edgeRenumber[1] = reversePointMap_[edgeRenumber[1]];
                    else
                        edgeRenumber[1] = addedPointRenumbering_[edgeRenumber[1]];
                }

                // Update the maps
                if (curFaces[faceI] < nOldFaces_)
                {
                    faceMap_[bFaceIndex] = curFaces[faceI];
                    reverseFaceMap_[curFaces[faceI]] = bFaceIndex;
                }
                else
                {
                    addedFaceRenumbering.insert(curFaces[faceI],bFaceIndex);
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
                    addedFaceRenumbering.insert(curFaces[nextNei],faceInOrder);
                }

                // Renumber the point labels in this face
                face& faceRenumber = oldFaces[curFaces[nextNei]];
                forAll(faceRenumber, pointI)
                {
                    if (faceRenumber[pointI] < nOldPoints_)
                        faceRenumber[pointI] = reversePointMap_[faceRenumber[pointI]];
                    else
                        faceRenumber[pointI] = addedPointRenumbering_[faceRenumber[pointI]];
                }

                // Renumber owner/neighbour
                label oldOwn = oldOwner[curFaces[nextNei]];
                label oldNei = oldNeighbour[curFaces[nextNei]];
                label ownerRenumber =   oldOwner[curFaces[nextNei]] < nOldCells_
                                      ? reverseCellMap_[oldOwn] : addedCellRenumbering_[oldOwn];
                label neighbourRenumber =   oldNeighbour[curFaces[nextNei]] < nOldCells_
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
                if (edgeModification_ && twoDMotion_)
                {
                    edge& edgeRenumber = oldEdgeToWatch[curFaces[nextNei]];
                    if (edgeRenumber[0] < nOldPoints_)
                        edgeRenumber[0] = reversePointMap_[edgeRenumber[0]];
                    else
                        edgeRenumber[0] = addedPointRenumbering_[edgeRenumber[0]];
                    if (edgeRenumber[1] < nOldPoints_)
                        edgeRenumber[1] = reversePointMap_[edgeRenumber[1]];
                    else
                        edgeRenumber[1] = addedPointRenumbering_[edgeRenumber[1]];
                    edgeToWatch_.append(edgeRenumber);
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
        if (edgeModification_ && twoDMotion_)
            edgeToWatch_.append(oldEdgeToWatch[oldIndex]);

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
                cellFaces[faceI] = addedFaceRenumbering[cellFaces[faceI]];
            }
        }
    }

    // Final check to ensure everything went okay
    if (debug)
    {
        if (sum(visited) != 0)
            FatalErrorIn("Foam::dynamicTopoFvMesh::reOrderFaces()") << nl
                    << " Algorithm did not visit every face in the mesh."
                    << " Something's messed up." << nl
                    << abort(FatalError);         
    }
}

// Reorder & renumber cells with bandwidth reduction after a topology change
void Foam::dynamicTopoFvMesh::reOrderCells()
{
    // *** Cell renumbering *** //
    // If cells were deleted during topology change, the numerical order ceases 
    // to be continuous. Also, cells are always added at the end of the list by 
    // virtue of the HashList append method. Thus, cells would now have to be 
    // reordered so that bandwidth is reduced and renumbered to be sequential.

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
        oldCells[cIter.index()] = cIter();
    cells_.clear();
    if (edgeModification_)
    {
        oldLengthScale.setSize(allCells);
        forAllIter(HashList<scalar>::iterator, lengthScale_, cIter)
            oldLengthScale[cIter.index()] = cIter();
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
        if (ncc[cellI] == 0) visited[cellI] = 1;
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
                        lengthScale_.append(oldLengthScale[currentCell]);

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
            FatalErrorIn("Foam::dynamicTopoFvMesh::reOrderCells()") << nl
                    << " Algorithm did not visit every cell in the mesh."
                    << " Something's messed up." << nl
                    << abort(FatalError);
    }
}

// Reorder the faces in upper-triangular order, and generate mapping information
void Foam::dynamicTopoFvMesh::reOrderMesh
(
    pointField& points,
    faceList& faces,
    labelList& owner,
    labelList& neighbour
)
{
    // Allocate for the mapping information
    pointMap_.setSize(nPoints_, -1);
    faceMap_.setSize(nFaces_, -1);
    cellMap_.setSize(nCells_, -1);

    if (debug)
    {
        Info << endl;
        Info << "=================" << endl;
        Info << " Mesh reOrdering " << endl;
        Info << "=================" << endl;
        Info << "Mesh Info [n]:" << endl;
        Info << "Points: " << nOldPoints_ << endl;
        Info << "Faces: " << nOldFaces_ << endl;
        Info << "Cells: " << nOldCells_ << endl;
        Info << "Internal Faces: " << nOldInternalFaces_ << endl;
        Info << "Patch Starts: " << oldPatchStarts_ << endl;
        Info << "Patch Sizes: " << oldPatchSizes_ << endl;
        Info << "=================" << endl;        
        Info << "Mesh Info [n+1]:" << endl;
        Info << "Points: " << nPoints_ << endl;
        Info << "Faces: " << nFaces_ << endl;
        Info << "Cells: " << nCells_ << endl;
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
}

// Calculate the edge length-scale for the mesh
void Foam::dynamicTopoFvMesh::calculateLengthScale()
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
        const fvBoundaryMesh& bdy = boundary();
        const labelList& own = allOwner();
        const pointField& pList = points();
        forAll(bdy,patchI)
        {
            const polyPatch& bdyPatch = bdy[patchI].patch();
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
                        cellLevels[ownCell] = level;
                        lengthScale[ownCell] = fixedLengthScalePatches_[patchName][0].scalarToken();
                        visitedCells++;
                    }
                    fixed = true; break;
                }
            }
            if (
                    (!fixed)
                 && (bdyPatch.type() != "wedge")
                 && (bdyPatch.type() != "empty")
                 && (bdyPatch.type() != "symmetryPlane")
               )
            {
                label pStart = bdyPatch.start();
                forAll(bdyPatch,faceI)
                {
                    label ownCell = own[pStart+faceI];
                    edge& etw = edgeToWatch_[pStart+faceI];
                    cellLevels[ownCell] = level;
                    lengthScale[ownCell] = mag(pList[etw[0]] - pList[etw[1]]);
                    visitedCells++;
                }
            }
        }

        // Perform multiple sweeps through the mesh...
        while (visitedCells < nCells())
        {
            // Loop through cells, and increment neighbour cells of the current level
            forAll(cellLevels,cellI)
            {
                if (cellLevels[cellI] == level)
                {
                    // Obtain the cells neighbouring this one
                    const labelList& cList = cc[cellI];
                    forAll(cList, indexI) {
                        label& ngbLevel = cellLevels[cList[indexI]];
                        if (ngbLevel == 0)
                        {
                            ngbLevel = level+1;
                            // Compute the mean of the existing neighbour length-scales
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

        // Copy the most recent length-scale values
        forAllIter(HashList<scalar>, lengthScale_, lIter)
        {
            lIter() = lengthScale[lIter.index()];
        }
    }
}

// Return the appropriate length-scale for boundary face
scalar Foam::dynamicTopoFvMesh::boundaryLengthScale
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
        if(boundary()[bFacePatch].name() == patchName)
            return (fixedLengthScalePatches_[patchName][0].scalarToken());
    }

    return lengthScale_[owner_[faceIndex]];
}

// 2D Edge-swapping engine
void Foam::dynamicTopoFvMesh::swap2DEdges(const topoMeshStruct *thread)
{
    bool found, foundinner;
    face f;
    edge firstEdge(0,0);
    labelList otherPointIndex(4), nextToOtherPoint(4);
    labelList commonFaceIndex(4), commonIntFaceIndex(4);
    labelList c0BdyIndex(2), c0IntIndex(2), c1BdyIndex(2), c1IntIndex(2);
    faceList  c0BdyFace(2),  c0IntFace(2),  c1BdyFace(2),  c1IntFace(2);
    faceList  commonFaces(4), commonIntFaces(4);
    edgeList  commonEdges(2);

    // Loop through faces assigned to this thread
    HashList<face>::iterator fBegin = faces_.getIterator(thread->faceStart), fEnd;
    if (thread->threadID != thread->nThreads-1)
        fEnd = faces_.getIterator(thread->faceStart+thread->faceSize);
    
    for(HashList<face>::iterator fIter = fBegin; fIter != fEnd; fIter++) 
    {
        // Retrieve the index for this iterator
        label findex = fIter.index();

        // Reference to this face...
        face& thisFace = fIter();

        // Get the two cells on either side...
        label c0 = owner_[findex], c1 = neighbour_[findex];

        // Consider only internal faces..
        if (c1 == -1) continue;

        // Get cell references
        cell &cell_0 = cells_[c0], &cell_1 = cells_[c1];

        // Consider only triangle-prisms
        if (cell_0.nFaces() > 5 || cell_1.nFaces() > 5) continue;

        // Find the interior/boundary faces.
        findPrismFaces(findex,cell_0,c0BdyFace,c0BdyIndex,c0IntFace,c0IntIndex);
        findPrismFaces(findex,cell_1,c1BdyFace,c1BdyIndex,c1IntFace,c1IntIndex);

//      if (debug) 
//      {
//          Info << "Cell: " << c0 << endl;
//          Info << "Boundary faces: " << c0BdyIndex[0] << ": " << c0BdyFace[0] << endl;
//          Info << "Boundary faces: " << c0BdyIndex[1] << ": " << c0BdyFace[1] << endl;
//          Info << "Interior faces: " << c0IntIndex[0] << ": " << c0IntFace[0] << endl;
//          Info << "Interior faces: " << c0IntIndex[1] << ": " << c0IntFace[1] << endl;
//          Info << "Cell: " << c1 << endl;
//          Info << "Boundary faces: " << c1BdyIndex[0] << ": " << c1BdyFace[0] << endl;
//          Info << "Boundary faces: " << c1BdyIndex[1] << ": " << c1BdyFace[1] << endl;
//          Info << "Interior faces: " << c1IntIndex[0] << ": " << c1IntFace[0] << endl;
//          Info << "Interior faces: " << c1IntIndex[1] << ": " << c1IntFace[1] << endl;
//      }

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
        point& pointZero = meshPoints_[zeroIndex];
        point& pointOne  = meshPoints_[oneIndex];
        point& pointTwo  = meshPoints_[twoIndex];
        point  center    = circumCenter
                           (
                               pointZero, 
                               pointOne, 
                               pointTwo, 
                               zeroIndex, 
                               oneIndex, 
                               twoIndex
                           );
        scalar radius    = (pointZero - center)&(pointZero - center);

        // Find the isolated point on the other face, and the point next to it 
        findIsolatedPoint
        (
            commonFaces[1], 
            commonEdges[0], 
            otherPointIndex[1], 
            nextToOtherPoint[1]
        );

        // ...and determine whether it lies in this circle
        point otherPoint = meshPoints_[otherPointIndex[1]];
        if ( ((otherPoint - center)&(otherPoint - center)) < radius )
        {
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

            // Set the flag
            topoChangeFlag_ = true;

            // Find the other three points that don't lie on shared edges
            // and the points next to them (for orientation)
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
            bool flipOption = false;

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
                    flipOption = true;
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
                    flipOption = true;
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
            flipOption = false;
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
                    flipOption = true;
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
                    flipOption = true;
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
                firstParent = c0;
            else
                firstParent = cellParents_[c0];

            if (c1 < nOldCells_)
                secondParent = c1;
            else
                secondParent = cellParents_[c1];               
            
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
    }
}

// Initialize the edge-length field
void Foam::dynamicTopoFvMesh::initEdgeLengths()
{
    lengthScale_.setSize(nCells_, 0.0);

    if (twoDMotion_)
    {
        // Allocate fields
        edge nullEdge(0,0);  
        edgeToWatch_.setSize(nFaces_,nullEdge);

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

// Initialize mesh edges and related connectivity lists
void Foam::dynamicTopoFvMesh::initEdges()
{
    /*
    label edgeIndex = 0;
    
    // List that keeps track of added edges
    List<labelHashSet> pointPoints(nPoints());
    
    // Initialize faceEdges
    faceEdges_.setSize(nFaces_, SLList<label>());
    
    // Loop through all faces in the mesh and add edges in order
    const faceList& fList = faces();
    forAll(fList, faceI)
    {
        edgeList eList = fList[faceI].edges();
        forAll(eList, edgeI)
        {
            if 
            (
                (!pointPoints[eList[edgeI][0]].found(eList[edgeI][1]))
             && (!pointPoints[eList[edgeI][1]].found(eList[edgeI][0]))
            )
            {
                // Add this edge to the HashList
                edges_.append(eList[edgeI]);
                edgeFaces_.append(SLList<label>());
                
                // Insert edge-face and face-edge
                faceEdges_[faceI].insert(edgeIndex);
                edgeFaces_[edgeIndex].insert(faceI);
                
                // Update the pointPoints list
                pointPoints[eList[edgeI][0]].insert(eList[edgeI][1]);
                pointPoints[eList[edgeI][1]].insert(eList[edgeI][0]);
                edgeIndex++;
            }
        }
    }
    
    Info << "Mesh edges: " << edgeIndex << endl;
    Info << "Mesh edges (polyMesh): " << this->edges().size() << endl;
    */
}

// Does the mesh perform edge-modification?
bool Foam::dynamicTopoFvMesh::edgeModification()
{
    return edgeModification_;
}

// 2D Edge-bisection/collapse engine
void Foam::dynamicTopoFvMesh::edgeBisectCollapse2D
(
    const topoMeshStruct *thread
)
{
    // Loop through all quad-faces and bisect/collapse 
    // edges (quad-faces) by the criterion:
    // Bisect when boundary edge-length > ratioMax_*originalLength
    // Collapse when boundary edge-length < ratioMin_*originalLength
    
    // Loop through faces assigned to this thread
    HashList<face>::iterator fEnd;
    if (thread->threadID != thread->nThreads-1)
        fEnd = faces_(thread->faceStart+thread->faceSize);    

    HashList<face>::iterator fIter = faces_(thread->faceStart);
    while(fIter != fEnd)
    {
        // Retrieve the index for this iterator
        label findex = fIter.index();

        // Reference to this face...
        face& thisFace = fIter();

        // Select only quad-faces
        if (thisFace.size() == 4)
        {
            // Measure the boundary edge-length of the face in question
            edge& checkEdge = edgeToWatch_[findex];
            point& a = meshPoints_[checkEdge[0]];
            point& b = meshPoints_[checkEdge[1]];
            scalar length = mag(b-a);

            // Get the two cells on either side...
            label c0 = owner_[findex], c1 = neighbour_[findex];

            // Determine the length-scale at this face
            scalar scale=0;
            if (c1 == -1)
            {
                scale = boundaryLengthScale(findex);

                // Check if this boundary face is adjacent to a sliver-cell,
                // and remove it by a two-step bisection/collapse operation.
                bool sliverRemoved = remove2DSliver(findex, thisFace);

                if (sliverRemoved)
                {
                    // Set the flag
                    topoChangeFlag_ = true;

                    // Move on to the next face
                    fIter++; 

                    // Remove the temporary interior face
                    removeFace(bisectInteriorFace_);

                    continue;
                }
            }
            else
            {
                scale = 0.5*(lengthScale_[c0]+lengthScale_[c1]);
            }

            //== Edge Bisection ==//
            if(length > ratioMax_*scale)
            {
                // Consider only triangle-prisms
                cell &cell_0 = cells_[c0];
                if (cell_0.nFaces() > 5) continue;
                if (c1 != -1)
                {
                    cell &cell_1 = cells_[c1];
                    if (cell_1.nFaces() > 5) continue;
                }

                // Set the flag
                topoChangeFlag_ = true;

                // Bisect this face
                bisectQuadFace(findex, thisFace);

                // Move on to the next face
                fIter++;
            }
            else
            //== Edge Collapse ==//
            if(length < ratioMin_*scale)
            {
                // Consider only triangle-prisms
                cell &cell_0 = cells_[c0];
                if (cell_0.nFaces() > 5) continue;
                if (c1 != -1)
                {
                    cell &cell_1 = cells_[c1];
                    if (cell_1.nFaces() > 5) continue;
                }

                // Collapse this face
                bool success = collapseQuadFace(findex, thisFace);

                // Increment the iterator to move on to the next face...
                fIter++;

                // The face can safely be deleted, since the iterator points
                // to the next valid face on the face-list.
                if (success)
                {
                    removeFace(findex);

                    // Set the flag
                    topoChangeFlag_ = true;
                }
            }
            else
            {
                // Move on to the next face. Increments are done within the 
                // for-loop, since the face might actually be deleted (due to a collapse) 
                // within the loop-body.
                fIter++;
            }
        }
        else
        {
            // Not a quad-face. Move on to the next one.
            fIter++;
        }
    }
}

// Method for the bisection of a quad-face in 2D
void Foam::dynamicTopoFvMesh::bisectQuadFace
(
    const label findex,
    face& thisFace
)
{
    // Local variables
    bool found;
    label replaceFace, n0=-1, n1=-1;
    edgeList commonEdges(2);
    labelList otherPointIndex(4), nextToOtherPoint(4);
    labelList c0BdyIndex(2), c0IntIndex(2), c1BdyIndex(2), c1IntIndex(2);
    faceList  c0BdyFace(2),  c0IntFace(2),  c1BdyFace(2),  c1IntFace(2);
    edge tmpEdge(0,0), firstEdge(0,0), secondEdge(0,0);

    // Get the two cells on either side...
    label c0 = owner_[findex], c1 = neighbour_[findex];

    // Find the prism faces for cell[0].
    cell &cell_0 = cells_[c0];
    findPrismFaces
    (
        findex,
        cell_0,
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
    cell tmpPrismCell(5);
    label newCellIndex0 = cells_.append(tmpPrismCell);
    cell &newCell0 = cells_[newCellIndex0];

    // Generate mapping information for this new cell
    label firstParent;
    const labelListList& cc = cellCells();
    labelHashSet c0MasterObjects(6);

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
    replaceFaceLabel(c0BdyIndex[1],-1,cell_0);

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
    replacePointLabel
    (
        commonEdges[0].otherVertex(nextToOtherPoint[0]), 
        newPtIndex0, 
        thisFace
    );
    
    replacePointLabel
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
            replaceFaceLabel(c0IntIndex[0],-1,cell_0);
            replaceFace = c0IntIndex[0];
            found = true; break;
        }
    }
    if (!found)
    {
        // The edge was obviously not found before
        replaceFaceLabel(c0IntIndex[1],-1,cell_0);
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
            n0 = neighbour_[replaceFace];
            // This face has to be reversed
            faces_[replaceFace] = faces_[replaceFace].reverseFace();
            owner_[replaceFace] = neighbour_[replaceFace];
            neighbour_[replaceFace] = newCellIndex0;
        }
    }
    else
    {
        n0 = owner_[replaceFace];
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
    replaceFaceLabel(-1, newFaceIndex, newCell0);
    replaceFaceLabel(-1, newFaceIndex, cell_0);

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
    replaceFaceLabel(-1, newFaceIndex, newCell0);

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
    replaceFaceLabel(-1, newFaceIndex, cell_0);

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
        replaceFaceLabel(-1, newFaceIndex, newCell0);

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
            cell_1, 
            c1BdyFace, 
            c1BdyIndex, 
            c1IntFace, 
            c1IntIndex
        );

        // Add a new prism cell to the end of the list
        label newCellIndex1 = cells_.append(tmpPrismCell);
        cell &newCell1 = cells_[newCellIndex1];

        // Generate mapping information for this new cell
        label secondParent;
        labelHashSet c1MasterObjects(6);

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
                replaceFaceLabel(c1IntIndex[0], -1, cell_1);
                replaceFace = c1IntIndex[0];
                found = true; break;
            }
        }
        if (!found)
        {
            // The edge was obviously not found before
            replaceFaceLabel(c1IntIndex[1], -1, cell_1);
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
                n1 = neighbour_[replaceFace];
                // This face has to be reversed
                faces_[replaceFace] = faces_[replaceFace].reverseFace();
                owner_[replaceFace] = neighbour_[replaceFace];
                neighbour_[replaceFace] = newCellIndex1;
            }
        }
        else
        {
            n1 = owner_[replaceFace];
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
        replaceFaceLabel(-1, newFaceIndex, newCell0);
        replaceFaceLabel(-1, newFaceIndex, newCell1);
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
            replaceFaceLabel(c1BdyIndex[0], -1, cell_1);
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
            replaceFaceLabel(c1BdyIndex[1], -1, cell_1);
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
        replaceFaceLabel(-1, newFaceIndex, newCell1);
        replaceFaceLabel(-1, newFaceIndex, cell_1);

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
        replaceFaceLabel(-1, newFaceIndex, cell_1);

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
        replaceFaceLabel(-1, newFaceIndex, newCell1);

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
bool Foam::dynamicTopoFvMesh::collapseQuadFace
(
    const label findex,
    face& thisFace
)
{
    // This face is to be collapsed...
    if (debug) 
    {
        Info << nl << nl 
             << "Face: " << findex << ": " << thisFace 
             << " is to be collapsed. " << endl;
    }

    // Local variables
    labelList c0BdyIndex(2), c0IntIndex(2), c1BdyIndex(2), c1IntIndex(2);
    faceList  c0BdyFace(2),  c0IntFace(2),  c1BdyFace(2),  c1IntFace(2);
    edge  tmpEdge(0, 0), firstEdge, secondEdge;
    face  tmpTriFace(3);

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
    label edgePatch1 = -1, edgePatch2 = -1;
    DynamicList<label> firstEdgeHull(10), secondEdgeHull(10);
    DynamicList<label> firstCells(10),    secondCells(10);
    DynamicList<label> firstTriFaces(10), secondTriFaces(10);
    
    bool firstEdgeBoundary  = constructPrismHull
                              (
                                  firstEdge, 
                                  findex, 
                                  firstEdgeHull, 
                                  firstCells, 
                                  firstTriFaces, 
                                  true,
                                  edgePatch1
                              );

    bool secondEdgeBoundary = constructPrismHull
                              (
                                  secondEdge, 
                                  findex, 
                                  secondEdgeHull, 
                                  secondCells, 
                                  secondTriFaces, 
                                  true,
                                  edgePatch2            
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
    label faceToKeep[2] = {0,0}, faceToThrow[2] = {0,0};
    
    cell &cell_0 = cells_[c0];
    findPrismFaces
    (
        findex, 
        cell_0, 
        c0BdyFace, 
        c0BdyIndex, 
        c0IntFace, 
        c0IntIndex
    );

    if (c1 != -1)
    {
        cell &cell_1 = cells_[c1];
        findPrismFaces
        (
            findex, 
            cell_1, 
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
            << "the face preferentially towards it.\n" << endl
            << "Face: " << findex << ": " << thisFace << endl;
        
        if (boundaryMesh()[edgePatch1].type() == "symmetryPlane")
        {
            secondEdgeBoundary = false;
        }
        else if (boundaryMesh()[edgePatch2].type() == "symmetryPlane")
        {
            firstEdgeBoundary = false;
        }        
    }    

    if (!firstEdgeBoundary && secondEdgeBoundary)
    {
        // Check whether the collapse is possible. 
        forAll(firstTriFaces, indexI)
        {
            if 
            (
                firstTriFaces[indexI] == c0BdyIndex[0] 
             || firstTriFaces[indexI] == c0BdyIndex[1]
            ) continue;
            if (c1 != -1) 
            {
                if 
                (
                    firstTriFaces[indexI] == c1BdyIndex[0] 
                 || firstTriFaces[indexI] == c1BdyIndex[1]
                ) continue;
            }
            face &triFace = faces_[firstTriFaces[indexI]];
            forAll(triFace, pointI)
            {
                tmpTriFace[pointI] = triFace[pointI];
                if (triFace[pointI] == cv0) tmpTriFace[pointI] = cv2;
                if (triFace[pointI] == cv1) tmpTriFace[pointI] = cv3;
            }
            // Compute the area and check if it's zero/negative
            scalar origArea = triFaceArea(triFace);
            scalar newArea  = triFaceArea(tmpTriFace);
            if 
            (
                (Foam::sign(origArea) != Foam::sign(newArea)) 
             || mag(newArea) < VSMALL
            ) return false;
        }
        // Collapse to the second node...
        forAll(firstEdgeHull,faceI)
        {
            face& replacementFace = faces_[firstEdgeHull[faceI]];
            replacePointLabel(cv0,cv2,replacementFace);
            replacePointLabel(cv1,cv3,replacementFace);
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
                secondTriFaces[indexI] == c0BdyIndex[0] 
             || secondTriFaces[indexI] == c0BdyIndex[1]
            ) continue;
            if (c1 != -1)
            {
                if 
                (
                    secondTriFaces[indexI] == c1BdyIndex[0] 
                 || secondTriFaces[indexI] == c1BdyIndex[1]
                ) continue;
            }
            face &triFace = faces_[secondTriFaces[indexI]];
            forAll(triFace, pointI)
            {
                tmpTriFace[pointI] = triFace[pointI];
                if (triFace[pointI] == cv2) tmpTriFace[pointI] = cv0;
                if (triFace[pointI] == cv3) tmpTriFace[pointI] = cv1;
            }
            // Compute the area and check if it's zero/negative
            scalar origArea = triFaceArea(triFace);
            scalar newArea  = triFaceArea(tmpTriFace);
            if 
            (
                (Foam::sign(origArea) != Foam::sign(newArea)) 
             || mag(newArea) < VSMALL
            ) return false;
        }
        // Collapse to the first node by default...
        forAll(secondEdgeHull,faceI)
        {
            face& replacementFace = faces_[secondEdgeHull[faceI]];
            replacePointLabel(cv2, cv0, replacementFace);
            replacePointLabel(cv3, cv1, replacementFace);
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
    label cellCheck[2] = {0,0};
    if (owner_[faceToThrow[0]] == c0)
    {
        cellCheck[0] = neighbour_[faceToThrow[0]];
        if (owner_[faceToKeep[0]] == c0)
        {
            if (
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
                // Keep orientation intact, and update the owner
                owner_[faceToKeep[0]] = neighbour_[faceToThrow[0]];
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
                if (
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
                    // Keep orientation intact, and update the owner
                    owner_[faceToKeep[1]] = neighbour_[faceToThrow[1]];
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

    // Remove the unwanted faces in the cell(s) adjacent to this face,
    // and correct the cells that contain discarded faces
    forAll(cell_0,faceI)
    {
        if (cell_0[faceI] != findex && cell_0[faceI] != faceToKeep[0])
           removeFace(cell_0[faceI]);
    }
    cells_.remove(c0);
    lengthScale_.remove(c0);
    replaceFaceLabel(faceToThrow[0], faceToKeep[0], cells_[cellCheck[0]]);
    
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
        cell &cell_1 = cells_[c1];
        forAll(cell_1, faceI)
        {
            if (cell_1[faceI] != findex && cell_1[faceI] != faceToKeep[1])
               removeFace(cell_1[faceI]);
        }
        cells_.remove(c1);
        lengthScale_.remove(c1);
        replaceFaceLabel(faceToThrow[1], faceToKeep[1], cells_[cellCheck[1]]);

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

    // Return a successful collapse
    return true;
}

// Check if the boundary face is adjacent to a sliver-cell,
// and remove it by a two-step bisection/collapse operation.
bool Foam::dynamicTopoFvMesh::remove2DSliver
(
    const label findex,
    face& thisFace
)
{
    labelList c0BdyIndex(2), c0IntIndex(2);
    faceList  c0BdyFace(2),  c0IntFace(2);

    // Measure the boundary edge-length of the face in question
    edge& checkEdge = edgeToWatch_[findex];
    point& a = meshPoints_[checkEdge[0]];
    point& b = meshPoints_[checkEdge[1]];
    scalar length = mag(b-a);

    // Determine the neighbouring cell
    label c0 = owner_[findex];

    // Find the prism-faces
    cell &cell_0 = cells_[c0];
    findPrismFaces(findex, cell_0, c0BdyFace, c0BdyIndex, c0IntFace, c0IntIndex);

    // Determine the boundary triangular face area
    scalar area = triFaceArea(c0BdyFace[0]);

    // This cell has to be removed...
    if (mag(area) < (0.05*length*length))
    {
        // Step 1: Bisect the boundary quad face
        bisectInteriorFace_ = -1;
        bisectQuadFace(findex, thisFace);

        // Step 2: Collapse the newly created internal quad face
        face& newFace = faces_[bisectInteriorFace_];
        bool success = collapseQuadFace(bisectInteriorFace_, newFace);

        // If this failed, something is terribly wrong
        if (!success)
        {
            FatalErrorIn("Foam::dynamicTopoFvMesh::remove2DSliver(const label, face&)")
                << "Attempt to remove sliver cell: "
                << c0 << ": " << cell_0
                << " failed. Possibly a highly skewed mesh."
                << abort(FatalError);
        }
        return true;
    }

    return false;
}

// Prepare thread structures
void Foam::dynamicTopoFvMesh::prepareThreads()
{
    // Fill in required thread-info
    label numThreads = threader_().getNumThreads(); 
    label pointsPerBlock = (nPoints_/numThreads)+1;
    label facesPerBlock  = (nFaces_/numThreads)+1;
    label cellsPerBlock  = (nCells_/numThreads)+1;

    // Information for the master thread
    structPtr_[0].pointSize  = pointsPerBlock;
    structPtr_[0].faceSize   = facesPerBlock;
    structPtr_[0].cellSize   = cellsPerBlock;
    threader_().setData(0, &(structPtr_[0]));

    // Information for all subsequent threads
    label pointsLeft = nPoints_, facesLeft = nFaces_, cellsLeft = nCells_;
    for (label i = 1; i < numThreads; i++)
    {
        pointsLeft -= pointsPerBlock;
        facesLeft -= facesPerBlock;
        cellsLeft -= cellsPerBlock;
        structPtr_[i].pointStart = 
            structPtr_[i-1].pointStart + structPtr_[i-1].pointSize;
        structPtr_[i].faceStart = 
            structPtr_[i-1].faceStart + structPtr_[i-1].faceSize;
        structPtr_[i].cellStart =
            structPtr_[i-1].cellStart + structPtr_[i-1].cellSize;
        structPtr_[i].pointSize = 
            (pointsLeft < pointsPerBlock) ? pointsLeft : pointsPerBlock;
        structPtr_[i].faceSize =
            (facesLeft < facesPerBlock) ? facesLeft : facesPerBlock;
        structPtr_[i].cellSize =
            (cellsLeft < cellsPerBlock) ? cellsLeft : cellsPerBlock;
        threader_().setData(i, &(structPtr_[i]));
    }    
}

// Update mesh corresponding to the given map
void Foam::dynamicTopoFvMesh::updateMesh(const mapPolyMesh& mpm)
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
void Foam::dynamicTopoFvMesh::mapFields(const mapPolyMesh& meshMap)
{
    if (debug)
    {
        Info << "void dynamicTopoFvMesh::mapFields(const mapPolyMesh& meshmMap): "
             << "Mapping fvFields."
             << endl;
    }
   
    //- Field mapping class
    dynamicTopoFvMeshMapper fieldMapper(*this,meshMap);

    // Map all the volFields in the objectRegistry
    MapGeometricFields<scalar, fvPatchField, dynamicTopoFvMeshMapper, volMesh>
        (fieldMapper);
    MapGeometricFields<vector, fvPatchField, dynamicTopoFvMeshMapper, volMesh>
        (fieldMapper);

    // Map all the surfaceFields in the objectRegistry
    MapGeometricFields<scalar, fvsPatchField, dynamicTopoFvMeshMapper, surfaceMesh>
        (fieldMapper);
     
    // Old volumes are not mapped since interpolation is 
    // performed at the same time level.
}

// Update the mesh for motion
// This routine assumes that all boundary motions have been defined
// and incorporated into the mesh for the current time-step.
void Foam::dynamicTopoFvMesh::updateMotion()
{
    if (solveForMotion_)
    {    
        // Solve for motion   
        movePoints(motionPtr_->newPoints());
    }
}

// MultiThreaded topology modifier [2D]
threadReturnType Foam::dynamicTopoFvMesh::threadedTopoModifier2D(void *arg)
{
    // Typecast the argument into the required structure
    multiThreader::threadInfo *threadStruct =
        reinterpret_cast<multiThreader::threadInfo *>(arg);            
    
    topoMeshStruct *meshStructure = 
        reinterpret_cast<topoMeshStruct*>(threadStruct->methodArg);
    
    dynamicTopoFvMesh *mesh = meshStructure->mesh;
    
    // 2D Edge-bisection/collapse engine
    if (mesh->edgeModification()) mesh->edgeBisectCollapse2D(meshStructure);

    if (debug) Info << nl << "2D Edge Bisection/Collapse complete." << endl;

    // 2D Edge-swapping engine
    mesh->swap2DEdges(meshStructure);

    if (debug) Info << nl << "2D Edge Swapping complete." << endl;    
    
    return threadReturnValue;
}

// Update the mesh for topology changes
// Return true if changes have occurred
bool Foam::dynamicTopoFvMesh::updateTopology()
{
    // Calculate the edge length-scale for the mesh
    calculateLengthScale();

    // Keep a copy of existing sizes
    nOldPoints_ = nPoints_;
    nOldFaces_  = nFaces_;
    nOldCells_  = nCells_;
    nOldInternalFaces_ = nInternalFaces_;
    for(label i=0; i<numPatches_; i++)
    {
        oldPatchSizes_[i] = patchSizes_[i];
        oldPatchStarts_[i] = patchStarts_[i];
        oldPatchNMeshPoints_[i] = patchNMeshPoints_[i];
    }

    // Print out the mesh bandwidth
    if (debug)
    {
        label band=0;
        const labelList& oldOwner = owner();
        const labelList& oldNeighbour = neighbour();
        forAll(oldOwner, faceI)
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
        // Update the iterators
        pIter++;
    }

    //== Connectivity changes ==//

    // Reset the flag
    topoChangeFlag_ = false;

    if ( twoDMotion_ )
    {        
        // Set the mesh-modifier static method
        threader_().setMethod(threadedTopoModifier2D);

        // Prepare threads for execution
        prepareThreads();
        
        // Start all threads
        threader_().executeThreads();
    }
    else
    {
        if (debug) Info << nl << "3D Edge Bisection/Collapse complete." << endl;

        if (debug) Info << nl << "3D Edge Swapping complete." << endl;
    }

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
            oldMeshPointLabels[i] = boundaryMesh()[i].meshPoints();
        
        // Obtain cell-centre information as well
        cellCentresPtr_.set
        (
            new vectorField(polyMesh::cellCentres())
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
                            if (oldMeshPointLabels[i][pointJ] == meshPointLabels[pointI])
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
                oldPatchNMeshPoints_
            )
        );

        // Update the underlying mesh, and map all related fields
        updateMesh(mapper_);
        
        // Discard old cell-centre information after mapping
        cellCentresPtr_.clear();

        // Update the motion-solver, if necessary
        if (motionPtr_.valid()) motionPtr_().updateMesh(mapper_);

        // Print out the mesh bandwidth
        if (debug)
        {
            label band=0;
            forAll(owner, faceI)
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
        faceMap_.clear();
        cellMap_.clear();
        reversePointMap_.clear();
        reverseFaceMap_.clear();
        reverseCellMap_.clear();
        boundaryPatches_.clear();
        cellsFromCells_.clear();
        cellParents_.clear();

        // Set new sizes for the reverse maps
        reversePointMap_.setSize(nPoints_);
        reverseFaceMap_.setSize(nFaces_);
        reverseCellMap_.setSize(nCells_);

        // Basic checks for mesh-validity
        if (debug) checkMesh(true);
    }

    return topoChangeFlag_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::dynamicTopoFvMesh::operator=(const dynamicTopoFvMesh& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::dynamicTopoFvMesh::operator=(const Foam::dynamicTopoFvMesh&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
