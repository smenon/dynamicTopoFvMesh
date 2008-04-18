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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::dynamicTopoFvMesh::dynamicTopoFvMesh(const IOobject& io) 
:  
    fvMesh(io),
    numPatches_(this->boundaryMesh().size()),
    topoChangeFlag_(false),                 
    dict_(	    
        IOobject
        (
        "dynamicMeshDict",
        this->time().constant(),
        (*this),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        )
    ),            
    twoDMotion_(dict_.subDict("dynamicTopoFvMesh").lookup("twoD")), 
    edgeModification_(dict_.subDict("dynamicTopoFvMesh").lookup("edgeModification")),
    solveForMotion_(dict_.subDict("dynamicTopoFvMesh").lookup("solveForMotion")),
    fluxInterpolation_(dict_.subDict("dynamicTopoFvMesh").lookup("fluxInterpolation")),
    mapper_(NULL),
    meshPoints_(this->points()),
    pointsZeroVol_(this->points()),
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
    nInternalEdges_(0), 
    ratioMin_(0.0),
    ratioMax_(0.0),
    growthFactor_("growthFactor",dimensionSet(0,2,-1,0,0),0.0),
    debug(false)
{
    // Initialize the motion-solver, if it was requested
    if (solveForMotion_) {
        motionPtr_.set(motionSolver::New(*this).ptr());
    }
    
    // Obtain the name of the registered flux-field
    if (fluxInterpolation_) {        
        fluxFieldName_ = word(dict_.subDict("dynamicTopoFvMesh").lookup("fluxField"));
    }
    
    // For tetrahedral meshes...
    if (!twoDMotion_) {
        // Obtain the tetrahedral metric to be used.        
        word tetMetric(dict_.subDict("dynamicTopoFvMesh").lookup("tetMetric"));
        if (tetMetric == "Knupp") {
            tetMetric_.set(new Knupp);
        } else if (tetMetric == "Dihedral") {
            tetMetric_.set(new Dihedral);
        } else {
            FatalErrorIn("dynamicTopoFvMesh::dynamicTopoFvMesh(const IOobject& io) ") << nl
                    << " Unrecognized tet-quality metric: " << tetMetric
                    << abort(FatalError);
        }
        
        // Initialize internal edges for swapping
        initInternalEdges();
    }    
        
    // Define edgeModification options
    if (edgeModification_) {
        const dictionary& edgeOptionDict = dict_.subDict("dynamicTopoFvMesh").subDict("edgeOptions");
        ratioMax_ = readScalar(edgeOptionDict.lookup("bisectionRatio"));
        ratioMin_ = readScalar(edgeOptionDict.lookup("collapseRatio"));
        growthFactor_.value() = readScalar(edgeOptionDict.lookup("growthFactor"));
        initEdgeLengths();
    }
    
    // Set sizes for the reverse maps
    reversePointMap_.setSize(nPoints_);
    reverseFaceMap_.setSize(nFaces_);
    reverseCellMap_.setSize(nCells_);
        
    // Create a displacement field for all boundaries defined in the mesh
    // and initialize values to zero
    const polyBoundaryMesh& boundary = this->boundaryMesh();
    displacementPtr_.setSize(numPatches_);
    for(label i=0; i<numPatches_; i++) {
        displacementPtr_.set(
           i,
           new vectorField(boundary[i].nPoints(), vector::zero)
        ); 
        oldPatchSizes_[i]  = patchSizes_[i]  = boundary[i].size();
        oldPatchStarts_[i] = patchStarts_[i] = boundary[i].start();
        oldPatchNMeshPoints_[i] = patchNMeshPoints_[i] = boundary[i].meshPoints().size();
    }    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicTopoFvMesh::~dynamicTopoFvMesh()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Access
Foam::vectorField& Foam::dynamicTopoFvMesh::boundaryDisplacementPatch(label& index) {
    // Reallocate if the number of boundary points have changed
    if (displacementPtr_[index].size() != this->boundaryMesh()[index].nPoints()) {
        displacementPtr_[index].clear();
        displacementPtr_[index].setSize(this->boundaryMesh()[index].nPoints(), vector::zero);        
    }
    return displacementPtr_[index];
}

// Find the circumcenter, given three points
inline Foam::vector Foam::dynamicTopoFvMesh::circumCenter(point& a, point& b, point& c, label& one, label& two, label& three)
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

// Find the cell center
inline Foam::vector Foam::dynamicTopoFvMesh::cellCenter(const cell& checkCell)
{
    vector fC, cC = vector::zero;
    forAll(checkCell,faceI) {
        face& faceCheck = faces_[checkCell[faceI]];
        fC = vector::zero;
        forAll(faceCheck,pointI) {
            fC += meshPoints_[faceCheck[pointI]];
        }
        cC += (fC/faceCheck.size());
    }
    return (cC/checkCell.size());
}

// Find the area of a triangle face. This function also assumes face right-handedness
inline Foam::scalar Foam::dynamicTopoFvMesh::triFaceArea(const face& triFace)
{
    vector v = meshPoints_[triFace[1]] - meshPoints_[triFace[0]];
    vector w = meshPoints_[triFace[2]] - meshPoints_[triFace[0]];
    
    vector area = (v^w);
    
    // Return the cross-product of the vectors
    return 0.5 * (area.z());
}

// Method to determine the old boundary patch index for a given face
// Similar to the polyBoundaryMesh routine, but works on local information
inline Foam::label Foam::dynamicTopoFvMesh::whichPatch(const label& index) const
{    
    if (index < nInternalFaces_) return -1;
    
    for(label i=0; i<numPatches_; i++) {
        if
        (
            index >= patchStarts_[i]
         && index < patchStarts_[i] + patchSizes_[i]
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
            ) << "Cannot find patch information for face index " << index << ": " << faces_[index] 
            << nl << " It appears that face ordering is inconsistent with patch information."    
            << abort(FatalError);        
    
    return -2;    
}

// Utility method to find the interior/boundary faces
// for an input quad-face and adjacent triangle-prism cell.
void Foam::dynamicTopoFvMesh::findPrismFaces(
     const label& findex, 
     const cell& c, 
     face bdyf[], 
     label bidx[], 
     face intf[], 
     label iidx[]
)
{
    label indexO=0, indexI=0;
    
    forAll(c, i) {
       label faceIndex = c[i];
       face& fi=faces_[faceIndex];
       if (neighbour_[faceIndex] == -1) {
           if (fi.size() == 3) {
               // Triangular face on the boundary
               bidx[indexO] = faceIndex;
               bdyf[indexO++] = fi;
           } else {
               // This seems to be a non-triangular face on the boundary
               // Consider this as "interior" and move on
               // (Don't count the face under consideration)
               if (faceIndex != findex) {
                   iidx[indexI] = faceIndex;
                   intf[indexI++] = fi;
               }
           }
       } else {
           // Face on the interior (Don't count the face under consideration)
           if (faceIndex != findex) {
               iidx[indexI] = faceIndex;
               intf[indexI++] = fi;
           }
       }
    }        
}

// Utility method to find the common edge between two faces.
// If an edge is found, returns the common edge on the first face in the argument
bool Foam::dynamicTopoFvMesh::findCommonEdge(const face& first, const face& second, edge& common)
{
    bool found=false;
    edgeList efi = first.edges();
    edgeList efj = second.edges();
    forAll(efi, edgeI) {
        forAll(efj, edgeJ) {
            if (efi[edgeI] == efj[edgeJ]) {
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
void Foam::dynamicTopoFvMesh::findIsolatedPoint(
     const face& f, 
     const edge& e, 
     label& ptIndex, 
     label& nextPtIndex
)
{
    bool found = false;
    forAll(f, pointI) {
        if ( f[pointI] != e[0] && f[pointI] != e[1] ) {
            ptIndex = f[pointI];
            nextPtIndex = f[(pointI+1)%3];
            found = true;
            break;
        }
    }   
    if (!found) {
        FatalErrorIn
            (
            "label dynamicTopoFvMesh::findIsolatedPoint(const face&,const edge&,label&,label&)"
            )   << "Cannot find isolated point in face " << f << endl
            << " Using edge: " << e
            << abort(FatalError);
    }    
}

// Utility method to replace a face-label in a given cell
inline void Foam::dynamicTopoFvMesh::replaceFaceLabel(
     const label& original, 
     const label& replacement, 
     cell& c
)
{
    bool found = false;
    forAll(c, faceI) 
        if (c[faceI] == original) { 
            c[faceI] = replacement; 
            found = true;
            break; 
        }
    if (!found) {
        FatalErrorIn
            (
            "label dynamicTopoFvMesh::replaceFaceLabel(const label&,const label&,cell&)"
            )   << "Cannot find face " << original << " in cell: " << c << endl
            << " Face: " << replacement << " was not used in replacement."
            << abort(FatalError);
    }
}

// Utility method to replace a point-label in a given face
inline void Foam::dynamicTopoFvMesh::replacePointLabel(
     const label& original, 
     const label& replacement, 
     face& f
)
{
    bool found = false;
    forAll(f, pointI) 
        if (f[pointI] == original) { 
            f[pointI] = replacement; 
            found = true;
            break; 
        }
    if (!found) {
        FatalErrorIn
            (
            "label dynamicTopoFvMesh::replacePointLabel(const label&,const label&,face&)"
            )   << "Cannot find point " << original << " in face: " << f << endl
            << " Point: " << replacement << " was not used in replacement."
            << abort(FatalError);        
    }    
}

// Utility method for face-insertion
Foam::label Foam::dynamicTopoFvMesh::insertFace(
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
    if (edgeModification_ && twoDMotion_) {
        edgeToWatch_.append(edgeToWatch);
    }
    if (fluxInterpolation_) {
        localPhi_.append(0.0);
    }
    
    // Keep track of added boundary faces in a separate hash-table
    // This information will be required at the reordering stage
    if (newNeighbour == -1) {
        boundaryPatches_.insert(newFaceIndex,patch);
        // Modify patch information for this boundary face
        patchSizes_[patch]++;
        for(label i=patch+1; i<numPatches_; i++)
            patchStarts_[i]++;
    } else {
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
void Foam::dynamicTopoFvMesh::removeFace(const label index) {  
    
    if (debug) {
        Info << "Removed face: " << index << endl;
        Info << faces_[index] << endl;
    }
    
    faces_.remove(index);
    owner_.remove(index);
    if (neighbour_[index] == -1) {
        // Modify patch information for this boundary face
        label rmFacePatch = whichPatch(index);
        patchSizes_[rmFacePatch]--;
        for(label i=rmFacePatch+1; i<numPatches_; i++)
            patchStarts_[i]--;
    } else {
        // Decrement the internal face count, and subsequent patch-starts
        nInternalFaces_--;
        for(label i=0; i<numPatches_; i++)
            patchStarts_[i]--;
    }
    neighbour_.remove(index);
    if (edgeModification_ && twoDMotion_) {
        edgeToWatch_.remove(index);
    }  
    if (fluxInterpolation_) {
        localPhi_.remove(index);
    }
    
    // Update the reverse face-map, but only if this is a face that existed
    // at time [n]. Added faces which are deleted during the topology change
    // needn't be updated.
    if (index < nOldFaces_) reverseFaceMap_[index] = -1;
    
    // Decrement the total face-count
    nFaces_--;
}

// Utility method to build a hull of faces/cells that are connected to the edge
// This will also determine whether the edge lies on a boundary
bool Foam::dynamicTopoFvMesh::constructPrismHull(
    const edge& edgeToCheck, 
    const label startFaceIndex, 
    DynamicList<label>& hullFaces,
    DynamicList<label>& hullCells,
    DynamicList<label>& hullTriFaces,
    bool requiresTriFaces
)
{
    // Get the two cells on either side...
    label c0 = owner_[startFaceIndex], c1 = neighbour_[startFaceIndex];
    
    bool isBoundary=false, foundQuadFace, foundTriFace;
    label faceToExclude, cellIndex;
    
    // Start a search from cell[0] and add to the list as we go along
    faceToExclude=startFaceIndex, cellIndex=c0, hullCells.append(c0);
    do {
        cell& cellToCheck = cells_[cellIndex];
        foundQuadFace = false; foundTriFace = false;
        forAll(cellToCheck,faceI) {
            if (cellToCheck[faceI] != faceToExclude) {
                face& faceToCheck = faces_[cellToCheck[faceI]];
                if (faceToCheck.nEdges() == 4 && !foundQuadFace) {
                    edgeList indexEdges = faceToCheck.edges();
                    forAll(indexEdges,edgeI) {
                        if (indexEdges[edgeI] == edgeToCheck) {
                            // Bingo... We have a match. Add to the dynamic list
                            hullFaces.append(cellToCheck[faceI]);
                            faceToExclude = cellToCheck[faceI];
                            foundQuadFace=true; break;
                        }  
                    }
                }
                if (requiresTriFaces && faceToCheck.nEdges() == 3 && !foundTriFace) {
                    hullTriFaces.append(cellToCheck[faceI]);
                    foundTriFace=true;
                }
            }
            // Found the faces we were looking for, break-out
            if (requiresTriFaces) {
                if (foundQuadFace && foundTriFace) break;
            } else {
                if (foundQuadFace) break;
            }
        }
#       ifdef FULLDEBUG
        if (requiresTriFaces) {
            if (!foundQuadFace || !foundTriFace)
                FatalErrorIn("dynamicTopoFvMesh::constructPrismHull(...)") << nl
                        << " Failed to find a suitable quad/tri face. Possibly not a prismatic mesh. " << nl
                        << abort(FatalError);
        } else {
            if (!foundQuadFace)
                FatalErrorIn("dynamicTopoFvMesh::constructPrismHull(...)") << nl
                        << " Failed to find a suitable quad face. Possibly not a prismatic mesh. " << nl
                        << abort(FatalError);            
        }
#       endif
        // Decide which cell to check next
        if (owner_[faceToExclude] == cellIndex)
            cellIndex = neighbour_[faceToExclude];
        else
            cellIndex = owner_[faceToExclude];
        if (cellIndex == -1) {
            isBoundary=true; break; 
        } else {
            if (cellIndex != c0) hullCells.append(cellIndex);
        }
    } while ( faceToExclude != startFaceIndex );
    
    if (c1 == -1) 
        isBoundary = true;
    else {
        // Check if the previous search hit a boundary. 
        // If yes, start another search in the reverse direction.
        if (isBoundary) {
            // Start a search from cell[1] and add to the list as we go along
            faceToExclude=startFaceIndex, cellIndex=c1, hullCells.append(c1);
            do {
                cell& cellToCheck = cells_[cellIndex];
                foundQuadFace = false; foundTriFace = false;
                forAll(cellToCheck, faceI) {
                    if (cellToCheck[faceI] != faceToExclude) {
                        face& faceToCheck = faces_[cellToCheck[faceI]];
                        if (faceToCheck.nEdges() == 4 && !foundQuadFace) {
                            edgeList indexEdges = faceToCheck.edges();
                            forAll(indexEdges, edgeI) {
                                if (indexEdges[edgeI] == edgeToCheck) {
                                    // Bingo... We have a match. Add to the dynamic list
                                    hullFaces.append(cellToCheck[faceI]);
                                    faceToExclude = cellToCheck[faceI];
                                    foundQuadFace=true; break;
                                }
                            }
                        }
                        if (requiresTriFaces && faceToCheck.nEdges() == 3 && !foundTriFace) {
                            hullTriFaces.append(cellToCheck[faceI]);
                            foundTriFace=true;
                        }
                    }
                    // Found the faces we were looking for, break-out
                    if (requiresTriFaces) {
                        if (foundQuadFace && foundTriFace) break;
                    } else {
                        if (foundQuadFace) break;
                    }
                }
#               ifdef FULLDEBUG
                if (requiresTriFaces) {
                    if (!foundQuadFace || !foundTriFace)
                        FatalErrorIn("dynamicTopoFvMesh::constructPrismHull(...)") << nl
                                << " Failed to find a suitable quad/tri face. Possibly not a prismatic mesh. " << nl
                                << abort(FatalError);
                } else {
                    if (!foundQuadFace)
                        FatalErrorIn("dynamicTopoFvMesh::constructPrismHull(...)") << nl
                                << " Failed to find a suitable quad face. Possibly not a prismatic mesh. " << nl
                                << abort(FatalError);            
                }
#               endif                
                // Decide which cell to check next
                if (owner_[faceToExclude] == cellIndex)
                    cellIndex = neighbour_[faceToExclude];
                else
                    cellIndex = owner_[faceToExclude];
                if (cellIndex == -1) {
                    break;
                } else {
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
void Foam::dynamicTopoFvMesh::reOrderPoints(pointField& points, pointField& pointsZeroVolume)
{
    // *** Point renumbering *** //
    // If points were deleted during topology change, the numerical order ceases to be continuous.
    // Loop through all points and renumber sequentially. Possible scope for bandwidth-reduction
    // on the motion-solver.
    
    label pointRenum = 0;
    
    addedPointRenumbering_.clear();
    
    HashList<point>::iterator ptIter = meshPoints_.begin();
    HashList<point>::iterator pzvIter = pointsZeroVol_.begin();
    while(ptIter != meshPoints_.end()) {       
        // Obtain the index for this point
        label pIndex = ptIter.index();      
        // Update the point info
        points[pointRenum] = ptIter();
        pointsZeroVolume[pointRenum] = pzvIter();
        // Renumber the point index
        meshPoints_.reNumber(pointRenum, ptIter);
        pointsZeroVol_.reNumber(pointRenum,pzvIter);
        // Added points are always numbered after nOldPoints_ 
        // (by virtue of the HashList append method)
        if (pIndex < nOldPoints_) {
            pointMap_[pointRenum]    = pIndex;
            reversePointMap_[pIndex] = pointRenum;
        } else {
            addedPointRenumbering_.insert(pIndex,pointRenum);
        }
        // Update the counter
        pointRenum++;
        // Update the iterators
        ptIter++; pzvIter++;
    }    
}

// Reorder faces in upper-triangular order after a topology change
void Foam::dynamicTopoFvMesh::reOrderFaces(faceList& faces, labelList& owner, labelList& neighbour)
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
    scalarField oldPhi(0);
    
    HashTable<label,label> addedFaceRenumbering;
    HashTable<label,label> addedReverseFaceRenumbering;
    
    // Make a copy of the old face-based HashLists, and clear them
    HashList<face>::iterator fIter = faces_.begin();
    HashList<label>::iterator oIter = owner_.begin();
    HashList<label>::iterator nIter = neighbour_.begin();
    while(fIter != faces_.end()) {
        oldFaces[fIter.index()] = fIter();
        oldOwner[oIter.index()] = oIter();
        oldNeighbour[nIter.index()] = nIter();
        fIter++; oIter++; nIter++;
    }
    faces_.clear(); owner_.clear(); neighbour_.clear();
    if (edgeModification_ && twoDMotion_) {
        oldEdgeToWatch.setSize(allFaces);
        for(HashList<edge>::iterator eIter = edgeToWatch_.begin(); eIter != edgeToWatch_.end(); eIter++)
            oldEdgeToWatch[eIter.index()] = eIter();
        edgeToWatch_.clear();
    }
    if (fluxInterpolation_) {
        oldPhi.setSize(allFaces);
        for(HashList<scalar>::iterator pIter = localPhi_.begin(); pIter != localPhi_.end(); pIter++)
            oldPhi[pIter.index()] = pIter();
        localPhi_.clear();        
    }
    
    // Mark the internal faces with -2 so that they are inserted first
    for(HashList<cell>::iterator cIter = cells_.begin(); cIter != cells_.end(); cIter++) {
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
    for(HashList<cell>::iterator cIter = cells_.begin(); cIter != cells_.end(); cIter++) {
        
        // Record the neighbour cell
        label cellI = cIter.index();
        const cell& curFaces = cIter();
        labelList neiCells(curFaces.size(), -1);
        
        label nNeighbours = 0;
        
        forAll (curFaces, faceI) {
            
            if (visited[curFaces[faceI]] == -2) {
                
                // Face is internal and gets reordered
                label own =   oldOwner[curFaces[faceI]] < nOldCells_ 
                            ? reverseCellMap_[oldOwner[curFaces[faceI]]] 
                            : addedCellRenumbering_[oldOwner[curFaces[faceI]]];
                label nei =   oldNeighbour[curFaces[faceI]] < nOldCells_ 
                            ? reverseCellMap_[oldNeighbour[curFaces[faceI]]]
                            : addedCellRenumbering_[oldNeighbour[curFaces[faceI]]];
                
                if (cellI == own)
                    neiCells[faceI] = nei;                    
                else if (cellI == nei)
                    neiCells[faceI] = own;                    
                else
                    FatalErrorIn("Foam::dynamicTopoFvMesh::reOrderFaces()") << nl
                            << " Could not determine neighbour faces."
                            << " Something's messed up." << nl
                            << abort(FatalError);

                nNeighbours++;
            }
            
            // Boundary faces are inserted normally. Update maps for now.
            // Face insertion for boundaries will be done after internal faces.
            if (visited[curFaces[faceI]] == -1) {
                label patchID = whichPatch(curFaces[faceI]);
                label bFaceIndex = boundaryPatchIndices[patchID]++;
                // Renumber the point-labels for this boundary-face
                face& faceRenumber = oldFaces[curFaces[faceI]];
                forAll(faceRenumber,pointI) {
                    if (faceRenumber[pointI] < nOldPoints_)
                        faceRenumber[pointI] = reversePointMap_[faceRenumber[pointI]];
                    else
                        faceRenumber[pointI] = addedPointRenumbering_[faceRenumber[pointI]];
                }
                // Renumber the edges in edgeToWatch
                if (edgeModification_ && twoDMotion_) {
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
                if (curFaces[faceI] < nOldFaces_){
                    faceMap_[bFaceIndex] = curFaces[faceI];
                    reverseFaceMap_[curFaces[faceI]] = bFaceIndex;
                } else {
                    addedFaceRenumbering.insert(curFaces[faceI],bFaceIndex);
                    addedReverseFaceRenumbering.insert(bFaceIndex,curFaces[faceI]);
                }                
            }
        }
        
        // Add internal faces in the increasing order of neighbours
        for (label neiSearch = 0; neiSearch < nNeighbours; neiSearch++) {
            
            // Find the lowest neighbour which is still valid
            label nextNei = -1;
            label minNei = nCells_;

            forAll (neiCells, ncI) {
                if (neiCells[ncI] > -1 && neiCells[ncI] < minNei) {
                    nextNei = ncI;
                    minNei = neiCells[ncI];
                }
            }

            if (nextNei > -1) {
                // Face is internal and gets reordered
                faceMap_[faceInOrder] = curFaces[nextNei];
                if (curFaces[nextNei] < nOldFaces_)
                    reverseFaceMap_[curFaces[nextNei]] = faceInOrder;
                
                // Renumber the point labels in this face
                face& faceRenumber = oldFaces[curFaces[nextNei]];
                forAll(faceRenumber, pointI) {
                    if (faceRenumber[pointI] < nOldPoints_)
                        faceRenumber[pointI] = reversePointMap_[faceRenumber[pointI]];
                    else
                        faceRenumber[pointI] = addedPointRenumbering_[faceRenumber[pointI]];
                }                
                
                // Renumber owner/neighbour
                label ownerRenumber =   oldOwner[curFaces[nextNei]] < nOldCells_ 
                                      ? reverseCellMap_[oldOwner[curFaces[nextNei]]] 
                                      : addedCellRenumbering_[oldOwner[curFaces[nextNei]]];
                label neighbourRenumber =   oldNeighbour[curFaces[nextNei]] < nOldCells_ 
                                          ? reverseCellMap_[oldNeighbour[curFaces[nextNei]]]
                                          : addedCellRenumbering_[oldNeighbour[curFaces[nextNei]]];
                     
                // Cell-reordering may require face-flipping..
                if (neighbourRenumber > ownerRenumber) {
                    faceRenumber.reverseFace();
                    if (fluxInterpolation_)
                        oldPhi[curFaces[nextNei]] *= -1.0;
                }
                
                // Insert entities into HashLists...
                faces_.append(faceRenumber);                
                owner_.append(ownerRenumber);
                neighbour_.append(neighbourRenumber);
                if (edgeModification_ && twoDMotion_) {
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
                if (fluxInterpolation_)
                    localPhi_.append(oldPhi[curFaces[nextNei]]);
                
                // Insert entities into mesh-reset lists
                faces[faceInOrder] = faceRenumber;                
                owner[faceInOrder] = ownerRenumber;
                neighbour[faceInOrder] = neighbourRenumber;

                // Stop the neighbour from being used again
                neiCells[nextNei] = -1;
                
                // Mark this face as visited
                visited[curFaces[nextNei]] = 0;

                faceInOrder++;
            } else {
                FatalErrorIn("dynamicTopoFvMesh::reOrderFaces()") << nl
                    << "Error in internal face insertion" << nl
                    << abort(FatalError);
            }
        }        
    }
    
    // All internal faces have been inserted. Now insert boundary faces.
    label oldIndex;
    for(label i=nInternalFaces_; i<nFaces_; i++) {
        if (faceMap_[i] == -1) {
            // This boundary face was added during the topology change
            oldIndex = addedFaceRenumbering[i];
        } else {
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
        if (fluxInterpolation_)
            localPhi_.append(oldPhi[oldIndex]);        
        
        // Insert entities into mesh-reset lists
        faces[faceInOrder] = oldFaces[oldIndex];                
        owner[faceInOrder] = ownerRenumber;
        neighbour[faceInOrder] = -1;        
        
        // Mark this face as visited
        visited[oldIndex] = 0;
        
        faceInOrder++;
    }
    
    // Renumber all cells
    for(HashList<cell>::iterator cIter = cells_.begin(); cIter != cells_.end(); cIter++) {
        
        cell& cellFaces = cIter();
        
        forAll(cellFaces,faceI) {
            if (cellFaces[faceI] < nOldFaces_) {
                cellFaces[faceI] = reverseFaceMap_[cellFaces[faceI]];
            } else {
                cellFaces[faceI] = addedFaceRenumbering[cellFaces[faceI]];
            }
        }
    }
    
    // Final check to ensure everything went okay
    if (debug) {
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
    // If cells were deleted during topology change, the numerical order ceases to be continuous.
    // Also, cells are always added at the end of the list by virtue of the HashList append 
    // method. Thus, cells would now have to be reordered so that bandwidth is reduced and 
    // renumbered to be sequential again.
    
    label currentCell, cellInOrder = 0, allCells = cells_.lastIndex() + 1; 
    SLList<label> nextCell; 
    labelList ncc(allCells, 0);
    labelList visited(allCells, 0);
    labelListList cellCellAddr(allCells);
    cellList oldCells(allCells);
    scalarField oldLengthScale(0), oldPressure(0);
    
    addedCellRenumbering_.clear();
    
    // Make a copy of the old cell-based HashLists, and clear them
    for(HashList<cell>::iterator cIter = cells_.begin(); cIter != cells_.end(); cIter++)
        oldCells[cIter.index()] = cIter();
    cells_.clear();
    if (edgeModification_) {
        oldLengthScale.setSize(allCells);
        for(HashList<scalar>::iterator cIter = lengthScale_.begin(); cIter != lengthScale_.end(); cIter++)
            oldLengthScale[cIter.index()] = cIter();
        lengthScale_.clear();
    }
    if (fluxInterpolation_) {
        oldPressure.setSize(allCells);
        for(HashList<scalar>::iterator cIter = localp_.begin(); cIter != localp_.end(); cIter++)
            oldLengthScale[cIter.index()] = cIter(); 
        localp_.clear();
    }    
    
    // Build a cell-cell addressing list
    HashList<label>::iterator ownIter = owner_.begin();
    HashList<label>::iterator neiIter = neighbour_.begin();    
    while(ownIter != owner_.end()) {
        if (neiIter() != -1) {
            ncc[ownIter()]++;
            ncc[neiIter()]++;     
        }
        ownIter++; neiIter++;
    }     
    forAll (cellCellAddr, cellI) {
        cellCellAddr[cellI].setSize(ncc[cellI]); 
        // Mark off deleted cells as "visited"
        if (ncc[cellI] == 0) visited[cellI] = 1;
    }
    ncc = 0;
    ownIter = owner_.begin(); neiIter = neighbour_.begin();
    while(ownIter != owner_.end()) {
        if (neiIter() != -1) {
            cellCellAddr[ownIter()][ncc[ownIter()]++] = neiIter();
            cellCellAddr[neiIter()][ncc[neiIter()]++] = ownIter();
        }
        ownIter++; neiIter++;
    }
    
    // Let's get to the "business bit" of the band-compression
    forAll (visited, cellI) {
        
        // Find the first cell that has not been visited yet
        if (visited[cellI] == 0) {
            
            // Use this cell as a start            
            currentCell = cellI;
            nextCell.append(currentCell);

            // Loop through the nextCell list. Add the first cell into the
            // cell order if it has not already been visited and ask for its
            // neighbours. If the neighbour in question has not been visited,
            // add it to the end of the nextCell list
            while (nextCell.size() > 0) {
                
                currentCell = nextCell.removeHead();
                
                if (visited[currentCell] == 0) {
                    
                    // Mark as visited and update cell mapping info
                    visited[currentCell] = 1;
                    cellMap_[cellInOrder] = currentCell;
                    if (currentCell < nOldCells_)
                        reverseCellMap_[currentCell] = cellInOrder;
                    else
                        addedCellRenumbering_.insert(currentCell,cellInOrder);
                    
                    // Insert entities into HashLists...
                    cells_.append(oldCells[currentCell]);
                    if (edgeModification_) 
                        lengthScale_.append(oldLengthScale[currentCell]);
                    if (fluxInterpolation_) 
                        localp_.append(oldPressure[currentCell]);
                    
                    cellInOrder++;

                    // Find if the neighbours have been visited
                    const labelList& neighbours = cellCellAddr[currentCell];

                    forAll (neighbours, nI) {
                        if (visited[neighbours[nI]] == 0) {
                            // Not visited, add to the list
                            nextCell.append(neighbours[nI]);
                        }
                    }
                }
            }
        }
    }        
    
    if(debug) {
        if (sum(visited) != allCells)
            FatalErrorIn("Foam::dynamicTopoFvMesh::reOrderCells()") << nl
                    << " Algorithm did not visit every cell in the mesh."
                    << " Something's messed up." << nl
                    << abort(FatalError);         
    }
}

// Reorder the faces in upper-triangular order, and generate mapping information
void Foam::dynamicTopoFvMesh::reOrderMesh(
    pointField& points,
    pointField& pointsZeroVolume,
    faceList& faces, 
    labelList& owner, 
    labelList& neighbour
)
{    
    // Allocate for the mapping information
    pointMap_.setSize(nPoints_, -1);
    faceMap_.setSize(nFaces_, -1);
    cellMap_.setSize(nCells_, -1);
    
    if (debug) {
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
    reOrderPoints(points, pointsZeroVolume);
    
    // Reorder the cells
    reOrderCells();
    
    // Reorder the faces
    reOrderFaces(faces, owner, neighbour);
}

// Calculate the edge length-scale for the mesh
void Foam::dynamicTopoFvMesh::calculateLengthScale()
{
    if (edgeModification_) {
        
        // Initialize the length-density field for the current time-step
        volScalarField lengthDensity
        (
            IOobject
            (
                "lengthDensity",
                this->time().timeName(),
                (*this),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (*this),
            dimensionedScalar("length", dimensionSet(0,1,0,0,0), 0),
            fixedValueFvPatchScalarField::typeName
        );
        
        // Set the boundary conditions for the Laplace's equation for length-density 
        forAll(lengthDensity.boundaryField(), patchI) {  
            fvPatchField<scalar> &bPatch = lengthDensity.boundaryField()[patchI];
            if ((bPatch.type() == "wedge")) {
                bPatch = 0.0;
            } else {
                label start = patchStarts_[patchI];
                forAll(bPatch, faceI) {
                    edge& etw = edgeToWatch_[start+faceI];
                    bPatch[faceI] = 1.0/mag(meshPoints_[etw[0]] - meshPoints_[etw[1]]);
                }
            }
        }
        
        // Solve for the length-density function
        Foam::solve
        (
            fvm::laplacian
            (
                growthFactor_,
                lengthDensity,
                "laplacian(growthFactor,lengthDensity)"
            )
        );

        // Copy the most recent length-scale values
        for(HashList<scalar>::iterator l = lengthScale_.begin(); l != lengthScale_.end(); l++) {
            l() = 1.0/lengthDensity.internalField()[l.index()];
        }
    }
}

// 2D Edge-swapping engine (for wedge meshes one-cell thick, w/o collapseEdges)
void Foam::dynamicTopoFvMesh::swap2DEdges()
{
    bool found, foundinner;       
    label otherPointIndex[4], nextToOtherPoint[4], commonFaceIndex[4], commonIntFaceIndex[4];
    label c0BdyIndex[2], c0IntIndex[2], c1BdyIndex[2], c1IntIndex[2];
    face  c0BdyFace[2],  c0IntFace[2],  c1BdyFace[2],  c1IntFace[2];
    face f, commonFaces[4], commonIntFaces[4];       
    edge commonEdges[2], firstEdge(0,0);
    //vector xC0o = vector::zero, xC0n = vector::zero;
    //vector xC1o = vector::zero, xC1n = vector::zero;    

    for(HashList<face>::iterator fIter = faces_.begin(); fIter != faces_.end(); fIter++) {
        
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
        // Note: This will NOT work for pure-wedge meshes
        findPrismFaces(findex,cell_0,c0BdyFace,c0BdyIndex,c0IntFace,c0IntIndex);
        findPrismFaces(findex,cell_1,c1BdyFace,c1BdyIndex,c1IntFace,c1IntIndex);

        /*
        if (debug) {
            Info << "Cell: " << c0 << endl;                
            Info << "Boundary faces: " << c0BdyIndex[0] << ": " << c0BdyFace[0] << endl;
            Info << "Boundary faces: " << c0BdyIndex[1] << ": " << c0BdyFace[1] << endl; 
            Info << "Interior faces: " << c0IntIndex[0] << ": " << c0IntFace[0] << endl;
            Info << "Interior faces: " << c0IntIndex[1] << ": " << c0IntFace[1] << endl;                
            Info << "Cell: " << c1 << endl;
            Info << "Boundary faces: " << c1BdyIndex[0] << ": " << c1BdyFace[0] << endl;
            Info << "Boundary faces: " << c1BdyIndex[1] << ": " << c1BdyFace[1] << endl; 
            Info << "Interior faces: " << c1IntIndex[0] << ": " << c1IntFace[0] << endl;
            Info << "Interior faces: " << c1IntIndex[1] << ": " << c1IntFace[1] << endl;                
        }
        */

        // Find the common faces / edges on the boundary
        // At the end of this loop, commonFaces [0] & [1] share commonEdge [0]
        // and commonFaces [2] & [3] share commonEdge [1]
        // Also, commonFaces[0] & [2] are connected to cell[0],
        // and commonFaces[1] & [3] are connected to cell[1]
        found = false; foundinner = false;
        edgeList e1 = c0BdyFace[0].edges(), e2 = c1BdyFace[0].edges();
        edgeList e3 = c0BdyFace[1].edges(), e4 = c1BdyFace[1].edges();
        forAll(e1,edgeI) {
            forAll(e2,edgeJ) {
                if (e1[edgeI] == e2[edgeJ]) {
                    // These two faces share an edge, store for posterity
                    commonFaces[0] = c0BdyFace[0]; commonFaces[1] = c1BdyFace[0];
                    commonFaces[2] = c0BdyFace[1]; commonFaces[3] = c1BdyFace[1];  
                    commonFaceIndex[0] = c0BdyIndex[0]; commonFaceIndex[1] = c1BdyIndex[0];
                    commonFaceIndex[2] = c0BdyIndex[1]; commonFaceIndex[3] = c1BdyIndex[1];
                    commonEdges[0] = e1[edgeI];
                    forAll(e3,edgeK) {
                        forAll(e4,edgeL) {
                            if (e3[edgeK] == e4[edgeL]) {
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
        if (!found) {
            // A match was obviously not found before, but we know the common faces now.
            commonFaces[0] = c0BdyFace[0]; commonFaces[1] = c1BdyFace[1];
            commonFaces[2] = c0BdyFace[1]; commonFaces[3] = c1BdyFace[0]; 
            commonFaceIndex[0] = c0BdyIndex[0]; commonFaceIndex[1] = c1BdyIndex[1];
            commonFaceIndex[2] = c0BdyIndex[1]; commonFaceIndex[3] = c1BdyIndex[0];                
            // Start a new search for common edges
            forAll(e1,edgeI) {
                forAll(e4,edgeJ) {  
                    if (e1[edgeI] == e4[edgeJ]) {
                        commonEdges[0] = e1[edgeI]; 
                        forAll(e2,edgeK) {
                            forAll(e3,edgeL) {
                                if (e2[edgeK] == e3[edgeL]) {
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

        // Construct a circle passing through the three points of the first face, define its radius
        label zeroIndex  = commonFaces[0][0];
        label oneIndex   = commonFaces[0][1];
        label twoIndex   = commonFaces[0][2];
        point& pointZero = meshPoints_[zeroIndex];
        point& pointOne  = meshPoints_[oneIndex];
        point& pointTwo  = meshPoints_[twoIndex];
        point  center    = circumCenter(pointZero, pointOne, pointTwo, zeroIndex, oneIndex, twoIndex);
        scalar radius    = (pointZero - center)&(pointZero - center);

        // Find the isolated point on the other face, and the point next to it (for orientation)
        findIsolatedPoint(commonFaces[1], commonEdges[0], otherPointIndex[1], nextToOtherPoint[1]);

        // ...and determine whether it lies in this circle
        point otherPoint = meshPoints_[otherPointIndex[1]];
        if ( ((otherPoint - center)&(otherPoint - center)) < radius ) {

            // This face needs to be flipped... 
            if (debug) {
                Info << nl << nl << "Face: " << findex << " needs to be flipped. " << endl;
                Info << "Cell[0]: " << c0 << ": " << cell_0 << endl;                
                Info << "Cell[1]: " << c1 << ": " << cell_1 << endl;                 
                Info << "Common Faces: Set 1: " << commonFaceIndex[0] << ": " << commonFaces[0] << ", " 
                                                << commonFaceIndex[1] << ": " << commonFaces[1] << endl;
                Info << "Common Faces: Set 2: " << commonFaceIndex[2] << ": " << commonFaces[2] << ", " 
                                                << commonFaceIndex[3] << ": " << commonFaces[3] << endl; 
                Info << "Old face: " << faces_[findex] << endl;                
            }

            // Set the flag
            topoChangeFlag_ = true;

            // Find the other three points that don't lie on shared edges
            // and the points next to them (for orientation)
            findIsolatedPoint(commonFaces[0], commonEdges[0], otherPointIndex[0], nextToOtherPoint[0]);
            findIsolatedPoint(commonFaces[2], commonEdges[1], otherPointIndex[2], nextToOtherPoint[2]);  
            findIsolatedPoint(commonFaces[3], commonEdges[1], otherPointIndex[3], nextToOtherPoint[3]);

            // Find the other two edges on the face being flipped (besides the common edges)
            // First edge detected belongs to cell[0] by default
            edgeList eThis = thisFace.edges();
            forAll(eThis,edgeI) {
                if ( eThis[edgeI] != commonEdges[0] && eThis[edgeI] != commonEdges[1] ) {
                    if (eThis[edgeI][0] == nextToOtherPoint[0] || eThis[edgeI][1] == nextToOtherPoint[0])
                        firstEdge = eThis[edgeI];
                }
            }       

            // Find the interior faces that share the first edge
            // At the end of this loop, commonIntFaces [0] & [1] share firstEdge 
            // and commonIntFaces [2] & [3] share the secondEdge,
            // where [0],[2] lie on cell[0] and [1],[3] lie on cell[1]
            found = false;
            edgeList e1 = c0IntFace[0].edges();
            forAll(e1,edgeI) {
                if ( e1[edgeI] == firstEdge ) { 
                    commonIntFaces[0] = c0IntFace[0]; commonIntFaces[2] = c0IntFace[1];
                    commonIntFaceIndex[0] = c0IntIndex[0]; commonIntFaceIndex[2] = c0IntIndex[1];
                    found = true; break; 
                }
            }
            if (!found) {
                // The edge was obviously not found before
                commonIntFaces[0] = c0IntFace[1]; commonIntFaces[2] = c0IntFace[0];
                commonIntFaceIndex[0] = c0IntIndex[1]; commonIntFaceIndex[2] = c0IntIndex[0];
            }

            found = false;
            edgeList e3 = c1IntFace[0].edges(); 
            forAll(e3,edgeI) {
                if ( e3[edgeI] == firstEdge ) { 
                    commonIntFaces[1] = c1IntFace[0]; commonIntFaces[3] = c1IntFace[1];
                    commonIntFaceIndex[1] = c1IntIndex[0]; commonIntFaceIndex[3] = c1IntIndex[1];
                    found = true; break; 
                }
            }
            if (!found) {
                // The edge was obviously not found before
                commonIntFaces[1] = c1IntFace[1]; commonIntFaces[3] = c1IntFace[0];
                commonIntFaceIndex[1] = c1IntIndex[1]; commonIntFaceIndex[3] = c1IntIndex[0];
            }       

            /*
            // Obtain the cell centers for both cells before their faces are modified
            if (fluxInterpolation_) {
                xC0o = cellCenter(cell_0);
                xC1o = cellCenter(cell_1);
            }            
            */
            
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
            if (edgeModification_) {
                edge& edgeToModify = edgeToWatch_[findex];
                edgeToModify[0] = otherPointIndex[0];
                edgeToModify[1] = otherPointIndex[1];   
            }            
            
            // Four modified boundary faces need to be constructed, but right-handedness is also important
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

            if (debug) {
                Info << "New flipped face: " << newFace << endl;
                Info << "New boundary face[0]" << commonFaceIndex[0] << ": " << newBdyFace0 << endl;
                Info << "New boundary face[1]" << commonFaceIndex[1] << ": " << newBdyFace1 << endl;
                Info << "New boundary face[2]" << commonFaceIndex[2] << ": " << newBdyFace2 << endl;
                Info << "New boundary face[3]" << commonFaceIndex[3] << ": " << newBdyFace3 << endl;

                if (fluxInterpolation_) {
                    // Output fluxes
                    Info << "This face: " << findex << " Phi:" << localPhi_[findex] 
                         << " Owner: " << owner_[findex] << " Neighbour:" << neighbour_[findex] << endl;
                    Info << "commonIntFaces[0]:" << commonIntFaceIndex[0] << " Phi:" << localPhi_[commonIntFaceIndex[0]]
                         << " Owner: " << owner_[commonIntFaceIndex[0]] << " Neighbour:" << neighbour_[commonIntFaceIndex[0]] << endl; 
                    Info << "commonIntFaces[1]:" << commonIntFaceIndex[1] << " Phi:" << localPhi_[commonIntFaceIndex[1]]
                         << " Owner: " << owner_[commonIntFaceIndex[1]] << " Neighbour:" << neighbour_[commonIntFaceIndex[1]] << endl;
                    Info << "commonIntFaces[2]:" << commonIntFaceIndex[2] << " Phi:" << localPhi_[commonIntFaceIndex[2]]
                         << " Owner: " << owner_[commonIntFaceIndex[2]] << " Neighbour:" << neighbour_[commonIntFaceIndex[2]] << endl;
                    Info << "commonIntFaces[3]:" << commonIntFaceIndex[3] << " Phi:" << localPhi_[commonIntFaceIndex[3]]
                         << " Owner: " << owner_[commonIntFaceIndex[3]] << " Neighbour:" << neighbour_[commonIntFaceIndex[3]] << endl;
                }
            }
            
            // Check the orientation of the two quad faces, and modify as necessary
            label newOwn=0, newNei=0;
            bool flipOption;

            // The quad face belonging to cell[1] now becomes a part of cell[0]
            if ( neighbour_[commonIntFaceIndex[1]] == -1 ) {
                // Boundary face
                // Face doesn't need to be flipped, just update the owner
                f          = commonIntFaces[1];
                newOwn     = c0;
                newNei     = -1;
                flipOption = false;
            } else if ( owner_[commonIntFaceIndex[1]] == c1 ) {
                // This face is on the interior, check for previous owner 
                // Upper-triangular ordering has to be maintained, however...
                if ( c0 > neighbour_[commonIntFaceIndex[1]] ) {
                    // Flip is necessary
                    f          = commonIntFaces[1].reverseFace();
                    newOwn     = neighbour_[commonIntFaceIndex[1]];
                    newNei     = c0;
                    flipOption = true;
                } else {
                    // Flip isn't necessary, just change the owner
                    f          = commonIntFaces[1];
                    newOwn     = c0;
                    newNei     = neighbour_[commonIntFaceIndex[1]];
                    flipOption = false;
                }
            } else if ( neighbour_[commonIntFaceIndex[1]] == c1 ) {
                // This face is on the interior, check for previous neighbour
                // Upper-triangular ordering has to be maintained, however...
                if ( c0 < owner_[commonIntFaceIndex[1]] ) {
                    // Flip is necessary
                    f          = commonIntFaces[1].reverseFace();
                    newOwn     = c0;
                    newNei     = owner_[commonIntFaceIndex[1]];
                    flipOption = true;
                } else {
                    // Flip isn't necessary, just change the neighbour
                    f          = commonIntFaces[1];
                    newOwn     = owner_[commonIntFaceIndex[1]];
                    newNei     = c0;
                    flipOption = false;
                }                    
            }

            faces_[commonIntFaceIndex[1]] = f;
            cell_0[c0count++] = commonIntFaceIndex[0];
            cell_0[c0count++] = commonIntFaceIndex[1];
            owner_[commonIntFaceIndex[1]] = newOwn;
            neighbour_[commonIntFaceIndex[1]] = newNei;
            // Modify the local-flux field
            if (flipOption && fluxInterpolation_) {
                localPhi_[commonIntFaceIndex[1]] *= -1.0;
            }            

            // The quad face belonging to cell[0] now becomes a part of cell[1]
            if ( neighbour_[commonIntFaceIndex[2]] == -1 ) {
                // Boundary face
                // Face doesn't need to be flipped, just update the owner
                f          = commonIntFaces[2];
                newOwn     = c1;
                newNei     = -1;
                flipOption = false;
            } else if ( owner_[commonIntFaceIndex[2]] == c0 ) {
                // This face is on the interior, check for previous owner 
                // Upper-triangular ordering has to be maintained, however...
                if ( c1 > neighbour_[commonIntFaceIndex[2]] ) {
                    // Flip is necessary
                    f          = commonIntFaces[2].reverseFace();
                    newOwn     = neighbour_[commonIntFaceIndex[2]];
                    newNei     = c1;
                    flipOption = true;
                } else {
                    // Flip isn't necessary, just change the owner
                    f          = commonIntFaces[2];
                    newOwn     = c1;
                    newNei     = neighbour_[commonIntFaceIndex[2]];
                    flipOption = false;
                }
            } else if ( neighbour_[commonIntFaceIndex[2]] == c0 ) {
                // This face is on the interior, check for previous neighbour
                // Upper-triangular ordering has to be maintained, however...
                if ( c1 < owner_[commonIntFaceIndex[2]] ) {
                    // Flip is necessary
                    f          = commonIntFaces[2].reverseFace();
                    newOwn     = c1;
                    newNei     = owner_[commonIntFaceIndex[2]];
                    flipOption = true;
                } else {
                    // Flip isn't necessary, just change the neighbour
                    f          = commonIntFaces[2];
                    newOwn     = owner_[commonIntFaceIndex[2]];
                    newNei     = c1;
                    flipOption = false;
                }                    
            }

            faces_[commonIntFaceIndex[2]] = f;
            cell_1[c1count++] = commonIntFaceIndex[2];
            cell_1[c1count++] = commonIntFaceIndex[3];
            owner_[commonIntFaceIndex[2]] = newOwn;
            neighbour_[commonIntFaceIndex[2]] = newNei;
            // Modify the local-flux field
            if (flipOption && fluxInterpolation_) {
                localPhi_[commonIntFaceIndex[2]] *= -1.0;
            } 
            
            /*
            // Obtain the cell centers for both cells after modification and interpolate pressure
            if (fluxInterpolation_) {
                xC0n = cellCenter(cell_0);
                xC1n = cellCenter(cell_1);
                localp_[c0] += ((xC0n - xC0o)&localGradp_[c0]);
                localp_[c1] += ((xC1n - xC1o)&localGradp_[c1]);              
            }
            */            
            
            // Calculate flux for the flipped face
            if (fluxInterpolation_) {
                scalar sign1, sign2;                
                sign1 = (owner_[commonIntFaceIndex[0]] == c0) ? 1.0 : -1.0;
                sign2 = (owner_[commonIntFaceIndex[1]] == c0) ? 1.0 : -1.0;
                localPhi_[findex] = - (sign1*localPhi_[commonIntFaceIndex[0]]) 
                                    - (sign2*localPhi_[commonIntFaceIndex[1]]);               
            }
        }
    }   
}

// Initialize the edge-length field
void Foam::dynamicTopoFvMesh::initEdgeLengths()
{
    lengthScale_.setSize(nCells_, 0.0);
    
    if (twoDMotion_) { 
        
        // Allocate fields
        edge nullEdge(0,0);  
        edgeToWatch_.setSize(nFaces_,nullEdge);
        
        // Loop through all quad-faces and build initial edge-lengths
        bool found;
        for(label findex=0; findex<nFaces_; findex++) {
            face& fi = faces_[findex];
            if (fi.size() == 4) {
                label c0 = owner_[findex];
                cell& cell_0 = cells_[c0];
                edgeList efi = fi.edges();
                
                // Look for a triangular face on the boundary
                found=false;
                forAll(cell_0, i) {
                    label faceIndex = cell_0[i];
                    face& fj=faces_[faceIndex];
                    if (neighbour_[faceIndex] == -1 && fj.size() == 3) {
                        // Match an edge, and compute its length
                        edgeList efj = fj.edges();
                        forAll(efi, indexI) {
                            forAll(efj, indexJ) {
                                if (efi[indexI] == efj[indexJ]) {
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

// Initialize internal edges (for 3D tet-meshes)
void Foam::dynamicTopoFvMesh::initInternalEdges()
{
    labelList tmpInternalEdges(this->nEdges(), 1);
    labelList tmpStartFaceIndices(this->nEdges(), -1);
    const edgeList& meshEdges           = this->edges();
    const labelListList& faceToEdgeList = this->faceEdges();
    const polyBoundaryMesh& boundary    = this->boundaryMesh();
    
    // Blank-out boundary edges...
    for(label i=0; i<numPatches_; i++) {
        for(label j=boundary[i].start(); j<boundary[i].start()+boundary[i].size(); j++) {
            const labelList& faceToEdge = faceToEdgeList[j];
            forAll(faceToEdge, edgeI) 
                tmpInternalEdges[faceToEdge[edgeI]] = 0; 
        }
    }

    // Assign the number of internal-edges
    nInternalEdges_ = sum(tmpInternalEdges);
    
    // Allocate the internal edge list, and provide edge-labels 
    // for later use
    edge nullEdge(0,0);
    label edgeCounter = 0;
    internalEdges_.setSize(nInternalEdges_, nullEdge);
    forAll(tmpInternalEdges, edgeI) {
        if (tmpInternalEdges[edgeI]) {
            internalEdges_[edgeCounter] = meshEdges[edgeI];
            tmpInternalEdges[edgeI] = edgeCounter++;
        } else {
            tmpInternalEdges[edgeI] = -1;
        }
    }
    
    // Loop through all internal faces and 
    // add start-face indices for internal edges
    startFaceIndex_.setSize(nInternalEdges_, -1);
    for(label i=0; i<nInternalFaces_; i++) {
        const labelList& faceToEdge = faceToEdgeList[i];
        forAll(faceToEdge, edgeI) {
            label& indexToCheck = tmpInternalEdges[faceToEdge[edgeI]];
            if (indexToCheck >= 0) {
                startFaceIndex_[indexToCheck] = i;
                indexToCheck = -1;
            }
        }
    }
    
    // Perform a final-check to ensure that allocation is complete
#   ifdef FULLDEBUG
    if (sum(tmpInternalEdges) != -this->nEdges())
        FatalErrorIn("dynamicTopoFvMesh::initInternalEdges()") << nl
                << " Internal edge-allocation failed. " << nl
                << abort(FatalError);
#   endif                
}

// 2D Edge-bisection/collapse engine
void Foam::dynamicTopoFvMesh::edgeBisectCollapse2D()
{
    // Loop through all quad-faces and bisect/collapse 
    // edges (quad-faces) by the criterion:
    // Bisect when boundary edge-length > ratioMax_*originalLength
    // Collapse when boundary edge-length < ratioMin_*originalLength
    
    HashList<face>::iterator fIter = faces_.begin();
    while(fIter != faces_.end()) {
        
        // Retrieve the index for this iterator
        label findex = fIter.index();
        
        // Reference to this face...
        face& thisFace = fIter();
        
        // Select only quad-faces
        if (thisFace.size() == 4) {
            
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
                scale = lengthScale_[c0];
            else
                scale = 0.5*(lengthScale_[c0]+lengthScale_[c1]);

            /*
            //-- For testing
            if (this->time().value() > 1.0 && this->time().value() < 1.02 && findex == 62) {
                // Set the flag
                topoChangeFlag_ = true;
                
                // Bisect this face
                bisectQuadFace(findex, thisFace);
                
                // Move on to the next face
                fIter++;                
            }
            //-- For testing
            */
            
            //== Edge Bisection ==//
            if(length > ratioMax_*scale) {                                  
                
                // Consider only triangle-prisms
                cell &cell_0 = cells_[c0];                
                if (cell_0.nFaces() > 5) continue;  
                if (c1 != -1) {
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
            if(length < ratioMin_*scale) {
                
                // Consider only triangle-prisms
                cell &cell_0 = cells_[c0];
                if (cell_0.nFaces() > 5) continue;
                if (c1 != -1) {
                    cell &cell_1 = cells_[c1];
                    if (cell_1.nFaces() > 5) continue;
                } 
                
                // Collapse this face
                bool success = collapseQuadFace(findex, thisFace);
                
                // Increment the iterator to move on to the next face...
                fIter++;
                
                // The face can safely be deleted, since the iterator points
                // to the next valid face on the face-list.
                if (success) {
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
void Foam::dynamicTopoFvMesh::bisectQuadFace(const label findex, face& thisFace)
{
    // Local variables
    bool found, flipOption;
    label otherPointIndex[4], nextToOtherPoint[4], replaceFace, n0=-1, n1=-1;
    label c0BdyIndex[2], c0IntIndex[2], c1BdyIndex[2], c1IntIndex[2];
    face  c0BdyFace[2],  c0IntFace[2],  c1BdyFace[2],  c1IntFace[2];
    edge  tmpEdge(0,0), commonEdges[2], firstEdge(0,0), secondEdge(0,0);
    //vector xC0o = vector::zero, xC0n = vector::zero, xnewC0 = vector::zero;
    //vector xC1o = vector::zero, xC1n = vector::zero, xnewC1 = vector::zero;    
    
    // Get the two cells on either side...
    label c0 = owner_[findex], c1 = neighbour_[findex];
    
    // Find the prism faces for cell[0].
    cell &cell_0 = cells_[c0];    
    findPrismFaces(findex,cell_0,c0BdyFace,c0BdyIndex,c0IntFace,c0IntIndex);

    if (debug) {
        Info << nl << nl << "Face: " << findex << ": " << thisFace << " is to be bisected. " << endl;
        Info << "Cell[0]: " << c0 << ": " << cell_0 << endl;
        forAll(cell_0, faceI)
            Info << cell_0[faceI] << ": " << faces_[cell_0[faceI]] << endl;
    }

    // Find the common-edge between the triangular boundary faces
    // and the face under consideration.
    findCommonEdge(c0BdyFace[0],thisFace,commonEdges[0]);
    findCommonEdge(c0BdyFace[1],thisFace,commonEdges[1]);

    // Find the isolated point on both boundary faces of cell[0]
    findIsolatedPoint(c0BdyFace[0], commonEdges[0], otherPointIndex[0], nextToOtherPoint[0]);
    findIsolatedPoint(c0BdyFace[1], commonEdges[1], otherPointIndex[1], nextToOtherPoint[1]);

    // Obtain the cell center before the faces are modified,
    // and add a new cell for the local pressure field
    if (fluxInterpolation_) {
        //xC0o = cellCenter(cell_0);
        localp_.append(localp_[c0]);
        //localGradp_.append(localGradp_[c0]);
    }    
    
    // Add two new points to the end of the list
    label newPtIndex0 = meshPoints_.append(0.5*(meshPoints_[commonEdges[0][0]] + meshPoints_[commonEdges[0][1]]));
    label newPtIndex1 = meshPoints_.append(0.5*(meshPoints_[commonEdges[1][0]] + meshPoints_[commonEdges[1][1]]));
    nPoints_ += 2; 

    // Add a new prism cell to the end of the list
    cell tmpPrismCell(5);
    label newCellIndex0 = cells_.append(tmpPrismCell);
    cell &newCell0 = cells_[newCellIndex0];

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
    forAll(eThis, edgeI) {
        if (eThis[edgeI] != commonEdges[0] && eThis[edgeI] != commonEdges[1]) {
            if (eThis[edgeI][0] == nextToOtherPoint[0] || eThis[edgeI][1] == nextToOtherPoint[0])
                firstEdge = eThis[edgeI];
            if (eThis[edgeI][0] == nextToOtherPoint[1] || eThis[edgeI][1] == nextToOtherPoint[1])
                secondEdge = eThis[edgeI];
        }
    }

    /*
    // Obtain face-area from the cross-product
    scalar oldArea = mag( (meshPoints_[thisFace[1]]-meshPoints_[thisFace[0]])
                         ^(meshPoints_[thisFace[2]]-meshPoints_[thisFace[1]]) );
    */
    
    // Modify point-labels on the quad face under consideration
    replacePointLabel(commonEdges[0].otherVertex(nextToOtherPoint[0]), newPtIndex0, thisFace);
    replacePointLabel(nextToOtherPoint[1], newPtIndex1, thisFace);
    
    /*
    // Recalculate for new area after bisection
    scalar newArea = mag( (meshPoints_[thisFace[1]]-meshPoints_[thisFace[0]])
                         ^(meshPoints_[thisFace[2]]-meshPoints_[thisFace[1]]) );  
    
    // Modify the flux for this face...
    scalar newBisectFlux=0.0;
    if (fluxInterpolation_) {
        newBisectFlux = (1.0-(newArea/oldArea))*localPhi_[findex];
        localPhi_[findex] *= (newArea/oldArea);
    }
    */

    // Change the edge-length criteria for this face
    tmpEdge[0] = nextToOtherPoint[0];
    tmpEdge[1] = newPtIndex0;
    edgeToWatch_[findex] = tmpEdge;
    
    if (debug) Info << "Modified thisFace: " << findex << ": " << thisFace << endl;

    // Find the interior face that contains secondEdge
    found = false;
    edgeList e1 = c0IntFace[0].edges();
    forAll(e1, edgeI) {
        if ( e1[edgeI] == secondEdge ) {
            replaceFaceLabel(c0IntIndex[0],-1,cell_0);
            replaceFace = c0IntIndex[0];
            found = true; break;
        }
    }
    if (!found) {
        // The edge was obviously not found before
        replaceFaceLabel(c0IntIndex[1],-1,cell_0);
        replaceFace = c0IntIndex[1];
    }
    
    // Add two new points for the zero-volume point-list
    pointsZeroVol_.append(meshPoints_[secondEdge.commonVertex(commonEdges[0])]);
    pointsZeroVol_.append(meshPoints_[secondEdge.commonVertex(commonEdges[1])]);    

    // Check if face reversal is necessary for the replacement
    scalar sign = 0.0;
    if (owner_[replaceFace] == c0) {
        if (neighbour_[replaceFace] == -1) {
            // Change the owner (no flux-flipping)
            owner_[replaceFace] = newCellIndex0;
            sign = 1.0; flipOption = false;
        } else {
            n0 = neighbour_[replaceFace];
            // This face has to be reversed (and flux is to be flipped)
            faces_[replaceFace] = faces_[replaceFace].reverseFace();
            owner_[replaceFace] = neighbour_[replaceFace];
            neighbour_[replaceFace] = newCellIndex0;
            sign = -1.0; flipOption = true;
        }
    } else {
        n0 = owner_[replaceFace];
        // Keep owner, but change neighbour (no flux-flipping)
        neighbour_[replaceFace] = newCellIndex0;
        sign = -1.0; flipOption = false;
    }
    // Modify the local-flux field
    if (flipOption && fluxInterpolation_) {
        localPhi_[replaceFace] *= -1.0;
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
    newFaceIndex = insertFace(-1, tmpQuadFace, c0, newCellIndex0, edgeToWatch);
    replaceFaceLabel(-1, newFaceIndex, newCell0);
    replaceFaceLabel(-1, newFaceIndex, cell_0);

    // Calculate fluxes for this face to satisfy zero-divergence. 
    // Assumes that no fluxes are present on boundary triangle faces    
    if (fluxInterpolation_) {
        //localPhi_[newFaceIndex] = newBisectFlux + (sign*localPhi_[replaceFace]);
        localPhi_[newFaceIndex] = -1.0*localPhi_[replaceFace];
    }    

    // Second boundary face; Owner = newCell[0] & Neighbour = [-1]
    tmpTriFace[0] = otherPointIndex[0];
    tmpTriFace[1] = newPtIndex0;
    tmpTriFace[2] = commonEdges[0].otherVertex(nextToOtherPoint[0]);
    edgeToWatch[0] = edgeToWatch[1] = 0;
    newFaceIndex = insertFace(whichPatch(c0BdyIndex[0]), tmpTriFace, newCellIndex0, -1, edgeToWatch);
    replaceFaceLabel(-1, newFaceIndex, newCell0);

    // Third boundary face; Owner = c[0] & Neighbour = [-1] 
    tmpTriFace[0] = otherPointIndex[1];
    tmpTriFace[1] = newPtIndex1;
    tmpTriFace[2] = commonEdges[1].otherVertex(nextToOtherPoint[1]);
    edgeToWatch[0] = edgeToWatch[1] = 0;
    newFaceIndex = insertFace(whichPatch(c0BdyIndex[1]), tmpTriFace, c0, -1, edgeToWatch);
    replaceFaceLabel(-1, newFaceIndex, cell_0);

    if (c1 == -1) {

        // The quad boundary face resulting from bisection; Owner = newCell[0] & Neighbour = [-1]
        tmpQuadFace[0] = newPtIndex1;
        tmpQuadFace[1] = nextToOtherPoint[1];
        tmpQuadFace[2] = commonEdges[0].otherVertex(nextToOtherPoint[0]);
        tmpQuadFace[3] = newPtIndex0;
        edgeToWatch[0] = newPtIndex1;
        edgeToWatch[1] = nextToOtherPoint[1];
        newFaceIndex = insertFace(whichPatch(findex), tmpQuadFace, newCellIndex0, -1, edgeToWatch);
        replaceFaceLabel(-1, newFaceIndex, newCell0);
     
        /*
        if (fluxInterpolation_) {
            localPhi_[newFaceIndex] = newBisectFlux;
            // Obtain the new cell centers and interpolate for pressure
            xC0n   = cellCenter(cell_0); 
            xnewC0 = cellCenter(newCell0);
            localp_[newCellIndex0] = localp_[c0] + ((xnewC0 - xC0o)&localGradp_[c0]);
            localp_[c0] += ((xC0n - xC0o)&localGradp_[c0]);
        }
        */

        if (debug) {
            Info << "Modified Cell[0]: " << c0 << ": " << cell_0 << endl;
            forAll(cell_0, faceI)
                Info << cell_0[faceI] << ": " << faces_[cell_0[faceI]] << endl;
            Info << "New Cell[0]: " << newCellIndex0 << ": " << newCell0 << endl;
            forAll(newCell0, faceI)
                Info << newCell0[faceI] << ": " << faces_[newCell0[faceI]] << endl;
        }

    } else {

        cell &cell_1 = cells_[c1];

        // Find the prism faces for cell[1].
        findPrismFaces(findex, cell_1, c1BdyFace, c1BdyIndex, c1IntFace, c1IntIndex);

        // Add a new prism cell to the end of the list
        label newCellIndex1 = cells_.append(tmpPrismCell);
        cell &newCell1 = cells_[newCellIndex1];
        
        // Add a new element to the lengthScale field
        // (Currently the same as cell[1])
        lengthScale_.append(lengthScale_[c1]);

        if (debug) {
            Info << "Cell[1]: " << c1 << ": " << cell_1 << endl;
            forAll(cell_1, faceI)
                Info << cell_1[faceI] << ": " << faces_[cell_1[faceI]] << endl;
        }   
        
        // Obtain the cell center before the faces are modified,
        // and add a new cell for the local pressure field
        if (fluxInterpolation_) { 
            //xC1o = cellCenter(cell_1);
            localp_.append(localp_[c1]);
            //localGradp_.append(localGradp_[c1]);
        }        
        
        // Find the interior face that contains secondEdge
        found = false;
        edgeList e2 = c1IntFace[0].edges();
        forAll(e2, edgeI) {
            if ( e2[edgeI] == secondEdge ) {
                replaceFaceLabel(c1IntIndex[0], -1, cell_1);
                replaceFace = c1IntIndex[0];
                found = true; break;
            }
        }
        if (!found) {
            // The edge was obviously not found before
            replaceFaceLabel(c1IntIndex[1], -1, cell_1);
            replaceFace = c1IntIndex[1];
        } 

        // Check if face reversal is necessary for the replacement
        if (owner_[replaceFace] == c1) {
            if (neighbour_[replaceFace] == -1) {
                // Change the owner (no flux-flipping)
                owner_[replaceFace] = newCellIndex1;
                sign = 1.0; flipOption = false;
            } else {    
                n1 = neighbour_[replaceFace];
                // This face has to be reversed (and flux is to be flipped)
                faces_[replaceFace] = faces_[replaceFace].reverseFace();
                owner_[replaceFace] = neighbour_[replaceFace];
                neighbour_[replaceFace] = newCellIndex1;
                sign = -1.0; flipOption = true;
            }
        } else {
            n1 = owner_[replaceFace];
            // Keep owner, but change neighbour (no flux-flipping)
            neighbour_[replaceFace] = newCellIndex1;
            sign = -1.0; flipOption = false;
        }
        // Modify the local-flux field
        if (flipOption && fluxInterpolation_) {
            localPhi_[replaceFace] *= -1.0;
        }        
        
        // Define attributes for the new prism cell
        newCell1[0] = replaceFace;

        // The interior face resulting from bisection; Owner = newCell[0] & Neighbour = newCell[1]
        tmpQuadFace[0] = newPtIndex1;
        tmpQuadFace[1] = nextToOtherPoint[1];
        tmpQuadFace[2] = commonEdges[0].otherVertex(nextToOtherPoint[0]);
        tmpQuadFace[3] = newPtIndex0;
        edgeToWatch[0] = newPtIndex1;
        edgeToWatch[1] = nextToOtherPoint[1];
        newFaceIndex = insertFace(-1, tmpQuadFace, newCellIndex0, newCellIndex1, edgeToWatch);
        replaceFaceLabel(-1, newFaceIndex, newCell0);
        replaceFaceLabel(-1, newFaceIndex, newCell1);
        newCell1[1] = newFaceIndex;
        
        /*
        if (fluxInterpolation_) {
            localPhi_[newFaceIndex] = newBisectFlux;
        } 
        */

        // Check for common edges among the two boundary faces
        // Find the isolated point on both boundary faces of cell[1]
        if(findCommonEdge(c1BdyFace[0],c0BdyFace[0],commonEdges[0])) {

            findCommonEdge(c1BdyFace[1],c0BdyFace[1],commonEdges[1]);
            findIsolatedPoint(c1BdyFace[0], commonEdges[0], otherPointIndex[2], nextToOtherPoint[2]);
            findIsolatedPoint(c1BdyFace[1], commonEdges[1], otherPointIndex[3], nextToOtherPoint[3]);

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

        } else {

            findCommonEdge(c1BdyFace[0],c0BdyFace[1],commonEdges[1]);
            findCommonEdge(c1BdyFace[1],c0BdyFace[0],commonEdges[0]);
            findIsolatedPoint(c1BdyFace[0], commonEdges[1], otherPointIndex[3], nextToOtherPoint[3]);
            findIsolatedPoint(c1BdyFace[1], commonEdges[0], otherPointIndex[2], nextToOtherPoint[2]);

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

        } 

        // New interior face; Owner = cell[1] & Neighbour = newCell[1]
        tmpQuadFace[0] = newPtIndex0;
        tmpQuadFace[1] = otherPointIndex[2];
        tmpQuadFace[2] = otherPointIndex[3];
        tmpQuadFace[3] = newPtIndex1;
        edgeToWatch[0] = newPtIndex1;
        edgeToWatch[1] = otherPointIndex[3];
        newFaceIndex = insertFace(-1, tmpQuadFace, c1, newCellIndex1, edgeToWatch);
        replaceFaceLabel(-1, newFaceIndex, newCell1);
        replaceFaceLabel(-1, newFaceIndex, cell_1);
        
        // Calculate fluxes for this face to satisfy zero-divergence. 
        // Assumes that no fluxes are present on boundary triangle faces 
        if (fluxInterpolation_) {
            //localPhi_[newFaceIndex] = - newBisectFlux + (sign*localPhi_[replaceFace]);
            localPhi_[newFaceIndex] = -1.0*localPhi_[replaceFace];
            if (n0!=-1) localp_[newCellIndex0] = (0.5*localp_[c0]) + (0.5*localp_[n0]);
            if (n1!=-1) localp_[newCellIndex1] = (0.5*localp_[c1]) + (0.5*localp_[n1]);
        }        

        // Second boundary face; Owner = cell[1] & Neighbour [-1]
        tmpTriFace[0] = otherPointIndex[2];
        tmpTriFace[1] = newPtIndex0;
        tmpTriFace[2] = commonEdges[0].otherVertex(nextToOtherPoint[2]);
        edgeToWatch[0] = edgeToWatch[1] = 0;
        newFaceIndex = insertFace(whichPatch(c1BdyIndex[0]), tmpTriFace, c1, -1, edgeToWatch);
        replaceFaceLabel(-1, newFaceIndex, cell_1);

        // Third boundary face; Owner = newCell[1] & Neighbour [-1] 
        tmpTriFace[0] = otherPointIndex[3];
        tmpTriFace[1] = newPtIndex1;
        tmpTriFace[2] = commonEdges[1].otherVertex(nextToOtherPoint[3]);
        edgeToWatch[0] = edgeToWatch[1] = 0;
        newFaceIndex = insertFace(whichPatch(c1BdyIndex[1]), tmpTriFace, newCellIndex1, -1, edgeToWatch);
        replaceFaceLabel(-1, newFaceIndex, newCell1);
        
        // Obtain the new cell centers and interpolate for pressure
        /*
        if (fluxInterpolation_) {
            xC0n   = cellCenter(cell_0);             
            xC1n   = cellCenter(cell_1); 
            xnewC0 = cellCenter(newCell0); 
            xnewC1 = cellCenter(newCell1);
            localp_[newCellIndex0] = localp_[c0] + ((xnewC0 - xC0o)&localGradp_[c0]);
            localp_[newCellIndex1] = localp_[c1] + ((xnewC1 - xC1o)&localGradp_[c1]);
            localp_[c0] += ((xC0n - xC0o)&localGradp_[c0]);
            localp_[c1] += ((xC1n - xC1o)&localGradp_[c1]);
        }
        */
        
        if (debug) {
            Info << "Modified Cell[0]: " << c0 << ": " << cell_0 << endl;
            forAll(cell_0, faceI)
                Info << cell_0[faceI] << ": " << faces_[cell_0[faceI]] << endl;
            Info << "New Cell[0]: " << newCellIndex0 << ": " << newCell0 << endl;
            forAll(newCell0, faceI)
                Info << newCell0[faceI] << ": " << faces_[newCell0[faceI]] << endl;
            Info << "Modified Cell[1]: " << c1 << ": " << cell_1 << endl;
            forAll(cell_1, faceI)
                Info << cell_1[faceI] << ": " << faces_[cell_1[faceI]] << endl;
            Info << "New Cell[1]: " << newCellIndex1 << ": " << newCell1 << endl;
            forAll(newCell1, faceI)
                Info << newCell1[faceI] << ": " << faces_[newCell1[faceI]] << endl;
            /*
            if (fluxInterpolation_) {
                Info << "Pressure["<< c0 << "] = " << localp_[c0] << endl;
                Info << "Pressure["<< c1 << "] = " << localp_[c1] << endl;
                Info << "Pressure["<< newCellIndex0 << "] = " << localp_[newCellIndex0] << endl;
                Info << "Pressure["<< newCellIndex1 << "] = " << localp_[newCellIndex1] << endl;
            }
            */
        }
    }    
    
    // Update the number of cells
    nCells_++;
    if (c1 != -1) nCells_++;   
}

// Method for the collapse of a quad-face in 2D
// Returns a boolean value indicating whether the collapse was valid
bool Foam::dynamicTopoFvMesh::collapseQuadFace(const label findex, face& thisFace)
{
    // This face is to be collapsed...
    if (debug) Info << nl << nl << "Face: " << findex << ": " << thisFace << " is to be collapsed. " << endl;
    
    // Local variables
    bool flipOption;
    label c0BdyIndex[2], c0IntIndex[2], c1BdyIndex[2], c1IntIndex[2];
    face  c0BdyFace[2],  c0IntFace[2],  c1BdyFace[2],  c1IntFace[2];
    edge  tmpEdge(0, 0), firstEdge, secondEdge;
    face  tmpTriFace(3);
    
    // Find the two edges from checkEdge
    edge& checkEdge = edgeToWatch_[findex];
    edgeList thisEdges = thisFace.edges();
    forAll(thisEdges,edgeI) {
        if ( (checkEdge[0] == thisEdges[edgeI][0] || checkEdge[0] == thisEdges[edgeI][1]) 
             && !(checkEdge == thisEdges[edgeI]) ) {
            firstEdge = thisEdges[edgeI];
        } else
        if ( (checkEdge[1] == thisEdges[edgeI][0] || checkEdge[1] == thisEdges[edgeI][1])
             && !(checkEdge == thisEdges[edgeI]) ) {
            secondEdge = thisEdges[edgeI];
        } else 
        if (!(checkEdge == thisEdges[edgeI])) {
            // This is the fourth edge
            tmpEdge = thisEdges[edgeI];
        }
    }

    // Build a hull of faces that are connected to each edge
    // This will also determine whether the edge lies on a boundary
    DynamicList<label> firstEdgeHull(10), secondEdgeHull(10);
    DynamicList<label> firstCells(10),    secondCells(10);
    DynamicList<label> firstTriFaces(10), secondTriFaces(10);
    bool firstEdgeBoundary  = constructPrismHull(firstEdge, findex, firstEdgeHull, firstCells, firstTriFaces, true);
    bool secondEdgeBoundary = constructPrismHull(secondEdge, findex, secondEdgeHull, secondCells, secondTriFaces, true);

    if (debug) {
        Info << "-------------------------" << endl;
        Info << "Hulls before modification" << endl;
        Info << "-------------------------" << endl;
        Info << "Cells belonging to first Edge Hull: " << firstCells << endl;
        forAll(firstCells,cellI) {
            cell &firstCurCell = cells_[firstCells[cellI]];
            Info << "Cell:: " << firstCurCell << endl;
            forAll(firstCurCell,faceI)
                Info << firstCurCell[faceI] << ": " << faces_[firstCurCell[faceI]] << endl;
        }
        Info << "First Edge Hull: " << firstEdgeHull << endl;
        forAll(firstEdgeHull,indexI)
            Info << firstEdgeHull[indexI] << ": " << faces_[firstEdgeHull[indexI]] << endl;
        Info << "Cells belonging to second Edge Hull: " << secondCells << endl;
        forAll(secondCells, cellI) {
            cell &secondCurCell = cells_[secondCells[cellI]];
            Info << "Cell:: " << secondCurCell << endl;
            forAll(secondCurCell, faceI)
                Info << secondCurCell[faceI] << ": " << faces_[secondCurCell[faceI]] << endl;
        }
        Info << "Second Edge Hull: " << secondEdgeHull << endl;
        forAll(secondEdgeHull, indexI)
            Info << secondEdgeHull[indexI] << ": " << faces_[secondEdgeHull[indexI]] << endl;                    
    } 

    // Check if this is an edge (face) that lies on a boundary
    // Might cause pinching of the boundary, and so, disabled for the moment
    if (firstEdgeBoundary && secondEdgeBoundary && debug) {
        Info << "Collapsing a boundary edge..." << endl;
        return false;
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
    findPrismFaces(findex, cell_0, c0BdyFace, c0BdyIndex, c0IntFace, c0IntIndex);
    if (c1 != -1) {
        cell &cell_1 = cells_[c1];
        findPrismFaces(findex, cell_1, c1BdyFace, c1BdyIndex, c1IntFace, c1IntIndex);
    }

    if (!firstEdgeBoundary && secondEdgeBoundary) {
        // Check whether the collapse is possible. 
        forAll(firstTriFaces, indexI) {
            if (firstTriFaces[indexI] == c0BdyIndex[0] || firstTriFaces[indexI] == c0BdyIndex[1]) continue;
            if (c1 != -1) if (firstTriFaces[indexI] == c1BdyIndex[0] || firstTriFaces[indexI] == c1BdyIndex[1]) continue;
            face &triFace = faces_[firstTriFaces[indexI]];
            forAll(triFace, pointI) {
                tmpTriFace[pointI] = triFace[pointI];
                if (triFace[pointI] == cv0) tmpTriFace[pointI] = cv2;
                if (triFace[pointI] == cv1) tmpTriFace[pointI] = cv3;
            }
            // Compute the area and check if it's zero/negative
            scalar origArea = triFaceArea(triFace);
            scalar newArea  = triFaceArea(tmpTriFace);
            if ((Foam::sign(origArea) != Foam::sign(newArea)) || mag(newArea) < SMALL) return false;
        }
        // Collapse to the second node...
        forAll(firstEdgeHull,faceI) {
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
            if (firstEdgeHull[faceI] == c0IntIndex[0]) {
                faceToKeep[0]  = c0IntIndex[1];
                faceToThrow[0] = c0IntIndex[0];
            }
            if (firstEdgeHull[faceI] == c0IntIndex[1]) {
                faceToKeep[0]  = c0IntIndex[0];
                faceToThrow[0] = c0IntIndex[1];
            }
            if (c1 != -1) {
                if (firstEdgeHull[faceI] == c1IntIndex[0]) {
                    faceToKeep[1]  = c1IntIndex[1];
                    faceToThrow[1] = c1IntIndex[0];
                }
                if (firstEdgeHull[faceI] == c1IntIndex[1]) {
                    faceToKeep[1]  = c1IntIndex[0];
                    faceToThrow[1] = c1IntIndex[1];
                }
            }
        }   
        // All triangular boundary faces also need to have point labels replaced
        forAll(firstCells,cellI) {
            cell& cellToCheck = cells_[firstCells[cellI]];
            forAll(cellToCheck,faceI) {
                face& faceToCheck = faces_[cellToCheck[faceI]];
                if (faceToCheck.size() == 3) {
                    forAll(faceToCheck,pointI) {
                        if (faceToCheck[pointI] == cv0) faceToCheck[pointI] = cv2;
                        if (faceToCheck[pointI] == cv1) faceToCheck[pointI] = cv3;
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
    } else {
        // Check whether the collapse is possible.
        forAll(secondTriFaces, indexI) {
            if (secondTriFaces[indexI] == c0BdyIndex[0] || secondTriFaces[indexI] == c0BdyIndex[1]) continue;
            if (c1 != -1) if (secondTriFaces[indexI] == c1BdyIndex[0] || secondTriFaces[indexI] == c1BdyIndex[1]) continue;
            face &triFace = faces_[secondTriFaces[indexI]];
            forAll(triFace, pointI) {
                tmpTriFace[pointI] = triFace[pointI];
                if (triFace[pointI] == cv2) tmpTriFace[pointI] = cv0;
                if (triFace[pointI] == cv3) tmpTriFace[pointI] = cv1;
            }
            // Compute the area and check if it's zero/negative
            scalar origArea = triFaceArea(triFace);
            scalar newArea  = triFaceArea(tmpTriFace);
            if ((Foam::sign(origArea) != Foam::sign(newArea)) || mag(newArea) < SMALL) return false;
        }
        // Collapse to the first node by default...
        forAll(secondEdgeHull,faceI) {
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
            if (secondEdgeHull[faceI] == c0IntIndex[0]) {
                faceToKeep[0]  = c0IntIndex[1];
                faceToThrow[0] = c0IntIndex[0];
            }
            if (secondEdgeHull[faceI] == c0IntIndex[1]) {
                faceToKeep[0]  = c0IntIndex[0];
                faceToThrow[0] = c0IntIndex[1];
            }
            if (c1 != -1) {
                if (secondEdgeHull[faceI] == c1IntIndex[0]) {
                    faceToKeep[1]  = c1IntIndex[1];
                    faceToThrow[1] = c1IntIndex[0];
                }
                if (secondEdgeHull[faceI] == c1IntIndex[1]) {
                    faceToKeep[1]  = c1IntIndex[0];
                    faceToThrow[1] = c1IntIndex[1];
                }
            }
        }
        // All triangular boundary faces also need to have point labels replaced
        forAll(secondCells, cellI) {
            cell& cellToCheck = cells_[secondCells[cellI]];
            forAll(cellToCheck, faceI) {
                face& faceToCheck = faces_[cellToCheck[faceI]];
                if (faceToCheck.size() == 3) {
                    forAll(faceToCheck, pointI) {
                        if (faceToCheck[pointI] == cv2) faceToCheck[pointI] = cv0;
                        if (faceToCheck[pointI] == cv3) faceToCheck[pointI] = cv1;
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

    if (debug) {
        Info << "------------------------" << endl;
        Info << "Hulls after modification" << endl;
        Info << "------------------------" << endl;
        Info << "Cells belonging to first Edge Hull: " << firstCells << endl;
        forAll(firstCells, cellI) {
            cell &firstCurCell = cells_[firstCells[cellI]];
            Info << "Cell:: " << firstCurCell << endl;
            forAll(firstCurCell, faceI)
                Info << firstCurCell[faceI] << ": " << faces_[firstCurCell[faceI]] << endl;
        }
        Info << "First Edge Hull: " << firstEdgeHull << endl;
        forAll(firstEdgeHull, indexI)
            Info << firstEdgeHull[indexI] << ": " << faces_[firstEdgeHull[indexI]] << endl;
        Info << "Cells belonging to second Edge Hull: " << secondCells << endl;
        forAll(secondCells, cellI) {
            cell &secondCurCell = cells_[secondCells[cellI]];
            Info << "Cell:: " << secondCurCell << endl;
            forAll(secondCurCell, faceI)
                Info << secondCurCell[faceI] << ": " << faces_[secondCurCell[faceI]] << endl;
        }
        Info << "Second Edge Hull: " << secondEdgeHull << endl;
        forAll(secondEdgeHull, indexI)
            Info << secondEdgeHull[indexI] << ": " << faces_[secondEdgeHull[indexI]] << endl;

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
        if (c1 != -1) {
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
    if (owner_[faceToThrow[0]] == c0) {
        cellCheck[0] = neighbour_[faceToThrow[0]];
        if (owner_[faceToKeep[0]] == c0) {
            if (
                   (neighbour_[faceToThrow[0]] > neighbour_[faceToKeep[0]])
                && (neighbour_[faceToKeep[0]] != -1)
               ) 
            {
                // This face is to be flipped
                faces_[faceToKeep[0]] = faces_[faceToKeep[0]].reverseFace();
                owner_[faceToKeep[0]] = neighbour_[faceToKeep[0]];
                neighbour_[faceToKeep[0]] = neighbour_[faceToThrow[0]];
                flipOption = true;
            } else {
                // Keep orientation intact, and update the owner
                owner_[faceToKeep[0]] = neighbour_[faceToThrow[0]];
                flipOption = false;
            }
        } else {
            // Keep orientation intact, and update the neighbour
            neighbour_[faceToKeep[0]] = neighbour_[faceToThrow[0]];
            flipOption = false;
        }
    } else {
        cellCheck[0] = owner_[faceToThrow[0]];
        if (neighbour_[faceToKeep[0]] == c0) {
            if (owner_[faceToThrow[0]] < owner_[faceToKeep[0]]) {   
                // This face is to be flipped
                faces_[faceToKeep[0]] = faces_[faceToKeep[0]].reverseFace();
                neighbour_[faceToKeep[0]] = owner_[faceToKeep[0]];
                owner_[faceToKeep[0]] = owner_[faceToThrow[0]];
                flipOption = true;
            } else {
                // Keep orientation intact, and update the neighbour
                neighbour_[faceToKeep[0]] = owner_[faceToThrow[0]];
                flipOption = false;
            }                
        } else {
            // Keep orientation intact, and update the owner
            owner_[faceToKeep[0]] = owner_[faceToThrow[0]];
            flipOption = false;
        }
    }
    // Modify the local-flux field
    if (flipOption && fluxInterpolation_) {
        localPhi_[faceToKeep[0]] *= -1.0;
    }      
    
    if (c1 != -1) {
        if (owner_[faceToThrow[1]] == c1) {
            cellCheck[1] = neighbour_[faceToThrow[1]];
            if (owner_[faceToKeep[1]] == c1) {
                if (
                       (neighbour_[faceToThrow[1]] > neighbour_[faceToKeep[1]])
                    && (neighbour_[faceToKeep[1]] != -1)
                   ) 
                {
                    // This face is to be flipped
                    faces_[faceToKeep[1]] = faces_[faceToKeep[1]].reverseFace();
                    owner_[faceToKeep[1]] = neighbour_[faceToKeep[1]];
                    neighbour_[faceToKeep[1]] = neighbour_[faceToThrow[1]];
                    flipOption = true;
                } else {
                    // Keep orientation intact, and update the owner
                    owner_[faceToKeep[1]] = neighbour_[faceToThrow[1]];
                    flipOption = false;
                }
            } else {
                // Keep orientation intact, and update the neighbour
                neighbour_[faceToKeep[1]] = neighbour_[faceToThrow[1]];
                flipOption = false;
            }
        } else {
            cellCheck[1] = owner_[faceToThrow[1]];
            if (neighbour_[faceToKeep[1]] == c1) {
                if (owner_[faceToThrow[1]] < owner_[faceToKeep[1]]) {
                    // This face is to be flipped
                    faces_[faceToKeep[1]] = faces_[faceToKeep[1]].reverseFace();
                    neighbour_[faceToKeep[1]] = owner_[faceToKeep[1]];
                    owner_[faceToKeep[1]] = owner_[faceToThrow[1]];
                    flipOption = true;
                } else {
                    // Keep orientation intact, and update the neighbour
                    neighbour_[faceToKeep[1]] = owner_[faceToThrow[1]];
                    flipOption = false;
                }
            } else {
                // Keep orientation intact, and update the owner
                owner_[faceToKeep[1]] = owner_[faceToThrow[1]];
                flipOption = false;
            }
        }   
        // Modify the local-flux field
        if (flipOption && fluxInterpolation_) {
            localPhi_[faceToKeep[1]] *= -1.0;
        }
    }

    // Remove the unwanted faces in the cell(s) adjacent to this face,
    // and correct the cells that contain discarded faces
    forAll(cell_0,faceI){
        if (cell_0[faceI] != findex && cell_0[faceI] != faceToKeep[0]) 
           removeFace(cell_0[faceI]);
    }
    cells_.remove(c0);
    lengthScale_.remove(c0);
    replaceFaceLabel(faceToThrow[0], faceToKeep[0], cells_[cellCheck[0]]);
    // Update the number of cells, and the reverse map
    nCells_--;
    if (c0 < nOldCells_) reverseCellMap_[c0] = -1;

    if (c1 != -1) {
        cell &cell_1 = cells_[c1];
        forAll(cell_1, faceI){
            if (cell_1[faceI] != findex && cell_1[faceI] != faceToKeep[1])
               removeFace(cell_1[faceI]);
        }
        cells_.remove(c1);
        lengthScale_.remove(c1);
        replaceFaceLabel(faceToThrow[1], faceToKeep[1], cells_[cellCheck[1]]);
        // Update the number of cells, and the reverse map
        nCells_--;
        if (c1 < nOldCells_) reverseCellMap_[c1] = -1;
    } 
    
    // Return a successful collapse
    return true;
}

// Update the mesh for motion
// This routine assumes that all boundary motions have been defined 
// and incorporated into the mesh for the current time-step.
void Foam::dynamicTopoFvMesh::updateMotion()
{     
    tmp<pointField> newPoints;
           
    if (solveForMotion_) {
        
        // Determine the kind of motion solver in use
        word solverType(dict_.lookup("solver"));

        if 
        (    
            (solverType == "laplaceCellDecomposition")
         || (solverType == "laplaceFaceDecomposition")
         || (solverType == "pseudoSolidCellDecomposition")
         || (solverType == "pseudoSolidFaceDecomposition")
        ) 
        {
            // Boundary motion specified for the tetDecompositionMotionSolver
            tetPointVectorField& motionU = const_cast<tetPointVectorField&>
                    (this->objectRegistry::lookupObject<tetPointVectorField>("motionU"));    

            // Assign boundary conditions to the motion solver
            for(label i=0; i<numPatches_; i++) 
                motionU.boundaryField()[i] == displacementPtr_[i]/this->time().deltaT().value();
        }

        if 
        (    
            (solverType == "velocityLaplacian")
        ) 
        {
            // Boundary motion specified for the fvMotionSolver
            pointVectorField& motionU = const_cast<pointVectorField&>
                    (this->objectRegistry::lookupObject<pointVectorField>("pointMotionU"));    

            // Assign boundary conditions to the motion solver
            for(label i=0; i<numPatches_; i++)
                motionU.boundaryField()[i] == displacementPtr_[i]/this->time().deltaT().value();
        }

        // Solve for motion   
        newPoints = motionPtr_->newPoints();
        //movePoints(motionPtr_->newPoints());
    } else {
        
        // Keep the points at the present state for a static mesh
        newPoints = this->points();
        
    }

    if (topoChangeFlag_) {
        // Move mesh points to state containing zero-volume cells
        movePoints(pointsZeroVolume_);
        //movePoints(points());
        resetMotion();
        setV0();        
    }
    
    // Now move points to the present state, and solve for mesh motion fluxes
    movePoints(newPoints);
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
    for(label i=0; i<numPatches_; i++) {
        oldPatchSizes_[i] = patchSizes_[i];
        oldPatchStarts_[i] = patchStarts_[i]; 
        oldPatchNMeshPoints_[i] = patchNMeshPoints_[i]; 
    }   
    
    // Print out the mesh bandwidth
    if (debug) {
        label band=0;
        const labelList& oldOwner = owner();
        const labelList& oldNeighbour = neighbour();
        forAll(oldOwner, faceI) {
            label diff = oldNeighbour[faceI] - oldOwner[faceI];
            if (diff > band) band = diff;
        }
        Info << "Mesh size: " << nCells()
            << "    Bandwidth before renumbering: " << band << endl;            
    }    
    
    // Obtain the most recent point-positions
    pointsZeroVolume_.clear(); pointsZeroVolume_.setSize(nPoints_);
    const pointField& currentPoints = this->points();
    HashList<point>::iterator pIter = meshPoints_.begin();
    HashList<point>::iterator pzvIter = pointsZeroVol_.begin(); 
    while (pIter != meshPoints_.end()) {
        pIter() = currentPoints[pIter.index()];
        pzvIter() = currentPoints[pzvIter.index()];
        pointsZeroVolume_[pIter.index()] = currentPoints[pIter.index()];
        // Update the iterators
        pIter++; pzvIter++;
    }

    // Obtain recent fluxes from the object-registry and make a local-copy
    if (fluxInterpolation_) {
        if (localPhi_.empty()) localPhi_.setSize(nFaces_,0.0);
        if (localp_.empty()) localp_.setSize(nCells_,0.0);
        /*
        if (localGradp_.empty()) localGradp_.setSize(nCells_,vector::zero);
        */
        // Copy old conservative fluxes
        surfaceScalarField& phi = const_cast<surfaceScalarField&>
                (this->objectRegistry::lookupObject<surfaceScalarField>(fluxFieldName_));
        forAll(phi.internalField(),faceI) {
            localPhi_[faceI] = phi.internalField()[faceI];
        }  
        for(label i=0; i<numPatches_; i++) {
            label start=patchStarts_[i];
            forAll(phi.boundaryField()[i],faceI)
                localPhi_[start+faceI] = phi.boundaryField()[i][faceI];
        }
        // Copy old pressure and its gradient
        volScalarField& p = const_cast<volScalarField&>
                (this->objectRegistry::lookupObject<volScalarField>("p"));
        //volVectorField gradp = fvc::grad(p);
        forAll(p.internalField(),cellI) {
            localp_[cellI] = p.internalField()[cellI];
            //localGradp_[cellI] = gradp.internalField()[cellI];
        }
    }    
    
    //== Connectivity changes ==//
    
    // Reset the flag
    topoChangeFlag_ = false;    
    
    if ( twoDMotion_ ) {        

        // 2D Edge-swapping engine
        //swap2DEdges();
        
        if (debug) Info << nl << "2D Edge Swapping complete." << endl;        
        
        // 2D Edge-bisection/collapse engine
        if ( edgeModification_ ) edgeBisectCollapse2D();
        
        if (debug) Info << nl << "2D Edge Bisection/Collapse complete." << endl;
        
    } else {
        
        if (debug) Info << nl << "3D Edge Bisection/Collapse complete." << endl;
        
        if (debug) Info << nl << "3D Edge Swapping complete." << endl;        
        
    }
    
    // Apply all pending topology changes, if necessary
    if (topoChangeFlag_) {
        
        // Reset and allocate for zero-volume points
        pointsZeroVolume_.clear(); pointsZeroVolume_.setSize(nPoints_);
        
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
        List<objectMap> cellsFromCells(0);         
        labelHashSet flipFaceFlux(0);
        labelListList pointZoneMap(0);
        labelListList faceZonePointMap(0);
        labelListList faceZoneFaceMap(0);
        labelListList cellZoneMap(0);
        pointField preMotionPoints(0);
        
        // Reorder the mesh and obtain current topological information
        reOrderMesh(points, pointsZeroVolume_, faces, owner, neighbour);  
        
        // Obtain the patch-point labels for mapping before resetting the mesh
        labelListList oldMeshPointLabels(numPatches_);
        for(label i=0; i<numPatches_; i++)
            oldMeshPointLabels[i] = this->boundaryMesh()[i].meshPoints();
        
        // Reset the mesh
        this->resetPrimitives(
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
        for(label i=0; i<numPatches_; i++) {
            const labelList& meshPointLabels = this->boundaryMesh()[i].meshPoints();
            patchNMeshPoints_[i] = meshPointLabels.size();
            patchPointMap[i].setSize(meshPointLabels.size(), -1);
            forAll(meshPointLabels, pointI) {
                // Check if the position has been maintained (This saves a search).
                // Otherwise, perform a search for the old position in the patch.
                if (pointI < oldPatchNMeshPoints_[i]) {
                    if (meshPointLabels[pointI] == oldMeshPointLabels[i][pointI]) {
                        // Good. Position is maintained. Make an entry
                        patchPointMap[i][pointI] = pointI;
                    } else {
                        // Start a linear search for the old position
                        bool foundOldPos=false;
                        forAll(oldMeshPointLabels[i],pointJ) {
                            if (oldMeshPointLabels[i][pointJ] == meshPointLabels[pointI]) {
                                patchPointMap[i][pointI] = pointJ;
                                foundOldPos=true; break;
                            }
                        }
                        if (!foundOldPos) {
                            // Couldn't find a match. Must be a new label.
                            patchPointMap[i][pointI] = -1;
                        }
                    }
                }
                //patchPointMap[i][pointI] = this->boundaryMesh()[i].meshPointMap()[meshPointLabels[pointI]];
            }
        }
        
        if (debug) {
            Info << "Patch MeshPoints: " << patchNMeshPoints_ << endl;
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
        
        // Update the motion-solver, if necessary
        if (motionPtr_.valid()) motionPtr_().updateMesh(mapper_);
        
        // Update the underlying mesh, and map all related fields
        updateMesh(mapper_);
        
        // Print out the mesh bandwidth
        if (debug) {
            label band=0;
            forAll(owner, faceI) {
                label diff = neighbour[faceI] - owner[faceI];
                if (diff > band) band = diff;
            }
            Info << "Mesh size: " << nCells()
                << "    Bandwidth after renumbering: " << band << endl;            
        }
        
        // Clear the current and reverse maps         
        pointMap_.clear(); 
        faceMap_.clear(); 
        cellMap_.clear();
        reversePointMap_.clear(); 
        reverseFaceMap_.clear(); 
        reverseCellMap_.clear();
        boundaryPatches_.clear();
        
        // Set new sizes for the reverse maps
        reversePointMap_.setSize(nPoints_);
        reverseFaceMap_.setSize(nFaces_);
        reverseCellMap_.setSize(nCells_);
    }    
    
    // Basic checks for mesh-validity
    if (debug) checkMesh(true);
    
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
