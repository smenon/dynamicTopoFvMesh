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

#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Push items on to the stack
inline void dynamicTopoFvMesh::stack::push(const label index)
{
    stackMutex_.lock();

    if (findIndex(index) == -1)
    {
        stack_.append(index);
    }

    stackMutex_.unlock();
}

//- Insert item onto stack (no checking)
inline void dynamicTopoFvMesh::stack::insert(const label index)
{
    stack_.append(index);
}

// Pop an item off the stack
inline label dynamicTopoFvMesh::stack::pop()
{
    stackMutex_.lock();

    const label index = stack_.remove();

    stackMutex_.unlock();

    return index;
}

// Remove a specific item off the stack
inline void dynamicTopoFvMesh::stack::remove(const label index)
{
    stackMutex_.lock();

    label loc = findIndex(index);

    if (loc != -1)
    {
        // Create a new list
        labelList newList(stack_.size() - 1);

        label n = 0;

        // Copy items upto location
        for(label i = 0; i < loc; i++)
        {
            newList[n++] = stack_[i];
        }

        // Copy items from location
        for(label i = (loc + 1); i < stack_.size(); i++)
        {
            newList[n++] = stack_[i];
        }

        // Overwrite
        stack_ = newList;
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

//- Return if a certain index was found in the stack
inline label dynamicTopoFvMesh::stack::findIndex(const label index)
{
    forAll(stack_, itemI)
    {
        if (stack_[itemI] == index)
        {
            return itemI;
        }
    }

    return -1;
}

//- Access the stack as a list
inline label dynamicTopoFvMesh::stack::operator[]
(
    const label index
)
{
    return stack_[index];
}

// Obtain the tangent-to-edge vector
inline vector dynamicTopoFvMesh::tangentToEdge
(
    const label eIndex
)
{
    const edge& thisEdge = edges_[eIndex];

    return (points_[thisEdge[1]] - points_[thisEdge[0]]);
}

// Find the circumcenter of a face
inline point dynamicTopoFvMesh::circumCenter
(
    const label fIndex
)
{
    // Obtain references
    const face& thisFace = faces_[fIndex];
    const point a = points_[thisFace[0]];
    const point b = points_[thisFace[1]];
    const point c = points_[thisFace[2]];

    scalar d1 =  (c - a)&(b - a);
    scalar d2 = -(c - b)&(b - a);
    scalar d3 =  (c - a)&(c - b);

    scalar c1 = d2*d3;
    scalar c2 = d3*d1;
    scalar c3 = d1*d2;

    scalar cd = c1 + c2 + c3;

    if (cd < VSMALL)
    {
        FatalErrorIn
        (
            "dynamicTopoFvMesh::testDelaunay(const label fIndex) "
        ) << nl << " Encountered a co-linear set of points: " << nl
                << " Point a :: " << thisFace[0] << ": " << a << nl
                << " Point b :: " << thisFace[1] << ": " << b << nl
                << " Point c :: " << thisFace[2] << ": " << c << nl
                << abort(FatalError);
    }

    return ((c2 + c3)*a + (c3 + c1)*b + (c1 + c2)*c)/(2*cd);
}

// Find the center of a triangle face
inline point dynamicTopoFvMesh::triFaceCenter
(
    const label fIndex
)
{
    face& thisFace = faces_[fIndex];

    return
    (
        points_[thisFace[0]]
      + points_[thisFace[1]]
      + points_[thisFace[2]]
    )/3.0;
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
    vector v = points_[triFace[1]] - points_[triFace[0]];
    vector w = points_[triFace[2]] - points_[triFace[0]];

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

// Find the volume of a tetrahedron.
// The function assumes points (a-b-c)
// are in counter-clockwise fashion when viewed from d.
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

// Find the volume sign of a tetrahedron
// The function assumes points (a-b-c)
// are in counter-clockwise fashion when viewed from d;
// (to obtain positive orientation)
inline label dynamicTopoFvMesh::tetVolumeSign
(
    const point& a,
    const point& b,
    const point& c,
    const point& d,
    const scalar tolerance
)
{
    scalar det = ((b - a) ^ (c - a)) & (d - a);

    if (mag(det) < tolerance)
    {
        return 0;
    }

    if (Foam::sign(det) < 0.0)
    {
        return -1;
    }

    return 1;
}

// Return the appropriate length-scale for boundary face
inline scalar dynamicTopoFvMesh::boundaryLengthScale
(
    const label fIndex
)
{
    label bFacePatch = whichPatch(fIndex);

    const polyBoundaryMesh& boundary = polyMesh::boundaryMesh();

    // Check fixed length-scale patches
    // If the value is negative, average face length-scales.
    if (fixedPatches_.found(boundary[bFacePatch].name()))
    {
        scalar dictValue =
        (
            fixedPatches_[boundary[bFacePatch].name()][0].scalarToken()
        );

        if (dictValue > 0.0)
        {
            return dictValue;
        };
    }

    const polyPatch& bdyPatch = boundary[bFacePatch];

    if
    (
        (freePatches_.found(bdyPatch.name())) ||
        (bdyPatch.type() == "processor") ||
        (bdyPatch.type() == "cyclic") ||
        (bdyPatch.type() == "wedge") ||
        (bdyPatch.type() == "empty") ||
        (bdyPatch.type() == "symmetryPlane")
    )
    {
        return lengthScale_[owner_[fIndex]];
    }

    if (twoDMesh_)
    {
        label eIndex = getTriBoundaryEdge(fIndex);

        const edge& e = edges_[eIndex];

        return mag(points_[e[0]] - points_[e[1]]);
    }
    else
    {
        // Average edge-lengths for this face
        scalar edgeLength = 0.0;

        const labelList& fEdges = faceEdges_[fIndex];

        forAll(fEdges, edgeI)
        {
            const edge& e = edges_[fEdges[edgeI]];

            edgeLength += mag(points_[e[0]] - points_[e[1]]);
        }

        return (edgeLength/fEdges.size());
    }
}

// Return length-scale at an face-location in the mesh [2D]
inline scalar dynamicTopoFvMesh::faceLengthScale
(
  const label fIndex
)
{
    // Reset the scale first
    scalar scale = 0.0;

    // Determine whether the face is internal
    if (whichPatch(fIndex) < 0)
    {
        scale =
        (
            0.5 *
            (
                lengthScale_[owner_[fIndex]]
                + lengthScale_[neighbour_[fIndex]]
            )
        );
    }
    else
    {
        scale = boundaryLengthScale(fIndex);
    }

    return scale;
}

// Compute length-scale at an edge-location in the mesh [3D]
inline scalar dynamicTopoFvMesh::edgeLengthScale
(
    const label eIndex
)
{
    // Reset the scale first
    scalar scale = 0.0;

    const labelList& eFaces = edgeFaces_[eIndex];

    label edgePatch = whichEdgePatch(eIndex);

    // Determine whether the edge is internal
    if (edgePatch < 0)
    {
        forAll(eFaces, faceI)
        {
#           ifdef FULLDEBUG
            // Check whether neighbour is valid
            if (neighbour_[eFaces[faceI]] == -1)
            {
                FatalErrorIn("dynamicTopoFvMesh::edgeLengthScale()")
                    << "Face: " << eFaces[faceI]
                    << ": " << faces_[eFaces[faceI]]
                    << " is not internal, while edge: "
                    << eIndex << ": " << edges_[eIndex] << " is."
                    << abort(FatalError);
            }
#           endif

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

        // If proximity-based refinement is requested,
        // test the proximity to the nearest non-neighbour.
        if (proximityPatches_.found(boundaryMesh()[edgePatch].name()))
        {
            // Obtain the edge location

        }

        // If curvature-based refinement is requested,
        // test the variation in face-normal directions.
        if (curvaturePatches_.found(boundaryMesh()[edgePatch].name()))
        {
            // Obtain face-normals for both faces.
            label count = 0;
            FixedList<vector, 2> fNorm;

            forAll(eFaces, faceI)
            {
                if (neighbour_[eFaces[faceI]] == -1)
                {
                    // Obtain the normal.
                    fNorm[count] = triFaceNormal(faces_[eFaces[faceI]]);

                    // Normalize it.
                    fNorm[count] /= mag(fNorm[count]);

                    count++;
                }
            }

            scalar deviation = (fNorm[0] & fNorm[1]);

            if (mag(deviation) < curvatureDeviation_)
            {
                // Get the edge-length.
                scalar length = edgeLength(eIndex);

                if (debug > 2 && self() == 0)
                {
                    Info << "Deviation: " << deviation << nl
                         << "curvatureDeviation: " << curvatureDeviation_
                         << ", Edge: " << eIndex << ", Length: " << length
                         << ", Scale: " << scale << nl
                         << " Half-length: " << (0.5*length) << nl
                         << " MinRatio: " << (ratioMin_*scale)
                         << endl;
                }

                // Would dividing this edge make it too small?
                if
                (
                    ((0.5*length) > (ratioMin_*scale)) &&
                    ((0.5*length) > minLengthScale_)
                )
                {
                    // Smoosh the return value for bisection.
                    scalar newScale = (length - SMALL)/ratioMax_;

                    return newScale;
                }
            }
        }
    }

    if (coupledModification_)
    {
        // Check if this is a coupled processor edge.
        // - Locally coupled edges use length-scale information only
        //   from the master patch.
        if (processorCoupledEdge(eIndex))
        {
            // Check patchSubMeshes for length-scale info.

        }
    }

    return scale;
}

// Check for face bisection
inline bool dynamicTopoFvMesh::checkFaceBisection
(
    const label fIndex
)
{
    // Measure the boundary edge-length of the face in question
    scalar length = edgeLength(getTriBoundaryEdge(fIndex));

    // Determine the length-scale at this face
    scalar scale = faceLengthScale(fIndex);

    if (length > ratioMax_*scale)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// Check for face collapse
inline bool dynamicTopoFvMesh::checkFaceCollapse
(
    const label fIndex
)
{
    // Measure the boundary edge-length of the face in question
    scalar length = edgeLength(getTriBoundaryEdge(fIndex));

    // Determine the length-scale at this face
    scalar scale = faceLengthScale(fIndex);

    if (length < ratioMin_*scale)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// Check for edge bisection
inline bool dynamicTopoFvMesh::checkEdgeBisection
(
    const label eIndex
)
{
    // Measure the edge-length
    scalar length = edgeLength(eIndex);

    // Determine the length-scale at this point in the mesh
    scalar scale = edgeLengthScale(eIndex);

    if (length > ratioMax_*scale)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// Check for edge collapse
inline bool dynamicTopoFvMesh::checkEdgeCollapse
(
    const label eIndex
)
{
    // Measure the edge-length
    scalar length = edgeLength(eIndex);

    // Determine the length-scale at this point in the mesh
    scalar scale = edgeLengthScale(eIndex);

    if (length < ratioMin_*scale)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// Check if a given face is a quad
inline bool dynamicTopoFvMesh::checkQuadFace(const label fIndex)
{
    return (faces_[fIndex].size() == 4);
}

// Compute the length of an edge
inline scalar dynamicTopoFvMesh::edgeLength
(
    const label eIndex
)
{
    const edge& thisEdge = edges_[eIndex];

    return mag(points_[thisEdge[1]] - points_[thisEdge[0]]);
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

// Return the cell-stack for a particular thread
inline dynamicTopoFvMesh::stack& dynamicTopoFvMesh::cellStack
(
    const label threadID
)
{
    return cellStack_[threadID];
}

// Return the integer threadID for a given pthread
// Return zero for single-threaded operation
inline label dynamicTopoFvMesh::self()
{
    if (threader_->multiThreaded())
    {
        for (label i = 1; i <= threader_->getNumThreads(); i++)
        {
            if (pthread_equal(structPtr_[i].ID(), pthread_self()))
            {
                return i;
            }
        }
    }

    return 0;
}

// Initialize cell-stacks
inline void dynamicTopoFvMesh::initCellStacks()
{
    forAll(cellStack_, stackI)
    {
        cellStack_[stackI].clear();
    }

    if (threader_->multiThreaded())
    {
        label tID = 1;

        forAll(cells_, cellI)
        {
            if (cells_[cellI].size())
            {
                cellStack_[tID].insert(cellI);

                tID = (tID == cellStack_.size() - 1) ? 1 : tID + 1;
            }
        }
    }
    else
    {
        forAll(cells_, cellI)
        {
            if (cells_[cellI].size())
            {
                cellStack_[0].insert(cellI);
            }
        }
    }
}

// Initialize edge-stacks
inline void dynamicTopoFvMesh::initEdgeStacks()
{
    forAll(edgeStack_, stackI)
    {
        edgeStack_[stackI].clear();
    }

    if (threader_->multiThreaded())
    {
        label tID = 1;

        // For coupled meshes, avoid certain edges.
        if (patchCoupling_.size() || procIndices_.size())
        {
            forAll(edges_, edgeI)
            {
                if
                (
                   !entitiesToAvoid_.found(edgeI) &&
                    edgeFaces_[edgeI].size()
                )
                {
                    edgeStack_[tID].insert(edgeI);

                    tID = (tID == edgeStack_.size() - 1) ? 1 : tID + 1;
                }
            }
        }
        else
        {
            forAll(edges_, edgeI)
            {
                if (edgeFaces_[edgeI].size())
                {
                    edgeStack_[tID].insert(edgeI);

                    tID = (tID == edgeStack_.size() - 1) ? 1 : tID + 1;
                }
            }
        }
    }
    else
    {
        // For coupled meshes, avoid certain edges.
        if (patchCoupling_.size() || procIndices_.size())
        {
            forAll(edges_, edgeI)
            {
                if
                (
                   !entitiesToAvoid_.found(edgeI) &&
                   edgeFaces_[edgeI].size()
                )
                {
                    edgeStack_[0].insert(edgeI);
                }
            }
        }
        else
        {
            forAll(edges_, edgeI)
            {
                if (edgeFaces_[edgeI].size())
                {
                    edgeStack_[0].insert(edgeI);
                }
            }
        }
    }
}

// Initialize the coupled edge-stack
inline void dynamicTopoFvMesh::initCoupledEdgeStack()
{
    // Clear existing lists/stacks.
    edgeStack(0).clear();

    const polyBoundaryMesh& boundary = boundaryMesh();

    // Loop though boundary faces and check whether
    // they belong to master/slave coupled patches.
    for (label faceI = nOldInternalFaces_; faceI < faces_.size(); faceI++)
    {
        label pIndex = whichPatch(faceI);

        if (pIndex == -1)
        {
            continue;
        }

        // Check if this is a locally coupled master face.
        if (patchCoupling_.found(pIndex))
        {
            const labelList& mfEdges = faceEdges_[faceI];

            forAll(mfEdges, edgeI)
            {
                // Add this to the coupled modification stack.
                edgeStack(0).push(mfEdges[edgeI]);
            }
        }

        // Check if this is a processor patch.
        if (isA<processorPolyPatch>(boundary[pIndex]))
        {
            const labelList& mfEdges = faceEdges_[faceI];

            // Check if this is a master processor patch.
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[pIndex])
            );

            label neiProcID = pp.neighbProcNo();

            if (neiProcID > Pstream::myProcNo())
            {
                // Add this to the coupled modification stack.
                forAll(mfEdges, edgeI)
                {
                    edgeStack(0).push(mfEdges[edgeI]);
                }
            }
        }
    }
}

// Build a list of entities that need to be avoided
// by regular topo-changes.
inline void dynamicTopoFvMesh::buildEntitiesToAvoid()
{
    if (!twoDMesh_)
    {
        entitiesToAvoid_.clear();

        // Build a set of entities to avoid during regular modifications,
        // and build a master stack for coupled modifications.
        const polyBoundaryMesh& boundary = boundaryMesh();

        // Determine locally coupled slave patches.
        labelHashSet localSlavePatches;

        forAllIter(Map<coupledPatchInfo>, patchCoupling_, patchI)
        {
            localSlavePatches.insert(patchI().slaveIndex());
        }

        // Loop through boundary faces and check whether
        // they belong to master/slave coupled patches.
        for (label faceI = nOldInternalFaces_; faceI < faces_.size(); faceI++)
        {
            label pIndex = whichPatch(faceI);

            if (pIndex == -1)
            {
                continue;
            }

            // Check if this is a locally coupled master face.
            if (patchCoupling_.found(pIndex))
            {
                const labelList& mfEdges = faceEdges_[faceI];

                forAll(mfEdges, edgeI)
                {
                    if (!entitiesToAvoid_.found(mfEdges[edgeI]))
                    {
                        // Avoid this edge during regular modification.
                        entitiesToAvoid_.insert(mfEdges[edgeI]);
                    }
                }
            }

            // Check if this is a locally coupled slave face.
            if (localSlavePatches.found(pIndex))
            {
                const labelList& mfEdges = faceEdges_[faceI];

                forAll(mfEdges, edgeI)
                {
                    if (!entitiesToAvoid_.found(mfEdges[edgeI]))
                    {
                        // Avoid this edge during regular modification.
                        entitiesToAvoid_.insert(mfEdges[edgeI]);
                    }
                }
            }

            // Check if this is a processor patch.
            if (isA<processorPolyPatch>(boundary[pIndex]))
            {
                const labelList& mfEdges = faceEdges_[faceI];

                // Add to the stack.
                forAll(mfEdges, edgeI)
                {
                    if (!entitiesToAvoid_.found(mfEdges[edgeI]))
                    {
                        // Avoid this edge during regular modification.
                        entitiesToAvoid_.insert(mfEdges[edgeI]);
                    }
                }
            }
        }
    }
}

// Check whether the specified edge is a coupled master edge.
inline bool dynamicTopoFvMesh::isCoupledMaster
(
    const label eIndex
)
{
    if (!coupledModification_)
    {
        return true;
    }

    return locallyCoupledEdge(eIndex);
}

// Initialize face-stacks
inline void dynamicTopoFvMesh::initFaceStacks()
{
    forAll(faceStack_, stackI)
    {
        faceStack_[stackI].clear();
    }

    if (threader_->multiThreaded())
    {
        label tID = 1;

        forAll(faces_, faceI)
        {
            if (faces_[faceI].size())
            {
                faceStack_[tID].insert(faceI);

                tID = (tID == faceStack_.size() - 1) ? 1 : tID + 1;
            }
        }
    }
    else
    {
        forAll(faces_, faceI)
        {
            if (faces_[faceI].size())
            {
                faceStack_[0].insert(faceI);
            }
        }
    }
}

// Method to determine the old boundary patch index for a given face
// Similar to the polyBoundaryMesh routine, but works on local information
inline label dynamicTopoFvMesh::whichPatch
(
    const label index
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
    const label index
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
            (
                (list[itemI] == labelA && list[nextI] == labelB) ||
                (list[itemI] == labelB && list[nextI] == labelA)
            ) &&
           !found
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
    labelList newList(list.size() - 1);

#   ifdef FULLDEBUG
    // Check that the item actually exists.
    if (findIndex(list, item) == -1)
    {
        FatalErrorIn
        (
            "label dynamicTopoFvMesh::sizeDownList()"
        )   << "Cannot find " << item << " in list: " << list << nl
            << abort(FatalError);
    }
#   endif

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

// Utility method to compute the minimum quality of a vertex hull
inline scalar dynamicTopoFvMesh::computeMinQuality
(
    const label eIndex
)
{
    scalar minQuality = GREAT;
    scalar cQuality = 0.0;

    // Obtain a reference to this edge and corresponding edgePoints
    const edge& edgeToCheck = edges_[eIndex];
    const labelList& hullVertices = edgePoints_[eIndex];

    // Obtain point references
    const point& a = points_[edgeToCheck[0]];
    const point& c = points_[edgeToCheck[1]];

    if (whichEdgePatch(eIndex) < 0)
    {
        // Internal edge.
        forAll(hullVertices, indexI)
        {
            label prevIndex = hullVertices.rcIndex(indexI);

            // Pick vertices off the list
            const point& b = points_[hullVertices[prevIndex]];
            const point& d = points_[hullVertices[indexI]];

            // Compute the quality
            cQuality = (*tetMetric_)(a, b, c, d);

            // Check if the quality is worse
            minQuality = cQuality < minQuality ? cQuality : minQuality;
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

            // Compute the quality
            cQuality = (*tetMetric_)(a, b, c, d);

            // Check if the quality is worse
            minQuality = cQuality < minQuality ? cQuality : minQuality;
        }
    }

    if (coupledModification_)
    {
        if (locallyCoupledEdge(eIndex))
        {
            // Compute the minimum quality of the slave edge as well.
            label slaveIndex = -1;

            // Determine the slave index.
            forAllIter(Map<coupledPatchInfo>, patchCoupling_, patchI)
            {
                if ((slaveIndex = patchI().findSlaveIndex(eIndex)) > -1)
                {
                    break;
                }
            }

            // Temporarily turn off coupledModification
            coupledModification_ = false;

            scalar slaveQuality = computeMinQuality(slaveIndex);

            minQuality = slaveQuality < minQuality ? slaveQuality : minQuality;

            // Turn it back on.
            coupledModification_ = true;
        }
        else
        if (processorCoupledEdge(eIndex))
        {
            // Compute the minimum quality across patchSubMeshes.

        }
    }

    // Ensure that the mesh is valid
    if (minQuality < 0.0)
    {
        if (debug > 3)
        {
            // Write out cells for post processing.
            labelHashSet iCells;

            labelList& eFaces = edgeFaces_[eIndex];

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

            writeVTK(Foam::name(eIndex) + "_inverted", iCells.toc());
        }

        FatalErrorIn("dynamicTopoFvMesh::computeMinQuality()")
            << "Encountered negative cell-quality!" << nl
            << "Edge: " << eIndex << ": " << edgeToCheck << nl
            << "EdgePoints: " << hullVertices << nl
            << "Minimum Quality: " << minQuality
            << abort(FatalError);
    }

    return minQuality;
}

// Utility method to compute the quality of a vertex hull
// around an edge after bisection.
inline scalar dynamicTopoFvMesh::computeBisectionQuality
(
    const label eIndex
)
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

            // Compute the quality
            cQuality = (*tetMetric_)(a, b, midPoint, d);

            // Check if the quality is worse
            minQuality = cQuality < minQuality ? cQuality : minQuality;
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

            // Compute the quality
            cQuality = (*tetMetric_)(a, b, midPoint, d);

            // Check if the quality is worse
            minQuality = cQuality < minQuality ? cQuality : minQuality;
        }
    }

    // Ensure that the mesh is valid
    if (minQuality < 0.0)
    {
        if (debug > 3)
        {
            // Write out cells for post processing.
            labelHashSet iCells;

            labelList& eFaces = edgeFaces_[eIndex];

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

            writeVTK(Foam::name(eIndex) + "_inverted", iCells.toc());
        }

        FatalErrorIn("dynamicTopoFvMesh::computeBisectionQuality()")
            << "Negative cell-quality after bisection!" << nl
            << "Edge: " << eIndex << ": " << edgeToCheck << nl
            << "EdgePoints: " << hullVertices << nl
            << "Minimum Quality: " << minQuality << nl
            << "Mid point: " << midPoint
            << abort(FatalError);
    }

    return minQuality;
}

// Method to determine whether the master face is locally coupled
inline bool dynamicTopoFvMesh::locallyCoupledFace
(
    const label fIndex
)
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    label patch = whichPatch(fIndex);

    if (patch == -1)
    {
        FatalErrorIn("dynamicTopoFvMesh::locallyCoupledFace()")
            << "Face: " << fIndex
            << " is not on a boundary and therefore, "
            << " cannot possibly be coupled."
            << abort(FatalError);
    }

    // Processor checks receive priority.
    if (isA<processorPolyPatch>(boundary[patch]))
    {
        return false;
    }

    // Check coupled master patches.
    if (patchCoupling_.found(patch))
    {
        return true;
    }
    else
    {
        return false;
    }
}

// Method to determine whether the master edge is locally coupled.
inline bool dynamicTopoFvMesh::locallyCoupledEdge
(
    const label eIndex,
    bool checkSlaves
)
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    labelList& eFaces = edgeFaces_[eIndex];

    bool foundPatch = false;

    // Search for boundary faces, and determine boundary type.
    forAll(eFaces, faceI)
    {
        if (neighbour_[eFaces[faceI]] == -1)
        {
            label patch = whichPatch(eFaces[faceI]);

            // Processor checks receive priority.
            if (isA<processorPolyPatch>(boundary[patch]))
            {
                return false;
            }

            // Check coupled master patches.
            if (patchCoupling_.found(patch))
            {
                foundPatch = true;
            }

            if (checkSlaves)
            {
                // Check on slave patches as well.
                forAllIter
                (
                    Map<coupledPatchInfo>::iterator,
                    patchCoupling_,
                    iterI
                )
                {
                    if (iterI().slaveIndex() == patch)
                    {
                        foundPatch = true;
                    }
                }
            }
        }
    }

    // Could not find any faces on locally coupled patches.
    return foundPatch;
}

// Method to determine whether the master face is on a processor
inline bool dynamicTopoFvMesh::processorCoupledFace
(
    const label fIndex
)
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    label patch = whichPatch(fIndex);

    if (patch == -1)
    {
        FatalErrorIn("dynamicTopoFvMesh::locallyCoupledFace()")
            << "Face: " << fIndex
            << " is not on a boundary and therefore, "
            << " cannot possibly be coupled."
            << abort(FatalError);
    }

    if (isA<processorPolyPatch>(boundary[patch]))
    {
        return true;
    }
    else
    {
        return false;
    }
}

// Method to determine whether the master edge is coupled at a processor.
inline bool dynamicTopoFvMesh::processorCoupledEdge
(
    const label eIndex
)
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    labelList& eFaces = edgeFaces_[eIndex];

    // Search for boundary faces, and determine boundary type.
    forAll(eFaces, faceI)
    {
        if (neighbour_[eFaces[faceI]] == -1)
        {
            label patch = whichPatch(eFaces[faceI]);

            if (isA<processorPolyPatch>(boundary[patch]))
            {
                return true;
            }
        }
    }

    // Could not find any faces on processor patches.
    return false;
}

// Obtain the maximum coupled patch index
inline label dynamicTopoFvMesh::getMaxCouplingIndex()
{
    // Allocate a size for coupled patches
    label index = -1;

    forAllIter(Map<coupledPatchInfo>, patchCoupling_, patchI)
    {
        index = index > patchI.key() ? index : patchI.key();
    }

    return index;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
