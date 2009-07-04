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

    if (findIndex(stack_,index) == -1)
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

    label loc = findIndex(stack_,index);

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
    edge& thisEdge = edges_[eIndex];

    return (points_[thisEdge[1]] - points_[thisEdge[0]]);
}

// Find the circumcenter of a face
inline point dynamicTopoFvMesh::circumCenter
(
    const label fIndex
)
{
    // Obtain references
    face& thisFace = faces_[fIndex];
    const point& a = points_[thisFace[0]];
    const point& b = points_[thisFace[1]];
    const point& c = points_[thisFace[2]];

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

// Return length-scale at an face-location in the mesh [2D]
inline void dynamicTopoFvMesh::meshFaceLengthScale
(
    const label fIndex,
    scalar& scale
)
{
    // Reset the scale first
    scale = 0.0;

    // Determine whether the face is internal
    if (whichPatch(fIndex) < 0)
    {
        scale = 0.5*
                (
                    lengthScale_[owner_[fIndex]]
                  + lengthScale_[neighbour_[fIndex]]
                );
    }
    else
    {
        scale = boundaryLengthScale(fIndex);
    }
}

// Compute length-scale at an edge-location in the mesh [3D]
inline void dynamicTopoFvMesh::meshEdgeLengthScale
(
    const label eIndex,
    scalar& scale
)
{
    // Reset the scale first
    scale = 0.0;

    labelList& eFaces = edgeFaces_[eIndex];

    // Determine whether the edge is internal
    if (whichEdgePatch(eIndex) < 0)
    {
        forAll(eFaces, faceI)
        {
#           ifdef FULLDEBUG
            // Check whether neighbour is valid
            if (neighbour_[eFaces[faceI]] == -1)
            {
                FatalErrorIn("dynamicTopoFvMesh::meshEdgeLengthScale()")
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
    }
}

// Check if a given face is a quad
inline bool dynamicTopoFvMesh::checkQuadFace(const label fIndex)
{
    return (faces_[fIndex].size() == 4);
}

// Compute the length of an edge
inline void dynamicTopoFvMesh::edgeLength
(
    const label eIndex,
    scalar& length
)
{
    edge& thisEdge = edges_[eIndex];

    length = mag(points_[thisEdge[1]] - points_[thisEdge[0]]);
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

// Return the priority edge-stack for a particular thread
inline dynamicTopoFvMesh::stack& dynamicTopoFvMesh::priorityEdgeStack
(
    const label threadID
)
{
    return priorityEdgeStack_[threadID];
}

// Return the integer threadID for a given pthread
// Return zero for single-threaded operation
inline label dynamicTopoFvMesh::self()
{
    if (threader_->multiThreaded())
    {
        for (label i = 1; i <= threader_->getNumThreads(); i++)
        {
            if (pthread_equal(structPtr_[i].pthreadID_, pthread_self()))
            {
                return i;
            }
        }
    }

    return 0;
}

// Initialize edge-stacks
inline void dynamicTopoFvMesh::initEdgeStacks()
{
    forAll(edgeStack_, stackI)
    {
        edgeStack_[stackI].clear();
        priorityEdgeStack_[stackI].clear();
    }

    label tID = 1;
    labelHashSet priorityEdges;

    // Build a list of prioritized edges for coupled patches
    if (patchCoupling_.size())
    {
        for
        (
            HashList<face>::iterator faceI = faces_(faces_.lastIndex());
            faceI.index() >= nOldInternalFaces_;
            faceI--
        )
        {
            label pIndex = whichPatch(faceI.index());

            forAllIter(Map<label>, patchCoupling_, patchI)
            {
                label mPatch = patchI.key();
                label sPatch = patchI();

                if
                (
                    (pIndex == mPatch) ||
                    (pIndex == sPatch && sPatch != -1)
                )
                {
                    labelList& faceEdges = faceEdges_[faceI.index()];

                    forAll(faceEdges, edgeI)
                    {
                        if (!priorityEdges.found(faceEdges[edgeI]))
                        {
                            priorityEdges.insert(faceEdges[edgeI]);

                            if (threader_->multiThreaded())
                            {
                                priorityEdgeStack_[tID].insert
                                (
                                    faceEdges[edgeI]
                                );

                                tID = (tID == edgeStack_.size()-1) ? 1:tID + 1;
                            }
                            else
                            {
                                priorityEdgeStack_[0].insert
                                (
                                    faceEdges[edgeI]
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    if (threader_->multiThreaded())
    {
        tID = 1;

        for
        (
            HashList<edge>::iterator iter = edges_.begin();
            iter != edges_.end();
            iter++
        )
        {
            // Ensure that this edge doesn't exist in the priority stack.
            if (priorityEdges.found(iter.index()))
            {
                continue;
            }

            edgeStack_[tID].insert(iter.index());

            tID = (tID == edgeStack_.size()-1) ? 1 : tID + 1;
        }
    }
    else
    {
        for
        (
            HashList<edge>::iterator iter = edges_.begin();
            iter != edges_.end();
            iter++
        )
        {
            if (priorityEdges.found(iter.index()))
            {
                continue;
            }

            edgeStack_[0].insert(iter.index());
        }
    }
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

        for
        (
            HashList<face>::iterator iter = faces_.begin();
            iter != faces_.end();
            iter++
        )
        {
            faceStack_[tID].insert(iter.index());

            tID = (tID == faceStack_.size()-1) ? 1 : tID + 1;
        }
    }
    else
    {
        for
        (
            HashList<face>::iterator iter = faces_.begin();
            iter != faces_.end();
            iter++
        )
        {
            faceStack_[0].insert(iter.index());
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

// Method to determine whether the edge is locally coupled
inline bool dynamicTopoFvMesh::locallyCoupledEdge
(
    const label eIndex
)
{
    const polyBoundaryMesh& bdy = boundaryMesh();

    forAll(bdy, patchI)
    {
        if (isA<processorPolyPatch>(bdy[patchI]))
        {
            continue;
        }

        if (slaveToMaster_[patchI].found(eIndex))
        {
            return true;
        }

        if (masterToSlave_[patchI].found(eIndex))
        {
            return true;
        }
    }

    return false;
}

// Method to determine whether the edge is on a processor patch
inline bool dynamicTopoFvMesh::coupledProcessorEdge
(
    const label eIndex
)
{
    const polyBoundaryMesh& bdy = boundaryMesh();

    forAll(bdy, patchI)
    {
        if (!isA<processorPolyPatch>(bdy[patchI]))
        {
            continue;
        }

        if (slaveToMaster_[patchI].found(eIndex))
        {
            return true;
        }

        if (masterToSlave_[patchI].found(eIndex))
        {
            return true;
        }
    }

    return false;
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

// Check for the occurrence of a label in the list
inline bool dynamicTopoFvMesh::foundInList
(
    const label item,
    const labelList& list,
    label& index
)
{
    forAll(list, itemI)
    {
        if (list[itemI] == item)
        {
            index = itemI;

            return true;
        }
    }

    index = -1;

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

// Utility method to compute the minimum quality of a vertex hull
inline void dynamicTopoFvMesh::computeMinQuality
(
    const label eIndex,
    scalar& minQuality,
    bool priority,
    bool compareCoupled
)
{
    minQuality = GREAT;
    scalar cQuality = 0.0;

    // Obtain a reference to this edge and corresponding edgePoints
    edge& edgeToCheck = edges_[eIndex];
    labelList& hullVertices = edgePoints_[eIndex];

    // Assume that hullVertices is ordered
    // in a CCW order around edgeToCheck[0]
    forAll(hullVertices, indexI)
    {
        // Compute quality as well
        if (indexI != 0)
        {
            cQuality = (*tetMetric_)
            (
                points_[edgeToCheck[0]],
                points_[hullVertices[indexI-1]],
                points_[edgeToCheck[1]],
                points_[hullVertices[indexI]]
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

    // If this is a prioritized edge, fill the buffer with the quality value.
    if
    (
        priority &&
        patchCoupling_.size() &&
        whichEdgePatch(eIndex) > -1
    )
    {
        label i = -1;

        const polyBoundaryMesh& bdy = boundaryMesh();

        // Search for this edge on all coupled master patches
        forAll(mList_, patchI)
        {
            if (masterToSlave_[patchI].found(eIndex))
            {
                // Find the position in the list
                foundInList(eIndex, mList_[patchI], i);

                sendMinQualityBuffer_[patchI][i] = minQuality;
            }
        }

        // Search for this edge on non-processor coupled patches
        forAll(sList_, patchI)
        {
            if (isA<processorPolyPatch>(bdy[patchI]))
            {
                continue;
            }

            if (slaveToMaster_[patchI].found(eIndex))
            {
                // Find the position in the list
                foundInList(eIndex, sList_[patchI], i);

                recvMinQualityBuffer_[patchI][i] = minQuality;
            }
        }
    }

    // Check if a coupled comparison is required.
    if
    (
        compareCoupled &&
        patchCoupling_.size() &&
        whichEdgePatch(eIndex) > -1
    )
    {
        // Check the received buffer for quality
        scalar rQ = GREAT;

        label i = -1;

        const polyBoundaryMesh& bdy = boundaryMesh();

        // Search for this edge on all coupled master patches
        forAll(mList_, patchI)
        {
            if (masterToSlave_[patchI].found(eIndex))
            {
                // Find the position in the list
                foundInList(eIndex, mList_[patchI], i);

                rQ = recvMinQualityBuffer_[patchI][masterToSlave_[patchI][i]];

                minQuality = rQ < minQuality ? rQ : minQuality;
            }
        }

        // Search for this edge on non-processor coupled patches
        forAll(sList_, patchI)
        {
            if (isA<processorPolyPatch>(bdy[patchI]))
            {
                continue;
            }

            if (slaveToMaster_[patchI].found(eIndex))
            {
                // Find the position in the list
                foundInList(eIndex, sList_[patchI], i);

                rQ = sendMinQualityBuffer_[patchI][slaveToMaster_[patchI][i]];

                minQuality = rQ < minQuality ? rQ : minQuality;
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
