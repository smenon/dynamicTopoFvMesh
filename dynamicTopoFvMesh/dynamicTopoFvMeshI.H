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

\*---------------------------------------------------------------------------*/

#include "tetPointRef.H"
#include "processorPolyPatch.H"
#include "lengthScaleEstimator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return a non-const reference to the lengthScaleEstimator
inline lengthScaleEstimator&
dynamicTopoFvMesh::lengthEstimator()
{
    if (!lengthEstimator_.valid())
    {
        FatalErrorIn
        (
            "inline lengthScaleEstimator& "
            "dynamicTopoFvMesh::lengthEstimator()"
        ) << nl
          << " Invalid request. Length scale estimator was not loaded. "
          << abort(FatalError);
    }

    return lengthEstimator_();
}


// Return a const reference to the lengthScaleEstimator
inline const lengthScaleEstimator&
dynamicTopoFvMesh::lengthEstimator() const
{
    if (!lengthEstimator_.valid())
    {
        FatalErrorIn
        (
            "inline const lengthScaleEstimator& "
            "dynamicTopoFvMesh::lengthEstimator() const"
        ) << nl
          << " Invalid request. Length scale estimator was not loaded. "
          << abort(FatalError);
    }

    return lengthEstimator_();
}


// Return a const reference to the multiThreader
inline const multiThreader&
dynamicTopoFvMesh::threader() const
{
    if (!threader_.valid())
    {
        FatalErrorIn
        (
            "inline const multiThreader& "
            "dynamicTopoFvMesh::threader() const"
        ) << nl
          << " Invalid request. multiThreader was not loaded. "
          << abort(FatalError);
    }

    return threader_();
}


// Does the mesh perform edge refinement?
inline bool dynamicTopoFvMesh::edgeRefinement() const
{
    return edgeRefinement_;
}


// Return a reference to the entity mutexes.
// The index 'entity' ranges from 0 to 3 for point/edge/face/cell.
inline const Mutex& dynamicTopoFvMesh::entityMutex
(
    const label entity
) const
{
    return entityMutex_[entity];
}


// Return the edge index for a provided edge
inline label dynamicTopoFvMesh::getEdgeIndex
(
    const edge& edgeToCheck
) const
{
    if (twoDMesh_)
    {
        // No efficient search method. Loop through all edges.
        forAll(edges_, edgeI)
        {
            if (edges_[edgeI] == edgeToCheck)
            {
                return edgeI;
            }
        }
    }
    else
    {
        // Look througg pointEdges list
        const labelList& pEdges = pointEdges_[edgeToCheck.start()];

        forAll(pEdges, edgeI)
        {
            if (edges_[pEdges[edgeI]] == edgeToCheck)
            {
                return pEdges[edgeI];
            }
        }
    }

    // Probably not a valid edge.
    FatalErrorIn
    (
        "inline label dynamicTopoFvMesh::getEdgeIndex"
        "(const edge& edgeToCheck) const"
    )
        << "Could not find an appropriate edge index for edge:"
        << edgeToCheck
        << abort(FatalError);

    return -1;
}


// Return length-scale at an face-location in the mesh [2D]
inline scalar dynamicTopoFvMesh::faceLengthScale
(
    const label fIndex
) const
{
    // Reset the scale first
    scalar scale = 0.0;

    label facePatch = whichPatch(fIndex);

    // Determine whether the face is internal
    if (facePatch < 0)
    {
#       ifdef FULLDEBUG
        // Check whether neighbour is valid
        if (neighbour_[fIndex] == -1)
        {
            FatalErrorIn
            (
                "inline scalar dynamicTopoFvMesh::faceLengthScale"
                "(const label fIndex) const"
            )
                << nl << "Face: " << fIndex
                << ": " << faces_[fIndex]
                << " is not internal."
                << abort(FatalError);
        }
#       endif

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
        // Fetch the fixed-length scale
        scale = lengthEstimator().fixedLengthScale(fIndex, facePatch);

        // If this is a floating face, pick the owner length-scale
        if (lengthEstimator().isFreePatch(facePatch))
        {
            scale = lengthScale_[owner_[fIndex]];
        }

        // If proximity-based refinement is requested,
        // test the proximity to the nearest non-neighbour.
        if (lengthEstimator().isProximityPatch(facePatch))
        {
            label proximityFace = -1;

            // Perform a proximity-check.
            scalar distance = testProximity(fIndex, proximityFace);

            if (debug > 3 && self() == 0)
            {
                if
                (
                    (proximityFace > -1) &&
                    ((distance / 5.0) < scale)
                )
                {
                    Info << " Closest opposing face detected for face: " << nl
                         << '\t' << fIndex
                         << " :: " << faces_[fIndex]
                         << " was face:\n"
                         << '\t' << proximityFace
                         << " :: " << polyMesh::faces()[proximityFace] << nl
                         << " with distance: " << distance
                         << endl;
                }
            }

            scale =
            (
                Foam::min
                (
                    scale,
                    ((distance / 3.0) - SMALL)/lengthEstimator().ratioMax()
                )
            );
        }

        // Limit scales if necessary
        lengthEstimator().limitScale(scale);
    }

    if (coupledModification_)
    {
        // Check if this is a coupled processor face.
        // - Locally coupled faces use length-scale information only
        //   from the master patch.
        if (processorCoupledEntity(fIndex))
        {
            // Check patchSubMeshes for length-scale info.

        }
    }

    return scale;
}


// Compute length-scale at an edge-location in the mesh [3D]
inline scalar dynamicTopoFvMesh::edgeLengthScale
(
    const label eIndex
) const
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
                FatalErrorIn
                (
                    "inline scalar dynamicTopoFvMesh::edgeLengthScale"
                    "(const label eIndex) const"
                )
                    << nl << "Face: " << eFaces[faceI]
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
                scale +=
                (
                    lengthEstimator().fixedLengthScale
                    (
                        eFaces[faceI],
                        edgePatch
                    )
                );
            }
        }

        scale *= 0.5;

        // If proximity-based refinement is requested,
        // test the proximity to the nearest non-neighbour.
        if (lengthEstimator().isProximityPatch(edgePatch))
        {
            label proximityFace = -1;

            // Perform a proximity-check.
            scalar distance = testProximity(eIndex, proximityFace);

            if (debug > 3 && self() == 0)
            {
                if
                (
                    (proximityFace > -1) &&
                    ((distance / 5.0) < scale)
                )
                {
                    Info << " Closest opposing face detected for edge: " << nl
                         << '\t' << eIndex
                         << " :: " << edges_[eIndex]
                         << " was face:\n"
                         << '\t' << proximityFace
                         << " :: " << polyMesh::faces()[proximityFace] << nl
                         << " with distance: " << distance
                         << endl;
                }
            }

            scale =
            (
                Foam::min
                (
                    scale,
                    ((distance / 3.0) - SMALL)/lengthEstimator().ratioMax()
                )
            );
        }

        // If curvature-based refinement is requested,
        // test the variation in face-normal directions.
        if (lengthEstimator().isCurvaturePatch(edgePatch))
        {
            // Obtain face-normals for both faces.
            label count = 0;
            FixedList<vector, 2> fNorm;

            forAll(eFaces, faceI)
            {
                if (neighbour_[eFaces[faceI]] == -1)
                {
                    // Obtain the normal.
                    fNorm[count] =
                    (
                        meshOps::faceNormal(faces_[eFaces[faceI]], points_)
                    );

                    // Normalize it.
                    fNorm[count] /= mag(fNorm[count]);

                    count++;
                }
            }

            scalar deviation = (fNorm[0] & fNorm[1]);
            scalar refDeviation = lengthEstimator().curvatureDeviation();

            if (mag(deviation) < refDeviation)
            {
                // Get the edge-length.
                scalar length = meshOps::edgeLength(edges_[eIndex], points_);

                if (debug > 3 && self() == 0)
                {
                    Info << "Deviation: " << deviation << nl
                         << "curvatureDeviation: " << refDeviation
                         << ", Edge: " << eIndex << ", Length: " << length
                         << ", Scale: " << scale << nl
                         << " Half-length: " << (0.5*length) << nl
                         << " MinRatio: "
                         << (lengthEstimator().ratioMin()*scale)
                         << endl;
                }

                scale =
                (
                    Foam::min
                    (
                        scale,
                        ((length - SMALL)/lengthEstimator().ratioMax())
                    )
                );
            }
        }

        // Limit scales if necessary
        lengthEstimator().limitScale(scale);
    }

    return scale;
}


// Check for edge bisection
inline bool dynamicTopoFvMesh::checkBisection
(
    const label index
) const
{
    scalar length = 0.0, scale = 0.0;

    if (twoDMesh_)
    {
        // Measure the boundary edge-length of the face in question
        length =
        (
            meshOps::edgeLength
            (
                edges_[getTriBoundaryEdge(index)],
                points_
            )
        );

        // Determine the length-scale at this face
        scale = faceLengthScale(index);
    }
    else
    {
        // Measure the edge-length
        length =
        (
            meshOps::edgeLength
            (
                edges_[index],
                points_
            )
        );

        // Determine the length-scale at this point in the mesh
        scale = edgeLengthScale(index);
    }

    if (length > lengthEstimator().ratioMax() * scale)
    {
        return true;
    }
    else
    {
        return false;
    }
}


// Check for edge collapse
inline bool dynamicTopoFvMesh::checkCollapse
(
    const label index
) const
{
    scalar length = 0.0, scale = 0.0;

    if (twoDMesh_)
    {
        // Measure the boundary edge-length of the face in question
        length =
        (
            meshOps::edgeLength
            (
                edges_[getTriBoundaryEdge(index)],
                points_
            )
        );

        // Determine the length-scale at this face
        scale = faceLengthScale(index);
    }
    else
    {
        // Measure the edge-length
        length =
        (
            meshOps::edgeLength
            (
                edges_[index],
                points_
            )
        );

        // Determine the length-scale at this point in the mesh
        scale = edgeLengthScale(index);
    }

    if (length < lengthEstimator().ratioMin() * scale)
    {
        return true;
    }
    else
    {
        return false;
    }
}


// Return the entity stack
inline stack& dynamicTopoFvMesh::Stack
(
    const label threadID
)
{
    return entityStack_[threadID];
}


// Return the integer threadID for a given pthread
// Return zero for single-threaded operation
inline label dynamicTopoFvMesh::self() const
{
    if (threader_->multiThreaded())
    {
        for (label i = 1; i <= threader_->getNumThreads(); i++)
        {
            if (pthread_equal(handlerPtr_[i].ID(), pthread_self()))
            {
                return i;
            }
        }
    }

    return 0;
}


// Initialize edge-stacks
inline void dynamicTopoFvMesh::initStacks()
{
    forAll(entityStack_, stackI)
    {
        entityStack_[stackI].clear();
    }

    // Prepare a filling sequence based on threading operation
    label tIndex = 0;
    labelList tID(threader().getNumThreads());

    if (threader_->multiThreaded())
    {
        forAll(tID, tI)
        {
            tID[tI] = (tI + 1);
        }
    }
    else
    {
        tID = 0;
    }

    if (twoDMesh_)
    {
        forAll(faces_, faceI)
        {
            // For coupled meshes, avoid certain faces.
            if (patchCoupling_.size() || procIndices_.size())
            {
                if (entitiesToAvoid_.found(faceI))
                {
                    continue;
                }
            }

            if (faces_[faceI].size() == 4)
            {
                Stack(tID[tIndex]).insert(faceI);

                tIndex = tID.fcIndex(tIndex);
            }
        }
    }
    else
    {
        forAll(edges_, edgeI)
        {
            // For coupled meshes, avoid certain edges.
            if (patchCoupling_.size() || procIndices_.size())
            {
                if (entitiesToAvoid_.found(edgeI))
                {
                    continue;
                }
            }

            if (edgeFaces_[edgeI].size())
            {
                Stack(tID[tIndex]).insert(edgeI);

                tIndex = tID.fcIndex(tIndex);
            }
        }
    }
}


// Initialize the coupled face-stack
inline void dynamicTopoFvMesh::initCoupledStack()
{
    // Clear existing lists/stacks.
    Stack(0).clear();

    const polyBoundaryMesh& boundary = boundaryMesh();

    // Loop though boundary faces and check whether
    // they belong to master/slave coupled patches.
    for (label faceI = nOldInternalFaces_; faceI < faces_.size(); faceI++)
    {
        // Add only valid faces
        if (faces_[faceI].empty())
        {
            continue;
        }

        label pIndex = whichPatch(faceI);

        if (pIndex == -1)
        {
            continue;
        }

        // Check if this is a locally coupled master face.
        if (patchCoupling_(pIndex))
        {
            // Add this to the coupled modification stack.
            if (twoDMesh_)
            {
                Stack(0).push(faceI);
            }
            else
            {
                const labelList& mfEdges = faceEdges_[faceI];

                forAll(mfEdges, edgeI)
                {
                    // Add this to the coupled modification stack.
                    Stack(0).push(mfEdges[edgeI]);
                }
            }
        }

        // Check if this is a processor patch.
        if (isA<processorPolyPatch>(boundary[pIndex]))
        {
            // Check if this is a master processor patch.
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[pIndex])
            );

            label neiProcID = pp.neighbProcNo();

            if (neiProcID > Pstream::myProcNo())
            {
                // Add this to the coupled modification stack.
                if (twoDMesh_)
                {
                    Stack(0).push(faceI);
                }
                else
                {
                    const labelList& mfEdges = faceEdges_[faceI];

                    forAll(mfEdges, edgeI)
                    {
                        Stack(0).push(mfEdges[edgeI]);
                    }
                }
            }
        }
    }
}


// Build a list of entities that need to be avoided
// by regular topo-changes.
inline void dynamicTopoFvMesh::buildEntitiesToAvoid()
{
    entitiesToAvoid_.clear();

    // Build a set of entities to avoid during regular modifications,
    // and build a master stack for coupled modifications.
    const polyBoundaryMesh& boundary = boundaryMesh();

    // Determine locally coupled slave patches.
    labelHashSet localSlavePatches;

    forAll(patchCoupling_, patchI)
    {
        if (patchCoupling_(patchI))
        {
            const coupleMap& cMap = patchCoupling_[patchI].patchMap();

            localSlavePatches.insert(cMap.slaveIndex());
        }
    }

    // Loop through boundary faces and check whether
    // they belong to master/slave coupled patches.
    for (label faceI = nOldInternalFaces_; faceI < faces_.size(); faceI++)
    {
        // Add only valid faces
        if (faces_[faceI].empty())
        {
            continue;
        }

        label pIndex = whichPatch(faceI);

        if (pIndex == -1)
        {
            continue;
        }

        // Check if this is a coupled face.
        if
        (
            patchCoupling_(pIndex) ||
            localSlavePatches.found(pIndex) ||
            isA<processorPolyPatch>(boundary[pIndex])
        )
        {
            if (twoDMesh_)
            {
                // Avoid this face during regular modification.
                entitiesToAvoid_.insert(faceI, empty());
            }
            else
            {
                const labelList& fEdges = faceEdges_[faceI];

                forAll(fEdges, edgeI)
                {
                    // Avoid this edge during regular modification.
                    entitiesToAvoid_.insert(fEdges[edgeI], empty());
                }
            }
        }
    }
}


// Check whether the specified edge is a coupled master edge.
inline bool dynamicTopoFvMesh::isCoupledMaster
(
    const label eIndex
) const
{
    if (!coupledModification_)
    {
        return true;
    }

    return locallyCoupledEntity(eIndex);
}


// Set coupled modification
inline void dynamicTopoFvMesh::setCoupledModification() const
{
    coupledModification_ = true;
}


// Unset coupled modification
inline void dynamicTopoFvMesh::unsetCoupledModification() const
{
    coupledModification_ = false;
}


// Method to determine the old boundary patch index for a given face
// Similar to the polyBoundaryMesh routine, but works on local information
inline label dynamicTopoFvMesh::whichPatch
(
    const label index
) const
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
            "label dynamicTopoFvMesh::whichPatch"
            "(const label index) const"
        )
            << "Cannot find patch information for face index: "
            << index << nl
            << " It appears that face ordering is"
            << " inconsistent with patch information." << nl
            << " Mesh info: " << nl
            << "    nOldInternalFaces: " << nOldInternalFaces_ << nl
            << "    oldPatchStarts: " << oldPatchStarts_ << nl
            << "    oldPatchSizes: " << oldPatchSizes_ << nl
            << abort(FatalError);
    }

    return -2;
}


// Method to determine the old boundary patch index for a given edge
inline label dynamicTopoFvMesh::whichEdgePatch
(
    const label index
) const
{
    if (index < nOldInternalEdges_)
    {
        return -1;
    }

    for(label i = 0; i<numPatches_; i++)
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
            "label dynamicTopoFvMesh::whichEdgePatch"
            "(const label index) const"
        )
            << "Cannot find patch information for edge index: "
            << index << nl
            << " It appears that edge ordering is"
            << " inconsistent with patch information." << nl
            << " Mesh info: " << nl
            << "    nOldInternalEdges: " << nOldInternalEdges_ << nl
            << "    oldEdgePatchStarts: " << oldEdgePatchStarts_ << nl
            << "    oldEdgePatchSizes: " << oldEdgePatchSizes_ << nl
            << abort(FatalError);
    }

    return -2;
}


// Set a particular face index as flipped.
inline void dynamicTopoFvMesh::setFlip(const label fIndex)
{
    if (fIndex < nOldFaces_)
    {
        if (flipFaces_.found(fIndex))
        {
            flipFaces_.erase(fIndex);
        }
        else
        {
            flipFaces_.insert(fIndex);
        }
    }
}


// Utility method to compute the minimum quality of a vertex hull
inline scalar dynamicTopoFvMesh::computeMinQuality
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
            cQuality = tetMetric_(a, b, c, d);

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

            // Compute the quality
            cQuality = tetMetric_(a, b, c, d);

            // Check if the quality is worse
            minQuality = Foam::min(cQuality, minQuality);
        }
    }

    if (coupledModification_)
    {
        if (locallyCoupledEntity(eIndex))
        {
            // Compute the minimum quality of the slave edge as well.
            label sIndex = -1;

            // Determine the slave index.
            forAll(patchCoupling_, patchI)
            {
                if (patchCoupling_(patchI))
                {
                    const label edgeEnum  = coupleMap::EDGE;
                    const coupleMap& cMap = patchCoupling_[patchI].patchMap();

                    if ((sIndex = cMap.findSlaveIndex(edgeEnum, eIndex)) > -1)
                    {
                        break;
                    }
                }
            }

            if (sIndex == -1)
            {
                FatalErrorIn
                (
                    "inline scalar dynamicTopoFvMesh::computeMinQuality"
                    "(const label eIndex) const"
                )
                    << nl << "Coupled maps were improperly specified." << nl
                    << " Slave index not found for: " << nl
                    << " Edge: " << eIndex << nl
                    << abort(FatalError);
            }

            // Temporarily turn off coupledModification
            unsetCoupledModification();

            scalar slaveQuality = computeMinQuality(sIndex);

            minQuality = Foam::min(slaveQuality, minQuality);

            // Turn it back on.
            setCoupledModification();
        }
        else
        if (processorCoupledEntity(eIndex))
        {
            // Compute the minimum quality across patchSubMeshes.

        }
    }

    // Ensure that the mesh is valid
    if (minQuality < 0.0)
    {
        // if (debug > 3)
        {
            // Write out faces and cells for post processing.
            labelHashSet iFaces, iCells, bFaces;

            const labelList& eFaces = edgeFaces_[eIndex];

            forAll(eFaces, faceI)
            {
                iFaces.insert(eFaces[faceI]);

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
            writeVTK(Foam::name(eIndex) + "_iFaces", iFaces.toc(), 2);

            // Write out the boundary patches (for post-processing reference)
            for
            (
                label faceI = nOldInternalFaces_;
                faceI < faces_.size();
                faceI++
            )
            {
                if (faces_[faceI].empty())
                {
                    continue;
                }

                label pIndex = whichPatch(faceI);

                if (pIndex != -1)
                {
                    bFaces.insert(faceI);
                }
            }

            writeVTK(Foam::name(eIndex) + "_bFaces", bFaces.toc(), 2);
        }

        FatalErrorIn
        (
            "inline scalar dynamicTopoFvMesh::computeMinQuality"
            "(const label eIndex) const"
        )
            << "Encountered negative cell-quality!" << nl
            << "Edge: " << eIndex << ": " << edgeToCheck << nl
            << "EdgePoints: " << hullVertices << nl
            << "Minimum Quality: " << minQuality
            << abort(FatalError);
    }

    return minQuality;
}


// Method to determine whether the master face is locally coupled
inline bool dynamicTopoFvMesh::locallyCoupledEntity
(
    const label index,
    bool checkSlaves
) const
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    // Bail out if no patchCoupling is present
    if (!patchCoupling_.size())
    {
        return false;
    }

    if (twoDMesh_)
    {
        label patch = whichPatch(index);

        if (patch == -1)
        {
            return false;
        }

        // Processor checks receive priority.
        if (isA<processorPolyPatch>(boundary[patch]))
        {
            return false;
        }

        // Check coupled master patches.
        if (patchCoupling_(patch))
        {
            return true;
        }
        else
        if (checkSlaves)
        {
            // Check on slave patches as well.
            forAll(patchCoupling_, pI)
            {
                if (patchCoupling_(pI))
                {
                    const coupleMap& cMap = patchCoupling_[pI].patchMap();

                    if (cMap.slaveIndex() == patch)
                    {
                        return true;
                    }
                }
            }
        }
    }
    else
    {
        const labelList& eFaces = edgeFaces_[index];

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
                if (patchCoupling_(patch))
                {
                    return true;
                }

                if (checkSlaves)
                {
                    // Check on slave patches as well.
                    forAll(patchCoupling_, pI)
                    {
                        if (patchCoupling_(pI))
                        {
                            const coupleMap& cMap =
                            (
                                patchCoupling_[pI].patchMap()
                            );

                            if (cMap.slaveIndex() == patch)
                            {
                                return true;
                            }
                        }
                    }
                }
            }
        }
    }

    // Could not find any faces on locally coupled patches.
    return false;
}


// Method to determine the locally coupled patch index
inline label dynamicTopoFvMesh::locallyCoupledEdgePatch
(
    const label eIndex
) const
{
    const labelList& eFaces = edgeFaces_[eIndex];

    // Search for boundary faces, and determine boundary type.
    forAll(eFaces, faceI)
    {
        if (neighbour_[eFaces[faceI]] == -1)
        {
            label patch = whichPatch(eFaces[faceI]);

            // Check coupled master patches.
            if (patchCoupling_(patch))
            {
                return patch;
            }

            // Check on slave patches as well.
            forAll(patchCoupling_, pI)
            {
                if (patchCoupling_(pI))
                {
                    const coupleMap& cMap = patchCoupling_[pI].patchMap();

                    if (cMap.slaveIndex() == patch)
                    {
                        return patch;
                    }
                }
            }
        }
    }

    // Could not find any faces on locally coupled patches.
    FatalErrorIn
    (
        "inline label dynamicTopoFvMesh::locallyCoupledEdgePatch"
        "(const label cIndex) const"
    )
        << "Edge: " << eIndex << ":: " << edges_[eIndex]
        << " does not lie on any coupled patches."
        << abort(FatalError);

    return -1;
}


// Method to determine whether the master entity is on a processor
inline bool dynamicTopoFvMesh::processorCoupledEntity
(
    const label index
) const
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    label patch = -2;

    if (twoDMesh_)
    {
        patch = whichPatch(index);

        if (patch == -1)
        {
            return false;
        }

        if (isA<processorPolyPatch>(boundary[patch]))
        {
            return true;
        }
    }
    else
    {
        const labelList& eFaces = edgeFaces_[index];

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
    }

    // Could not find any faces on processor patches.
    return false;
}


// Obtain the maximum coupled patch index
inline label dynamicTopoFvMesh::getMaxCouplingIndex() const
{
    // Allocate a size for coupled patches
    label index = -1;

    forAll(patchCoupling_, patchI)
    {
        if (patchCoupling_(patchI))
        {
            index = index > patchI ? index : patchI;
        }
    }

    return index;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
