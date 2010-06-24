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

#include "triPointRef.H"
#include "tetPointRef.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// Return a reference to the multiThreader
inline const multiThreader& dynamicTopoFvMesh::threader() const
{
    return threader_();
}


// Does the mesh perform edge refinement?
inline bool dynamicTopoFvMesh::edgeRefinement() const
{
    return edgeRefinement_;
}


// Return the magnitude of an edge
inline scalar dynamicTopoFvMesh::edgeLength
(
    const edge& edgeToCheck,
    const UList<point>& points
) const
{
    return mag(points[edgeToCheck.end()] - points[edgeToCheck.start()]);
}


// Compute the centroid for a given face
inline vector dynamicTopoFvMesh::faceCentre
(
    const face& faceToCheck,
    const UList<point>& points
) const
{
    face localFace(faceToCheck.size());
    pointField localPoints(faceToCheck.size(), vector::zero);

    // Fill in local points
    forAll(localPoints, pointI)
    {
        localFace[pointI] = pointI;
        localPoints[pointI] = points[faceToCheck[pointI]];
    }

    return localFace.centre(localPoints);
}


// Compute the normal for a given face
inline vector dynamicTopoFvMesh::faceNormal
(
    const face& faceToCheck,
    const UList<point>& points
) const
{
    face localFace(faceToCheck.size());
    pointField localPoints(faceToCheck.size(), vector::zero);

    // Fill in local points
    forAll(localPoints, pointI)
    {
        localFace[pointI] = pointI;
        localPoints[pointI] = points[faceToCheck[pointI]];
    }

    return localFace.normal(localPoints);
}


// Given a face and cell index, find the apex point for a tet cell.
inline label dynamicTopoFvMesh::tetApexPoint
(
    const label cIndex,
    const label fIndex
) const
{
    label apexPoint = -1;
    bool foundApex = false;

    const cell& cellToCheck = cells_[cIndex];
    const face& baseFace = faces_[fIndex];

    forAll(cellToCheck, faceI)
    {
        const face& faceToCheck = faces_[cellToCheck[faceI]];

        apexPoint = findIsolatedPoint(baseFace, faceToCheck);

        if (foundApex > -1)
        {
            foundApex = true;
            break;
        }
    }

    if (!foundApex)
    {
        Info << "fIndex: " << fIndex << ":: " << faces_[fIndex] << endl;
        Info << "cIndex: " << cIndex << ":: " << cellToCheck << endl;

        forAll(cellToCheck, faceI)
        {
            Info << '\t' << cellToCheck[faceI] << ":: "
                 << faces_[cellToCheck[faceI]] << endl;
        }

        FatalErrorIn("dynamicTopoFvMesh::tetApexPoint()")
            << "Could not find an apex point in cell: " << cIndex
            << " given face: " << fIndex
            << abort(FatalError);
    }

    return apexPoint;
}


// Given a cell index, find the centroid / volume of a cell.
inline void dynamicTopoFvMesh::cellCentreAndVolume
(
    const label cIndex,
    const UList<point>& points,
    vector& centre,
    scalar& volume
) const
{
    // Reset inputs
    volume = 0.0;
    centre = vector::zero;

    const cell& cellToCheck = cells_[cIndex];

    // Average face-centres to get an estimate centroid
    vector cEst = vector::zero;

    forAll(cellToCheck, faceI)
    {
        cEst += faceCentre(faces_[cellToCheck[faceI]], points);
    }

    cEst /= cellToCheck.size();

    forAll(cellToCheck, faceI)
    {
        vector fArea = faceNormal(faces_[cellToCheck[faceI]], points);
        vector fCentre = faceCentre(faces_[cellToCheck[faceI]], points);

        if (owner_[cellToCheck[faceI]] != cIndex)
        {
            fArea *= -1.0;
        }

        // Calculate 3*face-pyramid volume
        scalar pyr3Vol = fArea & (fCentre - cEst);

        // Calculate face-pyramid centre
        vector pc = (3.0/4.0)*fCentre + (1.0/4.0)*cEst;

        // Accumulate volume-weighted face-pyramid centre
        centre += pyr3Vol*pc;

        // Accumulate face-pyramid volume
        volume += pyr3Vol;
    }

    centre /= volume + VSMALL;
    volume *= (1.0/3.0);
}


// Utility member function for separating axis tests
inline label dynamicTopoFvMesh::whichSide
(
    const labelList& cellPoints,
    const UList<point>& points,
    const point& dir,
    const point& p
) const
{
    label nP = 0, nN = 0;

    forAll(cellPoints, pointI)
    {
        scalar t = (dir & (points[cellPoints[pointI]] - p));

        if (t > 0.0)
        {
            nP++;
        }
        else
        if (t < 0.0)
        {
            nN++;
        }

        if (nP && nN)
        {
            return 0;
        }
    }

    return (nP ? 1 : -1);
}


// Determine whether a segment intersects a face
inline bool dynamicTopoFvMesh::segmentFaceIntersection
(
    const edge& edgeToCheck,
    const face& faceToCheck,
    const UList<point>& points,
    vector& intersectionPoint
) const
{
    // Specify unit face-normals
    vector n = faceNormal(faceToCheck, points);
    n /= mag(n) + VSMALL;

    point p1 = points[edgeToCheck[0]];
    point p2 = points[edgeToCheck[1]];
    point p3 = points[faceToCheck[0]];

    // Compute uValue
    scalar numerator = n & (p3 - p1);
    scalar denominator = n & (p2 - p1);

    // Check if the edge is parallel to the face
    scalar tolerance = (1e-4 * mag(p2 - p1));

    if (mag(denominator) < tolerance)
    {
        return false;
    }

    scalar u = (numerator / denominator);

    // Check for intersection along line.
    if ((u > tolerance) && (u < (1.0 - tolerance)))
    {
        // Compute point of intersection
        intersectionPoint = p1 + u*(p2 - p1);

        // Also make sure that intPoint lies within the face
        if (pointInFace(faceToCheck, points, intersectionPoint))
        {
            return true;
        }
    }

    return false;
}


// Determine whether the particular point lies
// inside the given face
inline bool dynamicTopoFvMesh::pointInFace
(
    const face& faceToCheck,
    const UList<point>& points,
    const point& checkPoint
) const
{
    vector nf = faceNormal(faceToCheck, points);
    nf /= mag(nf) + VSMALL;

    forAll(faceToCheck, pI)
    {
        vector p, q;

        label nI = faceToCheck.fcIndex(pI);

        // Identify vectors
        p = (points[faceToCheck[nI]] - points[faceToCheck[pI]]);
        q = (checkPoint - points[faceToCheck[nI]]);

        scalar area = 0.5 * ((p ^ q) & nf);

        if (area < 0.0)
        {
            return false;
        }
    }

    // Passed test with all edges
    return true;
}


// Determine whether the particular point lies
// inside the given cell
inline bool dynamicTopoFvMesh::pointInCell
(
    const label cIndex,
    const cell& cellToCheck,
    const UList<face>& faces,
    const UList<label>& owner,
    const UList<point>& points,
    const point& checkPoint
) const
{
    vector xf(vector::zero), nf(vector::zero);

    forAll(cellToCheck, faceI)
    {
        xf = faceCentre(faces[cellToCheck[faceI]], points);
        nf = faceNormal(faces[cellToCheck[faceI]], points);

        if (((xf - checkPoint) & nf) > 0.0)
        {
            if (owner[cellToCheck[faceI]] != cIndex)
            {
                return false;
            }
        }
        else
        {
            if (owner[cellToCheck[faceI]] == cIndex)
            {
                return false;
            }
        }
    }

    // Passed test with all faces
    return true;
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
    FatalErrorIn("dynamicTopoFvMesh::getEdgeIndex")
        << "Could not find an appropriate edge index for edge:"
        << edgeToCheck
        << abort(FatalError);

    return -1;
}


// Return the appropriate length-scale for boundary face
inline scalar dynamicTopoFvMesh::boundaryLengthScale
(
    const label fIndex
) const
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
        }
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
            FatalErrorIn("dynamicTopoFvMesh::faceLengthScale()")
                << "Face: " << fIndex
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
        scale = boundaryLengthScale(fIndex);

        // If proximity-based refinement is requested,
        // test the proximity to the nearest non-neighbour.
        if (proximityPatches_.found(boundaryMesh()[facePatch].name()))
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

            scalar tol = 1e-6;

            scale = Foam::min(scale, ((distance / 3.0) - tol)/ratioMax_);
        }

        // If curvature-based refinement is requested,
        // test the variation in face-normal directions.
        if (curvaturePatches_.found(boundaryMesh()[facePatch].name()))
        {
            notImplemented("curvatureRefinement2D");
        }

        // Check if we're limited by min/max lengthScales
        scale = Foam::max(scale, minLengthScale_);
        scale = Foam::min(scale, maxLengthScale_);
    }

    if (coupledModification_)
    {
        // Check if this is a coupled processor face.
        // - Locally coupled faces use length-scale information only
        //   from the master patch.
        if (processorCoupledEdge(fIndex))
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

        scale *= 0.5;

        // If proximity-based refinement is requested,
        // test the proximity to the nearest non-neighbour.
        if (proximityPatches_.found(boundaryMesh()[edgePatch].name()))
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

            scalar tol = 1e-6;

            scale = Foam::min(scale, ((distance / 3.0) - tol)/ratioMax_);
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
                    fNorm[count] = faceNormal(faces_[eFaces[faceI]], points_);

                    // Normalize it.
                    fNorm[count] /= mag(fNorm[count]);

                    count++;
                }
            }

            scalar deviation = (fNorm[0] & fNorm[1]);

            if (mag(deviation) < curvatureDeviation_)
            {
                // Get the edge-length.
                scalar length = edgeLength(edges_[eIndex], points_);

                if (debug > 3 && self() == 0)
                {
                    Info << "Deviation: " << deviation << nl
                         << "curvatureDeviation: " << curvatureDeviation_
                         << ", Edge: " << eIndex << ", Length: " << length
                         << ", Scale: " << scale << nl
                         << " Half-length: " << (0.5*length) << nl
                         << " MinRatio: " << (ratioMin_*scale)
                         << endl;
                }

                scalar tol = 1e-6;

                scale = Foam::min(scale, ((length - tol)/ratioMax_));
            }
        }

        // Check if we're limited by min/max lengthScales
        scale = Foam::max(scale, minLengthScale_);
        scale = Foam::min(scale, maxLengthScale_);
    }

    return scale;
}


// Check for face bisection
inline bool dynamicTopoFvMesh::checkFaceBisection
(
    const label fIndex
) const
{
    // Measure the boundary edge-length of the face in question
    scalar length = edgeLength(edges_[getTriBoundaryEdge(fIndex)], points_);

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
) const
{
    // Measure the boundary edge-length of the face in question
    scalar length = edgeLength(edges_[getTriBoundaryEdge(fIndex)], points_);

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
) const
{
    // Measure the edge-length
    scalar length = edgeLength(edges_[eIndex], points_);

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
) const
{
    // Measure the edge-length
    scalar length = edgeLength(edges_[eIndex], points_);

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


// Return the face-stack
inline stack& dynamicTopoFvMesh::faceStack
(
    const label threadID
)
{
    return faceStack_[threadID];
}


// Return the edge-stack for a particular thread
inline stack& dynamicTopoFvMesh::edgeStack
(
    const label threadID
)
{
    return edgeStack_[threadID];
}


// Synchronize all processes
inline void dynamicTopoFvMesh::synchronize()
{
    bool reached = true;
    reduce(reached, andOp<bool>());
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
                    (findIndex(entitiesToAvoid_, edgeI) == -1) &&
                    (edgeFaces_[edgeI].size())
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
                    (findIndex(entitiesToAvoid_, edgeI) == -1) &&
                    (edgeFaces_[edgeI].size())
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


// Initialize the coupled face-stack
inline void dynamicTopoFvMesh::initCoupledFaceStack()
{
    // Clear existing lists/stacks.
    faceStack(0).clear();

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
            faceStack(0).push(faceI);
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
                faceStack(0).push(faceI);
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
        // Add only valid faces
        if (!faces_[faceI].size())
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
                if (findIndex(entitiesToAvoid_, faceI) == -1)
                {
                    // Avoid this face during regular modification.
                    entitiesToAvoid_.append(faceI);
                }
            }
            else
            {
                const labelList& fEdges = faceEdges_[faceI];

                forAll(fEdges, edgeI)
                {
                    if (findIndex(entitiesToAvoid_, fEdges[edgeI]) == -1)
                    {
                        // Avoid this edge during regular modification.
                        entitiesToAvoid_.append(fEdges[edgeI]);
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
) const
{
    if (!coupledModification_)
    {
        return true;
    }

    return locallyCoupledEdge(eIndex);
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

        // For coupled meshes, avoid certain edges.
        if (patchCoupling_.size() || procIndices_.size())
        {
            forAll(faces_, faceI)
            {
                if
                (
                    (findIndex(entitiesToAvoid_, faceI) == -1) &&
                    (faces_[faceI].size() == 4)
                )
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
                if (faces_[faceI].size() == 4)
                {
                    faceStack_[tID].insert(faceI);

                    tID = (tID == faceStack_.size() - 1) ? 1 : tID + 1;
                }
            }
        }
    }
    else
    {
        // For coupled meshes, avoid certain edges.
        if (patchCoupling_.size() || procIndices_.size())
        {
            forAll(faces_, faceI)
            {
                if
                (
                    (findIndex(entitiesToAvoid_, faceI) == -1) &&
                    (faces_[faceI].size() == 4)
                )
                {
                    faceStack_[0].insert(faceI);
                }
            }
        }
        else
        {
            forAll(faces_, faceI)
            {
                if (faces_[faceI].size() == 4)
                {
                    faceStack_[0].insert(faceI);
                }
            }
        }
    }
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
            "label dynamicTopoFvMesh::whichEdgePatch(const label index) const"
        )
            << "Cannot find patch information for edge index: " << index << nl
            << " It appears that edge ordering is"
            << " inconsistent with patch information."
            << abort(FatalError);
    }

    return -2;
}


// Utility method to find the isolated point given two triangular faces.
//  - Returns the point on checkFace that does not belong to baseFace.
inline label dynamicTopoFvMesh::findIsolatedPoint
(
    const face& baseFace,
    const face& checkFace
) const
{
    // Get the fourth point
    forAll(checkFace, pointI)
    {
        if
        (
            checkFace[pointI] != baseFace[0] &&
            checkFace[pointI] != baseFace[1] &&
            checkFace[pointI] != baseFace[2]
        )
        {
            return checkFace[pointI];
        }
    }

    return -1;
}


// Utility method to find the isolated point on a triangular face
// that doesn't lie on the specified edge. Also returns the point next to it.
inline void dynamicTopoFvMesh::findIsolatedPoint
(
    const face& f,
    const edge& e,
    label& ptIndex,
    label& nextPtIndex
) const
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
) const
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


// Method to insert labels in a face, so that
// right-handedness is preserved.
inline void dynamicTopoFvMesh::insertPointLabels
(
    const vector& refNorm,
    const vectorField& points,
    const labelHashSet& pLabels,
    face& modFace
) const
{
    // Need to configure a new face.
    face newFace(modFace);

    forAllConstIter(labelHashSet, pLabels, pIter)
    {
        forAll(newFace, pI)
        {
            label nI = newFace.fcIndex(pI);

            // Compute the normal.
            vector newNorm =
            (
                triPointRef
                (
                    points[newFace[pI]],
                    points[pIter.key()],
                    points[newFace[nI]]
                ).normal()
            );

            if ((refNorm & newNorm) > 0.0)
            {
                // Insert the point.
                insertLabel
                (
                    pIter.key(),
                    newFace[pI],
                    newFace[nI],
                    newFace
                );

                break;
            }
        }
    }

    // Take over storage
    modFace.transfer(newFace);
}


// Method to insert a label between two labels in a list
// Assumes that all labels are unique.
inline void dynamicTopoFvMesh::insertLabel
(
    const label newLabel,
    const label labelA,
    const label labelB,
    labelList& list
) const
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
template <class Type>
inline void dynamicTopoFvMesh::sizeUpList
(
    const Type item,
    List<Type>& list
) const
{
    // Create a new list
    List<Type> newList(list.size() + 1);

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
template <class Type>
inline void dynamicTopoFvMesh::sizeDownList
(
    const Type item,
    List<Type>& list
) const
{
    label index = -1;

    if ((index = findIndex(list, item)) > -1)
    {
        removeIndex(index, list);
    }
    else
    {
        FatalErrorIn("dynamicTopoFvMesh::sizeDownList")
            << "Item: " << item
            << " was not found in list. " << nl
            << " List: " << nl << list
            << abort(FatalError);
    }
}


// Remove an item at a particular index in the list
template <class Type>
inline void dynamicTopoFvMesh::removeIndex
(
    const label index,
    List<Type>& list
) const
{
    // Create a new list
    List<Type> newList(list.size() - 1);

    // Copy individual items
    label n = 0;

    forAll(list, itemI)
    {
        if (itemI == index)
        {
            continue;
        }

        newList[n++] = list[itemI];
    }

    // Overwrite
    list.transfer(newList);
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
        if (locallyCoupledEdge(eIndex))
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
                FatalErrorIn("dynamicTopoFvMesh::computeMinQuality")
                    << "Coupled maps were improperly specified." << nl
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
        if (processorCoupledEdge(eIndex))
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

        FatalErrorIn("dynamicTopoFvMesh::computeMinQuality()")
            << "Encountered negative cell-quality!" << nl
            << "Edge: " << eIndex << ": " << edgeToCheck << nl
            << "EdgePoints: " << hullVertices << nl
            << "Minimum Quality: " << minQuality
            << abort(FatalError);
    }

    return minQuality;
}


// Utility method to compute the quality of a tetrahedral cell
inline scalar dynamicTopoFvMesh::tetQuality
(
    const label cIndex
) const
{
    const cell& cellToCheck = cells_[cIndex];

    const face& baseFace = faces_[cellToCheck[0]];
    const face& checkFace = faces_[cellToCheck[1]];

    // Get the fourth point
    label apexPoint = findIsolatedPoint(baseFace, checkFace);

    // Something's wrong with connectivity.
    if (apexPoint == -1)
    {
        FatalErrorIn("dynamicTopoFvMesh::tetQuality()")
            << "Cell: " << cIndex
            << " has inconsistent connectivity."
            << abort(FatalError);
    }

    // Compute cell-quality and write-out
    if (owner_[cellToCheck[0]] == cIndex)
    {
        return
        (
            tetMetric_
            (
                points_[baseFace[2]],
                points_[baseFace[1]],
                points_[baseFace[0]],
                points_[apexPoint]
            )
        );
    }
    else
    {
        return
        (
            tetMetric_
            (
                points_[baseFace[0]],
                points_[baseFace[1]],
                points_[baseFace[2]],
                points_[apexPoint]
            )
        );
    }

    // This should never happen.
    FatalErrorIn("dynamicTopoFvMesh::tetQuality()")
        << nl << " Illegal connectivity."
        << abort(FatalError);

    return 0.0;
}


// Check whether the given face should not be bisected/collapsed
inline bool dynamicTopoFvMesh::checkFaceModification
(
    const label fIndex
) const
{
    // Internal faces don't count
    label facePatch = -1;

    if ((facePatch = whichPatch(fIndex)) < 0)
    {
        return false;
    }
    else
    {
        // Check whether this face shouldn't be modified
        if (findIndex(noModPatchIDs_, facePatch) > -1)
        {
            return true;
        }
    }

    return false;
}


// Check whether the given edge should not be bisected/collapsed
inline bool dynamicTopoFvMesh::checkEdgeModification
(
    const label eIndex
) const
{
    // Internal edges don't count
    label edgePatch = -1;

    if ((edgePatch = whichEdgePatch(eIndex)) < 0)
    {
        return false;
    }
    else
    {
        // Check whether this edge shouldn't be modified
        if (findIndex(noModPatchIDs_, edgePatch) > -1)
        {
            return true;
        }
    }

    return false;
}


// Method to determine whether the master face is locally coupled
inline bool dynamicTopoFvMesh::locallyCoupledFace
(
    const label fIndex,
    bool checkSlaves
) const
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    label patch = whichPatch(fIndex);

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

    return false;
}


// Method to determine whether the master edge is locally coupled.
inline bool dynamicTopoFvMesh::locallyCoupledEdge
(
    const label eIndex,
    bool checkSlaves
) const
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    const labelList& eFaces = edgeFaces_[eIndex];

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
                        const coupleMap& cMap = patchCoupling_[pI].patchMap();

                        if (cMap.slaveIndex() == patch)
                        {
                            return true;
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
    FatalErrorIn("dynamicTopoFvMesh::locallyCoupledEdgeIndex()")
        << "Edge: " << eIndex << ":: " << edges_[eIndex]
        << " does not lie on any coupled patches."
        << abort(FatalError);

    return -1;
}


// Method to determine whether the master face is on a processor
inline bool dynamicTopoFvMesh::processorCoupledFace
(
    const label fIndex
) const
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    label patch = whichPatch(fIndex);

    if (patch == -1)
    {
        return false;
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
) const
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    const labelList& eFaces = edgeFaces_[eIndex];

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
