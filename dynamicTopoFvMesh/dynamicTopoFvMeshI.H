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

// Obtain the tangent-to-edge vector
inline vector dynamicTopoFvMesh::edgeCentre
(
    const label eIndex
) const
{
    const edge& thisEdge = edges_[eIndex];

    return 0.5 * (points_[thisEdge[0]] + points_[thisEdge[1]]);
}

// Obtain the tangent-to-edge vector
inline vector dynamicTopoFvMesh::tangentToEdge
(
    const label eIndex
) const
{
    const edge& thisEdge = edges_[eIndex];

    return (points_[thisEdge[1]] - points_[thisEdge[0]]);
}

// Find the circumcentre of a face
inline point dynamicTopoFvMesh::circumCentre
(
    const label fIndex
) const
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

// Find the centre of a triangle face
inline point dynamicTopoFvMesh::quadFaceCentre
(
    const label fIndex
) const
{
    const face& thisFace = faces_[fIndex];

    return 0.25 *
    (
        points_[thisFace[0]]
      + points_[thisFace[1]]
      + points_[thisFace[2]]
      + points_[thisFace[3]]
    );
}

// Find the centre of a triangle face
inline point dynamicTopoFvMesh::triFaceCentre
(
    const label fIndex
) const
{
    const face& thisFace = faces_[fIndex];

    return
    (
        points_[thisFace[0]]
      + points_[thisFace[1]]
      + points_[thisFace[2]]
    )/3.0;
}

// Find the area of a quad face.
inline scalar dynamicTopoFvMesh::quadFaceArea
(
    const face& quadFace
) const
{
    return Foam::mag(quadFaceNormal(quadFace));
}

// Find the area of a triangle face.
// This function also assumes face right-handedness
inline scalar dynamicTopoFvMesh::triFaceArea
(
    const face& triFace
) const
{
    return Foam::mag(triFaceNormal(triFace));
}

// Find the normal of a quad face.
// This function also assumes face right-handedness
inline vector dynamicTopoFvMesh::quadFaceNormal
(
    const face& quadFace
) const
{
    vector n = vector::zero, v = vector::zero, w = vector::zero;

    point centrePoint = 0.25 *
    (
        points_[quadFace[0]]
      + points_[quadFace[1]]
      + points_[quadFace[2]]
      + points_[quadFace[3]]
    );

    // Add normals from four sub-triangles.
    v = points_[quadFace[1]] - points_[quadFace[0]];
    w = centrePoint - points_[quadFace[1]];

    n += 0.5 * (v ^ w);

    v = points_[quadFace[2]] - points_[quadFace[1]];
    w = centrePoint - points_[quadFace[2]];

    n += 0.5 * (v ^ w);

    v = points_[quadFace[3]] - points_[quadFace[2]];
    w = centrePoint - points_[quadFace[3]];

    n += 0.5 * (v ^ w);

    v = points_[quadFace[0]] - points_[quadFace[3]];
    w = centrePoint - points_[quadFace[0]];

    n += 0.5 * (v ^ w);

    return n;
}

// Find the normal of a triangle face.
// This function also assumes face right-handedness
inline vector dynamicTopoFvMesh::triFaceNormal
(
    const face& triFace
) const
{
    vector v = points_[triFace[1]] - points_[triFace[0]];
    vector w = points_[triFace[2]] - points_[triFace[0]];

    return 0.5 * (v ^ w);
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
    const face& checkFace = faces_[fIndex];

    forAll(cellToCheck, faceI)
    {
        const face& faceToCheck = faces_[cellToCheck[faceI]];

        forAll(faceToCheck, pointI)
        {
            if
            (
                faceToCheck[pointI] != checkFace[0] &&
                faceToCheck[pointI] != checkFace[1] &&
                faceToCheck[pointI] != checkFace[2]
            )
            {
                apexPoint = faceToCheck[pointI];

                foundApex = true;
                break;
            }
        }

        if (foundApex)
        {
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

// Given a cell index, find the centroid of a prism cell.
inline vector dynamicTopoFvMesh::prismCellCentre
(
    const label cIndex
) const
{
    vector centre = vector::zero;

    const cell& cellToCheck = cells_[cIndex];

    forAll(cellToCheck, faceI)
    {
        const face& faceToCheck = faces_[cellToCheck[faceI]];

        if (faceToCheck.size() == 3)
        {
            centre += triFaceCentre(cellToCheck[faceI]);
        }

        if (faceToCheck.size() == 4)
        {
            centre += quadFaceCentre(cellToCheck[faceI]);
        }
    }

    centre /= cellToCheck.size();

    return centre;
}

// Given a cell index, find the centroid of a tet cell.
inline vector dynamicTopoFvMesh::tetCellCentre
(
    const label cIndex
) const
{
    vector centre = vector::zero;

    const cell& cellToCheck = cells_[cIndex];

    forAll(cellToCheck, faceI)
    {
        centre += triFaceCentre(cellToCheck[faceI]);
    }

    centre /= cellToCheck.size();

    return centre;
}

// Compare two triangular faces.
// Identical to triFace definition.
inline label dynamicTopoFvMesh::compare
(
    const face& a,
    const face& b
) const
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
) const
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
    const scalar tolerance,
    scalar& vol
) const
{
    vol = (1.0/6.0)*(((b - a) ^ (c - a)) & (d - a));

    if (mag(vol) < tolerance)
    {
        return 0;
    }

    if (Foam::sign(vol) < 0.0)
    {
        return -1;
    }

    return 1;
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
) const
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
) const
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
) const
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
) const
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
inline bool dynamicTopoFvMesh::checkQuadFace
(
    const label fIndex
) const
{
    return (faces_[fIndex].size() == 4);
}

// Compute the length of an edge
inline scalar dynamicTopoFvMesh::edgeLength
(
    const label eIndex
) const
{
    const edge& thisEdge = edges_[eIndex];

    return mag(points_[thisEdge[1]] - points_[thisEdge[0]]);
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
        if (patchCoupling_.found(pIndex))
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
    entitiesToAvoid_.clear();

    // Build a set of entities to avoid during regular modifications,
    // and build a master stack for coupled modifications.
    const polyBoundaryMesh& boundary = boundaryMesh();

    // Determine locally coupled slave patches.
    labelHashSet localSlavePatches;

    forAllIter(Map<coupledPatchInfo>, patchCoupling_, patchI)
    {
        localSlavePatches.insert(patchI().patchMap().slaveIndex());
    }

    // Loop through boundary faces and check whether
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

        // Check if this is a coupled face.
        if
        (
            patchCoupling_.found(pIndex) ||
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

// Set slave modification
inline void dynamicTopoFvMesh::setSlaveModification() const
{
    slaveModification_ = true;
}

// Unset slave modification
inline void dynamicTopoFvMesh::unsetSlaveModification() const
{
    slaveModification_ = false;
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
inline void dynamicTopoFvMesh::sizeUpList
(
    const label item,
    labelList& list
) const
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
) const
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
            cQuality = (*tetMetric_)(a, b, c, d);

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
            cQuality = (*tetMetric_)(a, b, c, d);

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
            forAllConstIter(Map<coupledPatchInfo>, patchCoupling_, patchI)
            {
                if ((sIndex = patchI().patchMap().findSlaveIndex(eIndex)) > -1)
                {
                    break;
                }
            }

            if (sIndex == -1)
            {
                FatalErrorIn
                (
                    "dynamicTopoFvMesh::computeMinQuality"
                ) << "Coupled maps were improperly specified." << nl
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
    const cell& curCell = cells_[cIndex];

    const face& currFace = faces_[curCell[0]];
    const face& nextFace = faces_[curCell[1]];

    // Get the fourth point
    forAll(nextFace, pointI)
    {
        if
        (
            nextFace[pointI] != currFace[0] &&
            nextFace[pointI] != currFace[1] &&
            nextFace[pointI] != currFace[2]
        )
        {
            // Compute cell-quality and write-out
            if (owner_[curCell[0]] == cIndex)
            {
                return
                (
                    (*tetMetric_)
                    (
                        points_[currFace[2]],
                        points_[currFace[1]],
                        points_[currFace[0]],
                        points_[nextFace[pointI]]
                    )
                );
            }
            else
            {
                return
                (
                    (*tetMetric_)
                    (
                        points_[currFace[0]],
                        points_[currFace[1]],
                        points_[currFace[2]],
                        points_[nextFace[pointI]]
                    )
                );
            }
        }
    }

    // This should never happen.
    FatalErrorIn
    (
        "dynamicTopoFvMesh::tetQuality() "
    )
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
    const label fIndex
) const
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
) const
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    const labelList& eFaces = edgeFaces_[eIndex];

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
                forAllConstIter
                (
                    Map<coupledPatchInfo>,
                    patchCoupling_,
                    iterI
                )
                {
                    if (iterI().patchMap().slaveIndex() == patch)
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
) const
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    label patch = whichPatch(fIndex);

    if (patch == -1)
    {
        FatalErrorIn("dynamicTopoFvMesh::processorCoupledFace()")
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

    forAllConstIter(Map<coupledPatchInfo>, patchCoupling_, patchI)
    {
        index = index > patchI.key() ? index : patchI.key();
    }

    return index;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
