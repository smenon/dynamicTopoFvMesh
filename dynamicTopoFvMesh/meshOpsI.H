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
    meshOps

Description
    Various utility functions that perform mesh-related operations.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "triPointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace meshOps
{

// Return the magnitude of an edge, using UList
inline scalar edgeLength
(
    const edge& edgeToCheck,
    const UList<point>& points
)
{
    return mag(points[edgeToCheck.end()] - points[edgeToCheck.start()]);
}


// Compute the centroid for a given face, using UList
inline vector faceCentre
(
    const face& faceToCheck,
    const UList<point>& points
)
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


// Compute the normal for a given face, using UList
inline vector faceNormal
(
    const face& faceToCheck,
    const UList<point>& points
)
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


// Utility method to find the common edge between two faces.
inline bool findCommonEdge
(
    const label first,
    const label second,
    const UList<labelList>& faceEdges,
    label& common
)
{
    bool found = false;

    const labelList& fEi = faceEdges[first];
    const labelList& fEj = faceEdges[second];

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


// Utility method to find the interior (quad) / boundary (tri) faces
// for an input quad-face and adjacent triangle-prism cell.
inline void findPrismFaces
(
    const label fIndex,
    const label cIndex,
    const UList<face>& faces,
    const UList<cell>& cells,
    const UList<label>& neighbour,
    FixedList<face,2>& bdyf,
    FixedList<label,2>& bidx,
    FixedList<face,2>& intf,
    FixedList<label,2>& iidx
)
{
    label indexO = 0, indexI = 0;

    const cell& c = cells[cIndex];

    forAll(c, i)
    {
        label faceIndex = c[i];

        // Don't count the face under consideration
        if (faceIndex != fIndex)
        {
            const face& fi = faces[faceIndex];

            if (neighbour[faceIndex] == -1)
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


// Utility method to find the isolated point given two triangular faces.
//  - Returns the point on checkFace that does not belong to baseFace.
inline label findIsolatedPoint
(
    const face& baseFace,
    const face& checkFace
)
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
inline void findIsolatedPoint
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
        "inline void meshOps::findIsolatedPoint"
        "(const face&, const edge&, label&, label&)"
    )
        << "Cannot find isolated point in face " << f << endl
        << " Using edge: " << e
        << abort(FatalError);
}


// Given a face and cell index, find the apex point for a tet cell.
inline label tetApexPoint
(
    const label cIndex,
    const label fIndex,
    const UList<face> faces,
    const UList<cell> cells
)
{
    label apexPoint = -1;
    bool foundApex = false;

    const cell& cellToCheck = cells[cIndex];
    const face& baseFace = faces[fIndex];

    forAll(cellToCheck, faceI)
    {
        const face& faceToCheck = faces[cellToCheck[faceI]];

        apexPoint = findIsolatedPoint(baseFace, faceToCheck);

        if (foundApex > -1)
        {
            foundApex = true;
            break;
        }
    }

    if (!foundApex)
    {
        Info << "fIndex: " << fIndex << ":: " << faces[fIndex] << endl;
        Info << "cIndex: " << cIndex << ":: " << cellToCheck << endl;

        forAll(cellToCheck, faceI)
        {
            Info << '\t' << cellToCheck[faceI] << ":: "
                 << faces[cellToCheck[faceI]] << endl;
        }

        FatalErrorIn
        (
            "\n\n"
            "inline label meshOps::tetApexPoint\n"
            "(\n"
            "    const label cIndex,\n"
            "    const label fIndex,\n"
            "    const UList<face> faces,\n"
            "    const UList<cell> cells\n"
            ")\n"
        )
            << "Could not find an apex point in cell: " << cIndex
            << " given face: " << fIndex
            << abort(FatalError);
    }

    return apexPoint;
}


// Given a cell index, find the centroid / volume of a cell.
inline void cellCentreAndVolume
(
    const label cIndex,
    const UList<point>& points,
    const UList<face>& faces,
    const UList<cell>& cells,
    const UList<label>& owner,
    vector& centre,
    scalar& volume
)
{
    // Reset inputs
    volume = 0.0;
    centre = vector::zero;

    const cell& cellToCheck = cells[cIndex];

    // Average face-centres to get an estimate centroid
    vector cEst = vector::zero;

    forAll(cellToCheck, faceI)
    {
        cEst +=
        (
            meshOps::faceCentre
            (
                faces[cellToCheck[faceI]],
                points
            )
        );
    }

    cEst /= cellToCheck.size();

    forAll(cellToCheck, faceI)
    {
        vector fArea =
        (
            meshOps::faceNormal
            (
                faces[cellToCheck[faceI]],
                points
            )
        );

        vector fCentre =
        (
            meshOps::faceCentre
            (
                faces[cellToCheck[faceI]],
                points
            )
        );

        if (owner[cellToCheck[faceI]] != cIndex)
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
inline label whichSide
(
    const labelList& cellPoints,
    const UList<point>& points,
    const point& dir,
    const point& p
)
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


// Determine whether a point lies on a segment
inline bool pointSegmentIntersection
(
    const edge& segment,
    const UList<point>& points,
    const point& checkPoint
)
{
    vector u = points[segment.end()] - points[segment.start()];
    vector v = checkPoint - points[segment.start()];

    scalar magU = mag(u) + VSMALL;
    scalar magV = mag(v) + VSMALL;

    scalar tolerance = (meshOps::matchTol_ * magU);

    // Compare dot-products
    if ( 1.0 - ((u/magU) & (v/magV)) > tolerance )
    {
        return false;
    }

    // Compute uValue
    scalar uValue = (u & v) / ((u & u) + VSMALL);

    // Does point fall-off ends?
    // Note: Does not check if end-points coincide
    //       with checkPoint, so be careful here.
    if (uValue < 0.0 || uValue > 1.0)
    {
        return false;
    }

    // Intersects segment
    return true;
}


// Determine whether a two line segments intersect
inline bool segmentSegmentIntersection
(
    const edge& fromSegment,
    const edge& toSegment,
    const UList<point>& points,
    point& intPoint
)
{
    vector u = points[fromSegment.end()] - points[fromSegment.start()];
    vector v = points[toSegment.end()] - points[toSegment.start()];

    vector p = points[fromSegment.start()];
    vector q = points[toSegment.start()];

    vector w = (p - q);

    // Get parameters
    scalar a = (u & u), b = (u & v), c = (v & v), d = (u & w), e = (v & w);

    // Parallel / collinear check
    scalar denom = (a * c) - (b * b);

    if (mag(denom) < VSMALL)
    {
        return false;
    }

    scalar s = ((b * e) - (c * d))/(denom);
    scalar t = ((a * e) - (b * d))/(denom);

    // Out-of-bounds check
    if (s < 0.0 || t < 0.0 || s > 1.0 || t > 1.0)
    {
        return false;
    }

    // Proximity check
    scalar dist = mag( w + (s * u) - (t * v) );

    scalar tolerance = (meshOps::matchTol_ * Foam::min(mag(u), mag(v)));

    if (dist > tolerance)
    {
        return false;
    }

    // Calculate the intersection point
    intPoint = points[fromSegment.start()] + (s * u);

    return true;
}


// Determine whether a segment intersects a face
inline bool segmentFaceIntersection
(
    const edge& edgeToCheck,
    const face& faceToCheck,
    const UList<point>& points,
    vector& intersectionPoint
)
{
    // Specify unit face-normals
    vector n = meshOps::faceNormal(faceToCheck, points);
    n /= mag(n) + VSMALL;

    point p1 = points[edgeToCheck[0]];
    point p2 = points[edgeToCheck[1]];
    point p3 = points[faceToCheck[0]];

    // Compute uValue
    scalar numerator = n & (p3 - p1);
    scalar denominator = n & (p2 - p1);

    // Check if the edge is parallel to the face
    if (mag(denominator) < VSMALL)
    {
        return false;
    }

    scalar u = (numerator / denominator);

    scalar tolerance = (meshOps::matchTol_ * mag(p2 - p1));

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

    // Failed to fall within edge-bounds, or within face
    return false;
}


// Determine whether the particular point lies
// inside the given face
inline bool pointInFace
(
    const face& faceToCheck,
    const UList<point>& points,
    const point& checkPoint
)
{
    vector nf = meshOps::faceNormal(faceToCheck, points);
    nf /= mag(nf) + VSMALL;

    forAll(faceToCheck, pI)
    {
        label nI = faceToCheck.fcIndex(pI);

        // Identify vectors
        vector tn =
        (
            triPointRef
            (
                points[faceToCheck[pI]],
                points[faceToCheck[nI]],
                checkPoint
            ).normal()
        );

        scalar area = (tn & nf);

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
inline bool pointInCell
(
    const label cIndex,
    const cell& cellToCheck,
    const UList<face>& faces,
    const UList<label>& owner,
    const UList<point>& points,
    const point& checkPoint,
    const scalar& tolerance
)
{
    vector xf(vector::zero), nf(vector::zero);

    forAll(cellToCheck, faceI)
    {
        xf = meshOps::faceCentre(faces[cellToCheck[faceI]], points);
        nf = meshOps::faceNormal(faces[cellToCheck[faceI]], points);

        vector faceToPoint = (xf - checkPoint);

        if ((faceToPoint & nf) > (tolerance * magSqr(faceToPoint)))
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


// Method to insert labels in a face, so that
// right-handedness is preserved.
inline void insertPointLabels
(
    const vector& refNorm,
    const vectorField& points,
    const labelHashSet& pLabels,
    face& modFace
)
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
                meshOps::insertLabel
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
inline void insertLabel
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
            "inline void meshOps::insertLabel"
            "(const label, const label, const label, labelList&)"
        )   << nl << "Cannot insert " << newLabel
            << " in list: " << list << nl
            << " Labels: "
            << labelA << " and " << labelB
            << " were not found in sequence."
            << abort(FatalError);
    }

    // Transfer the list
    list.transfer(newList);
}


// Utility method to replace a label in a given list
inline void replaceLabel
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
            "inline void label meshOps::replaceLabel"
            "(const label, const label, labelList&)"
        )   << nl << "Cannot find " << original
            << " in list: " << list << nl
            << " Label: " << replacement
            << " was not used in replacement."
            << abort(FatalError);
    }
}


// Utility method to size-up the list to include an item
template <class Type>
inline void sizeUpList
(
    const Type item,
    List<Type>& list
)
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
inline void sizeDownList
(
    const Type item,
    List<Type>& list
)
{
    label index = -1;

    if ((index = findIndex(list, item)) > -1)
    {
        meshOps::removeIndex(index, list);
    }
    else
    {
        FatalErrorIn
        (
            "inline void meshOps::sizeDownList"
            "(const Type item, List<Type>& list)"
        )
            << nl << "Item: " << item
            << " was not found in list. " << nl
            << " List: " << nl << list
            << abort(FatalError);
    }
}


// Remove an item at a particular index in the list
template <class Type>
inline void removeIndex
(
    const label index,
    List<Type>& list
)
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


} // End namespace meshOps


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
