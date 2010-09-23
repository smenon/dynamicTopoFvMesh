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

// Compute the centroid for a given face, using UList
template <class scal, class scal2>
inline void faceCentre
(
    const face& faceToCheck,
    const UList<Vector<scal2> >& points,
    Vector<scal>& centre
)
{
    Vector<scal> a, b, c;

    // If the face is a triangle, do a direct calculation
    // to avoid round-off error-related problems
    if (faceToCheck.size() == 3)
    {
        a = points[faceToCheck[0]];
        b = points[faceToCheck[1]];
        c = points[faceToCheck[2]];

        centre = (pTraits<scal>::one / scal(3.0)) * (a + b + c);

        return;
    }

    // Reset to zero
    centre = Vector<scal>::zero;
    label nPoints = faceToCheck.size();

    register label pI;

    // Store the centre point in c
    for (pI = 0; pI < nPoints; pI++)
    {
        c += points[faceToCheck[pI]];
    }

    c /= nPoints;

    scal ta = pTraits<scal>::zero;
    scal sumA = pTraits<scal>::zero;
    Vector<scal> ttc = Vector<scal>::zero;
    Vector<scal> sumAc = Vector<scal>::zero;

    for (pI = 0; pI < nPoints; pI++)
    {
        a = points[faceToCheck[pI]];
        b = points[faceToCheck[(pI + 1) % nPoints]];

        // Calculate 3*triangle centre
        ttc = (a + b + c);

        // Calculate 2*triangle area
        ta = mag((a - c)^(b - c));

        sumA += ta;
        sumAc += ta*ttc;
    }

    if (sumA > VSMALL)
    {
        centre = (sumAc/(3.0 * sumA));
    }
    else
    {
        centre = c;
    }
}


// Compute the normal for a given face, using UList
template <class scal, class scal2>
inline void faceNormal
(
    const face& faceToCheck,
    const UList<Vector<scal2> >& points,
    Vector<scal>& normal
)
{
    Vector<scal> a, b, c;

    // If the face is a triangle, do a direct calculation
    // to avoid round-off error-related problems
    if (faceToCheck.size() == 3)
    {
        a = points[faceToCheck[0]];
        b = points[faceToCheck[1]];
        c = points[faceToCheck[2]];

        normal = 0.5 * ((b - a)^(c - a));

        return;
    }

    // Reset to zero
    normal = Vector<scal>::zero;
    label nPoints = faceToCheck.size();

    register label pI;

    // Store the centre point in c
    for (pI = 0; pI < nPoints; pI++)
    {
        c += points[faceToCheck[pI]];
    }

    c /= nPoints;

    for (pI = 0; pI < nPoints; pI++)
    {
        // Define this point
        a = points[faceToCheck[pI]];
        b = points[faceToCheck[(pI + 1) % nPoints]];

        normal += (0.5 * ((b - a)^(c - a)));
    }
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
template <class scal, class scal2>
inline void cellCentreAndVolume
(
    const label cIndex,
    const UList<Vector<scal2> >& points,
    const UList<face>& faces,
    const UList<cell>& cells,
    const UList<label>& owner,
    Vector<scal>& centre,
    scal& volume
)
{
    // Reset inputs
    volume = pTraits<scal>::zero;
    centre = Vector<scal>::zero;

    const cell& cellToCheck = cells[cIndex];

    // Average face-centres to get an estimate centroid
    Vector<scal> cEst(Vector<scal>::zero);
    Vector<scal> fCentre(Vector<scal>::zero);
    Vector<scal> fArea(Vector<scal>::zero);

    scal oneThird = (pTraits<scal>::one / scal(3.0));
    scal oneFourth = (pTraits<scal>::one / scal(4.0));
    scal threeFourths = (scal(3.0) / scal(4.0));

    forAll(cellToCheck, faceI)
    {
        meshOps::faceCentre
        (
            faces[cellToCheck[faceI]],
            points,
            fCentre
        );

        cEst += fCentre;
    }

    cEst /= scal(cellToCheck.size());

    forAll(cellToCheck, faceI)
    {
        const face& checkFace = faces[cellToCheck[faceI]];

        if (checkFace.empty())
        {
            continue;
        }

        meshOps::faceNormal
        (
            faces[cellToCheck[faceI]],
            points,
            fArea
        );

        meshOps::faceCentre
        (
            faces[cellToCheck[faceI]],
            points,
            fCentre
        );

        // Flip if necessary
        if (owner[cellToCheck[faceI]] != cIndex)
        {
            fArea *= -1.0;
        }

        // Calculate 3*face-pyramid volume
        scal pyr3Vol = fArea & (fCentre - cEst);

        // Calculate face-pyramid centre
        Vector<scal> pc = (threeFourths*fCentre) + (oneFourth*cEst);

        // Accumulate volume-weighted face-pyramid centre
        centre += pyr3Vol*pc;

        // Accumulate face-pyramid volume
        volume += pyr3Vol;
    }

    centre /= volume + VSMALL;
    volume *= oneThird;
}


// Determine whether a point lies on a segment
template <class scal>
inline bool pointSegmentIntersection
(
    const linePointRef& segment,
    const point& cP,
    const scal& matchTol
)
{
    // Fetch references
    const point& sS = segment.start();
    const point& sE = segment.end();

    // Define custom vectors
    Vector<scal> u = (Vector<scal>(sE) - Vector<scal>(sS));
    Vector<scal> v = (Vector<scal>(cP) - Vector<scal>(sS));

    scal magU = mag(u);
    scal magV = mag(v);

    // Stabilise around zero
    magU = (magU < VSMALL) ? (magU + VSMALL) : magU;
    magV = (magV < VSMALL) ? (magV + VSMALL) : magV;

    // Compare dot-products
    scal collinearity = mag(pTraits<scal>::one - ((u/magU) & (v/magV)));
    scal tolerance = (matchTol * magU);

    if (collinearity > tolerance)
    {
        return false;
    }

    // Compute uValue
    scal uValue = (u & v) / ((u & u) + VSMALL);

    // Does point fall-off ends?
    // Note: Does not check if end-points coincide
    //       with checkPoint, so be careful here.
    if (uValue < pTraits<scal>::zero || uValue > pTraits<scal>::one)
    {
        return false;
    }

    // Intersects segment
    return true;
}


// Determine whether a two line segments intersect
template <class scal>
inline bool segmentSegmentIntersection
(
    const linePointRef& fromSegment,
    const linePointRef& toSegment,
    const scal& matchTol,
    Vector<scal>& intPoint
)
{
    // Fetch references
    const point& fS = fromSegment.start(), fE = fromSegment.end();
    const point& tS = toSegment.start(), tE = toSegment.end();

    // Define custom vectors
    Vector<scal> u = (Vector<scal>(fE) - Vector<scal>(fS));
    Vector<scal> v = (Vector<scal>(tE) - Vector<scal>(tS));

    Vector<scal> p = Vector<scal>(fS);
    Vector<scal> q = Vector<scal>(tS);

    Vector<scal> w = (p - q);

    // Get parameters
    scal a = (u & u), b = (u & v), c = (v & v), d = (u & w), e = (v & w);

    // Parallel / collinear check
    scal denom = (a * c) - (b * b);

    if (mag(denom) < VSMALL)
    {
        return false;
    }

    scal s = ((b * e) - (c * d))/denom;
    scal t = ((a * e) - (b * d))/denom;

    // Out-of-bounds check
    // Note: Does not check if end-points coincide,
    //       so be careful here.
    if
    (
        s < pTraits<scal>::zero ||
        t < pTraits<scal>::zero ||
        s > pTraits<scal>::one ||
        t > pTraits<scal>::one
    )
    {
        return false;
    }

    // Proximity check
    scal dist = mag( w + (s * u) - (t * v) );
    scal tolerance = (scal(matchTol) * Foam::min(mag(u), mag(v)));

    if (dist > tolerance)
    {
        return false;
    }

    // Calculate the intersection point
    intPoint = p + (s * u);

    return true;
}


// Determine whether a segment intersects a triangular face
template <class scal>
inline bool segmentTriFaceIntersection
(
    const triPointRef& faceToCheck,
    const linePointRef& edgeToCheck,
    const scal& matchTol,
    Vector<scal>& intPoint
)
{
    // Define custom vectors
    const Vector<scal> p1(edgeToCheck.start());
    const Vector<scal> p2(edgeToCheck.end());

    const Vector<scal> pfA(faceToCheck.a());
    const Vector<scal> pfB(faceToCheck.b());
    const Vector<scal> pfC(faceToCheck.c());

    // Define custom face-normal
    Vector<scal> n(0.5 * ((pfB - pfA)^(pfC - pfA)));
    n /= mag(n) + VSMALL;

    // Compute uValue
    scal numerator = n & (pfA - p1);
    scal denominator = n & (p2 - p1);

    // Check if the edge is parallel to the face
    if (mag(denominator) < VSMALL)
    {
        return false;
    }

    scal u = (numerator / denominator);
    scal tolerance = (scal(matchTol) * mag(p2 - p1));

    // Check for intersection along line.
    if ((u > tolerance) && (u < (pTraits<scal>::one - tolerance)))
    {
        // Compute point of intersection
        intPoint = p1 + u*(p2 - p1);

        // Also make sure that intPoint lies within the face
        if
        (
            pointInTriFace
            (
                triangle<Vector<scal>, const Vector<scal>&>(pfA, pfB, pfC),
                intPoint
            )
        )
        {
            return true;
        }
    }

    // Failed to fall within edge-bounds, or within face
    return false;
}


// Determine whether the particular point lies
// inside the given triangular face
template <class scal>
inline bool pointInTriFace
(
    const triangle<Vector<scal>, const Vector<scal>&>& faceToCheck,
    const Vector<scal>& cP
)
{
    // Fetch references
    const Vector<scal>& a = faceToCheck.a();
    const Vector<scal>& b = faceToCheck.b();
    const Vector<scal>& c = faceToCheck.c();

    // Compute the normal
    Vector<scal> nf = (0.5 * ((b - a)^(c - a)));

    if ( ((0.5 * ((b - a)^(cP - a))) & nf) < pTraits<scal>::zero)
    {
        return false;
    }

    if ( ((0.5 * ((c - b)^(cP - b))) & nf) < pTraits<scal>::zero)
    {
        return false;
    }

    if ( ((0.5 * ((a - c)^(cP - c))) & nf) < pTraits<scal>::zero)
    {
        return false;
    }

    // Passed test with all edges
    return true;
}


// Determine whether the particular point lies
// inside the given cell
template <class scal>
inline bool pointInCell
(
    const label cIndex,
    const cell& cellToCheck,
    const UList<face>& faces,
    const UList<label>& owner,
    const UList<point>& points,
    const scal& matchTol,
    const point& cP
)
{
    Vector<scal> xf, nf;

    forAll(cellToCheck, faceI)
    {
        meshOps::faceCentre(faces[cellToCheck[faceI]], points, xf);
        meshOps::faceNormal(faces[cellToCheck[faceI]], points, nf);

        Vector<scal> faceToPoint(xf - Vector<scal>(cP));

        if ((faceToPoint & nf) > matchTol)
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
