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

#include "Pair.H"
#include "meshOps.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class T1, class T2>
inline const Vector<T1> convert(const Vector<T2>& v)
{
    Vector<T1> tv;

    tv.x() = v.x();
    tv.y() = v.y();
    tv.z() = v.z();

    return tv;
}


template<class T1, class T2>
inline const Field<Vector<T1> > convert(const Field<Vector<T2> >& f)
{
    Field<Vector<T1> > tf(f.size());

    forAll(tf, indexI)
    {
        tf[indexI] = convert<T1>(f[indexI]);
    }

    return tf;
}


namespace meshOps
{

// Compute the centroid for a given face, using UList
template <class T1, class T2>
inline void faceCentre
(
    const face& faceToCheck,
    const UList<Vector<T1> >& points,
    Vector<T2>& centre
)
{
    // Reset to zero
    centre = Vector<T2>::zero;

    const T2 three(3.0);
    Vector<T2> a(Vector<T2>::zero);
    Vector<T2> b(Vector<T2>::zero);
    Vector<T2> c(Vector<T2>::zero);

    // If the face is a triangle, do a direct calculation
    // to avoid round-off error-related problems
    if (faceToCheck.size() == 3)
    {
        a = convert<T2>(points[faceToCheck[0]]);
        b = convert<T2>(points[faceToCheck[1]]);
        c = convert<T2>(points[faceToCheck[2]]);

        centre = ((pTraits<T2>::one / three) * (a + b + c));

        return;
    }

    label nPoints = faceToCheck.size();

    register label pI;

    // Store the centre point in c
    for (pI = 0; pI < nPoints; pI++)
    {
        c += convert<T2>(points[faceToCheck[pI]]);
    }

    c /= nPoints;

    T2 ta = pTraits<T2>::zero;
    T2 sumA = pTraits<T2>::zero;
    Vector<T2> ttc = Vector<T2>::zero;
    Vector<T2> sumAc = Vector<T2>::zero;

    for (pI = 0; pI < nPoints; pI++)
    {
        a = convert<T2>(points[faceToCheck[pI]]);
        b = convert<T2>(points[faceToCheck[(pI + 1) % nPoints]]);

        // Calculate 3*triangle centre
        ttc = (a + b + c);

        // Calculate 2*triangle area
        ta = mag((a - c)^(b - c));

        sumA += ta;
        sumAc += ta*ttc;
    }

    if (sumA > VSMALL)
    {
        centre = (sumAc/(three * sumA));
    }
    else
    {
        centre = c;
    }
}


// Compute the normal for a given face, using UList
template <class T1, class T2>
inline void faceNormal
(
    const face& faceToCheck,
    const UList<Vector<T1> >& points,
    Vector<T2>& normal
)
{
    // Reset to zero
    normal = Vector<T2>::zero;

    const T2 pointFive(0.5);
    Vector<T2> a(Vector<T2>::zero);
    Vector<T2> b(Vector<T2>::zero);
    Vector<T2> c(Vector<T2>::zero);

    // If the face is a triangle, do a direct calculation
    // to avoid round-off error-related problems
    if (faceToCheck.size() == 3)
    {
        a = convert<T2>(points[faceToCheck[0]]);
        b = convert<T2>(points[faceToCheck[1]]);
        c = convert<T2>(points[faceToCheck[2]]);

        normal = (pointFive * ((b - a)^(c - a)));

        return;
    }

    label nPoints = faceToCheck.size();

    register label pI;

    // Store the centre point in c
    for (pI = 0; pI < nPoints; pI++)
    {
        c += convert<T2>(points[faceToCheck[pI]]);
    }

    c /= nPoints;

    for (pI = 0; pI < nPoints; pI++)
    {
        a = convert<T2>(points[faceToCheck[pI]]);
        b = convert<T2>(points[faceToCheck[(pI + 1) % nPoints]]);

        normal += (pointFive * ((b - a)^(c - a)));
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
//  - If checking is enabled, return whether the cell is closed
inline bool cellCentreAndVolume
(
    const label cIndex,
    const vectorField& points,
    const UList<face>& faces,
    const UList<cell>& cells,
    const UList<label>& owner,
    vector& centre,
    scalar& volume,
    bool checkClosed
)
{
    // Reset inputs
    volume = 0.0;
    centre = vector::zero;

    const cell& cellToCheck = cells[cIndex];

    // Average face-centres to get an estimate centroid
    vector cEst(vector::zero), fCentre(vector::zero), fArea(vector::zero);
    vector sumClosed(vector::zero), sumMagClosed(vector::zero);

    forAll(cellToCheck, faceI)
    {
        const face& checkFace = faces[cellToCheck[faceI]];

        if (checkFace.empty())
        {
            continue;
        }

        cEst += checkFace.centre(points);
    }

    cEst /= cellToCheck.size();

    forAll(cellToCheck, faceI)
    {
        const face& checkFace = faces[cellToCheck[faceI]];

        if (checkFace.empty())
        {
            continue;
        }

        fArea = checkFace.normal(points);
        fCentre = checkFace.centre(points);

        // Flip if necessary
        if (owner[cellToCheck[faceI]] != cIndex)
        {
            fArea *= -1.0;
        }

        // Calculate 3*face-pyramid volume
        scalar pyr3Vol = fArea & (fCentre - cEst);

        // Calculate face-pyramid centre
        vector pc = ((3.0 / 4.0) * fCentre) + ((1.0 / 4.0) * cEst);

        // Accumulate volume-weighted face-pyramid centre
        centre += pyr3Vol*pc;

        // Accumulate face-pyramid volume
        volume += pyr3Vol;

        // Accumulate areas, if necessary
        if (checkClosed)
        {
            sumClosed += fArea;
            sumMagClosed += cmptMag(fArea);
        }
    }

    centre /= volume + VSMALL;
    volume *= (1.0 / 3.0);

    if (checkClosed)
    {
        // Check the sum across components
        scalar closed = 0.0;

        for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
        {
            closed =
            (
                mag(sumClosed[cmpt]) / (sumMagClosed[cmpt] + VSMALL)
            );

            if (closed > 1e-06)
            {
                // Volume is invalid
                return false;
            }
        }
    }

    // Volume is valid
    return true;
}


// Given a cell index, find the centroid / volume of a cell.
template <class T1, class T2>
inline bool cellCentreAndVolumeT
(
    const label cIndex,
    const UList<Vector<T1> >& points,
    const UList<face>& faces,
    const UList<cell>& cells,
    const UList<label>& owner,
    Vector<T2>& centre,
    T2& volume,
    bool checkClosed = false
)
{
    // Reset inputs
    volume = pTraits<T2>::zero;
    centre = Vector<T2>::zero;

    const cell& cellToCheck = cells[cIndex];

    // Average face-centres to get an estimate centroid
    Vector<T2> cEst(Vector<T2>::zero);
    Vector<T2> fCentre(Vector<T2>::zero);
    Vector<T2> fArea(Vector<T2>::zero);

    const T2 three(3.0), four(4.0);
    const T2 oneThird = (pTraits<T2>::one / three);
    const T2 oneFourth = (pTraits<T2>::one / four);
    const T2 threeFourths = (three / four);

    forAll(cellToCheck, faceI)
    {
        const face& checkFace = faces[cellToCheck[faceI]];

        if (checkFace.empty())
        {
            continue;
        }

        meshOps::faceCentre
        (
            checkFace,
            points,
            fCentre
        );

        cEst += fCentre;
    }

    cEst /= T2(cellToCheck.size());

    forAll(cellToCheck, faceI)
    {
        const face& checkFace = faces[cellToCheck[faceI]];

        if (checkFace.empty())
        {
            continue;
        }

        meshOps::faceNormal
        (
            checkFace,
            points,
            fArea
        );

        meshOps::faceCentre
        (
            checkFace,
            points,
            fCentre
        );

        // Flip if necessary
        if (owner[cellToCheck[faceI]] != cIndex)
        {
            fArea *= -1.0;
        }

        // Calculate 3*face-pyramid volume
        T2 pyr3Vol = fArea & (fCentre - cEst);

        // Calculate face-pyramid centre
        Vector<T2> pc = (threeFourths*fCentre) + (oneFourth*cEst);

        // Accumulate volume-weighted face-pyramid centre
        centre += pyr3Vol*pc;

        // Accumulate face-pyramid volume
        volume += pyr3Vol;
    }

    centre /= volume + VSMALL;
    volume *= oneThird;

    // Volume is valid
    return true;
}


// Determine whether a point lies on a segment
template <class T>
inline bool pointSegmentIntersection
(
    const line<Vector<T>, const Vector<T>&>& segment,
    const Vector<T>& cP,
    const T& matchTol
)
{
    // Fetch references
    const Vector<T>& sS = segment.start();
    const Vector<T>& sE = segment.end();

    // Define custom vectors
    Vector<T> u = (sE - sS);
    Vector<T> v = (cP - sS);

    T magU = mag(u);
    T magV = mag(v);

    // Stabilise around zero
    magU = (magU < VSMALL) ? (magU + VSMALL) : magU;
    magV = (magV < VSMALL) ? (magV + VSMALL) : magV;

    // Compare dot-products
    T collinearity = mag(pTraits<T>::one - ((u/magU) & (v/magV)));
    T tolerance = (matchTol * magU);

    if (collinearity > tolerance)
    {
        return false;
    }

    // Compute uValue
    T uValue = (u & v) / ((u & u) + VSMALL);

    // Does point fall-off ends?
    // Note: Does not check if end-points coincide
    //       with checkPoint, so be careful here.
    if (uValue < pTraits<T>::zero || uValue > pTraits<T>::one)
    {
        return false;
    }

    // Intersects segment
    return true;
}


// Determine whether a two line segments intersect
template <class T>
inline bool segmentSegmentIntersection
(
    const line<Vector<T>, const Vector<T>&>& oldSegment,
    const line<Vector<T>, const Vector<T>&>& newSegment,
    const T& matchTol,
    Vector<T>& intPoint
)
{
    // Fetch references
    const Vector<T>& oS = oldSegment.start(), oE = oldSegment.end();
    const Vector<T>& nS = newSegment.start(), nE = newSegment.end();

    // Define custom vectors
    Vector<T> u = (oE - oS);
    Vector<T> v = (nE - nS);

    const Vector<T>& p = oS;
    const Vector<T>& q = nS;

    Vector<T> w = (p - q);

    // Get parameters
    T a = (u & u), b = (u & v), c = (v & v), d = (u & w), e = (v & w);

    // Parallel / collinear check
    T denom = (a * c) - (b * b);

    if (mag(denom) < VSMALL)
    {
        return false;
    }

    T s = ((b * e) - (c * d))/denom;
    T t = ((a * e) - (b * d))/denom;

    // Out-of-bounds check
    // Note: Does not check if end-points coincide,
    //       so be careful here.
    if
    (
        s < pTraits<T>::zero || t < pTraits<T>::zero ||
        s > pTraits<T>::one  || t > pTraits<T>::one
    )
    {
        return false;
    }

    // Proximity check
    T dist = mag( w + (s * u) - (t * v) );
    T tolerance = (matchTol * Foam::min(mag(u), mag(v)));

    if (dist > tolerance)
    {
        return false;
    }

    // Calculate the intersection point
    intPoint = p + (s * u);

    return true;
}


// Determine whether a segment intersects a triangular face
template <class T>
inline bool segmentTriFaceIntersection
(
    const triangle<Vector<T>, const Vector<T>&>& faceToCheck,
    const line<Vector<T>, const Vector<T>&>& edgeToCheck,
    const T& matchTol,
    Vector<T>& intPoint
)
{
    // Fetch references
    const Vector<T>& p1 = edgeToCheck.start();
    const Vector<T>& p2 = edgeToCheck.end();

    const Vector<T>& a = faceToCheck.a();
    const Vector<T>& b = faceToCheck.b();
    const Vector<T>& c = faceToCheck.c();

    // Define custom face-normal
    Vector<T> n = (T(0.5) * ((b - a)^(c - a)));
    n /= mag(n) + VSMALL;

    // Compute uValue
    T numerator = n & (a - p1);
    T denominator = n & (p2 - p1);

    // Check if the edge is parallel to the face
    if (mag(denominator) < VSMALL)
    {
        return false;
    }

    T u = (numerator / denominator);
    T tolerance = (matchTol * mag(p2 - p1));

    // Check for intersection along line.
    if ((u > tolerance) && (u < (pTraits<T>::one - tolerance)))
    {
        // Compute point of intersection
        intPoint = p1 + u*(p2 - p1);

        // Also make sure that intPoint lies within the face
        if (pointInTriFace(faceToCheck, intPoint, matchTol, false))
        {
            return true;
        }
    }

    // Failed to fall within edge-bounds, or within face
    return false;
}


// Determine whether the particular point lies
// inside the given triangular face
template <class T>
inline bool pointInTriFace
(
    const triangle<Vector<T>, const Vector<T>&>& faceToCheck,
    const Vector<T>& cP,
    const T& matchTol,
    bool testCoplanarity
)
{
    // Copy inputs
    const Vector<T>& a = faceToCheck.a();
    const Vector<T>& b = faceToCheck.b();
    const Vector<T>& c = faceToCheck.c();

    const T pointFive(0.5);

    // Compute the normal
    Vector<T> nf = (pointFive * ((b - a)^(c - a)));

    if ( ((pointFive * ((b - a)^(cP - a))) & nf) < pTraits<T>::zero)
    {
        return false;
    }

    if ( ((pointFive * ((c - b)^(cP - b))) & nf) < pTraits<T>::zero)
    {
        return false;
    }

    if ( ((pointFive * ((a - c)^(cP - c))) & nf) < pTraits<T>::zero)
    {
        return false;
    }

    if (testCoplanarity)
    {
        Vector<T> ftp(a - cP);

        // Normalize vectors
        nf /= mag(nf) + VSMALL;
        ftp /= mag(ftp) + VSMALL;

        if (mag(ftp & nf) > matchTol)
        {
            return false;
        }
    }

    // Passed test with all edges
    return true;
}


// Determine whether the particular point lies
// inside the given cell
template <class T>
inline bool pointInCell
(
    const label cIndex,
    const cell& cellToCheck,
    const UList<face>& faces,
    const UList<label>& owner,
    const UList<point>& points,
    const T& matchTol,
    const Vector<T>& cP
)
{
    Vector<T> xf, nf;

    forAll(cellToCheck, faceI)
    {
        meshOps::faceCentre(faces[cellToCheck[faceI]], points, xf);
        meshOps::faceNormal(faces[cellToCheck[faceI]], points, nf);

        Vector<T> faceToPoint(xf - cP);

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
template <class T>
inline void insertPointLabels
(
    const Vector<T>& refNorm,
    const Field<Vector<T> >& points,
    const labelHashSet& pLabels,
    face& modFace
)
{
    // Need to configure a new face.
    face newFace(modFace);

    const T pointFive(0.5);

    forAllConstIter(labelHashSet, pLabels, pIter)
    {
        forAll(newFace, pI)
        {
            label nI = newFace.fcIndex(pI);

            const Vector<T>& a = points[newFace[pI]];
            const Vector<T>& b = points[pIter.key()];
            const Vector<T>& c = points[newFace[nI]];

            // Compute the normal.
            Vector<T> newNorm = (pointFive * ((b - a)^(c - a)));

            if ((refNorm & newNorm) > pTraits<T>::zero)
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
    label index = -1;

    if ((index = findIndex(list, original)) > -1)
    {
        list[index] = replacement;
    }
    else
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
    const Type& item,
    List<Type>& list
)
{
    list.setSize(list.size() + 1, item);
}


// Utility method to size-down the list to remove an item
template <class Type>
inline void sizeDownList
(
    const Type& item,
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
            "(const Type& item, List<Type>& list)"
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


} // End namespace Foam

// ************************************************************************* //
