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
    faceSetAlgorithm

Description
    Class for convexSetAlgorithms pertaining to faces

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "faceSetAlgorithm.H"

#include "meshOps.H"
#include "triFace.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
faceSetAlgorithm::faceSetAlgorithm
(
    const polyMesh& mesh,
    const pointField& newPoints,
    const UList<edge>& newEdges,
    const UList<face>& newFaces,
    const UList<cell>& newCells,
    const UList<label>& newOwner,
    const UList<label>& newNeighbour,
    const List<objectMap>& pointsFromPoints,
    const Map<labelList>& modPoints
)
:
    convexSetAlgorithm
    (
        mesh,
        newPoints,
        newEdges,
        newFaces,
        newCells,
        newOwner,
        newNeighbour,
        pointsFromPoints,
        modPoints
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void faceSetAlgorithm::computeNormFactor(const label index) const
{
    // Clear existing fields
    parents_.clear();

    centres_.clear();
    weights_.clear();

#   if USE_MPFR

    mpCentres_.clear();
    mpWeights_.clear();

    if (highPrecision())
    {
        meshOps::faceNormal(newFaces_[index], newPoints_, mpRefNorm_);

        mpNormFactor_ = mag(mpRefNorm_);

        // Normalize for later use
        mpRefNorm_ /= mpNormFactor_ + VSMALL;
    }
    else
#   endif
    {
        meshOps::faceNormal(newFaces_[index], newPoints_, refNorm_);

        normFactor_ = mag(refNorm_);

        // Normalize for later use
        refNorm_ /= normFactor_ + VSMALL;
    }
}


// Compute intersections
bool faceSetAlgorithm::computeIntersection
(
    const label newIndex,
    const label oldIndex,
    const scalar& matchTol,
    bool output
) const
{
    bool intersects = false;

#   if USE_MPFR
    if (convexSetAlgorithm::highPrecision())
    {
        Field<mpVector> mpIntPoints;
        const mpScalar mpMatchTol(matchTol);

        // Invoke the high-precision variant
        intersects =
        (
            faceIntersection
            (
                newIndex,
                oldIndex,
                mpMatchTol,
                mpIntPoints,
                output
            )
        );

        if (intersects)
        {
            mpScalar area;
            mpVector centre;

            convexSetArea
            (
                newIndex,
                oldIndex,
                mpIntPoints,
                mpRefNorm_,
                area,
                centre,
                output
            );

            // Size-up the internal lists
            if (!output)
            {
                meshOps::sizeUpList(oldIndex, parents_);
                meshOps::sizeUpList(area, mpWeights_);
                meshOps::sizeUpList(centre, mpCentres_);
            }
        }
    }
    else
#   endif
    {
        vectorField intPoints;

        // Invoke the conventional variant
        intersects =
        (
            faceIntersection
            (
                newIndex,
                oldIndex,
                matchTol,
                intPoints,
                output
            )
        );

        if (intersects)
        {
            scalar area;
            vector centre;

            convexSetArea
            (
                newIndex,
                oldIndex,
                intPoints,
                refNorm_,
                area,
                centre,
                output
            );

            // Size-up the internal lists
            if (!output)
            {
                meshOps::sizeUpList(oldIndex, parents_);
                meshOps::sizeUpList(area, weights_);
                meshOps::sizeUpList(centre, centres_);
            }
        }
    }

    return intersects;
}


template <class T>
bool faceSetAlgorithm::faceIntersection
(
    const label newIndex,
    const label oldIndex,
    const T& matchTol,
    Field<Vector<T> >& intPoints,
    bool output
) const
{
    // Reset inputs
    intPoints.clear();

    // Fetch face references for each mesh
    const face& newFace = this->newFaces_[newIndex];
    const face& oldFace = this->mesh_.faces()[oldIndex];

    // Alias references
    const UList<point>& newPoints = this->newPoints_;
    const UList<point>& oldPoints = this->mesh_.points();

    const Map<labelList>& modPoints = this->modPoints_;
    const List<objectMap>& pfp = this->pointsFromPoints_;

    // Obtain face centre and projection normal
    Vector<T> xf, nf;

    meshOps::faceCentre(newFace, newPoints, xf);
    meshOps::faceNormal(newFace, newPoints, nf);

    nf /= mag(nf) + VSMALL;

    // Track all possible intersections from here on.
    label nInts = 0;
    Map<Vector<T> > intersections;
    Vector<T> intPoint = Vector<T>::zero;

    // Topologically check for common points,
    // and project uniques ones on to the face plane
    Map<label> projPoints;
    Map<labelList> commonPoints;
    Field<Vector<T> > projections(oldFace.size(), Vector<T>::zero);

    forAll(oldFace, pointI)
    {
        label oldPoint = oldFace[pointI];
        label pIndex = findIndex(newFace, oldPoint);

        Vector<T> r = convert<T>(oldPoints[oldPoint]);

        if (pIndex == -1)
        {
            // Project this point on to the newFace plane.
            projections[pointI] = xf + ((r - xf) - ((r - xf) & nf)*nf);

            projPoints.insert(oldPoint, pointI);
        }
        else
        {
            // If this point was modified by a collapse
            // to an edge mid-point, it can't be a common point.
            if (modPoints.found(newFace[pIndex]))
            {
                // Project this point on to the newFace plane.
                projections[pointI] = xf + ((r - xf) - ((r - xf) & nf)*nf);

                projPoints.insert(oldPoint, pointI);
            }
            else
            {
                commonPoints.insert(newFace[pIndex], labelList(0));

                projections[pointI] = r;

                intersections.set(++nInts, r);
            }
        }
    }

    // Add points if they resulted from
    // bisections of old face edges.
    forAll(newFace, pointI)
    {
        label pIndex = newFace[pointI];

        if (pIndex >= nOldPoints_)
        {
            // Check pointsFromPoints info
            label index = -1;

            forAll(pfp, indexI)
            {
                if (pfp[indexI].index() == pIndex)
                {
                    index = indexI;
                    break;
                }
            }

            const labelList& mObj = pfp[index].masterObjects();

            // Check if the old face contains all master points
            bool allMaster = true;

            forAll(mObj, pointJ)
            {
                if (findIndex(oldFace, mObj[pointJ]) == -1)
                {
                    allMaster = false;
                    break;
                }
            }

            if (allMaster)
            {
                commonPoints.insert(newFace[pointI], mObj);

                intersections.set
                (
                    ++nInts,
                    convert<T>(newPoints[newFace[pointI]])
                );
            }
        }
    }

    // If all points are common, this is identical to the old face.
    if (nInts == oldFace.size())
    {
        // Copy intersections
        intPoints.setSize(nInts, Vector<T>::zero);

        nInts = 0;

        forAllConstIter(typename Map<Vector<T> >, intersections, pI)
        {
            intPoints[nInts++] = pI();
        }

        if (output)
        {
            // Check for point nearness
            meshOps::checkPointNearness(intPoints, T(1e-20));

            meshOps::writeVTK
            (
                this->mesh_,
                "ccSet_" + Foam::name(newIndex)
              + '<' + Foam::name(oldIndex) + '>',
                intPoints.size(),
                intPoints.size(),
                intPoints.size(),
                convert<scalar>(intPoints)
            );
        }

        return true;
    }

    // Check for point-segment intersections
    bool pointIntersections = false;

    forAll(oldFace, pointI)
    {
        label oldPoint = oldFace[pointI];

        if (commonPoints.found(oldPoint))
        {
            continue;
        }

        const Vector<T>& checkPoint = projections[pointI];

        // Loop through all new edges, and find possible intersections
        // with (projections of) old face points,
        forAll(newFace, pointJ)
        {
            edge newEdge = newFace.faceEdge(pointJ);

            if
            (
                meshOps::pointSegmentIntersection
                (
                    line<Vector<T>, const Vector<T>&>
                    (
                        convert<T>(newPoints[newEdge.start()]),
                        convert<T>(newPoints[newEdge.end()])
                    ),
                    checkPoint,
                    matchTol
                )
            )
            {
                commonPoints.insert
                (
                    oldPoint,
                    labelList(newEdge)
                );

                intersections.set(++nInts, checkPoint);

                pointIntersections = true;
            }
        }
    }

    forAll(newFace, pointI)
    {
        label newPoint = newFace[pointI];

        if (commonPoints.found(newPoint))
        {
            continue;
        }

        // Check if this point was modified during a collapse.
        // We might be able to avoid a calculation if it is.
        if (modPoints.found(newPoint))
        {
            // Fetch the two points that it originated from
            const labelList& mP = modPoints[newPoint];

            bool foundPoint = false;

            forAll(oldFace, pointJ)
            {
                const edge edgeToCheck = oldFace.faceEdge(pointJ);

                if (edgeToCheck == edge(mP[0], mP[1]))
                {
                    commonPoints.insert
                    (
                        newPoint,
                        labelList(edgeToCheck)
                    );

                    intersections.set
                    (
                        ++nInts,
                        convert<T>(newPoints[newPoint])
                    );

                    foundPoint = true;
                    pointIntersections = true;

                    break;
                }
            }

            if (foundPoint)
            {
                continue;
            }
        }

        const Vector<T> checkPoint = convert<T>(newPoints[newPoint]);

        forAll(oldFace, pointJ)
        {
            label nextJ = oldFace.fcIndex(pointJ);
            edge oldEdge = oldFace.faceEdge(pointJ);

            if
            (
                meshOps::pointSegmentIntersection
                (
                    line<Vector<T>, const Vector<T>&>
                    (
                        projections[pointJ],
                        projections[nextJ]
                    ),
                    checkPoint,
                    matchTol
                )
            )
            {
                commonPoints.insert
                (
                    newPoint,
                    labelList(oldEdge)
                );

                intersections.set(++nInts, checkPoint);

                pointIntersections = true;
            }
        }
    }

    if (pointIntersections && output)
    {
        Info << "Point Intersections exist: " << nl
             << " newFaceIndex: " << newIndex
             << " oldFaceIndex: " << oldIndex
             << endl;
    }

    if (oldFace.size() == 3 && newFace.size() == 3)
    {
        // Perform tests specific to triangular faces

        // Check whether any old projections are within
        // the new face. Count these as 'intersections'.
        forAll(oldFace, pointI)
        {
            label oldPoint = oldFace[pointI];

            if (commonPoints.found(oldPoint))
            {
                // Only skip for shared-points.
                // If the point-position was modified
                // due to a collapse, then this point
                // could be inside the new face.
                if (commonPoints[oldPoint].empty())
                {
                    continue;
                }
                else
                {
                    // Fetch master objects
                    const labelList& mObj = commonPoints[oldPoint];

                    // Check if the new face
                    // contains all master points
                    bool allMaster = true;

                    forAll(mObj, pointJ)
                    {
                        if (findIndex(newFace, mObj[pointJ]) == -1)
                        {
                            allMaster = false;
                            break;
                        }
                    }

                    if (allMaster)
                    {
                        continue;
                    }
                }
            }

            const Vector<T>& checkPoint = projections[pointI];

            if
            (
                meshOps::pointInTriFace
                (
                    triangle<Vector<T>, const Vector<T>&>
                    (
                        convert<T>(newPoints[newFace[0]]),
                        convert<T>(newPoints[newFace[1]]),
                        convert<T>(newPoints[newFace[2]])
                    ),
                    checkPoint,
                    matchTol,
                    false
                )
            )
            {
                intersections.set(++nInts, checkPoint);
            }
        }

        // Check whether and new points are within
        // projected old faces. Count these as 'intersections'.
        forAll(newFace, pointI)
        {
            label newPoint = newFace[pointI];

            if (commonPoints.found(newPoint))
            {
                continue;
            }

            const Vector<T> checkPoint = convert<T>(newPoints[newPoint]);

            if
            (
                meshOps::pointInTriFace
                (
                    triangle<Vector<T>, const Vector<T>&>
                    (
                        projections[0],
                        projections[1],
                        projections[2]
                    ),
                    checkPoint,
                    matchTol,
                    false
                )
            )
            {
                intersections.set(++nInts, checkPoint);
            }
        }

        // Loop through all new edges, and find possible intersections
        // with (projections of) old face edges,
        forAll(newFace, pointI)
        {
            edge newEdge = newFace.faceEdge(pointI);
            label nextLabel = newFace.nextLabel(pointI);

            forAll(oldFace, pointJ)
            {
                label nextJ = oldFace.fcIndex(pointJ);
                edge oldEdge = oldFace.faceEdge(pointJ);

                // Form an edge-pair
                Pair<edge> edgePair(oldEdge, newEdge);

                bool disableCheck = false;

                // Check edges topologically
                if (edgePair.first() == edgePair.second())
                {
                    const edge& checkEdge = edgePair.first();

                    // Check if points were modified by a collapse.
                    // If both were modified, continue with check.
                    if
                    (
                        !modPoints.found(checkEdge.start()) &&
                        !modPoints.found(checkEdge.end())
                    )
                    {
                        disableCheck = true;
                    }

                    // Skip shared points
                    if (commonPoints.found(checkEdge.start()))
                    {
                        if (commonPoints[checkEdge.start()].empty())
                        {
                            disableCheck = true;
                        }
                    }

                    if (commonPoints.found(checkEdge.end()))
                    {
                        if (commonPoints[checkEdge.end()].empty())
                        {
                            disableCheck = true;
                        }
                    }
                }
                else
                {
                    // Check for common vertices
                    label cV = edgePair.first().commonVertex(edgePair.second());

                    if (cV > -1)
                    {
                        // If this point was modified by a collapse
                        // to an edge mid-point, it can't be a common point.
                        // So, allow the check to continue.
                        if (!modPoints.found(cV))
                        {
                            disableCheck = true;
                        }

                        // Skip shared points
                        if (commonPoints.found(cV))
                        {
                            if (commonPoints[cV].empty())
                            {
                                disableCheck = true;
                            }
                        }
                    }
                }

                if (disableCheck)
                {
                    continue;
                }

                // Also check for bisection / point-on-edge cases
                bool foundPointOnEdge = false;

                forAll(edgePair, indexI)
                {
                    const edge thisEdge = edgePair[indexI];

                    const edge otherEdge =
                    (
                        (thisEdge == edgePair.first()) ?
                        edgePair.second() : edgePair.first()
                    );

                    forAll(otherEdge, pointI)
                    {
                        label pIndex = otherEdge[pointI];

                        if (commonPoints.found(pIndex))
                        {
                            // Fetch masterObjects
                            const labelList& mObj = commonPoints[pIndex];

                            // Skip shared-points.
                            if (mObj.size())
                            {
                                // Check if the old edge
                                // contains all master points
                                bool allMaster = true;

                                forAll(mObj, pointJ)
                                {
                                    if (findIndex(thisEdge, mObj[pointJ]) == -1)
                                    {
                                        allMaster = false;
                                        break;
                                    }
                                }

                                if (allMaster)
                                {
                                    foundPointOnEdge = true;
                                }
                            }
                        }

                        if (foundPointOnEdge)
                        {
                            break;
                        }
                    }

                    if (foundPointOnEdge)
                    {
                        break;
                    }
                }

                if (foundPointOnEdge)
                {
                    continue;
                }

                bool foundIntersection = false;

                foundIntersection =
                (
                    meshOps::segmentSegmentIntersection
                    (
                        line<Vector<T>, const Vector<T>&>
                        (
                            projections[pointJ],
                            projections[nextJ]
                        ),
                        line<Vector<T>, const Vector<T>&>
                        (
                            convert<T>(newPoints[newFace[pointI]]),
                            convert<T>(newPoints[nextLabel])
                        ),
                        matchTol,
                        intPoint
                    )
                );

                if (foundIntersection)
                {
                    intersections.set(++nInts, intPoint);
                }
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "\n\n"
            "bool faceSetAlgorithm::faceIntersection\n"
            "(\n"
            "    const label,\n"
            "    const label,\n"
            "    const T&,\n"
            "    Field<Vector<T> >&\n"
            ") const\n"
        )
            << " Invalid face pair: " << nl
            << " Old face: " << oldIndex << "::" << oldFace << nl
            << " New face: " << newIndex << "::" << newFace << nl
            << abort(FatalError);
    }

    // Copy intersections
    intPoints.setSize(nInts, Vector<T>::zero);

    nInts = 0;

    forAllConstIter(typename Map<Vector<T> >, intersections, pI)
    {
        intPoints[nInts++] = pI();
    }

    // Check for concurrent points.
    if (output)
    {
        meshOps::checkPointNearness(intPoints, T(1e-20));

        meshOps::writeVTK
        (
            this->mesh_,
            "ccSet_" + Foam::name(newIndex)
          + '<' + Foam::name(oldIndex) + '>',
            intPoints.size(),
            intPoints.size(),
            intPoints.size(),
            convert<scalar>(intPoints)
        );
    }

    // Found a convex set of points.
    if (nInts >= 3)
    {
        return true;
    }

    // Does not intersect
    return false;
}


// Compute the area / centre of a polygon
// formed by a convex set of points.
template <class T>
void faceSetAlgorithm::convexSetArea
(
    const label newFaceIndex,
    const label oldFaceIndex,
    const Field<Vector<T> >& cvxSet,
    const Vector<T>& refNorm,
    T& fArea,
    Vector<T>& fCentre,
    bool output
) const
{
    // Reset inputs
    fArea = pTraits<T>::zero;
    fCentre = Vector<T>::zero;

    // Try the trivial case for a triangle.
    if (cvxSet.size() == 3)
    {
        const Vector<T>& a = cvxSet[0];
        const Vector<T>& b = cvxSet[1];
        const Vector<T>& c = cvxSet[2];

        fArea = mag(0.5 * ((b - a)^(c - a)));
        fCentre = (pTraits<T>::one / T(3.0)) * (a + b + c);

        if (output)
        {
            Info << " newFaceIndex: " << newFaceIndex
                 << " oldFaceIndex: " << oldFaceIndex << nl
                 << " Area: " << fArea << nl
                 << " Centre: " << fCentre << nl
                 << endl;
        }

        return;
    }

    // Track edges
    label nEdges = 0;
    edgeList testEdges(0);

    // Loop through all points, and build edges with every
    // other point in the set
    forAll(cvxSet, i)
    {
        forAll(cvxSet, j)
        {
            // Skip duplicates.
            if (j == i)
            {
                continue;
            }

            // Define the edge
            edge tmpEdge(i, j);

            // If this is an existing edge, skip it.
            bool foundExisting = false;

            forAll(testEdges, edgeI)
            {
                if (testEdges[edgeI] == tmpEdge)
                {
                    foundExisting = true;
                    break;
                }
            }

            if (foundExisting)
            {
                continue;
            }

            // Specify a tolerance for collinearity
            T tolerance(1e-14);

            // Compute the normal to this edge
            Vector<T> n;

            n = ((cvxSet[tmpEdge.end()] - cvxSet[tmpEdge.start()]) ^ refNorm);
            n /= mag(n) + VSMALL;

            label curEdgeSign = 0;
            bool foundInternalEdge = false;

            // Quick-reject test:
            //   Check all other points in the set,
            //   and decide if all points lie on one side.
            forAll(cvxSet, k)
            {
                // Skip duplicates.
                if (tmpEdge[0] == k || tmpEdge[1] == k)
                {
                    continue;
                }

                Vector<T> rfVec = (cvxSet[k] - cvxSet[i]);
                T dotProd = (rfVec/(mag(rfVec) + VSMALL)) & n;

                // Skip nearly collinear points.
                if (mag(dotProd) < tolerance)
                {
                    continue;
                }

                // Obtain the sign of this point.
                label eSign = Foam::sign(dotProd);

                // Update the current sign if necessary.
                if (curEdgeSign == 0)
                {
                    curEdgeSign = eSign;
                }
                else
                if (curEdgeSign != eSign)
                {
                    // Interior edge. Bail out.
                    foundInternalEdge = true;
                    break;
                }
            }

            if (foundInternalEdge)
            {
                continue;
            }

            // Looks like we found an edge on the boundary.
            // Check its sign to ensure that it points outward.
            if (curEdgeSign == 1)
            {
                n *= -1.0;
                tmpEdge = tmpEdge.reverseEdge();
            }

            // Add to the list of edges.
            testEdges.setSize(++nEdges, tmpEdge);
        }
    }

    // Sanity check - do points match edges?
    if (testEdges.size() != cvxSet.size())
    {
        WarningIn
        (
            "\n"
            "void meshOps::convexSetArea\n"
            "(\n"
            "    const label newFaceIndex,\n"
            "    const label oldFaceIndex,\n"
            "    const vectorField& cvxSet,\n"
            "    const vector& refNorm,\n"
            "    scalar& fArea,\n"
            "    vector& fCentre,\n"
            "    bool output\n"
            ") const"
        )
            << " Points do not match edges. " << nl
            << " newFaceIndex: " << newFaceIndex
            << " oldFaceIndex: " << oldFaceIndex
            << " nPoints: " << cvxSet.size() << nl
            << " nEdges: " << testEdges.size() << nl
            << " Edge list: " << testEdges << nl
            << " Set: " << cvxSet << nl
            << endl;
    }

    // Find an approximate face-centroid
    T sumA = 0.0;
    Vector<T> sumAc = Vector<T>::zero;
    Vector<T> xC = average(cvxSet);

    forAll(testEdges, edgeI)
    {
        const edge& e = testEdges[edgeI];

        Vector<T> c = cvxSet[e[0]] + cvxSet[e[1]] + xC;
        T a = mag((cvxSet[e[1]] - cvxSet[e[0]]) ^ (xC - cvxSet[e[0]]));

        sumA += a;
        sumAc += a*c;
    }

    fCentre = (1.0/3.0)*sumAc/(sumA + VSMALL);
    fArea = 0.5*sumA;

    if (output)
    {
        Info << " newFaceIndex: " << newFaceIndex
             << " oldFaceIndex: " << oldFaceIndex << nl
             << " Edges: " << testEdges << nl
             << " Area: " << fArea << nl
             << " Centre: " << fCentre << nl
             << endl;
    }
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
