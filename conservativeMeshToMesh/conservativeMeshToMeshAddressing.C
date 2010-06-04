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

Description
    Private member of conservativeMeshToMesh.

    Calculates mesh to mesh addressing pattern (for each cell from one mesh
    find the closest cells in the other mesh).

\*---------------------------------------------------------------------------*/

#include "conservativeMeshToMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void conservativeMeshToMesh::calcAddressingAndWeights
(
    const label cellStart,
    const label cellSize
)
{
    if (debug)
    {
        Info << "conservativeMeshToMesh::calculateIntersectionAddressing() : "
             << "calculating mesh-to-mesh cell addressing" << endl;
    }

    // Fetch nearest-cell addressing from meshToMesh
    const labelList& cAddr = meshToMesh::cellAddressing();

    for (label cellI = cellStart; cellI < (cellStart + cellSize); cellI++)
    {
        // Fetch the volume of the new cell
        scalar newVol = toMesh().cellVolumes()[cellI];

        // Obtain weighting factors for this cell.
        // Perform several attempts for robustness.
        bool consistent = false;
        scalar searchFactor = 1.0;

        // Fetch references
        labelList& parents = addressing_[cellI];
        scalarField& weights = weights_[cellI];
        vectorField& centres = centres_[cellI];

        for (label attempt = 0; attempt < 10; attempt++)
        {
            consistent =
            (
                computeCellWeights
                (
                    cellI,
                    newVol,
                    cAddr[cellI],
                    searchFactor,
                    parents,
                    weights,
                    centres
                )
            );

            if (consistent)
            {
                break;
            }
            else
            {
                // Expand the search radius and try again.
                searchFactor *= 1.2;
            }
        }

        if (!consistent)
        {
            // Write out for post-processing
            label uIdx = 0;
            labelList candid = cellParents(cellI, searchFactor, cAddr[cellI]);
            labelList unMatch(candid.size() - parents.size(), -1);

            forAll(candid, cI)
            {
                if (findIndex(parents, candid[cI]) == -1)
                {
                    unMatch[uIdx++] = candid[cI];
                }
            }

            writeVTK("nCell_" + Foam::name(cellI), toMesh(), cellI, 3);
            writeVTK("oCell_" + Foam::name(cellI), fromMesh(), candid, 3);
            writeVTK("mCell_" + Foam::name(cellI), fromMesh(), parents, 3);
            writeVTK("uCell_" + Foam::name(cellI), fromMesh(), unMatch, 3);

            FatalErrorIn
            (
                "conservativeMeshToMesh::calcIntersectionAddressing()"
            )
                << "Encountered non-conservative weighting factors." << nl
                << " Cell: " << cellI << nl
                << " New cell volume: " << newVol << nl
                << " Candidate parent: " << cAddr[cellI] << nl
                << " Parents: " << parents << nl
                << " Weights: " << weights << nl
                << " Sum(Weights): " << sum(weights/newVol) << nl
                << " Error: " << mag(1.0 - sum(weights/newVol))
                << abort(FatalError);
        }
    }
}


// Compute cell weighting factors for a particular cell
bool conservativeMeshToMesh::computeCellWeights
(
    const label newCellIndex,
    const scalar newCellVolume,
    const label oldCandidate,
    const scalar searchFactor,
    labelList& parents,
    scalarField& weights,
    vectorField& centres
) const
{
    // Obtain candidate parents for this cell
    labelList candidates =
    (
        cellParents
        (
            newCellIndex,
            searchFactor,
            oldCandidate
        )
    );

    // Track actual intersections
    label nIntersects = 0;

    // Compute intersection volumes with candidates
    scalarField intVolumes(candidates.size(), 0.0);
    vectorField intCentres(candidates.size(), vector::zero);

    forAll(candidates, indexI)
    {
        bool intersects =
        (
            cellIntersection
            (
                newCellIndex,
                candidates[indexI],
                intVolumes[indexI],
                intCentres[indexI]
            )
        );

        if (intersects)
        {
            nIntersects++;
        }
    }

    // Now copy only valid intersections.
    parents.setSize(nIntersects, -1);
    weights.setSize(nIntersects, 0.0);
    centres.setSize(nIntersects, vector::zero);

    // Reset counter
    nIntersects = 0;

    forAll(intVolumes, indexI)
    {
        if (intVolumes[indexI] > 0.0)
        {
            parents[nIntersects] = candidates[indexI];
            weights[nIntersects] = intVolumes[indexI];
            centres[nIntersects] = intCentres[indexI];

            nIntersects++;
        }
    }

    // Test weights for consistency
    if (mag(1.0 - (sum(intVolumes)/newCellVolume)) > (1e-2*newCellVolume))
    {
        return false;
    }

    return true;
}


// Obtain a list of possible parent cells from the old mesh.
labelList conservativeMeshToMesh::cellParents
(
    const label newCellIndex,
    const scalar searchFactor,
    const label oldCandidate
) const
{
    labelHashSet masterCells, finalCells;

    // Fetch connectivity from the old mesh.
    const cellList& fromCells = fromMesh().cells();
    const labelList& fromOwner = fromMesh().faceOwner();
    const labelList& fromNeighbour = fromMesh().faceNeighbour();

    // Insert the old candidate first
    masterCells.insert(oldCandidate);

    for (label attempt = 0; attempt < 10; attempt++)
    {
        // Fetch the initial set of candidates
        labelList initList = masterCells.toc();

        // Accumulate a larger stencil of cell neighbours
        forAll(initList, indexI)
        {
            const cell& cellToCheck = fromCells[initList[indexI]];

            forAll(cellToCheck, faceI)
            {
                // Add owner to the list.
                masterCells.set
                (
                    fromOwner[cellToCheck[faceI]],
                    empty()
                );

                if (fromMesh().isInternalFace(cellToCheck[faceI]))
                {
                    // Add the neighbour.
                    masterCells.set
                    (
                        fromNeighbour[cellToCheck[faceI]],
                        empty()
                    );
                }
            }
        }
    }

    // Fetch the new cell, and determine its bounds.
    const cell& toCell = toMesh().cells()[newCellIndex];

    // Prepare an axis-aligned bounding box around the cell,
    // and add cells according to cell-centre positions.
    boundBox cellBox(toCell.points(toMesh().faces(), toMesh().points()), false);

    // Define a search radius
    vector bC = cellBox.midpoint();
    vector bMax = searchFactor * (cellBox.max() - bC);

    const vectorField& cellCentres = fromMesh().cellCentres();

    forAllIter(labelHashSet, masterCells, mIter)
    {
        vector xC = (cellCentres[mIter.key()] - bC);

        if ((xC & xC) < (bMax & bMax))
        {
            finalCells.insert(mIter.key());
        }
    }

    if (debug)
    {
        Info << " Cell: " << newCellIndex
             << " No. of parent candidates: "
             << finalCells.size()
             << " searchFactor: "
             << searchFactor
             << endl;
    }

    return finalCells.toc();
}


// Return the intersection volume between cells in old/new meshes
bool conservativeMeshToMesh::cellIntersection
(
    const label newCellIndex,
    const label oldCellIndex,
    scalar& intVolume,
    vector& intCentre
) const
{
    scalar tolFactor = 1e-8;

    // Reset inputs
    intVolume = 0.0;
    intCentre = vector::zero;

    // Fetch references for each mesh
    const edgeList& fromEdges = fromMesh().edges();
    const faceList& fromFaces = fromMesh().faces();
    const pointField& fromPoints = fromMesh().points();
    const cell& fromCell = fromMesh().cells()[oldCellIndex];
    const labelList& fromCellEdges = fromMesh().cellEdges()[oldCellIndex];
    const labelList& fromCellPoints = fromMesh().cellPoints()[oldCellIndex];

    const edgeList& toEdges = toMesh().edges();
    const faceList& toFaces = toMesh().faces();
    const pointField& toPoints = toMesh().points();
    const cell& toCell = toMesh().cells()[newCellIndex];
    const labelList& toCellEdges = toMesh().cellEdges()[newCellIndex];
    const labelList& toCellPoints = toMesh().cellPoints()[newCellIndex];

    // Track all possible intersections from here on.
    label nInts = 0;
    vectorField tP(0);
    vector intPoint = vector::zero;
    FixedList<vector,2> segment(vector::zero);

    // Check if any points are coincident.
    Map<label> FtoT, TtoF;

    forAll(fromCellPoints, pointI)
    {
        forAll(toCellPoints, pointJ)
        {
            if
            (
                magSqr
                (
                    fromPoints[fromCellPoints[pointI]]
                  - toPoints[toCellPoints[pointJ]]
                ) < tolFactor
            )
            {
                FtoT.insert
                (
                    fromCellPoints[pointI],
                    toCellPoints[pointJ]
                );

                TtoF.insert
                (
                    toCellPoints[pointJ],
                    fromCellPoints[pointI]
                );

                tP.setSize(++nInts, toPoints[toCellPoints[pointJ]]);
            }
        }
    }

    // If all points are common, this is identical to the old cell.
    if (FtoT.size() == fromCellPoints.size())
    {
        convexSetVolume
        (
            newCellIndex,
            oldCellIndex,
            tolFactor,
            tP,
            intVolume,
            intCentre
        );

        return true;
    }

    // Check whether any old points are within
    // the new cell. Count these as 'intersections'.
    forAll(fromCellPoints, pointI)
    {
        if (FtoT.found(fromCellPoints[pointI]))
        {
            continue;
        }

        const point& checkPoint = fromMesh().points()[fromCellPoints[pointI]];

        if
        (
            pointInCell
            (
                newCellIndex,
                checkPoint,
                false
            )
        )
        {
            tP.setSize(++nInts, checkPoint);
        }
    }

    // Check whether any new points are within
    // the old cell. Count these as 'intersections'.
    forAll(toCellPoints, pointI)
    {
        if (TtoF.found(toCellPoints[pointI]))
        {
            continue;
        }

        const point& checkPoint = toMesh().points()[toCellPoints[pointI]];

        if
        (
            pointInCell
            (
                oldCellIndex,
                checkPoint,
                true
            )
        )
        {
            tP.setSize(++nInts, checkPoint);
        }
    }

    bool foundIntersection = false;

    // Loop through all old edges, and find possible
    // intersections with faces of the new cell.
    forAll(fromCellEdges, edgeI)
    {
        const edge& edgeToCheck = fromEdges[fromCellEdges[edgeI]];

        forAll(toCell, faceI)
        {
            // Avoid common points, since this implies that
            // the edge intersects at a face point
            const face& faceToCheck = toFaces[toCell[faceI]];

            bool foundCommon = false;

            forAllConstIter(Map<label>, TtoF, pIter)
            {
                if (faceToCheck.which(pIter.key()) > -1)
                {
                    if
                    (
                        (edgeToCheck[0] == pIter()) ||
                        (edgeToCheck[1] == pIter())
                    )
                    {
                        foundCommon = true;
                        break;
                    }
                }
            }

            if (foundCommon)
            {
                continue;
            }

            segment[0] = fromPoints[edgeToCheck[0]];
            segment[1] = fromPoints[edgeToCheck[1]];

            foundIntersection = false;

            foundIntersection =
            (
                segmentFaceIntersection
                (
                    segment,
                    edgeToCheck,
                    toCell[faceI],
                    tolFactor,
                    intPoint,
                    false
                )
            );

            if (foundIntersection)
            {
                // Add to the list.
                tP.setSize(++nInts, intPoint);
            }
        }
    }

    // Loop through all new edges, and find possible
    // intersections with faces of the old cell.
    forAll(toCellEdges, edgeI)
    {
        const edge& edgeToCheck = toEdges[toCellEdges[edgeI]];

        forAll(fromCell, faceI)
        {
            // Avoid common points, since this implies that
            // the edge intersects at a face point
            const face& faceToCheck = fromFaces[fromCell[faceI]];

            bool foundCommon = false;

            forAllConstIter(Map<label>, FtoT, pIter)
            {
                if (faceToCheck.which(pIter.key()) > -1)
                {
                    if
                    (
                        (edgeToCheck[0] == pIter()) ||
                        (edgeToCheck[1] == pIter())
                    )
                    {
                        foundCommon = true;
                        break;
                    }
                }
            }

            if (foundCommon)
            {
                continue;
            }

            segment[0] = toPoints[edgeToCheck[0]];
            segment[1] = toPoints[edgeToCheck[1]];

            foundIntersection = false;

            foundIntersection =
            (
                segmentFaceIntersection
                (
                    segment,
                    edgeToCheck,
                    fromCell[faceI],
                    tolFactor,
                    intPoint,
                    true
                )
            );

            if (foundIntersection)
            {
                // Add to the list.
                tP.setSize(++nInts, intPoint);
            }
        }
    }

    // Found a polyhedral intersecting volume.
    // Compute the volume from points and return.
    if (nInts >= 4)
    {
        convexSetVolume
        (
            newCellIndex,
            oldCellIndex,
            tolFactor,
            tP,
            intVolume,
            intCentre
        );

        return true;
    }

    // Does not intersect
    return false;
}


// Determine whether the particular point lies
// inside the given cell
bool conservativeMeshToMesh::pointInCell
(
    const label cellIndex,
    const point& checkPoint,
    const bool useFromMesh
) const
{
    const cell* cellPtr = NULL;

    if (useFromMesh)
    {
        cellPtr = &(fromMesh().cells()[cellIndex]);
    }
    else
    {
        cellPtr = &(toMesh().cells()[cellIndex]);
    }

    // Fetch the reference
    const cell& cellToCheck = *cellPtr;

    label own = -1;
    vector xf(vector::zero), nf(vector::zero);

    forAll(cellToCheck, faceI)
    {
        if (useFromMesh)
        {
            own = fromMesh().faceOwner()[cellToCheck[faceI]];
            xf = fromMesh().faceCentres()[cellToCheck[faceI]];
            nf = fromMesh().faceAreas()[cellToCheck[faceI]];
        }
        else
        {
            own = toMesh().faceOwner()[cellToCheck[faceI]];
            xf = toMesh().faceCentres()[cellToCheck[faceI]];
            nf = toMesh().faceAreas()[cellToCheck[faceI]];
        }

        if (((xf - checkPoint) & nf) > 0.0)
        {
            if (own != cellIndex)
            {
                return false;
            }
        }
        else
        {
            if (own == cellIndex)
            {
                return false;
            }
        }
    }

    // Passed test with all faces
    return true;
}


// Determine whether a particular line segment
// intersects with a given face
bool conservativeMeshToMesh::segmentFaceIntersection
(
    const FixedList<vector,2>& segment,
    const edge& segmentIndices,
    const label faceIndex,
    const scalar tolFactor,
    point& intPoint,
    const bool useFromMesh
) const
{
    vector nf(vector::zero), fp(vector::zero);

    if (useFromMesh)
    {
        nf = fromMesh().faceAreas()[faceIndex];
        fp = fromMesh().points()[fromMesh().faces()[faceIndex][0]];
    }
    else
    {
        nf = toMesh().faceAreas()[faceIndex];
        fp = toMesh().points()[toMesh().faces()[faceIndex][0]];
    }

    // Normalize
    nf /= mag(nf) + VSMALL;

    // Compute uValue
    scalar numerator = nf & (fp - segment[0]);
    scalar denominator = nf & (segment[1] - segment[0]);

    // Specify a tolerance
    scalar tolerance = tolFactor*mag(segment[1] - segment[0]);

    // Check if the edge is parallel to the face
    if (mag(denominator) < tolerance)
    {
        return false;
    }

    scalar u = (numerator / denominator);

    // Check for intersection along line.
    if ((u >= 0.0) && (u <= 1.0))
    {
        // Compute point of intersection
        intPoint = segment[0] + u*(segment[1] - segment[0]);

        // Also make sure that intPoint lies within triFace
        if (pointInFace(faceIndex, intPoint, useFromMesh))
        {
            return true;
        }
    }

    // Failed to fall within edge-bounds, or within face
    return false;
}


// Determine whether the particular point lies
// inside the given face
bool conservativeMeshToMesh::pointInFace
(
    const label faceIndex,
    const point& checkPoint,
    const bool useFromMesh
) const
{
    FixedList<vector, 2> segment(vector::zero);

    vector nf(vector::zero);
    const face* facePtr = NULL;

    if (useFromMesh)
    {
        nf = fromMesh().faceAreas()[faceIndex];
        facePtr = &(fromMesh().faces()[faceIndex]);
    }
    else
    {
        nf = toMesh().faceAreas()[faceIndex];
        facePtr = &(toMesh().faces()[faceIndex]);
    }

    // Normalize
    nf /= mag(nf) + VSMALL;

    const face& faceToCheck = *facePtr;

    forAll(faceToCheck, pI)
    {
        vector p, q;

        label nI = faceToCheck.fcIndex(pI);

        if (useFromMesh)
        {
            segment[0] = fromMesh().points()[faceToCheck[pI]];
            segment[1] = fromMesh().points()[faceToCheck[nI]];
        }
        else
        {
            segment[0] = toMesh().points()[faceToCheck[pI]];
            segment[1] = toMesh().points()[faceToCheck[nI]];
        }

        // Identify vectors
        p = (segment[1] - segment[0]);
        q = (checkPoint - segment[1]);

        scalar area = 0.5 * ((p ^ q) & nf);

        if (area < 0.0)
        {
            return false;
        }
    }

    // Passed test with all edges
    return true;
}


// Compute the volume of a polyhedron
// formed by a convex set of points.
void conservativeMeshToMesh::convexSetVolume
(
    const label newCellIndex,
    const label oldCellIndex,
    const scalar tolFraction,
    const vectorField& cvxSet,
    scalar& cVolume,
    vector& cCentre
) const
{
    // Reset inputs
    cVolume = 0.0;
    cCentre = vector::zero;

    // Try the trivial case for a tetrahedron.
    // No checking for orientation here.
    if (cvxSet.size() == 4)
    {
        cCentre = average(cvxSet);

        cVolume =
        (
            mag
            (
                (1.0/6.0) *
                (
                    ((cvxSet[1] - cvxSet[0]) ^ (cvxSet[2] - cvxSet[0])) &
                     (cvxSet[3] - cvxSet[0])
                )
            )
        );

        return;
    }

    // Track faces
    face tmpFace(3);
    DynamicList<face> testFaces(10);

    // Loop through all points, and build faces with every
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

            forAll(cvxSet, k)
            {
                // Skip duplicates.
                if (k == i || k == j)
                {
                    continue;
                }

                // Configure the face.
                tmpFace[0] = i;
                tmpFace[1] = j;
                tmpFace[2] = k;

                // Specify a tolerance
                scalar tolerance = tolFraction*magSqr(cvxSet[j] - cvxSet[i]);

                // Compute the normal to this face
                vector n = tmpFace.normal(cvxSet);

                n /= mag(n) + VSMALL;

                label curFaceSign = 0;
                bool foundInternalFace = false, foundCoPlanar = false;

                // Quick-reject test:
                //   Check all other points in the set,
                //   and decide if all points lie on one side.
                forAll(cvxSet, l)
                {
                    // Skip duplicates.
                    if (tmpFace[0] == l || tmpFace[1] == l || tmpFace[2] == l)
                    {
                        continue;
                    }

                    vector rfVec = (cvxSet[l] - cvxSet[i]);
                    scalar dotProd = (rfVec & n);

                    // Skip co-planar points for now,
                    // but keep note for later stage
                    if (mag(dotProd) < tolerance)
                    {
                        foundCoPlanar = true;
                        continue;
                    }

                    // Obtain the sign of this point.
                    label fSign = Foam::sign(dotProd);

                    // Update the current sign if necessary.
                    if (curFaceSign == 0)
                    {
                        curFaceSign = fSign;
                    }
                    else
                    if (curFaceSign != fSign)
                    {
                        // Interior face. Bail out.
                        foundInternalFace = true;
                        break;
                    }
                }

                if (foundInternalFace)
                {
                    continue;
                }

                // Include all other co-planar points,
                // if any were found in the quick-reject test
                if (foundCoPlanar)
                {
                    forAll(cvxSet, l)
                    {
                        // Skip duplicates.
                        if (findIndex(tmpFace, l) > -1)
                        {
                            continue;
                        }

                        vector rfVec = (cvxSet[l] - cvxSet[i]);
                        scalar dotProd = (rfVec & n);

                        if (mag(dotProd) < tolerance)
                        {
                            // Need to configure a new face.
                            face newFace(3);
                            bool foundLocation = false;

                            forAll(tmpFace, pI)
                            {
                                label nI = tmpFace.fcIndex(pI);

                                newFace[0] = tmpFace[pI];
                                newFace[1] = l;
                                newFace[2] = tmpFace[nI];

                                // Compute the normal.
                                vector nNew = newFace.normal(cvxSet);

                                if ((n & nNew) > 0.0)
                                {
                                    // Insert the point.
                                    insertLabel
                                    (
                                        l,
                                        tmpFace[pI],
                                        tmpFace[nI],
                                        tmpFace
                                    );

                                    foundLocation = true;

                                    break;
                                }
                            }

                            if (!foundLocation)
                            {
                                FatalErrorIn
                                (
                                    "conservativeMeshToMesh::convexSetVolume()"
                                )   << "Cannot find appropriate configuration."
                                    << nl << " New: " << newCellIndex
                                    << nl << " Old: " << oldCellIndex
                                    << nl << " Face: " << tmpFace.points(cvxSet)
                                    << nl << " with point: " << cvxSet[l]
                                    << nl << " Set: " << cvxSet
                                    << abort(FatalError);
                            }
                        }
                    }
                }

                // Looks like we found a face on the boundary.
                // Check its sign to ensure that it points outward.
                if (curFaceSign == 1)
                {
                    tmpFace = tmpFace.reverseFace();
                }

                // Ensure that the face wasn't checked in.
                bool alreadyCheckedIn = false;

                forAll(testFaces, faceI)
                {
                    const face& checkFace = testFaces[faceI];

                    if (face::compare(tmpFace, checkFace))
                    {
                        alreadyCheckedIn = true;
                        break;
                    }
                }

                // Add this face to the list of faces.
                if (!alreadyCheckedIn)
                {
                    testFaces.append(tmpFace);
                }

                // Reset the face size.
                tmpFace.setSize(3, -1);
            }
        }
    }

    // Weed-out faces that are sub-sets of larger faces
    label nFaces = 0;
    faceList cellFaces(testFaces.size());

    forAll(testFaces, faceI)
    {
        bool subset = false;

        const face& checkFace = testFaces[faceI];

        forAll(testFaces, faceJ)
        {
            const face& testFace = testFaces[faceJ];

            if (testFace.size() > checkFace.size())
            {
                bool foundUniquePoint = false;

                forAll(checkFace, pI)
                {
                    if (findIndex(testFace, checkFace[pI]) == -1)
                    {
                        foundUniquePoint = true;
                        break;
                    }
                }

                if (!foundUniquePoint)
                {
                    subset = true;
                    break;
                }
            }
        }

        if (!subset)
        {
            cellFaces[nFaces++] = testFaces[faceI];
        }
    }

    // Set to actual size
    cellFaces.setSize(nFaces);

    // Find an approximate cell-centroid
    vector xC = average(cvxSet);

    // Calculate volume from all accumulated faces.
    forAll(cellFaces, faceI)
    {
        const face& checkFace = cellFaces[faceI];

        vector xF = checkFace.centre(cvxSet);
        vector Sf = checkFace.normal(cvxSet);

        // Calculate 3*face-pyramid volume
        scalar pyr3Vol = Sf & (xF - xC);

        // Calculate face-pyramid centre
        vector pc = (3.0/4.0)*xF + (1.0/4.0)*xC;

        // Accumulate volume-weighted face-pyramid centre
        cCentre += pyr3Vol*pc;

        // Accumulate face-pyramid volume
        cVolume += pyr3Vol;
    }

    cCentre /= cVolume + VSMALL;
    cVolume *= (1.0/3.0);

    if (debug)
    {
        // Write out points for post-processing
        labelListList cpList(cvxSet.size(), labelList(1));

        forAll(cpList, i)
        {
            cpList[i][0] = i;
        }

        // For post-processing purposes, define a name
        word cvxSetName
        (
            "cvxSet_"
          + Foam::name(newCellIndex)
          + '_'
          + Foam::name(oldCellIndex)
        );

        writeVTK
        (
            cvxSetName,
            cvxSet.size(),
            cvxSet.size(),
            cvxSet.size(),
            cvxSet,
            cpList,
            0
        );

        Info << nl
             << " Convex set: " << cvxSetName << nl
             << " cellFaces: " << cellFaces
             << " Volume: " << cVolume << nl
             << endl;
    }
}


// Method to insert a label between two labels in a face
// Assumes that all labels are unique.
void conservativeMeshToMesh::insertLabel
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
        FatalErrorIn("label conservativeMeshToMesh::insertLabel()")
            << "Cannot insert " << newLabel << " in list: " << list << endl
            << " Labels: "
            << labelA << " and " << labelB << " were not found in sequence."
            << abort(FatalError);
    }

    // Transfer the list
    list.transfer(newList);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
