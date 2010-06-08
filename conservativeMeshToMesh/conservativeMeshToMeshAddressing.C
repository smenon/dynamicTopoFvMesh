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
#include "IOmanip.H"
#include "SortableList.H"

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
        // Fetch references
        labelList& parents = addressing_[cellI];
        scalarField& weights = weights_[cellI];
        vectorField& centres = centres_[cellI];

        // Obtain weighting factors for this cell.
        computeCellWeights
        (
            cellI,
            cAddr[cellI],
            parents,
            weights,
            centres
        );
    }
}


// Compute cell weighting factors for a particular cell
void conservativeMeshToMesh::computeCellWeights
(
    const label newCellIndex,
    const label oldCandidate,
    labelList& parents,
    scalarField& weights,
    vectorField& centres
) const
{
    scalar searchFactor = 1.0;

    label nOldIntersects = 0, nIntersects = 0, nAttempts = 0;

    // Maintain a list of candidates and intersection points
    labelList candidates;
    List<vectorField> tP;

    // Fetch the volume of the new cell
    scalar newCellVolume = toMesh().cellVolumes()[newCellIndex];

    while (nAttempts < 10)
    {
        // Obtain candidate parents for this cell
        candidates =
        (
            cellParents
            (
                newCellIndex,
                searchFactor,
                oldCandidate
            )
        );

        // Set sizes
        boolList intersects(candidates.size(), false);
        tP.setSize(candidates.size());

        // Test for intersections
        forAll(candidates, indexI)
        {
            intersects[indexI] =
            (
                cellIntersection
                (
                    newCellIndex,
                    candidates[indexI],
                    tP[indexI]
                )
            );

            if (intersects[indexI])
            {
                nIntersects++;
            }
        }

        if (nIntersects == nOldIntersects)
        {
            // Set sizes
            parents.setSize(nIntersects, -1);
            weights.setSize(nIntersects, 0.0);
            centres.setSize(nIntersects, vector::zero);

            // Reset counter
            nIntersects = 0;

            forAll(intersects, indexI)
            {
                if (intersects[indexI])
                {
                    parents[nIntersects] = candidates[indexI];

                    // Compute weights
                    convexSetVolume
                    (
                        newCellIndex,
                        parents[nIntersects],
                        tP[indexI],
                        weights[nIntersects],
                        centres[nIntersects]
                    );

                    nIntersects++;
                }
            }

            break;
        }
        else
        {
            nAttempts++;
            nOldIntersects = nIntersects;
            nIntersects = 0;

            // Expand the search radius and try again.
            searchFactor *= 1.6;
        }
    }

    // Test weights for consistency
    if (mag(newCellVolume - sum(weights)) > 1e-15)
    {
        // Write out for post-processing
        label uIdx = 0, cellI = newCellIndex;
        labelList unMatch(candidates.size() - parents.size(), -1);

        forAll(candidates, cI)
        {
            if (findIndex(parents, candidates[cI]) == -1)
            {
                unMatch[uIdx++] = candidates[cI];
            }
        }

        writeVTK("nCell_" + Foam::name(cellI), toMesh(), cellI, 3);
        writeVTK("oCell_" + Foam::name(cellI), fromMesh(), candidates, 3);
        writeVTK("mCell_" + Foam::name(cellI), fromMesh(), parents, 3);
        writeVTK("uCell_" + Foam::name(cellI), fromMesh(), unMatch, 3);

        // Write out intersection points
        forAll(tP, indexI)
        {
            if (tP[indexI].size() >= 4)
            {
                writeVTK(cellI, candidates[indexI], tP[indexI]);

                // Write out to screen
                scalar dummyWeight = 0.0;
                vector dummyCentre = vector::zero;

                convexSetVolume
                (
                    cellI,
                    candidates[indexI],
                    tP[indexI],
                    dummyWeight,
                    dummyCentre,
                    true
                );
            }
        }

        // Sort list for easier post-processing
        SortableList<label> sortedParents(parents);

        // In-place reorder weights
        const labelList& indices = sortedParents.indices();

        forAll(indices, indexI)
        {
            Info << sortedParents[indexI] << ": "
                 << setprecision(16)
                 << weights[indices[indexI]]
                 << endl;
        }

        FatalErrorIn
        (
            "conservativeMeshToMesh::calcIntersectionAddressing()"
        )
            << "Encountered non-conservative weighting factors." << nl
            << " Cell: " << newCellIndex << nl
            << " Candidate parent: " << oldCandidate << nl
            << " nAttempts: " << nAttempts << nl
            << setprecision(16)
            << " New cell volume: " << newCellVolume << nl
            << " Sum(Weights): " << sum(weights) << nl
            << " Error: " << (newCellVolume - sum(weights)) << nl
            << " Norm Sum(Weights): " << sum(weights/newCellVolume) << nl
            << " Norm Error: " << mag(1.0 - sum(weights/newCellVolume))
            << abort(FatalError);
    }
}


// Obtain a list of possible parent cells from the old mesh.
labelList conservativeMeshToMesh::cellParents
(
    const label newCellIndex,
    const scalar searchFactor,
    const label oldCandidate
) const
{
    labelList finalCells;
    DynamicList<label> masterCells(100);

    // Fetch connectivity from the old mesh.
    const labelListList& fromCellCells = fromMesh().cellCells();

    // Insert the old candidate first
    masterCells.append(oldCandidate);

    for (label attempt = 0; attempt < 10; attempt++)
    {
        // Fetch the initial set of candidates
        DynamicList<label> initList(masterCells);

        // Accumulate a larger stencil of cell neighbours
        forAll(initList, indexI)
        {
            const labelList& cc = fromCellCells[initList[indexI]];

            forAll(cc, cellI)
            {
                if (findIndex(masterCells, cc[cellI]) == -1)
                {
                    masterCells.append(cc[cellI]);
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

    label nEntries = 0;

    // Count the number of entries
    forAll(masterCells, cellI)
    {
        vector xC = (cellCentres[masterCells[cellI]] - bC);

        if ((xC & xC) < (bMax & bMax))
        {
            nEntries++;
        }
    }

    // Set size and reset counter
    finalCells.setSize(nEntries, -1);
    nEntries = 0;

    forAll(masterCells, cellI)
    {
        vector xC = (cellCentres[masterCells[cellI]] - bC);

        if ((xC & xC) < (bMax & bMax))
        {
            finalCells[nEntries++] = masterCells[cellI];
        }
    }

    if (debug)
    {
        Info << " Cell: " << newCellIndex
             << " No. of parent candidates: "
             << nEntries
             << " searchFactor: "
             << searchFactor
             << endl;
    }

    return finalCells;
}


// Return the intersection volume between cells in old/new meshes
bool conservativeMeshToMesh::cellIntersection
(
    const label newCellIndex,
    const label oldCellIndex,
    vectorField& tP
) const
{
    scalar tolFactor = 1e-12;

    // Reset inputs
    tP.clear();

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

    // Found a convex set of points.
    if (nInts >= 4)
    {
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
    const vectorField& cvxSet,
    scalar& cVolume,
    vector& cCentre,
    bool output
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
    label nFaces = 0;
    faceList testFaces(0);
    labelHashSet uniquePts;

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

                // Quick-reject test:
                //   If this is a subset of an existing face, skip it.
                bool foundSubSet = false;

                forAll(testFaces, faceI)
                {
                    const face& checkFace = testFaces[faceI];

                    if (checkFace.size() >= tmpFace.size())
                    {
                        bool foundUniquePoint = false;

                        forAll(tmpFace, pI)
                        {
                            if (findIndex(checkFace, tmpFace[pI]) == -1)
                            {
                                foundUniquePoint = true;
                                break;
                            }
                        }

                        if (!foundUniquePoint)
                        {
                            foundSubSet = true;
                            break;
                        }
                    }
                }

                if (foundSubSet)
                {
                    continue;
                }

                // Specify a tolerance for planarity
                scalar tolerance = 1e-14;
                //tolFraction*magSqr(cvxSet[j] - cvxSet[i]);

                // Compute the normal to this face
                vector n = tmpFace.normal(cvxSet);

                n /= mag(n);

                label curFaceSign = 0;
                bool foundInternalFace = false;

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
                    scalar dotProd = (rfVec/mag(rfVec)) & n;

                    // Skip nearly co-planar points.
                    if (mag(dotProd) < tolerance)
                    {
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

                // Looks like we found a face on the boundary.
                // Check its sign to ensure that it points outward.
                if (curFaceSign == 1)
                {
                    n *= -1.0;
                    tmpFace = tmpFace.reverseFace();
                }

                // Ensure that the face wasn't checked in.
                bool alreadyCheckedIn = false;

                forAll(testFaces, faceI)
                {
                    // Fetch a non-const reference, since this face
                    // might be modified in this loop.
                    face& checkFace = testFaces[faceI];

                    label nCommon = 0;

                    uniquePts.clear();

                    forAll(tmpFace, pI)
                    {
                        if (findIndex(checkFace, tmpFace[pI]) > -1)
                        {
                            nCommon++;
                        }
                        else
                        {
                            uniquePts.insert(tmpFace[pI]);
                        }
                    }

                    if (nCommon >= 2)
                    {
                        if (checkFace.size() >= tmpFace.size())
                        {
                            // Check for unique points
                            if (uniquePts.size() > 0)
                            {
                                // Compute the existing normal
                                vector eNorm = checkFace.normal(cvxSet);

                                scalar dotProd = (n & (eNorm/mag(eNorm)));

                                if
                                (
                                    (mag(1.0 - dotProd) < tolerance) &&
                                    (dotProd > 0.0)
                                )
                                {
                                    // Add all unique points to checkFace
                                    insertPointLabels
                                    (
                                        n,
                                        cvxSet,
                                        uniquePts,
                                        checkFace
                                    );

                                    alreadyCheckedIn = true;
                                    break;
                                }
                            }
                            else
                            {
                                // Subset face
                                alreadyCheckedIn = true;
                                break;
                            }
                        }
                        else
                        {
                            // checkFace is a subset. Replace it.
                            checkFace = tmpFace;

                            alreadyCheckedIn = true;
                            break;
                        }
                    }
                }

                // Add this face to the list of faces.
                if (!alreadyCheckedIn)
                {
                    testFaces.setSize(++nFaces, tmpFace);
                }

                // Reset the face size.
                tmpFace.setSize(3, -1);
            }
        }
    }

    // Account for planarity test failure.
    //  - Check for subsets.
    forAll(testFaces, faceI)
    {
        // Fetch a non-const reference, since this face
        // might be modified in this loop.
        face& checkFace = testFaces[faceI];

        // Account for deleted testFaces
        if (checkFace.empty())
        {
            continue;
        }

        // Compute the normal to this face
        vector n = checkFace.normal(cvxSet);

        forAll(testFaces, faceJ)
        {
            if (faceI == faceJ)
            {
                continue;
            }

            // Fetch a non-const reference, since this face
            // might be modified in this loop.
            face& testFace = testFaces[faceJ];

            if (checkFace.size() >= testFace.size())
            {
                label nCommon = 0;

                uniquePts.clear();

                forAll(testFace, pI)
                {
                    if (findIndex(checkFace, testFace[pI]) > -1)
                    {
                        nCommon++;
                    }
                    else
                    {
                        uniquePts.insert(testFace[pI]);
                    }
                }

                if (nCommon >= 3)
                {
                    // Delete the test face
                    testFace.clear();

                    // Add all unique points to checkFace
                    // Failed the tolerance test before,
                    // so don't check for it now
                    if (uniquePts.size())
                    {
                        insertPointLabels
                        (
                            n,
                            cvxSet,
                            uniquePts,
                            checkFace
                        );
                    }
                }
            }
        }
    }

    // Find an approximate cell-centroid
    vector xC = average(cvxSet);

    // Calculate volume from all accumulated faces.
    forAll(testFaces, faceI)
    {
        const face& checkFace = testFaces[faceI];

        if (checkFace.empty())
        {
            continue;
        }

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

    if (output)
    {
        Info << " newCellIndex: " << newCellIndex
             << " oldCellIndex: " << oldCellIndex << nl
             << " Faces: " << testFaces << nl
             // << " Volume: " << cVolume << nl
             // << " Centre: " << cCentre << nl
             << endl;
    }
}


// Method to insert labels in a face, so that
// right-handedness is preserved.
void conservativeMeshToMesh::insertPointLabels
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
