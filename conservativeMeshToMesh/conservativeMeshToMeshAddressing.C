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
    Private members of conservativeMeshToMesh.

    Calculates mesh to mesh addressing pattern (for each cell from one mesh
    find the intersecting cells in the other mesh).

\*---------------------------------------------------------------------------*/

#include "conservativeMeshToMesh.H"

#include "triFace.H"
#include "IOmanip.H"
#include "ListOps.H"
#include "clockTime.H"
#include "tetPointRef.H"

#include "tetIntersection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void conservativeMeshToMesh::calcAddressingAndWeights
(
    const label cellStart,
    const label cellSize,
    bool report
)
{
    if (debug)
    {
        Info<< "conservativeMeshToMesh::calculateIntersectionAddressing() : "
            << "calculating mesh-to-mesh cell addressing" << endl;
    }

    // Fetch nearest-cell addressing
    const labelList& cAddr = conservativeMeshToMesh::cellAddressing();

    clockTime sTimer;
    bool reported = false;
    scalar maxError = 0.0;
    label count = 0, oldStart = cellStart, nInconsistencies = 0;
    scalar interval = 0.5, oIndex = 0.0, nIndex = 0.0, matchTol = 1e-4;

    oIndex = ::floor(sTimer.elapsedTime() / interval);

    for (label cellI = cellStart; cellI < (cellStart + cellSize); cellI++)
    {
        count++;

        // Update the index, if its changed
        nIndex = ::floor(sTimer.elapsedTime() / interval);

        if ((nIndex - oIndex) > VSMALL)
        {
            oIndex = nIndex;

            scalar percent = 100.0 * (double(count) / (cellSize + VSMALL));

            // Report progress
            if (report)
            {
                Info<< "  Progress: " << percent << "% : "
                    << "  Cells processed: " << count
                    << "  out of " << cellSize << " total."
                    << "             \r"
                    << flush;

                reported = true;
            }
        }

        // Fetch references
        labelList& parents = addressing_[cellI];
        scalarField& weights = weights_[cellI];
        vectorField& centres = centres_[cellI];

        // Obtain weighting factors for this cell.
        bool consistent =
        (
            computeWeights
            (
                cellI,
                cAddr[cellI],
                srcMesh().cellCells(),
                matchTol,
                parents,
                weights,
                centres
            )
        );

        if (!consistent)
        {
            maxError = Foam::max(maxError, mag(1.0 - sum(weights)));

            nInconsistencies++;
        }

        if ((cellI - oldStart) > 50)
        {
            ctrMutex_.lock();

            counter_ += (cellI - oldStart);

            ctrMutex_.unlock();

            // Reset start index
            oldStart = cellI;
        }
    }

    // Final addition to counter
    ctrMutex_.lock();

    counter_ += ((cellStart + cellSize) - oldStart);

    ctrMutex_.unlock();

    if (reported && report)
    {
        Info<< "  Progress: 100%"
            << "  Entities processed: " << count
            << "  out of " << cellSize << " total."
            << "             \r"
            << endl;
    }
}


// Invert addressing from source to target
bool conservativeMeshToMesh::invertAddressing()
{
    // Read in addressing arrays from the source mesh
    IOList<labelList> srcAddressing
    (
        IOobject
        (
            "addressing",
            srcMesh().time().timeName(),
            srcMesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Check for compatibility
    bool compatible = true;
    label targetCells = tgtMesh().nCells();
    labelList nCellsPerCell(targetCells, 0);

    forAll(srcAddressing, cellI)
    {
        const labelList& srcAddr = srcAddressing[cellI];

        label maxIndex = findMax(srcAddr);

        if (srcAddressing[cellI][maxIndex] >= targetCells)
        {
            Info<< " Addressing is not compatible. " << endl;

            compatible = false;
            break;
        }

        forAll(srcAddr, j)
        {
            nCellsPerCell[srcAddr[j]]++;
        }
    }

    if (!compatible)
    {
        addressing_.clear();
        addressing_.setSize(targetCells);

        return false;
    }

    // Read weights
    IOList<scalarField> srcWeights
    (
        IOobject
        (
            "weights",
            srcMesh().time().timeName(),
            srcMesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read centres
    IOList<vectorField> srcCentres
    (
        IOobject
        (
            "centres",
            srcMesh().time().timeName(),
            srcMesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Set sizes
    forAll(nCellsPerCell, cellI)
    {
        addressing_[cellI].setSize(nCellsPerCell[cellI]);
        weights_[cellI].setSize(nCellsPerCell[cellI]);
        centres_[cellI].setSize(nCellsPerCell[cellI]);
    }

    nCellsPerCell = 0;

    // Invert addressing
    forAll(srcAddressing, cellI)
    {
        const labelList& srcAddr = srcAddressing[cellI];
        const scalarField& srcVl = srcWeights[cellI];
        const vectorField& srcCt = srcCentres[cellI];

        forAll(srcAddr, j)
        {
            label cellJ = srcAddr[j];

            addressing_[cellJ][nCellsPerCell[cellJ]] = cellI;
            weights_[cellJ][nCellsPerCell[cellJ]] = srcVl[j];
            centres_[cellJ][nCellsPerCell[cellJ]] = srcCt[j];

            nCellsPerCell[cellJ]++;
        }
    }

    // Check weights for consistency
    forAll(weights_, cellI)
    {
        if (mag(1.0 - sum(weights_[cellI])) > 5e-14)
        {
            Info<< " Weights are not compatible. " << nl
                << " Cell: " << cellI << nl
                << " Sum(weights): " << sum(weights_[cellI]) << nl
                << " Error: " << mag(1.0 - sum(weights_[cellI]))
                << endl;

            compatible = false;
            break;
        }
    }

    if (!compatible)
    {
        addressing_.clear();
        weights_.clear();
        centres_.clear();

        addressing_.setSize(targetCells);
        weights_.setSize(targetCells);
        centres_.setSize(targetCells);

        return false;
    }

    // Return a successful inversion
    return true;
}


// Decompose the input cell using face-centre
void conservativeMeshToMesh::decomposeCell
(
    const cell& polyCell,
    const point& cCentre,
    const faceList& faces,
    const pointField& meshPoints,
    const pointField& faceCentres,
    DynamicList<TetPoints>& decompTets
)
{
    TetPoints tPoints;

    decompTets.clear();

    forAll(polyCell, faceI)
    {
        label faceIndex = polyCell[faceI];

        const face& polyFace = faces[faceIndex];
        const vector& fCentre = faceCentres[faceIndex];

        if (polyFace.size() == 3)
        {
            tPoints[0] = meshPoints[polyFace[0]];
            tPoints[1] = meshPoints[polyFace[1]];
            tPoints[2] = meshPoints[polyFace[2]];
            tPoints[3] = cCentre;

            // Add a new entry
            decompTets.append(tPoints);
        }
        else
        {
            // Point ordering is irrelevant, since the
            // clipping algorithm is expected to work
            // regardless of orientation
            forAll(polyFace, pI)
            {
                tPoints[0] = meshPoints[polyFace[pI]];
                tPoints[1] = meshPoints[polyFace.nextLabel(pI)];
                tPoints[2] = fCentre;
                tPoints[3] = cCentre;

                // Add a new entry
                decompTets.append(tPoints);
            }
        }
    }
}


// Compute weighting factors for a particular cell
bool conservativeMeshToMesh::computeWeights
(
    const label index,
    const label srcCandidate,
    const labelListList& srcNeighbourList,
    const scalar mTol,
    labelList& parents,
    scalarField& weights,
    vectorField& centres
) const
{
    if (parents.size() || weights.size() || centres.size())
    {
        FatalErrorIn
        (
            "\n\n"
            "bool conservativeMeshToMesh::computeWeights\n"
            "(\n"
            "    const label index,\n"
            "    const label srcCandidate,\n"
            "    const labelListList& srcNeighbourList,\n"
            "    const scalar mTol,\n"
            "    labelList& parents,\n"
            "    scalarField& weights,\n"
            "    vectorField& centres\n"
            ") const\n"
        )
            << " Addressing has already been calculated." << nl
            << " Index: " << index << nl
            << " srcCandidate: " << srcCandidate << nl
            << " oldNeighborList: " << srcNeighbourList << nl
            << " Parents: " << parents << nl
            << " Weights: " << weights << nl
            << " Centres: " << centres << nl
            << abort(FatalError);
    }

    bool changed;
    label nAttempts = 0, mapCandidate = -1;

    const faceList& tgtFaces = tgtMesh().faces();
    const cellList& tgtCells = tgtMesh().cells();
    const pointField& tgtPoints = tgtMesh().points();
    const vectorField& tgtFaceCentres = tgtMesh().faceCentres();
    const vectorField& tgtCellCentres = tgtMesh().cellCentres();
    const scalarField& tgtCellVolumes = tgtMesh().cellVolumes();

    const faceList& srcFaces = srcMesh().faces();
    const cellList& srcCells = srcMesh().cells();
    const pointField& srcPoints = srcMesh().points();
    const vectorField& srcFaceCentres = srcMesh().faceCentres();
    const vectorField& srcCellCentres = srcMesh().cellCentres();

    const cell& tgtCell = tgtCells[index];
    const vector& tgtCentre = tgtCellCentres[index];
    const scalar& tgtVolume = tgtCellVolumes[index];

    if (srcCandidate < 0)
    {
        // Loop through all cells, and find minimum
        label minCell = -1;
        scalar minDist = GREAT;

        forAll(srcCellCentres, cellI)
        {
            const vector& srcCentre = srcCellCentres[cellI];

            if (magSqr(tgtCentre - srcCentre) < minDist)
            {
                minDist = magSqr(tgtCentre - srcCentre);
                minCell = cellI;
            }
        }

        // Reassign old candidate
        mapCandidate = minCell;
    }
    else
    {
        mapCandidate = srcCandidate;
    }

    // Maintain a check-list
    labelHashSet checked, skipped;

    // Decompose source / targets
    DynamicList<TetPoints> srcTets(15);
    DynamicList<TetPoints> tgtTets(15);

    // Configure the target cell
    decomposeCell
    (
        tgtCell,
        tgtCentre,
        tgtFaces,
        tgtPoints,
        tgtFaceCentres,
        tgtTets
    );

    // Loop and add intersections until nothing changes
    do
    {
        // Reset flag
        changed = false;

        // Fetch the set of candidates
        labelList checkList;

        if (nAttempts == 0)
        {
            checkList = labelList(1, mapCandidate);
        }
        else
        {
            checkList = checked.toc();
        }

        forAll(checkList, indexI)
        {
            labelList checkEntities;

            if (nAttempts == 0)
            {
                checkEntities = labelList(1, checkList[indexI]);
            }
            else
            {
                checkEntities = srcNeighbourList[checkList[indexI]];
            }

            forAll(checkEntities, entityI)
            {
                label checkEntity = checkEntities[entityI];

                // Skip if this is already
                // on the checked / skipped list
                if
                (
                    (checked.found(checkEntity)) ||
                    (skipped.found(checkEntity))
                )
                {
                    continue;
                }

                const cell& srcCell = srcCells[checkEntity];
                const vector& srcCentre = srcCellCentres[checkEntity];

                // Decompose this source cell
                decomposeCell
                (
                    srcCell,
                    srcCentre,
                    srcFaces,
                    srcPoints,
                    srcFaceCentres,
                    srcTets
                );

                // Accumulate volume / centroid over all intersections
                bool foundIntersect = false;

                scalar totalVolume = 0.0;
                vector totalCentre = vector::zero;

                forAll(tgtTets, tetI)
                {
                    // Initialize the intersection object
                    tetIntersection intersector(tgtTets[tetI]);

                    forAll(srcTets, tetJ)
                    {
                        // Evaluate for intersection
                        bool intersects = intersector.evaluate(srcTets[tetJ]);

                        if (intersects)
                        {
                            scalar volume = 0.0;
                            vector centre = vector::zero;

                            // Get volume / centroid
                            intersector.getVolumeAndCentre(volume, centre);

                            // Accumulate volume / centroid
                            totalVolume += volume;
                            totalCentre += (volume * centre);

                            foundIntersect = true;
                        }
                    }
                }

                if (foundIntersect)
                {
                    // Normalize centre
                    totalCentre /= totalVolume + VSMALL;

                    label oldSize = parents.size();

                    parents.setSize(oldSize + 1, checkEntity);
                    weights.setSize(oldSize + 1, totalVolume);
                    centres.setSize(oldSize + 1, totalCentre);

                    if (!checked.found(checkEntity))
                    {
                        checked.insert(checkEntity);
                    }

                    changed = true;
                }
                else
                {
                    // Add to the skipped list
                    if (!skipped.found(checkEntity))
                    {
                        skipped.insert(checkEntity);
                    }
                }
            }
        }

        nAttempts++;

        // Break out if we're taking too long
        if (nAttempts > 20)
        {
            break;
        }

    } while (changed);

    bool consistent = (mag(1.0 - (sum(weights)/tgtVolume)) < 1e-13);

    // Normalize weights
    weights /= tgtVolume;

    return consistent;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
