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

    // Fetch nearest-cell addressing from meshToMesh
    const labelList& cAddr = meshToMeshPtr_().cellAddressing();

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
        scalarField& volumes = volumes_[cellI];
        vectorField& centres = centres_[cellI];

        label precisionAttempts = 0;

        // Obtain weighting factors for this cell.
        bool consistent =
        (
            computeWeights
            (
                cellI,
                cAddr[cellI],
                srcMesh().cellCells(),
                matchTol,
                precisionAttempts,
                parents,
                weights,
                volumes,
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
            fromMesh().time().timeName(),
            fromMesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Check for compatibility
    bool compatible = true;
    label targetCells = toMesh().nCells();
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
    IOList<scalarField> srcVolumes
    (
        IOobject
        (
            "volumes",
            fromMesh().time().timeName(),
            fromMesh(),
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
            fromMesh().time().timeName(),
            fromMesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Set sizes
    forAll(nCellsPerCell, cellI)
    {
        addressing_[cellI].setSize(nCellsPerCell[cellI]);
        volumes_[cellI].setSize(nCellsPerCell[cellI]);
        centres_[cellI].setSize(nCellsPerCell[cellI]);
    }

    nCellsPerCell = 0;

    // Invert addressing
    forAll(srcAddressing, cellI)
    {
        const labelList& srcAddr = srcAddressing[cellI];
        const scalarField& srcVl = srcVolumes[cellI];
        const vectorField& srcCt = srcCentres[cellI];

        forAll(srcAddr, j)
        {
            label cellJ = srcAddr[j];

            addressing_[cellJ][nCellsPerCell[cellJ]] = cellI;
            volumes_[cellJ][nCellsPerCell[cellJ]] = srcVl[j];
            centres_[cellJ][nCellsPerCell[cellJ]] = srcCt[j];

            nCellsPerCell[cellJ]++;
        }
    }

    // Check weights for consistency
    const scalarField& V = toMesh().cellVolumes();

    forAll(V, cellI)
    {
        if (mag((V[cellI] - sum(volumes_[cellI])) / V[cellI]) > 5e-14)
        {
            Info<< " Weights are not compatible. " << nl
                << " Cell: " << cellI << nl
                << " Volume: " << V[cellI] << nl
                << " Sum(volumes): " << sum(volumes_[cellI]) << nl
                << " Error: "
                << mag((V[cellI] - sum(volumes_[cellI])) / V[cellI])
                << endl;

            compatible = false;
            break;
        }
    }

    if (!compatible)
    {
        addressing_.clear();
        weights_.clear();
        volumes_.clear();
        centres_.clear();

        addressing_.setSize(targetCells);
        weights_.setSize(targetCells);
        volumes_.setSize(targetCells);
        centres_.setSize(targetCells);

        return false;
    }

    // Return a successful inversion
    return true;
}


// Compute weighting factors for a particular cell
bool conservativeMeshToMesh::computeWeights
(
    const label index,
    const label oldCandidate,
    const labelListList& oldNeighbourList,
    const scalar mTol,
    label& precisionAttempts,
    labelList& parents,
    scalarField& weights,
    scalarField& volumes,
    vectorField& centres,
    bool highPrecision
) const
{
    if (parents.size() || weights.size() || volumes.size() || centres.size())
    {
        FatalErrorIn
        (
            "\n\n"
            "bool conservativeMeshToMesh::computeWeights\n"
            "(\n"
            "    const label index,\n"
            "    const label oldCandidate,\n"
            "    const labelListList& oldNeighbourList,\n"
            "    const scalar mTol,\n"
            "    label& precisionAttempts,\n"
            "    labelList& parents,\n"
            "    scalarField& weights,\n"
            "    scalarField& volumes,\n"
            "    vectorField& centres,\n"
            "    bool highPrecision\n"
            ") const\n"
        )
            << " Addressing has already been calculated." << nl
            << " Index: " << index << nl
            << " oldCandidate: " << oldCandidate << nl
            << " Parents: " << parents << nl
            << " Weights: " << weights << nl
            << " Volumes: " << volumes << nl
            << " Centres: " << centres << nl
            << abort(FatalError);
    }

    bool changed;
    label nAttempts = 0, nIntersects = 0;

    label mapCandidate = -1;

    if (oldCandidate < 0)
    {
        // Loop through all cells, and find minimum
        scalar minDist = GREAT;
        label minCell = -1;

        const vectorField& oldCentres = fromMesh().cellCentres();
        const vector newCentre = toMesh().cellCentres()[index];

        forAll(oldCentres, cellI)
        {
            if (magSqr(newCentre - oldCentres[cellI]) < minDist)
            {
                minDist = magSqr(newCentre - oldCentres[cellI]);
                minCell = cellI;
            }
        }

        // Reassign old candidate
        mapCandidate = minCell;
    }
    else
    {
        mapCandidate = oldCandidate;
    }

    // Fetch the volume of the new cell
    scalar newCellVolume = toMesh().cellVolumes()[index];

    // Maintain a check-list
    labelHashSet checked, skipped;

    // Configure the input tet
    FixedList<point, 4> tgtTetPoints, srcTetPoints;

    const cell& tgtCell = tgtMesh().cells()[index];

    tgtTetPoints =
    (
        tgtCell.points
        (
            tgtMesh().faces(),
            tgtMesh().points()
        )
    );

    // Initialize the intersection object
    tetIntersection tI(tgtTetPoints);

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
                checkEntities = oldNeighbourList[checkList[indexI]];
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

                const cell& srcCell = srcMesh().cells()[checkEntity];

                srcTetPoints =
                (
                    srcCell.points
                    (
                        srcMesh().faces(),
                        srcMesh().points()
                    )
                );

                // Evaluate for intersection
                bool intersect = tI.evaluate(srcTetPoints);

                if (intersect)
                {
                    scalar volume = 0.0;
                    vector centre = vector::zero;

                    // Get volume / centroid
                    tI.getVolumeAndCentre(volume, centre);

                    label oldSize = parents.size();

                    parents.setSize(oldSize + 1, checkEntity);
                    weights.setSize(oldSize + 1, volume);
                    volumes.setSize(oldSize + 1, volume);
                    centres.setSize(oldSize + 1, centre);

                    nIntersects++;

                    if (!checked.found(checkEntity))
                    {
                        checked.insert(checkEntity);
                    }

                    changed = true;
                }
                else
                if (nAttempts == 0)
                {
                    FatalErrorIn
                    (
                        "\n\n"
                        "bool conservativeMeshToMesh::computeWeights\n"
                        "(\n"
                        "    const label index,\n"
                        "    const label oldCandidate,\n"
                        "    const labelListList& oldNeighbourList,\n"
                        "    const scalar mTol,\n"
                        "    label& precisionAttempts,\n"
                        "    labelList& parents,\n"
                        "    scalarField& weights,\n"
                        "    scalarField& volumes,\n"
                        "    vectorField& centres,\n"
                        "    bool highPrecision\n"
                        ") const\n"
                    )
                        << " First intersection was not found." << nl
                        << " Index: " << index << nl
                        << " oldCandidate: " << oldCandidate << nl
                        << " mapCandidate: " << mapCandidate << nl
                        << abort(FatalError);
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

    bool consistent = (mag(1.0 - (sum(weights)/newCellVolume)) < 1e-13);

    // Normalize weights
    weights /= newCellVolume;

    return consistent;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
