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
    cellSetAlgorithm

Description
    Class for convexSetAlgorithms pertaining to cells

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "cellSetAlgorithm.H"

#include "meshOps.H"
#include "triFace.H"
#include "polyMesh.H"
#include "tetIntersection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
cellSetAlgorithm::cellSetAlgorithm
(
    const polyMesh& mesh,
    const bool isTwoDMesh,
    const pointField& newPoints,
    const UList<edge>& newEdges,
    const UList<face>& newFaces,
    const UList<cell>& newCells,
    const UList<label>& newOwner,
    const UList<label>& newNeighbour
)
:
    convexSetAlgorithm
    (
        mesh,
        isTwoDMesh,
        newPoints,
        newEdges,
        newFaces,
        newCells,
        newOwner,
        newNeighbour
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void cellSetAlgorithm::computeNormFactor(const label index) const
{
    // Clear existing fields
    parents_.clear();

    centres_.clear();
    weights_.clear();

    // Compute volume / centre (using refNorm_ as centre)
    meshOps::cellCentreAndVolume
    (
        index,
        newPoints_,
        newFaces_,
        newCells_,
        newOwner_,
        refNorm_,
        normFactor_
    );
}


bool cellSetAlgorithm::computeIntersection
(
    const label newIndex,
    const label oldIndex,
    const label offset,
    bool output
) const
{
    bool intersects = false;

    const pointField& newPoints = newPoints_;
    const pointField& oldPoints = mesh_.points();

    const cell& newCell = newCells_[newIndex];
    const cell& oldCell = mesh_.cells()[oldIndex];

    if (twoDMesh_)
    {
        // Decompose new / old prism cells into 14 tets each
        FixedList<FixedList<point, 4>, 14> clippingTets
        (
            FixedList<point, 4>(vector::zero)
        );

        FixedList<FixedList<point, 4>, 14> subjectTets
        (
            FixedList<point, 4>(vector::zero)
        );

        label ntOld = 0, ntNew = 0;
        vector oldCentre = vector::zero, newCentre = vector::zero;

        // Configure tets from oldCell
        forAll(oldCell, faceI)
        {
            const face& oldFace = mesh_.faces()[oldCell[faceI]];

            vector fCentre = oldFace.centre(oldPoints);

            if (oldFace.size() == 3)
            {
                subjectTets[ntOld][0] = oldPoints[oldFace[0]];
                subjectTets[ntOld][1] = oldPoints[oldFace[1]];
                subjectTets[ntOld][2] = oldPoints[oldFace[2]];

                ntOld++;
            }
            else
            {
                forAll(oldFace, pI)
                {
                    subjectTets[ntOld][0] = oldPoints[oldFace[pI]];
                    subjectTets[ntOld][1] = oldPoints[oldFace.nextLabel(pI)];
                    subjectTets[ntOld][2] = fCentre;

                    ntOld++;
                }
            }

            oldCentre += fCentre;
        }

        // Configure tets from newCell
        forAll(newCell, faceI)
        {
            const face& newFace = newFaces_[newCell[faceI]];

            vector fCentre = newFace.centre(newPoints);

            if (newFace.size() == 3)
            {
                clippingTets[ntNew][0] = newPoints[newFace[0]];
                clippingTets[ntNew][1] = newPoints[newFace[1]];
                clippingTets[ntNew][2] = newPoints[newFace[2]];

                ntNew++;
            }
            else
            {
                forAll(newFace, pI)
                {
                    clippingTets[ntNew][0] = newPoints[newFace[pI]];
                    clippingTets[ntNew][1] = newPoints[newFace.nextLabel(pI)];
                    clippingTets[ntNew][2] = fCentre;

                    ntNew++;
                }
            }

            newCentre += fCentre;
        }

        oldCentre /= 5.0;
        newCentre /= 5.0;

        // Fill-in last points for all tets
        for (label i = 0; i < 14; i++)
        {
            subjectTets[i][3] = oldCentre;
            clippingTets[i][3] = newCentre;
        }

        // Accumulate volume / centroid over all intersections
        bool foundIntersect = false;

        scalar totalVolume = 0.0;
        vector totalCentre = vector::zero;

        // Loop through all clipping tets
        for (label i = 0; i < 14; i++)
        {
            // Initialize the intersector
            tetIntersection intersector(clippingTets[i]);

            // Test for intersection and evaluate
            // against all subject tets
            for (label j = 0; j < 14; j++)
            {
                intersects = intersector.evaluate(subjectTets[j]);

                if (intersects)
                {
                    scalar volume = 0.0;
                    vector centre = vector::zero;

                    // Fetch volume and centre
                    intersector.getVolumeAndCentre(volume, centre);

                    // Accumulate volume / centroid
                    totalVolume += volume;
                    totalCentre += (volume * centre);

                    foundIntersect = true;
                }
            }
        }

        // Size-up the internal lists
        if (foundIntersect && !output)
        {
            // Normalize centre
            totalCentre /= totalVolume + VSMALL;

            // Normalize and check if this is worth it
            if (mag(totalVolume/normFactor_) > SMALL)
            {
                meshOps::sizeUpList((oldIndex - offset), parents_);
                meshOps::sizeUpList(totalVolume, weights_);
                meshOps::sizeUpList(totalCentre, centres_);

                intersects = true;
            }
            else
            {
                intersects = false;
            }
        }
    }
    else
    {
        // Configure points for clipping tetrahedron
        FixedList<point, 4> clippingTet(vector::zero);

        // Configure the clipping tetrahedron.
        const face& firstNewFace = newFaces_[newCell[0]];
        const face& secondNewFace = newFaces_[newCell[1]];

        // Find the isolated point
        label fourthNewPoint =
        (
            meshOps::findIsolatedPoint
            (
                firstNewFace,
                secondNewFace
            )
        );

        // Fill in points
        clippingTet[0] = newPoints[firstNewFace[0]];
        clippingTet[1] = newPoints[firstNewFace[1]];
        clippingTet[2] = newPoints[firstNewFace[2]];
        clippingTet[3] = newPoints[fourthNewPoint];

        // Configure points for subject tetrahedron
        FixedList<point, 4> subjectTet(vector::zero);

        // Configure the subject tetrahedron.
        const face& firstOldFace = mesh_.faces()[oldCell[0]];
        const face& secondOldFace = mesh_.faces()[oldCell[1]];

        // Find the isolated point
        label fourthOldPoint =
        (
            meshOps::findIsolatedPoint
            (
                firstOldFace,
                secondOldFace
            )
        );

        // Fill in points
        subjectTet[0] = oldPoints[firstOldFace[0]];
        subjectTet[1] = oldPoints[firstOldFace[1]];
        subjectTet[2] = oldPoints[firstOldFace[2]];
        subjectTet[3] = oldPoints[fourthOldPoint];

        // Initialize the intersector
        tetIntersection intersector(clippingTet);

        // Test for intersection and evaluate
        intersects = intersector.evaluate(subjectTet);

        if (intersects)
        {
            scalar volume = 0.0;
            vector centre = vector::zero;

            // Fetch volume and centre
            intersector.getVolumeAndCentre(volume, centre);

            // Normalize and check if this is worth it
            if (mag(volume/normFactor_) > SMALL)
            {
                if (output)
                {
                    writeVTK
                    (
                        "tetIntersectNew_"
                      + Foam::name(newIndex),
                        List<FixedList<point, 4> >(1, clippingTet)
                    );

                    writeVTK
                    (
                        "tetIntersectOld_"
                      + Foam::name(newIndex)
                      + '_'
                      + Foam::name(oldIndex),
                        List<FixedList<point, 4> >(1, subjectTet)
                    );

                    writeVTK
                    (
                        "tetIntersect_"
                      + Foam::name(newIndex)
                      + '_'
                      + Foam::name(oldIndex),
                        intersector.getIntersection()
                    );
                }

                // Size-up the internal lists
                meshOps::sizeUpList((oldIndex - offset), parents_);
                meshOps::sizeUpList(volume, weights_);
                meshOps::sizeUpList(centre, centres_);
            }
            else
            {
                intersects = false;
            }
        }
    }

    return intersects;
}


//- Write out tets as a VTK
void cellSetAlgorithm::writeVTK
(
    const word& name,
    const List<FixedList<point, 4> >& tetList
) const
{
    // Fill up all points
    label pI = 0;

    List<cell> allTets(tetList.size(), cell(4));
    pointField allPoints(4 * tetList.size());

    forAll(tetList, tetI)
    {
        allTets[tetI][0] = pI;
        allPoints[pI++] = tetList[tetI][0];

        allTets[tetI][1] = pI;
        allPoints[pI++] = tetList[tetI][1];

        allTets[tetI][2] = pI;
        allPoints[pI++] = tetList[tetI][2];

        allTets[tetI][3] = pI;
        allPoints[pI++] = tetList[tetI][3];
    }

    // Write out in cell-to-node addressing
    meshOps::writeVTK
    (
        mesh_,
        name,
        identity(tetList.size()),
        3,
        allPoints,
        List<edge>(0),
        List<face>(0),
        allTets,
        List<label>(0)
    );
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
