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
#include "triIntersection.H"

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
    const UList<label>& newNeighbour
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
        newNeighbour
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void faceSetAlgorithm::computeNormFactor(const label index) const
{
    // Clear existing fields
    parents_.clear();

    centres_.clear();
    weights_.clear();

    // Compute refNorm and normFactor
    refNorm_ = newFaces_[index].normal(newPoints_);
    normFactor_ = mag(refNorm_);

    // Normalize for later use
    refNorm_ /= normFactor_ + VSMALL;
}


// Compute intersections
bool faceSetAlgorithm::computeIntersection
(
    const label newIndex,
    const label oldIndex,
    bool output
) const
{
    bool intersects = false;

    const pointField& newPoints = newPoints_;
    const pointField& oldPoints = mesh_.points();

    const face& newFace = newFaces_[newIndex];
    const face& oldFace = mesh_.faces()[oldIndex];

    if (twoDMesh_ && newFace.size() == 4 && oldFace.size() == 4)
    {
        // Decompose new / old quad faces into 4 tris each
        FixedList<FixedList<point, 3>, 4> clippingTris
        (
            FixedList<point, 3>(vector::zero)
        );

        FixedList<FixedList<point, 3>, 4> subjectTris
        (
            FixedList<point, 3>(vector::zero)
        );

        label ntOld = 0, ntNew = 0;
        vector oldCentre = vector::zero, newCentre = vector::zero;

        // Configure tris from oldFace
        vector ofCentre = oldFace.centre(oldPoints);

        forAll(oldFace, pI)
        {
            subjectTris[ntOld][0] = oldPoints[oldFace[pI]];
            subjectTris[ntOld][1] = oldPoints[oldFace.nextLabel(pI)];
            subjectTris[ntOld][2] = ofCentre;

            ntOld++;
        }

        // Configure tris from newFace
        vector nfCentre = newFace.centre(newPoints);

        forAll(newFace, pI)
        {
            clippingTris[ntNew][0] = newPoints[newFace[pI]];
            clippingTris[ntNew][1] = newPoints[newFace.nextLabel(pI)];
            clippingTris[ntNew][2] = nfCentre;

            ntNew++;
        }

        // Accumulate area / centroid over all intersections
        bool foundIntersect = false;

        scalar totalArea = 0.0;
        vector totalCentre = vector::zero;

        // Loop through all clipping tris
        for (label i = 0; i < 4; i++)
        {
            // Initialize the intersector
            triIntersection intersector(clippingTris[i]);

            // Test for intersection and evaluate
            // against all subject tris
            for (label j = 0; j < 4; j++)
            {
                intersects = intersector.evaluate(subjectTris[j]);

                if (intersects)
                {
                    scalar area = 0.0;
                    vector centre = vector::zero;

                    // Fetch area and centre
                    intersector.getAreaAndCentre(area, centre);

                    // Accumulate area / centroid
                    totalArea += area;
                    totalCentre += (area * centre);

                    foundIntersect = true;
                }
            }
        }

        // Size-up the internal lists
        if (foundIntersect && !output)
        {
            // Normalize centre
            totalCentre /= totalArea + VSMALL;

            // Normalize and check if this is worth it
            if (mag(totalArea/normFactor_) > SMALL)
            {
                meshOps::sizeUpList(oldIndex, parents_);
                meshOps::sizeUpList(totalArea, weights_);
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
        // Configure points for clipping triangle
        FixedList<point, 3> clippingTri(vector::zero);

        // Fill in points
        clippingTri[0] = newPoints[newFace[0]];
        clippingTri[1] = newPoints[newFace[1]];
        clippingTri[2] = newPoints[newFace[2]];

        // Configure points for subject triangle
        FixedList<point, 3> subjectTri(vector::zero);

        // Fill in points
        subjectTri[0] = oldPoints[oldFace[0]];
        subjectTri[1] = oldPoints[oldFace[1]];
        subjectTri[2] = oldPoints[oldFace[2]];

        // Initialize the intersector
        triIntersection intersector(clippingTri);

        // Test for intersection and evaluate
        intersects = intersector.evaluate(subjectTri);

        if (intersects)
        {
            scalar area;
            vector centre;

            // Fetch area and centre
            intersector.getAreaAndCentre(area, centre);

            // Normalize and check if this is worth it
            if (mag(area/normFactor_) > SMALL)
            {
                // Size-up the internal lists
                if (!output)
                {
                    meshOps::sizeUpList(oldIndex, parents_);
                    meshOps::sizeUpList(area, weights_);
                    meshOps::sizeUpList(centre, centres_);
                }
            }
            else
            {
                intersects = false;
            }
        }
    }

    return intersects;
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
