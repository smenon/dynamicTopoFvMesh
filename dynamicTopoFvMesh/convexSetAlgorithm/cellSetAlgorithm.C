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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
cellSetAlgorithm::cellSetAlgorithm
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

void cellSetAlgorithm::computeNormFactor(const label index) const
{
    // Clear existing fields
    parents_.clear();

    centres_.clear();
    weights_.clear();

    // Compute volume / centre
    {
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
}


bool cellSetAlgorithm::computeIntersection
(
    const label newIndex,
    const label oldIndex,
    const scalar& matchTol,
    bool output
) const
{
    bool intersects = false;

    {
        // Invoke the conventional variant

        if (intersects)
        {
            scalar volume;
            vector centre;

            // Size-up the internal lists
            if (!output)
            {
                meshOps::sizeUpList(oldIndex, parents_);
                meshOps::sizeUpList(volume, weights_);
                meshOps::sizeUpList(centre, centres_);
            }
        }
    }

    return intersects;
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
