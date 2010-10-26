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
    convexSetAlgorithm

Description
    Base class for convexSetAlgorithms

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "meshOps.H"
#include "polyMesh.H"

#include "convexSetAlgorithm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
convexSetAlgorithm::convexSetAlgorithm
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
    twoDMesh_(mesh.nGeometricD() == 2),
    nOldPoints_(mesh.nPoints()),
    mesh_(mesh),
    newPoints_(newPoints),
    newEdges_(newEdges),
    newFaces_(newFaces),
    newCells_(newCells),
    newOwner_(newOwner),
    newNeighbour_(newNeighbour),
    pointsFromPoints_(pointsFromPoints),
    modPoints_(modPoints),
    highPrecision_(false)
{}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

// Output an entity as a VTK file
void convexSetAlgorithm::writeVTK
(
    const word& name,
    const label entity,
    const label primitiveType,
    const bool useOldConnectivity
) const
{
    writeVTK
    (
        name,
        labelList(1, entity),
        primitiveType,
        useOldConnectivity
    );
}


// Output a list of entities as a VTK file
void convexSetAlgorithm::writeVTK
(
    const word& name,
    const labelList& cList,
    const label primitiveType,
    const bool useOldConnectivity
) const
{
    if (useOldConnectivity)
    {
        const polyMesh& mesh = this->mesh_;

        meshOps::writeVTK
        (
            mesh,
            name,
            cList,
            primitiveType,
            mesh.points(),
            mesh.edges(),
            mesh.faces(),
            mesh.cells(),
            mesh.faceOwner()
        );
    }
    else
    {
        meshOps::writeVTK
        (
            this->mesh_,
            name,
            cList,
            primitiveType,
            newPoints_,
            newEdges_,
            newFaces_,
            newCells_,
            newOwner_
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool convexSetAlgorithm::consistent(const scalar tolerance) const
{
    if (weights_.size())
    {
        scalar normError = mag(1.0 - (sum(weights_)/normFactor_));

        if (normError < tolerance)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

#   if USE_MPFR

    if (weights_.size() && mpWeights_.size())
    {
        FatalErrorIn
        (
            "bool convexSetAlgorithm::consistent"
            "(const scalar tolerance) const"
        )
            << nl << " Inconsistent internal lists: " << nl
            << " weights: " << weights_ << nl
            << " mpWeights: " << mpWeights_ << nl
            << abort(FatalError);
    }

    if (mpWeights_.size())
    {
        scalar sumNormWeights = 0.0;

        forAll(mpWeights_, indexI)
        {
            sumNormWeights += (mpWeights_[indexI]/mpNormFactor_);
        }

        scalar normError = mag(1.0 - sumNormWeights);

        if (normError < tolerance)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

#   endif

    return false;
}


// Return the normFactor
scalar convexSetAlgorithm::normFactor() const
{
#   if USE_MPFR
    if (highPrecision())
    {
        return mpNormFactor_;
    }
    else
#   endif
    {
        return normFactor_;
    }

    return 0.0;
}


// Normalize stored weights
void convexSetAlgorithm::normalize(bool normSum) const
{
    if (normSum)
    {
        if (weights_.size())
        {
            weights_ /= sum(weights_);
        }
    }
    else
    {
        if (weights_.size())
        {
            weights_ /= normFactor_;
        }
    }

#   if USE_MPFR

    if (weights_.size() && mpWeights_.size())
    {
        FatalErrorIn
        (
            "void convexSetAlgorithm::normalize(bool normSum) const"
        )
            << nl << " Inconsistent internal lists: " << nl
            << " weights: " << weights_ << nl
            << " mpWeights: " << mpWeights_ << nl
            << abort(FatalError);
    }

    if (normSum)
    {
        if (mpWeights_.size())
        {
            mpScalar sumWeights = sum(mpWeights_);

            forAll(mpWeights_, indexI)
            {
                mpWeights_[indexI] /= sumWeights;
            }
        }
    }
    else
    {
        if (mpWeights_.size())
        {
            forAll(mpWeights_, indexI)
            {
                mpWeights_[indexI] /= mpNormFactor_;
            }
        }
    }

#   endif
}


// Extract weights and centres to lists
void convexSetAlgorithm::populateLists
(
    labelList& parents,
    vectorField& centres,
    scalarField& weights
) const
{
    // Clear inputs
    parents.clear();
    centres.clear();
    weights.clear();

    if (weights_.size())
    {
        parents = parents_;
        centres = centres_;
        weights = weights_;
    }

#   if USE_MPFR

    if (weights_.size() && mpWeights_.size())
    {
        FatalErrorIn
        (
            "\n\n"
            "void convexSetAlgorithm::populateLists\n"
            "(\n"
            "    labelList& parents,\n"
            "    vectorField& centres,\n"
            "    scalarField& weights\n"
            ") const\n"
        )
            << nl << " Inconsistent internal lists: " << nl
            << " weights: " << weights_ << nl
            << " mpWeights: " << mpWeights_ << nl
            << abort(FatalError);
    }

    if (mpWeights_.size())
    {
        parents = parents_;

        centres.setSize(parents.size());
        weights.setSize(parents.size());

        forAll(parents, i)
        {
            centres[i] = convert<scalar>(mpCentres_[i]);
            weights[i] = mpWeights_[i];
        }
    }

#   endif
}


// Return the highPrecision switch
bool convexSetAlgorithm::highPrecision() const
{
    return highPrecision_;
}


// Set higher precision mode
void convexSetAlgorithm::setHighPrecision() const
{
    highPrecision_ = true;
}


// Reset higher precision mode
void convexSetAlgorithm::unsetHighPrecision() const
{
    highPrecision_ = false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
