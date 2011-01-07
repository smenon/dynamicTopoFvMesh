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

#include "Time.H"
#include "IOMap.H"
#include "meshOps.H"
#include "polyMesh.H"
#include "objectMap.H"
#include "edgeIOList.H"
#include "cellIOList.H"

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
    const UList<label>& newNeighbour
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
    newNeighbour_(newNeighbour)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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

    return false;
}


// Return the normFactor
scalar convexSetAlgorithm::normFactor() const
{
    return normFactor_;
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
}


// Write out connectivity information to disk
bool convexSetAlgorithm::write() const
{
    pointIOField
    (
        IOobject
        (
            "newPoints",
            mesh_.time().timeName(),
            "convexSetAlgorithm",
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        newPoints_
    ).writeObject
    (
        IOstream::BINARY,
        IOstream::currentVersion,
        IOstream::COMPRESSED
    );

    edgeIOList
    (
        IOobject
        (
            "newEdges",
            mesh_.time().timeName(),
            "convexSetAlgorithm",
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        newEdges_
    ).writeObject
    (
        IOstream::BINARY,
        IOstream::currentVersion,
        IOstream::COMPRESSED
    );

    faceIOList
    (
        IOobject
        (
            "newFaces",
            mesh_.time().timeName(),
            "convexSetAlgorithm",
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        newFaces_
    ).writeObject
    (
        IOstream::BINARY,
        IOstream::currentVersion,
        IOstream::COMPRESSED
    );

    cellIOList
    (
        IOobject
        (
            "newCells",
            mesh_.time().timeName(),
            "convexSetAlgorithm",
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        newCells_
    ).writeObject
    (
        IOstream::BINARY,
        IOstream::currentVersion,
        IOstream::COMPRESSED
    );

    labelIOList
    (
        IOobject
        (
            "newOwner",
            mesh_.time().timeName(),
            "convexSetAlgorithm",
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        newOwner_
    ).writeObject
    (
        IOstream::BINARY,
        IOstream::currentVersion,
        IOstream::COMPRESSED
    );

    labelIOList
    (
        IOobject
        (
            "newNeighbour",
            mesh_.time().timeName(),
            "convexSetAlgorithm",
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        newNeighbour_
    ).writeObject
    (
        IOstream::BINARY,
        IOstream::currentVersion,
        IOstream::COMPRESSED
    );

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
