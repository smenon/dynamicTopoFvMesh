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
    interpolator

Description
    Helper class to perform field interpolation during topo-changes

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*----------------------------------------------------------------------------*/

#include "interpolator.H"
#include "mapPolyMesh.H"
#include "dynamicTopoFvMesh.H"
#include "GeometricField.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "fvPatchFieldsFwd.H"
#include "volFields.H"
#include "fvsPatchFieldsFwd.H"
#include "surfaceFields.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(interpolator,0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

interpolator::interpolator(dynamicTopoFvMesh& mesh)
:
    mesh_(mesh),
    fieldsRegistered_(true),
    phiPtr_(NULL)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

interpolator::~interpolator()
{
    clearOut();
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void interpolator::clearOut()
{
    deleteDemandDrivenData(phiPtr_);
}

//- Print out mapped fields
void interpolator::printFields() const
{
    // Output all fields that will be used for mapping.

    Info << "volScalarFields: " << endl;

    HashTable<const volScalarField*> vsf(mesh_.lookupClass<volScalarField>());

    forAllIter(HashTable<const volScalarField*>, vsf, fIter)
    {
        Info << "\t" << fIter()->name() << endl;
    }

    Info << "volVectorFields: " << endl;

    HashTable<const volVectorField*> vvf(mesh_.lookupClass<volVectorField>());

    forAllIter(HashTable<const volVectorField*>, vvf, fIter)
    {
        Info << "\t" << fIter()->name() << endl;
    }

    Info << "surfaceScalarFields: " << endl;

    HashTable<const surfaceScalarField*> ssf
    (
        mesh_.lookupClass<surfaceScalarField>()
    );

    forAllIter(HashTable<const surfaceScalarField*>, ssf, fIter)
    {
        Info << "\t" << fIter()->name() << endl;
    }
}

void interpolator::makePhi() const
{
    if (debug)
    {
        Info << "void interpolator::makePhi() const : "
             << "Assembling volume fluxes"
             << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (phiPtr_)
    {
        FatalErrorIn("interpolator::makePhi() const")
            << "Volume fluxes already exist"
            << abort(FatalError);
    }

    phiPtr_ = new resizableList<scalar>(mesh_.nFaces(), 0.0);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Add a new face entry
void interpolator::insertFace()
{
    // (*phiPtr_).append(0.0);
}

// Add a new cell entry.
void interpolator::insertCell
(
    const labelList& mapCells,
    const scalarList& mapWeights
)
{
    // Now loop through all volFields
    // registered for interpolation,
    // and map values for the new cell.

}

// Set the volume-flux for an existing face
void interpolator::setVolFlux
(
    const label faceIndex,
    const scalar volFlux
)
{
    // (*phiPtr_)[faceIndex] = volFlux;
}

// Register fields for interpolation
void interpolator::registerFields()
{
    // If fields have already been registered, bail out.
    if (fieldsRegistered_)
    {
        return;
    }

    bool registerSelectFields = false;

    // First check to see if only certain
    // fields need to be mapped.
    if (mesh_.dict_.subDict("dynamicTopoFvMesh").found("mapFields"))
    {
        registerSelectFields = true;
    }

    if (registerSelectFields)
    {
        const dictionary& mapDict =
        (
            mesh_.dict_.subDict("dynamicTopoFvMesh").subDict("mapFields")
        );

        if (mapDict.size() == 0)
        {
            FatalErrorIn("interpolator::registerFields()")
                << " Empty mapFields entry found." << nl
                << " Either specify fields required for mapping," << nl
                << " or remove the mapFields entry." << nl
                << abort(FatalError);
        }
    }
    else
    {
        // Register all fields from the registry.
        FatalErrorIn("interpolator::registerFields()")
            << " No mapFields entry found"
            << " in the dynamicTopoFvMesh sub-dictionary." << nl
            << " Mapping all registry fields is not currently available."
            << abort(FatalError);
    }

    // Set the flag.
    fieldsRegistered_ = true;
}

// Update fields after a topo-change operation
void interpolator::updateMesh(const mapPolyMesh& mpm)
{
    if (debug)
    {
        printFields();
    }

    // Check if phi exists
    bool foundPhi = mesh_.foundObject<surfaceScalarField>("phi");

    if (!foundPhi)
    {
        return;
    }

    label nOldFaces = mpm.nOldFaces(), newIndex = -1;
    const labelList& reverseFaceMap = mpm.reverseFaceMap();

    // Obtain references
    resizableList<scalar>& phi = *phiPtr_;

    // Fetch the volume-flux field from the registry.
    surfaceScalarField& mapPhi =
    (
        const_cast<surfaceScalarField&>
        (
            mesh_.lookupObject<surfaceScalarField>("phi")
        )
    );

    // Update mesh fluxes
    forAll(phi, faceI)
    {
        if (faceI < nOldFaces)
        {
            newIndex = reverseFaceMap[faceI];

            // Ensure that the index doesn't
            // belong to a deleted face.
            if (newIndex < 0)
            {
                continue;
            }
        }
        else
        {
            if (mesh_.deletedFaces_.found(faceI))
            {
                continue;
            }

            newIndex = mesh_.addedFaceRenumbering_[faceI];
        }

        // Figure out which patch the face belongs to,
        // and add accordingly to internal or boundary.
        label patch = mpm.mesh().boundaryMesh().whichPatch(newIndex);

        if (patch == -1)
        {
            mapPhi.internalField()[newIndex] = phi[faceI];
        }
        else
        {
            label i = mpm.mesh().boundaryMesh()[patch].whichFace(newIndex);

            mapPhi.boundaryField()[patch][i] = phi[faceI];
        }
    }

    // Clear-out demand-driven data.
    clearOut();
}

} // End namespace Foam

// ************************************************************************* //
