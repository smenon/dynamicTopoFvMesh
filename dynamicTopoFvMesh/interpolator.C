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
    fieldsRegistered_(false),
    meshPhiPtr_(NULL),
    V0Ptr_(NULL)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

interpolator::~interpolator()
{
    clearOut();
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void interpolator::clearOut()
{
    deleteDemandDrivenData(meshPhiPtr_);
    deleteDemandDrivenData(V0Ptr_);
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

void interpolator::makeMeshPhi() const
{
    if (debug)
    {
        Info << "void interpolator::makeMeshPhi() const : "
             << "Assembling mesh fluxes"
             << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (meshPhiPtr_)
    {
        FatalErrorIn("interpolator::makeMeshPhi() const")
            << "Mesh fluxes already exist"
            << abort(FatalError);
    }

    meshPhiPtr_ = new resizableList<scalar>(mesh_.nFaces(), 0.0);
}

void interpolator::makeV0() const
{
    if (debug)
    {
        Info << "void interpolator::makeV0() const : "
             << "Assembling old cell volumes"
             << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (V0Ptr_)
    {
        FatalErrorIn("interpolator::makeV0() const")
            << "Old cell volumes already exist"
            << abort(FatalError);
    }

    V0Ptr_ = new resizableList<scalar>(mesh_.nCells(), 0.0);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Add a new face entry
void interpolator::insertFace
(
    const scalar sweptVol
)
{
    scalar rDeltaT = 1.0/mesh_.time().deltaT().value();

    (*meshPhiPtr_).append(rDeltaT*sweptVol);
}

// Add a new cell entry.
void interpolator::insertCell
(
    const labelList& mapCells,
    const scalarList& mapWeights
)
{
    // Append a null value, which will be set later
    (*V0Ptr_).append(0.0);

    // Now loop through all volFields
    // registered for interpolation,
    // and map values for the new cell.

}

// Set the old-volume value for an existing cell
void interpolator::setOldVolume
(
    const label cellIndex,
    const scalar oldVolume
)
{
    if (oldVolume < 0.0)
    {
        WarningIn("interpolator::setOldVolume()") << nl
            << "Negative value prescribed for old-volume." << nl
            << "  cellIndex: " << cellIndex << nl
            << "  oldVolume: " << oldVolume << nl
            << endl;
    }

    (*V0Ptr_)[cellIndex] = oldVolume;
}

// Set the mesh-flux for an existing face
void interpolator::setMeshFlux
(
    const label faceIndex,
    const scalar sweptVol
)
{
    scalar rDeltaT = 1.0/mesh_.time().deltaT().value();

    (*meshPhiPtr_)[faceIndex] = (rDeltaT*sweptVol);
}

// Get the mesh-flux for an existing face.
scalar interpolator::getMeshFlux
(
    const label faceIndex
) const
{
    return (*meshPhiPtr_)[faceIndex];
}

// Flip the face-flux for an existing face
void interpolator::flipFaceFlux(const label faceIndex)
{
    (*meshPhiPtr_)[faceIndex] *= -1.0;

    // Loop through registered surfaceScalarFields and flip fluxes.
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

    }

    // Set the flag.
    fieldsRegistered_ = true;
}

// Update geometry for a repositioned mesh
void interpolator::movePoints(const scalarField& sweptVols) const
{
    // Update mesh-fluxes from the mesh.
    if (!meshPhiPtr_)
    {
        makeMeshPhi();
    }

    scalar rDeltaT = 1.0/mesh_.time().deltaT().value();

    resizableList<scalar>& meshPhi = *meshPhiPtr_;

    meshPhi = sweptVols;

    forAll(meshPhi, faceI)
    {
        meshPhi[faceI] *= rDeltaT;
    }

    // Update old-volumes from the mesh.
    if (!V0Ptr_)
    {
        makeV0();
    }

    resizableList<scalar>& V0 = *V0Ptr_;

    V0 = mesh_.V0().field();
}

// Update fields after a topo-change operation
void interpolator::updateMesh(const mapPolyMesh& mpm)
{
    if (debug)
    {
        printFields();
    }

    label nOldCells = mpm.nOldCells();
    const labelList& reverseCellMap = mpm.reverseCellMap();
    resizableList<scalar>& V0 = *V0Ptr_;

    // The setV0 member wipes out meshPhi, so this has to be done first.
    DimensionedField<scalar, volMesh>& mapV0 = mesh_.setV0();

    // Update old-volumes.
    forAll(V0, cellI)
    {
        if (cellI < nOldCells)
        {
            label newIndex = reverseCellMap[cellI];

            if (newIndex < 0)
            {
                continue;
            }

            mapV0[newIndex] = V0[cellI];
        }
        else
        {
            if (mesh_.deletedCells_.found(cellI))
            {
                continue;
            }

            mapV0[mesh_.addedCellRenumbering_[cellI]] = V0[cellI];
        }
    }

    label nOldFaces = mpm.nOldFaces(), newIndex = -1;
    const labelList& reverseFaceMap = mpm.reverseFaceMap();
    resizableList<scalar>& meshPhi = *meshPhiPtr_;

    // Check if meshPhi was stored, and remove if it is,
    // since sizes probably don't match after a topo-change.
    scalar t0 = mesh_.time().value() - mesh_.time().deltaT().value();

    IOobject meshPhiHeader
    (
        "meshPhi",
        mesh_.time().timeName(t0),
        mesh_,
        IOobject::NO_READ
    );

    // Wipe it out
    if (meshPhiHeader.headerOk())
    {
        rm(mesh_.time().path()/mesh_.time().timeName(t0)/"meshPhi");
    }

    // Set a new empty meshPhi.
    surfaceScalarField& mapPhi = mesh_.setPhi();

    // Update mesh fluxes
    forAll(meshPhi, faceI)
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
            mapPhi.internalField()[newIndex] = meshPhi[faceI];
        }
        else
        {
            label i = mpm.mesh().boundaryMesh()[patch].whichFace(newIndex);

            mapPhi.boundaryField()[patch][i] = meshPhi[faceI];
        }
    }

    // Clear-out demand-driven data.
    clearOut();
}

} // End namespace Foam

// ************************************************************************* //
