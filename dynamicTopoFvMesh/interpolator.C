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
#include "GeometricFields.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "fvPatchFieldsFwd.H"
#include "volFields.H"
#include "fvsPatchFieldsFwd.H"
#include "surfaceFields.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(interpolator,1);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

interpolator::interpolator(dynamicTopoFvMesh& mesh)
:
    mesh_(mesh),
    fieldsRegistered_(false)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

interpolator::~interpolator()
{
    clearOut();

    fieldsRegistered_ = false;
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void interpolator::clearOut()
{
    // Clear vol and surf maps (but not hash-table keys)
    clearInterpolationMaps(scalar, Scalar);
    clearInterpolationMaps(vector, Vector);
    clearInterpolationMaps(sphericalTensor, SphericalTensor);
    clearInterpolationMaps(symmTensor, SymmTensor);
    clearInterpolationMaps(tensor, Tensor);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Add a new face entry
void interpolator::insertFace
(
    const label patch,
    const label newFaceIndex,
    const labelList& mapFaces,
    const scalarField& mapWeights
)
{
    // Loop through all volMaps / surfaceMaps, and perform a weighted mapping.
    // If the key exists, overwrite it.
    label nOldFaces = mesh_.nOldFaces_;

    const polyBoundaryMesh& boundary = mesh_.boundaryMesh();

    // Map for each primitive type
    mapInternalFace(scalar, Scalar);
    mapInternalFace(vector, Vector);
    mapInternalFace(sphericalTensor, SphericalTensor);
    mapInternalFace(symmTensor, SymmTensor);
    mapInternalFace(tensor, Tensor);

    // Map boundary volBoundaryFields only if this is a boundary face.
    if (patch == -1)
    {
        return;
    }

    // Map for each primitive type
    mapBoundaryFace(scalar, Scalar);
    mapBoundaryFace(vector, Vector);
    mapBoundaryFace(sphericalTensor, SphericalTensor);
    mapBoundaryFace(symmTensor, SymmTensor);
    mapBoundaryFace(tensor, Tensor);
}

// Remove the face, if it exists in the map.
void interpolator::removeFace
(
    const label index
)
{
    // Loop through all surfMaps
    removeEntity(scalar, surfScalar);
    removeEntity(vector, surfVector);
    removeEntity(sphericalTensor, surfSphericalTensor);
    removeEntity(symmTensor, surfSymmTensor);
    removeEntity(tensor, surfTensor);
}

// Add a new cell entry.
void interpolator::insertCell
(
    const label newCellIndex,
    const labelList& mapCells,
    const scalarField& mapWeights
)
{
    // Loop through all volMaps, and perform a weighted mapping.
    // If the key exists, overwrite it.
    label nOldCells = mesh_.nOldCells_;

    // Map for each primitive type
    mapCell(scalar, Scalar);
    mapCell(vector, Vector);
    mapCell(sphericalTensor, SphericalTensor);
    mapCell(symmTensor, SymmTensor);
    mapCell(tensor, Tensor);
}

// Remove the cell, if it exists in the map.
void interpolator::removeCell
(
    const label index
)
{
    // Loop through all volMaps
    removeEntity(scalar, volScalar);
    removeEntity(vector, volVector);
    removeEntity(sphericalTensor, volSphericalTensor);
    removeEntity(symmTensor, volSymmTensor);
    removeEntity(tensor, volTensor);
}

// Set the volume-flux for an existing face
void interpolator::setPhi
(
    const label faceIndex,
    const scalar facePhi
)
{

}

// Register fields for interpolation
void interpolator::registerFields()
{
    // If fields have already been registered, bail out.
    if (fieldsRegistered_)
    {
        return;
    }

    // Fetch all volFields and surfaceFields from the registry.

    // volFields need surface-mapping as well (for boundaries),
    // so add them to surface maps.
    registerVolumeField(scalar, Scalar);
    registerVolumeField(vector, Vector);
    registerVolumeField(sphericalTensor, SphericalTensor);
    registerVolumeField(symmTensor, SymmTensor);
    registerVolumeField(tensor, Tensor);

    registerSurfaceField(scalar, Scalar);
    registerSurfaceField(vector, Vector);
    registerSurfaceField(sphericalTensor, SphericalTensor);
    registerSurfaceField(symmTensor, SymmTensor);
    registerSurfaceField(tensor, Tensor);

    // Set the flag.
    fieldsRegistered_ = true;
}

// Update fields after a topo-change operation
void interpolator::updateMesh(const mapPolyMesh& mpm)
{
    // First re-number maps after re-ordering

    // Now loop through maps, and fill-in values
    // for inserted elements.

    // Clear-out demand-driven data.
    clearOut();
}

} // End namespace Foam

// ************************************************************************* //
