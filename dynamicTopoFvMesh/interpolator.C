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

defineTypeNameAndDebug(interpolator,0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

interpolator::interpolator(dynamicTopoFvMesh& mesh)
:
    mesh_(mesh),
    fieldsRegistered_(false),
    Uname_("U"),
    phiName_("phi")
{
    // Check if alternatives names for velocity and flux are defined
    const dictionary& meshDict = mesh.dict_.subDict("dynamicTopoFvMesh");

    if (meshDict.found("Uname"))
    {
        Uname_ = word(meshDict.lookup("Uname"));
    }

    if (meshDict.found("phiName"))
    {
        phiName_ = word(meshDict.lookup("phiName"));
    }
}

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

    flipFaces_.clear();
}

//- Post-processing
void interpolator::writeFluxes
(
    const word& name
)
{
    // Fetch phi from the registry.
    const surfaceScalarField& pF =
    (
        mesh_.lookupObject<surfaceScalarField>(phiName_)
    );

    // Fetch face-centres from the mesh.
    const vectorField& fC = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();

    const polyBoundaryMesh& boundary = mesh_.boundaryMesh();

    // Make the directory
    fileName dirName(mesh_.time().path()/mesh_.time().timeName());
    //fileName dirName(mesh_.time().path()/"VTK"/mesh_.time().timeName());

    mkDir(dirName);

    // Open stream for output
    OFstream file(dirName/name + ".vtk");

    label numFaces = fC.size();

    forAll(boundary, patchI)
    {
        if
        (
            (boundary[patchI].type() == "wedge") ||
            (boundary[patchI].type() == "empty")
        )
        {
            numFaces -= boundary[patchI].size();
        }
    }

    // Write out the header
    file << "# vtk DataFile Version 2.0" << nl
         << name << ".vtk" << nl
         << "ASCII" << nl
         << "DATASET UNSTRUCTURED_GRID" << nl
         << "POINTS " << numFaces << " double" << nl;

    forAll(fC, faceI)
    {
        label patch = boundary.whichPatch(faceI);

        if (patch == -1)
        {
            file << fC[faceI].x() << ' '
                 << fC[faceI].y() << ' '
                 << fC[faceI].z() << ' '
                 << nl;
        }
        else
        if
        (
            (boundary[patch].type() != "wedge") &&
            (boundary[patch].type() != "empty")
        )
        {
            file << fC[faceI].x() << ' '
                 << fC[faceI].y() << ' '
                 << fC[faceI].z() << ' '
                 << nl;
        }
    }

    file << "CELLS " << numFaces << " " << (2 * numFaces) << endl;

    for(label i = 0; i < numFaces; i++)
    {
        file << 1 << ' ' << i << nl;
    }

    file << "CELL_TYPES " << numFaces << endl;

    for(label i = 0; i < numFaces; i++)
    {
        file << 1 << nl;
    }

    file << "POINT_DATA " << numFaces << endl;

    file << "FIELD PointFields 1" << endl;

    file << "Fluxes 3 " << numFaces << " double" << endl;

    forAll(fC, faceI)
    {
        vector n = (Sf[faceI]/(mag(Sf[faceI]) + VSMALL)), v = vector::zero;

        label patch = boundary.whichPatch(faceI);

        if (patch == -1)
        {
            v = (pF.internalField()[faceI]*n);
        }
        else
        {
            if
            (
                (boundary[patch].type() == "wedge") ||
                (boundary[patch].type() == "empty")
            )
            {
                continue;
            }

            label i = boundary[patch].whichFace(faceI);

            v = (pF.boundaryField()[patch][i]*n);
        }

        file << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
    }
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

    // Map volBoundaryFaces for each primitive type
    mapVolBoundaryFace(scalar, Scalar);
    mapVolBoundaryFace(vector, Vector);
    mapVolBoundaryFace(sphericalTensor, SphericalTensor);
    mapVolBoundaryFace(symmTensor, SymmTensor);
    mapVolBoundaryFace(tensor, Tensor);
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

    if (flipFaces_.found(index))
    {
        flipFaces_.erase(index);
    }
}

// Add a new cell entry.
void interpolator::insertCell
(
    const label newCellIndex,
    const labelList& mapCells,
    const scalarField& cellWeights
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

// Get the volume-flux for an existing face
scalar interpolator::getPhi
(
    const label faceIndex
)
{
    scalar phiVal = 0.0;

    if (surfScalarMap_.found(phiName_))
    {
        if
        (
            (faceIndex < mesh_.nOldFaces_) &&
            (!surfScalarMap_[phiName_].found(faceIndex))
        )
        {
            // Fetch phi from the registry.
            const surfaceScalarField& oF =
            (
                mesh_.lookupObject<surfaceScalarField>(phiName_)
            );

            const polyBoundaryMesh& boundary = mesh_.boundaryMesh();

            label oPatch = boundary.whichPatch(faceIndex);

            if (oPatch == -1)
            {
                phiVal = oF.internalField()[faceIndex];
            }
            else
            {
                label i = boundary[oPatch].whichFace(faceIndex);

                phiVal = oF.boundaryField()[oPatch][i];
            }
        }
        else
        {
            phiVal = surfScalarMap_[phiName_][faceIndex];
        }
    }

    // Check if the flux needs to be flipped.
    if (flipFaces_.found(faceIndex))
    {
        phiVal *= -1.0;
    }

    return phiVal;
}

// Set the volume-flux for an existing face
void interpolator::setPhi
(
    const label faceIndex,
    const scalar facePhi
)
{
    if (surfScalarMap_.found(phiName_))
    {
        surfScalarMap_[phiName_].set(faceIndex, facePhi);
    }
}

// Interpolate flux for an existing face
void interpolator::interpolatePhi
(
    const label owner,
    const label neighbour,
    const label faceIndex,
    const vector& Sf
)
{
    FixedList<label, 2> c(-1);
    FixedList<vector,2> U(vector::zero);
    FixedList<scalar, 2> w(0.0);

    if (neighbour == -1)
    {
        w[0] = 1.0;
        w[1] = 0.0;
    }
    else
    {
        w[0] = 0.5;
        w[1] = 0.5;
    }

    // Set check indices
    c[0] = owner; c[1] = neighbour;

    // Check volVectorMaps for an entry named 'U'.
    // If it exists, make an entry (or over-write, if necessary).
    if (volVectorMap_.found(Uname_))
    {
        // Looks like an entry exists.
        // Check if a new value was set for this cell.
        forAll(c, indexI)
        {
            if (c[indexI] == -1)
            {
                continue;
            }

            if (volVectorMap_[Uname_].found(c[indexI]))
            {
                U[indexI] = volVectorMap_[Uname_][c[indexI]];
            }
            else
            if (c[indexI] < mesh_.nOldCells_)
            {
                // Fetch old value from the registry.
                const volVectorField& oF =
                (
                    mesh_.lookupObject<volVectorField>(Uname_)
                );

                U[indexI] = oF.internalField()[c[indexI]];
            }
            else
            {
                // Couldn't find the cell anywhere.
                FatalErrorIn("interpolator::interpolatePhi()")
                    << " Looking for cell: " << c[indexI]
                    << " in maps, but couldn't find it." << nl
                    << " faceIndex: " << faceIndex
                    << " Sf: " << Sf
                    << abort(FatalError);
            }
        }

        if (dynamicTopoFvMesh::debug > 3)
        {
            Info << nl
                 << "Flux for face: " << faceIndex
                 << " :: " << ((w[0]*U[0] + w[1]*U[1]) & Sf)
                 << " Using cells: " << c
                 << " with U: " << U
                 << " w: " << w
                 << " and Sf: " << Sf
                 << endl;
        }

        // Take the dot-product and assign.
        setPhi
        (
            faceIndex,
            ((w[0]*U[0] + w[1]*U[1]) & Sf)
        );
    }
}

// Set a particular face index as flipped.
void interpolator::setFlip(const label fIndex)
{
    if (flipFaces_.found(fIndex))
    {
        if (dynamicTopoFvMesh::debug > 3)
        {
            Info << "UnFlipping face: " << fIndex << endl;
        }

        flipFaces_.erase(fIndex);
    }
    else
    {
        if (dynamicTopoFvMesh::debug > 3)
        {
            Info << "Flipping face: " << fIndex << endl;
        }

        flipFaces_.insert(fIndex);
    }
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
    // Mesh should've been reset when this function is called,
    // so the boundaryMesh contains up-to-date patch information.
    const polyBoundaryMesh& boundary = mpm.mesh().boundaryMesh();

    // Fill internal fields
    fillVolumeMaps(scalar, Scalar);
    fillVolumeMaps(vector, Vector);
    fillVolumeMaps(sphericalTensor, SphericalTensor);
    fillVolumeMaps(symmTensor, SymmTensor);
    fillVolumeMaps(tensor, Tensor);

    // Fill volBoundary fields
    fillVolBoundaryMaps(scalar, Scalar);
    fillVolBoundaryMaps(vector, Vector);
    fillVolBoundaryMaps(sphericalTensor, SphericalTensor);
    fillVolBoundaryMaps(symmTensor, SymmTensor);
    fillVolBoundaryMaps(tensor, Tensor);

    // Fill surface fields
    fillSurfaceMaps(scalar, Scalar);
    fillSurfaceMaps(vector, Vector);
    fillSurfaceMaps(sphericalTensor, SphericalTensor);
    fillSurfaceMaps(symmTensor, SymmTensor);
    fillSurfaceMaps(tensor, Tensor);

    // Flip fluxes for specific internal faces.
    forAllConstIter(labelHashSet, flipFaces_, fIter)
    {
        label newIndex = -1;

        if (fIter.key() < mpm.nOldFaces())
        {
            newIndex = mpm.reverseFaceMap()[fIter.key()];
        }
        else
        {
            newIndex = mesh_.addedFaceRenumbering_[fIter.key()];
        }

        forAllIter(HashTable<Map<scalar> >, surfScalarMap_, hIter)
        {
            if (volScalarMap_.found(hIter.key()))
            {
                continue;
            }

            if (!mesh_.foundObject<surfaceScalarField>(hIter.key()))
            {
                continue;
            }

            surfaceScalarField& sf =
            (
                const_cast<surfaceScalarField&>
                (
                    mesh_.lookupObject<surfaceScalarField>
                    (
                        hIter.key()
                    )
                )
            );

            // Flip fux for the internal face.
            sf.internalField()[newIndex] *= -1.0;
        }
    }

    // Clear-out demand-driven data.
    clearOut();
}

} // End namespace Foam

// ************************************************************************* //
