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
#include "fvCFD.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(interpolator,0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

interpolator::interpolator(dynamicTopoFvMesh& mesh)
:
    mesh_(mesh),
    fieldsRegistered_(false),
    uName_("U"),
    phiName_("phi")
{
    // Check if alternatives names for velocity and flux are defined
    const dictionary& meshDict = mesh.dict_.subDict("dynamicTopoFvMesh");

    if (meshDict.found("uName"))
    {
        uName_ = word(meshDict.lookup("uName"));
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

// Correct fluxes after topo-changes
void interpolator::correctFluxes()
{
    bool correctPhi = false;

    // Fetch the solution dictionary entries
    dictionary piso = mesh_.solutionDict().subDict("PISO");

    if (piso.found("correctPhi"))
    {
        correctPhi = Switch(piso.lookup("correctPhi"));
    }

    if (!correctPhi)
    {
        return;
    }

    // Fetch references
    dynamicTopoFvMesh& mesh = mesh_;
    const Time& runTime = mesh.time();

    label pRefCell = 0;
    scalar pRefValue = 0.0;

    int nNonOrthCorr = 0;

    if (piso.found("nNonOrthogonalCorrectors"))
    {
        nNonOrthCorr = readInt(piso.lookup("nNonOrthogonalCorrectors"));
    }

    word pName("p"), rUAName("rUA");

    // Check if alternatives names are defined
    const dictionary& meshDict = mesh.dict_.subDict("dynamicTopoFvMesh");

    if (meshDict.found("pName"))
    {
        pName = word(meshDict.lookup("pName"));
    }

    if (meshDict.found("rUAName"))
    {
        rUAName = word(meshDict.lookup("rUAName"));
    }

    // Fetch fields from the registry
    surfaceScalarField& phi = const_cast<surfaceScalarField&>
    (
        mesh_.lookupObject<surfaceScalarField>(phiName_)
    );

    volScalarField& p = const_cast<volScalarField&>
    (
        mesh_.lookupObject<volScalarField>(pName)
    );

    volScalarField& rUA = const_cast<volScalarField&>
    (
        mesh_.lookupObject<volScalarField>(rUAName)
    );

#   include "correctPhi.H"
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

    if (surfScalarMap_.found(phiName_))
    {
        surfaceScalarField& phi = const_cast<surfaceScalarField&>
        (
            mesh_.lookupObject<surfaceScalarField>(phiName_)
        );

        volVectorField& U = const_cast<volVectorField&>
        (
            mesh_.lookupObject<volVectorField>(uName_)
        );

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

            // Flip flux for the internal face.
            phi.internalField()[newIndex] *= -1.0;
        }

        // Interpolate fluxes for new / modified faces.
        surfaceScalarField phiU = fvc::interpolate(U) & mesh_.Sf();

        // Correct fluxes, if necessary
        correctFluxes();
    }

    // Clear-out demand-driven data.
    clearOut();
}

} // End namespace Foam

// ************************************************************************* //
