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
    topoMapper

Description
    Implementation of topoMapper

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*----------------------------------------------------------------------------*/

#include "topoMapper.H"
#include "fluxCorrector.H"
#include "topoCellMapper.H"
#include "topoPointMapper.H"
#include "topoSurfaceMapper.H"
#include "topoBoundaryMeshMapper.H"
#include "topoPointBoundaryMapper.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

List<pointScalarField*> topoMapper::psFieldPtrs_;
List<pointVectorField*> topoMapper::pvFieldPtrs_;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Store gradients prior to mesh reset
void topoMapper::storeGradients() const
{
    // Go in order of highest to lowest rank,
    // to avoid double storage
    storeGradients<vector>(gradTable_, vGradPtrs_);
    storeGradients<scalar>(gradTable_, sGradPtrs_);

    if (fvMesh::debug)
    {
        Info<< "Registered gradients: " << gradientTable() << endl;
    }
}


//- Store geometric information
void topoMapper::storeGeometry() const
{
    typedef volVectorField::PatchFieldType PatchFieldType;
    typedef volVectorField::GeometricBoundaryField GeomBdyFieldType;
    typedef volVectorField::DimensionedInternalField DimInternalField;

    // Wipe out existing information
    deleteDemandDrivenData(cellCentresPtr_);

    vectorField Cv(mesh_.cellCentres());
    vectorField Cf(mesh_.faceCentres());

    // Create and map the patch field values
    label nPatches = mesh_.boundary().size();

    // Create field parts
    PtrList<PatchFieldType> volCentrePatches(nPatches);

    // Define patch type names
    word emptyType("empty");
    word fixedValueType("fixedValue");

    // Create dummy types for initial field creation
    forAll(volCentrePatches, patchI)
    {
        if (mesh_.boundary()[patchI].type() == emptyType)
        {
            volCentrePatches.set
            (
                patchI,
                PatchFieldType::New
                (
                    emptyType,
                    mesh_.boundary()[patchI],
                    DimInternalField::null()
                )
            );
        }
        else
        {
            volCentrePatches.set
            (
                patchI,
                PatchFieldType::New
                (
                    fixedValueType,
                    mesh_.boundary()[patchI],
                    DimInternalField::null()
                )
            );
        }
    }

    // Set the cell-centres pointer.
    cellCentresPtr_ =
    (
        new volVectorField
        (
            IOobject
            (
                "cellCentres",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true
            ),
            mesh_,
            dimLength,
            SubField<vector>(Cv, mesh_.nCells()),
            volCentrePatches
        )
    );

    // Alias for convenience
    volVectorField& centres = *cellCentresPtr_;

    // Set correct references for patch internal fields
    GeomBdyFieldType& bf = centres.boundaryField();

    forAll(bf, patchI)
    {
        if (mesh_.boundary()[patchI].type() == emptyType)
        {
            bf.set
            (
                patchI,
                PatchFieldType::New
                (
                    emptyType,
                    mesh_.boundary()[patchI],
                    centres.dimensionedInternalField()
                )
            );
        }
        else
        {
            bf.set
            (
                patchI,
                PatchFieldType::New
                (
                    fixedValueType,
                    mesh_.boundary()[patchI],
                    centres.dimensionedInternalField()
                )
            );

            // Slice field to patch (forced assignment)
            bf[patchI] == mesh_.boundaryMesh()[patchI].patchSlice(Cf);
        }
    }

    // Set the cell-volumes pointer
    cellVolumesPtr_ = new scalarField(mesh_.cellVolumes());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from mesh and dictionary
topoMapper::topoMapper
(
    const fvMesh& mesh,
    const dictionary& dict,
    const bool disableGradients
)
:
    mesh_(mesh),
    dict_(dict),
    cellMap_(NULL),
    pointMap_(NULL),
    surfaceMap_(NULL),
    boundaryMap_(NULL),
    pointBoundaryMap_(NULL),
    fluxCorrector_(fluxCorrector::New(mesh, dict)),
    cellVolumesPtr_(NULL),
    cellCentresPtr_(NULL),
    disableGradients_(disableGradients)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * *  //

topoMapper::~topoMapper()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- De-register all pointFields
void topoMapper::deregisterPointFields(const objectRegistry& registry)
{
    topoMapper::deregisterPointFields(registry, psFieldPtrs_);
    topoMapper::deregisterPointFields(registry, pvFieldPtrs_);
}


//- Re-register all pointFields
void topoMapper::reregisterPointFields(const objectRegistry& registry)
{
    topoMapper::reregisterPointFields(registry, psFieldPtrs_);
    topoMapper::reregisterPointFields(registry, pvFieldPtrs_);
}


//- Return reference to the mesh
const fvMesh& topoMapper::mesh() const
{
    return mesh_;
}


//- Return reference to objectRegistry storing fields.
const objectRegistry& topoMapper::thisDb() const
{
    return mesh_;
}


//- Set mapping information
void topoMapper::setMapper(const mapPolyMesh& mpm) const
{
    if
    (
        cellMap_.valid() ||
        pointMap_.valid() ||
        surfaceMap_.valid() ||
        boundaryMap_.valid() ||
        pointBoundaryMap_.valid()
    )
    {
        FatalErrorIn
        (
            "void topoMapper::setMapper() const"
        ) << nl << " Mapper has already been set. "
          << abort(FatalError);
    }

    // Set pointers
    cellMap_.set(new topoCellMapper(mpm, *this));
    surfaceMap_.set(new topoSurfaceMapper(mpm, *this));
    boundaryMap_.set(new topoBoundaryMeshMapper(mesh(), mpm, *this));

    const pointMesh& pMesh = pointMesh::New(mpm.mesh());

    // Set point mappers
    pointMap_.set(new topoPointMapper(pMesh, mpm, *this));
    pointBoundaryMap_.set(new topoPointBoundaryMapper(pMesh, mpm, *this));
}


//- Set point weighting information
void topoMapper::setPointWeights
(
    const Xfer<scalarFieldList>& weights,
    const Xfer<vectorFieldList>& centres
) const
{
    pointWeights_.transfer(weights());
    pointCentres_.transfer(centres());
}


//- Set face weighting information
void topoMapper::setFaceWeights
(
    const Xfer<scalarFieldList>& weights,
    const Xfer<vectorFieldList>& centres
) const
{
    faceWeights_.transfer(weights());
    faceCentres_.transfer(centres());
}


//- Set cell weighting information
void topoMapper::setCellWeights
(
    const Xfer<scalarFieldList>& weights,
    const Xfer<vectorFieldList>& centres
) const
{
    cellWeights_.transfer(weights());
    cellCentres_.transfer(centres());
}


//- Set old patch mesh points information
void topoMapper::setOldPatchMeshPoints
(
    const Xfer<labelListList>& patchMeshPoints
) const
{
    oldPatchMeshPoints_.transfer(patchMeshPoints());
}


//- Set sub mesh map points information
void topoMapper::setSubMeshMapPointList
(
    const Xfer<MapPointList>& subMeshPoints
) const
{
    subMeshMapPointList_.transfer(subMeshPoints());
}


//- Set sub mesh patch map information
void topoMapper::setSubMeshPatchMaps
(
    const Xfer<labelMapLists>& subMeshPatchMaps
) const
{
    subMeshPatchMaps_.transfer(subMeshPatchMaps());
}


//- Renumber map points after re-ordering
void topoMapper::renumberMapPoints(const labelMap& map) const
{
    forAll(subMeshMapPointList_, mpI)
    {
        MapPoint& mp = subMeshMapPointList_[mpI];

        mp.first() = map[mp.first()];
    }
}


//- Set cell / patch offset information
void topoMapper::setOffsets
(
    const labelList& cellSizes,
    const labelList& cellStarts,
    const labelList& faceSizes,
    const labelList& faceStarts,
    const labelList& pointSizes,
    const labelList& pointStarts,
    const labelListList& patchSizes,
    const labelListList& patchStarts,
    const labelListList& pointPatchSizes,
    const labelListList& pointPatchStarts
) const
{
    cellSizes_ = cellSizes;
    cellStarts_ = cellStarts;
    faceSizes_ = faceSizes;
    faceStarts_ = faceStarts;
    pointSizes_ = pointSizes;
    pointStarts_ = pointStarts;
    patchSizes_ = patchSizes;
    patchStarts_ = patchStarts;
    pointPatchSizes_ = pointPatchSizes;
    pointPatchStarts_ = pointPatchStarts;
}


//- Fetch point weights
const topoMapper::scalarFieldList& topoMapper::pointWeights() const
{
    return pointWeights_;
}


//- Fetch face weights
const topoMapper::scalarFieldList& topoMapper::faceWeights() const
{
    return faceWeights_;
}


//- Fetch cell weights
const topoMapper::scalarFieldList& topoMapper::cellWeights() const
{
    return cellWeights_;
}


//- Fetch point centres
const topoMapper::vectorFieldList& topoMapper::pointCentres() const
{
    return pointCentres_;
}


//- Fetch face centres
const topoMapper::vectorFieldList& topoMapper::faceCentres() const
{
    return faceCentres_;
}


//- Fetch cell centres
const topoMapper::vectorFieldList& topoMapper::cellCentres() const
{
    return cellCentres_;
}


//- Fetch old patch mesh points information
const labelListList& topoMapper::oldPatchMeshPoints() const
{
    return oldPatchMeshPoints_;
}


//- Fetch cell sizes
const labelList& topoMapper::cellSizes() const
{
    return cellSizes_;
}


//- Fetch face sizes
const labelList& topoMapper::faceSizes() const
{
    return faceSizes_;
}


//- Fetch point sizes
const labelList& topoMapper::pointSizes() const
{
    return pointSizes_;
}


//- Fetch patch sizes
const labelListList& topoMapper::patchSizes() const
{
    return patchSizes_;
}


//- Fetch point patch sizes
const labelListList& topoMapper::pointPatchSizes() const
{
    return pointPatchSizes_;
}


//- Fetch cell starts
const labelList& topoMapper::cellStarts() const
{
    return cellStarts_;
}


//- Fetch face starts
const labelList& topoMapper::faceStarts() const
{
    return faceStarts_;
}


//- Fetch point starts
const labelList& topoMapper::pointStarts() const
{
    return pointStarts_;
}


//- Fetch patch starts
const labelListList& topoMapper::patchStarts() const
{
    return patchStarts_;
}


//- Fetch point patch starts
const labelListList& topoMapper::pointPatchStarts() const
{
    return pointPatchStarts_;
}


//- Fetch subMesh map point list
const topoMapper::MapPointList& topoMapper::subMeshMapPointList() const
{
    return subMeshMapPointList_;
}


//- Fetch subMesh patch maps
const topoMapper::labelMapLists& topoMapper::subMeshPatchMaps() const
{
    return subMeshPatchMaps_;
}


// Store old-time information for all registered fields
void topoMapper::storeOldTimes() const
{
    storeOldTimes<volScalarField>();
    storeOldTimes<volVectorField>();
    storeOldTimes<volSphericalTensorField>();
    storeOldTimes<volSymmTensorField>();
    storeOldTimes<volTensorField>();

    storeOldTimes<surfaceScalarField>();
    storeOldTimes<surfaceVectorField>();
    storeOldTimes<surfaceSphericalTensorField>();
    storeOldTimes<surfaceSymmTensorField>();
    storeOldTimes<surfaceTensorField>();
}


//- Store mesh information for the mapping stage
void topoMapper::storeMeshInformation() const
{
    // Store field-gradients
    storeGradients();

    // Store geometry
    storeGeometry();
}


//- Return non-const access to cell centres
volVectorField& topoMapper::volCentres() const
{
    if (!cellCentresPtr_)
    {
        FatalErrorIn
        (
            "const vectorField& topoMapper::volCentres() const"
        ) << nl << " Pointer has not been set. "
          << abort(FatalError);
    }

    return *cellCentresPtr_;
}


//- Return stored cell centre information
const vectorField& topoMapper::internalCentres() const
{
    if (!cellCentresPtr_)
    {
        FatalErrorIn
        (
            "const vectorField& topoMapper::internalCentres() const"
        ) << nl << " Pointer has not been set. "
          << abort(FatalError);
    }

    return *cellCentresPtr_;
}


//- Return non-const access to cell volumes
scalarField& topoMapper::internalVolumes() const
{
    if (!cellVolumesPtr_)
    {
        FatalErrorIn
        (
            "scalarField& topoMapper::internalVolumes() const"
        ) << nl << " Pointer has not been set. "
          << abort(FatalError);
    }

    return *cellVolumesPtr_;
}


//- Return stored patch centre information
const vectorField& topoMapper::patchCentres(const label i) const
{
    if (!cellCentresPtr_)
    {
        FatalErrorIn
        (
            "const vectorField& topoMapper::patchCentres"
            "(const label i) const"
        ) << nl << " Pointer has not been set. index: " << i
          << abort(FatalError);
    }

    return (*cellCentresPtr_).boundaryField()[i];
}


//- Return names of stored gradients
const wordList topoMapper::gradientTable() const
{
    return gradTable_.toc();
}


//- Fetch the gradient field (template specialisation)
template <>
const volVectorField& topoMapper::gradient(const word& name) const
{
    if (!gradTable_.found(name))
    {
        FatalErrorIn
        (
            "const volVectorField& topoMapper::gradient(const word& name) const"
        ) << nl << " Gradient for: " << name
          << " has not been stored."
          << abort(FatalError);
    }

    return sGradPtrs_[gradTable_[name].second()];
}


//- Fetch the gradient field (template specialisation)
template <>
const volTensorField& topoMapper::gradient(const word& name) const
{
    if (!gradTable_.found(name))
    {
        FatalErrorIn
        (
            "const volTensorField& topoMapper::gradient(const word& name) const"
        ) << nl << " Gradient for: " << name
          << " has not been stored."
          << abort(FatalError);
    }

    return vGradPtrs_[gradTable_[name].second()];
}


//- Fetch the pointField list (template specialisation)
template <>
const List<pointScalarField*>& topoMapper::pointFieldList() const
{
    return psFieldPtrs_;
}


//- Fetch the pointField list (template specialisation)
template <>
const List<pointVectorField*>& topoMapper::pointFieldList() const
{
    return pvFieldPtrs_;
}


//- Deregister gradient fields and centres,
//  but retain for mapping
void topoMapper::deregisterMeshInformation() const
{
    // Check out scalar gradients
    forAll(sGradPtrs_, fieldI)
    {
        mesh_.objectRegistry::checkOut(sGradPtrs_[fieldI]);
    }

    // Check out vector gradients
    forAll(vGradPtrs_, fieldI)
    {
        mesh_.objectRegistry::checkOut(vGradPtrs_[fieldI]);
    }

    // Check out cell centres
    mesh_.objectRegistry::checkOut(*cellCentresPtr_);
}


//- Correct fluxes after topology changes, if required
void topoMapper::correctFluxes() const
{
    if (surfaceFluxCorrector().required())
    {
        // Supply a list of inserted faces for interpolation
        surfaceFluxCorrector().interpolateFluxes
        (
            surfaceMap().insertedObjectLabels()
        );

        // Update fluxes
        surfaceFluxCorrector().updateFluxes();
    }
}


//- Return volume mapper
const topoCellMapper& topoMapper::volMap() const
{
    if (!cellMap_.valid())
    {
        FatalErrorIn
        (
            "const topoCellMapper& topoMapper::volMap() const"
        ) << nl << " Volume mapper has not been set. "
          << abort(FatalError);
    }

    return cellMap_();
}


//- Return point mapper
const topoPointMapper& topoMapper::pointMap() const
{
    if (!pointMap_.valid())
    {
        FatalErrorIn
        (
            "const topoPointMapper& topoMapper::pointMap() const"
        ) << nl << " Point mapper has not been set. "
          << abort(FatalError);
    }

    return pointMap_();
}


//- Return surface mapper
const topoSurfaceMapper& topoMapper::surfaceMap() const
{
    if (!surfaceMap_.valid())
    {
        FatalErrorIn
        (
            "const topoSurfaceMapper& topoMapper::surfaceMap() const"
        ) << nl << " Surface mapper has not been set. "
          << abort(FatalError);
    }

    return surfaceMap_();
}


//- Return boundary mapper
const topoBoundaryMeshMapper& topoMapper::boundaryMap() const
{
    if (!boundaryMap_.valid())
    {
        FatalErrorIn
        (
            "const topoBoundaryMeshMapper& topoMapper::boundaryMap() const"
        ) << nl << " Boundary mapper has not been set. "
          << abort(FatalError);
    }

    return boundaryMap_();
}


//- Return point boundary mapper
const topoPointBoundaryMapper& topoMapper::pointBoundaryMap() const
{
    if (!pointBoundaryMap_.valid())
    {
        FatalErrorIn
        (
            "const topoPointBoundaryMapper& "
            "topoMapper::pointBoundaryMap() const"
        ) << nl << " Point boundary mapper has not been set. "
          << abort(FatalError);
    }

    return pointBoundaryMap_();
}


//- Return flux-corrector
const fluxCorrector& topoMapper::surfaceFluxCorrector() const
{
    if (!fluxCorrector_.valid())
    {
        FatalErrorIn
        (
            "const fluxCorrector& topoMapper::surfaceFluxCorrector() const"
        ) << nl << " fluxCorrector has not been set. "
          << abort(FatalError);
    }

    return fluxCorrector_();
}


//- Clear out member data
void topoMapper::clear() const
{
    // Clear out mappers
    cellMap_.clear();
    pointMap_.clear();
    surfaceMap_.clear();
    boundaryMap_.clear();
    pointBoundaryMap_.clear();

    // Clear index maps
    gradTable_.clear();

    // Clear stored gradients
    sGradPtrs_.clear();
    vGradPtrs_.clear();

    // Wipe out geometry information
    deleteDemandDrivenData(cellVolumesPtr_);
    deleteDemandDrivenData(cellCentresPtr_);

    // Clear maps
    pointWeights_.clear();
    faceWeights_.clear();
    cellWeights_.clear();

    pointCentres_.clear();
    faceCentres_.clear();
    cellCentres_.clear();

    // Clear sizes / offsets
    cellSizes_.clear();
    cellStarts_.clear();

    faceSizes_.clear();
    faceStarts_.clear();

    pointSizes_.clear();
    pointStarts_.clear();

    patchSizes_.clear();
    patchStarts_.clear();

    pointPatchSizes_.clear();
    pointPatchStarts_.clear();

    // Clear old patch mesh points
    oldPatchMeshPoints_.clear();

    // Clear sub mesh information
    subMeshMapPointList_.clear();
    subMeshPatchMaps_.clear();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void topoMapper::operator=(const topoMapper& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "topoMapper::operator=(const topoMapper&)"
        )
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}


} // End namespace Foam

// ************************************************************************* //
