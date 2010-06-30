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
#include "topoCellMapper.H"
#include "faceMapper.H"
#include "fvSurfaceMapper.H"
#include "fvBoundaryMeshMapper.H"
#include "volFields.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Store gradients of fields on the mesh prior to topology changes
template <class Type, class gradType>
void topoMapper::storeGradients
(
    HashTable<autoPtr<gradType> >& gradTable
)
{
    // Fetch all fields from registry
    HashTable<const Type*> fields(mesh_.objectRegistry::lookupClass<Type>());

    forAllConstIter(typename HashTable<const Type*>, fields, fIter)
    {
        const Type& field = *fIter();

        // Compute the gradient.
        tmp<gradType> tGrad = fvc::grad(field);

        // Make a new entry, but don't register the field.
        gradTable.insert
        (
            field.name(),
            autoPtr<gradType>
            (
                new gradType
                (
                    IOobject
                    (
                        tGrad().name(),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    tGrad()
                )
            )
        );
    }
}


void topoMapper::clearOut()
{
    deleteDemandDrivenData(oldCellCentresPtr_);

    // Clear stored gradients
    sGrads_.clear();
    vGrads_.clear();

    // Clear maps
    faceWeights_.clear();
    cellWeights_.clear();

    faceCentres_.clear();
    cellCentres_.clear();
}

// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * *  //

topoMapper::~topoMapper()
{
    clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return reference to the mesh
const fvMesh&
topoMapper::mesh() const
{
    return mesh_;
}


//- Return reference to objectRegistry storing fields.
const objectRegistry&
topoMapper::db() const
{
    return mesh_;
}


//- Return mapping method
label topoMapper::method
(
    const word& typeName
) const
{
    return -1;
}


//- Set mapping information
void topoMapper::setMapper(const mapPolyMesh& mpm)
{
    if
    (
        faceMap_.valid() ||
        cellMap_.valid() ||
        surfaceMap_.valid() ||
        boundaryMap_.valid()
    )
    {
        FatalErrorIn
        (
            "void topoMapper::setMapper()"
        ) << nl << " Mapper has already been set. "
          << abort(FatalError);
    }

    // Set pointers
    faceMap_.set(new faceMapper(mpm));
    cellMap_.set(new topoCellMapper(mpm, *this));
    surfaceMap_.set(new fvSurfaceMapper(mesh(), faceMap_()));
    boundaryMap_.set(new fvBoundaryMeshMapper(mesh(), faceMap_()));
}


//- Set face weighting information
void topoMapper::setFaceWeights
(
    Map<scalarField>& weights,
    Map<vectorField>& centres
)
{
    faceWeights_.transfer(weights);
    faceCentres_.transfer(centres);
}


//- Set cell weighting information
void topoMapper::setCellWeights
(
    Map<scalarField>& weights,
    Map<vectorField>& centres
)
{
    cellWeights_.transfer(weights);
    cellCentres_.transfer(centres);
}


//- Fetch face weights
const Map<scalarField>&
topoMapper::faceWeights() const
{
    return faceWeights_;
}


//- Fetch cell weights
const Map<scalarField>&
topoMapper::cellWeights() const
{
    return cellWeights_;
}


//- Fetch face centres
const Map<vectorField>&
topoMapper::faceCentres() const
{
    return faceCentres_;
}


//- Fetch cell centres
const Map<vectorField>&
topoMapper::cellCentres() const
{
    return cellCentres_;
}


//- Store gradients prior to mesh reset
void topoMapper::storeGradients()
{
    storeGradients<volScalarField>(sGrads_);
    storeGradients<volVectorField>(vGrads_);
}


//- Set old cell-centre information
void topoMapper::setOldCellCentres
(
    const volVectorField& oldCentres
)
{
    if (oldCellCentresPtr_)
    {
        deleteDemandDrivenData(oldCellCentresPtr_);
    }

    // Set the pointer.
    // Only copy values, but don't register the field,
    // since we don't want it to be mapped like the others
    oldCellCentresPtr_ =
    (
        new volVectorField
        (
            IOobject
            (
                "OldCellCentres",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            oldCentres
        )
    );
}


//- Return old cell-centre information
const volVectorField&
topoMapper::oldCentres() const
{
    if (!oldCellCentresPtr_)
    {
        FatalErrorIn
        (
            "void topoMapper::oldCentres()"
        ) << nl << " Pointer has not been set. "
          << abort(FatalError);
    }

    return *oldCellCentresPtr_;
}


//- Fetch the gradient field (template specialisation)
template <>
const volVectorField&
topoMapper::gradient(const word& name) const
{
    return sGrads_[name]();
}


//- Fetch the gradient field (template specialisation)
template <>
const volTensorField&
topoMapper::gradient(const word& name) const
{
    return vGrads_[name]();
}


//- Correct fluxes after topology change
void topoMapper::correctFluxes()
{

}


//- Return volume mapper
const topoCellMapper&
topoMapper::volMap() const
{
    if (!cellMap_.valid())
    {
        FatalErrorIn
        (
            "const topoCellMapper& topoMapper::volMap()"
        ) << nl << " Volume mapper has not been set. "
          << abort(FatalError);
    }

    return cellMap_();
}


//- Return surface mapper
const fvSurfaceMapper&
topoMapper::surfaceMap() const
{
    if (!surfaceMap_.valid())
    {
        FatalErrorIn
        (
            "const fvSurfaceMapper& topoMapper::surfaceMap()"
        ) << nl << " Surface mapper has not been set. "
          << abort(FatalError);
    }

    return surfaceMap_();
}


//- Return boundary mapper
const fvBoundaryMeshMapper&
topoMapper::boundaryMap() const
{
    if (!boundaryMap_.valid())
    {
        FatalErrorIn
        (
            "const fvBoundaryMapper& topoMapper::boundaryMap()"
        ) << nl << " Boundary mapper has not been set. "
          << abort(FatalError);
    }

    return boundaryMap_();
}


//- Clear out member data
void topoMapper::clear()
{
    // Clear out pointers
    faceMap_.clear();
    cellMap_.clear();
    surfaceMap_.clear();
    boundaryMap_.clear();
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
