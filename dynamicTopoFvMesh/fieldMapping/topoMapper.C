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
#include "volFields.H"

// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * *  //

Foam::topoMapper::~topoMapper()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return reference to the mesh
const Foam::fvMesh&
Foam::topoMapper::mesh() const
{
    return mesh_;
}


//- Return reference to objectRegistry storing fields.
const Foam::objectRegistry&
Foam::topoMapper::db() const
{
    return mesh_;
}


//- Set mapping information
void Foam::topoMapper::setMapper(const mapPolyMesh& mpm)
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
    cellMap_.set(new topoCellMapper(mpm));
    surfaceMap_.set(new fvSurfaceMapper(mesh(), faceMap_()));
    boundaryMap_.set(new fvBoundaryMeshMapper(mesh(), faceMap_()));
}


//- Set face weighting information
void Foam::topoMapper::setFaceWeights
(
    Map<scalarField>& weights
)
{

}


//- Set cell weighting information
void Foam::topoMapper::setCellWeights
(
    Map<scalarField>& weights
)
{

}


//- Set old cell-centre information
void Foam::topoMapper::setOldCellCentres
(
    const volVectorField& oldCentres
)
{
    if (!cellMap_.valid())
    {
        FatalErrorIn
        (
            "void topoMapper::setOldCellCentres()"
        ) << nl << " Cell map has not been set. "
          << abort(FatalError);
    }

    cellMap_().setOldCellCentres(oldCentres);
}


//- Correct fluxes after topology change
void Foam::topoMapper::correctFluxes()
{

}


//- Return volume mapper
const Foam::FieldMapper&
Foam::topoMapper::volMap() const
{
    if (!cellMap_.valid())
    {
        FatalErrorIn
        (
            "const FieldMapper& topoMapper::volMap()"
        ) << nl << " Volume mapper has not been set. "
          << abort(FatalError);
    }

    return cellMap_();
}


//- Return surface mapper
const Foam::fvSurfaceMapper&
Foam::topoMapper::surfaceMap() const
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
const Foam::fvBoundaryMeshMapper&
Foam::topoMapper::boundaryMap() const
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
void Foam::topoMapper::clear()
{
    // Clear out pointers
    faceMap_.clear();
    cellMap_.clear();
    surfaceMap_.clear();
    boundaryMap_.clear();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //
void Foam::topoMapper::operator=(const topoMapper& rhs)
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

// ************************************************************************* //
