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
    mesquitePatchMapper

Description
    Implementation of the mesquitePatchMapper class

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "mapPolyMesh.H"
#include "mesquiteMapper.H"
#include "demandDrivenData.H"
#include "mesquitePatchMapper.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Calculate addressing for interpolative mapping
void mesquitePatchMapper::calcAddressing() const
{
    if (directAddrPtr_ || interpolationAddrPtr_)
    {
        FatalErrorIn("void mesquitePatchMapper::calcAddressing() const")
            << "Addressing already calculated."
            << abort(FatalError);
    }
}


//- Calculate weights for interpolative mapping
void mesquitePatchMapper::calcWeights() const
{
    if (weightsPtr_)
    {
        FatalErrorIn("void mesquitePatchMapper::calcWeights() const")
            << "Weights already calculated."
            << abort(FatalError);
    }
}


//- Clear out local storage
void mesquitePatchMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
    deleteDemandDrivenData(interpolationAddrPtr_);
    deleteDemandDrivenData(weightsPtr_);
    deleteDemandDrivenData(insertedPointLabelsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
mesquitePatchMapper::mesquitePatchMapper
(
    const pointPatch& patch,
    const mapPolyMesh& mpm,
    const mesquiteMapper& mapper
)
:
    patch_(patch),
    mpm_(mpm),
    tMapper_(mapper),
    direct_(false),
    sizeBeforeMapping_
    (
        patch_.index() < mpm.oldPatchNMeshPoints().size()
      ? mpm.oldPatchNMeshPoints()[patch_.index()] : 0
    ),
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL),
    insertedPointLabelsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mesquitePatchMapper::~mesquitePatchMapper()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return size
label mesquitePatchMapper::size() const
{
    return patch_.size();
}


//- Return size before mapping
label mesquitePatchMapper::sizeBeforeMapping() const
{
    return sizeBeforeMapping_;
}


//- Is the mapping direct
bool mesquitePatchMapper::direct() const
{
    return direct_;
}


//- Has unmapped elements
bool mesquitePatchMapper::hasUnmapped() const
{
    return false;
}


//- Return direct addressing
const unallocLabelList& mesquitePatchMapper::directAddressing() const
{
    if (!direct())
    {
        FatalErrorIn
        (
            "const unallocLabelList& "
            "mesquitePatchMapper::directAddressing() const"
        )   << "Requested direct addressing for an interpolative mapper."
            << abort(FatalError);
    }

    if (!directAddrPtr_)
    {
        calcAddressing();
    }

    return *directAddrPtr_;
}


//- Return interpolation addressing
const labelListList& mesquitePatchMapper::addressing() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const labelListList& "
            "mesquitePatchMapper::addressing() const"
        )   << "Requested interpolative addressing for a direct mapper."
            << abort(FatalError);
    }

    if (!interpolationAddrPtr_)
    {
        calcAddressing();
    }

    return *interpolationAddrPtr_;
}


//- Return weights
const scalarListList& mesquitePatchMapper::weights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const scalarListList& "
            "mesquitePatchMapper::weights() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!weightsPtr_)
    {
        calcWeights();
    }

    return *weightsPtr_;
}


//- Are there any inserted cells
bool mesquitePatchMapper::insertedObjects() const
{
    return insertedObjectLabels().size();
}


//- Return list of inserted cells
const labelList& mesquitePatchMapper::insertedObjectLabels() const
{
    if (!insertedPointLabelsPtr_)
    {
        calcAddressing();
    }

    return *insertedPointLabelsPtr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void mesquitePatchMapper::operator=(const mesquitePatchMapper& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("void mesquitePatchMapper::operator=")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}


} // End namespace Foam

// ************************************************************************* //
