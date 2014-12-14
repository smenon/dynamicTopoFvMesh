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
    mesquiteInternalMapper

Description
    Implementation of the mesquiteInternalMapper class

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "mapPolyMesh.H"
#include "mesquiteMapper.H"
#include "demandDrivenData.H"
#include "mesquiteInternalMapper.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Calculate addressing for interpolative mapping
void mesquiteInternalMapper::calcAddressing() const
{
    if (directAddrPtr_ || interpolationAddrPtr_)
    {
        FatalErrorIn("void mesquiteInternalMapper::calcAddressing() const")
            << "Addressing already calculated."
            << abort(FatalError);
    }
}


//- Calculate weights for interpolative mapping
void mesquiteInternalMapper::calcWeights() const
{
    if (weightsPtr_)
    {
        FatalErrorIn("void mesquiteInternalMapper::calcWeights() const")
            << "Weights already calculated."
            << abort(FatalError);
    }
}


//- Clear out local storage
void mesquiteInternalMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
    deleteDemandDrivenData(interpolationAddrPtr_);
    deleteDemandDrivenData(weightsPtr_);
    deleteDemandDrivenData(insertedPointLabelsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
mesquiteInternalMapper::mesquiteInternalMapper
(
    const pointMesh& mesh,
    const mapPolyMesh& mpm,
    const mesquiteMapper& mapper
)
:
    mpm_(mpm),
    tMapper_(mapper),
    direct_(false),
    sizeBeforeMapping_(mpm.nOldPoints()),
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL),
    insertedPointLabelsPtr_(NULL)
{
    // Check for possibility of direct mapping
    if ((min(mpm_.pointMap()) > -1) && mpm_.pointsFromPointsMap().empty())
    {
        direct_ = true;
    }
    else
    {
        direct_ = false;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mesquiteInternalMapper::~mesquiteInternalMapper()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return size
label mesquiteInternalMapper::size() const
{
    return mpm_.pointMap().size();
}


//- Return size before mapping
label mesquiteInternalMapper::sizeBeforeMapping() const
{
    return sizeBeforeMapping_;
}


//- Is the mapping direct
bool mesquiteInternalMapper::direct() const
{
    return direct_;
}


//- Has unmapped elements
bool mesquiteInternalMapper::hasUnmapped() const
{
    return false;
}


//- Return direct addressing
const unallocLabelList& mesquiteInternalMapper::directAddressing() const
{
    if (!direct())
    {
        FatalErrorIn
        (
            "const unallocLabelList& "
            "mesquiteInternalMapper::directAddressing() const"
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
const labelListList& mesquiteInternalMapper::addressing() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const labelListList& "
            "mesquiteInternalMapper::addressing() const"
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
const scalarListList& mesquiteInternalMapper::weights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const scalarListList& "
            "mesquiteInternalMapper::weights() const"
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
bool mesquiteInternalMapper::insertedObjects() const
{
    return insertedObjectLabels().size();
}


//- Return list of inserted cells
const labelList& mesquiteInternalMapper::insertedObjectLabels() const
{
    if (!insertedPointLabelsPtr_)
    {
        calcAddressing();
    }

    return *insertedPointLabelsPtr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void mesquiteInternalMapper::operator=(const mesquiteInternalMapper& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("void mesquiteInternalMapper::operator=")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}


} // End namespace Foam

// ************************************************************************* //
