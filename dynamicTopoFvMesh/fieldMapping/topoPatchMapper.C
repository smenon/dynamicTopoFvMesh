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
    topoPatchMapper

Description
    Implementation of the topoPatchMapper class

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*----------------------------------------------------------------------------*/

#include "topoMapper.H"
#include "mapPolyMesh.H"
#include "topoPatchMapper.H"
#include "demandDrivenData.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Clear out local storage
void topoPatchMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
    deleteDemandDrivenData(interpolationAddrPtr_);
    deleteDemandDrivenData(weightsPtr_);
}


//- Calculate addressing for mapping
void topoPatchMapper::calcAddressing() const
{

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
topoPatchMapper::topoPatchMapper
(
    const polyPatch& patch,
    const mapPolyMesh& mpm,
    const topoMapper& mapper
)
:
    patch_(patch),
    mpm_(mpm),
    tMapper_(mapper),
    direct_(false),
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

topoPatchMapper::~topoPatchMapper()
{
    clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return size
label topoPatchMapper::size() const
{
    return patch_.size();
}


//- Return size before mapping
label topoPatchMapper::sizeBeforeMapping() const
{
    return mpm_.oldPatchSizes()[patch_.index()];
}


//- Is the mapping direct
bool topoPatchMapper::direct() const
{
    return direct_;
}


//- Return direct addressing
const unallocLabelList& topoPatchMapper::directAddressing() const
{
    if (!direct())
    {
        FatalErrorIn
        (
            "const unallocLabelList& "
            "topoPatchMapper::directAddressing() const"
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
const labelListList& topoPatchMapper::addressing() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const labelListList& "
            "topoPatchMapper::addressing() const"
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
const scalarListList& topoPatchMapper::weights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const scalarListList& "
            "topoPatchMapper::weights() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!weightsPtr_)
    {
        calcAddressing();
    }

    return *weightsPtr_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void topoPatchMapper::operator=(const topoPatchMapper& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("void topoPatchMapper::operator=")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
