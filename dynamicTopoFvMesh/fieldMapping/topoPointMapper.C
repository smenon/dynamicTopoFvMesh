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
    topoPointMapper

Description
    Implementation of the topoPointMapper class

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "topoMapper.H"
#include "mapPolyMesh.H"
#include "topoPointMapper.H"
#include "demandDrivenData.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Calculate addressing for interpolative mapping
void topoPointMapper::calcAddressing() const
{
    if (directAddrPtr_ || interpolationAddrPtr_)
    {
        FatalErrorIn("void topoPointMapper::calcAddressing() const")
            << "Addressing already calculated."
            << abort(FatalError);
    }

    // Allocate for inserted point labels
    label nInsertedPoints = 0;

    insertedPointLabelsPtr_ = new labelList(mpm_.mesh().nPoints(), -1);
    labelList& insertedPoints = *insertedPointLabelsPtr_;

    if (direct())
    {
        // Direct addressing, no weights
        directAddrPtr_ = new labelList(mpm_.pointMap());
    }
    else
    {
        // Interpolative addressing
        interpolationAddrPtr_ = new labelListList(mpm_.mesh().nPoints());
        labelListList& addr = *interpolationAddrPtr_;

        const List<objectMap>& pfp = mpm_.pointsFromPointsMap();

        forAll(pfp, pfpI)
        {
            // Get addressing
            const labelList& mo = pfp[pfpI].masterObjects();

            label pointI = pfp[pfpI].index();

            if (addr[pointI].size() > 0)
            {
                FatalErrorIn
                (
                    "void topoPointMapper::calcAddressing() const"
                )
                    << "Master point " << pointI
                    << " mapped from points " << mo
                    << " is already mapping from: " << addr[pointI]
                    << abort(FatalError);
            }

            // Set master objects
            addr[pointI] = mo;
        }

        // Do mapped points.
        // Note that this can already be set by pointsFromPoints,
        // so check if addressing size still zero.
        const labelList& pm = mpm_.pointMap();

        forAll(pm, pointI)
        {
            // Mapped from a single point
            if (pm[pointI] > -1 && addr[pointI].empty())
            {
                addr[pointI] = labelList(1, pm[pointI]);
            }

            // Check for inserted points without any addressing
            if (pm[pointI] < 0 && addr[pointI].empty())
            {
                insertedPoints[nInsertedPoints++] = pointI;
            }
        }
    }

    // Shorten inserted points to actual size
    insertedPoints.setSize(nInsertedPoints);

    if (nInsertedPoints)
    {
        FatalErrorIn("void topoPointMapper::calcAddressing() const")
            << " Found " << nInsertedPoints << " points which are"
            << " not mapped from any parent points." << nl
            << " List: " << nl
            << insertedPoints
            << abort(FatalError);
    }
}


//- Calculate weights for interpolative mapping
void topoPointMapper::calcWeights() const
{
    if (weightsPtr_)
    {
        FatalErrorIn("void topoPointMapper::calcWeights() const")
            << "Weights already calculated."
            << abort(FatalError);
    }

    // Fetch interpolative addressing
    const labelListList& addr = addressing();

    // Allocate memory
    weightsPtr_ = new scalarListList(size());
    scalarListList& weights = *weightsPtr_;

    // Fetch maps
    const List<objectMap>& pfp = mpm_.pointsFromPointsMap();
    const List<scalarField>& mapPointWeights = tMapper_.pointWeights();

    // Fill in maps first
    forAll(pfp, indexI)
    {
        weights[pfp[indexI].index()] = mapPointWeights[indexI];
    }

    // Now do mapped points
    forAll(addr, pointI)
    {
        const labelList& mo = addr[pointI];

        // Check if this is indeed a mapped point
        if (mo.size() == 1 && weights[pointI].empty())
        {
            weights[pointI] = scalarList(1, 1.0);
        }
    }
}


//- Clear out local storage
void topoPointMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
    deleteDemandDrivenData(interpolationAddrPtr_);
    deleteDemandDrivenData(weightsPtr_);
    deleteDemandDrivenData(insertedPointLabelsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
topoPointMapper::topoPointMapper
(
    const pointMesh& mesh,
    const mapPolyMesh& mpm,
    const topoMapper& mapper
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

topoPointMapper::~topoPointMapper()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return size
label topoPointMapper::size() const
{
    return mpm_.pointMap().size();
}


//- Return size before mapping
label topoPointMapper::sizeBeforeMapping() const
{
    return sizeBeforeMapping_;
}


//- Is the mapping direct
bool topoPointMapper::direct() const
{
    return direct_;
}


//- Has unmapped elements
bool topoPointMapper::hasUnmapped() const
{
    return false;
}


//- Return direct addressing
const unallocLabelList& topoPointMapper::directAddressing() const
{
    if (!direct())
    {
        FatalErrorIn
        (
            "const unallocLabelList& "
            "topoPointMapper::directAddressing() const"
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
const labelListList& topoPointMapper::addressing() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const labelListList& "
            "topoPointMapper::addressing() const"
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
const scalarListList& topoPointMapper::weights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const scalarListList& "
            "topoPointMapper::weights() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!weightsPtr_)
    {
        calcWeights();
    }

    return *weightsPtr_;
}


//- Are there any inserted objects
bool topoPointMapper::insertedObjects() const
{
    return insertedObjectLabels().size();
}


//- Return list of inserted objects
const labelList& topoPointMapper::insertedObjectLabels() const
{
    if (!insertedPointLabelsPtr_)
    {
        calcAddressing();
    }

    return *insertedPointLabelsPtr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void topoPointMapper::operator=(const topoPointMapper& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("void topoPointMapper::operator=")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}


} // End namespace Foam

// ************************************************************************* //
