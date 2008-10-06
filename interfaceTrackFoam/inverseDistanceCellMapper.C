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
    inverseDistanceCellMapper

Description
    Implementation of the inverseDistanceCellMapper class

Author
    Sandeep Menon

\*----------------------------------------------------------------------------*/

#include "inverseDistanceCellMapper.H"
#include "demandDrivenData.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::inverseDistanceCellMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
    deleteDemandDrivenData(interpolationAddrPtr_);
    deleteDemandDrivenData(weightsPtr_);
}

void Foam::inverseDistanceCellMapper::calcAddressing() const
{
    if
    (
        directAddrPtr_
     || interpolationAddrPtr_
     || weightsPtr_
    )
    {
        FatalErrorIn("void inverseDistanceCellMapper::calcAddressing() const")
            << "Addressing already calculated."
            << abort(FatalError);
    }
    
    if (direct())
    {
        // Direct addressing, no weights

        directAddrPtr_ = new labelList(mpm_.cellMap());
    }
    else
    {
        // Interpolative addressing

        interpolationAddrPtr_ = new labelListList(mesh_.nCells());
        labelListList& addr = *interpolationAddrPtr_;

        weightsPtr_ = new scalarListList(mesh_.nCells());
        scalarListList& w = *weightsPtr_;    
    
        // Obtain cell-centre information from old/new meshes
        const vectorField& oldCentres = mesh_.oldCellCentres();
        const vectorField& newCentres = mesh_.cellCentres();
        
        const List<objectMap>& cfc = mpm_.cellsFromCellsMap();

        forAll (cfc, cfcI)
        {
            // Get addressing
            const labelList& mo = cfc[cfcI].masterObjects();

            label cellI = cfc[cfcI].index();

            if (addr[cellI].size() > 0)
            {
                FatalErrorIn("inverseDistanceCellMapper::calcAddressing()")
                    << "Master cell " << cellI
                    << " mapped from cell " << mo
                    << " already destination of mapping." << abort(FatalError);
            }

            // Map from masters, inverse-distance weights
            addr[cellI] = mo;
            scalar totalWeight = 0.0;
            w[cellI] = scalarList(mo.size(), 0.0);
            
            forAll (mo, oldCellI)
            {
                /*
                Info << cellI << ": " << newCentres[cellI] 
                     << ": " << mo[oldCellI]<< ": " 
                     << oldCentres[mo[oldCellI]] << endl;
                */
                w[cellI][oldCellI] = 
                   1.0/magSqr(newCentres[cellI] - oldCentres[mo[oldCellI]]);
                totalWeight += w[cellI][oldCellI];
            }
            
            // Normalize weights
            scalar normFactor = (1.0/totalWeight);
            forAll (mo, oldCellI)
            {            
                w[cellI][oldCellI] *= normFactor;
            }
        }
        
        // Do mapped cells. Note that this can already be set by cellsFromCells
        // so check if addressing size still zero.

        const labelList& cm = mpm_.cellMap();

        forAll (cm, cellI)
        {
            if (cm[cellI] > -1 && addr[cellI].size() == 0)
            {
                // Mapped from a single cell
                addr[cellI] = labelList(1, cm[cellI]);
                w[cellI] = scalarList(1, 1.0);
            }
        }  
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::inverseDistanceCellMapper::inverseDistanceCellMapper
(
    const dynamicTopoFvMesh& mesh,         
    const mapPolyMesh& mpm
)
:
    mesh_(mesh),
    mpm_(mpm),
    direct_(false),
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::inverseDistanceCellMapper::~inverseDistanceCellMapper()
{
    clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::inverseDistanceCellMapper::size() const
{
    return mpm_.cellMap().size();
}

Foam::label Foam::inverseDistanceCellMapper::sizeBeforeMapping() const
{
    return mpm_.nOldCells();
}

const Foam::unallocLabelList& 
Foam::inverseDistanceCellMapper::directAddressing() const
{
    if (!direct())
    {
        FatalErrorIn
        (
            "const unallocLabelList& inverseDistanceCellMapper::directAddressing() const"
        )   << "Requested direct addressing for an interpolative mapper."
            << abort(FatalError);
    }

    if (!directAddrPtr_)
    {
        calcAddressing();
    }

    return *directAddrPtr_;
}

const Foam::labelListList& Foam::inverseDistanceCellMapper::addressing() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const labelListList& inverseDistanceCellMapper::addressing() const"
        )   << "Requested interpolative addressing for a direct mapper."
            << abort(FatalError);
    }

    if (!interpolationAddrPtr_)
    {
        calcAddressing();
    }

    return *interpolationAddrPtr_;
}

const Foam::scalarListList& Foam::inverseDistanceCellMapper::weights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const scalarListList& cellMapper::weights() const"
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

void Foam::inverseDistanceCellMapper::operator=(const inverseDistanceCellMapper& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("inverseDistanceCellMapper::operator=(const inverseDistanceCellMapper&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
