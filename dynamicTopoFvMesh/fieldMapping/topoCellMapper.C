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
    topoCellMapper

Description
    Implementation of the topoCellMapper class

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*----------------------------------------------------------------------------*/

#include "topoCellMapper.H"
#include "demandDrivenData.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::topoCellMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
    deleteDemandDrivenData(interpolationAddrPtr_);
    deleteDemandDrivenData(weightsPtr_);
    deleteDemandDrivenData(oldCellCentresPtr_);
}

void Foam::topoCellMapper::calcAddressing() const
{
    if
    (
        directAddrPtr_
     || interpolationAddrPtr_
     || weightsPtr_
    )
    {
        FatalErrorIn("void topoCellMapper::calcAddressing() const")
            << "Addressing already calculated."
            << abort(FatalError);
    }
    
    if (!oldCellCentresPtr_)
    {
        FatalErrorIn("void topoCellMapper::calcAddressing() const")
            << "Cell centres has not been set."
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
        const vectorField& oldCentres = oldCellCentresPtr_->internalField();
        const vectorField& newCentres = mesh_.cellCentres();
        
        const List<objectMap>& cfc = mpm_.cellsFromCellsMap();

        forAll (cfc, cfcI)
        {
            // Get addressing
            const labelList& mo = cfc[cfcI].masterObjects();

            label cellI = cfc[cfcI].index();

            if (addr[cellI].size() > 0)
            {
                FatalErrorIn("topoCellMapper::calcAddressing()")
                    << "Master cell " << cellI
                    << " mapped from cell " << mo
                    << " already destination of mapping."
                    << abort(FatalError);
            }

            // Map from masters, inverse-distance weights
            addr[cellI] = mo;
            scalar totalWeight = 0.0;
            w[cellI] = scalarList(mo.size(), 0.0);
            
            forAll (mo, oldCellI)
            {
                w[cellI][oldCellI] =
                (
                    1.0/stabilise
                    (
                        magSqr
                        (
                            newCentres[cellI]
                          - oldCentres[mo[oldCellI]]
                        ),
                        VSMALL
                    )
                );

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
Foam::topoCellMapper::topoCellMapper
(         
    const mapPolyMesh& mpm
)
:
    mesh_(mpm.mesh()),
    mpm_(mpm),
    direct_(false),
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL),
    oldCellCentresPtr_(NULL)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::topoCellMapper::~topoCellMapper()
{
    clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Set old cell-centre information
void Foam::topoCellMapper::setOldCellCentres
(
    const volVectorField& oldCentres
) const
{
    if (oldCellCentresPtr_)
    {
        FatalErrorIn
        (
            "void topoCellMapper::setOldCellCentres()"
        ) << nl << " Pointer has already been set. "
          << abort(FatalError);
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


Foam::label Foam::topoCellMapper::size() const
{
    return mpm_.cellMap().size();
}


Foam::label Foam::topoCellMapper::sizeBeforeMapping() const
{
    return mpm_.nOldCells();
}


const Foam::unallocLabelList& 
Foam::topoCellMapper::directAddressing() const
{
    if (!direct())
    {
        FatalErrorIn
        (
            "const unallocLabelList& topoCellMapper::directAddressing() const"
        )   << "Requested direct addressing for an interpolative mapper."
            << abort(FatalError);
    }

    if (!directAddrPtr_)
    {
        calcAddressing();
    }

    return *directAddrPtr_;
}


const Foam::labelListList& Foam::topoCellMapper::addressing() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const labelListList& topoCellMapper::addressing() const"
        )   << "Requested interpolative addressing for a direct mapper."
            << abort(FatalError);
    }

    if (!interpolationAddrPtr_)
    {
        calcAddressing();
    }

    return *interpolationAddrPtr_;
}


const Foam::scalarListList& Foam::topoCellMapper::weights() const
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

void Foam::topoCellMapper::operator=
(
    const topoCellMapper& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("topoCellMapper::operator=")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
