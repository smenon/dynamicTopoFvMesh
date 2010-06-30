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

#include "topoMapper.H"
#include "mapPolyMesh.H"
#include "topoCellMapper.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Clear out local storage
void Foam::topoCellMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
    deleteDemandDrivenData(interpolationAddrPtr_);
    deleteDemandDrivenData(weightsPtr_);
    deleteDemandDrivenData(volumesPtr_);
    deleteDemandDrivenData(centresPtr_);
}


//- Calculate addressing for interpolative mapping
void Foam::topoCellMapper::calcAddressing() const
{
    if
    (
        directAddrPtr_
     || interpolationAddrPtr_
    )
    {
        FatalErrorIn("void topoCellMapper::calcAddressing() const")
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

        const List<objectMap>& cfc = mpm_.cellsFromCellsMap();

        forAll (cfc, cfcI)
        {
            // Get addressing
            const labelList& mo = cfc[cfcI].masterObjects();

            label cellI = cfc[cfcI].index();

            if (addr[cellI].size() > 0)
            {
                FatalErrorIn("void topoCellMapper::calcAddressing() const")
                    << "Master cell " << cellI
                    << " mapped from cell " << mo
                    << " already destination of mapping."
                    << abort(FatalError);
            }

            // Set master objects
            addr[cellI] = mo;
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
            }

            if (cm[cellI] < 0)
            {
                FatalErrorIn("void topoCellMapper::calcAddressing() const")
                    << "Cell " << cellI
                    << " is not mapped from any parent cell."
                    << abort(FatalError);
            }
        }
    }
}


//- Calculate inverse-distance weights for interpolative mapping
void Foam::topoCellMapper::calcInverseDistanceWeights() const
{
    if (weightsPtr_)
    {
        FatalErrorIn("void topoCellMapper::calcInverseDistanceWeights() const")
            << "Weights already calculated."
            << abort(FatalError);
    }

    // Fetch interpolative addressing
    const labelListList& addr = addressing();

    // Allocate memory
    weightsPtr_ = new scalarListList(mesh_.nCells());
    scalarListList& w = *weightsPtr_;

    // Obtain cell-centre information from old/new meshes
    const vectorField& oldCentres = tMapper_.oldCentres().internalField();
    const vectorField& newCentres = mesh_.cellCentres();

    forAll(addr, cellI)
    {
        const labelList& mo = addr[cellI];

        // Do mapped cells
        if (mo.size() == 1)
        {
            w[cellI] = scalarList(1, 1.0);
        }
        else
        {
            // Map from masters, inverse-distance weights
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
    }
}


//- Calculate intersection weights for conservative mapping
void Foam::topoCellMapper::calcIntersectionWeightsAndCentres() const
{
    if (volumesPtr_ || centresPtr_)
    {
        FatalErrorIn
        (
            "void topoCellMapper::"
            "calcIntersectionWeightsAndCentres() const"
        )
            << "Weights already calculated."
            << abort(FatalError);
    }

    // Fetch interpolative addressing
    const labelListList& addr = addressing();

    // Allocate memory
    volumesPtr_ = new List<scalarField>(mesh_.nCells());
    List<scalarField>& v = *volumesPtr_;

    centresPtr_ = new List<vectorField>(mesh_.nCells());
    List<vectorField>& x = *centresPtr_;

    // Obtain cell-centre / volume information from the mesh.
    //  - Point positions must be set to old locations
    //    for mapping to be valid.
    const vectorField& cellCentres = mesh_.cellCentres();
    const scalarField& cellVolumes = mesh_.cellVolumes();

    // Fetch maps
    const Map<vectorField>& mapCellCentres = tMapper_.cellCentres();
    const Map<scalarField>& mapCellWeights = tMapper_.cellWeights();

    forAll(addr, cellI)
    {
        const labelList& mo = addr[cellI];

        // Do mapped cells
        if (mo.size() == 1)
        {
            x[cellI] = vectorField(1, cellCentres[cellI]);
            v[cellI] = scalarField(1, cellVolumes[cellI]);
        }
        else
        {
            // Map from masters, intersection weights
            x[cellI] = mapCellCentres[cellI];
            v[cellI] = mapCellWeights[cellI];

            if (mag(sum(v[cellI]) - cellVolumes[cellI]) > 1e-16)
            {
                FatalErrorIn
                (
                    "void topoCellMapper::"
                    "calcIntersectionWeightsAndCentres() const"
                )
                    << "Weights are inconsistent." << nl
                    << " Cell volume: " << cellVolumes[cellI] << nl
                    << " sum(weights): " << sum(v[cellI]) << nl
                    << abort(FatalError);
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::topoCellMapper::topoCellMapper
(
    const mapPolyMesh& mpm,
    const topoMapper& mapper
)
:
    mesh_(mpm.mesh()),
    mpm_(mpm),
    tMapper_(mapper),
    direct_(false),
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL),
    volumesPtr_(NULL),
    centresPtr_(NULL)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::topoCellMapper::~topoCellMapper()
{
    clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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
            "const unallocLabelList&"
            "topoCellMapper::directAddressing() const"
        )   << "Requested direct addressing for an interpolative mapper."
            << abort(FatalError);
    }

    if (!directAddrPtr_)
    {
        calcAddressing();
    }

    return *directAddrPtr_;
}


const Foam::labelListList&
Foam::topoCellMapper::addressing() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const labelListList&"
            "topoCellMapper::addressing() const"
        )   << "Requested interpolative addressing for a direct mapper."
            << abort(FatalError);
    }

    if (!interpolationAddrPtr_)
    {
        calcAddressing();
    }

    return *interpolationAddrPtr_;
}


const Foam::scalarListList&
Foam::topoCellMapper::weights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const scalarListList&"
            "topoCellMapper::weights() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!weightsPtr_)
    {
        calcInverseDistanceWeights();
    }

    return *weightsPtr_;
}


const Foam::List<Foam::scalarField>&
Foam::topoCellMapper::intersectionWeights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const List<scalarField>&"
            "topoCellMapper::intersectionWeights() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!volumesPtr_)
    {
        calcIntersectionWeightsAndCentres();
    }

    return *volumesPtr_;
}


const Foam::List<Foam::vectorField>&
Foam::topoCellMapper::intersectionCentres() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const List<vectorField>&"
            "topoCellMapper::intersectionCentres() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!centresPtr_)
    {
        calcIntersectionWeightsAndCentres();
    }

    return *centresPtr_;
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
