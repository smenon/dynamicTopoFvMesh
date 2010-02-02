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

\*---------------------------------------------------------------------------*/

#include "changeMap.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Constructor
changeMap::changeMap()
:
    type_(-1),
    firstEdge_(-1),
    secondEdge_(-1),
    apexPoint_(-1),
    opposingFace_(-1)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Type
label& changeMap::type()
{
    return type_;
}

// For sliver-type cells, opposite edges
// are identified for removal.
label& changeMap::firstEdge()
{
    return firstEdge_;
}

label& changeMap::secondEdge()
{
    return secondEdge_;
}

// For cap-type cells, the face requiring splitting
// is identified for removal.
label& changeMap::apexPoint()
{
    return apexPoint_;
}

label& changeMap::opposingFace()
{
    return opposingFace_;
}

// Added entities
void changeMap::addPoint
(
    const label pIndex,
    const label master
)
{
    addedPoints_.insert(pIndex, master);
}

void changeMap::addEdge
(
    const label eIndex,
    const label master
)
{
    addedEdges_.insert(eIndex, master);
}

void changeMap::addFace
(
    const label fIndex,
    const label master
)
{
    addedFaces_.insert(fIndex, master);
}

void changeMap::addCell
(
    const label cIndex,
    const label master
)
{
    addedCells_.insert(cIndex, master);
}

// Return an added point
const Map<label>& changeMap::addedPointList() const
{
    return addedPoints_;
}

// Return the list of added entities
const Map<label>& changeMap::addedEdgeList() const
{
    return addedEdges_;
}

const Map<label>& changeMap::addedFaceList() const
{
    return addedFaces_;
}

const Map<label>& changeMap::addedCellList() const
{
    return addedCells_;
}

void changeMap::operator=(const changeMap& rhs)
{
    type_ = rhs.type_;

    firstEdge_    = rhs.firstEdge_;
    secondEdge_   = rhs.secondEdge_;
    apexPoint_    = rhs.apexPoint_;
    opposingFace_ = rhs.opposingFace_;

    addedPoints_  = rhs.addedPoints_;
    addedEdges_   = rhs.addedEdges_;
    addedFaces_   = rhs.addedFaces_;
    addedCells_   = rhs.addedCells_;
}

} // End namespace Foam

// ************************************************************************* //
