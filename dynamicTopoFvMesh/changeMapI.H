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
    changeMap

Description
    Accumulate topology change statistics

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Entity index
inline label& changeMap::index()
{
    return index_;
}


inline label changeMap::index() const
{
    return index_;
}


// Coupled patch index
inline label& changeMap::patchIndex()
{
    return pIndex_;
}


inline label changeMap::patchIndex() const
{
    return pIndex_;
}


// Type
inline label& changeMap::type()
{
    return type_;
}


inline label changeMap::type() const
{
    return type_;
}


// For sliver-type cells, opposite edges
// are identified for removal.
inline label& changeMap::firstEdge()
{
    return firstEdge_;
}


inline label& changeMap::secondEdge()
{
    return secondEdge_;
}


// For cap-type cells, the face requiring splitting
// is identified for removal.
inline label& changeMap::apexPoint()
{
    return apexPoint_;
}


inline label& changeMap::opposingFace()
{
    return opposingFace_;
}


// Added entities
inline void changeMap::addPoint
(
    const label pIndex,
    const labelList& master
)
{
    label curSize = addedPoints_.size();

    addedPoints_.setSize(curSize + 1);

    addedPoints_[curSize].first() = pIndex;
    addedPoints_[curSize].second() = master;
}


inline void changeMap::addEdge
(
    const label eIndex,
    const labelList& master
)
{
    label curSize = addedEdges_.size();

    addedEdges_.setSize(curSize + 1);

    addedEdges_[curSize].first() = eIndex;
    addedEdges_[curSize].second() = master;
}


inline void changeMap::addFace
(
    const label fIndex,
    const labelList& master
)
{
    label curSize = addedFaces_.size();

    addedFaces_.setSize(curSize + 1);

    addedFaces_[curSize].first() = fIndex;
    addedFaces_[curSize].second() = master;
}


inline void changeMap::addCell
(
    const label cIndex,
    const labelList& master
)
{
    label curSize = addedCells_.size();

    addedCells_.setSize(curSize + 1);

    addedCells_[curSize].first() = cIndex;
    addedCells_[curSize].second() = master;
}


// Return an added point
inline const List<Tuple2<label, labelList> >&
changeMap::addedPointList() const
{
    return addedPoints_;
}


// Return the list of added entities
inline const List<Tuple2<label, labelList> >&
changeMap::addedEdgeList() const
{
    return addedEdges_;
}


inline const List<Tuple2<label, labelList> >&
changeMap::addedFaceList() const
{
    return addedFaces_;
}


inline const List<Tuple2<label, labelList> >&
changeMap::addedCellList() const
{
    return addedCells_;
}


inline void changeMap::operator=(const changeMap& rhs)
{
    index_ = rhs.index_;
    pIndex_ = rhs.pIndex_;

    type_ = rhs.type_;

    // Copy sliver-specific information
    firstEdge_ = rhs.firstEdge_;
    secondEdge_ = rhs.secondEdge_;
    apexPoint_ = rhs.apexPoint_;
    opposingFace_ = rhs.opposingFace_;

    // Copy maps
    addedPoints_.setSize(rhs.addedPoints_.size());

    forAll(addedPoints_, indexI)
    {
        addedPoints_[indexI].first() = rhs.addedPoints_[indexI].first();
        addedPoints_[indexI].second() = rhs.addedPoints_[indexI].second();
    }

    addedEdges_.setSize(rhs.addedEdges_.size());

    forAll(addedEdges_, indexI)
    {
        addedEdges_[indexI].first() = rhs.addedEdges_[indexI].first();
        addedEdges_[indexI].second() = rhs.addedEdges_[indexI].second();
    }

    addedFaces_.setSize(rhs.addedFaces_.size());

    forAll(addedFaces_, indexI)
    {
        addedFaces_[indexI].first() = rhs.addedFaces_[indexI].first();
        addedFaces_[indexI].second() = rhs.addedFaces_[indexI].second();
    }

    addedCells_.setSize(rhs.addedCells_.size());

    forAll(addedCells_, indexI)
    {
        addedCells_[indexI].first() = rhs.addedCells_[indexI].first();
        addedCells_[indexI].second() = rhs.addedCells_[indexI].second();
    }
}


inline Ostream& operator<<(Ostream& os, changeMap& cm)
{
    // Write indices
    os << "Index: " << cm.index_ << endl;
    os << "patchIndex: " << cm.pIndex_ << endl;

    // Write out type
    os << " Type: " << cm.type_ << endl;

    // Write out sliver-specific information
    os << " firstEdge: " << cm.firstEdge_ << endl;
    os << " secondEdge: " << cm.secondEdge_ << endl;
    os << " apexPoint: " << cm.apexPoint_ << endl;
    os << " opposingFace: " << cm.opposingFace_ << endl;

    // Write out added entity information
    os << " addedPoints: " << cm.addedPoints_ << endl;
    os << " addedEdges: " << cm.addedEdges_ << endl;
    os << " addedFaces: " << cm.addedFaces_ << endl;
    os << " addedCells: " << cm.addedCells_ << endl;

    // Check state of IOstream
    os.check("Ostream& operator<<(Ostream&, changeMap&)");

    return os;
}


} // End namespace Foam

// ************************************************************************* //
