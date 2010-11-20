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

SourceFiles
    changeMapI.H

\*---------------------------------------------------------------------------*/

#ifndef changeMap_H
#define changeMap_H

#include "objectMap.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class changeMap;

// * * * * * * * * Forward declaration of friend fuctions  * * * * * * * * * //

Ostream& operator<<(Ostream&, changeMap&);


/*---------------------------------------------------------------------------*\
                           Class changeMap Declaration
\*---------------------------------------------------------------------------*/

class changeMap
{
    // Entity indices
    label index_;

    // Coupled patch index
    label pIndex_;

    // Sliver type index.
    //  Type 1: Sliver
    //  Type 2: Cap
    //  Type 3: Spade
    //  Type 4: Wedge
    label type_;

    // Data specific to sliver-type cells
    label firstEdge_;
    label secondEdge_;

    // Data specific to cap-type cells
    label apexPoint_;
    label opposingFace_;

    // Entities that were added during the operation.
    List<objectMap> addedPoints_;
    List<objectMap> addedEdges_;
    List<objectMap> addedFaces_;
    List<objectMap> addedCells_;

public:

    // Constructor
    changeMap()
    :
        index_(-1),
        pIndex_(-1),
        type_(-1),
        firstEdge_(-1),
        secondEdge_(-1),
        apexPoint_(-1),
        opposingFace_(-1),
        addedPoints_(0),
        addedEdges_(0),
        addedFaces_(0),
        addedCells_(0)
    {}

    //- Access

        // Entity index
        inline label& index();
        inline label index() const;

        // Coupled patch index
        inline label& patchIndex();
        inline label patchIndex() const;

        // Type
        inline label& type();
        inline label type() const;

        // For sliver-type cells, opposite edges
        // are identified for removal.
        inline label& firstEdge();
        inline label& secondEdge();

        // For cap-type cells, the face requiring splitting
        // is identified for removal.
        inline label& apexPoint();
        inline label& opposingFace();

        // Added entities
        inline void addPoint
        (
            const label pIndex,
            const labelList& master = labelList()
        );

        inline void addEdge
        (
            const label eIndex,
            const labelList& master = labelList()
        );

        inline void addFace
        (
            const label fIndex,
            const labelList& master = labelList()
        );

        inline void addCell
        (
            const label cIndex,
            const labelList& master = labelList()
        );

        // Return the list of added entities
        inline const List<objectMap>& addedPointList() const;
        inline const List<objectMap>& addedEdgeList() const;
        inline const List<objectMap>& addedFaceList() const;
        inline const List<objectMap>& addedCellList() const;

    //- Operators

        inline void operator=(const changeMap& rhs);

    //- IOstream Operators

        inline friend Ostream& operator<<(Ostream&, changeMap&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "changeMapI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
