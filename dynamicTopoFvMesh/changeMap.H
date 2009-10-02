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
    changeMap.C

\*---------------------------------------------------------------------------*/

#ifndef changeMap_H
#define changeMap_H

#include "labelList.H"
#include "labelHashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class changeMap Declaration
\*---------------------------------------------------------------------------*/

class changeMap
{
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

    // For bisections.
    label addedPoint_;

    // Entities that were added during the operation.
    labelHashSet addedEdges_;
    labelHashSet addedFaces_;
    labelHashSet addedCells_;

public:

    // Constructor
    changeMap();

    //- Access

        // Type
        label& type();

        // For sliver-type cells, opposite edges
        // are identified for removal.
        label& firstEdge();
        label& secondEdge();

        // For cap-type cells, the face requiring splitting
        // is identified for removal.
        label& apexPoint();
        label& opposingFace();

        // Added entities
        void addPoint(const label pIndex);

        void addEdge(const label eIndex);

        void addFace(const label fIndex);

        void addCell(const label cIndex);

        // Return an added point
        label addedPoint();

        // Return the list of added entities
        const labelList addedEdgeList();

        const labelList addedFaceList();

        const labelList addedCellList();

    //- Operators

        void operator=(const changeMap& rhs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
