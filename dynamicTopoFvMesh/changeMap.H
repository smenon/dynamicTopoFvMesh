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
#include "FixedList.H"

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

    // Entities that were added during the operation.
    //  - Add a master mapping for slaves.
    //  - Masters get a map of -1.
    List<FixedList<label,2> > addedPoints_;
    List<FixedList<label,2> > addedEdges_;
    List<FixedList<label,2> > addedFaces_;
    List<FixedList<label,2> > addedCells_;

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
        void addPoint
        (
            const label pIndex,
            const label master = -1
        );

        void addEdge
        (
            const label eIndex,
            const label master = -1
        );

        void addFace
        (
            const label fIndex,
            const label master = -1
        );

        void addCell
        (
            const label cIndex,
            const label master = -1
        );

        // Return the list of added entities
        const List<FixedList<label,2> >& addedPointList() const;
        const List<FixedList<label,2> >& addedEdgeList() const;
        const List<FixedList<label,2> >& addedFaceList() const;
        const List<FixedList<label,2> >& addedCellList() const;

    //- Operators

        void operator=(const changeMap& rhs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
