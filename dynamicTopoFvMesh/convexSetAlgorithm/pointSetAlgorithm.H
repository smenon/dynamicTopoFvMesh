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
    pointSetAlgorithm

Description
    Class for convexSetAlgorithms pertaining to points

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    pointSetAlgorithm.C

\*---------------------------------------------------------------------------*/

#ifndef pointSetAlgorithm_H
#define pointSetAlgorithm_H

#include "convexSetAlgorithm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class pointSetAlgorithm Declaration
\*---------------------------------------------------------------------------*/

class pointSetAlgorithm
:
    public convexSetAlgorithm
{

    // Private member data

        //- Tree for candidate search
        indexedOctree<mappingTreeData> searchTree_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointSetAlgorithm(const pointSetAlgorithm&);

        //- Disallow default bitwise assignment
        void operator=(const pointSetAlgorithm&);

public:

    //- Constructor

        pointSetAlgorithm
        (
            const polyMesh& mesh,
            const pointField& newPoints,
            const UList<edge>& newEdges,
            const UList<face>& newFaces,
            const UList<cell>& newCells,
            const UList<label>& newOwner,
            const UList<label>& newNeighbour
        );

    //- Destructor

        virtual ~pointSetAlgorithm()
        {}

    //- Virtual member functions

        // Dimensions of the algorithm
        virtual label dimension() const
        {
            return 0;
        }

        // Compute normFactor
        virtual void computeNormFactor(const label index) const;

        // Find the nearest mapping candidates
        virtual void findMappingCandidates(labelList& mapCandidates) const;

        // Write out mapping candidates
        virtual void writeMappingCandidates() const;

        // Check whether the bounding box contains the entity
        virtual bool contains(const label index) const;

        // Compute intersections
        virtual bool computeIntersection
        (
            const label newIndex,
            const label oldIndex,
            const label offset,
            bool output
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
