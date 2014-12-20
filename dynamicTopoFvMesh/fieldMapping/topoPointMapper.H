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
    This object provides mapping and fill-in information for internal
    point data between the two meshes after topological changes.
    It is constructed from mapPolyMesh.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    topoPointMapper.C

\*---------------------------------------------------------------------------*/

#ifndef topoPointMapper_H
#define topoPointMapper_H

#include "morphFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class pointMesh;
class topoMapper;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                        Class topoPointMapper Declaration
\*---------------------------------------------------------------------------*/

class topoPointMapper
:
    public morphFieldMapper
{
    // Private data

        //- Reference to mapPolyMesh
        const mapPolyMesh& mpm_;

        //- Reference to the topoMapper
        const topoMapper& tMapper_;

        //- Is the mapping direct
        bool direct_;

        //- Size before mapping
        mutable label sizeBeforeMapping_;

    // Demand-driven private data

        //- Direct addressing
        mutable labelList* directAddrPtr_;

        //- Interpolated addressing
        mutable labelListList* interpolationAddrPtr_;

        //- Interpolation weights
        mutable scalarListList* weightsPtr_;

        //- Inserted points
        mutable labelList* insertedPointLabelsPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        topoPointMapper(const topoPointMapper&);

        //- Disallow default bitwise assignment
        void operator=(const topoPointMapper&);

        //- Calculate addressing for interpolative mapping
        void calcAddressing() const;

        //- Calculate weights for interpolative mapping
        void calcWeights() const;

        //- Clear out local storage
        void clearOut();

public:

    // Constructors

        //- Construct from components
        topoPointMapper
        (
            const pointMesh& mesh,
            const mapPolyMesh& mpm,
            const topoMapper& mapper
        );

    // Destructor

        virtual ~topoPointMapper();

    // Member Functions

        //- Return size
        virtual label size() const;

        //- Return size before mapping
        virtual label sizeBeforeMapping() const;

        //- Is the mapping direct
        virtual bool direct() const;

        //- Has unmapped elements
        virtual bool hasUnmapped() const;

        //- Return direct addressing
        virtual const unallocLabelList& directAddressing() const;

        //- Return interpolation addressing
        virtual const labelListList& addressing() const;

        //- Return weights
        virtual const scalarListList& weights() const;

        //- Are there any inserted objects
        virtual bool insertedObjects() const;

        //- Return list of inserted objects
        virtual const labelList& insertedObjectLabels() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
