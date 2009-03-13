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
    This object provides mapping and fill-in information for cell data
    between the two meshes after the topological change.  It is
    constructed from mapPolyMesh, and uses inverse-distance weighting.    

Author
    Sandeep Menon

SourceFiles
    inverseDistanceCellMapper.C

\*---------------------------------------------------------------------------*/

#ifndef inverseDistanceCellMapper_H
#define inverseDistanceCellMapper_H

#include "FieldMapper.H"
#include "dynamicTopoFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class inverseDistanceCellMapper Declaration
\*---------------------------------------------------------------------------*/

class inverseDistanceCellMapper
:
    public FieldMapper
{
    // Private data

        //- Reference to polyMesh
        const dynamicTopoFvMesh& mesh_;

        //- Reference to mapPolyMesh
        const mapPolyMesh& mpm_;
        
        //- Is the mapping direct
        bool direct_;        

    // Demand-driven private data

        //- Direct addressing (only one for of addressing is used)
        mutable labelList* directAddrPtr_;

        //- Interpolated addressing (only one for of addressing is used)
        mutable labelListList* interpolationAddrPtr_;

        //- Interpolation weights
        mutable scalarListList* weightsPtr_;        
        

    // Private Member Functions

        //- Disallow default bitwise copy construct
        inverseDistanceCellMapper(const inverseDistanceCellMapper&);

        //- Disallow default bitwise assignment
        void operator=(const inverseDistanceCellMapper&);

        //- Calculate addressing for interpolative mapping
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();        

public:

    // Constructors    
    
        //- Construct from mapPolyMesh
        inverseDistanceCellMapper(const mapPolyMesh& mpm);

    // Destructor

        virtual ~inverseDistanceCellMapper();        

    // Member Functions        
        
        //- Return size
        virtual label size() const;

        //- Return size before mapping
        virtual label sizeBeforeMapping() const;

        //- Is the mapping direct
        virtual bool direct() const
        {
            return direct_;
        }

        //- Return direct addressing
        virtual const unallocLabelList& directAddressing() const;

        //- Return interpolated addressing
        virtual const labelListList& addressing() const;

        //- Return interpolaion weights
        virtual const scalarListList& weights() const;        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#    include "inverseDistanceCellMapper.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
