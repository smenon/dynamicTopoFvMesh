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
    dynamicTopoFvMeshMapper

Description
    Class holds all necessary information for mapping fields associated with
    dynamicTopoFvMesh and fvMesh.

Author
    Sandeep Menon

SourceFiles
    dynamicTopoFvMeshMapper.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicTopoFvMeshMapper_H
#define dynamicTopoFvMeshMapper_H

#include "faceMapper.H"
#include "inverseDistanceCellMapper.H"
#include "fvSurfaceMapper.H"
#include "fvBoundaryMeshMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class dynamicTopoFvMesh;    
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class dynamicTopoFvMeshMapper Declaration
\*---------------------------------------------------------------------------*/

class dynamicTopoFvMeshMapper
{
    // Private data

        //- Reference to dynamicTopoFvMesh
	const dynamicTopoFvMesh& mesh_;
        
        //- Face mapper
        faceMapper faceMap_;

        //- Inverse-distance based interpolative cell mapper
        inverseDistanceCellMapper cellMap_;

        //- Surface mapper
        fvSurfaceMapper surfaceMap_;

        //- Boundary mapper
        fvBoundaryMeshMapper boundaryMap_;        

    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicTopoFvMeshMapper(const dynamicTopoFvMeshMapper&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicTopoFvMeshMapper&);


public:

    // Constructors

        //- Construct from dynamicTopoFvMesh
        dynamicTopoFvMeshMapper
        (
            const dynamicTopoFvMesh& mesh,
            const mapPolyMesh& mpm
        )
	:
	    mesh_(mesh),                   
            faceMap_(mpm),
            cellMap_(mesh, mpm),
            surfaceMap_(dynamicCast<const fvMesh>(mesh), faceMap_),
            boundaryMap_(dynamicCast<const fvMesh>(mesh), faceMap_)            
	{}

    // Destructor

	~dynamicTopoFvMeshMapper();

    // Member Functions

	//- Return reference to the mesh
	const dynamicTopoFvMesh& mesh() const
	{
	    return mesh_;
	}
        
        //- Return volume mapper
        const FieldMapper& volMap() const
        {
            return cellMap_;
        }

        //- Return surface mapper
        const fvSurfaceMapper& surfaceMap() const
        {
            return surfaceMap_;
        }

        //- Return boundary mapper
        const fvBoundaryMeshMapper& boundaryMap() const
        {
            return boundaryMap_;
        }        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "dynamicTopoFvMeshMapper.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
