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
    mesquiteMapper

Description
    Class holds all necessary information for mapping fields associated with
    mesquiteMotionSolver.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef mesquiteMapper_H
#define mesquiteMapper_H

#include "mesquiteInternalMapper.H"
#include "mesquiteBoundaryMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class pointMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                      Class mesquiteMapper Declaration
\*---------------------------------------------------------------------------*/

class mesquiteMapper
{
    // Private data

        //- Reference to mesh
        const pointMesh& mesh_;

        //- Internal field mapper
        mesquiteInternalMapper internalMap_;

        //- Boundary field mapper
        mesquiteBoundaryMapper boundaryMap_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        mesquiteMapper(const mesquiteMapper&);

        //- Disallow default bitwise assignment
        void operator=(const mesquiteMapper&);

public:

    // Constructors

        //- Construct from pointMesh and mapPolyMesh
        mesquiteMapper(const pointMesh& mesh, const mapPolyMesh& mpm)
        :
            mesh_(mesh),
            internalMap_(mesh, mpm, *this),
            boundaryMap_(mesh, mpm, *this)
        {}

    // Member Functions

        //- Return reference to the pointMesh
        const pointMesh& mesh() const
        {
            return mesh_;
        }

        //- Return reference to objectRegistry storing fields
        const objectRegistry& thisDb() const
        {
            return mesh_();
        }

        //- Return the internal field mapper
        const mesquiteInternalMapper& internalMap() const
        {
            return internalMap_;
        }

        //- Return the boundary field mapper
        const mesquiteBoundaryMapper& boundaryMap() const
        {
            return boundaryMap_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
