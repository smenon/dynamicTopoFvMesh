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
    topoPointBoundaryMapper

Description
    This object provides mapping and fill-in information for boundary
    point data between the two meshes after topological changes.
    It is constructed from mapPolyMesh.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef topoPointBoundaryMapper_H
#define topoPointBoundaryMapper_H

#include "PtrList.H"
#include "pointMesh.H"
#include "topoPointPatchMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class topoPointBoundaryMapper Declaration
\*---------------------------------------------------------------------------*/

class topoPointBoundaryMapper
:
    public PtrList<topoPointPatchMapper>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        topoPointBoundaryMapper(const topoPointBoundaryMapper&);

        //- Disallow default bitwise assignment
        void operator=(const topoPointBoundaryMapper&);

public:

    // Constructors

        //- Construct from components
        topoPointBoundaryMapper
        (
            const pointMesh& mesh,
            const mapPolyMesh& mpm,
            const topoMapper& mapper
        )
        :
            PtrList<topoPointPatchMapper>(mesh.boundary().size())
        {
            const pointBoundaryMesh& patches = mesh.boundary();

            forAll(patches, patchI)
            {
                set
                (
                    patchI,
                    new topoPointPatchMapper
                    (
                        patches[patchI],
                        mpm,
                        mapper
                    )
                );
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
