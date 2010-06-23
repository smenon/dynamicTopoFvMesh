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
    topoMapper

Description
    Class holds all necessary information for mapping fields associated with
    dynamicTopoFvMesh and fvMesh.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    topoMapper.C

\*---------------------------------------------------------------------------*/

#ifndef topoMapper_H
#define topoMapper_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class faceMapper;
class topoCellMapper;
class fvSurfaceMapper;
class fvBoundaryMeshMapper;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                         Class topoMapper Declaration
\*---------------------------------------------------------------------------*/

class topoMapper
{
    // Private data

        //- Reference to fvMesh
        const fvMesh& mesh_;

        //- Face mapper
        autoPtr<faceMapper> faceMap_;

        //- Cell mapper
        autoPtr<topoCellMapper> cellMap_;

        //- Surface mapper
        autoPtr<fvSurfaceMapper> surfaceMap_;

        //- Boundary mapper
        autoPtr<fvBoundaryMeshMapper> boundaryMap_;

    // Demand-driven private data

        //- Old cell-centre information
        mutable volVectorField* oldCellCentresPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        topoMapper(const topoMapper&);

        //- Disallow default bitwise assignment
        void operator=(const topoMapper&);

        //- Clear out local storage
        void clearOut();

public:

    // Constructors

        //- Construct from dynamicTopoFvMesh
        topoMapper(const fvMesh& mesh)
        :
            mesh_(mesh),
            faceMap_(NULL),
            cellMap_(NULL),
            surfaceMap_(NULL),
            boundaryMap_(NULL),
            oldCellCentresPtr_(NULL)
        {}

    // Destructor

        ~topoMapper();

    // Member Functions

        //- Return reference to the mesh
        const fvMesh& mesh() const;

        //- Return reference to objectRegistry storing fields.
        const objectRegistry& db() const;

        //- Set mapping information
        void setMapper(const mapPolyMesh& mpm);

        //- Set face weighting information
        void setFaceWeights
        (
			Map<scalarField>& weights,
			Map<vectorField>& centres
		);

        //- Set cell weighting information
        void setCellWeights
        (
			Map<scalarField>& weights,
			Map<vectorField>& centres
		);

        //- Set old cell-centre information
        void setOldCellCentres(const volVectorField& oldCentres);

        //- Return old cell-centre information
        const volVectorField& oldCentres() const;

        //- Correct fluxes after topology change
        void correctFluxes();

        //- Return volume mapper
        const FieldMapper& volMap() const;

        //- Return surface mapper
        const fvSurfaceMapper& surfaceMap() const;

        //- Return boundary mapper
        const fvBoundaryMeshMapper& boundaryMap() const;

        //- Clear out member data
        void clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#include "topoMapper.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
