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
class topoCellMapper;
class topoSurfaceMapper;
class topoBoundaryMeshMapper;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                         Class topoMapper Declaration
\*---------------------------------------------------------------------------*/

class topoMapper
{
    // Private data

        //- Reference to fvMesh
        const fvMesh& mesh_;

        //- Cell mapper
        autoPtr<topoCellMapper> cellMap_;

        //- Surface mapper
        autoPtr<topoSurfaceMapper> surfaceMap_;

        //- Boundary mapper
        autoPtr<topoBoundaryMeshMapper> boundaryMap_;

    // Demand-driven private data

        // Stored gradients for mapping
        HashTable<autoPtr<volVectorField> > sGrads_;
        HashTable<autoPtr<volTensorField> > vGrads_;

        //- Old cell-centre information
        mutable volVectorField* oldCellCentresPtr_;

        //- Intersection volume weights
        Map<scalarField> faceWeights_;
        Map<scalarField> cellWeights_;

        //- Intersection centre weights
        Map<vectorField> faceCentres_;
        Map<vectorField> cellCentres_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        topoMapper(const topoMapper&);

        //- Disallow default bitwise assignment
        void operator=(const topoMapper&);

        // Store gradients of volFields on the mesh
        // prior to topology changes
        template <class Type, class gradType>
        void storeGradients
        (
            HashTable<autoPtr<gradType> >& gradTable
        );

        //- Clear out local storage
        void clearOut();

public:

    // Constructors

        //- Construct from dynamicTopoFvMesh
        topoMapper(const fvMesh& mesh)
        :
            mesh_(mesh),
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

        //- Fetch face weights
        const Map<scalarField>& faceWeights() const;

        //- Fetch cell weights
        const Map<scalarField>& cellWeights() const;

        //- Fetch face centres
        const Map<vectorField>& faceCentres() const;

        //- Fetch cell centres
        const Map<vectorField>& cellCentres() const;

        //- Store gradients prior to mesh reset
        void storeGradients();

        //- Set old cell-centre information
        void setOldCellCentres(const volVectorField& oldCentres);

        //- Return old cell-centre information
        const volVectorField& oldCentres() const;

        //- Fetch the gradient field
        template <class Type>
        const Type& gradient(const word& name) const;

        //- Correct fluxes after topology change
        void correctFluxes();

        //- Return volume mapper
        const topoCellMapper& volMap() const;

        //- Return surface mapper
        const topoSurfaceMapper& surfaceMap() const;

        //- Return boundary mapper
        const topoBoundaryMeshMapper& boundaryMap() const;

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
