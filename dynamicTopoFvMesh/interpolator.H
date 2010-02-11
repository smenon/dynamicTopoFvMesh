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
    interpolator

Description
    Helper class to perform field interpolation during topo-changes

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    interpolator.C

\*---------------------------------------------------------------------------*/

#ifndef interpolator_H
#define interpolator_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "typeInfo.H"
#include "scalarField.H"
#include "resizableList.H"

namespace Foam
{

// Class forward declarations
class mapPolyMesh;
class dynamicTopoFvMesh;

/*---------------------------------------------------------------------------*\
                        Class interpolator Declaration
\*---------------------------------------------------------------------------*/

class interpolator
{
    // Private data

        //- Reference to the base mesh
        dynamicTopoFvMesh& mesh_;

        //- Flag to determine whether
        //  fields have been registered
        bool fieldsRegistered_;

        //- Mesh motion fluxes
        mutable resizableList<scalar>* meshPhiPtr_;

        //- Old cell volume field
        mutable resizableList<scalar>* V0Ptr_;

    // Private Member Functions

        //- Make mesh motion fluxes
        void makeMeshPhi() const;

        //- Make old cell volumes
        void makeV0() const;

        //- Clear out data
        void clearOut();

        //- Print out mapped fields
        void printFields() const;

public:

    // Declare the name of the class and its debug switch
        TypeName("interpolator");

    // Constructors

        //- Construct from existing IOobject
        explicit interpolator(dynamicTopoFvMesh& mesh);

    // Destructor

        virtual ~interpolator();

    // Member Functions

        // Add a new face entry
        void insertFace
        (
            const scalar sweptVol
        );

        // Add a new cell entry
        void insertCell
        (
            const labelList& mapCells,
            const scalarList& mapWeights
        );

        // Set the old-volume value for an existing cell
        void setOldVolume
        (
            const label cellIndex,
            const scalar oldVolume
        );

        // Set the mesh-flux for an existing face
        void setMeshFlux
        (
            const label faceIndex,
            const scalar sweptVol
        );

        // Get the mesh-flux for an existing face.
        scalar getMeshFlux(const label faceIndex) const;

        // Flip the face-flux for an existing face
        void flipFaceFlux(const label faceIndex);

        // Register fields for interpolation
        void registerFields();

        // Update geometry for a repositioned mesh
        void movePoints(const scalarField& sweptVols) const;

        // Update fields after a topo-change operation
        void updateMesh(const mapPolyMesh& mpm);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
