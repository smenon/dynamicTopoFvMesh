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
#include "Map.H"
#include "fieldTypes.H"
#include "scalarField.H"
#include "labelHashSet.H"

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

        // List of maps for each primitive-type
        HashTable<Map<scalar> >          volScalarMap_;
        HashTable<Map<vector> >          volVectorMap_;
        HashTable<Map<sphericalTensor> > volSphericalTensorMap_;
        HashTable<Map<symmTensor> >      volSymmTensorMap_;
        HashTable<Map<tensor> >          volTensorMap_;

        HashTable<Map<scalar> >          surfScalarMap_;
        HashTable<Map<vector> >          surfVectorMap_;
        HashTable<Map<sphericalTensor> > surfSphericalTensorMap_;
        HashTable<Map<symmTensor> >      surfSymmTensorMap_;
        HashTable<Map<tensor> >          surfTensorMap_;

        labelHashSet flipFaces_;

    // Private Member Functions

        //- Clear out data
        void clearOut();

        //- Post-processing
        void writeVTK(const word& name);

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
            const label patch,
            const label newFaceIndex,
            const labelList& mapFaces,
            const scalarField& mapWeights
        );

        // Remove the face, if it exists in the map.
        void removeFace(const label index);

        // Add a new cell entry
        void insertCell
        (
            const label newCellIndex,
            const labelList& mapCells,
            const scalarField& mapWeights
        );

        // Remove the cell, if it exists in the map.
        void removeCell(const label index);

        // Set the volume-flux for an existing face
        void setPhi
        (
            const label faceIndex,
            const scalar facePhi
        );

        // Interpolate flux for an existing face
        void interpolatePhi
        (
            const label cellIndex,
            const label faceIndex,
            const vector& Sf
        );

        //- Post-processing
        void writeFluxes(const word& name);

        // Set a particular face index as flipped.
        void setFlip(const label fIndex);

        // Register fields for interpolation
        void registerFields();

        // Update fields after a topo-change operation
        void updateMesh(const mapPolyMesh& mpm);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Define a few macros for code brevity.
#include "interpolatorM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
