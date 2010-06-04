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
    conservativeMeshToMesh

Description
    Conservative mesh to mesh interpolation class.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    conservativeMeshToMesh.C

\*---------------------------------------------------------------------------*/

#ifndef conservativeMeshToMesh_H
#define conservativeMeshToMesh_H

#include "meshToMesh.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class conservativeMeshToMesh Declaration
\*---------------------------------------------------------------------------*/

class conservativeMeshToMesh
:
    public meshToMesh
{
    // Private data

        //- Interpolation cells
        IOList<labelList> addressing_;

        //- Interpolation weights
        IOList<scalarField> weights_;

        //- Interpolation centres
        IOList<vectorField> centres_;

    // Private Member Functions

        void calcAddressingAndWeights
        (
            const label cellStart,
            const label cellSize
        );

        // Compute cell weighting factors for a particular cell
        bool computeCellWeights
        (
            const label newCellIndex,
            const scalar newCellVolume,
            const label oldCandidate,
            const scalar searchFactor,
            labelList& parents,
            scalarField& weights,
            vectorField& centres
        ) const;

        // Obtain a list of possible parent cells from the old mesh.
        labelList cellParents
        (
            const label newCellIndex,
            const scalar searchFactor,
            const label oldCandidate
        ) const;

        // Return the intersection volume between cells in old/new meshes
        bool cellIntersection
        (
            const label newCellIndex,
            const label oldCellIndex,
            scalar& intVolume,
            vector& intCentre
        ) const;

        // Determine whether the particular point lies
        // inside the given cell
        bool pointInCell
        (
            const label cellIndex,
            const point& checkPoint,
            const bool useFromMesh
        ) const;

        // Determine whether a particular line segment
        // intersects with a given face
        bool segmentFaceIntersection
        (
            const FixedList<vector,2>& segment,
            const edge& segmentIndices,
            const label faceIndex,
            const scalar tolFactor,
            point& intPoint,
            const bool useFromMesh
        ) const;

        // Determine whether the particular point lies
        // inside the given face
        bool pointInFace
        (
            const label faceIndex,
            const point& checkPoint,
            const bool useFromMesh
        ) const;

        // Compute the volume / centre of a polyhedron
        // formed by a convex set of points.
        void convexSetVolume
        (
            const label newCellIndex,
            const label oldCellIndex,
            const scalar tolFraction,
            const vectorField& cvxSet,
            scalar& cVolume,
            vector& cCentre
        ) const;

        // Method to insert a label between two labels in a face
        void insertLabel
        (
            const label newLabel,
            const label labelA,
            const label labelB,
            labelList& list
        ) const;

    // Private debugging routines

        void writeVTK
        (
            const word& name,
            const fvMesh& mesh,
            const label entity,
            const label primitiveType
        ) const;

        void writeVTK
        (
            const word& name,
            const fvMesh& mesh,
            const labelList& cList,
            const label primitiveType
        ) const;

        void writeVTK
        (
            const word& name,
            const label nPoints,
            const label nCells,
            const label nTotalCells,
            const pointField& points,
            const labelListList& cpList,
            const label primitiveType,
            const Map<label>& reversePointMap = Map<label>(),
            const Map<label>& reverseCellMap = Map<label>()
        ) const;

    // Private static members

        static void calcAddressingAndWeightsThreaded(void *argument);

public:

    // Declare name of the class and its debug switch
    ClassName("conservativeMeshToMesh");

    // Constructors

        //- Construct from the two meshes assuming there is
        //  an exact mapping between all patches,
        //  with an additional option of being multi-threaded
        conservativeMeshToMesh
        (
            const fvMesh& fromMesh,
            const fvMesh& toMesh,
            const label nThreads = 1,
            const bool forceRecalculation = false,
            const bool writeAddressing = false
        );

    // Destructor

        ~conservativeMeshToMesh();


    //- Patch-field interpolation class
    class conservativePatchFieldInterpolator
    :
        public fvPatchFieldMapper
    {

    public:

        // Constructors

            //- Construct given addressing
            conservativePatchFieldInterpolator()
            {}

        // Destructor

            virtual ~conservativePatchFieldInterpolator()
            {}

        // Member Functions

            virtual label size() const
            {
                return 0;
            }

            virtual label sizeBeforeMapping() const
            {
                return 0;
            }

            virtual bool direct() const
            {
                return false;
            }
    };


    // Member Functions

        //- Interpolate internal field values
        template<class Type>
        void interpolateInternalField
        (
            Field<Type>& toF,
            const GeometricField<Type, fvPatchField, volMesh>& fromVf
        ) const;

        //- Interpolate volume field
        template<class Type>
        void interpolate
        (
            GeometricField<Type, fvPatchField, volMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        template<class Type>
        void interpolate
        (
            GeometricField<Type, fvPatchField, volMesh>&,
            const tmp<GeometricField<Type, fvPatchField, volMesh> >&
        ) const;

        //- Interpolate volume field
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh> > interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh> > interpolate
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh> >&
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "conservativeMeshToMeshInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
