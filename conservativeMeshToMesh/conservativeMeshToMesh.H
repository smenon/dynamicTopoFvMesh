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
    conservativeMeshToMeshAddressing.C
    conservativeMeshToMeshInterpolate.C
    conservativeMeshToMeshWriteVTK.C

\*---------------------------------------------------------------------------*/

#ifndef conservativeMeshToMesh_H
#define conservativeMeshToMesh_H

#include "meshOps.H"
#include "className.H"
#include "meshToMesh.H"
#include "tetPolyMesh.H"
#include "multiThreader.H"
#include "threadHandler.H"

// Use the MPFR library, if necessary
#define USE_MPFR 1

#if USE_MPFR
#include "multiPrecision.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class conservativeMeshToMesh Declaration
\*---------------------------------------------------------------------------*/

class conservativeMeshToMesh
{
    // Private data

        //- Pointers to meshToMesh interpolation
        autoPtr<meshToMesh> meshToMeshPtr_;

        //- Hold references to the original mesh
        const fvMesh& meshFrom_;
        const fvMesh& meshTo_;

        //- Pointers for optional decomposition
        autoPtr<fvMesh> srcTetFvMesh_, tgtTetFvMesh_;
        autoPtr<tetPolyMesh> srcTetMesh_, tgtTetMesh_;

        //- Addressing from polyhedra to tets
        labelList srcTetStarts_, srcTetSizes_;
        labelList tgtTetStarts_, tgtTetSizes_;

        //- Interpolation cells
        IOList<labelList> addressing_;

        //- Interpolation weights
        IOList<scalarField> weights_;

        //- Interpolation volumes
        IOList<scalarField> volumes_;

        //- Interpolation centres
        IOList<vectorField> centres_;

        //- Mutex for the progress counter
        Mutex ctrMutex_;

        //- Progress counter
        label counter_;

        //- 2D / 3D mesh characterstics
        bool twoDMesh_;

        //- Typedef for convenience
        typedef threadHandler<conservativeMeshToMesh> handler;

        //- Boundary addressing
        labelListList boundaryAddressing_;

    // Private Member Functions

        // Decompose a given mesh, using tetDecomposition
        static const fvMesh& decomposeMesh
        (
            const fvMesh& initialMesh,
            const bool decomp,
            autoPtr<tetPolyMesh>& tetDecomp,
            autoPtr<fvMesh>& tetDecompFvMesh,
            labelList& tetStarts,
            labelList& tetSizes
        );

        void calcAddressingAndWeights
        (
            const label cellStart,
            const label cellSize,
            bool report = false
        );

        // Invert addressing from source to target
        bool invertAddressing();

        // Compute weighting factors for a particular cell
        bool computeWeights
        (
            const label index,
            const label oldCandidate,
            const labelListList& oldNeighbourList,
            const scalar mTol,
            label& precisionAttempts,
            labelList& parents,
            scalarField& weights,
            scalarField& volumes,
            vectorField& centres,
            bool highPrecision = false
        ) const;

        // Return the intersection points between cells in old/new meshes
        template <class T>
        bool cellIntersection
        (
            const label newIndex,
            const label oldIndex,
            const T& matchTol,
            Field<Vector<T> >& intPoints,
            bool output = false
        ) const;

        // Find all indices in edge
        bool findInEdge
        (
            const labelList& cList,
            const edge& edgeToCheck
        ) const;

        // Compute the volume / centre of a polyhedron
        // formed by a convex set of points.
        void convexSetVolume
        (
            const label newCellIndex,
            const label oldCellIndex,
            const vectorField& cvxSet,
            scalar& cVolume,
            vector& cCentre,
            faceList& vFaces,
            bool output = false
        ) const;

        //- Interpolate internal field values from decomposition
        template<class Type>
        void interpolateInternalFieldDecomp
        (
            Field<Type>& toF,
            const GeometricField<Type, fvPatchField, volMesh>& fromVf
        ) const;

        //- Interpolate internal field values (conservative first-order)
        template<class Type>
        void interpolateInternalFieldConserveFirstOrder
        (
            Field<Type>& toF,
            const GeometricField<Type, fvPatchField, volMesh>& fromVf
        ) const;

        //- Interpolate internal field values (conservative)
        //  with supplied gradient
        template<class Type>
        void interpolateInternalFieldConserve
        (
            Field<Type>& toF,
            const GeometricField
            <
                Type,
                fvPatchField,
                volMesh
            >& fromVf,
            const GeometricField
            <
                typename outerProduct<vector, Type>::type,
                fvPatchField,
                volMesh
            >& fromgVf
        ) const;

        //- Interpolate internal field values (conservative)
        template<class Type>
        void interpolateInternalFieldConserve
        (
            Field<Type>& toF,
            const GeometricField<Type, fvPatchField, volMesh>& fromVf
        ) const;

        //- Interpolate internal field values (inverse-distance)
        template<class Type>
        void interpolateInternalFieldInvDist
        (
            Field<Type>& toF,
            const GeometricField<Type, fvPatchField, volMesh>& fromVf
        ) const;

        // Output an entity as a VTK file
        void writeVTK
        (
            const word& name,
            const label entity,
            const label primitiveType = 3,
            const bool useOldConnectivity = false
        ) const;

        // Output a list of entities as a VTK file
        void writeVTK
        (
            const word& name,
            const labelList& cList,
            const label primitiveType = 3,
            const bool useOldConnectivity = false,
            const UList<scalar>& field = UList<scalar>()
        ) const;

    // Private static members

        static void calcAddressingAndWeightsThreaded(void *argument);

public:

    // Declare name of the class and its debug switch
    ClassName("conservativeMeshToMesh");

    //- Enumeration specifying mapping type
    enum method
    {
        CONSERVATIVE,
        INVERSE_DISTANCE,
        CONSERVATIVE_FIRST_ORDER,
        DECOMP
    };

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
            const bool writeAddressing = false,
            const bool decompSource = false,
            const bool decompTarget = false
        );

    // Destructor

        ~conservativeMeshToMesh();


    // Member Functions

        // Access

            // Return (decomposed) source mesh
            const fvMesh& fromMesh() const;

            // Return (decomposed) target mesh
            const fvMesh& toMesh() const;

            // Return basic source mesh (undecomposed)
            const fvMesh& origSrcMesh() const;

            // Return basic target mesh (undecomposed)
            const fvMesh& origTgtMesh() const;

            // Return (decomposed) source mesh
            const fvMesh& srcMesh() const;

            // Return (decomposed) target mesh
            const fvMesh& tgtMesh() const;

        // Interpolation

            template<class Type>
            void interpolateInternalField
            (
                Field<Type>&,
                const GeometricField<Type, fvPatchField, volMesh>&,
                const label method
            ) const;

            template<class Type>
            void interpolateInternalField
            (
                Field<Type>&,
                const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
                const label method
            ) const;

            //- Interpolate volume field with a supplied gradient
            template<class Type>
            void interpolate
            (
                GeometricField<Type, fvPatchField, volMesh>&,
                const GeometricField<Type, fvPatchField, volMesh>&,
                const GeometricField
                <
                    typename outerProduct<vector, Type>::type,
                    fvPatchField,
                    volMesh
                >&,
                const label method
            ) const;

            //- Interpolate volume field
            template<class Type>
            void interpolate
            (
                GeometricField<Type, fvPatchField, volMesh>&,
                const GeometricField<Type, fvPatchField, volMesh>&,
                const label method
            ) const;

            template<class Type>
            void interpolate
            (
                GeometricField<Type, fvPatchField, volMesh>&,
                const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
                const label method
            ) const;

            //- Interpolate volume field
            template<class Type>
            tmp<GeometricField<Type, fvPatchField, volMesh> > interpolate
            (
                const GeometricField<Type, fvPatchField, volMesh>&,
                const label method
            ) const;

            template<class Type>
            tmp<GeometricField<Type, fvPatchField, volMesh> > interpolate
            (
                const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
                const label method
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
