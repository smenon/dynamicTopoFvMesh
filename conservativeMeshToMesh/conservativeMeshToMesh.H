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

\*---------------------------------------------------------------------------*/

#ifndef conservativeMeshToMesh_H
#define conservativeMeshToMesh_H

#include "fvMesh.H"
#include "meshOps.H"
#include "className.H"
#include "multiThreader.H"
#include "threadHandler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class conservativeMeshToMesh Declaration
\*---------------------------------------------------------------------------*/

class conservativeMeshToMesh
{
    // Private data

        //- Hold references to original meshes
        const fvMesh& meshSrc_;
        const fvMesh& meshTgt_;

        // Source cell indices for post-processing
        Map<label> srcCellIndex_;

        //- Interpolation cells
        IOList<labelList> addressing_;

        //- Interpolation weights
        IOList<scalarField> weights_;

        //- Interpolation centres
        IOList<vectorField> centres_;

        //- Cell addressing
        labelList cellAddressing_;

        //- Mutex for the progress counter
        Mutex ctrMutex_;

        //- Progress counter
        label counter_;

        //- Typedef for convenience
        typedef threadHandler<conservativeMeshToMesh> handler;

        //- Boundary addressing
        labelListList boundaryAddressing_;

    // Private Member Functions

        // Calculate nearest cell addressing
        void calcCellAddressing();

        // Calculate weighting factors
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
            labelList& parents,
            scalarField& weights,
            vectorField& centres
        ) const;

        //- Interpolate internal field values (conservative first-order)
        template<class Type>
        void interpolateInternalFieldConserveFirstOrder
        (
            Field<Type>& tgtF,
            const GeometricField<Type, fvPatchField, volMesh>& srcVf
        ) const;

        //- Interpolate internal field values (conservative)
        //  with supplied gradient
        template<class Type>
        void interpolateInternalFieldConserve
        (
            Field<Type>& tgtF,
            const GeometricField
            <
                Type,
                fvPatchField,
                volMesh
            >& srcVf,
            const GeometricField
            <
                typename outerProduct<vector, Type>::type,
                fvPatchField,
                volMesh
            >& srcVfGrad
        ) const;

        //- Interpolate internal field values (conservative)
        template<class Type>
        void interpolateInternalFieldConserve
        (
            Field<Type>& tgtF,
            const GeometricField<Type, fvPatchField, volMesh>& srcVf
        ) const;

        //- Interpolate internal field values (inverse-distance)
        template<class Type>
        void interpolateInternalFieldInvDist
        (
            Field<Type>& tgtF,
            const GeometricField<Type, fvPatchField, volMesh>& srcVf
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

        typedef FixedList<point, 4> TetPoints;

        // Decompose the input cell using face-centre
        static void decomposeCell
        (
            const cell& polyCell,
            const point& cCentre,
            const faceList& faces,
            const pointField& meshPoints,
            const pointField& faceCentres,
            DynamicList<TetPoints>& decompTets
        );

        // Multi-threaded weighting factor computation
        static void calcAddressingAndWeightsThreaded(void *argument);

public:

    // Declare name of the class and its debug switch
    ClassName("conservativeMeshToMesh");

    //- Enumeration specifying mapping type
    enum method
    {
        CONSERVATIVE,
        INVERSE_DISTANCE,
        CONSERVATIVE_FIRST_ORDER
    };

    // Constructors

        //- Construct from the two meshes assuming there is
        //  an exact mapping between all patches,
        //  with an additional option of being multi-threaded
        conservativeMeshToMesh
        (
            const fvMesh& srcMesh,
            const fvMesh& tgtMesh,
            const label nThreads = 1,
            const bool forceRecalculation = false,
            const bool writeAddressing = false
        );

    // Destructor

        ~conservativeMeshToMesh();

    //- Patch-field interpolation class
    class patchFieldInterpolator
    :
        public fvPatchFieldMapper
    {
        const labelList& directAddressing_;

    public:

        // Constructors

            //- Construct given addressing
            patchFieldInterpolator(const labelList& addr)
            :
                directAddressing_(addr)
            {}


        // Destructor

            virtual ~patchFieldInterpolator()
            {}


        // Member Functions

            virtual label size() const
            {
                return directAddressing_.size();
            }

            virtual label sizeBeforeMapping() const
            {
                return directAddressing_.size();
            }

            virtual bool direct() const
            {
                return true;
            }

            const labelList& directAddressing() const
            {
                return directAddressing_;
            }
    };

    // Member Functions

        // Access

            // Return source mesh
            const fvMesh& srcMesh() const;

            // Return target mesh
            const fvMesh& tgtMesh() const;

            // Fetch cell addressing
            const labelList& cellAddressing() const;

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

#if defined(__INTEL_COMPILER)

// Template specialisation
template<>
class token::Compound<Field<scalar> >
:
    public token::compound,
    public Field<scalar>
{
public:

    //- Runtime type information
    TypeName("Compound<T>");

    Compound(Istream& is)
    :
        Field<scalar>(is)
    {}

    label size() const
    {
        return Field<scalar>::size();
    }

    void write(Ostream& os) const
    {
        operator<<(os, static_cast<const Field<scalar>&>(*this));
    }
};

// Template specialisation
template<>
class token::Compound<Field<vector> >
:
    public token::compound,
    public Field<vector>
{
public:

    //- Runtime type information
    TypeName("Compound<T>");

    Compound(Istream& is)
    :
        Field<vector>(is)
    {}

    label size() const
    {
        return Field<vector>::size();
    }

    void write(Ostream& os) const
    {
        operator<<(os, static_cast<const Field<vector>&>(*this));
    }
};

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "conservativeMeshToMeshInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
