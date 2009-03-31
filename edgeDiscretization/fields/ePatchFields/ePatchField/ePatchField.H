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
    ePatchField

Description
    ePatchField<Type> abstract base class.  This class gives a fat-interface
    to all derived classes covering all possible ways in which they might be
    used.

SourceFiles
    ePatchField.C

\*---------------------------------------------------------------------------*/

#ifndef ePatchField_H
#define ePatchField_H

#include "ePatch.H"
#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

class objectRegistry;
class dictionary;
// class ePatchFieldMapper;
class EdgeMesh;


// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class ePatchField;

template<class Type>
Ostream& operator<<(Ostream&, const ePatchField<Type>&);

/*---------------------------------------------------------------------------*\
                         Class ePatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ePatchField
:
    public Field<Type>
{
    // Private data

        //- Reference to a patch
        const ePatch& patch_;

        //- Reference to internal field
        const DimensionedField<Type, EdgeMesh>& internalField_;

        //- Update index used so that updateCoeffs is called only once during
        //  the construction of the matrix
        bool updated_;


public:

    typedef ePatch Patch;


    //- Runtime type information
    TypeName("ePatchField");

    //- Debug switch to disallow the use of 
    static int disallowDefaultEPatchField;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            ePatchField,
            patch,
            (
                const ePatch& p,
                const DimensionedField<Type, EdgeMesh>& iF
            ),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            ePatchField,
            dictionary,
            (
                const ePatch& p,
                const DimensionedField<Type, EdgeMesh>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        ePatchField
        (
            const ePatch&,
            const DimensionedField<Type, EdgeMesh>&
        );

        //- Construct from patch and internal field and patch field
        ePatchField
        (
            const ePatch&,
            const DimensionedField<Type, EdgeMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        ePatchField
        (
            const ePatch&,
            const DimensionedField<Type, EdgeMesh>&,
            const dictionary&
        );

        //- Construct as copy
        ePatchField(const ePatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<ePatchField<Type> > clone() const
        {
            return tmp<ePatchField<Type> >(new ePatchField<Type>(*this));
        }

        //- Construct as copy setting internal field reference
        ePatchField
        (
            const ePatchField<Type>&,
            const DimensionedField<Type, EdgeMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<ePatchField<Type> > clone
        (
            const DimensionedField<Type, EdgeMesh>& iF
        ) const
        {
            return tmp<ePatchField<Type> >(new ePatchField<Type>(*this, iF));
        }


    // Selectors

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static tmp<ePatchField<Type> > New
        (
            const word&,
            const ePatch&,
            const DimensionedField<Type, EdgeMesh>&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static tmp<ePatchField<Type> > New
        (
            const ePatch&,
            const DimensionedField<Type, EdgeMesh>&,
            const dictionary&
        );

        //- Return a pointer to a new calculatedEPatchField created on
        //  freestore without setting patchField values
        /*
        template<class Type2>
        static tmp<ePatchField<Type> > NewCalculatedType
        (
            const ePatchField<Type2>&
        );
        */

    // Destructor

        virtual ~ePatchField<Type>()
        {}


    // Member functions

        // Access

            //- Return local objectRegistry
            const objectRegistry& db() const;

            //- Return patch
            const ePatch& patch() const
            {
                return patch_;
            }

            //- Return dimensioned internal field reference
            const DimensionedField<Type, EdgeMesh>&
            dimensionedInternalField() const
            {
                return internalField_;
            }

            //- Return internal field reference
            const Field<Type>& internalField() const
            {
                return internalField_;
            }

            //- Return the type of the calculated for of ePatchField
            static const word& calculatedType();

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Return true if this patch field is coupled
            virtual bool coupled() const
            {
                return false;
            }

            //- Return true if the boundary condition has already been updated
            bool updated() const
            {
                return updated_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            /*
            virtual void autoMap
            (
                const ePatchFieldMapper&
            );

            //- Reverse map the given ePatchField onto this ePatchField
            virtual void rmap
            (
                const ePatchField<Type>&,
                const labelList&
            );
            */

        // Evaluation functions

            //- Return patch-normal gradient
            // virtual tmp<Field<Type> > snGrad() const;

            //- Update the coefficients associated with the patch field
            //  Sets Updated to true
            virtual void updateCoeffs()
            {
                updated_ = true;
            }

            //- Return internal field next to patch as patch field
            // virtual tmp<Field<Type> > patchInternalField() const;

            //- Return patchField on the opposite patch of a coupled patch
            virtual tmp<Field<Type> > patchNeighbourField() const
            {
                notImplemented(type() + "patchNeighbourField()");
                return *this;
            }

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            )
            {}

            //- Evaluate the patch field, sets Updated to false
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );


            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueInternalCoeffs
            (
                const tmp<Field<scalar> >&
            ) const
            {
                notImplemented
                (
                    type()
                  + "::valueInternalCoeffs(const tmp<Field<scalar> >&)"
                );
                return *this;
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueBoundaryCoeffs
            (
                const tmp<Field<scalar> >&
            ) const
            {
                notImplemented
                (
                    type()
                  + "::valueBoundaryCoeffs(const tmp<Field<scalar> >&)"
                );
                return *this;
            }

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type> > gradientInternalCoeffs() const
            {
                notImplemented(type() + "::gradientInternalCoeffs()");
                return *this;
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type> > gradientBoundaryCoeffs() const
            {
                notImplemented(type() + "::gradientBoundaryCoeffs()");
                return *this;
            }


        //- Write
        virtual void write(Ostream&) const;


        // Check

            //- Check ePatchField<Type> against given ePatchField<Type>
            void check(const ePatchField<Type>&) const;


    // Member operators

        virtual void operator=(const UList<Type>&);

        virtual void operator=(const ePatchField<Type>&);
        virtual void operator+=(const ePatchField<Type>&);
        virtual void operator-=(const ePatchField<Type>&);
        virtual void operator*=(const ePatchField<scalar>&);
        virtual void operator/=(const ePatchField<scalar>&);

        virtual void operator+=(const Field<Type>&);
        virtual void operator-=(const Field<Type>&);

        virtual void operator*=(const Field<scalar>&);
        virtual void operator/=(const Field<scalar>&);

        virtual void operator=(const Type&);
        virtual void operator+=(const Type&);
        virtual void operator-=(const Type&);
        virtual void operator*=(const scalar);
        virtual void operator/=(const scalar);


        // Force an assignment irrespective of form of patch

        virtual void operator==(const ePatchField<Type>&);
        virtual void operator==(const Field<Type>&);
        virtual void operator==(const Type&);


    // Ostream operator

        friend Ostream& operator<< <Type>(Ostream&, const ePatchField<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ePatchField.C"
//#   include "calculatedEPatchField.H"
#endif


#define makeEPatchTypeFieldTypeName(typePatchTypeField)                    \
                                                                           \
defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);

#define makeEPatchFieldsTypeName(typePatchField)                           \
                                                                           \
makeEPatchTypeFieldTypeName(typePatchField##EPatchScalarField);            \
makeEPatchTypeFieldTypeName(typePatchField##EPatchVectorField);            \
makeEPatchTypeFieldTypeName(typePatchField##EPatchSphericalTensorField);   \
makeEPatchTypeFieldTypeName(typePatchField##EPatchSymmTensorField);        \
makeEPatchTypeFieldTypeName(typePatchField##EPatchTensorField);

#define makeEPatchTypeField(PatchTypeField, typePatchTypeField)            \
                                                                           \
defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);                \
                                                                           \
addToRunTimeSelectionTable                                                 \
(                                                                          \
    PatchTypeField, typePatchTypeField, patch                              \
);                                                                         \
                                                                           \
addToRunTimeSelectionTable                                                 \
(                                                                          \
    PatchTypeField, typePatchTypeField, dictionary                         \
);


#define makeEPatchFields(type)                                             \
                                                                           \
makeEPatchTypeField(ePatchScalarField, type##EPatchScalarField);           \
makeEPatchTypeField(ePatchVectorField, type##EPatchVectorField);           \
makeEPatchTypeField                                                        \
(                                                                          \
    ePatchSphericalTensorField,                                            \
    type##EPatchSphericalTensorField                                       \
);                                                                         \
makeEPatchTypeField(ePatchSymmTensorField, type##EPatchSymmTensorField);   \
makeEPatchTypeField(ePatchTensorField, type##EPatchTensorField);


#define makeEPatchTypeFieldTypedefs(type)                                  \
                                                                           \
typedef type##EPatchField<scalar> type##EPatchScalarField;                 \
typedef type##EPatchField<vector> type##EPatchVectorField;                 \
typedef type##EPatchField<sphericalTensor>                                 \
     type##EPatchSphericalTensorField;                                     \
typedef type##EPatchField<symmTensor> type##EPatchSymmTensorField;         \
typedef type##EPatchField<tensor> type##EPatchTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
