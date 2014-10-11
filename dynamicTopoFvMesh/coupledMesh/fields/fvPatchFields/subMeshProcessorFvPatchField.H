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
    subMeshProcessorFvPatchField

Description
    Customized processor type for subMeshes.

    The intent of this class is to provide processor patch field information
    on subMeshes, but avoiding all forms of communication.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef subMeshProcessorFvPatchField_H
#define subMeshProcessorFvPatchField_H

#include "coupledFvPatchField.H"
#include "subMeshProcessorFvPatch.H"
#include "processorLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class subMeshProcessorFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class subMeshProcessorFvPatchField
:
    public processorLduInterfaceField,
    public coupledFvPatchField<Type>
{
    // Private data

        //- Local reference cast into the processor patch
        const subMeshProcessorFvPatch& procPatch_;

public:

    //- Runtime type information
    TypeName(subMeshProcessorFvPatch::typeName_());

    // Constructors

        //- Construct from patch and internal field
        subMeshProcessorFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        subMeshProcessorFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        subMeshProcessorFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField onto a new patch
        subMeshProcessorFvPatchField
        (
            const subMeshProcessorFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        subMeshProcessorFvPatchField
        (
            const subMeshProcessorFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new subMeshProcessorFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        subMeshProcessorFvPatchField
        (
            const subMeshProcessorFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new subMeshProcessorFvPatchField<Type>(*this, iF)
            );
        }

    //- Destructor
    ~subMeshProcessorFvPatchField();

    // Member functions

        // Access

            //- Return true if running parallel
            virtual bool coupled() const
            {
                if (Pstream::parRun())
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }

            //- Return neighbour field given internal field
            virtual tmp<Field<Type> > patchNeighbourField() const;

        // Evaluation functions

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

            //- Evaluate the patch field
            virtual void evaluate(const Pstream::commsTypes commsType);

            //- Return patch-normal gradient
            virtual tmp<Field<Type> > snGrad
            (
                const scalarField& deltaCoeffs
            ) const;

            //- Is all data available
            virtual bool ready() const;

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                scalarField& result,
                const scalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                scalarField& result,
                const scalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                Field<Type>& result,
                const Field<Type>& psiInternal,
                const scalarField& coeffs,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>& result,
                const Field<Type>& psiInternal,
                const scalarField& coeffs,
                const Pstream::commsTypes commsType
            ) const;

        //- Processor coupled interface functions

            //- Return communicator used for comms
            virtual int comm() const
            {
                return procPatch_.comm();
            }

            //- Return processor number
            virtual int myProcNo() const
            {
                return procPatch_.myProcNo();
            }

            //- Return neigbour processor number
            virtual int neighbProcNo() const
            {
                return procPatch_.neighbProcNo();
            }

            //- Does the patch field perform the transformation
            virtual bool doTransform() const
            {
                return !(procPatch_.parallel() || pTraits<Type>::rank == 0);
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return procPatch_.forwardT();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return pTraits<Type>::rank;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "subMeshProcessorFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
