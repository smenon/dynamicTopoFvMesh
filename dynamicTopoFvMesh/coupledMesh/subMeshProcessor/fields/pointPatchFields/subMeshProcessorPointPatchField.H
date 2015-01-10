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
    subMeshProcessorPointPatchField

Description
    Customized processor type for subMeshes.

    The intent of this class is to provide processor patch field information
    on subMeshes, but avoiding all forms of communication.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef subMeshProcessorPointPatchField_H
#define subMeshProcessorPointPatchField_H

#include "coupledPointPatchField.H"
#include "subMeshProcessorPointPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class subMeshProcessorPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class subMeshProcessorPointPatchField
:
    public coupledPointPatchField<Type>
{
    // Private data

        //- Local reference cast into the processor patch
        const subMeshProcessorPointPatch& procPatch_;

public:

    //- Runtime type information
    TypeName(subMeshProcessorPointPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        subMeshProcessorPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch and internal field and patch field
        subMeshProcessorPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        subMeshProcessorPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given PointPatchField onto a new patch
        subMeshProcessorPointPatchField
        (
            const subMeshProcessorPointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy
        subMeshProcessorPointPatchField
        (
            const subMeshProcessorPointPatchField<Type>&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<Type> > clone() const
        {
            return autoPtr<pointPatchField<Type> >
            (
                new subMeshProcessorPointPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        subMeshProcessorPointPatchField
        (
            const subMeshProcessorPointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type> > clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<Type> >
            (
                new subMeshProcessorPointPatchField<Type>(*this, iF)
            );
        }

    //- Destructor
    ~subMeshProcessorPointPatchField();

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


        //- Constraint handling

            //- Return the constraint type this pointPatchField implements
            virtual const word& constraintType() const
            {
                return type();
            }


        // Evaluation functions

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            )
            {}

            //- Assume processor patch always collocated
            virtual void swapAddSeparated
            (
                const Pstream::commsTypes commsType,
                Field<Type>&
            ) const
            {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "subMeshProcessorPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
