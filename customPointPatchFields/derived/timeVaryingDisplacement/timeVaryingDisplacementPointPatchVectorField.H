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
    timeVaryingDisplacementPointPatchVectorField

Description
    Provide displacement based on profile file data.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    timeVaryingDisplacementPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingDisplacementPointPatchVectorField_H
#define timeVaryingDisplacementPointPatchVectorField_H

#include "fixedValuePointPatchFields.H"
#include "pointPatchFieldsFwd.H"
#include "interpolationTable.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class timeVaryingDisplacementPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingDisplacementPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    // Private data

        //- Axis of displacement
        vector axis_;

        //- Profile file containing displacement data
        fileName profileFile_;

        //- Interpolation table
        interpolationTable<scalar> table_;

        //- Engine time option
        Switch engineTime_;

        //- RPM for engine time
        scalar rpm_;

public:

    //- Runtime type information
    TypeName("timeVaryingDisplacement");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        timeVaryingDisplacementPointPatchVectorField
        (
            const timeVaryingDisplacementPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchVectorField> clone() const
        {
            return autoPtr<pointPatchVectorField>
            (
                new timeVaryingDisplacementPointPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingDisplacementPointPatchVectorField
        (
            const timeVaryingDisplacementPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchVectorField> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchVectorField>
            (
                new timeVaryingDisplacementPointPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
