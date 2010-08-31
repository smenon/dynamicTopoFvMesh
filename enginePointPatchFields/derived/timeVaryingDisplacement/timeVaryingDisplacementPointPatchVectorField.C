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

\*---------------------------------------------------------------------------*/

#include "timeVaryingDisplacementPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

timeVaryingDisplacementPointPatchVectorField::
timeVaryingDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(p, iF),
    axis_(vector::zero),
    profileFile_("NoFileSpecified")
{}


timeVaryingDisplacementPointPatchVectorField::
timeVaryingDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchVectorField(p, iF, dict),
    axis_(dict.lookup("axis")),
    profileFile_(dict.lookup("profileFile")),
    table_(profileFile_)
{
    // Normalize the axis
    axis_ /= mag(axis_) + VSMALL;

    // Set outOfBounds handling to clamp
    table_.outOfBounds(interpolationTable<scalar>::CLAMP);
}


timeVaryingDisplacementPointPatchVectorField::
timeVaryingDisplacementPointPatchVectorField
(
    const timeVaryingDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const PointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchVectorField(ptf, p, iF, mapper),
    axis_(ptf.axis_),
    profileFile_(ptf.profileFile_),
    table_(ptf.table_)
{}


timeVaryingDisplacementPointPatchVectorField::
timeVaryingDisplacementPointPatchVectorField
(
    const timeVaryingDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(ptf, iF),
    axis_(ptf.axis_),
    profileFile_(ptf.profileFile_),
    table_(ptf.table_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void timeVaryingDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->dimensionedInternalField().mesh()();
    const Time& t = mesh.time();

    // Fetch old / new displacement values
    scalar newVal = table_(t.value());
    scalar oldVal = table_(t.value() - t.deltaT().value());

    Field<vector>::operator= ( axis_ * (newVal - oldVal) );

    fixedValuePointPatchVectorField::updateCoeffs();
}


void timeVaryingDisplacementPointPatchVectorField::write(Ostream& os) const
{
    os.writeKeyword("axis")
        << axis_ << token::END_STATEMENT << nl;

    os.writeKeyword("profileFile")
        << profileFile_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    timeVaryingDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
