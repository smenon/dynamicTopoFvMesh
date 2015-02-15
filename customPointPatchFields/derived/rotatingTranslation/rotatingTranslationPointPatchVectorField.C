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
    rotatingTranslationPointPatchVectorField

Description
    Provide rotating / translating displacement for point fields.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "rotatingTranslationPointPatchVectorField.H"
#include "pointPatchFields.H"
//#include "RodriguesRotation.H"
#include "transformField.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Rotate / translate points
void rotatingTranslationPointPatchVectorField::rotatePoints
(
    pointField& tPoints
)
{
    vector p(vector::zero), q(vector::zero), r(rotAxis_);

    // Convert to radians
    scalar theta = rotMag_ * (constant::mathematical::pi/180.0);

    // Fetch sines / cosines
    scalar costheta = Foam::cos(theta);
    scalar sintheta = Foam::sin(theta);

    forAll(tPoints, pointI)
    {
        q = vector::zero;

        // Fetch the old point and translate it
        p = tPoints[pointI];

        // Translate to the origin
        p -= rotPoint_;

        // Apply the rotation matrix
        q.x() += (costheta + (1 - costheta) * r.x() * r.x()) * p.x();
        q.x() += ((1 - costheta) * r.x() * r.y() - r.z() * sintheta) * p.y();
        q.x() += ((1 - costheta) * r.x() * r.z() + r.y() * sintheta) * p.z();

        q.y() += ((1 - costheta) * r.x() * r.y() + r.z() * sintheta) * p.x();
        q.y() += (costheta + (1 - costheta) * r.y() * r.y()) * p.y();
        q.y() += ((1 - costheta) * r.y() * r.z() - r.x() * sintheta) * p.z();

        q.z() += ((1 - costheta) * r.x() * r.z() - r.y() * sintheta) * p.x();
        q.z() += ((1 - costheta) * r.y() * r.z() + r.x() * sintheta) * p.y();
        q.z() += (costheta + (1 - costheta) * r.z() * r.z()) * p.z();

        // Translate back to original location
        q += rotPoint_;

        // Assign to the mesh
        tPoints[pointI] = q + transVec_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rotatingTranslationPointPatchVectorField::
rotatingTranslationPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(p, iF),
    rotPoint_(vector::zero),
    rotAxis_(vector::zero),
    transVec_(vector::zero),
    rotMag_(0.0)
{}


rotatingTranslationPointPatchVectorField::
rotatingTranslationPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchVectorField(p, iF, dict),
    rotPoint_(dict.lookup("rotPoint")),
    rotAxis_(dict.lookup("rotAxis")),
    transVec_(dict.lookup("transVec")),
    rotMag_(readScalar(dict.lookup("rotMag")))
{
    // Normalize the rotation axis
    rotAxis_ /= mag(rotAxis_) + VSMALL;
}


rotatingTranslationPointPatchVectorField::
rotatingTranslationPointPatchVectorField
(
    const rotatingTranslationPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchVectorField(ptf, p, iF, mapper),
    rotPoint_(ptf.rotPoint_),
    rotAxis_(ptf.rotAxis_),
    transVec_(ptf.transVec_),
    rotMag_(ptf.rotMag_)
{}


rotatingTranslationPointPatchVectorField::
rotatingTranslationPointPatchVectorField
(
    const rotatingTranslationPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(ptf, iF),
    rotPoint_(ptf.rotPoint_),
    rotAxis_(ptf.rotAxis_),
    transVec_(ptf.transVec_),
    rotMag_(ptf.rotMag_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void rotatingTranslationPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Compute transform tensor by Rodrigues rotation
    // tensor rrTensor = RodriguesRotation(rotAxis_, rotMag_);

    pointField oldPoints = this->patch().localPoints();
    pointField transformPoints = this->patch().localPoints();

    // Rotate points
    rotatePoints(transformPoints);

    // Assign displacement vector
    Field<vector>::operator= (transformPoints - oldPoints);

    // Update rotation point
    rotPoint_ += transVec_;

    fixedValuePointPatchVectorField::updateCoeffs();
}


void rotatingTranslationPointPatchVectorField::write(Ostream& os) const
{
    os.writeKeyword("rotPoint")
        << rotPoint_ << token::END_STATEMENT << nl;

    os.writeKeyword("rotAxis")
        << rotAxis_ << token::END_STATEMENT << nl;

    os.writeKeyword("transVec")
        << transVec_ << token::END_STATEMENT << nl;

    os.writeKeyword("rotMag")
        << rotMag_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    rotatingTranslationPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
