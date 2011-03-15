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

\*---------------------------------------------------------------------------*/

#include "mapPolyMesh.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Map the patch field
template <class Type>
void topoPatchMapper::mapPatchField
(
    const word& fieldName,
    Field<Type>& pF
) const
{
    // Specify that mapping is conservative
    conservative_ = true;

    // Compute the integral of the source field
    Type intSource = pTraits<Type>::zero;
    Type intTarget = pTraits<Type>::zero;

    if (pF.size())
    {
        intSource = sum(pF * tMapper_.patchAreas(patch_.index()));
    }

    // Map patchField onto itself
    pF.autoMap(*this);

    // Compute the integral of the target field
    const polyPatch& ppI = mpm_.mesh().boundaryMesh()[patch_.index()];

    if (pF.size())
    {
        intTarget = sum(pF * mag(ppI.faceAreas()));
    }

    if (polyMesh::debug)
    {
        int oldP = Info().precision();

        // Compare the global integral
        Pout<< " Field : " << fieldName
            << " Patch : " << ppI.name()
            << " integral errors : " << setprecision(10)
            << " source : " << mag(intSource)
            << " target : " << mag(intTarget)
            << " norm : "
            << (mag(intTarget - intSource) / (mag(intSource) + VSMALL))
            << setprecision(oldP)
            << endl;
    }
}


} // End namespace Foam

// ************************************************************************* //
