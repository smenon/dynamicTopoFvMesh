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

#include "conservativeMeshToMesh.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void conservativeMeshToMesh::interpolateInternalField
(
    Field<Type>& toF,
    const GeometricField<Type, fvPatchField, volMesh>& fromVf
) const
{
    if (fromVf.mesh() != fromMesh())
    {
        FatalErrorIn
        (
            "conservativeMeshToMesh::interpolateInternalField(Field<Type>& toF,"
            "const GeometricField<Type, fvPatchField, volMesh>& fromVf) const"
        )   << "the argument field does not correspond to the right mesh. "
            << "Field size: " << fromVf.size()
            << " mesh size: " << fromMesh().nCells()
            << exit(FatalError);
    }

    if (toF.size() != toMesh().nCells())
    {
        FatalErrorIn
        (
            "conservativeMeshToMesh::interpolateInternalField(Field<Type>& toF,"
            "const GeometricField<Type, fvPatchField, volMesh>& fromVf) const"
        )   << "the argument field does not correspond to the right mesh. "
            << "Field size: " << toF.size()
            << " mesh size: " << toMesh().nCells()
            << exit(FatalError);
    }

    // Evaluate the boundary condition
    const_cast<GeometricField<Type, fvPatchField, volMesh>&>
    (fromVf).boundaryField().evaluate();

    // Get the gradient of the original field.
    typedef typename outerProduct<vector, Type>::type GradCmptType;

    tmp<GeometricField<GradCmptType, fvPatchField, volMesh> > gVf =
    (
        fvc::grad(fromVf)
    );

    // Fetch geometry
    const scalarField& toCellVols = toMesh().cellVolumes();
    const vectorField& fromCellCentres = fromMesh().cellCentres();

    forAll (toF, celli)
    {
        // Initialize to zero
        toF[celli] = pTraits<Type>::zero;

        // Fetch addressing and weights for this cell
        const labelList& addr = addressing_[celli];
        const scalarField& w = weights_[celli];
        const vectorField& x = centres_[celli];

        forAll(addr, cellj)
        {
            vector xCo = fromCellCentres[addr[cellj]];

            // Accumulate volume-weighted Taylor-series interpolate
            toF[celli] +=
            (
                w[cellj] *
                (
                    fromVf[addr[cellj]]
                  + (gVf()[addr[cellj]] & (x[cellj] - xCo))
                )
            );
        }

        // Divide by current volume
        toF[celli] /= toCellVols[celli];
    }
}


template<class Type>
void conservativeMeshToMesh::interpolate
(
    GeometricField<Type, fvPatchField, volMesh>& toVf,
    const GeometricField<Type, fvPatchField, volMesh>& fromVf
) const
{
    interpolateInternalField(toVf, fromVf);
}


template<class Type>
void conservativeMeshToMesh::interpolate
(
    GeometricField<Type, fvPatchField, volMesh>& toVf,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tfromVf
) const
{
    interpolate(toVf, tfromVf());
    tfromVf.clear();
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
conservativeMeshToMesh::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& fromVf
) const
{
    // Create and map the internal-field values
    Field<Type> internalField(toMesh().nCells(), pTraits<Type>::zero);

    interpolateInternalField(internalField, fromVf);

    // Check whether both meshes have got the same number
    // of boundary patches
    if (fromMesh().boundary().size() != toMesh().boundary().size())
    {
        FatalErrorIn
        (
            "conservativeMeshToMesh::interpolate"
            "(const GeometricField<Type, fvPatchField, volMesh>& fromVf) const"
        )   << "Incompatible meshes: different number of boundaries, "
               "only internal field may be interpolated"
            << exit(FatalError);
    }

    // Create and map the patch field values
    PtrList<fvPatchField<Type> > patchFields
    (
        fromVf.boundaryField().size()
    );

    forAll(fromVf.boundaryField(), patchI)
    {
        patchFields.set
        (
            patchI,
            fvPatchField<Type>::New
            (
                fromVf.boundaryField()[patchI],
                toMesh().boundary()[patchI],
                DimensionedField<Type, volMesh>::null(),
                conservativePatchFieldInterpolator()
            )
        );
    }

    // Create the complete field from the pieces
    tmp<GeometricField<Type, fvPatchField, volMesh> > ttoF
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "interpolated(" + fromVf.name() + ')',
                toMesh().time().timeName(),
                toMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            toMesh(),
            fromVf.dimensions(),
            internalField,
            patchFields
        )
    );

    return ttoF;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
conservativeMeshToMesh::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tfromVf
) const
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tint =
    (
        interpolate(tfromVf())
    );

    tfromVf.clear();

    return tint;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
