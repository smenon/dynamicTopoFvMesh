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
void conservativeMeshToMesh::interpolateInternalFieldConserveFirstOrder
(
    Field<Type>& toF,
    const GeometricField<Type, fvPatchField, volMesh>& fromVf
) const
{
    if (fromVf.mesh() != fromMesh())
    {
        FatalErrorIn
        (
            "\n\n"
            "void conservativeMeshToMesh::"
            "interpolateInternalFieldConserveFirstOrder\n"
            "(\n"
            "    Field<Type>& toF,\n"
            "    const GeometricField<Type, fvPatchField, volMesh>& fromVf\n"
            ") const\n"
        )   << "the argument field does not correspond to the right mesh. "
            << "Field size: " << fromVf.size()
            << " mesh size: " << fromMesh().nCells()
            << exit(FatalError);
    }

    if (toF.size() != toMesh().nCells())
    {
        FatalErrorIn
        (
            "\n\n"
            "void conservativeMeshToMesh::"
            "interpolateInternalFieldConserveFirstOrder\n"
            "(\n"
            "    Field<Type>& toF,\n"
            "    const GeometricField<Type, fvPatchField, volMesh>& fromVf\n"
            ") const\n"
        )   << "the argument field does not correspond to the right mesh. "
            << "Field size: " << toF.size()
            << " mesh size: " << toMesh().nCells()
            << exit(FatalError);
    }

    // Fetch geometry
    const scalarField& toCellVols = toMesh().cellVolumes();

    forAll (toF, celli)
    {
        // Initialize to zero
        toF[celli] = pTraits<Type>::zero;

        // Fetch addressing and weights for this cell
        const labelList& addr = addressing_[celli];
        const scalarField& w = weights_[celli];

        // Accumulate volume-weighted interpolate
        forAll(addr, cellj)
        {
            toF[celli] += (w[cellj] * fromVf[addr[cellj]]);
        }

        // Divide by current volume
        toF[celli] /= toCellVols[celli];
    }
}


template<class Type>
void conservativeMeshToMesh::interpolateInternalFieldConserve
(
    Field<Type>& toF,
    const GeometricField<Type, fvPatchField, volMesh>& fromVf
) const
{
    if (fromVf.mesh() != fromMesh())
    {
        FatalErrorIn
        (
            "\n\n"
            "void conservativeMeshToMesh::interpolateInternalFieldConserve\n"
            "(\n"
            "    Field<Type>& toF,\n"
            "    const GeometricField<Type, fvPatchField, volMesh>& fromVf\n"
            ") const\n"
        )   << "the argument field does not correspond to the right mesh. "
            << "Field size: " << fromVf.size()
            << " mesh size: " << fromMesh().nCells()
            << exit(FatalError);
    }

    if (toF.size() != toMesh().nCells())
    {
        FatalErrorIn
        (
            "\n\n"
            "void conservativeMeshToMesh::interpolateInternalFieldConserve\n"
            "(\n"
            "    Field<Type>& toF,\n"
            "    const GeometricField<Type, fvPatchField, volMesh>& fromVf\n"
            ") const\n"
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


//- Interpolate internal field values (inverse-distance)
template<class Type>
void conservativeMeshToMesh::interpolateInternalFieldInvDist
(
    Field<Type>& toF,
    const GeometricField<Type, fvPatchField, volMesh>& fromVf
) const
{
    if (fromVf.mesh() != fromMesh())
    {
        FatalErrorIn
        (
            "\n\n"
            "void conservativeMeshToMesh::interpolateInternalFieldInvDist\n"
            "(\n"
            "    Field<Type>& toF,\n"
            "    const GeometricField<Type, fvPatchField, volMesh>& fromVf\n"
            ") const\n"
        )   << "the argument field does not correspond to the right mesh. "
            << "Field size: " << fromVf.size()
            << " mesh size: " << fromMesh().nCells()
            << exit(FatalError);
    }

    if (toF.size() != toMesh().nCells())
    {
        FatalErrorIn
        (
            "\n\n"
            "void conservativeMeshToMesh::interpolateInternalFieldInvDist\n"
            "(\n"
            "    Field<Type>& toF,\n"
            "    const GeometricField<Type, fvPatchField, volMesh>& fromVf\n"
            ") const\n"
        )   << "the argument field does not correspond to the right mesh. "
            << "Field size: " << toF.size()
            << " mesh size: " << toMesh().nCells()
            << exit(FatalError);
    }

    // Fetch geometry
    const vectorField& newCentres = toMesh().cellCentres();
    const vectorField& oldCentres = fromMesh().cellCentres();

    forAll (toF, celli)
    {
        // Initialize to zero
        toF[celli] = pTraits<Type>::zero;

        scalar weight = 0.0, totalWeight = 0.0;
        const labelList& addr = addressing_[celli];

        forAll(addr, oldCellI)
        {
            weight =
            (
                1.0/stabilise
                (
                    magSqr
                    (
                        newCentres[celli]
                      - oldCentres[addr[oldCellI]]
                    ),
                    VSMALL
                )
            );

            // Accumulate field value
            toF[celli] += (fromVf[addr[oldCellI]] * weight);

            // Accumulate weights
            totalWeight += weight;
        }

        toF[celli] *= (1.0 / totalWeight);
    }
}


template<class Type>
void conservativeMeshToMesh::interpolateInternalField
(
    Field<Type>& toF,
    const GeometricField<Type, fvPatchField, volMesh>& fromVf,
    const label meth
) const
{
    switch (meth)
    {
        case CONSERVATIVE:
        {
            interpolateInternalFieldConserve
            (
                toF,
                fromVf
            );

            break;
        }

        case INVERSE_DISTANCE:
        {
            interpolateInternalFieldInvDist
            (
                toF,
                fromVf
            );

            break;
        }

        case CONSERVATIVE_FIRST_ORDER:
        {
            interpolateInternalFieldConserveFirstOrder
            (
                toF,
                fromVf
            );

            break;
        }

        default:
        {
            FatalErrorIn
            (
                "\n\n"
                "void conservativeMeshToMesh::interpolateInternalField\n"
            )   << "unknown interpolation scheme " << meth
                << exit(FatalError);
        }
    }
}


template<class Type>
void conservativeMeshToMesh::interpolateInternalField
(
    Field<Type>& toF,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tfromVf,
    const label meth
) const
{
    interpolateInternalField(toF, tfromVf(), meth);
    tfromVf.clear();
}


template<class Type>
void conservativeMeshToMesh::interpolate
(
    GeometricField<Type, fvPatchField, volMesh>& toVf,
    const GeometricField<Type, fvPatchField, volMesh>& fromVf,
    const label meth
) const
{
    interpolateInternalField(toVf, fromVf, meth);
}


template<class Type>
void conservativeMeshToMesh::interpolate
(
    GeometricField<Type, fvPatchField, volMesh>& toVf,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tfromVf,
    const label meth
) const
{
    interpolate(toVf, tfromVf(), meth);
    tfromVf.clear();
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
conservativeMeshToMesh::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& fromVf,
    const label meth
) const
{
    // Create and map the internal-field values
    Field<Type> internalField(toMesh().nCells(), pTraits<Type>::zero);

    interpolateInternalField(internalField, fromVf, meth);

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
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tfromVf,
    const label meth
) const
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tint =
    (
        interpolate(tfromVf(), meth)
    );

    tfromVf.clear();

    return tint;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
