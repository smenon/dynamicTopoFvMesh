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

//- Interpolate internal field values (conservative first-order)
template<class Type>
void conservativeMeshToMesh::interpolateInternalFieldConserveFirstOrder
(
    Field<Type>& toF,
    const GeometricField<Type, fvPatchField, volMesh>& fromVf
) const
{
    if (fromVf.mesh() != origSrcMesh())
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
            << " mesh size: " << origSrcMesh().nCells()
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
            << " mesh size: " << origTgtMesh().nCells()
            << exit(FatalError);
    }

    // Fetch geometry
    const scalarField& toCellVols = origTgtMesh().cellVolumes();

    forAll (toF, celli)
    {
        // Initialize to zero
        toF[celli] = pTraits<Type>::zero;

        // Fetch addressing and weights for this cell
        const labelList& addr = addressing_[celli];
        const scalarField& w = volumes_[celli];

        // Accumulate volume-weighted interpolate
        forAll(addr, cellj)
        {
            toF[celli] += (w[cellj] * fromVf[addr[cellj]]);
        }

        // Divide by current volume
        toF[celli] /= toCellVols[celli];
    }
}


//- Interpolate internal field values (conservative)
//  with supplied gradient
template<class Type>
void conservativeMeshToMesh::interpolateInternalFieldConserve
(
    Field<Type>& toF,
    const GeometricField
    <
        Type,
        fvPatchField,
        volMesh
    >& fromVf,
    const GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >& fromgVf
) const
{
    if (fromVf.mesh() != origSrcMesh())
    {
        FatalErrorIn
        (
            "\n\n"
            "void conservativeMeshToMesh::interpolateInternalFieldConserve\n"
            "(\n"
            "    Field<Type>& toF,\n"
            "    const GeometricField<Type, fvPatchField, volMesh>& fromVf,\n"
            "    const GeometricField\n"
            "    <\n"
            "        typename outerProduct<vector, Type>::type,\n"
            "        fvPatchField,\n"
            "        volMesh\n"
            "    >& fromgVf\n"
            ") const\n"
        )   << "The argument field does not correspond to the right mesh. "
            << "Field size: " << fromVf.size()
            << " mesh size: " << origSrcMesh().nCells()
            << exit(FatalError);
    }

    if (toF.size() != origTgtMesh().nCells())
    {
        FatalErrorIn
        (
            "\n\n"
            "void conservativeMeshToMesh::interpolateInternalFieldConserve\n"
            "(\n"
            "    Field<Type>& toF,\n"
            "    const GeometricField<Type, fvPatchField, volMesh>& fromVf,\n"
            "    const GeometricField\n"
            "    <\n"
            "        typename outerProduct<vector, Type>::type,\n"
            "        fvPatchField,\n"
            "        volMesh\n"
            "    >& fromgVf\n"
            ") const\n"
        )   << "The argument field does not correspond to the right mesh. "
            << "Field size: " << toF.size()
            << " mesh size: " << origTgtMesh().nCells()
            << exit(FatalError);
    }

    // Fetch geometry
    const scalarField& toCellVols = origTgtMesh().cellVolumes();
    const vectorField& fromCellCentres = origSrcMesh().cellCentres();

    forAll (toF, celli)
    {
        // Initialize to zero
        toF[celli] = pTraits<Type>::zero;

        // Fetch addressing and weights for this cell
        const labelList& addr = addressing_[celli];
        const scalarField& w = volumes_[celli];
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
                  + (fromgVf[addr[cellj]] & (x[cellj] - xCo))
                )
            );
        }

        // Divide by current volume
        toF[celli] /= toCellVols[celli];
    }
}


//- Interpolate internal field values (conservative)
template<class Type>
void conservativeMeshToMesh::interpolateInternalFieldConserve
(
    Field<Type>& toF,
    const GeometricField<Type, fvPatchField, volMesh>& fromVf
) const
{
    // Get the gradient of the original field.
    typedef typename outerProduct<vector, Type>::type GradCmptType;

    tmp<GeometricField<GradCmptType, fvPatchField, volMesh> > gVf =
    (
        fvc::grad(fromVf)
    );

    interpolateInternalFieldConserve(toF, fromVf, gVf());
}


//- Interpolate internal field values (inverse-distance)
template<class Type>
void conservativeMeshToMesh::interpolateInternalFieldInvDist
(
    Field<Type>& toF,
    const GeometricField<Type, fvPatchField, volMesh>& fromVf
) const
{
    if (fromVf.mesh() != origSrcMesh())
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
            << " mesh size: " << origSrcMesh().nCells()
            << exit(FatalError);
    }

    if (toF.size() != origTgtMesh().nCells())
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
            << " mesh size: " << origTgtMesh().nCells()
            << exit(FatalError);
    }

    // Fetch geometry
    const vectorField& newCentres = origTgtMesh().cellCentres();
    const vectorField& oldCentres = origSrcMesh().cellCentres();

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
    Field<Type>& toVf,
    const GeometricField<Type, fvPatchField, volMesh>& fromVf,
    const label method
) const
{
    switch (method)
    {
        case CONSERVATIVE:
        {
            interpolateInternalFieldConserve(toVf, fromVf);

            break;
        }

        case INVERSE_DISTANCE:
        {
            interpolateInternalFieldInvDist
            (
                toVf,
                fromVf
            );

            break;
        }

        case CONSERVATIVE_FIRST_ORDER:
        {
            interpolateInternalFieldConserveFirstOrder(toVf, fromVf);

            break;
        }

        default:
        {
            FatalErrorIn
            (
                "\n\n"
                "void conservativeMeshToMesh::interpolateInternalField\n"
            )   << "unknown interpolation scheme " << method
                << exit(FatalError);
        }
    }
}


template<class Type>
void conservativeMeshToMesh::interpolateInternalField
(
    Field<Type>& toF,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tfromVf,
    const label method
) const
{
    interpolateInternalField(toF, tfromVf(), method);
    tfromVf.clear();
}


//- Interpolate volume field with a supplied gradient
template<class Type>
void conservativeMeshToMesh::interpolate
(
    GeometricField<Type, fvPatchField, volMesh>& toVf,
    const GeometricField<Type, fvPatchField, volMesh>& fromVf,
    const GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >& fromgVf,
    const label method
) const
{
    switch (method)
    {
        case CONSERVATIVE:
        {
            interpolateInternalFieldConserve(toVf, fromVf, fromgVf);

            break;
        }

        case INVERSE_DISTANCE:
        {
            interpolateInternalFieldInvDist(toVf, fromVf);

            break;
        }

        case CONSERVATIVE_FIRST_ORDER:
        {
            interpolateInternalFieldConserveFirstOrder(toVf, fromVf);

            break;
        }

        default:
        {
            FatalErrorIn
            (
                "\n\n"
                "void conservativeMeshToMesh::interpolateInternalField()\n"
            )   << "unknown interpolation scheme " << method
                << exit(FatalError);
        }
    }
}


template<class Type>
void conservativeMeshToMesh::interpolate
(
    GeometricField<Type, fvPatchField, volMesh>& toVf,
    const GeometricField<Type, fvPatchField, volMesh>& fromVf,
    const label method
) const
{
    interpolateInternalField(toVf, fromVf, method);

    // Temporary: Not mapping boundary fields at this point.
    forAll(toVf.boundaryField(), patchI)
    {
        toVf.boundaryField()[patchI] = pTraits<Type>::zero;
    }
}


template<class Type>
void conservativeMeshToMesh::interpolate
(
    GeometricField<Type, fvPatchField, volMesh>& toVf,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tfromVf,
    const label method
) const
{
    interpolate(toVf, tfromVf(), method);
    tfromVf.clear();
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
conservativeMeshToMesh::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& fromVf,
    const label method
) const
{
    // Create and map the internal-field values
    Field<Type> internalField(origTgtMesh().nCells(), pTraits<Type>::zero);

    interpolateInternalField(internalField, fromVf, method);

    // Check whether both meshes have got the same number
    // of boundary patches
    if (origSrcMesh().boundary().size() != origTgtMesh().boundary().size())
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
                origTgtMesh().boundary()[patchI],
                DimensionedField<Type, volMesh>::null(),
                meshToMesh::patchFieldInterpolator(boundaryAddressing_[patchI])
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
                origTgtMesh().time().timeName(),
                origTgtMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            origTgtMesh(),
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
    const label method
) const
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tint =
    (
        interpolate(tfromVf(), method)
    );

    tfromVf.clear();

    return tint;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
