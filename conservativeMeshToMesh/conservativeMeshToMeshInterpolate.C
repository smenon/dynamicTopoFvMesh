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
    Field<Type>& tgtField,
    const GeometricField<Type, fvPatchField, volMesh>& srcVf
) const
{
    if (srcVf.mesh() != srcMesh())
    {
        FatalErrorIn
        (
            "\n\n"
            "void conservativeMeshToMesh::"
            "interpolateInternalFieldConserveFirstOrder\n"
            "(\n"
            "    Field<Type>& tgtField,\n"
            "    const GeometricField<Type, fvPatchField, volMesh>& srcVf\n"
            ") const\n"
        )   << " The argument field does not correspond to the right mesh.\n"
            << " Field size: " << srcVf.size() << nl
            << " Mesh size: " << srcMesh().nCells() << nl
            << exit(FatalError);
    }

    if (tgtField.size() != tgtMesh().nCells())
    {
        FatalErrorIn
        (
            "\n\n"
            "void conservativeMeshToMesh::"
            "interpolateInternalFieldConserveFirstOrder\n"
            "(\n"
            "    Field<Type>& tgtField,\n"
            "    const GeometricField<Type, fvPatchField, volMesh>& srcVf\n"
            ") const\n"
        )   << " The argument field does not correspond to the right mesh.\n"
            << " Field size: " << tgtField.size() << nl
            << " mesh size: " << tgtMesh().nCells() << nl
            << exit(FatalError);
    }

    forAll(tgtField, cellI)
    {
        // Initialize to zero
        tgtField[cellI] = pTraits<Type>::zero;

        // Fetch addressing and weights for this cell
        const labelList& addr = addressing_[cellI];
        const scalarField& w = weights_[cellI];

        // Accumulate volume-weighted interpolate
        forAll(addr, cellJ)
        {
            tgtField[cellI] += (w[cellJ] * srcVf[addr[cellJ]]);
        }
    }
}


//- Interpolate internal field values (conservative)
//  with supplied gradient
template<class Type>
void conservativeMeshToMesh::interpolateInternalFieldConserve
(
    Field<Type>& tgtField,
    const GeometricField
    <
        Type,
        fvPatchField,
        volMesh
    >& srcVf,
    const GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >& srcVfGrad
) const
{
    if (srcVf.mesh() != srcMesh())
    {
        FatalErrorIn
        (
            "\n\n"
            "void conservativeMeshToMesh::interpolateInternalFieldConserve\n"
            "(\n"
            "    Field<Type>& tgtField,\n"
            "    const GeometricField<Type, fvPatchField, volMesh>& srcVf,\n"
            "    const GeometricField\n"
            "    <\n"
            "        typename outerProduct<vector, Type>::type,\n"
            "        fvPatchField,\n"
            "        volMesh\n"
            "    >& srcVfGrad\n"
            ") const\n"
        )   << " The argument field does not correspond to the right mesh.\n"
            << " Field size: " << srcVf.size() << nl
            << " mesh size: " << srcMesh().nCells() << nl
            << exit(FatalError);
    }

    if (tgtField.size() != tgtMesh().nCells())
    {
        FatalErrorIn
        (
            "\n\n"
            "void conservativeMeshToMesh::interpolateInternalFieldConserve\n"
            "(\n"
            "    Field<Type>& tgtField,\n"
            "    const GeometricField<Type, fvPatchField, volMesh>& srcVf,\n"
            "    const GeometricField\n"
            "    <\n"
            "        typename outerProduct<vector, Type>::type,\n"
            "        fvPatchField,\n"
            "        volMesh\n"
            "    >& srcVfGrad\n"
            ") const\n"
        )   << " The argument field does not correspond to the right mesh.\n"
            << " Field size: " << tgtField.size() << nl
            << " mesh size: " << tgtMesh().nCells() << nl
            << exit(FatalError);
    }

    // Fetch geometry
    const vectorField& srcCellCentres = srcMesh().cellCentres();

    forAll(tgtField, cellI)
    {
        // Initialize to zero
        tgtField[cellI] = pTraits<Type>::zero;

        // Fetch addressing and weights for this cell
        const labelList& addr = addressing_[cellI];
        const scalarField& w = weights_[cellI];
        const vectorField& x = centres_[cellI];

        forAll(addr, cellJ)
        {
            vector xCo = srcCellCentres[addr[cellJ]];

            // Accumulate volume-weighted Taylor-series interpolate
            tgtField[cellI] +=
            (
                w[cellJ] *
                (
                    srcVf[addr[cellJ]]
                  + (srcVfGrad[addr[cellJ]] & (x[cellJ] - xCo))
                )
            );
        }
    }
}


//- Interpolate internal field values (conservative)
template<class Type>
void conservativeMeshToMesh::interpolateInternalFieldConserve
(
    Field<Type>& tgtField,
    const GeometricField<Type, fvPatchField, volMesh>& srcVf
) const
{
    // Get the gradient of the original field.
    typedef typename outerProduct<vector, Type>::type GradCmptType;

    tmp<GeometricField<GradCmptType, fvPatchField, volMesh> > srcVfGrad =
    (
        fvc::grad(srcVf)
    );

    interpolateInternalFieldConserve(tgtField, srcVf, srcVfGrad());
}


//- Interpolate internal field values (inverse-distance)
template<class Type>
void conservativeMeshToMesh::interpolateInternalFieldInvDist
(
    Field<Type>& tgtField,
    const GeometricField<Type, fvPatchField, volMesh>& srcVf
) const
{
    if (srcVf.mesh() != srcMesh())
    {
        FatalErrorIn
        (
            "\n\n"
            "void conservativeMeshToMesh::interpolateInternalFieldInvDist\n"
            "(\n"
            "    Field<Type>& tgtField,\n"
            "    const GeometricField<Type, fvPatchField, volMesh>& srcVf\n"
            ") const\n"
        )   << " The argument field does not correspond to the right mesh.\n"
            << " Field size: " << srcVf.size() << nl
            << " mesh size: " << srcMesh().nCells() << nl
            << exit(FatalError);
    }

    if (tgtField.size() != tgtMesh().nCells())
    {
        FatalErrorIn
        (
            "\n\n"
            "void conservativeMeshToMesh::interpolateInternalFieldInvDist\n"
            "(\n"
            "    Field<Type>& tgtField,\n"
            "    const GeometricField<Type, fvPatchField, volMesh>& srcVf\n"
            ") const\n"
        )   << " The argument field does not correspond to the right mesh.\n"
            << " Field size: " << tgtField.size() << nl
            << " mesh size: " << tgtMesh().nCells() << nl
            << exit(FatalError);
    }

    // Fetch geometry
    const vectorField& newCentres = tgtMesh().cellCentres();
    const vectorField& oldCentres = srcMesh().cellCentres();

    forAll(tgtField, cellI)
    {
        // Initialize to zero
        tgtField[cellI] = pTraits<Type>::zero;

        scalar weight = 0.0, totalWeight = 0.0;
        const labelList& addr = addressing_[cellI];

        forAll(addr, oldCellI)
        {
            weight =
            (
                1.0/stabilise
                (
                    magSqr
                    (
                        newCentres[cellI]
                      - oldCentres[addr[oldCellI]]
                    ),
                    VSMALL
                )
            );

            // Accumulate field value
            tgtField[cellI] += (srcVf[addr[oldCellI]] * weight);

            // Accumulate weights
            totalWeight += weight;
        }

        tgtField[cellI] *= (1.0 / totalWeight);
    }
}


template<class Type>
void conservativeMeshToMesh::interpolateInternalField
(
    Field<Type>& tgtVf,
    const GeometricField<Type, fvPatchField, volMesh>& srcVf,
    const label method
) const
{
    switch (method)
    {
        case CONSERVATIVE:
        {
            interpolateInternalFieldConserve(tgtVf, srcVf);

            break;
        }

        case INVERSE_DISTANCE:
        {
            interpolateInternalFieldInvDist(tgtVf, srcVf);

            break;
        }

        case CONSERVATIVE_FIRST_ORDER:
        {
            interpolateInternalFieldConserveFirstOrder(tgtVf, srcVf);

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

            break;
        }
    }
}


template<class Type>
void conservativeMeshToMesh::interpolateInternalField
(
    Field<Type>& tgtField,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsrcVf,
    const label method
) const
{
    interpolateInternalField(tgtField, tsrcVf(), method);
    tsrcVf.clear();
}


//- Interpolate volume field with a supplied gradient
template<class Type>
void conservativeMeshToMesh::interpolate
(
    GeometricField<Type, fvPatchField, volMesh>& tgtVf,
    const GeometricField<Type, fvPatchField, volMesh>& srcVf,
    const GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >& srcVfGrad,
    const label method
) const
{
    switch (method)
    {
        case CONSERVATIVE:
        {
            interpolateInternalFieldConserve(tgtVf, srcVf, srcVfGrad);

            break;
        }

        case INVERSE_DISTANCE:
        {
            interpolateInternalFieldInvDist(tgtVf, srcVf);

            break;
        }

        case CONSERVATIVE_FIRST_ORDER:
        {
            interpolateInternalFieldConserveFirstOrder(tgtVf, srcVf);

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

            break;
        }
    }
}


template<class Type>
void conservativeMeshToMesh::interpolate
(
    GeometricField<Type, fvPatchField, volMesh>& tgtVf,
    const GeometricField<Type, fvPatchField, volMesh>& srcVf,
    const label method
) const
{
    interpolateInternalField(tgtVf, srcVf, method);

    // Temporary: Not mapping boundary fields at this point.
    forAll(tgtVf.boundaryField(), patchI)
    {
        tgtVf.boundaryField()[patchI] = pTraits<Type>::zero;
    }
}


template<class Type>
void conservativeMeshToMesh::interpolate
(
    GeometricField<Type, fvPatchField, volMesh>& tgtVf,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsrcVf,
    const label method
) const
{
    interpolate(tgtVf, tsrcVf(), method);
    tsrcVf.clear();
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
conservativeMeshToMesh::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& srcVf,
    const label method
) const
{
    // Create and map the internal-field values
    Field<Type> internalField(tgtMesh().nCells(), pTraits<Type>::zero);

    interpolateInternalField(internalField, srcVf, method);

    // Check whether both meshes have got the same number
    // of boundary patches
    if (srcMesh().boundary().size() != tgtMesh().boundary().size())
    {
        FatalErrorIn
        (
            "conservativeMeshToMesh::interpolate"
            "(const GeometricField<Type, fvPatchField, volMesh>& srcVf) const"
        )   << "Incompatible meshes: different number of boundaries, "
               "only internal field may be interpolated"
            << exit(FatalError);
    }

    // Create and map the patch field values
    PtrList<fvPatchField<Type> > patchFields
    (
        srcVf.boundaryField().size()
    );

    forAll(srcVf.boundaryField(), patchI)
    {
        patchFields.set
        (
            patchI,
            fvPatchField<Type>::New
            (
                srcVf.boundaryField()[patchI],
                tgtMesh().boundary()[patchI],
                DimensionedField<Type, volMesh>::null(),
                patchFieldInterpolator(boundaryAddressing_[patchI])
            )
        );
    }

    // Create the complete field from the pieces
    tmp<GeometricField<Type, fvPatchField, volMesh> > ttgtF
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "interpolated(" + srcVf.name() + ')',
                tgtMesh().time().timeName(),
                tgtMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tgtMesh(),
            srcVf.dimensions(),
            internalField,
            patchFields
        )
    );

    return ttgtF;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
conservativeMeshToMesh::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsrcVf,
    const label method
) const
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tint =
    (
        interpolate(tsrcVf(), method)
    );

    tsrcVf.clear();

    return tint;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
