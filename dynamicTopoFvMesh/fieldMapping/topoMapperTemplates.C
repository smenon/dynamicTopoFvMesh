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

#include "fvc.H"
#include "leastSquaresGrad.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Store old-time information for all fields of the specified type
template <class GeomField>
void topoMapper::storeOldTimes() const
{
    typedef HashTable<const GeomField*> GeomFieldTable;

    // Fetch all fields from registry
    GeomFieldTable fields(mesh_.objectRegistry::lookupClass<GeomField>());

    // Store old-times before gradient computation
    for
    (
        typename GeomFieldTable::iterator fIter = fields.begin();
        fIter != fields.end();
        ++fIter
    )
    {
        fIter()->storeOldTimes();
    }
}


// Store gradients of fields on the mesh prior to topology changes
template <class Type, class gradType>
void topoMapper::storeGradients
(
    GradientTable& gradTable,
    PtrList<gradType>& gradList
) const
{
    // Define a few typedefs for convenience
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<Type, fvPatchField, volMesh> volType;
    typedef const GeometricField<Type, fvPatchField, volMesh> constVolType;

    typedef HashTable<constVolType*> volTypeTable;

    // Fetch all fields from registry
    volTypeTable fields(mesh_.objectRegistry::lookupClass<volType>());

    // Track field count
    label nFields = 0;

    // Store old-times before gradient computation
    for
    (
        typename volTypeTable::iterator fIter = fields.begin();
        fIter != fields.end();
        ++fIter
    )
    {
        fIter()->storeOldTimes();
        nFields++;
    }

    // Size up the list
    gradList.setSize(nFields);

    label fieldIndex = 0;

    for
    (
        typename volTypeTable::const_iterator fIter = fields.begin();
        fIter != fields.end();
        ++fIter
    )
    {
        const volType& field = *fIter();

        // Compute the gradient.

        // If the fvSolution dictionary contains an entry,
        // use that, otherwise, default to leastSquares
        word gradName("grad(" + field.name() + ')');

        // Register field under a name that's unique
        word registerName("remapGradient(" + field.name() + ')');

        // Make a new entry
        if (disableGradients_)
        {
            gradList.set
            (
                fieldIndex,
                new gradType
                (
                    IOobject
                    (
                        registerName,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        true
                    ),
                    mesh_,
                    dimensioned<GradType>
                    (
                        "0",
                        field.dimensions() / dimLength,
                        pTraits<GradType>::zero
                    )
                )
            );
        }
        else
        if (mesh_.schemesDict().subDict("gradSchemes").found(gradName))
        {
            gradList.set
            (
                fieldIndex,
                new gradType
                (
                    IOobject
                    (
                        registerName,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        true
                    ),
                    fvc::grad(field, gradName)()
                )
            );
        }
        else
        {
            gradList.set
            (
                fieldIndex,
                new gradType
                (
                    IOobject
                    (
                        registerName,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        true
                    ),
                    fv::leastSquaresGrad<Type>(mesh_).grad(field)()
                )
            );
        }

        // Add a map entry
        gradTable.insert
        (
            field.name(),
            GradientMap(registerName, fieldIndex++)
        );
    }
}


} // End namespace Foam

// ************************************************************************* //
