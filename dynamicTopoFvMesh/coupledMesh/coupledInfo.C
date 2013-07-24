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

#include "Time.H"
#include "coupledInfo.H"
#include "dynamicTopoFvMesh.H"
#include "emptyFvPatchFields.H"
#include "emptyFvsPatchFields.H"
#include "processorFvPatchFields.H"
#include "processorFvsPatchFields.H"
#include "fixedValueFvPatchFields.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Constructor for coupledInfo
coupledInfo::coupledInfo
(
    const dynamicTopoFvMesh& mesh,
    const coupleMap& cMap,
    const label mfzIndex,
    const label sfzIndex
)
:
    mesh_(mesh),
    builtMaps_(false),
    map_(cMap),
    masterFaceZone_(mfzIndex),
    slaveFaceZone_(sfzIndex)
{}


coupledInfo::coupledInfo
(
    const dynamicTopoFvMesh& mesh,
    const bool isTwoDMesh,
    const bool isLocal,
    const bool isSend,
    const label patchIndex,
    const label mPatch,
    const label sPatch,
    const label mfzIndex,
    const label sfzIndex
)
:
    mesh_(mesh),
    builtMaps_(false),
    map_
    (
        IOobject
        (
            "coupleMap_"
          + Foam::name(mPatch)
          + "_To_"
          + Foam::name(sPatch)
          + word(isLocal ? "_Local" : "_Proc")
          + word(isSend ? "_Send" : "_Recv"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        isTwoDMesh,
        isLocal,
        isSend,
        patchIndex,
        mPatch,
        sPatch
    ),
    masterFaceZone_(mfzIndex),
    slaveFaceZone_(sfzIndex)
{}


//- Construct given addressing
coupledInfo::subMeshMapper::subMeshMapper
(
    const coupledInfo& cInfo,
    const label patchI
)
:
    sizeBeforeMapping_(cInfo.baseMesh().boundary()[patchI].size()),
    directAddressing_
    (
        SubList<label>
        (
            cInfo.map().faceMap(),
            cInfo.subMesh().boundary()[patchI].size(),
            cInfo.subMesh().boundary()[patchI].patch().start()
        )
    )
{
    // Offset indices
    label pStart = cInfo.baseMesh().boundary()[patchI].patch().start();

    forAll(directAddressing_, faceI)
    {
        directAddressing_[faceI] -= pStart;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const dynamicTopoFvMesh& coupledInfo::baseMesh() const
{
    return mesh_;
}


void coupledInfo::setMesh
(
    label index,
    dynamicTopoFvMesh* mesh
)
{
    subMesh_.set(mesh);
}


dynamicTopoFvMesh& coupledInfo::subMesh()
{
    if (!subMesh_.valid())
    {
        FatalErrorIn("dynamicTopoFvMesh& coupledInfo::subMesh()")
            << " Sub-mesh pointer has not been set."
            << abort(FatalError);
    }

    return subMesh_();
}


const dynamicTopoFvMesh& coupledInfo::subMesh() const
{
    if (!subMesh_.valid())
    {
        FatalErrorIn("const dynamicTopoFvMesh& coupledInfo::subMesh() const")
            << " Sub-mesh pointer has not been set."
            << abort(FatalError);
    }

    return subMesh_();
}


bool coupledInfo::builtMaps() const
{
    return builtMaps_;
}


void coupledInfo::setBuiltMaps()
{
    builtMaps_ = true;
}


coupleMap& coupledInfo::map()
{
    return map_;
}


const coupleMap& coupledInfo::map() const
{
    return map_;
}


label coupledInfo::masterFaceZone() const
{
    return masterFaceZone_;
}


label coupledInfo::slaveFaceZone() const
{
    return slaveFaceZone_;
}


// Set subMesh centres
void coupledInfo::setCentres(PtrList<volVectorField>& centres) const
{
    // Fetch reference to subMesh
    const dynamicTopoFvMesh& mesh = subMesh();

    // Set size
    centres.setSize(1);

    vectorField Cv(mesh.cellCentres());
    vectorField Cf(mesh.faceCentres());

    // Create and map the patch field values
    label nPatches = mesh.boundary().size();

    // Create field parts
    PtrList<fvPatchField<vector> > volCentrePatches(nPatches);

    // Over-ride and set all patches to fixedValue
    for (label patchI = 0; patchI < nPatches; patchI++)
    {
        volCentrePatches.set
        (
            patchI,
            new fixedValueFvPatchField<vector>
            (
                mesh.boundary()[patchI],
                DimensionedField<vector, volMesh>::null()
            )
        );

        // Slice field to patch (forced assignment)
        volCentrePatches[patchI] ==
        (
            mesh.boundaryMesh()[patchI].patchSlice(Cf)
        );
    }

    // Set the cell-centres pointer.
    centres.set
    (
        0,
        new volVectorField
        (
            IOobject
            (
                "cellCentres",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimLength,
            SubField<vector>(Cv, mesh.nCells()),
            volCentrePatches
        )
    );
}


// Subset geometric field
template<class Type, template<class> class PatchField, class Mesh>
tmp<GeometricField<Type, PatchField, Mesh> >
coupledInfo::subSetField(const GeometricField<Type, PatchField, Mesh>& f) const
{
    typedef PatchField<Type> PatchFieldType;
    typedef GeometricField<Type, PatchField, Mesh> GeomFieldType;

    // Create and map the internal-field values
    Field<Type> internalField(f.internalField(), map().cellMap());

    // Create and map the patch field values
    label nPatches = subMesh().boundary().size();
    PtrList<PatchFieldType> patchFields(nPatches);

    // Define patch type names, assumed to be
    // common for volume and surface fields
    word emptyType(emptyPolyPatch::typeName);
    word processorType(processorPolyPatch::typeName);

    forAll(patchFields, patchI)
    {
        if (patchI == (nPatches - 1))
        {
            // Artificially set last patch
            patchFields.set
            (
                patchI,
                PatchFieldType::New
                (
                    emptyType,
                    subMesh().boundary()[patchI],
                    DimensionedField<Type, Mesh>::null()
                )
            );
        }
        else
        if (isA<processorPolyPatch>(subMesh().boundary()[patchI].patch()))
        {
            patchFields.set
            (
                patchI,
                PatchFieldType::New
                (
                    processorType,
                    subMesh().boundary()[patchI],
                    DimensionedField<Type, Mesh>::null()
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchI,
                PatchFieldType::New
                (
                    f.boundaryField()[patchI],
                    subMesh().boundary()[patchI],
                    DimensionedField<Type, Mesh>::null(),
                    subMeshMapper(*this, patchI)
                )
            );
        }
    }

    // Create new field from pieces
    tmp<GeomFieldType> subFld
    (
        new GeomFieldType
        (
            IOobject
            (
                "subField_" + f.name(),
                subMesh().time().timeName(),
                subMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            subMesh(),
            f.dimensions(),
            internalField,
            patchFields
        )
    );

    return subFld;
}


// Subset geometric fields from registry to output stream
template<class GeomField>
void coupledInfo::send
(
    const wordList& fieldNames,
    const word& fieldType,
    OSstream& strStream
) const
{
    strStream
        << fieldType << token::NL
        << token::BEGIN_BLOCK << token::NL;

    forAll(fieldNames, i)
    {
        // Fetch object from registry
        const GeomField& fld = mesh_.lookupObject<GeomField>(fieldNames[i]);

        // Subset the field
        tmp<GeomField> tsubFld = subSetField(fld);

        // Send field subset through stream
        strStream
            << fieldNames[i]
            << token::NL << token::BEGIN_BLOCK
            << tsubFld
            << token::NL << token::END_BLOCK
            << token::NL;
    }

    strStream
        << token::END_BLOCK << token::NL;
}


// Set geometric field pointer from input dictionary
template<class Type, template<class> class PatchField, class Mesh>
void coupledInfo::setField
(
    const wordList& fieldNames,
    const dictionary& fieldDicts,
    PtrList<GeometricField<Type, PatchField, Mesh> >& fields
) const
{
    typedef PatchField<Type> PatchFieldType;
    typedef GeometricField<Type, PatchField, Mesh> GeomFieldType;

    // Size up the pointer list
    fields.setSize(fieldNames.size());

    // Define patch type names, assumed to be
    // common for volume and surface fields
    word emptyType(emptyPolyPatch::typeName);
    word processorType(processorPolyPatch::typeName);

    forAll(fieldNames, i)
    {
        // Create and map the patch field values
        label nPatches = subMesh().boundary().size();

        // Create field parts
        PtrList<PatchFieldType> patchFields(nPatches);

        // Read dimensions
        dimensionSet dimSet
        (
            fieldDicts.subDict(fieldNames[i]).lookup("dimensions")
        );

        // Read the internal field
        Field<Type> internalField
        (
            "internalField",
            fieldDicts.subDict(fieldNames[i]),
            Mesh::size(subMesh())
        );

        // Create a temporary DimensionedField for patch-evaluation
        IOobject io
        (
            fieldNames[i],
            subMesh().time().timeName(),
            subMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        DimensionedField<Type, Mesh> tmpInternal
        (
            io,
            subMesh(),
            dimSet,
            internalField
        );

        forAll(patchFields, patchI)
        {
            if (patchI == (nPatches - 1))
            {
                // Artificially set last patch
                patchFields.set
                (
                    patchI,
                    PatchFieldType::New
                    (
                        emptyType,
                        subMesh().boundary()[patchI],
                        tmpInternal
                    )
                );
            }
            else
            if (isA<processorPolyPatch>(subMesh().boundary()[patchI].patch()))
            {
                patchFields.set
                (
                    patchI,
                    PatchFieldType::New
                    (
                        processorType,
                        subMesh().boundary()[patchI],
                        tmpInternal
                    )
                );
            }
            else
            {
                patchFields.set
                (
                    patchI,
                    PatchFieldType::New
                    (
                        subMesh().boundary()[patchI],
                        tmpInternal,
                        fieldDicts.subDict
                        (
                            fieldNames[i]
                        ).subDict("boundaryField").subDict
                        (
                            subMesh().boundary()[patchI].name()
                        )
                    )
                );
            }
        }

        fields.set
        (
            i,
            new GeomFieldType
            (
                io,
                subMesh(),
                dimSet,
                internalField,
                patchFields
            )
        );
    }
}


template <class GeomField>
void coupledInfo::resizeMap
(
    const label srcIndex,
    const subMeshMapper& internalMapper,
    const List<labelList>& internalReverseMaps,
    const PtrList<subMeshMapper>& boundaryMapper,
    const List<labelListList>& boundaryReverseMaps,
    const List<PtrList<GeomField> >& srcFields,
    GeomField& field
)
{
    // autoMap the internal field
    field.internalField().autoMap(internalMapper);

    // Reverse map for additional cells
    forAll(srcFields, pI)
    {
        // Fetch field for this processor
        const GeomField& srcField = srcFields[pI][srcIndex];

        field.internalField().rmap
        (
            srcField.internalField(),
            internalReverseMaps[pI]
        );
    }

    // Map physical boundary-fields
    forAll(boundaryMapper, patchI)
    {
        // autoMap the patchField
        field.boundaryField()[patchI].autoMap(boundaryMapper[patchI]);

        // Reverse map for additional patch faces
        forAll(srcFields, pI)
        {
            // Fetch field for this processor
            const GeomField& srcField = srcFields[pI][srcIndex];

            field.boundaryField()[patchI].rmap
            (
                srcField.boundaryField()[patchI],
                boundaryReverseMaps[pI][patchI]
            );
        }
    }
}


// Resize all fields in registry
template <class GeomField>
void coupledInfo::resizeMap
(
    const wordList& names,
    const objectRegistry& mesh,
    const subMeshMapper& internalMapper,
    const List<labelList>& internalReverseMaps,
    const PtrList<subMeshMapper>& boundaryMapper,
    const List<labelListList>& boundaryReverseMaps,
    const List<PtrList<GeomField> >& srcFields
)
{
    forAll(names, indexI)
    {
        // Fetch field from registry
        GeomField& field =
        (
            const_cast<GeomField&>
            (
                mesh.lookupObject<GeomField>(names[indexI])
            )
        );

        // Map the field
        coupledInfo::resizeMap
        (
            indexI,
            internalMapper,
            internalReverseMaps,
            boundaryMapper,
            boundaryReverseMaps,
            srcFields,
            field
        );
    }
}


// Resize boundaryFields for all fields in the registry
template<class Type, template<class> class PatchField, class Mesh>
void coupledInfo::resizeBoundaries
(
    const wordList& names,
    const label nOldPatches,
    const objectRegistry& mesh,
    const fvBoundaryMesh& boundary,
    const List<PtrList<GeometricField<Type, PatchField, Mesh> > >& srcFields
)
{
    typedef GeometricField<Type, PatchField, Mesh> GeomFieldType;
    typedef typename GeomFieldType::GeometricBoundaryField GeomBoundaryType;

    forAll(names, indexI)
    {
        // Fetch field from registry
        GeomFieldType& field =
        (
            const_cast<GeomFieldType&>
            (
                mesh.lookupObject<GeomFieldType>(names[indexI])
            )
        );

        GeomBoundaryType& bf = field.boundaryField();

        // Resize boundary
        label nPatches = boundary.size();

        // Existing fields are simply cloned here
        bf.setSize(nPatches);

        // Size up new patches
        for (label patchI = nOldPatches; patchI < nPatches; patchI++)
        {
            bf.set
            (
                patchI,
                PatchField<Type>::New
                (
                    boundary[patchI].type(),
                    boundary[patchI],
                    field
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void coupledInfo::operator=(const coupledInfo& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "void coupledInfo::operator=(const Foam::coupledInfo&)"
        )
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
