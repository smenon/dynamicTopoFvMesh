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
#include "pointMesh.H"
#include "coupledInfo.H"
#include "emptyPolyPatch.H"
#include "fvPatchFieldMapper.H"
#include "pointPatchFieldMapper.H"

#include "processorPolyPatch.H"
#include "fixedValueFvPatchFields.H"

#include "faceSetAlgorithm.H"
#include "cellSetAlgorithm.H"
#include "pointSetAlgorithm.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Sub-classes * * * * * * * * * * * * * * * //

//- Generic subMesh mapper
template <class MeshType>
class coupledInfo<MeshType>::subMeshPatchMapper
:
    public fvPatchFieldMapper
{
    label sizeBeforeMapping_;

    labelField directAddressing_;

public:

    // Constructors

        //- Construct from components
        inline subMeshPatchMapper
        (
            const label sbm,
            const labelList& da
        )
        :
            sizeBeforeMapping_(sbm),
            directAddressing_(da)
        {}

        //- Construct given addressing
        inline subMeshPatchMapper
        (
            const coupledInfo& cInfo,
            const label patchI
        );

    // Destructor

        virtual ~subMeshPatchMapper()
        {}

    // Member Functions

        label size() const
        {
            return directAddressing_.size();
        }

        label sizeBeforeMapping() const
        {
            return sizeBeforeMapping_;
        }

        bool direct() const
        {
            return true;
        }

        const unallocLabelList& directAddressing() const
        {
            return directAddressing_;
        }

        //- Has unmapped elements
        virtual bool hasUnmapped() const
        {
            return false;
        }
};


//- Generic subMesh point mapper
template <class MeshType>
class coupledInfo<MeshType>::subMeshPointMapper
:
    public pointPatchFieldMapper
{
    label sizeBeforeMapping_;

    labelField directAddressing_;

public:

    // Constructors

        //- Construct from components
        inline subMeshPointMapper
        (
            const label sbm,
            const labelList& da
        )
        :
            sizeBeforeMapping_(sbm),
            directAddressing_(da)
        {}

        //- Construct given addressing
        inline subMeshPointMapper
        (
            const coupledInfo& cInfo,
            const label patchI
        );

    // Destructor

        virtual ~subMeshPointMapper()
        {}

    // Member Functions

        label size() const
        {
            return directAddressing_.size();
        }

        label sizeBeforeMapping() const
        {
            return sizeBeforeMapping_;
        }

        bool direct() const
        {
            return true;
        }

        const unallocLabelList& directAddressing() const
        {
            return directAddressing_;
        }

        //- Has unmapped elements
        virtual bool hasUnmapped() const
        {
            return false;
        }
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct given mesh, coupleMap and master / slave indices
template <class MeshType>
coupledInfo<MeshType>::coupledInfo
(
    const MeshType& mesh,
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


// Construct from components
template <class MeshType>
coupledInfo<MeshType>::coupledInfo
(
    const MeshType& mesh,
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
template <class MeshType>
coupledInfo<MeshType>::subMeshPatchMapper::subMeshPatchMapper
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


//- Construct given addressing
template <class MeshType>
coupledInfo<MeshType>::subMeshPointMapper::subMeshPointMapper
(
    const coupledInfo& cInfo,
    const label patchI
)
:
    sizeBeforeMapping_(cInfo.baseMesh().boundary()[patchI].patch().nPoints())
{
    // Fetch the subMesh / baseMesh patches
    const polyPatch& sPatch = cInfo.subMesh().boundary()[patchI].patch();
    const polyPatch& bPatch = cInfo.baseMesh().boundary()[patchI].patch();

    // Fetch meshPoints for this patch, and the pointMap for the subMesh
    const labelList& sMeshPoints = sPatch.meshPoints();
    const labelList& pointMap = cInfo.map().pointMap();

    // Prepare direct addressing
    directAddressing_.setSize(sPatch.nPoints());

    forAll(sMeshPoints, pointI)
    {
        const label sMeshIndex = sMeshPoints[pointI];
        const label bMeshIndex = pointMap[sMeshIndex];
        const label localIndex = bPatch.whichPoint(bMeshIndex);

        directAddressing_[pointI] = localIndex;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return a const reference to the parent mesh
template <class MeshType>
const MeshType&
coupledInfo<MeshType>::baseMesh() const
{
    return mesh_;
}


// Set a new subMesh
template <class MeshType>
void coupledInfo<MeshType>::setMesh
(
    label index,
    MeshType* mesh
)
{
    subMesh_.set(mesh);
}


// Set the pointAlgorithm
template <class MeshType>
inline void
coupledInfo<MeshType>::setPointAlgorithm(pointSetAlgorithm* algorithm)
{
    pointAlgorithm_.set(algorithm);
}


// Set the faceAlgorithm
template <class MeshType>
inline void
coupledInfo<MeshType>::setFaceAlgorithm(faceSetAlgorithm* algorithm)
{
    faceAlgorithm_.set(algorithm);
}


// Set the cellAlgorithm
template <class MeshType>
inline void
coupledInfo<MeshType>::setCellAlgorithm(cellSetAlgorithm* algorithm)
{
    cellAlgorithm_.set(algorithm);
}


// Return a reference to the subMesh
template <class MeshType>
MeshType& coupledInfo<MeshType>::subMesh()
{
    if (!subMesh_.valid())
    {
        FatalErrorIn("MeshType& coupledInfo::subMesh()")
            << " Sub-mesh pointer has not been set."
            << abort(FatalError);
    }

    return subMesh_();
}


// Return a const reference to the subMesh
template <class MeshType>
const MeshType& coupledInfo<MeshType>::subMesh() const
{
    if (!subMesh_.valid())
    {
        FatalErrorIn("const MeshType& coupledInfo::subMesh() const")
            << " Sub-mesh pointer has not been set."
            << abort(FatalError);
    }

    return subMesh_();
}


// Return if maps have been built
template <class MeshType>
bool coupledInfo<MeshType>::builtMaps() const
{
    return builtMaps_;
}


// Set internal state of maps as built
template <class MeshType>
void coupledInfo<MeshType>::setBuiltMaps()
{
    builtMaps_ = true;
}


// Return a reference to the coupleMap
template <class MeshType>
coupleMap& coupledInfo<MeshType>::map()
{
    return map_;
}


// Return a const reference to the coupleMap
template <class MeshType>
const coupleMap& coupledInfo<MeshType>::map() const
{
    return map_;
}


// Return a const reference to the pointAlgorithm
template <class MeshType>
inline const pointSetAlgorithm&
coupledInfo<MeshType>::pointAlgorithm() const
{
    if (!pointAlgorithm_.valid())
    {
        FatalErrorIn("const pointSetAlgorithm& coupledInfo::pointAlgorithm()")
            << " Algorithm pointer has not been set."
            << abort(FatalError);
    }

    return pointAlgorithm_();
}


// Return a const reference to the faceAlgorithm
template <class MeshType>
inline const faceSetAlgorithm&
coupledInfo<MeshType>::faceAlgorithm() const
{
    if (!faceAlgorithm_.valid())
    {
        FatalErrorIn("const faceSetAlgorithm& coupledInfo::faceAlgorithm()")
            << " Algorithm pointer has not been set."
            << abort(FatalError);
    }

    return faceAlgorithm_();
}


// Return a const reference to the cellAlgorithm
template <class MeshType>
inline const cellSetAlgorithm&
coupledInfo<MeshType>::cellAlgorithm() const
{
    if (!cellAlgorithm_.valid())
    {
        FatalErrorIn("const cellSetAlgorithm& coupledInfo::cellAlgorithm()")
            << " Algorithm pointer has not been set."
            << abort(FatalError);
    }

    return cellAlgorithm_();
}


// Return the master face zone ID
template <class MeshType>
label coupledInfo<MeshType>::masterFaceZone() const
{
    return masterFaceZone_;
}


// Return the slave face zone ID
template <class MeshType>
label coupledInfo<MeshType>::slaveFaceZone() const
{
    return slaveFaceZone_;
}


// Subset point field
template <class MeshType>
template <class GeomField, class ZeroType>
tmp<GeomField>
coupledInfo<MeshType>::subSetPointField
(
    const GeomField& f,
    const ZeroType& zeroValue,
    const labelList& mapper
) const
{
    typedef typename GeomField::InternalField InternalField;
    typedef typename GeomField::PatchFieldType PatchFieldType;
    typedef typename GeomField::GeometricBoundaryField GeomBdyFieldType;
    typedef typename GeomField::DimensionedInternalField DimInternalField;

    // Fetch the reference to the pointMesh
    const pointMesh& pMesh = pointMesh::New(subMesh());

    // Create and map the internal-field values
    InternalField internalField(f.internalField(), mapper);

    // Create and map the patch field values
    label nPatches = pMesh.boundary().size();
    PtrList<PatchFieldType> patchFields(nPatches);

    // Define patch type names, assumed to be
    // common for point / volume / surface fields
    word emptyType(emptyPolyPatch::typeName);
    word processorType(processorPolyPatch::typeName);

    // Create dummy types for initial field creation
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
                    pMesh.boundary()[patchI],
                    DimInternalField::null()
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
                    PatchFieldType::calculatedType(),
                    pMesh.boundary()[patchI],
                    DimInternalField::null()
                )
            );
        }
    }

    // Create new field from pieces
    tmp<GeomField> subFld
    (
        new GeomField
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
            pMesh,
            f.dimensions(),
            internalField,
            patchFields
        )
    );

    // Set correct references for patch internal fields,
    // and map values from the supplied geometric field
    GeomBdyFieldType& bf = subFld().boundaryField();

    forAll(bf, patchI)
    {
        if (patchI == (nPatches - 1))
        {
            // Artificially set last patch
            bf.set
            (
                patchI,
                PatchFieldType::New
                (
                    emptyType,
                    pMesh.boundary()[patchI],
                    subFld().dimensionedInternalField()
                )
            );
        }
        else
        if (isA<processorPolyPatch>(subMesh().boundary()[patchI].patch()))
        {
            bf.set
            (
                patchI,
                PatchFieldType::New
                (
                    processorType,
                    pMesh.boundary()[patchI],
                    subFld().dimensionedInternalField()
                )
            );

            // Avoid dealing with uninitialised values
            // by artificially assigning to zero
            bf[patchI] == zeroValue;
        }
        else
        {
            bf.set
            (
                patchI,
                PatchFieldType::New
                (
                    f.boundaryField()[patchI],
                    pMesh.boundary()[patchI],
                    subFld().dimensionedInternalField(),
                    subMeshPointMapper(*this, patchI)
                )
            );
        }
    }

    return subFld;
}


// Subset geometric field
template <class MeshType>
template <class GeomField, class ZeroType>
tmp<GeomField>
coupledInfo<MeshType>::subSetField
(
    const GeomField& f,
    const ZeroType& zeroValue,
    const labelList& mapper
) const
{
    typedef typename GeomField::InternalField InternalField;
    typedef typename GeomField::PatchFieldType PatchFieldType;
    typedef typename GeomField::GeometricBoundaryField GeomBdyFieldType;
    typedef typename GeomField::DimensionedInternalField DimInternalField;

    // Create and map the internal-field values
    InternalField internalField(f.internalField(), mapper);

    // Create and map the patch field values
    label nPatches = subMesh().boundary().size();
    PtrList<PatchFieldType> patchFields(nPatches);

    // Define patch type names, assumed to be
    // common for point / volume / surface fields
    word emptyType(emptyPolyPatch::typeName);
    word processorType(processorPolyPatch::typeName);

    // Create dummy types for initial field creation
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
                    DimInternalField::null()
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
                    PatchFieldType::calculatedType(),
                    subMesh().boundary()[patchI],
                    DimInternalField::null()
                )
            );
        }
    }

    // Create new field from pieces
    tmp<GeomField> subFld
    (
        new GeomField
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

    // Set correct references for patch internal fields,
    // and map values from the supplied geometric field
    GeomBdyFieldType& bf = subFld().boundaryField();

    forAll(bf, patchI)
    {
        if (patchI == (nPatches - 1))
        {
            // Artificially set last patch
            bf.set
            (
                patchI,
                PatchFieldType::New
                (
                    emptyType,
                    subMesh().boundary()[patchI],
                    subFld().dimensionedInternalField()
                )
            );
        }
        else
        if (isA<processorPolyPatch>(subMesh().boundary()[patchI].patch()))
        {
            bf.set
            (
                patchI,
                PatchFieldType::New
                (
                    processorType,
                    subMesh().boundary()[patchI],
                    subFld().dimensionedInternalField()
                )
            );

            // Avoid dealing with uninitialised values
            // by artificially assigning to zero
            bf[patchI] == zeroValue;
        }
        else
        {
            bf.set
            (
                patchI,
                PatchFieldType::New
                (
                    f.boundaryField()[patchI],
                    subMesh().boundary()[patchI],
                    subFld().dimensionedInternalField(),
                    subMeshPatchMapper(*this, patchI)
                )
            );
        }
    }

    return subFld;
}


// Subset point fields from registry to output stream
template <class MeshType>
template <class GeomField, class ZeroType>
void coupledInfo<MeshType>::sendPointFields
(
    const wordList& fieldNames,
    const word& fieldType,
    const ZeroType& zeroValue,
    const labelList& mapper,
    OSstream& strStream
) const
{
    strStream
        << fieldType << token::NL
        << token::BEGIN_BLOCK << token::NL;

    forAll(fieldNames, i)
    {
        // Fetch object from registry
        const objectRegistry& db = mesh_.thisDb();

        const GeomField& fld = db.lookupObject<GeomField>(fieldNames[i]);

        // Subset the field
        tmp<GeomField> tsubFld = subSetPointField(fld, zeroValue, mapper);

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


// Subset geometric fields from registry to output stream
template <class MeshType>
template <class GeomField, class ZeroType>
void coupledInfo<MeshType>::sendFields
(
    const wordList& fieldNames,
    const word& fieldType,
    const ZeroType& zeroValue,
    const labelList& mapper,
    OSstream& strStream
) const
{
    strStream
        << fieldType << token::NL
        << token::BEGIN_BLOCK << token::NL;

    forAll(fieldNames, i)
    {
        // Fetch object from registry
        const objectRegistry& db = mesh_.thisDb();

        const GeomField& fld = db.lookupObject<GeomField>(fieldNames[i]);

        // Subset the field
        tmp<GeomField> tsubFld = subSetField(fld, zeroValue, mapper);

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


// Set point field pointers from input dictionary
template <class MeshType>
template <class GeomField>
void coupledInfo<MeshType>::setPointField
(
    const wordList& fieldNames,
    const dictionary& fieldDicts,
    const label internalSize,
    PtrList<GeomField>& fields
) const
{
    typedef typename GeomField::InternalField InternalField;
    typedef typename GeomField::PatchFieldType PatchFieldType;
    typedef typename GeomField::GeometricBoundaryField GeomBdyFieldType;
    typedef typename GeomField::DimensionedInternalField DimInternalField;

    // Fetch the reference to the pointMesh
    const pointMesh& pMesh = pointMesh::New(subMesh());

    // Size up the pointer list
    fields.setSize(fieldNames.size());

    // Define patch type names, assumed to be
    // common for point / volume / surface fields
    word emptyType(emptyPolyPatch::typeName);
    word processorType(processorPolyPatch::typeName);

    forAll(fieldNames, i)
    {
        // Create and map the patch field values
        label nPatches = pMesh.boundary().size();

        // Create field parts
        PtrList<PatchFieldType> patchFields(nPatches);

        // Read dimensions
        dimensionSet dimSet
        (
            fieldDicts.subDict(fieldNames[i]).lookup("dimensions")
        );

        // Read the internal field
        InternalField internalField
        (
            "internalField",
            fieldDicts.subDict(fieldNames[i]),
            internalSize
        );

        // Create dummy types for initial field creation
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
                        pMesh.boundary()[patchI],
                        DimInternalField::null()
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
                        PatchFieldType::calculatedType(),
                        pMesh.boundary()[patchI],
                        DimInternalField::null()
                    )
                );
            }
        }

        // Create field with dummy patches
        fields.set
        (
            i,
            new GeomField
            (
                IOobject
                (
                    fieldNames[i],
                    subMesh().time().timeName(),
                    subMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                pMesh,
                dimSet,
                internalField,
                patchFields
            )
        );

        // Set correct references for patch internal fields,
        // and fetch values from the supplied geometric field dictionaries
        GeomBdyFieldType& bf = fields[i].boundaryField();

        forAll(bf, patchI)
        {
            if (patchI == (nPatches - 1))
            {
                // Artificially set last patch
                bf.set
                (
                    patchI,
                    PatchFieldType::New
                    (
                        emptyType,
                        pMesh.boundary()[patchI],
                        fields[i].dimensionedInternalField()
                    )
                );
            }
            else
            if (isA<processorPolyPatch>(subMesh().boundary()[patchI].patch()))
            {
                bf.set
                (
                    patchI,
                    PatchFieldType::New
                    (
                        processorType,
                        pMesh.boundary()[patchI],
                        fields[i].dimensionedInternalField()
                    )
                );
            }
            else
            {
                bf.set
                (
                    patchI,
                    PatchFieldType::New
                    (
                        pMesh.boundary()[patchI],
                        fields[i].dimensionedInternalField(),
                        fieldDicts.subDict
                        (
                            fieldNames[i]
                        ).subDict("boundaryField").subDict
                        (
                            pMesh.boundary()[patchI].name()
                        )
                    )
                );
            }
        }
    }
}


// Set geometric field pointers from input dictionary
template <class MeshType>
template <class GeomField>
void coupledInfo<MeshType>::setField
(
    const wordList& fieldNames,
    const dictionary& fieldDicts,
    const label internalSize,
    PtrList<GeomField>& fields
) const
{
    typedef typename GeomField::InternalField InternalField;
    typedef typename GeomField::PatchFieldType PatchFieldType;
    typedef typename GeomField::GeometricBoundaryField GeomBdyFieldType;
    typedef typename GeomField::DimensionedInternalField DimInternalField;

    // Size up the pointer list
    fields.setSize(fieldNames.size());

    // Define patch type names, assumed to be
    // common for point / volume / surface fields
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
        InternalField internalField
        (
            "internalField",
            fieldDicts.subDict(fieldNames[i]),
            internalSize
        );

        // Create dummy types for initial field creation
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
                        DimInternalField::null()
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
                        PatchFieldType::calculatedType(),
                        subMesh().boundary()[patchI],
                        DimInternalField::null()
                    )
                );
            }
        }

        // Create field with dummy patches
        fields.set
        (
            i,
            new GeomField
            (
                IOobject
                (
                    fieldNames[i],
                    subMesh().time().timeName(),
                    subMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                subMesh(),
                dimSet,
                internalField,
                patchFields
            )
        );

        // Set correct references for patch internal fields,
        // and fetch values from the supplied geometric field dictionaries
        GeomBdyFieldType& bf = fields[i].boundaryField();

        forAll(bf, patchI)
        {
            if (patchI == (nPatches - 1))
            {
                // Artificially set last patch
                bf.set
                (
                    patchI,
                    PatchFieldType::New
                    (
                        emptyType,
                        subMesh().boundary()[patchI],
                        fields[i].dimensionedInternalField()
                    )
                );
            }
            else
            if (isA<processorPolyPatch>(subMesh().boundary()[patchI].patch()))
            {
                bf.set
                (
                    patchI,
                    PatchFieldType::New
                    (
                        processorType,
                        subMesh().boundary()[patchI],
                        fields[i].dimensionedInternalField()
                    )
                );
            }
            else
            {
                bf.set
                (
                    patchI,
                    PatchFieldType::New
                    (
                        subMesh().boundary()[patchI],
                        fields[i].dimensionedInternalField(),
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
    }
}


// Temporarily resize empty patchFields
template <class Type, template<class> class PatchField>
struct EmptyResize
{
    void operator()(const labelList& addr, PatchField<Type>& field)
    {}
};


// Partial template specialization for fvPatchFields
template <class Type>
struct EmptyResize<Type, fvPatchField>
{
    void operator()(const labelList& addr, fvPatchField<Type>& field)
    {
        const word& emptyType = emptyPolyPatch::typeName;

        if (field.empty() && addr.size() && field.type() != emptyType)
        {
            // Artificially set the size prior to remap,
            // since fvPatchField::autoMap appears to be
            // assigning the field to patchInternalField
            // (which is empty, since the patch is zero-sized)
            field.setSize(1);
        }
    }
};


// Resize map for individual field
template <class MeshType>
template <class Type, template<class> class Patch, class Mesh, class Mapper>
void coupledInfo<MeshType>::resizeMap
(
    const label srcIndex,
    const Mapper& internalMapper,
    const List<labelList>& internalReverseMaps,
    const PtrList<Mapper>& boundaryMapper,
    const List<labelListList>& boundaryReverseMaps,
    const List<PtrList<GeometricField<Type, Patch, Mesh> > >& srcFields,
    GeometricField<Type, Patch, Mesh>& field
)
{
    typedef GeometricField<Type, Patch, Mesh> GeomFieldType;
    typedef typename GeomFieldType::PatchFieldType PatchFieldType;

    // autoMap the internal field
    field.internalField().autoMap(internalMapper);

    // Reverse map for additional cells
    forAll(srcFields, pI)
    {
        // Fetch field for this processor
        const GeomFieldType& srcField = srcFields[pI][srcIndex];

        field.internalField().rmap
        (
            srcField.internalField(),
            internalReverseMaps[pI]
        );
    }

    // Map physical boundary-fields
    forAll(boundaryMapper, patchI)
    {
        PatchFieldType& patchField = field.boundaryField()[patchI];

        // Optionally resize empty fields temporarily
        EmptyResize<Type, Patch>()
        (
            boundaryMapper[patchI].directAddressing(),
            patchField
        );

        // autoMap the patchField
        patchField.autoMap(boundaryMapper[patchI]);

        // Reverse map for additional patch faces
        forAll(srcFields, pI)
        {
            // Fetch field for this processor
            const GeomFieldType& srcField = srcFields[pI][srcIndex];

            patchField.rmap
            (
                srcField.boundaryField()[patchI],
                boundaryReverseMaps[pI][patchI]
            );
        }
    }
}


// Resize all fields in registry
template <class MeshType>
template <class GeomField, class Mapper>
void coupledInfo<MeshType>::resizeMap
(
    const wordList& names,
    const objectRegistry& mesh,
    const Mapper& internalMapper,
    const List<labelList>& internalReverseMaps,
    const PtrList<Mapper>& boundaryMapper,
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
        coupledInfo<MeshType>::resizeMap
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
template <class MeshType>
template <class GeomField>
void coupledInfo<MeshType>::resizeBoundaries
(
    const objectRegistry& mesh,
    const fvBoundaryMesh& boundary
)
{
    typedef typename GeomField::PatchFieldType PatchFieldType;
    typedef typename GeomField::GeometricBoundaryField GeomBoundaryType;

    HashTable<const GeomField*> fields(mesh.lookupClass<GeomField>());

    forAllIter(typename HashTable<const GeomField*>, fields, fIter)
    {
        // Fetch field from registry
        GeomField& field = const_cast<GeomField&>(*fIter());

        GeomBoundaryType& bf = field.boundaryField();

        // Resize boundary
        label nPatches = boundary.size();
        label nOldPatches = field.boundaryField().size();

        // Create a new list of boundaries
        PtrList<PatchFieldType> newbf(nPatches);

        // Existing fields are mapped with new fvBoundaryMesh references
        for (label patchI = 0; patchI < nOldPatches; patchI++)
        {
            label oldPatchSize = bf[patchI].size();

            newbf.set
            (
                patchI,
                PatchFieldType::New
                (
                    bf[patchI],
                    boundary[patchI],
                    field,
                    subMeshPatchMapper(oldPatchSize, identity(oldPatchSize))
                )
            );
        }

        // Size up new patches
        for (label patchI = nOldPatches; patchI < nPatches; patchI++)
        {
            newbf.set
            (
                patchI,
                PatchFieldType::New
                (
                    boundary[patchI].type(),
                    boundary[patchI],
                    field
                )
            );
        }

        // Transfer contents with new patches
        bf.transfer(newbf);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// Disallow default bitwise assignment
template <class MeshType>
void coupledInfo<MeshType>::operator=(const coupledInfo& rhs)
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
