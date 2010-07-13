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
    topoPatchMapper

Description
    Implementation of the topoPatchMapper class

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*----------------------------------------------------------------------------*/

#include "topoMapper.H"
#include "mapPolyMesh.H"
#include "topoPatchMapper.H"
#include "demandDrivenData.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Clear out local storage
void topoPatchMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
    deleteDemandDrivenData(interpolationAddrPtr_);
    deleteDemandDrivenData(weightsPtr_);
    deleteDemandDrivenData(insertedFaceLabelsPtr_);
    deleteDemandDrivenData(insertedFaceAddressingPtr_);
}


//- Calculate the insertedFaceLabels list
void topoPatchMapper::calcInsertedFaceAddressing() const
{
    if (insertedFaceLabelsPtr_ || insertedFaceAddressingPtr_)
    {
        FatalErrorIn
        (
            "void topoPatchMapper::calcInsertedFaceAddressing() const"
        )   << " Inserted labels has already been calculated."
            << abort(FatalError);
    }

    // Information from the old patch
    const label oldPatchSize = mpm_.oldPatchSizes()[patch_.index()];
    const label oldPatchStart = mpm_.oldPatchStarts()[patch_.index()];
    const label oldPatchEnd = oldPatchStart + oldPatchSize;

    // Allocate for inserted face labels and addressing
    label nInsertedFaces = 0;

    insertedFaceLabelsPtr_ = new labelList(size(), -1);
    labelList& insertedFaces = *insertedFaceLabelsPtr_;

    insertedFaceAddressingPtr_ = new labelListList(size(), labelList(0));
    labelListList& insertedAddressing = *insertedFaceAddressingPtr_;

    // Fetch the current boundary
    const polyBoundaryMesh& boundary = mpm_.mesh().boundaryMesh();

    // Loop through the facesFromFaces map, and ensure that
    // inserted faces are only mapped from faces on the same patch.
    const List<objectMap>& fff = mpm_.facesFromFacesMap();

    forAll(fff, objectI)
    {
        const objectMap& fffI = fff[objectI];

        // Only pick boundary faces in this patch
        if (boundary.whichPatch(fffI.index()) == patch_.index())
        {
            if (fffI.masterObjects().empty())
            {
                FatalErrorIn
                (
                    "void topoPatchMapper::"
                    "calcInsertedFaceAddressing() const"
                )   << " Mapping for inserted boundary face is incorrect."
                    << " Found an empty masterObjects list."
                    << nl << " Face: " << fffI.index()
                    << nl << " Patch: " << patch_.name()
                    << abort(FatalError);
            }
            else
            {
                // Make an entry for the inserted label,
                // and renumber addressing to patch.
                insertedFaces[nInsertedFaces] =
                (
                    fffI.index() - patch_.patch().start()
                );

                // Make an entry for addressing
                labelList& addr = insertedAddressing[nInsertedFaces];

                // Renumber addressing to patch.
                // Also, check mapping for hits into
                // other patches / internal faces.
                addr = fffI.masterObjects();

                forAll(addr, faceI)
                {
                    if
                    (
                        addr[faceI] >= oldPatchStart
                     && addr[faceI] < oldPatchEnd
                    )
                    {
                        addr[faceI] -= oldPatchStart;
                    }
                    else
                    {
                        FatalErrorIn
                        (
                            "void topoPatchMapper::"
                            "calcInsertedFaceAddressing() const"
                        )
                            << "Addressing into another patch is not allowed."
                            << nl << " Patch face index: " << faceI
                            << nl << " addr[faceI]: " << addr[faceI]
                            << nl << " oldPatchStart: " << oldPatchStart
                            << nl << " oldPatchSize: " << oldPatchSize
                            << nl << " oldPatchEnd: " << oldPatchEnd
                            << abort(FatalError);
                    }
                }

                nInsertedFaces++;
            }
        }
    }

    // Shorten inserted faces to actual size
    insertedFaces.setSize(nInsertedFaces);
    insertedAddressing.setSize(nInsertedFaces);
}


//- Calculate addressing for mapping
void topoPatchMapper::calcAddressing() const
{
    if (directAddrPtr_ || interpolationAddrPtr_)
    {
        FatalErrorIn("void topoPatchMapper::calcAddressing() const")
            << "Addressing already calculated."
            << abort(FatalError);
    }

    // Information from the old patch
    const label oldPatchSize = mpm_.oldPatchSizes()[patch_.index()];
    const label oldPatchStart = mpm_.oldPatchStarts()[patch_.index()];
    const label oldPatchEnd = oldPatchStart + oldPatchSize;

    // Assemble the maps: slice to patch
    if (direct())
    {
        // Direct mapping - slice to size
        directAddrPtr_ = new labelList(patch_.patchSlice(mpm_.faceMap()));

        labelList& addr = *directAddrPtr_;

        // Shift to local patch indices.
        // Also, check mapping for hits into other patches / internal faces.
        forAll (addr, faceI)
        {
            if
            (
                addr[faceI] >= oldPatchStart
             && addr[faceI] < oldPatchEnd
            )
            {
                addr[faceI] -= oldPatchStart;
            }
            else
            {
                FatalErrorIn
                (
                    "void topoPatchMapper::calcAddressing() const"
                )
                    << "Addressing into another patch is not allowed."
                    << nl << " Patch face index: " << faceI
                    << nl << " addr[faceI]: " << addr[faceI]
                    << nl << " oldPatchStart: " << oldPatchStart
                    << nl << " oldPatchSize: " << oldPatchSize
                    << nl << " oldPatchEnd: " << oldPatchEnd
                    << abort(FatalError);
            }
        }
    }
    else
    {
        // Interpolative addressing
        interpolationAddrPtr_ = new labelListList(size(), labelList(0));
        labelListList& addr = *interpolationAddrPtr_;

        // Fetch the list of inserted faces / addressing
        const labelList& insertedFaces = insertedObjectLabels();
        const labelListList& insertedAddressing = insertedFaceAddressing();

        // Make entries
        forAll(insertedFaces, faceI)
        {
            addr[insertedFaces[faceI]] = insertedAddressing[faceI];
        }

        // Do mapped faces. Note that this can already be set by insertedFaces
        // so check if addressing size still zero.
        const labelList& fm = patch_.patchSlice(mpm_.faceMap());

        forAll(fm, faceI)
        {
            if (fm[faceI] > -1 && addr[faceI].size() == 0)
            {
                // Mapped from a single face
                label oldFace = fm[faceI];

                if
                (
                    oldFace >= oldPatchStart
                 && oldFace < oldPatchEnd
                )
                {
                    oldFace -= oldPatchStart;
                }
                else
                {
                    FatalErrorIn
                    (
                        "void topoPatchMapper::calcAddressing() const"
                    )
                        << "Addressing into another patch is not allowed."
                        << nl << " Patch face index: " << faceI
                        << nl << " faceMap[faceI]: " << oldFace
                        << nl << " oldPatchStart: " << oldPatchStart
                        << nl << " oldPatchSize: " << oldPatchSize
                        << nl << " oldPatchEnd: " << oldPatchEnd
                        << abort(FatalError);
                }

                addr[faceI] = labelList(1, oldFace);
            }
        }

        // Check if we missed anything
        forAll(addr, faceI)
        {
            if (addr[faceI].empty())
            {
                FatalErrorIn
                (
                    "void topoPatchMapper::calcAddressing() const"
                )
                    << "Addressing is missing."
                    << nl << " Patch face index: " << faceI
                    << abort(FatalError);
            }
        }
    }
}


//- Calculate inverse-distance weights for interpolative mapping
void topoPatchMapper::calcInverseDistanceWeights() const
{
    if (weightsPtr_)
    {
        FatalErrorIn
        (
            "void topoPatchMapper::calcInverseDistanceWeights() const"
        )
            << "Weights already calculated."
            << abort(FatalError);
    }

    // Fetch interpolative addressing
    const labelListList& addr = addressing();

    // Allocate memory
    weightsPtr_ = new scalarListList(size());
    scalarListList& w = *weightsPtr_;

    // Obtain cell-centre information from old/new meshes
    const vectorField& oldCentres = tMapper_.patchCentres(patch_.index());
    const vectorField& newCentres = patch_.patch().faceCentres();

    forAll(addr, faceI)
    {
        const labelList& mo = addr[faceI];

        // Do mapped faces
        if (mo.size() == 1)
        {
            w[faceI] = scalarList(1, 1.0);
        }
        else
        {
            // Map from masters, inverse-distance weights
            scalar totalWeight = 0.0;
            w[faceI] = scalarList(mo.size(), 0.0);

            forAll (mo, oldFaceI)
            {
                w[faceI][oldFaceI] =
                (
                    1.0/stabilise
                    (
                        magSqr
                        (
                            newCentres[faceI]
                          - oldCentres[mo[oldFaceI]]
                        ),
                        VSMALL
                    )
                );

                totalWeight += w[faceI][oldFaceI];
            }

            // Normalize weights
            scalar normFactor = (1.0/totalWeight);

            forAll (mo, oldFaceI)
            {
                w[faceI][oldFaceI] *= normFactor;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
topoPatchMapper::topoPatchMapper
(
    const fvPatch& patch,
    const mapPolyMesh& mpm,
    const topoMapper& mapper
)
:
    patch_(patch),
    mpm_(mpm),
    tMapper_(mapper),
    direct_(false),
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL),
    insertedFaceLabelsPtr_(NULL),
    insertedFaceAddressingPtr_(NULL)
{
    // Check for the possibility of direct mapping
    if (insertedObjects())
    {
        direct_ = false;
    }
    else
    {
        direct_ = true;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

topoPatchMapper::~topoPatchMapper()
{
    clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return size
label topoPatchMapper::size() const
{
    return patch_.size();
}


//- Return size before mapping
label topoPatchMapper::sizeBeforeMapping() const
{
    return mpm_.oldPatchSizes()[patch_.index()];
}


//- Is the mapping direct
bool topoPatchMapper::direct() const
{
    return direct_;
}


//- Return direct addressing
const unallocLabelList& topoPatchMapper::directAddressing() const
{
    if (!direct())
    {
        FatalErrorIn
        (
            "const unallocLabelList& "
            "topoPatchMapper::directAddressing() const"
        )   << "Requested direct addressing for an interpolative mapper."
            << abort(FatalError);
    }

    if (!directAddrPtr_)
    {
        calcAddressing();
    }

    return *directAddrPtr_;
}


//- Return interpolation addressing
const labelListList& topoPatchMapper::addressing() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const labelListList& "
            "topoPatchMapper::addressing() const"
        )   << "Requested interpolative addressing for a direct mapper."
            << abort(FatalError);
    }

    if (!interpolationAddrPtr_)
    {
        calcAddressing();
    }

    return *interpolationAddrPtr_;
}


//- Return weights
const scalarListList& topoPatchMapper::weights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const scalarListList& "
            "topoPatchMapper::weights() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!weightsPtr_)
    {
        calcInverseDistanceWeights();
    }

    return *weightsPtr_;
}


//- Are there any inserted faces
bool topoPatchMapper::insertedObjects() const
{
    return insertedObjectLabels().size();
}


//- Return list of inserted faces
const labelList& topoPatchMapper::insertedObjectLabels() const
{
    if (!insertedFaceLabelsPtr_)
    {
        calcInsertedFaceAddressing();
    }

    return *insertedFaceLabelsPtr_;
}


//- Return addressing for inserted faces
const labelListList& topoPatchMapper::insertedFaceAddressing() const
{
    if (!insertedFaceAddressingPtr_)
    {
        calcInsertedFaceAddressing();
    }

    return *insertedFaceAddressingPtr_;
}


//- Map the patch field
template <class Type>
void topoPatchMapper::mapPatchField(Field<Type>& pF) const
{
    pF.autoMap(*this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void topoPatchMapper::operator=(const topoPatchMapper& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("void topoPatchMapper::operator=")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
