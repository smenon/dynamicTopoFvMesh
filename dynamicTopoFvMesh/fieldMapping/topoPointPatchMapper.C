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
    topoPointPatchMapper

Description
    Implementation of the topoPointPatchMapper class

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "meshOps.H"
#include "topoMapper.H"
#include "mapPolyMesh.H"
#include "facePointPatch.H"
#include "demandDrivenData.H"
#include "topoPointPatchMapper.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Calculate the inserted point addressing list
void topoPointPatchMapper::calcInsertedAddressing() const
{
    if (insertedPointLabelsPtr_ || insertedPointAddressingPtr_)
    {
        FatalErrorIn
        (
            "void topoPointPatchMapper::calcInsertedAddressing() const"
        )   << " Inserted labels has already been calculated."
            << abort(FatalError);
    }

    // Allocate for inserted point labels and addressing
    label nInsertedPoints = 0;

    insertedPointLabelsPtr_ = new labelList(size(), -1);
    labelList& insertedPoints = *insertedPointLabelsPtr_;

    insertedPointAddressingPtr_ = new labelListList(size(), labelList(0));
    labelListList& insertedAddressing = *insertedPointAddressingPtr_;

    // Fetch the point maps for this patch
    const Map<label>& oldMeshPointMap = patchMeshPointMap();
    const Map<label>& newMeshPointMap = patch_.patch().meshPointMap();

    // Loop through the pointsFromPoints map,
    // and pick inserted points on the same patch.
    const List<objectMap>& pfp = mpm_.pointsFromPointsMap();

    Map<label>::const_iterator poIter, pnIter;

    forAll(pfp, pfpI)
    {
        const objectMap& pObj = pfp[pfpI];
        const label pIndex = pObj.index();

        // Check if the index belongs to this patch
        pnIter = newMeshPointMap.find(pIndex);

        if (pnIter == newMeshPointMap.end())
        {
            continue;
        }

        // Ensure that point is mapped
        const labelList& mo = pObj.masterObjects();

        if (mo.empty())
        {
            FatalErrorIn
            (
                "void topoPointPatchMapper::calcInsertedAddressing() const"
            )   << " Mapping for inserted boundary point is incorrect."
                << " Found an empty masterObjects list."
                << nl << " Point: " << pIndex
                << nl << " Patch: " << patch_.name()
                << abort(FatalError);
        }

        // Make an entry for the inserted label,
        // and renumber addressing to patch.
        insertedPoints[nInsertedPoints] = pnIter();

        // Make an entry for addressing
        labelList& addr = insertedAddressing[nInsertedPoints];

        addr.setSize(mo.size());

        forAll(mo, indexI)
        {
            const label oldIndex = mo[indexI];

            // Check if the index belongs to this patch
            poIter = oldMeshPointMap.find(oldIndex);

            if (poIter == oldMeshPointMap.end())
            {
                // Write out for post-processing
                meshOps::writeVTK
                (
                    mpm_.mesh(),
                    "patchPointPatchError_"
                  + Foam::name(pObj.index()),
                    labelList(1, pObj.index()),
                    0, mpm_.mesh().points()
                );

                FatalErrorIn
                (
                    "void topoPointPatchMapper::calcInsertedAddressing() const"
                )   << " Mapping for inserted boundary point is incorrect."
                    << " Found an master point that doesn't belong to the patch"
                    << nl << " Point: " << pIndex
                    << nl << " Master: " << oldIndex
                    << nl << " MasterObjects: " << mo
                    << nl << " Patch: " << patch_.name()
                    << abort(FatalError);
            }

            // Add renumbered entity
            addr[indexI] = poIter();
        }

        nInsertedPoints++;
    }

    // Shorten inserted points to actual size
    insertedPoints.setSize(nInsertedPoints);
    insertedAddressing.setSize(nInsertedPoints);
}


//- Calculate the patch mesh point map
void topoPointPatchMapper::calcPatchMeshPointMap() const
{
    if (patchMeshPointMapPtr_)
    {
        FatalErrorIn("void topoPointPatchMapper::calcPatchMeshPointMap() const")
            << "Mesh point map already calculated."
            << abort(FatalError);
    }

    // Fetch the mesh points for this patch
    const label patchIndex = patch_.index();
    const labelList& pMeshPoints = tMapper_.oldPatchMeshPoints()[patchIndex];

    // Allocate the map
    patchMeshPointMapPtr_ = new Map<label>(2 * pMeshPoints.size());

    // Insert point addressing
    Map<label>& mpMap = *patchMeshPointMapPtr_;

    forAll(pMeshPoints, pointI)
    {
        mpMap.insert(pMeshPoints[pointI], pointI);
    }
}


//- Calculate addressing for interpolative mapping
void topoPointPatchMapper::calcAddressing() const
{
    if (directAddrPtr_ || interpolationAddrPtr_)
    {
        FatalErrorIn("void topoPointPatchMapper::calcAddressing() const")
            << "Addressing already calculated."
            << abort(FatalError);
    }

    // Track unmapped points
    label nUnmappedPoints = 0;
    labelList unmappedPoints(size());

    // Fetch the current point maps for this patch
    const labelList& patchMap = mpm_.patchPointMap()[patch_.index()];

    if (direct())
    {
        // Direct addressing
        directAddrPtr_ = new labelList(size());

        labelList& addr = *directAddrPtr_;

        // Map all indices
        forAll(addr, pointI)
        {
            const label oldIndex = patchMap[pointI];

            if (oldIndex > -1)
            {
                addr[pointI] = oldIndex;
            }
            else
            {
                FatalErrorIn
                (
                    "void topoPointPatchMapper::calcAddressing() const"
                )   << " Mapping for inserted boundary point is incorrect."
                    << " Found a point which is not mapped from anything."
                    << nl << " Point: " << pointI
                    << nl << " Patch: " << patch_.name()
                    << abort(FatalError);
            }
        }
    }
    else
    {
        // Interpolative addressing
        interpolationAddrPtr_ = new labelListList(size());
        labelListList& addr = *interpolationAddrPtr_;

        // Fetch the list of inserted points / addressing
        const labelList& insertedPoints = insertedObjectLabels();
        const labelListList& insertedAddressing = insertedPointAddressing();

        // Make entries for inserted points
        forAll(insertedPoints, pointI)
        {
            addr[insertedPoints[pointI]] = insertedAddressing[pointI];
        }

        // Do mapped points.
        // Note that this can already be set by insertedPoints,
        // so check if addressing size still zero.
        forAll(patchMap, pointI)
        {
            // Mapped from a single point
            if (patchMap[pointI] > -1)
            {
                if (addr[pointI].empty())
                {
                    addr[pointI] = labelList(1, patchMap[pointI]);
                }
                else
                {
                    FatalErrorIn
                    (
                        "void topoPointPatchMapper::calcAddressing() const"
                    )
                        << "Master point " << pointI
                        << " to be mapped from: " << patchMap[pointI]
                        << " is already mapped from points: " << addr[pointI]
                        << abort(FatalError);
                }
            }

            // Check for unmapped points without any addressing
            if (patchMap[pointI] < 0 && addr[pointI].empty())
            {
                unmappedPoints[nUnmappedPoints++] = pointI;
            }
        }
    }

    // Shorten inserted points to actual size
    unmappedPoints.setSize(nUnmappedPoints);

    if (nUnmappedPoints)
    {
        FatalErrorIn("void topoPointPatchMapper::calcAddressing() const")
            << " Found " << nUnmappedPoints << " points which are"
            << " not mapped from any parent points." << nl
            << " Patch: " << patch_.name() << nl
            << " List: " << nl
            << unmappedPoints
            << abort(FatalError);
    }
}


//- Calculate weights for interpolative mapping
void topoPointPatchMapper::calcWeights() const
{
    if (weightsPtr_)
    {
        FatalErrorIn("void topoPointPatchMapper::calcWeights() const")
            << "Weights already calculated."
            << abort(FatalError);
    }

    // Fetch interpolative addressing
    const labelListList& addr = addressing();

    // Allocate memory
    weightsPtr_ = new scalarListList(size());
    scalarListList& weights = *weightsPtr_;

    forAll(addr, pointI)
    {
        const labelList& mo = addr[pointI];

        // Specify equal weights
        const scalar weight = (1.0 / scalar(mo.size()));

        weights[pointI] = scalarList(mo.size(), weight);
    }
}


//- Clear out local storage
void topoPointPatchMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
    deleteDemandDrivenData(interpolationAddrPtr_);
    deleteDemandDrivenData(weightsPtr_);
    deleteDemandDrivenData(insertedPointLabelsPtr_);
    deleteDemandDrivenData(insertedPointAddressingPtr_);
    deleteDemandDrivenData(patchMeshPointMapPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
topoPointPatchMapper::topoPointPatchMapper
(
    const pointPatch& patch,
    const mapPolyMesh& mpm,
    const topoMapper& mapper
)
:
    patch_(refCast<const facePointPatch>(patch)),
    mpm_(mpm),
    tMapper_(mapper),
    direct_(false),
    sizeBeforeMapping_(0),
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL),
    insertedPointLabelsPtr_(NULL),
    insertedPointAddressingPtr_(NULL),
    patchMeshPointMapPtr_(NULL)
{
    // Compute sizeBeforeMapping
    {
        label patchIndex = patch_.index();
        label totalSize = mpm_.oldPatchNMeshPoints()[patchIndex];

        // Fetch offset sizes from topoMapper
        const labelListList& sizes = tMapper_.pointPatchSizes();

        // Add offset sizes
        if (sizes.size())
        {
            // Fetch number of physical patches
            label nPhysical = sizes[0].size();

            if (patchIndex < nPhysical)
            {
                forAll(sizes, pI)
                {
                    totalSize += sizes[pI][patchIndex];
                }
            }
        }

        sizeBeforeMapping_ = totalSize;
    }

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

topoPointPatchMapper::~topoPointPatchMapper()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return size
label topoPointPatchMapper::size() const
{
    return patch_.size();
}


//- Return size before mapping
label topoPointPatchMapper::sizeBeforeMapping() const
{
    return sizeBeforeMapping_;
}


//- Is the mapping direct
bool topoPointPatchMapper::direct() const
{
    return direct_;
}


//- Has unmapped elements
bool topoPointPatchMapper::hasUnmapped() const
{
    return false;
}


//- Return direct addressing
const unallocLabelList& topoPointPatchMapper::directAddressing() const
{
    if (!direct())
    {
        FatalErrorIn
        (
            "const unallocLabelList& "
            "topoPointPatchMapper::directAddressing() const"
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
const labelListList& topoPointPatchMapper::addressing() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const labelListList& "
            "topoPointPatchMapper::addressing() const"
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
const scalarListList& topoPointPatchMapper::weights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const scalarListList& "
            "topoPointPatchMapper::weights() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!weightsPtr_)
    {
        calcWeights();
    }

    return *weightsPtr_;
}


//- Are there any inserted objects
bool topoPointPatchMapper::insertedObjects() const
{
    return insertedObjectLabels().size();
}


//- Return list of inserted objects
const labelList& topoPointPatchMapper::insertedObjectLabels() const
{
    if (!insertedPointLabelsPtr_)
    {
        calcInsertedAddressing();
    }

    return *insertedPointLabelsPtr_;
}


//- Return addressing for inserted points
const labelListList& topoPointPatchMapper::insertedPointAddressing() const
{
    if (!insertedPointAddressingPtr_)
    {
        calcInsertedAddressing();
    }

    return *insertedPointAddressingPtr_;
}


//- Return the patch mesh point map
const Map<label>& topoPointPatchMapper::patchMeshPointMap() const
{
    if (!patchMeshPointMapPtr_)
    {
        calcPatchMeshPointMap();
    }

    return *patchMeshPointMapPtr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void topoPointPatchMapper::operator=(const topoPointPatchMapper& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("void topoPointPatchMapper::operator=")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}


} // End namespace Foam

// ************************************************************************* //
