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
    interpolator

Description
    Macro definitions for code brevity

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#define mapCell(type, capsType)                                            \
    forAllIter(HashTable<Map<type> >, vol##capsType##Map_, hIter)          \
    {                                                                      \
        type mapVal = pTraits<type>::zero;                                 \
                                                                           \
        if (dynamicTopoFvMesh::debug > 3)                                  \
        {                                                                  \
            Info << nl << "Mapping: " << hIter.key() << endl;              \
        }                                                                  \
                                                                           \
        forAll(mapCells, cellI)                                            \
        {                                                                  \
            type oVal = pTraits<type>::zero;                               \
                                                                           \
            if                                                             \
            (                                                              \
                (mapCells[cellI] < nOldCells) &&                           \
                (!hIter().found(mapCells[cellI]))                          \
            )                                                              \
            {                                                              \
                const vol##capsType##Field& oF =                           \
                (                                                          \
                    mesh_.lookupObject<vol##capsType##Field>(hIter.key())  \
                );                                                         \
                                                                           \
                oVal = oF.internalField()[mapCells[cellI]];                \
            }                                                              \
            else                                                           \
            {                                                              \
                oVal = hIter()[mapCells[cellI]];                           \
            }                                                              \
                                                                           \
            if (dynamicTopoFvMesh::debug > 3)                              \
            {                                                              \
                Info << "\tCell:" << mapCells[cellI]                       \
                     << " Value: " << oVal << endl;                        \
            }                                                              \
                                                                           \
            mapVal += (mapWeights[cellI]*oVal);                            \
        }                                                                  \
                                                                           \
        if (dynamicTopoFvMesh::debug > 3)                                  \
        {                                                                  \
            Info << "\tNew value: " << mapVal << endl;                     \
        }                                                                  \
                                                                           \
        hIter().set(newCellIndex, mapVal);                                 \
    }

#define mapBoundaryFace(type, capsType)                                    \
    forAllIter(HashTable<Map<type> >, surf##capsType##Map_, hIter)         \
    {                                                                      \
        type mapVal = pTraits<type>::zero;                                 \
                                                                           \
        if (dynamicTopoFvMesh::debug > 3)                                  \
        {                                                                  \
            Info << nl << "Mapping: " << hIter.key() << endl;              \
        }                                                                  \
                                                                           \
        forAll(mapFaces, faceI)                                            \
        {                                                                  \
            type oVal = pTraits<type>::zero;                               \
                                                                           \
            if                                                             \
            (                                                              \
                (mapFaces[faceI] < nOldFaces) &&                           \
                (!hIter().found(mapFaces[faceI]))                          \
            )                                                              \
            {                                                              \
                const vol##capsType##Field& oF =                           \
                (                                                          \
                    mesh_.lookupObject<vol##capsType##Field>(hIter.key())  \
                );                                                         \
                                                                           \
                label oPatch = boundary.whichPatch(mapFaces[faceI]);       \
                                                                           \
                if (oPatch == -1)                                          \
                {                                                          \
                    FatalErrorIn("interpolator::insertFace()")             \
                        << nl                                              \
                        << " Cannot map from internal to boundary face."   \
                        << abort(FatalError);                              \
                }                                                          \
                                                                           \
                label i = boundary[oPatch].whichFace(mapFaces[faceI]);     \
                                                                           \
                oVal = oF.boundaryField()[oPatch][i];                      \
            }                                                              \
            else                                                           \
            {                                                              \
                oVal = hIter()[mapFaces[faceI]];                           \
            }                                                              \
                                                                           \
            if (dynamicTopoFvMesh::debug > 3)                              \
            {                                                              \
                Info << "\tFace:" << mapFaces[faceI]                       \
                     << " Value: " << oVal << endl;                        \
            }                                                              \
                                                                           \
            mapVal += (mapWeights[faceI]*oVal);                            \
        }                                                                  \
                                                                           \
        if (dynamicTopoFvMesh::debug > 3)                                  \
        {                                                                  \
            Info << "\tNew value: " << mapVal << endl;                     \
        }                                                                  \
                                                                           \
        hIter().set(newFaceIndex, mapVal);                                 \
    }

#define mapInternalFace(type, capsType)                                    \
    forAllIter(HashTable<Map<type> >, surf##capsType##Map_, hIter)         \
    {                                                                      \
        type mapVal = pTraits<type>::zero;                                 \
                                                                           \
        if (vol##capsType##Map_.found(hIter.key()))                        \
        {                                                                  \
            continue;                                                      \
        }                                                                  \
                                                                           \
        if (dynamicTopoFvMesh::debug > 3)                                  \
        {                                                                  \
            Info << nl << "Mapping: " << hIter.key() << endl;              \
        }                                                                  \
                                                                           \
        forAll(mapFaces, faceI)                                            \
        {                                                                  \
            type oVal = pTraits<type>::zero;                               \
                                                                           \
            if                                                             \
            (                                                              \
                (mapFaces[faceI] < nOldFaces) &&                           \
                (!hIter().found(mapFaces[faceI]))                          \
            )                                                              \
            {                                                              \
                const surface##capsType##Field& oF =                       \
                (                                                          \
                    mesh_.lookupObject<surface##capsType##Field>           \
                    (                                                      \
                        hIter.key()                                        \
                    )                                                      \
                );                                                         \
                                                                           \
                label oPatch = boundary.whichPatch(mapFaces[faceI]);       \
                                                                           \
                if (oPatch == -1)                                          \
                {                                                          \
                    if (patch != -1)                                       \
                    {                                                      \
                        FatalErrorIn("interpolator::insertFace()")         \
                            << nl                                          \
                            << " Cannot map from internal to boundary."    \
                            << abort(FatalError);                          \
                    }                                                      \
                                                                           \
                    oVal = oF.internalField()[mapFaces[faceI]];            \
                }                                                          \
                else                                                       \
                {                                                          \
                    if (patch == -1)                                       \
                    {                                                      \
                        FatalErrorIn("interpolator::insertFace()")         \
                            << nl                                          \
                            << " Cannot map from boundary to internal."    \
                            << abort(FatalError);                          \
                    }                                                      \
                                                                           \
                    label i = boundary[oPatch].whichFace(mapFaces[faceI]); \
                                                                           \
                    oVal = oF.boundaryField()[oPatch][i];                  \
                }                                                          \
            }                                                              \
            else                                                           \
            {                                                              \
                oVal = hIter()[mapFaces[faceI]];                           \
            }                                                              \
                                                                           \
            if (dynamicTopoFvMesh::debug > 3)                              \
            {                                                              \
                Info << "\tFace:" << mapFaces[faceI]                       \
                     << " Value: " << oVal << endl;                        \
            }                                                              \
                                                                           \
            mapVal += (mapWeights[faceI]*oVal);                            \
        }                                                                  \
                                                                           \
        if (dynamicTopoFvMesh::debug > 3)                                  \
        {                                                                  \
            Info << "\tNew value: " << mapVal << endl;                     \
        }                                                                  \
                                                                           \
        hIter().set(newFaceIndex, mapVal);                                 \
    }

#define registerVolumeField(type, capsType)                                \
    HashTable<const vol##capsType##Field*> v##type##f                      \
    (                                                                      \
        mesh_.lookupClass<vol##capsType##Field>()                          \
    );                                                                     \
                                                                           \
    forAllIter(HashTable<const vol##capsType##Field*>, v##type##f, fIter)  \
    {                                                                      \
        if (dynamicTopoFvMesh::debug)                                      \
        {                                                                  \
            Info << "Registering: " << fIter()->name() << endl;            \
        }                                                                  \
                                                                           \
        vol##capsType##Map_.insert(fIter()->name(), Map<type>());          \
        surf##capsType##Map_.insert(fIter()->name(), Map<type>());         \
    }

#define registerSurfaceField(type, capsType)                               \
    HashTable<const surface##capsType##Field*> s##type##f                  \
    (                                                                      \
        mesh_.lookupClass<surface##capsType##Field>()                      \
    );                                                                     \
                                                                           \
    forAllIter(HashTable<const surface##capsType##Field*>, s##type##f, f)  \
    {                                                                      \
        if (dynamicTopoFvMesh::debug)                                      \
        {                                                                  \
            Info << "Registering: " << f()->name() << endl;                \
        }                                                                  \
                                                                           \
        surf##capsType##Map_.insert(f()->name(), Map<type>());             \
    }

#define removeEntity(type, capsType)                                       \
    forAllIter(HashTable<Map<type> >, capsType##Map_, hIter)               \
    {                                                                      \
        if (hIter().found(index))                                          \
        {                                                                  \
            hIter().erase(index);                                          \
        }                                                                  \
    }

#define clearInterpolationMaps(type, capsType)                             \
    forAllIter(HashTable<Map<type> >, vol##capsType##Map_, hIter)          \
    {                                                                      \
        hIter().clear();                                                   \
    }                                                                      \
                                                                           \
    forAllIter(HashTable<Map<type> >, surf##capsType##Map_, hIter)         \
    {                                                                      \
        hIter().clear();                                                   \
    }                                                                      \

#define fillVolumeMaps(type, capsType)                                     \
    forAllIter(HashTable<Map<type> >, vol##capsType##Map_, hIter)          \
    {                                                                      \
        vol##capsType##Field& vf =                                         \
        (                                                                  \
            const_cast<vol##capsType##Field&>                              \
            (                                                              \
                mesh_.lookupObject<vol##capsType##Field>(hIter.key())      \
            )                                                              \
        );                                                                 \
                                                                           \
        if (dynamicTopoFvMesh::debug > 2)                                  \
        {                                                                  \
            Info << nl << "Updating: " << hIter.key() << endl;             \
        }                                                                  \
                                                                           \
        const labelList& reverseMap = mpm.reverseCellMap();                \
                                                                           \
        forAllConstIter(Map<type>, hIter(), mIter)                         \
        {                                                                  \
            label newIndex = -1;                                           \
                                                                           \
            if (mIter.key() < mpm.nOldCells())                             \
            {                                                              \
                newIndex = reverseMap[mIter.key()];                        \
            }                                                              \
            else                                                           \
            {                                                              \
                newIndex = mesh_.addedCellRenumbering_[mIter.key()];       \
            }                                                              \
                                                                           \
            vf.internalField()[newIndex] = mIter();                        \
        }                                                                  \
    }

#define fillBoundaryMaps(type, capsType)                                   \
    forAllIter(HashTable<Map<type> >, surf##capsType##Map_, hIter)         \
    {                                                                      \
        bool found =                                                       \
        (                                                                  \
            mesh_.foundObject<vol##capsType##Field>(hIter.key())           \
        );                                                                 \
                                                                           \
        if (!found)                                                        \
        {                                                                  \
            continue;                                                      \
        }                                                                  \
                                                                           \
        if (dynamicTopoFvMesh::debug > 2)                                  \
        {                                                                  \
            Info << nl << "Updating: " << hIter.key() << endl;             \
        }                                                                  \
                                                                           \
        vol##capsType##Field& vf =                                         \
        (                                                                  \
            const_cast<vol##capsType##Field&>                              \
            (                                                              \
                mesh_.lookupObject<vol##capsType##Field>(hIter.key())      \
            )                                                              \
        );                                                                 \
                                                                           \
        const labelList& reverseMap = mpm.reverseFaceMap();                \
                                                                           \
        forAllConstIter(Map<type>, hIter(), mIter)                         \
        {                                                                  \
            label newIndex = -1;                                           \
                                                                           \
            if (mIter.key() < mpm.nOldFaces())                             \
            {                                                              \
                newIndex = reverseMap[mIter.key()];                        \
            }                                                              \
            else                                                           \
            {                                                              \
                newIndex = mesh_.addedFaceRenumbering_[mIter.key()];       \
            }                                                              \
                                                                           \
            label patch = boundary.whichPatch(newIndex);                   \
                                                                           \
            if (patch == -1)                                               \
            {                                                              \
                FatalErrorIn("interpolator::updateMesh()")                 \
                    << nl                                                  \
                    << " Encountered an internal face"                     \
                    << " when mapping for the boundary." << nl             \
                    << " faceIndex: " << mIter.key() << nl                 \
                    << " newIndex: " << newIndex << nl                     \
                    << abort(FatalError);                                  \
            }                                                              \
            else                                                           \
            {                                                              \
                label i = boundary[patch].whichFace(newIndex);             \
                                                                           \
                vf.boundaryField()[patch][i] = mIter();                    \
            }                                                              \
        }                                                                  \
    }

#define fillSurfaceMaps(type, capsType)                                    \
    forAllIter(HashTable<Map<type> >, surf##capsType##Map_, hIter)         \
    {                                                                      \
        bool found =                                                       \
        (                                                                  \
            mesh_.foundObject<surface##capsType##Field>(hIter.key())       \
        );                                                                 \
                                                                           \
        if (!found)                                                        \
        {                                                                  \
            continue;                                                      \
        }                                                                  \
                                                                           \
        if (dynamicTopoFvMesh::debug > 2)                                  \
        {                                                                  \
            Info << nl << "Updating: " << hIter.key() << endl;             \
        }                                                                  \
                                                                           \
        surface##capsType##Field& sf =                                     \
        (                                                                  \
            const_cast<surface##capsType##Field&>                          \
            (                                                              \
                mesh_.lookupObject<surface##capsType##Field>(hIter.key())  \
            )                                                              \
        );                                                                 \
                                                                           \
        const labelList& reverseMap = mpm.reverseFaceMap();                \
                                                                           \
        forAllConstIter(Map<type>, hIter(), mIter)                         \
        {                                                                  \
            label newIndex = -1;                                           \
                                                                           \
            if (mIter.key() < mpm.nOldFaces())                             \
            {                                                              \
                newIndex = reverseMap[mIter.key()];                        \
            }                                                              \
            else                                                           \
            {                                                              \
                newIndex = mesh_.addedFaceRenumbering_[mIter.key()];       \
            }                                                              \
                                                                           \
            label patch = boundary.whichPatch(newIndex);                   \
                                                                           \
            if (patch == -1)                                               \
            {                                                              \
                sf.internalField()[newIndex] = mIter();                    \
            }                                                              \
            else                                                           \
            {                                                              \
                label i = boundary[patch].whichFace(newIndex);             \
                                                                           \
                sf.boundaryField()[patch][i] = mIter();                    \
            }                                                              \
        }                                                                  \
    }

// ************************************************************************* //
