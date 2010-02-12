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
        forAll(mapCells, cellI)                                            \
        {                                                                  \
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
                mapVal +=                                                  \
                (                                                          \
                    mapWeights[cellI]*oF.internalField()[mapCells[cellI]]  \
                );                                                         \
            }                                                              \
            else                                                           \
            {                                                              \
                mapVal += (mapWeights[cellI]*hIter()[mapCells[cellI]]);    \
            }                                                              \
        }                                                                  \
                                                                           \
        hIter().set(newCellIndex, mapVal);                                 \
    }

#define mapBoundaryFace(type, capsType)                                    \
    forAllIter(HashTable<Map<type> >, vol##capsType##Map_, hIter)          \
    {                                                                      \
        type mapVal = pTraits<type>::zero;                                 \
                                                                           \
        forAll(mapFaces, faceI)                                            \
        {                                                                  \
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
                mapVal +=                                                  \
                (                                                          \
                    mapWeights[faceI]*oF.boundaryField()[oPatch][i]        \
                );                                                         \
            }                                                              \
            else                                                           \
            {                                                              \
                mapVal += (mapWeights[faceI]*hIter()[mapFaces[faceI]]);    \
            }                                                              \
        }                                                                  \
                                                                           \
        hIter().set(newFaceIndex, mapVal);                                 \
    }

#define mapInternalFace(type, capsType)                                    \
    forAllIter(HashTable<Map<type> >, surf##capsType##Map_, hIter)         \
    {                                                                      \
        type mapVal = pTraits<type>::zero;                                 \
                                                                           \
        forAll(mapFaces, faceI)                                            \
        {                                                                  \
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
                    mapVal +=                                              \
                    (                                                      \
                        mapWeights[faceI] *                                \
                        (                                                  \
                            oF.internalField()[mapFaces[faceI]]            \
                        )                                                  \
                    );                                                     \
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
                    mapVal +=                                              \
                    (                                                      \
                        mapWeights[faceI]*oF.boundaryField()[oPatch][i]    \
                    );                                                     \
                }                                                          \
            }                                                              \
            else                                                           \
            {                                                              \
                mapVal += (mapWeights[faceI]*hIter()[mapFaces[faceI]]);    \
            }                                                              \
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
        if (debug)                                                         \
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
        if (debug)                                                         \
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

// ************************************************************************* //
