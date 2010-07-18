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
    coupledPatchInfo

Description
    An interface class that provides patch coupling functionality

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    coupledPatchInfo.C

\*---------------------------------------------------------------------------*/

#ifndef coupledPatchInfo_H
#define coupledPatchInfo_H

#include "autoPtr.H"
#include "coupleMap.H"

namespace Foam
{

// Class forward declarations
class dynamicTopoFvMesh;

/*---------------------------------------------------------------------------*\
                         Class coupledPatchInfo Declaration
\*---------------------------------------------------------------------------*/

class coupledPatchInfo
{
    // Reference to the parent mesh
    const dynamicTopoFvMesh& mesh_;

    // Auto pointer to a subMesh
    autoPtr<dynamicTopoFvMesh> subMesh_;

    // Flag to determine whether maps have been built.
    bool builtMaps_;

    // For locally coupled patches,
    // specify a master / slave index
    coupleMap patchMap_;

    // Zone IDs for patches associated with faceZones
    label masterFaceZone_;
    label slaveFaceZone_;

    //- Disallow default bitwise assignment
    void operator=(const coupledPatchInfo&);

public:

    // Constructor
    coupledPatchInfo
    (
        const dynamicTopoFvMesh& mesh,
        const coupleMap& cMap,
        const label mfzIndex = -1,
        const label sfzIndex = -1
    );

    coupledPatchInfo
    (
        const dynamicTopoFvMesh& mesh,
        const bool isLocal,
        const bool isSend,
        const label mPatch,
        const label sPatch,
        const label mfzIndex = -1,
        const label sfzIndex = -1
    );

    //- Access

        // Set a new subMesh
        void setMesh(label index, dynamicTopoFvMesh* mesh);

        // Return a reference to the subMesh
        dynamicTopoFvMesh& subMesh();

        // Have maps been built?
        bool builtMaps() const;

        // Change the flag
        void setBuiltMaps();

        // Return the coupleMap
        const coupleMap& patchMap() const;

        // Return the master / slave face zone IDs
        label masterFaceZone() const;
        label slaveFaceZone() const;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
