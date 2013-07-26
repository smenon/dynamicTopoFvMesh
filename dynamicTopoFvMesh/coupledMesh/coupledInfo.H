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
    coupledInfo

Description
    An interface class that provides patch coupling functionality

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    coupledInfo.C

\*---------------------------------------------------------------------------*/

#ifndef coupledInfo_H
#define coupledInfo_H

#include "autoPtr.H"
#include "coupleMap.H"
#include "volFieldsFwd.H"
#include "fvPatchFieldMapper.H"

namespace Foam
{

// Class forward declarations
class fvBoundaryMesh;
class dynamicTopoFvMesh;

/*---------------------------------------------------------------------------*\
                         Class coupledInfo Declaration
\*---------------------------------------------------------------------------*/

class coupledInfo
{
    // Reference to the parent mesh
    const dynamicTopoFvMesh& mesh_;

    // Auto pointer to a subMesh
    autoPtr<dynamicTopoFvMesh> subMesh_;

    // Flag to determine whether maps have been built.
    bool builtMaps_;

    // For locally coupled patches,
    // specify a master / slave index
    coupleMap map_;

    // Zone IDs for patches associated with faceZones
    label masterFaceZone_;
    label slaveFaceZone_;

    //- Disallow default bitwise assignment
    inline void operator=(const coupledInfo&);

public:

    // Constructor
    inline coupledInfo
    (
        const dynamicTopoFvMesh& mesh,
        const coupleMap& cMap,
        const label mfzIndex = -1,
        const label sfzIndex = -1
    );

    inline coupledInfo
    (
        const dynamicTopoFvMesh& mesh,
        const bool isTwoDMesh,
        const bool isLocal,
        const bool isSend,
        const label patchIndex,
        const label mPatch,
        const label sPatch,
        const label mfzIndex = -1,
        const label sfzIndex = -1
    );

    //- Access

        // Return a const reference to the parent mesh
        inline const dynamicTopoFvMesh& baseMesh() const;

        // Set a new subMesh
        inline void setMesh(label index, dynamicTopoFvMesh* mesh);

        // Return a reference to the subMesh
        inline dynamicTopoFvMesh& subMesh();

        // Return a const reference to the subMesh
        inline const dynamicTopoFvMesh& subMesh() const;

        // Have maps been built?
        inline bool builtMaps() const;

        // Change the flag
        inline void setBuiltMaps();

        // Return a reference to the coupleMap
        inline coupleMap& map();

        // Return a const reference to the coupleMap
        inline const coupleMap& map() const;

        // Return the master / slave face zone IDs
        inline label masterFaceZone() const;
        inline label slaveFaceZone() const;

    //- Interpolation

        //- Generic subMesh mapper
        class subMeshMapper
        :
            public fvPatchFieldMapper
        {
            label sizeBeforeMapping_;

            labelField directAddressing_;

        public:

            // Constructors

                //- Construct from components
                inline subMeshMapper
                (
                    const label sbm,
                    const labelList& da
                )
                :
                    sizeBeforeMapping_(sbm),
                    directAddressing_(da)
                {}

                //- Construct given addressing
                inline subMeshMapper
                (
                    const coupledInfo& cInfo,
                    const label patchI
                );

            // Destructor

                virtual ~subMeshMapper()
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
        };

        // Subset geometric field
        template<class GeomField>
        tmp<GeomField>
        subSetField(const GeomField& field) const;

        // Subset geometric fields from registry to output stream
        template<class GeomField>
        void send
        (
            const wordList& fieldNames,
            const word& fieldType,
            OSstream& strStream
        ) const;

        // Set geometric field pointer from input dictionary
        template<class Type, template<class> class PatchField, class Mesh>
        void setField
        (
            const wordList& fieldNames,
            const dictionary& fieldDicts,
            PtrList<GeometricField<Type, PatchField, Mesh> >& fields
        ) const;

        // Resize map for individual field
        template <class GeomField>
        static void resizeMap
        (
            const label srcIndex,
            const subMeshMapper& internalMapper,
            const List<labelList>& internalReverseMaps,
            const PtrList<subMeshMapper>& boundaryMapper,
            const List<labelListList>& boundaryReverseMaps,
            const List<PtrList<GeomField> >& srcFields,
            GeomField& field
        );

        // Resize map for all fields in registry
        template <class GeomField>
        static void resizeMap
        (
            const wordList& fieldNames,
            const objectRegistry& mesh,
            const subMeshMapper& internalMapper,
            const List<labelList>& internalReverseMaps,
            const PtrList<subMeshMapper>& boundaryMapper,
            const List<labelListList>& boundaryReverseMaps,
            const List<PtrList<GeomField> >& srcFields
        );

        // Resize boundaryFields for all fields in the registry
        template <class GeomField>
        static void resizeBoundaries
        (
            const objectRegistry& mesh,
            const fvBoundaryMesh& boundary
        );
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "coupledInfo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
