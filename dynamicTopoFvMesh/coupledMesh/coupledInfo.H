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

namespace Foam
{

// Class forward declarations
class fvBoundaryMesh;
class faceSetAlgorithm;
class cellSetAlgorithm;
class pointSetAlgorithm;

/*---------------------------------------------------------------------------*\
                         Class coupledInfo Declaration
\*---------------------------------------------------------------------------*/

template <class MeshType>
class coupledInfo
{

    // Private data

        // Reference to the parent mesh
        const MeshType& mesh_;

        // Auto pointer to a subMesh
        autoPtr<MeshType> subMesh_;

        // Auto pointer to the cellSetAlgorithm
        autoPtr<cellSetAlgorithm> cellAlgorithm_;

        // Auto pointer to the faceSetAlgorithm
        autoPtr<faceSetAlgorithm> faceAlgorithm_;

        // Auto pointer to the pointSetAlgorithm
        autoPtr<pointSetAlgorithm> pointAlgorithm_;

        // Flag to determine whether maps have been built.
        bool builtMaps_;

        // For locally coupled patches,
        // specify a master / slave index
        coupleMap map_;

        // Zone IDs for patches associated with faceZones
        label masterFaceZone_;
        label slaveFaceZone_;

    // Private member functions

        // Disallow default bitwise assignment
        inline void operator=(const coupledInfo&);

public:

    //- Forward declarations

        class subMeshMapper;
        class subMeshPointMapper;

    //- Constructors

        // Construct given mesh, coupleMap and master / slave indices
        inline coupledInfo
        (
            const MeshType& mesh,
            const coupleMap& cMap,
            const label mfzIndex = -1,
            const label sfzIndex = -1
        );

        // Construct from components
        inline coupledInfo
        (
            const MeshType& mesh,
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
        inline const MeshType& baseMesh() const;

        // Set a new subMesh
        inline void setMesh(label index, MeshType* mesh);

        // Set the pointAlgorithm
        inline void setPointAlgorithm(pointSetAlgorithm* algorithm);

        // Set the faceAlgorithm
        inline void setFaceAlgorithm(faceSetAlgorithm* algorithm);

        // Set the cellAlgorithm
        inline void setCellAlgorithm(cellSetAlgorithm* algorithm);

        // Return a reference to the subMesh
        inline MeshType& subMesh();

        // Return a const reference to the subMesh
        inline const MeshType& subMesh() const;

        // Return if maps have been built
        inline bool builtMaps() const;

        // Set internal state of maps as built
        inline void setBuiltMaps();

        // Return a reference to the coupleMap
        inline coupleMap& map();

        // Return a const reference to the coupleMap
        inline const coupleMap& map() const;

        // Return a reference to the pointAlgorithm
        inline const pointSetAlgorithm& pointAlgorithm() const;

        // Return a reference to the faceAlgorithm
        inline const faceSetAlgorithm& faceAlgorithm() const;

        // Return a const reference to the cellAlgorithm
        inline const cellSetAlgorithm& cellAlgorithm() const;

        // Return the master face zone ID
        inline label masterFaceZone() const;

        // Return the slave face zone ID
        inline label slaveFaceZone() const;

    //- Interpolation

        // Subset point field
        template <class GeomField, class ZeroType>
        tmp<GeomField>
        subSetPointField
        (
            const GeomField& f,
            const ZeroType& zeroValue,
            const labelList& mapper
        ) const;

        // Subset geometric field
        template <class GeomField, class ZeroType>
        tmp<GeomField>
        subSetField
        (
            const GeomField& f,
            const ZeroType& zeroValue,
            const labelList& mapper
        ) const;

        // Subset point fields from registry to output stream
        template <class GeomField, class ZeroType>
        void sendPointFields
        (
            const wordList& fieldNames,
            const word& fieldType,
            const ZeroType& zeroValue,
            const labelList& mapper,
            OSstream& strStream
        ) const;

        // Subset geometric fields from registry to output stream
        template <class GeomField, class ZeroType>
        void sendFields
        (
            const wordList& fieldNames,
            const word& fieldType,
            const ZeroType& zeroValue,
            const labelList& mapper,
            OSstream& strStream
        ) const;

        // Set point field pointers from input dictionary
        template <class GeomField>
        void setPointField
        (
            const wordList& fieldNames,
            const dictionary& fieldDicts,
            const label internalSize,
            PtrList<GeomField>& fields
        ) const;

        // Set geometric field pointers from input dictionary
        template <class GeomField>
        void setField
        (
            const wordList& fieldNames,
            const dictionary& fieldDicts,
            const label internalSize,
            PtrList<GeomField>& fields
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
