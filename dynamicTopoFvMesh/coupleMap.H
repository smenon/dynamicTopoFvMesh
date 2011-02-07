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
    coupleMap

Description
    Coupled patch information with registry support.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    coupleMapI.H
    coupleMap.C

\*---------------------------------------------------------------------------*/

#ifndef coupleMap_H
#define coupleMap_H

#include "regIOobject.H"
#include "labelList.H"
#include "pointField.H"
#include "fieldTypes.H"
#include "edgeList.H"
#include "faceList.H"
#include "cellList.H"
#include "labelPair.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class coupleMap Declaration
\*---------------------------------------------------------------------------*/

class coupleMap
:
    public regIOobject
{
    // Private data
    mutable bool twoDMesh_;

    mutable bool isLocal_;
    mutable bool isSend_;

    mutable label masterIndex_;
    mutable label slaveIndex_;

    mutable pointField pointBuffer_;
    mutable pointField oldPointBuffer_;

    // List of points that are required
    // during patch sub-mesh creation
    mutable labelList subMeshPoints_;
    mutable List<labelPair> globalProcPoints_;

    // Entity sizes (as specified by the entityType enumerant)
    mutable FixedList<label,10> nEntities_;

    // Maps for entities
    mutable FixedList<Map<label>,4> entityMap_;
    mutable FixedList<Map<label>,4> reverseEntityMap_;

    // Entity Buffers (as specified by the entityType enumerant)
    mutable FixedList<labelList,11> entityBuffer_;

    // List of entity indices with topological operations
    mutable labelList entityIndices_;

    // List of operations performed on entities
    mutable labelList entityOperations_;

    // List of point-locations to move points to
    mutable pointField moveNewPoints_;
    mutable pointField moveOldPoints_;

    //- Demand-driven connectivity data.
        mutable label nInternalFaces_;
        mutable labelList* ownerPtr_;
        mutable labelList* neighbourPtr_;
        mutable edgeList* edgesPtr_;
        mutable faceList* facesPtr_;
        mutable cellList* cellsPtr_;
        mutable labelListList* faceEdgesPtr_;

        inline void makeEdges() const;
        inline void makeFaces() const;
        inline void makeCells() const;

        inline void makeAddressing() const;

        inline void clearAddressing() const;

    //- Disallow default bitwise assignment
        inline void operator=(const coupleMap&);

public:

    // Enumerants for entity-types
    enum entityType
    {
        // Buffers and sizes
        POINT = 0,
        EDGE = 1,
        FACE = 2,
        CELL = 3,
        // Sizes only
        INTERNAL_EDGE = 4,
        INTERNAL_FACE = 5,
        SHARED_POINT = 6,
        GLOBAL_POINT = 7,
        NFE_SIZE = 8,
        NBDY = 9,
        // Buffers only
        FACE_EDGE = 4,
        NFE_BUFFER = 5,
        FACE_STARTS = 6,
        FACE_SIZES = 7,
        EDGE_STARTS = 8,
        EDGE_SIZES = 9,
        PATCH_ID = 10
    };

    enum operationType
    {
        BISECTION = 0,
        COLLAPSE_FIRST = 1,
        COLLAPSE_SECOND = 2,
        COLLAPSE_MIDPOINT = 3,
        REMOVE_CELL = 4,
        MOVE_POINT = 5,
        CONVERT_PATCH = 6
    };

    //- Runtime type information
    TypeName("coupleMap");

    // Constructor
    inline coupleMap
    (
        const IOobject& io,
        const bool twoDMesh,
        const bool isLocal,
        const bool isSend,
        const label masterIndex,
        const label slaveIndex
    );

    // Construct as copy
    inline coupleMap(const coupleMap&);

    // Destructor

        inline virtual ~coupleMap();

    // Interpolation functions

        //- Interpolate point field
        template<class Type>
        tmp<Field<Type> > pointInterpolate
        (
            const Map<label>& mPointMap,
            const Map<label>& sPointMap,
            const Field<Type>& pf,
            bool reverse = false
        ) const;

        template<class Type>
        tmp<Field<Type> > pointInterpolate
        (
            const Map<label>& mPointMap,
            const Map<label>& sPointMap,
            const tmp<Field<Type> >& tpf,
            bool reverse = false
        ) const;

        //- Interpolate face field
        template<class Type>
        tmp<Field<Type> > faceInterpolate
        (
            const label mStart,
            const label sStart,
            const Field<Type>& pf,
            bool reverse = false
        ) const;

        template<class Type>
        tmp<Field<Type> > faceInterpolate
        (
            const label mStart,
            const label sStart,
            const tmp<Field<Type> >& tpf,
            bool reverse = false
        ) const;

    //- Access

        inline label masterIndex() const;
        inline label slaveIndex() const;

        inline bool isTwoDMesh() const;
        inline bool isLocal() const;
        inline bool isProcessor() const;
        inline bool isSend() const;
        inline bool isRecv() const;

        inline pointField& pointBuffer() const;
        inline pointField& oldPointBuffer() const;

        inline labelList& subMeshPoints() const;
        inline List<labelPair>& globalProcPoints() const;

        inline void allocateBuffers() const;

        inline label findSlave
        (
            const label eType,
            const label Index
        ) const;

        inline label findMaster
        (
            const label eType,
            const label Index
        ) const;

        inline void removeSlave
        (
            const label eType,
            const label Index
        ) const;

        inline void removeMaster
        (
            const label eType,
            const label Index
        ) const;

        inline void mapSlave
        (
            const label eType,
            const label master,
            const label slave
        ) const;

        inline void mapMaster
        (
            const label eType,
            const label slave,
            const label master
        ) const;

        inline FixedList<label,10>& nEntities() const;
        inline label& nEntities(const label eType) const;

        inline Map<label>& entityMap(const label eType) const;
        inline Map<label>& reverseEntityMap(const label eType) const;

        inline FixedList<labelList,11>& entityBuffer() const;
        inline labelList& entityBuffer(const label eType) const;

        inline labelList& entityIndices() const;
        inline labelList& entityOperations() const;

        inline pointField& moveNewPoints() const;
        inline pointField& moveOldPoints() const;

        inline void pushOperation
        (
            const label index,
            const label opType,
            const point& newPoint = vector::zero,
            const point& oldPoint = vector::zero
        ) const;

        //- Demand-driven data
        inline label nInternalFaces() const;
        inline const labelList& owner() const;
        inline const labelList& neighbour() const;
        inline const edgeList& edges() const;
        inline const faceList& faces() const;
        inline const cellList& cells() const;
        inline const labelListList& faceEdges() const;

        inline void transferMaps
        (
            const label eType,
            Map<label>& newEntityMap,
            Map<label>& newReverseEntityMap
        ) const;

        inline void clearMaps() const;

        inline void clearBuffers() const;

        inline bool readData(Istream&);

        inline bool writeData(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "coupleMapI.H"

#ifdef NoRepository
#    include "coupleMap.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
