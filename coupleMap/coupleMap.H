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

\*---------------------------------------------------------------------------*/

#ifndef coupleMap_H
#define coupleMap_H

#include "regIOobject.H"
#include "labelList.H"
#include "pointField.H"
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
    mutable bool isLocal_;

    mutable label masterIndex_;
    mutable label slaveIndex_;

    mutable pointField pointBuffer_;

    // List of points that are required
    // during patch sub-mesh creation
    mutable labelList subMeshPoints_;

    // Entity sizes (as specified by the entityType enumerant)
    mutable FixedList<label,6> nEntities_;

    // Maps for entities
    mutable FixedList<Map<label>,4> entityMap_;
    mutable FixedList<Map<label>,4> reverseEntityMap_;

    // Entity Buffers (as specified by the entityType enumerant)
    mutable FixedList<labelList,5> entityBuffer_;

    //- Demand-driven connectivity data.
        mutable label nInternalFaces_;
        mutable labelList* ownerPtr_;
        mutable labelList* neighbourPtr_;

        void makeAddressing() const;

    //- Disallow default bitwise assignment
    void operator=(const coupleMap&);

public:

    // Enumerants for entity-types
    enum entityType
    {
        POINT = 0,
        EDGE = 1,
        FACE = 2,
        CELL = 3,
        INTERNAL_EDGE = 4,
        FACE_EDGE = 4,
        SHARED_POINT = 5
    };

    //- Runtime type information
    TypeName("coupleMap");

    // Constructor
    coupleMap
    (
        const IOobject& io,
        const bool isLocal,
        const label masterIndex,
        const label slaveIndex
    );

    // Construct as copy
    coupleMap(const coupleMap&);

    // Destructor

        virtual ~coupleMap();

    //- Access

        label masterIndex() const;
        label slaveIndex() const;

        bool isLocal() const;

        pointField& pointBuffer() const;

        labelList& subMeshPoints() const;

        label findSlaveIndex
        (
            const label eType,
            const label Index
        ) const;

        label findMasterIndex
        (
            const label eType,
            const label Index
        ) const;

        void removeSlaveIndex
        (
            const label eType,
            const label Index
        ) const;

        void removeMasterIndex
        (
            const label eType,
            const label Index
        ) const;

        void mapSlave
        (
            const label eType,
            const label master,
            const label slave
        ) const;

        void mapMaster
        (
            const label eType,
            const label slave,
            const label master
        ) const;

        FixedList<label,6>& nEntities() const;
        label& nEntities(const label eType) const;

        Map<label>& entityMap(const label eType) const;
        Map<label>& reverseEntityMap(const label eType) const;

        FixedList<labelList,5>& entityBuffer() const;
        labelList& entityBuffer(const label eType) const;

        label nInternalFaces() const;
        const labelList& owner() const;
        const labelList& neighbour() const;

        void transferMaps
        (
            const label eType,
            Map<label>& newEntityMap,
            Map<label>& newReverseEntityMap
        ) const;

        void clearMaps() const;

        bool writeData(Ostream&) const
        {
            return true;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
