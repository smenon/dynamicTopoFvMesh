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
    meshOps

Description
    Various utility functions that perform mesh-related operations.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    meshOpsI.H
    meshOps.C

\*---------------------------------------------------------------------------*/

#ifndef meshOps_H
#define meshOps_H

#include "Map.H"
#include "line.H"
#include "point.H"
#include "label.H"
#include "scalar.H"
#include "cellList.H"
#include "edgeList.H"
#include "faceList.H"
#include "triangle.H"
#include "objectMap.H"
#include "vectorField.H"
#include "labelHashSet.H"

// Use the MPFR library, if necessary
#define USE_MPFR 0

#if USE_MPFR
#include "multiPrecision.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                        Namespace meshOps Declaration
\*---------------------------------------------------------------------------*/

namespace meshOps
{
    // Compute the centre for a given face, using UList
    template <class T1, class T2>
    inline void faceCentre
    (
        const face& faceToCheck,
        const UList<Vector<T1> >& points,
        Vector<T2>& centre
    );

    // Compute the normal for a given face, using UList
    template <class T1, class T2>
    inline void faceNormal
    (
        const face& faceToCheck,
        const UList<Vector<T1> >& points,
        Vector<T2>& normal
    );

    // Method to find the common-edge between two faces.
    inline bool findCommonEdge
    (
        const label first,
        const label second,
        const UList<labelList>& faceEdges,
        label& common
    );

    // Method to find the interior/boundary faces
    // for an input quad-face and adjacent triangle-prism cell.
    inline void findPrismFaces
    (
        const label fIndex,
        const label cIndex,
        const UList<face>& faces,
        const UList<cell>& cells,
        const UList<label>& neighbour,
        FixedList<face,2>& bdyf,
        FixedList<label,2>& bidx,
        FixedList<face,2>& intf,
        FixedList<label,2>& iidx
    );

    // Utility method to build a hull of cells
    // connected to the edge (for 2D simplical meshes)
    void constructPrismHull
    (
        const label eIndex,
        const UList<face>& faces,
        const UList<cell>& cells,
        const UList<label>& owner,
        const UList<label>& neighbour,
        const UList<labelList>& edgeFaces,
        labelList& hullTriFaces,
        labelList& hullCells
    );

    // Utility method to build a hull of cells (and faces)
    // around an edge (for 3D simplical meshes)
    void constructHull
    (
        const label eIndex,
        const UList<face>& faces,
        const UList<edge>& edges,
        const UList<cell>& cells,
        const UList<label>& owner,
        const UList<label>& neighbour,
        const UList<labelList>& faceEdges,
        const UList<labelList>& edgeFaces,
        const UList<labelList>& edgePoints,
        labelList& hullEdges,
        labelList& hullFaces,
        labelList& hullCells,
        labelListList& ringEntities
    );

    // Utility method to find the isolated point
    // given two triangular faces.
    inline label findIsolatedPoint
    (
        const face& baseFace,
        const face& checkFace
    );

    // Method to find the isolated point on a triangular face
    // that doesn't lie on the specified edge
    inline void findIsolatedPoint
    (
        const face& f,
        const edge& e,
        label& ptIndex,
        label& nextPtIndex
    );

    // Given a face and cell index, find the apex point for a tet cell.
    inline label tetApexPoint
    (
        const label cIndex,
        const label fIndex,
        const UList<face> faces,
        const UList<cell> cells
    );

    // Given a cell index, find the centroid / volume of a cell.
    template <class T1, class T2>
    inline void cellCentreAndVolume
    (
        const label cIndex,
        const UList<Vector<T1> >& points,
        const UList<face>& faces,
        const UList<cell>& cells,
        const UList<label>& owner,
        Vector<T2>& centre,
        T2& volume
    );

    // Determine whether a point lies on a segment
    template <class T>
    inline bool pointSegmentIntersection
    (
        const line<Vector<T>, const Vector<T>&>& segment,
        const Vector<T>& cP,
        const T& matchTol
    );

    // Determine whether a two line segments intersect
    template <class T>
    inline bool segmentSegmentIntersection
    (
        const line<Vector<T>, const Vector<T>&>& oldSegment,
        const line<Vector<T>, const Vector<T>&>& newSegment,
        const T& matchTol,
        Vector<T>& intPoint
    );

    // Determine whether a segment intersects a triangular face
    template <class T>
    inline bool segmentTriFaceIntersection
    (
        const triangle<Vector<T>, const Vector<T>&>& faceToCheck,
        const line<Vector<T>, const Vector<T>&>& edgeToCheck,
        const T& matchTol,
        Vector<T>& intPoint
    );

    // Determine whether the particular point lies
    // inside the given triangular face
    template <class T>
    inline bool pointInTriFace
    (
        const triangle<Vector<T>, const Vector<T>&>& faceToCheck,
        const Vector<T>& checkPoint
    );

    // Determine whether the particular point lies
    // inside the given cell
    template <class T>
    inline bool pointInCell
    (
        const label cIndex,
        const cell& cellToCheck,
        const UList<face>& faces,
        const UList<label>& owner,
        const UList<point>& points,
        const T& matchTol,
        const Vector<T>& cP
    );

    // Dijkstra's algorithm for the shortest path problem
    bool Dijkstra
    (
        const Map<point>& points,
        const Map<edge>& edges,
        const label startPoint,
        const label endPoint,
        Map<label>& pi
    );

    // Method to insert labels in a face, so that
    // right-handedness is preserved.
    template <class T>
    inline void insertPointLabels
    (
        const Vector<T>& refNorm,
        const Field<Vector<T> >& points,
        const labelHashSet& pLabels,
        face& modFace
    );

    // Method to insert a label between two labels in a list
    inline void insertLabel
    (
        const label newLabel,
        const label labelA,
        const label labelB,
        labelList& list
    );

    // Utility method to replace a label in a given list
    inline void replaceLabel
    (
         const label original,
         const label replacement,
         labelList& list
    );

    // Utility method to size-up the list to include an item
    template <class Type>
    inline void sizeUpList
    (
        const Type item,
        List<Type>& list
    );

    // Utility method to size-down the list to remove an item
    template <class Type>
    inline void sizeDownList
    (
        const Type item,
        List<Type>& list
    );

    // Remove an item at a particular index in the list
    template <class Type>
    inline void removeIndex
    (
        const label index,
        List<Type>& list
    );

    // Utility method to check for concurrent points.
    template <class T>
    bool checkPointNearness
    (
        const Field<Vector<T> >& points,
        const T& magSqrTol
    );

    // Parallel send
    void pWrite
    (
        const label toID,
        const label& data
    );

    // Parallel send (for fixed lists)
    template <class Type, label Size>
    void pWrite
    (
        const label toID,
        const FixedList<Type,Size>& data
    );

    // Parallel send (for lists)
    template <class Type>
    void pWrite
    (
        const label toID,
        const UList<Type>& data
    );

    // Parallel receive
    void pRead
    (
        const label fromID,
        label& data
    );

    // Parallel receive (for fixed lists)
    template <class Type, label Size>
    void pRead
    (
        const label fromID,
        FixedList<Type,Size>& data
    );

    // Parallel receive (for lists)
    template <class Type>
    void pRead
    (
        const label fromID,
        UList<Type>& data
    );

    // Wait for buffer transfer completion.
    void waitForBuffers();

    // Select a list of elements from connectivity,
    // and output to a VTK format
    void writeVTK
    (
        const polyMesh& mesh,
        const word& name,
        const labelList& cList,
        const label primitiveType,
        const UList<point>& meshPoints,
        const UList<edge>& edges,
        const UList<face>& faces,
        const UList<cell>& cells,
        const UList<label>& owner
    );

    // Actual routine to write out the VTK file
    void writeVTK
    (
        const polyMesh& mesh,
        const word& name,
        const label nPoints,
        const label nCells,
        const label nTotalCells,
        const vectorField& points,
        const labelListList& cpList = labelListList(0),
        const label primitiveType = 3,
        const Map<label>& reversePointMap = Map<label>(),
        const Map<label>& reverseCellMap = Map<label>()
    );

} // End namespace meshOps


/*---------------------------------------------------------------------------*\
                    Class convexSetAlgorithm Declaration
\*---------------------------------------------------------------------------*/

class convexSetAlgorithm
{

protected:

    // Protected data

        const bool twoDMesh_;
        const label nOldPoints_;

        //- References to old-level connectivity
        //  [Before topo-changes, at old point-positions]
        const polyMesh& mesh_;

        //- References to new-level connectivity
        //  [After topo-changes, at old point-positions]
        const UList<point>& newPoints_;
        const UList<edge>& newEdges_;
        const UList<face>& newFaces_;
        const UList<cell>& newCells_;
        const UList<label>& newOwner_;
        const UList<label>& newNeighbour_;

        //- Reference to pointsFromPoints information
        const List<objectMap>& pointsFromPoints_;

        //- Reference to modified-points information
        const Map<labelList>& modPoints_;

        //- Switch for higher precision
        mutable bool highPrecision_;

        //- Entity parents
        mutable labelList parents_;

        //- Low-precision fields
        mutable vector refNorm_;
        mutable scalar normFactor_;
        mutable vectorField centres_;
        mutable scalarField weights_;

#       if USE_MPFR

        //- High-precision fields
        mutable mpVector mpRefNorm_;
        mutable mpScalar mpNormFactor_;
        mutable Field<mpVector> mpCentres_;
        mutable Field<mpScalar> mpWeights_;

#       endif

    //- Output connectivity to VTK

        // Output an entity as a VTK file
        void writeVTK
        (
            const word& name,
            const label entity,
            const label primitiveType = 3,
            const bool useOldConnectivity = false
        ) const;

        // Output a list of entities as a VTK file
        void writeVTK
        (
            const word& name,
            const labelList& cList,
            const label primitiveType = 3,
            const bool useOldConnectivity = false
        ) const;

public:

    //- Constructor

        // Construct from components
        convexSetAlgorithm
        (
            const polyMesh& mesh,
            const UList<point>& newPoints,
            const UList<edge>& newEdges,
            const UList<face>& newFaces,
            const UList<cell>& newCells,
            const UList<label>& newOwner,
            const UList<label>& newNeighbour,
            const List<objectMap>& pointsFromPoints,
            const Map<labelList>& modPoints
        );

    //- Destructor

        virtual ~convexSetAlgorithm()
        {}

    //- Member functions

        // Dimensions of the algorithm
        virtual label dimension() const = 0;

        // Return true if accumulated weights are consistent
        virtual bool consistent(const scalar tolerance) const;

        // Return the normFactor
        virtual scalar normFactor() const;

        // Normalize stored weights
        virtual void normalize(bool normSum) const;

        // Extract weights and centres to lists
        virtual void populateLists
        (
            labelList& parents,
            vectorField& centres,
            scalarField& weights
        ) const;

        // Return the highPrecision switch
        bool highPrecision() const
        {
            return highPrecision_;
        }

        // Set / Reset higher precision mode
        void setHighPrecision() const
        {
            highPrecision_ = true;
        }

        void unsetHighPrecision() const
        {
            highPrecision_ = false;
        }

        // Compute normFactor
        virtual void computeNormFactor(const label index) const = 0;

        // Compute intersections
        virtual bool computeInsersection
        (
            const label newIndex,
            const label oldIndex,
            const scalar& matchTol,
            bool output
        ) const = 0;
};


/*---------------------------------------------------------------------------*\
                     Class faceSetAlgorithm Declaration
\*---------------------------------------------------------------------------*/

class faceSetAlgorithm
:
    public convexSetAlgorithm
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        faceSetAlgorithm(const faceSetAlgorithm&);

        //- Disallow default bitwise assignment
        void operator=(const faceSetAlgorithm&);

        // Compute the area / centre of a polygon
        // formed by a convex set of points.
        template <class T>
        void convexSetArea
        (
            const label newIndex,
            const label oldIndex,
            const Field<Vector<T> >& cvxSet,
            const Vector<T>& refNorm,
            T& fArea,
            Vector<T>& fCentre,
            bool output
        ) const;

        //- Intersection member function
        //  Templated for precision
        template <class T>
        bool faceIntersection
        (
            const label newIndex,
            const label oldIndex,
            const T& matchTol,
            Field<Vector<T> >& intPoints,
            bool output
        ) const;

public:

    //- Constructor

        faceSetAlgorithm
        (
            const polyMesh& mesh,
            const UList<point>& newPoints,
            const UList<edge>& newEdges,
            const UList<face>& newFaces,
            const UList<cell>& newCells,
            const UList<label>& newOwner,
            const UList<label>& newNeighbour,
            const List<objectMap>& pointsFromPoints,
            const Map<labelList>& modPoints
        );

    //- Destructor

        virtual ~faceSetAlgorithm()
        {}

    //- Virtual member functions

        // Dimensions of the algorithm
        virtual label dimension() const
        {
            return 2;
        }

        // Compute normFactor
        virtual void computeNormFactor(const label index) const;

        // Compute intersections
        virtual bool computeInsersection
        (
            const label newIndex,
            const label oldIndex,
            const scalar& matchTol,
            bool output
        ) const;
};


/*---------------------------------------------------------------------------*\
                     Class cellSetAlgorithm Declaration
\*---------------------------------------------------------------------------*/

class cellSetAlgorithm
:
    public convexSetAlgorithm
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        cellSetAlgorithm(const cellSetAlgorithm&);

        //- Disallow default bitwise assignment
        void operator=(const cellSetAlgorithm&);

        // Compute the volume / centre of a polyhedron
        // formed by a convex set of points.
        template <class T>
        void convexSetVolume
        (
            const label newCellIndex,
            const label oldCellIndex,
            const Field<Vector<T> >& cvxSet,
            T& cVolume,
            Vector<T>& cCentre,
            bool output
        ) const;

        //- Intersection member function
        //  Templated for precision
        template <class T>
        bool cellIntersection
        (
            const label newIndex,
            const label oldIndex,
            const T& matchTol,
            Field<Vector<T> >& intPoints,
            bool output
        ) const;

public:

    //- Constructor

        cellSetAlgorithm
        (
            const polyMesh& mesh,
            const UList<point>& newPoints,
            const UList<edge>& newEdges,
            const UList<face>& newFaces,
            const UList<cell>& newCells,
            const UList<label>& newOwner,
            const UList<label>& newNeighbour,
            const List<objectMap>& pointsFromPoints,
            const Map<labelList>& modPoints
        );

    //- Destructor

        virtual ~cellSetAlgorithm()
        {}

    //- Virtual member functions

        // Dimensions of the algorithm
        virtual label dimension() const
        {
            return 3;
        }

        // Compute normFactor
        virtual void computeNormFactor(const label index) const;

        // Compute intersections
        virtual bool computeInsersection
        (
            const label newIndex,
            const label oldIndex,
            const scalar& matchTol,
            bool output
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshOpsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
