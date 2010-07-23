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
#include "point.H"
#include "label.H"
#include "scalar.H"
#include "cellList.H"
#include "edgeList.H"
#include "faceList.H"
#include "vectorField.H"
#include "labelHashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                        Namespace meshOps Declaration
\*---------------------------------------------------------------------------*/

namespace meshOps
{

    // Return the magnitude of an edge, using UList
    inline scalar edgeLength
    (
        const edge& edgeToCheck,
        const UList<point>& points
    );

    // Compute the centre for a given face, using UList
    inline vector faceCentre
    (
        const face& faceToCheck,
        const UList<point>& points
    );

    // Compute the normal for a given face, using UList
    inline vector faceNormal
    (
        const face& faceToCheck,
        const UList<point>& points
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
    inline void cellCentreAndVolume
    (
        const label cIndex,
        const UList<point>& points,
        const UList<face>& faces,
        const UList<cell>& cells,
        const UList<label>& owner,
        vector& centre,
        scalar& volume
    );

    // Utility member function for separating axis tests
    inline label whichSide
    (
        const labelList& cellPoints,
        const UList<point>& points,
        const point& dir,
        const point& p
    );

    // Determine whether a point lies on a segment
    inline bool pointSegmentIntersection
    (
        const edge& segment,
        const UList<point>& points,
        const point& checkPoint
    );

    // Determine whether a two line segments intersect
    inline bool segmentSegmentIntersection
    (
        const edge& fromSegment,
        const edge& toSegment,
        const UList<point>& points,
        point& intPoint
    );

    // Determine whether a segment intersects a face
    inline bool segmentFaceIntersection
    (
        const edge& edgeToCheck,
        const face& faceToCheck,
        const UList<point>& points,
        vector& intersectionPoint
    );

    // Determine whether the particular point lies
    // inside the given face
    inline bool pointInFace
    (
        const face& faceToCheck,
        const UList<point>& points,
        const point& checkPoint
    );

    // Determine whether the particular point lies
    // inside the given cell
    inline bool pointInCell
    (
        const label cIndex,
        const cell& cellToCheck,
        const UList<face>& faces,
        const UList<label>& owner,
        const UList<point>& points,
        const point& checkPoint,
        const scalar& tolerance
    );

    // Compute the area / centre of a polygon
    // formed by a convex set of points.
    void convexSetArea
    (
        const label newFaceIndex,
        const label oldFaceIndex,
        const vectorField& cvxSet,
        const vector& refNorm,
        scalar& fArea,
        vector& fCentre,
        bool output = false
    );

    // Compute the volume / centre of a polyhedron
    // formed by a convex set of points.
    void convexSetVolume
    (
        const label newCellIndex,
        const label oldCellIndex,
        const vectorField& cvxSet,
        scalar& cVolume,
        vector& cCentre,
        bool output = false
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
    inline void insertPointLabels
    (
        const vector& refNorm,
        const vectorField& points,
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
    bool checkPointNearness
    (
        const pointField& points,
        const scalar magSqrTol
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

    // Actual routine to write out the VTK file
    void writeVTK
    (
        const polyMesh& mesh,
        const word& name,
        const label nPoints,
        const label nCells,
        const label nTotalCells,
        const pointField& points,
        const labelListList& cpList = labelListList(0),
        const label primitiveType = 3,
        const Map<label>& reversePointMap = Map<label>(),
        const Map<label>& reverseCellMap = Map<label>()
    );

} // End namespace meshOps


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshOpsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
