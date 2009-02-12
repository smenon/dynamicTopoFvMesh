/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    dynamicTopoFvMesh

Description
    An implementation of dynamic changes to mesh-topology (w/ smoothing)

SourceFiles
    dynamicTopoFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicTopoFvMesh_H
#define dynamicTopoFvMesh_H

#include <sys/time.h>
#include "HashList.H"
#include "multiThreader.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class mapPolyMesh;
class objectMap;
class motionSolver;

/*---------------------------------------------------------------------------*\
                      Class dynamicTopoFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicTopoFvMesh
:
    public fvMesh
{
    // Private data

        //- Number of boundary-patches
        label numPatches_;

        //- Topology change flag
        bool topoChangeFlag_;

        //- Dynamic mesh dictionary
        IOdictionary dict_;

        //- Mesh characteristics [2D/3D]
        Switch twoDMesh_;

        //- Edge bisection/collapse option
        Switch edgeModification_;

        //- Sub-dictionary which specifies fixed length-scales for patches
        dictionary fixedLengthScalePatches_;

        //- List of noSwap patch IDs (which specifies patches on
        //  which 2-2 swaps are to be avoided)
        labelHashSet noSwapPatchIDs_;

        //- Should the mesh solve for motion on update?
        Switch solveForMotion_;

        //- The motion solver pointer
        autoPtr<motionSolver> motionPtr_;

        //- PolyMesh mapping structure
        autoPtr<mapPolyMesh> mapper_;

        //- Local copies of connectivity information
        HashList<point> meshPoints_;
        HashList<face>  faces_;
        HashList<label> owner_;
        HashList<label> neighbour_;
        HashList<cell>  cells_;
        HashList<edge>  edges_;
        HashList<labelList> pointEdges_;
        HashList<labelList> edgeFaces_;
        HashList<labelList> faceEdges_;

        //- Up-to-date edge connectivity information (with IO support)
        IOList<edge>        IOedges_;
        IOList<labelList>   IOpointEdges_;
        IOList<labelList>   IOedgeFaces_;
        IOList<labelList>   IOfaceEdges_;
        IOList<label>       IOedgePatchStarts_;
        IOList<label>       IOedgePatchSizes_;

        //- Size information
        labelList       oldPatchSizes_, patchSizes_;
        labelList       oldPatchStarts_, patchStarts_;
        labelList       oldEdgePatchSizes_, edgePatchSizes_;
        labelList       oldEdgePatchStarts_, edgePatchStarts_;
        labelList       oldPatchNMeshPoints_, patchNMeshPoints_;
        label           nOldPoints_, nPoints_;
        label           nOldEdges_,  nEdges_;
        label           nOldFaces_,  nFaces_;
        label           nOldCells_,  nCells_;
        label           nOldInternalFaces_, nInternalFaces_;
        label           nOldInternalEdges_, nInternalEdges_;

        //- Old-to-New
        labelList       reversePointMap_;
        labelList       reverseEdgeMap_;
        labelList       reverseFaceMap_;
        labelList       reverseCellMap_;

        //- New-to-Old
        labelList       pointMap_;
        labelList       edgeMap_;
        labelList       faceMap_;
        labelList       cellMap_;

        //- Maps for the renumbering of added entities
        Map<objectMap>  cellsFromCells_;
        Map<label>      addedPointRenumbering_;
        Map<label>      addedEdgeRenumbering_;
        Map<label>      addedFaceRenumbering_;
        Map<label>      addedCellRenumbering_;
        Map<label>      addedFacePatches_;
        Map<label>      addedEdgePatches_;
        Map<label>      cellParents_;

        // Cell-centre field
        autoPtr<vectorField> cellCentresPtr_;

        //- Edge length cut-off criteria
        //- ( for all faces in the mesh [2D] )
        HashList<scalar>    lengthScale_;
        HashList<edge>      edgeToWatch_;
        scalar              ratioMin_;
        scalar              ratioMax_;
        scalar              growthFactor_;
        scalar              maxLengthScale_;
        scalar              sliverThreshold_;
        label               bisectInteriorFace_;

        //- Maximum number of tets around an edge
        //  for 3D mesh flipping
        label               maxTetsPerEdge_;
        Switch              allowTableResize_;

        // Private timer sub-class
        class timer
        {
            // Private data
            
            // Defined in time.h
            struct timeval start_;
            struct timeval end_;

        public:

            // Constructor
            timer(){}

            // Member functions

                //- Start the timer
                void start();

                //- Stop the timer
                void stop();

                //- Report time between start/stop
                scalar reportTime();
        };

        // Private Stack sub-class
        class stack
        {
            // Private data
            labelHashSet    stack_;

            Mutex           stackMutex_;

        public:

            // Constructor
            stack(){}

            // Member functions for access to the stack

                //- Initialize the stack with a prescribed size
                void initStack(const label size);

                //- Initialize the stack with a list
                template <class Type>
                void initStack(HashList<Type>& list);

                //- Push items on to the stack
                inline void push(const label index);

                //- Pop an item off the stack
                inline label pop();

                //- Remove a specific index off the stack
                inline void remove(const label index);

                //- Return if a stack is empty or not
                inline bool empty();

                //- Mutex lock (for multiThreading)
                inline void lock();

                //- Mutex unlock (for multiThreading)
                inline void unlock();

                //- Print out the stack
                inline void print();
        };

        // Stack of faces to be checked for topo-changes [2D]
        stack faceStack_;

        // Stack of edges to be checked for topo-changes [3D]
        stack edgeStack_;

        //- Support for multithreading
        autoPtr<multiThreader> threader_;

        //- Structure used for multiThreading
        struct topoMeshStruct
        {
            // Pointer to this class
            dynamicTopoFvMesh *mesh_;

            // Total number of threads and thread ID
            label nThreads_, threadID_;

            // Default constructor
            topoMeshStruct();
            
        } *structPtr_;

        //- Mutexes for all entities in the mesh
        Mutex           pMutex_;
        Mutex           eMutex_;
        Mutex           fMutex_;
        Mutex           cMutex_;

        //- Mutex lists for all entities in the mesh
        HashList<Mutex> pointMutex_;
        HashList<Mutex> edgeMutex_;
        HashList<Mutex> faceMutex_;
        HashList<Mutex> cellMutex_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicTopoFvMesh(const dynamicTopoFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicTopoFvMesh&);

        // Find the circumcenter, given three points
        inline vector circumCenter
        (
            const point& a,
            const point& b,
            const point& c,
            const label& one,
            const label& two,
            const label& three
        );

        // Typedef for tetrahedral metrics
        typedef scalar (*tetMetricReturnType)
        (
            const point& p0,
            const point& p1,
            const point& p2,
            const point& p3
        );

        //- Quality metric for tetrahedra in 3D
        tetMetricReturnType tetMetric_;

        // Find the area of a triangle face.
        inline scalar triFaceArea(const face& triFace);

        // Find the normal of a triangle face.
        inline vector triFaceNormal(const face& triFace);

        // Find the volume of a tetrahedron.
        inline scalar tetVolume
        (
            const point& a,
            const point& b,
            const point& c,
            const point& d
        );

        // Initialize the length-scale field
        void initLengthScale();

        // Initialize mesh edges and related connectivity lists
        void initEdges();

        // Initialize mutex lists
        void initMutexLists();

        // Unlock mutex lists with HashSets
        void unlockMutexLists
        (
            labelHashSet& pLocks,
            labelHashSet& eLocks,
            labelHashSet& fLocks,
            labelHashSet& cLocks
        );

        // Return a reference to the multiThreader
        const multiThreader& threader() const;

        // Return length-scale at an edge-location in the mesh [3D]
        inline scalar meshEdgeLengthScale(const label eIndex);

        // Check if a given face is a quad
        inline bool checkQuadFace(const label fIndex);

        // Return the length of an edge
        inline scalar edgeLength(const label eIndex);

        // Return cell-labels on either side of the face
        inline void faceCells
        (
            const label fIndex,
            FixedList<label,2> cellLabels
        );

        // Return length-scale at an cell-location in the mesh
        inline scalar meshCellLengthScale(const label cIndex);

        // MultiThreaded topology modifier [2D]
        void threadedTopoModifier2D();

        // MultiThreaded topology modifier [3D]
        void threadedTopoModifier3D();

        // 2D Edge-swapping engine
        static void swap2DEdges(void *argument);

        // 2D Edge-bisection/collapse engine
        static void edgeBisectCollapse2D(void *argument);

        // 3D Edge-swapping engine
        static void swap3DEdges(void *argument);

        // 3D Edge-bisection/collapse engine
        static void edgeBisectCollapse3D(void *argument);

        // Return the face-stack
        inline stack& faceStack();

        // Return the edge-stack
        inline stack& edgeStack();

        // Attempt to lock the edge and its two points.
        // Return zero on success.
        inline bool lockEdge
        (
            const label eIndex,
            const label threadID
        );

        // Method to determine the boundary patch index for a given face
        inline label whichPatch(const label& index) const;

        // Method to determine the boundary patch index for a given edge
        inline label whichEdgePatch(const label& index) const;

        // Method for the swapping of a quad-face in 2D
        void swapQuadFace
        (
            const label findex,
            const FixedList<label,4>& commonFaceIndex,
            const FixedList<face,4>&  commonFaces,
            const FixedList<edge,2>&  commonEdges
        );

        // Method for the bisection of a quad-face in 2D
        void bisectQuadFace
        (
            const label findex
        );

        // Method for the collapse of a quad-face in 2D
        // Returns a boolean value indicating whether the collapse was valid
        bool collapseQuadFace
        (
            const label findex
        );

        // Method for the bisection of an edge in 3D
        void bisectEdge
        (
            const label eIndex,
            const label threadID
        );

        // Method for the collapse of an edge in 3D
        // Returns a boolean value indicating whether the collapse was valid
        bool collapseEdge
        (
            const label eIndex,
            const label threadID
        );

        // Utility to determine ring edges and the hull faces connected to them.
        void constructEdgeRing
        (
            const label eIndex,
            const edge& edgeToCheck,
            const labelList& hullVertices,
            const labelList& hullFaces,
            const labelList& hullCells,
            labelList& ringEdges,
            labelListList& hullEdgesAndFaces,
            FixedList<bool,2>& edgeBoundary
        );

        // Utility method to check for invalid edge-collapse.
        bool checkCollapse
        (
            const point& newPoint,
            const label pointIndex,
            const label cellIndex,
            const FixedList<bool,2>& edgeBoundary,
            labelHashSet& cellsChecked
        );

        // Check if the boundary face is adjacent to a sliver-cell,
        // and remove it by a two-step bisection/collapse operation.
        bool remove2DSliver
        (
            const label findex
        );

        // Method to find the interior/boundary faces
        // for an input quad-face and adjacent triangle-prism cell.
        inline void findPrismFaces
        (
            const label findex,
            const label cindex,
            FixedList<face,2>& bdyf,
            FixedList<label,2>& bidx,
            FixedList<face,2>& intf,
            FixedList<label,2>& iidx
        );

        // Method to find the common-edge between two faces.
        // If found, returns the common edge on the first face in the argument
        bool findCommonEdge
        (
            const face& first,
            const face& second,
            edge& common
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

        // Method to replace a label in a given list
        inline void replaceLabel
        (
            const label& original,
            const label& replacement,
            labelList& list
        );

        // Insert the specified face to the mesh
        label insertFace
        (
            const label patch,
            const face& newFace,
            const label newOwner,
            const label newNeighbour,
            const edge& edgeToWatch
        );

        // Remove the specified face from the mesh
        void removeFace(const label index);

        // Insert the specified edge to the mesh
        label insertEdge
        (
            const label patch,
            const edge& newEdge,
            const labelList& edgeFaces
        );

        // Remove the specified edge from the mesh
        void removeEdge(const label index);

        // Utility method to size-up the list to include an item
        inline void sizeUpList
        (
            const label item,
            labelList& list
        );

        // Utility method to size-down the list to remove an item
        inline void sizeDownList
        (
            const label item,
            labelList& list
        );

        // Utility method to build a hull of faces connected to the edge
        // This will also determine whether the edge lies on a boundary
        bool constructPrismHull
        (
            const edge& edgeToCheck,
            const label startFaceIndex,
            DynamicList<label>& hullFaces,
            DynamicList<label>& hullCells,
            DynamicList<label>& hullTriFaces,
            DynamicList<label>& edgePatches,
            bool requiresTriFaces
        );

        // Utility method to build a ring of vertices around the edge
        inline bool constructVertexRing
        (
            const label eIndex,
            const label threadID,
            DynamicList<label>& hullCells,
            DynamicList<label>& hullFaces,
            DynamicList<label>& hullVertices,
            scalar& minQuality,
            labelHashSet& pLocks,
            labelHashSet& eLocks,
            labelHashSet& fLocks,
            labelHashSet& cLocks,
            bool requiresQuality
        );

        // Check whether the given edge is on a bounding curve
        inline bool checkBoundingCurve(label eIndex);

        // Allocate dynamic programming tables
        inline void initTables
        (
            const label mMax,
            scalarListList& Q,
            labelListList& K,
            labelListList& triangulations
        );

        // Utility method to fill the dynamic programming tables
        // Returns the number of triangulations
        inline label fillTables
        (
            const label eIndex,
            const DynamicList<label>& hullVertices,
            const scalar minQuality,
            scalarListList& Q,
            labelListList& K
        );

        // Remove the edge according to the swap sequence
        void removeEdgeFlips
        (
            const label m,
            const label eIndex,
            const labelListList& K,
            const DynamicList<label>& hullCells,
            const DynamicList<label>& hullFaces,
            const DynamicList<label>& hullVertices,
            labelListList& triangulations
        );

        // Extract triangulations from the programming table
        void extractTriangulation
        (
            const label i,
            const label j,
            const labelListList& K,
            label& numTriangulations,
            labelListList& triangulations
        );

        // Identify the 3-2 swap from the triangulation sequence
        label identify32Swap
        (
            const label m,
            const edge& edgeToCheck,
            const DynamicList<label>& hullVertices,
            const labelListList& triangulations
        );

        // Routine to check whether the triangulation at the
        // index lies on the boundary of the vertex ring.
        bool boundaryTriangulation
        (
            const label index,
            label& isolatedVertex,
            labelListList& triangulations
        );

        // Routine to perform 2-3 swaps
        void swap23
        (
            const label isolatedVertex,
            const label edgeToCheckIndex,
            const edge& edgeToCheck,
            const label triangulationIndex,
            const labelListList& triangulations,
            const DynamicList<label>& hullCells,
            const DynamicList<label>& hullFaces,
            const DynamicList<label>& hullVertices
        );

        // Routine to perform 3-2 swaps
        void swap32
        (
            const label edgeToCheckIndex,
            const edge& edgeToCheck,
            const label triangulationIndex,
            const labelListList& triangulations,
            const DynamicList<label>& hullCells,
            const DynamicList<label>& hullFaces,
            const DynamicList<label>& hullVertices
        );

        // Reorder points after a topology change
        void reOrderPoints(pointField& points);

        // Reorder edges after a topology change
        void reOrderEdges();

        // Reorder faces in upper-triangular order after a topology change
        void reOrderFaces
        (
            faceList& faces,
            labelList& owner,
            labelList& neighbour
        );

        // Reorder & renumber cells with bandwidth
        // reduction after a topology change
        void reOrderCells();

        // Reorder the mesh in upper-triangular order,
        // and generate mapping information
        void reOrderMesh
        (
            pointField& points,
            faceList& faces,
            labelList& owner,
            labelList& neighbour
        );

        // Copy edge-based connectivity from HashLists
        void setEdgeConnectivity();

        // Check the state of edge-based connectivity HashLists
        void checkEdgeConnectivity();

        // Check tet-specific connectivity structures
        void checkTetConnectivity();

        // Calculate the edge length-scale for the mesh
        void calculateLengthScale();

        // Return the appropriate length-scale for boundary face
        scalar boundaryLengthScale(const label faceIndex);

        // Does the mesh perform edge-modification?
        bool edgeModification();

        // Private internal access to member data for static members
        HashList<point>& meshPoints() { return meshPoints_; }
        HashList<edge>& edgeToWatch() { return edgeToWatch_; }
        bool& topoChangeFlag() { return topoChangeFlag_; }
        scalar ratioMin() { return ratioMin_; }
        scalar ratioMax() { return ratioMax_; }
        scalar growthFactor() { return growthFactor_; }
        label maxTetsPerEdge() { return maxTetsPerEdge_; }
        bool allowTableResize() { return allowTableResize_; }

public:

    // Declare the name of the class and its debug switch
        ClassName("dynamicTopoFvMesh");

    // Constructors

        //- Construct from existing IOobject
        explicit dynamicTopoFvMesh(const IOobject&);

    // Destructor

        virtual ~dynamicTopoFvMesh();

    // Member Functions

        // Return the number of edges in the mesh.
        // Override of primitiveMesh member function
        label nEdges() const;

        // Return the number of internal edges in the mesh.
        // Override of primitiveMesh member function
        label nInternalEdges() const;

        // Return the list of edges in the mesh.
        // Override of primitiveMesh member function.
        const edgeList& edges() const;

        // Update the mesh for topology changes
        // Return true if changes have occurred
        bool updateTopology();

        // Update the mesh for motion
        void updateMotion();

        //- Update mesh corresponding to the given map
        virtual void updateMesh(const mapPolyMesh& mpm);

        // Map all fields in time using the given map.
        virtual void mapFields(const mapPolyMesh& mpm);

        // Set mesh-motion displacement boundary conditions
        void setMotionBC
        (
            const label& index,
            const vectorField& dispField
        );

        // Return the mesh-mapper
        const mapPolyMesh& meshMap();

        // Return old cell-centre information (prior to a topology change)
        const vectorField& oldCellCentres() const;

        // Return mesh length-scale values
        tmp<volScalarField> lengthScale();

        // Return mesh cell-quality values
        tmp<volScalarField> meshQuality();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
