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
    dynamicTopoFvMesh

Description
    An implementation of dynamic changes to mesh-topology

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    dynamicTopoFvMesh.C
    dynamicTopoFvMeshI.H
    dynamicTopoFvMeshCheck.C
    dynamicTopoFvMeshCoupled.C
    dynamicTopoFvMeshReOrder.C
    edgeBisect.C
    edgeCollapse.C
    edgeSwap.C
    coupledPatchInfo.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicTopoFvMesh_H
#define dynamicTopoFvMesh_H

#include "Time.H"
#include "resizableList.H"
#include "multiThreader.H"
#include "threadHandler.H"
#include "eMesh.H"
#include "stack.H"
#include "coupleMap.H"
#include "changeMap.H"
#include "tetMetric.H"
#include "dynamicFvMesh.H"
#include "StaticHashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class motionSolver;
class topoMapper;
class objectMap;
class lengthScaleEstimator;

/*---------------------------------------------------------------------------*\
                      Class dynamicTopoFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicTopoFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Number of boundary-patches
        label numPatches_;

        //- Topology change flag
        bool topoChangeFlag_;

        //- Is this a subMesh?
        bool isSubMesh_;

        //- Dynamic mesh dictionary
        IOdictionary dict_;

        //- Should all options be mandatorily specified?
        //  Handy for first-time use.
        Switch mandatory_;

        //- Mesh characteristics [2D/3D]
        Switch twoDMesh_;

        //- Edge refinement switch
        Switch edgeRefinement_;

        //- Switch for cell-bandwidth reduction
        Switch bandWidthReduction_;

        //- Coupled modification switch
        mutable Switch coupledModification_;

        //- Specify the re-meshing interval
        label interval_;

        //- Field mapper
        autoPtr<topoMapper> mapper_;

        //- Mesh motion solver
        autoPtr<motionSolver> mPtr_;

        //- Length scale estimator
        autoPtr<lengthScaleEstimator> lsPtr_;

        //- Edge-mesh to store up-to-date connectivity information
        autoPtr<eMesh> eMeshPtr_;

        //- Lists that dynamically resize during topo-changes
        resizableList<point> oldPoints_, points_;
        resizableList<face> faces_;
        resizableList<label> owner_, neighbour_;
        resizableList<cell> cells_;
        resizableList<edge> edges_;
        resizableList<labelList> pointEdges_, edgePoints_;
        resizableList<labelList> edgeFaces_, faceEdges_;
        resizableList<scalar> lengthScale_;

        //- Size information
        labelList oldPatchSizes_, patchSizes_;
        labelList oldPatchStarts_, patchStarts_;
        labelList oldEdgePatchSizes_, edgePatchSizes_;
        labelList oldEdgePatchStarts_, edgePatchStarts_;
        labelList oldPatchNMeshPoints_, patchNMeshPoints_;
        label nOldPoints_, nPoints_;
        label nOldEdges_, nEdges_;
        label nOldFaces_, nFaces_;
        label nOldCells_, nCells_;
        label nOldInternalFaces_, nInternalFaces_;
        label nOldInternalEdges_, nInternalEdges_;

        //- Maps [Old-to-New]
        labelList reversePointMap_;
        labelList reverseEdgeMap_;
        labelList reverseFaceMap_;
        labelList reverseCellMap_;

        //- Maps [New-to-Old]
        labelList pointMap_;
        labelList edgeMap_;
        labelList faceMap_;
        labelList cellMap_;

        //- Maps for the renumbering of added entities
        Map<label> addedPointRenumbering_;
        Map<label> addedEdgeRenumbering_;
        Map<label> addedFaceRenumbering_;
        Map<label> addedCellRenumbering_;
        Map<label> addedFacePatches_;
        Map<label> addedEdgePatches_;
        Map<label> addedPointZones_;
        Map<label> addedFaceZones_;
        Map<label> addedCellZones_;

        // Information for field-mapping
        Map<labelList> faceParents_;
        Map<scalarField> faceWeights_;
        Map<vectorField> faceCentres_;
        Map<labelList> cellParents_;
        Map<scalarField> cellWeights_;
        Map<vectorField> cellCentres_;

        // Information for mapPolyMesh
        List<objectMap> pointsFromPoints_;
        List<objectMap> facesFromPoints_;
        List<objectMap> facesFromEdges_;
        List<objectMap> facesFromFaces_;
        List<objectMap> cellsFromPoints_;
        List<objectMap> cellsFromEdges_;
        List<objectMap> cellsFromFaces_;
        List<objectMap> cellsFromCells_;

        //- Maps to keep track of entities deleted after addition
        labelHashSet deletedPoints_;
        labelHashSet deletedEdges_;
        labelHashSet deletedFaces_;
        labelHashSet deletedCells_;

        //- List of flipped faces
        labelHashSet flipFaces_;

        //- Run-time statistics
        label nModifications_, maxModifications_;
        FixedList<label, 2> nBisections_;
        FixedList<label, 2> nCollapses_;
        FixedList<label, 2> nSwaps_;

        //- Sliver exudation
        scalar sliverThreshold_;
        Map<scalar> thresholdSlivers_;

        //- Specific to proximity-based refinement
        List<labelPair> slicePairs_;

        //- 3D Mesh Flipping data
        label maxTetsPerEdge_;
        Switch allowTableResize_;
        labelList noSwapPatchIDs_;

        //- Stack-list of entities to be checked for topo-changes.
        List<stack> entityStack_;

        //- Support for multithreading
        autoPtr<IOmultiThreader> threader_;

        // Pointer-list of thread-handlers
        PtrList<threadHandler<dynamicTopoFvMesh> > handlerPtr_;

        // Scheduling sequences for thread locking
        labelList topoSeq_;

        // Entity mutexes used for synchronization
        // in multi-threaded reOrdering
        FixedList<Mutex, 4> entityMutex_;

        // Private sub-class for coupled patches (local and processor)
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

        // Local coupled patch information
        PtrList<coupledPatchInfo> patchCoupling_;

        typedef StaticHashTable<empty, label, Hash<label> > labelStaticHashSet;

        // Coupled entities to avoid during normal modification
        labelStaticHashSet entitiesToAvoid_;

        // List of processors that this sub-domain talks to
        labelList procIndices_;

        // Sub-Mesh pointers
        PtrList<coupledPatchInfo> sendPatchMeshes_;
        PtrList<coupledPatchInfo> recvPatchMeshes_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicTopoFvMesh(const dynamicTopoFvMesh&);

        //- Construct from components
        //- Used for subMeshes only.
        dynamicTopoFvMesh
        (
            const dynamicTopoFvMesh& mesh,
            const IOobject& io,
            const pointField& points,
            const label nInternalEdges,
            const edgeList& edges,
            const faceList& faces,
            const labelListList& faceEdges,
            const labelList& owner,
            const labelList& neighbour,
            const cellList& cells
        );

        //- Disallow default bitwise assignment
        void operator=(const dynamicTopoFvMesh&);

        // Perform a Delaunay test on an internal face
        bool testDelaunay(const label fIndex) const;

        //- Quality metric for tetrahedra in 3D
        tetMetric::tetMetricReturnType tetMetric_;

        // Return the magnitude of an edge
        inline scalar edgeLength
        (
            const edge& edgeToCheck,
            const UList<point>& points
        ) const;

        // Compute the centre for a given face
        inline vector faceCentre
        (
            const face& faceToCheck,
            const UList<point>& points
        ) const;

        // Compute the normal for a given face
        inline vector faceNormal
        (
            const face& faceToCheck,
            const UList<point>& points
        ) const;

        // Given a face and cell index, find the apex point for a tet cell.
        inline label tetApexPoint
        (
            const label cIndex,
            const label fIndex
        ) const;

        // Given a cell index, find the centroid / volume of a cell.
        inline void cellCentreAndVolume
        (
            const label cIndex,
            const UList<point>& points,
            vector& centre,
            scalar& volume
        ) const;

        // Utility member function for separating axis tests
        inline label whichSide
        (
            const labelList& cellPoints,
            const UList<point>& points,
            const point& dir,
            const point& p
        ) const;

        // Determine whether a segment intersects a face
        inline bool segmentFaceIntersection
        (
            const edge& edgeToCheck,
            const face& faceToCheck,
            const UList<point>& points,
            vector& intersectionPoint
        ) const;

        // Determine whether the particular point lies
        // inside the given face
        inline bool pointInFace
        (
            const face& faceToCheck,
            const UList<point>& points,
            const point& checkPoint
        ) const;

        // Determine whether the particular point lies
        // inside the given cell
        inline bool pointInCell
        (
            const label cIndex,
            const cell& cellToCheck,
            const UList<face>& faces,
            const UList<label>& owner,
            const UList<point>& points,
            const point& checkPoint
        ) const;

        // Obtain map weighting factors for a face
        void computeFaceWeights
        (
            const label fIndex,
            const labelList& mapCandidates,
            labelList& parents,
            scalarField& weights,
            vectorField& centres
        ) const;

        // Obtain map weighting factors for a cell
        void computeCellWeights
        (
            const label cIndex,
            const labelList& mapCandidates,
            labelList& parents,
            scalarField& weights,
            vectorField& centres
        ) const;

        // Test for intersection between faces in old/new meshes
        bool testFaceIntersection
        (
            const label newFaceIndex,
            const label oldFaceIndex
        ) const;

        // Test for intersection between cells in old/new meshes
        bool testCellIntersection
        (
            const label newCellIndex,
            const label oldCellIndex
        ) const;

        // Return the intersection points between faces in old/new meshes
        bool faceIntersection
        (
            const label newFaceIndex,
            const label oldFaceIndex,
            vectorField& intPoints
        ) const;

        // Return the intersection points between cells in old/new meshes
        bool cellIntersection
        (
            const label newCellIndex,
            const label oldCellIndex,
            vectorField& intPoints
        ) const;

        // Compute the area / centre of a polygon
        // formed by a convex set of points.
        void convexSetArea
        (
            const label newFaceIndex,
            const label oldFaceIndex,
            const vectorField& cvxSet,
            scalar& fArea,
            vector& fCentre,
            bool output = false
        ) const;

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
        ) const;

        // Obtain a list of possible parent faces from the old mesh
        labelList faceParents
        (
            const label fIndex,
            const scalar searchFactor,
            const labelList& oldCandidates
        ) const;

        // Obtain a list of possible parent cells from the old mesh
        labelList cellParents
        (
            const label cIndex,
            const scalar searchFactor,
            const labelList& oldCandidates
        ) const;

        // Initialize mesh edges and related connectivity lists
        void initEdges();

        // Load the mesh-quality metric from the library
        void loadMetric();

        // Load the mesh-motion solver
        void loadMotionSolver();

        // Load the field mapper
        void loadFieldMapper();

        // Load the length scale estimator
        void loadLengthScaleEstimator();

        // Initialize the threading environment
        void initializeThreadingEnvironment(const label specThreads = -1);

        // Return a non-const reference to the lengthScaleEstimator
        inline lengthScaleEstimator& lengthEstimator();

        // Return a const reference to the lengthScaleEstimator
        inline const lengthScaleEstimator& lengthEstimator() const;

        // Return a const reference to the multiThreader
        inline const multiThreader& threader() const;

        // Return a reference to the entity mutexes
        inline const Mutex& entityMutex(const label entity) const;

        // Return the edge index for a provided edge
        inline label getEdgeIndex(const edge& edgeToCheck) const;

        // Given a boundary quad face, return a boundary triangular face.
        label getTriBoundaryFace(const label fIndex) const;

        // Given a boundary face, pick out a boundary edge that
        // contains a triangular face.
        label getTriBoundaryEdge(const label fIndex) const;

        // Return length-scale at an face-location in the mesh [2D]
        inline scalar faceLengthScale(const label fIndex) const;

        // Return length-scale at an edge-location in the mesh [3D]
        inline scalar edgeLengthScale(const label eIndex) const;

        // Check for bisection
        inline bool checkBisection(const label index) const;

        // Check for collapse
        inline bool checkCollapse(const label index) const;

        // MultiThreaded topology modifier [2D]
        void threadedTopoModifier2D();

        // MultiThreaded topology modifier [3D]
        void threadedTopoModifier3D();

        // Lock all slave threads
        template <class Type>
        void lockSlaveThreads
        (
            const labelList& sequence,
            PtrList<threadHandler<Type> >& handler
        );

        // Synchronize all slave threads
        template <class Type>
        void synchronizeThreads
        (
            const labelList& sequence,
            PtrList<threadHandler<Type> >& handler
        );

        // 2D Edge-swapping engine
        static void swap2DEdges(void *argument);

        // 2D Edge-bisection/collapse engine
        static void edgeBisectCollapse2D(void *argument);

        // 3D Edge-swapping engine
        static void swap3DEdges(void *argument);

        // 3D Edge-bisection/collapse engine
        static void edgeBisectCollapse3D(void *argument);

        // Return the entity stack for a particular thread
        inline stack& Stack(const label threadID);

        // Synchronize all processes
        inline void synchronize();

        // Return the integer threadID for a given pthread
        inline label self() const;

        // Initialize stacks
        inline void initStacks();

        // Initialize the coupled stack
        inline void initCoupledStack();

        // Method to determine the boundary patch index for a given face
        inline label whichPatch(const label index) const;

        // Method to determine the boundary patch index for a given edge
        inline label whichEdgePatch(const label index) const;

        // Method for the swapping of a quad-face in 2D
        void swapQuadFace
        (
            const label fIndex
        );

        // Method for the bisection of a quad-face in 2D
        const changeMap
        bisectQuadFace
        (
            const label fIndex,
            bool checkOnly = false,
            const changeMap& masterMap = changeMap()
        );

        // Method for the collapse of a quad-face in 2D
        const changeMap
        collapseQuadFace
        (
            const label fIndex,
            label overRideCase = -1,
            bool checkOnly = false
        );

        // Method for the swapping of an edge in 3D
        void swapEdge
        (
            const label eIndex,
            bool forceOp = false
        );

        // Method for the bisection of an edge in 3D
        const changeMap
        bisectEdge
        (
            const label eIndex,
            bool checkOnly = false,
            bool forceOp = false
        );

        // Method for the trisection of a face in 3D
        const changeMap
        trisectFace
        (
            const label fIndex,
            bool checkOnly = false,
            bool forceOp = false
        );

        // Method for the collapse of an edge in 3D
        const changeMap
        collapseEdge
        (
            const label eIndex,
            label overRideCase = -1,
            bool checkOnly = false,
            bool forceOp = false
        );

        // Slice the mesh at a particular location
        void sliceMesh(const labelPair& pointPair);

        // Dijkstra's algorithm for the shortest path problem
        bool Dijkstra
        (
            const Map<point>& points,
            const Map<edge>& edges,
            const label startPoint,
            const label endPoint,
            Map<label>& pi
        ) const;

        // Utility method to check for invalid face-collapse.
        bool checkCollapse
        (
            const labelList& triFaces,
            const FixedList<label,2>& c0BdyIndex,
            const FixedList<label,2>& c1BdyIndex,
            const FixedList<label,2>& original,
            const FixedList<label,2>& replacement,
            const bool checkNeighbour
        ) const;

        // Utility method to check for invalid edge-collapse.
        bool checkCollapse
        (
            const point& newPoint,
            const point& oldPoint,
            const label pointIndex,
            const label cellIndex,
            labelHashSet& cellsChecked,
            bool forceOp = false
        ) const;

        // Utility method to check for concurrent points.
        bool checkPointNearness
        (
            const pointField& points,
            const scalar magSqrTol
        ) const;

        // Check if the boundary face is adjacent to a sliver-cell,
        // and remove it by a two-step bisection/collapse operation.
        void remove2DSliver(const label fIndex);

        // Indentify the sliver type in 3D
        const changeMap identifySliverType(const label cIndex) const;

        // Remove sliver cells in 3D
        void removeSlivers();

        // Method to find the interior/boundary faces
        // for an input quad-face and adjacent triangle-prism cell.
        void findPrismFaces
        (
            const label fIndex,
            const label cIndex,
            FixedList<face,2>& bdyf,
            FixedList<label,2>& bidx,
            FixedList<face,2>& intf,
            FixedList<label,2>& iidx
        ) const;

        // Method to find the common-edge between two faces.
        bool findCommonEdge
        (
            const label first,
            const label second,
            label& common
        ) const;

        // Method to find the isolated point on a triangular face
        // that doesn't lie on the specified edge
        inline void findIsolatedPoint
        (
            const face& f,
            const edge& e,
            label& ptIndex,
            label& nextPtIndex
        ) const;

        // Utility method to find the isolated point
        // given two triangular faces.
        inline label findIsolatedPoint
        (
            const face& baseFace,
            const face& checkFace
        ) const;

        // Method to replace a label in a given list
        inline void replaceLabel
        (
            const label original,
            const label replacement,
            labelList& list
        ) const;

        // Method to insert labels in a face, so that
        // right-handedness is preserved.
        inline void insertPointLabels
        (
            const vector& refNorm,
            const vectorField& points,
            const labelHashSet& pLabels,
            face& modFace
        ) const;

        // Method to insert a label between two labels in a list
        inline void insertLabel
        (
            const label newLabel,
            const label labelA,
            const label labelB,
            labelList& list
        ) const;

        // Insert the specified cell to the mesh
        label insertCell
        (
            const cell& newCell,
            const scalar lengthScale,
            const label zoneID = -1
        );

        // Set fill-in mapping information for a particular cell
        void setCellMapping
        (
            const label cIndex,
            const labelList& mapCells,
            const scalarField& mapWeights,
            const vectorField& mapCentres
        );

        // Set fill-in mapping information for a particular face
        void setFaceMapping
        (
            const label fIndex,
            const labelList& mapFaces = labelList(0),
            const scalarField& mapWeights = scalarField(0),
            const vectorField& mapCentres = vectorField(0)
        );

        // Remove the specified cells from the mesh
        const changeMap
        removeCells
        (
            const labelList& cList,
            const label patch
        );

        // Remove the specified cell from the mesh
        void removeCell
        (
            const label cIndex
        );

        // Insert the specified face to the mesh
        label insertFace
        (
            const label patch,
            const face& newFace,
            const label newOwner,
            const label newNeighbour,
            const label zoneID = -1
        );

        // Remove the specified face from the mesh
        void removeFace
        (
            const label fIndex
        );

        // Split a set of internal faces into boundary faces,
        // and add them to the specified patch.
        void splitInternalFaces
        (
            const label patchIndex,
            const labelList& internalFaces,
            const Map<bool>& cellColors
        );

        // Merge two boundary faces into one internal face
        void mergeBoundaryFaces
        (
            const label firstFace,
            const label secondFace
        );

        // Insert the specified edge to the mesh
        label insertEdge
        (
            const label patch,
            const edge& newEdge,
            const labelList& edgeFaces,
            const labelList& edgePoints = labelList::null()
        );

        // Remove the specified edge from the mesh
        void removeEdge
        (
            const label eIndex
        );

        // Insert the specified point to the mesh
        label insertPoint
        (
            const point& newPoint,
            const point& oldPoint,
            const labelList& mapPoints,
            const label zoneID = -1
        );

        // Remove the specified point from the mesh
        void removePoint
        (
            const label pIndex
        );

        // Utility method to size-up the list to include an item
        template <class Type>
        inline void sizeUpList
        (
            const Type item,
            List<Type>& list
        ) const;

        // Utility method to size-down the list to remove an item
        template <class Type>
        inline void sizeDownList
        (
            const Type item,
            List<Type>& list
        ) const;

        // Remove an item at a particular index in the list
        template <class Type>
        inline void removeIndex
        (
            const label index,
            List<Type>& list
        ) const;

        // Utility method to build a hull of cells connected to the edge [2D]
        void constructPrismHull
        (
            const label eIndex,
            labelHashSet& hullTriFaces,
            labelHashSet& hullCells
        ) const;

        // Utility method to build a hull of cells around an edge [3D]
        void constructHull
        (
            const label eIndex,
            labelList& hullEdges,
            labelList& hullFaces,
            labelList& hullCells,
            labelListList& ringEntities
        ) const;

        // Utility method to build edgePoints for an edge [3D]
        void buildEdgePoints
        (
            const label eIndex,
            const label checkIndex = 0
        );

        // Utility to check whether points of an edge lie on a boundary
        const FixedList<bool,2> checkEdgeBoundary(const label eIndex) const;

        // Set a particular face index as flipped.
        inline void setFlip(const label fIndex);

        // Utility method to compute the minimum quality of a vertex hull
        inline scalar computeMinQuality(const label eIndex) const;

        // Utility method to compute the quality of a tetrahedral cell
        inline scalar tetQuality(const label cIndex) const;

        // Utility method to compute the quality of a vertex hull
        // around an edge after bisection.
        scalar computeBisectionQuality(const label eIndex) const;

        // Utility method to compute the quality of cells
        // around a face after trisection.
        scalar computeTrisectionQuality(const label fIndex) const;

        // Check whether the given edge is on a bounding curve
        bool checkBoundingCurve(const label eIndex) const;

        // Allocate dynamic programming tables
        void initTables
        (
            labelList& m,
            PtrList<scalarListList>& Q,
            PtrList<labelListList>& K,
            PtrList<labelListList>& triangulations,
            const label checkIndex = -1
        ) const;

        // Check triangulation quality for an edge index
        bool checkQuality
        (
            const label eIndex,
            const labelList& m,
            const PtrList<scalarListList>& Q,
            const scalar minQuality,
            const label checkIndex = 0
        ) const;

        // Utility method to fill the dynamic programming tables
        bool fillTables
        (
            const label eIndex,
            const scalar minQuality,
            labelList& m,
            PtrList<scalarListList>& Q,
            PtrList<labelListList>& K,
            PtrList<labelListList>& triangulations,
            const label checkIndex = 0
        ) const;

        // Print out tables for debugging
        void printTables
        (
            const labelList& m,
            const PtrList<scalarListList>& Q,
            const PtrList<labelListList>& K,
            const label checkIndex = 0
        ) const;

        // Remove the edge according to the swap sequence
        const changeMap
        removeEdgeFlips
        (
            const label eIndex,
            const scalar minQuality,
            const PtrList<labelListList>& K,
            PtrList<labelListList>& triangulations,
            const label checkIndex = 0
        );

        // Extract triangulations from the programming table
        void extractTriangulation
        (
            const label i,
            const label j,
            const labelListList& K,
            label& numTriangulations,
            labelListList& triangulations
        ) const;

        // Identify the 3-2 swap from the triangulation sequence
        label identify32Swap
        (
            const label eIndex,
            const labelList& hullVertices,
            const labelListList& triangulations
        ) const;

        // Check old-volumes for the configuration.
        bool checkTriangulationVolumes
        (
            const label eIndex,
            const labelList& hullVertices,
            const labelListList& triangulations
        ) const;

        // Routine to check whether the triangulation at the
        // index lies on the boundary of the vertex ring.
        bool boundaryTriangulation
        (
            const label index,
            label& isolatedVertex,
            labelListList& triangulations
        ) const;

        // Routine to perform 2-3 swaps
        const changeMap
        swap23
        (
            const label isolatedVertex,
            const label eIndex,
            const label triangulationIndex,
            const label numTriangulations,
            const labelListList& triangulations,
            const labelList& hullVertices,
            const labelList& hullFaces,
            const labelList& hullCells
        );

        // Routine to perform 3-2 or 2-2 swaps
        const changeMap
        swap32
        (
            const label eIndex,
            const label triangulationIndex,
            const label numTriangulations,
            const labelListList& triangulations,
            const labelList& hullVertices,
            const labelList& hullFaces,
            const labelList& hullCells
        );

        // Write out mesh bandwidth
        void printBandwidth() const;

        // Reorder points after a topology change
        void reOrderPoints
        (
            pointField& points,
            pointField& preMotionPoints,
            labelListList& pointZoneMap,
            bool threaded = false
        );

        // Static equivalent for multi-threading
        static void reOrderPointsThread(void *argument);

        // Reorder edges after a topology change
        void reOrderEdges
        (
            edgeList& edges,
            labelListList& edgeFaces,
            labelListList& faceEdges,
            bool threaded = false
        );

        // Static equivalent for multi-threading
        static void reOrderEdgesThread(void *argument);

        // Reorder faces in upper-triangular order after a topology change
        void reOrderFaces
        (
            faceList& faces,
            labelList& owner,
            labelList& neighbour,
            labelListList& faceEdges,
            labelListList& faceZoneFaceMap,
            bool threaded = false
        );

        // Static equivalent for multi-threading
        static void reOrderFacesThread(void *argument);

        // Reorder & renumber cells with bandwidth
        // reduction after a topology change
        void reOrderCells
        (
            labelListList& cellZoneMap,
            bool threaded = false
        );

        // Static equivalent for multi-threading
        static void reOrderCellsThread(void *argument);

        // Reorder the mesh in upper-triangular order,
        // and generate mapping information
        void reOrderMesh
        (
            pointField& points,
            pointField& preMotionPoints,
            edgeList& edges,
            faceList& faces,
            labelList& owner,
            labelList& neighbour,
            labelListList& faceEdges,
            labelListList& edgeFaces,
            labelListList& pointZoneMap,
            labelListList& faceZoneFaceMap,
            labelListList& cellZoneMap
        );

        // Invoke reOrdering with multiple threads
        void threadedMeshReOrdering
        (
            pointField& points,
            pointField& preMotionPoints,
            edgeList& edges,
            faceList& faces,
            labelList& owner,
            labelList& neighbour,
            labelListList& faceEdges,
            labelListList& edgeFaces,
            labelListList& pointZoneMap,
            labelListList& faceZoneFaceMap,
            labelListList& cellZoneMap
        );

        // Reset the mesh and generate mapping information
        void resetMesh();

        // Output an entity as a VTK file
        void writeVTK
        (
            const word& name,
            const label entity,
            const label primitiveType = 3,
            const bool useOldConnectivity = false,
            const bool useOldPoints = false
        ) const;

        // Output a list of entities as a VTK file
        void writeVTK
        (
            const word& name,
            const labelList& cList,
            const label primitiveType = 3,
            const bool useOldConnectivity = false,
            const bool useOldPoints = false
        ) const;

        // Actual routine to write out the VTK file
        void writeVTK
        (
            const word& name,
            const label nPoints,
            const label nCells,
            const label nTotalCells,
            const pointField& points,
            const labelListList& cpList = labelListList(0),
            const label primitiveType = 3,
            const Map<label>& reversePointMap = Map<label>(),
            const Map<label>& reverseCellMap = Map<label>()
        ) const;

        // Check the state of local connectivity lists
        void checkConnectivity(const label maxErrors = 0) const;

        // Handle mesh slicing events.
        void handleMeshSlicing();

        // Test an edge / face for proximity with other non-neighbouring faces.
        // Return the scalar distance to the nearest face.
        scalar testProximity
        (
            const label index,
            label& proximityFace
        ) const;

        // Calculate the edge length-scale for the mesh
        void calculateLengthScale(bool dump = false);

        // Fill buffers with length-scale info
        // and exchange across processors.
        void exchangeLengthBuffers();

        // Parallel send
        void pWrite(const label toID, const label& data) const;

        // Parallel send (for fixed lists)
        template <class Type, label Size>
        void pWrite(const label toID, const FixedList<Type,Size>& data) const;

        // Parallel send (for lists)
        template <class Type>
        void pWrite(const label toID, const List<Type>& data) const;

        // Parallel receive
        void pRead(const label fromID, label& data) const;

        // Parallel receive (for fixed lists)
        template <class Type, label Size>
        void pRead(const label fromID, FixedList<Type,Size>& data) const;

        // Parallel receive (for lists)
        template <class Type>
        void pRead(const label fromID, List<Type>& data) const;

        // Method to determine whether the master entity is locally coupled
        inline bool locallyCoupledEntity
        (
            const label index,
            bool checkSlaves = false
        ) const;

        // Method to determine the locally coupled patch index
        inline label locallyCoupledEdgePatch
        (
            const label eIndex
        ) const;

        // Method to determine whether the master entity is on a processor
        inline bool processorCoupledEntity(const label index) const;

        // Obtain the maximum coupled patch index
        inline label getMaxCouplingIndex() const;

        // Read optional dictionary parameters
        void readOptionalParameters();

        // Does the mesh perform edge refinement?
        inline bool edgeRefinement() const;

        // Identify coupled patches.
        bool identifyCoupledPatches();

        // Read coupled patch information from dictionary.
        void readCoupledPatches();

        // Initialize coupled patch connectivity for topology modifications.
        static void initCoupledConnectivity(void *argument);

        // Build a list of entities that need to be avoided
        // by regular topo-changes.
        inline void buildEntitiesToAvoid();

        // Check whether the specified edge is a coupled master edge.
        inline bool isCoupledMaster(const label eIndex) const;

        // Set coupled modification
        inline void setCoupledModification() const;

        // Unset coupled modification
        inline void unsetCoupledModification() const;

        // Handle topology changes for coupled patches
        void handleCoupledPatches();

        // Build patch sub-meshes for processor patches
        void buildProcessorPatchMeshes();

        // Build patch sub-mesh for a specified processor
        void buildProcessorPatchMesh
        (
            coupledPatchInfo& subMesh,
            labelHashSet& commonCells
        );

        // Build coupled maps for locally coupled patches
        void buildLocalCoupledMaps();

        // Build coupled maps for coupled processor patches
        void buildProcessorCoupledMaps();

        // Wait for buffer transfer completion.
        void waitForBuffers() const;

        // Write out proc IDs for post-processing
        void writeProcIDs() const;

        // Dump cell-quality statistics
        bool meshQuality(bool outputOption);

public:

    // Declare the name of the class and its debug switch
        TypeName("dynamicTopoFvMesh");

    // Constructors

        //- Construct from existing IOobject
        explicit dynamicTopoFvMesh(const IOobject& io);

    // Destructor

        virtual ~dynamicTopoFvMesh();

    // Member Functions

        // Map all fields in time using a customized mapper
        virtual void mapFields(const mapPolyMesh& meshMap);

        // Update the mesh for motion / topology changes
        // Return true if topology changes have occurred
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dynamicTopoFvMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
