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
#include "changeMap.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class mapPolyMesh;
class objectMap;

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

        //- Is this a subMesh?
        bool isSubMesh_;

        //- Dynamic mesh dictionary
        IOdictionary dict_;

        //- Should all options be mandatorily specified?
        //  Useful for first-time use.
        Switch mandatory_;

        //- Mesh characteristics [2D/3D]
        Switch twoDMesh_;

        //- Edge bisection/collapse option
        Switch edgeModification_;

        //- Coupled modification switches
        Switch coupledModification_;
        Switch slaveModification_;

        //- Specify the re-meshing interval
        label interval_;

        //- Sub-dictionary which specifies
        //  fixed length-scales for patches
        dictionary fixedPatches_;

        //- Sub-dictionary which specifies
        //  floating length-scales for patches
        dictionary freePatches_;

        //- Sub-dictionary which specifies
        //  patches with curvature-based length-scale
        dictionary curvaturePatches_;

        //- Sub-dictionary which specifies
        //  patches with proximity-based length-scale
        dictionary proximityPatches_;

        //- PolyMesh mapping structure
        autoPtr<mapPolyMesh> mapper_;

        //- Local copies of connectivity information
        resizableList<point> points_;
        resizableList<face>  faces_;
        resizableList<label> owner_;
        resizableList<label> neighbour_;
        resizableList<cell>  cells_;
        resizableList<edge>  edges_;
        resizableList<labelList> pointEdges_;
        resizableList<labelList> edgePoints_;
        resizableList<labelList> edgeFaces_;
        resizableList<labelList> faceEdges_;

        // Edge-mesh to store up-to-date connectivity information
        autoPtr<eMesh>  eMeshPtr_;

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

        //- Maps to keep track of entities deleted after addition
        labelHashSet    deletedPoints_;
        labelHashSet    deletedEdges_;
        labelHashSet    deletedFaces_;
        labelHashSet    deletedCells_;

        // Geometry-fields for mapping
        autoPtr<vectorField> pointPositionsPtr_;
        autoPtr<vectorField> cellCentresPtr_;
        autoPtr<vectorField> faceCentresPtr_;

        //- Edge bisection/collapse criteria
        resizableList<scalar> lengthScale_;
        scalar                ratioMin_;
        scalar                ratioMax_;
        scalar                growthFactor_;
        scalar                curvatureDeviation_;
        scalar                minLengthScale_;
        scalar                maxLengthScale_;
        scalar                sliverThreshold_;
        label                 nModifications_;
        label                 nBisections_;
        label                 nCollapses_;
        label                 nSwaps_;
        label                 maxModifications_;
        Map<scalar>           thresholdSlivers_;

        // Specific to proximity-based refinement
        label               spatialRes_;
        boundBox            proxBoundBox_;
        labelListList       proximityBins_;

        // Buffers for length-scale calculations
        labelListList       sendLblBuffer_;
        labelListList       recvLblBuffer_;
        scalarListList      sendSclBuffer_;
        scalarListList      recvSclBuffer_;

        //- 3D Mesh Flipping data
        label               maxTetsPerEdge_;
        Switch              allowTableResize_;
        labelList           noSwapPatchIDs_;

        // Stack-list of faces to be checked for topo-changes [2D].
        List<stack>         faceStack_;

        // Stack-list of edges to be checked for topo-changes [3D].
        List<stack>         edgeStack_;

        //- Support for multithreading
        autoPtr<multiThreader> threader_;

        // Pointer-list of thread-handlers
        PtrList<threadHandler<dynamicTopoFvMesh> > handlerPtr_;
        PtrList<threadHandler<dynamicTopoFvMesh> > reOrderPtr_;

        // Scheduling sequences for thread locking
        labelList           topoSeq_;
        labelList           reOrderSeq_;

        // Entity mutexes used for synchronization
        // in multi-threaded reOrdering
        FixedList<Mutex, 4> entityMutex_;

        // Private sub-class for coupled patches (local and processor)
        class coupledPatchInfo
        {
            // Auto pointer to a subMesh
            autoPtr<dynamicTopoFvMesh> subMesh_;

            // For locally coupled patches, specify a slave index
            label    slaveIndex_;

            // Entity sizes (as specified by the entityType enumerant)
            FixedList<label,6> nEntities_;

            // Maps for entities
            // - For every local index,
            //   provide global index
            FixedList<Map<label>,4> entityMap_;

            // - For every global index,
            //   provide local index
            FixedList<Map<label>,4> reverseEntityMap_;

            // Master-Slave mapping
            Map<label> masterToSlave_;
            Map<label> slaveToMaster_;

            // Entity Buffers (as specified by the entityType enumerant)
            FixedList<labelList,5> entityBuffer_;

            // Auxiliary buffers
            pointField  pointBuffer_;
            scalarList  lengthBuffer_;

            //- Disallow default bitwise assignment
            void operator=(const coupledPatchInfo&);

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

            // Constructor
            coupledPatchInfo(const label slaveIndex = -1);

            //- Access

                // Set a new subMesh
                void setMesh(label index, dynamicTopoFvMesh* mesh);

                // Return the slave index
                label slaveIndex();

                // Access to entity sizes
                FixedList<label,6>& nEntities();
                label& nEntities(const label eType);

                // Access to maps
                Map<label>& entityMap(const label eType);
                Map<label>& reverseEntityMap(const label eType);

                // Master-slave mapping
                label findSlaveIndex(const label Index);
                label findMasterIndex(const label Index);

                // Master/slave entries
                void mapSlave(const label master, const label slave);
                void mapMaster(const label slave, const label master);

                // Access to buffers
                FixedList<labelList,5>& entityBuffer();
                labelList& entityBuffer(const label eType);

                // Access to auxiliary buffers
                pointField& pointBuffer();
                scalarList& lengthBuffer();

                // Clear maps
                void clearMaps();
        };

        // Coupled patch information
        Map<coupledPatchInfo> patchCoupling_;

        // Coupled entities to avoid during normal modification
        labelHashSet entitiesToAvoid_;

        // List of processors that this domain talks to
        labelList procIndices_;

        // List of points that are required during patch sub-mesh creation
        labelListList subMeshPoints_;

        // Sub-Mesh pointers
        PtrList<coupledPatchInfo> sendPatchMeshes_;
        PtrList<coupledPatchInfo> recvPatchMeshes_;

        // Scalar tolerance for geometric matching
        scalar gTol_;

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
            const cellList& cells
        );

        //- Disallow default bitwise assignment
        void operator=(const dynamicTopoFvMesh&);

        // Perform a Delaunay test on an internal face
        bool testDelaunay(const label fIndex);

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

        // Obtain the tangent-to-edge vector
        inline vector tangentToEdge(const label eIndex);

        // Find the circumcenter of three points
        inline point circumCenter(const label fIndex);

        // Find the center of a triangle face
        inline point triFaceCenter(const label fIndex);

        // Find the area of a triangle face.
        inline scalar triFaceArea(const face& triFace);

        // Find the normal of a triangle face.
        inline vector triFaceNormal(const face& triFace);

        // Given a face and cell index, find the apex point for a tet cell.
        inline label tetApexPoint(const label cIndex, const label fIndex);

        // Compare two triangular faces
        inline label compare(const face& a, const face& b);

        // Find the volume of a tetrahedron.
        inline scalar tetVolume
        (
            const point& a,
            const point& b,
            const point& c,
            const point& d
        );

        // Find the volume sign of a tetrahedron
        inline label tetVolumeSign
        (
            const point& a,
            const point& b,
            const point& c,
            const point& d,
            const scalar tolerance
        );

        // Initialize mesh edges and related connectivity lists
        void initEdges();

        // Load the mesh-quality metric library
        void loadMetricLibrary();

        // Initialize the threading environment
        void initializeThreadingEnvironment
        (
            const label specThreads = -1
        );

        // Return a reference to the multiThreader
        const multiThreader& threader() const;

        // Return a reference to the entity mutexes
        inline Mutex& entityMutex(const label entity);

        // Given a boundary quad face, return a boundary triangular face.
        label getTriBoundaryFace(const label fIndex);

        // Given a boundary face, pick out a boundary edge that
        // contains a triangular face.
        label getTriBoundaryEdge(const label fIndex);

        // Return the appropriate length-scale for boundary face
        inline scalar boundaryLengthScale(const label fIndex);

        // Return length-scale at an face-location in the mesh [2D]
        inline scalar faceLengthScale(const label fIndex);

        // Return length-scale at an edge-location in the mesh [3D]
        inline scalar edgeLengthScale(const label eIndex);

        // Check for face bisection
        inline bool checkFaceBisection(const label fIndex);

        // Check for face collapse
        inline bool checkFaceCollapse(const label fIndex);

        // Check for edge bisection
        inline bool checkEdgeBisection(const label eIndex);

        // Check for edge collapse
        inline bool checkEdgeCollapse(const label eIndex);

        // Check if a given face is a quad
        inline bool checkQuadFace(const label fIndex);

        // Return the length of an edge
        inline scalar edgeLength(const label eIndex);

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

        // Return the face-stack for a particular thread
        inline stack& faceStack(const label threadID);

        // Return the edge-stack for a particular thread
        inline stack& edgeStack(const label threadID);

        // Return the integer threadID for a given pthread
        inline label self();

        // Initialize face-stacks
        inline void initFaceStacks();

        // Initialize edge-stacks
        inline void initEdgeStacks();

        // Initialize the coupled edge-stack
        inline void initCoupledEdgeStack();

        // Method to determine the boundary patch index for a given face
        inline label whichPatch(const label index);

        // Method to determine the boundary patch index for a given edge
        inline label whichEdgePatch(const label index);

        // Method for the swapping of a quad-face in 2D
        void swapQuadFace
        (
            const label fIndex
        );

        // Method for the bisection of a quad-face in 2D
        void bisectQuadFace
        (
            const label fIndex
        );

        // Method for the collapse of a quad-face in 2D
        void collapseQuadFace
        (
            const label fIndex
        );

        // Method for the bisection of an edge in 3D
        const changeMap
        bisectEdge
        (
            const label eIndex
        );

        // Method for the trisection of a face in 3D
        const changeMap
        trisectFace
        (
            const label fIndex
        );

        // Method for the collapse of an edge in 3D
        const changeMap
        collapseEdge
        (
            const label eIndex,
            label overRideCase = -1,
            bool checkOnly = false
        );

        // Slice the mesh at a particular edge in 3D
        void sliceMesh
        (
            const label eIndex
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

        // Utility method to check for invalid edge-collapse.
        bool checkCollapse
        (
            const point& newPoint,
            const label pointIndex,
            const label cellIndex,
            labelHashSet& cellsChecked
        );

        // Check if the boundary face is adjacent to a sliver-cell,
        // and remove it by a two-step bisection/collapse operation.
        void remove2DSliver
        (
            const label fIndex
        );

        // Indentify the sliver type in 3D
        const changeMap
        identifySliverType
        (
            const label cIndex
        );

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
        );

        // Method to find the common-edge between two faces.
        bool findCommonEdge
        (
            const label first,
            const label second,
            label& common
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
            const label original,
            const label replacement,
            labelList& list
        );

        // Method to insert a label between two labels in a list
        inline void insertLabel
        (
            const label newLabel,
            const label labelA,
            const label labelB,
            labelList& list
        );

        // Insert the specified cell to the mesh
        label insertCell
        (
            const cell& newCell,
            const scalar lengthScale = -1.0,
            const label mappingCell = -1
        );

        // Remove the specified cell from the mesh,
        // (including lower entities like points, edges, faces, etc.)
        // and add exposed internal faces to the specified patch
        void removeCell
        (
            const label index,
            const label patch
        );

        // Remove the specified cell from the mesh
        void removeCell
        (
            const label index
        );

        // Insert the specified face to the mesh
        label insertFace
        (
            const label patch,
            const face& newFace,
            const label newOwner,
            const label newNeighbour
        );

        // Remove the specified face from the mesh
        void removeFace
        (
            const label index
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
            const label index
        );

        // Insert the specified point to the mesh
        label insertPoint
        (
            const point& newPoint
        );

        // Remove the specified point from the mesh
        void removePoint
        (
            const label index
        );

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

        // Utility method to build a hull of cells connected to the edge [2D]
        void constructPrismHull
        (
            const label eIndex,
            labelHashSet& hullTriFaces,
            labelHashSet& hullCells
        );

        // Utility method to build a hull of cells around an edge [3D]
        void constructHull
        (
            const label eIndex,
            labelList& hullEdges,
            labelList& hullFaces,
            labelList& hullCells,
            labelListList& ringEntities
        );

        // Utility method to build edgePoints for an edge [3D]
        void buildEdgePoints
        (
            const label eIndex,
            const label checkIndex = 0
        );

        // Utility to invert a connectivity list (pointEdges from edges, etc)
        template<class InList, class OutList>
        void invertConnectivity
        (
            const label nEntities,
            const resizableList<InList>& inEntities,
            resizableList<OutList>& outEntities
        );

        // Utility to check whether points of an edge lie on a boundary
        const FixedList<bool,2>
        checkEdgeBoundary
        (
            const label eIndex
        );

        // Utility method to compute the minimum quality of a vertex hull
        inline scalar computeMinQuality(const label eIndex);

        // Utility method to compute the quality of a vertex hull
        // around an edge after bisection.
        scalar computeBisectionQuality(const label eIndex);

        // Utility method to compute the quality of cells
        // around a face after trisection.
        scalar computeTrisectionQuality(const label fIndex);

        // Check whether the given edge is on a bounding curve
        bool checkBoundingCurve(const label eIndex);

        // Allocate dynamic programming tables
        void initTables
        (
            labelList& m,
            PtrList<scalarListList>& Q,
            PtrList<labelListList>& K,
            PtrList<labelListList>& triangulations,
            const label checkIndex = -1
        );

        // Check triangulation quality for an edge index
        bool checkQuality
        (
            const label eIndex,
            const labelList& m,
            const PtrList<scalarListList>& Q,
            const scalar minQuality,
            const label checkIndex = 0
        );

        // Utility method to fill the dynamic programming tables
        bool fillTables
        (
            const label eIndex,
            labelList& m,
            PtrList<scalarListList>& Q,
            PtrList<labelListList>& K,
            PtrList<labelListList>& triangulations,
            const label checkIndex = 0
        );

        // Remove the edge according to the swap sequence
        bool removeEdgeFlips
        (
            const label eIndex,
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
        );

        // Identify the 3-2 swap from the triangulation sequence
        label identify32Swap
        (
            const label eIndex,
            const labelList& hullVertices,
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
            const label eIndex,
            const label triangulationIndex,
            const label numTriangulations,
            const labelListList& triangulations,
            const labelList& hullVertices,
            const labelList& hullFaces,
            const labelList& hullCells
        );

        // Routine to perform 3-2 or 2-2 swaps
        void swap32
        (
            const label eIndex,
            const label triangulationIndex,
            const label numTriangulations,
            const labelListList& triangulations,
            const labelList& hullVertices,
            const labelList& hullFaces,
            const labelList& hullCells
        );

        // Reorder points after a topology change
        void reOrderPoints
        (
            pointField& points,
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
            bool threaded = false
        );

        // Static equivalent for multi-threading
        static void reOrderFacesThread(void *argument);

        // Reorder & renumber cells with bandwidth
        // reduction after a topology change
        void reOrderCells(bool threaded = false);

        // Static equivalent for multi-threading
        static void reOrderCellsThread(void *argument);

        // Reorder the mesh in upper-triangular order,
        // and generate mapping information
        void reOrderMesh
        (
            pointField& points,
            edgeList& edges,
            faceList& faces,
            labelList& owner,
            labelList& neighbour,
            labelListList& faceEdges,
            labelListList& edgeFaces
        );

        // Invoke reOrdering with multiple threads
        void threadedMeshReOrdering
        (
            pointField& points,
            edgeList& edges,
            faceList& faces,
            labelList& owner,
            labelList& neighbour,
            labelListList& faceEdges,
            labelListList& edgeFaces
        );

        // Output a list of cells as a VTK file
        void writeVTK
        (
            const word& name,
            const labelList& cList,
            const label primitiveType = 3
        );

        // Check the state of connectivity local lists
        void checkConnectivity();

        // Perform spatial hashing on a set of points
        void spatialHash
        (
            const pointField& pointLocations,
            const labelList& pointIndices,
            const boundBox& box,
            const label resolution,
            labelListList& bins,
            label removeIndex = false
        );

        // Prepare for proximity-based refinement, if necessary
        void prepareProximityPatches();

        // Test an edge for proximity with other non-neighbouring faces.
        // Return the scalar distance to the nearest face.
        scalar testProximity
        (
            const label eIndex,
            label& proximityFace
        );

        // Calculate the edge length-scale for the mesh
        void calculateLengthScale();

        // Compute the growth factor of an existing mesh
        scalar computeGrowthFactor();

        // Send length-scale info across processors
        void writeLengthScaleInfo
        (
            const labelList& cellLevels,
            const resizableList<scalar>& lengthScale
        );

        // Receive length-scale info across processors
        void readLengthScaleInfo
        (
            const label level,
            label& visitedCells,
            labelList& cellLevels,
            resizableList<scalar>& lengthScale,
            labelHashSet& levelCells
        );

        // Parallel send
        void pWrite(const label toID, const label& data);

        // Parallel send (for fixed lists)
        template <class Type, label Size>
        void pWrite(const label toID, const FixedList<Type,Size>& data);

        // Parallel send (for lists)
        template <class Type>
        void pWrite(const label toID, const List<Type>& data);

        // Parallel receive
        void pRead(const label fromID, label& data);

        // Parallel receive (for fixed lists)
        template <class Type, label Size>
        void pRead(const label fromID, FixedList<Type,Size>& data);

        // Parallel receive (for lists)
        template <class Type>
        void pRead(const label fromID, List<Type>& data);

        // Method to determine whether the master face is locally coupled
        inline bool locallyCoupledFace(const label fIndex);

        // Method to determine whether the master edge is locally coupled
        inline bool locallyCoupledEdge
        (
            const label eIndex,
            bool checkSlaves = false
        );

        // Method to determine whether the master face is on a processor
        inline bool processorCoupledFace(const label fIndex);

        // Method to determine whether the master edge is on a processor
        inline bool processorCoupledEdge(const label eIndex);

        // Obtain the maximum coupled patch index
        inline label getMaxCouplingIndex();

        // Read optional dictionary parameters
        void readOptionalParameters();

        // Read edge-modification options from the dictionary
        void readEdgeOptions(bool reRead = false);

        // Check for legitimacy of patches
        void checkPatches(const wordList& patchList);

        // Does the mesh perform edge-modification?
        bool edgeModification();

        // Build the global shared-point list
        void buildGlobalPoints();

        // Identify coupled patches.
        bool identifyCoupledPatches();

        // Read coupled patch information from dictionary.
        void readCoupledPatches();

        // Initialize coupled patches for topology modifications.
        void initCoupledPatches();

        // Build a list of entities that need to be avoided
        // by regular topo-changes.
        inline void buildEntitiesToAvoid();

        // Check whether the specified edge is a coupled master edge.
        inline bool isCoupledMaster(const label eIndex);

        // Handle topology changes for coupled patches
        void handleCoupledPatches();

        // Build patch sub-meshes for processor patches
        void buildCoupledPatchMeshes();

        // Build patch sub-mesh for a specified processor
        void buildCoupledPatchMesh
        (
            const label proc,
            coupledPatchInfo& subMesh,
            labelHashSet& cellsToAvoid
        );

        // Build coupled maps for locally coupled patches
        void buildLocalCoupledMaps();

        // Build coupled maps for coupled processor patches
        void buildProcCoupledMaps();

        // Wait for buffer transfer completion.
        void waitForBuffers();

public:

    // Declare the name of the class and its debug switch
        ClassName("dynamicTopoFvMesh");

    // Constructors

        //- Construct from existing IOobject
        explicit dynamicTopoFvMesh(const IOobject& io);

    // Destructor

        virtual ~dynamicTopoFvMesh();

    // Member Functions

        // Return reference to the dictionary
        const IOdictionary& dynamicMeshDict() const;

        // Update the mesh for topology changes
        // Return true if changes have occurred
        bool updateTopology();

        //- Update mesh corresponding to the given map
        virtual void updateMesh(const mapPolyMesh& mpm);

        // Map all fields in time using the given map.
        virtual void mapFields(const mapPolyMesh& mpm);

        // Return the mesh-mapper
        const mapPolyMesh& meshMap();

        // Return old cell-centre information (prior to a topology change)
        const vectorField& oldCellCentres() const;

        // Return old face-centre information (prior to a topology change)
        const vectorField& oldFaceCentres() const;

        // Return old face-centre information (prior to a topology change)
        const vectorField& oldPointPositions() const;

        // Return mesh length-scale values
        tmp<scalarField> lengthScale();

        // Return mesh cell-quality values
        tmp<scalarField> meshQuality(bool outputOption);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dynamicTopoFvMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
