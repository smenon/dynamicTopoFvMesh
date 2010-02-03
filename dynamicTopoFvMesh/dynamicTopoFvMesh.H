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

        //- Edge refinement switch
        Switch edgeRefinement_;

        //- Coupled modification switches
        mutable Switch coupledModification_;
        mutable Switch slaveModification_;

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
        Map<objectMap>  pointsFromPoints_;
        Map<objectMap>  cellsFromCells_;
        Map<label>      addedPointRenumbering_;
        Map<label>      addedEdgeRenumbering_;
        Map<label>      addedFaceRenumbering_;
        Map<label>      addedCellRenumbering_;
        Map<label>      addedFacePatches_;
        Map<label>      addedEdgePatches_;
        Map<label>      addedPointZones_;
        Map<label>      addedFaceZones_;
        Map<label>      addedCellZones_;
        Map<label>      pointParents_;
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
        labelList             noModPatchIDs_;

        // Specific to proximity-based refinement
        label                   spatialRes_;
        boundBox                proxBoundBox_;
        labelListList           proximityBins_;
        scalar                  sliceThreshold_;
        label                   sliceHoldOff_;
        mutable List<boundBox>  sliceBoxes_;
        mutable List<labelPair> slicePairs_;

        // Buffers for length-scale calculations
        labelListList       sendLblBuffer_;
        labelListList       recvLblBuffer_;
        scalarListList      sendSclBuffer_;
        scalarListList      recvSclBuffer_;

        //- 3D Mesh Flipping data
        mutable label       maxTetsPerEdge_;
        Switch              allowTableResize_;
        labelList           noSwapPatchIDs_;

        // Stack-list of faces to be checked for topo-changes [2D].
        List<stack>         faceStack_;

        // Stack-list of edges to be checked for topo-changes [3D].
        List<stack>         edgeStack_;

        //- Support for multithreading
        autoPtr<IOmultiThreader> threader_;

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
            // Reference to the parent mesh
            const dynamicTopoFvMesh& mesh_;

            // Auto pointer to a subMesh
            autoPtr<dynamicTopoFvMesh> subMesh_;

            // Flag to determine whether maps have been built.
            bool builtMaps_;

            // For locally coupled patches,
            // specify a master / slave index
            label masterIndex_;
            label slaveIndex_;

            // Zone IDs for patches associated with faceZones
            label masterFaceZone_;
            label slaveFaceZone_;

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
            coupledPatchInfo
            (
                const dynamicTopoFvMesh& mesh,
                const label masterIndex = -1,
                const label slaveIndex = -1,
                const label mfzIndex = -1,
                const label sfzIndex = -1
            );

            //- Access

                // Set a new subMesh
                void setMesh(label index, dynamicTopoFvMesh* mesh);

                // Have maps been built?
                bool builtMaps() const;

                // Change the flag
                void setBuiltMaps();

                // Return the slave index
                label slaveIndex() const;

                // Return the master / slave face zone IDs
                label masterFaceZone() const;
                label slaveFaceZone() const;

                // Access to entity sizes
                FixedList<label,6>& nEntities();
                label& nEntities(const label eType);

                // Access to maps
                Map<label>& entityMap(const label eType);
                Map<label>& reverseEntityMap(const label eType);

                // Access to master / slave maps
                const Map<label>& masterToSlaveMap();
                const Map<label>& slaveToMasterMap();

                // Master-slave mapping
                label findSlaveIndex(const label Index) const;
                label findMasterIndex(const label Index) const;

                // Remove a master/slave entry, if found.
                void removeSlaveIndex(const label Index);
                void removeMasterIndex(const label Index);

                // Master/slave entries
                void mapSlave(const label master, const label slave);
                void mapMaster(const label slave, const label master);

                // Access to buffers
                FixedList<labelList,5>& entityBuffer();
                labelList& entityBuffer(const label eType);

                // Access to auxiliary buffers
                pointField& pointBuffer();
                scalarList& lengthBuffer();

                // Transfer maps, replacing old ones
                void transferMaps
                (
                    Map<label>& newMasterToSlave,
                    Map<label>& newSlaveToMaster
                );

                // Clear maps
                void clearMaps();
        };

        // Coupled patch information
        Map<coupledPatchInfo> patchCoupling_;

        // Coupled entities to avoid during normal modification
        DynamicList<label> entitiesToAvoid_;

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
        bool testDelaunay(const label fIndex) const;

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

        // Obtain the edge centre
        inline point edgeCentre(const label eIndex) const;

        // Obtain the tangent-to-edge vector
        inline vector tangentToEdge(const label eIndex) const;

        // Find the circumcentre of three points
        inline point circumCentre(const label fIndex) const;

        // Find the centre of a quad face
        inline point quadFaceCentre(const label fIndex) const;

        // Find the centre of a triangle face
        inline point triFaceCentre(const label fIndex) const;

        // Find the area of a quad face.
        inline scalar quadFaceArea(const face& quadFace) const;

        // Find the area of a triangle face.
        inline scalar triFaceArea(const face& triFace) const;

        // Find the normal of a quad face.
        inline vector quadFaceNormal(const face& quadFace) const;

        // Find the normal of a triangle face.
        inline vector triFaceNormal(const face& triFace) const;

        // Given a face and cell index, find the apex point for a tet cell.
        inline label tetApexPoint
        (
            const label cIndex,
            const label fIndex
        ) const;

        // Given a cell index, find the centroid of a prism cell.
        inline vector prismCellCentre(const label cIndex) const;

        // Given a cell index, find the centroid of a tet cell.
        inline vector tetCellCentre(const label cIndex) const;

        // Compare two triangular faces
        inline label compare(const face& a, const face& b) const;

        // Find the volume of a tetrahedron.
        inline scalar tetVolume
        (
            const point& a,
            const point& b,
            const point& c,
            const point& d
        ) const;

        // Find the volume sign of a tetrahedron
        inline label tetVolumeSign
        (
            const point& a,
            const point& b,
            const point& c,
            const point& d,
            const scalar tolerance,
            scalar& vol
        ) const;

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
        inline const Mutex& entityMutex(const label entity) const;

        // Given a boundary quad face, return a boundary triangular face.
        label getTriBoundaryFace(const label fIndex) const;

        // Given a boundary face, pick out a boundary edge that
        // contains a triangular face.
        label getTriBoundaryEdge(const label fIndex) const;

        // Return the appropriate length-scale for boundary face
        inline scalar boundaryLengthScale(const label fIndex) const;

        // Return length-scale at an face-location in the mesh [2D]
        inline scalar faceLengthScale(const label fIndex) const;

        // Return length-scale at an edge-location in the mesh [3D]
        inline scalar edgeLengthScale(const label eIndex) const;

        // Check for face bisection
        inline bool checkFaceBisection(const label fIndex) const;

        // Check for face collapse
        inline bool checkFaceCollapse(const label fIndex) const;

        // Check for edge bisection
        inline bool checkEdgeBisection(const label eIndex) const;

        // Check for edge collapse
        inline bool checkEdgeCollapse(const label eIndex) const;

        // Check if a given face is a quad
        inline bool checkQuadFace(const label fIndex) const;

        // Return the length of an edge
        inline scalar edgeLength(const label eIndex) const;

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
        inline label self() const;

        // Initialize face-stacks
        inline void initFaceStacks();

        // Initialize edge-stacks
        inline void initEdgeStacks();

        // Initialize the coupled face-stack
        inline void initCoupledFaceStack();

        // Initialize the coupled edge-stack
        inline void initCoupledEdgeStack();

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

        // Method for the bisection of an edge in 3D
        const changeMap
        bisectEdge
        (
            const label eIndex,
            bool checkOnly = false,
            bool forceOp = false,
            const changeMap& masterMap = changeMap()
        );

        // Method for the trisection of a face in 3D
        const changeMap
        trisectFace
        (
            const label fIndex,
            bool checkOnly = false,
            bool forceOp = false,
            const changeMap& masterMap = changeMap()
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
            const label pointIndex,
            const label cellIndex,
            labelHashSet& cellsChecked,
            bool forceOp = false
        ) const;

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

        // Agglomerate tetrahedral cells around a certain cell index
        void agglomerateTetCells
        (
            const label cIndex,
            DynamicList<label>& aggCells,
            DynamicList<label>& polyCell,
            DynamicList<face>& polyCellFaces
        ) const;

        // Small polyhedron reconnection for sliver removal
        bool smallPolyhedronReconnection
        (
            const scalar minQuality,
            DynamicList<face>& polyCellFaces
        );

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

        // Method to replace a label in a given list
        inline void replaceLabel
        (
            const label original,
            const label replacement,
            labelList& list
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
            const scalar lengthScale = -1.0,
            const label mappingCell = -1,
            const label zoneID = -1
        );

        // Remove the specified cell from the mesh,
        // (including lower entities like points, edges, faces, etc.)
        // and add exposed internal faces to the specified patch
        void removeCell
        (
            const label cIndex,
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
            const labelList& mappingPoints,
            const label zoneID = -1
        );

        // Remove the specified point from the mesh
        void removePoint
        (
            const label pIndex
        );

        // Utility method to size-up the list to include an item
        inline void sizeUpList
        (
            const label item,
            labelList& list
        ) const;

        // Utility method to size-down the list to remove an item
        inline void sizeDownList
        (
            const label item,
            labelList& list
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

        // Utility to invert a connectivity list (pointEdges from edges, etc)
        template<class InList, class OutList>
        void invertConnectivity
        (
            const label nEntities,
            const resizableList<InList>& inEntities,
            resizableList<OutList>& outEntities
        ) const;

        // Utility to check whether points of an edge lie on a boundary
        const FixedList<bool,2>
        checkEdgeBoundary
        (
            const label eIndex
        ) const;

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

        // Check whether the given face should not be bisected/collapsed
        inline bool checkFaceModification(const label eIndex) const;

        // Check whether the given edge should not be bisected/collapsed
        inline bool checkEdgeModification(const label eIndex) const;

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
            const labelListList& triangulations,
            const scalar tolFraction
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

        // Reorder points after a topology change
        void reOrderPoints
        (
            pointField& points,
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

        // Output a list of cells as a VTK file
        void writeVTK
        (
            const word& name,
            const labelList& cList,
            const label primitiveType = 3
        ) const;

        // Check the state of local connectivity lists
        void checkConnectivity(label maxErrors = 0) const;

        // Perform spatial hashing on a set of points
        void spatialHash
        (
            const pointField& pointLocations,
            const labelList& pointIndices,
            const boundBox& box,
            const label resolution,
            labelListList& bins,
            label removeIndex = false
        ) const;

        // Prepare for proximity-based refinement, if necessary
        void prepareProximityPatches();

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
        void calculateLengthScale();

        // Calculate geometric tolerance for the mesh
        void calculateGeometricTolerance();

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

        // Method to determine whether the master face is locally coupled
        inline bool locallyCoupledFace(const label fIndex) const;

        // Method to determine whether the master edge is locally coupled
        inline bool locallyCoupledEdge
        (
            const label eIndex,
            bool checkSlaves = false
        ) const;

        // Method to determine whether the master face is on a processor
        inline bool processorCoupledFace(const label fIndex) const;

        // Method to determine whether the master edge is on a processor
        inline bool processorCoupledEdge(const label eIndex) const;

        // Obtain the maximum coupled patch index
        inline label getMaxCouplingIndex() const;

        // Read optional dictionary parameters
        void readOptionalParameters();

        // Read edge refinement options from the dictionary
        void readRefinementOptions(bool reRead = false);

        // Check for legitimacy of patches
        void checkPatches(const wordList& patchList) const;

        // Does the mesh perform edge refinement?
        bool edgeRefinement() const;

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
        inline bool isCoupledMaster(const label eIndex) const;

        // Set coupled modification
        inline void setCoupledModification() const;

        // Unset coupled modification
        inline void unsetCoupledModification() const;

        // Set slave modification
        inline void setSlaveModification() const;

        // Unset slave modification
        inline void unsetSlaveModification() const;

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
        void waitForBuffers() const;

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
        const mapPolyMesh& meshMap() const;

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
