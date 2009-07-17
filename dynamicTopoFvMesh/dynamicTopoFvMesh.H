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

\*---------------------------------------------------------------------------*/

#ifndef dynamicTopoFvMesh_H
#define dynamicTopoFvMesh_H

#include "Time.H"
#include "HashList.H"
#include "multiThreader.H"
#include "eMesh.H"
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

        //- Dynamic mesh dictionary
        IOdictionary dict_;

        //- Mesh characteristics [2D/3D]
        Switch twoDMesh_;

        //- Edge bisection/collapse option
        Switch edgeModification_;

        //- Specify the re-meshing interval
        label interval_;

        //- Sub-dictionary which specifies
        //  fixed length-scales for patches
        dictionary fixedPatches_;

        //- Sub-dictionary which specifies
        //  patches with curvature-based length-scale
        dictionary curvaturePatches_;

        //- List of noSwap patch IDs (which specifies patches on
        //  which 2-2 swaps are to be avoided)
        labelHashSet noSwapPatchIDs_;

        //- PolyMesh mapping structure
        autoPtr<mapPolyMesh> mapper_;

        //- Local copies of connectivity information
        HashList<point> points_;
        HashList<face>  faces_;
        HashList<label> owner_;
        HashList<label> neighbour_;
        HashList<cell>  cells_;
        HashList<edge>  edges_;
        HashList<labelList> pointEdges_;
        HashList<labelList> edgePoints_;
        HashList<labelList> edgeFaces_;
        HashList<labelList> faceEdges_;

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
        Map<label>      faceParents_;

        // Geometry-fields for mapping
        autoPtr<vectorField> pointPositionsPtr_;
        autoPtr<vectorField> cellCentresPtr_;
        autoPtr<vectorField> faceCentresPtr_;

        //- Edge bisection/collapse criteria
        HashList<scalar>    lengthScale_;
        scalar              ratioMin_;
        scalar              ratioMax_;
        scalar              growthFactor_;
        scalar              maxLengthScale_;
        scalar              sliverThreshold_;
        scalar              curvatureRatio_;
        label               nModifications_;
        label               nBisections_;
        label               nCollapses_;
        label               nSwaps_;
        label               maxModifications_;
        label               bisectInterior_;
        labelHashSet        thresholdSlivers_;
        scalarListList      curvatureFields_;

        // Coupled patch information
        Map<label>          patchCoupling_;

        // Buffers for length-scale calculations
        labelListList       sendLblBuffer_;
        labelListList       recvLblBuffer_;
        scalarListList      sendSclBuffer_;
        scalarListList      recvSclBuffer_;

        //- 3D Mesh Flipping data
        label               maxTetsPerEdge_;
        Switch              allowTableResize_;

        // Private Stack sub-class
        class stack
        {
            // Private data
            DynamicList<label> stack_;

            Mutex              stackMutex_;

        public:

            // Constructor
            stack(){}

            // Member functions for access to the stack

                //- Push items on to the stack
                inline void push(const label index);

                //- Insert item onto stack
                inline void insert(const label index);

                //- Pop an item off the stack
                inline label pop();

                //- Remove a specific index off the stack
                inline void remove(const label index);

                //- Return if a stack is empty or not
                inline bool empty();

                //- Return the size of the stack
                inline label size();

                //- Clear out the stack
                inline void clear();

                //- Print out the stack
                inline void print();

                //- Access the stack as a list
                inline label operator[](const label index);
        };

        // Stack-list of faces to be checked for topo-changes [2D].
        // One for each thread.
        List<stack> faceStack_;

        // Stack-list of edges to be checked for topo-changes [3D].
        // One for each thread.
        List<stack> edgeStack_;

        //- Support for multithreading
        autoPtr<multiThreader> threader_;

        //- Private sub-class used for multiThreading
        class topoMeshStruct
        {
            // Pointer to this class
            dynamicTopoFvMesh *mesh_;

            // Total number of threads
            label nThreads_;

            // ID generated by the pthreads API
            pthread_t pthreadID_;

            // Is this a master/slave thread
            bool master_;

            // Is coupled modifications required?
            bool isCoupled_;

            // Synchronization mutex
            Mutex syncMutex_;

        public:

            // Constructor
            topoMeshStruct
            (
                dynamicTopoFvMesh *mesh,
                label nThreads
            );

            // Return a reference to the mesh
            dynamicTopoFvMesh& mesh() { return *mesh_; }

            // Return the number of threads
            inline label nThreads() { return nThreads_; }

            // Designate as master thread
            inline void setMaster() { master_ = true; }

            // Designate as slave thread
            inline void setSlave() { master_ = false; }

            // Is this a master thread?
            inline bool master() { return (master_ == true); }

            // Is this a slave thread?
            inline bool slave() { return !master(); }

            // Lock this thread
            inline void lock() { syncMutex_.lock(); }

            // Unlock this thread
            inline void unlock() { syncMutex_.unlock(); }

            // Return the ID
            inline pthread_t& ID() { return pthreadID_; }

            // Is coupled modifications required?
            inline bool& isCoupled() { return isCoupled_; }
        };

        PtrList<topoMeshStruct> structPtr_;

        // Private sub-class for coupled patches
        class patchSubMesh
        {
            // Entity sizes
            label    nPoints_;
            label    nEdges_;
            label    nFaces_;
            label    nCells_;

            // Maps for entities
            // - For every local index,
            //   provide global index
            Map<label> pointMap_;
            Map<label> edgeMap_;
            Map<label> faceMap_;
            Map<label> cellMap_;

            // - For every global index,
            //   provide local index
            Map<label> reversePointMap_;
            Map<label> reverseEdgeMap_;
            Map<label> reverseFaceMap_;
            Map<label> reverseCellMap_;

            // Buffers
            pointField  pointBuffer_;
            labelList   edgeBuffer_;
            labelList   faceBuffer_;
            labelList   cellBuffer_;
            scalarList  lengthBuffer_;

        public:

            // Constructor
            patchSubMesh();

            //- Access

                // Access to entity sizes
                label& nPoints() { return nPoints_; }
                label& nEdges() { return nEdges_; }
                label& nFaces() { return nFaces_; }
                label& nCells() { return nCells_; }

                // Access to maps
                Map<label>& pointMap() { return pointMap_; }
                Map<label>& edgeMap() { return edgeMap_; }
                Map<label>& faceMap() { return faceMap_; }
                Map<label>& cellMap() { return cellMap_; }

                Map<label>& reversePointMap() { return reversePointMap_; }
                Map<label>& reverseEdgeMap() { return reverseEdgeMap_; }
                Map<label>& reverseFaceMap() { return reverseFaceMap_; }
                Map<label>& reverseCellMap() { return reverseCellMap_; }

                // Access to buffers
                pointField& pointBuffer() { return pointBuffer_; }
                labelList& edgeBuffer() { return edgeBuffer_; }
                labelList& faceBuffer() { return faceBuffer_; }
                labelList& cellBuffer() { return cellBuffer_; }
                scalarList& lengthBuffer() { return lengthBuffer_; }
        };

        PtrList<patchSubMesh> sendPatchMeshes_;
        PtrList<patchSubMesh> recvPatchMeshes_;

        // Mapping for coupled patches
        List<Map<label> > masterToSlave_;
        List<Map<label> > slaveToMaster_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicTopoFvMesh(const dynamicTopoFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicTopoFvMesh&);

        // Perform a Delaunay test on an internal face
        void testDelaunay
        (
            const label fIndex,
            bool& failed
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

        // Return a reference to the multiThreader
        const multiThreader& threader() const;

        // Given a boundary quad face, return a boundary triangular face.
        label getTriBoundaryFace(const label faceIndex);

        // Given a boundary face, pick out a boundary edge that
        // contains a triangular face.
        label getTriBoundaryEdge(const label faceIndex);

        // Return the appropriate length-scale for boundary face
        scalar boundaryLengthScale(const label faceIndex);

        // Compute the face-curvature for a given boundary face
        scalar faceCurvature(const label fIndex);

        // Return length-scale at an face-location in the mesh [2D]
        inline void meshFaceLengthScale
        (
            const label fIndex,
            scalar& scale
        );

        // Return length-scale at an edge-location in the mesh [3D]
        inline void meshEdgeLengthScale
        (
            const label eIndex,
            scalar& scale
        );

        // Check if a given face is a quad
        inline bool checkQuadFace(const label fIndex);

        // Return the length of an edge
        inline void edgeLength(const label eIndex, scalar& length);

        // MultiThreaded topology modifier [2D]
        void threadedTopoModifier2D();

        // MultiThreaded topology modifier [3D]
        void threadedTopoModifier3D();

        // Lock all slave threads
        void lockSlaveThreads();

        // Unlock all slave threads
        void unlockSlaveThreads();

        // Synchronize all slave threads
        void synchronizeThreads();

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
        void bisectEdge
        (
            const label eIndex
        );

        // Method for the collapse of an edge in 3D
        void collapseEdge
        (
            const label eIndex
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

        // Remove the specified cell from the mesh,
        // and add internal faces to the specified patch
        void removeCell
        (
            const label index,
            const label patch
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

        // Split an internal face into two boundary faces
        void splitInternalFace
        (
            const label internalFace
        );

        // Merge two boundary faces into one internal face
        label mergeBoundaryFaces
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

        // Check for the occurrence of a label in the list
        inline bool foundInList
        (
            const label item,
            const labelList& list,
            label& index
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
            const label eIndex
        );

        // Utility to check whether points of an edge lie on a boundary
        void checkEdgeBoundary
        (
            const label eIndex,
            FixedList<bool,2>& edgeBoundary
        );

        // Utility method to compute the minimum quality of a vertex hull
        inline void computeMinQuality
        (
            const label eIndex,
            scalar& minQuality
        );

        // Check whether the given edge is on a bounding curve
        bool checkBoundingCurve(label eIndex);

        // Allocate dynamic programming tables
        void initTables
        (
            scalarListList& Q,
            labelListList& K,
            labelListList& triangulations
        );

        // Check triangulation quality for an edge index
        bool checkQuality
        (
            const label eIndex,
            const label m,
            const scalarListList& Q,
            const scalar minQuality
        );

        // Utility method to fill the dynamic programming tables
        bool fillTables
        (
            const label eIndex,
            label& m,
            scalarListList& Q,
            labelListList& K,
            labelListList& triangulations
        );

        // Remove the edge according to the swap sequence
        void removeEdgeFlips
        (
            const label eIndex,
            const labelListList& K,
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
        void reOrderPoints(pointField& points);

        // Reorder edges after a topology change
        void reOrderEdges(edgeList& edges, labelListList& edgeFaces);

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
            edgeList& edges,
            faceList& faces,
            labelList& owner,
            labelList& neighbour,
            labelListList& edgeFaces
        );

        // Output a list of cells as a VTK file
        void writeVTK(const word& name, const labelList& cList);

        // Copy edge-based connectivity to HashLists
        void setEdgeConnectivity();

        // Check the state of connectivity HashLists
        void checkConnectivity();

        // Calculate the edge length-scale for the mesh
        void calculateLengthScale();

        // Send length-scale info across processors
        void writeLengthScaleInfo
        (
            const labelList& cellLevels,
            const scalarList& lengthScale
        );

        // Receive length-scale info across processors
        void readLengthScaleInfo
        (
            const label level,
            label& visitedCells,
            labelList& cellLevels,
            scalarList& lengthScale,
            labelHashSet& levelCells
        );

        // Parallel send
        void pWrite(const label toID, const label& data);

        template <class Type>
        void pWrite(const label toID, const List<Type>& data);

        // Parallel receive
        void pRead(const label fromID, label& data);

        template <class Type>
        void pRead(const label fromID, List<Type>& data);

        // Read edge-modification options from the dictionary
        void readEdgeOptions();

        // Does the mesh perform edge-modification?
        bool edgeModification();

        // Send and receive patch sub-meshes for coupled patches
        void sendAndRecvCoupledMeshes();

        // Handle topology changes for coupled patches
        void handleCoupledPatches();

        // Build patch sub-meshes for coupled patches
        void buildPatchSubMeshes();

        // Build patch sub-mesh for a specified coupled patch
        void buildPatchSubMesh
        (
            const label patchID,
            patchSubMesh& subMesh
        );

        // Build coupled maps
        void buildCoupledMaps(bool localOnly);

        // Wait for buffer transfer completion.
        void waitForBuffers();

        // Synchronize and exit for parallel runs.
        void synchronizeAndExit();

        // Private internal access to member data for static members
        scalar ratioMin() { return ratioMin_; }
        scalar ratioMax() { return ratioMax_; }
        scalar growthFactor() { return growthFactor_; }

public:

    // Declare the name of the class and its debug switch
        ClassName("dynamicTopoFvMesh");

    // Constructors

        //- Construct from existing IOobject
        explicit dynamicTopoFvMesh(const IOobject&);

    // Destructor

        virtual ~dynamicTopoFvMesh();

    // Member Functions

        // Return reference to the dictionary
        const dictionary& dynamicMeshDict() const;

        // Return reference to the edge mesh
        eMesh& EdgeMesh();

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

        // Set curvature information for a particular patch
        void setCurvatureField
        (
            const label pID,
            const scalarField& field
        );

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
