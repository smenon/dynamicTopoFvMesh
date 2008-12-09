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

#include <iomanip>
#include "HashList.H"
#include "fvCFD.H"
#include "motionSolver.H"
#include "mapPolyMesh.H"
#include "MapFvFields.H"
#include "tetDecompositionMotionSolver.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"
#include "MeshObject.H"

// Multithreading support
#include "multiThreader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

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
        
        //- Motion characteristics [2D/3D]
        Switch twoDMotion_;        
                
        //- Edge bisection/collapse option
        Switch edgeModification_;
        
        //- Sub-dictionary which specifies fixed length-scales for patches
        dictionary fixedLengthScalePatches_;
        
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
        
        //- Size information
        labelList       oldPatchSizes_, patchSizes_;
        labelList       oldPatchStarts_, patchStarts_;
        labelList       oldPatchNMeshPoints_, patchNMeshPoints_;
        label           nOldPoints_, nPoints_;
        label           nOldEdges_,  nEdges_;
        label           nOldFaces_,  nFaces_;
        label           nOldCells_,  nCells_;
        label           nOldInternalFaces_, nInternalFaces_;
        
        //- Old-to-New
        labelList       reversePointMap_;
        labelList       reverseFaceMap_;
        labelList       reverseCellMap_;
        
        //- New-to-Old
        labelList       pointMap_;
        labelList       faceMap_;
        labelList       cellMap_;
        
        //- Maps for the renumbering of added entities
        Map<objectMap>  cellsFromCells_;
        Map<label>      addedPointRenumbering_; 
        Map<label>      addedCellRenumbering_;
        Map<label>      boundaryPatches_;
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
        label               bisectInteriorFace_;

        //- Edge-face connectivity list
        HashList<label>     edgeStartFace_;

        //- Maximum number of tets around an edge
        //  for 3D mesh flipping
        label               maxTetsPerEdge_;
        
        //- Support for multithreading
        autoPtr<multiThreader> threader_;
        
        //- Structure used for multiThreading
        struct topoMeshStruct
        {
            // Pointer to this class
            dynamicTopoFvMesh *mesh;
            
            // Total number of threads and thread ID
            label nThreads, threadID;
            
            // Entity start positions & sizes
            label pointStart, edgeStart, faceStart, cellStart;
            label pointSize, edgeSize, faceSize, cellSize;
            
            // Default constructor
            topoMeshStruct()
            {
                mesh = 0;
                nThreads = threadID = -1;
                pointStart = edgeStart = faceStart = cellStart = 0;
                pointSize  = edgeSize = faceSize  = cellSize  = 0;
            }
        } *structPtr_;
                
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
        
        // Initialize the edge-length field
        void initEdgeLengths();
        
        // Initialize mesh edges and related connectivity lists
        void initEdges();
        
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
        
        // Method to determine the boundary patch index for a given face
        inline label whichPatch(const label& index) const;
        
        // Method for the bisection of a quad-face in 2D
        void bisectQuadFace
        (
            const label findex, 
            face& thisFace
        );
        
        // Method for the collapse of a quad-face in 2D
        // Returns a boolean value indicating whether the collapse was valid
        bool collapseQuadFace
        (
            const label findex, 
            face& thisFace
        );
        
        // Check if the boundary face is adjacent to a sliver-cell,
        // and remove it by a two-step bisection/collapse operation.
        bool remove2DSliver
        (
            const label findex, 
            face& thisFace
        );        
        
        // Method to find the interior/boundary faces
        // for an input quad-face and adjacent triangle-prism cell.
        inline void findPrismFaces
        (
            const label& findex,
            const cell& c,
            faceList& bdyf,
            labelList& bidx,
            faceList& intf,
            labelList& iidx
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
        void findIsolatedPoint
        (
            const face& f, 
            const edge& e, 
            label& ptIndex, 
            label& nextPtIndex
        );

        // Method to determine whether the face contains the edge
        inline bool edgeOnFace
        (
            const face& f,
            const edge& e
        );
        
        // Method to replace a face-label in a given cell
        inline void replaceFaceLabel
        (
            const label& original, 
            const label& replacement, 
            cell& c
        );
        
        // Method to replace a point-label in a given face
        inline void replacePointLabel
        (
            const label& original, 
            const label& replacement, 
            face& f
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

        // Remove the specified edge from the mesh
        void removeEdge(const label index);
        
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
        // Returns whether the edge lies on a boundary 
        inline bool constructVertexRing
        (
            const label eIndex,
            const edge& edgeToCheck,
            DynamicList<label>& hullVertices,
            scalar& minQuality
        );

        // Utility method to fill the dynamic programming tables
        // Returns the number of triangulations
        inline label fillTables
        (
            const edge& edgeToCheck,
            const DynamicList<label>& hullVertices,
            const scalar minQuality,
            scalarListList& Q,
            labelListList& K
        );

        // Print out tables for debugging
        void printTables
        (
            const label m,
            const scalarListList& Q,
            const labelListList& K
        );

        // Remove the edge according to the swap sequence
        void removeEdgeFlips
        (
            const label m,
            const edge& edgeToCheck,
            const labelListList& K,
            const DynamicList<label>& hullVertices,
            labelListList& triangulations
        );

        // Extract triangulations from the programming table
        void extractTriangulation
        (
            const label i,
            const label j,
            const labelListList& K,
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
        
        // Reorder points after a topology change
        void reOrderPoints(pointField& points);
        
        // Reorder faces in upper-triangular order after a topology change
        void reOrderFaces
        (
            faceList& faces, 
            labelList& owner, 
            labelList& neighbour
        );
        
        // Reorder & renumber cells with bandwidth reduction after a topology change
        void reOrderCells();
        
        // Reorder the mesh in upper-triangular order, and generate mapping information
        void reOrderMesh
        (
            pointField& points, 
            faceList& faces, 
            labelList& owner, 
            labelList& neighbour
        );
        
        // Calculate the edge length-scale for the mesh
        void calculateLengthScale();
        
        // Return the appropriate length-scale for boundary face
        scalar boundaryLengthScale(const label faceIndex);
        
        // Prepare thread structures
        void prepareThreads(const label numThreads);
        
        // Does the mesh perform edge-modification?
        bool edgeModification();
        
        // Private internal access to member data for static members
        HashList<point>& meshPoints() { return meshPoints_; }
        HashList<edge>&  meshEdges() { return edges_; }
        HashList<face>&  meshFaces() { return faces_; }
        HashList<cell>&  meshCells() { return cells_; }
        HashList<label>& meshOwner() { return owner_; }
        HashList<label>& meshNeighbour() { return neighbour_; }
        bool& topoChangeFlag() { return topoChangeFlag_; }
        
        HashList<scalar>& meshLengthScale() { return lengthScale_; }
        HashList<edge>& edgeToWatch() { return edgeToWatch_; }
        const scalar ratioMin() { return ratioMin_; }
        const scalar ratioMax() { return ratioMax_; }
        const scalar growthFactor() { return growthFactor_; }
        const label bisectInteriorFace() { return bisectInteriorFace_; }
        const label nOldCells() { return nOldCells_; }
        const label maxTetsPerEdge() { return maxTetsPerEdge_; }
        
        Map<objectMap>& cellsFromCells() { return cellsFromCells_; }
        Map<label>& cellParents() { return cellParents_; }
        
public:
    
    // Declare the name of the class and its debug switch
        ClassName("dynamicTopoFvMesh");
    
    // Constructors

        //- Construct from existing IOobject
        explicit dynamicTopoFvMesh(const IOobject&);        

    // Destructor

        virtual ~dynamicTopoFvMesh();

    // Member Functions

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
        const autoPtr<mapPolyMesh> meshMap();
        
        // Return old cell-centre information (prior to a topology change)
        const vectorField& oldCellCentres() const;
        
        // Return mesh length-scale values
        tmp<volScalarField> lengthScale();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "dynamicTopoFvMesh.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
