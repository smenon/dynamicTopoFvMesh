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

// Tetrahedral mesh-quality metrics
#include "tetQuality.H"
#include "Knupp.H"
#include "Dihedral.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
    
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
        
        //- Quality metric for tetrahedra in 3D
        autoPtr<tetQuality> tetMetric_;        
        
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
        
        //- Edge-face connectivity lists
        HashList<SLList<label> > faceEdges_;
        HashList<SLList<label> > edgeFaces_;
        
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
            label pointStart, faceStart, cellStart;
            label pointSize, faceSize, cellSize;
            
            // Default constructor
            topoMeshStruct()
            {
                mesh = 0;
                nThreads = threadID = -1;
                pointStart = faceStart = cellStart = 0;
                pointSize  = faceSize  = cellSize  = 0;
            }
        };
        
        // Structure array pointer
        topoMeshStruct *structPtr_;
        
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
        
        // Find the area of a triangle face.
        inline scalar triFaceArea(const face&);
        
        // Initialize the edge-length field
        void initEdgeLengths();
        
        // Initialize mesh edges and related connectivity lists
        void initEdges();
        
        // MultiThreaded topology modifier [2D]
        static threadReturnType threadedTopoModifier2D(void *arg);
        
        // 2D Edge-swapping engine
        void swap2DEdges(const topoMeshStruct *);
        
        // 2D Edge-bisection/collapse engine
        void edgeBisectCollapse2D(const topoMeshStruct *);
        
        // Method to determine the boundary patch index for a given face
        inline label whichPatch(const label&) const;
        
        // Method for the bisection of a quad-face in 2D
        void bisectQuadFace(const label, face&);
        
        // Method for the collapse of a quad-face in 2D
        // Returns a boolean value indicating whether the collapse was valid
        bool collapseQuadFace(const label, face&);
        
        // Check if the boundary face is adjacent to a sliver-cell,
        // and remove it by a two-step bisection/collapse operation.
        bool remove2DSliver(const label, face&);        
        
        // Method to find the interior/boundary faces
        // for an input quad-face and adjacent triangle-prism cell.
        void findPrismFaces
        (
            const label& findex,
            const cell& c,
            faceList& bdyf,
            labelList& bidx,
            faceList& intf,
            labelList& iidx
        );
        
        // Method to find the interior/boundary faces
        // for an input tri-face and adjacent tet/pyramid cell.
        void findTetPyramidFaces
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
        bool findCommonEdge(const face&, const face&, edge&);
        
        // Method to find the isolated point on a triangular face
        // that doesn't lie on the specified edge
        void findIsolatedPoint(const face&, const edge&, label&, label&);
        
        // Method to replace a face-label in a given cell
        inline void replaceFaceLabel(const label&, const label&, cell&);
        
        // Method to replace a point-label in a given face
        inline void replacePointLabel(const label&, const label&, face&);
        
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
        void removeFace(const label);
        
        // Utility method to build a hull of faces connected to the edge
        // This will also determine whether the edge lies on a boundary
        bool constructPrismHull
        (
            const edge& edgeToCheck,
            const label startFaceIndex,
            DynamicList<label>& hullFaces,
            DynamicList<label>& hullCells,
            DynamicList<label>& hullTriFaces,
            bool requiresTriFaces
        );
        
        // Reorder points after a topology change
        void reOrderPoints(pointField&);
        
        // Reorder faces in upper-triangular order after a topology change
        void reOrderFaces(faceList&, labelList&, labelList&);
        
        // Reorder & renumber cells with bandwidth reduction after a topology change
        void reOrderCells();
        
        // Reorder the mesh in upper-triangular order, and generate mapping information
        void reOrderMesh(pointField&, faceList&, labelList&, labelList&);
        
        // Calculate the edge length-scale for the mesh
        void calculateLengthScale();
        
        // Return the appropriate length-scale for boundary face
        scalar boundaryLengthScale(const label);
        
        // Prepare thread structures
        void prepareThreads();
        
        // Does the mesh perform edge-modification?
        bool edgeModification();
        
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
        
        // Access
        void setMotionBC(const label&, const vectorField&);
        
        // Return the mesh-mapper
        const autoPtr<mapPolyMesh> meshMap();
        
        // Return old cell-centre information (prior to a topology change)
        const vectorField& oldCellCentres() const;
        
        HashList<scalar>& lengthDensity()
        {
            return lengthScale_;
        }           
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "dynamicTopoFvMesh.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
