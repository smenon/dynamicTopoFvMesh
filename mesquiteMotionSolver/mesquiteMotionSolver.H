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
    mesquiteMotionSolver

Description
    Thin interface to the Mesquite Mesh Improvement library.
    Also provides surface-mesh smoothing capabilities.

SourceFiles
    mesquiteMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef mesquiteMotionSolver_H
#define mesquiteMotionSolver_H

#include "motionSolver.H"
#include "Time.H"
#include "Map.H"
#include "Switch.H"
#include "edgeList.H"
#include "labelPair.H"
#include "pointFields.H"
#include "MeshObject.H"
#include "HashSet.H"

#include "mesquiteHeaders.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class mesquiteMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class mesquiteMotionSolver
:
    public motionSolver,
    public MeshObject
    <
        polyMesh,
        Foam::UpdateableMeshObject,
        mesquiteMotionSolver
    >
{

    // Private data

        //- The mesh to be moved
        const polyMesh& Mesh_;

        //- Mesh characteristics [2D/3D]
        Switch twoDMesh_;

        //- Flag to determine if arrays have been initialized
        bool arraysInitialized_;

        //- Number of points
        unsigned long nPoints_;

        //- Number of elements
        unsigned long nCells_;

        //- Polyhedra decomposition type
        //  1: for cell, 2: for face
        label decompType_;

        //- Number of points from poly decomposition
        unsigned long nDecompPoints_;

        //- List of decomposition centres
        List<labelPair> decompCellCentres_;
        List<labelPair> decompFaceCentres_;

        //- Number of auxiliary points
        labelList nAuxPoints_;

        //- Number of auxiliary elements
        labelList nAuxCells_;

        //- Switch to toggle surface smoothing
        Switch surfaceSmoothing_;

        //- Switch to toggle volume correction
        Switch volumeCorrection_;

        //- Specify tolerance for volume correction
        scalar volCorrTolerance_;

        //- Specify max volume correction iterations
        label volCorrMaxIter_;

        //- Specify tolerance for the CG solver
        scalar tolerance_;

        //- Specify multiple mesh-motion sweeps
        label nSweeps_;

        //- Specify the surface-smoothing interval
        label surfInterval_;

        //- Specify a relaxation factor for surface-smoothing
        scalar relax_;

        //- Vertex coordinate array passed into Mesquite
        mutable List<double> vtxCoords_;

        //- Connectivity array passed into mesquite
        mutable List<unsigned long> cellToNode_;

        //- Flag array for vertices (fixed/free)
        mutable List<int> fixFlags_;

        //- Array with element types
        mutable List<Mesquite::EntityTopology> mixedTypes_;

        //- Demand-driven pointCells addressing
        mutable labelListList pointCells_;

        //- Demand-drive cell offset addressing
        mutable labelList cellOffset_;

        //- Reference points
        mutable pointIOField refPoints_;

        //- Original points prior to surface-smoothing
        mutable pointField origPoints_;

        //- Dictionary for fixedValue patches
        mutable autoPtr<IOdictionary> fvpDict_;

        //- Base points for pointDisplacement boundary conditions
        mutable autoPtr<pointVectorField> basePoints_;

        //- Boundary conditions
        mutable autoPtr<pointVectorField> boundaryConditions_;

        //- Quality metric list
        wordList qMetrics_;

        //- Pointers to base element metric class
        HashTable<autoPtr<Mesquite::QualityMetric> > qMetricTable_;

        //- Pointer to base objective function
        autoPtr<Mesquite::ObjectiveFunction> objFunction_;

        //- Pointer to base optimization algorithm (vertex mover)
        autoPtr<Mesquite::VertexMover> optAlgorithm_;

        //- Termination criteria
        Mesquite::TerminationCriterion tcInner_;
        Mesquite::TerminationCriterion tcOuter_;

        //- Data specific to Laplacian surface smoothing
        labelList pIDs_;
        vectorField bdy_;
        vectorField gradEdge_;
        scalarField edgeMarker_;
        scalarField pointMarker_;
        scalarField edgeConstant_;
        List<vectorField> pNormals_;

        //- Data for the auxiliary entities
        labelList procIndices_;
        scalarField pointFractions_;
        labelListList procPointLabels_;
        List<List<labelPair> > globalProcPoints_;

        List<Map<label> > sendPointMap_, recvPointMap_;
        List<labelList> sendCellToNode_, recvCellToNode_;
        List<vectorField> sendPointBuffer_, recvPointBuffer_;

        //- Explicit patch coupling for surface smoothing
        Map<label> patchCoupling_;

        //- Old domain volume
        scalar oldVolume_;

    // Private Member Functions

        // Sparse Matrix multiply
        void A(const vectorField& p, vectorField& w);

        // Dot-product
        scalar dot(const vectorField& f1, const vectorField& f2);

        // CG solver
        label CG
        (
            const vectorField& b,
            vectorField& p,
            vectorField& r,
            vectorField& w,
            vectorField& x
        );

        // Compute the normalization factor for the matrix
        scalar normFactor
        (
            const vectorField& x,
            const vectorField& b,
            const vectorField& w,
            vectorField& tmpField
        );

        // Component-wise sumMag
        scalar cmptSumMag(const vectorField& field);

        // Transfer enumerant
        enum TransferType { Send, Recv };

        // Transfer buffers for surface point fields
        void transferBuffers
        (
            vectorField& field,
            TransferType type,
            bool fix = false
        );

        // Apply boundary conditions
        void applyBCs(vectorField &field);

        // Private member function to compute connectivity data
        void initArrays();

        // Private member function to construct parallel connectivity data
        void initParallelConnectivity();

        // Parallel send
        void parWrite(const label toID, const label& data);

        // Parallel send (for lists)
        template <class Type>
        void parWrite(const label toID, const UList<Type>& data);

        // Parallel receive
        void parRead(const label fromID, label& data);

        // Parallel receive (for lists)
        template <class Type>
        void parRead(const label fromID, UList<Type>& data);

        // Identify coupled patches
        void identifyCoupledPatches();

        // Prepare mesquite connectivity for parallel runs
        void initMesquiteParallelArrays();

        // Prepare flags for boundary layer zones
        void initBoundaryLayerZones();

        // Prepare flags for fixed zones
        void initFixedZones();

        // Fix decomposed centroids
        void fixDecomposedCentroids();

        // Compute and return cell typeMap
        labelList cellTypeMap() const;

        // Initialize pointCells addressing
        void initPointCells() const;

        // Compute and return pointCells addressing
        const labelListList& pointCells() const;

        // Compute and return cellOffset addressing
        const labelList& cellOffset() const;

        // Clear demand-driven addressing
        void clearDemandDrivenAddressing();

        // Dump current mesh for debugging
        void writeDecomposition() const;

        // Compute centroids based on up-to-date point positions
        void computeCentroids
        (
            const pointField& points,
            pointField& faceCentres,
            pointField& cellCentres
        );

        // Apply centroid positions
        void applyCentroids();

        // Copy auxiliary points to/from buffers
        void copyAuxiliaryPoints(bool firstCopy);

        // Private member functions to read dictionary options
        void readOptions();

        // Apply fixed-value boundary conditions, if any.
        void applyFixedValuePatches();

        // Apply boundary layer displacement
        void applyBoundaryLayerPatches(pointField& displacement);

        // Private member function to perform Laplacian surface smoothing
        void smoothSurfaces();

        // Compute quality of a tetrahedral cell
        scalar tetQuality
        (
            const label cIndex,
            const pointField& pField,
            bool returnVolume = false
        );

        // Private member function to check for invalid
        // cells and correct if necessary.
        void correctInvalidCells();

        // Correct global volume
        void correctGlobalVolume();

        // Enforce cylindrical constraints for slip-patches
        void enforceCylindricalConstraints();

        // Prepare point-normals with updated point positions
        void preparePointNormals();

        // Prepare non-uniform edge constants with updated point positions
        void prepareEdgeConstants(const vectorField& p);

        // Utility method to check validity of cells connected to a point.
        bool checkValidity
        (
            const vector& x,
            const labelList& jList,
            scalar& beta
        );

        // Update on topology change
        void update(const mapPolyMesh&);

public:

    //- Runtime type information
    TypeName("mesquiteMotionSolver");


    // Constructors

        //- Construct from polyMesh
        mesquiteMotionSolver(const polyMesh& mesh);

        //- Construct from polyMesh and IOdictionary
        mesquiteMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict
        );

    // Destructor

        virtual ~mesquiteMotionSolver();


    // Member Functions

        //- Return reference to the mesh to be moved
        const polyMesh& mesh() const
        {
            return Mesh_;
        }

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Move points
        virtual bool movePoints();

        //- Move points
        virtual void movePoints(const pointField&);

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
