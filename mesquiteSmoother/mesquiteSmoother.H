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
    mesquiteSmoother

Description
    Thin interface to the Mesquite Mesh Improvement library.

SourceFiles
    mesquiteSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef mesquiteSmoother_H
#define mesquiteSmoother_H

#include "motionSolver.H"
#include "Time.H"
#include "Map.H"
#include "Switch.H"
#include "pointIOField.H"
#include "labelHashSet.H"

#include "lbfgs.h"

// Have gcc ignore certain warnings while including mesquite headers
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wempty-body"
#include "Mesquite_all_headers.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class fvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class mesquiteSmoother
:
    public motionSolver
{

protected:

    // Protected data

    //- The mesh to be moved
    const polyMesh& Mesh_;

    //- Mesh characteristics [2D/3D]
    Switch twoDMesh_;

    //- Number of points
    unsigned long nPoints_;

    //- Number of elements
    unsigned long nCells_;

    //- Number of auxiliary points
    unsigned long nAuxPoints_;

    //- Number of auxiliary elements
    unsigned long nAuxCells_;

    //- Switch to toggle surface smoothing
    Switch surfaceSmoothing_;

    //- Switch to toggle volume correction
    Switch volumeCorrection_;

    //- Specify tolerance for volume correction
    scalar volCorrTolerance_;

    //- Specify max volume correction iterations
    label volCorrMaxIter_;

    // Specify tolerance for the CG solver
    scalar tolerance_;

    //- Specify multiple mesh-motion sweeps
    label nSweeps_;

    //- Specify the surface-smoothing interval
    label surfInterval_;

    //- Specify a relaxation factor for surface-smoothing
    scalar relax_;

    //- Vertex coordinate array passed into Mesquite
    mutable double* vtxCoords_;

    //- Connectivity array passed into mesquite
    mutable unsigned long* cellToNode_;

    //- Flag array for vertices (fixed/free)
    mutable int* fixFlags_;

    //- Reference points
    mutable pointIOField refPoints_;

    //- Original points prior to surface-smoothing
    mutable pointField origPoints_;

    //- The quality metric
    word qMetric_;

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
    labelList offsets_;
    List<vectorField> pNormals_;
    List<vectorField> gradEdgeV_;
    List<vectorField> localPts_;

    //- Explicit patch coupling for surface smoothing
    Map<label> patchCoupling_;

    vectorField bV_;
    vectorField xV_;
    vectorField pV_;
    vectorField rV_;
    vectorField wV_;

    scalar oldVolume_;

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

    // Apply boundary conditions
    void applyBCs(vectorField &field);

    // Private member function to compute connectivity data
    void initArrays();

    // Private member functions to read dictionary options
    void readOptions();

    // Apply fixed-value boundary conditions, if any.
    void applyFixedValuePatches();

    // Private member function to perform Laplacian surface smoothing
    void smoothSurfaces();

    // Compute volume of a tetrahedral cell
    scalar tetVolume
    (
        const label cIndex,
        const pointField& pField
    );

    // Compute quality of a tetrahedral cell
    scalar tetQuality
    (
        const label cIndex,
        const pointField& pField
    );

    // Private member function to check for invalid
    // cells and correct if necessary.
    void correctInvalidCells();

    // Private member function to correct for volume changes
    // from surface smoothing.
    void correctGlobalVolume();

    // Enforce cylindrical constraints for slip-patches
    void enforceCylindricalConstraints();

    // Write a particular point out for post-processing
    void writePoint
    (
        const label pointIndex,
        const labelList& jList
    );

    // Output a list of points / cells as a VTK file.
    void writeVTK
    (
        const word& name,
        const polyMesh& mesh,
        const pointField& points,
        const labelListList& cellPoints
    );

    // Prepare point-normals with updated point positions
    void preparePointNormals();

    // Class used in callback
    class optInfo
    {
        // Reference to the smoother
        const mesquiteSmoother& msRef_;

        // Beta value for optimization
        scalar beta_;

        // Jacobian connectivity list
        const labelList& jList_;

        // Points for reference
        const pointField& refPts_;

    public:

        optInfo
        (
            const mesquiteSmoother& msRef,
            const scalar& beta,
            const labelList& jList,
            const pointField& refPoints
        )
        :
            msRef_(msRef),
            beta_(beta),
            jList_(jList),
            refPts_(refPoints)
        {}

        const mesquiteSmoother& motionSolver()
        {
            return msRef_;
        }

        const labelList& orderedPoints()
        {
            return jList_;
        }

        const pointField& referencePoints()
        {
            return refPts_;
        }

        scalar beta()
        {
            return beta_;
        }
    };

    // Static function for callback
    static lbfgsfloatval_t _evaluate
    (
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
    );

    // Actual function evaulation routine
    lbfgsfloatval_t evaluate
    (
        const labelList& jList,
        const pointField& refPoints,
        const scalar beta,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g
    ) const;

    // Utility method to check validity of cells connected to a point.
    bool checkValidity
    (
        const vector& x,
        const labelList& jList,
        scalar& beta
    );

    // Clear out addressing
    void clearOut();

public:

    //- Runtime type information
    ClassName("mesquiteSmoother");


    // Constructors

        //- Construct from polyMesh
        mesquiteSmoother(const polyMesh& mesh);

        mesquiteSmoother(const polyMesh& mesh, Istream& msData);

    // Destructor

        virtual ~mesquiteSmoother();


    // Member Functions

        //- Return reference to the mesh to be moved
        const polyMesh& mesh() const
        {
            return Mesh_;
        }

        //- Provide new points for motion.  Solves for motion
        virtual tmp<pointField> newPoints();

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
