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
#include "boolList.H"
#include "pointIOField.H"
#include "scalarMatrix.H"
#include "labelHashSet.H"
#include "tensor2D.H"

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
                       Class mesquiteSmoother Declaration
\*---------------------------------------------------------------------------*/

class mesquiteSmoother
:
    public motionSolver
{

protected:

    // Protected data

    //- The mesh to be moved
    const polyMesh& Mesh_;

    //- Number of points
    unsigned long nPoints_;

    //- Number of elements
    unsigned long nCells_;

    //- Switch to toggle surface smoothing
    Switch surfaceSmoothing_;

    //- Switch to prioritize surface smoothing
    //  for certain vertices
    Switch prioritizeVertices_;

    // Specify tolerance for the solver
    scalar tolerance_;

    //- Specify multiple mesh-motion sweeps
    label nSweeps_;

    //- Specify the surface-smoothing interval
    label surfInterval_;

    //- Vertex coordinate array passed into Mesquite
    mutable double* vtxCoords_;

    //- Connectivity array passed into mesquite
    mutable unsigned long* cellToNode_;

    //- Flag array for vertices (fixed/free)
    mutable int* fixFlags_;

    //- Reference points
    mutable pointIOField refPoints_;

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

    // Specify patches on which surface smoothing is applied
    labelHashSet slipPatchIDs_;

    //- Data specific to surface smoothing
    label patchIndex_;
    label pointIndex_;
    labelList pIDs_;
    List<vectorField> pNormals_;
    List<vectorField> localPts_;
    List<labelListList> jPoints_;
    List<tensorField> T_;
    List<List<scalarField> > a_;
    boolListList pFlags_;
    List<scalar> meanFnValues_;
    List<scalarField> pFnValues_;

    //- Explicit patch coupling for surface smoothing
    Map<label> patchCoupling_;

    // Mapping for locally coupled patches
    List<Map<label> > masterToSlave_;

    // Private member function to compute connectivity data
    void initArrays();

    // Private member functions to read dictionary options
    void readOptions();

    // Private member function to perform surface smoothing
    void smoothSurfaces();

    // Solve for surface point positions
    void solveForSurfacePoints();

    // Optimize the point-position in local coordinates.
    void optimizePoint();

    // Evaluate the function at a point.
    void evaluateFunction
    (
        const label type,
        const vector& x,
        scalar& fn,
        vector& fnGrad
    );

    // Evaluate the constraint for a given point.
    void evaluateConstraint
    (
        const label type,
        const vector& x,
        scalar& cstr,
        vector& grad
    );

    // Solve the QP sub-problem.
    void QPsolver
    (
        const tensor& H,
        const vector& f,
        const vector& Ac,
        const scalar& bc,
        vector& x,
        scalar& lambda,
        const bool normalize = true
    );

    // QR Decomposition using LAPACK (for vectors)
    void qrDecompose
    (
        const vector& a,
        tensor& Q,
        vector& R,
        tensor& Z
    );

    // Pseudo-inverse using LAPACK (for vectors)
    vector pInv(const vector& a);

    // Pseudo-inverse using LAPACK (for 2D tensors)
    tensor2D pInv(const tensor2D& a);

    // Compute a search direction vector
    vector computeSearchDirection
    (
        const tensor& Z,
        const tensor& H,
        const vector& gf,
        label& type
    );

    // Prepare a priority list of surface vertices for smoothing.
    void preparePriorityList();

    // Prepare point-normals with updated point positions
    void preparePointNormals();

    // Utility method to prepare the least-squares matrix
    void prepareLeastSquaresMatrix
    (
        const label nPoints,
        const scalarListList& lsqCoeffs,
        const scalarField& lsqWeights,
        const scalarField& lsqSource,
        scalarMatrix& lsqMatrix,
        scalarField& source
    );

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
