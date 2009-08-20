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

    //- Number of points
    unsigned long nPoints_;

    //- Number of elements
    unsigned long nCells_;

    //- Switch to toggle surface smoothing
    Switch surfaceSmoothing_;

    // Specify tolerance for the CG solver
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
    //Mesquite::TerminationCriterion untInner_;
    //Mesquite::TerminationCriterion untOuter_;

    //- Specific to mesh untangling
    //autoPtr<Mesquite::UntangleBetaQualityMetric> untangleMetric_;
    //autoPtr<Mesquite::LPtoPTemplate> untangleFunc_;
    //autoPtr<Mesquite::ConjugateGradient> untangleGlobal_;

    // Specify patches on which slip BCs are applied
    labelHashSet slipPatchIDs_;

    //- Data specific to Laplacian surface smoothing
    labelList pIDs_;
    labelList offsets_;
    List<vectorField> pNormals_;
    List<vectorField> gradEdgeV_;
    List<vectorField> localPts_;

    //- Explicit patch coupling for surface smoothing
    Map<label> patchCoupling_;

    // Mapping for locally coupled patches
    List<Map<label> > masterToSlave_;

    vectorField bV_;
    vectorField xV_;
    vectorField pV_;
    vectorField rV_;
    vectorField wV_;

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

    // Private member function to perform Laplacian surface smoothing
    void smoothSurfaces();

    // Prepare point-normals with updated point positions
    void preparePointNormals();

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
