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
    springMotionSolver

Description
    Class for the spring-analogy [2D] / ball-vertex [3D] mesh motion solver.

SourceFiles
    springMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef springMotionSolver_H
#define springMotionSolver_H

#include "motionSolver.H"
#include "faceList.H"
#include "labelHashSet.H"
#include "pointIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{    

// Class forward declarations
class polyMesh;
class dynamicTopoFvMesh;

/*---------------------------------------------------------------------------*\
                           Class fvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class springMotionSolver
:
    public motionSolver
{

protected:

    // Protected data

    //- The mesh to be moved
    const dynamicTopoFvMesh& Mesh_;
    
    // Reference points
    mutable pointIOField refPoints_;
    
    // Edge gradient field
    scalarField gradEdge_;

    // Edge stiffness field
    scalarField stiffness_;

    // Point-to-opposing face connectivity
    List<faceList> pointFaces_;

    // Stiffness of edge between point and virtual-point
    scalarListList pfStiffness_;

    // Interpolation coefficients for the virtual point
    scalarListList xi_;
    scalarListList eta_;

    // Specify patches on which slip BCs are applied
    labelHashSet slipPatchIDs_;
    
    // Component being solved for (0,1 or 2)
    label cmpt_;
    
    // ID of the patch being smoothed
    label pID_;

    // Y-direction constraint
    scalar fixedY_;
    
    // Variables used in the CG solver
    scalarField b_;
    scalarField x_;
    scalarField p_;
    scalarField r_;
    scalarField w_;
    
    // Protected member functions
    
    // Sparse matrix-vector multiply
    void A(const scalarField& p, scalarField& w);

    // Matrix multiply [2D]
    void computeA2D(const scalarField& p, scalarField& w);

    // Matrix multiply [3D]
    void computeA3D(const scalarField& p, scalarField& w);
    
    // Basic CG solver
    label CG(const scalarField&, scalarField&);
    
    // Initialize fields for the CG solver
    void initCG(label);
    
    // Set Dirichlet conditions for the solution field (if any)
    void setDirichlet(scalarField &x);

    // Apply boundary conditions
    void applyBCs(scalarField &field);

    // Compute point-to-opposing face connectivity
    void computePointFaces();
    
public:

    //- Runtime type information
    ClassName("springMotionSolver");


    // Constructors

        //- Construct from polyMesh
        springMotionSolver(const polyMesh& mesh);

        springMotionSolver(const polyMesh& mesh, Istream& msData);

    // Destructor

        virtual ~springMotionSolver();


    // Member Functions

        //- Return reference to the mesh to be moved
        const dynamicTopoFvMesh& mesh() const
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
