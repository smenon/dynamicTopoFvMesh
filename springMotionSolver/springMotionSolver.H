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
    fvMotionSolver

Description
    Virtual base class for spring-analogy mesh motion solver.

SourceFiles
    springMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef springMotionSolver_H
#define springMotionSolver_H

#include "motionSolver.H"

#include "edgeList.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class springMotionSolver
:
    public motionSolver
{

protected:

    // Protected data

    //- The polyMesh to be moved
    const polyMesh& polyMesh_;
    
    // Reference points
    mutable pointField refPoints_;
    
    // Edge gradient field
    scalarField gradEdge_;
    
    // Component being solved for (0,1 or 2)
    label cmpt_;
    
    // ID of the patch being smoothed
    label pID_;
    
    // Variables used in the CG solver
    scalarField b_;
    scalarField x_;
    scalarField p_;
    scalarField r_;
    scalarField w_;
    
    // Protected member functions
    
    // Sparse matrix-vector multiply
    void A(const scalarField&, scalarField&);
    
    // Basic CG solver
    label CG(const scalarField&, scalarField&);
    
    // Initialize fields for the CG solver
    void initCG(label);
    
    // Apply boundary conditions
    void applyBCs(scalarField &field);

public:

    //- Runtime type information
    TypeName("springMotionSolver");


    // Constructors

        //- Construct from polyMesh
        springMotionSolver(const polyMesh& mesh);

        springMotionSolver(const polyMesh& mesh, Istream& msData);

    // Destructor

        virtual ~springMotionSolver();


    // Member Functions

        //- Return reference to the fvMesh to be moved
        const polyMesh& mesh() const
        {
            return polyMesh_;
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
