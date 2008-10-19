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

\*---------------------------------------------------------------------------*/

#include "springMotionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(springMotionSolver, 0);
    addToRunTimeSelectionTable
    (
        motionSolver,
        springMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::springMotionSolver::springMotionSolver(const polyMesh& mesh)
:
    motionSolver(mesh),
    polyMesh_(mesh),
    refPoints_
    (
        IOobject
        (
            "refPoints",
            mesh.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),            
        mesh.points()
    ),       
    cmpt_(-1),
    pID_(-1),
    fixedY_(0.0)
{
    // Check if points on certain axes need to be maintained
    if (found("fixedY")) {
        fixedY_ = readScalar(lookup("fixedY"));
    }
}


Foam::springMotionSolver::springMotionSolver(const polyMesh& mesh, Istream& msData)
:
    motionSolver(mesh),
    polyMesh_(mesh),
    refPoints_
    (
        IOobject
        (
            "refPoints",
            mesh.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),            
        mesh.points()
    ), 
    cmpt_(-1),
    pID_(-1),
    fixedY_(0.0)
{
    // Check if points on certain axes need to be maintained
    if (found("fixedY")) 
    {
        fixedY_ = readScalar(lookup("fixedY"));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::springMotionSolver::~springMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Basic CG solver
Foam::label Foam::springMotionSolver::CG(const scalarField& b, scalarField& x)
{
    // Local variables
    scalar alpha, beta, delta_old, delta_new, tolerance = 1e-16;
    label maxIter = x.size(), iter = 0;

    // Set Dirichlet conditions on the solution field (if any)
    setDirichlet(x);
    
    A(x,w_);
    r_ = b - w_; 
    p_ = r_;    
    delta_new = gSumProd(r_,r_);  
    Info << " Initial residual: " << delta_new;
    while ( (iter < maxIter) && (delta_new > tolerance) ) 
    {
        A(p_,w_);
        alpha = delta_new / gSumProd(p_,w_);
        x  += (alpha*p_);
        r_ -= (alpha*w_);
        delta_old = delta_new;
        delta_new = gSumProd(r_,r_);
        beta = delta_new / delta_old;
        p_ = r_ + (beta*p_);
        iter++;
    }
    Info << " Final residual: " << delta_new;
    
    return iter;
}

// Sparse matrix-vector multiply
void Foam::springMotionSolver::A(const scalarField& p, scalarField& w)
{
    w = 0.0;
    
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();
    
    // Gradient (n2e)
    forAll(boundary[pID_].edges(),edgeI) 
    {
        gradEdge_[edgeI] =   p[boundary[pID_].edges()[edgeI][1]] 
                           - p[boundary[pID_].edges()[edgeI][0]];
    }
    
    // Divergence (e2n)
    forAll(boundary[pID_].edges(),edgeI) 
    {
        w[boundary[pID_].edges()[edgeI][0]] += gradEdge_[edgeI];
        w[boundary[pID_].edges()[edgeI][1]] -= gradEdge_[edgeI];
    }
    
    // Apply boundary conditions
    applyBCs(w);    
}

// Set Dirichlet conditions on the solution field (if any)
void Foam::springMotionSolver::setDirichlet(scalarField &x)
{
    // For wedges, loop through boundaries, and fix y-direction on the axis
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();
    if (boundary[pID_].type() == "wedge" && cmpt_ == 1) 
    {
        for(label i = boundary[pID_].nInternalEdges(); 
            i < boundary[pID_].edges().size(); 
            i++)
        {
            if (    
                    boundary[pID_].localPoints()[boundary[pID_].edges()[i][0]][1] < (fixedY_+SMALL)
                 && boundary[pID_].localPoints()[boundary[pID_].edges()[i][1]][1] < (fixedY_+SMALL)
               )
            {
                // This edge lies on the axis. 
                x[boundary[pID_].edges()[i][0]] = fixedY_;
                x[boundary[pID_].edges()[i][1]] = fixedY_; 
            }	    
        }    
    }
}

// Apply boundary conditions
void Foam::springMotionSolver::applyBCs(scalarField &field)
{
    // Loop through boundary edges and blank-out residuals
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();
    for(label i = boundary[pID_].nInternalEdges(); 
        i < boundary[pID_].edges().size(); 
        i++)
    {
        // If wedge patches are present, set axis nodes to slip
        if (boundary[pID_].type() == "wedge") 
        {            
            if (    
                    boundary[pID_].localPoints()[boundary[pID_].edges()[i][0]][1] < (fixedY_+SMALL)
                 && boundary[pID_].localPoints()[boundary[pID_].edges()[i][1]][1] < (fixedY_+SMALL)
               )
            {
                // This edge lies on the axis. If the Y component is being solved for,
                // blank out the residual.
                if (cmpt_ == 1) {
                    field[boundary[pID_].edges()[i][0]] = 0.0;
                    field[boundary[pID_].edges()[i][1]] = 0.0;
                }
		continue;
            }
        }
        field[boundary[pID_].edges()[i][0]] = 0.0;
        field[boundary[pID_].edges()[i][1]] = 0.0;
    }
}

// Initialize fields for the CG solver
void Foam::springMotionSolver::initCG(label nUnknowns)
{
    if (r_.size() != nUnknowns) 
    {
        b_.setSize(nUnknowns, 0.0);
        x_.setSize(nUnknowns, 0.0);
        p_.setSize(nUnknowns, 0.0);
        r_.setSize(nUnknowns, 0.0);
        w_.setSize(nUnknowns, 0.0);
    }

    const polyBoundaryMesh& boundary = mesh().boundaryMesh();    
    if (gradEdge_.size() != boundary[pID_].edges().size()) 
    {
        gradEdge_.setSize(boundary[pID_].edges().size(), 0.0);
    }
}

Foam::tmp<Foam::pointField> Foam::springMotionSolver::newPoints()
{
    solve();
    return curPoints();
}

//- Return point location obtained from the current motion field
Foam::tmp<Foam::pointField> 
Foam::springMotionSolver::curPoints() const
{      
    tmp<pointField> tcurPoints(refPoints_);    
    
    return tcurPoints;
}

void Foam::springMotionSolver::solve()
{
    if (twoDMotion()) {
        // Loop through patches, check for empty/wedge, and smooth them
        const polyBoundaryMesh& boundary = mesh().boundaryMesh();
        for(label i=0; i<boundary.size(); i++) 
        {
            if(boundary[i].type() == "wedge" || boundary[i].type() == "empty") 
            { 
                pID_ = i;
                
                // Initialize the solver variables
                const labelList& meshPts = boundary[pID_].meshPoints();
                
                initCG(meshPts.size());
                
                for(cmpt_=0; cmpt_<2; cmpt_++) 
                {
                    forAll(meshPts,pointI) 
                    {
                        x_[pointI] = refPoints_[meshPts[pointI]][cmpt_];
                    }
                    
                    Info << "Solving for component: " << cmpt_;
                    label iters = CG(b_,x_);
                    Info << " No Iterations: " << iters << endl;                    

                    forAll(meshPts,pointI) 
                    {
                        refPoints_[meshPts[pointI]][cmpt_] = x_[pointI];
                    }
                }
                
                // Correct z-direction for wedge-patches
                if(boundary[i].type() == "wedge")
                {
                    const wedgePolyPatch& wedgePatch = 
                        refCast<const wedgePolyPatch>(boundary[i]);                                                
                    
                    vector centrePoint = vector::zero;
                    forAll(meshPts,pointI) 
                    {
                        centrePoint.x() = refPoints_[meshPts[pointI]][0];
                        centrePoint.y() = refPoints_[meshPts[pointI]][1];
                        // Transform using the wedge transform tensor
                        refPoints_[meshPts[pointI]][2] 
                            = (wedgePatch.faceT()&centrePoint)[2];
                    }                    
                }
            }
        }
    }    
}

void Foam::springMotionSolver::updateMesh(const mapPolyMesh& mpm)
{    
    motionSolver::updateMesh(mpm);
    
    refPoints_.clear();
    b_.clear();
    x_.clear();
    p_.clear();
    r_.clear();
    w_.clear();  
    gradEdge_.clear();
    
    refPoints_ = polyMesh_.points();
}

// ************************************************************************* //
