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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    topoInterTrackFoam

Description
    Incompressible laminar CFD code for interface between fluid phases using 
    a dynamic mesh.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicTopoFvMesh.H"
#include "freeSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    // For debug purposes, assign <root> and <case>
    /*
    argc = 3; 
    argv[1] = new char[100];
    argv[2] = new char[100];
    strcpy(argv[1],"/home/smenon/OpenFOAM/smenon-1.4.1-dev/run");
    strcpy(argv[2],"smaller");  
    */
    
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicMesh.H"
#   include "initContinuityErrs.H"
#   include "initTotalVolume.H"    
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    //polyMesh::debug = true;
    //mesh.debug = true;
    
    while (runTime.run())
    {
#       include "readPISOControls.H"
#       include "readFreeSurfaceControls.H"        
#       include "readTimeControls.H"
#       include "checkTotalVolume.H"
#       include "CourantNo.H"

#       include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;       

        // Update free-surface displacement directions
        interface.updateDisplacementDirections();        

        // Set boundary conditions for the motionSolver and solve for mesh-motion
        interface.moveMeshPointsForOldFreeSurfDisplacement();
        mesh.setMotionBC(interface.aPatchID()) = interface.totalDisplacement();
        mesh.updateMotion();  
        
#       include "volContinuity.H"
        
        if (smooth) interface.smooth();        
        
        // Update the free-surface
        interface.movePoints();

        Info<< "\nMax surface Courant Number = " 
            << interface.maxCourantNumber() << endl << endl;        
        
        for (int corr=0; corr<nOuterCorr; corr++)
        {        
            // Update interface BCs
            interface.updateBoundaryConditions();            
            
            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, rho, U);

#           include "UEqn.H"

            rUA = 1.0/UEqn.A();            
            
            // --- PISO loop

            for (int corr=0; corr<nCorr; corr++)
            {
                U = rUA*UEqn.H();
                
                phi = (fvc::interpolate(U) & mesh.Sf());
                     //+ fvc::ddtPhiCorr(rUA, U, phi);

                //adjustPhi(phi, U, p);

                for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                {
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rUA, p) == fvc::div(phi)
                    );

                    pEqn.setReference(pRefCell, pRefValue);

                    if (corr == nCorr - 1 && nonOrth == nNonOrthCorr)
                    {
                        pEqn.solve(mesh.solver(p.name() + "Final"));
                    }
                    else
                    {
                        pEqn.solve(mesh.solver(p.name()));
                    }

                    if (nonOrth == nNonOrthCorr)
                    {
                        phi -= pEqn.flux();
                    }
                }

#               include "continuityErrs.H"

                U -= rUA*fvc::grad(p);
                U.correctBoundaryConditions();
            }
            
            // Move the free-surface
            interface.movePoints();
            
            // Update motion fluxes - fluxes are still relative at this point
            phiNet = fvc::interpolate(rho)*phi;            

#           include "hEqn.H"            
#           include "freeSurfaceContinuityErrs.H"            
            
            Info << endl;
        }       

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        
        // Make the fluxes absolute
        fvc::makeAbsolute(phi, rho, U);  
        
        bool meshChanged = mesh.updateTopology(); 
        
        if (meshChanged)
        {
#           include "checkTotalVolume.H"
            /*
            // Obtain interpolated fluxes from the mesh, and reconstruct U
            forAll(phi.internalField(),faceI) {
                phi.internalField()[faceI] = mesh.interpolatedPhi()[faceI];
            }
            forAll(mesh.boundaryMesh(),patchI) {
                label start=mesh.boundaryMesh()[patchI].start();
                forAll(phi.boundaryField()[patchI],faceI) {
                    phi.boundaryField()[patchI][faceI] = mesh.interpolatedPhi()[start+faceI];
                }
            }
            U = fvc::reconstruct(phi); 
            */
            interface.updateMesh(mesh.meshMap());
#           include "correctPhi.H"            
#           include "CourantNo.H"
        }        

        runTime.write(); 
#       include "meshInfo.H"        
    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
