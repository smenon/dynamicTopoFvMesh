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
    a dynamic mesh, including non-Newtonian effects.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "viscosityModel.H"
#include "dynamicTopoFvMesh.H"
#include "fluidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    // For debug purposes, assign <root> and <case>
    argc = 3; 
    argv[1] = new char[100];
    argv[2] = new char[100];
    strcpy(argv[1],"/home/smenon/OpenFOAM/smenon-1.4.1-dev/run");
    strcpy(argv[2],"ligament");  

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicMesh.H"
#   include "initContinuityErrs.H"
#   include "initTotalVolume.H"    
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    //polyMesh::debug = true;
    mesh.debug = true;
    
    while (runTime.run())
    {
#       include "readPISOControls.H"       
#       include "readTimeControls.H"
#       include "checkTotalVolume.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;       

        // Update free-surface displacement directions
        interface.updateDisplacementDirections();      

        // Set boundary conditions for the motionSolver and solve for mesh-motion
        interface.restorePosition();
        mesh.setMotionBC(interface.patchID(), interface.displacement());
        mesh.updateMotion();  
        
#       include "volContinuity.H"         
        
        for (int corr=0; corr<nOuterCorr; corr++)
        {        
            // Update boundary conditions on velocity and pressure
            interface.updateBoundaryConditions();         
            
            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, U);

#           include "UEqn.H"

            volScalarField rUA = 1.0/UEqn.A();            

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

            // Update the interface with the fluid velocity
            interface.movePoints();            
            
#           include "freeSurfaceContinuityErrs.H"  
            
            Info << endl;
        }
        
        bool meshChanged = mesh.updateTopology(); 
        
        if (meshChanged)
        {      
#           include "checkTotalVolume.H"
            phi = linearInterpolate(U) & mesh.Sf();
#           include "correctPhi.H"            
#           include "CourantNo.H"
            interface.updateMesh(mesh.meshMap());            
        }     
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;        
 
        runTime.write(); 
#       include "meshInfo.H"        
    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
