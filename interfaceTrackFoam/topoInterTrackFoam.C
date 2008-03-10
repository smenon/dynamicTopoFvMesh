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
    argc = 3; 
    argv[1] = new char[100];
    argv[2] = new char[100];
    strcpy(argv[1],"/home/smenon/OpenFOAM/smenon-1.4.1-dev/run");
    strcpy(argv[2],"inkJet");  
    
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicMesh.H"
#   include "initContinuityErrs.H"
#   include "initTotalVolume.H"    
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    mesh.debug = true;
    
    while (runTime.run())
    {
#       include "readPISOControls.H"
#       include "readFreeSurfaceControls.H"        
#       include "readTimeControls.H"
#       include "checkTotalVolume.H"
#       include "CourantNo.H"

        // Make the fluxes absolute
        fvc::makeAbsolute(phi, rho, U);

#       include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;       
        
        bool meshChanged = mesh.updateTopology();

        if (meshChanged)
        {
            interface.updateMesh(mesh.meshMap()); 
#           include "checkTotalVolume.H"
#           include "correctPhi.H"
#           include "CourantNo.H"
        }

        // Update the free-surface
        interface.updateDisplacementDirections();

        interface.moveMeshPointsForOldFreeSurfDisplacement();
        
#       include "volContinuity.H"        

        if (smooth) interface.smooth();

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
            
#           include "freeSurfaceContinuityErrs.H"            
            
            Info << endl;
        }
        
#       include "meshInfo.H"        
        
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
