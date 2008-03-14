/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-2007 Zeljko Tukovic
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

Application
    icoDyMFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids
    with dynamic mesh.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicTopoFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // For debug purposes, assign <root> and <case>
    argc = 3; 
    argv[1] = new char[100];
    argv[2] = new char[100];
    strcpy(argv[1],"/home/smenon/OpenFOAM/smenon-1.4.1-dev/run");
    strcpy(argv[2],"cylinder");     
    
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicMesh.H"
#   include "initContinuityErrs.H"
#   include "initTotalVolume.H"
#   include "createNuFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    mesh.debug = true;
    //polyMesh::debug = true;
    
    while (runTime.run())
    {
#       include "readPISOControls.H"
#       include "readTimeControls.H"
#       include "checkTotalVolume.H"
#       include "CourantNo.H"
        
        // Make the fluxes absolute
        fvc::makeAbsolute(phi, U);

#       include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Assign boundary conditions to the motion solver
        label patchID = -1;
        forAll (mesh.boundary(), patchI)
        {
            if(mesh.boundary()[patchI].name() == "topWall") {
                patchID = patchI;
                break;
            }
        }         
        mesh.boundaryDisplacementPatch(patchID) = vector(0,-0.05,0)*mesh.time().deltaT().value();          
        
        bool meshChanged = mesh.updateTopology(); 
        
        if (meshChanged)
        {
#           include "checkTotalVolume.H"
            // Obtain interpolated fluxes from the mesh, and reconstruct U
            forAll(phi.internalField(),faceI) {
                phi.internalField()[faceI] = mesh.interpolatedFluxes()[faceI];
            }
            forAll(mesh.boundaryMesh(),patchI) {
                label start=mesh.boundaryMesh()[patchI].start();
                forAll(phi.boundaryField()[patchI],faceI)
                    phi.boundaryField()[patchI][faceI] = mesh.interpolatedFluxes()[start+faceI];
            }
            U = fvc::reconstruct(phi);           
//#           include "correctPhi.H"
#           include "CourantNo.H"
        }
        
        volScalarField divPhi = fvc::div(phi);
        divPhi.write();

        // Solve for mesh-motion
        mesh.updateMotion();         
        
#       include "volContinuity.H"        
        
        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);

#       include "NuUEqn.H"

        // --- PISO loop

        for (int corr=0; corr<nCorr; corr++)
        {
            rUA = 1.0/UEqn.A();

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

#           include "continuityErrs.H"            

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, U);

            U -= rUA*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
