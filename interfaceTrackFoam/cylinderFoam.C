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

Application
    cylinderFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids
    with dynamic mesh.

Author
    Sandeep Menon

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidInterface.H"
#include "dynamicTopoFvMesh.H"

// Mesh motion solvers
#include "motionSolver.H"
#include "tetDecompositionMotionSolver.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"
#include "setMotionBC.H"
#include "rotatePoints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicMesh.H"
#   include "initContinuityErrs.H"
#   include "initTotalVolume.H"
#   include "createFields.H"

    // Initialize the motion solver
    autoPtr<motionSolver> mPtr = motionSolver::New(mesh);

    // Define the rotation axis and angle from the dictionary
    IOdictionary rotationParams
                 (
                    IOobject
                    (
                        "flippingFoamDict",
                        runTime.findInstance
                        (
                            "",
                            "flippingFoamDict"
                        ),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    )
                 );

    dictionary patchNames(rotationParams.subDict("patchNames"));
    wordList patches = patchNames.toc();

    vector p1(rotationParams.lookup("axisPointStart"));
    vector p2(rotationParams.lookup("axisPointEnd"));
    vector t(rotationParams.lookup("translation"));
    doubleScalar angle = readScalar(rotationParams.lookup("angle"));

    // Convert angle to radians
    angle *= (3.14159/180.0);

    Info<< "\nStarting time loop\n" << endl;
    
    while (runTime.run())
    {
#       include "readPISOControls.H"
#       include "readTimeControls.H"
#       include "checkTotalVolume.H"
#       include "CourantNo.H"       

#       include "setDeltaT.H"       
        
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;        
        
        // Assign boundary conditions to the motion solver
        rotatePoints(mesh, patches, angle, p1, p2, t);
        
        // Solve for mesh-motion
        mesh.movePoints(mPtr->newPoints());
        
#       include "volContinuity.H"        

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);

#       include "UEqn.H"

        // --- PISO loop

        for (int corr=0; corr<nCorr; corr++)
        {
            volScalarField rUA = 1.0/UEqn.A();

            U = rUA*UEqn.H();
            phi = (fvc::interpolate(U) & mesh.Sf());
            //     + fvc::ddtPhiCorr(rUA, U, phi);

            // adjustPhi(phi, U, p);

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

            U -= rUA*fvc::grad(p);
            U.correctBoundaryConditions();
        }  
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;           

        bool meshChanged = mesh.updateTopology(); 
  
        if (meshChanged)
        {
#           include "checkTotalVolume.H"

            // Update the motion solver
            mPtr->updateMesh(mesh.meshMap());

            // Obtain flux from mapped velocity        
            phi = (fvc::interpolate(U) & mesh.Sf());              
#           include "correctPhi.H"           
#           include "CourantNo.H"
        }

        // Write out current parameters
        rotationParams.instance() = runTime.timeName();
        rotationParams.add("patchNames", patchNames, true);
        rotationParams.add("axisPointStart", p1, true);
        rotationParams.add("axisPointEnd", p2, true);
        rotationParams.add("translation", t, true);
        rotationParams.add("angle", angle*(180/3.14159), true);
        
        runTime.write();

        if (runTime.outputTime())
        {
            // Write out mesh quality
            volScalarField meshQuality
            (
                IOobject
                (
                    "meshQuality",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("scalar", dimless, 0.0),
                "zeroGradient"
            );

            meshQuality.internalField() = mesh.meshQuality(true);
            meshQuality.write();

            // Write out the mesh length scales
            volScalarField lengthScale
            (
                IOobject
                (
                    "lengthScale",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("scalar", dimLength, 0.0),
                "zeroGradient"
            );

            lengthScale.internalField() = mesh.lengthScale();
            lengthScale.write();

            // Write out divergence-free fluxes
            volScalarField divPhi = fvc::div(phi);
            divPhi.write();
        }
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
