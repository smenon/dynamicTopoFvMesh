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
#include "fluidInterface.H"

// Mesh motion solvers
#include "motionSolver.H"
#include "tetDecompositionMotionSolver.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"
#include "setMotionBC.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicMesh.H"
#   include "initContinuityErrs.H"
#   include "initTotalVolume.H"
#   include "createFields.H"

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimMass/dimVolume, 0)
    );

    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimPressure*dimTime, 0)
    );

    fluidInterface interface(mesh, rho, U, p, phi);

    mu =
    (
        interface.fluidIndicator()*
        (
            interface.muFluidA()
          - interface.muFluidB()
        )
      + interface.muFluidB()
    );

    rho =
    (
        interface.fluidIndicator()*
        (
              interface.rhoFluidA()
            - interface.rhoFluidB()
        )
      + interface.rhoFluidB()
    );

    Info<< "Reading field rUA if present\n" << endl;
    volScalarField rUA
    (
        IOobject
        (
            "rUA",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        runTime.deltaT()/rho
    );

    if (interface.twoFluids())
    {
        // Set initial velocity only for the liquid phase
        U = interface.fluidIndicator()*U;
    }

    Info << "\nStarting time loop\n" << endl;

    // Initialize the motion solver
    autoPtr<motionSolver> mPtr = motionSolver::New(mesh);

    while (runTime.run())
    {
#       include "readPISOControls.H"
#       include "readTimeControls.H"
#       include "checkTotalVolume.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Update the interface
        interface.updateInterface();

        setMotionBC(mesh, interface.aPatchID(), interface.displacement());

        if (interface.twoFluids())
        {
            // Interpolate displacement to the shadow patch
            pointField dispB = interface.interpolatorAB().pointInterpolate
                               (
                                   interface.displacement()
                               );

            setMotionBC(mesh, interface.bPatchID(), dispB);
        }

        // Solve for motion
        mesh.movePoints(mPtr->newPoints());

#       include "volContinuity.H"

        for (int corr=0; corr<nOuterCorr; corr++)
        {
            // Update boundary conditions on velocity and pressure
            interface.updateBoundaryConditions();

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, rho, U);

            fvVectorMatrix UEqn
            (
                fvm::ddt(rho, U)
              + fvm::div(fvc::interpolate(rho)*phi, U)
              - fvm::laplacian(mu, U)
            );

            solve(UEqn == -fvc::grad(p));

            rUA = 1.0/UEqn.A();

            // --- PISO loop

            for (int corr=0; corr<nCorr; corr++)
            {
                U = rUA*UEqn.H();

                phi = (fvc::interpolate(U) & mesh.Sf());

                adjustPhi(phi, U, p);

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

                // Make the fluxes relative
                fvc::makeRelative(phi, rho, U);

                U -= rUA*fvc::grad(p);
                U.correctBoundaryConditions();
            }

#           include "freeSurfaceContinuityErrs.H"

            Info << endl;
        }

        // Make the fluxes absolute
        fvc::makeAbsolute(phi, rho, U);

        bool meshChanged = mesh.updateTopology();

        if (meshChanged)
        {
#           include "checkTotalVolume.H"
            phi = (fvc::interpolate(U) & mesh.Sf());
#           include "correctPhi.H"
#           include "CourantNo.H"

            // Update the interface
            interface.updateMesh(mesh.meshMap());

            // Update the motion solver
            mPtr->updateMesh(mesh.meshMap());

            // Update viscosity with the new fluid indicator
            mu =
            (
                interface.fluidIndicator()*
                (
                    interface.muFluidA()
                  - interface.muFluidB()
                )
              + interface.muFluidB()
            );

            // Update density with the new fluid indicator
            rho =
            (
                interface.fluidIndicator()*
                (
                      interface.rhoFluidA()
                    - interface.rhoFluidB()
                )
              + interface.rhoFluidB()
            );
        }

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;

        runTime.write();

#       include "meshInfo.H"
    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
