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
#include "dynamicFvMesh.H"
#include "freeSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "initContinuityErrs.H"
#   include "initTotalVolume.H"
#   include "createFields.H"

    freeSurface interface(mesh, U, p, phi);

    Info<< "Reading field rAU if present\n" << endl;
    volScalarField rAU
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        runTime.deltaT()/interface.rho()
    );

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readPISOControls.H"
#       include "readTimeControls.H"
#       include "checkTotalVolume.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        interface.updateDisplacementDirections();

        bool meshChanged = interface.restorePosition();

        if (meshChanged)
        {
#           include "correctPhi.H"
        }

#       include "volContinuity.H"

        interface.moveSurfacePoints();

        Info << "\nMax surface Courant Number = "
             << interface.maxCourantNumber() << endl << endl;

        for (int corr=0; corr<nOuterCorr; corr++)
        {
            // Update boundary conditions on velocity and pressure
            interface.updateBoundaryConditions();

#           include "CourantNo.H"

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, interface.rho(), U);

            fvVectorMatrix UEqn
            (
                fvm::ddt(interface.rho(), U)
              + fvm::div(fvc::interpolate(interface.rho())*phi, U)
              - fvm::laplacian(interface.mu(), U)
            );

            solve(UEqn == -fvc::grad(p));

            rAU = 1.0/UEqn.A();

            // --- PISO loop

            for (int corr=0; corr<nCorr; corr++)
            {
                U = rAU*UEqn.H();

                phi = (fvc::interpolate(U) & mesh.Sf());

                for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                {
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rAU, p) == fvc::div(phi)
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

                U -= rAU*fvc::grad(p);
                U.correctBoundaryConditions();
            }

            interface.moveSurfacePoints();

#           include "freeSurfaceContinuityErrs.H"

            Info << endl;
        }

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;

        runTime.write();
    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
