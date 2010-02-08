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
    inCylinderFoam

Description
    Solver for cold-flow in-cylinder engine simulations with dynamic mesh.

Author
    Sandeep Menon

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "initContinuityErrs.H"
#   include "initTotalVolume.H"
#   include "createFields.H"

#   include "readEngineOptions.H"

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

#       include "setEngineOptions.H"

#       include "volContinuity.H"

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);

#       include "UEqn.H"

        // --- PISO loop

        for (int corr=0; corr<nCorr; corr++)
        {
            rUA = 1.0/UEqn.A();

            U = rUA*UEqn.H();
            phi = (fvc::interpolate(U) & mesh.Sf());

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

            // Some boundary conditions require fluxes to be relative
            fvc::makeRelative(phi, U);

            U -= rUA*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        // Make the fluxes absolute before manipulating the mesh.
        fvc::makeAbsolute(phi, U);

        // Update the mesh
        mesh.update();

        runTime.write();
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
