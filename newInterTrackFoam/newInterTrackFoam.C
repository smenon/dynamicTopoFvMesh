/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    interfaceTrackinFoam

Description
    Incompressible laminar CFD code for simulation of a single bubble rising 
    in a stil liquid. Interface between fluid phases is tracked using moving
    mesh.

\*---------------------------------------------------------------------------*/

#include "OFstream.H"
#include "fvCFD.H"
#include "inletOutletFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "motionSolver.H"
#include "freeSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"

#   include "createMeshNoClear.H"

#   include "createFields.H"

#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info << "Time = " << runTime.value() << endl << endl;

#       include "readPISOControls.H"
#       include "readFreeSurfaceControls.H"

        interface.updateDisplacementDirections();

        interface.moveMeshPointsForOldFreeSurfDisplacement();

        if (smooth)
        {
            interface.smooth();
        }

        interface.movePoints();

        Info<< "\nMax surface Courant Number = " 
            << interface.maxCourantNumber() << endl << endl;

        for (int corr=0; corr<nOuterCorr; corr++)
        {
            // Update interface bc
            interface.updateBoundaryConditions();

            // Make the flux relative
            phi -= fvc::meshPhi(rho, U);

#           include "CourantNo.H"

            fvVectorMatrix UEqn
            (
                fvm::ddt(rho, U)
              + fvm::div(fvc::interpolate(rho)*phi, U)
              - fvm::laplacian(mu, U)
            );

            solve(UEqn == -fvc::grad(p));

            volScalarField rUA = 1.0/UEqn.A();

            // --- PISO loop
            for (int i=0; i<nCorr; i++)
            {
                U = UEqn.H()*rUA;
                
                phi = (fvc::interpolate(U) & mesh.Sf());

                // Non-orthogonal pressure corrector loop
                for (label nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                {
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rUA, p) 
                     == fvc::div(phi)
                    );

                    pEqn.setReference(pRefCell, pRefValue);
                    pEqn.solve();

                    if (nonOrth == nNonOrthCorr)
                    {
                        phi -= pEqn.flux();
                    }
                }

#               include "continuityErrs.H"

                // Momentum corrector
                U -= rUA*fvc::grad(p);
                U.correctBoundaryConditions();
            }

            // Move mesh
            interface.movePoints();

#           include "freeSurfaceContinuityErrs.H"
        }

        runTime.write();

        Info << "ExecutionTime = "
            << scalar(runTime.elapsedCpuTime())
            << " s\n" << endl << endl;
    }

    Info << "End\n" << endl;

    return(0);
}

// ************************************************************************* //
