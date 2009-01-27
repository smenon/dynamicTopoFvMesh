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

// Included for point-normals post-processing
#include "pointMesh.H"
#include "pointFields.H"
#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Adjust fluid viscosity for temperature variations
// according to the Vogel-Fulcher-Tammann (VFT) model
void adjustViscosity
(
    const dictionary& coeffDict,
    const volScalarField& T,
    const scalar nuRef,
    volScalarField& nu
)
{
    // Read coeffiecients for the VFT model
    dictionary vftCoeffs(coeffDict.subDict("vftCoeffs"));

    scalar Tref = readScalar(vftCoeffs.lookup("Tref"));
    scalar Tv   = readScalar(vftCoeffs.lookup("Tv"));
    scalar Ta   = readScalar(vftCoeffs.lookup("Ta"));

    // Compute constants for the model
    scalar a = (Ta / Tref);
    scalar b = Tv;
    scalar c = (Tref / (Tref - Tv));

    // Correct the internalField
    forAll(nu.internalField(), cellI)
    {
        scalar Tcell = T.internalField()[cellI];

        nu.internalField()[cellI] =
            nuRef*::exp(a*((Tref / (Tcell - b)) - c));
    }
}

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicMesh.H"
#   include "initContinuityErrs.H"
#   include "initTotalVolume.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    bool nonNewtonian = false;
    if (interface.found("nonNewtonian"))
    {
        nonNewtonian = Switch(interface.lookup("nonNewtonian"));
    }

    bool adjustNuForTemperature = false;
    if (interface.found("adjustNuForTemperature"))
    {
        adjustNuForTemperature =
            Switch(interface.lookup("adjustNuForTemperature"));
    }

    bool adjustSigmaForTemperature = false;
    if (interface.found("adjustSigmaForTemperature"))
    {
        adjustSigmaForTemperature =
            Switch(interface.lookup("adjustSigmaForTemperature"));
    }

    // Introduce the non-Newtonian transport model
    autoPtr<viscosityModel> nuModel(NULL);
    if (nonNewtonian)
    {
        nuModel = viscosityModel::New("nu",interface,U,phi);
    }

    //polyMesh::debug = true;
    //mesh.debug = true;

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

        if (adjustSigmaForTemperature)
        {
            interface.adjustSurfaceTension(T);
        }

        if (adjustNuForTemperature)
        {
            adjustViscosity
            (
                interface,
                T,
                (interface.muFluidA()/interface.rhoFluidA()).value(),
                nu
            );
        }

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

        if (nonNewtonian)
        {
            nuModel->correct();
            nu = nuModel->nu();
        }

        // Make the fluxes relative
        fvc::makeRelative(phi, U);

        dimensionedScalar DT
        (
            "DT",
            interface.condFluidA()/
            (interface.CpFluidA()*interface.rhoFluidA())
        );

        // Passive heat-transfer
        solve
        (
            fvm::ddt(T)
          + fvm::div(phi, T)
          - fvm::laplacian(DT, T)
        );

        // Make the fluxes absolute
        fvc::makeAbsolute(phi, U);

        runTime.write();
#       include "meshInfo.H"

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

    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
