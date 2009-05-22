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

// Mesh motion solvers
#include "motionSolver.H"
#include "tetDecompositionMotionSolver.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"
#include "setMotionBC.H"

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

    fluidInterface interface(mesh, U, p, phi);

    // Obtain fluid indicator from the interface
    volScalarField* fluidIndicatorPtr = NULL;

    fluidIndicatorPtr = new volScalarField
    (
        IOobject
        (
            "fluidIndicator",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        interface.fluidIndicator()
    );

    volScalarField& fluidIndicator = *fluidIndicatorPtr;

    volScalarField nu
    (
        IOobject
        (
            "nu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluidIndicator*
        (
            (interface.muFluidA()/interface.rhoFluidA())
          - (interface.muFluidB()/interface.rhoFluidB())
        )
      + (interface.muFluidB()/interface.rhoFluidB())
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    // Initialize the motion solver
    autoPtr<motionSolver> mPtr = motionSolver::New(mesh);

    bool nonNewtonian = false;
    bool solveForTemperature = false;
    bool adjustNuForTemperature = false;
    bool adjustSigmaForTemperature = false;

    // Read in flags from the dictionary
    if (interface.found("nonNewtonian"))
    {
        nonNewtonian = Switch(interface.lookup("nonNewtonian"));
    }

    if (interface.found("solveForTemperature"))
    {
        solveForTemperature =
            Switch(interface.lookup("solveForTemperature"));
    }

    // Maintain the temperature field as a pointer
    autoPtr<volScalarField> TPtr(NULL);
    if (solveForTemperature)
    {
        Info<< "Reading field T\n" << endl << flush;

        TPtr.set
        (
            new volScalarField
            (
                IOobject
                (
                    "T",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );

        if (interface.found("adjustNuForTemperature"))
        {
            adjustNuForTemperature =
                Switch(interface.lookup("adjustNuForTemperature"));
        }

        if (interface.found("adjustSigmaForTemperature"))
        {
            adjustSigmaForTemperature =
                Switch(interface.lookup("adjustSigmaForTemperature"));
        }
    }

    // Introduce the non-Newtonian transport model
    autoPtr<viscosityModel> nuModel(NULL);
    if (nonNewtonian)
    {
        nuModel = viscosityModel::New("nu",interface,U,phi);
    }

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

        if (solveForTemperature)
        {
            if (adjustSigmaForTemperature)
            {
                Info << "Adjusting surface tension for temperature." << endl;
                interface.adjustSurfaceTension(TPtr());
            }

            if (adjustNuForTemperature)
            {
                Info << "Adjusting viscosity for temperature." << endl;
                adjustViscosity
                (
                    interface,
                    TPtr(),
                    (interface.muFluidA()/interface.rhoFluidA()).value(),
                    nu
                );
            }
        }

        // Set boundary conditions for the motionSolver and solve for mesh-motion
        interface.restorePosition();
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

                // Make the fluxes relative
                fvc::makeRelative(phi, U);

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

        if (solveForTemperature)
        {
            dimensionedScalar DT
            (
                "DT",
                interface.condFluidA()/
                (interface.CpFluidA()*interface.rhoFluidA())
            );

            // Passive heat-transfer
            solve
            (
                fvm::ddt(TPtr())
              + fvm::div(phi, TPtr())
              - fvm::laplacian(DT, TPtr())
            );
        }

        // Make the fluxes absolute
        fvc::makeAbsolute(phi, U);

        runTime.write();
#       include "meshInfo.H"

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
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
