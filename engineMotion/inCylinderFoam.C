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
#include "dynamicTopoFvMesh.H"
#include "interpolationTable.H"
#include "motionSolver.H"

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

    // Define the engine parameters
    IOdictionary engineDict
    (
        IOobject
        (
            "engineDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Obtain the piston patch name
    dictionary piston(engineDict.subDict("piston"));
    wordList pistonName = piston.toc();

    // Read in necessary information from the dictionary
    dictionary pistonDict = engineDict.subDict(pistonName[0]);

    // Read the profile
    interpolationTable<scalar> pistonTable
    (
        fileName(pistonDict.lookup("profileFile"))
    );

    // Set outOfBounds handling to clamp
    pistonTable.outOfBounds(interpolationTable<scalar>::CLAMP);

    // Read in the piston axis
    vector pistonAxis = pistonDict.lookup("axis");

    pistonAxis /= mag(pistonAxis) + VSMALL;

    scalar oldStroke = pistonTable(runTime.value()), currentStroke = 0.0;

    // Obtain the number of valves in the system
    dictionary valves(engineDict.subDict("valves"));
    wordList valveList = valves.toc();

    label numValves = valveList.size();

    PtrList<interpolationTable<scalar> > valveLiftTables(numValves);
    List<scalar> oldLift(numValves, 0.0), currentLift(numValves, 0.0);
    List<vector> valveAxes(numValves, vector::zero);

    // Read in all necessary information from the dictionary
    forAll(valveList, valveI)
    {
        dictionary valveDict = engineDict.subDict(valveList[valveI]);

        // Read in the lift profile table
        valveLiftTables.set
        (
            valveI,
            new interpolationTable<scalar>
            (
                fileName(valveDict.lookup("profileFile"))
            )
        );

        // Set outOfBounds handling to clamp
        valveLiftTables[valveI].outOfBounds
        (
            interpolationTable<scalar>::CLAMP
        );

        // Read in the valve axis
        valveAxes[valveI] = valveDict.lookup("axis");

        // Normalize the axis
        valveAxes[valveI] /= mag(valveAxes[valveI]) + VSMALL;

        // Interpolate and obtain the lift value
        // for the current time-step
        oldLift[valveI] = valveLiftTables[valveI](runTime.value());
    }

    const polyBoundaryMesh& boundary = mesh.boundaryMesh();

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

        // Copy existing point locations
        pointField meshPoints(mesh.points());

        // Update the piston
        currentStroke = pistonTable(runTime.value());

        const labelList& pistonPoints =
        (
            boundary[boundary.findPatchID(pistonName[0])].meshPoints()
        );

        vector pD = (currentStroke - oldStroke)*pistonAxis;

        Info << "Piston " << pistonName[0]
             << ": Current Stroke value: " << currentStroke << endl;
        Info << "Piston " << pistonName[0]
             << ": Old Stroke value: " << oldStroke << endl;

        forAll(pistonPoints, index)
        {
            meshPoints[pistonPoints[index]] += pD;
        }

        oldStroke = currentStroke;

        vector vD = vector::zero;

        // Update valves
        forAll(valveList, valveI)
        {
            // Obtain the current lift value
            currentLift[valveI] = valveLiftTables[valveI](runTime.value());

            // Move the valve patch by the difference in lift
            const labelList& valvePoints =
            (
                boundary[boundary.findPatchID(valveList[valveI])].meshPoints()
            );

            // Define displacement
            vD = (currentLift[valveI] - oldLift[valveI])*valveAxes[valveI];

            Info << "Valve " << valveList[valveI]
                 << ": Current Lift value: " << currentLift[valveI] << endl;
            Info << "Valve " << valveList[valveI]
                 << ": Old Lift value: " << oldLift[valveI] << endl;

            dictionary valveDict = engineDict.subDict(valveList[valveI]);

            forAll(valvePoints, index)
            {
                meshPoints[valvePoints[index]] += vD;
            }

            oldLift[valveI] = currentLift[valveI];
        }

        // Obtain the field from the registry
        pointField& refPoints = const_cast<pointField&>
        (
            mesh.lookupObject<pointField>("refPoints")
        );

        // Assign boundary conditions to the motion solver
        refPoints = meshPoints;

        // Solve for mesh-motion
        mesh.movePoints(mPtr->newPoints());

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

        runTime.write();
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
