/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright held by original author
    \\/      M anipulation   |
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
    engineMotion

Description
    Driver routine to test engine valve and piston motion.

Author
    Sandeep Menon

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicTopoFvMesh.H"
#include "interpolationTable.H"
#include "motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicMesh.H"

    // Initialize the motion solver
    autoPtr<motionSolver> mPtr = motionSolver::New(mesh);

    // Define the rotation axis and angle from the dictionary
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

    // Read in the piston axis
    vector pistonAxis = pistonDict.lookup("axis");

    pistonAxis /= mag(pistonAxis) + VSMALL;

    scalar oldStroke = pistonTable(runTime.value()), currentStroke = 0.0;

    // Obtain the number of valves in the system
    label numValves = readLabel(engineDict.lookup("numValves"));

    PtrList<interpolationTable<scalar> > valveLiftTables(numValves);
    List<scalar> oldLift(numValves, 0.0), currentLift(numValves, 0.0);
    List<vector> valveAxes(numValves, vector::zero);

    // Obtain the valve patch names
    dictionary valves(engineDict.subDict("valves"));
    wordList valveList = valves.toc();

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

        // Read in the valve axis
        valveAxes[valveI] = valveDict.lookup("axis");

        // Normalize the axis
        valveAxes[valveI] /= mag(valveAxes[valveI]) + VSMALL;

        // Interpolate and obtain the lift value
        // for the current time-step
        oldLift[valveI] = valveLiftTables[valveI](runTime.value());
    }

    const polyBoundaryMesh& boundary = mesh.boundaryMesh();

    for (runTime++; !runTime.end(); runTime++)
    {
        Info << "Time = " << runTime.value() << endl << endl;

        // Copy existing point locations
        pointField meshPoints(mesh.points());

        // Update the piston
        currentStroke = pistonTable(runTime.value());

        const labelList& pistonPoints =
        (
            boundary[boundary.findPatchID(pistonName[0])].meshPoints()
        );

        vector p = (currentStroke - oldStroke)*pistonAxis;

        Info << "Piston " << pistonName[0]
             << ": Current Stroke value: " << currentStroke << endl;
        Info << "Piston " << pistonName[0]
             << ": Old Stroke value: " << oldStroke << endl;

        forAll(pistonPoints, index)
        {
            meshPoints[pistonPoints[index]] += p;
        }

        oldStroke = currentStroke;

        vector v = vector::zero;

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
            v = (currentLift[valveI] - oldLift[valveI])*valveAxes[valveI];

            Info << "Valve " << valveList[valveI]
                 << ": Current Lift value: " << currentLift[valveI] << endl;
            Info << "Valve " << valveList[valveI]
                 << ": Old Lift value: " << oldLift[valveI] << endl;

            forAll(valvePoints, index)
            {
                meshPoints[valvePoints[index]] += v;
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

        // Update mesh motion
        mesh.movePoints(mPtr->newPoints());

        // Update mesh for topology changes
        bool meshChanged = mesh.updateTopology();

        if (meshChanged)
        {
            // Update the motion solver
            mPtr->updateMesh(mesh.meshMap());
        }

        runTime.write();

#       include "meshInfo.H"
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
