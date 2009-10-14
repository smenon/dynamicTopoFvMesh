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
    flippingFoam

Description
    Driver routine to test mesh-motion and topology changes.

Author
    Sandeep Menon

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicTopoFvMesh.H"

// Mesh motion solvers
#include "motionSolver.H"
#include "tetDecompositionMotionSolver.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"
#include "setMotionBC.H"
#include "rotatePoints.H"

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
    IOdictionary rotationParams
                 (
                    IOobject
                    (
                        "flippingFoamDict",
                        mesh.time().constant(),
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
    bool solveForMotion = readBool(rotationParams.lookup("solveForMotion"));

    // Convert angle to radians
    angle *= (3.14159/180.0);

    for (runTime++; !runTime.end(); runTime++)
    {
        Info << "Time = " << runTime.value() << endl << endl;

        // Translate the axis
        p1 += t; p2 += t;

        // Update boundary points and solve for mesh-motion
        rotatePoints(mesh, patches, angle, p1, p2, t, solveForMotion);

        // Update mesh motion
        if (solveForMotion)
        {
            mesh.movePoints(mPtr->newPoints());
        }

        // Update mesh for topology changes
        bool meshChanged = mesh.updateTopology();

        if (meshChanged)
        {
            // Update the motion solver
            mPtr->updateMesh(mesh.meshMap());
        }

        // Write out current parameters
        rotationParams.instance() = runTime.timeName();
        rotationParams.add("patchNames", patchNames, true);
        rotationParams.add("axisPointStart", p1, true);
        rotationParams.add("axisPointEnd", p2, true);
        rotationParams.add("translation", t, true);
        rotationParams.add("angle", angle*(180/3.14159), true);
        rotationParams.add("solveForMotion", solveForMotion, true);

        runTime.write();

#       include "meshInfo.H"
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
