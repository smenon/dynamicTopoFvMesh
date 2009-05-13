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
    // Set mutex debug option before mesh is created.
    // Mutex::debug = true;

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

    for (runTime++; !runTime.end(); runTime++)
    {
        Info << "Time = " << runTime.value() << endl << endl;

        // Translate the axis
        p1 += t; p2 += t;

        // Update boundary points and solve for mesh-motion
        rotatePoints(mesh, patches, angle, p1, p2, t);

        // Update mesh motion
        mesh.movePoints(mPtr->newPoints());

        // Obtain mesh stats before topo-changes
        mesh.meshQuality(true);

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
        }
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
