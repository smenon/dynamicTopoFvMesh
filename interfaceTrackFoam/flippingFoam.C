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

// Rotate all points belonging to the 'patchNames' by 'theta'
// about an axis p1-p2, translate by a distance defined by 't'
void rotatePoints
(
    dynamicTopoFvMesh& mesh,
    wordList& patchNames,
    doubleScalar theta,
    vector p1,
    vector p2,
    vector t
)
{
    label nP = 0;
    vector p, q;
    doubleScalar costheta = Foam::cos(theta);
    doubleScalar sintheta = Foam::sin(theta);

    labelList patchID(patchNames.size(), -1);
    labelHashSet pointSet;

    // Define the rotation axis and normalize it
    vector r = (p2-p1)/mag(p2-p1);

    // Fetch the mesh-points
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
    const pointField& oldMeshPoints(mesh.points());

    // Copy existing point locations
    pointField meshPoints(oldMeshPoints);

    // Match patch names and add to a HashSet to avoid moving
    // patch points twice
    forAll(patchNames, wordI)
    {
        forAll(mesh.boundaryMesh(), patchI)
        {
            if (bMesh[patchI].name() == patchNames[wordI])
            {
                patchID[nP++] = patchI;

                // Add all points of this patch
                const labelList& patchPoints = bMesh[patchI].meshPoints();

                forAll(patchPoints, index)
                {
                    if (!pointSet.found(patchPoints[index]))
                    {
                        pointSet.insert(patchPoints[index]);
                    }
                }

                break;
            }
        }
    }

    // Move all patch points cumulatively
    labelList allPatchPoints = pointSet.toc();

    forAll(allPatchPoints, index)
    {
        q = vector::zero;

        // Fetch the old point and translate it
        p = oldMeshPoints[allPatchPoints[index]] + t;

        // Translate to the origin
        p -= p1;

        // Apply the rotation matrix
        q.x() += (costheta + (1 - costheta) * r.x() * r.x()) * p.x();
        q.x() += ((1 - costheta) * r.x() * r.y() - r.z() * sintheta) * p.y();
        q.x() += ((1 - costheta) * r.x() * r.z() + r.y() * sintheta) * p.z();

        q.y() += ((1 - costheta) * r.x() * r.y() + r.z() * sintheta) * p.x();
        q.y() += (costheta + (1 - costheta) * r.y() * r.y()) * p.y();
        q.y() += ((1 - costheta) * r.y() * r.z() - r.x() * sintheta) * p.z();

        q.z() += ((1 - costheta) * r.x() * r.z() - r.y() * sintheta) * p.x();
        q.z() += ((1 - costheta) * r.y() * r.z() + r.x() * sintheta) * p.y();
        q.z() += (costheta + (1 - costheta) * r.z() * r.z()) * p.z();

        // Translate back to original location
        q += p1;

        // Assign to the mesh
        meshPoints[allPatchPoints[index]] = q;
    }

    // Update the displacement BCs for mesh motion
    setMotionBC(mesh, meshPoints);
}

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
                        runTime.constant(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
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

    // Enable/disable debugging
    mesh.debug = true;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info << "Time = " << runTime.value() << endl << endl;

        p1 += t; p2 += t;

        // Update boundary points and solve for mesh-motion
        rotatePoints(mesh, patches, angle, p1, p2, t);

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

            meshQuality.internalField() = mesh.meshQuality();
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
