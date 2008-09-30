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

// Rotate all points belonging to the "box" patch by 'theta' 
// about an axis p1-p2, translate by a distance defined by 't', and update the mesh
void rotatePoints(dynamicTopoFvMesh& mesh, doubleScalar theta, vector p1, vector p2, vector t)
{
    label patchID = -1;
    vector p_orig, p, q;
    doubleScalar costheta = Foam::cos(theta), sintheta = Foam::sin(theta);
    
    // Define the rotation axis and normalize it
    vector r = (p2-p1)/mag(p2-p1);
        
    // Fetch the mesh-points
    pointField meshPoints = mesh.points();
    
    // Find the box patch
    forAll (mesh.boundary(), patchI)
    {
        if(mesh.boundary()[patchI].name() == "box") {
            patchID = patchI;
            break;
        }
    }   
    
    // Move points that lie on the box patch
    const labelList& boxPoints = mesh.boundaryMesh()[patchID].meshPoints();  
    vectorField Displacement(boxPoints.size(),vector::zero);
    
    forAll (boxPoints, index)
    {
        q = vector::zero;
        // Fetch the point
        p_orig = p = meshPoints[boxPoints[index]] + t;
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
        meshPoints[boxPoints[index]] = q;
        // Change displacement conditions as well
        Displacement[index] = q-p_orig;
    }   
    
    // Update the displacement
    mesh.setMotionBC(patchID, Displacement);
    
    // Move the boundary points
    mesh.movePoints(meshPoints);      
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // For debug purposes, assign <root> and <case>
    argc = 3; 
    argv[1] = new char[100];
    argv[2] = new char[100];
    strcpy(argv[1],"/home/smenon/OpenFOAM/smenon-1.4.1-dev/run");
    strcpy(argv[2],"holeinsqr");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicMesh.H"

    // Define the rotation axis and angle (in radians)
    //vector p1(1.0,0.5,0.0), p2(1.0,0.5,1.0), t(0.0,0.0,0.0); // holeinslab
    //vector p1(0.75,1.0,0.0), p2(0.75,1.0,1.0), t(0.005,0.0,0.0); // Street
    vector p1(0.5,0.5,0.0), p2(0.5,0.5,1.0), t(0.0,0.0,0.0); // OffsetSqr
    doubleScalar angle = (0.72*3.14159)/180.0;
    
    // Enable/disable debugging
    mesh.debug = true;

    for (runTime++; !runTime.end(); runTime++)
    {    
        Info << "Time = " << runTime.value() << endl << endl;
        
        // Update boundary points and move them
        rotatePoints(mesh, angle, p1, p2, t);
        p1 += t; p2 += t;
        
        // Update mesh (Solve for motion and topology)        
        mesh.updateMotion();        
        mesh.updateTopology();
        
        runTime.write();
    }

    Info << "End\n" << endl;
    
    delete [] argv[1];
    delete [] argv[2];

    return 0;
}


// ************************************************************************* //
