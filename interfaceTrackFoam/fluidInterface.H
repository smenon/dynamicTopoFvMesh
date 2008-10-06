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

Class
    fluidInterface

Description
    Lagrangian interface tracking class

Author
    Zeljko Tukovic / Sandeep Menon

SourceFiles
    fluidInterface.C

\*---------------------------------------------------------------------------*/

#ifndef fluidInterface_H
#define fluidInterface_H

#include "fvCFD.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fluidInterface Declaration
\*---------------------------------------------------------------------------*/

class fluidInterface
:
    public IOdictionary
{
    // Private data
        
        //- Reference to the underlying mesh
        fvMesh& mesh_;
        
        //- Reference to the density field
        volScalarField& rho_;

        //- Reference to the velocity field
        volVectorField& U_;

        //- Reference to the pressure field
        volScalarField& p_;

        //- Reference to the flux field
        surfaceScalarField& phi_;    
        
        //- Name of the interface patch
        word interfacePatch_;
        
        //- Patch ID of the interface
        label patchID_;
        
        //- Time index for the interface
        label curTimeIndex_;
        
        //- Viscosity of fluid A
        dimensionedScalar muA_;
        
        //- Density of fluid A
        dimensionedScalar rhoA_;
        
        //- Surface tension
        dimensionedScalar surfaceTension_;
        
        //- Interface displacement field
        mutable vectorField* displacementPtr_;
        
        //- Field which determines the motion of interface points
        mutable scalarField* motionPointsMaskPtr_;        
        
        //- Points which are attached to the interface, defining its shape
        mutable vectorIOField* controlPointsPtr_;        
        
        //- Displacement direction of interface points
        mutable vectorField* pointsDisplacementDirPtr_;

        //- Displacement direction of interface control points
        mutable vectorField* facesDisplacementDirPtr_;        
        
        //- Finite area mesh which consists
        mutable faMesh* aMeshPtr_;    
        
        //- Free surface patches which do not move
        wordList fixedFreeSurfacePatches_;        
        
    // Private Member Functions
        
        // Make demand-driven data
        
            void makeControlPoints();
            void makeMotionPointsMask();
            void makeDisplacement();
            void makeDirections();            
            void makeFaMesh();
            
        //- Clear all demand-driven data
        void clearOut();   
        
        //- Update pressure boundary conditions
        void updatePressure();
        
        //- Update velocity boundary conditions
        void updateVelocity();
        
        //- Initialize control points 
        void initializeControlPointsPosition();
        
        //- Return reference to the control points
        vectorField& controlPoints();        
        
        //- Return reference to point displacement direction field
        vectorField& pointsDisplacementDir();

        //- Return reference to control points displacement direction field
        vectorField& facesDisplacementDir();   
        
        //- Return reference to the point-mask
        scalarField& motionPointsMask();
        
        //- Move control points by deltaH and calculate interface 
        //  points displacement for the new control points position
        tmp<vectorField> pointDisplacement(const scalarField& deltaH);        

        //- Disallow default bitwise copy construct
        fluidInterface(const fluidInterface&);

        //- Disallow default bitwise assignment
        void operator=(const fluidInterface&);


public:

    // Declare name of the class
    ClassName("fluidInterface");

    // Constructors

        fluidInterface
        (
            fvMesh& m,
            volScalarField& rho,
            volVectorField& U,
            volScalarField& p,
            surfaceScalarField& phi
        );


    // Destructor

        ~fluidInterface();


    // Member Functions

        // Return constant reference to the mesh
        fvMesh& mesh();
        
        // Return constant reference to the area mesh
        faMesh& areaMesh();

        // Return the interface patchID
        const label patchID();

        // Return reference to the displacement field
        vectorField& displacement();

        // Update the interface with the fluid velocity
        void movePoints();
        
        // Update control points end displacement directions
        void updateDisplacementDirections();        
        
        // Restore interface position
        void restorePosition();
        
        // Update boundary conditions on velocity and pressure
        void updateBoundaryConditions();

        // Update on topology change
        void updateMesh(const mapPolyMesh& mpm);
        
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
