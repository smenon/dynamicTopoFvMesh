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
#include "patchToPatchInterpolation.H"

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

        //- Reference to the velocity field
        volVectorField& U_;

        //- Reference to the pressure field
        volScalarField& p_;

        //- Reference to the flux field
        surfaceScalarField& phi_;

        //- Name of the interface patch
        word interfacePatch_;

        //- Name of the interface shadow patch
        word shadowPatch_;

        //- Switch for one/two fluids
        Switch twoFluids_;

        //- Patch ID of the interface (side A)
        label aPatchID_;

        //- Patch ID of the interface (side B)
        label bPatchID_;

        //- Time index for the interface
        label curTimeIndex_;

        //- Viscosity of fluid A
        dimensionedScalar muFluidA_;

        //- Density of fluid A
        dimensionedScalar rhoFluidA_;

        //- Conductivity of fluid A
        dimensionedScalar condFluidA_;

        //- Specific Heat of fluid A
        dimensionedScalar CpFluidA_;

        //- Viscosity of fluid B
        dimensionedScalar muFluidB_;

        //- Density of fluid B
        dimensionedScalar rhoFluidB_;

        //- Conductivity of fluid B
        dimensionedScalar condFluidB_;

        //- Specific Heat of fluid B
        dimensionedScalar CpFluidB_;

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

        //- Surface tension field
        mutable areaScalarField* surfaceTensionPtr_;

        //- Finite area mesh
        mutable faMesh* aMeshPtr_;

        //- Free-surface velocity field
        mutable areaVectorField* UsPtr_;

        //- Free surface patches which do not move
        wordList fixedFreeSurfacePatches_;

        //- Patch to patch interpolation (A to B)
        mutable patchToPatchInterpolation* interpolatorABPtr_;

        //- Patch to patch interpolation (B to A)
        mutable patchToPatchInterpolation* interpolatorBAPtr_;

        //- Fluid indicator
        mutable volScalarField* fluidIndicatorPtr_;

    // Private Member Functions

        // Make demand-driven data

            void makeFluidIndicator();
            void makeInterpolators();
            void makeControlPoints();
            void makeMotionPointsMask();
            void makeDisplacement();
            void makeDirections();
            void makeFaMesh();
            void makeSurfaceTension();
            void makeUs();

        //- Clear all demand-driven data
        void clearOut();

        //- Update pressure boundary conditions
        void updatePressure();

        //- Update velocity boundary conditions
        void updateVelocity();

        //- Compute and return the gradient of surface-tension
        tmp<areaVectorField> surfaceTensionGrad();

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

        //- Return reference to the free-surface velocity field
        areaVectorField& Us();

        //- Correct boundary conditions for Us
        void correctUsBoundaryConditions();

        //- Move correctedFvPatchField fvSubMeshes
        void moveFvSubMeshes();

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

        // Return true if two fluids are present
        bool twoFluids() const
        {
            return twoFluids_;
        }

        // Return the interface patchID
        label aPatchID();

        // Return the shadow patchID
        label bPatchID();

        //- Viscosity of fluid A
        const dimensionedScalar& muFluidA() const
        {
            return muFluidA_;
        }

        //- Density of fluid A
        const dimensionedScalar& rhoFluidA() const
        {
            return rhoFluidA_;
        }

        //- Conductivity of fluid A
        const dimensionedScalar& condFluidA() const
        {
            return condFluidA_;
        }

        //- Specific Heat of fluid A
        const dimensionedScalar& CpFluidA() const
        {
            return CpFluidA_;
        }

        //- Viscosity of fluid B
        const dimensionedScalar& muFluidB() const
        {
            return muFluidB_;
        }

        //- Density of fluid B
        const dimensionedScalar& rhoFluidB() const
        {
            return rhoFluidB_;
        }

        //- Conductivity of fluid B
        const dimensionedScalar& condFluidB() const
        {
            return condFluidB_;
        }

        //- Specific Heat of fluid B
        const dimensionedScalar& CpFluidB() const
        {
            return CpFluidB_;
        }

        //- Return the surface tension field
        areaScalarField& surfaceTension();

        //- Return the face curvature field
        const scalarField& faceCurvatures();

        //- Return reference to interpolator (A to B)
        const patchToPatchInterpolation& interpolatorAB();

        //- Return reference to interpolator (B to A)
        const patchToPatchInterpolation& interpolatorBA();

        //- Return an indicator field for two fluids
        const volScalarField& fluidIndicator();

        //- Adjust the surface-tension for temperature
        void adjustSurfaceTension(const volScalarField& T);

        // Return reference to the displacement field
        vectorField& displacement();

        // Update the interface with fluid velocity
        void updateInterface();

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
