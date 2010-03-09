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
    freeSurface

Description
    Lagrangian interface tracking class

Author
    Zeljko Tukovic / Sandeep Menon

SourceFiles
    freeSurface.C

\*---------------------------------------------------------------------------*/

#ifndef freeSurface_H
#define freeSurface_H

#include "fvCFD.H"
#include "faCFD.H"
#include "patchToPatchInterpolation.H"
#include "dynamicFvMesh.H"
#include "surfactantProperties.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class freeSurface Declaration
\*---------------------------------------------------------------------------*/

class freeSurface
:
    public IOdictionary,
    public MeshObject<fvMesh, freeSurface>
{
    // Private data

        //- Reference to the underlying mesh
        dynamicFvMesh& mesh_;

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

        //- Time index for the interface
        label curTimeIndex_;

        //- Switch for one/two fluids
        Switch twoFluids_;

        //- Free-surface points displacement direction
        //  parallel to free-surface point normals
        Switch normalMotionDir_;

        //- Devoid of surfactants?
        Switch cleanInterface_;

        //- Patch ID of the interface (side A)
        label aPatchID_;

        //- Patch ID of the interface (side B)
        label bPatchID_;

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

        //- Reference to gravitational acceleration
        dimensionedVector g_;

        //- Surface tension
        dimensionedScalar cleanInterfaceSurfTension_;

        //- Free surface patches which do not move
        wordList fixedFreeSurfacePatches_;

        //- Free surface patches which require
        //  point normal correction
        wordList pointNormalsCorrectionPatches_;

    // Demand driven data

        //- Patch to patch interpolation (A to B)
        mutable patchToPatchInterpolation* interpolatorABPtr_;

        //- Patch to patch interpolation (B to A)
        mutable patchToPatchInterpolation* interpolatorBAPtr_;

        //- Points which are attached to the interface, defining its shape
        mutable vectorIOField* controlPointsPtr_;

        //- Field which determines the motion of interface points
        mutable scalarField* motionPointsMaskPtr_;

        //- Displacement direction of interface points
        mutable vectorField* pointsDisplacementDirPtr_;

        //- Displacement direction of interface control points
        mutable vectorField* facesDisplacementDirPtr_;

        //- Face centres of interface faces
        mutable vectorField* areaCentresPtr_;

        //- Total displacement of free-surface points
        //  in one time step
        mutable vectorIOField* totalDisplacementPtr_;

        //- Finite area mesh
        mutable faMesh* aMeshPtr_;

        //- Free-surface velocity field
        mutable areaVectorField* UsPtr_;

        //- Free-surface fluid flux
        mutable edgeScalarField* phisPtr_;

        //- Free-surface surfactant concetration
        mutable areaScalarField* surfactConcPtr_;

        //- Surface tension field
        mutable areaScalarField* surfaceTensionPtr_;

        //- Surfactant properties
        mutable surfactantProperties* surfactantPtr_;

        //- Fluid indicator
        mutable volScalarField* fluidIndicatorPtr_;

        //- Density
        mutable volScalarField* rhoPtr_;

        //- Viscosity
        mutable volScalarField* muPtr_;

    // Private Member Functions

        // Make demand-driven data

            void makeInterpolators() const;
            void makeControlPoints() const;
            void makeMotionPointsMask() const;
            void makeDirections() const;
            void makeTotalDisplacement() const;
            void readTotalDisplacement() const;
            void makeFaMesh() const;
            void makeUs() const;
            void readUs() const;
            void makePhis() const;
            void makeSurfactConc() const;
            void makeSurfaceTension() const;
            void makeSurfactant() const;
            void makeFluidIndicator() const;
            void makeRho() const;
            void makeMu() const;

        //- Clear all demand-driven data
        void clearOut() const;

        //- Update pressure boundary conditions
        void updatePressure();

        //- Update velocity boundary conditions
        void updateVelocity();

        //- Initialize control points
        void initializeControlPointsPosition() const;

        //- Return reference to interpolator (A to B)
        const patchToPatchInterpolation& interpolatorAB() const;

        //- Return reference to interpolator (B to A)
        const patchToPatchInterpolation& interpolatorBA() const;

        //- Return reference to the control points
        vectorField& controlPoints() const;

        //- Return reference to the point-mask
        scalarField& motionPointsMask() const;

        //- Return reference to point displacement direction field
        vectorField& pointsDisplacementDir() const;

        //- Return reference to control points displacement direction field
        vectorField& facesDisplacementDir() const;

        //- Return reference to the area centres field
        vectorField& areaCentrePositions() const;

        //- Return total points displacement
        vectorField& totalDisplacement() const;

        // Return constant reference to the area mesh
        faMesh& aMesh() const;

        //- Return reference to the free-surface velocity field
        areaVectorField& Us() const;

        //- Return free-surface fluid flux field
        edgeScalarField& Phis() const;

        //- Return free-surface surfactant concentration field
        areaScalarField& surfactantConcentration() const;

        //- Return the surface tension field
        areaScalarField& surfaceTension() const;

        //- Return surface tension gradient
        tmp<areaVectorField> surfaceTensionGrad();

        //- Correct boundary conditions for Us
        void correctUsBoundaryConditions() const;

        //- Return surfactant properties
        const surfactantProperties& surfactant() const;

        //- Return the face curvature field
        const scalarField& faceCurvatures();

        //- Return an indicator field for two fluids
        const volScalarField& fluidIndicator() const;

        //- Move correctedFvPatchField fvSubMeshes
        void moveFvSubMeshes();

        //- Move control points by deltaH and calculate interface
        //  points displacement for the new control points position
        tmp<vectorField> pointDisplacement
        (
            const scalarField& deltaH
        ) const;

        //- Disallow default bitwise copy construct
        freeSurface(const freeSurface&);

        //- Disallow default bitwise assignment
        void operator=(const freeSurface&);


public:

    // Declare name of the class
    ClassName("freeSurface");

    // Constructors

        freeSurface
        (
            dynamicFvMesh& m,
            volVectorField& U,
            volScalarField& p,
            surfaceScalarField& phi
        );


    // Destructor

        ~freeSurface();


    // Member Functions

        // Return constant reference to the mesh
        dynamicFvMesh& mesh() const;

        //- Return reference to velocity field
        volVectorField& U() const
        {
            return U_;
        }

        //- Return reference to pressure field
        volScalarField& p() const
        {
            return p_;
        }

        //- Return reference to flux field
        const surfaceScalarField& phi() const
        {
            return phi_;
        }

        // Return true if two fluids are present
        bool twoFluids() const
        {
            return twoFluids_;
        }

        //- Motion direction switch
        bool normalMotionDir() const
        {
            return normalMotionDir_;
        }

        //- Clean interface switch
        bool cleanInterface() const
        {
            return cleanInterface_;
        }

        // Return the interface patchID
        label aPatchID() const;

        // Return the shadow patchID
        label bPatchID() const;

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

        //- Gravity
        const dimensionedVector& g() const
        {
            return g_;
        }

        //- Surface tension for clean interface
        const dimensionedScalar& cleanInterfaceSurfTension() const
        {
            return cleanInterfaceSurfTension_;
        }

        //- Return the combined density field for two fluids
        volScalarField& rho() const;

        //- Return the combined viscosity field for two fluids
        volScalarField& mu() const;

        //- Adjust the surface-tension for temperature
        void adjustSurfaceTension(const volScalarField& T);

        //- Maximal surface tension based Courant number
        scalar maxCourantNumber();

        // Update the interface with the fluid velocity
        void moveSurfacePoints();

        // Update control points end displacement directions
        void updateDisplacementDirections() const;

        // Restore interface position
        bool restorePosition();

        // Update boundary conditions on velocity and pressure
        void updateBoundaryConditions();

        // Update for mesh motion
        bool movePoints() const;

        // Update on topology change
        bool updateMesh(const mapPolyMesh&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
