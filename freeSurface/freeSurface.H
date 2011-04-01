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
    freeSurfacePointDisplacement.C
    makeFreeSurfaceData.C

\*---------------------------------------------------------------------------*/

#ifndef freeSurface_H
#define freeSurface_H

#include "fvCFD.H"
#include "faCFD.H"
#include "IOpatchToPatchInterpolation.H"
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

        //- Reference to fluid flux field
        const surfaceScalarField& phi_;

        //- Name of the interface patch
        word interfacePatch_;

        //- Name of the interface shadow patch
        word shadowPatch_;

        //- Are patches ordered for interpolation
        Switch orderedPatches_;

        //- Time index for the interface
        mutable label curTimeIndex_;

        //- Switch for one/two fluids
        Switch twoFluids_;

        //- Free-surface points displacement direction
        //  parallel to free-surface point normals
        Switch normalMotionDir_;

        //- Free-surface points displacement direction
        //  if not normal motion direction
        vector motionDir_;

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

        //- Number of free-surface correctors
        label nFreeSurfCorr_;

    // Demand-driven data

        //- Patch to patch interpolation object which deals with
        //  interpolation of flow variables between two free-surface
        //  (interface) patches
        mutable IOpatchToPatchInterpolation* interpolatorABPtr_;

        //- Patch to patch interpolation object which deals with
        //  interpolation of flow variables between two free-surface
        //  (interface) patches
        mutable IOpatchToPatchInterpolation* interpolatorBAPtr_;

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

        //- Free-surface flux
        mutable areaScalarField* phiPtr_;

        //- Free-surface flux
        mutable areaScalarField* ddtPhiPtr_;

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
            void makePhis() const;
            void makePhi() const;
            void makeDdtPhi() const;
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

        //- Initialize control points position
        void initializeControlPointsPosition() const;

        //- Return reference to interpolator from free surface side A to
        //  free surface side B
        const IOpatchToPatchInterpolation& interpolatorAB() const;

        //- Return reference to interpolator from free surface side B to
        //  free surface side A
        const IOpatchToPatchInterpolation& interpolatorBA() const;

        // Interpolate from A to B using coupleMap
        template<class Type>
        tmp<Field<Type> > pointInterpolateAB(const Field<Type>& field);

        // Interpolate from B to A using coupleMap
        template<class Type>
        tmp<Field<Type> > pointInterpolateBA(const Field<Type>& field);

        // Interpolate from A to B using coupleMap
        template<class Type>
        tmp<Field<Type> > faceInterpolateAB(const Field<Type>& field);

        // Interpolate from B to A using coupleMap
        template<class Type>
        tmp<Field<Type> > faceInterpolateBA(const Field<Type>& field);

        //- Return control points
        vectorField& controlPoints() const;

        //- Return reference to motion points mask field
        scalarField& motionPointsMask() const;

        //- Return reference to point displacement direction field
        vectorField& pointsDisplacementDir() const;

        //- Return reference to control points displacement direction field
        vectorField& facesDisplacementDir() const;

        //- Return total points displacement
        vectorField& totalDisplacement() const;

        //- Return reference to finite area mesh
        faMesh& aMesh();

        //- Return const reference to finite area mesh
        const faMesh& aMesh() const;

        //- Return reference to free-surface velocity field
        areaVectorField& Us();

        //- Return const reference to free-surface velocity field
        const areaVectorField& Us() const;

        //- Return reference to surface tension field
        areaScalarField& surfaceTension();

        //- Return const reference to surface tension field
        const areaScalarField& surfaceTension() const;

        //- Return surface tension gradient
        tmp<areaVectorField> surfaceTensionGrad();

        //- Return surfactant properties
        const surfactantProperties& surfactant() const;

        //- Return reference to surfactant concentration field
        areaScalarField& surfactantConcentration();

        //- Return const reference to surfactant concentration field
        const areaScalarField& surfactantConcentration() const;

        //- Return free-surface fluid flux field
        edgeScalarField& Phis();

        //- Return reference to free-surface fluid flux field
        areaScalarField& phi();

        //- Return const reference to free-surface fluid flux field
        const areaScalarField& phi() const;

        //- Return free-surface fluid flux field
        areaScalarField& ddtPhi();

        //- Move correctedFvPatchField fvSubMeshes
        void moveFvSubMeshes() const;

        //- Move control points by deltaH and calculate interface
        //  points displacement for the new control points position
        tmp<vectorField> pointDisplacement(const scalarField& deltaH) const;

        //- Update free-surface flux
        void updateSurfaceFlux();

        //- Update free-surface surfactant concentration
        void updateSurfactantConcentration();

        //- Correct surface velocity boundary conditions
        void correctUsBoundaryConditions();

        // Correct points for cyclics
        void correctCyclics(pointField& points) const;

        //- Disallow default bitwise copy construct
        freeSurface(const freeSurface&);

        //- Disallow default bitwise assignment
        void operator=(const freeSurface&);

public:

    // Declare name of the class and it's debug switch
    ClassName("freeSurface");


    // Constructors

    freeSurface
    (
        dynamicFvMesh& m,
        volVectorField& Ub,
        volScalarField& Pb,
        const surfaceScalarField& phi
    );


    // Destructor

    ~freeSurface();


    // Member Functions


        //- Return reference to fvMesh
        dynamicFvMesh& mesh()
        {
            return mesh_;
        }

        //- Return reference to fvMesh
        const dynamicFvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return reference to DB
        const Time& DB() const
        {
            return U_.mesh().time();
        }

        //- Return reference to velocity field
        volVectorField& U()
        {
            return U_;
        }

        //- Return reference to velocity field
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return reference to pressure field
        volScalarField& p()
        {
            return p_;
        }

        //- Return reference to pressure field
        const volScalarField& p() const
        {
            return p_;
        }

        //- Are there exist two fluid
        //  in the free-surface model
        bool twoFluids() const
        {
            return twoFluids_;
        }

        //- Motion direction swithc
        bool normalMotionDir() const
        {
            return normalMotionDir_;
        }

        //- Is it interface clean
        bool cleanInterface() const
        {
            return cleanInterface_;
        }

        //- ID of patch which represent A side
        //  of free-surface
        const label& aPatchID() const
        {
            return aPatchID_;
        }

        //- ID of patch which represent B side
        //  of free-surface
        const label& bPatchID() const
        {
            return bPatchID_;
        }

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

        //- Densiti of fluid B
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

        //- Return referenct to fluid indicator
        const volScalarField& fluidIndicator() const;

        //- Return the combined density field for two fluids
        const volScalarField& rho() const;

        //- Return the combined viscosity field for two fluids
        const volScalarField& mu() const;

        //- Update control points end displacement directions
        void updateDisplacementDirections() const;

        //- Predict interface points position
        bool predictPoints();

        //- Correct interface points position
        bool correctPoints();

        //- Move only free-surface points
        bool movePoints(const scalarField& interfacePhi);

        //- Move mesh points for old time free-surface points displacement
        bool moveMeshPointsForOldFreeSurfDisplacement();

        //- Move mesh points
        bool moveMeshPoints();

        //- Update free-surface boundary conditions
        void updateBoundaryConditions();

        //- Calculate free-surface undulation indicator
        tmp<scalarField> divPhiSUndulationIndicator();
        tmp<scalarField> gradGradDivPhiS();

        //- Calculate total pressure force
        vector totalPressureForce() const;

        //- Calculate total viscous force
        vector totalViscousForce() const;

        //- Calculate total surface tension force
        vector totalSurfaceTensionForce() const;

        //- Calculate free-surface undulation indicator
        tmp<scalarField> undulationIndicator();

        //- Chech free-surface faces flatness
        tmp<scalarField> checkFaceFlatness();

        //- Maximal surface tension based Courant number
        scalar maxCourantNumber();

        //- Update properties
        void updateProperties();

        // Update for mesh motion
        bool movePoints() const;

        // Update on topology change
        bool updateMesh(const mapPolyMesh&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
