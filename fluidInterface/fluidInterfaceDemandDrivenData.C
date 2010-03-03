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
    Demand-driven data for fluidInterface

Author
    Sandeep Menon

\*----------------------------------------------------------------------------*/

#include "fluidInterface.H"

#include "wedgeFaPatch.H"
#include "wedgeFaPatchFields.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void fluidInterface::makeInterpolators() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if
    (
        interpolatorBAPtr_
     || interpolatorABPtr_
    )
    {
        FatalErrorIn("fluidInterface::makeInterpolators()")
            << "patch to patch interpolators already exist."
            << abort(FatalError);
    }

    if (aPatchID() == -1)
    {
        FatalErrorIn("fluidInterface::makeInterpolators()")
            << "Free surface patch A not defined."
            << abort(FatalError);
    }

    if (bPatchID() == -1)
    {
        FatalErrorIn("fluidInterface::makeInterpolators()")
            << "Free surface patch B not defined."
            << abort(FatalError);
    }

    interpolatorBAPtr_ = new patchToPatchInterpolation
    (
        mesh().boundaryMesh()[bPatchID()],
        mesh().boundaryMesh()[aPatchID()],
        intersection::VISIBLE
    );

    const scalarField& faceDistBA =
        interpolatorBAPtr_->faceDistanceToIntersection();

    forAll(faceDistBA, faceI)
    {
        if (mag(faceDistBA[faceI] - GREAT) < SMALL)
        {
            FatalErrorIn("fluidInterface::makeInterpolators()")
                << "Error in B-to-A face patchToPatchInterpolation."
                << abort(FatalError);
        }
    }

    const scalarField& pointDistBA =
        interpolatorBAPtr_->pointDistanceToIntersection();

    forAll(pointDistBA, pointI)
    {
        if (mag(pointDistBA[pointI] - GREAT) < SMALL)
        {
            FatalErrorIn("fluidInterface::makeInterpolators()")
                << "Error in B-to-A point patchToPatchInterpolation."
                << abort(FatalError);
        }
    }

    interpolatorABPtr_ = new patchToPatchInterpolation
    (
        mesh().boundaryMesh()[aPatchID()],
        mesh().boundaryMesh()[bPatchID()],
        intersection::VISIBLE
    );

    const scalarField& faceDistAB =
        interpolatorABPtr_->faceDistanceToIntersection();

    forAll(faceDistAB, faceI)
    {
        if (mag(faceDistAB[faceI] - GREAT) < SMALL)
        {
            FatalErrorIn("fluidInterface::makeInterpolators()")
                << "Error in A-to-B face patchToPatchInterpolation."
                << abort(FatalError);
        }
    }

    const scalarField& pointDistAB =
        interpolatorABPtr_->pointDistanceToIntersection();

    forAll(pointDistAB, pointI)
    {
        if (mag(pointDistAB[pointI] - GREAT)<SMALL)
        {
            FatalErrorIn("fluidInterface::makeInterpolators()")
                << "Error in A-to-B point patchToPatchInterpolation."
                << abort(FatalError);
        }
    }

    Info << "Checking A-to-B and B-to-A interpolators..." << endl;

    scalar maxDist = max
    (
        mag
        (
            interpolatorABPtr_->faceInterpolate
            (
                vectorField
                (
                    mesh().boundaryMesh()[aPatchID()].faceCentres()
                )
            )
          - mesh().boundaryMesh()[bPatchID()].faceCentres()
        )
    );

    scalar maxDistPt = max
    (
        mag
        (
            interpolatorABPtr_->pointInterpolate
            (
                vectorField
                (
                    mesh().boundaryMesh()[aPatchID()].localPoints()
                )
            )
          - mesh().boundaryMesh()[bPatchID()].localPoints()
        )
    );

    Info << "A-to-B interpolation error, face: " << maxDist
         << ", point: " << maxDistPt << endl;

    maxDist = max
    (
        mag
        (
            interpolatorBAPtr_->faceInterpolate
            (
                vectorField
                (
                    mesh().boundaryMesh()[bPatchID()].faceCentres()
                )
            )
          - mesh().boundaryMesh()[aPatchID()].faceCentres()
        )
    );

    maxDistPt = max
    (
        mag
        (
            interpolatorBAPtr_->pointInterpolate
            (
                vectorField
                (
                    mesh().boundaryMesh()[bPatchID()].localPoints()
                )
            )
          - mesh().boundaryMesh()[aPatchID()].localPoints()
        )
    );

    Info << "B-to-A interpolation error, face: " << maxDist
         << ", point: " << maxDistPt << endl;
}

void fluidInterface::makeControlPoints() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (controlPointsPtr_)
    {
        FatalErrorIn("fluidInterface::makeControlPoints()")
            << "Control points already exist"
            << abort(FatalError);
    }

    IOobject controlPointsHeader
    (
        "controlPoints",
        mesh().time().timeName(),
        mesh(),
        IOobject::MUST_READ
    );

    if (controlPointsHeader.headerOk())
    {
        controlPointsPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "controlPoints",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
    }
    else
    {
        controlPointsPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "controlPoints",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                aMesh().areaCentres().internalField()
            );

        initializeControlPointsPosition();
    }
}

void fluidInterface::makeMotionPointsMask() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (motionPointsMaskPtr_)
    {
        FatalErrorIn("fluidInterface::motionPointsMask()")
            << "Motion points mask already exists"
            << abort(FatalError);
    }

    if(aPatchID() == -1)
    {
        FatalErrorIn("fluidInterface::makeMotionPointsMask()")
            << "Interface patch is not defined."
            << abort(FatalError);
    }

    motionPointsMaskPtr_ = new scalarField
    (
        mesh().boundaryMesh()[aPatchID()].nPoints(),
        1.0
    );

    // Mark free surface boundary points
    // which do not belong to empty and wedge patches
    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            (
                aMesh().boundary()[patchI].type()
             != "empty"
            )
         && (
                aMesh().boundary()[patchI].type()
             != "wedge"
            )
        )
        {
            labelList patchPoints =
                aMesh().boundary()[patchI].pointLabels();

            forAll(patchPoints, pointI)
            {
                motionPointsMask()[patchPoints[pointI]] = 0.0;
            }
        }
    }

    // Mark free-surface boundary point
    // at the axis of 2-D axisymmetic cases
    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == "wedge"
        )
        {
            const wedgeFaPatch& wedgePatch =
                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

            if (wedgePatch.axisPoint() > -1)
            {
                motionPointsMask()[wedgePatch.axisPoint()] = 0;
            }
        }
    }
}

void fluidInterface::makeDirections() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if
    (
        pointsDisplacementDirPtr_ ||
        facesDisplacementDirPtr_
    )
    {
        FatalErrorIn("fluidInterface::makeDirections()")
            << "Point and face displacement directions "
            << "already exists"
            << abort(FatalError);
    }

    if(aPatchID() == -1)
    {
        FatalErrorIn("fluidInterface::makeDirections()")
            << "Interface patch is not defined."
            << abort(FatalError);
    }

    pointsDisplacementDirPtr_ =
        new vectorField
        (
            mesh().boundaryMesh()[aPatchID()].nPoints(),
            vector::zero
        );

    facesDisplacementDirPtr_ =
        new vectorField
        (
            mesh().boundaryMesh()[aPatchID()].size(),
            vector::zero
        );

    areaCentresPtr_ =
        new vectorField
        (
            mesh().boundaryMesh()[aPatchID()].size(),
            vector::zero
        );

    if (!normalMotionDir())
    {
        if (mag(g_).value() < SMALL)
        {
            FatalErrorIn("fluidInterface::makeDirections()")
                << "Zero gravity"
                << abort(FatalError);
        }

        facesDisplacementDir() = -(g_/mag(g_)).value();
        pointsDisplacementDir() = -(g_/mag(g_)).value();
    }

    updateDisplacementDirections();
}

void fluidInterface::makeTotalDisplacement() const
{
    if (debug)
    {
        Info << "fluidInterface::makeTotalDisplacement() : "
             << "Making zero total points displacement"
             << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (totalDisplacementPtr_)
    {
        FatalErrorIn("fluidInterface::makeTotalDisplacement()")
            << "Total points displacement already exists"
            << abort(FatalError);
    }

    totalDisplacementPtr_ =
        new vectorIOField
        (
            IOobject
            (
                "totalDisplacement",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            vectorField
            (
                mesh().boundaryMesh()[aPatchID()].nPoints(),
                vector::zero
            )
        );
}

void fluidInterface::readTotalDisplacement() const
{
    if (debug)
    {
        Info<< "fluidInterface::readTotalDisplacement() : "
            << "Reading total points displacement if present"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (totalDisplacementPtr_)
    {
        FatalErrorIn("fluidInterface::readTotalDisplacement()")
            << "Total points displacement already exists"
            << abort(FatalError);
    }

    if
    (
        IOobject
        (
            "totalDisplacement",
            mesh().time().timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ).headerOk()
    )
    {
        totalDisplacementPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "totalDisplacement",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
    }
}

// Make the finite-area mesh
void fluidInterface::makeFaMesh() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (aMeshPtr_)
    {
        FatalErrorIn("fluidInterface::makeFaMesh()")
            << "Area mesh already exists."
            << abort(FatalError);
    }

    // Creating faMesh
    aMeshPtr_ = new faMesh
    (
        mesh(),
        "faMeshDefinition"
    );

    // Set point normal correction patches
    boolList& correction = aMeshPtr_->correctPatchPointNormals();

    forAll(pointNormalsCorrectionPatches_, patchI)
    {
        word patchName = pointNormalsCorrectionPatches_[patchI];

        label patchID = aMeshPtr_->boundary().findPatchID(patchName);

        if (patchID == -1)
        {
            FatalErrorIn
            (
                "fluidInterface::makeFaMesh()"
            )   << "Patch name for point normals correction don't exist"
                << abort(FatalError);
        }

        correction[patchID] = true;
    }
}

void fluidInterface::makeUs() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (UsPtr_)
    {
        FatalErrorIn("fluidInterface::makeUs()")
            << "Free-surface velocity field already exists"
            << abort(FatalError);
    }

    wordList patchFieldTypes
    (
        aMesh().boundary().size(),
        calculatedFaPatchVectorField::typeName
    );

    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == "wedge"
        )
        {
            patchFieldTypes[patchI] = "wedge";
        }
    }

    UsPtr_ = new areaVectorField
    (
        IOobject
        (
            "Us",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh(),
        dimensioned<vector>("Us", dimVelocity, vector::zero),
        patchFieldTypes
    );
}

void fluidInterface::readUs() const
{
    if (debug)
    {
        Info<< "fluidInterface::readUs() : "
            << "Reading free-surface velocity field if present"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (UsPtr_)
    {
        FatalErrorIn("fluidInterface::readUs()")
            << "Free-surface velocity field already exists"
            << abort(FatalError);
    }

    if
    (
        IOobject
        (
            "Us",
            mesh().time().timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ).headerOk()
    )
    {
        UsPtr_ =
        (
            new areaVectorField
            (
                IOobject
                (
                    "Us",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                aMesh()
            )
        );

        correctUsBoundaryConditions();
    }
}

void fluidInterface::makePhis() const
{
    if (debug)
    {
        Info<< "fluidInterface::makePhis() : "
            << "Making free-surface fluid flux"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (phisPtr_)
    {
        FatalErrorIn("fluidInterface::makePhis()")
            << "Free-surface fluid flux already exists"
            << abort(FatalError);
    }

    phisPtr_ = new edgeScalarField
    (
        IOobject
        (
            "phis",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearEdgeInterpolate(Us()) & aMesh().Le()
    );
}

void fluidInterface::makeSurfactConc() const
{
    if (debug)
    {
        Info<< "fluidInterface::makeSurfactConc() : "
            << "Making free-surface surfactant concentration field"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (surfactConcPtr_)
    {
        FatalErrorIn("fluidInterface::makeSurfactConc()")
            << "Free-surface surfactant concentration field already exists"
            << abort(FatalError);
    }


    surfactConcPtr_ = new areaScalarField
    (
        IOobject
        (
            "Cs",
            mesh().time().timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh()
    );
}

void fluidInterface::makeSurfaceTension() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (surfaceTensionPtr_)
    {
        FatalErrorIn("fluidInterface::makeSurfaceTension()")
            << "Surface tension field already exists."
            << abort(FatalError);
    }

    surfaceTensionPtr_ = new areaScalarField
    (
        IOobject
        (
            "surfaceTension",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cleanInterfaceSurfTension()
      + surfactant().surfactR()*
        surfactant().surfactT()*
        surfactant().surfactSaturatedConc()*
        log(1.0 - surfactantConcentration()/
        surfactant().surfactSaturatedConc())
    );
}

void fluidInterface::makeSurfactant() const
{
    if (debug)
    {
        Info<< "fluidInterface::makeSurfactant() : "
            << "Making surfactant properties"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (surfactantPtr_)
    {
        FatalErrorIn("fluidInterface::makeSurfactant()")
            << "Surfactant properties already exists"
            << abort(FatalError);
    }


    const dictionary& surfactProp = subDict("surfactantProperties");

    surfactantPtr_ = new surfactantProperties(surfactProp);
}

void fluidInterface::makeFluidIndicator() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (fluidIndicatorPtr_)
    {
        FatalErrorIn("fluidInterface::makeFluidIndicator()")
            << "fluid indicator already exists"
            << abort(FatalError);
    }

    fluidIndicatorPtr_ = new volScalarField
    (
        IOobject
        (
            "fluidIndicator",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("1", dimless, 1.0),
        "zeroGradient"
    );

    volScalarField& fluidIndicator = *fluidIndicatorPtr_;

    if (twoFluids_)
    {
        // Find start cell
        label pointOnShadowPatch = mesh().boundaryMesh()[bPatchID()][0][0];

        label startCell = mesh().pointCells()[pointOnShadowPatch][0];

        // Get cell-cells addressing
        const labelListList& cellCells = mesh().cellCells();

        SLList<label> slList(startCell);

        while (slList.size())
        {
            label curCell = slList.removeHead();

            if (fluidIndicator[curCell] == 1)
            {
                fluidIndicator[curCell] = 0.0;

                for (int i = 0; i < cellCells[curCell].size(); i++)
                {
                    slList.append(cellCells[curCell][i]);
                }
            }
        }
    }

    fluidIndicator.correctBoundaryConditions();
}

void fluidInterface::makeRho() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (rhoPtr_)
    {
        FatalErrorIn("fluidInterface::makeRho()")
            << "rho already exists"
            << abort(FatalError);
    }

    rhoPtr_ = new volScalarField
    (
        IOobject
        (
            "rho",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("0", dimMass/dimVolume, 0)
    );

    volScalarField& rho = *rhoPtr_;

    rho = (fluidIndicator()*(rhoFluidA() - rhoFluidB()) + rhoFluidB());

    rho.correctBoundaryConditions();
}

void fluidInterface::makeMu() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (muPtr_)
    {
        FatalErrorIn("fluidInterface::makeMu()")
            << "mu already exists"
            << abort(FatalError);
    }

    muPtr_ = new volScalarField
    (
        IOobject
        (
            "mu",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("0", dimPressure*dimTime, 0)
    );

    volScalarField& mu = *muPtr_;

    mu = (fluidIndicator()*(muFluidA() - muFluidB()) + muFluidB());

    mu.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return reference to interpolator (A to B)
const patchToPatchInterpolation& fluidInterface::interpolatorAB() const
{
    if (!interpolatorABPtr_)
    {
        makeInterpolators();
    }

    return *interpolatorABPtr_;
}

//- Return reference to interpolator (B to A)
const patchToPatchInterpolation& fluidInterface::interpolatorBA() const
{
    if (!interpolatorBAPtr_)
    {
        makeInterpolators();
    }

    return *interpolatorBAPtr_;
}

//- Return reference to the control points
vectorField& fluidInterface::controlPoints() const
{
    if (!controlPointsPtr_)
    {
        makeControlPoints();
    }

    return *controlPointsPtr_;
}

scalarField& fluidInterface::motionPointsMask() const
{
    if (!motionPointsMaskPtr_)
    {
        makeMotionPointsMask();
    }

    return *motionPointsMaskPtr_;
}

// Return reference to point displacement direction field
vectorField& fluidInterface::pointsDisplacementDir() const
{
    if (!pointsDisplacementDirPtr_)
    {
        makeDirections();
    }

    return *pointsDisplacementDirPtr_;
}

vectorField& fluidInterface::areaCentrePositions() const
{
    if (!areaCentresPtr_)
    {
        makeDirections();
    }

    return *areaCentresPtr_;
}

// Return reference to control points displacement direction field
vectorField& fluidInterface::facesDisplacementDir() const
{
    if (!facesDisplacementDirPtr_)
    {
        makeDirections();
    }

    return *facesDisplacementDirPtr_;
}

vectorField& fluidInterface::totalDisplacement() const
{
    if (!totalDisplacementPtr_)
    {
        makeTotalDisplacement();
    }

    return *totalDisplacementPtr_;
}

faMesh& fluidInterface::aMesh() const
{
    if (!aMeshPtr_)
    {
        makeFaMesh();
    }

    return *aMeshPtr_;
}

areaVectorField& fluidInterface::Us() const
{
    if (!UsPtr_)
    {
        makeUs();
    }

    return *UsPtr_;
}

edgeScalarField& fluidInterface::Phis() const
{
    if (!phisPtr_)
    {
        makePhis();
    }

    return *phisPtr_;
}

const volScalarField& fluidInterface::fluidIndicator() const
{
    if (!fluidIndicatorPtr_)
    {
        makeFluidIndicator();
    }

    return *fluidIndicatorPtr_;
}

areaScalarField& fluidInterface::surfactantConcentration() const
{
    if (!surfactConcPtr_)
    {
        makeSurfactConc();
    }

    return *surfactConcPtr_;
}

areaScalarField& fluidInterface::surfaceTension() const
{
    if (!surfaceTensionPtr_)
    {
        makeSurfaceTension();
    }

    return *surfaceTensionPtr_;
}

const surfactantProperties& fluidInterface::surfactant() const
{
    if (!surfactantPtr_)
    {
        makeSurfactant();
    }

    return *surfactantPtr_;
}

tmp<areaVectorField> fluidInterface::surfaceTensionGrad()
{
    tmp<areaVectorField> tgrad
    (
        new areaVectorField
        (
            IOobject
            (
                "surfaceTensionGrad",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (-fac::grad(surfactantConcentration())*
            surfactant().surfactR()*surfactant().surfactT()/
            (1.0 - surfactantConcentration()/
            surfactant().surfactSaturatedConc()))()
        )
    );

    return tgrad;
}

volScalarField& fluidInterface::rho() const
{
    if (!rhoPtr_)
    {
        makeRho();
    }

    return *rhoPtr_;
}

volScalarField& fluidInterface::mu() const
{
    if (!muPtr_)
    {
        makeMu();
    }

    return *muPtr_;
}

} // End namespace Foam

// ************************************************************************* //
