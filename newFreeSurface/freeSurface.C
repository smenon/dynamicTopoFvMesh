/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description

\*---------------------------------------------------------------------------*/

#include "freeSurface.H"
#include "volFields.H"
#include "primitivePatchInterpolation.H"
#include "fixedGradientFvPatchFields.H"
#include "zeroGradientCorrectedFvPatchFields.H"
#include "fixedGradientCorrectedFvPatchFields.H"
#include "fixedValueCorrectedFvPatchFields.H"
#include "emptyFaPatch.H"
#include "wedgeFaPatch.H"
#include "demandDrivenData.H"
#include "EulerDdtScheme.H"
#include "CrankNicholsonDdtScheme.H"
#include "backwardDdtScheme.H"
#include "facLnGrad.H"

#include "wallFvPatch.H"

#include "fixedGradientFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"
#include "transformField.H"
#include "twoDPointCorrector.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(freeSurface, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void freeSurface::clearOut()
{
    deleteDemandDrivenData(interpolatorABPtr_);
    deleteDemandDrivenData(interpolatorBAPtr_);
    deleteDemandDrivenData(controlPointsPtr_);
    deleteDemandDrivenData(motionPointsMaskPtr_);
    deleteDemandDrivenData(pointsDisplacementDirPtr_);
    deleteDemandDrivenData(facesDisplacementDirPtr_);
    deleteDemandDrivenData(totalDisplacementPtr_);
    deleteDemandDrivenData(aMeshPtr_);
    deleteDemandDrivenData(UsPtr_);
    deleteDemandDrivenData(phisPtr_);
    deleteDemandDrivenData(surfactConcPtr_);
    deleteDemandDrivenData(surfaceTensionPtr_);
    deleteDemandDrivenData(surfactantPtr_);
    deleteDemandDrivenData(fluidIndicatorPtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

freeSurface::freeSurface
(
    fvMesh& m,
    const volScalarField& rho,
    volVectorField& Ub, 
    volScalarField& Pb, 
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "freeSurfaceProperties",
            Ub.mesh().time().constant(),
            Ub.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(m),
    rho_(rho),
    U_(Ub),
    p_(Pb),
    phi_(phi),
    curTimeIndex_(Ub.mesh().time().timeIndex()),
    twoFluids_
    (
        this->lookup("twoFluids")
    ),
    normalMotionDir_
    (
        this->lookup("normalMotionDir")
    ),
    cleanInterface_
    (
        this->lookup("cleanInterface")
    ),
    aPatchID_(-1),
    bPatchID_(-1),
    muFluidA_
    (
        this->lookup("muFluidA")
    ),
    muFluidB_
    (
        this->lookup("muFluidB")
    ),
    rhoFluidA_
    (
        this->lookup("rhoFluidA")
    ),
    rhoFluidB_
    (
        this->lookup("rhoFluidB")
    ),
    condFluidA_
    (
        this->lookup("condFluidA")
    ),
    condFluidB_
    (
        this->lookup("condFluidB")
    ),
    CpFluidA_
    (
        this->lookup("CpFluidA")
    ),
    CpFluidB_
    (
        this->lookup("CpFluidB")
    ),            
    g_(this->lookup("g")),
    cleanInterfaceSurfTension_
    (
        this->lookup("surfaceTension")
    ),
    fixedFreeSurfacePatches_
    (
        this->lookup("fixedFreeSurfacePatches")
    ),
    pointNormalsCorrectionPatches_
    (
        this->lookup("pointNormalsCorrectionPatches")
    ),
    interpolatorABPtr_(NULL),
    interpolatorBAPtr_(NULL),
    controlPointsPtr_(NULL),
    motionPointsMaskPtr_(NULL),
    pointsDisplacementDirPtr_(NULL),
    facesDisplacementDirPtr_(NULL),
    totalDisplacementPtr_(NULL),
    aMeshPtr_(NULL),
    UsPtr_(NULL),
    phisPtr_(NULL),
    surfactConcPtr_(NULL),
    surfaceTensionPtr_(NULL),
    surfactantPtr_(NULL),
    fluidIndicatorPtr_(NULL)
{
    // Make finite area mesh
    if (!aMeshPtr_)
    {
        makeFaMesh();
    }

    // Create free-surface velocity field
    if (!UsPtr_)
    {
        makeUs();
    }

    // Create surfactant concentration and fluid flux field
    if(!cleanInterface())
    {
        if (!surfactConcPtr_)
        {
            makeSurfactConc();
        }

        if (!phisPtr_)
        {
            makePhis();
        }
    }

    // Detect the free surface patch
    forAll (mesh().boundary(), patchI)
    {
        if(mesh().boundary()[patchI].name() == "freeSurface")
        {
            aPatchID_ = patchI;
                
            Info<< "Found free surface patch. ID: " << aPatchID_
                << endl;
        }
    }

    if(aPatchID() == -1)
    {
        FatalErrorIn("freeSurface::freeSurface(...)")
            << "Free surface patch not defined.  Please make sure that "
                << " the free surface patches is named as freeSurface"
                << abort(FatalError);
    }

    // Detect the free surface shadow patch
    if (twoFluids())
    {
        forAll (mesh().boundary(), patchI)
        {
            if(mesh().boundary()[patchI].name() == "freeSurfaceShadow")
            {
                bPatchID_ = patchI;
                    
                Info<< "Found free surface shadow patch. ID: "
                    << bPatchID_ << endl;
            }
        }

        if(bPatchID() == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Free surface shadow patch not defined. "
                    << "Please make sure that the free surface shadow patch "
                    << "is named as freeSurfaceShadow."
                    << abort(FatalError);
        }
    }

    // Define free-surface points displacement directions
    if (!pointsDisplacementDirPtr_)
    {
        makeDirections();
    }

    // Define motion points displacement mask
    if (!motionPointsMaskPtr_)
    {
        makeMotionPointsMask();
    }
        
    // Mark free surface boundary points 
    // which do not belonge to empty and wedge patch
    forAll(aMesh().boundary(), patchI)
    {
        if 
        (
            (
                aMesh().boundary()[patchI].type()
             != emptyFaPatch::typeName
            )
         && (
                aMesh().boundary()[patchI].type()
             != wedgeFaPatch::typeName
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
         == wedgeFaPatch::typeName
        )
        {
            const wedgeFaPatch& wedgePatch =
                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

            if(wedgePatch.axisPoint() > -1)
            {
                motionPointsMask()[wedgePatch.axisPoint()] = 0;

                Info << "Axis point: " 
                    << wedgePatch.axisPoint()
                    << "vector: " 
                    << aMesh().points()[wedgePatch.axisPoint()] << endl;
            }
        }
    }

    // Define control points position
    if (!controlPointsPtr_)
    {
        makeControlPoints();
    }

    // Read free-surface points total displacement if present
    if (!totalDisplacementPtr_)
    {
        readTotalDisplacement();
    }

    // Define variable surface tension calculation
    if (!cleanInterface())
    {
        // Read and create surfactant properties

        if (!surfactantPtr_)
        {
            makeSurfactant();
        }
            

        if (!surfaceTensionPtr_)
        {
            makeSurfaceTension();
        }
    }

    // Set point normal correction patches
    boolList& correction = aMesh().correctPatchPointNormals();

    forAll(pointNormalsCorrectionPatches_, patchI)
    {
        word patchName = pointNormalsCorrectionPatches_[patchI];

        label patchID = aMesh().boundary().findPatchID(patchName);

        if(patchID == -1)
        {
            FatalErrorIn
            (
                "freeSurface::freeSurface(...)"
            )   << "patch name for point normals correction don't exist"
                << abort(FatalError);
        }

        correction[patchID] = true;
    }
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

freeSurface::~freeSurface()
{
    clearOut();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void freeSurface::updateDisplacementDirections()
{
    if(normalMotionDir())
    {
        // Update point displacement correction
        pointsDisplacementDir() = aMesh().pointAreaNormals();

        // TEMPORARY FIX:
        // Correct point displacement direction 
        // at the "right" symmetryPlane which represents the axis
        // of an axisymmetric case
        forAll(aMesh().boundary(), patchI)
        {
            label centerLinePatchID = 
                aMesh().boundary().findPatchID("right");

            if(centerLinePatchID != -1)
            {
		const labelList& pointLabels = 
                    aMesh().boundary()[centerLinePatchID].pointLabels();

                forAll(pointLabels, pointI)
                {
		    pointsDisplacementDir()[pointLabels[pointI]] = vector(1,0,0);
		}

	    } else {
                    Info << "Warning: right polyPatch does not exist. " 
                        << "Free surface points displacement directions "
                        << "will not be corrected at the axis (centerline)" 
                        << endl;
	    }	    
	}	

        forAll(aMesh().boundary(), patchI)
        {
            if(aMesh().boundary()[patchI].type() == wedgeFaPatch::typeName)
            {
                const wedgeFaPatch& wedgePatch =
                    refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);
		
                vector axis = wedgePatch.axis();

                label centerLinePatchID = 
                    aMesh().boundary().findPatchID("right");
	
                if(centerLinePatchID != -1)
                {
                    const labelList& pointLabels = 
                        aMesh().boundary()[centerLinePatchID].pointLabels();
                    
                    forAll(pointLabels, pointI)
                    {
                        vector dir = 
                            pointsDisplacementDir()[pointLabels[pointI]];
                    
                        dir = (dir&axis)*axis;
                        dir /= mag(dir);
                
                        pointsDisplacementDir()[pointLabels[pointI]] = dir;
                    }
                }
                else
                {
                    Info << "Warning: centerline polyPatch does not exist. " 
                        << "Free surface points displacement directions "
                        << "will not be corrected at the axis (right)" 
                        << endl; 
                }
            
                break;   
            }
        }

	// Update face displacement direction
        facesDisplacementDir() =
            aMesh().faceAreaNormals().internalField();

        // Correction of control points postion
        const vectorField& Cf = aMesh().centres().internalField();

        controlPoints() = 
            facesDisplacementDir()
           *(facesDisplacementDir()&(controlPoints() - Cf))
          + Cf;
    }    
}


bool freeSurface::movePoints()
{
    pointField newMeshPoints = mesh().points();

    scalarField sweptVolCorr = 
        phi().boundaryField()[aPatchID()]
      - fvc::meshPhi(rho(),U())().boundaryField()[aPatchID()];

    word ddtScheme
    (
        mesh().ddtScheme("ddt(" + rho().name() + ',' + U().name()+')')
    );

    if 
    (
        ddtScheme
     == fv::CrankNicholsonDdtScheme<vector>::typeName
    )
    {
        sweptVolCorr *= (1.0/2.0)*DB().deltaT().value();
    }
    else if
    (
        ddtScheme
     == fv::EulerDdtScheme<vector>::typeName
    )
    {
        sweptVolCorr *= DB().deltaT().value();
    }
    else if
    (
        ddtScheme
     == fv::backwardDdtScheme<vector>::typeName
    )
    {
        sweptVolCorr *= (2.0/3.0)*DB().deltaT().value();
    }   
    else
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Unsupported temporal differencing scheme : "
                << ddtScheme
                << abort(FatalError);
    }


    const scalarField& Sf = aMesh().S();
    const vectorField& Nf = aMesh().faceAreaNormals().internalField();

    scalarField deltaH =
        sweptVolCorr/(Sf*(Nf & facesDisplacementDir()));
        

    pointField displacement = pointDisplacement(deltaH);

    // Move only free-surface points

    const labelList& meshPointsA = 
        mesh().boundaryMesh()[aPatchID()].meshPoints();

    forAll (displacement, pointI)
    {
        newMeshPoints[meshPointsA[pointI]] += displacement[pointI];
    }

    if(twoFluids_)
    {
        const labelList& meshPointsB = 
            mesh().boundaryMesh()[bPatchID_].meshPoints();

        pointField displacementB =
            interpolatorAB().pointInterpolate
            (
                displacement
            );
        
        forAll (displacementB, pointI)
        {
            newMeshPoints[meshPointsB[pointI]] += displacementB[pointI]; 
        }
    }

    // Update total displacement field
    /*
    if(totalDisplacementPtr_ && (curTimeIndex_ < DB().timeIndex()))
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Total displacement of free surface points "
                << "from previous time step is not absorbed by the mesh."
                << abort(FatalError);
    }
    else*/ if (curTimeIndex_ < DB().timeIndex())
    {
        totalDisplacement() = displacement;

        curTimeIndex_ = DB().timeIndex();
    }
    else
    {
        totalDisplacement() += displacement;
    }

    twoDPointCorrector twoDPointCorr(mesh());

    twoDPointCorr.correctPoints(newMeshPoints);
    
    mesh().movePoints(newMeshPoints);

    aMesh().movePoints(mesh().points());

    // Move correctedFvPatchField fvSubMeshes
    moveFvSubMeshes();

    return true;
}

//- Move correctedFvPatchField fvSubMeshes
void freeSurface::moveFvSubMeshes()
{
    // Move correctedFvPatchField fvSubMeshes
    forAll(U().boundaryField(), patchI)
    {
        if
        (
            (
                U().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<vector>::typeName
            )
        )
        {
            correctedFvPatchField<vector>& pU =
                refCast<correctedFvPatchField<vector> >
                (
                    U().boundaryField()[patchI]
                );

            pU.movePatchSubMesh();
        }
    }

    forAll(p().boundaryField(), patchI)
    {
        if
        (
            (
                p().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<scalar>::typeName
            )
        )
        {
            correctedFvPatchField<scalar>& pP =
                refCast<correctedFvPatchField<scalar> >
                (
                    p().boundaryField()[patchI]
                );

            pP.movePatchSubMesh();
        }
    }    
}


bool freeSurface::moveMeshPointsForOldFreeSurfDisplacement()
{
    if(totalDisplacementPtr_)
    {
        
        pointField newPoints = mesh().points();

        const labelList& meshPointsA = 
            mesh().boundaryMesh()[aPatchID()].meshPoints();

        forAll (totalDisplacement(), pointI)
        {
            newPoints[meshPointsA[pointI]] -= totalDisplacement()[pointI]; 
        }

        mesh().movePoints(newPoints);

        aMesh().movePoints(mesh().points());

        // Move correctedFvPatchField fvSubMeshes
        moveFvSubMeshes();        
    }      

    return true;
}

//- Update mesh corresponding to the given map
void freeSurface::updateMesh(const mapPolyMesh& mpm)
{
    // For the moment, wipe-out existing data and re-create.
    
    // Re-make faMesh.
    deleteDemandDrivenData(aMeshPtr_);
    makeFaMesh();
    
    // Re-make the freeSurface velocity field
    deleteDemandDrivenData(UsPtr_);
    makeUs();
    
    // Re-make the motionPointsMask
    deleteDemandDrivenData(motionPointsMaskPtr_);
    makeMotionPointsMask();    

    // Mark free surface boundary points 
    // which do not belong to empty and wedge patches
    forAll(aMesh().boundary(), patchI)
    {
        if 
        (
            (
                aMesh().boundary()[patchI].type()
             != emptyFaPatch::typeName
            )
         && (
                aMesh().boundary()[patchI].type()
             != wedgeFaPatch::typeName
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
    
    // Re-make the displacement directions
    deleteDemandDrivenData(controlPointsPtr_);
    deleteDemandDrivenData(pointsDisplacementDirPtr_);
    deleteDemandDrivenData(facesDisplacementDirPtr_);    
    makeDirections();

    // Mark free-surface boundary point 
    // at the axis of 2-D axisymmetic cases
    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == wedgeFaPatch::typeName
        )
        {
            const wedgeFaPatch& wedgePatch =
                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

            if(wedgePatch.axisPoint() > -1)
            {
                motionPointsMask()[wedgePatch.axisPoint()] = 0;

                Info << "Axis point: " 
                    << wedgePatch.axisPoint()
                    << "vector: " 
                    << aMesh().points()[wedgePatch.axisPoint()] << endl;
            }
        }
    }
  
    // Set point normal correction patches
    boolList& correction = aMesh().correctPatchPointNormals();

    forAll(pointNormalsCorrectionPatches_, patchI)
    {
        word patchName = pointNormalsCorrectionPatches_[patchI];

        label patchID = aMesh().boundary().findPatchID(patchName);

        if(patchID == -1)
        {
            FatalErrorIn
            (
                "freeSurface::freeSurface(...)"
            )   << "patch name for point normals correction don't exist"
                << abort(FatalError);
        }

        correction[patchID] = true;
    }          
}

void freeSurface::updateBoundaryConditions()
{
    updateVelocity();
    updateSurfactantConcentration();
    updatePressure();
}


void freeSurface::updateVelocity()
{
    if(twoFluids())
    {
        vectorField nA = mesh().boundary()[aPatchID()].nf();

        vectorField nB = mesh().boundary()[bPatchID()].nf();

        // Update fixedValue boundary condition on patch B

        scalarField DnB = interpolatorBA().faceInterpolate
        (
            mesh().boundary()[bPatchID()].deltaCoeffs()
        );

        scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();


        vectorField UtPA = 
            U().boundaryField()[aPatchID()].patchInternalField();

        if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedGradientCorrectedFvPatchField<vector>& aU =
                refCast<fixedGradientCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            UtPA += aU.corrVecGrad();
        }


        UtPA -= nA*(nA & UtPA);


        vectorField UtPB = interpolatorBA().faceInterpolate
        (
            U().boundaryField()[bPatchID()].patchInternalField()
        );

        if
        (
            U().boundaryField()[bPatchID()].type()
         == fixedValueCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedValueCorrectedFvPatchField<vector>& bU =
                refCast<fixedValueCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[bPatchID()]
                );

            UtPB += interpolatorBA().faceInterpolate(bU.corrVecGrad());
        }

        UtPB -= nA*(nA & UtPB);


        vectorField UtFs = muFluidA().value()*DnA*UtPA 
            + muFluidB().value()*DnB*UtPB;

        vectorField UnFs = nA*phi().boundaryField()[aPatchID()]
            /mesh().boundary()[aPatchID()].magSf();

        Us().internalField() += UnFs - nA*(nA&Us().internalField());
        correctUsBoundaryConditions();


        UtFs -= (muFluidA().value() - muFluidB().value())*
            (fac::grad(Us())&aMesh().faceAreaNormals())().internalField();

        if(!cleanInterface())
        {
            UtFs += surfaceTensionGrad()().internalField();
        }

        UtFs /= muFluidA().value()*DnA + muFluidB().value()*DnB + VSMALL;

        
        Us().internalField() = UnFs + UtFs;
        correctUsBoundaryConditions();


        U().boundaryField()[bPatchID()] == 
            interpolatorAB().faceInterpolate(UtFs)
          + nB*fvc::meshPhi(rho(),U())().boundaryField()[bPatchID()]/
            mesh().boundary()[bPatchID()].magSf();


        // Update fixedGradient boundary condition on patch A

        vectorField nGradU = 
          - nA*(muFluidA().value() - muFluidB().value())*
            fac::div(Us())().internalField();
            
        nGradU += muFluidA().value()*(UtFs - UtPA)
            *mesh().boundary()[aPatchID()].deltaCoeffs();

        nGradU /= muFluidA().value() + VSMALL;

        if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedGradientCorrectedFvPatchField<vector>& aU =
                refCast<fixedGradientCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientFvPatchField<vector>::typeName
        )
        {
            fixedGradientFvPatchField<vector>& aU =
                refCast<fixedGradientFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );
            
            aU.gradient() = nGradU;
        }
        else
        {
            FatalErrorIn("freeSurface::updateVelocity()")
                << "Bounary condition on " << U().name() 
                    <<  " for freeSurface patch is " 
                    << U().boundaryField()[aPatchID()].type() 
                    << ", instead " 
                    << fixedGradientCorrectedFvPatchField<vector>::typeName 
                    << " or "
                    << fixedGradientFvPatchField<vector>::typeName 
                    << abort(FatalError);
        }
    }
    else
    {
        // Update fixedValue boundary condition on patch B

        vectorField nA = aMesh().faceAreaNormals().internalField();

        scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();

        vectorField UtPA = 
            U().boundaryField()[aPatchID()].patchInternalField();
        UtPA -= nA*(nA & UtPA);

        vectorField UtFs = 
            muFluidA().value()*DnA*UtPA;

        if(!cleanInterface())
        {
            UtFs += surfaceTensionGrad()().internalField();
        }

        vectorField UnFs = 
            nA*phi().boundaryField()[aPatchID()]
           /mesh().boundary()[aPatchID()].magSf();

        Us().internalField() += UnFs - nA*(nA&Us().internalField());
        correctUsBoundaryConditions();


        UtFs -= muFluidA().value()*
            (fac::grad(Us())&aMesh().faceAreaNormals())().internalField();

        UtFs /= muFluidA().value()*DnA + VSMALL;


        Us().internalField() = UtFs + UnFs;
        correctUsBoundaryConditions();

        
        // Update fixedGradient boundary condition on patch A

        vectorField nGradU = 
            -nA*muFluidA().value()*fac::div(Us())().internalField();

        nGradU += 
            muFluidA().value()
           *(UtFs - UtPA)*mesh().boundary()[aPatchID()].deltaCoeffs();


        nGradU /= muFluidA().value() + VSMALL;

        if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedGradientCorrectedFvPatchField<vector>& aU =
                refCast<fixedGradientCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientFvPatchField<vector>::typeName
        )
        {
            fixedGradientFvPatchField<vector>& aU =
                refCast<fixedGradientFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else
        {
            FatalErrorIn("freeSurface::updateVelocity()")
                << "Bounary condition on " << U().name() 
                    <<  " for freeSurface patch is " 
                    << U().boundaryField()[aPatchID()].type() 
                    << ", instead " 
                    << fixedGradientCorrectedFvPatchField<vector>::typeName 
                    << " or "
                    << fixedGradientFvPatchField<vector>::typeName 
                    << abort(FatalError);
        }
    }
}


void freeSurface::updatePressure()
{
    // Correct pressure boundary condition at the free-surface
        
    vectorField nA = mesh().boundary()[aPatchID()].nf();

    if(twoFluids())
    {
        scalarField pA =
            interpolatorBA().faceInterpolate
            (
                p().boundaryField()[bPatchID()]
            );

        const scalarField& K = aMesh().faceCurvatures().internalField();

        if(Pstream::master())
        {
            Info << "Free surface curvature: min = " << min(K)
                << ", max = " << max(K)
                << ", average = " << average(K) << endl << flush;


            if(cleanInterface())
            {   
                pA -= cleanInterfaceSurfTension().value()*(K - average(K));
            }
            else
            {
                scalarField surfTensionK =
                    surfaceTension().internalField()*K;
                
                pA -= surfTensionK - average(surfTensionK);
            }
        }


        if(true)
        {
            scalarField nnGradU = 
                nA&U().boundaryField()[aPatchID()].snGrad();
            
            pA += 2.0*(muFluidA().value() - muFluidB().value())*nnGradU;
        }


        vector R0 = vector(0, 0, 0);

        pA -= (rhoFluidA().value() - rhoFluidB().value())*
            (
                g_.value()
              & (
                    mesh().C().boundaryField()[aPatchID()]
                  - R0
                )
            );


        p().boundaryField()[aPatchID()] == pA;
    }
    else
    {
        vector R0 = vector(0, 0, 0);

        scalarField pA =
          - rhoFluidA().value()*
            (
                g_.value()
              & (
                    mesh().C().boundaryField()[aPatchID()]
                  - R0
                )
            );

        const scalarField& K = aMesh().faceCurvatures().internalField();
        
        if(Pstream::master())
        {
            Info << "Free surface curvature: min = " << min(K)
                << ", max = " << max(K) << ", average = " << average(K) 
                << endl;

            if(cleanInterface())
            {
                pA -= cleanInterfaceSurfTension().value()*(K - average(K));
            }
            else
            {
                scalarField surfTensionK =
                    surfaceTension().internalField()*K;
                
                pA -= surfTensionK - average(surfTensionK);
            }
        }

        if(true)
        {
            scalarField nnGradU =
                nA&U().boundaryField()[aPatchID()].snGrad();

            pA += 2.0*muFluidA().value()*nnGradU;
        }

        p().boundaryField()[aPatchID()] == pA;
    }


    // Set modified pressure at patches with fixed apsolute
    // pressure

    vector R0 = vector(0, 0, 0);

    for (int patchI=0; patchI < p().boundaryField().size(); patchI++)
    {
        if 
        (
            p().boundaryField()[patchI].type()
         == fixedValueFvPatchScalarField::typeName
        )
        {
            if (patchI != aPatchID())
            {
                p().boundaryField()[patchI] ==
                    p().boundaryField()[patchI]
                  - rho().boundaryField()[patchI]*
                    (g_.value()&(mesh().C().boundaryField()[patchI] - R0));
            }
        }
    }
}


void freeSurface::updateSurfaceFlux()
{
    Phis() = fac::interpolate(Us()) & aMesh().Le();
}


void freeSurface::updateSurfactantConcentration()
{
    if(!cleanInterface())
    {
        Info << "Correct surfactant concentration" << endl << flush;
        
        updateSurfaceFlux();        

        // Crate and solve the surfactanta transport equation
        faScalarMatrix CsEqn
        (
            fam::ddt(surfactantConcentration())
          + fam::div(Phis(), surfactantConcentration())
          - fam::laplacian
            (
                surfactant().surfactDiffusion(),
                surfactantConcentration()
            )
        );


        if(surfactant().soluble())
        {
            const scalarField& C =
                mesh().boundary()[aPatchID()]
               .lookupPatchField<volScalarField, scalar>("C");

            areaScalarField Cb
            (
                IOobject
                (
                    "Cb",
                    DB().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                aMesh(),
                dimensioned<scalar>("Cb", dimMoles/dimVolume, 0),
                zeroGradientFaPatchScalarField::typeName
            );

            Cb.internalField() = C;
            Cb.correctBoundaryConditions();

            CsEqn += 
                fam::Sp
                (
                    surfactant().surfactAdsorptionCoeff()*Cb
                  + surfactant().surfactAdsorptionCoeff()
                   *surfactant().surfactDesorptionCoeff(),
                    surfactantConcentration()
                )
              - surfactant().surfactAdsorptionCoeff()
               *Cb*surfactant().surfactSaturatedConc();
        }

        CsEqn.solve();

        Info << "Correct surface tension" << endl;

        surfaceTension() =
            cleanInterfaceSurfTension()
          + surfactant().surfactR()
           *surfactant().surfactT()
           *surfactant().surfactSaturatedConc()
           *log(1.0 - surfactantConcentration()
           /surfactant().surfactSaturatedConc());
        
        if(neg(min(surfaceTension().internalField())))
        {
            FatalErrorIn
            (
                "void freeSurface::correctSurfactantConcentration()"
            ) 
                << "Surface tension has negative value" 
                    << abort(FatalError);
        }
    }
}


void freeSurface::correctUsBoundaryConditions()
{
    forAll(Us().boundaryField(), patchI)
    {
        if
        (
            Us().boundaryField()[patchI].type()
         == calculatedFaPatchVectorField::typeName
        )
        {
            vectorField& pUs = Us().boundaryField()[patchI];

            pUs = Us().boundaryField()[patchI].patchInternalField();

            label ngbPolyPatchID = 
                aMesh().boundary()[patchI].ngbPolyPatchIndex();

            if(ngbPolyPatchID != -1)
            {
                if
                (
                    (
                        U().boundaryField()[ngbPolyPatchID].type() 
                     == slipFvPatchVectorField::typeName
                    )
                 ||
                    (
                        U().boundaryField()[ngbPolyPatchID].type() 
                     == symmetryFvPatchVectorField::typeName
                    )
                )
                {
                    vectorField N = 
                        aMesh().boundary()[patchI].ngbPolyPatchFaceNormals();

                    pUs -= N*(N&pUs);
                }
            }
        }
    }

    Us().correctBoundaryConditions();
}


tmp<scalarField> freeSurface::divPhiSUndulationIndicator()
{
    tmp<scalarField> tUndulation
    (
        new scalarField
        (
            aMesh().nInternalEdges(),
            0.0
        )
    );

    scalarField& undulation = tUndulation();

    areaScalarField div = fac::div(Phis());
    scalarField lnGrad = fac::lnGrad(div)().internalField();
    vectorField grad = fac::grad(div)().internalField();
    tensorField gradGrad = fac::grad(fac::grad(div))().internalField();

    const unallocLabelList& P = aMesh().owner();
    const unallocLabelList& N = aMesh().neighbour();

    const vectorField& Le = aMesh().Le().internalField();
    const scalarField& magLe = aMesh().magLe().internalField();


    labelList infl(aMesh().nFaces(), 0);

    forAll(P, eI)
    {
        const tensorField& curT = aMesh().edgeTransformTensors()[eI];

        scalar min1 = 
            ((Le[eI]*Le[eI])&&transform(curT[0].T(), transform(curT[1], gradGrad[P[eI]])));
            
        scalar max1 =
            ((Le[eI]*Le[eI])&&transform(curT[0].T(), transform(curT[1], gradGrad[N[eI]])));
        
        if(min1*max1<0)
        {
            infl[P[eI]] = 1;
            infl[N[eI]] = 1;
        }
    }


    forAll(undulation, eI)
    {
        if( (infl[P[eI]]==0) && (infl[N[eI]]==0) )
        {

        const tensorField& curT = aMesh().edgeTransformTensors()[eI];

        scalar min = 
            (Le[eI]&transform(curT[0].T(), transform(curT[1], grad[P[eI]])))
            /magLe[eI];
            
        scalar max =
            (Le[eI]&transform(curT[0].T(), transform(curT[1], grad[N[eI]])))
            /magLe[eI];


            if(min < max)
            {
                if(lnGrad[eI] < min)
                {
                    undulation[eI] = lnGrad[eI] - min;
                }
                else if(lnGrad[eI] > max)
                {
                    undulation[eI] = lnGrad[eI] - max;
                }
            }
            else
            {
                if(lnGrad[eI] > min)
                {
                    undulation[eI] = lnGrad[eI] - min;
                }
                else if(lnGrad[eI] < max)
                {
                    undulation[eI] = lnGrad[eI] - max;
                }
            }
        }
    }
        
    return tUndulation;
}


tmp<scalarField> freeSurface::gradGradDivPhiS()
{
    tmp<scalarField> tUndulation
    (
        new scalarField
        (
            aMesh().nInternalEdges(),
            0.0
        )
    );

    scalarField& undulation = tUndulation();

    areaScalarField div = fac::div(Phis());
    scalarField lnGrad = fac::lnGrad(div)().internalField();
    vectorField grad = fac::grad(div)().internalField();
    areaTensorField gradGrad = fac::grad(fac::grad(div));

    edgeVectorField me = aMesh().Le()/aMesh().magLe();

    undulation = ((me*me)&&fac::interpolate(gradGrad))().internalField();
        
    return tUndulation;
}


vector freeSurface::totalPressureForce()
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    const scalarField& P = p().boundaryField()[aPatchID()];

    vectorField pressureForces = S*P*n;

    return sum(pressureForces);
}


vector freeSurface::totalViscousForce()
{
    const scalarField& S = aMesh().S();
    const vectorField& n = aMesh().faceAreaNormals().internalField();

    vectorField nGradU =
        U().boundaryField()[aPatchID()].snGrad();

    vectorField viscousForces = 
      - muFluidA().value()*S*
        (
            transform(I + n*n, nGradU) 
          + fac::grad(Us()&aMesh().faceAreaNormals())().internalField()
        );

    return sum(viscousForces);
}


vector freeSurface::totalSurfaceTensionForce()
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    const scalarField& K = aMesh().faceCurvatures().internalField();


    vectorField surfTensionForces = S*n;


    if(cleanInterface())
    {   
        surfTensionForces *= cleanInterfaceSurfTension().value()*K;
    }
    else
    {
        surfTensionForces *= surfaceTension().internalField()*K;
    }


    return sum(surfTensionForces);
}


tmp<scalarField> freeSurface::undulationIndicator()
{
    tmp<scalarField> tUndulation
    (
        new scalarField
        (
            aMesh().nFaces(),
            0.0
        )
    );

    scalarField& undulation = tUndulation();

    primitivePatchInterpolation patchInterpolator
        (
            mesh().boundaryMesh()[aPatchID()]
        );

    undulation = 
        asin(mag(
            aMesh().faceAreaNormals().internalField()^
            patchInterpolator.pointToFaceInterpolate
            (
                aMesh().pointAreaNormals()
            )
        ))*180.0/M_PI;

    return tUndulation;
}


void freeSurface::smooth()
{
    if
    (
        Pstream::master()
    )
    {
        Info << "Smoothing free-surface faces...";

        const faceList& faces = aMesh().faces();

        scalarField deltaH = facesDisplacementDir()&
        (
            aMesh().centres().internalField()
          - controlPoints()
        );

        pointField displacement = pointDisplacement(deltaH);

        const pointField& points = aMesh().points();


        pointField newPoints = points + displacement;

        scalarField sweptVol(faces.size(), 0.0);

        forAll(faces, faceI)
        {
            sweptVol[faceI] = -faces[faceI].sweptVol(points, newPoints);
        }

        vectorField faceArea(faces.size(), vector::zero);

        forAll (faceArea, faceI)
        {
            faceArea[faceI] = faces[faceI].normal(newPoints);
        }

        forAll(deltaH, faceI)
        {
            deltaH[faceI] = sweptVol[faceI]/
                (faceArea[faceI] & facesDisplacementDir()[faceI]);
        }


        displacement = pointDisplacement(deltaH);

        // Move only free-surface points

        pointField newMeshPoints = mesh().points();

        const labelList& meshPointsA = 
            mesh().boundaryMesh()[aPatchID()].meshPoints();

        forAll (displacement, pointI)
        {
            newMeshPoints[meshPointsA[pointI]] += displacement[pointI]; 
        }

        if(twoFluids_)
        {
            const labelList& meshPointsB = 
                mesh().boundaryMesh()[bPatchID_].meshPoints();
            
            pointField displacementB =
                interpolatorAB().pointInterpolate
                (
                    displacement
                );
        
            forAll (displacementB, pointI)
            {
                newMeshPoints[meshPointsB[pointI]] += displacementB[pointI]; 
            }
        }

        mesh().movePoints(newMeshPoints);
       
        aMesh().movePoints(mesh().points());


        // Update total displacement field

        /*if(totalDisplacementPtr_ && (curTimeIndex_ < DB().timeIndex()))
        {
            FatalErrorIn("freeSurface::movePoints()")
                << "Total displacement of free surface points "
                << "from previous time step is not absorbed by the mesh."
                << abort(FatalError);
        }
        else*/ if (curTimeIndex_ < DB().timeIndex())
        {
            totalDisplacement() = displacement;
        
            curTimeIndex_ = DB().timeIndex();
        }
        else
        {
            totalDisplacement() += displacement;
        }

        Info << "done" << endl << flush;
    }
}


void freeSurface::initializeControlPointsPosition()
{
    scalarField deltaH = scalarField(aMesh().nFaces(), 0.0);

    pointField displacement = pointDisplacement(deltaH);


    const faceList& faces = aMesh().faces();
    const pointField& points = aMesh().points();


    pointField newPoints = points + displacement;

    scalarField sweptVol(faces.size(), 0.0);

    forAll(faces, faceI)
    {
        sweptVol[faceI] = -faces[faceI].sweptVol(points, newPoints);
    }

    vectorField faceArea(faces.size(), vector::zero);

    forAll (faceArea, faceI)
    {
        faceArea[faceI] = faces[faceI].normal(newPoints);
    }

    forAll(deltaH, faceI)
    {
        deltaH[faceI] = sweptVol[faceI]/
            (faceArea[faceI] & facesDisplacementDir()[faceI]);
    }

    pointDisplacement(deltaH);
}


scalar freeSurface::maxCourantNumber()
{
    scalar CoNum = 0;

    if(Pstream::master())
    {
        if(cleanInterface())
        {
            const scalarField& dE =aMesh().lPN();

            CoNum = max
            (
                DB().deltaT().value()/
                sqrt
                (
                    rhoFluidA().value()*dE*dE*dE/
                    2.0/M_PI/(cleanInterfaceSurfTension().value() + SMALL)
                )
            );
        }
        else
        {
            scalarField sigmaE = 
                linearEdgeInterpolate(surfaceTension())().internalField()
              + SMALL;

            const scalarField& dE =aMesh().lPN();

            CoNum = max
            (
                DB().deltaT().value()/
                sqrt
                (
                    rhoFluidA().value()*dE*dE*dE/
                    2.0/M_PI/sigmaE
                )
            );
        }
    }

    return CoNum;
}


tmp<scalarField> freeSurface::checkFaceFlatness()
{
    const pointField& points = aMesh().points();
    const faceList& faces = aMesh().faces();
    const pointField& Cf = aMesh().centres().internalField();
    const scalarField& Sf = aMesh().S();

    tmp<scalarField> tFlatness
    (
        new scalarField
        (
            aMesh().nFaces(),
            1.0
        )
    );
    scalarField& flatness = tFlatness();


    forAll(flatness, faceI)
    {
        const face& f  = faces[faceI];        

        if (f.size() > 3 && Sf[faceI] > VSMALL)
        {
            const point& fc = Cf[faceI];

            // Calculate the sum of magnitude of areas and compare to magnitude
            // of sum of areas.
            scalar sumA = 0.0;

            forAll(f, pointI)
            {
                const point& thisPoint = points[f[pointI]];
                const point& nextPoint = points[f.nextLabel(pointI)];

                // Triangle around Cf.
                vector n = 0.5*((nextPoint - thisPoint)^(fc - thisPoint));
                sumA += mag(n);
            }
            
            flatness[faceI] = Sf[faceI] / (sumA+VSMALL);
        }
    }

    Info << "Free-surface faces flatness, min: " << min(flatness)
        << ", average: " << average(flatness) 
        << ", sd: " 
        << std::sqrt(sum(sqr(flatness-1))/sum(sqr(flatness))) << endl;

    return tFlatness;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
