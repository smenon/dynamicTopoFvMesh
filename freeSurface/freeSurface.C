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
#include "transformField.H"

#include "emptyFaPatch.H"
#include "wedgeFaPatch.H"
#include "wallFvPatch.H"

#include "EulerDdtScheme.H"
#include "CrankNicholsonDdtScheme.H"
#include "backwardDdtScheme.H"

#include "tetFemMatrices.H"
#include "tetPointFields.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"
#include "fixedValueTetPolyPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "twoDPointCorrector.H"

#include "cyclicPolyPatch.H"
#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "zeroGradientCorrectedFvPatchFields.H"
#include "fixedGradientCorrectedFvPatchFields.H"
#include "fixedValueCorrectedFvPatchFields.H"

#include "mapPolyMesh.H"
#include "setMotionBC.H"

#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(freeSurface, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void freeSurface::clearOut() const
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
    deleteDemandDrivenData(phiPtr_);
    deleteDemandDrivenData(ddtPhiPtr_);
    deleteDemandDrivenData(surfactConcPtr_);
    deleteDemandDrivenData(surfaceTensionPtr_);
    deleteDemandDrivenData(surfactantPtr_);
    deleteDemandDrivenData(fluidIndicatorPtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

freeSurface::freeSurface
(
    dynamicFvMesh& m,
    volVectorField& Ub,
    volScalarField& Pb,
    const surfaceScalarField& sfPhi
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
            IOobject::NO_WRITE,
            false
        )
    ),
    MeshObject<fvMesh, freeSurface>(m),
    mesh_(m),
    U_(Ub),
    p_(Pb),
    phi_(sfPhi),
    interfacePatch_(word(lookup("interfacePatch"))),
    orderedPatches_(true),
    curTimeIndex_(Ub.mesh().time().timeIndex()),
    twoFluids_(lookup("twoFluids")),
    normalMotionDir_(lookup("normalMotionDir")),
    motionDir_(0, 0, 0),
    cleanInterface_(lookup("cleanInterface")),
    aPatchID_(-1),
    bPatchID_(-1),
    muFluidA_(lookup("muFluidA")),
    rhoFluidA_(lookup("rhoFluidA")),
    condFluidA_(lookup("condFluidA")),
    CpFluidA_(lookup("CpFluidA")),
    muFluidB_(lookup("muFluidB")),
    rhoFluidB_(lookup("rhoFluidB")),
    condFluidB_(lookup("condFluidB")),
    CpFluidB_(lookup("CpFluidB")),
    g_(lookup("g")),
    cleanInterfaceSurfTension_(lookup("surfaceTension")),
    fixedFreeSurfacePatches_(lookup("fixedFreeSurfacePatches")),
    pointNormalsCorrectionPatches_(lookup("pointNormalsCorrectionPatches")),
    nFreeSurfCorr_(readInt(lookup("nFreeSurfaceCorrectors"))),
    interpolatorABPtr_(NULL),
    interpolatorBAPtr_(NULL),
    controlPointsPtr_(NULL),
    motionPointsMaskPtr_(NULL),
    pointsDisplacementDirPtr_(NULL),
    facesDisplacementDirPtr_(NULL),
    areaCentresPtr_(NULL),
    totalDisplacementPtr_(NULL),
    aMeshPtr_(NULL),
    UsPtr_(NULL),
    phisPtr_(NULL),
    phiPtr_(NULL),
    ddtPhiPtr_(NULL),
    surfactConcPtr_(NULL),
    surfaceTensionPtr_(NULL),
    surfactantPtr_(NULL),
    fluidIndicatorPtr_(NULL),
    rhoPtr_(NULL),
    muPtr_(NULL)
{
    if (twoFluids_)
    {
        // Obtain the shadow patch name
        shadowPatch_ = word(lookup("shadowPatch"));
    }

    // Read motion direction
    if (!normalMotionDir_)
    {
        motionDir_ = vector(lookup("motionDir"));
        motionDir_ /= mag(motionDir_) + SMALL;
    }

    // Loop through all boundary patches and determine the interface patchID
    const fvBoundaryMesh& bdy = mesh().boundary();

    forAll (bdy, patchI)
    {
        if (bdy[patchI].name() == interfacePatch_)
        {
            aPatchID_ = patchI;
        }

        if (twoFluids_)
        {
            if (bdy[patchI].name() == shadowPatch_)
            {
                bPatchID_ = patchI;
            }
        }
    }

    if (aPatchID_ == -1)
    {
        FatalErrorIn("freeSurface::freeSurface()")
            << "Patch with name " << interfacePatch_ << " is not defined."
            << abort(FatalError);
    }

    if (twoFluids_)
    {
        if (bPatchID_ == -1)
        {
            FatalErrorIn("freeSurface::freeSurface()")
                << "Patch with name " << shadowPatch_ << " is not defined."
                << abort(FatalError);
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
            )   << "Patch name for point normals correction does not exist"
                << abort(FatalError);
        }

        correction[patchID] = true;
    }

    // Read free-surface points total displacement if present
    readTotalDisplacement();

    // Read control points positions if present
    controlPoints();

    // Set interface flux
    phi();
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

freeSurface::~freeSurface()
{
    clearOut();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Interpolate from A to B using either interpolator or direct copy
template<class Type>
tmp<Field<Type> > Foam::freeSurface::pointInterpolateAB
(
    const Field<Type>& field
)
{
    if (orderedPatches_)
    {
        label srcSize = mesh_.boundaryMesh()[aPatchID_].nPoints();
        label patchSize = mesh_.boundaryMesh()[bPatchID_].nPoints();

        if (patchSize != field.size() || patchSize != srcSize)
        {
            FatalErrorIn("freeSurface::pointInterpolateAB()")
                << "Sizes don't match: " << nl
                << "  srcSize: " << srcSize << nl
                << "  patchSize: " << patchSize << nl
                << "  Field size: " << field.size() << nl
                << "  Master: " << aPatchID() << nl
                << "  Slave: " << bPatchID() << nl
                << abort(FatalError);
        }

        tmp<Field<Type> > tresult
        (
            new Field<Type>
            (
                patchSize,
                pTraits<Type>::zero
            )
        );

        Field<Type>& result = tresult();

        // Fetch face-lists
        const faceList& mFaces = mesh_.boundaryMesh()[aPatchID_].localFaces();
        const faceList& sFaces = mesh_.boundaryMesh()[bPatchID_].localFaces();

        // Assume coupled face-point ordering
        forAll(mFaces, faceI)
        {
            const face& mFace = mFaces[faceI];
            const face& sFace = sFaces[faceI];

            label fS = sFace.size();

            forAll(mFace, pointI)
            {
                label masterIndex = mFace[pointI];
                label slaveIndex = sFace[(fS - pointI) % fS];

                // Copy value
                result[slaveIndex] = field[masterIndex];
            }
        }

        return tresult;
    }

    // Return conventional interpolation
    return interpolatorAB().pointInterpolate(field);
}


// Interpolate from B to A using either interpolator or direct copy
template<class Type>
tmp<Field<Type> > Foam::freeSurface::pointInterpolateBA
(
    const Field<Type>& field
)
{
    if (orderedPatches_)
    {
        label srcSize = mesh_.boundaryMesh()[bPatchID_].nPoints();
        label patchSize = mesh_.boundaryMesh()[aPatchID_].nPoints();

        if (patchSize != field.size() || patchSize != srcSize)
        {
            FatalErrorIn("freeSurface::pointInterpolateBA()")
                << "Sizes don't match: " << nl
                << "  srcSize: " << srcSize << nl
                << "  patchSize: " << patchSize << nl
                << "  Field size: " << field.size() << nl
                << "  Master: " << aPatchID() << nl
                << "  Slave: " << bPatchID() << nl
                << abort(FatalError);
        }

        tmp<Field<Type> > tresult
        (
            new Field<Type>
            (
                patchSize,
                pTraits<Type>::zero
            )
        );

        Field<Type>& result = tresult();

        // Fetch face-lists
        const faceList& mFaces = mesh_.boundaryMesh()[aPatchID_].localFaces();
        const faceList& sFaces = mesh_.boundaryMesh()[bPatchID_].localFaces();

        // Assume coupled face-point ordering
        forAll(mFaces, faceI)
        {
            const face& mFace = mFaces[faceI];
            const face& sFace = sFaces[faceI];

            label fS = sFace.size();

            forAll(mFace, pointI)
            {
                label masterIndex = mFace[pointI];
                label slaveIndex = sFace[(fS - pointI) % fS];

                // Copy value
                result[masterIndex] = field[slaveIndex];
            }
        }

        return tresult;
    }

    // Return conventional interpolation
    return interpolatorBA().pointInterpolate(field);
}


// Interpolate from A to B using either interpolator or direct copy
template<class Type>
tmp<Field<Type> > Foam::freeSurface::faceInterpolateAB
(
    const Field<Type>& field
)
{
    if (orderedPatches_)
    {
        label srcSize = mesh_.boundaryMesh()[aPatchID_].size();
        label patchSize = mesh_.boundaryMesh()[bPatchID_].size();

        if (patchSize != field.size() || patchSize != srcSize)
        {
            FatalErrorIn("freeSurface::faceInterpolateAB()")
                << "Sizes don't match: " << nl
                << "  srcSize: " << srcSize << nl
                << "  patchSize: " << patchSize << nl
                << "  Field size: " << field.size() << nl
                << "  Master: " << aPatchID() << nl
                << "  Slave: " << bPatchID() << nl
                << abort(FatalError);
        }

        tmp<Field<Type> > tresult
        (
            new Field<Type>
            (
                patchSize,
                pTraits<Type>::zero
            )
        );

        Field<Type>& result = tresult();

        // Assume ordered faces, direct copy
        result = field;

        return tresult;
    }

    // Return conventional interpolation
    return interpolatorAB().faceInterpolate(field);
}


// Interpolate from B to A using either interpolator or direct copy
template<class Type>
tmp<Field<Type> > Foam::freeSurface::faceInterpolateBA
(
    const Field<Type>& field
)
{
    if (orderedPatches_)
    {
        label srcSize = mesh_.boundaryMesh()[bPatchID_].size();
        label patchSize = mesh_.boundaryMesh()[aPatchID_].size();

        if (patchSize != field.size() || patchSize != srcSize)
        {
            FatalErrorIn("freeSurface::faceInterpolateBA()")
                << "Sizes don't match: " << nl
                << "  srcSize: " << srcSize << nl
                << "  patchSize: " << patchSize << nl
                << "  Field size: " << field.size() << nl
                << "  Master: " << aPatchID() << nl
                << "  Slave: " << bPatchID() << nl
                << abort(FatalError);
        }

        tmp<Field<Type> > tresult
        (
            new Field<Type>
            (
                patchSize,
                pTraits<Type>::zero
            )
        );

        Field<Type>& result = tresult();

        // Assume ordered faces, direct copy
        result = field;

        return tresult;
    }

    // Return conventional interpolation
    return interpolatorBA().faceInterpolate(field);;
}


void freeSurface::updateDisplacementDirections() const
{
    if (normalMotionDir())
    {
        // Update point displacement correction
        pointsDisplacementDir() = aMesh().pointAreaNormals();

        // Update face displacement direction
        facesDisplacementDir() = aMesh().faceAreaNormals().internalField();

        // Correction of control points postion
        const vectorField& Cf = aMesh().areaCentres().internalField();

        controlPoints() =
            facesDisplacementDir()
          * (facesDisplacementDir() & (controlPoints() - Cf))
          + Cf;
    }
}


bool freeSurface::predictPoints()
{
    for
    (
        int freeSurfCorr=0;
        freeSurfCorr<nFreeSurfCorr_;
        freeSurfCorr++
    )
    {
        movePoints(phi_.boundaryField()[aPatchID()]);
    }

    return true;
}


bool freeSurface::correctPoints()
{
    for
    (
        int freeSurfCorr=0;
        freeSurfCorr<nFreeSurfCorr_;
        freeSurfCorr++
    )
    {
        movePoints(phi_.boundaryField()[aPatchID()]);
    }

    return true;
}


bool freeSurface::movePoints(const scalarField& interfacePhi)
{
    pointField newMeshPoints = mesh().points();

    scalarField sweptVolCorr =
    (
        interfacePhi
      - fvc::meshPhi(rho(),U())().boundaryField()[aPatchID()]
    );

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
        if (DB().timeIndex() == 1)
        {
            sweptVolCorr *= DB().deltaT().value();
        }
        else
        {
            sweptVolCorr *= (2.0/3.0)*DB().deltaT().value();
        }
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

    scalarField deltaH = sweptVolCorr/(Sf*(Nf & facesDisplacementDir()));

    pointField displacement = pointDisplacement(deltaH);

    // Move only free-surface points
    const labelList& meshPointsA =
    (
        mesh().boundaryMesh()[aPatchID()].meshPoints()
    );

    forAll (displacement, pointI)
    {
        newMeshPoints[meshPointsA[pointI]] += displacement[pointI];
    }

    if (twoFluids_)
    {
        // Interpolate directly or using patchToPatchInterpolation
        pointField displacementB = pointInterpolateAB(displacement);

        const labelList& meshPtsB =
        (
            mesh().boundaryMesh()[bPatchID_].meshPoints()
        );

        forAll (displacementB, pointI)
        {
            newMeshPoints[meshPtsB[pointI]] += displacementB[pointI];
        }
    }

    // Update total displacement field
    if (totalDisplacementPtr_ && (curTimeIndex_ < DB().timeIndex()))
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Total displacement of free surface points "
            << "from previous time step is not absorbed by the mesh."
            << abort(FatalError);
    }
    else if (curTimeIndex_ < DB().timeIndex())
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

    // Correct points for cyclics
    correctCyclics(newMeshPoints);

    mesh().movePoints(newMeshPoints);
    aMesh().movePoints(mesh().points());

    // Move correctedFvPatchField fvSubMeshes
    moveFvSubMeshes();

    return true;
}


// Correct points for cyclics
void freeSurface::correctCyclics(pointField& points) const
{
    // Correct points on cyclics
    const polyBoundaryMesh& boundary = mesh().boundaryMesh();

    forAll(boundary, patchI)
    {
        if (!isA<cyclicPolyPatch>(boundary[patchI]))
        {
            continue;
        }

        // Cast to cyclic
        const cyclicPolyPatch& cyclicPatch =
        (
            refCast<const cyclicPolyPatch>(boundary[patchI])
        );

        bool translate =
        (
            cyclicPatch.transform()
         == cyclicPolyPatch::TRANSLATIONAL
        );

        label patchStart = boundary[patchI].start();
        label halfSize = (boundary[patchI].size() / 2);

        for (label faceI = 0; faceI < halfSize; faceI++)
        {
            label half0Index = (patchStart + faceI);
            label half1Index = (patchStart + halfSize + faceI);

            const face& half0Face = mesh().faces()[half0Index];
            const face& half1Face = mesh().faces()[half1Index];

            label fS = half0Face.size();

            forAll(half0Face, pointI)
            {
                label masterIndex = half0Face[pointI];
                label slaveIndex = half1Face[(fS - pointI) % fS];

                if (translate)
                {
                    points[slaveIndex] =
                    (
                        points[masterIndex]
                      + cyclicPatch.separationVector()
                    );
                }
                else
                {
                    points[slaveIndex] =
                    (
                        cyclicPatch.transform
                        (
                            points[masterIndex],
                            faceI
                        )
                    );
                }
            }
        }
    }
}


bool freeSurface::moveMeshPointsForOldFreeSurfDisplacement()
{
    if (!totalDisplacementPtr_)
    {
        return false;
    }

    pointField newPoints = mesh().points();

    const labelList& meshPtsA = mesh().boundaryMesh()[aPatchID()].meshPoints();

    forAll (meshPtsA, pointI)
    {
        newPoints[meshPtsA[pointI]] -= totalDisplacement()[pointI];
    }

    if (twoFluids())
    {
        // Interpolate displacement to the shadow patch
        pointField dispB = pointInterpolateAB(totalDisplacement());

        const labelList& meshPtsB =
        (
            mesh().boundaryMesh()[bPatchID()].meshPoints()
        );

        forAll(meshPtsB,pointI)
        {
            newPoints[meshPtsB[pointI]] -= dispB[pointI];
        }
    }

    // Now set boundary conditions for the motion-solver
    setMotionBC(mesh(), aPatchID(), totalDisplacement());

    if (twoFluids())
    {
        // Interpolate displacement to the shadow patch
        pointField dispB = pointInterpolateAB(totalDisplacement());

        setMotionBC(mesh(), bPatchID(), dispB);
    }

    twoDPointCorrector twoDPointCorr(mesh());
    twoDPointCorr.correctPoints(newPoints);

    // Correct points for cyclics
    correctCyclics(newPoints);

    // Move mesh points to old positions
    mesh().movePoints(newPoints);

    deleteDemandDrivenData(totalDisplacementPtr_);

    // Now update the mesh with motion / topology changes.
    bool meshChanged = mesh().update();

    aMesh().movePoints(mesh().points());

    // Move correctedFvPatchField fvSubMeshes
    moveFvSubMeshes();

    if (meshChanged)
    {
        //reInitializeControlPointsPosition();
    }

    return meshChanged;
}


//- Move correctedFvPatchField fvSubMeshes
void freeSurface::moveFvSubMeshes() const
{
    // Move correctedFvPatchField fvSubMeshes
    forAll(U().boundaryField(), patchI)
    {
        if
        (
            (
                U().boundaryField()[patchI].type()
             == "fixedGradientCorrected"
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == "fixedValueCorrected"
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == "zeroGradientCorrected"
            )
        )
        {
            correctedFvPatchField<vector>& pU =
            (
                refCast<correctedFvPatchField<vector> >
                (
                    U_.boundaryField()[patchI]
                )
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
             == "fixedGradientCorrected"
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == "fixedValueCorrected"
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == "zeroGradientCorrected"
            )
        )
        {
            correctedFvPatchField<scalar>& pP =
            (
                refCast<correctedFvPatchField<scalar> >
                (
                    p_.boundaryField()[patchI]
                )
            );

            pP.movePatchSubMesh();
        }
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
    if (twoFluids())
    {
        vectorField nA = mesh().boundary()[aPatchID()].nf();
        vectorField nB = mesh().boundary()[bPatchID()].nf();

        scalarField DnB = faceInterpolateBA
        (
            mesh().boundary()[bPatchID()].deltaCoeffs()
        );

        scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();

        vectorField UtPA = U().boundaryField()[aPatchID()].patchInternalField();

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

        vectorField UtPB = faceInterpolateBA
        (
            U().boundaryField()[bPatchID()].patchInternalField()()
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

            UtPB += faceInterpolateBA(bU.corrVecGrad());
        }

        UtPB -= nA*(nA & UtPB);

        vectorField UtFs =
        (
            muFluidA().value()*DnA*UtPA
          + muFluidB().value()*DnB*UtPB
        );

        vectorField UnFs =
            nA*phi_.boundaryField()[aPatchID()]
           /mesh().boundary()[aPatchID()].magSf();

        Us().internalField() += UnFs - nA*(nA&Us().internalField());
        correctUsBoundaryConditions();

        UtFs -= (muFluidA().value() - muFluidB().value())*
            (fac::grad(Us())&aMesh().faceAreaNormals())().internalField();

        vectorField tangentialSurfaceTensionForce(nA.size(), vector::zero);

        if (!cleanInterface())
        {
            tangentialSurfaceTensionForce =
                surfaceTensionGrad()().internalField();
        }
        else
        {
            vectorField surfaceTensionForce =
                cleanInterfaceSurfTension().value()
               *fac::edgeIntegrate
                (
                    aMesh().Le()*aMesh().edgeLengthCorrection()
                )().internalField();

            tangentialSurfaceTensionForce =
                surfaceTensionForce
              - cleanInterfaceSurfTension().value()
               *aMesh().faceCurvatures().internalField()*nA;
        }

        UtFs += tangentialSurfaceTensionForce;

        UtFs /= muFluidA().value()*DnA + muFluidB().value()*DnB + VSMALL;

        Us().internalField() = UnFs + UtFs;
        correctUsBoundaryConditions();

        // Store old-time velocity field U()
        U().oldTime();

        U().boundaryField()[bPatchID()] ==
            faceInterpolateAB(UtFs)
          + nB*fvc::meshPhi(rho(),U())().boundaryField()[bPatchID()]/
            mesh().boundary()[bPatchID()].magSf();

        if
        (
            p().boundaryField()[bPatchID()].type()
         == fixedGradientFvPatchField<scalar>::typeName
        )
        {
            fixedGradientFvPatchField<scalar>& pB =
                refCast<fixedGradientFvPatchField<scalar> >
                (
                    p().boundaryField()[bPatchID()]
                );

            pB.gradient() =
               - rhoFluidB().value()
                *(
                     nB&fvc::ddt(U())().boundaryField()[bPatchID()]
                 );
        }

        // Update fixedGradient boundary condition on patch A
        vectorField nGradU =
            muFluidB().value()*(UtPB - UtFs)*DnA
          + tangentialSurfaceTensionForce
          - muFluidA().value()*nA*fac::div(Us())().internalField()
          + (muFluidB().value() - muFluidA().value())
           *(fac::grad(Us())().internalField()&nA);

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
                << "Boundary condition on " << U().name()
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
        vectorField nA = aMesh().faceAreaNormals().internalField();

        vectorField UnFs =
            nA*phi_.boundaryField()[aPatchID()]
           /mesh().boundary()[aPatchID()].magSf();

        // Correct normal component of free-surface velocity
        Us().internalField() += UnFs - nA*(nA&Us().internalField());
        correctUsBoundaryConditions();

        vectorField tangentialSurfaceTensionForce(nA.size(), vector::zero);

        if (!cleanInterface())
        {
             tangentialSurfaceTensionForce =
                 surfaceTensionGrad()().internalField();
        }
        else
        {
            vectorField surfaceTensionForce =
                cleanInterfaceSurfTension().value()
               *fac::edgeIntegrate
                (
                    aMesh().Le()*aMesh().edgeLengthCorrection()
                )().internalField();

            tangentialSurfaceTensionForce =
                surfaceTensionForce
              - cleanInterfaceSurfTension().value()
               *aMesh().faceCurvatures().internalField()*nA;

            if (muFluidA().value() < SMALL)
            {
                tangentialSurfaceTensionForce = vector::zero;
            }
        }

        vectorField tnGradU =
            tangentialSurfaceTensionForce/(muFluidA().value() + VSMALL)
          - (fac::grad(Us())&aMesh().faceAreaNormals())().internalField();

        vectorField UtPA =
            U().boundaryField()[aPatchID()].patchInternalField();
        UtPA -= nA*(nA & UtPA);

        scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();

        vectorField UtFs = UtPA + tnGradU/DnA;

        Us().internalField() = UtFs + UnFs;
        correctUsBoundaryConditions();

        vectorField nGradU =
            tangentialSurfaceTensionForce/(muFluidA().value() + VSMALL)
          - nA*fac::div(Us())().internalField()
          - (fac::grad(Us())().internalField()&nA);

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
                << "Boundary condition on " << U().name()
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

    if (twoFluids())
    {
        scalarField pA = faceInterpolateBA(p().boundaryField()[bPatchID()]);

        const scalarField& K = aMesh().faceCurvatures().internalField();

        Info << "Free surface curvature: min = " << gMin(K)
             << ", max = " << gMax(K)
             << ", average = " << gAverage(K) << endl << flush;

        if (cleanInterface())
        {
//             pA -= cleanInterfaceSurfTension().value()*(K - gAverage(K));
            pA -= cleanInterfaceSurfTension().value()*K;
        }
        else
        {
            scalarField surfTensionK = surfaceTension().internalField()*K;

            pA -= surfTensionK - gAverage(surfTensionK);
        }

        pA -= 2.0*(muFluidA().value() - muFluidB().value())
            *fac::div(Us())().internalField();

        vector R0 = gAverage(mesh().C().boundaryField()[aPatchID()]);

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
        vector R0 = gAverage(mesh().C().boundaryField()[aPatchID()]);

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

        Info << "Free surface curvature: min = " << gMin(K)
             << ", max = " << gMax(K) << ", average = " << gAverage(K)
             << endl;

        if (cleanInterface())
        {
//             pA -= cleanInterfaceSurfTension().value()*(K - gAverage(K));
            pA -= cleanInterfaceSurfTension().value()*K;
        }
        else
        {
            scalarField surfTensionK = surfaceTension().internalField()*K;

            pA -= surfTensionK - gAverage(surfTensionK);
        }

        pA -= 2.0*muFluidA().value()*fac::div(Us())().internalField();

        p().boundaryField()[aPatchID()] == pA;
    }


    // Set modified pressure at patches with fixed absolute pressure
    vector R0 = gAverage(mesh().C().boundaryField()[aPatchID()]);

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
                    // 0
                  - rho().boundaryField()[patchI]
                   *(g_.value()&(mesh().C().boundaryField()[patchI] - R0));
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
    if (!cleanInterface())
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


        if (surfactant().soluble())
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

        if (neg(min(surfaceTension().internalField())))
        {
            FatalErrorIn
            (
                "void freeSurface::correctSurfactantConcentration()"
            )
                << "Surface tension is negative"
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

            if (ngbPolyPatchID != -1)
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


vector freeSurface::totalPressureForce() const
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    const scalarField& P = p().boundaryField()[aPatchID()];

    vectorField pressureForces = S*P*n;

    return gSum(pressureForces);
}


vector freeSurface::totalViscousForce() const
{
    const scalarField& S = aMesh().S();
    const vectorField& n = aMesh().faceAreaNormals().internalField();

    vectorField nGradU =
        U().boundaryField()[aPatchID()].snGrad();

    vectorField viscousForces =
      - muFluidA().value()*S
       *(
            nGradU
          + (fac::grad(Us())().internalField()&n)
          - (n*fac::div(Us())().internalField())
        );

    return gSum(viscousForces);
}


vector freeSurface::totalSurfaceTensionForce() const
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    const scalarField& K = aMesh().faceCurvatures().internalField();

    vectorField surfTensionForces(n.size(), vector::zero);

    if(cleanInterface())
    {
        surfTensionForces =
            S*cleanInterfaceSurfTension().value()
           *fac::edgeIntegrate
            (
                aMesh().Le()*aMesh().edgeLengthCorrection()
            )().internalField();
    }
    else
    {
        surfTensionForces *= surfaceTension().internalField()*K;
    }

    return gSum(surfTensionForces);
}


void freeSurface::initializeControlPointsPosition() const
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

//     const scalarField& S = aMesh().S().field();
//     const vectorField& N = aMesh().faceAreaNormals().internalField();

//     forAll(deltaH, faceI)
//     {
//         deltaH[faceI] = sweptVol[faceI]/
//             (S[faceI]*(N[faceI] & facesDisplacementDir()[faceI]));
//     }

//     displacement = pointDisplacement(deltaH);
}


scalar freeSurface::maxCourantNumber()
{
    scalar CoNum = 0;

    if (cleanInterface())
    {
        const scalarField& dE =aMesh().lPN();

        CoNum = gMax
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

        CoNum = gMax
        (
            DB().deltaT().value()/
            sqrt
            (
                rhoFluidA().value()*dE*dE*dE/
                2.0/M_PI/sigmaE
            )
        );
    }

    return CoNum;
}


void freeSurface::updateProperties()
{
    muFluidA_ = dimensionedScalar(this->lookup("muFluidA"));
    muFluidB_ = dimensionedScalar(this->lookup("muFluidB"));
    rhoFluidA_ = dimensionedScalar(this->lookup("rhoFluidA"));
    rhoFluidB_ = dimensionedScalar(this->lookup("rhoFluidB"));

    g_ = dimensionedVector(this->lookup("g"));

    cleanInterfaceSurfTension_ =
        dimensionedScalar(this->lookup("surfaceTension"));
}


// Update for mesh motion
bool Foam::freeSurface::movePoints() const
{
    return true;
}


// Update on topology change
bool Foam::freeSurface::updateMesh(const mapPolyMesh& mpm) const
{
    if (debug)
    {
        Info << "Clearing out freeSurface after topology change" << endl;
    }

    // Wipe out demand-driven data
    clearOut();

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
