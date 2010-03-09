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

\*----------------------------------------------------------------------------*/

#include "freeSurface.H"

#include "wedgeFaPatch.H"
#include "fixedGradientFvPatchFields.H"
#include "fixedGradientCorrectedFvPatchFields.H"
#include "fixedValueCorrectedFvPatchFields.H"
#include "mapPolyMesh.H"
#include "setMotionBC.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

defineTypeNameAndDebug(freeSurface, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::freeSurface::freeSurface
(
    dynamicFvMesh& m,
    volVectorField& U,
    volScalarField& p,
    surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "freeSurfaceProperties",
            U.mesh().time().constant(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    MeshObject<fvMesh, freeSurface>(m),
    mesh_(m),
    U_(U),
    p_(p),
    phi_(phi),
    interfacePatch_(word(lookup("interfacePatch"))),
    curTimeIndex_(U.mesh().time().timeIndex()),
    twoFluids_(lookup("twoFluids")),
    normalMotionDir_(lookup("normalMotionDir")),
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

    // Force creation of demand-driven data
    // (and read from disk if this is a restart)
    //makeControlPoints();
    //readTotalDisplacement();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::freeSurface::~freeSurface()
{
    // Clear out other demand-driven data
    clearOut();
}

void Foam::freeSurface::clearOut() const
{
    deleteDemandDrivenData(interpolatorABPtr_);
    deleteDemandDrivenData(interpolatorBAPtr_);
    deleteDemandDrivenData(controlPointsPtr_);
    deleteDemandDrivenData(motionPointsMaskPtr_);
    deleteDemandDrivenData(pointsDisplacementDirPtr_);
    deleteDemandDrivenData(facesDisplacementDirPtr_);
    deleteDemandDrivenData(areaCentresPtr_);
    deleteDemandDrivenData(totalDisplacementPtr_);
    deleteDemandDrivenData(aMeshPtr_);
    deleteDemandDrivenData(UsPtr_);
    deleteDemandDrivenData(phisPtr_);
    deleteDemandDrivenData(surfactConcPtr_);
    deleteDemandDrivenData(surfaceTensionPtr_);
    deleteDemandDrivenData(surfactantPtr_);
    deleteDemandDrivenData(fluidIndicatorPtr_);
    deleteDemandDrivenData(rhoPtr_);
    deleteDemandDrivenData(muPtr_);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return constant reference to the mesh
dynamicFvMesh& Foam::freeSurface::mesh() const
{
    return mesh_;
}

// Return the interface patchID
label Foam::freeSurface::aPatchID() const
{
    return aPatchID_;
}

// Return the shadow patchID
label Foam::freeSurface::bPatchID() const
{
    return bPatchID_;
}

void Foam::freeSurface::initializeControlPointsPosition() const
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

scalar Foam::freeSurface::maxCourantNumber()
{
    scalar CoNum = 0;

    if (cleanInterface())
    {
        const scalarField& dE = aMesh().lPN();

        CoNum = max
        (
            mesh().time().deltaT().value()/
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
            mesh().time().deltaT().value()/
            sqrt
            (
                rhoFluidA().value()*dE*dE*dE/
                2.0/M_PI/sigmaE
            )
        );
    }

    return CoNum;
}

//- Update control points end displacement directions
void Foam::freeSurface::updateDisplacementDirections() const
{
    // Update point displacement correction
    pointsDisplacementDir() = aMesh().pointAreaNormals();

    // Update face displacement direction
    facesDisplacementDir() = aMesh().faceAreaNormals().internalField();

    // Correction of control points position
    areaCentrePositions() = aMesh().areaCentres().internalField();

    controlPoints() =
        areaCentrePositions()
      + facesDisplacementDir()
      * (facesDisplacementDir()&(controlPoints() - areaCentrePositions()));
}

//- Move control points by deltaH and calculate interface
//  points displacement for the new control points position
tmp<vectorField> Foam::freeSurface::pointDisplacement
(
    const scalarField& deltaH
) const
{
    const pointField& points = aMesh().patch().localPoints();
    const labelListList& pointFaces = aMesh().patch().pointFaces();

    // Update control points
    controlPoints() += facesDisplacementDir()*deltaH;

    tmp<vectorField> tdisplacement
    (
        new vectorField
        (
            points.size(),
            vector::zero
        )
    );

    vectorField& displacement = tdisplacement();

    // Calculate displacement of internal points

    const edgeList& edges = aMesh().patch().edges();

    labelList internalPoints = aMesh().internalPoints();

    forAll (internalPoints, pointI)
    {
        label curPoint = internalPoints[pointI];

        const labelList& curPointFaces = pointFaces[curPoint];

        scalarField w(curPointFaces.size(), 0.0);

        forAll (curPointFaces, faceI)
        {
            label curFace = curPointFaces[faceI];

            scalar magDistance =
            (
                 stabilise
                 (
                     mag
                     (
                        controlPoints()[curFace]
                      - points[curPoint]
                     ), VSMALL
                 )
            );

            w[faceI] = 1.0/magDistance;
        }

        w /= sum(w);

        vector Q = vector::zero;

        forAll (curPointFaces, faceI)
        {
            label curFace = curPointFaces[faceI];

            Q += w[faceI]*controlPoints()[curFace];
        }

        displacement[curPoint] =
        (
           -1.0*
           (
               pointsDisplacementDir()[curPoint] & (points[curPoint] - Q)
              *pointsDisplacementDir()[curPoint]
           )
        );
    }

    // Calculate displacement of points
    // which belonge to empty and wedge patches

    labelList boundaryPoints = aMesh().boundaryPoints();

    const labelListList& edgeFaces = aMesh().patch().edgeFaces();
    const labelListList& pointEdges = aMesh().patch().pointEdges();

    vectorField pointNormals = aMesh().pointAreaNormals();

    forAll (boundaryPoints, pointI)
    {
        label curPoint = boundaryPoints[pointI];

        if (motionPointsMask()[curPoint])
        {
            // Calculating mirror points
            const labelList& curPointEdges = pointEdges[curPoint];

            vectorField mirrorPoints(2, vector::zero);

            label counter = -1;

            forAll (curPointEdges, edgeI)
            {
                label curEdge = curPointEdges[edgeI];

                if(edgeFaces[curEdge].size() == 1)
                {
                    label patchID = -1;

                    forAll(aMesh().boundary(), patchI)
                    {
                        forAll(aMesh().boundary()[patchI], eI)
                        {
                            if(aMesh().boundary()[patchI][eI] == curEdge)
                            {
                                patchID = patchI;
                                break;
                            }
                        }
                    }

                    if
                    (
                        patchID > -1
                     && (
                            aMesh().boundary()[patchID].type()
                         == "wedge"
                        )
                    )
                    {
                        const wedgeFaPatch& wedgePatch =
                            refCast<const wedgeFaPatch>
                            (
                                aMesh().boundary()[patchID]
                            );

                        mirrorPoints[++counter] =
                            transform
                            (
                                wedgePatch.faceT(),
                                controlPoints()[edgeFaces[curEdge][0]]
                            );
                    }
                    else
                    {
                        vector nE =
                            pointNormals[edges[curEdge].start()]
                          + pointNormals[edges[curEdge].end()];

                        nE /= mag(nE);

                        vector eP =
                            controlPoints()[edgeFaces[curEdge][0]]
                          - edges[curEdge].centre(points);

                        mirrorPoints[++counter] =
                            edges[curEdge].centre(points)
                          + ((2.0*nE*nE - I)&eP);
                    }
                }
            }


            // Calculating LS plane interpolation
            const labelList& curPointFaces = pointFaces[curPoint];

            symmTensor M = symmTensor::zero;

            vector S = vector::zero;

            scalarField w(curPointFaces.size() + 2, 0.0);

            forAll (curPointFaces, faceI)
            {
                label curFace = curPointFaces[faceI];

                w[faceI] = 1.0/mag
                    (
                        controlPoints()[curFace]
                      - points[curPoint]
                    );
            }

            forAll (mirrorPoints, pI)
            {
                w[curPointFaces.size() + pI] = 1.0/mag
                    (
                        mirrorPoints[pI]
                      - points[curPoint]
                    );
            }

            w /= sum(w);


            forAll (curPointFaces, faceI)
            {
                label curFace = curPointFaces[faceI];

                M = M + sqr(w[faceI])*sqr(controlPoints()[curFace]);

                S += sqr(w[faceI])*controlPoints()[curFace];
            }

            forAll (mirrorPoints, pI)
            {
                M = M + sqr(w[curPointFaces.size()+pI])*sqr(mirrorPoints[pI]);

                S += sqr(w[curPointFaces.size()+pI])*mirrorPoints[pI];
            }


            vector N = inv(M)&S;

            N /= mag(N);

            scalar p = (S&N)/sum(sqr(w));

            displacement[curPoint] =
                pointsDisplacementDir()[curPoint]*
                (p - (points[curPoint]&N))/
                (pointsDisplacementDir()[curPoint]&N);
        }
    }


    forAll(aMesh().boundary(), patchI)
    {
        bool fixedPatch = false;

        forAll(fixedFreeSurfacePatches_, fpI)
        {
            label fixedPatchID = aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[fpI]
            );

            if (fixedPatchID == -1)
            {
                FatalErrorIn
                (
                    "freeSurface::pointDisplacement(...)"
                )
                    << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
            }

            if (fixedPatchID == patchI)
            {
                fixedPatch = true;
            }
        }

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
         && !fixedPatch
        )
        {
            labelList patchPoints =
                aMesh().boundary()[patchI].pointLabels();

            labelListList patchPointEdges =
                aMesh().boundary()[patchI].pointEdges();

            unallocLabelList patchEdgeFaces =
                aMesh().boundary()[patchI].edgeFaces();

            forAll(patchPoints, pointI)
            {
                forAll(patchPointEdges[pointI], edgeI)
                {
                    label curEdge = patchPointEdges[pointI][edgeI];

                    displacement[patchPoints[pointI]] +=
                        pointsDisplacementDir()[patchPoints[pointI]]
                       *(
                            pointsDisplacementDir()[patchPoints[pointI]]
                           &(
                                controlPoints()[patchEdgeFaces[curEdge]]
                              - points[patchPoints[pointI]]
                            )
                        );
                }

                displacement[patchPoints[pointI]] /=
                    patchPointEdges[pointI].size();
            }
        }
    }


    // Calculate displacement of axis point
    forAll (aMesh().boundary(), patchI)
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
                label axisPoint = wedgePatch.axisPoint();

                displacement[axisPoint] =
                    pointsDisplacementDir()[axisPoint]
                   *(
                        pointsDisplacementDir()[axisPoint]
                       &(
                            controlPoints()[pointFaces[axisPoint][0]]
                          - points[axisPoint]
                        )
                    );
            }
        }
    }

    return tdisplacement;
}

//- Return the face curvature field
const scalarField& Foam::freeSurface::faceCurvatures()
{
    return aMesh().faceCurvatures().internalField();
}

//- Adjust the surface-tension for temperature
void Foam::freeSurface::adjustSurfaceTension(const volScalarField& T)
{
    const labelList& fCells = mesh().boundaryMesh()[aPatchID()].faceCells();

    scalar adjustedValue = 0.0;

    // Set the surface tension field according to temperature
    forAll(fCells, faceI)
    {
        adjustedValue =
            cleanInterfaceSurfTension().value()
            - (0.00018*(T.internalField()[fCells[faceI]] - 298));

        surfaceTension().internalField()[faceI] =
        (
            adjustedValue < 0 ?
            cleanInterfaceSurfTension().value() : adjustedValue
        );
    }

    Info << "Surface Tension: "
         << " Min: " << min(surfaceTension().internalField())
         << " Max: " << max(surfaceTension().internalField())
         << " Average: " << average(surfaceTension().internalField())
         << endl;

    if (min(surfaceTension().internalField()) < 0)
    {
        FatalErrorIn("freeSurface::adjustSurfaceTension()")
            << "Negative surface-tension!"
            << abort(FatalError);
    }
}

// Update the interface with the fluid velocity
void Foam::freeSurface::moveSurfacePoints()
{
    pointField newMeshPoints = mesh().points();

    scalarField sweptVolCorr =
    (
        phi().boundaryField()[aPatchID()]
      - fvc::meshPhi(rho(),U())().boundaryField()[aPatchID()]
    );

    word ddtScheme
    (
        mesh().ddtScheme("ddt(" + rho().name() + ',' + U().name()+')')
    );

    if (ddtScheme == "CrankNicholson")
    {
        sweptVolCorr *= (1.0/2.0)*mesh().time().deltaT().value();
    }
    else
    if (ddtScheme == "Euler")
    {
        sweptVolCorr *= mesh().time().deltaT().value();
    }
    else
    if (ddtScheme == "backward")
    {
        sweptVolCorr *= (2.0/3.0)*mesh().time().deltaT().value();
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
    const labelList& meshPtsA =
    (
        mesh().boundaryMesh()[aPatchID()].meshPoints()
    );

    forAll (displacement, pointI)
    {
        newMeshPoints[meshPtsA[pointI]] += displacement[pointI];
    }

    if (twoFluids_)
    {
        const labelList& meshPtsB =
        (
            mesh().boundaryMesh()[bPatchID_].meshPoints()
        );

        pointField displacementB =
        (
            interpolatorAB().pointInterpolate(displacement)
        );

        forAll (displacementB, pointI)
        {
            newMeshPoints[meshPtsB[pointI]] += displacementB[pointI];
        }
    }

    // Update total displacement field
    if (totalDisplacementPtr_ && (curTimeIndex_ < mesh().time().timeIndex()))
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Total displacement of free surface points "
            << "from previous time step is not absorbed by the mesh."
            << abort(FatalError);
    }
    else
    if (curTimeIndex_ < mesh().time().timeIndex())
    {
        totalDisplacement() = displacement;

        curTimeIndex_ = mesh().time().timeIndex();
    }
    else
    {
        totalDisplacement() += displacement;
    }

    // Move interface points
    mesh().movePoints(newMeshPoints);
    aMesh().movePoints(newMeshPoints);
    moveFvSubMeshes();
}

// Restore interface position
bool Foam::freeSurface::restorePosition()
{
    if (!totalDisplacementPtr_)
    {
        return false;
    }

    pointField newMeshPoints = mesh().points();

    const labelList& meshPtsA = mesh().boundaryMesh()[aPatchID()].meshPoints();

    forAll (meshPtsA, pointI)
    {
        newMeshPoints[meshPtsA[pointI]] -= totalDisplacement()[pointI];
    }

    if (twoFluids())
    {
        const labelList& meshPtsB =
        (
            mesh().boundaryMesh()[bPatchID()].meshPoints()
        );

        // Interpolate displacement to the shadow patch
        pointField dispB =
        (
            interpolatorAB().pointInterpolate(totalDisplacement())
        );

        forAll(meshPtsB,pointI)
        {
            newMeshPoints[meshPtsB[pointI]] -= dispB[pointI];
        }
    }

    mesh().movePoints(newMeshPoints);

    // Now set boundary conditions for the motion-solver
    setMotionBC(mesh(), aPatchID(), totalDisplacement());

    if (twoFluids())
    {
        // Interpolate displacement to the shadow patch
        pointField dispB =
        (
            interpolatorAB().pointInterpolate(totalDisplacement())
        );

        setMotionBC(mesh(), bPatchID(), dispB);
    }

    deleteDemandDrivenData(totalDisplacementPtr_);

    // Now update the mesh with motion / topology changes.
    bool meshChanged = mesh().update();

    aMesh().movePoints(mesh().points());
    moveFvSubMeshes();

    return meshChanged;
}

void Foam::freeSurface::correctUsBoundaryConditions() const
{
    forAll(Us().boundaryField(), patchI)
    {
        if
        (
            Us().boundaryField()[patchI].type()
         == "calculated"
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
                     == "slip"
                    )
                 ||
                    (
                        U().boundaryField()[ngbPolyPatchID].type()
                     == "symmetryPlane"
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

//- Move correctedFvPatchField fvSubMeshes
void Foam::freeSurface::moveFvSubMeshes()
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
                    U().boundaryField()[patchI]
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
                    p().boundaryField()[patchI]
                )
            );

            pP.movePatchSubMesh();
        }
    }
}

// Update boundary conditions on velocity and pressure
void Foam::freeSurface::updateBoundaryConditions()
{
    // Update boundary conditions
    updateVelocity();
    updatePressure();
}

//- Update pressure boundary conditions
void Foam::freeSurface::updatePressure()
{
    // Correct pressure boundary condition at the free-surface
    vectorField nA = mesh().boundary()[aPatchID()].nf();

    if (twoFluids())
    {
        scalarField pA =
        (
            interpolatorBA().faceInterpolate
            (
                p().boundaryField()[bPatchID()]
            )
        );

        const scalarField& K = aMesh().faceCurvatures().internalField();

        Info << "Free surface curvature: min = " << min(K)
            << ", max = " << max(K)
            << ", average = " << average(K) << endl << flush;

        if (cleanInterface())
        {
            pA -= cleanInterfaceSurfTension().value()*(K - average(K));
        }
        else
        {
            scalarField surfTensionK =
            (
                surfaceTension().internalField()*K
            );

            pA -= surfTensionK - average(surfTensionK);
        }

        scalarField nnGradU = (nA & U().boundaryField()[aPatchID()].snGrad());

        pA += 2.0*(muFluidA().value() - muFluidB().value())*nnGradU;

        vector R0 = vector(0, 0, 0);

        pA -=
        (
            (rhoFluidA().value() - rhoFluidB().value())*
            (
                g_.value()
              & (
                    mesh().C().boundaryField()[aPatchID()]
                  - R0
                )
            )
        );

        p().boundaryField()[aPatchID()] == pA;
    }
    else
    {
        vector R0 = vector(0, 0, 0);

        scalarField pA =
        (
          - rhoFluidA().value()*
            (
                g_.value()
              & (
                    mesh().C().boundaryField()[aPatchID()]
                  - R0
                )
            )
        );

        const scalarField& K = aMesh().faceCurvatures().internalField();

        Info << "Free surface curvature: min = " << min(K)
            << ", max = " << max(K) << ", average = " << average(K)
            << endl;

        if (cleanInterface())
        {
            pA -= cleanInterfaceSurfTension().value()*(K - average(K));
        }
        else
        {
            scalarField surfTensionK =
            (
                surfaceTension().internalField()*K
            );

            pA -= surfTensionK - average(surfTensionK);
        }

        scalarField nnGradU = (nA & U().boundaryField()[aPatchID()].snGrad());

        pA += 2.0*muFluidA().value()*nnGradU;

        p().boundaryField()[aPatchID()] == pA;
    }

    // Set modified pressure at patches with fixed absolute pressure
    vector R0 = vector(0, 0, 0);

    forAll(p().boundaryField(), patchI)
    {
        if (p().boundaryField()[patchI].type() == "fixedValue")
        {
            if (patchI != aPatchID())
            {
                p().boundaryField()[patchI] ==
                (
                    p().boundaryField()[patchI]
                  - rho().boundaryField()[patchI]*
                    (g_.value()&(mesh().C().boundaryField()[patchI] - R0))
                );
            }
        }
    }
}

//- Update velocity boundary conditions
void Foam::freeSurface::updateVelocity()
{
    if (twoFluids_)
    {
        vectorField nA = mesh().boundary()[aPatchID()].nf();
        vectorField nB = mesh().boundary()[bPatchID()].nf();

        // Update fixedValue boundary condition on patch B
        scalarField DnB = interpolatorBA().faceInterpolate
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

        vectorField UtFs =
        (
            muFluidA().value()*DnA*UtPA
          + muFluidB().value()*DnB*UtPB
        );

        vectorField UnFs = nA*phi().boundaryField()[aPatchID()]
            /mesh().boundary()[aPatchID()].magSf();

        Us().internalField() += UnFs - nA*(nA&Us().internalField());
        correctUsBoundaryConditions();

        UtFs -= (muFluidA().value() - muFluidB().value())*
            (fac::grad(Us()) & aMesh().faceAreaNormals())().internalField();

        if (!cleanInterface())
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

        scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();

        vectorField UtPA = U().boundaryField()[aPatchID()].patchInternalField();

        // Obtain the tangential component of velocity
        UtPA -= nA*(nA & UtPA);

        // Tangential force
        vectorField UtFs = muFluidA().value()*DnA*UtPA;

        // Add effects of surface-tension gradient
        if (!cleanInterface())
        {
            UtFs += surfaceTensionGrad()().internalField();
        }

        // Normal
        vectorField UnFs =
            nA*phi().boundaryField()[aPatchID()]
          / mesh().boundary()[aPatchID()].magSf();

        Us().internalField() += UnFs - nA*(nA&Us().internalField());
        correctUsBoundaryConditions();

        UtFs -= muFluidA().value()*
            (fac::grad(Us())&aMesh().faceAreaNormals())().internalField();

        UtFs /= muFluidA().value()*DnA + VSMALL;

        // Update free-surface velocity
        Us().internalField() = UtFs + UnFs;
        correctUsBoundaryConditions();

        // Update fixedGradient boundary condition on the free-surface
        vectorField nGradU =
            -nA*muFluidA().value()*fac::div(Us())().internalField();

        nGradU += muFluidA().value()*
            (UtFs - UtPA)*mesh().boundary()[aPatchID()].deltaCoeffs();

        nGradU /= muFluidA().value() + VSMALL;

        if
        (
            U().boundaryField()[aPatchID()].type()
         == "fixedGradientCorrected"
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
         == "fixedGradient"
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
                << " is " << U().boundaryField()[aPatchID()].type()
                << ", instead of "
                << fixedGradientFvPatchField<vector>::typeName
                << " or "
                << fixedGradientCorrectedFvPatchField<vector>::typeName
                << abort(FatalError);
        }
    }
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

    // Copy old data
    /*
    vectorField oldXf(areaCentrePositions());
    vectorField oldCp(controlPoints());
    */

    // Wipe out demand-driven data
    clearOut();

    /*
    updateDisplacementDirections();

    vectorField& Cp = controlPoints();

    const labelList& fMap = mpm.faceMap();
    const vectorField& XfNew = areaCentrePositions();

    // Map old data to new field sizes
    vectorField XfOld(XfNew.size(), vector::zero);
    vectorField CpOld(XfNew.size(), vector::zero);

    label opStart = mpm.oldPatchStarts()[aPatchID()];
    label npStart = mpm.mesh().boundaryMesh()[aPatchID()].start();

    forAll(XfNew, pI)
    {
        label newIndex = (npStart + pI);
        label oldIndex = fMap[newIndex];

        if (oldIndex > -1)
        {
            XfOld[pI] = oldXf[oldIndex - opStart];
            CpOld[pI] = oldCp[oldIndex - opStart];
        }
        else
        {
            FatalErrorIn("freeSurface::updateMesh()")
                << nl << " Failed to determine mapping face. "
                << abort(FatalError);
        }
    }

    vectorField corrVec(Cp.size(), vector::zero);

    forAll(Cp, pI)
    {
        vector oCv = (CpOld[pI] - XfOld[pI]);
        vector nCv = (CpOld[pI] - XfNew[pI]);

        oCv /= (mag(oCv) + VSMALL);

        vector corr = (nCv - (nCv & oCv)*oCv);

        Cp[pI] = ((XfNew[pI] + nCv) - corr);

        corrVec[pI] = corr;
    }

    if (debug)
    {
        Info << "corrVec: " << endl;
        Info << corrVec << endl;

        Info << "oldCp: " << endl;
        Info << oldCp << endl;

        Info << "Cp: " << endl;
        Info << Cp << endl;
    }
    */

    return true;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::freeSurface::operator=(const freeSurface& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("freeSurface::operator=(const freeSurface&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}


// ************************************************************************* //
