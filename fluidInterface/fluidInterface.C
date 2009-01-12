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
    Sandeep Menon

\*----------------------------------------------------------------------------*/

#include "fluidInterface.H"

#include "wedgeFaPatch.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidInterface::fluidInterface
(
    fvMesh& m,
    volVectorField& U,
    volScalarField& p,
    surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "fluidInterfaceProperties",
            U.mesh().time().constant(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(m),
    U_(U),
    p_(p),
    phi_(phi),
    interfacePatch_(word(lookup("interfacePatch"))),
    patchID_(-1),
    curTimeIndex_(U.mesh().time().timeIndex()),
    muFluidA_(lookup("muFluidA")),
    rhoFluidA_(lookup("rhoFluidA")),
    condFluidA_(lookup("condFluidA")),
    CpFluidA_(lookup("CpFluidA")),
    surfaceTension_(lookup("surfaceTension")),
    displacementPtr_(NULL),
    motionPointsMaskPtr_(NULL),
    controlPointsPtr_(NULL),            
    pointsDisplacementDirPtr_(NULL),
    facesDisplacementDirPtr_(NULL),
    surfaceTensionPtr_(NULL),
    aMeshPtr_(NULL),
    UsPtr_(NULL),
    fixedFreeSurfacePatches_
    (
        lookup("fixedFreeSurfacePatches")
    )
{
    // Loop through all boundary patches and determine the patchID
    const fvBoundaryMesh& bdy = mesh().boundary();
    forAll (bdy, patchI)
    {
        if(bdy[patchI].name() == interfacePatch_)
        {
            patchID_ = patchI;
        }
    }

    if(patchID_ == -1)
    {
        FatalErrorIn("fluidInterface::fluidInterface()")
            << "Patch with name " << interfacePatch_ << " is not defined."
            << abort(FatalError);
    }
    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidInterface::~fluidInterface()
{}

void Foam::fluidInterface::clearOut()
{
    deleteDemandDrivenData(displacementPtr_);
    deleteDemandDrivenData(controlPointsPtr_);
    deleteDemandDrivenData(motionPointsMaskPtr_);
    deleteDemandDrivenData(pointsDisplacementDirPtr_);
    deleteDemandDrivenData(facesDisplacementDirPtr_);
    deleteDemandDrivenData(aMeshPtr_);
    deleteDemandDrivenData(surfaceTensionPtr_);
    deleteDemandDrivenData(UsPtr_);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return constant reference to the mesh
fvMesh& Foam::fluidInterface::mesh()
{
    return mesh_;
}

// Return constant reference to the area mesh
faMesh& Foam::fluidInterface::areaMesh()
{
    if (!aMeshPtr_)
    {
        makeFaMesh();
    }

    return *aMeshPtr_;     
}

// Return the interface patchID
const label Foam::fluidInterface::patchID()
{
    return patchID_;
}

// Make the displacement field
void Foam::fluidInterface::makeDisplacement()
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (displacementPtr_)
    {
        FatalErrorIn("fluidInterface::makeDisplacement()")
            << "Displacement field already exists."
            << abort(FatalError);
    }    
    
    displacementPtr_ =
        new vectorField
        (
            mesh().boundaryMesh()[patchID()].nPoints(),
            vector::zero
        );
}

// Make the finite-area mesh
void Foam::fluidInterface::makeFaMesh()
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
    
    // Lookup for point-normals correction
    wordList pointNormalsCorrectionPatches
    (
        lookup("pointNormalsCorrectionPatches")
    );
    
    // Set point normal correction patches
    boolList& correction = aMeshPtr_->correctPatchPointNormals();

    forAll(pointNormalsCorrectionPatches, patchI)
    {
        word patchName = pointNormalsCorrectionPatches[patchI];

        label patchID = aMeshPtr_->boundary().findPatchID(patchName);

        if(patchID == -1)
        {
            FatalErrorIn
            (
                "fluidInterface::makeFaMesh()"
            )   << "patch name for point normals correction don't exist"
                << abort(FatalError);
        }

        correction[patchID] = true;
    }    
}

void Foam::fluidInterface::makeSurfaceTension()
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
        areaMesh(),
        surfaceTension_
    );
}

void Foam::fluidInterface::makeUs()
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (UsPtr_)
    {
        FatalErrorIn("Foam::fluidInterface::makeUs()")
            << "Free-surface velocity field already exists"
            << abort(FatalError);
    }

    wordList patchFieldTypes
    (
        areaMesh().boundary().size(),
        calculatedFaPatchVectorField::typeName
    );

    forAll(areaMesh().boundary(), patchI)
    {
        if
        (
            areaMesh().boundary()[patchI].type()
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
        areaMesh(),
        dimensioned<vector>("Us", dimVelocity, vector::zero),
        patchFieldTypes
    );
}

void Foam::fluidInterface::makeMotionPointsMask()
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (motionPointsMaskPtr_)
    {
        FatalErrorIn("Foam::fluidInterface::motionPointsMask()")
            << "Motion points mask already exists"
            << abort(FatalError);
    }

    if(patchID() == -1)
    {
        FatalErrorIn("Foam::fluidInterface::makeMotionPointsMask()")
            << "Interface patch is not defined."
            << abort(FatalError);
    }

    motionPointsMaskPtr_ = new scalarField
    (
        mesh().boundaryMesh()[patchID()].nPoints(),
        1.0
    );
    
    // Mark free surface boundary points
    // which do not belong to empty and wedge patches
    forAll(areaMesh().boundary(), patchI)
    {
        if
        (
            (
                areaMesh().boundary()[patchI].type()
             != "empty"
            )
         && (
                areaMesh().boundary()[patchI].type()
             != "wedge"
            )
        )
        {
            labelList patchPoints =
                areaMesh().boundary()[patchI].pointLabels();

            forAll(patchPoints, pointI)
            {
                motionPointsMask()[patchPoints[pointI]] = 0.0;
            }
        }
    }  
    
    // Mark free-surface boundary point
    // at the axis of 2-D axisymmetic cases
    forAll(areaMesh().boundary(), patchI)
    {
        if
        (
            areaMesh().boundary()[patchI].type()
         == "wedge"
        )
        {
            const wedgeFaPatch& wedgePatch =
                refCast<const wedgeFaPatch>(areaMesh().boundary()[patchI]);

            if(wedgePatch.axisPoint() > -1)
            {
                motionPointsMask()[wedgePatch.axisPoint()] = 0;
            }
        }
    }    
}

void Foam::fluidInterface::makeControlPoints()
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (controlPointsPtr_)
    {
        FatalErrorIn("Foam::fluidInterface::makeControlPoints()")
            << "Control points already exist"
            << abort(FatalError);
    }

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
            areaMesh().centres().internalField()
        );

    initializeControlPointsPosition();
}

void Foam::fluidInterface::initializeControlPointsPosition()
{
    scalarField deltaH = scalarField(areaMesh().nFaces(), 0.0);

    pointField displacement = pointDisplacement(deltaH);

    const faceList& faces = areaMesh().faces();
    const pointField& points = areaMesh().points();

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


void Foam::fluidInterface::makeDirections()
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if
    (
        pointsDisplacementDirPtr_ ||
        facesDisplacementDirPtr_
    )
    {
        FatalErrorIn("Foam::fluidInterface::makeDirections()")
            << "Point and face displacement directions "
            << "already exists"
            << abort(FatalError);
    }
    
    if(patchID() == -1)
    {
        FatalErrorIn("freeSurface::makeDirections()")
            << "Interface patch is not defined."
            << abort(FatalError);
    }

    pointsDisplacementDirPtr_ =
        new vectorField
        (
            mesh().boundaryMesh()[patchID()].nPoints(),
            vector::zero
        );

    facesDisplacementDirPtr_ =
        new vectorField
        (
            mesh().boundaryMesh()[patchID()].size(),
            vector::zero
        );
    
    updateDisplacementDirections();    
}

//- Update control points end displacement directions
void Foam::fluidInterface::updateDisplacementDirections()
{
    // Update point displacement correction
    pointsDisplacementDir() = 
        areaMesh().pointAreaNormals();    
    
    // Update face displacement direction
    facesDisplacementDir() =
        areaMesh().faceAreaNormals().internalField();

    // Correction of control points postion
    const vectorField& Cf = 
        areaMesh().centres().internalField();
    
    controlPoints() =
        facesDisplacementDir()
       *(facesDisplacementDir()&(controlPoints() - Cf))
      + Cf;
}

//- Move control points by deltaH and calculate interface 
//  points displacement for the new control points position
tmp<vectorField> Foam::fluidInterface::pointDisplacement(const scalarField& deltaH)
{
    const pointField& points = areaMesh().patch().localPoints();
    const labelListList& pointFaces = areaMesh().patch().pointFaces();
    
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

    const edgeList& edges = areaMesh().patch().edges();

    labelList internalPoints = areaMesh().internalPoints();

    forAll (internalPoints, pointI)
    {
        label curPoint = internalPoints[pointI];

        const labelList& curPointFaces = pointFaces[curPoint];

        symmTensor M = symmTensor::zero;

        vector S = vector::zero;


        scalarField w(curPointFaces.size(), 0.0);

        forAll (curPointFaces, faceI)
        {
            label curFace = curPointFaces[faceI];

            w[faceI] = 1.0/mag
            (
                controlPoints()[curFace]
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

        vector N = inv(M)&S;

        N /= mag(N);

        scalar p = (S&N)/sum(sqr(w));

        displacement[curPoint] = 
            pointsDisplacementDir()[curPoint]*
            (p - (points[curPoint]&N))/
            (pointsDisplacementDir()[curPoint]&N);
    }



    // Calculate displacement of points 
    // which belonge to empty and wedge patches

    labelList boundaryPoints = areaMesh().boundaryPoints();

    const labelListList& edgeFaces = areaMesh().patch().edgeFaces();
    const labelListList& pointEdges = areaMesh().patch().pointEdges();

    vectorField pointNormals = areaMesh().pointAreaNormals();
            
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

                    forAll(areaMesh().boundary(), patchI)
                    {
                        forAll(areaMesh().boundary()[patchI], eI)
                        {
                            if(areaMesh().boundary()[patchI][eI] == curEdge)
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
                            areaMesh().boundary()[patchID].type()
                         == "wedge"
                        )
                    )
                    {
                        const wedgeFaPatch& wedgePatch =
                            refCast<const wedgeFaPatch>
                            (
                                areaMesh().boundary()[patchID]
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


    forAll(areaMesh().boundary(), patchI)
    {
        bool fixedPatch = false;

        forAll(fixedFreeSurfacePatches_, fpI)
        {
            label fixedPatchID = areaMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[fpI]
            );

            if(fixedPatchID == -1)
            {
                FatalErrorIn("fluidInterface::pointDisplacement(...)")
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
                areaMesh().boundary()[patchI].type()
             != "empty"
            )
         && ( 
                areaMesh().boundary()[patchI].type()
             != "wedge"
            )            
         && !fixedPatch
        )
        {
            labelList patchPoints =
                areaMesh().boundary()[patchI].pointLabels();

            labelListList patchPointEdges = 
                areaMesh().boundary()[patchI].pointEdges();

            unallocLabelList patchEdgeFaces = 
                areaMesh().boundary()[patchI].edgeFaces();

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
    forAll (areaMesh().boundary(), patchI)
    {
        if
        (
            areaMesh().boundary()[patchI].type() 
         == "wedge"
        )
        {
            const wedgeFaPatch& wedgePatch =
                refCast<const wedgeFaPatch>(areaMesh().boundary()[patchI]);

            if(wedgePatch.axisPoint() > -1)
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

// Return reference to the displacement field
vectorField& Foam::fluidInterface::displacement()
{
    if (!displacementPtr_)
    {
        makeDisplacement();
    }

    return *displacementPtr_;    
}

//- Return reference to the control points
vectorField& Foam::fluidInterface::controlPoints()
{
    if (!controlPointsPtr_)
    {
        makeControlPoints();
    }

    return *controlPointsPtr_;    
}

// Return reference to point displacement direction field
vectorField& Foam::fluidInterface::pointsDisplacementDir()
{
    if (!pointsDisplacementDirPtr_)
    {
        makeDirections();
    }

    return *pointsDisplacementDirPtr_;    
}

// Return reference to control points displacement direction field
vectorField& Foam::fluidInterface::facesDisplacementDir()
{
    if (!facesDisplacementDirPtr_)
    {
        makeDirections();
    }

    return *facesDisplacementDirPtr_;    
}

scalarField& Foam::fluidInterface::motionPointsMask()
{
    if (!motionPointsMaskPtr_)
    {
        makeMotionPointsMask();
    }

    return *motionPointsMaskPtr_;
}

areaScalarField& Foam::fluidInterface::surfaceTension()
{
    if (!surfaceTensionPtr_)
    {
        makeSurfaceTension();
    }
    
    return *surfaceTensionPtr_;
}

areaVectorField& Foam::fluidInterface::Us()
{
    if (!UsPtr_)
    {
        makeUs();
    }

    return *UsPtr_;
}

//- Adjust the surface-tension for temperature
void Foam::fluidInterface::adjustSurfaceTension(const volScalarField& T)
{
    const labelList& fCells = mesh().boundaryMesh()[patchID()].faceCells();

    scalar adjustedValue = 0.0;

    // Set the surface tension field according to temperature
    forAll(fCells, faceI)
    {
        adjustedValue =
            surfaceTension_.value()
            - (0.00018*(T.internalField()[fCells[faceI]] - 298));

        surfaceTension().internalField()[faceI] =
            adjustedValue < 0 ? surfaceTension_.value() : adjustedValue;
    }
    
    Info << "Surface Tension: " 
         << " Min: " << min(surfaceTension().internalField())
         << " Max: " << max(surfaceTension().internalField())
         << " Average: " << average(surfaceTension().internalField())
         << endl;

    if (min(surfaceTension().internalField()) < 0)
    {
        FatalErrorIn("fluidInterface::adjustSurfaceTension()")
            << "Negative surface-tension!"
            << abort(FatalError);
    }
}

// Update the interface with the fluid velocity
void Foam::fluidInterface::movePoints()
{
    // Obtain the current time-step
    scalar dt = mesh().time().deltaT().value();
    
    // Swept-volume correction for the explicit-Euler scheme
    scalarField sweptVolCorr = 
    (
        phi_.boundaryField()[patchID()]
      - fvc::meshPhi(U_)().boundaryField()[patchID()]
    )*dt;

    const scalarField& Sf = areaMesh().S();
    const vectorField& Nf = areaMesh().faceAreaNormals().internalField();

    scalarField deltaH = sweptVolCorr/(Sf*(Nf & facesDisplacementDir()));    
    
    // Obtain the interface displacement vectors
    vectorField disp = pointDisplacement(deltaH);
    
    if (curTimeIndex_ < mesh().time().timeIndex())
    {
        displacement() = disp;

        curTimeIndex_ = mesh().time().timeIndex();
    }
    else
    {
        displacement() += disp;
    }    
    
    // Move the areaMesh to recalculate surface curvature
    pointField newMeshPoints = mesh().points();
    
    const labelList& meshPts = mesh().boundaryMesh()[patchID()].meshPoints();
    
    forAll(meshPts,pointI)
    {
        newMeshPoints[meshPts[pointI]] += disp[pointI];
    }
    
    mesh().movePoints(newMeshPoints);
    
    areaMesh().movePoints(newMeshPoints);
}

// Restore interface position
void Foam::fluidInterface::restorePosition()
{
    vectorField& disp = displacement();
    
    pointField newMeshPoints = mesh().points();
    
    const labelList& meshPts = mesh().boundaryMesh()[patchID()].meshPoints();
    
    forAll (meshPts, pointI)
    {
        newMeshPoints[meshPts[pointI]] -= disp[pointI];
    }

    mesh().movePoints(newMeshPoints);

    areaMesh().movePoints(newMeshPoints);    
}

void Foam::fluidInterface::correctUsBoundaryConditions()
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
                areaMesh().boundary()[patchI].ngbPolyPatchIndex();

            if(ngbPolyPatchID != -1)
            {
                if
                (
                    (
                        U_.boundaryField()[ngbPolyPatchID].type()
                     == "slip"
                    )
                 ||
                    (
                        U_.boundaryField()[ngbPolyPatchID].type()
                     == "symmetryPlane"
                    )
                )
                {
                    vectorField N =
                        areaMesh().boundary()[patchI].ngbPolyPatchFaceNormals();

                    pUs -= N*(N&pUs);
                }
            }
        }
    }

    Us().correctBoundaryConditions();
}

// Update boundary conditions on velocity and pressure
void Foam::fluidInterface::updateBoundaryConditions()
{
    updateVelocity();
    updatePressure();
}

//- Update pressure boundary conditions
void Foam::fluidInterface::updatePressure()
{
    const scalarField& K = areaMesh().faceCurvatures().internalField();

    scalarField pA(K.size(),0.0);    
    
    Info << "Free surface curvature: min = " << min(K)
        << ", max = " << max(K) << ", average = " << average(K)
        << endl;

    scalarField surfTensionK = surfaceTension().internalField()*K;

    dimensionSet dimSurfaceTension = (dimForce/dimLength);
    
    // Check the dimensions of pressure and 
    // adjust surface-tension accordingly        
    if (p_.dimensions() == (dimPressure/dimDensity))
    {
        if (surfaceTension_.dimensions() == dimSurfaceTension)
        {
            surfTensionK /= rhoFluidA().value();
        }
        else 
        if 
        (
            surfaceTension_.dimensions() != (dimSurfaceTension/dimDensity)
        )
        {
            FatalErrorIn("fluidInterface::updatePressure()")
                << "Dimensions for surface tension are inconsistent."
                << abort(FatalError);            
        }
    }
    
    if 
    (
         (p_.dimensions() == dimPressure)
      && (surfaceTension_.dimensions() != dimSurfaceTension)
    )
    {
        FatalErrorIn("fluidInterface::updatePressure()")
            << "Dimensions for surface tension are inconsistent."
            << abort(FatalError);         
    }
    
    pA -= (surfTensionK - average(surfTensionK));
        
    p_.boundaryField()[patchID()] == pA;
}

//- Update velocity boundary conditions
void Foam::fluidInterface::updateVelocity()
{
    vectorField nA = areaMesh().faceAreaNormals().internalField();

    scalarField DnA = mesh().boundary()[patchID()].deltaCoeffs();

    vectorField UtPA = U_.boundaryField()[patchID()].patchInternalField();

    // Obtain the tangential component of velocity
    UtPA -= nA*(nA & UtPA);

    // Tangential force
    vectorField UtFs = muFluidA().value()*DnA*UtPA;

    // Add effects of surface-tension gradient
    UtFs += surfaceTensionGrad()().internalField();

    // Normal
    vectorField UnFs =
        nA*phi_.boundaryField()[patchID()]/mesh().boundary()[patchID()].magSf();

    Us().internalField() += UnFs - nA*(nA&Us().internalField());
    correctUsBoundaryConditions();

    UtFs -= muFluidA().value()*
        (fac::grad(Us())&areaMesh().faceAreaNormals())().internalField();

    UtFs /= muFluidA().value()*DnA + VSMALL;

    // Update free-surface velocity
    Us().internalField() = UtFs + UnFs;
    correctUsBoundaryConditions();

    // Update fixedGradient boundary condition on the free-surface
    vectorField nGradU =
        -nA*muFluidA().value()*fac::div(Us())().internalField();

    nGradU += muFluidA().value()*
        (UtFs - UtPA)*mesh().boundary()[patchID()].deltaCoeffs();

    nGradU /= muFluidA().value() + VSMALL;

    if
    (
        U_.boundaryField()[patchID()].type()
     == "fixedGradient"
    )
    {
        fixedGradientFvPatchField<vector>& aU =
            refCast<fixedGradientFvPatchField<vector> >
            (
                U_.boundaryField()[patchID()]
            );

        aU.gradient() == nGradU;
    }
    else
    {
        FatalErrorIn("fluidInterface::updateVelocity()")
            << "Boundary condition on " << U_.name()
            << " is " << U_.boundaryField()[patchID()].type()
            << ", instead of "
            << fixedGradientFvPatchField<vector>::typeName
            << abort(FatalError);
    }
    
}

tmp<areaVectorField> Foam::fluidInterface::surfaceTensionGrad()
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
            -fac::grad(surfaceTension())()
        )
    );

    return tgrad;
}

// Update on topology change
void Foam::fluidInterface::updateMesh(const mapPolyMesh& mpm)
{
    // Wipe out demand-driven data
    clearOut();
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::fluidInterface::operator=(const fluidInterface& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("fluidInterface::operator=(const fluidInterface&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}


// ************************************************************************* //
