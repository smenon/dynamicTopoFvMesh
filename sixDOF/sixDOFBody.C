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
    sixDOFBody

Description
    Six DOF Body

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "sixDOFBody.H"
#include "ODESolver.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"
#include "uniformDimensionedFields.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Calculate transformation tensor
//   - Assumes rot is in radians
tensor sixDOFBody::calcTensor(const vector& rot) const
{
    scalar crX = Foam::cos(rot.x()), srX = Foam::sin(rot.x());
    scalar crY = Foam::cos(rot.y()), srY = Foam::sin(rot.y());
    scalar crZ = Foam::cos(rot.x()), srZ = Foam::sin(rot.x());

    tensor Rx
    (
        1.0, 0.0, 0.0,
        0.0, crX,-srX,
        0.0, srX, crX
    );

    tensor Ry
    (
        crY, 0.0, srY,
        0.0, 1.0, 0.0,
       -srY, 0.0, crY
    );

    tensor Rz
    (
        crZ,-srZ, 0.0,
        srZ, crZ, 0.0,
        0.0, 0.0, 1.0
    );

    // Return rotation matrix
    return (Rz & Ry & Rx);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
sixDOFBody::sixDOFBody
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& bodyDict,
    const label index
)
:
    mesh_(mesh),
    dict_(dict),
    bodyDict_(bodyDict),
    index_(index),
    patchIndex_
    (
        mesh.boundaryMesh().findPatchID
        (
            bodyDict.lookup("patchName")
        )
    ),
    odeSolver_
    (
        ODESolver::New
        (
            dict.subDict("ODECoeffs").lookup("ODESolver"),
            ode_
        )
    ),
    eps_(readScalar(dict.subDict("ODECoeffs").lookup("eps"))),
    hEst_(readScalar(dict.subDict("ODECoeffs").lookup("hEst"))),
    weightFactor_(bodyDict.lookup("weightFactor")),
    localVelocity_(vector::zero),
    localOldVelocity_
    (
        bodyDict.lookupOrDefault
        (
            "initialVelocity",
            vector::zero
        )
    ),
    globalVelocity_(vector::zero),
    localOmega_(vector::zero),
    localOldOmega_
    (
        bodyDict.lookupOrDefault
        (
            "initialAngularVelocity",
            vector::zero
        )
    ),
    globalOmega_(vector::zero),
    displacement_(vector::zero),
    totalDisplacement_(vector::zero),
    rotation_(vector::zero),
    totalRotation_(vector::zero),
    fAvg_(vector::zero),
    mAvg_(vector::zero),
    fTotal_(vector::zero),
    Fs_(vector::zero),
    Ms_(vector::zero),
    aMin_(bodyDict.lookupOrDefault("minAcceleration", -GREAT)),
    aMax_(bodyDict.lookupOrDefault("maxAcceleration", GREAT)),
    Cg_(bodyDict.lookup("CentreOfGravity")),
    Rot_(bodyDict.lookup("Rotation")),
    calcTransDOF_(bodyDict.lookup("calcTranslationDOF")),
    calcRotDOF_(bodyDict.lookup("calcRotationDOF")),
    fConst_
    (
        bodyDict.lookupOrDefault
        (
            "constantForce",
            vector::zero
        )
    ),
    mConst_
    (
        bodyDict.lookupOrDefault
        (
            "constantMoment",
            vector::zero
        )
    ),
    mass_(readScalar(bodyDict.lookup("mass"))),
    J_(bodyDict.lookup("MomentOfInertia")),
    kTrans_(bodyDict.lookup("translationSpringConstant")),
    dTrans_(bodyDict.lookup("translationDampingConstant")),
    kRot_(bodyDict.lookup("rotationSpringConstant")),
    dRot_(bodyDict.lookup("rotationDampingConstant")),
    pName_(bodyDict.lookupOrDefault("pName", word("p"))),
    gName_(bodyDict.lookupOrDefault("gName", word("g")))
{
    // Convert rotation to radians, if necessary
    {
        Rot_ *= (mathematicalConstant::pi / 180.0);
    }

    // Initialize transformation tensors
    globalToLocal_ = calcTensor( 1.0 * Rot_);
    localToGlobal_ = calcTensor(-1.0 * Rot_);

    // Attempt to read restart info from disk
    IOdictionary sixDOFDict
    (
        IOobject
        (
            "sixDOFDict",
            mesh_.time().timeName(),
            "uniform",
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    if (sixDOFDict.headerOk())
    {
        Info<< "Reading sixDOF information from dictionary" << endl;

        totalDisplacement_ = sixDOFDict.lookup("totalDisplacement");
        totalRotation_ = sixDOFDict.lookup("totalRotation");
        localOldVelocity_ = sixDOFDict.lookup("localOldVelocity");
        localOldOmega_ = sixDOFDict.lookup("localOldOmega");
        Fs_[0] = sixDOFDict.lookup("Fs0");
        Fs_[1] = sixDOFDict.lookup("Fs1");
        Ms_[0] = sixDOFDict.lookup("Ms0");
        Ms_[1] = sixDOFDict.lookup("Ms1");

        // Update CentreOfGravity
        Cg_ += totalDisplacement_;
    }
}


sixDOFBody::odeSixDOF::odeSixDOF()
:
    coeffs_(6, 0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

sixDOFBody::~sixDOFBody()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Calculate forces
void sixDOFBody::calculateForces()
{
    // Fetch pressure from objectRegistry
    const volScalarField& p = mesh_.lookupObject<volScalarField>(pName_);

    // Fetch gravity from objectRegistry
    vector g =
    (
        mesh_.lookupObject<uniformDimensionedVectorField>(gName_).value()
    );

    // Store forces / moments
    for (label i = 2; i > 0; i--)
    {
        Fs_[i] = Fs_[i-1];
        Ms_[i] = Ms_[i-1];
    }

    scalar sumWeightFactor = cmptSum(weightFactor_);

    // Update force
    Fs_[0] =
    (
        gSum
        (
            p.boundaryField()[patchIndex_]
          * (
                globalToLocal_
              & mesh_.Sf().boundaryField()[patchIndex_]
            )
        )
    );

    // Update centre-of-gravity
    Cg_ += (mesh_.time().deltaT().value() * localVelocity_);

    // Update moment
    Ms_[0] =
    (
        gSum
        (
            (
                globalToLocal_
              & (mesh_.Cf().boundaryField()[patchIndex_] - Cg_)
            )
          ^ (
                p.boundaryField()[patchIndex_]
              * (
                    globalToLocal_
                  & mesh_.Sf().boundaryField()[patchIndex_]
                )
            )
        )
    );

    // Calculate average forces / moments
    fAvg_ = vector::zero;
    mAvg_ = vector::zero;

    for (label j = 0; j < 3; j++)
    {
        for (label i = 0; i < 3; i++)
        {
            fAvg_[j] += (weightFactor_[i] * Fs_[i][j] / sumWeightFactor);
            mAvg_[j] += (weightFactor_[i] * Ms_[i][j] / sumWeightFactor);
        }
    }

    // Accumulate force
    fTotal_ = (fAvg_ + (globalToLocal_ & (mass_ * g)));

    // Write out a bit
    Info<< " Force on body: "
        << mesh_.boundaryMesh()[patchIndex_].name() << nl
        << " Pressure force: " << Fs_[0] << nl
        << " Moment: " << Ms_[0] << nl
        << endl;
}


//- Balance forces
void sixDOFBody::balanceForces()
{
    // Fetch time-step
    scalar deltaT = mesh_.time().deltaT().value();

    // Update total displacement / rotation
    totalDisplacement_ += (deltaT * localOldVelocity_);
    totalRotation_ += (rotation_ * 180.0 / mathematicalConstant::pi);

    // Solve for translation
    forAll(calcTransDOF_, dofI)
    {
        if (calcTransDOF_[dofI])
        {
            // Fetch coefficients
            scalarField& y = ode_.coeffs();

            scalar magTrans = mag(totalDisplacement_[dofI]);

            // Set coefficients
            y[0] = 0.0;
            y[1] = localOldVelocity_[dofI];
            y[2] = mass_;
            y[3] = dTrans_[dofI];
            y[4] = magTrans * kTrans_[dofI];
            y[5] = fTotal_[dofI] + fConst_[dofI];

            // Solve the ODE
            odeSolver_->solve
            (
                0.0,
                deltaT,
                eps_,
                hEst_
            );

            // Update velocity / displacement
            localVelocity_[dofI] = y[1];
            displacement_[dofI] = y[0];
        }
    }

    // Write out a bit
    Info<< " Velocity CG: " << localVelocity_ << nl
        << " Total Displacement: " << totalDisplacement_ << nl
        << endl;

    // Update old information
    localOldVelocity_ = localVelocity_;
    globalVelocity_ = (localToGlobal_ & localVelocity_);

    // Solve for rotation
    forAll(calcRotDOF_, dofI)
    {
        if (calcRotDOF_[dofI])
        {
            // Fetch coefficients
            scalarField& y = ode_.coeffs();

            scalar magRot =
            (
                mag
                (
                    totalRotation_[dofI]
                  * (mathematicalConstant::pi / 180.0)
                )
            );

            // Set coefficients
            y[0] = 0.0;
            y[1] = localOldOmega_[dofI];
            y[2] = J_[dofI];
            y[3] = dRot_[dofI];
            y[4] = magRot * kRot_[dofI];
            y[5] = mAvg_[dofI] + mConst_[dofI];

            // Solve the ODE
            odeSolver_->solve
            (
                0.0,
                deltaT,
                eps_,
                hEst_
            );

            // Update angular velocity / rotation
            localOmega_[dofI] = y[1];
            rotation_[dofI] = y[0];
        }
    }

    // Write out a bit
    Info<< " Angular Velocity CG: "
        << localOmega_ * 180.0 / mathematicalConstant::pi << nl
        << " Total Rotation: " << totalRotation_ << nl
        << endl;

    // Update old information
    localOldOmega_ = localOmega_;
    globalOmega_ = (localToGlobal_ & localOmega_);

    // Now set displacement BCs for mesh motion
    vectorField points(mesh_.boundaryMesh()[patchIndex_].localPoints() - Cg_);

    vectorField Ux(((points ^ vector(1,0,0)) & globalOmega_) * vector(1,0,0));
    vectorField Uy(((points ^ vector(0,1,0)) & globalOmega_) * vector(0,1,0));
    vectorField Uz(((points ^ vector(0,0,1)) & globalOmega_) * vector(0,0,1));

    // Calculate displacement for motionSolver
    vectorField motionDisplacement
    (
        (Ux + Uy + Uz + globalVelocity_) * deltaT
    );

    pointField& refPoints = const_cast<pointField&>
    (
        mesh_.lookupObject<pointField>("refPoints")
    );

    // Assign boundary conditions to the motion solver
    const labelList& meshPts = mesh_.boundaryMesh()[patchIndex_].meshPoints();

    forAll(meshPts,pointI)
    {
        refPoints[meshPts[pointI]] += motionDisplacement[pointI];
    }

    // If this is an output time-step, write out
    if (mesh_.time().outputTime())
    {
        IOdictionary sixDOFDict
        (
            IOobject
            (
                "sixDOFDict",
                mesh_.time().timeName(),
                "uniform",
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        // Add relevant info
        sixDOFDict.add("totalDisplacement", totalDisplacement_);
        sixDOFDict.add("totalRotation", totalRotation_);
        sixDOFDict.add("localOldVelocity", localOldVelocity_);
        sixDOFDict.add("localOldOmega", localOldOmega_);
        sixDOFDict.add("Fs0", Fs_[0]);
        sixDOFDict.add("Fs1", Fs_[1]);
        sixDOFDict.add("Ms0", Ms_[0]);
        sixDOFDict.add("Ms1", Ms_[1]);

        // Write it out
        sixDOFDict.regIOobject::write();
    }
}


void sixDOFBody::odeSixDOF::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    dydx[0] = y[1];
    dydx[1] = (y[5]-y[3]*y[1]-y[4]*y[0]) / y[2];
    dydx[2] = 0.0;
    dydx[3] = 0.0;
    dydx[4] = 0.0;
    dydx[5] = 0.0;
}


void sixDOFBody::odeSixDOF::jacobian
(
    const scalar x,
    const scalarField& y,
    scalarField& dfdx,
    scalarSquareMatrix& dfdy
) const
{
    notImplemented("odeSixDOF::jacobian() const");
}


void sixDOFBody::odeSixDOF::update(const scalar delta)
{

}


} // End namespace Foam

// ************************************************************************* //
