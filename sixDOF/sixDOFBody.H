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
    Six-DOF Body

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    sixDOFBody.C

\*---------------------------------------------------------------------------*/

#ifndef sixDOFBody_H
#define sixDOFBody_H

#include "ODE.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class ODESolver;

/*---------------------------------------------------------------------------*\
                           Class sixDOFBody Declaration
\*---------------------------------------------------------------------------*/

class sixDOFBody
{
    // Private data

        //- ODE class for six-DOF
        class odeSixDOF
        :
            public ODE
        {
            // Private data

                //- ODE coefficients
                scalarField coeffs_;

        public:

            // Constructors

                odeSixDOF();

            // Destructor

                virtual ~odeSixDOF()
                {}

            // Member Functions

                //- Return number of equations
                virtual label nEqns() const
                {
                    return 6;
                }

                //- Return reference to interpolation coefficients
                virtual scalarField& coeffs()
                {
                    return coeffs_;
                }

                //- Return reference to interpolation coefficients
                virtual const scalarField& coeffs() const
                {
                    return coeffs_;
                }

                //- Return derivatives
                virtual void derivatives
                (
                    const scalar,
                    const scalarField&,
                    scalarField&
                ) const;

                //- Return Jacobian
                virtual void jacobian
                (
                    const scalar x,
                    const scalarField& y,
                    scalarField& dfdx,
                    scalarSquareMatrix& dfdy
                ) const;

                //- Update ODE after the solution, advancing by delta
                virtual void update(const scalar delta);
        };

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to base dictionary
        const dictionary& dict_;

        //- Reference to body sub-dictionary
        const dictionary& bodyDict_;

        //- Index of body
        label index_;

        //- Patch index of body
        label patchIndex_;

        //- Six DOF ODE
        odeSixDOF ode_;

        //- ODE solver
        autoPtr<ODESolver> odeSolver_;

        //- ODE solvers options
        scalar eps_, hEst_;

        //- Transform tensors
        tensor globalToLocal_, localToGlobal_;

        vector weightFactor_;
        vector localVelocity_, localOldVelocity_, globalVelocity_;
        vector localOmega_, localOldOmega_, globalOmega_;
        vector displacement_, totalDisplacement_;
        vector rotation_, totalRotation_;
        vector fAvg_, mAvg_, fTotal_;

        //- Stored forces / moments
        FixedList<vector, 3> Fs_, Ms_;

        //- Min / max acceleration
        scalar aMin_, aMax_;

        //- Initial centre-of-gravity / rotation
        vector Cg_, Rot_;

        //- DOFs to calculate
        FixedList<bool, 3> calcTransDOF_, calcRotDOF_;

        //- Constant external force / moment
        vector fConst_, mConst_;

        //- Body mass / moment-of-inertia
        scalar mass_;
        vector J_;

        // Linear spring / damping constants
        vector kTrans_, dTrans_;
        vector kRot_, dRot_;

        // Names for relevant fields in objectRegistry
        word pName_, gName_;

    // Private member functions

        //- Calculate transformation tensor
        tensor calcTensor(const vector& rot) const;

public:

    // Constructors

        //- Construct from components
        sixDOFBody
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& bodyDict,
            const label index
        );

    // Destructor

        virtual ~sixDOFBody();

    // Member Functions

        //- Calculate forces
        void calculateForces();

        //- Balance forces
        void balanceForces();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

