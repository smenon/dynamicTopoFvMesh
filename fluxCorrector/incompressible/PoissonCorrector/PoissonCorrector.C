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
    PoissonCorrector

Description
    Flux-correction after topo-changes, using a Poisson solver.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    PoissonCorrector.C

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "PoissonCorrector.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(PoissonCorrector, 0);
addToRunTimeSelectionTable(fluxCorrector, PoissonCorrector, mesh);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from fvMesh and dictionary
PoissonCorrector::PoissonCorrector
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fluxCorrector(mesh, dict),
    required_(dict.subDict("PoissonCorrector").lookup("correctFluxes"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Is flux-correction required?
bool PoissonCorrector::required() const
{
    return required_;
}


//- Optionally cache an old-flux field to maintain divergence
void PoissonCorrector::cacheFluxes() const
{
    if (required())
    {
        const dictionary& subDict = dict().subDict("PoissonCorrector");

        // Search the dictionary for field information
        word UName(subDict.lookup("U"));
        word phiName(subDict.lookup("phi"));

        // Lookup fields from the registry
        const volVectorField& U =
        (
            mesh().lookupObject<volVectorField>(UName)
        );

        const surfaceScalarField& phi =
        (
            mesh().lookupObject<surfaceScalarField>(phiName)
        );

        // Cache the old divergence of the velocity field
        divUPtr_.reset
        (
            new volScalarField
            (
                "divU_cached",
                fvc::div(fvc::absolute(phi, U))
            )
        );
    }
}


//- Interpolate fluxes to a specified list of faces
void PoissonCorrector::interpolateFluxes(const labelList& faces) const
{
    if (required())
    {
        const dictionary& subDict = dict().subDict("PoissonCorrector");

        // Search the dictionary for field information
        word UName(subDict.lookup("U"));
        word phiName(subDict.lookup("phi"));

        // Lookup fields from the registry
        volVectorField& U =
        (
            const_cast<volVectorField&>
            (
                mesh().lookupObject<volVectorField>(UName)
            )
        );

        surfaceScalarField& phi =
        (
            const_cast<surfaceScalarField&>
            (
                mesh().lookupObject<surfaceScalarField>(phiName)
            )
        );

        const surfaceVectorField& Sf = mesh().Sf();

        // Interpolate mapped velocity to faces
        tmp<surfaceScalarField> tphiU = fvc::interpolate(U) & Sf;

        // Alias for convenience
        surfaceScalarField& phiU = tphiU();

        phiU.rename("phiU");

        // Over-write phi entirely
        phi = phiU;

        // If the mesh is moving, correct the velocity BCs on the moving walls
        // to ensure the corrected fluxes and velocity are consistent
        if (mesh().moving())
        {
            forAll(U.boundaryField(), patchi)
            {
                if (U.boundaryField()[patchi].fixesValue())
                {
                    U.boundaryField()[patchi].initEvaluate();
                }
            }

            forAll(U.boundaryField(), patchi)
            {
                if (U.boundaryField()[patchi].fixesValue())
                {
                    U.boundaryField()[patchi].evaluate();

                    phi.boundaryField()[patchi] =
                    (
                        U.boundaryField()[patchi] & Sf.boundaryField()[patchi]
                    );
                }
            }
        }
    }
}


//- Update fluxes in the registry, if required
void PoissonCorrector::updateFluxes() const
{
    if (!required())
    {
        return;
    }

    const dictionary& subDict = dict().subDict("PoissonCorrector");

    // Search the dictionary for field information
    word pName(subDict.lookup("p"));
    word rAUName(subDict.lookup("rAU"));
    word phiName(subDict.lookup("phi"));
    word solName(subDict.lookup("solutionDict"));

    // Fetch references
    const fvMesh& mesh = fluxCorrector::mesh();
    const Time& runTime = mesh.time();

    // Lookup fields from the registry
    const volScalarField& p = mesh.lookupObject<volScalarField>(pName);
    const volScalarField& rAU = mesh.lookupObject<volScalarField>(rAUName);

    surfaceScalarField& phi =
    (
        const_cast<surfaceScalarField&>
        (
            mesh.lookupObject<surfaceScalarField>(phiName)
        )
    );

    // Initialize BCs list for pcorr to zero-gradient
    wordList pcorrTypes
    (
        p.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    // Set BCs of pcorr to fixed-value for patches at which p is fixed
    forAll(p.boundaryField(), patchI)
    {
        if (p.boundaryField()[patchI].fixesValue())
        {
            pcorrTypes[patchI] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pcorr", p.dimensions(), 0.0),
        pcorrTypes
    );

    // Fetch the relevant dictionary (PISO, SIMPLE, PIMPLE, etc)
    dictionary solDict = mesh.solutionDict().subDict(solName);

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    label nNonOrthCorr = 0;

    setRefCell(p, solDict, pRefCell, pRefValue);

    if (solDict.found("nNonOrthogonalCorrectors"))
    {
        nNonOrthCorr = readLabel(solDict.lookup("nNonOrthogonalCorrectors"));
    }

    for (label nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAU, pcorr) == fvc::div(phi) - divUPtr_()
        );

        pcorrEqn.setReference(pRefCell, pRefValue);
        pcorrEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pcorrEqn.flux();
        }
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //
void PoissonCorrector::operator=(const PoissonCorrector& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("PoissonCorrector::operator=(const PoissonCorrector&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
