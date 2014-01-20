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
        surfaceScalarField& phi =
        (
            const_cast<surfaceScalarField&>
            (
                mesh().lookupObject<surfaceScalarField>(phiName)
            )
        );

        const volVectorField& U = mesh().lookupObject<volVectorField>(UName);

        // Interpolate mapped velocity to faces
        tmp<surfaceScalarField> tphiU = fvc::interpolate(U) & mesh().Sf();

        // Alias for convenience
        surfaceScalarField& phiU = tphiU();

        phiU.rename("phiU");

        forAll(faces, faceI)
        {
            phi[faces[faceI]] = phiU[faces[faceI]];
//            phi[faces[faceI]] = 0.0;
        }

        // Over-write phi entirely
//        phi = phiU;

//        forAll(phi.internalField(), faceI)
//        {
//            phi.internalField()[faceI] = 0.0;
//        }
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
    word UName(subDict.lookup("U"));
    word rAUName(subDict.lookup("rAU"));
    word phiName(subDict.lookup("phi"));

    // Fetch references
    const fvMesh& mesh = fluxCorrector::mesh();
    const Time& runTime = mesh.time();

    const volScalarField& p = mesh.lookupObject<volScalarField>(pName);
    //const volVectorField& U = mesh.lookupObject<volVectorField>(UName);
    const volScalarField& rAU = mesh.lookupObject<volScalarField>(rAUName);

    surfaceScalarField& phi =
    (
        const_cast<surfaceScalarField&>
        (
            mesh.lookupObject<surfaceScalarField>(phiName)
        )
    );

    wordList pcorrTypes(p.boundaryField().types());

    for (label i=0; i<p.boundaryField().size(); i++)
    {
        if (p.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
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

    dictionary piso = mesh.solutionDict().subDict("PISO");

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, piso, pRefCell, pRefValue);

    int nNonOrthCorr = 0;
    if (piso.found("nNonOrthogonalCorrectors"))
    {
        nNonOrthCorr = readInt(piso.lookup("nNonOrthogonalCorrectors"));
    }

    //dimensionedScalar rAUf("(1|A(U))", dimTime, 1.0);

    // Adjust dimensions for cases involving
    // density-normalised pressure
    //if (p.dimensions() == dimPressure)
    //{
    //    rAU.dimensions() /= dimDensity;
    //}

    //adjustPhi(phi, U, pcorr);

#   include "initContinuityErrs.H"
#   include "continuityErrs.H"

    {
#       include "CourantNo.H"
    }

    // Write out phi prior to correction
//    surfaceVectorField prePhi = phi * (mesh.Sf() / mesh.magSf());
//    prePhi.rename("prePhi");

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAU, pcorr) == fvc::div(phi)
        );

        pcorrEqn.setReference(pRefCell, pRefValue);
        pcorrEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pcorrEqn.flux();
        }
    }

#   include "continuityErrs.H"

    {
#       include "CourantNo.H"
    }

    // Write out phi prior to correction
//    surfaceVectorField postPhi = phi * (mesh.Sf() / mesh.magSf());
//    postPhi.rename("postPhi");
//
//    if (mesh.time().outputTime())
//    {
//        prePhi.write();
//        postPhi.write();
//    }
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
