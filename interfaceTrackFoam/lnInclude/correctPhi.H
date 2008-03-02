{
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

#   include "continuityErrs.H"

    for(int nonOrth=0; nonOrth<=nPCorr; nonOrth++)
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rUA, pcorr) == fvc::div(phi)
        );

        pcorrEqn.setReference(pRefCell, pRefValue);
        pcorrEqn.solve();

        if (nonOrth == nPCorr)
        {
            phi -= pcorrEqn.flux();
        }
    }  
}


