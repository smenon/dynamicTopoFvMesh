Info<< "\nReading field p\n" << endl << flush;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl << flush;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


#   include "createPhi.H"


volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimMass/dimVolume, 0)
);


freeSurface interface(mesh, rho, U, p, phi);


volScalarField* fluidIndicatorPtr = NULL;

if(args.options().found("parallel"))
{
    fluidIndicatorPtr = new volScalarField
    (
        IOobject
        (
            "fluidIndicator",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
}
else
{
    fluidIndicatorPtr = new volScalarField
    (
        IOobject
        (
            "fluidIndicator",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        interface.fluidIndicator()
    );
}

volScalarField& fluidIndicator = *fluidIndicatorPtr;


rho = fluidIndicator*
    (
        interface.rhoFluidA() 
      - interface.rhoFluidB()
    ) 
  + interface.rhoFluidB();


volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fluidIndicator*
    (
        interface.muFluidA() 
      - interface.muFluidB()
    ) 
  + interface.muFluidB()
);


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

