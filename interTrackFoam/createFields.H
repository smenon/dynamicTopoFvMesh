
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

Info<< "Creating the interface\n" << endl;
freeSurface interface(mesh, U, p, phi);

Info<< "Reading field rAU if present\n" << endl;
volScalarField rAU
(
    IOobject
    (
        "rAU",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    runTime.deltaT()/interface.rho()
);

#   include "setRefCell.H"

