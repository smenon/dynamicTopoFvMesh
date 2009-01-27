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
    
    Info<< "Reading field T\n" << endl << flush;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );    

#   include "createPhi.H"     

    fluidInterface interface(mesh, U, p, phi);

    volScalarField nu
    (
        IOobject
        (
            "nu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        interface.muFluidA()/interface.rhoFluidA()
    );
        
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p,mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

