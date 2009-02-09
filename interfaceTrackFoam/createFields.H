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

    // Obtain fluid indicator from the interface
    volScalarField* fluidIndicatorPtr = NULL;

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

    volScalarField& fluidIndicator = *fluidIndicatorPtr;

    volScalarField nu
    (
        IOobject
        (
            "nu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluidIndicator*
        (
            (interface.muFluidA()/interface.rhoFluidA())
          - (interface.muFluidB()/interface.rhoFluidB())
        )
      + (interface.muFluidB()/interface.rhoFluidB())
    );
        
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p,mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

