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

    fluidInterface interface(mesh, rho, U, p, phi);   
    
    // Introduce the non-Newtonian transport model
    // autoPtr<viscosityModel> nuModel = viscosityModel::New("nu",interface,U,phi);
    
    dimensionedScalar nu
    (
        IOdictionary
        (
            IOobject
            (
                "fluidInterfaceProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).lookup("nu")
    );    

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p,mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

