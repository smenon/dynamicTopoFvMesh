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
    
    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );

#   include "createPhi.H"

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p,mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

    Info<< "Reading field rUA if present\n" << endl;
    volScalarField rUA
    (
        IOobject
        (
            "rUA",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        runTime.deltaT(),
        zeroGradientFvPatchScalarField::typeName
    );
    
    {
        const dictionary& piso = mesh.solutionDict().subDict("PISO");

        bool momentumPredictor = true;
        if (piso.found("momentumPredictor"))
        {
            momentumPredictor = Switch(piso.lookup("momentumPredictor"));
        }

        bool momentumPredictorSave = momentumPredictor;
        momentumPredictor = false;
#       include "UEqn.H"
        momentumPredictor = momentumPredictorSave;

        rUA = 1.0/UEqn.A();
    }