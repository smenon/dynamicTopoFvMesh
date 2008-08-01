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
    
    // Introduce the non-Newtonian transport model
    autoPtr<viscosityModel> nuModel = viscosityModel::New("nu",interface,U,phi);
    
    volScalarField cond
    (
        IOobject
        (
            "cond",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluidIndicator*
        (
            interface.condFluidA()
          - interface.condFluidB()
        )
      + interface.condFluidB()
    );

    volScalarField Cp
    (
        IOobject
        (
            "Cp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluidIndicator*
        (
            interface.CpFluidA()
          - interface.CpFluidB()
        )
      + interface.CpFluidB()
    );

    surfaceScalarField phiNet = fvc::interpolate(rho)*phi;
    phiNet.rename("phiNet");    

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p,mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

    volScalarField rUA
    (
        IOobject
        (
            "rUA",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("rUA", dimTime/dimDensity, 1),
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

