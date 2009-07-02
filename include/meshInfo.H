    
    if (runTime.outputTime()) 
    {
        // Write out the mesh length scales
        volScalarField lengthScale
        (
            IOobject
            (
                "lengthScale",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("scalar", dimLength, 0.0),
            "zeroGradient"
        );

        lengthScale.internalField() = mesh.lengthScale();
        lengthScale.write();

        // Write out mesh quality
        volScalarField meshQuality
        (
            IOobject
            (
                "meshQuality",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("scalar", dimless, 0.0),
            "zeroGradient"
        );

        meshQuality.internalField() = mesh.meshQuality(true);
        meshQuality.write();
    }
