    
    if (runTime.outputTime()) 
    {       
        // Instantiate a pointMesh object
        pointMesh pMesh(mesh);

        pointVectorField pNormals
        (
            IOobject
            (
                "pointNormals",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedVector("vector", dimless, vector::zero),
            "fixedValue"
        );

        pNormals.boundaryField()[interface.aPatchID()] ==
            interface.areaMesh().pointAreaNormals();
    
        pNormals.write();

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
        
        // Write out the pressure gradient
        volVectorField gradP = fvc::grad(p);
        gradP *= -1.0; gradP.write();
    }
