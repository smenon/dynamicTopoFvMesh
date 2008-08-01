    
    // Mesh length-density function
    volScalarField lDensity
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
        dimensionedScalar("lScale", dimless, 0),
        zeroGradientFvPatchScalarField::typeName            
    );
    
    // Obtain length-density values from the mesh
    HashList<scalar> &meshLengthDensity = mesh.lengthDensity();
    for(HashList<scalar>::iterator ldIter = meshLengthDensity.begin(); ldIter != meshLengthDensity.end(); ldIter++) {
        lDensity.internalField()[ldIter.index()] = ldIter();
    }
    
    if (runTime.outputTime()) lDensity.write();
