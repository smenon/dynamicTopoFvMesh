
    // Cell-skewness values
    /*
    volScalarField skewness
    (
        IOobject
        (
            "skewness",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("cellSkewness", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    );        
    scalarField meshskewness = mesh.skewness()();
    forAll(skewness,index) skewness[index] = meshskewness[index];
    */
    
    // Mesh length-density function
    volScalarField lDensity
    (
        IOobject
        (
            "lengthDensity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("lDensity", dimless, 0),
        zeroGradientFvPatchScalarField::typeName            
    );
    
    // Obtain length-density values from the mesh
    HashList<scalar> &meshLengthDensity = mesh.lengthDensity();
    for(HashList<scalar>::iterator ldIter = meshLengthDensity.begin(); ldIter != meshLengthDensity.end(); ldIter++) {
        lDensity.internalField()[ldIter.index()] = ldIter();
    }
