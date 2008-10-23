    
    if (runTime.outputTime()) 
    {        
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
            fixedValuePointPatchVectorField::typeName          
        );

        pNormals.boundaryField()[interface.patchID()] == interface.areaMesh().pointAreaNormals();
    
        pNormals.write();
        lDensity.write();
    }
