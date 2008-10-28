    
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
            fixedValuePointPatchVectorField::typeName          
        );

        pNormals.boundaryField()[interface.patchID()] == interface.areaMesh().pointAreaNormals();
    
        pNormals.write();

        // Write out the mesh length scales
        mesh.lengthScale()().write();
        
        // Write out the pressure gradient
        volVectorField gradP = fvc::grad(p);
        gradP *= -1.0; gradP.write();
    }
