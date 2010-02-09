        const polyBoundaryMesh& boundary = mesh.boundaryMesh();

        // Copy existing point locations
        pointField meshPoints(mesh.points());

        // Update the piston
        currentStroke = pistonTable(runTime.value());

        const labelList& pistonPoints =
        (
            boundary[boundary.findPatchID(pistonName[0])].meshPoints()
        );

        vector pD = (currentStroke - oldStroke)*pistonAxis;

        Info << "Piston " << pistonName[0]
             << ": Current Stroke value: " << currentStroke << endl;
        Info << "Piston " << pistonName[0]
             << ": Old Stroke value: " << oldStroke << endl;

        // Check if a fixed-constant over-ride displacement
        // has been specified.
        if (pistonDict.found("displacement"))
        {
            scalar pDisp = readScalar(pistonDict.lookup("displacement"));

            pD = (pDisp*pistonAxis);
        }

        forAll(pistonPoints, index)
        {
            meshPoints[pistonPoints[index]] += pD;
        }

        oldStroke = currentStroke;

        vector vD = vector::zero;

        // Update valves
        forAll(valveList, valveI)
        {
            // Obtain the current lift value
            currentLift[valveI] = valveLiftTables[valveI](runTime.value());

            // Move the valve patch by the difference in lift
            const labelList& valvePoints =
            (
                boundary[boundary.findPatchID(valveList[valveI])].meshPoints()
            );

            // Define displacement
            vD = (currentLift[valveI] - oldLift[valveI])*valveAxes[valveI];

            Info << "Valve " << valveList[valveI]
                 << ": Current Lift value: " << currentLift[valveI] << endl;
            Info << "Valve " << valveList[valveI]
                 << ": Old Lift value: " << oldLift[valveI] << endl;

            dictionary valveDict = engineDict.subDict(valveList[valveI]);

            // Check if a fixed-constant over-ride displacement
            // has been specified.
            if (valveDict.found("displacement"))
            {
                scalar vDisp = readScalar(valveDict.lookup("displacement"));

                vD = (vDisp*valveAxes[valveI]);
            }

            forAll(valvePoints, index)
            {
                meshPoints[valvePoints[index]] += vD;
            }

            oldLift[valveI] = currentLift[valveI];
        }

        // Obtain the field from the registry
        pointField& refPoints = const_cast<pointField&>
        (
            mesh.lookupObject<pointField>("refPoints")
        );

        // Assign boundary conditions to the motion solver
        refPoints = meshPoints;
