    // Define the engine parameters
    IOdictionary engineDict
    (
        IOobject
        (
            "engineDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Obtain the piston patch name
    dictionary piston(engineDict.subDict("piston"));
    wordList pistonName = piston.toc();

    // Read in necessary information from the dictionary
    dictionary pistonDict = engineDict.subDict(pistonName[0]);

    // Read the profile
    interpolationTable<scalar> pistonTable
    (
        fileName(pistonDict.lookup("profileFile"))
    );

    // Set outOfBounds handling to clamp
    pistonTable.outOfBounds(interpolationTable<scalar>::CLAMP);

    // Read in the piston axis
    vector pistonAxis = pistonDict.lookup("axis");

    pistonAxis /= mag(pistonAxis) + VSMALL;

    scalar oldStroke = pistonTable(runTime.value()), currentStroke = 0.0;

    // Obtain the number of valves in the system
    dictionary valves(engineDict.subDict("valves"));
    wordList valveList = valves.toc();

    label numValves = valveList.size();

    PtrList<interpolationTable<scalar> > valveLiftTables(numValves);
    List<scalar> oldLift(numValves, 0.0), currentLift(numValves, 0.0);
    List<vector> valveAxes(numValves, vector::zero);

    // Read in all necessary information from the dictionary
    forAll(valveList, valveI)
    {
        dictionary valveDict = engineDict.subDict(valveList[valveI]);

        // Read in the lift profile table
        valveLiftTables.set
        (
            valveI,
            new interpolationTable<scalar>
            (
                fileName(valveDict.lookup("profileFile"))
            )
        );

        // Set outOfBounds handling to clamp
        valveLiftTables[valveI].outOfBounds
        (
            interpolationTable<scalar>::CLAMP
        );

        // Read in the valve axis
        valveAxes[valveI] = valveDict.lookup("axis");

        // Normalize the axis
        valveAxes[valveI] /= mag(valveAxes[valveI]) + VSMALL;

        // Interpolate and obtain the lift value
        // for the current time-step
        oldLift[valveI] = valveLiftTables[valveI](runTime.value());
    }
