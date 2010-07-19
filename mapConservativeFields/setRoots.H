    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("source dir");

    argList::validOptions.insert("sourceTime", "scalar");
    argList::validOptions.insert("method", "label");
    argList::validOptions.insert("nThreads", "label");
    argList::validOptions.insert("forceRecalc", "");
    argList::validOptions.insert("writeAddr", "");

    argList args(argc, argv);

    if (!args.check())
    {
         FatalError.exit();
    }

    fileName rootDirTarget(args.rootPath());
    fileName caseDirTarget(args.globalCaseName());

    fileName casePath(args.additionalArgs()[0]);
    fileName rootDirSource = casePath.path();
    fileName caseDirSource = casePath.name();

    Info<< "Source: " << rootDirSource << " " << caseDirSource << nl
        << "Target: " << rootDirTarget << " " << caseDirTarget << endl;

    label nThreads = 1;

    if (args.options().found("nThreads"))
    {
        nThreads = readLabel(IStringStream(args.options()["nThreads"])());
    }

    label method = 0;

    if (args.options().found("method"))
    {
        method = readLabel(IStringStream(args.options()["method"])());
    }

    bool forceRecalc = false;

    if (args.options().found("forceRecalc"))
    {
        forceRecalc = true;
    }

    bool writeAddr = false;

    if (args.options().found("writeAddr"))
    {
        writeAddr = true;
    }
