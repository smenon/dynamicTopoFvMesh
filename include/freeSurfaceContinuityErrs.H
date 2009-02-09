Info<< "Free surface flux: sum local = "
    << sum(mag(phi.boundaryField()[interface.aPatchID()]))
    << ", global = "
    << sum(phi.boundaryField()[interface.aPatchID()]) << endl;

scalarField netPhiA =
    phi.boundaryField()[interface.aPatchID()]
  - fvc::meshPhi(U)().boundaryField()[interface.aPatchID()];

Info<< "Free surface continuity error : sum local = " 
    << sum(mag(netPhiA)) << ", global = " << sum(netPhiA) 
    << endl << endl;

