Info<< "Free surface flux: sum local = "
    << sum(mag(phi.boundaryField()[interface.patchID()])) 
    << ", global = "
    << sum(phi.boundaryField()[interface.patchID()]) << endl;

scalarField netPhiA =
    phi.boundaryField()[interface.patchID()]
  - fvc::meshPhi(U)().boundaryField()[interface.patchID()];

Info<< "Free surface continuity error : sum local = " 
    << sum(mag(netPhiA)) << ", global = " << sum(netPhiA) 
    << endl << endl;

