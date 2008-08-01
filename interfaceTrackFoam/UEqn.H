    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(fvc::interpolate(rho)*phi, U)
      - fvm::laplacian(nuModel->nu()*rho, U)
    );

    if (momentumPredictor)
    {
        solve(UEqn == -fvc::grad(p));
    }
