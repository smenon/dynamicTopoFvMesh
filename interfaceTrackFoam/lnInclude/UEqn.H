    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(fvc::interpolate(rho)*phi, U)
      - fvm::laplacian(mu, U)
    );

    if (momentumPredictor)
    {
        solve(UEqn == -fvc::grad(p));
    }
