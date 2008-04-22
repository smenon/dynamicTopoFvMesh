{
volScalarField DpDt = fvc::DDt(phi,p); // here, phi is volumetric flux
       
        solve
        (
            Cp*fvm::ddt(rho, T)
          + Cp*fvm::div(phiNet, T)
          - fvm::laplacian(cond,T)
         ==
             DpDt 
        );

}
