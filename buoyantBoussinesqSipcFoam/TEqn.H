alphat = turbulence->nut()/Prt;
alphat.correctBoundaryConditions();

volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

fvScalarMatrix TEqn
(
    fvm::ddt(T)
  + fvm::div(phi, T)
//      	  - fvm::Sp(fvc::div(phi), T)
  - fvm::laplacian(alphaEff, T)
  ==
    fvOptions(T)
);

TEqn.relax();
TEqn.solve();
//	T.correctBoundaryConditions();
rhok = - beta*(T - TRef);
