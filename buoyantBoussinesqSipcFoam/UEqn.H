// Convection-diffusion matrix
fvVectorMatrix HUEqn
(
  	fvm::div(phi, U) 
      + turbulence->divDevReff(U)
);

// Time derivative matrix
fvVectorMatrix ddtUEqn(fvm::ddt(U));

if (piso.momentumPredictor())
{
	solve
	(
	    ddtUEqn + HUEqn 
	    == 
	    fvc::reconstruct
            (
                (
                  fvc::interpolate(rhok)*ghf
		- fvc::snGrad(p)*mesh.magSf()
                )
            )
	);
}
