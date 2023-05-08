volScalarField rUA("rUA", 1.0/ddtUEqn.A());

phi = (fvc::interpolate(U) & mesh.Sf());
adjustPhi(phi, U, p);

while (piso.correctNonOrthogonal())
{
	fvScalarMatrix pEqn
	(
		fvm::laplacian(rUA, dp) == fvc::div(phi)
	);

        pEqn.setReference(dpRefCell, dpRefValue);
        pEqn.solve
        (
        	mesh.solver(dp.select(piso.finalInnerIter()))
        );

        if (piso.finalNonOrthogonalIter())
        {
		phi -= pEqn.flux();
		p += dp;
	}

	#include "continuityErrs.H"

	U -= rUA*fvc::grad(dp);
	U.correctBoundaryConditions();
}
