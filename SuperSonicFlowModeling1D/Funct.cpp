#include "Funct.h"
#include "PhysLawsFunct.h"
#include "Schemes.h"

std::vector<Funct*> ConservativeVariablesFactory::create(Unknowns& vars,
	Funct* TotalEnergy)
{
	return { new W1Funct(vars),
				new W2Funct(vars),
					new W3Funct(vars, TotalEnergy) };
}

std::vector<Funct*> InviscidFluxesFactory::create(Unknowns& vars,
	Funct* TotalEnthalpy) 
{
	return { new F1Funct(vars),
				new F2Funct(vars),
					new F3Funct(vars, TotalEnthalpy) };
}

namespace RoeFunct
{
	Array Dissipation::operator()(uint i)
	{
		const double del_p(vars.p[i + 1] - vars.p[i]);
		const double del_u(vars.u[i + 1] - vars.u[i]);
		const double del_rho(vars.rho[i + 1] - vars.rho[i]);
		const double coeff1(abs(u(i) + c(i)) *
			(del_p + rho(i) * c(i) * del_u) / (2. * pow(c(i), 2.0)));
		const double coeff2(abs(u(i) - c(i)) *
			(del_p - rho(i) * c(i) * del_u) / (2. * pow(c(i), 2.0)));
		const double coeff3(abs(u(i)) * (del_rho - del_p / pow(c(i), 2.0)));
		return 
			coeff1 * Array{ 1., u(i) + c(i), H(i) + c(i) * u(i) } +
			coeff2 * Array{ 1., u(i) - c(i), H(i) - c(i) * u(i) } +
			coeff3 * Array{ 1., u(i), pow(u(i), 2) / 2. };
	}
}