#include "UnknownsProperties.h"
#include "EulerEquations.h"
#include "Funct.h"

void Unknowns::RecalculateVariables(const ConservativeVariables& vars, 
	const Properties& prop)
{
	for (auto i = 1; i < u.size() - 1; ++i)
	{
		const double w1((*vars(0))[i]), w2((*vars(1))[i]), w3((*vars(2))[i]);
		p[i] = w1 * (prop.gamma - 1.) * (w3 / w1 - pow(w2 / w1, 2.) / 2.);
		u[i] = w2 / w1;
		rho[i] = w1;
	}
}