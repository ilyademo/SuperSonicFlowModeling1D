#include "Schemes.h"
#include "EulerEquations.h"
#include "MathOperations.h"
#include "GodunovFunctions.h"
#include "PhysLawsFunct.h"

bool FluxProcessor::CheckConvFluxDirection(Euler& eq, uint i)
{
	return (eq.vars.u[i + 1] + eq.vars.u[i]) / 2. >= 0;
}

Array FluxProcessor::FirstOrderUpdwind(Euler& eq, uint i)
{
	Array deltaF(3, 0);
	bool front = CheckConvFluxDirection(eq, i);
	bool back = CheckConvFluxDirection(eq, i - 1);
	Array unitSec{ 0, 1.0, 0 };
	if (front)
		deltaF += ((*(eq.F))[i] + eq.vars.p[i + 1] * unitSec);
	else
		deltaF += ((*(eq.F))[i + 1] + eq.vars.p[i] * unitSec);
	if (back)
		deltaF -= ((*eq.F)[i - 1] + eq.vars.p[i] * unitSec);
	else
		deltaF -= ((*eq.F)[i] + eq.vars.p[i - 1] * unitSec);
	return deltaF;
}

Array FluxProcessor::FortranGOGO(double pl, double pr, double ul, double ur,
	double rhol, double rhor, Properties& prop)
{
	double pbig(0.), ubig(0.), rbig1(0.), rbig2(0.), dl1(0.), dl2(0.),
		dp1(0.), dp2(0.), udot(0.);
	double pf(0.), uf(0.), rhof(0.);
	raspad(&(prop.gamma), &pl, &rhol, &ul,
		&pr, &rhor, &ur, &pbig, &ubig, &rbig1, &rbig2,
		&dl1, &dl2, &dp1, &dp2);
	potok(&(prop.gamma), &pl, &rhol, &ul,
		&pr, &rhor, &ur, &pbig, &ubig, &rbig1, &rbig2,
		&dl1, &dl2, &dp1, &dp2, &udot, &pf, &uf, &rhof);
	double T(PhysLaws::IdealGasLaw::CalcTemp(pf, rhof, prop.R));
	double Hf(prop.Cp * T + uf * uf / 2.0);
	return Array{ rhof * uf, rhof * uf * uf + pf, rhof * uf * Hf };
}

Array FluxProcessor::GodunovFirstOrder(Euler& eq, uint i)
{
	Array deltaF(3, 0);
	double pl, pr, ul, ur, rhol, rhor;

	pl = eq.vars.p[i]; pr = eq.vars.p[i + 1];
	ul = eq.vars.u[i]; ur = eq.vars.u[i + 1];
	rhol = eq.vars.rho[i]; rhor = eq.vars.rho[i + 1];
	deltaF += FortranGOGO(pl, pr, ul, ur, rhol, rhor, eq.prop);

	pl = eq.vars.p[i - 1]; pr = eq.vars.p[i];
	ul = eq.vars.u[i - 1]; ur = eq.vars.u[i];
	rhol = eq.vars.rho[i - 1]; rhor = eq.vars.rho[i];
	deltaF -= FortranGOGO(pl, pr, ul, ur, rhol, rhor, eq.prop);

	return deltaF;
}

double FluxProcessor::gLimiter(const Array& arr, uint currentIndex)
{
	int i = currentIndex;
	double denominator(arr[i] - arr[i - 1]);
	if (denominator == 0)
		return 1;
	double r = (arr[i + 1] - arr[i]) / denominator;
	if (r < 0)
		return 0;
	else
		return (r * r + r) / (r * r + 1);
}

Array FluxProcessor::GodunovSecondOrder(Euler& eq, uint i)
{
	using namespace std;
	using PairDbl = pair<double&, double&>;
	Array deltaF(3, 0);
	double pl, pr, ul, ur, rhol, rhor;
	vector<PairDbl> quantities{ { pl, pr }, { ul, ur }, { rhol, rhor } };
	vector<vector<double>> vars{ eq.vars.p, eq.vars.u, eq.vars.rho };
	auto limiterCoeff{ [](vector<double> vec, uint j) {
		return FluxProcessor::gLimiter(vec, j); } };

	//Forward flux
	for (size_t k = 0; k < 3; ++k)
	{
		if (i != vars[k].size() - 2)
		{
			quantities[k].first = vars[k][i] + 1. / 2. * limiterCoeff(vars[k], i) * (vars[k][i] - vars[k][i - 1]);
			quantities[k].second = vars[k][i + 1] - 1. / 2. * limiterCoeff(vars[k], i + 1) * (vars[k][i + 1] - vars[k][i]);
		}
		else
		{
			quantities[k].first = vars[k][i];
			quantities[k].second = vars[k][i + 1];
		}
	}
	deltaF += FortranGOGO(pl, pr, ul, ur, rhol, rhor, eq.prop);

	//Backward flux
	for (size_t k = 0; k < 3; ++k)
	{
		if (i != 1)
		{
			quantities[k].first = vars[k][i - 1] + 1. / 2. * limiterCoeff(vars[k], i - 1) * (vars[k][i - 1] - vars[k][i - 2]);
			quantities[k].second = vars[k][i] - 1. / 2. * limiterCoeff(vars[k], i) * (vars[k][i] - vars[k][i - 1]);
		}
		else
		{
			quantities[k].first = vars[k][i - 1];
			quantities[k].second = vars[k][i];
		}
	}
	deltaF -= FortranGOGO(pl, pr, ul, ur, rhol, rhor, eq.prop);

	return deltaF;
}

Array FluxProcessor::Roe(Euler& eq, uint i)
{
	Array deltaF(3, 0);
	Funct* rho = new RoeFunct::Rho(eq.vars);
	Funct* u = new RoeFunct::U(eq.vars);
	Funct* H = new RoeFunct::H(eq.vars, eq.enrg_fnctrs->H);
	Funct* p = new RoeFunct::P(rho, H, u, eq.prop);
	Funct* c = new RoeFunct::C(H, u, eq.prop);
	RoeFunct::Dissipation D(eq.vars, *u, *c, *rho, *H, *p);
	Array unitSec{ 0, 1.0, 0 };
	auto& flux(*(eq.F));

	auto rightF = 0.5 * ((flux[i] + eq.vars.p[i] * unitSec) +
		(flux[i + 1] + eq.vars.p[i + 1] * unitSec)) - 0.5 * D(i);
	auto leftF = 0.5 * ((flux[i - 1] + eq.vars.p[i - 1] * unitSec) +
		(flux[i] + eq.vars.p[i] * unitSec)) - 0.5 * D(i - 1);
	deltaF = rightF - leftF;
	delete rho, u, H, c, p;
	return deltaF;
}

double FluxProcessor::GetDeltaForRoe(const std::vector<double>& vec, uint i,
	const RoeFunct::Direction& dir)
{
	auto limiterCoeff{ [](std::vector<double> vec, uint j) {
		return FluxProcessor::gLimiter(vec, j); } };
	double result(0);
	switch (dir) 
	{
	case RoeFunct::Direction::Forward:
		if (i != vec.size() - 2)
		{
			result -= vec[i] + 1. / 2. * limiterCoeff(vec, i) * (vec[i] - vec[i - 1]);
			result += vec[i + 1] - 1. / 2. * limiterCoeff(vec, i + 1) * (vec[i + 1] - vec[i]);
		}
		else
		{
			result -= vec[i];
			result += vec[i + 1];
		}
		break;
	case RoeFunct::Direction::Backward:
		if (i != 1)
		{
			result -= vec[i - 1] + 1. / 2. * limiterCoeff(vec, i - 1) * (vec[i - 1] - vec[i - 2]);
			result += vec[i] - 1. / 2. * limiterCoeff(vec, i) * (vec[i] - vec[i - 1]);
		}
		else
		{
			result -= vec[i - 1];
			result += vec[i];
		}
		break;
	}
	return result;
}