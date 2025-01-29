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
	//do not use this, fix fluxes before
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

auto FluxProcessor::CalculateFaceQuantities(Euler& eq, uint i, const FaceLocation& dir)
{
	using namespace std;
	using PairDbl = pair<double&, double&>;
	Array deltaF(3, 0);
	FaceValues f;
	vector<PairDbl> quantities{ { f.pl, f.pr }, { f.ul, f.ur }, { f.rhol, f.rhor } };
	vector<vector<double>> vars{ eq.vars.p, eq.vars.u, eq.vars.rho };
	auto limiterCoeff{ [](vector<double> vec, uint j) {
		return FluxProcessor::gLimiter(vec, j); } };
	switch (dir)
	{
	case FaceLocation::Front:
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
		break;

	case FaceLocation::Back:
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
		break;
	}
	return f;
}

Array FluxProcessor::GodunovSecondOrder(Euler& eq, uint i)
{
	auto f = FluxProcessor::CalculateFaceQuantities(eq, i, FaceLocation::Front);
	auto rightF = FortranGOGO(f.pl, f.pr, f.ul, f.ur, f.rhol, f.rhor, eq.prop);

	f = FluxProcessor::CalculateFaceQuantities(eq, i, FaceLocation::Back);
	auto leftF = FortranGOGO(f.pl, f.pr, f.ul, f.ur, f.rhol, f.rhor, eq.prop);

	return rightF - leftF;
}

Array FluxProcessor::Roe(Euler& eq, uint i)
{
	FaceValues f;
	Array rightF(3, 0), leftF(3, 0);
	std::unique_ptr<Funct> rho = std::make_unique<RoeFunct::Rho>(eq.vars);
	std::unique_ptr<Funct> u = std::make_unique<RoeFunct::U>(eq.vars);
	std::unique_ptr<Funct> H = std::make_unique<RoeFunct::H>(eq.vars, *eq.enrg_fnctrs->H);
	std::unique_ptr<Funct> p = std::make_unique<RoeFunct::P>(*rho, *H, *u, eq.prop);
	std::unique_ptr<Funct> c = std::make_unique<RoeFunct::C>(*H, *u, eq.prop);
	RoeFunct::Dissipation D(eq.vars, *u, *c, *rho, *H, *p);
	auto schemeOrder(static_cast<Order>(std::stoi(eq.parameters["RoeOrder"])));
	auto& F(*(eq.F));
	switch (schemeOrder)
	{
	case Order::High:
		f = FluxProcessor::CalculateFaceQuantities(eq, i, FaceLocation::Front);
		rightF = 0.5 * F.GetFluxByValuesOnFace(f) - 0.5 * D(f);
		f = FluxProcessor::CalculateFaceQuantities(eq, i, FaceLocation::Back);
		leftF = 0.5 * F.GetFluxByValuesOnFace(f) - 0.5 * D(f);
		break;
	case Order::Low:
		rightF = 0.5 * (F[i] + F[i + 1]) - 0.5 * D(i);
		leftF = 0.5 * (F[i - 1] + F[i]) - 0.5 * D(i - 1);
		break;
	default:
		throw "Unrecognized Roe scheme order";
	}	
	return rightF - leftF;
}

Array FluxProcessor::StegerWorming(Euler& eq, uint i)
{
	using namespace StegerWorming;
	FaceValues f;
	Array rightF(3, 0), leftF(3, 0);
	std::unique_ptr<Funct> c = std::make_unique<PhysLaws::SoundSpeed>(eq.vars, eq.prop);
	std::unique_ptr<FunctArray> eigenValues = std::make_unique<EigenValues>(eq.vars, *c);
	std::unique_ptr<FunctArray> lambdaPlus = std::make_unique<EigenValuesPositive>(*eigenValues);
	std::unique_ptr<FunctArray> lambdaMinus = std::make_unique<EigenValuesNegative>(*eigenValues);
	std::unique_ptr<FunctArray> fluxPlus = std::make_unique<Flux>(eq.vars, eq.prop, *lambdaPlus, *c);
	std::unique_ptr<FunctArray> fluxMinus = std::make_unique<Flux>(eq.vars, eq.prop, *lambdaMinus, *c);
	auto& fPlus(*fluxPlus);
	auto& fMinus(*fluxMinus);

	auto schemeOrder(static_cast<Order>(std::stoi(eq.parameters["StegerWormingOrder"])));
	switch (schemeOrder)
	{
	case Order::High:
		f = FluxProcessor::CalculateFaceQuantities(eq, i, FaceLocation::Front);
		rightF = fPlus(f, FaceSide::Left) + fMinus(f, FaceSide::Right);
		f = FluxProcessor::CalculateFaceQuantities(eq, i, FaceLocation::Back);
		leftF = fPlus(f, FaceSide::Left) + fMinus(f, FaceSide::Right);
		break;
	case Order::Low:
		rightF = fPlus(i) + fMinus(i + 1);
		leftF = fPlus(i - 1) + fMinus(i);
		break;
	default:
		throw "Unrecognized StegerWorming scheme order";
	}
	return rightF - leftF;
}