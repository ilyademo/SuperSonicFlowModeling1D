#include "EulerEquations.h"
#include "PhysLawsFunct.h"
#include "Funct.h"
#include "EquationBase.h"
#include "MathOperations.h"
#include "GodunovFunctions.h"
#include "Schemes.h"

Array EulerVariables::operator[](uint i)
{
	return { (*vars[0])(i), (*vars[1])(i), (*vars[2])(i) };
}

void ConservativeVariables::Recalc(uint i)
{
	for (auto var : vars)
		var->RecalcData(i);
}

void ConservativeVariables::Set(const Array& vec, uint i)
{
	for (auto j = 0; j < vars.size(); ++j)
		(*vars[j])[i] = vec[j];
}

Euler::Euler(const std::string& inputFile)
{
	Read(inputFile);
}

void Euler::Initialize()
{
	const int numOfCells = std::stoi(parameters["NX"]);
	const double L = std::stod(parameters["L"]);
	const double diaphr_loc = std::stoi(parameters["d_loc"]);
	const double dx = L / numOfCells;
	const double pl = std::stod(parameters["pl"]);
	const double rhol = std::stod(parameters["rhol"]);
	const double pr = std::stod(parameters["pr"]);
	const double rhor = std::stod(parameters["rhor"]);
	const double u1 = std::stod(parameters["u1"]);
	const double u2 = std::stod(parameters["u2"]);
	vars.u.resize(numOfCells + 2, 0);
	vars.rho.resize(numOfCells + 2, 0);
	vars.p.resize(numOfCells + 2, 0);
	x.resize(numOfCells + 2, 0);
	prop.T = new PhysLaws::IdealGasLaw(vars, prop);
	for (auto i = 1; i < numOfCells + 1; ++i)
	{
		x[i] = i * dx - dx / 2;
		if (i <= diaphr_loc)
		{
			vars.u[i] = u1;
			vars.p[i] = pl;
			vars.rho[i] = rhol; 
		}
		else
		{
			vars.u[i] = u2;
			vars.p[i] = pr;
			vars.rho[i] = rhor;
		}
	}
	prop.gamma = std::stod(parameters["gamma"]);
	prop.Cp = std::stod(parameters["Cp"]);
	prop.R = prop.Cp * (prop.gamma - 1) / prop.gamma;
	prop.Cv = prop.Cp - prop.R;
	enrg_fnctrs = new EnergyFunctors(vars, prop);
	
	EulerQuantitiesFactory* factory = new ConservativeVariablesFactory();
	w = new ConservativeVariables(factory->create(vars, enrg_fnctrs->E));

	factory = new InviscidFluxesFactory();
	F = new InviscidFluxes(factory->create(vars, enrg_fnctrs->H));
	SetWallBoundaryCondition();
}

void Euler::SetWallBoundaryCondition()
{
	vars.u[0] = -vars.u[1];
	vars.u[x.size() - 1] = -vars.u[x.size() - 2];
	vars.rho[0] = vars.rho[1];
	vars.rho[x.size() - 1] = vars.rho[x.size() - 2];
	vars.p[0] = vars.p[1];
	vars.p[x.size() - 1] = vars.p[x.size() - 2];
}

bool Euler::CheckConvFluxDirection(uint i)
{
	return (vars.u[i + 1] + vars.u[i]) / 2. >= 0;
}

Array Euler::CalculateDeltaFlux(uint i)
{
	const auto scheme(std::stoi(parameters["scheme"]));
	switch (scheme)
	{
	case uint(SchemeType::FOU):
		return FluxProcessor::FirstOrderUpdwind(*this, i);
	case uint(SchemeType::GodunovFirstOrder):
		return FluxProcessor::GodunovFirstOrder(*this, i);
	case uint(SchemeType::GodunovSecondOrder):
		return FluxProcessor::GodunovSecondOrder(*this, i);
	case uint(SchemeType::Roe):
		return FluxProcessor::Roe(*this, i);
	default:
		throw "Unrecognized scheme order";
		break;
	}
}

void Euler::UpdateTimeStep(const Array& u,
	const double CFL, const double dx, double& dt)
{
	double dt_local(dt);
	for (auto i = 1; i < u.size() - 1; ++i)
	{
		if (u[i] > 0)
			dt_local = CFL * dx / u[i];
		if (dt_local < dt)
		{
			std::cout << "Update time step from " <<
				dt << " to " << dt_local << '\n';
			dt = dt_local;
		}
	}
}

void Euler::Solve()
{
	int numOfCells = std::stoi(parameters["NX"]);
	const double L = std::stod(parameters["L"]);
	const double dx = L / numOfCells;
	const double time = std::stod(parameters["time"]);
	const double CFL = std::stod(parameters["CFL"]);
	double dt(0.), t(0.);
	dt = 1e-7;
	do
	{
		Visitor* printVisitor = new PrintVisitor();
		for (int j = 1; j < std::stoi(parameters["NX"]) + 1; ++j)
		{
			Array deltaF(CalculateDeltaFlux(j));
			auto w_old = (*w)[j];
			Array w_new(3, 0);
			w_new = w_old - dt / dx * deltaF;
			w->Set(w_new, j);
		}
		vars.RecalculateVariables(*w, prop);
		//if (!vars.u > 0)
		//	UpdateTimeStep(vars.u, CFL, dx, dt);
		
		this->SetWallBoundaryCondition();
		t += dt;
		//std::cout << t << '\n';
	} while (t < time);
	std::cout << "Actual t " << t << ", neseccery time " << time << '\n';
}

void Euler::MakeOutput(const std::string& outFile)
{
	std::ofstream outStream(outFile);
	const std::string tab{ "		" };
	outStream << "Variables = X,Ro,U,P,T" << '\n';
	outStream << "Zone i = " << x.size() - 2 << '\n';
	outStream << std::fixed << std::setprecision(6);
	for (auto i = 1; i < x.size() - 1; ++i)
	{
		outStream << std::setw(16);
		outStream << x[i] << tab << vars.rho[i] << tab <<
			vars.u[i] << tab << vars.p[i] << tab << (*prop.T)(i) << '\n';
	}
}

void PrintVisitor::print(EulerVariables* vars, int i)
{
	std::cout << vars->getKindVar() << " at " << i << ":\n";
	for (int j = 0; j < vars->size(); j++)
		std::cout << (*(*vars)(j))(i) << " ";
	std::cout << '\n';
}