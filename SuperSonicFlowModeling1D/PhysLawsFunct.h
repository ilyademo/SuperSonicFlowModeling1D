#ifndef PHYS_LAWS_PUNCT
#define PHYS_LAWS_PUNCT
#include "UnknownsProperties.h"
#include "Funct.h"

namespace PhysLaws
{
	class TotalEnergy : public Funct
	{
	private:
		Unknowns& vars;
		Properties& prop;
	public:
		TotalEnergy(Unknowns& vars, Properties& prop) :
			vars(vars), prop(prop) {};
		double operator()(uint i)
		{
			auto p(vars.p[i]);
			auto rho(vars.rho[i]);
			auto e(1. / (prop.gamma - 1) * p / rho);
			return e + pow(vars.u[i], 2) / 2.;
		}
	};

	class IdealGasLaw : public Funct
	{
	private:
		Unknowns& vars;
		Properties& prop;
	public:
		IdealGasLaw(Unknowns& vars, Properties& prop) : vars(vars),
			prop(prop) {};

		double operator()(uint i)
		{
			return  vars.p[i] / (vars.rho[i] * prop.R);
		}
		static double CalcTemp(double pf, double rhof, double R)
		{
			return pf / (rhof * R);
		}
	};

	class TotalEnthalpy : public Funct
	{
	private:
		Unknowns& vars;
		Properties& prop;
	public:
		TotalEnthalpy(Unknowns& vars, Properties& prop) :
			vars(vars), prop(prop) {};
		double operator()(uint i)
		{
			return prop.Cp * (*prop.T)(i) + pow(vars.u[i], 2) / 2.;
		}
		double Calculate(double rho, double u, double p)
		{
			return prop.Cp * IdealGasLaw::CalcTemp(p, rho, prop.R) + u * u / 2;
		}
	};

	class SoundSpeed : public Funct
	{
	private:
		Unknowns& vars;
		Properties& prop;
	public:
		SoundSpeed(Unknowns& vars, Properties& prop) :
			vars(vars), prop(prop) {};
		double operator()(uint i)
		{
			return sqrt(prop.gamma * vars.p[i] / vars.rho[i]);
		}
		double Calculate(const double& rho, const double& p)
		{
			return sqrt(prop.gamma * p / rho);
		}
	};
}
struct EnergyFunctors
{
	Funct* E;
	Funct* H;
	EnergyFunctors() : E(nullptr), H(nullptr) {};
	EnergyFunctors(Unknowns& vars, Properties& prop) :
		E(new PhysLaws::TotalEnergy(vars, prop)),
		H(new PhysLaws::TotalEnthalpy(vars, prop)) {};
};
#endif // !PHYS_LAWS_PUNCT

