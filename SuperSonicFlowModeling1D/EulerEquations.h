#ifndef EULER_EQUATIONS_H
#define EULER_EQUATIONS_H
#include <map>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>

#include "EquationBase.h"
#include "UnknownsProperties.h"

using uint = unsigned int;
using Array = std::vector<double>;

class Funct;
struct EnergyFunctors;

class EulerVariables;

struct Visitor
{
	virtual void print(EulerVariables* ref, int i) = 0;
	virtual ~Visitor() = default;
};

class EulerVariables
{
protected:
	std::vector<Funct*> vars;
public:
	virtual void print(Visitor& visitor, int i) = 0;
	virtual std::string getKindVar() = 0;
	virtual size_t size() = 0;
	EulerVariables(std::vector<Funct*> vars) : vars(vars) {};
	virtual ~EulerVariables() {}
	Array operator[](uint i);
	Funct* operator()(int i) const { return vars[i]; }
};

class ConservativeVariables : public EulerVariables
{
public:
	void print(Visitor& visitor, int i) override {visitor.print(this, i);}
	std::string getKindVar() override { return "W"; }
	size_t size() override { return vars.size(); }
	ConservativeVariables(std::vector<Funct*> W) : EulerVariables(W) {};
	~ConservativeVariables() {}
	void Set(const Array& vec, uint i);
	void Recalc(uint i);
};

class InviscidFluxes : public EulerVariables
{
public:
	void print(Visitor& visitor, int i) override { visitor.print(this, i); }
	std::string getKindVar() override { return "F"; }
	size_t size() override { return vars.size(); }
	InviscidFluxes(std::vector<Funct*> F) : EulerVariables(F) {};
	~InviscidFluxes() {}	
};

class Euler : public Equation
{
private:
	ConservativeVariables* w;
	InviscidFluxes* F;
	EnergyFunctors* enrg_fnctrs;
	Unknowns vars;
	Properties prop;

	friend struct FluxProcessor;
	bool CheckConvFluxDirection(uint);
	Array CalculateDeltaFlux(uint);
	void UpdateTimeStep(const Array&, const double,
		const double, double&);
public:
	Euler(const std::string&);
	void Initialize();
	void SetWallBoundaryCondition();
	void Solve();
	void MakeOutput(const std::string&);
};

class PrintVisitor : public Visitor
{
public:
	void print(EulerVariables* vars, int i) override;
};

#endif // !EULER_EQUATIONS_H