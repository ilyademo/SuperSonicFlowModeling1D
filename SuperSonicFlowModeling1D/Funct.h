#ifndef FUNCT_H
#define FUNCT_H
#include "EquationBase.h"
#include "UnknownsProperties.h"
#include "MathOperations.h"
struct EnergyFunctors;

class Funct
{
public:
	Funct() {};
	virtual double operator()(uint) = 0;
	virtual double& operator[](uint) {double a; return a; };
	virtual ~Funct() {};
	virtual void RecalcData(uint i) {};
};

class FData : public Funct
{
protected:
	Array data;
public:
	FData() {};
	FData(const Unknowns& vars) { data.resize(vars.u.size(), 0); };
	void RecalcData(uint i)
	{
		data[i] = this->operator()(i);
	}
	double& operator[](uint i) final
	{
		return data[i];
	}
	virtual double operator()(uint) = 0;
	virtual ~FData() {};
};

class W1Funct : public FData
{
	Unknowns& vars;
public:
	W1Funct(Unknowns& vars) : FData(vars), vars(vars)
	{};
	double operator()(uint i) override
	{
		return vars.rho[i];
	}
	~W1Funct() {};
};

class W2Funct : public FData
{
	Unknowns& vars;
public:
	W2Funct(Unknowns& vars) : FData(vars), vars(vars)
	{};
	double operator()(uint i) override
	{
		return vars.rho[i] * vars.u[i];
	}
	~W2Funct() {};
};

class W3Funct : public FData
{
	Funct* E;
	Unknowns& vars;
public:
	W3Funct(Unknowns& vars, Funct* E) : FData(vars), 
		vars(vars), E(E) {};
	double operator()(uint i) override
	{
		return vars.rho[i] * (*E)(i);
	}
	~W3Funct() {};
};

class F1Funct : public Funct
{
	Unknowns& vars;
public:
	F1Funct(Unknowns& vars) : vars(vars) {};
	double operator()(uint i) override
	{
		return vars.rho[i] * vars.u[i];
	}
	~F1Funct() {};
};

class F2Funct : public Funct
{
	Unknowns& vars;
public:
	F2Funct(Unknowns& vars) : vars(vars) {};
	double operator()(uint i) override
	{
		return vars.rho[i] * pow(vars.u[i], 2);
	}
	~F2Funct() {};
};

class F3Funct : public Funct
{
	Unknowns& vars;
	Funct* H;
public:
	F3Funct(Unknowns& vars, Funct* H) : vars(vars),
		H(H) {};
	double operator()(uint i) override
	{
		return vars.rho[i] * vars.u[i] * (*H)(i);
	}
	~F3Funct() {};
};

class EulerQuantitiesFactory
{
public:
	virtual std::vector<Funct*> create(Unknowns&, Funct*) = 0;
};

class ConservativeVariablesFactory : public EulerQuantitiesFactory
{
public:
	std::vector<Funct*> create(Unknowns&, Funct*) override; 
};

class InviscidFluxesFactory : public EulerQuantitiesFactory
{
public:
	std::vector<Funct*> create(Unknowns&, Funct*) override;
};

namespace RoeFunct
{
	enum class Direction{Forward, Backward};
	class Rho : public Funct
	{
		Unknowns & vars;
	public:
		Rho(Unknowns& vars) : vars(vars) {};
		double operator()(uint i) override
		{
			return std::sqrt(vars.rho[i] * vars.rho[i + 1]);
		}
		~Rho() {};
	};
	class U : public Funct
	{
		Unknowns& vars;
	public:
		U(Unknowns& vars) : vars(vars) {};
		double operator()(uint i) override
		{
			return (vars.u[i] * sqrt(vars.rho[i]) +
				vars.u[i + 1] * sqrt(vars.rho[i + 1])) /
				(sqrt(vars.rho[i]) + sqrt(vars.rho[i + 1]));
		}
		~U() {};
	};
	class H : public Funct
	{
		Unknowns& vars;
		Funct* m_H;
	public:
		H(Unknowns& vars, Funct* H) : vars(vars),
			m_H(H) {};
		double operator()(uint i) override
		{
			return ((*m_H)(i) * sqrt(vars.rho[i]) +
				(*m_H)(i + 1) * sqrt(vars.rho[i + 1])) /
				(sqrt(vars.rho[i]) + sqrt(vars.rho[i + 1]));
		}
		~H() {};
	};
	class P : public Funct
	{
		Funct* m_rho;
		Funct* m_H;
		Funct* m_u;
		Properties& m_prop;
	public:
		P(Funct* rho, Funct* H, Funct* u, Properties& prop) : m_rho(rho),
			m_H(H), m_u(u), m_prop(prop) {};
		double operator()(uint i) override
		{
			double g(m_prop.gamma);
			return (g - 1.) / g * (*m_rho)(i) *
				((*m_H)(i) - 1. / 2. * pow((*m_u)(i), 2.0));
		}
		~P() {};
	};
	class C : public Funct
	{
		Funct* m_H;
		Funct* m_u;
		Properties& m_prop;
	public:
		C(Funct* H, Funct* u, Properties& prop) : m_H(H), m_u(u),
			m_prop(prop) {};
		double operator()(uint i) override
		{
			double g(m_prop.gamma);
			return sqrt((g - 1.) * ((*m_H)(i) -
				1. / 2. * pow((*m_u)(i), 2)));
		}
		~C() {};
	};
	class Dissipation
	{
		Unknowns& vars;
		Funct& u;
		Funct& c;
		Funct& rho;
		Funct& H;
		Funct& p;
	public:
		Dissipation(Unknowns& vars, Funct& u, Funct& c,
			Funct& rho, Funct& H, Funct& p) : vars(vars), u(u), c(c),
			rho(rho), H(H), p(p) {};
		Array operator() (uint);
		~Dissipation() {};
	};
}

#endif // !FUNCT_H
