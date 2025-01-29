#ifndef FUNCT_H
#define FUNCT_H
#include "EquationBase.h"
#include "UnknownsProperties.h"
#include "MathOperations.h"
struct EnergyFunctors;
struct FaceValues;
enum class Order;

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
	double Calculate(double rho, double u)
	{
		return rho * u;
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
		return vars.rho[i] * pow(vars.u[i], 2) + vars.p[i];
	}
	double Calculate(double rho, double u, double p)
	{
		return rho * u * u + p;
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
	double Calculate(double, double, double);
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
	class InterfaceVar : public Funct
	{
	public:
		InterfaceVar() {};
		virtual double operator()(uint) = 0;
		virtual double Calculate(const FaceValues&) = 0;
	};
	class Rho : public InterfaceVar
	{
		Unknowns & vars;
	public:
		Rho(Unknowns& vars) : vars(vars) {};
		double operator()(uint i) override
		{
			return std::sqrt(vars.rho[i] * vars.rho[i + 1]);
		}
		double Calculate(const FaceValues& f) override;
		~Rho() {};
	};
	class U : public InterfaceVar
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
		double Calculate(const FaceValues& f) override;
		~U() {};
	};
	class H : public InterfaceVar
	{
		Unknowns& vars;
		Funct& m_H;
	public:
		H(Unknowns& vars, Funct& H) : vars(vars),
			m_H(H) {};
		double operator()(uint i) override
		{
			return (m_H(i) * sqrt(vars.rho[i]) +
				m_H(i + 1) * sqrt(vars.rho[i + 1])) /
				(sqrt(vars.rho[i]) + sqrt(vars.rho[i + 1]));
		}
		double Calculate(const FaceValues& f) override;
		~H() {};
	};
	class P : public InterfaceVar
	{
		Funct& m_rho;
		Funct& m_H;
		Funct& m_u;
		Properties& m_prop;
	public:
		P(Funct& rho, Funct& H, Funct& u, Properties& prop) : m_rho(rho),
			m_H(H), m_u(u), m_prop(prop) {};
		double operator()(uint i) override
		{
			double g(m_prop.gamma);
			return (g - 1.) / g * m_rho(i) *
				(m_H(i) - 1. / 2. * pow(m_u(i), 2.0));
		}
		double Calculate(const FaceValues& f) override;
		~P() {};
	};
	class C : public InterfaceVar
	{
		Funct& m_H;
		Funct& m_u;
		Properties& m_prop;
	public:
		C(Funct& H, Funct& u, Properties& prop) : m_H(H), m_u(u),
			m_prop(prop) {};
		double operator()(uint i) override
		{
			double g(m_prop.gamma);
			return sqrt((g - 1.) * (m_H(i) -
				1. / 2. * pow(m_u(i), 2)));
		}
		double Calculate(const FaceValues& f) override;
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
		Array operator()(const FaceValues&);
		~Dissipation() {};
	private:
		double del_p = 0, del_u = 0, del_rho = 0;
		double H_ = 0, rho_ = 0, u_ = 0, c_ = 0;
		double coeff1 = 0, coeff2 = 0, coeff3 = 0;
		Array arr1 = { 0, 0, 0 }, arr2 = { 0, 0, 0 },
			arr3 = { 0, 0, 0 };
		Array Calculate();
	};
}

namespace StegerWorming
{
	enum class FaceSide
	{
		Left,
		Right
	};

	class FunctArray
	{
	public:
		FunctArray() {};
		virtual Array operator()(uint) = 0;
		virtual Array operator()(FaceValues&, const FaceSide&) = 0;
		virtual ~FunctArray() {};
	};

	class EigenValues : public FunctArray
	{
	private:
		Unknowns& vars;
		Funct& c;
	public:
		EigenValues(Unknowns& vars, Funct& soundSpeed) :
			vars(vars), c(soundSpeed) {};
		Array operator()(uint i) override
		{
			return Array{ vars.u[i], vars.u[i] + c(i), vars.u[i] - c(i) };
		}
		Array operator()(FaceValues&, const FaceSide&) override;
	};

	class EigenValuesPositive : public FunctArray
	{
	private:
		FunctArray& lambda;
	public:
		EigenValuesPositive(FunctArray& lambda) : lambda(lambda) {};
		Array operator()(uint i)
		{
			return (lambda(i) + mod(lambda(i))) * 0.5;
		}
		Array operator()(FaceValues&, const FaceSide&) override;
	};

	class EigenValuesNegative : public FunctArray
	{
	private:
		FunctArray& lambda;
	public:
		EigenValuesNegative(FunctArray& lambda) : lambda(lambda) {};
		Array operator()(uint i)
		{
			return (lambda(i) - mod(lambda(i))) * 0.5;
		}
		Array operator()(FaceValues&, const FaceSide&) override;
	};
	class Flux : public FunctArray
	{
	private:
		Unknowns& vars;
		Properties& prop;
		FunctArray& lambda;
		Funct& soundSpeed;
	public:
		Flux(Unknowns& vars, Properties& prop, FunctArray& lambda, Funct& c) :
			vars(vars), prop(prop), lambda(lambda), soundSpeed(c) {};
		Array Calculate()
		{
			coeff1 = 2. * (g - 1.) * l1 + l2 + l3;
			coeff2 = 2. * (g - 1.) * l1 * u + l2 * (u + c) + l3 * (u - c);
			coeff3 = (g - 1.) * l1 * u * u + 0.5 * l2 * pow(u + c, 2) +
				0.5 * l3 * (u - c, 2) + (3. - g) * c * c / (2. * (g - 1)) * (l2 + l3);
			return rho / (2. * g) * Array { coeff1, coeff2, coeff3 };
		}
		Array operator()(uint i)
		{
			l1 = lambda(i)[0]; l2 = lambda(i)[1]; l3 = lambda(i)[2];
			g = prop.gamma; rho = vars.rho[i]; c = soundSpeed(i); u = vars.u[i];
			return this->Calculate();
		}
		Array operator()(FaceValues&, const FaceSide&);
	private:
		double l1 = 0, l2 = 0, l3 = 0,
			g = 0, rho = 0, c = 0, u = 0,
			coeff1 = 0, coeff2 = 0, coeff3 = 0;
	};
}

#endif // !FUNCT_H
