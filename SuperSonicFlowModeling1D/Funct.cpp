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

double F3Funct::Calculate(double rho, double u, double p)
{
	auto& Enthalpy = *dynamic_cast<PhysLaws::TotalEnthalpy*>(H);
	return rho * u * Enthalpy.Calculate(rho, u, p);
}


namespace RoeFunct
{
	double Rho::Calculate(const FaceValues& f)
	{
		return sqrt(f.rhol * f.rhor);
	}
	double U::Calculate(const FaceValues& f)
	{
		return (f.ul * sqrt(f.rhol) +
			f.ur * f.rhor) /
			(sqrt(f.rhol) + sqrt(f.rhor));
	}
	double H::Calculate(const FaceValues& f)
	{
		auto& Enthalpy = *dynamic_cast<PhysLaws::TotalEnthalpy*>(&m_H);
		auto Hl(Enthalpy.Calculate(f.rhol, f.ul, f.pl));
		auto Hr(Enthalpy.Calculate(f.rhor, f.ur, f.pr));
		return (Hl * sqrt(f.rhol) +
			Hr * sqrt(f.rhor)) /
			(sqrt(f.rhol) + sqrt(f.rhor));
	}
	double P::Calculate(const FaceValues& f)
	{
		auto& H_ = *dynamic_cast<RoeFunct::H*>(&m_H);
		auto& rho_ = *dynamic_cast<RoeFunct::Rho*>(&m_rho);
		auto& u_ = *dynamic_cast<RoeFunct::U*>(&m_u);
		double g(m_prop.gamma);
		return (g - 1.) / g * rho_.Calculate(f) *
			(H_.Calculate(f) - 1. / 2. * pow(u_.Calculate(f), 2.0));
	}
	double C::Calculate(const FaceValues& f)
	{
		auto& H_ = *dynamic_cast<RoeFunct::H*>(&m_H);
		auto& u_ = *dynamic_cast<RoeFunct::U*>(&m_u);
		double g(m_prop.gamma);
		return sqrt((g - 1.) * (H_.Calculate(f) -
			1. / 2. * pow(u_.Calculate(f), 2)));
	}

	Array Dissipation::Calculate()
	{
		coeff1 = abs(u_ + c_) *
			(del_p + rho_ * c_ * del_u) / (2. * pow(c_, 2.0));
		coeff2 = abs(u_ - c_) *
			(del_p - rho_ * c_ * del_u) / (2. * pow(c_, 2.0));
		coeff3 = abs(u_) * (del_rho - del_p / pow(c_, 2.0));
		arr1 = Array{ 1., u_ + c_, H_ + c_ * u_ };
		arr2 = Array{ 1., u_ - c_, H_ - c_ * u_ };
		arr3 = Array{ 1., u_, pow(u_, 2) / 2. };
		return 
			coeff1 * arr1 +
			coeff2 * arr2 +
			coeff3 * arr3;
	}

	Array Dissipation::operator()(const FaceValues& f)
	{
		del_p = f.pr - f.pl;
		del_u = f.ur - f.ul;
		del_rho = f.rhor - f.rhol;
		H_ = (*dynamic_cast<RoeFunct::H*>(&H)).Calculate(f);
		rho_ = (*dynamic_cast<RoeFunct::Rho*>(&rho)).Calculate(f);
		u_ = (*dynamic_cast<RoeFunct::U*>(&u)).Calculate(f);
		c_ = (*dynamic_cast<RoeFunct::C*>(&c)).Calculate(f);
		return this->Calculate();
	}

	Array Dissipation::operator()(uint i)
	{
		del_p = vars.p[i + 1] - vars.p[i];
		del_u = vars.u[i + 1] - vars.u[i];
		del_rho = vars.rho[i + 1] - vars.rho[i];
		H_ = H(i);
		rho_ = rho(i);
		u_ = u(i);
		c_ = c(i);
		return this->Calculate();
	}
}

namespace StegerWorming
{
	Array EigenValues::operator()(FaceValues& f, const FaceSide& side)
	{
		auto& c_local = (*dynamic_cast<PhysLaws::SoundSpeed*>(&c));
		double u(0.), c(0.);
		switch (side)
		{
		case FaceSide::Left:
			c = c_local.Calculate(f.rhol, f.pl);
			u = f.ul;
			break;
		case FaceSide::Right:
			c = c_local.Calculate(f.rhor, f.pr);
			u = f.ur;
			break;
		default:
			throw "Unrecognized face side";
		}
		return Array{ u, u + c, u - c };
	}

	Array EigenValuesPositive::operator()(FaceValues& f, const FaceSide& side)
	{
		Array lambda_arr(lambda(f, side));
		return (lambda_arr + mod(lambda_arr)) * 0.5;
	}

	Array EigenValuesNegative::operator()(FaceValues& f, const FaceSide& side)
	{
		Array lambda_arr(lambda(f, side));
		return (lambda_arr - mod(lambda_arr)) * 0.5;
	}

	Array Flux::operator()(FaceValues& f, const FaceSide& side)
	{
		Array lambda_arr(lambda(f, side));
		auto& c_local = (*dynamic_cast<PhysLaws::SoundSpeed*>(&soundSpeed));
		l1 = lambda_arr[0]; l2 = lambda_arr[1]; l3 = lambda_arr[2];
		g = prop.gamma;
		switch (side)
		{
		case FaceSide::Left:
			rho = f.rhol;
			c = c_local.Calculate(f.rhol, f.pl);
			u = f.ul;
			break;
		case FaceSide::Right:
			rho = f.rhor;
			c = c_local.Calculate(f.rhor, f.pr);
			u = f.ur;
			break;
		default:
			throw "Unrecognized face side";
		}
		return this->Calculate();
	}
}