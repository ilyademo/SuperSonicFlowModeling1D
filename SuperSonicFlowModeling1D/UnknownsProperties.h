#ifndef UNKNOWNS_PROPERTIES_H
#define UNKNOWNS_PROPERTIES_H
#include <iostream>
#include <vector>

class Funct;
class ConservativeVariables;
struct Properties;

using Array = std::vector<double>;
struct Unknowns
{
	Unknowns() {};
	Unknowns(const Unknowns& obj) : rho(obj.rho),
		p(obj.p), u(obj.u){};
	Array rho, u, p;
	void RecalculateVariables(const ConservativeVariables&,
		const Properties&);
};

struct Properties
{
	Properties() : T(nullptr), gamma(0), Cp(0) {};
	double gamma;
	double Cp;
	double R;
	double Cv;
	Funct* T;
};

#endif // !UNKNOWNS_PROPERTIES_H

