#pragma once
#include "EquationBase.h"

class ConvectionEquation : public Equation
{
public:
	ConvectionEquation() {};
	ConvectionEquation(const std::string& inputFile) 
	{
		Read(inputFile);
	}
	void Solve();

	std::vector<double> gGaussAnalytic();

	void GaussAnalytic(const std::string&);
	void SinAnalytic(const std::string&);
	void StepAnalytic(const std::string&);
	~ConvectionEquation() {};
};

class BurgersEquation : public Equation
{
public:
	BurgersEquation(const std::string& inputFile)
	{
		Read(inputFile);
	}
	void Solve();

	void StepAnalyticOutput(const std::string&);

	~BurgersEquation() {};
};