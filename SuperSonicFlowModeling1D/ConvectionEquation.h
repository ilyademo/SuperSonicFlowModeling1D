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

	std::vector<double> GaussAnalytic();
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

	void StepAnalyticOutput();

	~BurgersEquation() {};
};