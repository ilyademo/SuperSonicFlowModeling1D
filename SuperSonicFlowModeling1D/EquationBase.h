#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <map>
#define PI 3.1415926535

using uint = unsigned int;
using Array = std::vector<double>;

class Equation
{
protected:
	std::map<std::string, std::string> parameters;
	std::vector<double> U;
	std::vector<double> x;
	
	void SetPeriodicBoundaryCondition();
	void SetBoundaryCondition();
	
	void InitializeSinField();
	void InitializeForwardStepField();
	void InitializeGaussField();

	double gLimiter(int);
	
	void Read(const std::string& inputFile);
public:
	void MakeOutput(const std::string&);
	
	void Initialize();

	std::vector<double> gSolution();

	std::vector<double> gMesh();
};