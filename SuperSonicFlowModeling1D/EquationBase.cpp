#include "EquationBase.h"

std::vector<double> Equation::gMesh()
{
	return x;
}

std::vector<double> Equation::gSolution()
{
	return U;
}

void Equation::Read(const std::string& inputFile)
{
	std::ifstream input(inputFile);
	std::string var, value;
	while (input >> var >> value)
	{
		parameters[var] = value;
	}
}

void Equation::SetPeriodicBoundaryCondition()
{
	int numOfCells = std::stoi(parameters["NX"]);
	U[numOfCells + 1] = U[1];
	U[0] = U[numOfCells];
}

void Equation::SetBoundaryCondition()
{
	if (parameters["BCType"] == "Periodic")
	{
		return SetPeriodicBoundaryCondition();
	}
}
void Equation::Initialize()
{
	if (parameters["Init"] == "sin")
		return InitializeSinField();
	else if (parameters["Init"] == "forward_step")
		return InitializeForwardStepField();
	else if (parameters["Init"] == "gauss")
		return InitializeGaussField();
}
void Equation::InitializeSinField()
{
	int numOfCells = std::stoi(parameters["NX"]);
	double L = std::stod(parameters["L"]);
	double dx = L / numOfCells;
	for (auto i = 0; i < numOfCells + 2; ++i)
	{
		x.push_back(i * dx - dx / 2);
		U.push_back(sin(6 * PI * x[i]) + 3.0);
	}
}
void Equation::InitializeForwardStepField()
{
	int numOfCells = std::stoi(parameters["NX"]);
	double L = std::stod(parameters["L"]);
	double dx = L / numOfCells;
	for (auto i = 0; i < (numOfCells + 2) / 2; ++i)
	{
		x.push_back(i * dx - dx / 2);
		U.push_back(std::stod(parameters["Ul"]));
	}
	for (auto i = (numOfCells + 2) / 2; i < numOfCells + 2; ++i)
	{
		x.push_back(i * dx - dx / 2);
		U.push_back(std::stod(parameters["Ur"]));
	}
}

void Equation::InitializeGaussField()
{
	int numOfCells = std::stoi(parameters["NX"]);
	double L = std::stod(parameters["L"]);
	double a(std::stod(parameters["a"]));
	double dx = L / numOfCells;
	double CFL(std::stod(parameters["CFL"]));
	double dt(CFL * dx / a);
	for (auto i = 0; i < numOfCells + 2; ++i)
	{
		x.push_back(i * dx - dx / 2);
		U.push_back(exp(-pow((x[i] - 0.5 * L), 2) / 0.02));
	}
}
void Equation::MakeOutput(const std::string& outFile)
{
	std::ofstream outStream(outFile);
	for (auto i = 1; i < U.size() - 1; ++i)
	{
		outStream << x[i] << " " << U[i] << '\n';
	}
}

double Equation::gLimiter(int currentIndex)
{
	int i = currentIndex;
	double denominator(U[i] - U[i - 1]);
	if (denominator == 0)
		return 1;
	double r = (U[i + 1] - U[i]) / denominator;
	if (r < 0)
		return 0;
	else
		return (r * r + r) / (r * r + 1);
}