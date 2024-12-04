#include "ConvectionEquation.h"

void ConvectionEquation::Solve()
{
	int numOfTimeSteps(std::stoi(parameters["NT"]));
	int numOfCells(std::stoi(parameters["NX"]));
	double CFL(std::stod(parameters["CFL"]));
	double convCoeffcient(std::stod(parameters["a"]));
	std::vector<double> U_old(U);
	double dx(std::stod(parameters["L"]) / numOfCells);
	double dt(CFL * dx / convCoeffcient);
	const uint schemeOrder(std::stoi(parameters["schemeOrder"]));
	const std::string TVD_scheme(parameters["TVD"]);
	bool useTVD(TVD_scheme == "van_albada");
	double limiter(1.0);
	double U_face_forward, U_face_backward;
	for (int i = 1; i <= numOfTimeSteps; ++i)
	{
		for (int j = 1; j < numOfCells + 1; ++j)
		{
			if (useTVD)
				limiter = gLimiter(j);
			if (convCoeffcient >= 0)
			{
				switch (schemeOrder)
				{
				case 1:
					U[j] = U_old[j] - (convCoeffcient * dt / dx)
						* (U_old[j] - U_old[j - 1]);
					break;
				case 2:
					U_face_forward = U_old[j] + 1. / 2. * limiter * (U_old[j] - U_old[j - 1]);
					if (j != 1)
						U_face_backward = U_old[j - 1] + 1. / 2. * limiter * (U_old[j - 1] - U_old[j - 2]);
					else
						U_face_backward = U_old[j - 1] + 1. / 2. * limiter * (U_old[j - 1] - U_old[numOfCells - 1]);
					U[j] = U_old[j] - convCoeffcient * dt / dx * (U_face_forward - U_face_backward);
					break;
				default:
					throw "Unrecognized scheme order";
					break;
				}
			}
			else
			{
				switch (schemeOrder)
				{
				case 1:
					U[j] = U_old[j] - (convCoeffcient * dt / dx)
						* (U_old[j + 1] - U_old[j]);
					break;
				case 2:
					U_face_forward = U_old[j] + 1. / 2. * limiter * (U_old[j] - U_old[j + 1]);
					U_face_backward = U_old[j - 1] + 1. / 2. * limiter * (U_old[j - 1] - U_old[j]);
					U[j] = U_old[j] - convCoeffcient * dt / dx * (U_face_forward - U_face_backward);
					break;
				default:
					throw "Unrecognized scheme order";
					break;
				}
			}
		}
		SetBoundaryCondition();
		U_old = U;
	}
}

void BurgersEquation::Solve()
{
	int numOfTimeSteps(std::stoi(parameters["NT"]));
	int numOfCells(std::stoi(parameters["NX"]));
	double CFL(std::stod(parameters["CFL"]));
	double convCoeffcient(std::stod(parameters["a"]));
	std::vector<double> U_old(U);
	double dx(std::stod(parameters["L"]) / numOfCells);
	double dt(CFL * dx / convCoeffcient);
	const uint schemeOrder(std::stoi(parameters["schemeOrder"]));
	double U_face_forward, U_face_backward;
	for (int i = 1; i <= numOfTimeSteps; ++i)
	{
		for (int j = 1; j < numOfCells + 1; ++j)
		{
			if (U[j] >= 0)
			{
				switch (schemeOrder)
				{
				case 1:
					U[j] = U_old[j] - dt / (2.0 * dx) *
						(std::pow(U_old[j], 2) - std::pow(U_old[j - 1], 2));
					break;
				case 2:
					U_face_forward = U_old[j] + 1. / 2. * (U_old[j] - U_old[j - 1]);
					if (j != 1)
						U_face_backward = U_old[j - 1] + 1. / 2. * (U_old[j - 1] - U_old[j - 2]);
					else
						U_face_backward = U_old[j - 1] + 1. / 2. * (U_old[j - 1] - U_old[numOfCells - 1]);
					U[j] = U_old[j] - dt / (2.0 * dx) *
						(std::pow(U_face_forward, 2) - std::pow(U_face_backward, 2));
					break;
				default:
					throw "Unrecognized scheme order";
					break;
				}
			}
			else
			{
				U[j] = U_old[j] - dt / (2.0 * dx) *
					(std::pow(U_old[j + 1], 2) - std::pow(U_old[j], 2));
			}
		}
		SetBoundaryCondition();
		U_old = U;
	}
}

void BurgersEquation::StepAnalyticOutput()
{
	std::vector<double> u_an(U.size() - 2);
	double Ul(std::stod(parameters["Ul"]));
	double Ur(std::stod(parameters["Ur"]));
	int numOfTimeSteps(std::stoi(parameters["NT"]));
	int numOfCells(std::stoi(parameters["NX"]));
	double CFL(std::stod(parameters["CFL"]));
	double L(std::stod(parameters["L"]));
	double dx(L / numOfCells);
	double convCoeffcient(std::stod(parameters["a"]));;
	double dt(CFL * dx / convCoeffcient);
	double S = (Ul + Ur) / 2.0;
	double T = numOfTimeSteps * dt;
	if (Ul > Ur)
	{
		for (int i = 0; i < u_an.size(); ++i)
		{
			(x[i + 1] - 0.5 * L) / T <= S ?
				u_an[i] = Ul :
				u_an[i] = Ur;
		}
	}
	else
	{
		for (int i = 0; i < u_an.size(); ++i)
		{
			double tmp((x[i + 1] - 0.5 * L) / T);
			if (tmp < Ul)
				u_an[i] = Ul;
			else if ((tmp >= Ul) && (tmp <= Ur))
				u_an[i] = tmp;
			else if (tmp > Ur)
				u_an[i] = Ur;
		}
	}

	std::ofstream outStream("OutputFiles/BurgersAnalyticStep.plt");
	for (auto i = 1; i < U.size() - 1; ++i)
	{
		outStream << x[i] << " " << u_an[i - 1] << '\n';
	}
}


std::vector<double> ConvectionEquation::GaussAnalytic()
{
	int numOfCells = std::stoi(parameters["NX"]);
	double L = std::stod(parameters["L"]);
	double a(std::stod(parameters["a"]));
	double dx = L / numOfCells;
	double CFL(std::stod(parameters["CFL"]));
	int NT(std::stoi(parameters["NT"]));
	double dt(CFL * dx / a);

	std::vector<double> result;
	for (auto i = 1; i < numOfCells + 1; ++i)
	{
		result.push_back(std::exp(-pow((x[i] - 0.5 * L - a * NT * dt), 2) / 0.02));
	}
	return result;
}