#include "ConvectionEquation.h"
#include "EulerEquations.h"
#include "GodunovFunctions.h"

void Output(std::vector<double> x, std::vector<double> u, const std::string& outFile)
{
	std::ofstream outStream(outFile);
	for (auto i = 1; i < u.size() - 1; ++i)
	{
		outStream << x[i] << " " << u[i] << '\n';
	}
}

int main()
{
	try
	{
		//{
		//	std::string inp("InputFiles/inputConvEquation.txt");
		//	std::string out("OutputFiles/outputConvEquation_gauss_2rd.plt");
		//	auto eq_default = new ConvectionEquation(inp);
		//	eq_default->Initialize();
		//	eq_default->Solve();
		//	eq_default->MakeOutput(out);
		//	std::string out_analytical("OutputFiles/outputConvEquation_gauss_analytical.plt");
		//	//eq_default->SinAnalytic(out_analytical);
		//	//eq_default->StepAnalytic(out_analytical);
		//	eq_default->GaussAnalytic(out_analytical);
		//}

		//{
		//	std::string inp("InputFiles/inputBurgersEquation.txt");
		//	std::string out("OutputFiles/outputBurgersEquation_gauss_2rd_TVD.plt");
		//	auto eq_burgers = new BurgersEquation(inp);
		//	eq_burgers->Initialize();
		//	std::string out_start_field("OutputFiles/outputBurgersEquation_gauss_intial.plt");
		//	Output(eq_burgers->gMesh(), eq_burgers->gSolution(), out_start_field);
		//	eq_burgers->Solve();
		//	eq_burgers->MakeOutput(out);
		//	//eq_burgers->StepAnalyticOutput(out_analytical);
		//}

		/*{
			std::string inp("InputFiles/inputGauss.txt");
			std::string out("OutputFiles/outputGauss.plt");
			auto eq_new = new ConvectionEquation(inp);
			eq_new->Initialize();
			eq_new->Solve();
			auto result = eq_new->gSolution();
			auto mesh = eq_new->gMesh();
			eq_new->MakeOutput(out);
			auto analytical_result = eq_new->gGaussAnalytic();
			double max_value{ 0 };
			int max_index{ 0 };
			double x1{ 0 };
			for (int i = 1; i < analytical_result.size() - 1; ++i)
			{
				if (analytical_result[i] > max_value)
				{
					max_index = i;
					max_value = analytical_result[i];
				}
			}
			int x1_index(max_index - mesh.size() / 8);
			x1 = analytical_result[x1_index];
			double err2(abs(max_value - result[max_index]));
			double integral_error(0.);
			double dx(mesh[2] - mesh[1]);
			for (auto i = 1; i < analytical_result.size() - 1; ++i)
			{
				integral_error += abs(result[i] - analytical_result[i])*dx;
			}

			std::cout << integral_error << " " << err2 << '\n';

			Output(mesh, analytical_result, "OutputFiles/AnalyticGauss.plt");
		}*/
		double pbig(0.), ubig(0.), rbig1(0.), rbig2(0.), dl1(0.), dl2(0.),
			dp1(0.), dp2(0.);
		double gamma(1.4);
		//individual
		double pl(0.4), pr(0.4), rhol(1.0), rhor(1.0), ul(-2.0), ur(2.0);
		////osadchi
		//pl = 460.894; pr = 46.0950; rhol = 5.99924; rhor = 5.99242; ul = 19.5975; ur = -6.19633;
		raspad(&gamma, &pl, &rhol, &ul,
			&pr, &rhor, &ur, &pbig, &ubig, &rbig1, &rbig2,
			&dl1, &dl2, &dp1, &dp2);
		std::cout << "ddd";

		/*std::string inp("InputFiles/inputEuler.txt");
		std::string out("OutputFiles/Test5/SW-2.plt");
		std::ofstream outStream(out);
		auto eq_new = new Euler(inp);
		eq_new->Initialize();
		eq_new->Solve();

		eq_new->MakeOutput(out);*/
	}
	catch (const char* err_message)
	{
		std::cout << err_message << '\n';
	}
	std::cout << "End of program!\n";
}
