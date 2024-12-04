#include "ConvectionEquation.h"
#include "EulerEquations.h"

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
		/*{
			std::string inp("InputFiles/inputConvEquation.txt");
			std::string out("OutputFiles/outputConvEquation.plt");
			auto eq_default = new ConvectionEquation(inp);
			eq_default->Initialize();
			eq_default->Solve();
			eq_default->MakeOutput(out);
		}*/

		/*{
			std::string inp("inputfiles/inputburgersequation.txt");
			std::string out("outputfiles/outputburgersequation.plt");
			auto eq_burgers = new BurgersEquation(inp);
			eq_burgers->Initialize();
			eq_burgers->Solve();
			eq_burgers->MakeOutput(out);
			eq_burgers->StepAnalyticOutput();
		}*/

		//{
		//	std::string inp("InputFiles/inputGauss.txt");
		//	std::string out("OutputFiles/outputGauss.plt");
		//	auto eq_new = new ConvectionEquation(inp);
		//	eq_new->Initialize();
		//	eq_new->Solve();
		//	auto result = eq_new->gSolution();
		//	auto mesh = eq_new->gMesh();
		//	eq_new->MakeOutput(out);
		//	auto analytical_result = eq_new->GaussAnalytic();
		//	double max_value{ 0 };
		//	int max_index{ 0 };
		//	double x1{ 0 };
		//	for (int i = 1; i < analytical_result.size() - 1; ++i)
		//	{
		//		if (analytical_result[i] > max_value)
		//		{
		//			max_index = i;
		//			max_value = analytical_result[i];
		//		}
		//	}
		//	/*int x1_index(max_index - mesh.size() / 8);
		//	x1 = analytical_result[x1_index];*/
		//	//double err1(abs(x1 - result[x1_index]));
		//	double err2(abs(max_value - result[max_index]));

		//	std::cout << /*err1 << " " <<*/ err2 << '\n';

		//	Output(mesh, analytical_result, "OutputFiles/AnalyticGauss.plt");
		//}
		std::string inp("InputFiles/inputEuler.txt");
		std::string out("OutputFiles/outputEulerRow.plt");
		std::ofstream outStream(out);
		auto eq_new = new Euler(inp);
		eq_new->Initialize();
		eq_new->Solve();

		eq_new->MakeOutput(out);
	}
	catch (const char* err_message)
	{
		std::cout << err_message << '\n';
	}
	std::cout << "End of program!\n";
}
