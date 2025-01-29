#ifndef SCHEMES_H
#define SCHEMES_H
#include <vector>
using Array = std::vector<double>;
using uint = unsigned int;

class Euler;
struct Properties;
namespace RoeFunct { enum class Direction; }

enum class SchemeType
{
	FOU = 0,
	GodunovFirstOrder,
	GodunovSecondOrder,
	Roe,
	StegerWorming
};

enum class FaceLocation
{
	Front,
	Back
};

struct FaceValues
{
	double pl = 0, pr = 0, ul = 0, ur = 0, rhol = 0, rhor = 0;
	FaceValues() {};
};

enum class Order { Low = 1, High };

struct FluxProcessor
{
private:
	static bool CheckConvFluxDirection(Euler&, uint);
	static Array FortranGOGO(double, double, double,
		double, double, double, Properties&);
	static double gLimiter(const Array&, uint);
	static auto CalculateFaceQuantities(Euler&, uint, const FaceLocation&);
public:
	static Array FirstOrderUpdwind(Euler&, uint);
	static Array GodunovFirstOrder(Euler&, uint);
	static Array GodunovSecondOrder(Euler&, uint);
	static Array Roe(Euler&, uint);
	static Array StegerWorming(Euler&, uint);

};

#endif // !SCHEMES_H

