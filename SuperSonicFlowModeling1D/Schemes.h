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
	Roe
};

struct FluxProcessor
{
private:
	static bool CheckConvFluxDirection(Euler&, uint);
	static Array FortranGOGO(double, double, double,
		double, double, double, Properties&);
	static double gLimiter(const Array&, uint);
public:
	static Array FirstOrderUpdwind(Euler&, uint);
	static Array GodunovFirstOrder(Euler&, uint);
	static Array GodunovSecondOrder(Euler&, uint);
	static Array Roe(Euler&, uint);
	static double GetDeltaForRoe(const std::vector<double>&, uint, const RoeFunct::Direction&);

};

#endif // !SCHEMES_H

