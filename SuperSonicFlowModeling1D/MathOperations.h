#ifndef MATH_OPERATIONS_H
#define MATH_OPERATIONS_H
#include <vector>
#include <iostream>
#include "EquationBase.h"

inline Array operator*(const Array& lvalue, double rvalue)
{
	Array result(lvalue);
	for (auto i = 0; i < lvalue.size(); i++)
	{
		result[i] = lvalue[i] * rvalue;
	}
	return result;
}

inline Array operator*(double lvalue, const Array& rvalue)
{
	return operator*(rvalue, lvalue);
}

inline Array& operator*(const Array& lvalue, const Array& rvalue)
{
	if (lvalue.size() != rvalue.size())
		throw "different sizes of arrays in multiply";
	Array result(lvalue);
	for (int i = 0; i < lvalue.size(); ++i)
	{
		result[i] = lvalue[i] * rvalue[i];
	}
	return result;
}

inline Array& operator/(const Array& lvalue, const Array& rvalue)
{
	if (lvalue.size() != rvalue.size())
		throw "different sizes of arrays in division";
	Array result(lvalue);
	for (int i = 0; i < lvalue.size(); ++i)
	{
		result[i] = lvalue[i] / rvalue[i];
	}
	return result;
}

inline Array operator+(const Array& lvalue, const Array& rvalue)
{
	if (lvalue.size() != rvalue.size())
		throw "different sizes of arrays in sum";
	Array result(lvalue);
	for (int i = 0; i < lvalue.size(); ++i)
	{
		result[i] = lvalue[i] + rvalue[i];
	}
	return result;
}

inline Array& operator+=(Array& lvalue, const Array& rvalue)
{
	for (int i = 0; i < lvalue.size(); ++i)
	{
		lvalue[i] = lvalue[i] + rvalue[i];
	}
	return lvalue;
}

inline Array operator+(const Array& lvalue, double rvalue)
{
	Array result(lvalue);
	for (int i = 0; i < lvalue.size(); ++i)
	{
		result[i] = lvalue[i] + rvalue;
	}
	return result;
}
inline Array operator+(double lvalue, const Array& rvalue)
{
	return operator+(rvalue, lvalue);
}
inline Array operator-(const Array& lvalue, const Array& rvalue)
{
	Array result(lvalue);
	for (int i = 0; i < lvalue.size(); ++i)
	{
		result[i] = lvalue[i] - rvalue[i];
	}
	return result;
}
inline Array& operator-=(Array& lvalue, const Array& rvalue)
{
	for (int i = 0; i < lvalue.size(); ++i)
	{
		lvalue[i] = lvalue[i] - rvalue[i];
	}
	return lvalue;
}

inline double operator!(const Array& vec)
{
	double res(0.);
	for (auto i = 0; i < vec.size(); ++i)
	{
		res += vec[i] * vec[i];
	}
	return sqrt(res);
}

inline std::ostream& operator<<(std::ostream& os, const Array& rvalue)
{
	for (auto i : rvalue)
		os << i << " ";
	return os;
}

#endif // !MATH_OPERATIONS_H