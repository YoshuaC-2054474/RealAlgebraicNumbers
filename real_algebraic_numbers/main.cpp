#include "RealAlgebraicNumber.h"
#include <string>
#include <iostream>
//#include "Polynomial3.h"

int main()
{
	/*auto a = RealAlgebraicNumber({-2,0,1}, {1.4,1.5});
	auto b = RealAlgebraicNumber({-3,0,1}, { 1.7,1.8 });*/
	/*auto a = RealAlgebraicNumber({ -2,1 }, { 1.9,2.1 });
	auto b = RealAlgebraicNumber({ -3,1 }, { 2.9,3.1 });*/
	auto a = RealAlgebraicNumber({ -2,1,2 }, { 0.5,1.0 });
	auto b = RealAlgebraicNumber({ -3,1,3 }, { 0.5,1.0});
	auto c = a + b;
	std::string result = c.toString();
	std::cout << result << std::endl;

	return 0;
}