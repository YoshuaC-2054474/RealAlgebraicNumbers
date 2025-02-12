#include "RealAlgebraicNumber.h"
#include <string>
#include <iostream>
//#include "Polynomial3.h"

int main()
{
	/*auto a = RealAlgebraicNumber({ -2,1 }, { 1.9,2.1 });
	auto b = RealAlgebraicNumber({ -3,1 }, { 2.9,3.1 });*/
	/*auto a = RealAlgebraicNumber({-2,0,1}, {1.4,1.45});
	auto b = RealAlgebraicNumber({-3,0,1}, { 1.7,1.75 });*/
	/*auto a = RealAlgebraicNumber({ -2,1,2 }, { {75,100},{8,10} });
	auto b = RealAlgebraicNumber({ -3,1,3 }, { {8,10},{85,100} });*/
	auto a = RealAlgebraicNumber({ -60,42,-128,256,-12,10 }, { {50,100},{1,1} });
	auto b = RealAlgebraicNumber({ -3,1,3 }, { {8,10},{85,100} });
	auto c = a / b;
	std::cout << a.toString() << std::endl;
	std::cout << "/ " << std::endl;
	std::cout << b.toString() << std::endl;
	std::cout << "= " << std::endl;
	std::cout << c.toString() << std::endl;
	std::cout << "sqrt(a) = " << std::endl;
	auto d = a.sqrt();
	std::cout << d.toString() << std::endl;

	return 0;
}