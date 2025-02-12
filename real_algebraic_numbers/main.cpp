#include "RealAlgebraicNumber.h"
#include <string>
#include <iostream>
//#include "Polynomial3.h"

int main()
{
	/*auto a = RealAlgebraicNumber({-2,0,1}, {1.4,1.5});
	auto b = RealAlgebraicNumber({-3,0,1}, { 1.7,1.8 });*/
	auto a = RealAlgebraicNumber({ -2,1 }, { 1.9,2.1 });
	auto b = RealAlgebraicNumber({ -3,1 }, { 2.9,3.1 });
	/*auto a = RealAlgebraicNumber({ -2,1,2 }, { 0.75,0.8 });
	auto b = RealAlgebraicNumber({ -3,1,3 }, { 0.8,0.85});*/
	auto c = a - b;
	std::cout << a.toString() << std::endl;
	std::cout << "- " << std::endl;
	std::cout << b.toString() << std::endl;
	std::cout << "= " << std::endl;
	std::cout << c.toString() << std::endl;

	return 0;
}