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
	try
	{
		/*auto a = RealAlgebraicNumber({ -60,42,-128,256,-12,10 }, { {50,100},{1,1} });
		auto b = RealAlgebraicNumber({ -3,1,3 }, { {8,10},{85,100} });
		auto c = a / b;
		std::cout << a.toString() << std::endl;
		std::cout << "/ " << std::endl;
		std::cout << b.toString() << std::endl;
		std::cout << "= " << std::endl;
		std::cout << c.toString() << std::endl;
		std::cout << "sqrt(a) = " << std::endl;
		auto d = a.sqrt();
		std::cout << d.toString() << std::endl;*/

		/*auto a = RealAlgebraicNumber({ -2,0,1 }, { {14,10},{15,10} });
		std::cout << a.toString() << std::endl;
		std::cout << "a^2 = " << std::endl;
		auto b = a.pow(2);
		std::cout << b.toString() << std::endl;*/

		//auto a = RealAlgebraicNumber({ -16,0,0,0,1 }, { {14,10},{15,10} });

		/*auto a = RealAlgebraicNumber({ -16,1,2 }, { {258,100},{300,100} });
		auto b = a + a;
		std::cout << "a + a = " << b.toString() << std::endl;
		auto c = b - a;
		std::cout << "a + a - a = " << c.toString() << std::endl;
		std::cout << "a == c = " << (a == c) << std::endl;

		auto d = a.sqrt();
		std::cout << "sqrt(a) = " << d.toString() << std::endl;
		auto e = d.pow(2);
		std::cout << "sqrt(a)^2 = " << e.toString() << std::endl;
		std::cout << "a == e = " << (a == e) << std::endl;*/

		auto f = RealAlgebraicNumber({ -2,0,1 }, { {1,1},{2,1} }); // sqrt(2)
		auto g = RealAlgebraicNumber({ -3,0,1 }, { {1,1},{2,1} }); // sqrt(3)
		auto h = f * g;
		std::cout << "sqrt(2) * sqrt(3) = " << h.toString() << std::endl;
	}
	catch (const std::exception& e)
	{
		std::cout << e.what() << std::endl;
	}
	

	return 0;
}