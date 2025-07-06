#include "RealAlgebraicNumber.h"
#include <string>
#include <iostream>
//#include "Polynomial3.h"

#include "MyTimer.h"


void testingFunction() {
	PROFILE_FUNCTION
	/*auto a = RealAlgebraicNumber({ -2,1 }, { 1.9,2.1 });
	auto b = RealAlgebraicNumber({ -3,1 }, { 2.9,3.1 });*/
	/*auto a = RealAlgebraicNumber({-2,0,1}, {1.4,1.45});
	auto b = RealAlgebraicNumber({-3,0,1}, { 1.7,1.75 });*/
	/*auto a = RealAlgebraicNumber({ -2,1,2 }, { {75,100},{8,10} });
	auto b = RealAlgebraicNumber({ -3,1,3 }, { {8,10},{85,100} });*/
	try {
		RealAlgebraicNumber().testOperators();


		auto a = RealAlgebraicNumber({-60, 42, -128, 256, -12, 10}, {{50, 100}, {1, 1}});
		auto b = RealAlgebraicNumber({-3, 1, 3}, {{8, 10}, {85, 100}});
		auto c = a / b;
		/*std::cout << a.toString() << std::endl;
		std::cout << "/ " << std::endl;
		std::cout << b.toString() << std::endl;
		std::cout << "= " << std::endl;
		std::cout << c.toString() << std::endl;
		std::cout << "sqrt(a) = " << std::endl;*/
		auto d = a.sqrt();
		//std::cout << d.toString() << std::endl;

		a = RealAlgebraicNumber({-2, 0, 1}, {{14, 10}, {15, 10}});
		/*std::cout << a.toString() << std::endl;
		std::cout << "a^2 = " << std::endl;*/
		b = a.pow(2);
		//std::cout << b.toString() << std::endl;

		a = RealAlgebraicNumber({-16, 0, 0, 0, 1}, {{14, 10}, {15, 10}});

		a = RealAlgebraicNumber({-16, 1, 2}, {{258, 100}, {300, 100}});
		b = a + a;
		//std::cout << "a + a = " << b.toString() << std::endl;
		c = b - a;
		/*std::cout << "a + a - a = " << c.toString() << std::endl;
		std::cout << "a == c = " << (a == c) << std::endl;*/

		d = a.sqrt();
		/*std::cout << "sqrt(a) = " << d.toString() << std::endl;*/
		auto e = d.pow(2);
		/*std::cout << "sqrt(a)^2 = " << e.toString() << std::endl;
		std::cout << "a == e = " << (a == e) << std::endl;*/

		auto f = RealAlgebraicNumber({-2, 0, 1}, {{1, 1}, {2, 1}}); // sqrt(2)
		auto g = RealAlgebraicNumber({-3, 0, 1}, {{1, 1}, {2, 1}}); // sqrt(3)
		auto h = f * g;
		auto inv = f.inverse();
		/*std::cout << "inverse(sqrt(2)) = " << inv.toString() << std::endl;
		std::cout << "sqrt(2) * sqrt(3) = " << h.toString() << std::endl;*/
		auto l = f - g;
		//std::cout << "sqrt(2) - sqrt(3) = " << l.toString() << std::endl;
		auto i = f.pow(2);
		//std::cout << "sqrt(2)^2 = " << i.toString() << std::endl;
		auto j = g.pow(2);
		//std::cout << "sqrt(3)^2 = " << j.toString() << std::endl;
		auto k = i + j;
		//std::cout << "sqrt(2)^2 + sqrt(3)^2 = " << k.toString() << std::endl;
		RealAlgebraicNumber two = 2;
		/*if (i == two)
			std::cout << "sqrt(2)^2 == 2" << std::endl;*/
		RealAlgebraicNumber div = 2 / 3;
		/*std::cout << "2/3 = " << div.toString() << std::endl;
		std::cout << "inverse(2) = " << two.inverse().toString() << std::endl;*/

		a = 2;
		auto t1 = a.sqrt();
		//std::cout << "sqrt(2) = " << t1.toString() << std::endl;
		auto t2 = t1 + 2;
		//std::cout << "sqrt(2) + 2 = " << t2.toString() << std::endl;
		auto t3 = t2.sqrt();
		//std::cout << "sqrt(sqrt(2) + 2) = " << t3.toString() << std::endl;
		auto t4 = t3 + 2;
		//std::cout << "sqrt(sqrt(2) + 2) + 2 = " << t4.toString() << std::endl;
		auto t5 = t4.sqrt();
		//std::cout << "sqrt(sqrt(sqrt(2) + 2) + 2) = " << t5.toString() << std::endl;

		b = ((a.sqrt() + 2).sqrt() + 2).sqrt();
		//std::cout << "b = " << b.toDecimalString() << std::endl;

		auto temp1 = t5.pow(2);
		//std::cout << "sqrt(sqrt(2) + 2) + 2 = " << temp1.toString() << std::endl;
		auto temp2 = temp1 - 2;
		//std::cout << "sqrt(sqrt(2) + 2) = " << temp2.toString() << std::endl;
		auto temp3 = temp2.pow(2);
		//std::cout << "sqrt(2) + 2 = " << temp3.toString() << std::endl;
		auto temp4 = temp3 - 2;
		//std::cout << "sqrt(2) = " << temp4.toString() << std::endl;
		auto temp5 = temp4.pow(2);
		//std::cout << "2 = " << temp5.toString() << std::endl;

		c = ((b.pow(2) - 2).pow(2) - 2).pow(2);
		//std::cout << "c = " << c.toDecimalString() << std::endl;
		//std::cout << "a == c = " << (a == c) << std::endl;

		auto aa = (RealAlgebraicNumber(3504) / 100) - 35;
		auto bb = 9 / RealAlgebraicNumber(93).sqrt() + (RealAlgebraicNumber(5654) / 100);
		auto cc = (RealAlgebraicNumber(192) / 10) / 96;
		auto dd = (RealAlgebraicNumber(6586) / 100) - 27 - (RealAlgebraicNumber(7305) / 100) + (
			RealAlgebraicNumber(2797) / 100) / (RealAlgebraicNumber(818) / 10).pow(2);
		auto ee = (RealAlgebraicNumber(9758) / 100) / (RealAlgebraicNumber(6145) / 100);
		auto ff = RealAlgebraicNumber(28).pow(2) * 78 + RealAlgebraicNumber(42).sqrt();
		auto gg = 2 * 46 + RealAlgebraicNumber(3).pow(2);
		auto hh = (RealAlgebraicNumber(5431) / 100) - 32 / RealAlgebraicNumber(83).pow(4) / RealAlgebraicNumber(94).
			pow(3);
		auto ii = 55 - 32 + 76 - (RealAlgebraicNumber(597) / 10);
		auto jj = (RealAlgebraicNumber(2661) / 100).pow(2) + (RealAlgebraicNumber(1878) / 100) / 48;

		/*std::cout << "aa = " << aa.toString() << std::endl;
		std::cout << "bb = " << bb.toString() << std::endl;
		std::cout << "cc = " << cc.toString() << std::endl;
		std::cout << "dd = " << dd.toString() << std::endl;
		std::cout << "ee = " << ee.toString() << std::endl;
		std::cout << "ff = " << ff.toString() << std::endl;
		std::cout << "gg = " << gg.toString() << std::endl;
		std::cout << "hh = " << hh.toString() << std::endl;
		std::cout << "ii = " << ii.toString() << std::endl;
		std::cout << "jj = " << jj.toString() << std::endl;*/
	}
	catch (const std::exception& e) {
		std::cout << e.what() << std::endl;
	}
}

int main() {
	InitializePerformanceFrequency();
	for (int i = 0; i < 1; ++i) {
		std::cout << i;
		testingFunction();
	}
	atexit(PrintProfilingReport);

	return 0;
}
