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
		//RealAlgebraicNumber().testOperators();
		RealAlgebraicNumber().extensiveTest();

		RealAlgebraicNumber a = 2;
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
		if (a != temp5)
			std::cout << "Error: a != temp5\n";
		if (t1 != temp4)
			std::cout << "Error: t1 != temp4\n";
		if (t2 != temp3)
			std::cout << "Error: t2 != temp3\n";
		if (t3 != temp2)
			std::cout << "Error: t3 != temp2\n";
		if (t4 != temp1)
			std::cout << "Error: t4 != temp1\n";

		/*RealAlgebraicNumber b = ((a.sqrt() + 2).sqrt() + 2).sqrt();
		RealAlgebraicNumber c = ((b.pow(2) - 2).pow(2) - 2).pow(2);*/
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

		/*if (aa.toDecimalString(3) != "0.04")
			std::cout << "Error: aa != expected value\n";*/
		if (bb != RealAlgebraicNumber(Polynomial({247682299, -8763700, 77500}), Rational(574, 10), Rational(575, 10)))
			std::cout << "Error: bb != expected value\n";
		if (cc != RealAlgebraicNumber(Polynomial({-1, 5}), Rational(1, 10), Rational(3, 10)))
			std::cout << "Error: cc != expected value\n";
		if (dd != RealAlgebraicNumber(Polynomial({285931907, 8364050}), Rational(-3418, 100), Rational(-3419, 100)))
			std::cout << "Error: dd != expected value\n";
		if (ee != RealAlgebraicNumber(Polynomial({-9758, 6145}), Rational(158, 100), Rational(159, 100)))
			std::cout << "Error: ee != expected value\n";
		// ff polynomial overflows int
		/*if (hh.toDecimalString(20) != "54.30999999999918819065")
			std::cout << "Error: ff != expected value\n";*/


		/*std::cout << "aa = " << aa.toDecimalString(20) << std::endl;
		std::cout << "bb = " << bb.toDecimalString(20) << std::endl;
		std::cout << "cc = " << cc.toDecimalString(20) << std::endl;
		std::cout << "dd = " << dd.toDecimalString(20) << std::endl;
		std::cout << "ee = " << ee.toDecimalString(20) << std::endl;
		std::cout << "ff = " << ff.toDecimalString(20) << std::endl;
		std::cout << "gg = " << gg.toDecimalString(20) << std::endl;
		std::cout << "hh = " << hh.toDecimalString(20) << std::endl;
		std::cout << "ii = " << ii.toDecimalString(20) << std::endl;
		std::cout << "jj = " << jj.toDecimalString(20) << std::endl;*/
	}
	catch (const std::exception& e) {
		std::cout << e.what() << '\n';
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
