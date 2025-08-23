#include "RealAlgebraicNumber.h"
#include <string>
#include <iostream>
//#include "Polynomial3.h"
#include "MyTimer.h"

#ifdef NORMAL_CONFIG
	#include "TestOperators.h"
#endif

using Ran = RealAlgebraicNumber;

// Helper macro for testing assertions
#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        std::cerr << "FAILED: " << message /*<< " at " << __FILE__ << ":" << __LINE__ */ << std::endl; \
        failedTests++; \
    } /*else { \
        std::cout << "PASSED: " << message << std::endl; \
    }*/

// Helper function to print a Ran
static void printRAN(const std::string& name, const Ran& ran) {
	std::cout << name << " = " << ran.toString() << "\n"; //(Decimal: " << ran.toDecimalString(15) << ")\n";
}

void extensiveTest() {
	int failedTests = 0;

	/*Ran ranDefault;
	TEST_ASSERT(ranDefault.isZero(), "Default constructor creates zero")

	Ran ranIntPos(5);
	TEST_ASSERT(ranIntPos == Ran(Polynomial({ -5, 1 }), 5, 5), "Integer constructor (positive)")

	Ran ranIntNeg(-3);
	TEST_ASSERT(ranIntNeg == Ran(Polynomial({ 3, 1 }), -3, -3), "Integer constructor (negative)")

	Ran ranIntZero(0);
	TEST_ASSERT(ranIntZero.isZero(), "Integer constructor (zero)")

	Ran ranDoublePos(2.5);
	auto ranDoublePosTest = Ran({Rational(-5, 2), 1}, Rational(25, 10), Rational(25, 10));
	TEST_ASSERT(ranDoublePos == ranDoublePosTest,
	            "Double constructor (positive)\n\t" + ranDoublePos.toString() +
	            "\n\t" + ranDoublePosTest.toString())

	Ran ranDoubleNeg(-1.75);
	auto ranDoubleNegTest = Ran({Rational(7, 4), 1}, Rational(-175, 100), Rational(-175, 100));
	TEST_ASSERT(ranDoubleNeg == ranDoubleNegTest,
	            "Double constructor (positive)\n\t" + ranDoubleNeg.toString() +
	            "\n\t" + ranDoubleNegTest.toString())

	Ran ranDoubleZero(0.0);
	TEST_ASSERT(ranDoubleZero.isZero(), "Double constructor (zero)")*/

	/*Ran sqrt2PolyInterval(Polynomial({-2, 0, 1}), Rational(14, 10), Rational(15, 10));

	Ran cbrt3CoeffsInterval(Polynomial({-3, 0, 0, 1}), Rational(14, 10), Rational(15, 10));*/

	Ran two(2);
	Ran three(3);
	Ran minusFive(-5);
	Ran half(Rational(1,2));

	// Addition
	Ran sum1 = two + three;
	TEST_ASSERT(sum1 == Ran(5), "2 + 3 == 5")

	Ran sum2 = two + minusFive;
	TEST_ASSERT(sum2 == Ran(-3), "2 + (-5) == -3")

	Ran sum3 = minusFive + two;
	TEST_ASSERT(sum3 == Ran(-3), "(-5) + 2 == -3")

	Ran sumZero = two + Ran(-2);
	TEST_ASSERT(sumZero.isZero(), "2 + (-2) == 0")

	// Subtraction
	Ran diff1 = three - two;
	TEST_ASSERT(diff1 == Ran(1), "3 - 2 == 1")

	Ran diff2 = two - three;
	TEST_ASSERT(diff2 == Ran(-1), "2 - 3 == -1")

	Ran diffSelf = two - two;
	TEST_ASSERT(diffSelf.isZero(), "2 - 2 == 0")

	// Unary Negation
	Ran negTwo = -two;
	TEST_ASSERT(negTwo == Ran(-2), "-(2) == -2")

	Ran negMinusFive = -minusFive;
	TEST_ASSERT(negMinusFive == Ran(5), "-(-5) == 5")

	// Multiplication
	Ran prod1 = two * three;
	TEST_ASSERT(prod1 == Ran(6), "2 * 3 == 6")

	Ran prod2 = two * minusFive;
	TEST_ASSERT(prod2 == Ran(-10), "2 * (-5) == -10")

	Ran prodZero = two * Ran(0);
	TEST_ASSERT(prodZero.isZero(), "2 * 0 == 0")

	Ran prodOne = two * Ran(1);
	TEST_ASSERT(prodOne == two, "2 * 1 == 2")

	// Division
	Ran div1 = three / two;
	TEST_ASSERT(div1 == half + Ran(1), "3 / 2 == 1.5") // 1.5 is 1 + 0.5

	Ran div2 = minusFive / two;
	TEST_ASSERT(div2 == Ran(Rational(-5, 2)), "(-5) / 2 == -2.5")

	Ran divSelf = two / two;
	TEST_ASSERT(divSelf == Ran(1), "2 / 2 == 1")

	// Division by zero - should ideally be handled by Rational class or throw
	// For now, assuming Rational handles division by zero gracefully (e.g., infinity or error state)
	// Ran div_by_zero = two / Ran(0); // This will likely crash or produce invalid results

	Ran aVal(10);
	Ran bVal(4);

	Ran testA = aVal;
	testA += bVal;
	TEST_ASSERT(testA == Ran(14), "a += b")

	testA = aVal;
	testA -= bVal;
	TEST_ASSERT(testA == Ran(6), "a -= b")

	testA = aVal;
	testA *= bVal;
	TEST_ASSERT(testA == Ran(40), "a *= b")

	testA = aVal;
	testA /= bVal;
	TEST_ASSERT(testA == Ran(Rational(10, 4)), "a /= b")

	Ran x(10);
	Ran y(10);
	Ran z(12);
	Ran negX(-10);

	TEST_ASSERT(x == y, "x == y (equal values)")
	TEST_ASSERT(x != z, "x != z (unequal values)")
	TEST_ASSERT(x < z, "x < z")
	TEST_ASSERT(z > x, "z > x")
	TEST_ASSERT(x <= y, "x <= y (equal)")
	TEST_ASSERT(x <= z, "x <= z (less than)")
	TEST_ASSERT(x >= y, "x >= y (equal)")
	TEST_ASSERT(z >= x, "z >= x (greater than)")

	TEST_ASSERT(negX < x, "-x < x")
	TEST_ASSERT(x > negX, "x > -x")
	TEST_ASSERT(negX <= x, "-x <= x")
	TEST_ASSERT(x >= negX, "x >= -x")

	// Test equality for numbers with same polynomial but different initial intervals
	Ran sqrt2A(Polynomial({-2, 0, 1}), Rational(14, 10), Rational(15, 10)); // ~1.414
	//Ran sqrt2B(Polynomial({-2, 0, 1}), Rational(141, 100), Rational(142, 100)); // ~1.41
	//TEST_ASSERT(sqrt2A == sqrt2B, "sqrt(2) == sqrt(2) (different initial intervals)")

	Ran minusSqrt2(Polynomial({-2, 0, 1}), Rational(-15, 10), Rational(-14, 10)); // ~-1.414
	TEST_ASSERT(sqrt2A != minusSqrt2, "sqrt(2) != -sqrt(2)")
	TEST_ASSERT(minusSqrt2 < sqrt2A, "-sqrt(2) < sqrt(2)")

	// inverse()
	Ran invTwo = two.inverse();
	TEST_ASSERT(invTwo == half, "inverse(2) == 0.5")

	Ran invHalf = half.inverse();
	TEST_ASSERT(invHalf == two, "inverse(0.5) == 2")

	Ran invSqrt2 = sqrt2A.inverse();
	Ran expectedInvSqrt2(Polynomial({1, 0, -2}), Rational(7, 10), Rational(8, 10));
	TEST_ASSERT(invSqrt2 == expectedInvSqrt2, "inverse(sqrt(2)) == sqrt(2)/2")

	Ran doubleInvSqrt2 = invSqrt2.inverse();
	TEST_ASSERT(doubleInvSqrt2 == sqrt2A, "inverse(inverse(sqrt(2))) == sqrt(2)")

	// sqrt()
	Ran four(4);
	Ran sqrtFour = four.sqrt();
	TEST_ASSERT(sqrtFour == two, "sqrt(4) == 2")

	Ran sqrtTwo = two.sqrt();
	TEST_ASSERT(sqrtTwo == sqrt2A, "sqrt(2) == sqrt(2)")

	Ran eight(8);
	Ran cbrtEight = eight.sqrt(3); // Cube root of 8
	TEST_ASSERT(cbrtEight == two, "cbrt(8) == 2")

	Ran negEight(-8);
	Ran cbrtNegEight = negEight.sqrt(3); // Cube root of -8
	TEST_ASSERT(cbrtNegEight == Ran(-2), "cbrt(-8) == -2")

	// Test even root of negative number - should throw
	//try {
	//	Ran negFour(-4);
	//	Ran sqrtNegFour = negFour.sqrt();
	//	TEST_ASSERT(false, "sqrt(-4) did NOT throw exception (FAILURE)") // Should not reach here
	//}
	//catch (const std::invalid_argument& e) {
	//	// TEST_ASSERT(true, "sqrt(-4) threw expected exception (SUCCESS): " + std::string(e.what()))
	//}
	//catch (...) {
	//	TEST_ASSERT(false, "sqrt(-4) threw unexpected exception (FAILURE)")
	//}

	// pow()
	Ran twoPowThree = two.pow(3);
	TEST_ASSERT(twoPowThree == Ran(8), "2^3 == 8")

	Ran twoPowZero = two.pow(0);
	TEST_ASSERT(twoPowZero == Ran(1), "2^0 == 1")

	Ran minusTwoPowThree = negTwo.pow(3);
	TEST_ASSERT(minusTwoPowThree == Ran(-8), "(-2)^3 == -8")

	Ran minusTwoPowTwo = negTwo.pow(2);
	TEST_ASSERT(minusTwoPowTwo == Ran(4), "(-2)^2 == 4")

	Ran sqrt2PowTwo = sqrt2A.pow(2);
	TEST_ASSERT(sqrt2PowTwo == two, "sqrt(2)^2 == 2")

	Ran zeroVal(0);
	TEST_ASSERT(zeroVal.isZero(), "isZero() on 0")
	TEST_ASSERT(!two.isZero(), "isZero() on 2")
	TEST_ASSERT(!sqrt2A.isZero(), "isZero() on sqrt(2)")

	Ran computedZero = (sqrt2A * sqrt2A) - two;
	TEST_ASSERT(computedZero.isZero(), "(sqrt(2)*sqrt(2)) - 2 == 0")

	//Ran piApprox(Polynomial({-314159, 100000, 0, 0, 0, 1}), Rational(314159, 100000),
	//                             Rational(314160, 100000)); // x^5 - 3.14159 = 0, just for testing

	//Ran simpleSqrt2(Polynomial({-2, 0, 1}), Rational(1), Rational(2)); // sqrt(2)

	//Ran oneThird(Rational(1, 3));

	//Ran negOneThird(Rational(-1, 3));

	// (sqrt(2) + sqrt(3))^2
	Ran sqrt3Val(Polynomial({-3, 0, 1}), Rational(17, 10), Rational(18, 10));

	Ran sumSqrt2Sqrt3 = sqrt2A + sqrt3Val;
	Ran expectedVal = Ran(5) + (Ran(2) * Ran(6).
		sqrt());

	Ran chainedResult = sumSqrt2Sqrt3.pow(2);
	TEST_ASSERT(chainedResult == expectedVal, "(sqrt(2) + sqrt(3))^2 == 5 + 2*sqrt(6)")

	//Ran smallDiffA(Polynomial({-1000000001, 1000000000}), Rational(1000000001, 1000000000),
	//                               Rational(1000000001, 1000000000)); // 1.000000001
	//Ran smallDiffB(Polynomial({-1000000002, 1000000000}), Rational(1000000002, 1000000000),
	//                               Rational(1000000002, 1000000000)); // 1.000000002

	//TEST_ASSERT(smallDiffA < smallDiffB, "1.000000001 < 1.000000002")
	//TEST_ASSERT(smallDiffB > smallDiffA, "1.000000002 > 1.000000001")
	//TEST_ASSERT(smallDiffA != smallDiffB, "1.000000001 != 1.000000002")

	Ran cbrt2Val(Polynomial({-2, 0, 0, 1}), Rational(12, 10), Rational(13, 10)); // ~1.2599
	Ran cbrt2Cubed = cbrt2Val.pow(3);
	TEST_ASSERT(cbrt2Cubed == two, "cbrt(2)^3 == 2")

	if (failedTests > 0) {
		std::cerr << "SUMMARY: " << failedTests << " tests FAILED.\n";
	}
}

void testingFunction() {
	PROFILE_FUNCTION

	try {
		//testOperatorsRational();
		//testOperatorsPolynomial();
		//testOperatorsRAN();
		extensiveTest();

		Ran a = 2;
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

		/*Ran b = ((a.sqrt() + 2).sqrt() + 2).sqrt();
		Ran c = ((b.pow(2) - 2).pow(2) - 2).pow(2);*/
		//std::cout << "c = " << c.toDecimalString() << std::endl;
		//std::cout << "a == c = " << (a == c) << std::endl;

		auto aa = (Ran(3504) / 100) - 35;
		auto bb = 9 / Ran(93).sqrt() + (Ran(5654) / 100);
		auto cc = (Ran(192) / 10) / 96;
		auto dd = (Ran(6586) / 100) - 27 - (Ran(7305) / 100) + (
			Ran(2797) / 100) / (Ran(818) / 10).pow(2);
		auto ee = (Ran(9758) / 100) / (Ran(6145) / 100);
		auto ff = Ran(28).pow(2) * 78 + Ran(42).sqrt();
		auto gg = 2 * 46 + Ran(3).pow(2);
		auto hh = (Ran(5431) / 100) - 32 / Ran(83).pow(4) / Ran(94).
			pow(3);
		auto ii = 55 - 32 + 76 - (Ran(597) / 10);
		auto jj = (Ran(2661) / 100).pow(2) + (Ran(1878) / 100) / 48;

		/*if (aa.toDecimalString(3) != "0.04")
			std::cout << "Error: aa != expected value\n";*/
		/*if (bb != Ran(Polynomial({247682299, -8763700, 77500}), Rational(574, 10), Rational(575, 10)))
			std::cout << "Error: bb != expected value\n";
		if (cc != Ran(Polynomial({-1, 5}), Rational(1, 10), Rational(3, 10)))
			std::cout << "Error: cc != expected value\n";
		if (dd != Ran(Polynomial({285931907, 8364050}), Rational(-3418, 100), Rational(-3419, 100)))
			std::cout << "Error: dd != expected value\n";
		if (ee != Ran(Polynomial({-9758, 6145}), Rational(158, 100), Rational(159, 100)))
			std::cout << "Error: ee != expected value\n";*/
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
	for (int i = 0; i < 200; ++i) {
		std::cout << i;
		testingFunction();
	}
	atexit(PrintProfilingReport);

	return 0;
}
