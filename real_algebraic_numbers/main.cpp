#include "RealAlgebraicNumber.h"
#include <string>
#include <iostream>
//#include "Polynomial3.h"
#include "MyTimer.h"
#include "TestOperators.h"

// Helper macro for testing assertions
#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        std::cerr << "FAILED: " << message /*<< " at " << __FILE__ << ":" << __LINE__ */ << std::endl; \
        failedTests++; \
    } /*else { \
        std::cout << "PASSED: " << message << std::endl; \
    }*/

// Helper function to print a RealAlgebraicNumber
static void printRAN(const std::string& name, const RealAlgebraicNumber& ran) {
	std::cout << name << " = " << ran.toString() << "\n"; //(Decimal: " << ran.toDecimalString(15) << ")\n";
}

void extensiveTest() {
	//std::cout << "--- Starting Extensive RealAlgebraicNumber Tests ---" << std::endl;
	int failedTests = 0;

	// Call the existing basic operator tests first
	//std::cout << "\n--- Running Basic Operator Tests (from testOperators()) ---" << std::endl;
	//std::cout << "--- Basic Operator Tests Finished ---" << std::endl;

	//std::cout << "\n--- Testing Constructors ---" << std::endl;
	RealAlgebraicNumber ranDefault;
	TEST_ASSERT(ranDefault.isZero(), "Default constructor creates zero")
	//printRAN("ran_default", ran_default);

	RealAlgebraicNumber ranIntPos(5);
	TEST_ASSERT(ranIntPos == RealAlgebraicNumber(Polynomial({ -5, 1 }), 5, 5), "Integer constructor (positive)")
	//printRAN("ran_int_pos", ran_int_pos);

	RealAlgebraicNumber ranIntNeg(-3);
	TEST_ASSERT(ranIntNeg == RealAlgebraicNumber(Polynomial({ 3, 1 }), -3, -3), "Integer constructor (negative)")
	//printRAN("ran_int_neg", ran_int_neg);

	RealAlgebraicNumber ranIntZero(0);
	TEST_ASSERT(ranIntZero.isZero(), "Integer constructor (zero)")
	//printRAN("ran_int_zero", ran_int_zero);

	RealAlgebraicNumber ranDoublePos(2.5);
	auto ranDoublePosTest = RealAlgebraicNumber({Rational(-5, 2), 1}, Rational(25, 10), Rational(25, 10));
	TEST_ASSERT(ranDoublePos == ranDoublePosTest,
	            "Double constructor (positive)\n\t" + ranDoublePos.toString() +
	            "\n\t" + ranDoublePosTest.toString())
	//printRAN("ran_double_pos", ran_double_pos);

	RealAlgebraicNumber ranDoubleNeg(-1.75);
	auto ranDoubleNegTest = RealAlgebraicNumber({Rational(7, 4), 1}, Rational(-175, 100), Rational(-175, 100));
	TEST_ASSERT(ranDoubleNeg == ranDoubleNegTest,
	            "Double constructor (positive)\n\t" + ranDoubleNeg.toString() +
	            "\n\t" + ranDoubleNegTest.toString())
	//printRAN("ran_double_neg", ran_double_neg);

	RealAlgebraicNumber ranDoubleZero(0.0);
	TEST_ASSERT(ranDoubleZero.isZero(), "Double constructor (zero)")
	//printRAN("ran_double_zero", ran_double_zero);

	// Test a known irrational number: sqrt(2)
	// Polynomial: x^2 - 2 = 0, Interval: [1.4, 1.5]
	RealAlgebraicNumber sqrt2PolyInterval(Polynomial({-2, 0, 1}), Rational(14, 10), Rational(15, 10));
	//TEST_ASSERT(sqrt2PolyInterval.toDecimalString(5) == "1.41421", "Polynomial+Interval constructor (sqrt(2))")
	//printRAN("sqrt2_poly_interval", sqrt2_poly_interval);

	// Test a known irrational number: cube root of 3
	// Polynomial: x^3 - 3 = 0, Interval: [1.4, 1.5]
	RealAlgebraicNumber cbrt3CoeffsInterval(Polynomial({-3, 0, 0, 1}), Rational(14, 10), Rational(15, 10));
	//TEST_ASSERT(cbrt3CoeffsInterval.toDecimalString(5) == "1.44224", "Coeffs+Interval constructor (cbrt(3))")
	//printRAN("cbrt3_coeffs_interval", cbrt3_coeffs_interval);

	//std::cout << "\n--- Testing Arithmetic Operations ---" << std::endl;
	RealAlgebraicNumber two(2);
	RealAlgebraicNumber three(3);
	RealAlgebraicNumber minusFive(-5);
	RealAlgebraicNumber half(Rational(1, 2));

	// Addition
	RealAlgebraicNumber sum1 = two + three;
	TEST_ASSERT(sum1 == RealAlgebraicNumber(5), "2 + 3 == 5")
	//printRAN("sum1", sum1);

	RealAlgebraicNumber sum2 = two + minusFive;
	TEST_ASSERT(sum2 == RealAlgebraicNumber(-3), "2 + (-5) == -3")
	//printRAN("sum2", sum2);

	RealAlgebraicNumber sum3 = minusFive + two;
	TEST_ASSERT(sum3 == RealAlgebraicNumber(-3), "(-5) + 2 == -3")
	//printRAN("sum3", sum3);

	RealAlgebraicNumber sumZero = two + RealAlgebraicNumber(-2);
	TEST_ASSERT(sumZero.isZero(), "2 + (-2) == 0")
	//printRAN("sum_zero", sum_zero);

	// Subtraction
	RealAlgebraicNumber diff1 = three - two;
	TEST_ASSERT(diff1 == RealAlgebraicNumber(1), "3 - 2 == 1")
	//printRAN("diff1", diff1);

	RealAlgebraicNumber diff2 = two - three;
	TEST_ASSERT(diff2 == RealAlgebraicNumber(-1), "2 - 3 == -1")
	//printRAN("diff2", diff2);

	RealAlgebraicNumber diffSelf = two - two;
	TEST_ASSERT(diffSelf.isZero(), "2 - 2 == 0")
	//printRAN("diff_self", diff_self);

	// Unary Negation
	RealAlgebraicNumber negTwo = -two;
	TEST_ASSERT(negTwo == RealAlgebraicNumber(-2), "-(2) == -2")
	//printRAN("neg_two", neg_two);

	RealAlgebraicNumber negMinusFive = -minusFive;
	TEST_ASSERT(negMinusFive == RealAlgebraicNumber(5), "-(-5) == 5")
	//printRAN("neg_minus_five", neg_minus_five);

	// Multiplication
	RealAlgebraicNumber prod1 = two * three;
	TEST_ASSERT(prod1 == RealAlgebraicNumber(6), "2 * 3 == 6")
	//printRAN("prod1", prod1);

	RealAlgebraicNumber prod2 = two * minusFive;
	TEST_ASSERT(prod2 == RealAlgebraicNumber(-10), "2 * (-5) == -10")
	//printRAN("prod2", prod2);

	RealAlgebraicNumber prodZero = two * RealAlgebraicNumber(0);
	TEST_ASSERT(prodZero.isZero(), "2 * 0 == 0")
	//printRAN("prod_zero", prod_zero);

	RealAlgebraicNumber prodOne = two * RealAlgebraicNumber(1);
	TEST_ASSERT(prodOne == two, "2 * 1 == 2")
	//printRAN("prod_one", prod_one);

	// Division
	RealAlgebraicNumber div1 = three / two;
	/*std::cout << div1.toString() << std::endl;
	auto test = half + RealAlgebraicNumber(1);
	std::cout << test.toString() << std::endl;*/
	TEST_ASSERT(div1 == half + RealAlgebraicNumber(1), "3 / 2 == 1.5") // 1.5 is 1 + 0.5
	//printRAN("div1", div1);

	RealAlgebraicNumber div2 = minusFive / two;
	TEST_ASSERT(div2 == RealAlgebraicNumber(Rational(-5, 2)), "(-5) / 2 == -2.5")
	//printRAN("div2", div2);

	RealAlgebraicNumber divSelf = two / two;
	TEST_ASSERT(divSelf == RealAlgebraicNumber(1), "2 / 2 == 1")
	//printRAN("div_self", div_self);

	// Division by zero - should ideally be handled by Rational class or throw
	// For now, assuming Rational handles division by zero gracefully (e.g., infinity or error state)
	// RealAlgebraicNumber div_by_zero = two / RealAlgebraicNumber(0); // This will likely crash or produce invalid results

	//std::cout << "\n--- Testing Compound Assignment Operations ---" << std::endl;
	RealAlgebraicNumber aVal(10);
	RealAlgebraicNumber bVal(4);

	RealAlgebraicNumber testA = aVal;
	testA += bVal;
	TEST_ASSERT(testA == RealAlgebraicNumber(14), "a += b")
	//printRAN("test_a after +=", test_a);

	testA = aVal;
	testA -= bVal;
	TEST_ASSERT(testA == RealAlgebraicNumber(6), "a -= b")
	//printRAN("test_a after -=", test_a);

	testA = aVal;
	testA *= bVal;
	TEST_ASSERT(testA == RealAlgebraicNumber(40), "a *= b")
	//printRAN("test_a after *=", test_a);

	testA = aVal;
	testA /= bVal;
	TEST_ASSERT(testA == RealAlgebraicNumber(Rational(10, 4)), "a /= b")
	//printRAN("test_a after /=", test_a);

	//std::cout << "\n--- Testing Comparison Operators ---" << std::endl;
	RealAlgebraicNumber x(10);
	RealAlgebraicNumber y(10);
	RealAlgebraicNumber z(12);
	RealAlgebraicNumber negX(-10);

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
	RealAlgebraicNumber sqrt2A(Polynomial({-2, 0, 1}), Rational(14, 10), Rational(15, 10)); // ~1.414
	RealAlgebraicNumber sqrt2B(Polynomial({-2, 0, 1}), Rational(141, 100), Rational(142, 100)); // ~1.41
	TEST_ASSERT(sqrt2A == sqrt2B, "sqrt(2) == sqrt(2) (different initial intervals)")
	//printRAN("sqrt2_a", sqrt2_a);
	//printRAN("sqrt2_b", sqrt2_b);

	RealAlgebraicNumber minusSqrt2(Polynomial({-2, 0, 1}), Rational(-15, 10), Rational(-14, 10)); // ~-1.414
	TEST_ASSERT(sqrt2A != minusSqrt2, "sqrt(2) != -sqrt(2)")
	TEST_ASSERT(minusSqrt2 < sqrt2A, "-sqrt(2) < sqrt(2)")

	//std::cout << "\n--- Testing Special Functions ---" << std::endl;

	// inverse()
	RealAlgebraicNumber invTwo = two.inverse();
	TEST_ASSERT(invTwo == half, "inverse(2) == 0.5")
	//printRAN("inv_two", inv_two);

	RealAlgebraicNumber invHalf = half.inverse();
	TEST_ASSERT(invHalf == two, "inverse(0.5) == 2")
	//printRAN("inv_half", inv_half);

	RealAlgebraicNumber invSqrt2 = sqrt2A.inverse();
	RealAlgebraicNumber expectedInvSqrt2(Polynomial({1, 0, -2}), Rational(7, 10), Rational(8, 10));
	// 1/sqrt(2) = sqrt(2)/2
	TEST_ASSERT(invSqrt2 == expectedInvSqrt2, "inverse(sqrt(2)) == sqrt(2)/2")
	//printRAN("inv_sqrt2", inv_sqrt2);

	RealAlgebraicNumber doubleInvSqrt2 = invSqrt2.inverse();
	TEST_ASSERT(doubleInvSqrt2 == sqrt2A, "inverse(inverse(sqrt(2))) == sqrt(2)")
	//printRAN("double_inv_sqrt2", double_inv_sqrt2);

	// sqrt()
	RealAlgebraicNumber four(4);
	RealAlgebraicNumber sqrtFour = four.sqrt();
	TEST_ASSERT(sqrtFour == two, "sqrt(4) == 2")
	//printRAN("sqrt_four", sqrt_four);

	RealAlgebraicNumber sqrtTwo = two.sqrt();
	TEST_ASSERT(sqrtTwo == sqrt2A, "sqrt(2) == sqrt(2)")
	//printRAN("sqrt_two", sqrt_two);

	RealAlgebraicNumber eight(8);
	RealAlgebraicNumber cbrtEight = eight.sqrt(3); // Cube root of 8
	TEST_ASSERT(cbrtEight == two, "cbrt(8) == 2")
	//printRAN("cbrt_eight", cbrt_eight);

	RealAlgebraicNumber negEight(-8);
	RealAlgebraicNumber cbrtNegEight = negEight.sqrt(3); // Cube root of -8
	TEST_ASSERT(cbrtNegEight == RealAlgebraicNumber(-2), "cbrt(-8) == -2")
	//printRAN("cbrt_neg_eight", cbrt_neg_eight);

	// Test even root of negative number - should throw
	//std::cout << "Attempting sqrt(-4) (should throw): ";
	try {
		RealAlgebraicNumber negFour(-4);
		RealAlgebraicNumber sqrtNegFour = negFour.sqrt();
		TEST_ASSERT(false, "sqrt(-4) did NOT throw exception (FAILURE)") // Should not reach here
	}
	catch (const std::invalid_argument& e) {
		// TEST_ASSERT(true, "sqrt(-4) threw expected exception (SUCCESS): " + std::string(e.what()))
	}
	catch (...) {
		TEST_ASSERT(false, "sqrt(-4) threw unexpected exception (FAILURE)")
	}

	// pow()
	RealAlgebraicNumber twoPowThree = two.pow(3);
	TEST_ASSERT(twoPowThree == RealAlgebraicNumber(8), "2^3 == 8")
	//printRAN("two_pow_three", two_pow_three);

	RealAlgebraicNumber twoPowZero = two.pow(0);
	TEST_ASSERT(twoPowZero == RealAlgebraicNumber(1), "2^0 == 1")
	//printRAN("two_pow_zero", two_pow_zero);

	RealAlgebraicNumber minusTwoPowThree = negTwo.pow(3);
	TEST_ASSERT(minusTwoPowThree == RealAlgebraicNumber(-8), "(-2)^3 == -8")
	//printRAN("minus_two_pow_three", minus_two_pow_three);

	RealAlgebraicNumber minusTwoPowTwo = negTwo.pow(2);
	TEST_ASSERT(minusTwoPowTwo == RealAlgebraicNumber(4), "(-2)^2 == 4")
	//printRAN("minus_two_pow_two", minus_two_pow_two);

	RealAlgebraicNumber sqrt2PowTwo = sqrt2A.pow(2);
	TEST_ASSERT(sqrt2PowTwo == two, "sqrt(2)^2 == 2")
	//printRAN("sqrt2_pow_two", sqrt2_pow_two);

	//RealAlgebraicNumber two_pow_half = two.pow(Rational(1, 2)); // Assuming pow can take Rational exponent
	// NOTE: The current pow function only takes int exponent. This test would require a different pow overload.
	// For now, testing integer powers only.

	// isZero()
	RealAlgebraicNumber zeroVal(0);
	TEST_ASSERT(zeroVal.isZero(), "isZero() on 0")
	TEST_ASSERT(!two.isZero(), "isZero() on 2")
	TEST_ASSERT(!sqrt2A.isZero(), "isZero() on sqrt(2)")

	// Test a number that is computationally zero but might not be exactly 0
	// e.g., (sqrt(2) * sqrt(2)) - 2
	RealAlgebraicNumber computedZero = (sqrt2A * sqrt2A) - two;
	TEST_ASSERT(computedZero.isZero(), "(sqrt(2)*sqrt(2)) - 2 == 0")
	//printRAN("computed_zero", computed_zero);

	//std::cout << "\n--- Testing toDecimalString with various precisions ---" << std::endl;
	RealAlgebraicNumber piApprox(Polynomial({-314159, 100000, 0, 0, 0, 1}), Rational(314159, 100000),
	                             Rational(314160, 100000)); // x^5 - 3.14159 = 0, just for testing
	// This is not actually pi, but a number defined by a polynomial and interval to test string conversion.
	// Let's use a simpler irrational for `toDecimalString`
	RealAlgebraicNumber simpleSqrt2(Polynomial({-2, 0, 1}), Rational(1), Rational(2)); // sqrt(2)
	//printRAN("simple_sqrt2", simple_sqrt2);

	//std::string sSqrt25 = simpleSqrt2.toDecimalString(5);
	//std::cout << "sqrt(2) to 5 decimal places: " << s_sqrt2_5 << std::endl;
	//TEST_ASSERT(sSqrt25 == "1.41421", "toDecimalString(5) for sqrt(2)")

	//std::string sSqrt210 = simpleSqrt2.toDecimalString(10);
	//std::cout << "sqrt(2) to 10 decimal places: " << s_sqrt2_10 << std::endl;
	//TEST_ASSERT(sSqrt210 == "1.4142135623", "toDecimalString(10) for sqrt(2)")

	//std::string sSqrt21 = simpleSqrt2.toDecimalString(1);
	//std::cout << "sqrt(2) to 1 decimal place: " << s_sqrt2_1 << std::endl;
	//TEST_ASSERT(sSqrt21 == "1.4", "toDecimalString(1) for sqrt(2)")

	RealAlgebraicNumber oneThird(Rational(1, 3));
	//std::string sOneThird5 = oneThird.toDecimalString(5);
	//std::cout << "1/3 to 5 decimal places: " << s_one_third_5 << std::endl;
	//TEST_ASSERT(sOneThird5 == "0.33333", "toDecimalString(5) for 1/3")

	RealAlgebraicNumber negOneThird(Rational(-1, 3));
	//std::string sNegOneThird5 = negOneThird.toDecimalString(5);
	//std::cout << "-1/3 to 5 decimal places: " << s_neg_one_third_5 << std::endl;
	//TEST_ASSERT(sNegOneThird5 == "-0.33333", "toDecimalString(5) for -1/3")


	//std::cout << "\n--- Testing Complex Chained Operations ---" << std::endl;
	// (sqrt(2) + sqrt(3))^2
	RealAlgebraicNumber sqrt3Val(Polynomial({-3, 0, 1}), Rational(17, 10), Rational(18, 10)); // sqrt(3)
	//printRAN("sqrt3_val", sqrt3_val);

	RealAlgebraicNumber sumSqrt2Sqrt3 = sqrt2A + sqrt3Val;
	//printRAN("sum_sqrt2_sqrt3", sum_sqrt2_sqrt3);
	// (sqrt(2) + sqrt(3))^2 = 2 + 3 + 2*sqrt(6) = 5 + 2*sqrt(6)
	RealAlgebraicNumber expectedVal = RealAlgebraicNumber(5) + (RealAlgebraicNumber(2) * RealAlgebraicNumber(6).
		sqrt());
	//printRAN("expected_val", expected_val);

	RealAlgebraicNumber chainedResult = sumSqrt2Sqrt3.pow(2);
	//printRAN("chained_result", chained_result);
	TEST_ASSERT(chainedResult == expectedVal, "(sqrt(2) + sqrt(3))^2 == 5 + 2*sqrt(6)")


	// Test for very close numbers
	//std::cout << "\n--- Testing very close numbers ---" << std::endl;
	RealAlgebraicNumber smallDiffA(Polynomial({-1000000001, 1000000000}), Rational(1000000001, 1000000000),
	                               Rational(1000000001, 1000000000)); // 1.000000001
	RealAlgebraicNumber smallDiffB(Polynomial({-1000000002, 1000000000}), Rational(1000000002, 1000000000),
	                               Rational(1000000002, 1000000000)); // 1.000000002

	TEST_ASSERT(smallDiffA < smallDiffB, "1.000000001 < 1.000000002")
	TEST_ASSERT(smallDiffB > smallDiffA, "1.000000002 > 1.000000001")
	TEST_ASSERT(smallDiffA != smallDiffB, "1.000000001 != 1.000000002")
	//printRAN("small_diff_a", small_diff_a);
	//printRAN("small_diff_b", small_diff_b);

	// Test a number that is a root of a higher degree polynomial
	// e.g., a root of x^5 - 10x + 1 = 0
	// This requires finding an isolating interval for a specific root.
	// For demonstration, let's pick a polynomial with known rational roots for simplicity,
	// or a simple irrational one that's not just sqrt.
	// Let's test x^3 - 2 = 0 (cbrt(2))
	RealAlgebraicNumber cbrt2Val(Polynomial({-2, 0, 0, 1}), Rational(12, 10), Rational(13, 10)); // ~1.2599
	//printRAN("cbrt2_val", cbrt2_val);
	RealAlgebraicNumber cbrt2Cubed = cbrt2Val.pow(3);
	TEST_ASSERT(cbrt2Cubed == two, "cbrt(2)^3 == 2")
	//printRAN("cbrt2_cubed", cbrt2_cubed);

	//std::cout << "\n--- Extensive RealAlgebraicNumber Tests Finished ---" << std::endl;
	if (failedTests > 0) {
		std::cerr << "SUMMARY: " << failedTests << " tests FAILED.\n";
	}
	/*else {
		std::cout << "SUMMARY: All tests PASSED.\n";
	}*/
}

void testingFunction() {
	PROFILE_FUNCTION
	/*auto a = RealAlgebraicNumber({ -2,1 }, { 1.9,2.1 });
	auto b = RealAlgebraicNumber({ -3,1 }, { 2.9,3.1 });*/
	/*auto a = RealAlgebraicNumber({-2,0,1}, {1.4,1.45});
	auto b = RealAlgebraicNumber({-3,0,1}, { 1.7,1.75 });*/
	/*auto a = RealAlgebraicNumber({ -2,1,2 }, { {75,100},{8,10} });
	auto b = RealAlgebraicNumber({ -3,1,3 }, { {8,10},{85,100} });*/
	try {
		testOperatorsRational();
		testOperatorsPolynomial();
		testOperatorsRAN();
		//RealAlgebraicNumber().testOperators();
		extensiveTest();

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
	//InitializePerformanceFrequency();
	for (int i = 0; i < 1; ++i) {
		std::cout << i;
		testingFunction();
	}
	//atexit(PrintProfilingReport);

	return 0;
}
