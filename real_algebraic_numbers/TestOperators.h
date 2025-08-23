#ifndef TEST_OPERATORS_H
#define TEST_OPERATORS_H

#ifdef NORMAL_CONFIG
	#include "Rational.h"
#endif

#include "RealAlgebraicNumber.h"
#include "Polynomial.h"

void fail(const std::string& test_name) {
	std::cout << "FAIL: " << test_name << " not working!\n";
}

void testOperatorsRational()
{
	// =================================================================================
	// CONSTRUCTOR TESTS + SIMPLIFY
	// =================================================================================

	// Default constructor
	Rational a;
	if (a.numerator != 0 || a.denominator != 1)
		std::cout << "FAIL: Rational default constructor not working!\n";

	// cpp_int constructor with simplification
	auto b = Rational(cpp_int(797634), cpp_int(2658780));
	if (b.numerator != 3 || b.denominator != 10)
		std::cout << "FAIL: Rational cpp_int constructor not working!\n";

	// int constructor with simplification
	auto c = Rational(797634, 2658780);
	if (c.numerator != 3 || c.denominator != 10)
		std::cout << "FAIL: Rational int constructor not working!\n";

	// Single cpp_int constructor
	Rational d = cpp_int(100);
	if (d.numerator != 100 || d.denominator != 1)
		std::cout << "FAIL: Rational cpp_int single-param constructor not working!\n";

	// Single int constructor
	Rational e = 100;
	if (e.numerator != 100 || e.denominator != 1)
		std::cout << "FAIL: Rational int single-param constructor not working!\n";

	// Long long constructor
	Rational ll = static_cast<long long>(42);
	if (ll.numerator != 42 || ll.denominator != 1)
		std::cout << "FAIL: Rational long long constructor not working!\n";

	// Double constructor
	Rational f = static_cast<double>(3.14);
	if (f.numerator != 157 || f.denominator != 50)
		std::cout << "FAIL: Rational double constructor not working!\n";

	// Float constructor
	/*Rational g = static_cast<float>(3.14);
	if (g.numerator != 157 || g.denominator != 50)
		std::cout << "FAIL: Rational float constructor not working! " << g << "\n";*/

	// Copy constructor
	Rational h = c;
	if (h.numerator != 3 || h.denominator != 10)
		std::cout << "FAIL: Rational copy constructor not working!\n";

	// =================================================================================
	// ARITHMETIC OPERATORS
	// =================================================================================

	// Addition
	auto ab = Rational(3, 51) + Rational(26, 31);
	if (ab.numerator != 473 || ab.denominator != 527)
		std::cout << "FAIL: Rational + operator not working!\n";

	// Subtraction
	auto ac = Rational(3, 51) - Rational(26, 31);
	if (ac.numerator != -411 || ac.denominator != 527)
		std::cout << "FAIL: Rational - operator not working!\n";

	// Multiplication
	auto ad = Rational(3, 51) * Rational(26, 31);
	if (ad.numerator != 26 || ad.denominator != 527)
		std::cout << "FAIL: Rational * operator not working!\n";

	// Division
	auto ae = Rational(3, 51) / Rational(26, 31);
	if (ae.numerator != 31 || ae.denominator != 442)
		std::cout << "FAIL: Rational / operator not working!\n";

	// Unary minus
	auto af = -Rational(5, 20);
	if (af.numerator != -1 || af.denominator != 4)
		std::cout << "FAIL: Rational unary - operator not working!\n";

	// Modulo operator
	auto mod_test = Rational(7, 3) % Rational(2, 1);
	if (mod_test.numerator != 1 || mod_test.denominator != 3)
		std::cout << "FAIL: Rational % operator not working!\n";

	// =================================================================================
	// ASSIGNMENT OPERATORS
	// =================================================================================

	// Addition assignment
	auto bc = Rational(3, 51);
	bc += Rational(26, 31);
	if (bc.numerator != 473 || bc.denominator != 527)
		std::cout << "FAIL: Rational += operator not working!\n";

	// Subtraction assignment
	auto bd = Rational(3, 51);
	bd -= Rational(26, 31);
	if (bd.numerator != -411 || bd.denominator != 527)
		std::cout << "FAIL: Rational -= operator not working!\n";

	// Multiplication assignment
	auto be = Rational(3, 51);
	be *= Rational(26, 31);
	if (be.numerator != 26 || be.denominator != 527)
		std::cout << "FAIL: Rational *= operator not working!\n";

	// Division assignment
	auto bf = Rational(3, 51);
	bf /= Rational(26, 31);
	if (bf.numerator != 31 || bf.denominator != 442)
		std::cout << "FAIL: Rational /= operator not working!\n";

	// =================================================================================
	// COMPARISON OPERATORS
	// =================================================================================

	// Equality
	if ((Rational(3, 10) == Rational(30, 100)) == false)
		std::cout << "FAIL: Rational == operator not working (should be equal)!\n";
	if ((Rational(3, 10) == Rational(30, 10)) == true)
		std::cout << "FAIL: Rational == operator not working (should not be equal)!\n";

	// Inequality
	if ((Rational(3, 10) != Rational(30, 100)) == true)
		std::cout << "FAIL: Rational != operator not working (should be equal)!\n";
	if ((Rational(3, 10) != Rational(30, 10)) == false)
		std::cout << "FAIL: Rational != operator not working (should not be equal)!\n";

	// Less than
	if ((Rational(3, 10) < Rational(4, 10)) == false)
		std::cout << "FAIL: Rational < operator not working!\n";
	if ((Rational(4, 10) < Rational(3, 10)) == true)
		std::cout << "FAIL: Rational < operator not working (reverse case)!\n";

	// Greater than
	if ((Rational(4, 10) > Rational(3, 10)) == false)
		std::cout << "FAIL: Rational > operator not working!\n";
	if ((Rational(3, 10) > Rational(4, 10)) == true)
		std::cout << "FAIL: Rational > operator not working (reverse case)!\n";

	// Less than or equal
	if ((Rational(3, 10) <= Rational(4, 10)) == false)
		std::cout << "FAIL: Rational <= operator not working!\n";
	if ((Rational(3, 10) <= Rational(3, 10)) == false)
		std::cout << "FAIL: Rational <= operator not working (equal case)!\n";

	// Greater than or equal
	if ((Rational(4, 10) >= Rational(3, 10)) == false)
		std::cout << "FAIL: Rational >= operator not working!\n";
	if ((Rational(3, 10) >= Rational(3, 10)) == false)
		std::cout << "FAIL: Rational >= operator not working (equal case)!\n";

	// =================================================================================
	// CONVERSION OPERATORS AND METHODS
	// =================================================================================

	// toInt() method
	Rational int_test(15, 3);
	if (int res = int_test.toInt(); res != 5)
		std::cout << "FAIL: Rational toInt() method not working! " << res << "\n";

	Rational int_test_round_down(16, 3);
	if (int res = int_test_round_down.toInt(); res != 5)
		std::cout << "FAIL: Rational toInt() method not working! " << res << "\n";

	Rational int_test_round_up(14, 3);
	if (int res = int_test_round_up.toInt(); res != 5)
		std::cout << "FAIL: Rational toInt() method not working! " << res << "\n";

	// double conversion operator
	Rational conv_test(1, 2);
	double conv_result = static_cast<double>(conv_test);
	if (std::abs(conv_result - 0.5) > 1e-10)
		std::cout << "FAIL: Rational double conversion operator not working!\n";

	// =================================================================================
	// STRING METHODS
	// =================================================================================

	// toString() method
	Rational str_test(22, 7);
	std::string str_result = str_test.toString();
	if (str_result != "22/7")
		std::cout << "FAIL: Rational toString() method not working!\n";

	// toString() for integer
	Rational int_str_test(5, 1);
	std::string int_str_result = int_str_test.toString();
	if (int_str_result != "5")
		std::cout << "FAIL: Rational toString() method not working for integers!\n";

	// toDecimalString() method
	Rational dec_test(1, 4);
	std::string dec_result = dec_test.toDecimalString(3);
	if (dec_result != "0.250")
		std::cout << "FAIL: Rational toDecimalString() method not working!\n";

	// =================================================================================
	// UTILITY METHODS
	// =================================================================================

	// abs() method
	Rational abs_test(-3, 4);
	Rational abs_result = abs_test.abs();
	if (abs_result.numerator != 3 || abs_result.denominator != 4)
		std::cout << "FAIL: Rational abs() method not working!\n";

	// abs() method for positive
	Rational abs_pos_test(3, 4);
	Rational abs_pos_result = abs_pos_test.abs();
	if (abs_pos_result.numerator != 3 || abs_pos_result.denominator != 4)
		std::cout << "FAIL: Rational abs() method not working for positive numbers!\n";

	// inverse() method
	Rational inv_test(3, 4);
	Rational inv_result = inv_test.inverse();
	if (inv_result.numerator != 4 || inv_result.denominator != 3)
		std::cout << "FAIL: Rational inverse() method not working!\n";

	// pow() method
	Rational pow_test(2, 3);
	Rational pow_result = pow_test.pow(2);
	if (pow_result.numerator != 4 || pow_result.denominator != 9)
		std::cout << "FAIL: Rational pow() method not working!\n";

	// pow() method with negative exponent
	Rational pow_neg_test(2, 3);
	Rational pow_neg_result = pow_neg_test.pow(-2);
	if (pow_neg_result.numerator != 9 || pow_neg_result.denominator != 4)
		std::cout << "FAIL: Rational pow() method not working for negative exponents!\n";

	// pow() method with zero exponent
	Rational pow_zero_result = pow_test.pow(0);
	if (pow_zero_result.numerator != 1 || pow_zero_result.denominator != 1)
		std::cout << "FAIL: Rational pow() method not working for zero exponent!\n";

	// sqrt() method (returns double)
	Rational sqrt_test(9, 4);
	double sqrt_result = sqrt_test.sqrt(2);
	if (std::abs(sqrt_result - 1.5) > 1e-10)
		std::cout << "FAIL: Rational sqrt() method not working!\n";

	// sqrt() method with different root
	Rational cube_test(8, 1);
	double cube_result = cube_test.sqrt(3);
	if (std::abs(cube_result - 2.0) > 1e-10)
		std::cout << "FAIL: Rational sqrt() method not working for cube root!\n";

	// gcd() method
	Rational gcd_test1(6, 8);
	Rational gcd_test2(9, 12);
	Rational gcd_result = gcd_test1.gcd(gcd_test2);
	if (gcd_result.numerator != 3 || gcd_result.denominator != 4)
		std::cout << "FAIL: Rational gcd() method not working! " << gcd_result << "\n";

	// floor() method
	Rational floor_test(5, 2);
	if (auto fl = floor_test.floor(); fl.numerator != 2 || fl.denominator != 1)
		std::cout << "FAIL: Rational floor() method not working!\n";

	// ceil() method
	Rational ceil_test(5, 2);
	if (auto ce = ceil_test.ceil(); ce.numerator != 3 || ce.denominator != 1)
		std::cout << "FAIL: Rational ceil() method not working!\n";

	// =================================================================================
	// STREAM OUTPUT OPERATOR
	// =================================================================================

	// Test stream output (we can't easily test this without capturing output)
	std::ostringstream oss;
	oss << Rational(22, 7);
	if (oss.str() != "22/7")
		std::cout << "FAIL: Rational stream output operator not working!\n";

	// Test stream output for integer
	std::ostringstream oss2;
	oss2 << Rational(5, 1);
	if (oss2.str() != "5")
		std::cout << "FAIL: Rational stream output operator not working for integers!\n";

	// =================================================================================
	// EDGE CASES AND ERROR CONDITIONS
	// =================================================================================

	// Test zero operations
	Rational zero(0, 1);
	Rational non_zero(3, 4);

	// Zero addition
	Rational zero_add = zero + non_zero;
	if (zero_add.numerator != 3 || zero_add.denominator != 4)
		std::cout << "FAIL: Rational zero addition not working!\n";

	// Zero multiplication
	Rational zero_mult = zero * non_zero;
	if (zero_mult.numerator != 0 || zero_mult.denominator != 1)
		std::cout << "FAIL: Rational zero multiplication not working!\n";

	// Negative number tests
	Rational neg1(-2, 3);
	Rational neg2(2, -3);
	if (neg1.numerator != -2 || neg1.denominator != 3)
		std::cout << "FAIL: Rational negative numerator handling not working!\n";
	if (neg2.numerator != -2 || neg2.denominator != 3)
		std::cout << "FAIL: Rational negative denominator normalization not working!\n";
}

void testOperatorsPolynomial() {
	// =================================================================================
	// CONSTRUCTOR TESTS
	// =================================================================================

	// Default constructor
	Polynomial p_default;
	if (p_default.degree != -1 || !p_default.coefficients.empty()) {
		fail("Default constructor");
	}

	// initializer_list<int> constructor
	Polynomial p_int_list({ 1, 2, 3 }); // Represents 1 + 2x + 3x^2
	if (p_int_list.degree != 2 || p_int_list.coefficients.size() != 3 ||
		p_int_list.coefficients[0] != 1 || p_int_list.coefficients[1] != 2 || p_int_list.coefficients[2] != 3) {
		fail("Initializer list (int) constructor");
	}

	// vector<Rational> constructor
	std::vector<Rational> coeffs_rat = { Rational(1, 2), Rational(0), Rational(3, 4) }; // Represents 1/2 + 3/4x^2
	Polynomial p_rat_vec(coeffs_rat);
	if (p_rat_vec.degree != 2 || p_rat_vec.coefficients.size() != 3 ||
		p_rat_vec.coefficients[0] != Rational(1, 2) || p_rat_vec.coefficients[1] != Rational(0) || p_rat_vec.coefficients[2] != Rational(3, 4)) {
		fail("Vector of Rational constructor");
	}

	// Rational constructor
	Polynomial p_rational_c(Rational(5, 7)); // Represents 5/7
	if (p_rational_c.degree != 0 || p_rational_c.coefficients.size() != 1 || p_rational_c.coefficients[0] != Rational(5, 7)) {
		fail("Rational constructor");
	}

	// isZero() and coeff() tests
	Polynomial p_zero;
	if (!p_zero.isZero() || p_zero.coeff(0) != 0 || p_zero.coeff(5) != 0) {
		fail("isZero() or coeff() for zero polynomial");
	}

	if (p_int_list.isZero() || p_int_list.coeff(1) != 2 || p_int_list.coeff(3) != 0) {
		fail("isZero() or coeff() for non-zero polynomial");
	}

	// =================================================================================
	// ARITHMETIC OPERATORS
	// =================================================================================

	Polynomial p1({ 1, 2, 3 }); // 1 + 2x + 3x^2
	Polynomial p2({ 5, -1 });   // 5 - x

	// operator+
	Polynomial p_add = p1 + p2; // 6 + x + 3x^2
	Polynomial p_add_expected({ 6, 1, 3 });
	if (p_add != p_add_expected) {
		fail("+ operator");
	}

	// operator-
	Polynomial p_sub = p1 - p2; // -4 + 3x + 3x^2
	Polynomial p_sub_expected({ -4, 3, 3 });
	if (p_sub != p_sub_expected) {
		fail("- operator");
	}

	// operator*
	Polynomial p_mul = p1 * p2; // (1+2x+3x^2)(5-x) = 5 -x +10x -2x^2 +15x^2 -3x^3 = 5 + 9x + 13x^2 - 3x^3
	Polynomial p_mul_expected({ 5, 9, 13, -3 });
	if (p_mul != p_mul_expected) {
		fail("* operator");
	}

	// operator/ (integer division, should be based on polyDivide)
	Polynomial p_div = p_mul / p2; // (5 + 9x + 13x^2 - 3x^3) / (5 - x) = 1+2x+3x^2
	if (p_div != p1) {
		fail("/ operator (polynomial)");
	}

	// operator/ (scalar division)
	Polynomial p_div_scalar = p1 / Rational(2); // 1/2 + x + 3/2x^2
	Polynomial p_div_scalar_expected(std::vector<Rational>{Rational(1, 2), Rational(1), Rational(3, 2)});
	if (p_div_scalar != p_div_scalar_expected) {
		fail("/ operator (scalar)");
	}

	// operator- (unary)
	Polynomial p_unary_minus = -p1; // -1 - 2x - 3x^2
	Polynomial p_unary_minus_expected({ -1, -2, -3 });
	if (p_unary_minus != p_unary_minus_expected) {
		fail("Unary - operator");
	}

	// =================================================================================
	// ASSIGNMENT OPERATORS
	// =================================================================================

	Polynomial p_assign1({ 1, 1 }); // 1 + x
	Polynomial p_assign2({ -1, 1 }); // -1 + x

	// operator+=
	Polynomial p_add_assign = p_assign1;
	p_add_assign += p_assign2; // (1+x) + (-1+x) = 2x
	Polynomial p_add_assign_expected({ 0, 2 });
	if (p_add_assign != p_add_assign_expected) {
		fail("+= operator");
	}

	// operator-=
	Polynomial p_sub_assign = p_assign1;
	p_sub_assign -= p_assign2; // (1+x) - (-1+x) = 2
	Polynomial p_sub_assign_expected({ 2 });
	if (p_sub_assign != p_sub_assign_expected) {
		fail("-= operator");
	}

	// operator*=
	Polynomial p_mul_assign = p_assign1;
	p_mul_assign *= p_assign2; // (1+x)(-1+x) = x^2 - 1
	Polynomial p_mul_assign_expected({ -1, 0, 1 });
	if (p_mul_assign != p_mul_assign_expected) {
		fail("*= operator");
	}

	// operator/=
	Polynomial p_div_assign = p_mul_assign;
	p_div_assign /= p_assign2; // (x^2-1) / (-1+x) = 1+x
	if (p_div_assign != p_assign1) {
		fail("/= operator");
	}

	// =================================================================================
	// COMPARISON OPERATORS
	// =================================================================================

	Polynomial p_eq1({ 1, 2, 3 });
	Polynomial p_eq2({ 1, 2, 3 });
	Polynomial p_neq({ 1, 2, 4 });

	// operator==
	if (!(p_eq1 == p_eq2)) {
		fail("== operator (equal)");
	}
	if (p_eq1 == p_neq) {
		fail("== operator (not equal)");
	}

	// operator!=
	if (p_eq1 != p_eq2) {
		fail("!= operator (equal)");
	}
	if (!(p_eq1 != p_neq)) {
		fail("!= operator (not equal)");
	}

	// =================================================================================
	// UTILITY METHODS
	// =================================================================================

	Polynomial p_util({ 1, -2, 1 }); // (x-1)^2 = x^2 - 2x + 1

	// toString()
	if (p_util.toString() != "x^2 - 2x + 1") {
		fail("toString()");
		std::cout << p_util.toString() << "\n";
	}

	// reflectY()
	Polynomial p_reflect = p_util.reflectY(); // ((-x)-1)^2 = x^2+2x+1
	Polynomial p_reflect_expected({ 1, 2, 1 });
	if (p_reflect != p_reflect_expected) {
		fail("reflectY()");
	}

	// evaluate()
	if (p_util.evaluate(Rational(3)) != Rational(4)) { // (3-1)^2 = 4
		fail("evaluate()");
	}

	// getLargestCoeff()
	Polynomial p_largest({ -1, 5, -10, 100, 3 });
	if (p_largest.getLargestCoeff() != Rational(100)) {
		fail("getLargestCoeff()");
	}
}

void testOperatorsRAN() {
	////PROFILE_FUNCTION
	
	// Define some key algebraic numbers for testing
	// sqrt(2) is a root of x^2 - 2 = 0, in the interval [1, 2]
	RealAlgebraicNumber sqrt2(Polynomial({ -2, 0, 1 }), { Rational(1), Rational(2) });
	// sqrt(3) is a root of x^2 - 3 = 0, in the interval [1, 2]
	RealAlgebraicNumber sqrt3(Polynomial({ -3, 0, 1 }), { Rational(1), Rational(2) });
	// 2*sqrt(2) is a root of x^2 - 8 = 0, in the interval [2, 3]
	RealAlgebraicNumber two_sqrt2(Polynomial({ -8, 0, 1 }), { Rational(2), Rational(3) });
	// 5 is a root of x - 5 = 0
	RealAlgebraicNumber five(Polynomial({ -5, 1 }), { Rational(5), Rational(5) });
	// -5 is a root of x + 5 = 0
	RealAlgebraicNumber minus_five(Polynomial({ 5, 1 }), { Rational(-5), Rational(-5) });
	// 0.5 is a root of 2x-1=0
	RealAlgebraicNumber half(Polynomial({ -1, 2 }), { Rational(1, 2), Rational(1, 2) });

	// =================================================================================
	// CONSTRUCTOR TESTS
	// =================================================================================

	// int constructor
	RealAlgebraicNumber ran_int(10);
	if (ran_int.polynomial != Polynomial({ -10, 1 }) || ran_int.interval.lowerBound != 10 || ran_int.interval.upperBound != 10) {
		fail("int constructor");
	}

	// Rational constructor
	RealAlgebraicNumber ran_rational(Rational(1, 2));
	if (ran_rational.polynomial != Polynomial({ -1 , 2 }) || ran_rational.interval.lowerBound != Rational(1, 2) || ran_rational.interval.upperBound != Rational(1, 2)) {
		fail("Rational constructor");
	}

	// double constructor
	RealAlgebraicNumber ran_double(0.5);
	if (ran_double.polynomial != Polynomial({ -1, 2 }) || ran_double.interval.lowerBound != Rational(1, 2) || ran_double.interval.upperBound != Rational(1, 2)) {
		fail("double constructor");
	}

	// copy constructor (implicit)
	RealAlgebraicNumber ran_copy = sqrt2;
	if (ran_copy.polynomial != sqrt2.polynomial || ran_copy.interval.lowerBound != sqrt2.interval.lowerBound || ran_copy.interval.upperBound != sqrt2.interval.upperBound) {
		fail("copy constructor");
	}

	// =================================================================================
	// ARITHMETIC OPERATORS
	// =================================================================================

	// Addition
	RealAlgebraicNumber ran_add = five + half;
	if (ran_add != RealAlgebraicNumber(Rational(11, 2))) {
		fail("+ operator (RAN + RAN)");
	}
	ran_add = sqrt2 + sqrt3;
	if (ran_add.polynomial != Polynomial({ 1, 0, -10, 0, 1 }) || ran_add.interval.lowerBound < 1 || ran_add.interval.upperBound > 4) {
		fail("+ operator (RAN + RAN with sqrt2 and sqrt3)");
	}

	// Subtraction
	RealAlgebraicNumber ran_sub = five - half;
	if (ran_sub != RealAlgebraicNumber(Rational(9, 2))) {
		fail("- operator (RAN - RAN)");
	}
	ran_sub = sqrt2 - sqrt3;
	if (ran_sub.polynomial != Polynomial({ 1, 0, -10, 0, 1 }) || ran_sub.interval.lowerBound < -1 || ran_sub.interval.upperBound > 0) {
		fail("- operator (RAN - RAN with sqrt2 and sqrt3)");
	}

	// Multiplication
	RealAlgebraicNumber ran_mul = two_sqrt2 * sqrt2;
	if (ran_mul != RealAlgebraicNumber(Rational(4))) {
		fail("* operator (RAN * RAN)");
	}
	ran_mul = sqrt2 * sqrt3;
	if (ran_mul.polynomial != Polynomial({ -6, 0, 1 }) || ran_mul.interval.lowerBound < 0 || ran_mul.interval.upperBound > 5) {
		fail("* operator (RAN * RAN with sqrt2 and sqrt3)");
	}

	// Division
	RealAlgebraicNumber ran_div = two_sqrt2 / sqrt2;
	if (ran_div != RealAlgebraicNumber(Rational(2))) {
		fail("/ operator (RAN / RAN)");
	}
	ran_div = sqrt2 / sqrt3;
	if (ran_div.polynomial != Polynomial({ -2, 0, 3 }) || ran_div.interval.lowerBound < 0 || ran_div.interval.upperBound > 5) {
		fail("/ operator (RAN / RAN with sqrt2 and sqrt3)");
	}

	// Mixed type arithmetic (RAN + Rational)
	RealAlgebraicNumber ran_add_rat = five + Rational(0.5);
	if (ran_add_rat != RealAlgebraicNumber(Rational(11, 2))) {
		fail("+ operator (RAN + Rational)");
	}

	// Unary minus
	RealAlgebraicNumber ran_unary_minus = -five;
	if (ran_unary_minus != minus_five) {
		fail("Unary - operator");
	}

	auto orderRes = RealAlgebraicNumber({ 190, -48, 3 }, { .lowerBound = {88, 10}, .upperBound = {89, 10} });
	if (auto orderOfOperations = 2 + 3 * 2 + sqrt2 / sqrt3; orderOfOperations != orderRes)
		std::cout << "operator precedence not working!\n"
		<< "\ta + b * a + c / d = " << orderOfOperations << "\n";

	// =================================================================================
	// ASSIGNMENT OPERATORS
	// =================================================================================

	RealAlgebraicNumber ran_assign = two_sqrt2;
	ran_assign += sqrt2;
	if (ran_assign != RealAlgebraicNumber(Polynomial({ -18, 0, 1 }), { Rational(4), Rational(5) })) { // 3*sqrt(2) is a root of x^2-18=0
		fail("+= operator");
	}

	ran_assign = two_sqrt2;
	ran_assign -= sqrt2;
	if (ran_assign != sqrt2) {
		fail("-= operator");
	}

	ran_assign = sqrt2;
	ran_assign *= sqrt2;
	if (ran_assign != RealAlgebraicNumber(Rational(2))) {
		fail("*= operator");
	}

	ran_assign = two_sqrt2;
	ran_assign /= sqrt2;
	if (ran_assign != RealAlgebraicNumber(Rational(2))) {
		fail("/= operator");
	}

	// =================================================================================
	// COMPARISON OPERATORS
	// =================================================================================

	// Equality (==)
	if (!(sqrt2 == sqrt2)) {
		fail("== operator (equal)");
	}
	if (sqrt2 == sqrt3) {
		fail("== operator (not equal)");
	}

	// Inequality (!=)
	if (sqrt2 != sqrt2) {
		fail("!= operator (equal)");
	}
	if (!(sqrt2 != sqrt3)) {
		fail("!= operator (not equal)");
	}

	// Less than (<)
	if (!(sqrt2 < sqrt3)) {
		fail("< operator");
	}

	// Greater than (>)
	if (!(sqrt3 > sqrt2)) {
		fail("> operator");
	}

	// Less than or equal (<=)
	if (!(sqrt2 <= sqrt2)) {
		fail("<= operator (equal)");
	}

	// Greater than or equal (>=)
	if (!(sqrt3 >= sqrt3)) {
		fail(">= operator (equal)");
	}

	// =================================================================================
	// UTILITY METHODS
	// =================================================================================

	// inverse()
	RealAlgebraicNumber inv_sqrt2 = sqrt2.inverse();
	if (inv_sqrt2 != RealAlgebraicNumber(Polynomial({ -1, 0, 2 }), { Rational(0), Rational(1) })) {
		fail("inverse()");
	}
	if (inv_sqrt2.inverse() != sqrt2) {
		fail("double inverse()");
	}

	// sqrt()
	RealAlgebraicNumber sqrt_of_four = RealAlgebraicNumber(Rational(4)).sqrt();
	if (sqrt_of_four != RealAlgebraicNumber(Rational(2))) {
		fail("sqrt()");
	}
	if (sqrt_of_four.pow(2) != RealAlgebraicNumber(Rational(4))) {
		fail("pow() followed by sqrt() should give original value");
	}

	// pow()
	RealAlgebraicNumber p_sqr = sqrt2.pow(2);
	if (p_sqr != RealAlgebraicNumber(Rational(2))) {
		fail("pow()");
	}

	// isPositive() and isNegative()
	if (!five.isPositive() || five.isNegative()) {
		fail("isPositive()/isNegative() for positive number");
	}
	if (minus_five.isPositive() || !minus_five.isNegative()) {
		fail("isPositive()/isNegative() for negative number");
	}

	RealAlgebraicNumber a = 2;
	RealAlgebraicNumber b = 3;
	auto c = sqrt2;
	auto d = sqrt3;
	std::vector<Rational> ratCoeff = { cpp_int(-9007199254740992), 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };
	auto ratExp = RealAlgebraicNumber(ratCoeff, { .lowerBound = {39, 1}, .upperBound = {40, 1} });
	auto ratExp2 = RealAlgebraicNumber({ -32, 0,0,0,0,0,1 }, { .lowerBound = {17, 10}, .upperBound = {18, 10} });
	auto ratSq = RealAlgebraicNumber({ -8,0,0,0,0,1 }, { .lowerBound = { 15, 10 }, .upperBound = { 16, 10 } });
	auto ratSq2 = RealAlgebraicNumber({ -8,0,0,0,0,0,0,0,0,0,1 }, { .lowerBound = { 12, 10 }, .upperBound = { 13, 10 } });

	// chained expressions
	if (a.inverse().inverse() != a || c.inverse().inverse() != c)
		std::cout << "double inverse() not working!\n";
	if (b.pow(2).sqrt() != b || d.sqrt().pow(2) != d)
		std::cout << "sqrt() and pow() not working together! (1)\n";
	if (a.sqrt().pow(2) != a || c.pow(2).sqrt() != c)
		std::cout << "sqrt() and pow() not working together! (2)\n";

	// sqrt and pow with rationals
	if (a.sqrt({ 5,3 }) != ratSq || c.sqrt({ 5,3 }) != ratSq2)
		std::cout << "sqrt() not working with Rational!\n";
	if (a.pow({ 53,10 }) != ratExp || c.pow({ 5,3 }) != ratExp2)
		std::cout << "pow() not working with Rational!\n";

	// complex order of operations
	auto oooResult = RealAlgebraicNumber({ -63, -6, 1 }, { .lowerBound = { 114, 10 }, .upperBound = { 115, 10 } });
	if (auto complexOrderOfOperations = a + b * c.pow(3) + d / b.sqrt(); complexOrderOfOperations != oooResult)
	{
		std::cout << "operator precedence not working with complex operations!\n"
			<< "\ta + b * c.pow(3) + d / b.sqrt() = " << complexOrderOfOperations << "\n";
	}
}

#endif