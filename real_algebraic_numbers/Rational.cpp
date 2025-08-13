#include "Rational.h"
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <numeric>
#include <vector>

#include "MyTimer.h"

void Rational::simplify() {
	if (denominator == 0) {
		throw std::invalid_argument("Denominator cannot be zero 1");
	}

	// Early exit for already simplified cases
	if (numerator == 0) {
		denominator = 1;
		return;
	}

	// Handle sign normalization first
	if (denominator < 0) {
		numerator = -numerator;
		denominator = -denominator;
	}

	const cpp_int g = boost::multiprecision::gcd(numerator, denominator);
	if (g > 1) {
		numerator /= g;
		denominator /= g;
	}
}

Rational::Rational(const cpp_int& num, const cpp_int& den) : numerator(num), denominator(den) {
	if (denominator == 0) {
		throw std::invalid_argument("Denominator cannot be zero 2");
	}
	simplify();
}

Rational::Rational(const int num, const int den) : numerator(num), denominator(den) {
	if (denominator == 0) {
		throw std::invalid_argument("Denominator cannot be zero 3");
	}
	simplify();
}

Rational::Rational(const cpp_int& numerator) : numerator(numerator), denominator(1) {}

Rational::Rational(const int numerator) : numerator(numerator), denominator(1) {}

Rational::Rational(const long long numerator) : numerator(numerator), denominator(1) {}

Rational::Rational(const double numerator, const cpp_int& maxDenominator) {
	//PROFILE_FUNCTION
	/*std::cout << "In Double Constructor Rational: " << numerator << ", Max Denominator: " << maxDenominator <<
		std::endl;*/
	if (maxDenominator <= 0) {
		throw std::invalid_argument("Max denominator must be positive");
	}

	constexpr double epsilon = 1e-10; // Tolerance for approximation
	int h[3] = {0, 1, 0};
	int k[3] = {1, 0, 0};
	double x = std::abs(numerator);
	const int sign = numerator < 0 ? -1 : 1;

	// Continued fraction coefficients
	for (int i = 0; ; ++i) {
		const int a = static_cast<int>(std::floor(x));
		h[2] = a * h[1] + h[0];
		k[2] = a * k[1] + k[0];

		if (k[2] > maxDenominator) {
			// Compare previous two convergents to choose the best
			const double diff1 = std::abs(static_cast<double>(h[1]) / k[1] - x);
			const double diff2 = std::abs(static_cast<double>(h[2]) / k[2] - x);
			if (diff1 < diff2) {
				h[2] = h[1];
				k[2] = k[1];
			}
			break;
		}

		// Shift history
		h[0] = h[1];
		h[1] = h[2];
		k[0] = k[1];
		k[1] = k[2];

		if (std::abs(x - a) < epsilon) break; // Terminate if exact
		x = 1.0 / (x - a);
	}

	this->numerator = sign * h[2];
	if (k[2] == 0) {
		throw std::invalid_argument("Denominator cannot be zero 4");
	}
	this->denominator = k[2];

	simplify();
}

Rational::Rational(const float numerator) : Rational(static_cast<double>(numerator)) {}

Rational::Rational(const Rational& other) = default;

std::string Rational::toString() const {
	//PROFILE_FUNCTION
	if (denominator == 1) return numerator.str();

	// Pre-allocate string with estimated size to avoid reallocations
	std::string result;
	const std::string numStr = numerator.str();
	const std::string denStr = denominator.str();
	result.reserve(numStr.length() + denStr.length() + 1);

	result = numStr;
	result += '/';
	result += denStr;
	return result;
}

std::string Rational::toDecimalString(int precision) const {
	//PROFILE_FUNCTION
	if (precision < 0) {
		throw std::invalid_argument("Precision must be non-negative.");
	}

	cpp_int num = this->numerator;
	cpp_int den = this->denominator;

	// --- 1. Handle the sign ---
	// The final result is negative if numerator and denominator have opposite signs.
	// We then proceed with positive numbers for the calculation.
	std::string signStr;
	if (num < 0) {
		signStr = "-";
		num *= -1;
	}

	// --- 2. Calculate the integer part ---
	cpp_int integerPart = num / den;

	// --- 3. Calculate the initial remainder for the fractional part ---
	cpp_int remainder = num % den;

	// If there is no remainder, the number is an integer.
	if (remainder == 0) {
		return signStr + integerPart.str() + (precision > 0 ? "." + std::string(precision, '0') : "");
	}

	// --- 4. Calculate the fractional part ---
	std::string fractionalStr;

	// Loop for each decimal place required.
	for (int i = 0; i < precision; ++i) {
		// "Bring down a zero" by multiplying the remainder by 10.
		remainder *= 10;

		// The next digit is the integer result of dividing the new remainder by the denominator.
		cpp_int digit = remainder / den;
		fractionalStr += digit.str();

		// The new remainder is what's left over.
		remainder %= den;

		// Optimization: If remainder becomes 0, the decimal terminates.
		// We can just pad with zeros and finish early.
		if (remainder == 0) {
			fractionalStr += std::string(precision - i - 1, '0');
			break;
		}
	}

	// --- 5. Combine and return the result ---
	if (precision == 0) {
		// Note: Standard rounding would be applied here if needed.
		// This implementation simply truncates.
		return signStr + integerPart.str();
	}

	return signStr + integerPart.str() + "." + fractionalStr;
}

void Rational::print() const {
	std::cout << numerator << "/" << denominator;
}

Rational Rational::operator%(const Rational& other) const {
	//PROFILE_FUNCTION
	if (other.numerator == 0) throw std::invalid_argument("Division by zero");
	const cpp_int num = numerator * other.denominator;
	cpp_int den = denominator * other.numerator;
	if (den == 0) {
		throw std::invalid_argument("Denominator cannot be zero 4");
	}
	return {num % den, den};
}

// Assignment Operators
Rational& Rational::operator+=(const Rational& other) {
	//PROFILE_FUNCTION

	// Handle zero cases
	if (other.numerator == 0) return *this;
	if (numerator == 0) {
		*this = other;
		return *this;
	}

	// Optimize for same denominators
	if (denominator == other.denominator) {
		numerator += other.numerator;
		if (numerator != 0) {
			const cpp_int g = boost::multiprecision::gcd(numerator, denominator);
			if (g > 1) {
				numerator /= g;
				denominator /= g;
			}
		}
		else {
			denominator = 1;
		}
		return *this;
	}

	numerator = numerator * other.denominator + other.numerator * denominator;
	denominator = denominator * other.denominator;
	if (denominator == 0) {
		throw std::invalid_argument("Denominator cannot be zero 4");
	}
	simplify();
	return *this;
}

Rational& Rational::operator-=(const Rational& other) {
	// Handle zero cases
	if (other.numerator == 0) return *this;
	if (numerator == 0) {
		numerator = -other.numerator;
		denominator = other.denominator;
		return *this;
	}

	// Optimize for same denominators
	if (denominator == other.denominator) {
		numerator -= other.numerator;
		if (numerator != 0) {
			const cpp_int g = boost::multiprecision::gcd(numerator, denominator);
			if (g > 1) {
				numerator /= g;
				denominator /= g;
			}
		}
		else {
			denominator = 1;
		}
		return *this;
	}

	numerator = numerator * other.denominator - other.numerator * denominator;
	denominator = denominator * other.denominator;
	if (denominator == 0) {
		throw std::invalid_argument("Denominator cannot be zero 4");
	}
	simplify();
	return *this;
}

//Rational& Rational::operator*=(const Rational& other) {
//	PROFILE_FUNCTION
//	numerator = numerator * other.numerator;
//	denominator = denominator * other.denominator;
//	simplify();
//	return *this;
//}

Rational& Rational::operator*=(const Rational& other) {
	//PROFILE_FUNCTION

	// Handle zero cases early
	if (numerator == 0 || other.numerator == 0) {
		numerator = 0;
		denominator = 1;
		return *this;
	}

	// Handle unit cases
	if (other.numerator == other.denominator) {
		return *this; // Multiplying by 1
	}
	if (numerator == denominator) {
		*this = other; // This is 1, so result is other
		return *this;
	}

	// Cross-reduce to prevent intermediate overflow
	const cpp_int gcd1 = boost::multiprecision::gcd(numerator, other.denominator);
	const cpp_int gcd2 = boost::multiprecision::gcd(other.numerator, denominator);

	numerator = (numerator / gcd1) * (other.numerator / gcd2);
	denominator = (denominator / gcd2) * (other.denominator / gcd1);

	// Handle sign normalization
	if (denominator < 0) {
		numerator = -numerator;
		denominator = -denominator;
	}
	if (denominator == 0) {
		throw std::invalid_argument("Denominator cannot be zero 4");
	}

	return *this;
}

//Rational& Rational::operator/=(const Rational& other) {
//	PROFILE_FUNCTION
//	if (other.numerator == 0) {
//		throw std::overflow_error("Division by zero (rational number has zero numerator).");
//	}
//	numerator = numerator * other.denominator;
//	denominator = denominator * other.numerator;
//	simplify();
//	return *this;
//}

Rational& Rational::operator/=(const Rational& other) {
	if (denominator == 0 || other.denominator == 0) {
		throw std::invalid_argument("Denominator cannot be zero 4");
	}
	if (other.numerator == 0) {
		throw std::overflow_error("Division by zero (rational number has zero numerator).");
	}
	if (numerator == 0) return *this;

	// Handle unit cases
	if (other.numerator == other.denominator) {
		return *this; // Dividing by 1
	}

	// Cross-reduce to prevent intermediate overflow
	const cpp_int gcd1 = boost::multiprecision::gcd(numerator, other.numerator);
	const cpp_int gcd2 = boost::multiprecision::gcd(other.denominator, denominator);

	const cpp_int tempNumerator = (numerator / gcd1) * (other.denominator / gcd2);
	const cpp_int tempDenominator = (denominator / gcd2) * (other.numerator / gcd1);

	numerator = tempNumerator;
	denominator = tempDenominator;

	// Handle sign normalization
	if (denominator < 0) {
		numerator = -numerator;
		denominator = -denominator;
	}

	if (denominator == 0 || other.denominator == 0) {
		throw std::invalid_argument("Denominator cannot be zero 4");
	}

	return *this;
}

Rational operator+(const Rational& lhs, const Rational& rhs) {
	//PROFILE_FUNCTION
	Rational lhsCopy = lhs;
	lhsCopy += rhs;
	return lhsCopy;
}

Rational operator-(const Rational& lhs, const Rational& rhs) {
	//PROFILE_FUNCTION
	Rational lhsCopy = lhs;
	lhsCopy -= rhs;
	return lhsCopy;
}

Rational Rational::operator-() const {
	if (denominator == 0) {
		throw std::invalid_argument("Denominator cannot be zero 4");
	}
	return {-numerator, denominator};
}

Rational operator*(const Rational& lhs, const Rational& rhs) {
	//PROFILE_FUNCTION
	Rational lhsCopy = lhs;
	lhsCopy *= rhs;
	return lhsCopy;
}

Rational operator/(const Rational& lhs, const Rational& rhs) {
	//PROFILE_FUNCTION
	Rational lhsCopy = lhs;
	lhsCopy /= rhs;
	return lhsCopy;
}

// Comparison Operators
bool Rational::operator==(const Rational& other) const {
	// Fast path for identical denominators
	if (denominator == other.denominator) {
		return numerator == other.numerator;
	}
	return numerator * other.denominator == other.numerator * denominator;
}

bool Rational::operator!=(const Rational& other) const {
	return !(*this == other);
}

bool Rational::operator<(const Rational& other) const {
	// Fast path for identical denominators
	if (denominator == other.denominator) {
		return numerator < other.numerator;
	}

	// Check signs first for quick comparison
	const bool thisNeg = numerator < 0;
	const bool otherNeg = other.numerator < 0;

	if (thisNeg && !otherNeg) return true;
	if (!thisNeg && otherNeg) return false;

	return numerator * other.denominator < other.numerator * denominator;
}

bool Rational::operator>(const Rational& other) const {
	return other < *this;
}

bool Rational::operator<=(const Rational& other) const {
	return !(*this > other);
}

bool Rational::operator>=(const Rational& other) const {
	return !(*this < other);
}

// Utility Methods
Rational Rational::abs() const {
	if (denominator == 0) {
		throw std::invalid_argument("Denominator cannot be zero 4");
	}
	return numerator >= 0 ? *this : Rational(-numerator, denominator);
}

Rational Rational::inverse() const {
	if (numerator == 0) throw std::invalid_argument("Cannot invert zero");
	return {denominator, numerator};
}

double Rational::sqrt(const int n) const {
	if (n == 0) throw std::invalid_argument("Zeroth root is undefined");
	if (n % 2 == 0 && numerator < 0)
		throw std::invalid_argument("Even root of negative number");
	if (n == 1) return numerator.convert_to<double>() / denominator.convert_to<double>();

	const double val = numerator.convert_to<double>() / denominator.convert_to<double>();
	return std::pow(val, 1.0 / n);
}

Rational Rational::pow(const int n) const {
	if (n == 0) return {1};
	if (n == 1) return *this;
	if (n < 0) {
		if (numerator == 0) throw std::invalid_argument("Cannot raise zero to negative power");
		return Rational(denominator, numerator).pow(-n);
	}

	const cpp_int resNum = boost::multiprecision::pow(numerator, n);
	const cpp_int resDen = boost::multiprecision::pow(denominator, n);
	if (resDen == 0) {
		throw std::invalid_argument("Denominator cannot be zero 4");
	}
	return {resNum, resDen};
}

//Rational Rational::floor() const {
//	return {numerator / denominator, 1};
//}
//
//Rational Rational::ceil() const {
//	cpp_int q = numerator / denominator;
//	cpp_int r = numerator % denominator;
//
//	if (r < 0) {
//		q -= 1;
//		r += denominator;
//	}
//
//	cpp_int ceil = r == 0 ? q : q + 1;
//
//	return {ceil, 1};
//}

Rational Rational::gcd(const Rational& other) const {
	const cpp_int numGcd = boost::multiprecision::gcd(numerator * other.denominator, other.numerator * denominator);
	const cpp_int denLcm = denominator * other.denominator;
	if (denLcm == 0) {
		throw std::invalid_argument("Denominator cannot be zero 4");
	}
	return {numGcd, denLcm};
}

//bool Rational::isInteger() const {
//	return denominator == 1;
//}

//std::vector<cpp_int> Rational::factorNumerator() const {
//	std::vector<cpp_int> factors;
//	cpp_int n = boost::multiprecision::abs(numerator);
//
//	if (n == 0) {
//		factors.emplace_back(0);
//		return factors;
//	}
//
//	if (n == 1) {
//		factors.emplace_back(1);
//		return factors;
//	}
//
//	// Factor out 2s
//	while (n % 2 == 0) {
//		factors.emplace_back(2);
//		n /= 2;
//	}
//
//	// Check odd factors from 3 onwards
//	for (cpp_int i = 3; i * i <= n; i += 2) {
//		while (n % i == 0) {
//			factors.push_back(i);
//			n /= i;
//		}
//	}
//
//	// If n is still greater than 1, it's a prime factor
//	if (n > 1) {
//		factors.push_back(n);
//	}
//
//	return factors;
//}


void Rational::testOperators() const
{
	std::cout << "Starting comprehensive Rational tests...\n";

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
	Rational g = static_cast<float>(3.14);
	if (g.numerator != 157 || g.denominator != 50)
		std::cout << "FAIL: Rational float constructor not working!\n";

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
	if (ac.numerator != -441 || ac.denominator != 527)
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
	if (af.numerator != -1 || af.denominator != 4)  // Fixed: should check af, not ae
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
	if (bd.numerator != -441 || bd.denominator != 527)
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
	if (int_test.toInt() != 15)
		std::cout << "FAIL: Rational toInt() method not working!\n";

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
	// GCD of 6/8 and 9/12 should be 3/24 = 1/8
	if (gcd_result.numerator != 3 || gcd_result.denominator != 32)
		std::cout << "FAIL: Rational gcd() method not working!\n";

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

	std::cout << "Rational comprehensive tests completed.\n";
}