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
	const Rational fl = (*this / other).floor();
	Rational r = *this - fl * other;
	return r;
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

// returns closest integer (rounding up or down)
int Rational::toInt() const {
	if (denominator == 1)
		return numerator.convert_to<int>();
	const Rational fl = this->floor();
	const Rational cl = this->ceil();
	if (*this - fl < cl - *this) // Closer to floor
		return fl.numerator.convert_to<int>();
	return cl.numerator.convert_to<int>();
}

Rational Rational::floor() const {
	return {numerator / denominator, 1};
}

Rational Rational::ceil() const {
	cpp_int q = numerator / denominator;
	cpp_int r = numerator % denominator;

	if (r < 0) {
		q -= 1;
		r += denominator;
	}

	cpp_int ceil = r == 0 ? q : q + 1;

	return {ceil, 1};
}

Rational Rational::gcd(const Rational& other) const {
	const cpp_int numGcd = boost::multiprecision::gcd(numerator * other.denominator, other.numerator * denominator);
	const cpp_int denLcm = denominator * other.denominator;
	if (denLcm == 0) {
		throw std::invalid_argument("Denominator cannot be zero 4");
	}
	return {numGcd, denLcm};
}