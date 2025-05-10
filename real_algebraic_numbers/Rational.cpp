#include "Rational.h"
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <numeric>

Rational::Rational(const cpp_int& num, const cpp_int& den) {
	if (den == 0) throw std::invalid_argument("Denominator cannot be zero");
	const cpp_int g = boost::multiprecision::gcd(num, den);
	numerator = num / g;
	denominator = den / g;
	if (denominator < 0) {
		numerator = -numerator;
		denominator = -denominator;
	}
}

Rational::Rational(const int num, const int den) {
	if (den == 0) throw std::invalid_argument("Denominator cannot be zero");
	const int g = std::gcd(num, den);
	numerator = num / g;
	denominator = den / g;
	if (denominator < 0) {
		numerator = -numerator;
		denominator = -denominator;
	}
}

Rational::Rational(const cpp_int& numerator) {
	this->numerator = numerator;
	this->denominator = 1;
}

Rational::Rational(const int numerator) {
	this->numerator = cpp_int(numerator);
	this->denominator = 1;
}

Rational::Rational(const double numer, const cpp_int& maxDenominator) {
	if (maxDenominator <= 0) {
		throw std::invalid_argument("Max denominator must be positive");
	}

	constexpr double epsilon = 1e-10; // Tolerance for approximation
	int h[3] = {0, 1, 0};
	int k[3] = {1, 0, 0};
	double x = std::abs(numer);
	const int sign = numer < 0 ? -1 : 1;

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

	numerator = sign * h[2];
	denominator = k[2];

	// Simplify
	const cpp_int g = boost::multiprecision::gcd(numerator, denominator);
	numerator /= g;
	denominator /= g;
}

Rational::Rational(const float numer) : Rational(static_cast<double>(numer)) {}

Rational::Rational(const Rational& other) = default;

std::string Rational::toString() const {
	if (denominator == 1) return numerator.str();
	return numerator.str() + "/" + denominator.str();
}

void Rational::print() const {
	std::cout << numerator << "/" << denominator;
}

// Arithmetic Operators
Rational Rational::operator+(const Rational& other) const {
	cpp_int num = numerator * other.denominator + other.numerator * denominator;
	cpp_int den = denominator * other.denominator;
	return {num, den};
}

Rational Rational::operator-(const Rational& other) const {
	cpp_int num = numerator * other.denominator - other.numerator * denominator;
	cpp_int den = denominator * other.denominator;
	return {num, den};
}

Rational operator-(const cpp_int& lsh, const Rational& other) {
	return Rational(lsh) - other;
}

Rational Rational::operator-() const {
	return {-numerator, denominator};
}

Rational Rational::operator*(const Rational& other) const {
	cpp_int num = numerator * other.numerator;
	cpp_int den = denominator * other.denominator;
	return {num, den};
}

Rational operator*(const cpp_int& lsh, const Rational& other) {
	return Rational(lsh) * other;
}

Rational Rational::operator/(const Rational& other) const {
	if (other.numerator == 0) throw std::invalid_argument("Division by zero");
	cpp_int num = numerator * other.denominator;
	cpp_int den = denominator * other.numerator;
	return {num, den};
}

Rational operator/(const cpp_int& lsh, const Rational& other) {
	return Rational(lsh) / other;
}

Rational Rational::operator%(const Rational& other) const {
	if (other.numerator == 0) throw std::invalid_argument("Division by zero");
	const cpp_int num = numerator * other.denominator;
	cpp_int den = denominator * other.numerator;
	return {num % den, den};
}

// Assignment Operators
Rational& Rational::operator+=(const Rational& other) {
	*this = *this + other;
	return *this;
}

Rational& Rational::operator-=(const Rational& other) {
	*this = *this - other;
	return *this;
}

Rational& Rational::operator*=(const Rational& other) {
	*this = *this * other;
	return *this;
}

Rational& Rational::operator/=(const Rational& other) {
	*this = *this / other;
	return *this;
}

// Comparison Operators
bool Rational::operator==(const Rational& other) const {
	return numerator * other.denominator == other.numerator * denominator;
}

bool Rational::operator!=(const Rational& other) const {
	return !(*this == other);
}

bool Rational::operator<(const Rational& other) const {
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
	if (numerator >= 0) return *this;
	return {-numerator, denominator};
	//return { mpz_class::abs(numerator), denominator.abs()};
}

Rational Rational::inverse() const {
	if (numerator == 0) throw std::invalid_argument("Cannot invert zero");
	return {denominator, numerator};
}

double Rational::sqrt(const int n) const {
	if (n == 0) throw std::invalid_argument("Zeroth root is undefined");
	if (n % 2 == 0 && numerator < 0)
		throw std::invalid_argument("Even root of negative number");
	const double val = numerator.convert_to<double>() / denominator.convert_to<double>();
	return std::pow(val, 1.0 / n);
}

Rational Rational::gcd(const Rational& other) const {
	cpp_int numGcd = boost::multiprecision::gcd(numerator * other.denominator, other.numerator * denominator);
	cpp_int denLcm = (denominator * other.denominator); // / boost::multiprecision::gcd(denominator, other.denominator);
	return {numGcd, denLcm};
}

//cpp_int Rational::computeGcd(const cpp_int& a, const cpp_int& b) {
//    cpp_int g;
//    cpp_int(g.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
//    return g;
//}

std::vector<cpp_int> Rational::factorNumerator() const {
	std::vector<cpp_int> factors;
	//cpp_int num = numerator;
	//int numInt = num.convert_to<int>();
	if (numerator > 0) {
		for (cpp_int i = 1; i <= numerator; ++i) {
			if (numerator % i == 0) factors.push_back(i);
		}
	}
	else {
		for (cpp_int i = -1; i >= numerator; --i) {
			if (numerator % i == 0) factors.push_back(i);
		}
	}
	return factors;
}

bool Rational::isInteger() const {
	return denominator == 1;
}
