#include "Rational2.h"
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <numeric>

Rational::Rational(const long long num, const long long den) {
	if (den == 0) throw std::invalid_argument("Denominator cannot be zero");
	const long long g = std::gcd(num, den);
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

Rational::Rational(const long long numerator) {
	this->numerator = numerator;
	this->denominator = 1;
}

Rational::Rational(const int numerator) {
	this->numerator = numerator;
	this->denominator = 1;
}

Rational::Rational(const double numer, const long long maxDenominator) {
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
	const long long g = std::gcd(numerator, denominator);
	numerator /= g;
	denominator /= g;
}

Rational::Rational(const float numerator) : Rational(static_cast<double>(numerator)) {}

Rational::Rational(const Rational& other) = default;

std::string Rational::toString() const {
	if (denominator == 1) return std::to_string(numerator);
	return std::to_string(numerator) + "/" + std::to_string(denominator);
}

void Rational::print() const {
	std::cout << numerator << "/" << denominator;
}

// Arithmetic Operators
Rational Rational::operator+(const Rational& other) const {
	long long num = numerator * other.denominator + other.numerator * denominator;
	long long den = denominator * other.denominator;
	return {num, den};
}

Rational Rational::operator-(const Rational& other) const {
	long long num = numerator * other.denominator - other.numerator * denominator;
	long long den = denominator * other.denominator;
	return {num, den};
}

Rational operator-(const long long lsh, const Rational& other) {
	return Rational(lsh) - other;
}

Rational Rational::operator-() const {
	return {-numerator, denominator};
}

Rational Rational::operator*(const Rational& other) const {
	long long num = numerator * other.numerator;
	long long den = denominator * other.denominator;
	return {num, den};
}

Rational operator*(const long long lsh, const Rational& other) {
	return Rational(lsh) * other;
}

Rational Rational::operator/(const Rational& other) const {
	if (other.numerator == 0) throw std::invalid_argument("Division by zero");
	long long num = numerator * other.denominator;
	long long den = denominator * other.numerator;
	return {num, den};
}

Rational operator/(const long long lsh, const Rational& other) {
	return Rational(lsh) / other;
}

Rational Rational::operator%(const Rational& other) const {
	if (other.numerator == 0) throw std::invalid_argument("Division by zero");
	const long long num = numerator * other.denominator;
	long long den = denominator * other.numerator;
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
	const double val = static_cast<double>(numerator) / static_cast<double>(denominator);
	return std::pow(val, 1.0 / n);
}

Rational Rational::gcd(const Rational& other) const {
	long long numGcd = std::gcd(numerator * other.denominator, other.numerator * denominator);
	long long denLcm = (denominator * other.denominator);
	// / boost::multiprecision::gcd(denominator, other.denominator);
	return {numGcd, denLcm};
}

//cpp_int Rational::computeGcd(const cpp_int& a, const cpp_int& b) {
//    cpp_int g;
//    cpp_int(g.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
//    return g;
//}

std::vector<long long> Rational::factorNumerator() const {
	std::vector<long long> factors;
	//cpp_int num = numerator;
	//int numInt = num.convert_to<int>();
	if (numerator > 0) {
		for (long long i = 1; i <= numerator; ++i) {
			if (numerator % i == 0) factors.push_back(i);
		}
	}
	else {
		for (long long i = -1; i >= numerator; --i) {
			if (numerator % i == 0) factors.push_back(i);
		}
	}
	return factors;
}

bool Rational::isInteger() const {
	return denominator == 1;
}
