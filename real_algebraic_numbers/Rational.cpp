#include "Rational.h"
#include <stdexcept>
#include <numeric>

Rational::Rational(const int numerator, const int denominator)
{
	if (denominator == 0) {
		throw std::invalid_argument("Denominator cannot be zero");
	}
	const int g = std::gcd(numerator, denominator);
	this->numerator = numerator / g;
	this->denominator = denominator / g;
}

Rational::Rational(const long long numerator, const long long denominator)
{
	if (denominator == 0) {
		throw std::invalid_argument("Denominator cannot be zero");
	}
	const long long g = std::gcd(numerator, denominator);
	this->numerator = numerator / g;
	this->denominator = denominator / g;
}

Rational::Rational(const int numerator)
{
	this->numerator = numerator;
	this->denominator = 1;
}

Rational::Rational(const long long numerator)
{
	this->numerator = numerator;
	this->denominator = 1;
}

Rational::Rational(double numer)
{
	int denom = 1;
	while (numer != std::floor(numer))
	{
		numer *= 10;
		denom *= 10;
	}
	const int numerInt = static_cast<int>(std::floor(numer));
	const int g = std::gcd(numerInt, denom);
	this->numerator = numerInt / g;
	this->denominator = denom / g;
}

Rational::Rational(float numer)
{
	int denom = 1;
	while (numer != std::floor(numer))
	{
		numer *= 10;
		denom *= 10;
	}
	const int numerInt = static_cast<int>(std::floor(numer));
	const int g = std::gcd(numerInt, denom);
	this->numerator = numerInt / g;
	this->denominator = denom / g;
}

Rational::Rational(const Rational& other)
{
	this->numerator = other.numerator;
	this->denominator = other.denominator;
}

std::string Rational::toString() const
{
	if (denominator == 1)
	{
		return std::to_string(numerator);
	}
	else
	{
		return std::to_string(numerator) + "/" + std::to_string(denominator);
	}
}

void Rational::print() const
{
	if (denominator == 1)
	{
		std::cout << numerator;
	}
	else
	{
		std::cout << numerator << "/" << denominator;
	}
}

bool willAdditionOverflow(const Rational& a, const Rational& b) {
	if (b > 0 && a > std::numeric_limits<long long>::max() - b)
		return true;
	if (b < 0 && a < std::numeric_limits<long long>::min() - b)
		return true;
	return false;
}

Rational Rational::operator+(const Rational& other) const
{
	if (willAdditionOverflow(*this, other))
		throw std::overflow_error("Addition overflow error");

	if (other.numerator == 0 && other.denominator == 1)
	{
		return *this;
	}
	const long long num = static_cast<long long>(numerator) * other.denominator + static_cast<long long>(other.numerator) * denominator;
	const long long denom = static_cast<long long>(denominator) * other.denominator;
	return { num, denom };
}

Rational Rational::operator-(const Rational& other) const
{
	if (other.numerator == 0 && other.denominator == 1)
	{
		return *this;
	}
	const long long num = static_cast<long long>(numerator) * other.denominator - static_cast<long long>(other.numerator) * denominator;
	const long long denom = static_cast<long long>(denominator) * other.denominator;
	return { num, denom };
}

Rational operator-(const long long lsh, const Rational& other)
{
	if (other.numerator == 0 && other.denominator == 1)
	{
		return lsh;
	}
	const Rational temp = lsh;
	const long long num = static_cast<long long>(temp.numerator) * other.denominator - static_cast<long long>(other.numerator) * temp.denominator;
	const long long denom = static_cast<long long>(temp.denominator) * other.denominator;
	return { num, denom };
}

Rational Rational::operator-() const
{
	return { -numerator, denominator };
}

bool willMultiplicationOverflow(const Rational& a, const Rational& b) {
	// Handle zero cases
	if (a == 0 || b == 0) return false;
	// Check overflow depending on the signs
	if (a > 0) {
		if (b > 0)
			return a > std::numeric_limits<long long>::max() / b;
		else
			return b < std::numeric_limits<long long>::min() / a;
	}
	else {
		if (b > 0)
			return a < std::numeric_limits<long long>::min() / b;
		else
			return a != 0 && b < std::numeric_limits<long long>::max() / a;
	}
}

Rational Rational::operator*(const Rational& other) const
{
	if (willMultiplicationOverflow(*this, other))
		throw std::overflow_error("Multiplication overflow error");

	if (other.numerator == 0 && other.denominator == 1) return { 0, 1 };
	if (other.numerator == 1 && other.denominator == 1) return *this;
	if (other.numerator == -1 && other.denominator == 1) return -(*this);

	const long long num = static_cast<long long>(numerator) * other.numerator;
	const long long denom = static_cast<long long>(denominator) * other.denominator;
	return { num, denom };
}

Rational Rational::operator/(const Rational& other) const
{
	if (other.numerator == 0 && other.denominator == 1) throw std::invalid_argument("Cannot divide by zero");
	if (other.numerator == 1 && other.denominator == 1) return *this;
	if (other.numerator == -1 && other.denominator == 1) return -(*this);

	const long long num = static_cast<long long>(numerator) * other.denominator;
	const long long denom = static_cast<long long>(denominator) * other.numerator;
	return { num, denom };
}

Rational operator/(const long long lsh, const Rational& other)
{
	if (other.numerator == 0 && other.denominator == 1) throw std::invalid_argument("Cannot divide by zero");
	if (other.numerator == 1 && other.denominator == 1) return lsh;
	if (other.numerator == -1 && other.denominator == 1) return -lsh;

	const Rational temp = lsh;
	const long long num = static_cast<long long>(temp.numerator) * other.denominator;
	const long long denom = static_cast<long long>(temp.denominator) * other.numerator;
	return { num, denom };
}

//Rational Rational::operator%(const Rational& other) const
//{
//	long long numGcd = std::gcd(numerator, other.numerator);
//	long long denomGcd = denominator * other.denominator;
//	return { numGcd, denomGcd };
//}

Rational Rational::operator+=(const Rational& other)
{
	*this = *this + other;
	return *this;
}

Rational Rational::operator-=(const Rational& other)
{
	*this = *this - other;
	return *this;
}

Rational Rational::operator*=(const Rational& other)
{
	*this = *this * other;
	return *this;
}

Rational Rational::operator/=(const Rational& other)
{
	*this = *this / other;
	return *this;
}

bool Rational::operator==(const Rational& other) const
{
	const long long left = static_cast<long long>(numerator) * other.denominator;
	const long long right = static_cast<long long>(other.numerator) * denominator;
	return left == right;
}

bool Rational::operator!=(const Rational& other) const
{
	const long long left = static_cast<long long>(numerator) * other.denominator;
	const long long right = static_cast<long long>(other.numerator) * denominator;
	return left != right;
}

bool Rational::operator<(const Rational& other) const
{
	const long long left = static_cast<long long>(numerator) * other.denominator;
	const long long right = static_cast<long long>(other.numerator) * denominator;
	return left < right;
}

bool Rational::operator>(const Rational& other) const
{
	const long long left = static_cast<long long>(numerator) * other.denominator;
	const long long right = static_cast<long long>(other.numerator) * denominator;
	return left > right;
}

bool Rational::operator<=(const Rational& other) const
{
	const long long left = static_cast<long long>(numerator) * other.denominator;
	const long long right = static_cast<long long>(other.numerator) * denominator;
	return left <= right;
}

bool Rational::operator>=(const Rational& other) const
{
	const long long left = static_cast<long long>(numerator) * other.denominator;
	const long long right = static_cast<long long>(other.numerator) * denominator;
	return left >= right;
}

Rational Rational::abs() const
{
	return { std::abs(numerator), std::abs(denominator) };
}

Rational Rational::inverse() const
{
	if (numerator == 0)
	{
		throw std::invalid_argument("Cannot invert zero");
	}
	else if (numerator < 0)
	{
		return { -denominator, -numerator };
	}
	else
	{
		return { denominator, numerator };
	}
}

Rational Rational::sqrt(const int n) const
{
	if (n == 0)
	{
		throw std::invalid_argument("Cannot take zeroth root");
	}
	if (n % 2 == 0 && numerator < 0)
	{
		throw std::invalid_argument("Cannot take even root of negative number");
	}

	float result = static_cast<float>(numerator) / denominator;
	result = std::pow(result, 1.0 / n);

	return result;

	//long long numSqrt = static_cast<long long>(std::pow(std::abs(numerator), 1.0 / n));
	//long long denomSqrt = static_cast<long long>(std::pow(std::abs(denominator), 1.0 / n));

	//if (n % 2 != 0 && numerator < 0)
	//{
	//	return { -numSqrt, denomSqrt };
	//}

	////long long numSqrt = static_cast<long long>(std::pow(numerator, 1.0 / n));
	//return { numSqrt, denomSqrt };
}

Rational Rational::gcd(const Rational& other) const {
	long long numGcd = std::gcd(numerator, other.numerator);
	long long denomGcd = denominator * other.denominator;
	return { numGcd, denomGcd };
}