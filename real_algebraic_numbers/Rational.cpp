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

Rational Rational::operator+(const Rational& other) const
{
	const long long num = static_cast<long long>(numerator) * other.denominator + static_cast<long long>(other.numerator) * denominator;
	const long long denom = static_cast<long long>(denominator) * other.denominator;
	return {num, denom};
}

Rational Rational::operator-(const Rational& other) const
{
	const long long num = static_cast<long long>(numerator) * other.denominator - static_cast<long long>(other.numerator) * denominator;
	const long long denom = static_cast<long long>(denominator) * other.denominator;
	return {num, denom};
}

Rational Rational::operator-() const
{
	return { -numerator, denominator };
}

Rational Rational::operator*(const Rational& other) const
{
	const long long num = static_cast<long long>(numerator) * other.numerator;
	const long long denom = static_cast<long long>(denominator) * other.denominator;
	return { num, denom };
}

Rational Rational::operator/(const Rational& other) const
{
	const long long num = static_cast<long long>(numerator) * other.denominator;
	const long long denom = static_cast<long long>(denominator) * other.numerator;
	return { num, denom };
}

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