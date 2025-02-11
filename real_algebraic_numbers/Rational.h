#ifndef RATIONAL_H
#define RATIONAL_H

#include <numeric>
#include <algorithm>
#include <iostream>
#include <string>
#include <stdexcept>

class Rational
{
public:
	long long numerator;
	long long denominator;

	Rational() : numerator(0), denominator(1) {}
	Rational(int numerator, int denominator)
	{
		if (denominator == 0) {
			throw std::invalid_argument("Denominator cannot be zero");
		}
		int g = std::gcd(numerator, denominator);
		this->numerator = numerator / g;
		this->denominator = denominator / g;
	}
	Rational(long long numerator, long long denominator)
	{
		if (denominator == 0) {
			throw std::invalid_argument("Denominator cannot be zero");
		}
		int g = std::gcd(numerator, denominator);
		this->numerator = numerator / g;
		this->denominator = denominator / g;
	}
	Rational(int numerator)
	{
		this->numerator = numerator;
		this->denominator = 1;
	}
	Rational(long long numerator)
	{
		this->numerator = numerator;
		this->denominator = 1;
	}
	Rational(double numer)
	{
		int denom = 1;
		while (numer != std::floor(numer))
		{
			numer *= 10;
			denom *= 10;
		}
		int numerInt = static_cast<int>(std::floor(numer));
		int g = std::gcd(numerInt, denom);
		this->numerator = numerInt / g;
		this->denominator = denom / g;

	}
	Rational(const Rational& other)
	{
		this->numerator = other.numerator;
		this->denominator = other.denominator;
	}

	std::string toString() const
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

	void print() const
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

	friend std::ostream& operator<<(std::ostream& os, const Rational& rational) {
		os << rational.numerator << "/" << rational.denominator;
		return os;
	}

	explicit operator double() const
	{
		return static_cast<double>(numerator) / denominator;
	}

	Rational operator+(const Rational& other) const
	{
		long long num = static_cast<long long>(numerator) * other.denominator + static_cast<long long>(other.numerator) * denominator;
		long long denom = static_cast<long long>(denominator) * other.denominator;
		/*int g = std::gcd(num, denom);
		num /= g;
		denom /= g;*/
		return Rational(num, denom);
	}
	Rational operator+(const int other) const
	{
		return *this + other;
	}
	Rational operator-(const Rational& other) const
	{
		long long num = static_cast<long long>(numerator) * other.denominator - static_cast<long long>(other.numerator) * denominator;
		long long denom = static_cast<long long>(denominator) * other.denominator;
		/*int g = std::gcd(num, denom);
		num /= g;
		denom /= g;*/
		return Rational(num, denom);
	}
	Rational operator-() const
	{
		return Rational(-numerator, denominator);
	}
	Rational operator*(const Rational& other) const
	{
		long long num = static_cast<long long>(numerator) * other.numerator;
		long long denom = static_cast<long long>(denominator) * other.denominator;
		/*int g = std::gcd(num, denom);
		num /= g;
		denom /= g;*/
		return Rational(num, denom);
	}
	Rational operator/(const Rational& other) const
	{
		long long num = static_cast<long long>(numerator) * other.denominator;
		long long denom = static_cast<long long>(denominator) * other.numerator;
		/*int g = std::gcd(num, denom);
		num /= g;
		denom /= g;*/
		return Rational(num, denom);
	}

	Rational operator+=(const Rational& other)
	{
		*this = *this + other;
		return *this;
	}
	Rational operator-=(const Rational& other)
	{
		*this = *this - other;
		return *this;
	}
	Rational operator*=(const Rational& other)
	{
		*this = *this * other;
		return *this;
	}
	Rational operator/=(const Rational& other)
	{
		*this = *this / other;
		return *this;
	}

	bool operator==(const Rational& other) const
	{
		long long left = static_cast<long long>(numerator) * other.denominator;
		long long right = static_cast<long long>(other.numerator) * denominator;
		return left == right;
	}
	bool operator!=(const Rational& other) const
	{
		long long left = static_cast<long long>(numerator) * other.denominator;
		long long right = static_cast<long long>(other.numerator) * denominator;
		return left != right;
	}
	bool operator<(const Rational& other) const
	{
		long long left = static_cast<long long>(numerator) * other.denominator;
		long long right = static_cast<long long>(other.numerator) * denominator;
		return left < right;
	}
	bool operator>(const Rational& other) const
	{
		long long left = static_cast<long long>(numerator) * other.denominator;
		long long right = static_cast<long long>(other.numerator) * denominator;
		return left > right;
	}
	bool operator<=(const Rational& other) const
	{
		long long left = static_cast<long long>(numerator) * other.denominator;
		long long right = static_cast<long long>(other.numerator) * denominator;
		return left <= right;
	}
	bool operator>=(const Rational& other) const
	{
		long long left = static_cast<long long>(numerator) * other.denominator;
		long long right = static_cast<long long>(other.numerator) * denominator;
		return left >= right;
	}

	Rational abs() const
	{
		return Rational(std::abs(numerator), std::abs(denominator));
	}
};


#endif