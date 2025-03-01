#ifndef RATIONAL_H
#define RATIONAL_H

#include <iostream>
#include <string>

class Rational
{
public:
	Rational() : numerator(0), denominator(1) {}
	Rational(int numerator, int denominator);
	Rational(long long numerator, long long denominator);
	Rational(int numerator);
	Rational(long long numerator);
	Rational(double numer);
	Rational(float numer);
	Rational(const Rational& other);

	std::string toString() const;
	void print() const;

	friend std::ostream& operator<<(std::ostream& os, const Rational& rational) {
		os << rational.numerator << "/" << rational.denominator;
		return os;
	}

	explicit operator double() const
	{
		return static_cast<double>(numerator) / denominator;
	}

	Rational operator+(const Rational& other) const;
	Rational operator-(const Rational& other) const;
	friend Rational operator-(const long long lsh, const Rational& other);
	Rational operator-() const;
	Rational operator*(const Rational& other) const;
	Rational operator/(const Rational& other) const;
	friend Rational operator/(const long long lsh, const Rational& other);
	/*Rational operator%(const Rational& other) const;*/

	Rational operator+=(const Rational& other);
	Rational operator-=(const Rational& other);
	Rational operator*=(const Rational& other);
	Rational operator/=(const Rational& other);

	bool operator==(const Rational& other) const;
	bool operator!=(const Rational& other) const;
	bool operator<(const Rational& other) const;
	bool operator>(const Rational& other) const;
	bool operator<=(const Rational& other) const;
	bool operator>=(const Rational& other) const;

	Rational abs() const;
	Rational inverse() const;
	Rational sqrt(const int n) const;
	Rational gcd(const Rational& other) const;
private:
	long long numerator;
	long long denominator;
};

#endif