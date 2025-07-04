#ifndef RATIONAL_H
#define RATIONAL_H

#include <iostream>
#include <string>
#include <boost/multiprecision/cpp_int.hpp>
using namespace boost::multiprecision;

class Rational {
public:
	Rational() : numerator(0), denominator(1) {}
	Rational(const cpp_int& num, const cpp_int& den);
	Rational(int num, int den);
	Rational(const cpp_int& numerator);
	Rational(int numerator);
	Rational(long long numerator);
	Rational(double numerator, const cpp_int& maxDenominator = 1000000);
	Rational(float numerator);
	Rational(const Rational& other);

	std::string toString() const;
	std::string toDecimalString(int precision) const;
	void print() const;

	int toInt() const {
		return numerator.convert_to<int>();
		//return { numerator.convert_to<int>(), denominator.convert_to<int>() };
	}

	/*friend std::ostream& operator<<(std::ostream& os, const Rational& rational) {
	    os << rational.numerator << "/" << rational.denominator;
	    return os;
	}*/

	explicit operator double() const {
		return numerator.convert_to<double>() / denominator.convert_to<double>();
	}

	/*Rational operator+(const Rational& other) const;
	Rational operator-(const Rational& other) const;
	friend Rational operator-(const cpp_int& lsh, const Rational& other);
	
	Rational operator*(const Rational& other) const;
	friend Rational operator*(const cpp_int& lsh, const Rational& other);
	Rational operator/(const Rational& other) const;
	friend Rational operator/(const cpp_int& lsh, const Rational& other);*/
	Rational operator%(const Rational& other) const;

	Rational& operator+=(const Rational& other);
	Rational& operator-=(const Rational& other);
	Rational& operator*=(const Rational& other);
	Rational& operator/=(const Rational& other);

	friend Rational operator+(const Rational& lhs, const Rational& rhs);
	friend Rational operator-(const Rational& lhs, const Rational& rhs);
	Rational operator-() const;
	friend Rational operator*(const Rational& lhs, const Rational& rhs);
	friend Rational operator/(const Rational& lhs, const Rational& rhs);

	bool operator==(const Rational& other) const;
	bool operator!=(const Rational& other) const;
	bool operator<(const Rational& other) const;
	bool operator>(const Rational& other) const;
	bool operator<=(const Rational& other) const;
	bool operator>=(const Rational& other) const;

	Rational abs() const;
	Rational inverse() const;
	double sqrt(int n = 2) const;
	Rational pow(int n = 2) const;
	Rational gcd(const Rational& other) const;
	std::vector<cpp_int> factorNumerator() const;
	bool isInteger() const;

private:
	cpp_int numerator;
	cpp_int denominator;

	void simplify();

	//static cpp_int computeGcd(const cpp_int a, const cpp_int b);
};

#endif
