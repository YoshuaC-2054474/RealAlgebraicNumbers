#ifndef RATIONAL2_H
#define RATIONAL2_H

#include <iostream>
#include <string>
#include <vector>
//#include <boost/multiprecision/cpp_int.hpp>
//using namespace boost::multiprecision;

class Rational {
public:
	Rational() : numerator(0), denominator(1) {}
	Rational(long long num, long long den);
	Rational(int num, int den);
	Rational(long long numerator);
	Rational(int numerator);
	Rational(double numer, long long maxDenominator = 1000000);
	Rational(float numer);
	Rational(const Rational& other);

	std::string toString() const;
	void print() const;

	int toInt() const {
		return static_cast<int>(numerator / denominator);
		//return { numerator.convert_to<int>(), denominator.convert_to<int>() };
	}

	/*friend std::ostream& operator<<(std::ostream& os, const Rational& rational) {
		os << rational.numerator << "/" << rational.denominator;
		return os;
	}*/

	explicit operator double() const {
		return static_cast<double>(numerator) / static_cast<double>(denominator);
		//return numerator.convert_to<double>() / denominator.convert_to<double>();
	}

	Rational operator+(const Rational& other) const;
	Rational operator-(const Rational& other) const;
	friend Rational operator-(long long lsh, const Rational& other);
	Rational operator-() const;
	Rational operator*(const Rational& other) const;
	friend Rational operator*(long long lsh, const Rational& other);
	Rational operator/(const Rational& other) const;
	friend Rational operator/(long long lsh, const Rational& other);
	Rational operator%(const Rational& other) const;

	Rational& operator+=(const Rational& other);
	Rational& operator-=(const Rational& other);
	Rational& operator*=(const Rational& other);
	Rational& operator/=(const Rational& other);

	bool operator==(const Rational& other) const;
	bool operator!=(const Rational& other) const;
	bool operator<(const Rational& other) const;
	bool operator>(const Rational& other) const;
	bool operator<=(const Rational& other) const;
	bool operator>=(const Rational& other) const;

	Rational abs() const;
	Rational inverse() const;
	double sqrt(int n = 2) const;
	Rational gcd(const Rational& other) const;
	std::vector<long long> factorNumerator() const;
	bool isInteger() const;

private:
	long long numerator;
	long long denominator;

	//static cpp_int computeGcd(const cpp_int a, const cpp_int b);
};

#endif
