#ifndef REAL_ALGEBRAIC_NUMBER_H
#define REAL_ALGEBRAIC_NUMBER_H

#include <vector>
#include <string>
//#include <algorithm>
//#include <cmath>

#include "Polynomial.h"

struct Interval {
	Rational lower_bound;
	Rational upper_bound;
};

class RealAlgebraicNumber {
public:
	Polynomial polynomial;
	Interval interval = {.lower_bound = 0, .upper_bound = 0};

	RealAlgebraicNumber();
	RealAlgebraicNumber(const Polynomial& polynomial, const Interval& interval);
	RealAlgebraicNumber(const Polynomial& polynomial, const Rational& lowerBound, const Rational& upperBound);
	//RealAlgebraicNumber(const RealAlgebraicNumber& other);
	//RealAlgebraicNumber(const Polynomial& polynomial, int lowerBound, int upperBound);
	//RealAlgebraicNumber(const std::vector<int>& coefficients, const Interval& interval);
	RealAlgebraicNumber(const std::vector<Rational>& coefficients, const Rational& lowerBound,
	                    const Rational& upperBound);
	RealAlgebraicNumber(int value);
	RealAlgebraicNumber(double value);

	//RealAlgebraicNumber fromInteger(const int n);

	//void fromRational(const int p, int q)
	//{
	//	polynomial.coefficients = { q, -p };
	//	polynomial.degree = 1;
	//	interval.lower_bound = static_cast<double>(p) / q - 0.1;
	//	interval.upper_bound = static_cast<double>(p) / q + 0.1;
	//}

	//void fromRadical(const int n, const int p, const int q) // n-th root of p/q
	//{
	//	std::vector<int> coefficients(n + 1, 0);
	//	coefficients[0] = q;
	//	coefficients[n] = -p;
	//	polynomial.coefficients = coefficients;
	//	polynomial.degree = n;
	//	//interval.lower_bound = 0;
	//	//interval.upper_bound = 0;
	//}
	//~RealAlgebraicNumber();

	friend RealAlgebraicNumber operator+(const RealAlgebraicNumber& lhs, const RealAlgebraicNumber& rhs);
	friend RealAlgebraicNumber operator-(const RealAlgebraicNumber& lhs, const RealAlgebraicNumber& rhs);
	RealAlgebraicNumber operator-() const;
	friend RealAlgebraicNumber operator*(const RealAlgebraicNumber& lhs, const RealAlgebraicNumber& rhs);
	friend RealAlgebraicNumber operator/(const RealAlgebraicNumber& lhs, const RealAlgebraicNumber& rhs);

	RealAlgebraicNumber operator+=(const RealAlgebraicNumber& other);
	RealAlgebraicNumber operator-=(const RealAlgebraicNumber& other);
	RealAlgebraicNumber operator*=(const RealAlgebraicNumber& other);
	RealAlgebraicNumber operator/=(const RealAlgebraicNumber& other);

	bool operator==(const RealAlgebraicNumber& other) const;
	bool operator!=(const RealAlgebraicNumber& other) const;
	bool operator<(const RealAlgebraicNumber& other) const;
	bool operator>(const RealAlgebraicNumber& other) const;
	bool operator<=(const RealAlgebraicNumber& other) const;
	bool operator>=(const RealAlgebraicNumber& other) const;

	RealAlgebraicNumber inverse() const;
	RealAlgebraicNumber abs() const;
	RealAlgebraicNumber sqrt(int n = 2);
	RealAlgebraicNumber pow(int n) const;

	bool isZero() const;
	bool isPositive() const;
	bool isNegative() const;
	int sign() const;

	std::vector<int> getRoots() const;
	std::vector<int> getInterval() const;
	std::vector<int> getPolynomial() const;
	int getDegree() const;

	std::string toString() const;
	std::string toDecimalString(int precision = 10) const;

	static void testOperators();

private:
	// Function to count sign variations in a sequence
	static int countSignVariations(const std::vector<Rational>& sequence);

	// Function to evaluate a polynomial at a given x-value
	static Rational evaluatePoly(const std::vector<Rational>& sequence, const Rational& x);

	// Function to count variations in Sturm sequence at x
	static int variationCount(const std::vector<Polynomial>& sturm, const Rational& x);

	void normalize();

	void refine();
};

#endif
