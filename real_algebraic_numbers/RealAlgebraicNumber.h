#ifndef REAL_ALGEBRAIC_NUMBER_H
#define REAL_ALGEBRAIC_NUMBER_H

#include <vector>
#include <string>
#include "MyMatrix.h"
#include "Polynomial.h"

struct Interval {
	Rational lowerBound;
	Rational upperBound;
};

class RealAlgebraicNumber {
public:
	Polynomial polynomial;
	Interval interval = {.lowerBound = 0, .upperBound = 0};

	RealAlgebraicNumber();
	RealAlgebraicNumber(const Polynomial& polynomial, const Interval& interval);
	RealAlgebraicNumber(const Polynomial& polynomial, const Rational& lowerBound, const Rational& upperBound);
	RealAlgebraicNumber(const std::vector<Rational>& coefficients, const Rational& lowerBound,
	                    const Rational& upperBound);
	RealAlgebraicNumber(int value);
	RealAlgebraicNumber(const Rational& value);
	RealAlgebraicNumber(double value);

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
	RealAlgebraicNumber sqrt(const Rational& n = 2) const;
	RealAlgebraicNumber pow(const Rational& n) const;

	bool isZero() const;
	bool isPositive() const;
	bool isNegative() const;

	friend std::ostream& operator<<(std::ostream& os, const RealAlgebraicNumber& ran);
	std::string toString() const;
	std::string toDecimalString(int precision = 10) const;

private:
	friend long long binomialCoeff(int n, int k);
	friend MyMatrix<Polynomial> constructSylvesterMatrixForSum(const Polynomial& p, const Polynomial& q);
	friend MyMatrix<Polynomial> constructSylvesterMatrixForProduct(const Polynomial& p, const Polynomial& q);
	friend MyMatrix<Polynomial> constructSylvesterMatrixForPower(const Polynomial& p, int k);
	friend Rational minRational(const Rational& r1, const Rational& r2, const Rational& r3, const Rational& r4);
	friend Rational maxRational(const Rational& r1, const Rational& r2, const Rational& r3, const Rational& r4);

	int intervalToOrder();

	// Function to count variations in Sturm sequence at x
	static int variationCount(const std::vector<Polynomial>& sturm, const Rational& x);

	void normalize();
	void refine();
	void refineToTolerance();
};

#endif
