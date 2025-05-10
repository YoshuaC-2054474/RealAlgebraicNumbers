#ifndef POLYNOMIAL3_H
#define POLYNOMIAL3_H

#include "Rational.h"
#include <string>
#include <vector>


class Polynomial {
public:
	int degree;
	std::vector<Rational> coefficients;
	std::vector<Polynomial> sturm_sequence;
	bool is_normalized = false;

	Polynomial() : degree(-1), coefficients({}) {}
	Polynomial(std::initializer_list<int> coeffs);
	Polynomial(const std::vector<Rational>& coeffs);
	Polynomial(int zeroCoeff);

	void normalize(const Rational& lowerBound, const Rational& upperBound);
	void testNormalize();

	// Check if the polynomial is zero
	bool isZero() const;

	Polynomial operator+(const Polynomial& other) const;
	Polynomial operator-(const Polynomial& other) const;
	Polynomial operator*(const Polynomial& other) const;
	Polynomial operator/(const Polynomial& other) const;

	Polynomial& operator+=(const Polynomial& other) {
		*this = *this + other;
		return *this;
	}

	Polynomial& operator-=(const Polynomial& other) {
		*this = *this - other;
		return *this;
	}

	Polynomial& operator*=(const Polynomial& other) {
		*this = *this * other;
		return *this;
	}


	Polynomial& operator/=(const Polynomial& other) {
		*this = *this / other;
		return *this;
	}

	Polynomial operator-() const { return *this * Polynomial({-1}); }

	bool operator==(const Polynomial& other) const;
	bool operator!=(const Polynomial& other) const { return !(*this == other); }

	std::string toString() const;
	void print() const;

	Rational evaluate(const Rational& x) const;

	// Compute the derivative
	Polynomial derivative() const;

	static Polynomial polyTrim(const Polynomial& poly);
	static std::pair<std::vector<Rational>, std::vector<Rational>> polyDivide(
		const Polynomial& dividend, const Polynomial& divisor);
	static std::vector<Rational> polyNegate(const std::vector<Rational>& poly);

	Polynomial reflectY() const;

	std::vector<Polynomial> sturmSequence(const Polynomial& p);
};

#endif
