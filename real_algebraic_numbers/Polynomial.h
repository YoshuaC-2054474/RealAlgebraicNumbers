#ifndef POLYNOMIAL3_H
#define POLYNOMIAL3_H

#include "Rational.h"
#include <boost/rational.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <string>
#include <vector>

//using namespace boost::multiprecision;
//using Rational = boost::rational<cpp_int>;


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

	Polynomial(const Rational& c) {
		coefficients.push_back(c);
		degree = 0;
		/*trim_zeros();*/
	}

	void normalize(const Rational& lowerBound, const Rational& upperBound);

	// Check if the polynomial is zero
	bool isZero() const;

	Rational coeff(const int idx) const {
		if (idx < 0 || idx >= static_cast<int>(coefficients.size())) {
			return 0;
		}
		return coefficients[idx];
	}

	Polynomial operator+(const Polynomial& other) const;
	Polynomial operator-(const Polynomial& other) const;
	Polynomial operator*(const Polynomial& other) const;
	Polynomial operator/(const Polynomial& other) const;
	Polynomial operator/(const Rational& scalar) const;

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

private:
	friend Rational findGcd(const std::vector<Rational>& arr);
	friend std::vector<Polynomial> minimalPolynomialsNtl(const Polynomial& poly);
};

#endif
