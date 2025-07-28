#include "Polynomial.h"
#include <iostream>
#include <numeric>
#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <NTL/ZZX.h>
#include <NTL/ZZXFactoring.h>
#include "MyTimer.h"

// Use specific NTL functions to avoid namespace conflicts
using NTL::ZZX;
using NTL::ZZ;
using NTL::Vec;
using NTL::Pair;
using NTL::SetCoeff;
using NTL::to_ZZ;
using NTL::to_int;
using NTL::factor;

Polynomial::Polynomial(const std::initializer_list<int> coeffs) {
	coefficients.reserve(coeffs.size());
	for (const int coeff : coeffs) {
		coefficients.emplace_back(coeff);
	}

	degree = coefficients.empty() ? -1 : static_cast<int>(coefficients.size()) - 1;
}

Polynomial::Polynomial(const std::vector<Rational>& coeffs) : coefficients(coeffs) {
	degree = coefficients.empty() ? -1 : static_cast<int>(coefficients.size()) - 1;
}

Polynomial::Polynomial(const int zeroCoeff) {
	const std::vector<Rational> coeffs = {zeroCoeff};
	coefficients = coeffs;

	degree = coefficients.empty() ? -1 : static_cast<int>(coefficients.size()) - 1;
}

Rational findGcd(const std::vector<Rational>& arr) {
	if (arr.empty()) return {1};

	Rational res = arr[0];
	for (size_t i = 1; i < arr.size(); ++i) {
		res = boost::multiprecision::gcd(arr[i].numerator(), res.numerator()) / boost::multiprecision::lcm(
			arr[i].denominator(), res.denominator());
		if (res == 1) return {1};
	}
	return res;
}

std::vector<Polynomial> minimalPolynomialsNtl(const Polynomial& poly) {
	PROFILE_FUNCTION
	ZZX f;
	f.SetMaxLength(static_cast<long>(poly.coefficients.size()));

	for (long i = 0; i < static_cast<long>(poly.coefficients.size()); ++i) {
		SetCoeff(f, i, to_ZZ(boost::rational_cast<int>(poly.coefficients[i])));
	}

	Vec<Pair<ZZX, long>> factors;
	ZZ content;
	factor(content, factors, f);

	std::vector<Polynomial> result;
	result.reserve(factors.length());
	for (long i = 0; i < factors.length(); ++i) {
		const ZZX& factor = factors[i].a;
		std::vector<Rational> coeffs;
		coeffs.reserve(factor.rep.length());
		for (long j = 0; j < factor.rep.length(); ++j) {
			coeffs.emplace_back(to_int(factor.rep[j]));
		}
		result.emplace_back(std::move(coeffs));
	}

	return result;
}

void Polynomial::normalize(const Rational& lowerBound, const Rational& upperBound) {
	PROFILE_FUNCTION
	if (is_normalized) return;

	// Remove leading zeros
	while (coefficients.size() > 1 && coefficients.back() == 0) {
		coefficients.pop_back();
	}
	degree = coefficients.empty() ? -1 : static_cast<int>(coefficients.size()) - 1;

	if (coefficients.empty()) return;

	// turn x^2 + 3/2x into 2x^2 + 3x
	std::vector<cpp_int> allDenominators;
	for (const Rational& coeff : coefficients) {
		/*if (coeff.denominator != 1) {
			allDenominators.push_back(coeff.denominator);
		}*/
		if (coeff.denominator() != 1) {
			allDenominators.push_back(coeff.denominator());
		}
	}
	for (const cpp_int& denom : allDenominators) {
		for (Rational& coeff : coefficients) {
			coeff *= denom; // Multiply all coefficients by the denominator
		}
	}

	// Divide by GCD
	if (const Rational gcd = findGcd(coefficients); gcd > 1) {
		for (Rational& coeff : coefficients) {
			coeff /= gcd;
		}
	}

	std::vector<Polynomial> minimalPolys = minimalPolynomialsNtl(*this);

	// Find polynomial with root in the interval
	for (const Polynomial& poly : minimalPolys) {
		const bool firstSign = poly.evaluate(lowerBound) < 0;
		const bool secondSign = poly.evaluate(upperBound) < 0;
		if (firstSign != secondSign) {
			*this = poly;
			break;
		}
	}

	// Make leading coefficient positive
	if (!coefficients.empty() && coefficients.back() < 0) {
		for (Rational& coeff : coefficients) {
			coeff *= -1;
		}
	}

	is_normalized = true;
}

bool Polynomial::isZero() const {
	return coefficients.empty() || (coefficients.size() == 1 && coefficients[0] == 0);
}

Polynomial Polynomial::operator+(const Polynomial& other) const {
	const int maxDeg = (this->degree > other.degree) ? this->degree : other.degree;
	std::vector<Rational> resultCoeffs;
	resultCoeffs.reserve(maxDeg + 1);

	for (int i = 0; i <= maxDeg; ++i) {
		resultCoeffs.emplace_back(this->coeff(i) + other.coeff(i));
	}
	return {resultCoeffs};
}

Polynomial Polynomial::operator-(const Polynomial& other) const {
	PROFILE_FUNCTION
	const int maxDeg = (this->degree > other.degree) ? this->degree : other.degree;
	std::vector<Rational> resultCoeffs;
	resultCoeffs.reserve(maxDeg + 1);

	for (int i = 0; i <= maxDeg; ++i) {
		resultCoeffs.emplace_back(this->coeff(i) - other.coeff(i));
	}
	return {resultCoeffs};
}

Polynomial Polynomial::operator*(const Polynomial& other) const {
	PROFILE_FUNCTION
	if (this->isZero() || other.isZero()) {
		return {0};
	}

	const int newDeg = this->degree + other.degree;
	std::vector<Rational> resultCoeffs(newDeg + 1, 0);

	for (int i = 0; i <= this->degree; ++i) {
		if (this->coeff(i) == 0) continue; // Skip zero coefficients
		for (int j = 0; j <= other.degree; ++j) {
			resultCoeffs[i + j] += this->coeff(i) * other.coeff(j);
		}
	}
	return {resultCoeffs};
}

Polynomial Polynomial::operator/(const Polynomial& other) const {
	PROFILE_FUNCTION
	if (other.isZero()) {
		throw std::runtime_error("Division by zero polynomial!");
	}
	if (this->isZero()) {
		return {0};
	}

	auto [quotient, remainder] = polyDivide(*this, other);
	return polyTrim(Polynomial(quotient));
}

Polynomial Polynomial::operator/(const Rational& scalar) const {
	if (scalar == 0) {
		throw std::runtime_error("Division by zero scalar in polynomial.");
	}

	std::vector<Rational> resultCoeffs;
	resultCoeffs.reserve(coefficients.size());
	for (const Rational& coeff : coefficients) {
		resultCoeffs.emplace_back(coeff / scalar);
	}
	return {resultCoeffs};
}

bool Polynomial::operator==(const Polynomial& other) const {
	if (degree != other.degree) return false;

	return std::equal(coefficients.begin(), coefficients.end(), other.coefficients.begin());
}

std::string Polynomial::toString() const {
	PROFILE_FUNCTION
	if (coefficients.empty()) return "0";

	std::string output;
	output.reserve(coefficients.size() * 10); // Reserve space to avoid reallocations

	for (int i = static_cast<int>(coefficients.size()) - 1; i >= 0; --i) {
		if (coefficients[i] == 0) continue;

		if (coefficients[i] > 0 && !output.empty()) output += " + ";
		else if (coefficients[i] < 0) output += " - ";

		//const Rational absCoeff = coefficients[i].abs();
		const Rational absCoeff = abs(coefficients[i]);
		if (absCoeff != 1 || i == 0) {
			//output += absCoeff.toString();
			output += absCoeff.numerator().str(10, 10) + "/" + absCoeff.denominator().str(10, 10);
		}

		if (i > 1) {
			output += "x^" + std::to_string(i);
		}
		else if (i == 1) {
			output += "x";
		}
	}
	return output;
}

void Polynomial::print() const {
	std::cout << toString() << "\n";
}

Rational Polynomial::evaluate(const Rational& x) const {
	PROFILE_FUNCTION
	Rational result = 0;
	const int coeffLength = static_cast<int>(coefficients.size());
	for (int i = coeffLength - 1; i >= 0; --i) {
		result = result * x + coefficients[i];
	}
	return result;
}

Polynomial Polynomial::derivative() const {
	if (isZero() || degree <= 0) return {};

	std::vector<Rational> result;
	result.reserve(degree);

	for (int i = 1; i <= degree; ++i) {
		result.emplace_back(coefficients[i] * i);
	}
	return {result};
}

Polynomial Polynomial::polyTrim(const Polynomial& poly) {
	PROFILE_FUNCTION
	Polynomial result = poly;

	// Remove trailing zeros with a small epsilon tolerance
	/*while (result.coefficients.size() > 1 && result.coefficients.back().abs() < 1e-9) {
		result.coefficients.pop_back();
	}*/
	while (result.coefficients.size() > 1 && abs(result.coefficients.back()) < cpp_int(1e-9)) {
		result.coefficients.pop_back();
	}

	result.degree = result.coefficients.empty() ? -1 : static_cast<int>(result.coefficients.size()) - 1;
	return result;
}

std::pair<std::vector<Rational>, std::vector<Rational>> Polynomial::polyDivide(
	const Polynomial& dividend, const Polynomial& divisor) {
	PROFILE_FUNCTION
	const Polynomial a = polyTrim(dividend);
	const Polynomial b = polyTrim(divisor);

	if (b.coefficients.empty()) {
		throw std::runtime_error("Division by zero polynomial!");
	}

	if (a.degree < b.degree) {
		return {std::vector<Rational>(), a.coefficients};
	}

	const int degA = a.degree;
	const int degB = b.degree;

	std::vector<Rational> quotient(degA - degB + 1, Rational(0));
	std::vector<Rational> remainder = a.coefficients;

	// Perform polynomial long division
	while (static_cast<int>(remainder.size()) - 1 >= degB && !remainder.empty()) {
		const int degR = static_cast<int>(remainder.size()) - 1;
		const Rational factor = remainder.back() / b.coefficients.back();
		const int power = degR - degB;
		quotient[power] = factor;

		// Subtract factor * (B shifted by "power") from remainder
		for (int i = 0; i <= degB; ++i) {
			remainder[power + i] -= factor * b.coefficients[i];
		}

		// Remove leading zero
		/*if (!remainder.empty() && remainder.back().abs() < 1e-9) {
			remainder.pop_back();
		}*/
		if (!remainder.empty() && abs(remainder.back()) < cpp_int(1e-9)) {
			remainder.pop_back();
		}
	}

	return {polyTrim(Polynomial(quotient)).coefficients, polyTrim(Polynomial(remainder)).coefficients};
}

std::vector<Rational> Polynomial::polyNegate(const std::vector<Rational>& poly) {
	std::vector<Rational> neg;
	neg.reserve(poly.size());
	for (const Rational& coeff : poly) {
		neg.emplace_back(-coeff);
	}
	return neg;
}

Polynomial Polynomial::reflectY() const {
	std::vector<Rational> result;
	result.reserve(coefficients.size());

	for (int i = 0; i < static_cast<int>(coefficients.size()); ++i) {
		if (i % 2 == 1) {
			result.emplace_back(-coefficients[i]);
		}
		else {
			result.emplace_back(coefficients[i]);
		}
	}
	return {result};
}

std::vector<Polynomial> Polynomial::sturmSequence(const Polynomial& p) {
	PROFILE_FUNCTION
	std::vector<Polynomial> seq;
	seq.reserve(p.degree + 1); // Reserve space for efficiency

	const Polynomial s0 = polyTrim(p);
	if (s0.coefficients.empty()) {
		throw std::runtime_error("Zero polynomial provided.");
	}
	seq.push_back(s0);

	// S1 = p'(x)
	const Polynomial s1 = s0.derivative();
	seq.push_back(s1);

	// Compute subsequent sequence members until remainder becomes zero
	while (true) {
		const Polynomial& sPrev = seq[seq.size() - 2];
		const Polynomial& sCurr = seq.back();

		if (sCurr.degree == 0) {
			break;
		}

		// Compute polynomial remainder of S_prev divided by S_curr
		auto [quotient, remainder] = polyDivide(sPrev, sCurr);

		// If remainder is zero, stop
		/*if (remainder.empty() || (remainder.size() == 1 && remainder[0].abs() < 1e-9)) {
			break;
		}*/
		if (remainder.empty() || (remainder.size() == 1 && abs(remainder[0]) < cpp_int(1e-9))) {
			break;
		}

		// Next term in Sturm sequence is the negative of the remainder
		std::vector<Rational> next = polyNegate(remainder);
		seq.emplace_back(next);
	}

	sturm_sequence = seq;
	return seq;
}
