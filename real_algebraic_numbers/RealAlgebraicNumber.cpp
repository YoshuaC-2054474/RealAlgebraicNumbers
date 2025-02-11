#include "RealAlgebraicNumber.h"
#include <cmath>
#include <algorithm>
#include <iostream>

using Matrix = std::vector<std::vector<Polynomial>>;

Polynomial determinant(const Matrix& mat) {
	int n = mat.size();
	if (n == 1) return mat[0][0];

	auto det = Polynomial();
	for (int j = 0; j < n; ++j) {
		Matrix minor;
		for (int i = 1; i < n; ++i) {
			std::vector<Polynomial> row;
			for (int k = 0; k < n; ++k)
				if (k != j) row.push_back(mat[i][k]);
			minor.push_back(row);
		}
		Polynomial term = mat[0][j] * determinant(minor);
		if (j % 2) term = term * Polynomial({ -1 });
		det += term;
	}
	return det;
}

int binom(int n, int k) {
	int result = 1;
	if (k > n) return 0;
	k = std::min(k, n - k);
	for (int i = 1; i <= k; i++) {
		result = result * (n - i + 1) / i;
	}
	return result;
}

std::vector<Polynomial> poly_x_minus_y(const Polynomial& f) {
	int n = f.degree;  // degree of f(x)
	// result[j] will be a polynomial in x corresponding to the coefficient of y^j.
	std::vector<Polynomial> result(n + 1);

	// For each power j of y (from 0 to n)
	for (int j = 0; j <= n; j++) {
		// The polynomial coefficient for y^j has degree (n - j) in x.
		std::vector<Rational> poly(n - j + 1, 0);  // poly[k] will be the coefficient for x^k
		// Sum contributions from every term in f(x) that contributes to y^j.
		for (int i = j; i <= n; i++) {
			// The term from a_i * (x-y)^i:
			// In the binomial expansion, the y^j term is:
			//   a_i * (-1)^j * binom(i, j) * x^(i-j) * y^j.
			Rational coeff = f.coefficients[i] * ((j % 2 == 0) ? 1 : -1) * binom(i, j);
			poly[i - j] += coeff;
		}
		result[j] = poly;
	}
	return result;
}

Matrix createSylvesterMatrix(const std::vector<Polynomial>& f_sub, const Polynomial& g) {
	int m = f_sub.size() - 1;
	int n = g.degree;
	int matrixSize = m + n;
	Matrix sylvesterMatrix;

	for (int i = 0; i < n; i++)
	{
		std::vector<Polynomial> row(matrixSize, Polynomial({ 0 }));
		for (int j = 0; j <= m; j++)
		{
			if (i + j < matrixSize)
			{
				row[i + j] = f_sub[j];
			}
		}
		sylvesterMatrix.push_back(row);
	}

	for (int i = 0; i < m; i++)
	{
		std::vector<Polynomial> row(matrixSize, Polynomial({ 0 }));
		for (int j = 0; j <= n; j++)
		{
			if (i + j < matrixSize)
			{
				row[i + j] = (j <= g.degree) ? Polynomial({ g.coefficients[j] }) : Polynomial({ 0 });
			}
		}
		sylvesterMatrix.push_back(row);
	}

	return sylvesterMatrix;
}

RealAlgebraicNumber::RealAlgebraicNumber()
	: polynomial({ 0 }), interval({ 0.0, 0.0 })
{
	//normalize();
}

RealAlgebraicNumber::RealAlgebraicNumber(const Polynomial& polynomial, const Interval& interval)
	: polynomial(polynomial), interval(interval)
{
	//normalize();
}

RealAlgebraicNumber::RealAlgebraicNumber(const Polynomial& polynomial, const Rational lowerBound, const Rational upperBound)
	: polynomial(polynomial), interval{ lowerBound, upperBound }
{
	//normalize();
}

RealAlgebraicNumber::RealAlgebraicNumber(const std::vector<Rational>& coefficients, const Rational lowerBound, const Rational upperBound)
	: polynomial(coefficients), interval{ lowerBound, upperBound }
{
	//normalize();
}

void RealAlgebraicNumber::fromInteger(const int n)
{
	this->polynomial = Polynomial({ 1, -n });
	this->interval.lower_bound = n;
	this->interval.upper_bound = n;
}

RealAlgebraicNumber RealAlgebraicNumber::operator+(const RealAlgebraicNumber& other) const
{
	this->polynomial.print();
	this->interval.lower_bound.print();
	std::cout << " - ";
	this->interval.upper_bound.print();
	std::cout << "\n";
	other.polynomial.print();
	other.interval.lower_bound.print();
	std::cout << " - ";
	other.interval.upper_bound.print();
	std::cout << "\n";

	std::vector<Polynomial> f_sub = poly_x_minus_y(this->polynomial);

	Matrix sylvester = createSylvesterMatrix(f_sub, other.polynomial);

	Polynomial f3 = determinant(sylvester);
	f3.print();

	//Polynomial derivative = f3.derivative();
	//std::vector<double> coDouble(f3.coefficients.begin(), f3.coefficients.end());

	std::vector<Polynomial> sturm;
	if (f3.sturm_sequence.empty()) {
		sturm = f3.sturmSequence(f3);
	}
	else
	{
		sturm = f3.sturm_sequence;
	}

	//for (const auto& sequence : sturm) {
	//	/*for (double coeff : sequence) {
	//		std::cout << coeff << " ";
	//	
	//	}*/
	//	sequence.print();
	//	std::cout << "\n";
	//}

	// Compute initial interval bounds
	Rational l3 = interval.lower_bound + other.interval.lower_bound;
	Rational r3 = interval.upper_bound + other.interval.upper_bound;

	// Refine interval until exactly one root is isolated
	while (variationCount(sturm, l3) - variationCount(sturm, r3) > 1) {
		std::cout << l3 << " - " << r3 << "\n";
		auto current = RealAlgebraicNumber(f3, l3, r3);
		current.refine();
		l3 = current.interval.lower_bound;
		r3 = current.interval.upper_bound;
	}

	const Interval sumInterval = { l3, r3 };
	return { f3, sumInterval };
}

RealAlgebraicNumber RealAlgebraicNumber::operator-(const RealAlgebraicNumber& other) const
{
	return *this + (-other);
}

RealAlgebraicNumber RealAlgebraicNumber::operator-() const
{
	return RealAlgebraicNumber(polynomial.reflectY(), -interval.upper_bound, -interval.lower_bound);
}

int RealAlgebraicNumber::countSignVariations(const std::vector<Rational>& sequence) const {
	int variations = 0;
	int prevSign = 0;

	for (Rational value : sequence) {
		if (value != 0) {
			int currentSign = (value > 0) ? 1 : -1;
			if (prevSign != 0 && currentSign != prevSign) {
				variations++;
			}
			prevSign = currentSign;
		}
	}
	return variations;
}

Rational RealAlgebraicNumber::evaluatePoly(const std::vector<Rational>& sequence, const Rational& x) const {
	Rational result = 0;
	for (int i = sequence.size() - 1; i >= 0; i--) {
		result = result * x + sequence[i];
	}
	/*for (const Rational coeff : sequence) {
		result = result * x + coeff;
	}*/
	return result;
}

int RealAlgebraicNumber::variationCount(const std::vector<Polynomial>& sturm, const Rational& x) const {
	std::vector<Rational> evaluations;

	for (const auto& sequence : sturm) {
		evaluations.push_back(evaluatePoly(sequence.coefficients, x));
	}

	return countSignVariations(evaluations);
}

void RealAlgebraicNumber::normalize()
{
	Rational fInf = 0;
	for (const Rational coeff : polynomial.coefficients) {
		if (fInf < coeff.abs())
			fInf = coeff.abs();
		//fInf = std::max(fInf, coeff.abs());
	}

	const Rational p = Rational(1) / (Rational(1) + fInf);

	//if (polynomial.sturm_sequence.empty()) {
	//	Polynomial derivative = polynomial.derivative();
	//	//polynomial.sturm(derivative);
	//}

	std::vector<Polynomial> sturm;
	if (polynomial.sturm_sequence.empty()) {
		sturm = polynomial.sturmSequence(polynomial);
	}
	else
	{
		sturm = polynomial.sturm_sequence;
	}

	const int varL = variationCount(sturm, interval.lower_bound);
	const int varNegP = variationCount(sturm, -p);
	const int varP = variationCount(sturm, p);
	//int varR = variationCount(polynomial.coefficients, interval.upper_bound);

	if (varL > varNegP) {
		interval.upper_bound = -p; // Adjust interval to exclude zero
	}
	else if (varNegP > varP) {
		interval.lower_bound = 0;
		interval.upper_bound = 0; // Indicating zero root
	}
	else {
		interval.lower_bound = p; // Adjust to keep positive sign interval
	}
}

void RealAlgebraicNumber::refine()
{
	std::vector<Polynomial> sturm;
	if (polynomial.sturm_sequence.empty()) {
		sturm = polynomial.sturmSequence(polynomial);
	}
	else
	{
		sturm = polynomial.sturm_sequence;
	}

	const Rational m = (interval.lower_bound + interval.upper_bound) / 2;

	const int varL = variationCount(sturm, interval.lower_bound);
	const int varM = variationCount(sturm, m);
	if (varL > varM)
	{
		interval.upper_bound = m;
	}
	else
	{
		interval.lower_bound = m;
	}
}

std::string RealAlgebraicNumber::toString() const
{
	/*Rational prevLower = interval.lower_bound;
	Rational prevUpper = interval.upper_bound;
	while (true) {
		RealAlgebraicNumber current = {polynomial, prevLower, prevUpper};
		current.refine();
		if (current.interval.lower_bound == prevLower && current.interval.upper_bound == prevUpper) {
			break;
		}
		prevLower = current.interval.lower_bound;
		prevUpper = current.interval.upper_bound;
	}*/
	std::string output;
	output += "Real Algebraic Number:\n";
	output += polynomial.toString() + "\n";
	output += interval.lower_bound.toString() + " <= x <= " + interval.upper_bound.toString();

	return output;
}