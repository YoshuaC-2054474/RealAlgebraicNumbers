#include "RealAlgebraicNumber.h"
#include <cmath>
#include <algorithm>
#include <iostream>

using Matrix = std::vector<std::vector<Polynomial>>;
using Matrix2 = std::vector<std::vector<std::vector<Rational>>>;

Polynomial determinant(const Matrix& mat) {
	const int n = mat.size();
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

std::vector<Rational> determinant2(const Matrix2& mat) {
	const int n = mat.size();
	if (n == 1) return mat[0][0];

	std::vector<Rational> det(mat[0][0].size(), 0);
	for (int j = 0; j < n; ++j) {
		Matrix2 minor;
		for (int i = 1; i < n; ++i) {
			std::vector<std::vector<Rational>> row;
			for (int k = 0; k < n; ++k)
				if (k != j) row.push_back(mat[i][k]);
			minor.push_back(row);
		}
		std::vector<Rational> term = mat[0][j];
		std::vector<Rational> minorDet = determinant2(minor);
		for (size_t k = 0; k < term.size(); ++k) {
			term[k] *= minorDet[k];
		}
		if (j % 2) {
			for (Rational& coeff : term) {
				coeff = -coeff;
			}
		}
		for (size_t k = 0; k < det.size(); ++k) {
			det[k] += term[k];
		}
	}
	return det;
}

int binom(const int n, int k) {
	int result = 1;
	if (k > n) return 0;
	k = std::min(k, n - k);
	for (int i = 1; i <= k; i++) {
		result = result * (n - i + 1) / i;
	}
	return result;
}

std::vector<Polynomial> poly_x_minus_y(const Polynomial& f) {
	const int n = f.degree;  // degree of f(x)
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
			const Rational coeff = f.coefficients[i] * ((j % 2 == 0) ? 1 : -1) * binom(i, j);
			poly[i - j] += coeff;
		}
		result[j] = poly;
	}
	return result;
}

std::vector<Polynomial> poly_x_over_y(const Polynomial& f)
{
	int n = f.degree; // degree of f
	// We will have (n+1) polynomials: one for each power of y, from y^0 up to y^n.
	// For j = n-i, the corresponding polynomial in x is the monomial f[i]*x^i.
	std::vector<Polynomial> result(n + 1);
	for (int i = 0; i <= n; i++) {
		int j = n - i; // j is the power of y
		// Create a polynomial in x which is a monomial of degree i.
		// Represented as a vector of length (i+1) with the coefficient for x^i equal to f[i].
		std::vector<Rational> poly(i + 1, 0);
		poly[i] = f.coefficients[i];
		// Place this polynomial in the result corresponding to y^j.
		result[j] = poly;
	}
	return result;
}

Matrix createSylvesterMatrix(const std::vector<Polynomial>& f_sub, const Polynomial& g) {
	const int m = f_sub.size() - 1;
	const int n = g.degree;
	const int matrixSize = m + n;
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
	//this->polynomial.normalize();
}

RealAlgebraicNumber::RealAlgebraicNumber(const Polynomial& polynomial, const Interval& interval)
	: polynomial(polynomial), interval(interval)
{
	//normalize();
	while ((this->interval.lower_bound - this->interval.upper_bound).abs() > 0.01)
	{
		refine();
	}
	this->polynomial.normalize();
}

RealAlgebraicNumber::RealAlgebraicNumber(const Polynomial& polynomial, const Rational& lowerBound, const Rational& upperBound)
	: polynomial(polynomial), interval{ lowerBound, upperBound }
{
	//normalize();
	while ((this->interval.lower_bound - this->interval.upper_bound).abs() > 0.01)
	{
		refine();
	}
	this->polynomial.normalize();
}

RealAlgebraicNumber::RealAlgebraicNumber(const std::vector<Rational>& coefficients, const Rational& lowerBound, const Rational& upperBound)
	: polynomial(coefficients), interval{ lowerBound, upperBound }
{
	//normalize();
	while ((this->interval.lower_bound - this->interval.upper_bound).abs() > 0.01)
	{
		refine();
	}
	this->polynomial.normalize();
}

RealAlgebraicNumber RealAlgebraicNumber::fromInteger(const int n)
{
	this->polynomial = Polynomial({ -n, 1 });
	this->interval.lower_bound = n-0.01;
	this->interval.upper_bound = n+0.01;
	this->polynomial.normalize();
	return *this;
}

RealAlgebraicNumber RealAlgebraicNumber::operator+(RealAlgebraicNumber& other)
{
	const std::vector<Polynomial> fXminY = poly_x_minus_y(this->polynomial);
	const Matrix sylvester = createSylvesterMatrix(fXminY, other.polynomial);
	Polynomial f3 = determinant(sylvester);
	f3.normalize();

	std::vector<Polynomial> sturm;
	if (f3.sturm_sequence.empty()) {
		sturm = f3.sturmSequence(f3);
	}
	else
	{
		sturm = f3.sturm_sequence;
	}

	// Compute initial interval bounds
	Rational l3 = interval.lower_bound + other.interval.lower_bound;
	Rational r3 = interval.upper_bound + other.interval.upper_bound;

	// Refine interval until exactly one root is isolated
	while (variationCount(sturm, l3) - variationCount(sturm, r3) > 1) {
		this->refine();
		other.refine();
		l3 = interval.lower_bound + other.interval.lower_bound;
		r3 = interval.upper_bound + other.interval.upper_bound;
	}

	const Interval sumInterval = { l3, r3 };
	return { f3, sumInterval };
}

RealAlgebraicNumber RealAlgebraicNumber::operator-(const RealAlgebraicNumber& other)
{
	RealAlgebraicNumber temp = -other;
	return *this + temp;
}

static Rational minRational(const Rational& r1, const Rational& r2, const Rational& r3, const Rational& r4) {
	Rational min = r1;
	if (r2 < min) min = r2;
	if (r3 < min) min = r3;
	if (r4 < min) min = r4;
	return min;
}

static Rational maxRational(const Rational& r1, const Rational& r2, const Rational& r3, const Rational& r4) {
	Rational max = r1;
	if (r2 > max) max = r2;
	if (r3 > max) max = r3;
	if (r4 > max) max = r4;
	return max;
}

void printSylvester(const Matrix& sylvester)
{
	for (int i = 0; i < sylvester.size(); i++)
	{
		for (int j = 0; j < sylvester[i].size(); j++)
		{
			std::cout << sylvester[i][j].toString() << " ";
		}
		std::cout << std::endl;
	}
}

RealAlgebraicNumber RealAlgebraicNumber::operator*(RealAlgebraicNumber& other)
{
	const std::vector<Polynomial> fXoverY = poly_x_over_y(this->polynomial);
	const Matrix sylvester = createSylvesterMatrix(fXoverY, other.polynomial);
	//printSylvester(sylvester);
	Polynomial f3 = determinant(sylvester);
	/*std::vector<Rational> f3Coeff = determinant2(sylvester);
	Polynomial f3 = {f3Coeff};*/
	f3.normalize();

	Polynomial test = { 0,-20,1 };
	test.normalize();

	/*for (auto coef : f3.coefficients)
	{
		std::cout << coef.toString() << " ";
	}*/

	std::vector<Polynomial> sturm;
	if (f3.sturm_sequence.empty()) {
		sturm = f3.sturmSequence(f3);
	}
	else
	{
		sturm = f3.sturm_sequence;
	}

	/*for (int i = 0; i < sturm.size(); i++)
	{
		std::cout << sturm[i].toString() << std::endl;
	}*/

	Rational l3 = minRational(
		interval.lower_bound * other.interval.lower_bound,
		interval.lower_bound * other.interval.upper_bound,
		interval.upper_bound * other.interval.lower_bound,
		interval.upper_bound * other.interval.upper_bound);
	Rational r3 = maxRational(
		interval.lower_bound * other.interval.lower_bound,
		interval.lower_bound * other.interval.upper_bound,
		interval.upper_bound * other.interval.lower_bound,
		interval.upper_bound * other.interval.upper_bound);

	// Refine interval until exactly one root is isolated
	while (variationCount(sturm, l3) - variationCount(sturm, r3) > 1) {
		//auto f1 = RealAlgebraicNumber(polynomial, interval.lower_bound, r3);
		this->refine();
		other.refine();
		l3 = minRational(
			interval.lower_bound * other.interval.lower_bound,
			interval.lower_bound * other.interval.upper_bound,
			interval.upper_bound * other.interval.lower_bound,
			interval.upper_bound * other.interval.upper_bound);
		r3 = maxRational(
			interval.lower_bound * other.interval.lower_bound,
			interval.lower_bound * other.interval.upper_bound,
			interval.upper_bound * other.interval.lower_bound,
			interval.upper_bound * other.interval.upper_bound);
	}

	return { f3, { l3, r3 } };
}

RealAlgebraicNumber RealAlgebraicNumber::operator/(const RealAlgebraicNumber& other)
{
	RealAlgebraicNumber temp = other.inverse();
	return *this * temp;
}

RealAlgebraicNumber RealAlgebraicNumber::operator-() const
{
	return { polynomial.reflectY(), -interval.upper_bound, -interval.lower_bound };
}

RealAlgebraicNumber RealAlgebraicNumber::inverse() const
{
	std::vector<Rational> inverseCo(polynomial.coefficients.rbegin(), polynomial.coefficients.rend());
	Rational inverseLower = interval.upper_bound.inverse();
	Rational inverseUpper = interval.lower_bound.inverse();
	return { inverseCo, inverseLower, inverseUpper };
}

bool RealAlgebraicNumber::operator==(RealAlgebraicNumber& other)
{
	// make both polynomials minimal
	// if minimal polynomials are not equal, return false
	// if minimal polynomials are equal, check if intervals overlap
	//
	if (this->polynomial != other.polynomial)
	{
		return false;
	}
	
	for (int i = 0; i < 50; i ++)
	{
		if (this->interval.lower_bound == other.interval.lower_bound && this->interval.upper_bound == other.interval.upper_bound)
		{
			std::cout << "Equal after " << i << std::endl;
			return true;
		}

		if (this->interval.lower_bound > other.interval.upper_bound || this->interval.upper_bound < other.interval.lower_bound)
		{
			return false;
		}
		this->refine();
		other.refine();
	}
	std::cout << "Equal after 50" << std::endl;
	return true;
}

bool RealAlgebraicNumber::operator!=(RealAlgebraicNumber& other)
{
	return !(*this == other);
}

bool RealAlgebraicNumber::operator<(RealAlgebraicNumber& other)
{
	if (*this == other)
	{
		return false;
	}

	while (true)
	{
		if (this->interval.upper_bound < other.interval.lower_bound)
		{
			return true;
		}
		if (this->interval.lower_bound > other.interval.upper_bound)
		{
			return false;
		}
		this->refine();
		other.refine();
	}
}

bool RealAlgebraicNumber::operator>(RealAlgebraicNumber& other)
{
	if (*this == other)
	{
		return false;
	}

	while (true)
	{
		if (this->interval.lower_bound > other.interval.upper_bound)
		{
			return true;
		}
		if (this->interval.upper_bound < other.interval.lower_bound)
		{
			return false;
		}
		this->refine();
		other.refine();
	}
}

bool RealAlgebraicNumber::operator<=(RealAlgebraicNumber& other)
{
	return !(*this > other);
}

bool RealAlgebraicNumber::operator>=(RealAlgebraicNumber& other)
{
	return !(*this < other);
}

RealAlgebraicNumber RealAlgebraicNumber::sqrt(const int n)
{
	// 1. Check Non-Negativity (for even n): If n is even, ensure the algebraic number is non-negative.
	// If it is negative, the n-th root is not real and the operation is undefined.
	// 2. Construct the New Polynomial: For the given polynomial p(x) and the root a,
	// construct the polynomial q(x) = p(x^2). This polynomial will have a^(1/n) as one of its roots.
	// 3. Adjust the Interval : Compute the new interval for a^(1/n) by taking the n-th roots of the endpoints of the original interval:
	// If n is even or n is odd and a is non-negative, the new interval is (sqrt(l), sqrt(r)).
	// If n is odd and a is negative, the new interval is (-sqrt(r.abs()), -sqrt(l.abs())).
	// 4. Refine the Interval: Refine the interval until exactly one root is isolated.
	// 5. Return the Result: Return the algebraic number with the new polynomial and interval.
	if (n % 2 == 0 && interval.lower_bound < 0)
	{
		throw std::invalid_argument("Cannot compute even root of negative number");
	}

	const int newDegree = polynomial.degree * n;
	std::vector<Rational> newCoeffs(newDegree + 1, 0);

	for (int i = 0; i < polynomial.coefficients.size(); ++i) {
		newCoeffs[i * n] = polynomial.coefficients[i];
	}

	Polynomial f3 = { newCoeffs };
	//newPoly.print();

	/*std::vector<Polynomial> sturm;
	if (f3.sturm_sequence.empty()) {
		sturm = f3.sturmSequence(f3);
	}
	else
	{
		sturm = f3.sturm_sequence;
	}*/

	Rational l3;
	Rational r3;

	if (n % 2 != 0 && interval.lower_bound < 0)
	{
		l3 = -interval.upper_bound.abs().sqrt(n);
		r3 = -interval.lower_bound.abs().sqrt(n);
	}
	else
	{
		l3 = interval.lower_bound.sqrt(n);
		r3 = interval.upper_bound.sqrt(n);
	}

	//while (variationCount(sturm, l3) - variationCount(sturm, r3) > 1) {
	//	//auto f1 = RealAlgebraicNumber(polynomial, interval.lower_bound, r3);
	//	this->refine();
	//	if (n % 2 != 0 && interval.lower_bound < 0)
	//	{
	//		l3 = -interval.upper_bound.abs().sqrt(n);
	//		r3 = -interval.lower_bound.abs().sqrt(n);
	//	}
	//	else
	//	{
	//		l3 = interval.lower_bound.sqrt(n);
	//		r3 = interval.upper_bound.sqrt(n);
	//	}
	//}

	return { f3, { l3, r3 } };
}

RealAlgebraicNumber RealAlgebraicNumber::pow(int n) const
{
	if (n == 0) {
		std::vector<Rational> coefficients = { 1 };
		return { coefficients, 1, 1 };
	}
	//if (n == 1) return *this;
	RealAlgebraicNumber result = *this;
	for (int i = 1; i < n; i++) {
		RealAlgebraicNumber temp = *this;
		result *= temp;
	}
	return result;
}

int RealAlgebraicNumber::countSignVariations(const std::vector<Rational>& sequence) const {
	int variations = 0;
	int prevSign = 0;

	for (const Rational& value : sequence) {
		if (value != 0) {
			const int currentSign = (value > 0) ? 1 : -1;
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
	for (const Rational& coeff : polynomial.coefficients) {
		if (fInf < coeff.abs())
			fInf = coeff.abs();
		//fInf = std::max(fInf, coeff.abs());
	}

	const Rational p = Rational(1) / (Rational(1) + fInf);

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
	output += polynomial.toString() + " @ ";
	const double lower = static_cast<double>(interval.lower_bound);
	output += std::to_string(lower);
	output += " <= x <= ";
	const double upper = static_cast<double>(interval.upper_bound);
	output += std::to_string(upper);

	return output;
}