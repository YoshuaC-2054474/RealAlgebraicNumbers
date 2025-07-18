#include "RealAlgebraicNumber.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <utility>

#include "MyTimer.h"
#include "MyMatrix.h"

// Hash specialization for std::pair<int, int>
namespace std {
	template<>
	struct hash<std::pair<int, int>> {
		std::size_t operator()(const std::pair<int, int>& p) const {
			return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
		}
	};
}

// Optimized binomial coefficient with memoization
static std::unordered_map<std::pair<int, int>, long long> binomialCache;

long long binomialCoeff(const int n, int k) {
	//PROFILE_FUNCTION
	if (k < 0 || k > n) {
		return 0;
	}
	if (k == 0 || k == n) {
		return 1;
	}
	if (k > n / 2) {
		// Use symmetry (nCk = nC(n-k)) to reduce calculations
		k = n - k;
	}

	// Check cache first
	auto key = std::make_pair(n, k);
	auto it = binomialCache.find(key);
	if (it != binomialCache.end()) {
		return it->second;
	}

	long long res = 1;
	for (int i = 1; i <= k; ++i) {
		res = res * (n - i + 1) / i;
	}

	// Store in cache
	binomialCache[key] = res;
	return res;
}


MyMatrix<Polynomial> constructSylvesterMatrixForSum(const Polynomial& p, const Polynomial& q) {
	PROFILE_FUNCTION
		if (p.isZero() || q.isZero()) {
			return MyMatrix<Polynomial>(0, 0);
		}

	const int m = p.degree;
	const int n = q.degree;
	const int matSize = m + n;

	MyMatrix<Polynomial> sylvester(matSize, matSize);

	// Fill n rows for p(x) - optimized with fewer temporary objects
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j <= m; ++j) {
			sylvester.set_element(i, i + j, Polynomial(p.coeff(j)));
		}
	}

	// Pre-allocate coefficient vector and polynomial for G(x)
	std::vector<Polynomial> coeffs_of_G_in_x(n + 1);
	static const Polynomial zPoly({ Rational(0), Rational(1) });

	for (int s = 0; s <= n; ++s) {
		Polynomial g_s_z_poly;

		for (int j = s; j <= n; ++j) {
			const Rational q_j = q.coeff(j);
			if (q_j == 0) continue;

			const long long termBinomial = binomialCoeff(j, s);
			const Rational term_minus_one_power = ((s & 1) == 0) ? Rational(1) : Rational(-1);
			const Rational scalarCoeff = q_j * Rational(termBinomial) * term_minus_one_power;

			// Create z^(j-s) more efficiently
			const int power = j - s;
			if (power == 0) {
				g_s_z_poly += Polynomial(scalarCoeff);
			}
			else {
				std::vector<Rational> zCoeffs(power + 1, Rational(0));
				zCoeffs[power] = scalarCoeff;
				g_s_z_poly += Polynomial(std::move(zCoeffs));
			}
		}
		coeffs_of_G_in_x[s] = std::move(g_s_z_poly);
	}

	// Fill the m rows for G(x)
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j <= n; ++j) {
			sylvester.set_element(n + i, i + j, coeffs_of_G_in_x[j]);
		}
	}
	return sylvester;
}


MyMatrix<Polynomial> constructSylvesterMatrixForProduct(const Polynomial& p, const Polynomial& q) {
	//PROFILE_FUNCTION
	const int m = p.degree;
	const int n = q.degree;

	if (m == -1 || n == -1) {
		return MyMatrix<Polynomial>(0, 0);
	}

	const int matSize = m + n;
	MyMatrix<Polynomial> sylvester(matSize, matSize);

	// Fill rows for p(x) - optimized
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j <= m; ++j) {
			sylvester.set_element(i, i + j, Polynomial(p.coeff(j)));
		}
	}

	// Fill rows for Q_zx(x) - optimized with pre-allocated vectors
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j <= n; ++j) {
			const Rational q_j = q.coeff(j);

			if (q_j == 0) {
				sylvester.set_element(n + i, i + (n - j), Polynomial(Rational(0)));
			}
			else {
				std::vector<Rational> z_poly_coeffs(j + 1, Rational(0));
				z_poly_coeffs[j] = q_j;
				sylvester.set_element(n + i, i + (n - j), Polynomial(std::move(z_poly_coeffs)));
			}
		}
	}
	return sylvester;
}


MyMatrix<Polynomial> constructSylvesterMatrixForPower(const Polynomial& p, const int k) {
	//PROFILE_FUNCTION
	if (k <= 0) {
		throw std::invalid_argument("Exponent k must be positive for power resultant construction.");
	}
	if (p.isZero()) {
		return MyMatrix<Polynomial>(0, 0);
	}

	const int m = p.degree;
	const int n = k;
	const int matSize = m + n;
	MyMatrix<Polynomial> sylvester(matSize, matSize);

	// Fill n rows for p(x)
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j <= m; ++j) {
			sylvester.set_element(i, i + j, Polynomial(p.coeff(j)));
		}
	}

	// Pre-create polynomials for G(x) coefficients
	static const Polynomial zPoly({ Rational(0), Rational(1) });
	static const Polynomial negOne(Rational(-1));
	static const Polynomial zero(Rational(0));

	// Fill m rows for G(x) = -x^k + z
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j <= n; ++j) {
			if (j == k) {
				sylvester.set_element(n + i, i + j, negOne);
			}
			else if (j == 0) {
				sylvester.set_element(n + i, i + j, zPoly);
			}
			else {
				sylvester.set_element(n + i, i + j, zero);
			}
		}
	}
	return sylvester;
}

RealAlgebraicNumber::RealAlgebraicNumber()
	: polynomial({ 0 }), interval({ .lowerBound = 0.0, .upperBound = 0.0 }) {
}

RealAlgebraicNumber::RealAlgebraicNumber(const Polynomial& polynomial, const Interval& interval)
	: polynomial(polynomial), interval(interval) {
	refineToTolerance();
	this->polynomial.normalize(this->interval.lowerBound, this->interval.upperBound);
}

RealAlgebraicNumber::RealAlgebraicNumber(const Polynomial& polynomial, const Rational& lowerBound,
	const Rational& upperBound)
	: polynomial(polynomial), interval{ .lowerBound = lowerBound, .upperBound = upperBound } {
	refineToTolerance();
	this->polynomial.normalize(this->interval.lowerBound, this->interval.upperBound);
}

RealAlgebraicNumber::RealAlgebraicNumber(const std::vector<Rational>& coefficients, const Rational& lowerBound,
	const Rational& upperBound)
	: polynomial(coefficients), interval{ .lowerBound = lowerBound, .upperBound = upperBound } {
	refineToTolerance();
	this->polynomial.normalize(this->interval.lowerBound, this->interval.upperBound);
}

RealAlgebraicNumber::RealAlgebraicNumber(const int value)
	: polynomial({ -value, 1 }), interval{ .lowerBound = Rational(value), .upperBound = Rational(value) } {
	this->polynomial.normalize(this->interval.lowerBound, this->interval.upperBound);
}

RealAlgebraicNumber::RealAlgebraicNumber(const double value) {
	const Rational rationalValue(value);
	this->polynomial = Polynomial({ -rationalValue, 1 });
	this->interval.lowerBound = rationalValue;
	this->interval.upperBound = rationalValue;
	this->polynomial.normalize(this->interval.lowerBound, this->interval.upperBound);
}

void RealAlgebraicNumber::refineToTolerance() {
	static const Rational tolerance(1, 100); // 0.01
	while ((this->interval.lowerBound - this->interval.upperBound).abs() > tolerance) {
		refine();
	}
}

RealAlgebraicNumber operator+(const RealAlgebraicNumber& lhs, const RealAlgebraicNumber& rhs) {
	PROFILE_FUNCTION
	RealAlgebraicNumber lhsCopy = lhs;
	lhsCopy += rhs;
	return lhsCopy;
}

RealAlgebraicNumber operator-(const RealAlgebraicNumber& lhs, const RealAlgebraicNumber& rhs) {
	PROFILE_FUNCTION
	RealAlgebraicNumber lhsCopy = lhs;
	lhsCopy -= rhs;
	return lhsCopy;
}

RealAlgebraicNumber RealAlgebraicNumber::operator-() const {
	PROFILE_FUNCTION
	return { polynomial.reflectY(), -interval.upperBound, -interval.lowerBound };
}

RealAlgebraicNumber operator*(const RealAlgebraicNumber& lhs, const RealAlgebraicNumber& rhs) {
	PROFILE_FUNCTION
	RealAlgebraicNumber lhsCopy = lhs;
	lhsCopy *= rhs;
	return lhsCopy;
}

RealAlgebraicNumber operator/(const RealAlgebraicNumber& lhs, const RealAlgebraicNumber& rhs) {
	PROFILE_FUNCTION
	RealAlgebraicNumber lhsCopy = lhs;
	lhsCopy /= rhs;
	return lhsCopy;
}

RealAlgebraicNumber RealAlgebraicNumber::operator+=(const RealAlgebraicNumber& other) {
	PROFILE_FUNCTION

	// Construct Sylvester matrix and compute determinant
	MyMatrix<Polynomial> sylvester_mat = constructSylvesterMatrixForSum(this->polynomial, other.polynomial);
	Polynomial f3 = sylvester_mat.determinant();

	// Get or compute Sturm sequence
	std::vector<Polynomial> sturm;
	if (f3.sturm_sequence.empty()) {
		sturm = f3.sturmSequence(f3);
	}
	else {
		sturm = f3.sturm_sequence;
	}

	// Create working copies only if needed
	RealAlgebraicNumber thisCopy = *this;
	RealAlgebraicNumber otherCopy = other;

	// Compute initial interval bounds
	Rational l3 = thisCopy.interval.lowerBound + otherCopy.interval.lowerBound;
	Rational r3 = thisCopy.interval.upperBound + otherCopy.interval.upperBound;

	// Refine interval until exactly one root is isolated
	while (variationCount(sturm, l3) - variationCount(sturm, r3) > 1) {
		thisCopy.refine();
		otherCopy.refine();
		l3 = thisCopy.interval.lowerBound + otherCopy.interval.lowerBound;
		r3 = thisCopy.interval.upperBound + otherCopy.interval.upperBound;
	}

	f3.normalize(l3, r3);

	// Update this object
	this->polynomial = std::move(f3);
	this->interval = { .lowerBound = l3, .upperBound = r3 };

	return *this;
}

RealAlgebraicNumber RealAlgebraicNumber::operator-=(const RealAlgebraicNumber& other) {
	PROFILE_FUNCTION
	return *this += -other;
}

Rational minRational(const Rational& r1, const Rational& r2, const Rational& r3, const Rational& r4) {
	//PROFILE_FUNCTION
	Rational min = r1;
	if (r2 < min) min = r2;
	if (r3 < min) min = r3;
	if (r4 < min) min = r4;
	return min;
}

Rational maxRational(const Rational& r1, const Rational& r2, const Rational& r3, const Rational& r4) {
	//PROFILE_FUNCTION
	Rational max = r1;
	if (r2 > max) max = r2;
	if (r3 > max) max = r3;
	if (r4 > max) max = r4;
	return max;
}

RealAlgebraicNumber RealAlgebraicNumber::operator*=(const RealAlgebraicNumber& other) {
	PROFILE_FUNCTION

	// Create working copies and normalize
	RealAlgebraicNumber thisCopy = *this;
	RealAlgebraicNumber otherCopy = other;

	otherCopy.polynomial.normalize(otherCopy.interval.lowerBound, otherCopy.interval.upperBound);
	thisCopy.polynomial.normalize(thisCopy.interval.lowerBound, thisCopy.interval.upperBound);

	// Construct Sylvester matrix and compute determinant
	MyMatrix<Polynomial> sylvester_prod_mat = constructSylvesterMatrixForProduct(
		thisCopy.polynomial, otherCopy.polynomial);
	Polynomial f3 = sylvester_prod_mat.determinant();

	// Compute initial interval bounds using optimized min/max
	const auto& tLower = thisCopy.interval.lowerBound;
	const auto& tUpper = thisCopy.interval.upperBound;
	const auto& oLower = otherCopy.interval.lowerBound;
	const auto& oUpper = otherCopy.interval.upperBound;

	Rational l3 = minRational(tLower * oLower, tLower * oUpper, tUpper * oLower, tUpper * oUpper);
	Rational r3 = maxRational(tLower * oLower, tLower * oUpper, tUpper * oLower, tUpper * oUpper);

	f3.normalize(l3, r3);

	// Get or compute Sturm sequence
	std::vector<Polynomial> sturm;
	if (f3.sturm_sequence.empty()) {
		sturm = f3.sturmSequence(f3);
	}
	else {
		sturm = f3.sturm_sequence;
	}

	// Refine interval until exactly one root is isolated
	while (variationCount(sturm, l3) - variationCount(sturm, r3) > 1) {
		thisCopy.refine();
		otherCopy.refine();

		const auto& tLower = thisCopy.interval.lowerBound;
		const auto& tUpper = thisCopy.interval.upperBound;
		const auto& oLower = otherCopy.interval.lowerBound;
		const auto& oUpper = otherCopy.interval.upperBound;

		l3 = minRational(tLower * oLower, tLower * oUpper, tUpper * oLower, tUpper * oUpper);
		r3 = maxRational(tLower * oLower, tLower * oUpper, tUpper * oLower, tUpper * oUpper);
	}

	f3.normalize(l3, r3);

	// Update this object
	this->polynomial = std::move(f3);
	this->interval = { .lowerBound = l3, .upperBound = r3 };

	return *this;
}

RealAlgebraicNumber RealAlgebraicNumber::operator/=(const RealAlgebraicNumber& other) {
	PROFILE_FUNCTION
	return *this *= other.inverse();
}

bool RealAlgebraicNumber::operator==(const RealAlgebraicNumber& other) const {
	PROFILE_FUNCTION
		// make both polynomials minimal
		// if minimal polynomials are not equal, return false
		// if minimal polynomials are equal, check if intervals overlap

	RealAlgebraicNumber otherCopy = other;
	RealAlgebraicNumber thisCopy = *this;
	otherCopy.polynomial.normalize(otherCopy.interval.lowerBound, otherCopy.interval.upperBound);
	thisCopy.polynomial.normalize(thisCopy.interval.lowerBound, thisCopy.interval.upperBound);

	if (thisCopy.polynomial != otherCopy.polynomial) {
		return false;
	}

	//for (int i = 0; i < 50; i++) {
	//	// intervals are the same
	//	if (thisCopy.interval.lowerBound == otherCopy.interval.lowerBound && thisCopy.interval.upperBound ==
	//		otherCopy.interval.upperBound) {
	//		return true;
	//	}
	//	// one interval is contained in the other
	//	if (thisCopy.interval.lowerBound > otherCopy.interval.lowerBound && thisCopy.interval.upperBound < otherCopy.
	//		interval.upperBound) {
	//		return true;
	//	}
	//	if (thisCopy.interval.lowerBound < otherCopy.interval.lowerBound && thisCopy.interval.upperBound > otherCopy.
	//		interval.upperBound) {
	//		return true;
	//	}

	//	// intervals are separated
	//	if (thisCopy.interval.lowerBound > otherCopy.interval.upperBound || thisCopy.interval.upperBound < otherCopy.
	//		interval.lowerBound) {
	//		return false;
	//	}
	//	thisCopy.refine();
	//	otherCopy.refine();
	//}


	//std::cout << "Equal after 50" << '\n';
	/*return true;*/

	/*if (this->polynomial != other.polynomial) {
		return false;
	}

	return (*this - other).isZero();*/

	// intervals are the same
	if (thisCopy.interval.lowerBound == otherCopy.interval.lowerBound && thisCopy.interval.upperBound ==
		otherCopy.interval.upperBound) {
		return true;
	}
	// one interval is contained in the other
	if (thisCopy.interval.lowerBound > otherCopy.interval.lowerBound && thisCopy.interval.upperBound < otherCopy.
		interval.upperBound) {
		return true;
	}
	if (thisCopy.interval.lowerBound < otherCopy.interval.lowerBound && thisCopy.interval.upperBound > otherCopy.
		interval.upperBound) {
		return true;
	}

	// intervals are separated
	if (thisCopy.interval.lowerBound > otherCopy.interval.upperBound || thisCopy.interval.upperBound < otherCopy.
		interval.lowerBound) {
		return false;
	}

	// intervals overlap, use order representation
	int thisOrder = thisCopy.intervalToOrder();
	int otherOrder = otherCopy.intervalToOrder();

	return thisOrder == otherOrder;
}

bool RealAlgebraicNumber::operator!=(const RealAlgebraicNumber& other) const {
	PROFILE_FUNCTION
		return !(*this == other);
}

bool RealAlgebraicNumber::operator<(const RealAlgebraicNumber& other) const {
	//PROFILE_FUNCTION
	RealAlgebraicNumber otherCopy = other;
	RealAlgebraicNumber thisCopy = *this;

	if (thisCopy == otherCopy) {
		return false;
	}

	while (true) {
		if (thisCopy.interval.upperBound < otherCopy.interval.lowerBound) {
			return true;
		}
		if (thisCopy.interval.lowerBound > otherCopy.interval.upperBound) {
			return false;
		}
		thisCopy.refine();
		otherCopy.refine();
	}
}

bool RealAlgebraicNumber::operator>(const RealAlgebraicNumber& other) const {
	//PROFILE_FUNCTION
	return other < *this;
}

bool RealAlgebraicNumber::operator<=(const RealAlgebraicNumber& other) const {
	//PROFILE_FUNCTION
	return !(*this > other);
}

bool RealAlgebraicNumber::operator>=(const RealAlgebraicNumber& other) const {
	//PROFILE_FUNCTION
	return !(*this < other);
}


RealAlgebraicNumber RealAlgebraicNumber::inverse() const {
	PROFILE_FUNCTION
		std::vector<Rational> inverseCo;
	inverseCo.reserve(polynomial.coefficients.size());

	// Use reverse iterators for better performance
	std::reverse_copy(polynomial.coefficients.begin(), polynomial.coefficients.end(),
		std::back_inserter(inverseCo));

	return { std::move(inverseCo), interval.upperBound.inverse(), interval.lowerBound.inverse() };
}

RealAlgebraicNumber RealAlgebraicNumber::sqrt(const int n) const {
	PROFILE_FUNCTION
	// 1. Check Non-Negativity (for even n): If n is even, ensure the algebraic number is non-negative.
	// If it is negative, the n-th root is not real and the operation is undefined.
	// 2. Construct the New Polynomial: For the given polynomial p(x) and the root a,
	// construct the polynomial q(x) = p(x^2). This polynomial will have a^(1/n) as one of its roots.
	// 3. Adjust the Interval : Compute the new interval for a^(1/n) by taking the n-th roots of the endpoints of the original interval:
	// If n is even or n is odd and a is non-negative, the new interval is (sqrt(l), sqrt(r)).
	// If n is odd and a is negative, the new interval is (-sqrt(r.abs()), -sqrt(l.abs())).
	// 4. Refine the Interval: Refine the interval until exactly one root is isolated.
	// 5. Return the Result: Return the algebraic number with the new polynomial and interval.
	if (n % 2 == 0 && interval.lowerBound < 0) {
		throw std::invalid_argument("Cannot compute even root of negative number");
	}

	const int newDegree = polynomial.degree * n;
	std::vector<Rational> newCoeffs(newDegree + 1, 0);

	for (size_t i = 0; i < polynomial.coefficients.size(); ++i) {
		newCoeffs[i * n] = polynomial.coefficients[i];
	}

	Polynomial f3 = {newCoeffs};
	//newPoly.print();

	std::vector<Polynomial> sturm;
	if (f3.sturm_sequence.empty()) {
		sturm = f3.sturmSequence(f3);
	}
	else
	{
		sturm = f3.sturm_sequence;
	}

	RealAlgebraicNumber thisCopy = *this;

	Rational l3;
	Rational r3;

	if (n % 2 != 0 && interval.lowerBound < 0) {
		l3 = -thisCopy.interval.upperBound.abs().sqrt(n) - 0.01;
		r3 = -thisCopy.interval.lowerBound.abs().sqrt(n) + 0.01;
	}
	else {
		l3 = thisCopy.interval.lowerBound.sqrt(n) - 0.01;
		r3 = thisCopy.interval.upperBound.sqrt(n) + 0.01;
	}

	while (variationCount(sturm, l3) - variationCount(sturm, r3) > 1) {
		//auto f1 = RealAlgebraicNumber(polynomial, interval.lowerBound, r3);
		thisCopy.refine();
		if (n % 2 != 0 && thisCopy.interval.lowerBound < 0)
		{
			l3 = -thisCopy.interval.upperBound.abs().sqrt(n);
			r3 = -thisCopy.interval.lowerBound.abs().sqrt(n);
		}
		else
		{
			l3 = thisCopy.interval.lowerBound.sqrt(n);
			r3 = thisCopy.interval.upperBound.sqrt(n);
		}
	}

	f3.normalize(l3, r3);
	return {f3, {.lowerBound = l3, .upperBound = r3}};
}

RealAlgebraicNumber RealAlgebraicNumber::pow(int n) const {
	PROFILE_FUNCTION

	if (n == 0) {
		std::vector<Rational> coefficients = {-1, 1};
		return {coefficients, 0.99, 1.01}; // Represents the number 1
	}

	RealAlgebraicNumber thisCopy = *this;
	thisCopy.polynomial.normalize(thisCopy.interval.lowerBound, thisCopy.interval.upperBound);

	MyMatrix<Polynomial> sylvesterMat = constructSylvesterMatrixForPower(thisCopy.polynomial, n);
	Polynomial f3 = sylvesterMat.determinant();

	std::vector<Polynomial> sturm;
	if (f3.sturm_sequence.empty()) {
		sturm = f3.sturmSequence(f3);
	}
	else
	{
		sturm = f3.sturm_sequence;
	}

	Rational l3;
	Rational r3;

	if (thisCopy.interval.lowerBound >= 0) {
		// Both endpoints are non-negative
		l3 = thisCopy.interval.lowerBound.pow(n);
		r3 = thisCopy.interval.upperBound.pow(n);
	}
	else if (thisCopy.interval.upperBound <= 0) {
		// Both endpoints are non-positive
		if (n % 2 == 0) {
			// Even power: result is positive, order flips
			l3 = thisCopy.interval.upperBound.pow(n);
			r3 = thisCopy.interval.lowerBound.pow(n);
		}
		else {
			// Odd power: result stays negative, order same
			l3 = thisCopy.interval.lowerBound.pow(n);
			r3 = thisCopy.interval.upperBound.pow(n);
		}
	}

	while (variationCount(sturm, l3) - variationCount(sturm, r3) > 1) {
		//auto f1 = RealAlgebraicNumber(polynomial, interval.lowerBound, r3);
		thisCopy.refine();
		if (thisCopy.interval.lowerBound >= 0) {
			// Both endpoints are non-negative
			l3 = thisCopy.interval.lowerBound.pow(n);
			r3 = thisCopy.interval.upperBound.pow(n);
		}
		else if (thisCopy.interval.upperBound <= 0) {
			// Both endpoints are non-positive
			if (n % 2 == 0) {
				// Even power: result is positive, order flips
				l3 = thisCopy.interval.upperBound.pow(n);
				r3 = thisCopy.interval.lowerBound.pow(n);
			}
			else {
				// Odd power: result stays negative, order same
				l3 = thisCopy.interval.lowerBound.pow(n);
				r3 = thisCopy.interval.upperBound.pow(n);
			}
		}
	}

	f3.normalize(l3, r3);
	return {f3, {.lowerBound = l3, .upperBound = r3}};
}

bool RealAlgebraicNumber::isZero() const {
	if (this->polynomial == Polynomial({0, 1})) {
		return this->interval.lowerBound <= 0 && this->interval.upperBound >= 0;
	}
	return false;
}

void RealAlgebraicNumber::testOperators() {
	//PROFILE_FUNCTION
	RealAlgebraicNumber a = 2;
	RealAlgebraicNumber b = 3;
	auto c = RealAlgebraicNumber({-2, 0, 1}, {.lowerBound = {1, 1}, .upperBound = {2, 1}}); // sqrt(2)
	auto d = RealAlgebraicNumber({-3, 0, 1}, {.lowerBound = {1, 1}, .upperBound = {2, 1}}); // sqrt(3)

	RealAlgebraicNumber e = 5;
	RealAlgebraicNumber f = -1;
	auto g = RealAlgebraicNumber({1, 0, -10, 0, 1}, {.lowerBound = {3146, 1000}, .upperBound = {3147, 1000}});
	// x^4 - 10x^2 + 1
	auto h = RealAlgebraicNumber({1, 0, -10, 0, 1}, {.lowerBound = {-318, 1000}, .upperBound = {-317, 1000}});
	RealAlgebraicNumber i = 6;
	auto j = RealAlgebraicNumber({2, -3}, {.lowerBound = {666, 1000}, .upperBound = {667, 1000}});
	auto k = RealAlgebraicNumber({4, 0, -12, 0, 9}, {.lowerBound = {816, 1000}, .upperBound = {817, 1000}});
	// 9x^4 - 12x^2 + 4
	auto l = RealAlgebraicNumber({36, 0, -12, 0, 1}, {.lowerBound = {244, 100}, .upperBound = {245, 100}});
	// x^4 - 12x^2 + 36

	if (a + b != e || c + d != g)
		std::cout << "operator+ not working!\n";
	if (a - b != f || c - d != h)
		std::cout << "operator- not working!\n";
	if (a + b - b != a || c + d - c != d)
		std::cout << "operator+ and operator- not working together!\n";
	if (a * b != i || c * d != l)
		std::cout << "operator* not working!\n";
	if (a / b != j || c / d != k)
		std::cout << "operator/ not working!\n";
	if ((a * b) / b != a || (c * d) / c != d)
		std::cout << "operator* and operator/ not working together!\n";

	auto aa = a;
	auto cc = c;
	if (aa += b, cc += d; aa != e || cc != g)
		std::cout << "operator+= not working!\n";
	aa = a;
	cc = c;
	if (aa -= b, cc -= d; aa != f || cc != h)
		std::cout << "operator-= not working!\n";
	aa = a;
	cc = c;
	if (aa *= b, cc *= d; aa != i || cc != l)
		std::cout << "operator*= not working!\n";
	aa = a;
	cc = c;
	if (aa /= b, cc /= d; aa != j || cc != k)
		std::cout << "operator/= not working!\n";

	auto aAlt = RealAlgebraicNumber({-2, 1}, {.lowerBound = {1, 1}, .upperBound = {4, 1}});
	auto cAlt = RealAlgebraicNumber({-2, 0, 1}, {.lowerBound = {1414, 1000}, .upperBound = {1415, 1000}});

	if ((a == aAlt) == false || (c == cAlt) == false)
		std::cout << "operator== not working!\n";
	if ((a != b) == false || (c != d) == false)
		std::cout << "operator!= not working!\n";
	if ((a < b) == false || (c < d) == false)
		std::cout << "operator< not working!\n";
	if ((b > a) == false || (d > c) == false)
		std::cout << "operator> not working!\n";
	if ((a <= b) == false || (c <= d) == false || (a <= a) == false || (c <= c) == false)
		std::cout << "operator<= not working!\n";
	if ((b >= a) == false || (d >= c) == false || (b >= b) == false || (d >= d) == false)
		std::cout << "operator>= not working!\n";

	auto invA = RealAlgebraicNumber({-1, 2}, {.lowerBound = {3, 10}, .upperBound = {6, 10}});
	auto invC = RealAlgebraicNumber({-1, 0, 2}, {.lowerBound = {707, 1000}, .upperBound = {708, 1000}});

	if (a.inverse() != invA || c.inverse() != invC)
		std::cout << "inverse() not working!\n";
	if (a.inverse().inverse() != a || c.inverse().inverse() != c)
		std::cout << "double inverse() not working!\n";
	if (a.sqrt() != c || b.sqrt() != d)
		std::cout << "sqrt() not working!\n";
	if (c.pow(2) != a || d.pow(2) != b)
		std::cout << "pow() not working!\n";
	if (b.pow(2).sqrt() != b || d.sqrt().pow(2) != d)
		std::cout << "sqrt() and pow() not working together! (1)\n";
	if (a.sqrt().pow(2) != a || c.pow(2).sqrt() != c)
		std::cout << "sqrt() and pow() not working together! (2)\n";

	std::cout << "\n";
}

int RealAlgebraicNumber::countSignVariations(const std::vector<Rational>& sequence) {
	//PROFILE_FUNCTION
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

Rational RealAlgebraicNumber::evaluatePoly(const std::vector<Rational>& sequence, const Rational& x) {
	PROFILE_FUNCTION
		Rational result = 0;
	const int seqSize = static_cast<int>(sequence.size());
	for (int i = seqSize - 1; i >= 0; i--) {
		result = result * x + sequence[i];
	}
	return result;
}

int RealAlgebraicNumber::variationCount(const std::vector<Polynomial>& sturm, const Rational& x) {
	PROFILE_FUNCTION
	std::vector<Rational> evaluations;
	evaluations.reserve(sturm.size());

	for (const auto& poly : sturm) {
		evaluations.push_back(evaluatePoly(poly.coefficients, x));
	}

	return countSignVariations(evaluations);
}

int RealAlgebraicNumber::intervalToOrder() {
	PROFILE_FUNCTION

	// Find maximum absolute coefficient more efficiently
	Rational fInf = polynomial.coefficients.empty() ? Rational(0) : polynomial.coefficients[0].abs();
	for (size_t i = 1; i < polynomial.coefficients.size(); ++i) {
		const Rational abs_coeff = polynomial.coefficients[i].abs();
		if (abs_coeff > fInf) {
			fInf = abs_coeff;
		}
	}

	// Get or compute Sturm sequence
	std::vector<Polynomial> sturm;
	if (polynomial.sturm_sequence.empty()) {
		sturm = polynomial.sturmSequence(polynomial);
	}
	else {
		sturm = polynomial.sturm_sequence;
	}

	return variationCount(sturm, -1 - fInf) - variationCount(sturm, interval.upperBound);
}

void RealAlgebraicNumber::normalize() {
	PROFILE_FUNCTION

	// Find maximum absolute coefficient more efficiently
	Rational fInf = polynomial.coefficients.empty() ? Rational(0) : polynomial.coefficients[0].abs();
	for (size_t i = 1; i < polynomial.coefficients.size(); ++i) {
		const Rational abs_coeff = polynomial.coefficients[i].abs();
		if (abs_coeff > fInf) {
			fInf = abs_coeff;
		}
	}

	const Rational p = Rational(1) / (Rational(1) + fInf);

	// Get or compute Sturm sequence
	std::vector<Polynomial> sturm;
	if (polynomial.sturm_sequence.empty()) {
		sturm = polynomial.sturmSequence(polynomial);
	}
	else {
		sturm = polynomial.sturm_sequence;
	}

	const int varL = variationCount(sturm, interval.lowerBound);
	const int varNegP = variationCount(sturm, -p);
	const int varP = variationCount(sturm, p);

	if (varL > varNegP) {
		interval.upperBound = -p;
	}
	else if (varNegP > varP) {
		interval.lowerBound = Rational(0);
		interval.upperBound = Rational(0);
	}
	else {
		interval.lowerBound = p;
	}
}


void RealAlgebraicNumber::refine() {
	PROFILE_FUNCTION

		// Get or compute Sturm sequence
	std::vector<Polynomial> sturm;
	if (polynomial.sturm_sequence.empty()) {
		sturm = polynomial.sturmSequence(polynomial);
	}
	else {
		sturm = polynomial.sturm_sequence;
	}

	const Rational m = (interval.lowerBound + interval.upperBound) / 2;

	const int varL = variationCount(sturm, interval.lowerBound);
	const int varM = variationCount(sturm, m);

	if (varL > varM) {
		interval.upperBound = m;
	}
	else {
		interval.lowerBound = m;
	}
}

std::string RealAlgebraicNumber::toString() const {
	PROFILE_FUNCTION
		std::string output;
	output += polynomial.toString() + " @ ";
	const double lower = static_cast<double>(interval.lowerBound);
	output += std::to_string(lower);
	output += " <= x <= ";
	const double upper = static_cast<double>(interval.upperBound);
	output += std::to_string(upper);

	return output;
}

std::string RealAlgebraicNumber::toDecimalString(int precision) const {
	PROFILE_FUNCTION

	// Early exit for linear polynomials
	if (this->polynomial.degree == 1 && this->polynomial.coefficients[1] == 1) {
		return (this->polynomial.coefficients[0] * -1).toDecimalString(precision);
	}

	RealAlgebraicNumber temp = *this;

	while (true) {
		const std::string lowerString = temp.interval.lowerBound.toDecimalString(precision);
		const std::string upperString = temp.interval.upperBound.toDecimalString(precision);

		// Find decimal point positions
		const size_t lowerDotPos = lowerString.find('.');
		const size_t upperDotPos = upperString.find('.');

		if (lowerDotPos == std::string::npos || upperDotPos == std::string::npos) {
			temp.refine();
			continue;
		}

		// Extract decimal parts
		const std::string lowerDecimal = lowerString.substr(lowerDotPos + 1);
		const std::string upperDecimal = upperString.substr(upperDotPos + 1);

		// Count matching digits
		int matchingDigits = 0;
		const size_t minSize = lowerDecimal.size() < upperDecimal.size() ? lowerDecimal.size() : upperDecimal.size();
		for (size_t i = 0; i < minSize; ++i) {
			if (lowerDecimal[i] == upperDecimal[i]) {
				matchingDigits++;
			}
			else {
				break;
			}
		}

		if (matchingDigits >= precision) {
			// Return result with trailing zeros removed
			std::string result = lowerString;
			result.erase(result.find_last_not_of('0') + 1);
			if (result.back() == '.') {
				result.pop_back();
			}
			return result;
		}

		temp.refine();
	}
}

