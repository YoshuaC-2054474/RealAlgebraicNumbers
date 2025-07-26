#include "RealAlgebraicNumber.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <utility>

#include <limits> // For numeric_limits
#include <iomanip> // For std::setprecision

#include "MyTimer.h"
#include "MyMatrix.h"

// Helper macro for testing assertions
#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        std::cerr << "FAILED: " << message /*<< " at " << __FILE__ << ":" << __LINE__ */ << std::endl; \
        failed_tests++; \
    } /*else { \
        std::cout << "PASSED: " << message << std::endl; \
    }*/

// Helper function to print a RealAlgebraicNumber
void printRAN(const std::string& name, const RealAlgebraicNumber& ran) {
	std::cout << name << " = " << ran.toString() << " (Decimal: " << ran.toDecimalString(15) << ")" << std::endl;
}

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
	: polynomial({ 0, 1 }), interval({ .lowerBound = 0.0, .upperBound = 0.0 }) {
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

RealAlgebraicNumber::RealAlgebraicNumber(const Rational& value)
	: polynomial({ -value, 1 }), interval{ .lowerBound = value, .upperBound = value } {
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

RealAlgebraicNumber RealAlgebraicNumber::abs() const {
	PROFILE_FUNCTION
	// If the polynomial is zero, return zero
	if (this->isZero()) {
		return { Polynomial({0,1}), 0, 0};
	}
	// If the polynomial is positive, return it as is
	if (this->isPositive()) {
		return *this;
	}
	// If the polynomial is negative, reflect it across the x-axis
	return -(*this);
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

bool RealAlgebraicNumber::isPositive() const {
	// Check if the polynomial is positive in the interval
	RealAlgebraicNumber thisCopy = *this;
	thisCopy.normalize();

	if (thisCopy.interval.lowerBound > 0) {
		return true;
	}
	if (thisCopy.interval.upperBound < 0) {
		return false;
	}
	return false;
}

bool RealAlgebraicNumber::isNegative() const {
	// Check if the polynomial is negative in the interval
	RealAlgebraicNumber thisCopy = *this;
	thisCopy.normalize();
	if (thisCopy.interval.upperBound < 0) {
		return true;
	}
	if (thisCopy.interval.lowerBound > 0) {
		return false;
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
	auto adivb = a / b;
	std::cout << "a / b = " << adivb.toString() << '\n';
	std::cout << "j = " << j.toString() << '\n';
	auto cdivd = c / d;
	std::cout << "c / d = " << cdivd.toString() << '\n';
	std::cout << "k = " << k.toString() << '\n';
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

void RealAlgebraicNumber::extensiveTest() {
	//std::cout << "--- Starting Extensive RealAlgebraicNumber Tests ---" << std::endl;
	int failed_tests = 0;

	// Call the existing basic operator tests first
	//std::cout << "\n--- Running Basic Operator Tests (from testOperators()) ---" << std::endl;
	testOperators();
	//std::cout << "--- Basic Operator Tests Finished ---" << std::endl;

	//std::cout << "\n--- Testing Constructors ---" << std::endl;
	RealAlgebraicNumber ran_default;
	TEST_ASSERT(ran_default.isZero(), "Default constructor creates zero");
	//printRAN("ran_default", ran_default);

	RealAlgebraicNumber ran_int_pos(5);
	TEST_ASSERT(ran_int_pos == RealAlgebraicNumber(Polynomial({ -5, 1 }), 5, 5), "Integer constructor (positive)");
	//printRAN("ran_int_pos", ran_int_pos);

	RealAlgebraicNumber ran_int_neg(-3);
	TEST_ASSERT(ran_int_neg == RealAlgebraicNumber(Polynomial({ 3, 1 }), -3, -3), "Integer constructor (negative)");
	//printRAN("ran_int_neg", ran_int_neg);

	RealAlgebraicNumber ran_int_zero(0);
	TEST_ASSERT(ran_int_zero.isZero(), "Integer constructor (zero)");
	//printRAN("ran_int_zero", ran_int_zero);

	RealAlgebraicNumber ran_double_pos(2.5);
	TEST_ASSERT(ran_double_pos == RealAlgebraicNumber({ Rational(-5, 2), 1 }, Rational(25, 10), Rational(25, 10)), "Double constructor (positive)");
	//printRAN("ran_double_pos", ran_double_pos);

	RealAlgebraicNumber ran_double_neg(-1.75);
	TEST_ASSERT(ran_double_neg == RealAlgebraicNumber({ Rational(7, 4), 1 }, Rational(-175, 100), Rational(-175, 100)), "Double constructor (negative)");
	//printRAN("ran_double_neg", ran_double_neg);

	RealAlgebraicNumber ran_double_zero(0.0);
	TEST_ASSERT(ran_double_zero.isZero(), "Double constructor (zero)");
	//printRAN("ran_double_zero", ran_double_zero);

	// Test a known irrational number: sqrt(2)
	// Polynomial: x^2 - 2 = 0, Interval: [1.4, 1.5]
	RealAlgebraicNumber sqrt2_poly_interval(Polynomial({ -2, 0, 1 }), Rational(14, 10), Rational(15, 10));
	TEST_ASSERT(sqrt2_poly_interval.toDecimalString(5) == "1.41421", "Polynomial+Interval constructor (sqrt(2))");
	//printRAN("sqrt2_poly_interval", sqrt2_poly_interval);

	// Test a known irrational number: cube root of 3
	// Polynomial: x^3 - 3 = 0, Interval: [1.4, 1.5]
	RealAlgebraicNumber cbrt3_coeffs_interval(Polynomial({ -3, 0, 0, 1 }), Rational(14, 10), Rational(15, 10));
	TEST_ASSERT(cbrt3_coeffs_interval.toDecimalString(5) == "1.44224", "Coeffs+Interval constructor (cbrt(3))");
	//printRAN("cbrt3_coeffs_interval", cbrt3_coeffs_interval);

	//std::cout << "\n--- Testing Arithmetic Operations ---" << std::endl;
	RealAlgebraicNumber two(2);
	RealAlgebraicNumber three(3);
	RealAlgebraicNumber minus_five(-5);
	RealAlgebraicNumber half(Rational(1, 2));

	// Addition
	RealAlgebraicNumber sum1 = two + three;
	TEST_ASSERT(sum1 == RealAlgebraicNumber(5), "2 + 3 == 5");
	//printRAN("sum1", sum1);

	RealAlgebraicNumber sum2 = two + minus_five;
	TEST_ASSERT(sum2 == RealAlgebraicNumber(-3), "2 + (-5) == -3");
	//printRAN("sum2", sum2);

	RealAlgebraicNumber sum3 = minus_five + two;
	TEST_ASSERT(sum3 == RealAlgebraicNumber(-3), "(-5) + 2 == -3");
	//printRAN("sum3", sum3);

	RealAlgebraicNumber sum_zero = two + RealAlgebraicNumber(-2);
	TEST_ASSERT(sum_zero.isZero(), "2 + (-2) == 0");
	//printRAN("sum_zero", sum_zero);

	// Subtraction
	RealAlgebraicNumber diff1 = three - two;
	TEST_ASSERT(diff1 == RealAlgebraicNumber(1), "3 - 2 == 1");
	//printRAN("diff1", diff1);

	RealAlgebraicNumber diff2 = two - three;
	TEST_ASSERT(diff2 == RealAlgebraicNumber(-1), "2 - 3 == -1");
	//printRAN("diff2", diff2);

	RealAlgebraicNumber diff_self = two - two;
	TEST_ASSERT(diff_self.isZero(), "2 - 2 == 0");
	//printRAN("diff_self", diff_self);

	// Unary Negation
	RealAlgebraicNumber neg_two = -two;
	TEST_ASSERT(neg_two == RealAlgebraicNumber(-2), "-(2) == -2");
	//printRAN("neg_two", neg_two);

	RealAlgebraicNumber neg_minus_five = -minus_five;
	TEST_ASSERT(neg_minus_five == RealAlgebraicNumber(5), "-(-5) == 5");
	//printRAN("neg_minus_five", neg_minus_five);

	// Multiplication
	RealAlgebraicNumber prod1 = two * three;
	TEST_ASSERT(prod1 == RealAlgebraicNumber(6), "2 * 3 == 6");
	//printRAN("prod1", prod1);

	RealAlgebraicNumber prod2 = two * minus_five;
	TEST_ASSERT(prod2 == RealAlgebraicNumber(-10), "2 * (-5) == -10");
	//printRAN("prod2", prod2);

	RealAlgebraicNumber prod_zero = two * RealAlgebraicNumber(0);
	TEST_ASSERT(prod_zero.isZero(), "2 * 0 == 0");
	//printRAN("prod_zero", prod_zero);

	RealAlgebraicNumber prod_one = two * RealAlgebraicNumber(1);
	TEST_ASSERT(prod_one == two, "2 * 1 == 2");
	//printRAN("prod_one", prod_one);

	// Division
	RealAlgebraicNumber div1 = three / two;
	std::cout << div1.toString() << std::endl;
	auto test = half + RealAlgebraicNumber(1);
	std::cout << test.toString() << std::endl;
	TEST_ASSERT(div1 == half + RealAlgebraicNumber(1), "3 / 2 == 1.5"); // 1.5 is 1 + 0.5
	//printRAN("div1", div1);

	RealAlgebraicNumber div2 = minus_five / two;
	TEST_ASSERT(div2 == RealAlgebraicNumber(Rational(-5, 2)), "(-5) / 2 == -2.5");
	//printRAN("div2", div2);

	RealAlgebraicNumber div_self = two / two;
	TEST_ASSERT(div_self == RealAlgebraicNumber(1), "2 / 2 == 1");
	//printRAN("div_self", div_self);

	// Division by zero - should ideally be handled by Rational class or throw
	// For now, assuming Rational handles division by zero gracefully (e.g., infinity or error state)
	// RealAlgebraicNumber div_by_zero = two / RealAlgebraicNumber(0); // This will likely crash or produce invalid results

	//std::cout << "\n--- Testing Compound Assignment Operations ---" << std::endl;
	RealAlgebraicNumber a_val(10);
	RealAlgebraicNumber b_val(4);

	RealAlgebraicNumber test_a = a_val;
	test_a += b_val;
	TEST_ASSERT(test_a == RealAlgebraicNumber(14), "a += b");
	//printRAN("test_a after +=", test_a);

	test_a = a_val;
	test_a -= b_val;
	TEST_ASSERT(test_a == RealAlgebraicNumber(6), "a -= b");
	//printRAN("test_a after -=", test_a);

	test_a = a_val;
	test_a *= b_val;
	TEST_ASSERT(test_a == RealAlgebraicNumber(40), "a *= b");
	//printRAN("test_a after *=", test_a);

	test_a = a_val;
	test_a /= b_val;
	TEST_ASSERT(test_a == RealAlgebraicNumber(Rational(10, 4)), "a /= b");
	//printRAN("test_a after /=", test_a);

	//std::cout << "\n--- Testing Comparison Operators ---" << std::endl;
	RealAlgebraicNumber x(10);
	RealAlgebraicNumber y(10);
	RealAlgebraicNumber z(12);
	RealAlgebraicNumber neg_x(-10);

	TEST_ASSERT(x == y, "x == y (equal values)");
	TEST_ASSERT(x != z, "x != z (unequal values)");
	TEST_ASSERT(x < z, "x < z");
	TEST_ASSERT(z > x, "z > x");
	TEST_ASSERT(x <= y, "x <= y (equal)");
	TEST_ASSERT(x <= z, "x <= z (less than)");
	TEST_ASSERT(x >= y, "x >= y (equal)");
	TEST_ASSERT(z >= x, "z >= x (greater than)");

	TEST_ASSERT(neg_x < x, "-x < x");
	TEST_ASSERT(x > neg_x, "x > -x");
	TEST_ASSERT(neg_x <= x, "-x <= x");
	TEST_ASSERT(x >= neg_x, "x >= -x");

	// Test equality for numbers with same polynomial but different initial intervals
	RealAlgebraicNumber sqrt2_a(Polynomial({ -2, 0, 1 }), Rational(14, 10), Rational(15, 10)); // ~1.414
	RealAlgebraicNumber sqrt2_b(Polynomial({ -2, 0, 1 }), Rational(141, 100), Rational(142, 100)); // ~1.41
	TEST_ASSERT(sqrt2_a == sqrt2_b, "sqrt(2) == sqrt(2) (different initial intervals)");
	//printRAN("sqrt2_a", sqrt2_a);
	//printRAN("sqrt2_b", sqrt2_b);

	RealAlgebraicNumber minus_sqrt2(Polynomial({ -2, 0, 1 }), Rational(-15, 10), Rational(-14, 10)); // ~-1.414
	TEST_ASSERT(sqrt2_a != minus_sqrt2, "sqrt(2) != -sqrt(2)");
	TEST_ASSERT(minus_sqrt2 < sqrt2_a, "-sqrt(2) < sqrt(2)");

	//std::cout << "\n--- Testing Special Functions ---" << std::endl;

	// inverse()
	RealAlgebraicNumber inv_two = two.inverse();
	TEST_ASSERT(inv_two == half, "inverse(2) == 0.5");
	//printRAN("inv_two", inv_two);

	RealAlgebraicNumber inv_half = half.inverse();
	TEST_ASSERT(inv_half == two, "inverse(0.5) == 2");
	//printRAN("inv_half", inv_half);

	RealAlgebraicNumber inv_sqrt2 = sqrt2_a.inverse();
	RealAlgebraicNumber expected_inv_sqrt2(Polynomial({ 1, 0, -2 }), Rational(7, 10), Rational(8, 10)); // 1/sqrt(2) = sqrt(2)/2
	TEST_ASSERT(inv_sqrt2 == expected_inv_sqrt2, "inverse(sqrt(2)) == sqrt(2)/2");
	//printRAN("inv_sqrt2", inv_sqrt2);

	RealAlgebraicNumber double_inv_sqrt2 = inv_sqrt2.inverse();
	TEST_ASSERT(double_inv_sqrt2 == sqrt2_a, "inverse(inverse(sqrt(2))) == sqrt(2)");
	//printRAN("double_inv_sqrt2", double_inv_sqrt2);

	// sqrt()
	RealAlgebraicNumber four(4);
	RealAlgebraicNumber sqrt_four = four.sqrt();
	TEST_ASSERT(sqrt_four == two, "sqrt(4) == 2");
	//printRAN("sqrt_four", sqrt_four);

	RealAlgebraicNumber sqrt_two = two.sqrt();
	TEST_ASSERT(sqrt_two == sqrt2_a, "sqrt(2) == sqrt(2)");
	//printRAN("sqrt_two", sqrt_two);

	RealAlgebraicNumber eight(8);
	RealAlgebraicNumber cbrt_eight = eight.sqrt(3); // Cube root of 8
	TEST_ASSERT(cbrt_eight == two, "cbrt(8) == 2");
	//printRAN("cbrt_eight", cbrt_eight);

	RealAlgebraicNumber neg_eight(-8);
	RealAlgebraicNumber cbrt_neg_eight = neg_eight.sqrt(3); // Cube root of -8
	TEST_ASSERT(cbrt_neg_eight == RealAlgebraicNumber(-2), "cbrt(-8) == -2");
	//printRAN("cbrt_neg_eight", cbrt_neg_eight);

	// Test even root of negative number - should throw
	//std::cout << "Attempting sqrt(-4) (should throw): ";
	try {
		RealAlgebraicNumber neg_four(-4);
		RealAlgebraicNumber sqrt_neg_four = neg_four.sqrt();
		TEST_ASSERT(false, "sqrt(-4) did NOT throw exception (FAILURE)"); // Should not reach here
	}
	catch (const std::invalid_argument& e) {
		//std::cout << "PASSED: Caught expected exception: " << e.what() << std::endl;
	}
	catch (...) {
		TEST_ASSERT(false, "sqrt(-4) threw unexpected exception (FAILURE)");
	}

	// pow()
	RealAlgebraicNumber two_pow_three = two.pow(3);
	TEST_ASSERT(two_pow_three == RealAlgebraicNumber(8), "2^3 == 8");
	//printRAN("two_pow_three", two_pow_three);

	RealAlgebraicNumber two_pow_zero = two.pow(0);
	TEST_ASSERT(two_pow_zero == RealAlgebraicNumber(1), "2^0 == 1");
	//printRAN("two_pow_zero", two_pow_zero);

	RealAlgebraicNumber minus_two_pow_three = neg_two.pow(3);
	TEST_ASSERT(minus_two_pow_three == RealAlgebraicNumber(-8), "(-2)^3 == -8");
	//printRAN("minus_two_pow_three", minus_two_pow_three);

	RealAlgebraicNumber minus_two_pow_two = neg_two.pow(2);
	TEST_ASSERT(minus_two_pow_two == RealAlgebraicNumber(4), "(-2)^2 == 4");
	//printRAN("minus_two_pow_two", minus_two_pow_two);

	RealAlgebraicNumber sqrt2_pow_two = sqrt2_a.pow(2);
	TEST_ASSERT(sqrt2_pow_two == two, "sqrt(2)^2 == 2");
	//printRAN("sqrt2_pow_two", sqrt2_pow_two);

	//RealAlgebraicNumber two_pow_half = two.pow(Rational(1, 2)); // Assuming pow can take Rational exponent
	// NOTE: The current pow function only takes int exponent. This test would require a different pow overload.
	// For now, testing integer powers only.

	// isZero()
	RealAlgebraicNumber zero_val(0);
	TEST_ASSERT(zero_val.isZero(), "isZero() on 0");
	TEST_ASSERT(!two.isZero(), "isZero() on 2");
	TEST_ASSERT(!sqrt2_a.isZero(), "isZero() on sqrt(2)");

	// Test a number that is computationally zero but might not be exactly 0
	// e.g., (sqrt(2) * sqrt(2)) - 2
	RealAlgebraicNumber computed_zero = (sqrt2_a * sqrt2_a) - two;
	TEST_ASSERT(computed_zero.isZero(), "(sqrt(2)*sqrt(2)) - 2 == 0");
	//printRAN("computed_zero", computed_zero);

	//std::cout << "\n--- Testing toDecimalString with various precisions ---" << std::endl;
	RealAlgebraicNumber pi_approx(Polynomial({ -314159, 100000, 0, 0, 0, 1 }), Rational(314159, 100000), Rational(314160, 100000)); // x^5 - 3.14159 = 0, just for testing
	// This is not actually pi, but a number defined by a polynomial and interval to test string conversion.
	// Let's use a simpler irrational for `toDecimalString`
	RealAlgebraicNumber simple_sqrt2(Polynomial({ -2, 0, 1 }), Rational(1), Rational(2)); // sqrt(2)
	//printRAN("simple_sqrt2", simple_sqrt2);

	std::string s_sqrt2_5 = simple_sqrt2.toDecimalString(5);
	//std::cout << "sqrt(2) to 5 decimal places: " << s_sqrt2_5 << std::endl;
	TEST_ASSERT(s_sqrt2_5 == "1.41421", "toDecimalString(5) for sqrt(2)");

	std::string s_sqrt2_10 = simple_sqrt2.toDecimalString(10);
	//std::cout << "sqrt(2) to 10 decimal places: " << s_sqrt2_10 << std::endl;
	TEST_ASSERT(s_sqrt2_10 == "1.4142135623", "toDecimalString(10) for sqrt(2)");

	std::string s_sqrt2_1 = simple_sqrt2.toDecimalString(1);
	//std::cout << "sqrt(2) to 1 decimal place: " << s_sqrt2_1 << std::endl;
	TEST_ASSERT(s_sqrt2_1 == "1.4", "toDecimalString(1) for sqrt(2)");

	RealAlgebraicNumber one_third(Rational(1, 3));
	std::string s_one_third_5 = one_third.toDecimalString(5);
	//std::cout << "1/3 to 5 decimal places: " << s_one_third_5 << std::endl;
	TEST_ASSERT(s_one_third_5 == "0.33333", "toDecimalString(5) for 1/3");

	RealAlgebraicNumber neg_one_third(Rational(-1, 3));
	std::string s_neg_one_third_5 = neg_one_third.toDecimalString(5);
	std::cout << "-1/3 to 5 decimal places: " << s_neg_one_third_5 << std::endl;
	TEST_ASSERT(s_neg_one_third_5 == "-0.33333", "toDecimalString(5) for -1/3");


	//std::cout << "\n--- Testing Complex Chained Operations ---" << std::endl;
	// (sqrt(2) + sqrt(3))^2
	RealAlgebraicNumber sqrt3_val(Polynomial({ -3, 0, 1 }), Rational(17, 10), Rational(18, 10)); // sqrt(3)
	//printRAN("sqrt3_val", sqrt3_val);

	RealAlgebraicNumber sum_sqrt2_sqrt3 = sqrt2_a + sqrt3_val;
	//printRAN("sum_sqrt2_sqrt3", sum_sqrt2_sqrt3);
	// (sqrt(2) + sqrt(3))^2 = 2 + 3 + 2*sqrt(6) = 5 + 2*sqrt(6)
	RealAlgebraicNumber expected_val = RealAlgebraicNumber(5) + (RealAlgebraicNumber(2) * RealAlgebraicNumber(6).sqrt());
	//printRAN("expected_val", expected_val);

	RealAlgebraicNumber chained_result = sum_sqrt2_sqrt3.pow(2);
	//printRAN("chained_result", chained_result);
	TEST_ASSERT(chained_result == expected_val, "(sqrt(2) + sqrt(3))^2 == 5 + 2*sqrt(6)");


	// Test for very close numbers
	//std::cout << "\n--- Testing very close numbers ---" << std::endl;
	RealAlgebraicNumber small_diff_a(Polynomial({ -1000000001, 1000000000 }), Rational(1000000001, 1000000000), Rational(1000000001, 1000000000)); // 1.000000001
	RealAlgebraicNumber small_diff_b(Polynomial({ -1000000002, 1000000000 }), Rational(1000000002, 1000000000), Rational(1000000002, 1000000000)); // 1.000000002

	TEST_ASSERT(small_diff_a < small_diff_b, "1.000000001 < 1.000000002");
	TEST_ASSERT(small_diff_b > small_diff_a, "1.000000002 > 1.000000001");
	TEST_ASSERT(small_diff_a != small_diff_b, "1.000000001 != 1.000000002");
	//printRAN("small_diff_a", small_diff_a);
	//printRAN("small_diff_b", small_diff_b);

	// Test a number that is a root of a higher degree polynomial
	// e.g., a root of x^5 - 10x + 1 = 0
	// This requires finding an isolating interval for a specific root.
	// For demonstration, let's pick a polynomial with known rational roots for simplicity,
	// or a simple irrational one that's not just sqrt.
	// Let's test x^3 - 2 = 0 (cbrt(2))
	RealAlgebraicNumber cbrt2_val(Polynomial({ -2, 0, 0, 1 }), Rational(12, 10), Rational(13, 10)); // ~1.2599
	//printRAN("cbrt2_val", cbrt2_val);
	RealAlgebraicNumber cbrt2_cubed = cbrt2_val.pow(3);
	TEST_ASSERT(cbrt2_cubed == two, "cbrt(2)^3 == 2");
	//printRAN("cbrt2_cubed", cbrt2_cubed);

	//std::cout << "\n--- Extensive RealAlgebraicNumber Tests Finished ---" << std::endl;
	if (failed_tests > 0) {
		std::cerr << "SUMMARY: " << failed_tests << " tests FAILED." << std::endl;
	}
	else {
		std::cout << "SUMMARY: All tests PASSED." << std::endl;
	}
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

