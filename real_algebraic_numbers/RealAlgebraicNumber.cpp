#include "RealAlgebraicNumber.h"
#include <cmath>
#include <algorithm>
#include <iostream>

#include "MyTimer.h"
#include "MyMatrix.h"

long long binomialCoeff(const int n, int k) {
	PROFILE_FUNCTION
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
	long long res = 1;
	for (int i = 1; i <= k; ++i) {
		res = res * (n - i + 1) / i;
	}
	return res;
}

MyMatrix<Polynomial> constructSylvesterMatrixForSum(const Polynomial& p, const Polynomial& q) {
	PROFILE_FUNCTION
	if (p.isZero() || q.isZero()) {
		// Special handling for zero polynomials: sum would be the other poly or zero.
		// For resultant, result would typically be zero.
		return MyMatrix<Polynomial>(0, 0);
	}

	int m = p.degree;
	int n = q.degree;

	int matSize = m + n;
	MyMatrix<Polynomial> sylvester(matSize, matSize);

	// Fill n rows for p(x) (shifted coefficients of p)
	for (int i = 0; i < n; ++i) {
		// rows 0 to n-1
		for (int j = 0; j <= m; ++j) {
			// p.coeff(j) is a Rational (coefficient of x^j in p).
			// It becomes a degree 0 Polynomial in z for the matrix entry.
			sylvester.set_element(i, i + j, Polynomial(p.coeff(j)));
		}
	}

	// Fill m rows for G(x) = q(z-x) (shifted coefficients of G)
	// We need to compute the coefficients of G(x) as polynomials in z.
	// g_s(z) = Sum_{j=s to n} q_j * binomialCoeff(j, s) * (-1)^s * z^(j-s)

	// Polynomials in x that form rows of G(x) part
	// Each of these coeff_of_xs[s] is a Polynomial in z.
	std::vector<Polynomial> coeffs_of_G_in_x(n + 1); // coeffs_of_G_in_x[s] holds g_s(z)

	Polynomial zPoly({Rational(0), Rational(1)}); // Represents the variable 'z'

	for (int s = 0; s <= n; ++s) {
		// Iterate through powers of x (x^s) for G(x)
		Polynomial g_s_z_poly; // This will accumulate terms for the coefficient of x^s in G(x)

		for (int j = s; j <= n; ++j) {
			// Summation over j
			Rational q_j = q.coeff(j); // Coefficient of y^j in Q(y)
			if (q_j == 0) continue;

			long long termBinomial = binomialCoeff(j, s);
			Rational term_minus_one_power = ((s % 2 == 0) ? Rational(1) : Rational(-1));

			// Calculate the Rational coefficient part for this term: q_j * C(j,s) * (-1)^s
			Rational scalarCoeff = q_j * Rational(termBinomial) * term_minus_one_power;

			// Create the z-polynomial part: z^(j-s)
			Polynomial z_power_poly;
			if (j - s == 0) {
				z_power_poly = Polynomial(Rational(1)); // z^0 = 1
			}
			else {
				// Equivalent to z_poly.power(j-s) but Polynomial class doesn't have power yet.
				// Manually create polynomial representing z^(j-s)
				std::vector<Rational> zCoeffs(j - s + 1, Rational(0));
				zCoeffs[j - s] = Rational(1);
				z_power_poly = Polynomial(zCoeffs);
			}

			// Add to g_s_z_poly: scalar_coeff * z_power_poly
			g_s_z_poly += (z_power_poly * scalarCoeff);
		}
		coeffs_of_G_in_x[s] = g_s_z_poly;
	}

	// Now fill the m rows for G(x) in the Sylvester matrix
	for (int i = 0; i < m; ++i) {
		// rows n to n+m-1
		for (int j = 0; j <= n; ++j) {
			// coeffs_of_G_in_x[j] is the coefficient of x^j in G(x), which is a Polynomial in z.
			sylvester.set_element(n + i, i + j, coeffs_of_G_in_x[j]);
		}
	}
	return sylvester;
}

MyMatrix<Polynomial> constructSylvesterMatrixForProduct(const Polynomial& p, const Polynomial& q) {
	PROFILE_FUNCTION
	const int m = p.degree;
	const int n = q.degree;

	if (m == -1 || n == -1) {
		// If one of the polynomials is zero
		return MyMatrix<Polynomial>(0, 0); // Or handle as appropriate, usually det is 0
	}

	const int matSize = m + n;
	MyMatrix<Polynomial> sylvester(matSize, matSize);

	// Fill rows for p(x)
	// The coefficients of p(x) are Rational. We convert them to Polynomials of degree 0 in 'z'.
	for (int i = 0; i < n; ++i) {
		// n rows for p(x)
		for (int j = 0; j <= m; ++j) {
			// p.coeff(j) is the coefficient of x^j.
			// In the Sylvester matrix, terms are filled in descending power.
			// Example: for p(x) = p2 x^2 + p1 x + p0
			// Row 0: p2 p1 p0 0 0
			// Row 1: 0 p2 p1 p0 0
			// (If m=2, n=3, mat_size=5)
			// Column index for x^(m-k) in row i is i+k.
			// Coefficient for x^k is p.coeff(k)
			sylvester.set_element(i, i + j, Polynomial(p.coeff(j)));
		}
	}

	// Fill rows for Q_zx(x) = x^n * q(z/x)
	// Q_zx(x) = Sum_{j=0 to n} (q_j * z^j) * x^(n-j)
	// The coefficient of x^k in Q_zx(x) is q_{n-k} * z^{n-k}.
	// This requires a Polynomial in 'z' for the coefficient.
	// Let's assume a global 'z' variable or convention for simplicity.
	// For this example, we'll manually construct the polynomial terms for 'z'.
	// Here, `Polynomial` represents polynomials in `z` for the matrix elements.
	for (int i = 0; i < m; ++i) {
		// m rows for Q_zx(x)
		for (int j = 0; j <= n; ++j) {
			// j is the original power of y in q(y)
			// The term in Q_zx(x) is (q_j * z^j) * x^(n-j)
			// So, the coefficient of x^(n-j) is q_j * z^j.
			// The index of this coefficient in the `coeffs` vector of a Polynomial `P_in_x`
			// would be `n-j`.
			// We need to place this at column `i + (n-j)` in the Sylvester matrix.

			// Build the coefficient (a polynomial in 'z'): q_j * z^j
			Rational q_j = q.coeff(j); // The j-th coefficient of q(y)

			auto coeff_in_z = Polynomial(Rational(0)); // Start with 0
			if (q_j != 0) {
				std::vector<Rational> z_poly_coeffs(j + 1, Rational(0));
				z_poly_coeffs[j] = q_j; // Coefficient of z^j is q_j
				coeff_in_z = Polynomial(z_poly_coeffs);
			}

			// Place this polynomial in z into the matrix at the correct position.
			// The index for the power of x in Q_zx is (n-j).
			// This term contributes to column `i + (n-j)`.
			sylvester.set_element(n + i, i + (n - j), coeff_in_z);
		}
	}
	return sylvester;
}

MyMatrix<Polynomial> constructSylvesterMatrixForPower(const Polynomial& p, const int k) {
	PROFILE_FUNCTION
	if (k <= 0) {
		throw std::invalid_argument("Exponent k must be positive for power resultant construction.");
	}
	if (p.isZero()) {
		return MyMatrix<Polynomial>(0, 0); // Or handle as zero polynomial
	}

	const int m = p.degree;
	const int n = k; // Degree of G(x) = -x^k + z

	const int matSize = m + n;
	MyMatrix<Polynomial> sylvester(matSize, matSize);

	// Fill n rows for p(x) (shifted coefficients of p)
	for (int i = 0; i < n; ++i) {
		// rows 0 to n-1
		for (int j = 0; j <= m; ++j) {
			// p.coeff(j) is the coefficient of x^j.
			// Place p.coeff(j) at column i+j in matrix.
			sylvester.set_element(i, i + j, Polynomial(p.coeff(j)));
			// Polynomial(Rational) creates a degree 0 poly in 'z'
		}
	}

	// Fill m rows for G(x) = -x^k + z (shifted coefficients of G)
	// Coefficients of G(x) as a polynomial in x:
	// G.coeff(k) = -1 (constant, for x^k term)
	// G.coeff(0) = z (this is a Polynomial in z, so {Rational(0), Rational(1)} for z^1)
	// G.coeff(other_idx) = 0

	const Polynomial zPoly({Rational(0), Rational(1)}); // Represents the variable 'z' itself

	for (int i = 0; i < m; ++i) {
		// rows n to n+m-1
		for (int j = 0; j <= n; ++j) {
			Polynomial G_coeff_at_j; // This will be a Polynomial in 'z'

			if (j == k) {
				// Coefficient for x^k is -1
				G_coeff_at_j = Polynomial(Rational(-1));
			}
			else if (j == 0) {
				// Coefficient for x^0 is z
				G_coeff_at_j = zPoly;
			}
			else {
				// Other coefficients are 0
				G_coeff_at_j = Polynomial(Rational(0));
			}
			// Place G_coeff_at_j at column n+i+j in matrix (n is starting row for G-block)
			sylvester.set_element(n + i, i + j, G_coeff_at_j);
		}
	}
	return sylvester;
}

RealAlgebraicNumber::RealAlgebraicNumber()
	: polynomial({0}), interval({.lowerBound = 0.0, .upperBound = 0.0}) {}

RealAlgebraicNumber::RealAlgebraicNumber(const Polynomial& polynomial, const Interval& interval)
//: polynomial(std::move(polynomial)), interval(std::move(interval))
{
	this->polynomial = polynomial;
	this->interval = interval;
	//normalize();
	while ((this->interval.lowerBound - this->interval.upperBound).abs() > 0.01) {
		refine();
	}
	this->polynomial.normalize(this->interval.lowerBound, this->interval.upperBound);
}

RealAlgebraicNumber::RealAlgebraicNumber(const Polynomial& polynomial, const Rational& lowerBound,
                                         const Rational& upperBound)
//: polynomial(std::move(polynomial)), interval{.lowerBound = lowerBound, .upperBound = upperBound}
{
	this->polynomial = polynomial;
	this->interval.lowerBound = lowerBound;
	this->interval.upperBound = upperBound;
	//normalize();
	while ((this->interval.lowerBound - this->interval.upperBound).abs() > 0.01) {
		refine();
	}
	this->polynomial.normalize(this->interval.lowerBound, this->interval.upperBound);
}

RealAlgebraicNumber::RealAlgebraicNumber(const std::vector<Rational>& coefficients, const Rational& lowerBound,
                                         const Rational& upperBound)
	: polynomial(coefficients), interval{.lowerBound = lowerBound, .upperBound = upperBound} {
	//normalize();
	while ((this->interval.lowerBound - this->interval.upperBound).abs() > 0.01) {
		refine();
	}
	this->polynomial.normalize(this->interval.lowerBound, this->interval.upperBound);
}

RealAlgebraicNumber::RealAlgebraicNumber(const int value) {
	this->polynomial = Polynomial({-value, 1});
	this->interval.lowerBound = value;
	this->interval.upperBound = value;
	this->polynomial.normalize(this->interval.lowerBound, this->interval.upperBound);
}

RealAlgebraicNumber::RealAlgebraicNumber(const double value) {
	const Rational rationalValue = value;
	this->polynomial = Polynomial({-rationalValue, 1});
	this->interval.lowerBound = rationalValue;
	this->interval.upperBound = rationalValue;
	this->polynomial.normalize(this->interval.lowerBound, this->interval.upperBound);
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
	return {polynomial.reflectY(), -interval.upperBound, -interval.lowerBound};
}

Rational minRational(const Rational& r1, const Rational& r2, const Rational& r3, const Rational& r4) {
	PROFILE_FUNCTION
	Rational min = r1;
	if (r2 < min) min = r2;
	if (r3 < min) min = r3;
	if (r4 < min) min = r4;
	return min;
}

Rational maxRational(const Rational& r1, const Rational& r2, const Rational& r3, const Rational& r4) {
	PROFILE_FUNCTION
	Rational max = r1;
	if (r2 > max) max = r2;
	if (r3 > max) max = r3;
	if (r4 > max) max = r4;
	return max;
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
	RealAlgebraicNumber otherCopy = other;
	RealAlgebraicNumber thisCopy = *this;

	MyMatrix<Polynomial> sylvester_mat = constructSylvesterMatrixForSum(this->polynomial, other.polynomial);
	Polynomial f3 = sylvester_mat.determinant();

	std::vector<Polynomial> sturm;
	if (f3.sturm_sequence.empty()) {
		sturm = f3.sturmSequence(f3);
	}
	else {
		sturm = f3.sturm_sequence;
	}

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

	const Interval sumInterval = {.lowerBound = l3, .upperBound = r3};

	this->polynomial = f3;
	this->interval = sumInterval;

	return *this;
}

RealAlgebraicNumber RealAlgebraicNumber::operator-=(const RealAlgebraicNumber& other) {
	PROFILE_FUNCTION
	*this += -other;
	/*this->polynomial = result.polynomial;
	this->interval = result.interval;*/
	return *this;
}

RealAlgebraicNumber RealAlgebraicNumber::operator*=(const RealAlgebraicNumber& other) {
	PROFILE_FUNCTION
	RealAlgebraicNumber otherCopy = other;
	RealAlgebraicNumber thisCopy = *this;

	// make sure they are normalised
	otherCopy.polynomial.normalize(otherCopy.interval.lowerBound, otherCopy.interval.upperBound);
	thisCopy.polynomial.normalize(thisCopy.interval.lowerBound, thisCopy.interval.upperBound);

	MyMatrix<Polynomial> sylvester_prod_mat = constructSylvesterMatrixForProduct(
		thisCopy.polynomial, otherCopy.polynomial);

	Polynomial f3 = sylvester_prod_mat.determinant();

	/*std::vector<Polynomial> sturm;
	if (f3.sturm_sequence.empty()) {
		sturm = f3.sturmSequence(f3);
	}
	else {
		sturm = f3.sturm_sequence;
	}*/

	Rational l3 = minRational(
		thisCopy.interval.lowerBound * otherCopy.interval.lowerBound,
		thisCopy.interval.lowerBound * otherCopy.interval.upperBound,
		thisCopy.interval.upperBound * otherCopy.interval.lowerBound,
		thisCopy.interval.upperBound * otherCopy.interval.upperBound);
	Rational r3 = maxRational(
		thisCopy.interval.lowerBound * otherCopy.interval.lowerBound,
		thisCopy.interval.lowerBound * otherCopy.interval.upperBound,
		thisCopy.interval.upperBound * otherCopy.interval.lowerBound,
		thisCopy.interval.upperBound * otherCopy.interval.upperBound);

	f3.normalize(l3, r3);

	std::vector<Polynomial> sturm;
	if (f3.sturm_sequence.empty()) {
		sturm = f3.sturmSequence(f3);
	}
	else {
		sturm = f3.sturm_sequence;
	}

	// Refine interval until exactly one root is isolated
	while (variationCount(sturm, l3) - variationCount(sturm, r3) > 1) {
		//auto f1 = RealAlgebraicNumber(polynomial, interval.lowerBound, r3);
		thisCopy.refine();
		otherCopy.refine();
		l3 = minRational(
			thisCopy.interval.lowerBound * otherCopy.interval.lowerBound,
			thisCopy.interval.lowerBound * otherCopy.interval.upperBound,
			thisCopy.interval.upperBound * otherCopy.interval.lowerBound,
			thisCopy.interval.upperBound * otherCopy.interval.upperBound);
		r3 = maxRational(
			thisCopy.interval.lowerBound * otherCopy.interval.lowerBound,
			thisCopy.interval.lowerBound * otherCopy.interval.upperBound,
			thisCopy.interval.upperBound * otherCopy.interval.lowerBound,
			thisCopy.interval.upperBound * otherCopy.interval.upperBound);
	}

	f3.normalize(l3, r3);

	this->polynomial = f3;
	this->interval = {.lowerBound = l3, .upperBound = r3};

	return *this;
}

RealAlgebraicNumber RealAlgebraicNumber::operator/=(const RealAlgebraicNumber& other) {
	PROFILE_FUNCTION
	*this *= other.inverse();
	/*this->polynomial = result.polynomial;
	this->interval = result.interval;*/
	return *this;
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
	PROFILE_FUNCTION
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
	PROFILE_FUNCTION
	return other < *this;
}

bool RealAlgebraicNumber::operator<=(const RealAlgebraicNumber& other) const {
	PROFILE_FUNCTION
	return !(*this > other);
}

bool RealAlgebraicNumber::operator>=(const RealAlgebraicNumber& other) const {
	PROFILE_FUNCTION
	return !(*this < other);
}

RealAlgebraicNumber RealAlgebraicNumber::inverse() const {
	PROFILE_FUNCTION
	std::vector<Rational> inverseCo(polynomial.coefficients.rbegin(), polynomial.coefficients.rend());
	Rational inverseLower = interval.upperBound.inverse();
	Rational inverseUpper = interval.lowerBound.inverse();
	return {inverseCo, inverseLower, inverseUpper};
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

	if (n % 2 != 0 && interval.lowerBound < 0) {
		l3 = -interval.upperBound.abs().sqrt(n) - 0.01;
		r3 = -interval.lowerBound.abs().sqrt(n) + 0.01;
	}
	else {
		l3 = interval.lowerBound.sqrt(n) - 0.01;
		r3 = interval.upperBound.sqrt(n) + 0.01;
	}

	//while (variationCount(sturm, l3) - variationCount(sturm, r3) > 1) {
	//	//auto f1 = RealAlgebraicNumber(polynomial, interval.lowerBound, r3);
	//	this->refine();
	//	if (n % 2 != 0 && interval.lowerBound < 0)
	//	{
	//		l3 = -interval.upperBound.abs().sqrt(n);
	//		r3 = -interval.lowerBound.abs().sqrt(n);
	//	}
	//	else
	//	{
	//		l3 = interval.lowerBound.sqrt(n);
	//		r3 = interval.upperBound.sqrt(n);
	//	}
	//}

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

	Rational l3;
	Rational r3;

	// Handle k > 0
	if (this->interval.lowerBound >= 0) {
		// Both endpoints are non-negative
		l3 = this->interval.lowerBound.pow(n);
		r3 = this->interval.upperBound.pow(n);
	}
	else if (this->interval.upperBound <= Rational(0)) {
		// Both endpoints are non-positive
		if (n % 2 == 0) {
			// Even power: result is positive, order flips
			l3 = this->interval.upperBound.pow(n);
			r3 = this->interval.lowerBound.pow(n);
		}
		else {
			// Odd power: result stays negative, order same
			l3 = this->interval.lowerBound.pow(n);
			r3 = this->interval.upperBound.pow(n);
		}
	}
	else {
		// Interval spans zero (lower < 0 < upper)
		auto minVal = Rational(0);
		Rational lowerPowerN = this->interval.lowerBound.pow(n);
		Rational upperPowerN = this->interval.upperBound.pow(n);
		Rational maxVal = lowerPowerN > upperPowerN ? lowerPowerN : upperPowerN;

		// If k is even, the minimum value will be 0.
		// If k is odd, the minimum value will be lower.power(k).
		if (n % 2 != 0) {
			// Even power, like [-2, 3]^2 -> [0, 9]
			minVal = lowerPowerN;
			maxVal = upperPowerN;
		}

		l3 = minVal;
		r3 = maxVal;
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
	PROFILE_FUNCTION
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
	if (a + b - b != a || c + d - c != d) // TODO: c+d-c seems to be wrong, check operator+ ?
		std::cout << "operator+ and operator- not working together!\n";
	if (a * b != i || c * d != l)
		std::cout << "operator* not working!\n";
	if (a / b != j || c / d != k)
		std::cout << "operator/ not working!\n";
	if ((a * b) / b != a || (c * d) / c != d) // TODO: c*d/c seems to be wrong, check operator* ?
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
		// TODO: a.inverse().inverse() == a (factoring not complete)
		std::cout << "double inverse() not working!\n";
	if (a.sqrt() != c || b.sqrt() != d)
		std::cout << "sqrt() not working!\n";
	if (c.pow(2) != a || d.pow(2) != b)
		std::cout << "pow() not working!\n";
	if (b.pow(2).sqrt() != b || d.sqrt().pow(2) != d) // TODO: d.sqrt().pow(2) == d (pow issue? or sqrt interval)
		std::cout << "sqrt() and pow() not working together! (1)\n";
	if (a.sqrt().pow(2) != a || c.pow(2).sqrt() != c)
		std::cout << "sqrt() and pow() not working together! (2)\n";

	std::cout << "\n";
}

int RealAlgebraicNumber::countSignVariations(const std::vector<Rational>& sequence) {
	PROFILE_FUNCTION
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

	for (const auto& sequence : sturm) {
		evaluations.push_back(evaluatePoly(sequence.coefficients, x));
	}

	return countSignVariations(evaluations);
}

int RealAlgebraicNumber::intervalToOrder() {
	PROFILE_FUNCTION
	Rational fInf = 0;
	for (const Rational& coeff : polynomial.coefficients) {
		if (fInf < coeff.abs())
			fInf = coeff.abs();
		//fInf = std::max(fInf, coeff.abs());
	}

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
	else {
		sturm = polynomial.sturm_sequence;
	}

	const int varL = variationCount(sturm, interval.lowerBound);
	const int varNegP = variationCount(sturm, -p);
	const int varP = variationCount(sturm, p);
	//int varR = variationCount(polynomial.coefficients, interval.upperBound);

	if (varL > varNegP) {
		interval.upperBound = -p; // Adjust interval to exclude zero
	}
	else if (varNegP > varP) {
		interval.lowerBound = 0;
		interval.upperBound = 0; // Indicating zero root
	}
	else {
		interval.lowerBound = p; // Adjust to keep positive sign interval
	}
}

void RealAlgebraicNumber::refine() {
	PROFILE_FUNCTION
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
	if (this->polynomial.degree == 1 && this->polynomial.coefficients[1] == 1) {
		Rational zerothCoeff = this->polynomial.coefficients[0];
		return (zerothCoeff * -1).toDecimalString(0);
	}

	RealAlgebraicNumber temp = *this;

	while (true) {
		std::string lowerString = temp.interval.lowerBound.toDecimalString(precision);
		std::string upperString = temp.interval.upperBound.toDecimalString(precision);

		// only keeping part after '.'
		std::string lowerStringDecimal = lowerString.substr(lowerString.find('.') + 1);
		std::string upperStringDecimal = upperString.substr(upperString.find('.') + 1);

		int matchingNumbers = 0;
		for (size_t i = 0; i < lowerStringDecimal.size(); i++) {
			if (lowerStringDecimal[i] == upperStringDecimal[i]) {
				matchingNumbers++;
			}
			else {
				break;
			}
		}

		if (matchingNumbers >= precision) {
			// remove trailing zeros
			lowerString.erase(lowerString.find_last_not_of('0') + 1);
			return lowerString;
		}

		temp.refine();
	}
}
