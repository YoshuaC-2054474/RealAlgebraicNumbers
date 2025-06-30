#include "RealAlgebraicNumber.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <Eigen/Dense>

#include "MyTimer.h"
#include "MyMatrix.h"

using matrix = std::vector<std::vector<Polynomial>>;
using MatrixXp = Eigen::Matrix<Polynomial, Eigen::Dynamic, Eigen::Dynamic>;

MyMatrix<Polynomial> constructSylvesterMatrixForProduct(
	const Polynomial& P, // P(x)
	const Polynomial& Q // Q(y)
) {
	int m = P.degree;
	int n = Q.degree;

	if (m == -1 || n == -1) {
		// If one of the polynomials is zero
		return MyMatrix<Polynomial>(0, 0); // Or handle as appropriate, usually det is 0
	}

	int mat_size = m + n;
	MyMatrix<Polynomial> sylvester(mat_size, mat_size);

	// Fill rows for P(x)
	// The coefficients of P(x) are Rational. We convert them to Polynomials of degree 0 in 'z'.
	for (int i = 0; i < n; ++i) {
		// n rows for P(x)
		for (int j = 0; j <= m; ++j) {
			// P.coeff(j) is the coefficient of x^j.
			// In the Sylvester matrix, terms are filled in descending power.
			// Example: for P(x) = p2 x^2 + p1 x + p0
			// Row 0: p2 p1 p0 0 0
			// Row 1: 0 p2 p1 p0 0
			// (If m=2, n=3, mat_size=5)
			// Column index for x^(m-k) in row i is i+k.
			// Coefficient for x^k is P.coeff(k)
			sylvester.set_element(i, i + j, Polynomial(P.coeff(j)));
		}
	}

	// Fill rows for Q_zx(x) = x^n * Q(z/x)
	// Q_zx(x) = Sum_{j=0 to n} (q_j * z^j) * x^(n-j)
	// The coefficient of x^k in Q_zx(x) is q_{n-k} * z^{n-k}.
	// This requires a Polynomial in 'z' for the coefficient.
	// Let's assume a global 'z' variable or convention for simplicity.
	// For this example, we'll manually construct the polynomial terms for 'z'.
	// Here, `Polynomial` represents polynomials in `z` for the matrix elements.
	for (int i = 0; i < m; ++i) {
		// m rows for Q_zx(x)
		for (int j = 0; j <= n; ++j) {
			// j is the original power of y in Q(y)
			// The term in Q_zx(x) is (q_j * z^j) * x^(n-j)
			// So, the coefficient of x^(n-j) is q_j * z^j.
			// The index of this coefficient in the `coeffs` vector of a Polynomial `P_in_x`
			// would be `n-j`.
			// We need to place this at column `i + (n-j)` in the Sylvester matrix.

			// Build the coefficient (a polynomial in 'z'): q_j * z^j
			Rational q_j = Q.coeff(j); // The j-th coefficient of Q(y)

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

MyMatrix<Polynomial> constructSylvesterMatrixForPower(const Polynomial& P, int k) {
	if (k <= 0) {
		throw std::invalid_argument("Exponent k must be positive for power resultant construction.");
	}
	if (P.isZero()) {
		return MyMatrix<Polynomial>(0, 0); // Or handle as zero polynomial
	}

	int m = P.degree;
	int n = k; // Degree of G(x) = -x^k + z

	int mat_size = m + n;
	MyMatrix<Polynomial> sylvester(mat_size, mat_size);

	// Fill n rows for P(x) (shifted coefficients of P)
	for (int i = 0; i < n; ++i) {
		// rows 0 to n-1
		for (int j = 0; j <= m; ++j) {
			// P.coeff(j) is the coefficient of x^j.
			// Place P.coeff(j) at column i+j in matrix.
			sylvester.set_element(i, i + j, Polynomial(P.coeff(j)));
			// Polynomial(Rational) creates a degree 0 poly in 'z'
		}
	}

	// Fill m rows for G(x) = -x^k + z (shifted coefficients of G)
	// Coefficients of G(x) as a polynomial in x:
	// G.coeff(k) = -1 (constant, for x^k term)
	// G.coeff(0) = z (this is a Polynomial in z, so {Rational(0), Rational(1)} for z^1)
	// G.coeff(other_idx) = 0

	Polynomial z_poly({Rational(0), Rational(1)}); // Represents the variable 'z' itself

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
				G_coeff_at_j = z_poly;
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


MatrixXp convertToEigen(const matrix& mat) {
	PROFILE_FUNCTION
	int n = mat.size();
	MatrixXp eigenMat(n, n);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			eigenMat(i, j) = mat[i][j];
	return eigenMat;
}

static Polynomial determinantEigen(const MatrixXp& mat) {
	PROFILE_FUNCTION
	/*const int n = static_cast<int>(mat.size());
	if (n == 1) return mat[0][0];

	auto det = Polynomial();
	for (int j = 0; j < n; ++j) {
		matrix minor;
		for (int i = 1; i < n; ++i) {
			std::vector<Polynomial> row;
			for (int k = 0; k < n; ++k)
				if (k != j) row.push_back(mat[i][k]);
			minor.push_back(row);
		}
		Polynomial term = mat[0][j] * determinant(minor);
		if (j % 2) term = term * Polynomial({-1});
		det += term;
	}
	return det;*/

	const int n = mat.rows();
	if (n == 1) return mat(0, 0);

	Polynomial det;
	for (int j = 0; j < n; ++j) {
		// Create minor by excluding row 0 and column j
		MatrixXp minor(n - 1, n - 1);
		for (int row = 1; row < n; ++row) {
			int minor_col = 0;
			for (int col = 0; col < n; ++col) {
				if (col != j) {
					minor(row - 1, minor_col++) = mat(row, col);
				}
			}
		}

		// Recursive call and sign adjustment
		Polynomial term = mat(0, j) * determinantEigen(minor);
		if (j % 2) term *= Polynomial({-1}); // Alternate signs
		det += term;
	}
	return det;
}

Polynomial determinantBareiss(const matrix& input) {
	PROFILE_FUNCTION
	try {
		int n = input.size();
		if (n == 0) return Polynomial(0); // Handle empty matrix

		matrix mat = input; // Work on a copy
		Polynomial sign(1); // Track row swap signs

		Polynomial prev_pivot(1); // Initially 1 for division in step 0

		for (int k = 0; k < n; ++k) {
			// Find a non-zero pivot in column k
			int pivot_row = -1;
			for (int i = k; i < n; ++i) {
				if (mat[i][k] != Polynomial(0)) {
					pivot_row = i;
					break;
				}
			}
			if (pivot_row == -1) return Polynomial(0); // Singular matrix

			// Swap rows to bring pivot to diagonal
			if (pivot_row != k) {
				std::swap(mat[k], mat[pivot_row]);
				sign = sign * Polynomial(-1); // Invert sign
			}

			const Polynomial& pivot = mat[k][k];

			// Eliminate entries below and to the right using Bareiss' formula
			for (int i = k + 1; i < n; ++i) {
				for (int j = k + 1; j < n; ++j) {
					mat[i][j] = (mat[i][j] * pivot - mat[i][k] * mat[k][j]) / prev_pivot;
				}
			}

			prev_pivot = pivot; // Update for next iteration
		}

		// Determinant is the bottom-right element multiplied by the sign
		return mat[n - 1][n - 1] * sign;
	}
	catch (const std::exception& e) {
		std::cout << "Bareis error" << e.what() << std::endl;
	}
}

static int binom(const int n, int k) {
	PROFILE_FUNCTION
	int result = 1;
	if (k > n) return 0;
	k = k < (n - k) ? k : (n - k);
	//k = std::min(k, n - k);
	for (int i = 1; i <= k; i++) {
		result = result * (n - i + 1) / i;
	}
	return result;
}

static std::vector<Polynomial> poly_x_minus_y(const Polynomial& f) {
	PROFILE_FUNCTION
	const int n = f.degree; // degree of f(x)
	// result[j] will be a polynomial in x corresponding to the coefficient of y^j.
	std::vector<Polynomial> result(n + 1);

	// For each power j of y (from 0 to n)
	for (int j = 0; j <= n; j++) {
		// The polynomial coefficient for y^j has degree (n - j) in x.
		std::vector<Rational> poly(n - j + 1, 0); // poly[k] will be the coefficient for x^k
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

static std::vector<Polynomial> poly_x_over_y(const Polynomial& f) {
	PROFILE_FUNCTION
	const int n = f.degree; // degree of f
	// We will have (n+1) polynomials: one for each power of y, from y^0 up to y^n.
	// For j = n-i, the corresponding polynomial in x is the monomial f[i]*x^i.
	std::vector<Polynomial> result(n + 1);

	for (int i = 0; i <= n; i++) {
		const int j = n - i; // j is the power of y
		// Create a polynomial in x which is a monomial of degree i.
		// Represented as a vector of length (i+1) with the coefficient for x^i equal to f[i].
		std::vector<Rational> poly(i + 1, 0);
		poly[i] = f.coefficients[i];
		// Place this polynomial in the result corresponding to y^j.
		result[j] = poly;
	}
	return result;
}

static matrix createSylvesterMatrix(const std::vector<Polynomial>& fSub, const Polynomial& g) {
	PROFILE_FUNCTION
	const int m = static_cast<int>(fSub.size()) - 1;
	const int n = g.degree;
	const int matrixSize = m + n;
	matrix sylvesterMatrix;

	for (int i = 0; i < n; i++) {
		std::vector<Polynomial> row(matrixSize, Polynomial({0}));
		for (int j = 0; j <= m; j++) {
			if (i + j < matrixSize) {
				row[i + j] = fSub[j];
			}
		}
		sylvesterMatrix.push_back(row);
	}

	for (int i = 0; i < m; i++) {
		std::vector<Polynomial> row(matrixSize, Polynomial({0}));
		for (int j = 0; j <= n; j++) {
			if (i + j < matrixSize) {
				row[i + j] = (j <= g.degree) ? Polynomial({g.coefficients[j]}) : Polynomial({0});
			}
		}
		sylvesterMatrix.push_back(row);
	}

	return sylvesterMatrix;
}

RealAlgebraicNumber::RealAlgebraicNumber()
	: polynomial({0}), interval({.lower_bound = 0.0, .upper_bound = 0.0}) {}

RealAlgebraicNumber::RealAlgebraicNumber(const Polynomial& polynomial, const Interval& interval)
//: polynomial(std::move(polynomial)), interval(std::move(interval))
{
	this->polynomial = polynomial;
	this->interval = interval;
	//normalize();
	while ((this->interval.lower_bound - this->interval.upper_bound).abs() > 0.01) {
		refine();
	}
	this->polynomial.normalize(this->interval.lower_bound, this->interval.upper_bound);
}

RealAlgebraicNumber::RealAlgebraicNumber(const Polynomial& polynomial, const Rational& lowerBound,
                                         const Rational& upperBound)
//: polynomial(std::move(polynomial)), interval{.lower_bound = lowerBound, .upper_bound = upperBound}
{
	this->polynomial = polynomial;
	this->interval.lower_bound = lowerBound;
	this->interval.upper_bound = upperBound;
	//normalize();
	while ((this->interval.lower_bound - this->interval.upper_bound).abs() > 0.01) {
		refine();
	}
	this->polynomial.normalize(this->interval.lower_bound, this->interval.upper_bound);
}

RealAlgebraicNumber::RealAlgebraicNumber(const std::vector<Rational>& coefficients, const Rational& lowerBound,
                                         const Rational& upperBound)
	: polynomial(coefficients), interval{.lower_bound = lowerBound, .upper_bound = upperBound} {
	//normalize();
	while ((this->interval.lower_bound - this->interval.upper_bound).abs() > 0.01) {
		refine();
	}
	this->polynomial.normalize(this->interval.lower_bound, this->interval.upper_bound);
}

RealAlgebraicNumber::RealAlgebraicNumber(const int value) {
	this->polynomial = Polynomial({-value, 1});
	this->interval.lower_bound = value - 0.01;
	this->interval.upper_bound = value + 0.01;
	this->polynomial.normalize(this->interval.lower_bound, this->interval.upper_bound);
}

RealAlgebraicNumber::RealAlgebraicNumber(const double value) {
	const Rational rationalValue = value;
	this->polynomial = Polynomial({-rationalValue, 1});
	this->interval.lower_bound = rationalValue - 0.01;
	this->interval.upper_bound = rationalValue + 0.01;
	this->polynomial.normalize(this->interval.lower_bound, this->interval.upper_bound);
}

//RealAlgebraicNumber RealAlgebraicNumber::fromInteger(const int n)
//{
//	this->polynomial = Polynomial({ -n, 1 });
//	this->interval.lower_bound = n-0.01;
//	this->interval.upper_bound = n+0.01;
//	this->polynomial.normalize(this->interval.lower_bound, this->interval.upper_bound);
//	return *this;
//}

static void printSylvester(const matrix& sylvester) {
	PROFILE_FUNCTION
	for (const std::vector<Polynomial>& i : sylvester) {
		for (const Polynomial& j : i) {
			std::cout << j.toString() << " | ";
		}
		std::cout << '\n';
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
	return {polynomial.reflectY(), -interval.upper_bound, -interval.lower_bound};
}

static Rational minRational(const Rational& r1, const Rational& r2, const Rational& r3, const Rational& r4) {
	PROFILE_FUNCTION
	Rational min = r1;
	if (r2 < min) min = r2;
	if (r3 < min) min = r3;
	if (r4 < min) min = r4;
	return min;
}

static Rational maxRational(const Rational& r1, const Rational& r2, const Rational& r3, const Rational& r4) {
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

	const std::vector<Polynomial> fXminY = poly_x_minus_y(thisCopy.polynomial);
	const matrix sylvester = createSylvesterMatrix(fXminY, otherCopy.polynomial);
	//printSylvester(sylvester);
	//Polynomial f3 = determinant(sylvester);
	MatrixXp eigenMat = convertToEigen(sylvester);
	Polynomial f3 = determinantEigen(eigenMat);

	//Polynomial f3test = determinantBareiss(sylvester);
	//std::cout << (f3 == f3test) << std::endl;

	std::vector<Polynomial> sturm;
	if (f3.sturm_sequence.empty()) {
		sturm = f3.sturmSequence(f3);
	}
	else {
		sturm = f3.sturm_sequence;
	}

	// Compute initial interval bounds
	Rational l3 = thisCopy.interval.lower_bound + otherCopy.interval.lower_bound;
	Rational r3 = thisCopy.interval.upper_bound + otherCopy.interval.upper_bound;

	// Refine interval until exactly one root is isolated
	while (variationCount(sturm, l3) - variationCount(sturm, r3) > 1) {
		thisCopy.refine();
		otherCopy.refine();
		l3 = thisCopy.interval.lower_bound + otherCopy.interval.lower_bound;
		r3 = thisCopy.interval.upper_bound + otherCopy.interval.upper_bound;
	}

	f3.normalize(l3, r3);

	const Interval sumInterval = {.lower_bound = l3, .upper_bound = r3};

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
	otherCopy.polynomial.normalize(otherCopy.interval.lower_bound, otherCopy.interval.upper_bound);
	thisCopy.polynomial.normalize(thisCopy.interval.lower_bound, thisCopy.interval.upper_bound);

	//const std::vector<Polynomial> fXoverY = poly_x_over_y(thisCopy.polynomial);
	//const matrix sylvester = createSylvesterMatrix(fXoverY, otherCopy.polynomial);
	////printSylvester(sylvester);
	////Polynomial f3 = determinant(sylvester);
	////MatrixXp eigenMat = convertToEigen(sylvester);
	////Polynomial f3 = determinantEigen(eigenMat);

	//Polynomial f3 = determinantBareiss(sylvester);

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
		thisCopy.interval.lower_bound * otherCopy.interval.lower_bound,
		thisCopy.interval.lower_bound * otherCopy.interval.upper_bound,
		thisCopy.interval.upper_bound * otherCopy.interval.lower_bound,
		thisCopy.interval.upper_bound * otherCopy.interval.upper_bound);
	Rational r3 = maxRational(
		thisCopy.interval.lower_bound * otherCopy.interval.lower_bound,
		thisCopy.interval.lower_bound * otherCopy.interval.upper_bound,
		thisCopy.interval.upper_bound * otherCopy.interval.lower_bound,
		thisCopy.interval.upper_bound * otherCopy.interval.upper_bound);

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
		//auto f1 = RealAlgebraicNumber(polynomial, interval.lower_bound, r3);
		thisCopy.refine();
		otherCopy.refine();
		l3 = minRational(
			thisCopy.interval.lower_bound * otherCopy.interval.lower_bound,
			thisCopy.interval.lower_bound * otherCopy.interval.upper_bound,
			thisCopy.interval.upper_bound * otherCopy.interval.lower_bound,
			thisCopy.interval.upper_bound * otherCopy.interval.upper_bound);
		r3 = maxRational(
			thisCopy.interval.lower_bound * otherCopy.interval.lower_bound,
			thisCopy.interval.lower_bound * otherCopy.interval.upper_bound,
			thisCopy.interval.upper_bound * otherCopy.interval.lower_bound,
			thisCopy.interval.upper_bound * otherCopy.interval.upper_bound);
	}

	f3.normalize(l3, r3);

	this->polynomial = f3;
	this->interval = {.lower_bound = l3, .upper_bound = r3};

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

	if (thisCopy.polynomial != otherCopy.polynomial) {
		return false;
	}

	for (int i = 0; i < 50; i++) {
		if (thisCopy.interval.lower_bound == otherCopy.interval.lower_bound && thisCopy.interval.upper_bound ==
			otherCopy.interval.upper_bound) {
			//std::cout << "Equal after " << i << '\n';
			return true;
		}
		if (thisCopy.interval.lower_bound > otherCopy.interval.lower_bound && thisCopy.interval.upper_bound < otherCopy.
			interval.upper_bound) {
			return true;
		}
		if (thisCopy.interval.lower_bound < otherCopy.interval.lower_bound && thisCopy.interval.upper_bound > otherCopy.
			interval.upper_bound) {
			return true;
		}

		if (thisCopy.interval.lower_bound > otherCopy.interval.upper_bound || thisCopy.interval.upper_bound < otherCopy.
			interval.lower_bound) {
			return false;
		}
		thisCopy.refine();
		otherCopy.refine();
	}
	//std::cout << "Equal after 50" << '\n';
	return true;
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
		if (thisCopy.interval.upper_bound < otherCopy.interval.lower_bound) {
			return true;
		}
		if (thisCopy.interval.lower_bound > otherCopy.interval.upper_bound) {
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
	Rational inverseLower = interval.upper_bound.inverse();
	Rational inverseUpper = interval.lower_bound.inverse();
	return {inverseCo, inverseLower, inverseUpper};
}

RealAlgebraicNumber RealAlgebraicNumber::sqrt(const int n) {
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
	if (n % 2 == 0 && interval.lower_bound < 0) {
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

	if (n % 2 != 0 && interval.lower_bound < 0) {
		l3 = -interval.upper_bound.abs().sqrt(n) - 0.01;
		r3 = -interval.lower_bound.abs().sqrt(n) + 0.01;
	}
	else {
		l3 = interval.lower_bound.sqrt(n) - 0.01;
		r3 = interval.upper_bound.sqrt(n) + 0.01;
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

	return {f3, {.lower_bound = l3, .upper_bound = r3}};
}

RealAlgebraicNumber RealAlgebraicNumber::pow(int n) const {
	PROFILE_FUNCTION
	//if (n == 0) {
	//	std::vector<Rational> coefficients = {-1, 1};
	//	return {coefficients, 0.99, 1.01};
	//}
	////if (n == 1) return *this;
	//RealAlgebraicNumber result = *this;
	//for (int i = 1; i < n; i++) {
	//	//RealAlgebraicNumber temp = *this;
	//	result *= *this;
	//}
	//return result;

	if (n == 0) {
		std::vector<Rational> coefficients = {-1, 1};
		return {coefficients, 0.99, 1.01}; // Represents the number 1
	}

	RealAlgebraicNumber thisCopy = *this;
	thisCopy.polynomial.normalize(thisCopy.interval.lower_bound, thisCopy.interval.upper_bound);

	MyMatrix<Polynomial> sylvester_mat = constructSylvesterMatrixForPower(thisCopy.polynomial, n);
	Polynomial f3 = sylvester_mat.determinant();

	Rational l3;
	Rational r3;

	// Handle k > 0
	if (this->interval.lower_bound >= 0) {
		// Both endpoints are non-negative
		l3 = this->interval.lower_bound.pow(n);
		r3 = this->interval.upper_bound.pow(n);
	}
	else if (this->interval.upper_bound <= Rational(0)) {
		// Both endpoints are non-positive
		if (n % 2 == 0) {
			// Even power: result is positive, order flips
			l3 = this->interval.upper_bound.pow(n);
			r3 = this->interval.lower_bound.pow(n);
		}
		else {
			// Odd power: result stays negative, order same
			l3 = this->interval.lower_bound.pow(n);
			r3 = this->interval.upper_bound.pow(n);
		}
	}
	else {
		// Interval spans zero (lower < 0 < upper)
		auto min_val = Rational(0);
		Rational lowerPowerN = this->interval.lower_bound.pow(n);
		Rational upperPowerN = this->interval.upper_bound.pow(n);
		Rational max_val = lowerPowerN > upperPowerN ? lowerPowerN : upperPowerN;

		// If k is even, the minimum value will be 0.
		// If k is odd, the minimum value will be lower.power(k).
		if (n % 2 != 0) {
			// Even power, like [-2, 3]^2 -> [0, 9]
			min_val = lowerPowerN;
			max_val = upperPowerN;
		}

		l3 = min_val;
		r3 = max_val;
	}

	f3.normalize(l3, r3);
	return {f3, {.lower_bound = l3, .upper_bound = r3}};
}

void RealAlgebraicNumber::testOperators() {
	PROFILE_FUNCTION
	RealAlgebraicNumber a = 2;
	RealAlgebraicNumber b = 3;
	auto c = RealAlgebraicNumber({-2, 0, 1}, {.lower_bound = {1, 1}, .upper_bound = {2, 1}}); // sqrt(2)
	auto d = RealAlgebraicNumber({-3, 0, 1}, {.lower_bound = {1, 1}, .upper_bound = {2, 1}}); // sqrt(3)

	RealAlgebraicNumber e = 5;
	RealAlgebraicNumber f = -1;
	auto g = RealAlgebraicNumber({1, 0, -10, 0, 1}, {.lower_bound = {3146, 1000}, .upper_bound = {3147, 1000}});
	// x^4 - 10x^2 + 1
	auto h = RealAlgebraicNumber({1, 0, -10, 0, 1}, {.lower_bound = {-318, 1000}, .upper_bound = {-317, 1000}});
	RealAlgebraicNumber i = 6;
	auto j = RealAlgebraicNumber({2, -3}, {.lower_bound = {666, 1000}, .upper_bound = {667, 1000}});
	auto k = RealAlgebraicNumber({4, 0, -12, 0, 9}, {.lower_bound = {816, 1000}, .upper_bound = {817, 1000}});
	// 9x^4 - 12x^2 + 4
	auto l = RealAlgebraicNumber({36, 0, -12, 0, 1}, {.lower_bound = {244, 100}, .upper_bound = {245, 100}});
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

	auto aAlt = RealAlgebraicNumber({-2, 1}, {.lower_bound = {1, 1}, .upper_bound = {4, 1}});
	auto cAlt = RealAlgebraicNumber({-2, 0, 1}, {.lower_bound = {1414, 1000}, .upper_bound = {1415, 1000}});

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

	auto invA = RealAlgebraicNumber({-1, 2}, {.lower_bound = {3, 100}, .upper_bound = {6, 100}});
	auto invC = RealAlgebraicNumber({-1, 0, 2}, {.lower_bound = {707, 1000}, .upper_bound = {708, 1000}});

	std::cout << a.inverse().toString() << std::endl;
	std::cout << c.inverse().toString() << std::endl;

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

void RealAlgebraicNumber::refine() {
	PROFILE_FUNCTION
	std::vector<Polynomial> sturm;
	if (polynomial.sturm_sequence.empty()) {
		sturm = polynomial.sturmSequence(polynomial);
	}
	else {
		sturm = polynomial.sturm_sequence;
	}

	const Rational m = (interval.lower_bound + interval.upper_bound) / 2;

	const int varL = variationCount(sturm, interval.lower_bound);
	const int varM = variationCount(sturm, m);
	if (varL > varM) {
		interval.upper_bound = m;
	}
	else {
		interval.lower_bound = m;
	}
}

std::string RealAlgebraicNumber::toString() const {
	PROFILE_FUNCTION
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

std::string RealAlgebraicNumber::toDecimalString(int precision) const {
	PROFILE_FUNCTION
	if (this->polynomial.degree == 1) {
		Rational zeroth_coeff = this->polynomial.coefficients[0];
		return (zeroth_coeff * -1).toDecimalString(0);
	}

	RealAlgebraicNumber temp = *this;

	while (true) {
		std::string lower_string = temp.interval.lower_bound.toDecimalString(precision);
		std::string upper_string = temp.interval.upper_bound.toDecimalString(precision);

		// only keeping part after '.'
		std::string lower_string_decimal = lower_string.substr(lower_string.find('.') + 1);
		std::string upper_string_decimal = upper_string.substr(upper_string.find('.') + 1);

		int matchingNumbers = 0;
		for (size_t i = 0; i < lower_string_decimal.size(); i++) {
			if (lower_string_decimal[i] == upper_string_decimal[i]) {
				matchingNumbers++;
			}
			else {
				break;
			}
		}

		if (matchingNumbers >= precision) {
			// remove trailing zeros
			lower_string.erase(lower_string.find_last_not_of('0') + 1);
			return lower_string;
		}

		temp.refine();
	}
}
