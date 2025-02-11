#include "Polynomial.h"
#include <iostream>

Polynomial::Polynomial()
	: degree(-1), coefficients({}), derivative({})
{}

Polynomial::Polynomial(const std::vector<int>& coefficients)
	: degree(static_cast<int>(coefficients.size()) - 1), coefficients(coefficients)
{
	computeDerivative();
	computeSturmSequence();
}

void Polynomial::computeDerivative() {
	if (degree < 1) return;
	for (size_t i = 1; i < coefficients.size(); ++i) {
		derivative.push_back(coefficients[i] * i);
	}
}

void Polynomial::computeSturmSequence() {
	sturm_sequence.push_back(coefficients);
	sturm_sequence.push_back(derivative);

	while (!sturm_sequence.back().empty()) {
		std::vector<int> rem = polyMod(sturm_sequence[sturm_sequence.size() - 2], sturm_sequence.back());
		if (rem.empty()) break;
		sturm_sequence.push_back(rem);
	}
}

//void Polynomial::generate_sturm_sequence() {
//	sturm_sequence.clear();
//	if (degree == -1) return;
//
//	// Add P0 (original polynomial) and P1 (derivative)
//	sturm_sequence.push_back(coefficients);
//	Polynomial p1 = derivative;
//	sturm_sequence.push_back(p1.coefficients);
//
//	Polynomial p0 = *this;
//	while (!p1.degree == -1 && p1.degree > 0) {
//		Polynomial p_next = -(p0 % p1);
//		sturm_sequence.push_back(p_next.coefficients);
//		p0 = p1;
//		p1 = p_next;
//	}
//}

std::vector<int> Polynomial::polyMod(const std::vector<int>& dividend, const std::vector<int>& divisor) {
	if (divisor[0] == 0) return {}; // Division by zero
	std::vector<int> rem = dividend;
	int degDiv = static_cast<int>(divisor.size()) - 1;

	while (rem.size() >= divisor.size()) {
		const int coef = rem[0] / divisor[0];
		for (size_t i = 0; i < divisor.size(); i++) {
			rem[i] -= coef * divisor[i];
		}
		while (!rem.empty() && rem[0] == 0) {
			rem.erase(rem.begin());
		}
	}

	// Negate remainder for Sturm sequence property
	for (int& coeff : rem) {
		coeff = -coeff;
	}

	return rem;
}

std::vector<int> Polynomial::polyAdd(const std::vector<int>& p1, const std::vector<int>& p2) const {
	std::vector<int> result(std::max(p1.size(), p2.size()), 0);
	for (size_t i = 0; i < p1.size(); i++) result[i] += p1[i];
	for (size_t i = 0; i < p2.size(); i++) result[i] += p2[i];
	return result;
}

std::vector<int> Polynomial::polyMultiply(const std::vector<int>& p1, const std::vector<int>& p2) const {
	std::vector<int> result(p1.size() + p2.size() - 1, 0);
	for (size_t i = 0; i < p1.size(); i++) {
		for (size_t j = 0; j < p2.size(); j++) {
			result[i + j] += p1[i] * p2[j];
		}
	}
	return result;
}

Eigen::MatrixXd Polynomial::buildSylvesterMatrix(const Polynomial& a, const Polynomial& b) const {
	const std::vector<int>& coA = a.coefficients;
	const std::vector<int>& coB = b.coefficients;
	const int m = a.degree; // Degree of A
	const int n = b.degree; // Degree of B
	const int size = m + n;

	Eigen::MatrixXd s = Eigen::MatrixXd::Zero(size, size);

	// Fill the first n rows with shifted coefficients of A
	for (int row = 0; row < n; row++) {
		for (size_t col = 0; col < coA.size(); col++) {
			s(row, row + col) = coA[col]; // Shifted placement
		}
	}

	// Fill the last m rows with shifted coefficients of B
	for (int row = 0; row < m; row++) {
		for (size_t col = 0; col < coB.size(); col++) {
			s(n + row, row + col) = coB[col]; // Shifted placement
		}
	}

	return s;
}

std::vector<int> Polynomial::determinantAsPolynomial(Eigen::MatrixXd& s) const {
	const int n = s.rows();
	if (n == 1) return { static_cast<int>(s(0, 0)) }; // Base case: single value matrix

	std::vector<int> resultantPoly = { 0 }; // Initialize resultant polynomial

	// Expansion along the first row (Laplace Expansion)
	for (int col = 0; col < n; col++) {
		Eigen::MatrixXd subMatrix(n - 1, n - 1);

		// Build submatrix by excluding row 0 and column col
		for (int i = 1; i < n; i++) {
			int subCol = 0;
			for (int j = 0; j < n; j++) {
				if (j == col) continue;
				subMatrix(i - 1, subCol++) = s(i, j);
			}
		}

		// Recursively compute determinant of the submatrix
		std::vector<int> subDet = determinantAsPolynomial(subMatrix);

		// Multiply by (-1)^col * coefficient from first row
		if (static_cast<int>(s(0, col)) != 0) {
			std::vector<int> term = polyMultiply({ static_cast<int>(s(0, col)) }, subDet);
			if (col % 2 == 1) {
				for (int& coef : term) coef = -coef; // Apply sign
			}
			resultantPoly = polyAdd(resultantPoly, term);
		}
	}

	return resultantPoly;
}

Polynomial Polynomial::computeResultantPolynomial(const Polynomial& other) const {
	Eigen::MatrixXd s = buildSylvesterMatrix(coefficients, other.coefficients);
	std::cout << s << "\n";
	return determinantAsPolynomial(s);
}

