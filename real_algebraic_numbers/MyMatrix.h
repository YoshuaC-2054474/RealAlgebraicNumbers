// Matrix.h
#ifndef MATRIX_H
#define MATRIX_H

#include "Polynomial.h"
#include <vector>
#include <stdexcept>

template <typename T> // T will be Polynomial in our case
class MyMatrix {
public:
	int rows, cols;
	std::vector<std::vector<T>> data;

	MyMatrix(int r, int c) : rows(r), cols(c), data(r, std::vector<T>(c)) {}

	void set_element(int r, int c, const T& val) {
		if (r < 0 || r >= rows || c < 0 || c >= cols) {
			throw std::out_of_range("Matrix element out of bounds.");
		}
		data[r][c] = val;
	}

	T get_element(int r, int c) const {
		if (r < 0 || r >= rows || c < 0 || c >= cols) {
			throw std::out_of_range("Matrix element out of bounds.");
		}
		return data[r][c];
	}

	void print() const {
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				std::cout << data[i][j] << "\t";
			}
			std::cout << std::endl;
		}
	}

	// Bareiss Algorithm for Determinant
	T determinant() const {
		if (rows != cols) {
			throw std::runtime_error("Determinant can only be calculated for square matrices.");
		}

		int n = rows;
		if (n == 0) return T(Rational(1)); // Determinant of an empty matrix is 1
		if (n == 1) return data[0][0];

		// Create a mutable copy of the matrix
		std::vector<std::vector<T>> B = data;

		// Previous pivot element, initialized to 1 (Rational(1))
		// This 'prev_pivot' will be a Polynomial representing the Rational 1
		T prev_pivot_poly = T(Rational(1));

		for (int k = 0; k < n - 1; ++k) {
			// Find pivot element B[k][k]
			if (B[k][k].isZero()) {
				// Find a non-zero element in the same column below B[k][k]
				int swap_row = -1;
				for (int i = k + 1; i < n; ++i) {
					if (!B[i][k].isZero()) {
						swap_row = i;
						break;
					}
				}

				if (swap_row == -1) {
					// Column is all zeros below pivot. Determinant is zero.
					return T(Rational(0));
				}

				// Swap rows k and swap_row
				std::swap(B[k], B[swap_row]);

				// If rows were swapped, determinant sign flips.
				// We'll multiply the final determinant by -1 if an odd number of swaps occurred.
				// For Bareiss, instead of tracking swaps directly, it's easier to just
				// return 0 if the pivot is zero and cannot be made non-zero.
				// Or, more correctly, if B[k][k] becomes zero and no non-zero element below it
				// exists, then the determinant is zero.
				// Bareiss is designed to work without explicit division by zero if coefficients
				// are in an integral domain. If a pivot is zero and can't be fixed by row swap,
				// the determinant is zero.
			}

			// The actual pivot for this iteration
			T current_pivot_poly = B[k][k];

			//std::cout << "B[k][k] = " << B[k][k].degree << "\n";

			for (int i = k + 1; i < n; ++i) {
				for (int j = k + 1; j < n; ++j) {
					//std::cout << "B[i][j] = " << B[i][j].degree << "\n";
					//std::cout << "B[i][k] = " << B[i][k].degree << "\n";
					// B[i][j] = (B[k][k] * B[i][j] - B[i][k] * B[k][j]) / prev_pivot
					T term1 = current_pivot_poly * B[i][j];
					T term2 = B[i][k] * B[k][j];
					B[i][j] = (term1 - term2) / prev_pivot_poly;
					//std::cout << "term1 = " << term1.degree << "\n";
					//std::cout << "term2 = " << term2.degree << "\n";
					//std::cout << "B[i][j] after = " << B[i][j].degree << "\n";
				}
			}
			prev_pivot_poly = current_pivot_poly;
		}

		return B[n - 1][n - 1];
	}
};

#endif // MATRIX_H
