#ifndef MATRIX_H
#define MATRIX_H

#include "Polynomial.h"
#include <vector>
#include <stdexcept>
#include "MyTimer.h"

template <typename T>
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
		std::cout << "Matrix (" << rows << "x" << cols << "):" << std::endl;
		int longest = 0;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				int len = data[i][j].toString().length();
				if (len > longest) {
					longest = len;
				}
			}
		}
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				std::cout << data[i][j];
				for (int k = 0; k < longest - data[i][j].toString().length() + 1; ++k) {
					std::cout << " ";
				}
				std::cout << "| ";
			}
			std::cout << std::endl;
		}
	}

    // Bareiss Algorithm for Determinant (Optimized)
    T determinant() const {
        //PROFILE_FUNCTION
        if (rows != cols) {
            throw std::runtime_error("Determinant can only be calculated for square matrices.");
        }

        const int n = rows;
        if (n == 0) return T(Rational(1));
        if (n == 1) return data[0][0];

        std::vector<std::vector<T>> B = data;
        T prev_pivot_poly = T(Rational(1));
        int sign = 1;

        for (int k = 0; k < n - 1; ++k) {
            if (B[k][k].isZero()) {
                int swap_row = -1;
                for (int i = k + 1; i < n; ++i) {
                    if (!B[i][k].isZero()) {
                        swap_row = i;
                        break;
                    }
                }

                if (swap_row == -1) {
                    return T(Rational(0));
                }

                std::swap(B[k], B[swap_row]);
                sign *= -1;
            }

            const T& current_pivot_poly = B[k][k];

            for (int i = k + 1; i < n; ++i) {
                for (int j = k + 1; j < n; ++j) {
                    B[i][j] = (current_pivot_poly * B[i][j] - B[i][k] * B[k][j]) / prev_pivot_poly;
                }
            }
            prev_pivot_poly = current_pivot_poly;
        }

        T result = B[n - 1][n - 1];
        if (sign == -1) {
            return -result;
        }
        return result;
    }
};

#endif // MATRIX_H