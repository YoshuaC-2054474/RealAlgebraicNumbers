#ifndef POLYNOMIAL2_H
#define POLYNOMIAL2_H

#include <vector>
#include <stdexcept>
#include <iostream>

class Polynomial {
public:
    int degree;
    std::vector<int> coefficients;
    std::vector<std::vector<int>> sturm_sequence;

    Polynomial() : degree(-1), coefficients({}) {}
    Polynomial(const std::vector<int>& coeffs) {
        coefficients = coeffs;
        while (!coefficients.empty() && coefficients.back() == 0) {
            coefficients.pop_back();
        }
        degree = coefficients.empty() ? -1 : coefficients.size() - 1;
    }

    // Check if the polynomial is zero
    bool isZero() const {
        return degree == -1;
    }

    // Arithmetic operations
    Polynomial operator-() const {
        std::vector<int> result;
        for (int c : coefficients) result.push_back(-c);
        return Polynomial(result);
    }

    Polynomial operator+(const Polynomial& other) const {
        std::vector<int> result(std::max(coefficients.size(), other.coefficients.size()), 0);
        for (size_t i = 0; i < coefficients.size(); ++i) result[i] += coefficients[i];
        for (size_t i = 0; i < other.coefficients.size(); ++i) result[i] += other.coefficients[i];
        return Polynomial(result);
    }

    Polynomial operator-(const Polynomial& other) const {
        return *this + (-other);
    }

    Polynomial operator*(const Polynomial& other) const {
        if (isZero() || other.isZero()) return Polynomial();
        std::vector<int> result(degree + other.degree + 1, 0);
        for (int i = 0; i <= degree; ++i) {
            for (int j = 0; j <= other.degree; ++j) {
                result[i + j] += coefficients[i] * other.coefficients[j];
            }
        }
        return Polynomial(result);
    }

    // Compute the derivative
    Polynomial derivative() const {
        if (degree < 1) return Polynomial();
        std::vector<int> deriv_coeffs;
        for (size_t i = 1; i < coefficients.size(); ++i) {
            deriv_coeffs.push_back(coefficients[i] * i);
        }
        return Polynomial(deriv_coeffs);
    }

    // Polynomial division (returns quotient and remainder)
    friend std::pair<Polynomial, Polynomial> divide(const Polynomial& dividend, const Polynomial& divisor) {
        if (divisor.isZero()) throw std::invalid_argument("Division by zero polynomial");
        if (dividend.degree < divisor.degree) return { Polynomial(), dividend };

        Polynomial remainder = dividend;
        std::vector<int> quotient_coeffs(dividend.degree - divisor.degree + 1, 0);

        while (!remainder.isZero() && remainder.degree >= divisor.degree) {
            int deg_diff = remainder.degree - divisor.degree;
            int factor = remainder.coefficients.back() / divisor.coefficients.back();

            quotient_coeffs[deg_diff] = factor;
            std::vector<int> term_coeffs(deg_diff + 1, 0);
            term_coeffs[deg_diff] = factor;
            Polynomial term(term_coeffs);

            remainder = remainder - term * divisor;
        }

        return { Polynomial(quotient_coeffs), remainder };
    }

    // Remainder operator
    friend Polynomial operator%(const Polynomial& dividend, const Polynomial& divisor) {
        return divide(dividend, divisor).second;
    }

    // Generate the Sturm sequence
    void generate_sturm_sequence() {
        sturm_sequence.clear();
        if (isZero()) return;

        // Add P0 (original polynomial) and P1 (derivative)
        sturm_sequence.push_back(coefficients);
        Polynomial p1 = derivative();
        sturm_sequence.push_back(p1.coefficients);

        Polynomial p0 = *this;
        while (!p1.isZero() && p1.degree > 0) {
            Polynomial p_next = -(p0 % p1);
            sturm_sequence.push_back(p_next.coefficients);
            p0 = p1;
            p1 = p_next;
        }
    }
};

#endif