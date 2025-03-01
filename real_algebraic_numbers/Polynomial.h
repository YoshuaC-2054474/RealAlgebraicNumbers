#ifndef POLYNOMIAL3_H
#define POLYNOMIAL3_H

#include "Rational3.h"
#include <string>
#include <vector>

class Polynomial {
public:
    int degree;
    std::vector<Rational> coefficients;
    std::vector<Polynomial> sturm_sequence;

	Polynomial() : degree(-1), coefficients({}) {}
    Polynomial(const std::initializer_list<int> coeffs);
    Polynomial(const std::vector<Rational>& coeffs);

    void normalize();

    // Check if the polynomial is zero
    bool isZero() const;

    Polynomial operator+(const Polynomial& other);
    Polynomial operator-(const Polynomial& other);
    Polynomial operator*(const Polynomial& other) const;

    Polynomial& operator+=(const Polynomial& other) { *this = *this + other; return *this; }
    Polynomial& operator-=(const Polynomial& other) { *this = *this - other; return *this; }
    Polynomial& operator*=(const Polynomial& other) { *this = *this * other; return *this; }
    Polynomial operator-() const { return *this * Polynomial({-1}); }

    std::string toString() const;
    void print() const;

    // Compute the derivative
    Polynomial derivative() const;

    Polynomial polyTrim(const Polynomial& poly);
    std::pair<std::vector<Rational>, std::vector<Rational>> polyDivide(const Polynomial& dividend, const Polynomial& divisor);
    std::vector<Rational> polyNegate(const std::vector<Rational>& poly);

    Polynomial reflectY() const;

    std::vector<Polynomial> sturmSequence(const Polynomial& p);
};

#endif