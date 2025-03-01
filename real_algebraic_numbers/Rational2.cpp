#include "Rational2.h"
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <numeric>

Rational::Rational(const mpz_class& num, const mpz_class& den) {
    if (den == 0) throw std::invalid_argument("Denominator cannot be zero");
    mpz_class g = computeGcd(num, den);
    numerator = num / g;
    denominator = den / g;
    if (denominator < 0) {
        numerator = -numerator;
        denominator = -denominator;
    }
}

Rational::Rational(int num, int den) {
    if (den == 0) throw std::invalid_argument("Denominator cannot be zero");
    int g = std::gcd(num, den);
    numerator = num / g;
    denominator = den / g;
    if (denominator < 0) {
        numerator = -numerator;
        denominator = -denominator;
    }
}

Rational::Rational(const int numerator)
{
    this->numerator = numerator;
    this->denominator = 1;
}

Rational::Rational(const mpz_class& num) : numerator(num), denominator(1) {}

Rational::Rational(double numer) {
    mpz_class denom = 1;
    while (numer != std::floor(numer)) {
        numer *= 10;
        denom *= 10;
    }
    mpz_class numInt(numer);
    mpz_class g = computeGcd(numInt, denom);
    numerator = numInt / g;
    denominator = denom / g;
}

Rational::Rational(float numer) : Rational(static_cast<double>(numer)) {}

Rational::Rational(const Rational& other) = default;

std::string Rational::toString() const {
    if (denominator == 1) return numerator.get_str();
    return numerator.get_str() + "/" + denominator.get_str();
}

void Rational::print() const {
    std::cout << numerator << "/" << denominator;
}

// Arithmetic Operators
Rational Rational::operator+(const Rational& other) const {
    mpz_class num = numerator * other.denominator + other.numerator * denominator;
    mpz_class den = denominator * other.denominator;
    return { num, den };
}

Rational Rational::operator-(const Rational& other) const {
    mpz_class num = numerator * other.denominator - other.numerator * denominator;
    mpz_class den = denominator * other.denominator;
    return { num, den };
}

Rational operator-(const mpz_class& lsh, const Rational& other) {
    return Rational(lsh) - other;
}

Rational Rational::operator-() const {
    return { -numerator, denominator };
}

Rational Rational::operator*(const Rational& other) const {
    mpz_class num = numerator * other.numerator;
    mpz_class den = denominator * other.denominator;
    return { num, den };
}

Rational Rational::operator/(const Rational& other) const {
    if (other.numerator == 0) throw std::invalid_argument("Division by zero");
    mpz_class num = numerator * other.denominator;
    mpz_class den = denominator * other.numerator;
    return { num, den };
}

Rational operator/(const mpz_class& lsh, const Rational& other) {
    return Rational(lsh) / other;
}

// Assignment Operators
Rational& Rational::operator+=(const Rational& other) {
    *this = *this + other;
    return *this;
}

Rational& Rational::operator-=(const Rational& other) {
    *this = *this - other;
    return *this;
}

Rational& Rational::operator*=(const Rational& other) {
    *this = *this * other;
    return *this;
}

Rational& Rational::operator/=(const Rational& other) {
    *this = *this / other;
    return *this;
}

// Comparison Operators
bool Rational::operator==(const Rational& other) const {
    return numerator * other.denominator == other.numerator * denominator;
}

bool Rational::operator!=(const Rational& other) const {
    return !(*this == other);
}

bool Rational::operator<(const Rational& other) const {
    return numerator * other.denominator < other.numerator * denominator;
}

bool Rational::operator>(const Rational& other) const {
    return other < *this;
}

bool Rational::operator<=(const Rational& other) const {
    return !(*this > other);
}

bool Rational::operator>=(const Rational& other) const {
    return !(*this < other);
}

// Utility Methods
Rational Rational::abs() const {
	if (numerator >= 0) return *this;
    return { -numerator, denominator };
    //return { mpz_class::abs(numerator), denominator.abs()};
}

Rational Rational::inverse() const {
    if (numerator == 0) throw std::invalid_argument("Cannot invert zero");
    return { denominator, numerator };
}

Rational Rational::sqrt(int n) const {
    if (n == 0) throw std::invalid_argument("Zeroth root is undefined");
    if (n % 2 == 0 && numerator < 0)
        throw std::invalid_argument("Even root of negative number");
    double val = mpz_get_d(numerator.get_mpz_t()) / mpz_get_d(denominator.get_mpz_t());
    return Rational(std::pow(val, 1.0 / n));
}

Rational Rational::gcd(const Rational& other) const {
    mpz_class numGcd = computeGcd(numerator, other.numerator);
    mpz_class denLcm = (denominator * other.denominator) / computeGcd(denominator, other.denominator);
    return { numGcd, denLcm };
}

mpz_class Rational::computeGcd(const mpz_class& a, const mpz_class& b) {
    mpz_class g;
    mpz_gcd(g.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    return g;
}