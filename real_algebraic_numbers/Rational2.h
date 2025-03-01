#ifndef RATIONAL2_H
#define RATIONAL2_H

#include <gmpxx.h>
#include <iostream>
#include <string>

class Rational {
public:
    Rational() : numerator(0), denominator(1) {}
    Rational(const mpz_class& num, const mpz_class& den);
    Rational(int num, int den);
    Rational(const mpz_class& numerator);
    Rational(const int numerator);
    Rational(double numer);
    Rational(float numer);
    Rational(const Rational& other);

    std::string toString() const;
    void print() const;

    /*friend std::ostream& operator<<(std::ostream& os, const Rational& rational) {
        os << rational.numerator << "/" << rational.denominator;
        return os;
    }*/

    explicit operator double() const {
        return mpz_get_d(numerator.get_mpz_t()) / mpz_get_d(denominator.get_mpz_t());
    }

    Rational operator+(const Rational& other) const;
    Rational operator-(const Rational& other) const;
    friend Rational operator-(const mpz_class& lsh, const Rational& other);
    Rational operator-() const;
    Rational operator*(const Rational& other) const;
    Rational operator/(const Rational& other) const;
    friend Rational operator/(const mpz_class& lsh, const Rational& other);

    Rational& operator+=(const Rational& other);
    Rational& operator-=(const Rational& other);
    Rational& operator*=(const Rational& other);
    Rational& operator/=(const Rational& other);

    bool operator==(const Rational& other) const;
    bool operator!=(const Rational& other) const;
    bool operator<(const Rational& other) const;
    bool operator>(const Rational& other) const;
    bool operator<=(const Rational& other) const;
    bool operator>=(const Rational& other) const;

    Rational abs() const;
    Rational inverse() const;
    Rational sqrt(int n) const;
    Rational gcd(const Rational& other) const;

private:
    mpz_class numerator;
    mpz_class denominator;

    static mpz_class computeGcd(const mpz_class& a, const mpz_class& b);
};

#endif