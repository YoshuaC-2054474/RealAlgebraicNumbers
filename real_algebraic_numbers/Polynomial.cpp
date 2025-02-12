#include "Polynomial.h"
#include <iostream>

Polynomial::Polynomial(const std::initializer_list<int> coeffs)
{
    //int exponent = 0;
    for (auto it = coeffs.begin(); it < coeffs.end(); ++it) {
        coefficients.push_back(*it);
    }

    degree = coefficients.empty() ? -1 : coefficients.size() - 1;
}

Polynomial::Polynomial(const std::vector<Rational>& coeffs)
{
	coefficients = coeffs;
	while (!coefficients.empty() && coefficients.back() == 0) {
		coefficients.pop_back();
	}

	degree = coefficients.empty() ? -1 : coefficients.size() - 1;
}

bool Polynomial::isZero() const {
    return degree == -1;
}

Polynomial Polynomial::operator+(const Polynomial& other) {
    std::vector<Rational> result = coefficients;
    for (int i = 0; i < other.coefficients.size(); i++)
    {
        while (result.size() <= i) result.push_back(0);
        result[i] += other.coefficients[i];
    }
    return result;
}

Polynomial Polynomial::operator-(const Polynomial& other) {
    std::vector<Rational> result = coefficients;
    for (int i = 0; i < other.coefficients.size(); i++)
    {
        while (result.size() <= i) result.push_back(0);
        result[i] -= other.coefficients[i];
    }
    return result;
}

Polynomial Polynomial::operator*(const Polynomial& other) const {
    Polynomial result;
    for (int i = 0; i < coefficients.size(); i++)
    {
        for (int j = 0; j < other.coefficients.size(); j++)
        {
            while (result.coefficients.size() <= i + j) result.coefficients.push_back(0);
            result.coefficients[i + j] += coefficients[i] * other.coefficients[j];
        }
    }

    result.degree = result.coefficients.size() - 1;
    return result;
}

std::string Polynomial::toString() const {
    if (coefficients.empty()) return "0";
    std::string output;
    for (int i = coefficients.size() - 1; i >= 0; i--) {
        if (coefficients[i] == 0) continue;
        //if (coefficients[i] != 1 && coefficients[i] != -1)
        if (coefficients[i] > 0 && !output.empty()) output += " + ";
        else if (coefficients[i] < 0) output += " - ";
        if (coefficients[i].abs() != 1.0 || i == 0)
            output += coefficients[i].abs().toString();
        if (i > 1) output += "x^" + std::to_string(i);
        if (i == 1) output += "x";
    }
    return output;
}

void Polynomial::print() const {
    std::cout << toString() << "\n";
}

Polynomial Polynomial::derivative() const {
    if (isZero()) return {};
    std::vector<Rational> result(degree, 0);
    for (int i = 1; i <= degree; ++i) {
        result[i - 1] = coefficients[i] * i;
    }
    return { result };
}

Polynomial Polynomial::polyTrim(const Polynomial& poly) {
    Polynomial result = poly;
    while (!result.coefficients.empty() && result.coefficients.back().abs() < 1e-9) {
        result.coefficients.pop_back();
    }
    result.degree = result.coefficients.empty() ? -1 : result.coefficients.size() - 1;
    return result;
}

std::pair<std::vector<Rational>, std::vector<Rational>> Polynomial::polyDivide(const Polynomial& dividend, const Polynomial& divisor) {
    const Polynomial A = polyTrim(dividend);
    const Polynomial B = polyTrim(divisor);

    if (B.coefficients.empty()) {
        throw std::runtime_error("Division by zero polynomial!");
    }

    const int degA = A.degree;
    const int degB = B.degree;

    std::vector<Rational> quotient(degA - degB + 1, 0);
    std::vector<Rational> remainder = A.coefficients;

    // Perform polynomial long division
    while (remainder.size() >= B.coefficients.size() && !remainder.empty()) {
        const int degR = remainder.size() - 1;
        Rational factor = remainder.back() / B.coefficients.back();
        const int power = degR - degB;
        quotient[power] = factor;

        // Subtract factor * (B shifted by "power") from remainder.
        std::vector<Rational> sub(power, 0);  // zeros for lower-degree terms
        for (const Rational& coeff : B.coefficients) {
            sub.push_back(coeff * factor);
        }

        // Ensure remainder has enough size.
        if (remainder.size() < sub.size())
            remainder.resize(sub.size(), 0.0);
        for (size_t i = 0; i < sub.size(); i++) {
            remainder[i] -= sub[i];
        }
        remainder = polyTrim(remainder).coefficients;
    }

    quotient = polyTrim(quotient).coefficients;
    remainder = polyTrim(remainder).coefficients;
    return make_pair(quotient, remainder);
}

std::vector<Rational> Polynomial::polyNegate(const std::vector<Rational>& poly) {
    std::vector<Rational> neg(poly.size());
    for (size_t i = 0; i < poly.size(); i++) {
        neg[i] = -poly[i];
    }
    return neg;
}

Polynomial Polynomial::reflectY() const
{
    std::vector<Rational> result = coefficients;
    for (int i = result.size() - 1; i >= 0; i--)
    {
        if (i % 2 == 1)
        {
            result[i] = -result[i];
        }
    }
    return { result };
}

// Generate the Sturm sequence
std::vector<Polynomial> Polynomial::sturmSequence(const Polynomial& p) {
    std::vector<Polynomial> seq;

    const Polynomial S0 = polyTrim(p);
    if (S0.coefficients.empty()) {
        throw std::runtime_error("Zero polynomial provided.");
    }
    seq.push_back(S0);

    // S1 = p'(x)
    const Polynomial S1 = S0.derivative();
    seq.push_back(S1);

    // Compute subsequent sequence members until remainder becomes zero.
    while (true) {
        // Let S_prev = S_{i-1} and S_curr = S_i.
        const Polynomial& S_prev = seq[seq.size() - 2];
        const Polynomial& S_curr = seq.back();

        // Compute polynomial remainder of S_prev divided by S_curr.
        auto divRes = polyDivide(S_prev, S_curr);
        std::vector<Rational> rem = divRes.second;

        // If remainder is zero, stop.
        if (rem.empty() || (rem.size() == 1 && rem[0].abs() < 1e-9)) {
            break;
        }

        // Next term in Sturm sequence is the negative of the remainder.
        std::vector<Rational> next = polyNegate(rem);
        seq.push_back(next);
    }

    sturm_sequence = seq;
    return seq;
}