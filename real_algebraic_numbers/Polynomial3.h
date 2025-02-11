#ifndef POLYNOMIAL3_H
#define POLYNOMIAL3_H

#include "Rational.h"
#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>

using namespace std;

class Polynomial {
public:
    int degree;
    std::vector<Rational> coefficients;  // Maps exponent to coefficient
    std::vector<Polynomial> sturm_sequence;

	Polynomial() : degree(-1), coefficients({}) {}
    /*Polynomial(double coeff, int exp) {
        coefficients[exp] = coeff;
    }*/
    Polynomial(const std::initializer_list<int> coeffs)
	{
        //int exponent = 0;
        for (auto it = coeffs.begin();  it < coeffs.end(); ++it) {
            coefficients.push_back(*it);
        }

        degree = coefficients.empty() ? -1 : coefficients.size() - 1;
    }

	Polynomial(const std::vector<Rational>& coeffs)
	{
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

    //Polynomial(double constant) : Polynomial(constant, 0) {}

    //void addTerm(double coeff, int exp) {
    //    //if (coeff == 0.0) return;
    //    coefficients[exp] += coeff;
    //    //if (coefficients[exp] == 0.0) coefficients.erase(exp);
    //}

    Polynomial operator+(const Polynomial& other) {
        //Polynomial result = *this;
		std::vector<Rational> result = coefficients;
        for (int i = 0; i < other.coefficients.size(); i++)
        {
			while (result.size() <= i) result.push_back(0);
        	result[i] += other.coefficients[i];
        }
        //for (const auto& term : other.coefficients) result.addTerm(term, term.first);
        return result;
    }

    Polynomial operator-(const Polynomial& other) {
        //Polynomial result = *this;
        std::vector<Rational> result = coefficients;
        for (int i = 0; i < other.coefficients.size(); i++)
        {
            while (result.size() <= i) result.push_back(0);
            result[i] -= other.coefficients[i];
        }
        //for (const auto& term : other.coefficients) result.addTerm(-term.second, term.first);
        return result;
    }

    Polynomial operator*(const Polynomial& other) const {
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
        /*for (const auto& term1 : coefficients)
            for (const auto& term2 : other.coefficients)
                result.addTerm(term1.second * term2.second, term1.first + term2.first);*/
        return result;
    }

    // Polynomial division (returns quotient and remainder)
    //friend std::pair<Polynomial, Polynomial> divide(const Polynomial& dividend, const Polynomial& divisor) {
    //    if (divisor.isZero()) throw std::invalid_argument("Division by zero polynomial");
    //    if (dividend.degree < divisor.degree) return { Polynomial(), dividend };

    //    Polynomial remainder = dividend;
    //    std::vector<int> quotient_coeffs(dividend.degree - divisor.degree + 1, 0);

    //    while (!remainder.isZero() && remainder.degree >= divisor.degree) {
    //        int deg_diff = remainder.degree - divisor.degree;
    //        int factor = remainder.coefficients.back() / divisor.coefficients.back();

    //        quotient_coeffs[deg_diff] = factor;
    //        std::vector<int> term_coeffs(deg_diff + 1, 0);
    //        term_coeffs[deg_diff] = factor;
    //        Polynomial term(term_coeffs);

    //        remainder = remainder - term * divisor;
    //    }

    //    return { Polynomial(quotient_coeffs), remainder };
    //}

    //// Remainder operator
    //friend Polynomial operator%(const Polynomial& dividend, const Polynomial& divisor) {
    //    return divide(dividend, divisor).second;
    //}

    Polynomial& operator+=(const Polynomial& other) { *this = *this + other; return *this; }
    Polynomial& operator-=(const Polynomial& other) { *this = *this - other; return *this; }
    Polynomial& operator*=(const Polynomial& other) { *this = *this * other; return *this; }
    Polynomial operator-() const { return *this * Polynomial({-1}); }

    std::string toString() const {
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
			//output += std::to_string(i);
			
		}
        /*std::stringstream ss;
        bool first = true;
        for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it) {
            double coeff = it->second;
            int exp = it->first;
            if (!first) ss << (coeff > 0 ? " + " : " - ");
            else if (coeff < 0) ss << "-";
            if (first) first = false;
            coeff = std::abs(coeff);
            if (coeff != 1 || exp == 0) ss << coeff;
            if (exp > 0) ss << "x" << (exp > 1 ? "^" + std::to_string(exp) : "");
        }
        return ss.str();*/
        return output;
    }
		
    void print() const {
		std::cout << toString() << "\n";
	}

    // Compute the derivative
    Polynomial derivative() const {
		if (isZero()) return Polynomial();
		std::vector<Rational> result(degree, 0);
		for (int i = 1; i <= degree; ++i) {
			result[i - 1] = coefficients[i] * i;
		}
		return Polynomial(result);
    }

    Polynomial polyTrim(const Polynomial& poly) {
        Polynomial result = poly;
        while (!result.coefficients.empty() && result.coefficients.back().abs() < 1e-9) {
            result.coefficients.pop_back();
        }
		result.degree = result.coefficients.empty() ? -1 : result.coefficients.size() - 1;
        return result;
    }

    pair<vector<Rational>, vector<Rational>> polyDivide(const Polynomial& dividend, const Polynomial& divisor) {
        Polynomial A = polyTrim(dividend);
        Polynomial B = polyTrim(divisor);

        if (B.coefficients.empty()) {
            throw runtime_error("Division by zero polynomial!");
        }

        int degA = A.degree;
        int degB = B.degree;

        vector<Rational> quotient(degA - degB + 1, 0);
        vector<Rational> remainder = A.coefficients;

        // Perform polynomial long division
        while (remainder.size() >= B.coefficients.size() && !remainder.empty()) {
            int degR = remainder.size() - 1;
            Rational factor = remainder.back() / B.coefficients.back();
            int power = degR - degB;
            quotient[power] = factor;

            // Subtract factor * (B shifted by "power") from remainder.
            vector<Rational> sub(power, 0);  // zeros for lower-degree terms
            for (Rational coeff : B.coefficients) {
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

    vector<Rational> polyNegate(const vector<Rational>& poly) {
        vector<Rational> neg(poly.size());
        for (size_t i = 0; i < poly.size(); i++) {
            neg[i] = -poly[i];
        }
        return neg;
    }

    Polynomial reflectY() const
	{
		std::vector<Rational> result = coefficients;
		for (int i = result.size() -1 ; i >= 0; i--)
		{
			if (i % 2 == 1)
			{
				result[i] = -result[i];
			}
		}
		return Polynomial(result);
	}

    /*vector<double> polyDerivative(const vector<double>& poly) {
        vector<double> deriv;
        for (size_t i = 1; i < poly.size(); i++) {
            deriv.push_back(i * poly[i]);
        }
        return polyTrim(deriv);
    }*/

    // Generate the Sturm sequence
    vector<Polynomial> sturmSequence(const Polynomial& p) {
        vector<Polynomial> seq;

        // S0 = p(x)
        Polynomial S0 = polyTrim(p);
        if (S0.coefficients.empty()) {
            throw runtime_error("Zero polynomial provided.");
        }
        seq.push_back(S0);

        // S1 = p'(x)
        Polynomial S1 = S0.derivative();
        seq.push_back(S1);

        // Compute subsequent sequence members until remainder becomes zero.
        while (true) {
            // Let S_prev = S_{i-1} and S_curr = S_i.
            const Polynomial& S_prev = seq[seq.size() - 2];
            const Polynomial& S_curr = seq.back();

            // Compute polynomial remainder of S_prev divided by S_curr.
            auto divRes = polyDivide(S_prev, S_curr);
            vector<Rational> rem = divRes.second;

            // If remainder is zero, stop.
            if (rem.empty() || (rem.size() == 1 && rem[0].abs() < 1e-9)) {
                break;
            }

            // Next term in Sturm sequence is the negative of the remainder.
            vector<Rational> next = polyNegate(rem);
            seq.push_back(next);
        }

        sturm_sequence = seq;
        return seq;
    }
};

#endif