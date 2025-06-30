#include "Polynomial.h"
#include <iostream>
#include <numeric>
#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <NTL/ZZX.h>
#include <NTL/ZZXFactoring.h>
using namespace NTL;
#include "MyTimer.h"

Polynomial::Polynomial(const std::initializer_list<int> coeffs) {
	for (auto it = coeffs.begin(); it < coeffs.end(); ++it) {
		coefficients.emplace_back(*it);
	}

	degree = coefficients.empty() ? -1 : static_cast<int>(coefficients.size()) - 1;
}

Polynomial::Polynomial(const std::vector<Rational>& coeffs) {
	coefficients = coeffs;

	degree = coefficients.empty() ? -1 : static_cast<int>(coefficients.size()) - 1;
}

Polynomial::Polynomial(const int zeroCoeff) {
	const std::vector<Rational> coeffs = {zeroCoeff};
	coefficients = coeffs;

	degree = coefficients.empty() ? -1 : static_cast<int>(coefficients.size()) - 1;
}

Rational findGcd(const std::vector<Rational>& arr) {
	PROFILE_FUNCTION
	Rational res = arr[0];

	for (size_t i = 1; i < arr.size(); i++) {
		res = arr[i].gcd(res);
		if (res == 1)
			return 1;
	}

	return res;
}

std::vector<Polynomial> minimalPolynomialsNtl(const Polynomial& poly) {
	PROFILE_FUNCTION
	ZZX f;
	f.SetMaxLength(static_cast<long>(poly.coefficients.size()));

	for (long i = 0; i < static_cast<long>(poly.coefficients.size()); ++i) {
		SetCoeff(f, i, to_ZZ(poly.coefficients[i].toInt()));
	}

	Vec<Pair<ZZX, long>> factors;
	ZZ content;
	factor(content, factors, f);

	// convert 'factors' to vector of Polynomials
	std::vector<Polynomial> result;
	for (long i = 0; i < factors.length(); ++i) {
		ZZX factor = factors[i].a;
		std::vector<Rational> coeffs;
		coeffs.reserve(factor.rep.length());
		for (long j = 0; j < factor.rep.length(); ++j) {
			coeffs.emplace_back(to_int(factor.rep[j]));
		}
		result.emplace_back(coeffs);
	}

	return result;
}

void Polynomial::normalize(const Rational& lowerBound, const Rational& upperBound) {
	PROFILE_FUNCTION
	if (is_normalized) return;

	while (coefficients.size() > 1 && coefficients.back() == 0) {
		coefficients.pop_back();
	}
	degree = coefficients.empty() ? -1 : static_cast<int>(coefficients.size()) - 1;

	if (coefficients.empty()) return;

	if (const Rational gcd = findGcd(coefficients); gcd > 1) {
		for (size_t i = 0; i < coefficients.size(); i++) {
			coefficients[i] /= gcd;
		}
	}

	//const Rational highestCo = coefficients.back();
	//int highestCoInt = highestCo.toInts();
	////std::cout << "\n";
	//if (highestCo != 0)
	//{
	//    for (int i = 0; i < coefficients.size(); i++) {
	//        coefficients[i] /= highestCo;
	//    }
	//}

	std::vector<Polynomial> minimalPolys = minimalPolynomialsNtl(*this);

	// minimalPolys = get_minimal_polynomials(*this);
	//std::cout << "\n\nMinimal Polynomial of " << toString() << "\n";
	for (const Polynomial& poly : minimalPolys) {
		//std::cout << "\t" << poly.toString() << '\n';
		//std::cout << "\t" << poly.toString() << " -> " << poly.polyTrim(poly).toString() << '\n';
		const bool firstSign = poly.evaluate(lowerBound) < 0;
		const bool secondSign = poly.evaluate(upperBound) < 0;
		if (firstSign != secondSign) {
			// poly has root in interval
			*this = poly;
			break;
		}
	}
	//std::cout << "= " << toString() << "\n\n";

	// make first coefficient positive
	if (coefficients[coefficients.size() - 1] < 0) {
		for (int i = 0; i < coefficients.size(); i++) {
			coefficients[i] *= -1;
		}
	}

	//std::cout << "Normalized Polynomial: " << toString() << "\n";
	//std::cout << "Degree: " << degree << "\n";
	//std::cout << "Coefficients: ";
	//for (const auto& coeff : coefficients) {
	//    std::cout << coeff.toInts() << " ";
	//}
	//std::cout << "\n";

	//std::cout << toString() << "\n";

	is_normalized = true;
	//return minimalPolys;
}

//void Polynomial::testNormalize()
//{
//	Polynomial p1({ -4,0,1 });
//	std::vector<Polynomial> n1 = p1.normalize();
//    if (std::find(n1.begin(), n1.end(), Polynomial({2,1})) == n1.end()
//        || std::find(n1.begin(), n1.end(), Polynomial({ -2,1 })) == n1.end())
//    {
//	    std::cout << "!!!!!!!!!!!!!! Test 1 not passed !!!!!!!!!!!!!!\n";
//    }
//
//    Polynomial p2({ -16,0,0,0,9 });
//    std::vector<Polynomial> n2 = p2.normalize();
//    if (std::find(n2.begin(), n2.end(), Polynomial({ 4,0,3 })) == n2.end()
//        || std::find(n2.begin(), n2.end(), Polynomial({ -4,0,3 })) == n2.end())
//    {
//        std::cout << "!!!!!!!!!!!!!! Test 2 not passed !!!!!!!!!!!!!!\n";
//    }
//
//	Polynomial p3({ 8,0,0,1 });
//	std::vector<Polynomial> n3 = p3.normalize();
//	if (std::find(n3.begin(), n3.end(), Polynomial({2,1})) == n3.end()
//        || std::find(n3.begin(), n3.end(), Polynomial({ 4,-2,1 })) == n3.end())
//	{
//		std::cout << "!!!!!!!!!!!!!! Test 3 not passed !!!!!!!!!!!!!!\n";
//	}
//
//	Polynomial p4({ -1,0,0,27 });
//	std::vector<Polynomial> n4 = p4.normalize();
//	if (std::find(n4.begin(), n4.end(), Polynomial({ -1,3 })) == n4.end()
//        || std::find(n4.begin(), n4.end(), Polynomial({ 1,3,9 })) == n4.end())
//	{
//		std::cout << "!!!!!!!!!!!!!! Test 4 not passed !!!!!!!!!!!!!!\n";
//	}
//
//	Polynomial p5({ 6,5,1 });
//	std::vector<Polynomial> n5 = p5.normalize();
//    if (std::find(n5.begin(), n5.end(), Polynomial({ 2,1 })) == n5.end()
//        || std::find(n5.begin(), n5.end(), Polynomial({ 3,1 })) == n5.end())
//    {
//        std::cout << "!!!!!!!!!!!!!! Test 5 not passed !!!!!!!!!!!!!!\n";
//    }
//
//	Polynomial p6({ 3,-7,2 });
//	std::vector<Polynomial> n6 = p6.normalize();
//	if (std::find(n6.begin(), n6.end(), Polynomial({-1,2})) == n6.end()
//        || std::find(n6.begin(), n6.end(), Polynomial({ -3,1 })) == n6.end())
//	{
//		std::cout << "!!!!!!!!!!!!!! Test 6 not passed !!!!!!!!!!!!!!\n";
//	}
//
//	Polynomial p7({ -6,-1,1 });
//	std::vector<Polynomial> n7 = p7.normalize();
//	if (std::find(n7.begin(), n7.end(), Polynomial({ 2,1 })) == n7.end()
//        || std::find(n7.begin(), n7.end(), Polynomial({ -3,1 })) == n7.end())
//	{
//		std::cout << "!!!!!!!!!!!!!! Test 7 not passed !!!!!!!!!!!!!!\n";
//
//	}
//
//	Polynomial p8({ 4,0,-3,1 });
//	std::vector<Polynomial> n8 = p8.normalize();
//	if (std::find(n8.begin(), n8.end(), Polynomial({ 1,1 })) == n8.end()
//        || std::find(n8.begin(), n8.end(), Polynomial({ -2,1 })) == n8.end())
//	{
//		std::cout << "!!!!!!!!!!!!!! Test 8 not passed !!!!!!!!!!!!!!\n";
//	}
//
//	Polynomial p9({ 4,0,-5,0,1 });
//	std::vector<Polynomial> n9 = p9.normalize();
//	if (std::find(n9.begin(), n9.end(), Polynomial({ 1,1 })) == n9.end()
//        || std::find(n9.begin(), n9.end(), Polynomial({ -1,1 })) == n9.end()
//        || std::find(n9.begin(), n9.end(), Polynomial({ 2,1 })) == n9.end()
//        || std::find(n9.begin(), n9.end(), Polynomial({ -2,1 })) == n9.end())
//	{
//		std::cout << "!!!!!!!!!!!!!! Test 9 not passed !!!!!!!!!!!!!!\n";
//	}
//
//	Polynomial p10({ -12,-4,3,1 });
//	std::vector<Polynomial> n10 = p10.normalize();
//    if (std::find(n10.begin(), n10.end(), Polynomial({ 3,1 })) == n10.end()
//        || std::find(n10.begin(), n10.end(), Polynomial({ 2,1 })) == n10.end()
//        || std::find(n10.begin(), n10.end(), Polynomial({ -2,1 })) == n10.end())
//    {
//		std::cout << "!!!!!!!!!!!!!! Test 10 not passed !!!!!!!!!!!!!!\n";
//    }
//
//	Polynomial p11({ 0,3,11,6 });
//	std::vector<Polynomial> n11 = p11.normalize();
//	if (std::find(n11.begin(), n11.end(), Polynomial({ 0,1 })) == n11.end()
//        || std::find(n11.begin(), n11.end(), Polynomial({ 1,3 })) == n11.end()
//        || std::find(n11.begin(), n11.end(), Polynomial({ 3,2 })) == n11.end())
//	{
//		std::cout << "!!!!!!!!!!!!!! Test 11 not passed !!!!!!!!!!!!!!\n";
//	}
//
//	Polynomial p12({ 1,1,1 });
//	std::vector<Polynomial> n12 = p12.normalize();
//    if (std::find(n12.begin(), n12.end(), Polynomial({ 1,1,1 })) == n12.end())
//    {
//		std::cout << "!!!!!!!!!!!!!! Test 12 not passed !!!!!!!!!!!!!!\n";
//    }
//}

bool Polynomial::isZero() const {
	return coefficients.empty() || (coefficients.size() == 1 && coefficients[0] == 0);
}

Polynomial Polynomial::operator+(const Polynomial& other) const {
	PROFILE_FUNCTION
	/*std::vector<Rational> result = coefficients;
	for (size_t i = 0; i < other.coefficients.size(); i++) {
		while (result.size() <= i) result.emplace_back(0);
		result[i] += other.coefficients[i];
	}
	return result;*/
	const int maxDeg = this->degree > other.degree ? this->degree : other.degree;
	std::vector<Rational> resultCoeffs(maxDeg + 1);
	for (int i = 0; i <= maxDeg; ++i) {
		resultCoeffs[i] = this->coeff(i) + other.coeff(i);
	}
	return resultCoeffs;
}

Polynomial Polynomial::operator-(const Polynomial& other) const {
	PROFILE_FUNCTION
	/*std::vector<Rational> result = coefficients;
	for (size_t i = 0; i < other.coefficients.size(); i++) {
		while (result.size() <= i) result.emplace_back(0);
		result[i] -= other.coefficients[i];
	}
	return result;*/
	const int maxDeg = this->degree > other.degree ? this->degree : other.degree;
	std::vector<Rational> resultCoeffs(maxDeg + 1);
	for (int i = 0; i <= maxDeg; ++i) {
		resultCoeffs[i] = this->coeff(i) - other.coeff(i);
	}
	return resultCoeffs;
}

Polynomial Polynomial::operator*(const Polynomial& other) const {
	PROFILE_FUNCTION
	/*Polynomial result;
	for (size_t i = 0; i < coefficients.size(); i++) {
		for (size_t j = 0; j < other.coefficients.size(); j++) {
			while (result.coefficients.size() <= i + j) result.coefficients.emplace_back(0);
			result.coefficients[i + j] += coefficients[i] * other.coefficients[j];
		}
	}

	result.degree = static_cast<int>(result.coefficients.size()) - 1;

	result = polyTrim(result);

	return result;*/
	if (this->isZero() || other.isZero()) {
		return Rational(0);
	}
	const int newDeg = this->degree + other.degree;
	std::vector<Rational> resultCoeffs(newDeg + 1, Rational(0));
	for (int i = 0; i <= this->degree; ++i) {
		for (int j = 0; j <= other.degree; ++j) {
			resultCoeffs[i + j] += (this->coeff(i) * other.coeff(j));
		}
	}
	return resultCoeffs;
}

Polynomial Polynomial::operator/(const Polynomial& other) const {
	PROFILE_FUNCTION
	if (other.isZero()) {
		throw std::runtime_error("Division by zero polynomial!");
	}
	if (this->isZero()) {
		return Rational(0);
	}

	auto [quotient, remainder] = polyDivide(*this, other);

	Polynomial result = polyTrim(quotient);
	return result;
}

Polynomial Polynomial::operator/(const Rational& scalar) const {
	if (scalar == 0) {
		throw std::runtime_error("Division by zero scalar in polynomial.");
	}
	std::vector<Rational> resultCoeffs(coefficients.size());
	for (size_t i = 0; i < coefficients.size(); ++i) {
		resultCoeffs[i] = coefficients[i] / scalar;
	}
	return resultCoeffs;
}

bool Polynomial::operator==(const Polynomial& other) const {
	PROFILE_FUNCTION
	if (degree != other.degree) return false;

	for (size_t i = 0; i < coefficients.size(); i++) {
		if (coefficients[i] != other.coefficients[i]) return false;
	}
	return true;
}

std::string Polynomial::toString() const {
	PROFILE_FUNCTION
	if (coefficients.empty()) return "0";
	std::string output;
	const int coeffLength = static_cast<int>(coefficients.size());
	for (int i = coeffLength - 1; i >= 0; i--) {
		if (coefficients[i] == 0) continue;
		if (coefficients[i] >= 0 && !output.empty()) output += " + ";
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

Rational Polynomial::evaluate(const Rational& x) const {
	PROFILE_FUNCTION
	Rational result = 0;
	const int coeffLength = static_cast<int>(coefficients.size());
	for (int i = coeffLength - 1; i >= 0; --i) {
		result = result * x + coefficients[i];
	}
	return result;
}

Polynomial Polynomial::derivative() const {
	PROFILE_FUNCTION
	if (isZero()) return {};
	std::vector<Rational> result(degree, 0);
	for (int i = 1; i <= degree; ++i) {
		result[i - 1] = coefficients[i] * i;
	}
	return {result};
}

Polynomial Polynomial::polyTrim(const Polynomial& poly) {
	PROFILE_FUNCTION
	Polynomial result = poly;
	while (!result.coefficients.empty() && result.coefficients.back().abs() < 1e-9) {
		result.coefficients.pop_back();
	}
	result.degree = result.coefficients.empty() ? -1 : static_cast<int>(result.coefficients.size()) - 1;
	return result;
}

std::pair<std::vector<Rational>, std::vector<Rational>> Polynomial::polyDivide(
	const Polynomial& dividend, const Polynomial& divisor) {
	PROFILE_FUNCTION
	const Polynomial a = polyTrim(dividend);
	const Polynomial b = polyTrim(divisor);

	if (b.coefficients.empty()) {
		throw std::runtime_error("Division by zero polynomial!");
	}

	const int degA = a.degree;
	const int degB = b.degree;

	std::vector<Rational> quotient(degA - degB + 1, 0);
	std::vector<Rational> remainder = a.coefficients;

	// Perform polynomial long division
	while (remainder.size() >= b.coefficients.size() && !remainder.empty()) {
		const size_t degR = remainder.size() - 1;
		Rational factor = remainder.back() / b.coefficients.back();
		const size_t power = degR - degB;
		quotient[power] = factor;

		// Subtract factor * (B shifted by "power") from remainder.
		std::vector<Rational> sub(power, 0); // zeros for lower-degree terms
		for (const Rational& coeff : b.coefficients) {
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
	PROFILE_FUNCTION
	std::vector<Rational> neg(poly.size());
	for (size_t i = 0; i < poly.size(); i++) {
		neg[i] = -poly[i];
	}
	return neg;
}

Polynomial Polynomial::reflectY() const {
	PROFILE_FUNCTION
	std::vector<Rational> result = coefficients;
	const int coeffLength = static_cast<int>(result.size());
	for (int i = coeffLength - 1; i >= 0; i--) {
		if (i % 2 == 1) {
			result[i] = -result[i];
		}
	}
	return {result};
}

// Generate the Sturm sequence
std::vector<Polynomial> Polynomial::sturmSequence(const Polynomial& p) {
	PROFILE_FUNCTION
	std::vector<Polynomial> seq;

	const Polynomial s0 = polyTrim(p);
	if (s0.coefficients.empty()) {
		throw std::runtime_error("Zero polynomial provided.");
	}
	seq.push_back(s0);

	// S1 = p'(x)
	const Polynomial s1 = s0.derivative();
	seq.push_back(s1);

	// Compute subsequent sequence members until remainder becomes zero.
	while (true) {
		// Let S_prev = S_{i-1} and S_curr = S_i.
		const Polynomial& sPrev = seq[seq.size() - 2];
		const Polynomial& sCurr = seq.back();

		if (sCurr.degree == 0) {
			break;
		}

		// Compute polynomial remainder of S_prev divided by S_curr.
		auto divRes = polyDivide(sPrev, sCurr);
		std::vector<Rational> rem = divRes.second;

		// If remainder is zero, stop.
		if (rem.empty() || (rem.size() == 1 && rem[0].abs() < 1e-9)) {
			break;
		}

		// Next term in Sturm sequence is the negative of the remainder.
		std::vector<Rational> next = polyNegate(rem);
		seq.emplace_back(next);
	}

	sturm_sequence = seq;
	return seq;
}
