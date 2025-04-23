#include "Polynomial.h"
#include <iostream>
#include <numeric>
#include <boost/multiprecision/cpp_int.hpp>
//#include <fmpz_poly.h>
//#include <fmpz_poly_factor.h>

Polynomial::Polynomial(const std::initializer_list<int> coeffs)
{
    for (auto it = coeffs.begin(); it < coeffs.end(); ++it) {
        coefficients.emplace_back(*it);
    }

    degree = coefficients.empty() ? -1 : static_cast<int>(coefficients.size()) - 1;
}

Polynomial::Polynomial(const std::vector<Rational>& coeffs)
{
	coefficients = coeffs;

	degree = coefficients.empty() ? -1 : static_cast<int>(coefficients.size()) - 1;
}

static Rational findGCD(const std::vector<Rational>& arr) {
    Rational res = arr[0];

    for (int i = 1; i < arr.size(); i++) {
        res = arr[i].gcd(res);
        if (res == 1)
            return 1;
    }

    return res;
}

static bool isPerfectCube(const int n) {
    if (n == 0) {
        return true; 
    }
    const int absN = abs(n);
    int low = 1, high = absN;
    while (low <= high) {
        int mid = low + (high - low) / 2;
        if (const int cube = mid * mid * mid; cube == absN) {
            mid = (n < 0) ? -mid : mid;
            return (mid * mid * mid == n);
        }
        else if (cube < absN) {
            low = mid + 1;
        }
        else {
            high = mid - 1;
        }
    }
    return false;
}

// Helper to calculate integer roots via Rational Root Theorem
static std::vector<int> getPossibleRoots(const std::vector<Rational>& poly) {
    std::vector<int> roots;
    if (poly.empty()) return roots;

    const int constant = poly[0].toInts();
    const int leading = poly.back().toInts();

    std::vector<int> constantFactors, leadingFactors;

    // Generate factors of constant term (include negatives)
    for (int i = 1; i <= abs(constant); ++i) {
        if (constant % i == 0) {
            constantFactors.push_back(i);
            constantFactors.push_back(-i);
        }
    }

    // Generate factors of leading coefficient (include negatives)
    for (int i = 1; i <= abs(leading); ++i) {
        if (leading % i == 0) {
            leadingFactors.push_back(i);
            leadingFactors.push_back(-i);
        }
    }

    // Generate possible roots (p/q)
    for (const int p : constantFactors) {
        for (const int q : leadingFactors) {
            if (q != 0 && p % q == 0) {
                if (int root = p / q; std::ranges::find(roots.begin(), roots.end(), root) == roots.end()) {
                    roots.push_back(root);
                }
            }
        }
    }
    return roots;
}

// Helper to perform synthetic division
static std::vector<Rational> syntheticDivide(const std::vector<Rational>& poly, const int root) {
    std::vector<Rational> quotient;
    if (poly.empty()) return quotient;

    quotient.resize(poly.size() - 1);
    Rational carry = 0;

    for (size_t i = poly.size() - 1; i > 0; --i) {
        quotient[i - 1] = poly[i] + carry;
        carry = quotient[i - 1] * root;
    }

    // Check if remainder is zero
    if (poly[0] + carry != 0) {
        return {}; // Not a valid root
    }

    return quotient;
}

// Main function to get minimal polynomials
static std::vector<Polynomial> get_minimal_polynomials(const Polynomial& poly) {
    bool printFactoringSteps = false;

    std::vector<Polynomial> result;
    //if (poly.coefficients.size() <= 1) return result;

	// find lowest non-zero coefficient
    size_t lowestNz = 0;
	for (size_t i = 0; i < poly.coefficients.size(); i++) {
		if (poly.coefficients[i] != 0) {
			lowestNz = i;
			break;
		}
	}

    std::vector<std::pair<Rational, int>> nonZeroTerms;
    for (size_t exponent = 0; exponent < poly.coefficients.size(); ++exponent) {
        if (Rational coeff = poly.coefficients[exponent]; coeff != 0) {
            nonZeroTerms.emplace_back(coeff, exponent);
        }
    }

	if (lowestNz > 0) {
        std::vector<Rational> lCoeffs;
        lCoeffs.reserve(lowestNz+1);
        for (size_t i = 0; i < lowestNz; i++)
        {
			lCoeffs.emplace_back(0);
        }
		lCoeffs.emplace_back(1);
        result.emplace_back(lCoeffs);

		std::vector<Rational> newCoeffs;
        for (size_t i = lowestNz; i < poly.coefficients.size(); i++) {
			newCoeffs.push_back(poly.coefficients[i]);
		}
		Polynomial newPoly(newCoeffs);

        if (printFactoringSteps)
			std::cout << poly.toString() << " -> (" << result[0].toString() << ") (" << newPoly.toString() << ")\n";
		//result.extend(get_minimal_polynomials(newPoly));
		std::vector<Polynomial> subMinimals = get_minimal_polynomials(newPoly);
        result.reserve(result.size() + distance(subMinimals.begin(), subMinimals.end()));
        result.insert(result.end(), subMinimals.begin(), subMinimals.end());
	}
	else if (nonZeroTerms.size() == 2) {
		// Check if both exponents are even (check difference of squares)
        if (nonZeroTerms[0].second % 2 == 0 && nonZeroTerms[1].second % 2 == 0) {
            // Check if the terms have opposite signs
            if ((nonZeroTerms[0].first > 0) != (nonZeroTerms[1].first > 0)) {
                // Check if absolute values of coefficients are perfect squares
                Rational absCoeff1 = nonZeroTerms[0].first.abs();
                Rational absCoeff2 = nonZeroTerms[1].first.abs();

				Rational sqrtCoeff1 = absCoeff1.sqrt();
                bool coeff1IsSquare = sqrtCoeff1.isInteger();

				Rational sqrtCoeff2 = absCoeff2.sqrt();
				bool coeff2IsSquare = sqrtCoeff2.isInteger();

				if (coeff1IsSquare && coeff2IsSquare) {
					// Create the difference of squares polynomial
					int halvedExponent1 = nonZeroTerms[0].second / 2;
					int halvedExponent2 = nonZeroTerms[1].second / 2;
					std::vector<Rational> diffOfSquaresCoeffs1;
					diffOfSquaresCoeffs1.resize(std::max(halvedExponent1, halvedExponent2) + 1, 0);
					diffOfSquaresCoeffs1[halvedExponent1] = sqrtCoeff1;
					diffOfSquaresCoeffs1[halvedExponent2] = sqrtCoeff2;
                    std::vector<Rational> diffOfSquaresCoeffs2;
                    diffOfSquaresCoeffs2.resize(std::max(halvedExponent1, halvedExponent2) + 1, 0);
                    diffOfSquaresCoeffs2[halvedExponent1] = -sqrtCoeff1;
                    diffOfSquaresCoeffs2[halvedExponent2] = sqrtCoeff2;

					Polynomial newPoly1(diffOfSquaresCoeffs1);
					Polynomial newPoly2(diffOfSquaresCoeffs2);

                    if (printFactoringSteps)
						std::cout << poly.toString() << " -> (" << newPoly1.toString() << ") (" << newPoly2.toString() << ")\n";

                    std::vector<Polynomial> subMinimals1 = get_minimal_polynomials(newPoly1);
                    result.reserve(result.size() + distance(subMinimals1.begin(), subMinimals1.end()));
                    result.insert(result.end(), subMinimals1.begin(), subMinimals1.end());

                    std::vector<Polynomial> subMinimals2 = get_minimal_polynomials(newPoly2);
                    result.reserve(result.size() + distance(subMinimals2.begin(), subMinimals2.end()));
                    result.insert(result.end(), subMinimals2.begin(), subMinimals2.end());
				}
            }
        }
    	// Check if both exponents are divisible by 3 (check sum/difference of cubes)
    	else if (nonZeroTerms[0].second % 3 == 0 && nonZeroTerms[1].second % 3 == 0) {
            if (nonZeroTerms[0].first.isInteger() && nonZeroTerms[1].first.isInteger())
            {
				int coeff1 = nonZeroTerms[0].first.toInts();
				int coeff2 = nonZeroTerms[1].first.toInts();
				if (isPerfectCube(coeff1) && isPerfectCube(coeff2))
				{
					// Create the sum/difference of cubes polynomial
					int halvedExponent1 = nonZeroTerms[0].second / 3;
					int halvedExponent2 = nonZeroTerms[1].second / 3;
                    int b = static_cast<int>(std::cbrt(abs(coeff1)));
                    int a = static_cast<int>(std::cbrt(abs(coeff2)));
					b = (coeff1 < 0) ? -b : b;
					a = (coeff2 < 0) ? -a : a;

					std::vector<Rational> sumOfCubesCoeffs1;
					sumOfCubesCoeffs1.resize(std::max(halvedExponent1, halvedExponent2) + 1, 0);
					sumOfCubesCoeffs1[halvedExponent1] = b;
					sumOfCubesCoeffs1[halvedExponent2] = a;

					std::vector<Rational> sumOfCubesCoeffs2;
					sumOfCubesCoeffs2.resize(std::max(halvedExponent1*2, halvedExponent2*2) + 1, 0);
					sumOfCubesCoeffs2[halvedExponent1*2] = std::pow(b, 2);
					sumOfCubesCoeffs2[halvedExponent2*2] = std::pow(a, 2);
					int middleTerm = coeff1 < 0 ? a * b : -a * b;
					sumOfCubesCoeffs2[halvedExponent1+halvedExponent2] = -a * b;

					Polynomial newPoly1(sumOfCubesCoeffs1);
					Polynomial newPoly2(sumOfCubesCoeffs2);

                    if (printFactoringSteps)
						std::cout << poly.toString() << " -> (" << newPoly1.toString() << ") (" << newPoly2.toString() << ")\n";

					std::vector<Polynomial> subMinimals1 = get_minimal_polynomials(newPoly1);
					result.reserve(result.size() + distance(subMinimals1.begin(), subMinimals1.end()));
					result.insert(result.end(), subMinimals1.begin(), subMinimals1.end());

					std::vector<Polynomial> subMinimals2 = get_minimal_polynomials(newPoly2);
					result.reserve(result.size() + distance(subMinimals2.begin(), subMinimals2.end()));
					result.insert(result.end(), subMinimals2.begin(), subMinimals2.end());

				}
            }
        }
    }
    else if (nonZeroTerms.size() == 3 && poly.degree == 2) {
        Rational c = poly.coefficients[0];
        Rational b = poly.coefficients[1];
        Rational a = poly.coefficients[2];

        // Function to generate all divisors of x, including negatives
        auto getDivisors = [](const Rational& x) {
            std::vector<int> divisors;
            if (x == 0) return divisors;
            const int absX = abs(x.toInts());
            for (int i = 1; i <= absX; ++i) {
                if (absX % i == 0) {
                    divisors.push_back(i);
                    divisors.push_back(-i);
                }
            }
            return divisors;
        };

        std::vector<int> divisorsA = getDivisors(a);
        std::vector<int> divisorsC = getDivisors(c);

        //bool factorFound = false;

        for (int p : divisorsC) {
            if (p == 0) continue;
            for (int q : divisorsA) {
                if (q == 0) continue;
                // Check if a*p² + b*p*q + c*q² == 0
                if (a * p * p + b * p * q + c * q * q == 0) {
                    Rational k = a / q;
                    Rational m = -c / p;

                    Polynomial newPoly1({ -p, q });
                    Polynomial newPoly2({ m, k });

                    if (printFactoringSteps)
						std::cout << poly.toString() << " -> (" << newPoly1.toString() << ") (" << newPoly2.toString() << ")\n";

                    std::vector<Polynomial> subMinimals1 = get_minimal_polynomials(newPoly1);
                    result.reserve(result.size() + distance(subMinimals1.begin(), subMinimals1.end()));
                    result.insert(result.end(), subMinimals1.begin(), subMinimals1.end());

                    std::vector<Polynomial> subMinimals2 = get_minimal_polynomials(newPoly2);
                    result.reserve(result.size() + distance(subMinimals2.begin(), subMinimals2.end()));
                    result.insert(result.end(), subMinimals2.begin(), subMinimals2.end());
                    //factorFound = true;
                }
            }
        }
    }
    else if (nonZeroTerms.size() >= 3 /*&& poly.coefficients.size() < 4*/)
    {
        std::vector<int> possibleRoots = getPossibleRoots(poly.coefficients);
        for (int root : possibleRoots) {
            std::vector<Rational> quotient = syntheticDivide(poly.coefficients, root);
            if (!quotient.empty()) {
                Polynomial newPoly1({ -root, 1 }); // Factor (x - root)
				Polynomial newPoly2(quotient); // Quotient polynomial

                if (printFactoringSteps)
					std::cout << poly.toString() << " -> (" << newPoly1.toString() << ") (" << newPoly2.toString() << ")\n";

                std::vector<Polynomial> subMinimals1 = get_minimal_polynomials(newPoly1);
                result.reserve(result.size() + distance(subMinimals1.begin(), subMinimals1.end()));
                result.insert(result.end(), subMinimals1.begin(), subMinimals1.end());

                std::vector<Polynomial> subMinimals2 = get_minimal_polynomials(newPoly2);
                result.reserve(result.size() + distance(subMinimals2.begin(), subMinimals2.end()));
                result.insert(result.end(), subMinimals2.begin(), subMinimals2.end());
                break;
            }
        }
    }

	if (result.empty()) {
		// If no factors found, return the polynomial itself.
		result.push_back(poly);
	}

    return result;
}
void Polynomial::normalize(const Rational& lowerBound, const Rational& upperBound)
{
	if (is_normalized) return;

	while (coefficients.size() > 1 && coefficients.back() == 0) {
		coefficients.pop_back();
	}
	degree = coefficients.empty() ? -1 : static_cast<int>(coefficients.size()) - 1;

	if (coefficients.empty()) return;

	//int gcdInt = gcd.toInts();
    if (const Rational gcd = findGCD(coefficients); gcd > 1)
    {
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

	const std::vector<Polynomial> minimalPolys = get_minimal_polynomials(*this);
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
    return degree == -1;
}

Polynomial Polynomial::operator+(const Polynomial& other) const {
    std::vector<Rational> result = coefficients;
    for (size_t i = 0; i < other.coefficients.size(); i++)
    {
        while (result.size() <= i) result.emplace_back(0);
        result[i] += other.coefficients[i];
    }
    return result;
}

Polynomial Polynomial::operator-(const Polynomial& other) const {
    std::vector<Rational> result = coefficients;
    for (size_t i = 0; i < other.coefficients.size(); i++)
    {
        while (result.size() <= i) result.emplace_back(0);
        result[i] -= other.coefficients[i];
    }
    return result;
}

Polynomial Polynomial::operator*(const Polynomial& other) const {
    Polynomial result;
    for (size_t i = 0; i < coefficients.size(); i++)
    {
        for (size_t j = 0; j < other.coefficients.size(); j++)
        {
            while (result.coefficients.size() <= i + j) result.coefficients.emplace_back(0);
            result.coefficients[i + j] += coefficients[i] * other.coefficients[j];
        }
    }

    result.degree = static_cast<int>(result.coefficients.size()) - 1;
    return result;
}

bool Polynomial::operator==(const Polynomial& other) const {
	if (degree != other.degree) return false;

	for (size_t i = 0; i < coefficients.size(); i++)
	{
		if (coefficients[i] != other.coefficients[i]) return false;
	}
	return true;
}

std::string Polynomial::toString() const {
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

Rational Polynomial::evaluate(const Rational& x) const
{
	Rational result = 0;
	const int coeffLength = static_cast<int>(coefficients.size());
	for (int i = coeffLength - 1; i >= 0; --i) {
		result = result * x + coefficients[i];
	}
	return result;
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
    result.degree = result.coefficients.empty() ? -1 : static_cast<int>(result.coefficients.size()) - 1;
    return result;
}

std::pair<std::vector<Rational>, std::vector<Rational>> Polynomial::polyDivide(const Polynomial& dividend, const Polynomial& divisor) {
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
        std::vector<Rational> sub(power, 0);  // zeros for lower-degree terms
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
    std::vector<Rational> neg(poly.size());
    for (size_t i = 0; i < poly.size(); i++) {
        neg[i] = -poly[i];
    }
    return neg;
}

Polynomial Polynomial::reflectY() const
{
    std::vector<Rational> result = coefficients;
	const int coeffLength = static_cast<int>(result.size());
    for (int i = coeffLength - 1; i >= 0; i--)
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