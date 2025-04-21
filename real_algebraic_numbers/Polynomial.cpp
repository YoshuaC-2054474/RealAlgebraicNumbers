#include "Polynomial.h"
#include <iostream>
#include <numeric>
#include <boost/multiprecision/cpp_int.hpp>
//#include <fmpz_poly.h>
//#include <fmpz_poly_factor.h>

Polynomial::Polynomial(const std::initializer_list<int> coeffs)
{
    //int exponent = 0;
    for (auto it = coeffs.begin(); it < coeffs.end(); ++it) {
        coefficients.push_back(*it);
    }

    degree = coefficients.empty() ? -1 : coefficients.size() - 1;
    //normalize();
}

Polynomial::Polynomial(const std::vector<Rational>& coeffs)
{
	coefficients = coeffs;
	/*while (!coefficients.empty() && coefficients.back() == 0) {
		coefficients.pop_back();
	}*/

	degree = coefficients.empty() ? -1 : coefficients.size() - 1;
    //normalize();
}

Rational findGCD(const std::vector<Rational>& arr) {
    Rational res = arr[0];
     //res.toInts();
    for (int i = 1; i < arr.size(); i++) {
        res = arr[i].gcd(res);// gcd(arr[i], res);
        //arr[i].toInts();
        //res.toInts();
        // If res becomes 1 at any iteration then it remains 1
        // So no need to check further
        if (res == 1)
            return 1;
    }
    //res.toInts();
    //std::cout << "\n";
    return res;
}

bool isPerfectCube(const int n) {
    if (n == 0) {
        return true; 
    }
    const int abs_n = abs(n);
    int low = 1, high = abs_n;
    while (low <= high) {
        int mid = low + (high - low) / 2;
        int cube = mid * mid * mid;
        if (cube == abs_n) {
            mid = (n < 0) ? -mid : mid;
            return (mid * mid * mid == n);
        }
        else if (cube < abs_n) {
            low = mid + 1;
        }
        else {
            high = mid - 1;
        }
    }
    return false;
}

// Helper to calculate integer roots via Rational Root Theorem
std::vector<int> getPossibleRoots(const std::vector<Rational>& poly) {
    std::vector<int> roots;
    if (poly.empty()) return roots;

    const int constant = poly[0].toInts();
    const int leading = poly.back().toInts();

    std::vector<int> constant_factors, leading_factors;

    // Generate factors of constant term (include negatives)
    for (int i = 1; i <= abs(constant); ++i) {
        if (constant % i == 0) {
            constant_factors.push_back(i);
            constant_factors.push_back(-i);
        }
    }

    // Generate factors of leading coefficient (include negatives)
    for (int i = 1; i <= abs(leading); ++i) {
        if (leading % i == 0) {
            leading_factors.push_back(i);
            leading_factors.push_back(-i);
        }
    }

    // Generate possible roots (p/q)
    for (int p : constant_factors) {
        for (int q : leading_factors) {
            if (q != 0 && p % q == 0) {
                int root = p / q;
                if (find(roots.begin(), roots.end(), root) == roots.end()) {
                    roots.push_back(root);
                }
            }
        }
    }
    return roots;
}

// Helper to perform synthetic division
std::vector<Rational> syntheticDivide(const std::vector<Rational>& poly, int root) {
    std::vector<Rational> quotient;
    if (poly.empty()) return quotient;

    quotient.resize(poly.size() - 1);
    Rational carry = 0;

    for (int i = poly.size() - 1; i > 0; --i) {
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
std::vector<Polynomial> get_minimal_polynomials(const Polynomial& poly) {
    std::vector<Polynomial> result;
    //if (poly.coefficients.size() <= 1) return result;

	// find lowest non-zero coefficient
    int lowestNz = 0;
	for (int i = 0; i < poly.coefficients.size(); i++) {
		if (poly.coefficients[i] != 0) {
			lowestNz = i;
			break;
		}
	}

    std::vector<std::pair<Rational, int>> non_zero_terms;
    for (int exponent = 0; exponent < poly.coefficients.size(); ++exponent) {
        Rational coeff = poly.coefficients[exponent];
        if (coeff != 0) {
            non_zero_terms.emplace_back(coeff, exponent);
        }
    }

	if (lowestNz > 0) {
        std::vector<Rational> lCoeffs;
        for (size_t i = 0; i < lowestNz; i++)
        {
			lCoeffs.push_back(0);
        }
		lCoeffs.push_back(1);
        result.push_back(lCoeffs);

		std::vector<Rational> newCoeffs;
        for (int i = lowestNz; i < poly.coefficients.size(); i++) {
			newCoeffs.push_back(poly.coefficients[i]);
		}
		Polynomial newPoly(newCoeffs);

        std::cout << poly.toString() << " -> (" << result[0].toString() << ") (" << newPoly.toString() << ")\n";
		//result.extend(get_minimal_polynomials(newPoly));
		std::vector<Polynomial> subMinimals = get_minimal_polynomials(newPoly);
        result.reserve(result.size() + distance(subMinimals.begin(), subMinimals.end()));
        result.insert(result.end(), subMinimals.begin(), subMinimals.end());
	}
	else if (non_zero_terms.size() == 2) {
		// Check if both exponents are even (check difference of squares)
        if (non_zero_terms[0].second % 2 == 0 && non_zero_terms[1].second % 2 == 0) {
            // Check if the terms have opposite signs
            if ((non_zero_terms[0].first > 0) != (non_zero_terms[1].first > 0)) {
                // Check if absolute values of coefficients are perfect squares
                Rational abs_coeff1 = non_zero_terms[0].first.abs();
                Rational abs_coeff2 = non_zero_terms[1].first.abs();

				Rational sqrt_coeff1 = abs_coeff1.sqrt();
                bool coeff1_isSquare = sqrt_coeff1.isInteger();

				Rational sqrt_coeff2 = abs_coeff2.sqrt();
				bool coeff2_isSquare = sqrt_coeff2.isInteger();

				if (coeff1_isSquare && coeff2_isSquare) {
					// Create the difference of squares polynomial
					int halvedExponent1 = non_zero_terms[0].second / 2;
					int halvedExponent2 = non_zero_terms[1].second / 2;
					std::vector<Rational> diff_of_squares_coeffs1;
					diff_of_squares_coeffs1.resize(std::max(halvedExponent1, halvedExponent2) + 1, 0);
					diff_of_squares_coeffs1[halvedExponent1] = sqrt_coeff1;
					diff_of_squares_coeffs1[halvedExponent2] = sqrt_coeff2;
                    std::vector<Rational> diff_of_squares_coeffs2;
                    diff_of_squares_coeffs2.resize(std::max(halvedExponent1, halvedExponent2) + 1, 0);
                    diff_of_squares_coeffs2[halvedExponent1] = -sqrt_coeff1;
                    diff_of_squares_coeffs2[halvedExponent2] = sqrt_coeff2;

					Polynomial newPoly1(diff_of_squares_coeffs1);
					Polynomial newPoly2(diff_of_squares_coeffs2);
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
    	else if (non_zero_terms[0].second % 3 == 0 && non_zero_terms[1].second % 3 == 0) {
            if (non_zero_terms[0].first.isInteger() && non_zero_terms[1].first.isInteger())
            {
				int coeff1 = non_zero_terms[0].first.toInts();
				int coeff2 = non_zero_terms[1].first.toInts();
				if (isPerfectCube(coeff1) && isPerfectCube(coeff2))
				{
					// Create the sum/difference of cubes polynomial
					int halvedExponent1 = non_zero_terms[0].second / 3;
					int halvedExponent2 = non_zero_terms[1].second / 3;
                    int b = static_cast<int>(std::cbrt(abs(coeff1)));
                    int a = static_cast<int>(std::cbrt(abs(coeff2)));
					b = (coeff1 < 0) ? -b : b;
					a = (coeff2 < 0) ? -a : a;

					std::vector<Rational> sum_of_cubes_coeffs1;
					sum_of_cubes_coeffs1.resize(std::max(halvedExponent1, halvedExponent2) + 1, 0);
					sum_of_cubes_coeffs1[halvedExponent1] = b;
					sum_of_cubes_coeffs1[halvedExponent2] = a;

					std::vector<Rational> sum_of_cubes_coeffs2;
					sum_of_cubes_coeffs2.resize(std::max(halvedExponent1*2, halvedExponent2*2) + 1, 0);
					sum_of_cubes_coeffs2[halvedExponent1*2] = std::pow(b, 2);
					sum_of_cubes_coeffs2[halvedExponent2*2] = std::pow(a, 2);
					int middleTerm = coeff1 < 0 ? a * b : -a * b;
					sum_of_cubes_coeffs2[halvedExponent1+halvedExponent2] = -a*b;

					Polynomial newPoly1(sum_of_cubes_coeffs1);
					Polynomial newPoly2(sum_of_cubes_coeffs2);
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
    else if (non_zero_terms.size() == 3 && poly.degree == 2) {
        Rational c = poly.coefficients[0];
        Rational b = poly.coefficients[1];
        Rational a = poly.coefficients[2];

        // Function to generate all divisors of x, including negatives
        auto getDivisors = [](const Rational& x) {
            std::vector<int> divisors;
            if (x == 0) return divisors;
            const int abs_x = abs(x.toInts());
            for (int i = 1; i <= abs_x; ++i) {
                if (abs_x % i == 0) {
                    divisors.push_back(i);
                    divisors.push_back(-i);
                }
            }
            return divisors;
        };

        std::vector<int> divisors_a = getDivisors(a);
        std::vector<int> divisors_c = getDivisors(c);

        bool factorFound = false;

        // Iterate all possible pairs (m from divisors of a, n from divisors of c)
   //     for (int m : divisors_a) {
			//if (factorFound) break; // Stop if a factor has already been found

   //         if (m == 0) continue;
   //         Rational p = a / m; // m * p must equal a

   //         for (int n : divisors_c) {
   //             if (factorFound) break; // Stop if a factor has already been found

   //             if (n == 0) continue;
   //             Rational q = c / n; // n * q must equal c

   //             // Check if the combination satisfies mq + np = b
   //             if (m * q + n * p == b) {
   //                 // Construct the factors (mx + n) and (px + q)
   //                 std::vector<Rational> factor1 = { n, m }; // Represents mx + n
   //                 std::vector<Rational> factor2 = { q, p }; // Represents px + q

			//		Polynomial newPoly1(factor1);
			//		Polynomial newPoly2(factor2);

			//		std::cout << poly.toString() << " -> (" << newPoly1.toString() << ") (" << newPoly2.toString() << ")\n";

   //                 std::vector<Polynomial> subMinimals1 = get_minimal_polynomials(newPoly1);
   //                 result.reserve(result.size() + distance(subMinimals1.begin(), subMinimals1.end()));
   //                 result.insert(result.end(), subMinimals1.begin(), subMinimals1.end());

   //                 std::vector<Polynomial> subMinimals2 = get_minimal_polynomials(newPoly2);
   //                 result.reserve(result.size() + distance(subMinimals2.begin(), subMinimals2.end()));
   //                 result.insert(result.end(), subMinimals2.begin(), subMinimals2.end());

   //             	factorFound = true;
   //             }
   //         }
   //     }
        for (int p : divisors_c) {
            if (p == 0) continue;
            for (int q : divisors_a) {
                if (q == 0) continue;
                // Check if a*p² + b*p*q + c*q² == 0
                if (a * p * p + b * p * q + c * q * q == 0) {
                    Rational k = a / q;
                    Rational m = -c / p;

                    Polynomial newPoly1({ -p, q });
                    Polynomial newPoly2({ m, k });

                    std::cout << poly.toString() << " -> (" << newPoly1.toString() << ") (" << newPoly2.toString() << ")\n";

                    std::vector<Polynomial> subMinimals1 = get_minimal_polynomials(newPoly1);
                    result.reserve(result.size() + distance(subMinimals1.begin(), subMinimals1.end()));
                    result.insert(result.end(), subMinimals1.begin(), subMinimals1.end());

                    std::vector<Polynomial> subMinimals2 = get_minimal_polynomials(newPoly2);
                    result.reserve(result.size() + distance(subMinimals2.begin(), subMinimals2.end()));
                    result.insert(result.end(), subMinimals2.begin(), subMinimals2.end());
                    factorFound = true;
                }
            }
        }
    }
    else if (non_zero_terms.size() >= 3 /*&& poly.coefficients.size() < 4*/)
    {
        std::vector<int> possible_roots = getPossibleRoots(poly.coefficients);
        for (int root : possible_roots) {
            std::vector<Rational> quotient = syntheticDivide(poly.coefficients, root);
            if (!quotient.empty()) {
                Polynomial newPoly1({ -root, 1 }); // Factor (x - root)
				Polynomial newPoly2(quotient); // Quotient polynomial

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

std::vector<Polynomial> Polynomial::normalize()
{
	if (isNomalized) return {};

	while (coefficients.size() > 1 && coefficients.back() == 0) {
		coefficients.pop_back();
	}
	degree = coefficients.empty() ? -1 : coefficients.size() - 1;

	if (coefficients.empty()) return {};

    const Rational gcd = findGCD(coefficients);
	int gcdInt = gcd.toInts();
    if (gcd > 1)
    {
        for (int i = 0; i < coefficients.size(); i++) {
            coefficients[i] /= gcd;
        }
    }

    /*fmpz_poly_t f;
    fmpz_poly_init(f);

    for (int i = 0; i < coefficients.size(); i++)
    {
		fmpz_poly_set_coeff_si(f, i, coefficients[i].toInts());
    }
	std::cout << "\n\nNormalized Polynomial: " << toString() << "\n";
    fmpz_poly_print(f);*/

    //const Rational highestCo = coefficients.back();
    //int highestCoInt = highestCo.toInts();
    ////std::cout << "\n";
    //if (highestCo != 0)
    //{
    //    for (int i = 0; i < coefficients.size(); i++) {
    //        coefficients[i] /= highestCo;
    //    }
    //}

	std::vector<Polynomial> minimalPolys = get_minimal_polynomials(*this);
	std::cout << "\n\nMinimal Polynomials of " << toString() << "\n";
	for (const Polynomial& poly : minimalPolys) {
		std::cout << "\t" << poly.toString() << '\n';
	}
	std::cout << "\n\n";

	isNomalized = true;
    return minimalPolys;
}

void Polynomial::testNormalize()
{
	Polynomial p1({ -4,0,1 });
	std::vector<Polynomial> n1 = p1.normalize();
    if (std::find(n1.begin(), n1.end(), Polynomial({2,1})) == n1.end()
        || std::find(n1.begin(), n1.end(), Polynomial({ -2,1 })) == n1.end())
    {
	    std::cout << "!!!!!!!!!!!!!! Test 1 not passed !!!!!!!!!!!!!!\n";
    }

    Polynomial p2({ -16,0,0,0,9 });
    std::vector<Polynomial> n2 = p2.normalize();
    if (std::find(n2.begin(), n2.end(), Polynomial({ 4,0,3 })) == n2.end()
        || std::find(n2.begin(), n2.end(), Polynomial({ -4,0,3 })) == n2.end())
    {
        std::cout << "!!!!!!!!!!!!!! Test 2 not passed !!!!!!!!!!!!!!\n";
    }

	Polynomial p3({ 8,0,0,1 });
	std::vector<Polynomial> n3 = p3.normalize();
	if (std::find(n3.begin(), n3.end(), Polynomial({2,1})) == n3.end()
        || std::find(n3.begin(), n3.end(), Polynomial({ 4,-2,1 })) == n3.end())
	{
		std::cout << "!!!!!!!!!!!!!! Test 3 not passed !!!!!!!!!!!!!!\n";
	}

	Polynomial p4({ -1,0,0,27 });
	std::vector<Polynomial> n4 = p4.normalize();
	if (std::find(n4.begin(), n4.end(), Polynomial({ -1,3 })) == n4.end()
        || std::find(n4.begin(), n4.end(), Polynomial({ 1,3,9 })) == n4.end())
	{
		std::cout << "!!!!!!!!!!!!!! Test 4 not passed !!!!!!!!!!!!!!\n";
	}

	Polynomial p5({ 6,5,1 });
	std::vector<Polynomial> n5 = p5.normalize();
    if (std::find(n5.begin(), n5.end(), Polynomial({ 2,1 })) == n5.end()
        || std::find(n5.begin(), n5.end(), Polynomial({ 3,1 })) == n5.end())
    {
        std::cout << "!!!!!!!!!!!!!! Test 5 not passed !!!!!!!!!!!!!!\n";
    }

	Polynomial p6({ 3,-7,2 });
	std::vector<Polynomial> n6 = p6.normalize();
	if (std::find(n6.begin(), n6.end(), Polynomial({-1,2})) == n6.end()
        || std::find(n6.begin(), n6.end(), Polynomial({ -3,1 })) == n6.end())
	{
		std::cout << "!!!!!!!!!!!!!! Test 6 not passed !!!!!!!!!!!!!!\n";
	}

	Polynomial p7({ -6,-1,1 });
	std::vector<Polynomial> n7 = p7.normalize();
	if (std::find(n7.begin(), n7.end(), Polynomial({ 2,1 })) == n7.end()
        || std::find(n7.begin(), n7.end(), Polynomial({ -3,1 })) == n7.end())
	{
		std::cout << "!!!!!!!!!!!!!! Test 7 not passed !!!!!!!!!!!!!!\n";

	}

	Polynomial p8({ 4,0,-3,1 });
	std::vector<Polynomial> n8 = p8.normalize();
	if (std::find(n8.begin(), n8.end(), Polynomial({ 1,1 })) == n8.end()
        || std::find(n8.begin(), n8.end(), Polynomial({ -2,1 })) == n8.end())
	{
		std::cout << "!!!!!!!!!!!!!! Test 8 not passed !!!!!!!!!!!!!!\n";
	}

	Polynomial p9({ 4,0,-5,0,1 });
	std::vector<Polynomial> n9 = p9.normalize();
	if (std::find(n9.begin(), n9.end(), Polynomial({ 1,1 })) == n9.end()
        || std::find(n9.begin(), n9.end(), Polynomial({ -1,1 })) == n9.end()
        || std::find(n9.begin(), n9.end(), Polynomial({ 2,1 })) == n9.end()
        || std::find(n9.begin(), n9.end(), Polynomial({ -2,1 })) == n9.end())
	{
		std::cout << "!!!!!!!!!!!!!! Test 9 not passed !!!!!!!!!!!!!!\n";
	}

	Polynomial p10({ -12,-4,3,1 });
	std::vector<Polynomial> n10 = p10.normalize();
    if (std::find(n10.begin(), n10.end(), Polynomial({ 3,1 })) == n10.end()
        || std::find(n10.begin(), n10.end(), Polynomial({ 2,1 })) == n10.end()
        || std::find(n10.begin(), n10.end(), Polynomial({ -2,1 })) == n10.end())
    {
		std::cout << "!!!!!!!!!!!!!! Test 10 not passed !!!!!!!!!!!!!!\n";
    }

	Polynomial p11({ 0,3,11,6 });
	std::vector<Polynomial> n11 = p11.normalize();
	if (std::find(n11.begin(), n11.end(), Polynomial({ 0,1 })) == n11.end()
        || std::find(n11.begin(), n11.end(), Polynomial({ 1,3 })) == n11.end()
        || std::find(n11.begin(), n11.end(), Polynomial({ 3,2 })) == n11.end())
	{
		std::cout << "!!!!!!!!!!!!!! Test 11 not passed !!!!!!!!!!!!!!\n";
	}

	Polynomial p12({ 1,1,1 });
	std::vector<Polynomial> n12 = p12.normalize();
    if (std::find(n12.begin(), n12.end(), Polynomial({ 1,1,1 })) == n12.end())
    {
		std::cout << "!!!!!!!!!!!!!! Test 12 not passed !!!!!!!!!!!!!!\n";
    }
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

bool Polynomial::operator==(const Polynomial& other) const
{
	if (degree != other.degree) return false;

	for (int i = 0; i < coefficients.size(); i++)
	{
		if (coefficients[i] != other.coefficients[i]) return false;
	}
	return true;
}

std::string Polynomial::toString() const {
    if (coefficients.empty()) return "0";
    std::string output;
    for (int i = coefficients.size() - 1; i >= 0; i--) {
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

        if (S_curr.degree == 0) {
            break;
        }

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