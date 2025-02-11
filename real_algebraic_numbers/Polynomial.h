#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <Eigen/Dense>

class Polynomial
{
public:
	int degree;
	std::vector<int> coefficients;
	std::vector<std::vector<int>> sturm_sequence;
	std::vector<int> derivative;

	Polynomial();
	Polynomial(const std::vector<int>& coefficients);
	Polynomial computeResultantPolynomial(const Polynomial& other) const;

    Polynomial operator+(const Polynomial& other) const {
        std::vector<int> result(std::max(coefficients.size(), other.coefficients.size()), 0);
        for (size_t i = 0; i < coefficients.size(); ++i) result[i] += coefficients[i];
        for (size_t i = 0; i < other.coefficients.size(); ++i) result[i] += other.coefficients[i];
        return Polynomial(result);
    }

    Polynomial operator*(const Polynomial& other) const {
        if (degree == -1 || other.degree == -1) return Polynomial();
        std::vector<int> result(degree + other.degree + 1, 0);
        for (int i = 0; i <= degree; ++i) {
            for (int j = 0; j <= other.degree; ++j) {
                result[i + j] += coefficients[i] * other.coefficients[j];
            }
        }
        return Polynomial(result);
    }

    Polynomial operator-() const {
        std::vector<int> result;
        for (int c : coefficients) result.push_back(-c);
        return Polynomial(result);
    }

    Polynomial operator-(const Polynomial& other) const {
        return *this + (-other);
    }
private:
	// Function to compute the derivative of a polynomial
	void computeDerivative();

	// Function to compute Sturm sequence of a polynomial
	void computeSturmSequence();
    void generate_sturm_sequence(); //second

	// Polynomial division with remainder (for Sturm sequence)
	std::vector<int> polyMod(const std::vector<int>& dividend, const std::vector<int>& divisor);

    std::vector<int> polyAdd(const std::vector<int>& p1, const std::vector<int>& p2) const;

    std::vector<int> polyMultiply(const std::vector<int>& p1, const std::vector<int>& p2) const;

	Eigen::MatrixXd buildSylvesterMatrix(const Polynomial& a, const Polynomial& b) const;

	std::vector<int> determinantAsPolynomial(Eigen::MatrixXd& s) const;

    
};

#endif