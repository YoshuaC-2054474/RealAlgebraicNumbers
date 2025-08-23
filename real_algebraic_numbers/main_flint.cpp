#include <qqbar.h>
#include <fmpq.h>
#include <string>
#include <iostream>
#include <vector>
#include <chrono>
#include <stdexcept>
#include <memory>
#include "MyTimer.h"

// Using qqbar as our RealAlgebraicNumber
using Ran = qqbar_struct;
using RanPtr = qqbar_ptr;

// Helper function to initialize a qqbar
void initRan(Ran& ran) {
    qqbar_init(&ran);
}

// Helper function to clear a qqbar
void clearRan(Ran& ran) {
    qqbar_clear(&ran);
}

// Helper class for automatic memory management of qqbar
class AutoRan {
public:
    AutoRan() { initRan(value); }
    AutoRan(int n) {
        initRan(value);
        qqbar_set_si(&value, n);
    }
    AutoRan(const AutoRan& other) {
        initRan(value);
        qqbar_set(&value, &other.value);
    }

    // Add assignment operator
    AutoRan& operator=(const AutoRan& other) {
        if (this != &other) {
            qqbar_set(&value, &other.value);
        }
        return *this;
    }

    ~AutoRan() { clearRan(value); }

    Ran* operator->() { return &value; }
    Ran& get() { return value; }
    const Ran& get() const { return value; }

    // Conversion to double
    double toDouble(int prec = 53) const {
        acb_t approx;
        acb_init(approx);
        qqbar_get_acb(approx, &value, prec);
        double result = arf_get_d(arb_midref(acb_realref(approx)), ARF_RND_NEAR);
        acb_clear(approx);
        return result;
    }

    void printApprox() const {
        qqbar_printn(&value, 20);
    }

    // Arithmetic operations
    AutoRan operator+(const AutoRan& other) const {
        AutoRan result;
        qqbar_add(&result.value, &value, &other.value);
        return result;
    }

    AutoRan operator-(const AutoRan& other) const {
        AutoRan result;
        qqbar_sub(&result.value, &value, &other.value);
        return result;
    }

    AutoRan operator*(const AutoRan& other) const {
        AutoRan result;
        qqbar_mul(&result.value, &value, &other.value);
        return result;
    }

    AutoRan operator/(const AutoRan& other) const {
        AutoRan result;
        qqbar_div(&result.value, &value, &other.value);
        return result;
    }

    AutoRan operator-() const {
        AutoRan result;
        qqbar_neg(&result.value, &value);
        return result;
    }

    // Comparison operators
    bool operator==(const AutoRan& other) const {
        return qqbar_equal(&value, &other.value);
    }

    bool operator!=(const AutoRan& other) const {
        return !qqbar_equal(&value, &other.value);
    }

    bool operator<(const AutoRan& other) const {
        return qqbar_cmp_re(&value, &other.value) < 0;
    }

    bool operator<=(const AutoRan& other) const {
        return qqbar_cmp_re(&value, &other.value) <= 0;
    }

    bool operator>(const AutoRan& other) const {
        return qqbar_cmp_re(&value, &other.value) > 0;
    }

    bool operator>=(const AutoRan& other) const {
        return qqbar_cmp_re(&value, &other.value) >= 0;
    }

    // Compound assignment
    AutoRan& operator+=(const AutoRan& other) {
        qqbar_add(&value, &value, &other.value);
        return *this;
    }

    AutoRan& operator-=(const AutoRan& other) {
        qqbar_sub(&value, &value, &other.value);
        return *this;
    }

    AutoRan& operator*=(const AutoRan& other) {
        qqbar_mul(&value, &value, &other.value);
        return *this;
    }

    AutoRan& operator/=(const AutoRan& other) {
        qqbar_div(&value, &value, &other.value);
        return *this;
    }

    // Power function
    AutoRan pow(int exponent) const {
        AutoRan result;
        qqbar_pow_si(&result.value, &value, exponent);
        return result;
    }

    // Square root
    AutoRan sqrt() const {
        int sgn = *this < 0 ? -1 : 1;
        if (sgn < 0) {
            throw std::invalid_argument("Square root of negative number");
        }

        AutoRan result;
        qqbar_sqrt(&result.value, &value);
        return result;
    }

    // n-th root
    AutoRan root(unsigned int n) const {
        AutoRan result;

        // Handle negative numbers for odd roots
        if (n % 2 == 1) { // Odd root
            int sgn = *this < 0 ? -1 : 1;
            if (sgn < 0) {
                // For negative numbers and odd roots, we can compute the root
                AutoRan abs_val;
                qqbar_abs(&abs_val.get(), &value);
                qqbar_root_ui(&result.value, &abs_val.get(), n);
                qqbar_neg(&result.value, &result.value);
                return result;
            }
        }

        // For even roots or positive numbers, use the standard approach
        qqbar_root_ui(&result.value, &value, n);
        return result;
    }

    // Absolute value
    AutoRan abs() const {
        AutoRan result;
        qqbar_abs(&result.value, &value);
        return result;
    }

    // Inverse
    AutoRan inverse() const {
        AutoRan result;
        qqbar_inv(&result.value, &value);
        return result;
    }

    // Check if zero
    bool isZero() const {
        return qqbar_is_zero(&value);
    }

private:
    Ran value;
};

// Helper functions for creating AutoRan objects
AutoRan make_rational(int num, int den = 1) {
    AutoRan result;
    fmpq_t rational;
    fmpq_init(rational);
    fmpq_set_si(rational, num, den);
    qqbar_set_fmpq(&result.get(), rational);
    fmpq_clear(rational);
    return result;
}

AutoRan make_sqrt(int n) {
    AutoRan result;
    qqbar_sqrt_ui(&result.get(), n);
    return result;
}

// Timer class for profiling
//class MyTimer {
//private:
//    std::chrono::time_point<std::chrono::high_resolution_clock> start;
//    std::string name;
//
//public:
//    MyTimer(const std::string& name) : name(name) {
//        start = std::chrono::high_resolution_clock::now();
//    }
//
//    ~MyTimer() {
//        auto end = std::chrono::high_resolution_clock::now();
//        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//        std::cout << name << " took " << duration.count() << " microseconds\n";
//    }
//};

//#define //PROFILE_FUNCTION MyTimer timer(__FUNCTION__)

// Helper macro for testing assertions
#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        std::cerr << "FAILED: " << message << std::endl; \
        failedTests++; \
    }

// Helper function to print a Ran
static void printRAN(const std::string& name, const AutoRan& ran) {
    std::cout << name << " = ";
    ran.printApprox();
    std::cout << "\n";
}

void extensiveTest() {
    //std::cout << "--- Starting Extensive Ran Tests ---" << std::endl;
    int failedTests = 0;

    // Test constructors
    /*AutoRan ranDefault;
    TEST_ASSERT(ranDefault.isZero(), "Default constructor creates zero")

        AutoRan ranIntPos(5);
    TEST_ASSERT(ranIntPos == AutoRan(5), "Integer constructor (positive)")

        AutoRan ranIntNeg(-3);
    TEST_ASSERT(ranIntNeg == AutoRan(-3), "Integer constructor (negative)")

        AutoRan ranIntZero(0);
    TEST_ASSERT(ranIntZero.isZero(), "Integer constructor (zero)")*/

    // Test arithmetic operations
    AutoRan two(2);
    AutoRan three(3);
    AutoRan minusFive(-5);
	AutoRan half = make_rational(1, 2);

    // Addition
    AutoRan sum1 = two + three;
    TEST_ASSERT(sum1 == AutoRan(5), "2 + 3 == 5")

    AutoRan sum2 = two + minusFive;
    TEST_ASSERT(sum2 == AutoRan(-3), "2 + (-5) == -3")

    AutoRan sum3 = minusFive + two;
    TEST_ASSERT(sum3 == AutoRan(-3), "(-5) + 2 == -3")

    AutoRan sumZero = two + AutoRan(-2);
    TEST_ASSERT(sumZero.isZero(), "2 + (-2) == 0")

    // Subtraction
    AutoRan diff1 = three - two;
    TEST_ASSERT(diff1 == AutoRan(1), "3 - 2 == 1")

    AutoRan diff2 = two - three;
    TEST_ASSERT(diff2 == AutoRan(-1), "2 - 3 == -1")

    AutoRan diffSelf = two - two;
    TEST_ASSERT(diffSelf.isZero(), "2 - 2 == 0")

    // Unary Negation
    AutoRan negTwo = -two;
    TEST_ASSERT(negTwo == AutoRan(-2), "-(2) == -2")

    AutoRan negMinusFive = -minusFive;
    TEST_ASSERT(negMinusFive == AutoRan(5), "-(-5) == 5")

    // Multiplication
    AutoRan prod1 = two * three;
    TEST_ASSERT(prod1 == AutoRan(6), "2 * 3 == 6")

    AutoRan prod2 = two * minusFive;
    TEST_ASSERT(prod2 == AutoRan(-10), "2 * (-5) == -10")

    AutoRan prodZero = two * AutoRan(0);
    TEST_ASSERT(prodZero.isZero(), "2 * 0 == 0")

    AutoRan prodOne = two * AutoRan(1);
    TEST_ASSERT(prodOne == two, "2 * 1 == 2")

    // Division
    AutoRan div1 = three / two;
    TEST_ASSERT(div1 == half + AutoRan(1), "3 / 2 == 1.5") // 1.5 is 1 + 0.5

    AutoRan div2 = minusFive / two;
    TEST_ASSERT(div2 == make_rational(-5, 2), "(-5) / 2 == -2.5")

    AutoRan divSelf = two / two;
    TEST_ASSERT(divSelf == AutoRan(1), "2 / 2 == 1")

    // Compound
    AutoRan aVal(10);
    AutoRan bVal(4);

    AutoRan testA = aVal;
    testA += bVal;
    TEST_ASSERT(testA == AutoRan(14), "a += b")

    testA = aVal;
    testA -= bVal;
    TEST_ASSERT(testA == AutoRan(6), "a -= b")

    testA = aVal;
    testA *= bVal;
    TEST_ASSERT(testA == AutoRan(40), "a *= b")

    testA = aVal;
    testA /= bVal;
    TEST_ASSERT(testA == make_rational(10, 4), "a /= b")

    AutoRan x(10);
    AutoRan y(10);
    AutoRan z(12);
    AutoRan negX(-10);

    TEST_ASSERT(x == y, "x == y (equal values)")
    TEST_ASSERT(x != z, "x != z (unequal values)")
    TEST_ASSERT(x < z, "x < z")
    TEST_ASSERT(z > x, "z > x")
    TEST_ASSERT(x <= y, "x <= y (equal)")
    TEST_ASSERT(x <= z, "x <= z (less than)")
    TEST_ASSERT(x >= y, "x >= y (equal)")
    TEST_ASSERT(z >= x, "z >= x (greater than)")

    TEST_ASSERT(negX < x, "-x < x")
    TEST_ASSERT(x > negX, "x > -x")
    TEST_ASSERT(negX <= x, "-x <= x")
    TEST_ASSERT(x >= negX, "x >= -x")

    AutoRan sqrt2A = make_sqrt(2); // ~1.414

    AutoRan minusSqrt2 = -make_sqrt(2); // ~-1.414
    TEST_ASSERT(sqrt2A != minusSqrt2, "sqrt(2) != -sqrt(2)")
    TEST_ASSERT(minusSqrt2 < sqrt2A, "-sqrt(2) < sqrt(2)")

    // inverse()
    AutoRan invTwo = two.inverse();
    TEST_ASSERT(invTwo == half, "inverse(2) == 0.5")

    AutoRan invHalf = half.inverse();
    TEST_ASSERT(invHalf == two, "inverse(0.5) == 2")

    AutoRan invSqrt2 = sqrt2A.inverse();
    AutoRan expectedInvSqrt2 = make_rational(1) / sqrt2A;
    TEST_ASSERT(invSqrt2 == expectedInvSqrt2, "inverse(sqrt(2)) == sqrt(2)/2")

    AutoRan doubleInvSqrt2 = invSqrt2.inverse();
    TEST_ASSERT(doubleInvSqrt2 == sqrt2A, "inverse(inverse(sqrt(2))) == sqrt(2)")

    // sqrt()
    AutoRan four(4);
    AutoRan sqrtFour = four.sqrt();
    TEST_ASSERT(sqrtFour == two, "sqrt(4) == 2")

    AutoRan sqrtTwo = two.sqrt();
    TEST_ASSERT(sqrtTwo == sqrt2A, "sqrt(2) == sqrt(2)")

    AutoRan eight(8);
    AutoRan cbrtEight = eight.root(3); // Cube root of 8
    TEST_ASSERT(cbrtEight == two, "cbrt(8) == 2")

    AutoRan negEight(-8);
    AutoRan cbrtNegEight = negEight.root(3); // Cube root of -8
    TEST_ASSERT(cbrtNegEight == AutoRan(-2), "cbrt(-8) == -2")

    // pow()
    AutoRan twoPowThree = two.pow(3);
    TEST_ASSERT(twoPowThree == AutoRan(8), "2^3 == 8")

    AutoRan twoPowZero = two.pow(0);
    TEST_ASSERT(twoPowZero == AutoRan(1), "2^0 == 1")

    AutoRan minusTwoPowThree = negTwo.pow(3);
    TEST_ASSERT(minusTwoPowThree == AutoRan(-8), "(-2)^3 == -8")

    AutoRan minusTwoPowTwo = negTwo.pow(2);
    TEST_ASSERT(minusTwoPowTwo == AutoRan(4), "(-2)^2 == 4")

    AutoRan sqrt2PowTwo = sqrt2A.pow(2);
    TEST_ASSERT(sqrt2PowTwo == two, "sqrt(2)^2 == 2")

	// zero checks
    AutoRan zeroVal(0);
    TEST_ASSERT(zeroVal.isZero(), "isZero() on 0")
    TEST_ASSERT(!two.isZero(), "isZero() on 2")
    TEST_ASSERT(!sqrt2A.isZero(), "isZero() on sqrt(2)")

    AutoRan computedZero = (sqrt2A * sqrt2A) - two;
    TEST_ASSERT(computedZero.isZero(), "(sqrt(2)*sqrt(2)) - 2 == 0")

	// chained operations
    AutoRan sqrt3Val = make_sqrt(3);

    AutoRan sumSqrt2Sqrt3 = sqrt2A + sqrt3Val;
    AutoRan expectedVal = AutoRan(5) + AutoRan(2) * AutoRan(6).sqrt();

    AutoRan chainedResult = sumSqrt2Sqrt3.pow(2);
    TEST_ASSERT(chainedResult == expectedVal, "(sqrt(2) + sqrt(3))^2 == 5 + 2*sqrt(6)")

    AutoRan cbrt2Val = make_rational(2).root(3); // ~1.2599
    AutoRan cbrt2Cubed = cbrt2Val.pow(3);
    TEST_ASSERT(cbrt2Cubed == two, "cbrt(2)^3 == 2")

    //std::cout << "--- Extensive Ran Tests Finished ---" << std::endl;
    if (failedTests > 0) {
        std::cerr << "SUMMARY: " << failedTests << " tests FAILED.\n";
    }
}

void testingFunction() {
    PROFILE_FUNCTION

    try {
        extensiveTest();

		//std::cout << "--- Starting AutoRan Tests ---" << std::endl;

        AutoRan a = 2;
        AutoRan t1 = a.sqrt();
        //std::cout << "sqrt(2) = " << t1.toDouble() << std::endl;

        AutoRan t2 = t1 + 2;
        //std::cout << "sqrt(2) + 2 = " << t2.toDouble() << std::endl;

        AutoRan t3 = t2.sqrt();
        //std::cout << "sqrt(sqrt(2) + 2) = " << t3.toDouble() << std::endl;

        AutoRan t4 = t3 + 2;
        //std::cout << "sqrt(sqrt(2) + 2) + 2 = " << t4.toDouble() << std::endl;

        AutoRan t5 = t4.sqrt();
        //std::cout << "sqrt(sqrt(sqrt(2) + 2) + 2) = " << t5.toDouble() << std::endl;

		AutoRan temp1 = t5.pow(2);
		//std::cout << "temp1 = " << temp1.toDouble() << std::endl;
		AutoRan temp2 = temp1 - 2;
		//std::cout << "temp2 = " << temp2.toDouble() << std::endl;
		AutoRan temp3 = temp2.pow(2);
		//std::cout << "temp3 = " << temp3.toDouble() << std::endl;
		AutoRan temp4 = temp3 - 2;
		//std::cout << "temp4 = " << temp4.toDouble() << std::endl;
		AutoRan temp5 = temp4.pow(2);
		//std::cout << "temp5 = " << temp5.toDouble() << std::endl;

		if (a != temp5)
			std::cout << "Error: a != temp5\n";
		if (t1 != temp4)
			std::cout << "Error: t1 != temp4\n";
		if (t2 != temp3)
			std::cout << "Error: t2 != temp3\n";
		if (t3 != temp2)
			std::cout << "Error: t3 != temp2\n";
		if (t4 != temp1)
			std::cout << "Error: t4 != temp1\n";


        // Test complex expressions
        auto aa = (make_rational(3504, 100)) - 35;
        auto bb = make_rational(9) / make_sqrt(93) + (make_rational(5654, 100));
        auto cc = (make_rational(192, 10)) / 96;
        auto dd = make_rational(6586, 100) - 27 - make_rational(7305, 100) + make_rational(2797, 100) / make_rational(818, 10).pow(2);
		auto ee = make_rational(9758, 100) / make_rational(6145, 100);
		auto ff = make_rational(28).pow(2) * 78 + make_sqrt(42);
        auto gg = make_rational(2) * 46 + make_rational(3).pow(2);
		auto hh = make_rational(5431, 100) - make_rational(32) / make_rational(83).pow(4) / make_rational(94).pow(3);
		auto ii = make_rational(55) - 32 + 76 - make_rational(597, 10);
		auto jj = make_rational(2661, 100).pow(2) + (make_rational(1878, 100) / 48);

        /*std::cout << "aa = ";
    	aa.printApprox();
    	std::cout << std::endl;
		std::cout << "bb = ";
		bb.printApprox();
		std::cout << std::endl;
		std::cout << "cc = ";
		cc.printApprox();
		std::cout << std::endl;
		std::cout << "dd = ";
		dd.printApprox();
		std::cout << std::endl;
		std::cout << "ee = ";
		ee.printApprox();
		std::cout << std::endl;
		std::cout << "ff = ";
		ff.printApprox();
		std::cout << std::endl;
		std::cout << "gg = ";
		gg.printApprox();
		std::cout << std::endl;
		std::cout << "hh = ";
		hh.printApprox();
		std::cout << std::endl;
		std::cout << "ii = ";
		ii.printApprox();
		std::cout << std::endl;
		std::cout << "jj = ";
		jj.printApprox();
		std::cout << std::endl;*/
    }
    catch (const std::exception& e) {
        std::cout << e.what() << '\n';
    }
}

int main() {
    // Initialize FLINT
    flint_rand_t state;
    flint_randinit(state);

    InitializePerformanceFrequency();
    for (int i = 0; i < 200; ++i) {
        std::cout << i;
        testingFunction();
    }
    atexit(PrintProfilingReport);

    // Clean up FLINT
    flint_randclear(state);
    flint_cleanup();

    return 0;
}