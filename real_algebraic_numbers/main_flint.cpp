#include <iostream>
#include <string>
#include <chrono>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "qqbar.h"
#include "acb.h"
#include "acb_poly.h"

void qqbar_set_algebraic_from_poly(qqbar_t rop, const fmpz_poly_t poly, double approx) {
    if (fmpz_poly_degree(poly) <= 0) {
        // Handle constant polynomials
        qqbar_set_si(rop, fmpz_poly_get_coeff_si(poly, 0));
        return;
    }

    // Use Arb to find roots and their isolating intervals
    acb_poly_t acb_poly;
    acb_poly_init(acb_poly);
    acb_poly_set_fmpz_poly(acb_poly, poly, FLINT_DEFAULT_PREC);

    acb_struct* roots;
    long num_roots = acb_poly_find_roots(acb_poly, &roots, 0);

    slong selected_root_index = -1;
    double min_dist_sq = -1.0;

    // Iterate through the roots to find the one closest to our approximation
    for (slong i = 0; i < num_roots; ++i) {
        acb_t root_i;
        acb_init(root_i);
        acb_set(root_i, roots + i);

        double real_part = arf_get_d(acb_realref(root_i), ARF_RND_NEAR);
        double imag_part = arf_get_d(acb_imagref(root_i), ARF_RND_NEAR);

        // Check if the root is a real root
        if (fabs(imag_part) < 1e-9) {
            double dist_sq = (real_part - approx) * (real_part - approx);
            if (selected_root_index == -1 || dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
                selected_root_index = i;
            }
        }
        acb_clear(root_i);
    }

    if (selected_root_index != -1) {
        // We found the closest real root, now get its isolating interval
        acb_t selected_root;
        acb_init(selected_root);
        acb_set(selected_root, roots + selected_root_index);

        arb_t lo, hi;
        arb_init(lo);
        arb_init(hi);

        // The real part of the root provides the interval's midpoint
        arb_get_interval_arf(lo, hi, acb_realref(selected_root), FLINT_DEFAULT_PREC);

        // Pass the polynomial and the isolating interval to the core FLINT function
        qqbar_set_algebraic(rop, poly, lo, hi);

        arb_clear(lo);
        arb_clear(hi);
        acb_clear(selected_root);
    }

    acb_poly_clear(acb_poly);
    for (slong i = 0; i < num_roots; ++i) {
        acb_clear(roots + i);
    }
    flint_free(roots);
}

int main() {
    // Initialize FLINT's memory system
    flint_cleanup();
    //flint_init();

    // ----------------------------------------------------------------------
    // Define the algebraic numbers
    // ----------------------------------------------------------------------

    // The Golden Ratio: phi = (1 + sqrt(5))/2
    // Minimal polynomial: x^2 - x - 1 = 0
    qqbar_t phi;
    qqbar_init(phi);
    fmpz_poly_t poly_phi;
    fmpz_poly_init(poly_phi);
    fmpz_poly_set_coeff_si(poly_phi, 2, 1);
    fmpz_poly_set_coeff_si(poly_phi, 1, -1);
    fmpz_poly_set_coeff_si(poly_phi, 0, -1);
    qqbar_set_algebraic_from_poly(phi, poly_phi, 1.618); // Use a double to select the correct root
    fmpz_poly_clear(poly_phi);

    // The cube root of 2: alpha = 2^(1/3)
    // Minimal polynomial: x^3 - 2 = 0
    qqbar_t alpha;
    qqbar_init(alpha);
    fmpz_poly_t poly_alpha;
    fmpz_poly_init(poly_alpha);
    fmpz_poly_set_coeff_si(poly_alpha, 3, 1);
    fmpz_poly_set_coeff_si(poly_alpha, 0, -2);
    qqbar_set_algebraic_from_poly(alpha, poly_alpha, 1.259); // Use a double to select the correct root
    fmpz_poly_clear(poly_alpha);

    // The square root of 3
    // Minimal polynomial: x^2 - 3 = 0
    qqbar_t sqrt3;
    qqbar_init(sqrt3);
    fmpz_poly_t poly_sqrt3;
    fmpz_poly_init(poly_sqrt3);
    fmpz_poly_set_coeff_si(poly_sqrt3, 2, 1);
    fmpz_poly_set_coeff_si(poly_sqrt3, 0, -3);
    qqbar_set_algebraic_from_poly(sqrt3, poly_sqrt3, 1.732);
    fmpz_poly_clear(poly_sqrt3);

    // A rational number, 1/2
    qqbar_t half;
    qqbar_init(half);
    qqbar_set_d(half, 0.5);

    // ----------------------------------------------------------------------
    // Start the clock for profiling
    // ----------------------------------------------------------------------
    auto start = std::chrono::high_resolution_clock::now();

    // ----------------------------------------------------------------------
    // Perform complex calculations
    // ----------------------------------------------------------------------

    qqbar_t temp1, temp2, temp3, temp4, temp5, result;
    qqbar_init(temp1);
    qqbar_init(temp2);
    qqbar_init(temp3);
    qqbar_init(temp4);
    qqbar_init(temp5);
    qqbar_init(result);

    // temp1 = phi + alpha
    qqbar_add(temp1, phi, alpha);

    // temp2 = sqrt3 * half
    qqbar_mul(temp2, sqrt3, half);

    // temp3 = temp1 - temp2
    qqbar_sub(temp3, temp1, temp2); // (phi + alpha) - (sqrt3 * 1/2)

    // temp4 = temp3 * phi
    qqbar_mul(temp4, temp3, phi); // ((phi + alpha) - (sqrt3 * 1/2)) * phi

    // temp5 = alpha / sqrt3
    qqbar_div(temp5, alpha, sqrt3);

    // result = temp4 + temp5
    qqbar_add(result, temp4, temp5); // (((phi + alpha) - (sqrt3 * 1/2)) * phi) + (alpha / sqrt3)

    // A more complex operation: sqrt(result + phi)
    qqbar_add(result, result, phi);
    qqbar_sqrt(result, result);

    // ----------------------------------------------------------------------
    // Stop the clock
    // ----------------------------------------------------------------------
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    // ----------------------------------------------------------------------
    // Print the results and timing
    // ----------------------------------------------------------------------
    std::cout << "Final result (as minimal polynomial): " << std::endl;
    fmpz_poly_print(qqbar_poly(result));
    std::cout << std::endl;

    std::cout << "Approximate value: " << qqbar_get_d(result) << std::endl;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;

    // ----------------------------------------------------------------------
    // Clean up
    // ----------------------------------------------------------------------
    qqbar_clear(phi);
    qqbar_clear(alpha);
    qqbar_clear(sqrt3);
    qqbar_clear(half);
    qqbar_clear(temp1);
    qqbar_clear(temp2);
    qqbar_clear(temp3);
    qqbar_clear(temp4);
    qqbar_clear(temp5);
    qqbar_clear(result);

    flint_cleanup();

    return 0;
}