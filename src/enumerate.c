#define _POSIX_C_SOURCE 200809L

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <gmp.h>

#include "la.h"
#include "lp.h"
#include "enumerate.h"

static void search(long dimensions, mpq_t transform[dimensions][2 * dimensions], mpq_t transformOffset[dimensions], mpq_t A[2 * dimensions][2 * dimensions], mpq_t b[2 * dimensions], mpq_t x[2 * dimensions], mpq_t fixed[dimensions], long depth, long* count_out, mpq_t (**results_out)[dimensions]) {
    if (depth == dimensions) {
        printf("found: ");
        vec_print(dimensions, fixed, stdout);
        printf("\n");

        *count_out += 1;
        *results_out = realloc(*results_out, *count_out * sizeof(mpq_t[dimensions]));

        vec_init(dimensions, (*results_out)[*count_out - 1]);
        vec_copy(dimensions, 0, (*results_out)[*count_out - 1], 0, fixed);
    } else {
        mpq_t temp;
        mpz_t min;
        mpz_t max;

        mpq_init(temp);
        mpz_init(min);
        mpz_init(max);

        // lower bound
        vec_neg(2 * dimensions, transform[depth], transform[depth]);
        simplex_solve(dimensions + depth, 2 * dimensions, A, b, transform[depth], x);
        vec_neg(2 * dimensions, transform[depth], transform[depth]);

        vec_dot(2 * dimensions, temp, transform[depth], x);
        mpq_sub(temp, transformOffset[depth], temp);
        mpz_cdiv_q(min, mpq_numref(temp), mpq_denref(temp));

        // upper bound
        simplex_solve(dimensions + depth, 2 * dimensions, A, b, transform[depth], x);

        vec_dot(2 * dimensions, temp, transform[depth], x);
        mpq_sub(temp, transformOffset[depth], temp);
        mpz_fdiv_q(max, mpq_numref(temp), mpq_denref(temp));

        // rhs of new row of A
        mpq_set_z(temp, min);
        mpq_sub(temp, transformOffset[depth], temp);

        printf("%ld: ", depth);
        mpz_out_str(stdout, 10, min);
        printf(" -> ");
        mpz_out_str(stdout, 10, max);
        printf("\n");

        // min <= max, min += 1, rhs -= 1
        for (; mpz_cmp(min, max) <= 0; mpz_add_ui(min, min, 1), mpz_sub(mpq_numref(temp), mpq_numref(temp), mpq_denref(temp))) {
            if (mpq_sgn(temp) >= 0) {
                vec_copy(dimensions, 0, A[dimensions + depth], 0, transform[depth]);
                mpq_set(b[dimensions + depth], temp);
            } else {
                vec_neg(dimensions, A[dimensions + depth], transform[depth]);
                mpq_neg(b[dimensions + depth], temp);
            }

            mpq_set_z(fixed[depth], min);
            search(dimensions, transform, transformOffset, A, b, x, fixed, depth + 1, count_out, results_out);
        }

        mpq_clear(temp);
        mpz_clear(min);
        mpz_clear(max);
    }
}

void enumerate(long dimensions, mpq_t basis[dimensions][dimensions], mpq_t lower[dimensions], mpq_t upper[dimensions], long* count_out, mpq_t (**results_out)[dimensions]) {
    mpq_t transform[dimensions][2 * dimensions];
    mpq_t transformOffset[dimensions];
    mpq_t A[2 * dimensions][2 * dimensions];
    mpq_t b[2 * dimensions];
    mpq_t x[2 * dimensions];
    mpq_t fixed[dimensions];
    mpq_t lu[dimensions][dimensions];
    long pivots[dimensions];

    mat_init(dimensions, 2 * dimensions, transform);
    vec_init(dimensions, transformOffset);
    mat_init(2 * dimensions, 2 * dimensions, A);
    vec_init(2 * dimensions, b);
    vec_init(2 * dimensions, x);
    vec_init(dimensions, fixed);
    mat_init(dimensions, dimensions, lu);

    mat_copy(dimensions, dimensions, 0, 0, dimensions, lu, 0, 0, dimensions, basis);
    mat_lu(dimensions, lu, pivots);

    // inverse * identity
    for (long i = 0; i < dimensions; ++i) {
        mpq_set_ui(transform[i][i], 1, 1);
        solve_ptlu(dimensions, lu, pivots, transform[i]);
    }

    // inverse * upper
    vec_copy(dimensions, 0, transformOffset, 0, upper);
    solve_utltp(dimensions, lu, pivots, transformOffset);

    // A = [I|I]
    for (long i = 0; i < dimensions; ++i) {
        mpq_set_ui(A[i][i], 1, 1);
        mpq_set_ui(A[i][dimensions + i], 1, 1);
    }

    vec_sub(dimensions, b, upper, lower);

    search(dimensions, transform, transformOffset, A, b, x, fixed, 0, count_out, results_out);
}
