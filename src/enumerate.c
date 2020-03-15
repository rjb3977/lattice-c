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

static void search(long dimensions, long depth, mpq_t transform[dimensions][2 * dimensions], mpq_t offset[dimensions], mpq_t fixed[dimensions], mpq_t A[2 * dimensions][3 * dimensions], mpq_t b[2 * dimensions], mpq_t c[2][3 * dimensions], mpq_t λ[2 * dimensions], mpq_t s[dimensions], mpq_t d[2 * dimensions], mpq_t x[2 * dimensions], long B[2 * dimensions], long N[dimensions], mpq_t lu[2 * dimensions][2 * dimensions], long pivots[2 * dimensions], long *corrections, long **indices, mpq_t(**columns), long *count_out, mpq_t (**results_out)[dimensions]) {
    if (depth == dimensions) {
        printf("found: ");
        vec_print(dimensions, fixed, stdout);
        printf("\n");

        *count_out += 1;
        *results_out = realloc(*results_out, *count_out * sizeof(mpq_t[dimensions]));

        vec_init(dimensions, (*results_out)[*count_out - 1]);
        vec_set(dimensions, (*results_out)[*count_out - 1], fixed);
    } else {
        mpq_t temp;
        mpz_t min;
        mpz_t max;

        mpq_init(temp);
        mpz_init(min);
        mpz_init(max);

        // lower bound
        vec_neg(dimensions, c[0], transform[depth]);
        simplex_solve(dimensions, depth, A, b, c, λ, s, d, x, B, N, lu, pivots, corrections, indices, (void*) columns);

        vec_dot(2 * dimensions, temp, transform[depth], x);
        mpq_sub(temp, offset[depth], temp);
        mpz_cdiv_q(min, mpq_numref(temp), mpq_denref(temp));

        // upper bound
        vec_set(dimensions, c[0], transform[depth]);
        simplex_solve(dimensions, depth, A, b, c, λ, s, d, x, B, N, lu, pivots, corrections, indices, (void*) columns);

        vec_dot(2 * dimensions, temp, transform[depth], x);
        mpq_sub(temp, offset[depth], temp);
        mpz_fdiv_q(max, mpq_numref(temp), mpq_denref(temp));

        // rhs of new row of A
        mpq_set_z(temp, min);
        mpq_sub(temp, offset[depth], temp);

        // printf("%ld: ", depth);
        // mpz_out_str(stdout, 10, min);
        // printf(" -> ");
        // mpz_out_str(stdout, 10, max);
        // printf("\n");

        // min <= max, min += 1, rhs -= 1
        for (; mpz_cmp(min, max) <= 0; mpz_add_ui(min, min, 1), mpz_sub(mpq_numref(temp), mpq_numref(temp), mpq_denref(temp))) {
            if (mpq_sgn(temp) >= 0) {
                vec_set(dimensions, A[dimensions + depth], transform[depth]);
                mpq_set(b[dimensions + depth], temp);
            } else {
                vec_neg(dimensions, A[dimensions + depth], transform[depth]);
                mpq_neg(b[dimensions + depth], temp);
            }

            mpq_set_z(fixed[depth], min);
            search(dimensions, depth + 1, transform, offset, fixed, A, b, c, λ, s, d, x, B, N, lu, pivots, corrections, indices, columns, count_out, results_out);
        }

        mpq_clear(temp);
        mpz_clear(min);
        mpz_clear(max);
    }
}

void enumerate(long dimensions, mpq_t basis[dimensions][dimensions], mpq_t lower[dimensions], mpq_t upper[dimensions], long* count_out, mpq_t (**results_out)[dimensions]) {
    mpq_t transform[dimensions][2 * dimensions];
    mpq_t offset[dimensions];
    mpq_t fixed[dimensions];

    mpq_t A[2 * dimensions][3 * dimensions];
    mpq_t b[2 * dimensions];
    mpq_t c[2][3 * dimensions];
    mpq_t λ[2 * dimensions];
    mpq_t s[dimensions];
    mpq_t d[2 * dimensions];
    mpq_t x[2 * dimensions];
    long B[2 * dimensions];
    long N[dimensions];
    mpq_t lu[2 * dimensions][2 * dimensions];
    long pivots[2 * dimensions];
    long corrections = 0;
    long *indices = 0;
    mpq_t *columns = NULL;

    mat_init(dimensions, 2 * dimensions, transform);
    vec_init(dimensions, offset);
    vec_init(dimensions, fixed);

    mat_init(2 * dimensions, 3 * dimensions, A);
    vec_init(2 * dimensions, b);
    mat_init(2, 3 * dimensions, c);
    vec_init(2 * dimensions, λ);
    vec_init(dimensions, s);
    vec_init(2 * dimensions, d);
    vec_init(2 * dimensions, x);
    mat_init(2 * dimensions, 2 * dimensions, lu);

    mat_copy(dimensions, dimensions, 0, 0, 2 * dimensions, lu, 0, 0, dimensions, basis);
    mat_lu_2(dimensions, 2 * dimensions, lu, pivots);

    // inverse * identity
    for (long i = 0; i < dimensions; ++i) {
        mpq_set_ui(transform[i][i], 1, 1);
        solve_ptlu_2(dimensions, 2 * dimensions, lu, pivots, transform[i]);
    }

    // inverse * upper
    vec_set(dimensions, offset, upper);
    solve_utltp_2(dimensions, 2 * dimensions, lu, pivots, offset);

    // A = [I|I|0]
    //     [0|0|I]
    // c = [.|0|0]
    //     [0|0|1]
    for (long i = 0; i < dimensions; ++i) {
        mpq_set_ui(A[i][i], 1, 1);
        mpq_set_ui(A[i][dimensions + i], 1, 1);
        mpq_set_ui(A[dimensions + i][2 * dimensions + i], 1, 1);
        mpq_set_ui(c[1][2 * dimensions + i], 1, 1);
    }

    vec_sub(dimensions, b, upper, lower);

    search(dimensions, 0, transform, offset, fixed, A, b, c, λ, s, d, x, B, N, lu, pivots, &corrections, &indices, &columns, count_out, results_out);

    mat_clear(dimensions, 2 * dimensions, transform);
    vec_clear(dimensions, offset);
    vec_clear(dimensions, fixed);

    mat_clear(2 * dimensions, 3 * dimensions, A);
    vec_clear(2 * dimensions, b);
    mat_clear(2, 3 * dimensions, c);
    vec_clear(2 * dimensions, λ);
    vec_clear(dimensions, s);
    vec_clear(2 * dimensions, d);
    vec_clear(2 * dimensions, x);
    mat_clear(2 * dimensions, 2 * dimensions, lu);
}
