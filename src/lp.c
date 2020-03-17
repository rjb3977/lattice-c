#define _POSIX_C_SOURCE 200809L

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <gmp.h>

#include "la.h"
#include "lp.h"

static bool simplex_step(long rows, long cols, long stride, mpq_t A[rows][stride], mpq_t b[rows], mpq_t c[cols], mpq_t λ[rows], mpq_t s[cols - rows], mpq_t d[rows], mpq_t x[rows], long B[rows], long N[cols - rows], long lu_stride, mpq_t lu[rows][lu_stride], long pivots[rows], long *corrections, long **indices, mpq_t (**columns)[rows]) {
    mpq_t temp;
    mpq_init(temp);

    for (long i = 0; i < rows; ++i) {
        mpq_set(λ[i], c[B[i]]);
    }

    solve_utltp_corrected_2(rows, lu_stride, lu, pivots, *corrections, *indices, *columns, λ);

    for (long i = 0; i < cols - rows; ++i) {
        mpq_set_ui(s[i], 0, 1);

        for (long j = 0; j < rows; ++j) {
            mpq_mul(temp, λ[j], A[j][N[i]]);
            mpq_add(s[i], s[i], temp);
        }

        mpq_sub(s[i], c[N[i]], s[i]);
    }

    // select entering, make this better maybe?

    long q = -1;

    for (long i = 0; i < cols - rows; ++i) {
        if (mpq_sgn(s[i]) == -1 && (q == -1 || mpq_cmp(s[i], s[q]) < 0)) {
            q = i;
        }
    }

    if (q == -1) {
        mpq_clear(temp);
        return true;
    }

    for (long i = 0; i < rows; ++i) {
        mpq_set(d[i], A[i][N[q]]);
    }

    solve_ptlu_corrected_2(rows, lu_stride, lu, pivots, *corrections, *indices, *columns, d);

    // select exiting, this is straightforward and optimal

    long p = -1;
    mpq_t xq;
    mpq_init(xq);

    for (long i = 0; i < rows; ++i) {
        if (mpq_sgn(d[i]) == 1) {
            mpq_div(temp, x[i], d[i]);

            if (p == -1 || mpq_cmp(temp, xq) < 0) {
                mpq_set(xq, temp);
                p = i;
            }
        }
    }

    assert(p != -1 && "solution is unbounded");

    for (long i = 0; i < rows; ++i) {
        mpq_mul(temp, xq, d[i]);
        mpq_sub(x[i], x[i], temp);
    }

    mpq_set(x[p], xq);

    long Bp = B[p];
    long Nq = N[q];

    B[p] = Nq;
    N[q] = Bp;

    *corrections += 1;
    *indices = realloc(*indices, *corrections * sizeof(long));
    *columns = realloc(*columns, *corrections * sizeof(mpq_t[rows]));

    (*indices)[*corrections - 1] = p;

    vec_init(rows, (*columns)[*corrections - 1]);
    vec_set(rows, (*columns)[*corrections - 1], d);

    mpq_clear(xq);
    mpq_clear(temp);

    return false;
}

static void simplex_init(long size, long depth, mpq_t A[2 * size][3 * size], mpq_t b[2 * size], mpq_t c[2][3 * size], mpq_t λ[2 * size], mpq_t s[size], mpq_t d[2 * size], mpq_t x[2 * size], long B[2 * size], long N[size], mpq_t lu[2 * size][2 * size], long pivots[2 * size], long *corrections, long **indices, mpq_t (**columns)[size + depth]) {
    // printf("A\n");
    // mat_print(2 * size, 3 * size, A, stdout);
    // printf("\n\n");

    // printf("b\n");
    // mat_print(2 * size, 1, (void*) b, stdout);
    // printf("\n\n");

    for (long i = 0; i < size; ++i) {
        N[i] = i;
    }

    for (long i = 0; i < size + depth; ++i) {
        B[i] = size + i;
    }

    for (long i = 0; i < size + depth; ++i) {
        for (long j = 0; j < size + depth; ++j) {
            if (i == j) {
                mpq_set_ui(lu[i][j], 1, 1);
            } else {
                mpq_set_ui(lu[i][j], 0, 1);
            }
        }
    }

    for (long i = 0; i < size + depth; ++i) {
        pivots[i] = i;
    }

    vec_set(size + depth, x, b);

    while (!simplex_step(size + depth, 2 * size + depth, 3 * size, A, b, c[1], λ, s, d, x, B, N, 2 * size, lu, pivots, corrections, indices, columns));

    // Swap all artificial basic variables with real variables.
    // In effect, we're exiting the artificial variables and
    // entering the real variables along an edge of length zero.
    for (long i = 0; i < size + depth; ++i) {
        if (B[i] >= 2 * size) {
            assert(mpq_sgn(x[i]) == 0 && "artificial variable is nonzero");

            for (long j = 0; j < 2 * size; ++j) {
                if (N[j] < 2 * size && mpq_sgn(A[i][N[j]]) != 0) {
                    long temp = B[i];
                    B[i] = N[j];
                    N[j] = temp;
                }
            }
        }
    }

    // Find all real nonbasic variables.
    for (long i = 0, j = 0; i < size - depth; ++i, ++j) {
        for (;; ++j) {
            if (N[j] < 2 * size) {
                long temp = N[i];
                N[i] = N[j];
                N[j] = temp;
                break;
            }
        }
    }

    // Clear LU updates.
    mat_clear(*corrections, size + depth, *columns);
    free(*indices);
    free(*columns);

    *corrections = 0;
    *indices = NULL;
    *columns = NULL;

    // Reset basis LU.
    for (long row = 0; row < size + depth; ++row) {
        for (long col = 0; col < size + depth; ++col) {
            mpq_set(lu[row][col], A[row][B[col]]);
        }
    }

    mat_lu_2(size + depth, 2 * size, lu, pivots);
}

void simplex_solve(long size, long depth, mpq_t A[2 * size][3 * size], mpq_t b[2 * size], mpq_t c[2][3 * size], mpq_t λ[2 * size], mpq_t s[size], mpq_t d[2 * size], mpq_t x[2 * size], long B[2 * size], long N[size], mpq_t lu[2 * size][2 * size], long pivots[2 * size], long *corrections, long **indices, mpq_t (**columns)[2 * size]) {
    simplex_init(size, depth, A, b, c, λ, s, d, x, B, N, lu, pivots, corrections, indices, columns);

    while (!simplex_step(size + depth, 2 * size, 3 * size, A, b, c[0], λ, s, d, x, B, N, 2 * size, lu, pivots, corrections, indices, columns));

    // Set output
    for (long i = 0; i < 2 * size; ++i) {
        mpq_set_ui(d[i], 0, 1);
    }

    for (long i = 0; i < size + depth; ++i) {
        mpq_set(d[B[i]], x[i]);
    }

    vec_set(2 * size, x, d);

    // Clear LU updates
    mat_clear(*corrections, size + depth, *columns);
    free(*indices);
    free(*columns);

    *corrections = 0;
    *indices = NULL;
    *columns = NULL;
}
