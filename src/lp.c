#define _POSIX_C_SOURCE 200809L

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <gmp.h>

#include "la.h"
#include "lp.h"

void simplex_solve(long rows, long cols, mpq_t A[rows][cols], mpq_t b[rows], mpq_t c[cols], mpq_t x[cols]) {
    // declare necessary variables
    mpq_t λ[rows];
    mpq_t s[cols - rows];
    mpq_t d[rows];
    mpq_t _x[rows];

    long B[rows];
    long N[cols - rows];

    mpq_t lu[rows][cols];
    long pivots[rows];

    long corrections;
    long *indices;
    mpq_t (*columns)[rows];

    // initialize necessary variables
    vec_init(rows, λ);
    vec_init(cols - rows, s);
    vec_init(rows, d);
    vec_init(rows, _x);

    mat_init(rows, cols, lu);

    corrections = 0;
    indices = NULL;
    columns = NULL;

    // solve system
    simplex_init(rows, cols, A, b, c, λ, s, d, _x, B, N, lu, pivots, &corrections, &indices, &columns);

    long k = 1;
    while (!simplex_step(rows, cols, A, b, c, λ, s, d, _x, B, N, lu, pivots, &corrections, &indices, &columns)) ++k;

    printf("solve iterations:     %ld\n", k);

    // set output
    for (long i = 0; i < cols; ++i) {
        mpq_set_ui(x[i], 0, 1);
    }

    for (long i = 0; i < rows; ++i) {
        mpq_set(x[B[i]], _x[i]);
    }

    // clean up variables
    vec_clear(rows, λ);
    vec_clear(cols - rows, s);
    vec_clear(rows, d);
    vec_clear(rows, _x);

    mat_clear(rows, cols, lu);

    free(indices);
    free(columns);
}

void simplex_init(long rows, long cols, mpq_t A[rows][cols], mpq_t b[rows], mpq_t c[cols], mpq_t λ[rows], mpq_t s[cols - rows], mpq_t d[rows], mpq_t x[rows], long B[rows], long N[cols - rows], mpq_t lu[rows][rows], long pivots[rows], long *corrections, long **indices, mpq_t (**columns)[rows]) {
    // Set up variables for artificial system.
    mpq_t _A[rows][cols + rows];
    mpq_t _c[cols + rows];
    mpq_t _s[cols];
    long _N[cols];

    mat_init(rows, cols + rows, _A);
    vec_init(cols + rows, _c);
    vec_init(cols, _s);

    mat_copy(rows, cols, 0, 0, cols + rows, _A, 0, 0, cols, A);
    vec_copy(rows, 0, x, 0, b);

    // Set up identity sub matrix, gradient, and initial basic variables.
    for (long i = 0; i < rows; ++i) {
        mpq_set_ui(_A[i][cols + i], 1, 1);
        mpq_set_ui(_c[cols + i], 1, 1);
        B[i] = cols + i;
    }

    // Set up initial nonbasic variables.
    for (long i = 0; i < cols; ++i) {
        _N[i] = i;
    }

    // Set up initial basis LU.
    for (long row = 0; row < rows; ++row) {
        for (long col = 0; col < rows; ++col) {
            if (row == col) {
                mpq_set_ui(lu[row][col], 1, 1);
            } else {
                mpq_set_ui(lu[row][col], 0, 1);
            }
        }

        pivots[row] = row;
    }

    // Solve artificial system.
    long k = 1;
    while (!simplex_step(rows, cols + rows, _A, b, _c, λ, _s, d, x, B, _N, lu, pivots, corrections, indices, columns)) ++k;

    printf("pre-solve iterations: %ld\n", k);

    // Swap all artificial basic variables with real variables.
    // In effect, we're exiting the artificial variables and
    // entering the real variables along an edge of length zero.
    for (long i = 0; i < rows; ++i) {
        if (B[i] >= cols) {
            assert(mpq_sgn(x[i]) == 0 && "no valid solution");

            for (long j = 0; j < cols; ++j) {
                if (_N[j] < cols && mpq_sgn(_A[i][_N[j]]) != 0) {
                    long temp = B[i];
                    B[i] = _N[j];
                    _N[j] = temp;
                }
            }
        }
    }

    // Find all real nonbasic variables.
    for (long i = 0, j = 0; i < cols - rows; ++i, ++j) {
        for (;; ++j) {
            if (_N[j] < cols) {
                N[i] = _N[j];
                break;
            }
        }
    }

    // Clear LU updates.
    mat_clear(*corrections, rows, *columns);
    free(*indices);
    free(*columns);

    *corrections = 0;
    *indices = NULL;
    *columns = NULL;

    // Reset basis LU.
    for (long row = 0; row < rows; ++row) {
        for (long col = 0; col < rows; ++col) {
            mpq_set(lu[row][col], A[row][B[col]]);
        }
    }

    mat_lu(rows, lu, pivots);

    mat_clear(rows, cols + rows, _A);
    vec_clear(cols + rows, _c);
    vec_clear(cols, _s);
}

bool simplex_step(long rows, long cols, mpq_t A[rows][cols], mpq_t b[rows], mpq_t c[cols], mpq_t λ[rows], mpq_t s[cols - rows], mpq_t d[rows], mpq_t x[rows], long B[rows], long N[cols - rows], mpq_t lu[rows][rows], long pivots[rows], long *corrections, long **indices, mpq_t (**columns)[rows]) {
    assert(rows <= cols && "solution is overconstrained");

    mpq_t temp;
    mpq_init(temp);

    for (long i = 0; i < rows; ++i) {
        mpq_set(λ[i], c[B[i]]);
    }

    solve_utltp_corrected(rows, lu, pivots, *corrections, *indices, *columns, λ);

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
        return true;
    }

    for (long i = 0; i < rows; ++i) {
        mpq_set(d[i], A[i][N[q]]);
    }

    solve_ptlu_corrected(rows, lu, pivots, *corrections, *indices, *columns, d);

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

    assert(mpq_sgn(x[p]) == 0 && "exiting variable did not go to zero");

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
    vec_copy(rows, 0, (*columns)[*corrections - 1], 0, d);

    mpq_clear(xq);
    mpq_clear(temp);
    return false;
}
