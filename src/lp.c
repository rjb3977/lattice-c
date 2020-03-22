#define _POSIX_C_SOURCE 200809L

#include <assert.h>
#include <gmp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "la.h"

static void lp_pivot(matrix_t* table, long* B, long* N, long variables, long constraints, long entering, long exiting) {
    const long a = matrix_rows(table) - constraints;
    const long b = matrix_cols(table) - 1;

    assert(0 <= entering && entering < b);
    assert(0 <= exiting && exiting < constraints);

    mpq_t t0;
    mpq_init(t0);

    mpq_t t1;
    mpq_init(t1);

    for (long col = 0; col < matrix_cols(table); ++col) {
        if (col == entering) {
            continue;
        }

        mpq_ptr x = matrix_at(table, a + exiting, col);
        mpq_ptr y = matrix_at(table, a + exiting, entering);

        mpq_div(x, x, y);
    }

    for (long row = 0; row < matrix_rows(table); ++row) {
        if (row == a + exiting) {
            continue;
        }

        mpq_ptr x = matrix_at(table, row, entering);

        for (long col = 0; col < matrix_cols(table); ++col) {
            if (col == entering) {
                continue;
            }

            mpq_ptr y = matrix_at(table, row, col);

            mpq_mul(t0, x, matrix_at(table, a + exiting, col));
            mpq_sub(y, y, t0);
        }

        mpq_div(x, x, matrix_at(table, a + exiting, entering));
        mpq_neg(x, x);
    }

    mpq_inv(matrix_at(table, a + exiting, entering), matrix_at(table, a + exiting, entering));

    long _entering = N[entering];
    long _exiting = B[exiting];

    N[entering] = _exiting;
    B[exiting] = _entering;

    mpq_clear(t0);
    mpq_clear(t1);
}

// B: row    -> variable
// N: column -> variable
static bool lp_step(matrix_t* table, long* B, long* N, long variables, long constraints) {
    const long a = matrix_rows(table) - constraints; // first row of A
    const long b = matrix_cols(table) - 1;           // first col of b

    mpq_t t0;
    mpq_init(t0);

    mpq_t t1;
    mpq_init(t1);

    bool bland = false;

    for (long row = 0; row < constraints; ++row) {
        if (mpq_sgn(matrix_at(table, a + row, b)) == 0) {
            bland = true;
            break;
        }
    }

    // column, [0, b)
    long entering = -1;

    for (long col = 0; col < b; ++col) {
        mpq_ptr x = matrix_at(table, 0, col);

        if (mpq_sgn(x) > 0 && (entering == -1 || mpq_cmp(x, t0) > 0)) {
            entering = col;
            mpq_set(t0, x);

            if (bland) {
                break;
            }
        }
    }

    if (entering == -1) {
        mpq_clear(t0);
        mpq_clear(t1);

        return true;
    }

    // row, [0, constraints]
    long exiting = -1;

    for (long row = 0; row < constraints; ++row) {
        mpq_ptr x = matrix_at(table, a + row, entering);
        mpq_ptr y = matrix_at(table, a + row, b);

        if (mpq_sgn(x) > 0) {
            mpq_div(t1, y, x);

            if (exiting == -1 || mpq_cmp(t1, t0) < 0) {
                exiting = row;
                mpq_set(t0, t1);
            }
        }
    }

    assert(exiting != -1);

    mpq_clear(t0);
    mpq_clear(t1);

    // printf("entering, exiting: %ld, %ld\n", entering, exiting);

    lp_pivot(table, B, N, variables, constraints, entering, exiting);

    return false;
}

void lp_solve(matrix_t* dest, const matrix_t* initial_table, long dimensions, long depth) {
    assert(matrix_rows(dest) == dimensions);
    assert(matrix_cols(dest) == 1);

    matrix_t* table = matrix_alloc(2 + dimensions + depth, dimensions + 1);

    // row [1, 2 + size + depth) (Z and A)
    for (long row = 1; row < 2 + dimensions + depth; ++row) {
        for (long col = 0; col < dimensions + 1; ++col) {
            mpq_set(matrix_at(table, row, col), matrix_cat(initial_table, row - 1, col));
        }
    }

    // row 0 (for finding initial basis)
    for (long row = dimensions; row < dimensions + depth; ++row) {
        for (long col = 0; col < dimensions + 1; ++col) {
            mpq_add(matrix_at(table, 0, col), matrix_at(table, 0, col), matrix_at(table, 2 + row, col));
        }
    }

    long B[dimensions + depth];
    long N[dimensions];

    for (long i = 0; i < 2 * dimensions + depth; ++i) {
        if (i < dimensions) {
            N[i] = i;
        } else {
            B[i - dimensions] = i;
        }
    }

    while (!lp_step(table, B, N, 2 * dimensions + depth, dimensions + depth)) {
        //
    }

    for (long row = 0; row < dimensions + depth; ++row) {
        if (B[row] >= 2 * dimensions) {
            assert(mpq_sgn(matrix_at(table, 2 + row, dimensions)) == 0);

            for (long col = 0; col < dimensions; ++col) {
                if (N[col] < 2 * dimensions && mpq_sgn(matrix_at(table, 2 + row, col)) != 0) {
                    lp_pivot(table, B, N, 2 * dimensions + depth, dimensions + depth, col, row);
                }
            }
        }
    }

    for (long c0 = 0, c1 = dimensions - 1; c0 < dimensions - depth; ++c0) {
        if (N[c0] >= 2 * dimensions) {
            for (;; --c1) {
                if (N[c1] < 2 * dimensions) {
                    for (long row = 0; row < 2 + dimensions + depth; ++row) {
                        mpq_swap(matrix_at(table, row, c0), matrix_at(table, row, c1));
                    }

                    long temp = N[c0];
                    N[c0] = N[c1];
                    N[c1] = temp;

                    break;
                }
            }
        }
    }

    for (long row = 0; row < 2 + dimensions + depth; ++row) {
        mpq_swap(matrix_at(table, row, dimensions - depth), matrix_at(table, row, dimensions));
    }

    matrix_t* view = matrix_view(table, 1, 0, 1 + dimensions + depth, dimensions - depth + 1);

    while (!lp_step(view, B, N, 2 * dimensions, dimensions + depth)) {
        //
    }

    for (long i = 0; i < dimensions; ++i) {
        mpq_set_ui(matrix_at(dest, i, 0), 0, 1);
    }

    for (long i = 0; i < dimensions + depth; ++i) {
        if (B[i] < dimensions) {
            mpq_set(matrix_at(dest, B[i], 0), matrix_at(table, 2 + i, dimensions - depth));
        }
    }

    matrix_free(view);
    matrix_free(table);
}
