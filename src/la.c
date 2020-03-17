#define _POSIX_C_SOURCE 200809L

#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <gmp.h>

#include "la.h"

mpq_t *vec_malloc(long length) {
    mpq_t *dest = malloc(length * sizeof(mpq_t));
    vec_init(length, dest);

    return dest;
}

mpq_t *vec_dup(long length, const mpq_t *src) {
    mpq_t *dest = vec_malloc(length);
    vec_set(length, dest, src);

    return dest;
}

void vec_free(long length, mpq_t *src) {
    vec_clear(length, src);
    free(src);
}

void vec_init(long length, mpq_t dest[length]) {
    for (long i = 0; i < length; ++i) {
        mpq_init(dest[i]);
    }
}

void vec_clear(long length, mpq_t dest[length]) {
    for (long i = 0; i < length; ++i) {
        mpq_clear(dest[i]);
    }
}

void vec_set(long length, mpq_t dest[length], const mpq_t src[length]) {
    for (long i = 0; i < length; ++i) {
        mpq_set(dest[i], src[i]);
    }
}

void vec_swap(long length, mpq_t src1[length], mpq_t src2[length]) {
    for (long i = 0; i < length; ++i) {
        mpq_swap(src1[i], src2[i]);
    }
}

void vec_copy(long length, long dest_i, mpq_t dest[length], long src_i, const mpq_t src[length]) {
    for (long i = 0; i < length; ++i) {
        mpq_set(dest[dest_i + i], src[src_i + i]);
    }
}

void vec_add(long length, mpq_t dest[length], const mpq_t src1[length], const mpq_t src2[length]) {
    for (long i = 0; i < length; ++i) {
        mpq_add(dest[i], src1[i], src2[i]);
    }
}

void vec_sub(long length, mpq_t dest[length], const mpq_t src1[length], const mpq_t src2[length]) {
    for (long i = 0; i < length; ++i) {
        mpq_sub(dest[i], src1[i], src2[i]);
    }
}

void vec_neg(long length, mpq_t dest[length], const mpq_t src[length]) {
    for (long i = 0; i < length; ++i) {
        mpq_neg(dest[i], src[i]);
    }
}

void vec_dot(long length, mpq_t dest, const mpq_t src1[length], const mpq_t src2[length]) {
    mpq_t temp;
    mpq_init(temp);
    mpq_set_ui(dest, 0, 1);

    for (long i = 0; i < length; ++i) {
        mpq_mul(temp, src1[i], src2[i]);
        mpq_add(dest, dest, temp);
    }

    mpq_clear(temp);
}

void vec_print(long length, const mpq_t src[length], FILE* stream) {
    for (long i = 0; i < length; ++i) {
        if (i != 0) {
            printf(", ");
        }

        mpq_out_str(stream, 10, src[i]);
    }
}

mpq_t *mat_malloc(long rows, long cols) {
    mpq_t *dest = malloc(rows * cols * sizeof(mpq_t));
    mat_init(rows, cols, (void*) dest);

    return dest;
}

mpq_t *mat_dup(long rows, long cols, const mpq_t *src) {
    mpq_t *dest = mat_malloc(rows, cols);
    mat_set(rows, cols, (void*) dest, (const void*) src);

    return dest;
}

void mat_free(long rows, long cols, mpq_t *src) {
    mat_clear(rows, cols, (void*) src);
    free(src);
}

void mat_init(long rows, long cols, mpq_t dest[rows][cols]) {
    for (long row = 0; row < rows; ++row) {
        for (long col = 0; col < cols; ++col) {
            mpq_init(dest[row][col]);
        }
    }
}

void mat_clear(long rows, long cols, mpq_t dest[rows][cols]) {
    for (long row = 0; row < rows; ++row) {
        for (long col = 0; col < cols; ++col) {
            mpq_clear(dest[row][col]);
        }
    }
}

void mat_set(long rows, long cols, mpq_t dest[rows][cols], const mpq_t src[rows][cols]) {
    for (long row = 0; row < rows; ++row) {
        for (long col = 0; col < cols; ++col) {
            mpq_set(dest[row][col], src[row][col]);
        }
    }
}

void mat_swap(long rows, long cols, mpq_t src1[rows][cols], mpq_t src2[rows][cols]) {
    for (long row = 0; row < rows; ++row) {
        for (long col = 0; col < cols; ++col) {
            mpq_swap(src1[row][col], src2[row][col]);
        }
    }
}

void mat_copy(long rows, long cols, long dest_row, long dest_col, long dest_stride, mpq_t dest[][dest_stride], long src_row, long src_col, long src_stride, const mpq_t src[][src_stride]) {
    for (long row = 0; row < rows; ++row) {
        for (long col = 0; col < cols; ++col) {
            mpq_set(dest[dest_row + row][dest_col + col], src[src_row + row][src_col + col]);
        }
    }
}

long mat_rr(long rows, long cols, long limit, mpq_t matrix[rows][cols]) {
    mpq_t temp;
    mpq_init(temp);

    long map[limit];
    long rank = 0;

    for (long i = 0; i < limit && rank < rows; ++i) {
        long pivotRow = -1;

        for (long row = rank; row < rows; ++row) {
            if (mpq_sgn(matrix[row][i]) != 0) {
                pivotRow = row;
                break;
            }
        }

        if (pivotRow == -1) {
            map[i] = -1;
            continue;
        }

        if (pivotRow != rank) {
            for (long col = rank; col < cols; ++col) {
                mpq_swap(matrix[rank][col], matrix[pivotRow][col]);
            }
        }

        for (long col = rank; col < rows; ++col) {
            mpq_div(matrix[rank][col], matrix[rank][col], matrix[rank][i]);
        }

        for (long row = rank + 1; row < rows; ++row) {
            mpq_set_ui(matrix[row][rank], 0, 1);

            for (long col = rank + 1; col < cols; ++col) {
                mpq_mul(temp, matrix[rank][col], matrix[row][i]);
                mpq_sub(matrix[row][col], matrix[row][col], matrix[rank][col]);
            }
        }
    }

    mpq_clear(temp);
    return rank;
}

void mat_print(long rows, long cols, const mpq_t src[rows][cols], FILE* stream) {
    for (long row = 0; row < rows; ++row) {
        if (row != 0) {
            printf("\n");
        }

        for (long col = 0; col < cols; ++col) {
            if (col != 0) {
                printf(", ");
            }

            mpq_out_str(stream, 10, src[row][col]);
        }
    }
}

void mat_lu(long size, long stride, mpq_t matrix[size][stride], long pivots[size]) {
    mpq_t temp;
    mpq_init(temp);

    for (long i = 0; i < size; ++i) {
        pivots[i] = i;
    }

    for (long i = 0; i < size; ++i) {
        long pivotRow = -1;

        for (long row = i; row < size; ++row) {
            if (mpq_sgn(matrix[row][i]) != 0) {
                pivotRow = row;
                break;
            }
        }

        if (pivotRow == -1) {
            continue;
        }

        pivots[i] = pivotRow;

        if (pivotRow != i) {
            for (long col = 0; col < size; ++col) {
                mpq_swap(matrix[i][col], matrix[pivotRow][col]);
            }
        }

        for (long row = i + 1; row < size; ++row) {

            mpq_div(matrix[row][i], matrix[row][i], matrix[i][i]);
        }

        for (long row = i + 1; row < size; ++row) {
            for (long col = i + 1; col < size; ++col) {
                mpq_mul(temp, matrix[row][i], matrix[i][col]);
                mpq_sub(matrix[row][col], matrix[row][col], temp);
            }
        }
    }

    mpq_clear(temp);
}

void solve_p(long size, const long pivots[size], mpq_t x[size]) {
    for (long i = size - 1; i >= 0; --i) {
        mpq_swap(x[i], x[pivots[i]]);
    }
}

void solve_pt(long size, const long pivots[size], mpq_t x[size]) {
    for (long i = 0; i < size; ++i) {
        mpq_swap(x[i], x[pivots[i]]);
    }
}

void solve_l(long size, long stride, const mpq_t lu[size][stride], mpq_t x[size]) {
    mpq_t temp;
    mpq_init(temp);

    for (long row = 0; row < size; ++row) {
        for (long col = 0; col < row; ++col) {
            mpq_mul(temp, lu[row][col], x[col]);
            mpq_sub(x[row], x[row], temp);
        }
    }

    mpq_clear(temp);
}

void solve_lt(long size, long stride, const mpq_t lu[size][stride], mpq_t x[size]) {
    mpq_t temp;
    mpq_init(temp);

    for (long row = size - 1; row >= 0; --row) {
        for (long col = size - 1; col > row; --col) {
            mpq_mul(temp, lu[col][row], x[col]);
            mpq_sub(x[row], x[row], temp);
        }
    }

    mpq_clear(temp);
}

void solve_u(long size, long stride, const mpq_t lu[size][stride], mpq_t x[size]) {
    mpq_t temp;
    mpq_init(temp);

    for (long row = size - 1; row >= 0; --row) {
        for (long col = size - 1; col > row; --col) {
            mpq_mul(temp, lu[row][col], x[col]);
            mpq_sub(x[row], x[row], temp);
        }

        mpq_div(x[row], x[row], lu[row][row]);
    }

    mpq_clear(temp);
}

void solve_ut(long size, long stride, const mpq_t lu[size][stride], mpq_t x[size]) {
    mpq_t temp;
    mpq_init(temp);

    for (long row = 0; row < size; ++row) {
        for (long col = 0; col < row; ++col) {
            mpq_mul(temp, lu[col][row], x[col]);
            mpq_sub(x[row], x[row], temp);
        }

        mpq_div(x[row], x[row], lu[row][row]);
    }

    mpq_clear(temp);
}

void solve_ptlu(long size, long stride, const mpq_t lu[size][stride], const long pivots[size], mpq_t x[size]) {
    solve_pt(size, pivots, x);
    solve_l(size, stride, lu, x);
    solve_u(size, stride, lu, x);
}

void solve_utltp(long size, long stride, const mpq_t lu[size][stride], const long pivots[size], mpq_t x[size]) {
    solve_ut(size, stride, lu, x);
    solve_lt(size, stride, lu, x);
    solve_p(size, pivots, x);
}

void solve_ptlu_corrected(long size, long stride, const mpq_t lu[size][stride], const long pivots[size], long corrections, const long indices[corrections], const mpq_t columns[corrections][size], mpq_t x[size]) {
    mpq_t temp;
    mpq_init(temp);

    solve_ptlu(size, stride, lu, pivots, x);

    for (long i = 0; i < corrections; ++i) {
        long index = indices[i];
        const mpq_t* column = columns[i];
        mpq_ptr xp = x[index];

        mpq_div(xp, xp, column[index]);

        for (long row = 0; row < size; ++row) {
            if (row != index) {
                mpq_mul(temp, xp, column[row]);
                mpq_sub(x[row], x[row], temp);
            }
        }
    }

    mpq_clear(temp);
}

void solve_utltp_corrected(long size, long stride, const mpq_t lu[size][stride], const long pivots[size], long corrections, const long indices[corrections], const mpq_t columns[corrections][size], mpq_t x[size]) {
    mpq_t temp[2];
    mpq_init(temp[0]);
    mpq_init(temp[1]);

    for (long i = corrections - 1; i >= 0; --i) {
        long index = indices[i];
        const mpq_t* column = columns[i];

        mpq_set_ui(temp[0], 0, 1);

        for (long row = 0; row < size; ++row) {
            if (row != index) {
                mpq_mul(temp[1], x[row], column[row]);
                mpq_add(temp[0], temp[0], temp[1]);
            }
        }

        mpq_sub(x[index], x[index], temp[0]);
        mpq_div(x[index], x[index], column[index]);
    }

    solve_utltp(size, stride, lu, pivots, x);

    mpq_clear(temp[0]);
    mpq_clear(temp[1]);
}
