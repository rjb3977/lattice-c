#define _POSIX_C_SOURCE 200809L

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "la.h"

#ifndef NDEBUG
size_t max_size;
#endif

struct matrix_s {
    mpq_t* _data;
    long _stride;
    bool _view;

    long _row;
    long _col;
    long _rows;
    long _cols;
};

long matrix_rows(const matrix_t* src) {
    return src->_rows;
}

long matrix_cols(const matrix_t* src) {
    return src->_cols;
}

mpq_ptr matrix_at(matrix_t* src, long row, long col) {
    assert(0 <= row && row < src->_rows);
    assert(0 <= col && col < src->_cols);

    long _row = src->_row + row;
    long _col = src->_col + col;

    return src->_data[src->_stride * _row + _col];
}

mpq_srcptr matrix_cat(const matrix_t* src, long row, long col) {
    assert(0 <= row && row < src->_rows);
    assert(0 <= col && col < src->_cols);

    long _row = src->_row + row;
    long _col = src->_col + col;

    return src->_data[src->_stride * _row + _col];
}

matrix_t* matrix_alloc(long rows, long cols) {
    matrix_t* dest = malloc(sizeof(matrix_t));
    dest->_data = malloc(rows * cols * sizeof(mpq_t));
    dest->_stride = cols;
    dest->_view = false;
    dest->_row = 0;
    dest->_col = 0;
    dest->_rows = rows;
    dest->_cols = cols;

    for (long row = 0; row < dest->_rows; ++row) {
        for (long col = 0; col < dest->_cols; ++col) {
            mpq_init(matrix_at(dest, row, col));
        }
    }

    return dest;
}

matrix_t* matrix_dup(const matrix_t* src) {
    matrix_t* dest = matrix_alloc(matrix_rows(src), matrix_cols(src));
    matrix_set(dest, src);

    return dest;
}

matrix_t* matrix_view(matrix_t* src, long row, long col, long rows, long cols) {
    assert(0 <= row && row < src->_rows);
    assert(0 <= col && col < src->_rows);
    assert(0 <= rows && row + rows <= src->_rows);
    assert(0 <= cols && col + cols <= src->_cols);

    matrix_t* dest = malloc(sizeof(matrix_t));
    dest->_data = src->_data;
    dest->_stride = src->_stride;
    dest->_view = true;
    dest->_row = src->_row + row;
    dest->_col = src->_col + col;
    dest->_rows = rows;
    dest->_cols = cols;

    return dest;
}

void matrix_free(matrix_t* src) {
    if (!src->_view) {
        for (long row = 0; row < src->_rows; ++row) {
            for (long col = 0; col < src->_cols; ++col) {
                mpq_clear(matrix_at(src, row, col));
            }
        }

        free(src->_data);
    }

    free(src);
}

void matrix_set(matrix_t* dest, const matrix_t* src) {
    assert(dest->_rows == src->_rows);
    assert(dest->_cols == src->_cols);

    for (long row = 0; row < dest->_rows; ++row) {
        for (long col = 0; col < dest->_cols; ++col) {
            mpq_set(matrix_at(dest, row, col), matrix_cat(src, row, col));
        }
    }
}

void matrix_add(matrix_t* dest, const matrix_t* src1, const matrix_t* src2) {
    assert(dest->_rows == src1->_rows);
    assert(dest->_cols == src1->_cols);
    assert(dest->_rows == src2->_rows);
    assert(dest->_cols == src2->_cols);

    for (long row = 0; row < dest->_rows; ++row) {
        for (long col = 0; col < dest->_cols; ++col) {
            mpq_add(matrix_at(dest, row, col), matrix_cat(src1, row, col), matrix_cat(src2, row, col));
        }
    }
}

void matrix_sub(matrix_t* dest, const matrix_t* src1, const matrix_t* src2) {
    assert(dest->_rows == src1->_rows);
    assert(dest->_cols == src1->_cols);
    assert(dest->_rows == src2->_rows);
    assert(dest->_cols == src2->_cols);

    for (long row = 0; row < dest->_rows; ++row) {
        for (long col = 0; col < dest->_cols; ++col) {
            mpq_sub(matrix_at(dest, row, col), matrix_cat(src1, row, col), matrix_cat(src2, row, col));
        }
    }
}

void matrix_neg(matrix_t* dest, const matrix_t* src) {
    assert(dest->_rows == src->_rows);
    assert(dest->_cols == src->_cols);

    for (long row = 0; row < dest->_rows; ++row) {
        for (long col = 0; col < dest->_cols; ++col) {
            mpq_neg(matrix_at(dest, row, col), matrix_cat(src, row, col));
        }
    }
}

void matrix_dot(mpq_ptr dest, const matrix_t* src1, const matrix_t* src2) {
    assert(src1->_rows == src2->_rows);
    assert(src1->_cols == 1);
    assert(src2->_cols == 1);

    mpq_t temp;
    mpq_init(temp);

    mpq_set_ui(dest, 0, 1);

    for (long row = 0; row < src1->_rows; ++row) {
        mpq_mul(temp, matrix_cat(src1, row, 0), matrix_cat(src2, row, 0));
        mpq_add(dest, dest, temp);
    }

    mpq_clear(temp);
}

void matrix_lu(matrix_t* src, long* pivots) {
    assert(src->_rows == src->_cols);

    long size = src->_rows;
    mpq_t temp;
    mpq_init(temp);

    for (long i = 0; i < size; ++i) {
        pivots[i] = i;
    }

    for (long i = 0; i < size; ++i) {
        long pivotRow = -1;

        for (long row = i; row < size; ++row) {
            if (mpq_sgn(matrix_at(src, row, i)) != 0) {
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
                mpq_swap(matrix_at(src, i, col), matrix_at(src, pivotRow, col));
            }
        }

        for (long row = i + 1; row < size; ++row) {
            mpq_div(matrix_at(src, row, i), matrix_at(src, row, i), matrix_at(src, i, i));
        }

        for (long row = i + 1; row < size; ++row) {
            for (long col = i + 1; col < size; ++col) {
                mpq_mul(temp, matrix_at(src, row, i), matrix_at(src, i, col));
                mpq_sub(matrix_at(src, row, col), matrix_at(src, row, col), temp);
            }
        }
    }

    mpq_clear(temp);
}

void solve_p(matrix_t* dest, const long* pivots) {
    long size = dest->_rows;

    for (long col = 0; col < dest->_cols; ++col) {
        for (long row = size - 1; row >= 0; --row) {
            mpq_swap(matrix_at(dest, row, col), matrix_at(dest, pivots[row], col));
        }
    }
}

void solve_pt(matrix_t* dest, const long* pivots) {
    long size = dest->_rows;

    for (long col = 0; col < dest->_cols; ++col) {
        for (long row = 0; row < size; ++row) {
            mpq_swap(matrix_at(dest, row, col), matrix_at(dest, pivots[row], col));
        }
    }
}

void solve_l(matrix_t* dest, const matrix_t* src) {
    assert(src->_rows == src->_cols);
    assert(dest->_rows == src->_rows);

    long size = src->_rows;

    mpq_t temp;
    mpq_init(temp);

    for (long dcol = 0; dcol < dest->_cols; ++dcol) {
        for (long row = 0; row < size; ++row) {
            for (long col = 0; col < row; ++col) {
                mpq_mul(temp, matrix_cat(src, row, col), matrix_at(dest, col, dcol));
                mpq_sub(matrix_at(dest, row, dcol), matrix_at(dest, row, dcol), temp);
            }
        }
    }

    mpq_clear(temp);
}

void solve_lt(matrix_t* dest, const matrix_t* src) {
    assert(src->_rows == src->_cols);
    assert(dest->_rows == src->_rows);

    long size = src->_rows;

    mpq_t temp;
    mpq_init(temp);

    for (long dcol = 0; dcol < dest->_cols; ++dcol) {
        for (long row = size - 1; row >= 0; --row) {
            for (long col = size - 1; col > row; --col) {
                mpq_mul(temp, matrix_cat(src, col, row), matrix_at(dest, col, dcol));
                mpq_sub(matrix_at(dest, row, dcol), matrix_at(dest, row, dcol), temp);
            }
        }
    }

    mpq_clear(temp);
}

void solve_u(matrix_t* dest, const matrix_t* src) {
    assert(src->_rows == src->_cols);
    assert(dest->_rows == src->_rows);

    long size = src->_rows;

    mpq_t temp;
    mpq_init(temp);

    for (long dcol = 0; dcol < dest->_cols; ++dcol) {
        for (long row = size - 1; row >= 0; --row) {
            for (long col = size - 1; col > row; --col) {
                mpq_mul(temp, matrix_cat(src, row, col), matrix_at(dest, col, dcol));
                mpq_sub(matrix_at(dest, row, dcol), matrix_at(dest, row, dcol), temp);
            }

            mpq_div(matrix_at(dest, row, dcol), matrix_at(dest, row, dcol), matrix_cat(src, row, row));
        }
    }

    mpq_clear(temp);
}

void solve_ut(matrix_t* dest, const matrix_t* src) {
    assert(src->_rows == src->_cols);
    assert(dest->_rows == src->_rows);

    long size = src->_rows;

    mpq_t temp;
    mpq_init(temp);

    for (long dcol = 0; dcol < dest->_cols; ++dcol) {
        for (long row = 0; row < size; ++row) {
            for (long col = 0; col < row; ++col) {
                mpq_mul(temp, matrix_cat(src, col, row), matrix_at(dest, col, dcol));
                mpq_sub(matrix_at(dest, row, dcol), matrix_at(dest, row, dcol), temp);
            }

            mpq_div(matrix_at(dest, row, dcol), matrix_at(dest, row, dcol), matrix_cat(src, row, row));
        }
    }

    mpq_clear(temp);
}

void solve_ptlu(matrix_t* dest, const matrix_t* src, const long* pivots) {
    assert(src->_rows == src->_cols);
    assert(dest->_rows == src->_rows);

    solve_pt(dest, pivots);
    solve_l(dest, src);
    solve_u(dest, src);
}

void solve_utltp(matrix_t* dest, const matrix_t* src, const long* pivots) {
    assert(src->_rows == src->_cols);
    assert(dest->_rows == src->_rows);

    solve_ut(dest, src);
    solve_lt(dest, src);
    solve_p(dest, pivots);
}

void matrix_print(FILE* dest, const matrix_t* src) {
    for (long row = 0; row < src->_rows; ++row) {
        if (row != 0) {
            fprintf(dest, "\n");
        }

        for (long col = 0; col < src->_cols; ++col) {
            if (col != 0) {
                fprintf(dest, " ");
            }

            mpq_out_str(dest, 10, matrix_cat(src, row, col));
        }
    }
}

void matrix_print_t(FILE* dest, const matrix_t* src) {
    for (long col = 0; col < src->_cols; ++col) {
        if (col != 0) {
            fprintf(dest, "\n");
        }

        for (long row = 0; row < src->_rows; ++row) {
            if (row != 0) {
                fprintf(dest, " ");
            }

            mpq_out_str(dest, 10, matrix_cat(src, row, col));
        }
    }
}
