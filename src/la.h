#pragma once
#define _POSIX_C_SOURCE 200809L

#include <stdbool.h>
#include <stdio.h>
#include <gmp.h>

typedef struct matrix_s matrix_t;

long matrix_rows(const matrix_t* src);
long matrix_cols(const matrix_t* src);

mpq_ptr matrix_at(matrix_t* src, long row, long col);
mpq_srcptr matrix_cat(const matrix_t* src, long row, long col);

matrix_t* matrix_alloc(long rows, long cols);
matrix_t* matrix_dup(const matrix_t* src);
matrix_t* matrix_view(matrix_t* src, long row, long col, long rows, long cols);
void matrix_free(matrix_t* src);

void matrix_set(matrix_t* dest, const matrix_t* src);
void matrix_add(matrix_t* dest, const matrix_t* src1, const matrix_t* src2);
void matrix_sub(matrix_t* dest, const matrix_t* src1, const matrix_t* src2);
void matrix_neg(matrix_t* dest, const matrix_t* src);
void matrix_dot(mpq_ptr dest, const matrix_t* src1, const matrix_t* src2);

void matrix_lu(matrix_t* src, long* pivots);
void solve_p(matrix_t* dest, const long* pivots);
void solve_pt(matrix_t* dest, const long* pivots);
void solve_l(matrix_t* dest, const matrix_t* src);
void solve_lt(matrix_t* dest, const matrix_t* src);
void solve_u(matrix_t* dest, const matrix_t* src);
void solve_ut(matrix_t* dest, const matrix_t* src);
void solve_ptlu(matrix_t* dest, const matrix_t* src, const long* pivots);
void solve_utltp(matrix_t* dest, const matrix_t* src, const long* pivots);

void matrix_print(FILE* dest, const matrix_t* src);
void matrix_print_t(FILE* dest, const matrix_t* src);
