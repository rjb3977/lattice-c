#pragma once
#define _POSIX_C_SOURCE 200809L

#include <stdbool.h>
#include <stdio.h>
#include <gmp.h>

#ifndef NDEBUG
extern size_t max_size;

#undef mpq_add
#undef mpq_sub
#undef mpq_mul
#undef mpq_div

#define __mpq_max_size(a, b, c, f)  do { \
                                        size_t sbn = mpz_sizeinbase(mpq_numref(b), 2); size_t sbd = mpz_sizeinbase(mpq_denref(b), 2); \
                                        size_t scn = mpz_sizeinbase(mpq_numref(c), 2); size_t scd = mpz_sizeinbase(mpq_denref(c), 2); \
                                        f(a, b, c); \
                                        size_t san = mpz_sizeinbase(mpq_numref(a), 2); size_t sad = mpz_sizeinbase(mpq_denref(a), 2); \
                                        if (san > max_size) max_size = san; \
                                        if (sad > max_size) max_size = sad; \
                                        if (sbn > max_size) max_size = sbn; \
                                        if (sbd > max_size) max_size = sbd; \
                                        if (scn > max_size) max_size = scn; \
                                        if (scd > max_size) max_size = scd; \
                                    } while(0)

#define mpq_add(a, b, c) __mpq_max_size(a, b, c, __gmpq_add)
#define mpq_sub(a, b, c) __mpq_max_size(a, b, c, __gmpq_sub)
#define mpq_mul(a, b, c) __mpq_max_size(a, b, c, __gmpq_mul)
#define mpq_div(a, b, c) __mpq_max_size(a, b, c, __gmpq_div)

#endif

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
