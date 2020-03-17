#pragma once
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <gmp.h>

mpq_t *vec_malloc(long length);
mpq_t *vec_dup(long length, const mpq_t *src);
void vec_free(long length, mpq_t *src);
void vec_init(long length, mpq_t dest[length]);
void vec_clear(long length, mpq_t dest[length]);
void vec_set(long length, mpq_t dest[length], const mpq_t src[length]);
void vec_swap(long length, mpq_t src1[length], mpq_t src2[length]);
void vec_copy(long length, long dest_i, mpq_t dest[length], long src_i, const mpq_t src[length]);
void vec_add(long length, mpq_t dest[length], const mpq_t src1[length], const mpq_t src2[length]);
void vec_sub(long length, mpq_t dest[length], const mpq_t src1[length], const mpq_t src2[length]);
void vec_neg(long length, mpq_t dest[length], const mpq_t src[length]);
void vec_dot(long length, mpq_t dest, const mpq_t src1[length], const mpq_t src2[length]);
void vec_print(long length, const mpq_t src[length], FILE* stream);

mpq_t *mat_malloc(long rows, long cols);
mpq_t *mat_dup(long rows, long cols, const mpq_t *src);
void mat_free(long rows, long cols, mpq_t *src);
void mat_init(long rows, long cols, mpq_t dest[rows][cols]);
void mat_clear(long rows, long cols, mpq_t dest[rows][cols]);
void mat_set(long rows, long cols, mpq_t dest[rows][cols], const mpq_t src[rows][cols]);
void mat_swap(long rows, long cols, mpq_t src1[rows][cols], mpq_t src2[rows][cols]);
void mat_copy(long rows, long cols, long dest_row, long dest_col, long dest_stride, mpq_t dest[][dest_stride], long src_row, long src_col, long src_stride, const mpq_t src[][src_stride]);
long mat_rr(long rows, long cols, long limit, mpq_t matrix[rows][cols]);
void mat_lu(long size, long stride, mpq_t matrix[size][stride], long pivots[size]);
void mat_print(long rows, long cols, const mpq_t src[rows][cols], FILE* stream);

void solve_p(long size, const long pivots[size], mpq_t x[size]);
void solve_pt(long size, const long pivots[size], mpq_t x[size]);
void solve_l(long size, long stride, const mpq_t lu[size][stride], mpq_t x[size]);
void solve_lt(long size, long stride, const mpq_t lu[size][stride], mpq_t x[size]);
void solve_u(long size, long stride, const mpq_t lu[size][stride], mpq_t x[size]);
void solve_ut(long size, long stride, const mpq_t lu[size][stride], mpq_t x[size]);
void solve_ptlu(long size, long stride, const mpq_t lu[size][stride], const long pivots[size], mpq_t x[size]);
void solve_utltp(long size, long stride, const mpq_t lu[size][stride], const long pivots[size], mpq_t x[size]);
void solve_ptlu_corrected(long size, long stride, const mpq_t lu[size][stride], const long pivots[size], long corrections, const long indices[corrections], const mpq_t columns[corrections][size], mpq_t x[size]);
void solve_utltp_corrected(long size, long stride, const mpq_t lu[size][stride], const long pivots[size], long corrections, const long indices[corrections], const mpq_t columns[corrections][size], mpq_t x[size]);
