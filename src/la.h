#pragma once
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <gmp.h>

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

void mat_init(long rows, long cols, mpq_t dest[rows][cols]);
void mat_clear(long rows, long cols, mpq_t dest[rows][cols]);
void mat_set(long rows, long cols, mpq_t dest[rows][cols], const mpq_t src[rows][cols]);
void mat_swap(long rows, long cols, mpq_t src1[rows][cols], mpq_t src2[rows][cols]);
void mat_copy(long rows, long cols, long dest_row, long dest_col, long dest_stride, mpq_t dest[][dest_stride], long src_row, long src_col, long src_stride, const mpq_t src[][src_stride]);
long mat_rr(long rows, long cols, long limit, mpq_t matrix[rows][cols]);
void mat_lu_2(long size, long stride, mpq_t matrix[size][stride], long pivots[size]);
void mat_print(long rows, long cols, const mpq_t src[rows][cols], FILE* stream);

void solve_p_2(long size, const long pivots[size], mpq_t x[size]);
void solve_pt_2(long size, const long pivots[size], mpq_t x[size]);
void solve_l_2(long size, long stride, const mpq_t lu[size][stride], mpq_t x[size]);
void solve_lt_2(long size, long stride, const mpq_t lu[size][stride], mpq_t x[size]);
void solve_u_2(long size, long stride, const mpq_t lu[size][stride], mpq_t x[size]);
void solve_ut_2(long size, long stride, const mpq_t lu[size][stride], mpq_t x[size]);
void solve_ptlu_2(long size, long stride, const mpq_t lu[size][stride], const long pivots[size], mpq_t x[size]);
void solve_utltp_2(long size, long stride, const mpq_t lu[size][stride], const long pivots[size], mpq_t x[size]);
void solve_ptlu_corrected_2(long size, long stride, const mpq_t lu[size][stride], const long pivots[size], long corrections, const long indices[corrections], const mpq_t columns[corrections][size], mpq_t x[size]);
void solve_utltp_corrected_2(long size, long stride, const mpq_t lu[size][stride], const long pivots[size], long corrections, const long indices[corrections], const mpq_t columns[corrections][size], mpq_t x[size]);
