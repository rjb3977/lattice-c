#pragma once
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <gmp.h>

void vec_init(long length, mpq_t dest[length]);
void vec_clear(long length, mpq_t dest[length]);
void vec_copy(long length, long dest_i, mpq_t dest[length], long src_i, mpq_t src[length]);
void vec_add(long length, mpq_t dest[length], mpq_t src1[length], mpq_t src2[length]);
void vec_sub(long length, mpq_t dest[length], mpq_t src1[length], mpq_t src2[length]);
void vec_neg(long length, mpq_t dest[length], mpq_t src[length]);
void vec_dot(long length, mpq_t dest, mpq_t src1[length], mpq_t src2[length]);
void vec_print(long length, mpq_t src[length], FILE* stream);

void mat_init(long rows, long cols, mpq_t dest[rows][cols]);
void mat_clear(long rows, long cols, mpq_t dest[rows][cols]);
void mat_copy(long rows, long cols, long dest_row, long dest_col, long dest_stride, mpq_t dest[][dest_stride], long src_row, long src_col, long src_stride, mpq_t src[][src_stride]);
long mat_rr(long rows, long cols, long limit, mpq_t matrix[rows][cols]);
void mat_lu(long size, mpq_t matrix[size][size], long pivots[size]);
void mat_print(long rows, long cols, mpq_t src[rows][cols], FILE* stream);

void solve_p(long size, long pivots[size], mpq_t x[size]);
void solve_pt(long size, long pivots[size], mpq_t x[size]);
void solve_l(long size, mpq_t lu[size][size], mpq_t x[size]);
void solve_lt(long size, mpq_t lu[size][size], mpq_t x[size]);
void solve_u(long size, mpq_t lu[size][size], mpq_t x[size]);
void solve_ut(long size, mpq_t lu[size][size], mpq_t x[size]);
void solve_ptlu(long size, mpq_t lu[size][size], long pivots[size], mpq_t x[size]);
void solve_utltp(long size, mpq_t lu[size][size], long pivots[size], mpq_t x[size]);
void solve_ptlu_corrected(long size, mpq_t lu[size][size], long pivots[size], long corrections, long indices[corrections], mpq_t columns[corrections][size], mpq_t x[size]);
void solve_utltp_corrected(long size, mpq_t lu[size][size], long pivots[size], long corrections, long indices[corrections], mpq_t columns[corrections][size], mpq_t x[size]);
