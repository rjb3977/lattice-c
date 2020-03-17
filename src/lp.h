#pragma once
#define _POSIX_C_SOURCE 200809L

#include <gmp.h>

void simplex_solve(long size, long depth, mpq_t A[2 * size][3 * size], mpq_t b[2 * size], mpq_t c[2][3 * size], mpq_t λ[2 * size], mpq_t s[size], mpq_t d[2 * size], mpq_t x[2 * size], long B[2 * size], long N[size], mpq_t lu[2 * size][2 * size], long pivots[2 * size], long *corrections, long **indices, mpq_t (**columns)[2 * size]);
