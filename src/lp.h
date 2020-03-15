#pragma once

#include <stdbool.h>
#include <gmp.h>

void simplex_solve(long rows, long cols, mpq_t A[rows][cols], mpq_t b[rows], mpq_t c[cols], mpq_t x[cols]);
void simplex_init(long rows, long cols, mpq_t A[rows][cols], mpq_t b[rows], mpq_t c[cols], mpq_t λ[rows], mpq_t s[cols - rows], mpq_t d[rows], mpq_t x[rows], long B[rows], long N[cols - rows], mpq_t lu[rows][rows], long pivots[rows], long *corrections, long **indices, mpq_t (**columns)[rows]);
bool simplex_step(long rows, long cols, mpq_t A[rows][cols], mpq_t b[rows], mpq_t c[cols], mpq_t λ[rows], mpq_t s[cols - rows], mpq_t d[rows], mpq_t x[rows], long B[rows], long N[cols - rows], mpq_t lu[rows][rows], long pivots[rows], long *corrections, long **indices, mpq_t (**columns)[rows]);
