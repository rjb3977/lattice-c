#pragma once
#define _POSIX_C_SOURCE 200809L

#include <gmp.h>

void enumerate(long dimensions, mpq_t basis[dimensions][dimensions], mpq_t lower[dimensions], mpq_t upper[dimensions], long* count_out, mpq_t (**results_out)[dimensions]);
