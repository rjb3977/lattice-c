#pragma once
#define _POSIX_C_SOURCE 200809L

#include "la.h"

void enumerate(const matrix_t* basis, const matrix_t* lower, const matrix_t* upper, long* count_out, matrix_t*** results_out, long thread_max);
