#pragma once
#define _POSIX_C_SOURCE 200809L

#include <stdbool.h>
#include <stdio.h>
#include <gmp.h>

#include "la.h"

void lp_solve(matrix_t* dest, const matrix_t* initial_table, long dimensions, long depth);
