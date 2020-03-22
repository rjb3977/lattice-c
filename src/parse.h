#pragma once
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>

#include "la.h"

bool parse_data(FILE* stream, matrix_t** basis_out, matrix_t** lower_out, matrix_t** upper_out);
