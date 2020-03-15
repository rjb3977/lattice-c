#pragma once
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <gmp.h>

bool parse_data(FILE *stream, long *dimensions_out, mpq_t **basis_out, mpq_t **lower_out, mpq_t **upper_out);
