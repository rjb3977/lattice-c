#define _POSIX_C_SOURCE 200809L

#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "parse.h"
#include "la.h"

#define error_if(__x)      do { if (__x) { free(line); return false; } } while (false)
#define get_line()         error_if(getline(&line, &n, stream) == -1)
#define get_long(__x, __s) do { __x = strtol(__s, &end, 10); error_if(__s == end); } while (false)

static bool parse_row(FILE *stream, matrix_t* dest, long row) {
    char *line = NULL;
    char *end = NULL;
    size_t n = 0;

    get_line();
    char *pos = line;

    for (long col = 0; col < matrix_cols(dest); ++col) {
        char *pos_end = pos;
        char space;

        while (*pos_end && isspace(*pos_end)) ++pos_end;
        while (*pos_end && !isspace(*pos_end)) ++pos_end;

        space = *pos_end;
        *pos_end = 0;

        error_if(mpq_set_str(matrix_at(dest, row, col), pos, 10));
        mpq_canonicalize(matrix_at(dest, row, col));

        *pos_end = space;
        pos = pos_end;
    }

    char *line_end = pos;
    while(*line_end && isspace(*line_end)) ++line_end;

    error_if(*line_end);
    free(line);

    return true;
}

static bool parse_basis(FILE *stream, matrix_t* dest) {
    char *line = NULL;
    char *end = NULL;
    size_t n = 0;

    get_line();

    for (long row = 0; row < matrix_rows(dest); ++row) {
        error_if(!parse_row(stream, dest, row));
    }

    free(line);

    return true;
}

static bool parse_bounds(FILE *stream, matrix_t* lower, matrix_t* upper) {
    char *line = NULL;
    char *end = NULL;
    size_t n = 0;

    get_line();

    error_if(!parse_row(stream, lower, 0));
    error_if(!parse_row(stream, upper, 0));

    free(line);

    return true;
}

bool parse_data(FILE* stream, matrix_t** basis_out, matrix_t** lower_out, matrix_t** upper_out) {
    char *line = NULL;
    char *end = NULL;
    size_t n = 0;

    long dimensions;

    get_line();
    get_long(dimensions, line);

    error_if(dimensions < 1);

    *basis_out = matrix_alloc(dimensions, dimensions);
    *lower_out = matrix_alloc(1, dimensions);
    *upper_out = matrix_alloc(1, dimensions);

    if (!parse_basis(stream, *basis_out) || !parse_bounds(stream, *lower_out, *upper_out)) {
        matrix_free(*basis_out);
        matrix_free(*lower_out);
        matrix_free(*upper_out);

        free(*basis_out);
        free(*lower_out);
        free(*upper_out);

        *basis_out = NULL;
        *lower_out = NULL;
        *upper_out = NULL;

        error_if(true);
    }

    free(line);

    return true;
}
