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

static bool parse_row(FILE *stream, long dimensions, mpq_t row[dimensions]) {
    char *line = NULL;
    char *end = NULL;
    size_t n = 0;

    get_line();
    char *pos = line;

    for (long col = 0; col < dimensions; ++col) {
        char *pos_end = pos;
        char space;

        while (*pos_end && isspace(*pos_end)) ++pos_end;
        while (*pos_end && !isspace(*pos_end)) ++pos_end;

        space = *pos_end;
        *pos_end = 0;

        error_if(mpq_set_str(row[col], pos, 10));
        mpq_canonicalize(row[col]);

        *pos_end = space;
        pos = pos_end;
    }

    char *line_end = pos;
    while(*line_end && isspace(*line_end)) ++line_end;

    error_if(*line_end);
    free(line);

    return true;
}

static bool parse_basis(FILE *stream, long dimensions, mpq_t basis[dimensions][dimensions]) {
    char *line = NULL;
    char *end = NULL;
    size_t n = 0;

    get_line();

    for (long row = 0; row < dimensions; ++row) {
        error_if(!parse_row(stream, dimensions, basis[row]));
    }

    free(line);

    return true;
}

static bool parse_bounds(FILE *stream, long dimensions, mpq_t lower[dimensions], mpq_t upper[dimensions]) {
    char *line = NULL;
    char *end = NULL;
    size_t n = 0;

    get_line();

    error_if(!parse_row(stream, dimensions, lower));
    error_if(!parse_row(stream, dimensions, upper));

    free(line);

    return true;
}

bool parse_data(FILE *stream, long *dimensions_out, mpq_t **basis_out, mpq_t **lower_out, mpq_t **upper_out) {
    char *line = NULL;
    char *end = NULL;
    size_t n = 0;

    get_line();
    get_long(*dimensions_out, line);

    error_if(*dimensions_out < 1);

    *basis_out = malloc(*dimensions_out * *dimensions_out * sizeof(mpq_t));
    *lower_out = malloc(*dimensions_out * sizeof(mpq_t));
    *upper_out = malloc(*dimensions_out * sizeof(mpq_t));

    mat_init(*dimensions_out, *dimensions_out, (void*) *basis_out);
    vec_init(*dimensions_out, *lower_out);
    vec_init(*dimensions_out, *upper_out);

    if (!parse_basis(stream, *dimensions_out, (void*) *basis_out) || !parse_bounds(stream, *dimensions_out, *lower_out, *upper_out)) {
        mat_clear(*dimensions_out, *dimensions_out, (void*) *basis_out);
        vec_clear(*dimensions_out, *lower_out);
        vec_clear(*dimensions_out, *upper_out);

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
