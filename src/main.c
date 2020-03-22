#define _POSIX_C_SOURCE 200809L

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <unistd.h>

#include "parse.h"
#include "enumerate.h"
#include "la.h"
#include "lp.h"

static void get_duration(const struct timespec* start, const struct timespec* end, long* d_out, long* h_out, long* m_out, long* s_out, long* ms_out, long* us_out, long* ns_out) {
    long s = end->tv_sec - start->tv_sec;
    long ns = end->tv_nsec - start->tv_nsec;

    if (ns < 0)
    {
        ns += 1000000000;
        s -= 1;
    }

    if (ns_out != NULL) *ns_out = ns % 1000;
    if (us_out != NULL) *us_out = (ns / 1000) % 1000;
    if (ms_out != NULL) *ms_out = (ns / 1000 / 1000) % 1000;
    if (s_out != NULL) *s_out = s % 60;
    if (m_out != NULL) *m_out = (s / 60) % 60;
    if (h_out != NULL) *h_out = (s / 60 / 60) % 24;
    if (d_out != NULL) *d_out = s / 60 / 60 / 24;
}

static void test(void) {
    matrix_t* table = matrix_alloc(6, 5);
    matrix_t* x = matrix_alloc(4, 1);

    mpq_set_si(matrix_at(table, 0, 0), 1, 2);
    mpq_set_si(matrix_at(table, 0, 1), 0, 1);
    mpq_set_si(matrix_at(table, 0, 2), -1, 2);
    mpq_set_si(matrix_at(table, 0, 3), 1, 2);
    mpq_set_si(matrix_at(table, 0, 4), 0, 1);

    mpq_set_si(matrix_at(table, 1, 0), 1, 1);
    mpq_set_si(matrix_at(table, 1, 1), 0, 1);
    mpq_set_si(matrix_at(table, 1, 2), 0, 1);
    mpq_set_si(matrix_at(table, 1, 3), 0, 1);
    mpq_set_si(matrix_at(table, 1, 4), 4, 1);

    mpq_set_si(matrix_at(table, 2, 0), 0, 1);
    mpq_set_si(matrix_at(table, 2, 1), 1, 1);
    mpq_set_si(matrix_at(table, 2, 2), 0, 1);
    mpq_set_si(matrix_at(table, 2, 3), 0, 1);
    mpq_set_si(matrix_at(table, 2, 4), 4, 1);

    mpq_set_si(matrix_at(table, 3, 0), 0, 1);
    mpq_set_si(matrix_at(table, 3, 1), 0, 1);
    mpq_set_si(matrix_at(table, 3, 2), 1, 1);
    mpq_set_si(matrix_at(table, 3, 3), 0, 1);
    mpq_set_si(matrix_at(table, 3, 4), 4, 1);

    mpq_set_si(matrix_at(table, 4, 0), 0, 1);
    mpq_set_si(matrix_at(table, 4, 1), 0, 1);
    mpq_set_si(matrix_at(table, 4, 2), 0, 1);
    mpq_set_si(matrix_at(table, 4, 3), 1, 1);
    mpq_set_si(matrix_at(table, 4, 4), 4, 1);

    mpq_set_si(matrix_at(table, 5, 0), -1, 2);
    mpq_set_si(matrix_at(table, 5, 1), 0, 1);
    mpq_set_si(matrix_at(table, 5, 2), 1, 2);
    mpq_set_si(matrix_at(table, 5, 3), 0, 1);
    mpq_set_si(matrix_at(table, 5, 4), 1, 1);

    lp_solve(x, table, 4, 1);

    exit(0);
}

int main(int argc, char** argv) {
    FILE* stream = stdin;

    if (argc >= 2) {
        stream = fopen(argv[1], "r");

        if (!stream) {
            fprintf(stderr, "error opening file %s\n", argv[1]);
            exit(1);
        }
    }

    long dimensions;
    matrix_t* basis;
    matrix_t* lower;
    matrix_t* upper;

    if (!parse_data(stream, &basis, &lower, &upper)) {
        fprintf(stderr, "error parsing file %s\n", argc >= 2 ? argv[1] : "(stdin)");
        exit(1);
    }

    fclose(stream);

    long count = 0;
    matrix_t** results = NULL;

    long threads = sysconf(_SC_NPROCESSORS_ONLN);
    // long threads = 1;

    struct timespec start;
    struct timespec end;

    clock_gettime(CLOCK_MONOTONIC, &start);
    enumerate(basis, lower, upper, &count, &results, threads);
    clock_gettime(CLOCK_MONOTONIC, &end);

    long elapsed_h;
    long elapsed_m;
    long elapsed_s;
    long elapsed_ms;

    get_duration(&start, &end, NULL, &elapsed_h, &elapsed_m, &elapsed_s, &elapsed_ms, NULL, NULL);

    for (long i = 0; i < count; ++i) {
        matrix_print_t(stdout, results[i]);
        printf("\n");
    }

    printf("\n");
    printf("elapsed: %02ld:%02ld:%02ld.%03ld\n", elapsed_h, elapsed_m, elapsed_s, elapsed_ms);
    printf("count:   %ld\n", count);

    matrix_free(basis);
    matrix_free(lower);
    matrix_free(upper);

    for (long i = 0; i < count; ++i) {
        matrix_free(results[i]);
    }

    free(results);

    return 0;
}
