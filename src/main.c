#define _POSIX_C_SOURCE 200809L

#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
    mpq_t *basis;
    mpq_t *lower;
    mpq_t *upper;

    if (!parse_data(stream, &dimensions, &basis, &lower, &upper)) {
        fprintf(stderr, "error parsing file %s\n", argc >= 2 ? argv[1] : "(stdin)");
        exit(1);
    }

    fclose(stream);

    long count = 0;
    mpq_t (*results)[dimensions] = NULL;

    long threads = sysconf(_SC_NPROCESSORS_ONLN);

    struct timespec start;
    struct timespec end;

    clock_gettime(CLOCK_MONOTONIC, &start);
    enumerate(dimensions, (void*) basis, (void*) lower, (void*) upper, &count, &results, threads);
    clock_gettime(CLOCK_MONOTONIC, &end);

    long elapsed_h;
    long elapsed_m;
    long elapsed_s;
    long elapsed_ms;

    get_duration(&start, &end, NULL, &elapsed_h, &elapsed_m, &elapsed_s, &elapsed_ms, NULL, NULL);

    mat_print(count, dimensions, results, stdout);
    printf("\n");
    printf("elapsed: %02ld:%02ld:%02ld.%03ld\n", elapsed_h, elapsed_m, elapsed_s, elapsed_ms);
    printf("count:   %ld\n", count);

    mat_clear(dimensions, dimensions, (void*) basis);
    vec_clear(dimensions, lower);
    vec_clear(dimensions, upper);
    mat_clear(count, dimensions, results);

    free(basis);
    free(lower);
    free(upper);
    free(results);

    return 0;
}
