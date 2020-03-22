#define _POSIX_C_SOURCE 200809L

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <pthread.h>

#include "la.h"
#include "lp.h"

typedef struct {
    long dimensions;
    long depth;

    const matrix_t* transform; // mpq_t[dimensions][2 * dimensions]
    const matrix_t* offset;    // mpq_t[dimensions]
    matrix_t* fixed;           // mpq_t[dimensions]

    matrix_t* table;           // mpq_t[2 * dimensions][3 * dimensions]
    matrix_t* x;               // mpq_t[2 * dimensions]

    long* count_out;
    matrix_t*** results_out;  // mpq_t[*count_out][dimensions]

    pthread_mutex_t *mutex;
    pthread_cond_t *finished;
    volatile long* thread_count;
    long thread_max;
} search_info_t;

static search_info_t *search_info_dup(const search_info_t* src) {
    search_info_t* dest = malloc(sizeof(search_info_t));

    dest->dimensions = src->dimensions;
    dest->depth = src->depth;

    dest->transform = src->transform;
    dest->offset = src->offset;
    dest->fixed = matrix_dup(src->fixed);

    dest->table = matrix_dup(src->table);
    dest->x = matrix_dup(src->x);

    dest->count_out = src->count_out;
    dest->results_out = src->results_out;

    dest->mutex = src->mutex;
    dest->finished = src->finished;
    dest->thread_count = src->thread_count;
    dest->thread_max = src->thread_max;

    return dest;
}

static void search_info_free(search_info_t* src) {
    matrix_free(src->fixed);
    matrix_free(src->table);
    matrix_free(src->x);

    free(src);
}

static void* search_thread(void*);

static void search(search_info_t *info) {
    if (info->depth == info->dimensions) {
        pthread_mutex_lock(info->mutex);

        printf("found: ");
        matrix_print_t(stdout, info->fixed);
        printf("\n");

        *info->count_out += 1;
        *info->results_out = realloc(*info->results_out, *info->count_out * sizeof(matrix_t*));
        (*info->results_out)[*info->count_out - 1] = matrix_dup(info->fixed);

        pthread_mutex_unlock(info->mutex);
    } else {
        mpq_t t0;
        mpq_t t1;
        mpq_t t2;
        mpz_t min;
        mpz_t max;

        mpq_init(t0);
        mpq_init(t1);
        mpq_init(t2);
        mpz_init(min);
        mpz_init(max);

        // lower bound
        for (long col = 0; col < info->dimensions; ++col) {
            mpq_neg(matrix_at(info->table, 0, col), matrix_cat(info->transform, info->depth, col));
        }

        lp_solve(info->x, info->table, info->dimensions, info->depth);

        mpq_set(t0, matrix_cat(info->offset, info->depth, 0));

        for (long col = 0; col < info->dimensions; ++col) {
            mpq_mul(t1, matrix_cat(info->transform, info->depth, col), matrix_cat(info->x, col, 0));
            mpq_add(t0, t0, t1);
        }

        mpz_cdiv_q(min, mpq_numref(t0), mpq_denref(t0));
        mpq_set(t2, t0);

        // upper bound
        for (long col = 0; col < info->dimensions; ++col) {
            mpq_set(matrix_at(info->table, 0, col), matrix_cat(info->transform, info->depth, col));
        }

        lp_solve(info->x, info->table, info->dimensions, info->depth);

        mpq_set(t0, matrix_cat(info->offset, info->depth, 0));

        for (long col = 0; col < info->dimensions; ++col) {
            mpq_mul(t1, matrix_cat(info->transform, info->depth, col), matrix_cat(info->x, col, 0));
            mpq_add(t0, t0, t1);
        }

        mpz_fdiv_q(max, mpq_numref(t0), mpq_denref(t0));

        // rhs of new row = offset - x
        mpq_set_z(t0, min);
        mpq_sub(t0, matrix_cat(info->offset, info->depth, 0), t0);

        // min <= max, min += 1, rhs -= 1
        for (; mpz_cmp(min, max) <= 0; mpz_add_ui(min, min, 1), mpz_sub(mpq_numref(t0), mpq_numref(t0), mpq_denref(t0))) {
            if (mpq_sgn(t0) >= 0) {
                for (long col = 0; col < info->dimensions; ++col) {
                    mpq_neg(matrix_at(info->table, 1 + info->dimensions + info->depth, col), matrix_cat(info->transform, info->depth, col));
                }

                mpq_set(matrix_at(info->table, 1 + info->dimensions + info->depth, info->dimensions), t0);
            } else {
                for (long col = 0; col < info->dimensions; ++col) {
                    mpq_set(matrix_at(info->table, 1 + info->dimensions + info->depth, col), matrix_cat(info->transform, info->depth, col));
                }

                mpq_neg(matrix_at(info->table, 1 + info->dimensions + info->depth, info->dimensions), t0);
            }

            mpq_set_z(matrix_at(info->fixed, info->depth, 0), min);
            info->depth += 1;

            pthread_mutex_lock(info->mutex);

            if (*info->thread_count < info->thread_max) {
                *info->thread_count += 1;
                pthread_mutex_unlock(info->mutex);

                pthread_t thread;
                pthread_create(&thread, NULL, search_thread, search_info_dup(info));
                pthread_detach(thread);
            } else {
                pthread_mutex_unlock(info->mutex);
                search(info);
            }

            info->depth -= 1;
        }

        mpq_clear(t0);
        mpq_clear(t1);
        mpq_clear(t2);
        mpz_clear(min);
        mpz_clear(max);
    }
}

void* search_thread(void* data) {
    search_info_t* info = data;
    search(info);

    pthread_mutex_lock(info->mutex);
    *info->thread_count -= 1;
    pthread_cond_broadcast(info->finished);
    pthread_mutex_unlock(info->mutex);

    search_info_free(info);

    return NULL;
}

static void search_parallel(search_info_t *root) {
    *root->thread_count = 1;

    pthread_t thread;
    pthread_create(&thread, NULL, search_thread, search_info_dup(root));
    pthread_detach(thread);
    pthread_mutex_lock(root->mutex);

    while (*root->thread_count > 0) {
        pthread_cond_wait(root->finished, root->mutex);
    }

    pthread_mutex_unlock(root->mutex);
}

void enumerate(const matrix_t* basis, const matrix_t* lower, const matrix_t* upper, long* count_out, matrix_t*** results_out, long thread_max) {
    assert(matrix_rows(basis) == matrix_cols(basis));

    long dimensions = matrix_rows(basis);

    search_info_t* root = malloc(sizeof(search_info_t));

    matrix_t* transform = matrix_alloc(dimensions, dimensions);
    matrix_t* offset = matrix_alloc(dimensions, 1);

    for (long i = 0; i < dimensions; ++i) {
        mpq_set_ui(matrix_at(transform, i, i), 1, 1);
        mpq_set(matrix_at(offset, i, 0), matrix_cat(lower, 0, i));
    }

    matrix_t* lu = matrix_dup(basis);
    long pivots[dimensions];

    matrix_lu(lu, pivots);
    solve_utltp(transform, lu, pivots);
    solve_utltp(offset, lu, pivots);

    matrix_free(lu);

    root->dimensions = dimensions;
    root->depth = 0;

    root->transform = transform;
    root->offset = offset;
    root->fixed = matrix_alloc(dimensions, 1);

    root->table = matrix_alloc(2 * dimensions + 1, dimensions + 1);
    root->x = matrix_alloc(dimensions, 1);

    root->count_out = count_out;
    root->results_out = results_out;

    root->mutex = malloc(sizeof(pthread_mutex_t));
    root->finished = malloc(sizeof(pthread_cond_t));
    root->thread_count = malloc(sizeof(long));
    root->thread_max = thread_max;

    pthread_mutex_init(root->mutex, NULL);
    pthread_cond_init(root->finished, NULL);

    for (long i = 0; i < dimensions; ++i) {
        mpq_set_ui(matrix_at(root->table, i + 1, i), 1, 1);
        mpq_sub(matrix_at(root->table, i + 1, dimensions), matrix_cat(upper, 0, i), matrix_cat(lower, 0, i));
    }

    search_parallel(root);

    matrix_free(transform);
    matrix_free(offset);

    pthread_mutex_destroy(root->mutex);
    pthread_cond_destroy(root->finished);

    free(root->mutex);
    free(root->finished);
    free((void*) root->thread_count);

    search_info_free(root);
}
