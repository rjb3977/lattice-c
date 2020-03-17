#define _POSIX_C_SOURCE 200809L

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <gmp.h>
#include <pthread.h>

#include "la.h"
#include "lp.h"
#include "enumerate.h"

typedef struct {
    long dimensions;
    long depth;

    const mpq_t* transform; // mpq_t[dimensions][2 * dimensions]
    const mpq_t* offset;    // mpq_t[dimensions]
    mpq_t* fixed;           // mpq_t[dimensions]

    mpq_t *A;               // mpq_t[2 * dimensions][3 * dimensions]
    mpq_t *b;               // mpq_t[2 * dimensions]
    mpq_t *c;               // mpq_t[2][3 * dimensions]
    mpq_t *λ;               // mpq_t[2 * dimensions]
    mpq_t *s;               // mpq_t[dimensions]
    mpq_t *d;               // mpq_t[2 * dimensions]
    mpq_t *x;               // mpq_t[2 * dimensions]

    long *B;                // long[2 * dimensions]
    long *N;                // long[dimensions]

    mpq_t *lu;              // mpq_t[2 * dimensions][2 * dimensions]
    long *pivots;           // long[2 * dimensions]

    long corrections;       //
    long *indices;          // long[corrections]
    mpq_t (*columns);       // mpq_t[corrections][?]

    long *count_out;        // long[1]
    mpq_t (**results_out);  // mpq_t[*count_out][dimensions]

    pthread_mutex_t *mutex; // pthread_mutex_t[1]
    pthread_cond_t *finished;
    volatile long* thread_count;
    long thread_max;
} search_info_t;

static mpq_t *vec_malloc(long length) {
    mpq_t *dest = malloc(length * sizeof(mpq_t));
    vec_init(length, dest);

    return dest;
}

static mpq_t *vec_dup(long length, const mpq_t *src) {
    mpq_t *dest = vec_malloc(length);
    vec_set(length, dest, src);

    return dest;
}

static void vec_free(long length, mpq_t *src) {
    vec_clear(length, src);
    free(src);
}

static mpq_t *mat_malloc(long rows, long cols) {
    mpq_t *dest = malloc(rows * cols * sizeof(mpq_t));
    mat_init(rows, cols, (void*) dest);

    return dest;
}

static mpq_t *mat_dup(long rows, long cols, const mpq_t *src) {
    mpq_t *dest = mat_malloc(rows, cols);
    mat_set(rows, cols, (void*) dest, (const void*) src);

    return dest;
}

static void mat_free(long rows, long cols, mpq_t *src) {
    mat_clear(rows, cols, (void*) src);
    free(src);
}

static search_info_t *search_info_dup(const search_info_t* src) {
    search_info_t* dest = malloc(sizeof(search_info_t));

    dest->dimensions = src->dimensions;
    dest->depth = src->depth;

    dest->transform = src->transform;
    dest->offset = src->offset;
    dest->fixed = vec_dup(src->dimensions, src->fixed);

    dest->A = mat_dup(2 * src->dimensions, 3 * src->dimensions, src->A);
    dest->b = vec_dup(2 * src->dimensions, src->b);
    dest->c = mat_dup(2, 3 * src->dimensions, src->c);
    dest->λ = vec_malloc(2 * src->dimensions);
    dest->s = vec_malloc(src->dimensions);
    dest->d = vec_malloc(2 * src->dimensions);
    dest->x = vec_malloc(2 * src->dimensions);

    dest->B = malloc(2 * src->dimensions * sizeof(long));
    dest->N = malloc(src->dimensions * sizeof(long));

    dest->lu = mat_malloc(2 * src->dimensions, 2 * src->dimensions);
    dest->pivots = malloc(2 * src->dimensions * sizeof(long));

    dest->corrections = 0;
    dest->indices = NULL;
    dest->columns = NULL;

    dest->count_out = src->count_out;
    dest->results_out = src->results_out;

    dest->mutex = src->mutex;
    dest->finished = src->finished;
    dest->thread_count = src->thread_count;
    dest->thread_max = src->thread_max;

    return dest;
}

static void search_info_free(search_info_t* src) {
    vec_free(src->dimensions, src->fixed);

    mat_free(2 * src->dimensions, 3 * src->dimensions, src->A);
    vec_free(2 * src->dimensions, src->b);
    mat_free(2, 3 * src->dimensions, src->c);
    vec_free(2 * src->dimensions, src->λ);
    vec_free(src->dimensions, src->s);
    vec_free(2 * src->dimensions, src->d);
    vec_free(2 * src->dimensions, src->x);

    free(src->B);
    free(src->N);

    mat_free(2 * src->dimensions, 2 * src->dimensions, src->lu);
    free(src->pivots);

    free(src->indices);
    free(src->columns);

    free(src);
}

        // search_info_t* next = search_info_dup(start);

        // if (mpq_sgn(temp) >= 0) {
        //     vec_set(next->dimensions, next->A + (next->dimensions + next->depth) * 3 * next->dimensions, next->transform + next->depth * 2 * next->dimensions);
        //     mpq_set(next->b[start->dimensions + next->depth], temp);
        // } else {
        //     vec_neg(next->dimensions, next->A + (next->dimensions + next->depth) * 3 * next->dimensions, next->transform + next->depth * 2 * next->dimensions);
        //     mpq_neg(next->b[next->dimensions + next->depth], temp);
        // }

        // mpq_set_z(next->fixed[next->depth], min);
        // next->depth += 1;

        // *count += 1;
        // *nexts = realloc(*nexts, *count * sizeof(search_info_t));
        // (*nexts)[*count - 1] = next;

static void* search_thread(void*);

static void search(search_info_t *info) {
    if (info->depth == info->dimensions) {
        pthread_mutex_lock(info->mutex);

        printf("found: ");
        vec_print(info->dimensions, info->fixed, stdout);
        printf("\n");

        *info->count_out += 1;
        *info->results_out = realloc(*info->results_out, *info->count_out * sizeof(mpq_t[info->dimensions]));

        vec_init(info->dimensions, *info->results_out + (*info->count_out - 1) * info->dimensions);
        vec_set(info->dimensions, *info->results_out + (*info->count_out - 1) * info->dimensions, info->fixed);

        pthread_mutex_unlock(info->mutex);
    } else {
        mpq_t temp;
        mpz_t min;
        mpz_t max;

        mpq_init(temp);
        mpz_init(min);
        mpz_init(max);

        // lower bound
        vec_neg(info->dimensions, info->c, info->transform + info->depth * 2 * info->dimensions);
        simplex_solve(info->dimensions, info->depth, (void *)info->A, info->b, (void *)info->c, info->λ, info->s, info->d, info->x, info->B, info->N, (void *)info->lu, info->pivots, &info->corrections, &info->indices, (void *) &info->columns);

        vec_dot(2 * info->dimensions, temp, info->transform + info->depth * 2 * info->dimensions, info->x);
        mpq_sub(temp, info->offset[info->depth], temp);
        mpz_cdiv_q(min, mpq_numref(temp), mpq_denref(temp));

        // upper bound
        vec_set(info->dimensions, info->c, info->transform + info->depth * 2 * info->dimensions);
        simplex_solve(info->dimensions, info->depth, (void *)info->A, info->b, (void *)info->c, info->λ, info->s, info->d, info->x, info->B, info->N, (void *)info->lu, info->pivots, &info->corrections, &info->indices, (void *) &info->columns);

        vec_dot(2 * info->dimensions, temp, info->transform + info->depth * 2 * info->dimensions, info->x);
        mpq_sub(temp, info->offset[info->depth], temp);
        mpz_fdiv_q(max, mpq_numref(temp), mpq_denref(temp));

        // rhs of new row of A
        mpq_set_z(temp, min);
        mpq_sub(temp, info->offset[info->depth], temp);

        // min <= max, min += 1, rhs -= 1
        for (; mpz_cmp(min, max) <= 0; mpz_add_ui(min, min, 1), mpz_sub(mpq_numref(temp), mpq_numref(temp), mpq_denref(temp))) {
            if (mpq_sgn(temp) >= 0) {
                vec_set(info->dimensions, info->A + (info->dimensions + info->depth) * 3 * info->dimensions, info->transform + info->depth * 2 * info->dimensions);
                mpq_set(info->b[info->dimensions + info->depth], temp);
            } else {
                vec_neg(info->dimensions, info->A + (info->dimensions + info->depth) * 3 * info->dimensions, info->transform + info->depth * 2 * info->dimensions);
                mpq_neg(info->b[info->dimensions + info->depth], temp);
            }

            mpq_set_z(info->fixed[info->depth], min);
            info->depth += 1;

            pthread_mutex_lock(info->mutex);

            if (*info->thread_count < info->thread_max) {
                *info->thread_count += 1;
                pthread_mutex_unlock(info->mutex);

                pthread_t thread;
                pthread_create(&thread, NULL, search_thread, search_info_dup(info));
            } else {
                pthread_mutex_unlock(info->mutex);
                search(info);
            }

            info->depth -= 1;
        }

        mpq_clear(temp);
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
    pthread_mutex_lock(root->mutex);

    while (*root->thread_count > 0) {
        pthread_cond_wait(root->finished, root->mutex);
    }
}

void enumerate(long dimensions, mpq_t basis[dimensions][dimensions], mpq_t lower[dimensions], mpq_t upper[dimensions], long *count_out, mpq_t (**results_out)[dimensions], long thread_max) {
    search_info_t* root = malloc(sizeof(search_info_t));

    root->dimensions = dimensions;
    root->depth = 0;

    //root->transform = mat_malloc(dimensions, 2 * dimensions);
    //root->offset = vec_malloc(dimensions);
    root->fixed = vec_malloc(dimensions);

    root->A = mat_malloc(2 * dimensions, 3 * dimensions);
    root->b = vec_malloc(2 * dimensions);
    root->c = mat_malloc(2, 3 * dimensions);
    root->λ = vec_malloc(2 * dimensions);
    root->s = vec_malloc(dimensions);
    root->d = vec_malloc(2 * dimensions);
    root->x = vec_malloc(2 * dimensions);

    root->B = malloc(2 * dimensions * sizeof(long));
    root->N = malloc(dimensions * sizeof(long));

    root->lu = mat_malloc(2 * dimensions, 2 * dimensions);
    root->pivots = malloc(2 * dimensions * sizeof(long));
    root->corrections = 0;
    root->indices = NULL;
    root->columns = NULL;

    root->count_out = count_out;
    root->results_out = (void*) results_out;

    root->mutex = malloc(sizeof(pthread_mutex_t));
    root->finished = malloc(sizeof(pthread_cond_t));
    root->thread_count = malloc(sizeof(long));
    root->thread_max = thread_max;

    pthread_mutex_init(root->mutex, NULL);
    pthread_cond_init(root->finished, NULL);

    mpq_t *transform = mat_malloc(dimensions, 2 * dimensions);
    mpq_t *offset = vec_malloc(dimensions);

    root->transform = transform;
    root->offset = offset;

    mat_copy(dimensions, dimensions, 0, 0, 2 * dimensions, (void*) root->lu, 0, 0, dimensions, basis);
    mat_lu_2(dimensions, 2 * dimensions, (void*) root->lu, root->pivots);

    // inverse * identity
    for (long i = 0; i < dimensions; ++i) {
        mpq_set_ui(transform[i * 2 * dimensions + i], 1, 1);
        solve_ptlu_2(dimensions, 2 * dimensions, (void*) root->lu, root->pivots, transform + i * 2 * dimensions);
    }

    // inverse * upper
    vec_set(dimensions, offset, upper);
    solve_utltp_2(dimensions, 2 * dimensions, (void*) root->lu, root->pivots, offset);

    // A = [I|I|0]
    //     [0|0|I]
    // c = [.|0|0]
    //     [0|0|1]
    for (long i = 0; i < dimensions; ++i) {
        mpq_set_ui(root->A[i * 3 * dimensions + i], 1, 1);
        mpq_set_ui(root->A[i * 3 * dimensions + dimensions + i], 1, 1);
        mpq_set_ui(root->A[(dimensions + i) * 3 * dimensions + 2 * dimensions + i], 1, 1);
        mpq_set_ui(root->c[3 * dimensions + 2 * dimensions + i], 1, 1);
    }

    vec_sub(dimensions, root->b, upper, lower);

    search_parallel(root);

    mat_free(dimensions, 2 * dimensions, transform);
    vec_free(dimensions, offset);

    pthread_mutex_destroy(root->mutex);
    pthread_cond_destroy(root->finished);

    free(root->mutex);
    free(root->finished);
    free((void*) root->thread_count);

    search_info_free(root);
}

// static void search_old(long dimensions, long depth, mpq_t transform[dimensions][2 * dimensions], mpq_t offset[dimensions], mpq_t fixed[dimensions], mpq_t A[2 * dimensions][3 * dimensions], mpq_t b[2 * dimensions], mpq_t c[2][3 * dimensions], mpq_t λ[2 * dimensions], mpq_t s[dimensions], mpq_t d[2 * dimensions], mpq_t x[2 * dimensions], long B[2 * dimensions], long N[dimensions], mpq_t lu[2 * dimensions][2 * dimensions], long pivots[2 * dimensions], long *corrections, long **indices, mpq_t(**columns), long *count_out, mpq_t (**results_out)[dimensions]) {
//     if (depth == dimensions) {
//         printf("found: ");
//         vec_print(dimensions, fixed, stdout);
//         printf("\n");

//         *count_out += 1;
//         *results_out = realloc(*results_out, *count_out * sizeof(mpq_t[dimensions]));

//         vec_init(dimensions, (*results_out)[*count_out - 1]);
//         vec_set(dimensions, (*results_out)[*count_out - 1], fixed);
//     } else {
//         mpq_t temp;
//         mpz_t min;
//         mpz_t max;

//         mpq_init(temp);
//         mpz_init(min);
//         mpz_init(max);

//         // lower bound
//         vec_neg(dimensions, c[0], transform[depth]);
//         simplex_solve(dimensions, depth, A, b, c, λ, s, d, x, B, N, lu, pivots, corrections, indices, (void*) columns);

//         vec_dot(2 * dimensions, temp, transform[depth], x);
//         mpq_sub(temp, offset[depth], temp);
//         mpz_cdiv_q(min, mpq_numref(temp), mpq_denref(temp));

//         // upper bound
//         vec_set(dimensions, c[0], transform[depth]);
//         simplex_solve(dimensions, depth, A, b, c, λ, s, d, x, B, N, lu, pivots, corrections, indices, (void*) columns);

//         vec_dot(2 * dimensions, temp, transform[depth], x);
//         mpq_sub(temp, offset[depth], temp);
//         mpz_fdiv_q(max, mpq_numref(temp), mpq_denref(temp));

//         // rhs of new row of A
//         mpq_set_z(temp, min);
//         mpq_sub(temp, offset[depth], temp);

//         // printf("%ld: ", depth);
//         // mpz_out_str(stdout, 10, min);
//         // printf(" -> ");
//         // mpz_out_str(stdout, 10, max);
//         // printf("\n");

//         // min <= max, min += 1, rhs -= 1
//         for (; mpz_cmp(min, max) <= 0; mpz_add_ui(min, min, 1), mpz_sub(mpq_numref(temp), mpq_numref(temp), mpq_denref(temp))) {
//             if (mpq_sgn(temp) >= 0) {
//                 vec_set(dimensions, A[dimensions + depth], transform[depth]);
//                 mpq_set(b[dimensions + depth], temp);
//             } else {
//                 vec_neg(dimensions, A[dimensions + depth], transform[depth]);
//                 mpq_neg(b[dimensions + depth], temp);
//             }

//             mpq_set_z(fixed[depth], min);
//             search_old(dimensions, depth + 1, transform, offset, fixed, A, b, c, λ, s, d, x, B, N, lu, pivots, corrections, indices, columns, count_out, results_out);
//         }

//         mpq_clear(temp);
//         mpz_clear(min);
//         mpz_clear(max);
//     }
// }
