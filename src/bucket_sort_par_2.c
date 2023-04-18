#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>

typedef int bool;
#define true 1
#define false 0

#ifndef N
#define N 10
#endif /*N*/

#ifndef B
#define B 2
#endif /*B*/

#ifndef P
#define P 2
#endif /*P*/

#ifndef T
#define T 1
#endif /*T*/


float* init_array(int n) {
    float* result = (float*) malloc(n * sizeof(float));
    for (int i = 0; i < n; i++){
        result[i] = 0;
    }
    return result;
}

float** init_2d_array(int n, int m) {
    float** result = (float**) malloc(n * sizeof(float*));
    for (int i = 0; i < n; i++){
        result[i] = (float*) malloc(m * sizeof(float));
    }
    return result;
}

omp_lock_t* init_bucket_locks() {
    omp_lock_t* locks = (omp_lock_t*) malloc(B * sizeof(omp_lock_t));
    for (int i = 0; i < B; i++){
        omp_init_lock(&locks[i]);
    }
    return locks;
}

void free_bucket_locks(omp_lock_t* bucket_locks) {
    for (int i = 0; i < B; i++){
        omp_destroy_lock(&bucket_locks[i]);
    }
    free(bucket_locks);
}

void free_2d_array(float** arr, int n) {
    for (int i = 0; i < n; i++){
        free(arr[i]);
    }
    free(arr);
}

int* init_ints(int n) {
    int* result = (int*) malloc(n*sizeof(int));
    int i;
    for (i = 0; i < n; i++) {
        result[i] = 0;
    }
    return result;
}

int asc(const void * a, const void * b) {
    float va = *(const float*) a;
    float vb = *(const float*) b;
    return (va > vb) - (va < vb);
}

void fill_array_randomly(float* arr, int n, int tid) {
    int i;
    unsigned short xi[3];
    time_t sec = time(NULL);

    xi[0] = tid + sec;
    xi[1] = tid + 2 + sec;
    xi[2] = tid + 1 + sec;

#pragma omp for private(i) schedule(guided)
    for (i = 0; i < n; i++){
        arr[i] = erand48(xi);
    }
}

void distribute_to_buckets(const float* arr, float** buckets, int* bucket_ind,  omp_lock_t* bucket_locks) {
    int i;
    int number_of_elements_in_bucket;
#pragma omp for private(i, number_of_elements_in_bucket) schedule(guided)
    for (i = 0; i < N; i++) {
        int number_of_bucket = (int) B*arr[i];

        omp_set_lock(&bucket_locks[number_of_bucket]);
        number_of_elements_in_bucket = bucket_ind[number_of_bucket];
        buckets[number_of_bucket][number_of_elements_in_bucket] = arr[i];
        bucket_ind[number_of_bucket]++;
        omp_unset_lock(&bucket_locks[number_of_bucket]);
    }
}

void sort_buckets(float** buckets, int* bucket_ind) {
    int i;
#pragma omp for private(i) schedule(guided)
    for (i = 0; i < B; i++){
        qsort(buckets[i], bucket_ind[i], sizeof(*(buckets[i])), asc);
    }
}

void merge_buckets(float* arr, float** buckets, const int* bucket_ind, const int* cummulative_sum_arr) {
    int i, j, arr_i;
#pragma omp for private(i, j, arr_i) schedule(guided)
    for (i=0; i < B; i++) {
        arr_i = i == 0 ? 0 : cummulative_sum_arr[i-1];
        for (j=0; j < bucket_ind[i]; j++) {
            arr[arr_i] = buckets[i][j];
            arr_i ++;
        }
    }
}

bool is_sorted(const float* arr, int n) {
    for (int i = 1; i < n; i++){
        if (arr[i-1] > arr[i]){
            return false;
        }
    }
    return true;
}

double avg(const double* arr, int n) {
    double sum = 0;
    for (int i = 0; i < n; i++){
        sum += arr[i];
    }
    return sum / (float) n;
}


int* calculate_cumulative_sum(const int* arr, int n, int *cummulative_sum_arr) {
    int i;
    cummulative_sum_arr[0] = arr[0];
    for (i = 1; i < n; i++){
        cummulative_sum_arr[i] = cummulative_sum_arr[i-1] + arr[i];
    }
    return cummulative_sum_arr;
}

int main() {
    omp_set_num_threads(P);

    float* arr = init_array(N);
    float** buckets = init_2d_array(B, 3*(N/B));
    int* bucket_ind = init_ints(B);
    omp_lock_t* bucket_locks = init_bucket_locks();
    int *cummulative_sum_arr = init_ints(B);

    double* t_a_arr = malloc(P*sizeof(double));
    double* t_b_arr = malloc(P*sizeof(double));
    double* t_c_arr = malloc(P*sizeof(double));

    double t_start = omp_get_wtime();
#pragma omp parallel
    {
        int tid = omp_get_thread_num();

        fill_array_randomly(arr, N, tid);
        t_a_arr[tid] = omp_get_wtime();

        distribute_to_buckets(arr, buckets, bucket_ind, bucket_locks);
        t_b_arr[tid] = omp_get_wtime();

        sort_buckets(buckets, bucket_ind);
        t_c_arr[tid] = omp_get_wtime();

#pragma omp single
        cummulative_sum_arr = calculate_cumulative_sum(bucket_ind, B, cummulative_sum_arr);
        merge_buckets(arr, buckets, bucket_ind, cummulative_sum_arr);
    }
    double t_end = omp_get_wtime();

    bool sorted = is_sorted(arr, N);
    double t_a_delta = avg(t_a_arr, P) - t_start;
    double t_b_delta = avg(t_b_arr, P) - avg(t_a_arr, P);
    double t_c_delta = avg(t_c_arr, P) - avg(t_b_arr, P);
    double t_d_delta = t_end - avg(t_c_arr, P);
    double t_all = t_end - t_start;
    printf("%d,%d,%d,%d,%f,%f,%f,%f,%f,%s\n", N, B, P, T, t_a_delta, t_b_delta, t_c_delta, t_d_delta, t_all,
           sorted ? "true" : "false");

    free(arr);
    free(bucket_ind);
    free_2d_array(buckets, B);
    free_bucket_locks(bucket_locks);
    free(t_a_arr);
    free(t_b_arr);
    free(t_c_arr);
}
