#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>

#ifndef N
#define N 10000000
#endif /*N*/

#ifndef B
#define B 4
#endif /*B*/

#ifndef P
#define P 4
#endif /*P*/

#ifndef T
#define T 1
#endif /*T*/

int* init_int_array(int n){
    return (int*) malloc(n * sizeof(int));
}

float* init_array(int n) {
    float* result = (float*) malloc(n * sizeof(float));
    int i;
    for (i = 0; i < n; i++){
        result[i] = 0;
    }
    return result;
}

float** init_2d_array(int n, int m) {
    float** result = (float**) malloc(n * sizeof(float*));
    int i;
    for (i = 0; i < n; i++){
        result[i] = (float*) malloc(m * sizeof(float));
    }
    return result;
}

void free_2d_array(float** arr, int n) {
    int i;
    for (i = 0; i < n; i++){
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

void fill_array_rand(float* arr, int n, int tid) {
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

void add_to_bucket(float** buckets, int* bucket_ind, float val){
    int number_of_bucket = (int) B*val;
    #pragma omp critical
    {
        int number_of_elements_in_bucket = bucket_ind[number_of_bucket];
        buckets[number_of_bucket][number_of_elements_in_bucket] = val;
        bucket_ind[number_of_bucket]++;
    }
}


int compute_buck_start_ind(int thr){
    float step = (float) B / (float) P;
    int result = (int) ((float) thr * step);
    return result;
}

int compute_buck_end_ind(int thr){
    if (thr == P - 1) return B-1;

    return compute_buck_start_ind(thr + 1) - 1;
}


void distribute_to_buckets(float* arr, float** buckets, int* bucket_ind){
    int tid = omp_get_thread_num(); // Numer wątku
    int chunk_size = N / P; // Rozmiar fragmentu tablicy na jeden wątek
    int start = tid * chunk_size; // Początkowy indeks fragmentu tablicy dla danego wątku
    int end = start + chunk_size; // Końcowy indeks fragmentu tablicy dla danego wątku

    for (int i = start; i < end; i++) {
        add_to_bucket(buckets, bucket_ind, arr[i]);
    }
}

void sort_buckets(float** buckets, int* bucket_ind, int start_ind, int end_ind) {
    int i;
    for (i = start_ind; i <= end_ind; i++){
        qsort(buckets[i], bucket_ind[i], sizeof(*(buckets[i])), asc);
    }
}

void merge_buckets(float* arr, float** buckets, const int* bucket_ind, int bucket_start_ind, int from, int to){
    int j;
    int bucket_num = bucket_start_ind;
    int i = from;
    while (i < to) {
        for (j = 0; j < bucket_ind[bucket_num]; j++) {
            arr[i] = buckets[bucket_num][j];
            i++;
        }
        bucket_num++;
    }
}


void print_float_arr(float* arr, int n){
    int i;
    for (i = 0; i < n; i++) {
        printf("%4f, ", arr[i]);
    }
    printf("\n");
}

int is_sorted(const float* arr, int n){
    int i;
    for (i = 1; i < n; i++){
        if (arr[i-1] > arr[i]){
            return 0;
        }
    }
    return 1;
}

char* print_is_sorted(int sorted){
    return sorted ? "true" : "false";
}

double avg(const double* arr, int n){
    int i;
    double sum = 0;
    for (i = 0; i < n; i++){
        sum += arr[i];
    }
    return sum / (float) n;
}

int main() {
    omp_set_num_threads(P);

    float* arr = init_array(N);
    float** buckets = init_2d_array(B, 5*(N/B));
    int* bucket_ind = init_ints(B);

    double* t_a_arr = malloc(P*sizeof(double));
    double* t_b_arr = malloc(P*sizeof(double));
    double* t_c_arr = malloc(P*sizeof(double));

    double t_start = omp_get_wtime();
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        fill_array_rand(arr, N, tid);
        t_a_arr[tid] = omp_get_wtime();
        distribute_to_buckets(arr, buckets, bucket_ind);
        #pragma omp barrier
        t_b_arr[tid] = omp_get_wtime();
        int buck_start_ind = compute_buck_start_ind(tid);
        int buck_end_ind = compute_buck_end_ind(tid);
        sort_buckets(buckets, bucket_ind, buck_start_ind, buck_end_ind);
        #pragma omp barrier
        t_c_arr[tid] = omp_get_wtime();
        int merge_start_ind = 0;
        for (int i = 0; i < buck_start_ind; i++) {
            merge_start_ind += bucket_ind[i];
        }
        int merge_end_ind = merge_start_ind;
        for (int i = buck_start_ind; i <= buck_end_ind; i++) {
            merge_end_ind += bucket_ind[i];
        }
        merge_buckets(arr, buckets, bucket_ind, buck_start_ind, merge_start_ind, merge_end_ind);
    }
    double t_end = omp_get_wtime();

    int sorted = is_sorted(arr, N);

    double t_a_delta = avg(t_a_arr, P) - t_start;
    double t_b_delta = avg(t_b_arr, P) - avg(t_a_arr, P);
    double t_c_delta = avg(t_c_arr, P) - avg(t_b_arr, P);
    double t_d_delta = t_end - avg(t_c_arr, P);
    double t_all = t_end - t_start;
    printf("%d,%d,%d,%d,%f,%f,%f,%f,%f,%s\n", N, B, P, T, t_a_delta, t_b_delta, t_c_delta, t_d_delta, t_all, print_is_sorted(sorted));

    free(arr);
    free(bucket_ind);
    free_2d_array(buckets, B);
    free(t_a_arr);
    free(t_b_arr);
    free(t_c_arr);
    return 0;
}