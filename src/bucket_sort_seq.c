#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>

#ifndef N
#define N 15
#endif /*N*/

#ifndef B
#define B 5
#endif /*B*/

#ifndef T
#define T 1
#endif /*T*/


int* init_int_array(long long n){
	return (int*) malloc(n * sizeof(int));
}

float* init_array(long long n) {
	return (float*) malloc(n * sizeof(float));
}

float** init_2d_array(long long m, long long n) {
	float** result = (float**) malloc(n * sizeof(float*));
	int i;
	for (i = 0; i < m; i++){
		result[i] = (float*) malloc(m * sizeof(float));
	}
	return result;
}

void free_2d_array(float** arr, long long m) {
	int i;
	for (i = 0; i < m; i++){
		free(arr[i]);
	}
	free(arr);
}

int* init_bucket_indices(int n) {	
	int* bucket_ind = init_int_array(B);
	
	int i;
	for (i = 0; i < B; i++) {
		bucket_ind[i] = 0;
	}
	return bucket_ind;
}

int asc(const void * a, const void * b) {
	float va = *(const float*) a;
	float vb = *(const float*) b;
	return (va > vb) - (va < vb);
}

void fill_array_rand(float* arr, int n) {
	unsigned short xi[3];
	time_t t = time(NULL);

    xi[0] = t;
    xi[1] = t + 2;
    xi[2] = t + 1;

    for (int i = 0; i < n; i++){
        arr[i] = erand48(xi);
    }
}

void distribute_to_buckets(float* arr, float** buckets, int* bucket_ind, int n, int b){
	int i;
	int number_of_bucket;
	int number_of_elements_in_bucket;
	for (i = 0; i < n; i++){
		number_of_bucket = (int) b*arr[i];
		number_of_elements_in_bucket = bucket_ind[number_of_bucket];
		buckets[number_of_bucket][number_of_elements_in_bucket] = arr[i];
		bucket_ind[number_of_bucket]++;		
	}

}

void sort_buckets(float** buckets, int* bucket_ind, int b) {
	int i;
	for (i = 0; i < b; i++){
		qsort(buckets[i], bucket_ind[i], sizeof(int), asc);
	}
}

void merge_buckets(float* arr, float** buckets, int* bucket_ind, int n){	
	int j;
	int i = 0;
	int bucket_num = 0;
	while (i < n) {
		for (j = 0; j < bucket_ind[bucket_num]; j++) {
			arr[i] = buckets[bucket_num][j];
			i++;;
		}
		bucket_num++;
	}	
}

void print_arr(float* arr, int n){
	int i;
	for (i = 0; i < N; i++) {
		printf("%4f, ", arr[i]);
	} 
	printf("\n");
}

int is_sorted(const float* arr, int n) {
    for (int i = 1; i < n; i++){
        if (arr[i-1] > arr[i]){
            return 0;
        }
    }
    return 1;
}

int main() {
	float* arr = init_array(N);
	float** buckets = init_2d_array(B, 10*(N/B));
	int* bucket_ind = init_bucket_indices(B);

	double t_start = omp_get_wtime();
	fill_array_rand(arr, N);
    double t_a = omp_get_wtime();

	distribute_to_buckets(arr, buckets, bucket_ind, N, B);
    double t_b = omp_get_wtime();

	sort_buckets(buckets, bucket_ind, B);
    double t_c = omp_get_wtime();

	merge_buckets(arr, buckets, bucket_ind, N);
	double t_end = omp_get_wtime();

    int sorted = is_sorted(arr, N);
    double t_a_delta = t_a - t_start;
    double t_b_delta = t_b - t_a;
    double t_c_delta = t_c - t_b;
    double t_d_delta = t_end - t_c;
    double t_all = t_end - t_start;
    printf("%d,%d,%d,%d,%f,%f,%f,%f,%f,%s\n", N, B, 1, T, t_a_delta, t_b_delta, t_c_delta, t_d_delta, t_all,
           sorted ? "true" : "false");

    free(arr);
	free(bucket_ind);
	free_2d_array(buckets, B);
	return 0;
}
