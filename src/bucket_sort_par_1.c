#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>

#ifndef N
#define N 100 
#endif /*N*/

#ifndef B
#define B 4
#endif /*B*/

#ifndef P
#define P 2
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

	tid = omp_get_thread_num();
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
	int number_of_elements_in_bucket = bucket_ind[number_of_bucket];
	buckets[number_of_bucket][number_of_elements_in_bucket] = val;	
	bucket_ind[number_of_bucket]++;
}


int compute_offset(int thr, int n){
	return  (thr * ((float) n / (float) P));
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

float compute_min_val(int buck_start_ind){
	return  (float) buck_start_ind / (float) B; 
}

float compute_max_val(int buck_end_ind){
	return buck_end_ind != B ? (float) (buck_end_ind + 1) / (float) B : 1.0;
}
void distribute_to_buckets(float* arr, float** buckets, int* bucket_ind, float min_val, float max_val, int offset){
	int i;
	float val;
	int ind;
	for (i = 0; i < N; i++) {
		ind = (i + offset) % N;
		val = arr[ind];
		
		if ( min_val <= val && val < max_val){ 
			add_to_bucket(buckets, bucket_ind, val);
		}	
	} 
}

void sort_buckets(float** buckets, int* bucket_ind, int start_ind, int end_ind) {
	int i;
	for (i = start_ind; i <= end_ind; i++){
		qsort(buckets[i], bucket_ind[i], sizeof(int), asc);
	}
}

void merge_buckets(float* arr, float** buckets, int* bucket_ind, int bucket_start_ind, int from, int to){	
	int j;
	int bucket_num = bucket_start_ind;
	int i = from;
	while (i < to) {
		for (j = 0; j < bucket_ind[bucket_num]; j++) {
			arr[i] = buckets[bucket_num][j];
			i++;;
		}
		bucket_num++;
	}	
}

int* compute_cummulative_sum(int* arr, int* cum_arr, int n){
	int i;
	cum_arr[0] = arr[0];
	for (i = 1; i < n; i++){
		cum_arr[i] = cum_arr[i-1] + arr[i];
	}		
	return cum_arr;
}
 
void print_int_arr(int* arr, int n) {
	int i;
	for (i = 0; i < n; i++) {
		printf("%d, ", arr[i]);
	}
	printf("\n");
}


void print_float_arr(float* arr, int n){
	int i;
	for (i = 0; i < n; i++) {
		printf("%4f, ", arr[i]);
	} 
	printf("\n");
}

int is_sorted(float* arr, int n){
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

double avg(double* arr, int n){
	int i;
	double sum = 0;
	for (i = 0; i < n; i++){
		sum += arr[i];
	}
	return sum / (float) n; 
}

int main() {
	omp_set_num_threads(P);

	int buck_len = 10*(N/B);

	float* arr = init_array(N);
	float** buckets = init_2d_array(B, buck_len);
	int* bucket_ind = init_ints(B);
	int* cum_buck_ind = init_ints(B);
	double* t_a_arr = malloc(P*sizeof(double));
	double* t_b_arr = malloc(P*sizeof(double));
	double* t_c_arr = malloc(P*sizeof(double));
	
	double t_start = omp_get_wtime();	
	#pragma omp parallel
	{
	int thr = omp_get_thread_num();
	fill_array_rand(arr, N, thr);
	t_a_arr[thr] = omp_get_wtime(); 

	int offset = compute_offset(thr, N);
	int buck_start_ind = compute_buck_start_ind(thr); 
	int buck_end_ind = compute_buck_end_ind(thr);
	float min_val = compute_min_val(buck_start_ind);
	float max_val = compute_max_val(buck_end_ind);
	distribute_to_buckets(arr, buckets, bucket_ind, min_val, max_val, offset);
	t_b_arr[thr] = omp_get_wtime();	

	sort_buckets(buckets, bucket_ind, buck_start_ind, buck_end_ind);  
	t_c_arr[thr] = omp_get_wtime();

	#pragma omp barrier
	compute_cummulative_sum(bucket_ind, cum_buck_ind, B);
	int merge_start_ind = thr == 0 ? 0 : cum_buck_ind[buck_start_ind - 1]; //inclusive
	int merge_end_ind = cum_buck_ind[buck_end_ind];  //exclusive

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
	free(cum_buck_ind);
	free_2d_array(buckets, B);
	free(t_a_arr);
	free(t_b_arr);
	free(t_c_arr);
	return 0;
}
