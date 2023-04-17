#include <stdio.h>
#include <time.h>
#include <stdlib.h>


#ifndef N
#define N 1000000
#endif /*N*/

int main() {
    int tid = 1;
    unsigned short xi[3];
    float val;
    time_t sec = time(NULL);

    xi[0] = tid + sec;
    xi[1] = tid + 2 + sec;
    xi[2] = tid + 1 + sec;

    for (int i = 0; i < N; i++) {
        val = erand48(xi);
        printf("%f\n", val);
    }
}