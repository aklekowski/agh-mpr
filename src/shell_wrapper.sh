#!/bin/bash


echo "N,B,P,TEST_CASE_NUM,T_A,T_B,T_C,T_D,T_E,IS_SORTED"
for N in 10000 100000 1000000 10000000
do
        for B in 10 100 1000 10000
        do
                for P in {4..4}
                do
                        for test_case_num in {1..20}
                        do
                                gcc -Wall bucket_sort_par_2.c -o bucket_sort.o -DN=$N -DB=$B -DP=$P -DT=$test_case_num -fopenmp
                                ./bucket_sort.o
                        done
                done
        done
done
rm bucket_sort.o
