//
//  main.cpp
//  load_benchmark
//
//  Created by Tim Svensson on 2018-09-28.
//  Copyright © 2018 Tim Svensson. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
//#include <intrin.h>

#define L1_BYTES (256e3)
#define L2_BYTES (1e6)
#define L3_BYTES (1e6)

#define NUMBER_OF_RUNS 5
#define ITERATIONS_PER_RUN 10

#define ARRAY_SIZE L3_BYTES

struct timeval start_time, end_time;

int load_32_bit_int_seq(int* array, int array_size, int sum)
{
    __asm mov al, sum;
    for (int test=0; test<array_size; test++) {
        int crnt_int = array[test];
        __asm add al, crnt_int;
    }
    __asm mov al, sum;
    
    return sum;
}

int main(int argc, const char * argv[])
{
    // loop NUMBER_OF_RUNS times
    for (int run = 0; run < NUMBER_OF_RUNS; run++)
    {
        int* array = (int*)malloc(sizeof(int)*ARRAY_SIZE);
        
        gettimeofday(&start_time, NULL);
        
        int result = -1;
        for (int i = 0; i < ITERATIONS_PER_RUN; i++)
            result = load_32_bit_int_seq(array, ARRAY_SIZE, 1);
        
        gettimeofday(&end_time, NULL);
        
        // get time taken and Gops
        double time_in_sec  = (end_time.tv_sec+end_time.tv_usec/1000000.0)
        - (start_time.tv_sec+start_time.tv_usec/1000000.0);
        double GBps           = 32 * (double)ITERATIONS_PER_RUN * ITERATIONS_PER_RUN * ARRAY_SIZE
        / time_in_sec / 1e9;
        
        // print
        printf("Run: %d\n", run);
        printf("Result: %d, ", result);
        printf("completed %d arrays in %g seconds for %g GB/s.\n", ITERATIONS_PER_RUN, time_in_sec, GBps);
        
        free(array);
    }
    
    return 0;
}
