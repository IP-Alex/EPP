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

int load_32_bit_int_seq(int array_size, int sum)
{
    int add = 0;
    __asm mov al, sum; //move sum to the lower significant 8 bits of the 16 bits wide %ax register
    for (int test=0; test<array_size; test+=1) {
        __asm add al, add; //add resgister %al content with scale variable and storing the result in %al
        /*
         See http://locklessinc.com/articles/gcc_asm/
         Format: Instruction : constraint for %0 (variable) : constraint for %1 (variable)
         Constraints: +r input, r register
         */
    }
    __asm mov al, sum; //moving the result from %al back to the variable sum
    
    return 0;
}

int main(int argc, const char * argv[])
{
    // loop NUMBER_OF_RUNS times
    for (int run = 0; run < NUMBER_OF_RUNS; run++)
    {
        // start timer
        gettimeofday(&start_time, NULL);
        // loop ITERATIONS_PER_RUN times
        
        int result = -1;
        for (int i = 0; i < ITERATIONS_PER_RUN; i++)
            result = load_32_bit_int_seq(ARRAY_SIZE, 1);
        // stop timer
        gettimeofday(&end_time, NULL);
        
        // get time taken and Gops
        double time_in_sec  = (end_time.tv_sec+end_time.tv_usec/1000000.0)
                            - (start_time.tv_sec+start_time.tv_usec/1000000.0);
        double GB           = 32 * (double)ITERATIONS_PER_RUN * ITERATIONS_PER_RUN * ARRAY_SIZE
                            / time_in_sec / 1e9;
        
        // print
        printf("Run: %d\n", run);
        printf("Result: %d, ", result);
        printf("completed %d arrays in %g seconds for %g GB.\n", ITERATIONS_PER_RUN, time_in_sec, GB);
    }
    
    return 0;
}
