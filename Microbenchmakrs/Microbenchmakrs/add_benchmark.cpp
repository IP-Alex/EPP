//
//  add_benchmark.c
//  Add Microbenchmark Tutorial
//
//  Created by Ricardo Alves on 16/10/14.


#include <stdio.h>
#include <omp.h>
#include <time.h>
#include "gettimeofdayWin.h"
#include <intrin.h>
#include "opencl_helper.h"

#define NUMBER_OF_RUNS 10
#define TESTS_PER_RUN 1000*1000*1000*1000*1000L

struct timeval start_time, end_time;

void add_microbenchmark1(char scale);
void add_microbenchmark2(char scale);
void add_microbenchmark3(char scale);
void add_microbenchmark4(char scale);
void add_microbenchmark5(char scale);
void add_microbenchmark6(char scale);
void add_microbenchmark7(char scale);
void add_microbenchmark8(char scale);
void add_microbenchmark9(char scale);
void add_microbenchmark10(char scale);

void pause(void){
	printf("press Enter to exit\n");
	getchar();
}


int main(int argc, const char * argv[])
{
	int benchmark=10;

	switch(benchmark){
		case 1:
			printf("== Benchmark 1 == Baseline:\n");
			add_microbenchmark1(1);
			break;
		case 2:
		    printf("== Benchmark 2 == :\n");
		    add_microbenchmark2(1);
			break;
		case 3:
		    printf("== Benchmark 3 == :\n");
		    add_microbenchmark3(1);
			break;
		case 4:
		    printf("== Benchmark 4 == :\n");
		    add_microbenchmark4(1);
			break;
		case 5:
		    printf("== Benchmark 5 == :\n");
		    add_microbenchmark5(10);
			break;
		case 6:
			printf("== Benchmark 6 == :\n");
			add_microbenchmark6(50);
			break;
		case 7:
			printf("== Benchmark 7 == :\n");
			add_microbenchmark7(50);
			break;
		case 8:
			printf("== Benchmark 8 == :\n");
			add_microbenchmark8(1);
			break;
		case 9:
			printf("== Benchmark 9 == :\n");
			add_microbenchmark9(1);
			break;
		case 10:
			printf("== Benchmark 10 == :\n");
			add_microbenchmark10(1);
			break;
		default:
			printf("Invalid benchmark!\n");
	}
	pause();
    return 0;
}



/*
 Baseline add benchmark
 */
void add_microbenchmark1(char scale) {
    
    unsigned long LOCAL_TESTS_PER_RUN = TESTS_PER_RUN*scale;
    
    char * char_data = (char*)malloc(LOCAL_TESTS_PER_RUN*sizeof(char));
    char sum = 1;
    
    for (int run=0; run<NUMBER_OF_RUNS; run++) {
        
        gettimeofday(&start_time, NULL);
        
        for (unsigned long test=0; test<LOCAL_TESTS_PER_RUN; test+=1) {
            sum +=char_data[test];
        }
        
        gettimeofday(&end_time, NULL);
        
        double time_in_sec = (end_time.tv_sec+end_time.tv_usec/1000000.0) - (start_time.tv_sec+start_time.tv_usec/1000000.0);
        double GOPS = (LOCAL_TESTS_PER_RUN/time_in_sec)/1000000000;
        printf("Completed %ld adds in %g seconds for %g GOPS.\n", LOCAL_TESTS_PER_RUN, time_in_sec, GOPS);
    }
    
    free(char_data);
}

/*
	add benchmark 2
 */
void add_microbenchmark2(char scale) {
    
    unsigned long LOCAL_TESTS_PER_RUN = TESTS_PER_RUN*scale;
    
    char * char_data = (char*)malloc(LOCAL_TESTS_PER_RUN*sizeof(char));
    char sum = 1;
    
    for (int run=0; run<NUMBER_OF_RUNS; run++) {
        
        gettimeofday(&start_time, NULL);
        
        for (unsigned long test=0; test<LOCAL_TESTS_PER_RUN; test+=1) {
            sum +=char_data[test];
        }
        
        gettimeofday(&end_time, NULL);
        
        double time_in_sec = (end_time.tv_sec+end_time.tv_usec/1000000.0) - (start_time.tv_sec+start_time.tv_usec/1000000.0);
        double GOPS = (LOCAL_TESTS_PER_RUN/time_in_sec)/1000000000;
        printf("sum: %d  ", sum);
        printf("Completed %ld adds in %g seconds for %g GOPS.\n", LOCAL_TESTS_PER_RUN, time_in_sec, GOPS);
    }
    
    free(char_data);
}



/*
	add benchmark 3
 */
void add_microbenchmark3(char scale) {
    
    unsigned long LOCAL_TESTS_PER_RUN = TESTS_PER_RUN*scale;
    
    char sum = 1;
    
    for (unsigned long run=0; run<NUMBER_OF_RUNS; run++) {
        
        gettimeofday(&start_time, NULL);
        
        for (int test=0; test<LOCAL_TESTS_PER_RUN; test+=1) {
            sum += scale;
        }
        
        gettimeofday(&end_time, NULL);
        
        double time_in_sec = (end_time.tv_sec+end_time.tv_usec/1000000.0) - (start_time.tv_sec+start_time.tv_usec/1000000.0);
        double GOPS = (LOCAL_TESTS_PER_RUN/time_in_sec)/1000000000;
        printf("sum: %d  ", sum);
        printf("Completed %ld adds in %g seconds for %g GOPS.\n", LOCAL_TESTS_PER_RUN, time_in_sec, GOPS);
    }
}

/*
	add benchmark 4
 */
void add_microbenchmark4(char scale) {
    
    unsigned long LOCAL_TESTS_PER_RUN = TESTS_PER_RUN*scale;

    char sum = 1;
    
    for (unsigned long run=0; run<NUMBER_OF_RUNS; run++) {
        
        gettimeofday(&start_time, NULL);
        
		__asm mov al, sum; //move sum to the lower significant 8 bits of the 16 bits wide %ax register
        for (int test=0; test<LOCAL_TESTS_PER_RUN; test+=1) {
			__asm add al, scale; //add resgister %al content with scale variable and storing the result in %al
            /*
             See http://locklessinc.com/articles/gcc_asm/
             Format: Instruction : constraint for %0 (variable) : constraint for %1 (variable)
             Constraints: +r input, r register
            */
        }
		__asm mov al, sum; //moving the result from %al back to the variable sum
        
        gettimeofday(&end_time, NULL);
        
        double time_in_sec = (end_time.tv_sec+end_time.tv_usec/1000000.0) - (start_time.tv_sec+start_time.tv_usec/1000000.0);
        double GOPS = (LOCAL_TESTS_PER_RUN/time_in_sec)/1000000000;
        printf("sum: %d  ", sum);
        printf("Completed %ld adds in %g seconds for %g GOPS.\n", LOCAL_TESTS_PER_RUN, time_in_sec, GOPS);
    }
}

/*
	add benchmark 5
 */
void add_microbenchmark5(char scale) {
    
    unsigned long LOCAL_TESTS_PER_RUN = TESTS_PER_RUN*scale;

    char sum = 1;

    __m128i sum_v0 = _mm_set1_epi8(0);
    __m128i addthis_v = _mm_set1_epi8(scale);

    for (int run=0; run<NUMBER_OF_RUNS; run++) {
        
        gettimeofday(&start_time, NULL);

        for (long test=0; test<LOCAL_TESTS_PER_RUN; test+=16) {
            sum_v0 = _mm_add_epi8(addthis_v, sum_v0);
        }

        gettimeofday(&end_time, NULL);

        double time_in_sec = (end_time.tv_sec+end_time.tv_usec/1000000.0) - (start_time.tv_sec+start_time.tv_usec/1000000.0);
        double GOPS = (LOCAL_TESTS_PER_RUN/time_in_sec)/1000000000;
        for (int i=0; i<16;i++)
            sum+=((char*)(&sum_v0))[i];
        printf("sum: %d  ", sum);
        printf("Completed %lu adds in %g seconds for %g GOPS.\n", LOCAL_TESTS_PER_RUN, time_in_sec, GOPS);
    }
}


/*
	add benchmark 6
 */
void add_microbenchmark6(char scale) {
    
    unsigned long LOCAL_TESTS_PER_RUN = TESTS_PER_RUN*scale;

    char sum = 1;
    
    __m128i sum_v0 = _mm_set1_epi8(0);
    __m128i sum_v1 = _mm_set1_epi8(0);
    __m128i sum_v2 = _mm_set1_epi8(0);
    __m128i sum_v3 = _mm_set1_epi8(0);
    __m128i addthis_v = _mm_set1_epi8(scale);
    
    for (int run=0; run<NUMBER_OF_RUNS; run++) {
        
        gettimeofday(&start_time, NULL);
        
        for (long test=0; test<LOCAL_TESTS_PER_RUN; test+=16*4) {
            sum_v0 = _mm_add_epi8(addthis_v, sum_v0);
            sum_v1 = _mm_add_epi8(addthis_v, sum_v1);
            sum_v2 = _mm_add_epi8(addthis_v, sum_v2);
            sum_v3 = _mm_add_epi8(addthis_v, sum_v3);
        }
        
        gettimeofday(&end_time, NULL);
        
        double time_in_sec = (end_time.tv_sec+end_time.tv_usec/1000000.0) - (start_time.tv_sec+start_time.tv_usec/1000000.0);
        double GOPS = (LOCAL_TESTS_PER_RUN/time_in_sec)/1000000000;
        for (int i=0; i<16;i++)
            sum+=((char*)(&sum_v0))[i]+((char*)(&sum_v1))[i]+((char*)(&sum_v2))[i]+((char*)(&sum_v3))[i];
        printf("sum: %d\n", sum);
        printf("Completed %lu adds in %g seconds for %g GOPS.\n", LOCAL_TESTS_PER_RUN, time_in_sec, GOPS);
    }
}


char reduce(__m128i *sum, int array_size=4, int sse_width = 16) {
	char result = 0;
	for (int i = 0; i < array_size; i++) {
		for (int j = 0; j < sse_width; j++) 
			result += ((char*)(&sum[i]))[j];
	}
	return result;
}

/*
	add benchmark 7
 */
void add_microbenchmark7(char scale) {
	//TODO :)
	unsigned long LOCAL_TESTS_PER_RUN = TESTS_PER_RUN * scale;

	char sum = 1;
	const int num_threads = 4;
	int tid = 0;

	__m128i sum_v0[num_threads] = { _mm_set1_epi8(0) };
	__m128i sum_v1[num_threads] = { _mm_set1_epi8(0) };
	__m128i sum_v2[num_threads] = { _mm_set1_epi8(0) };
	__m128i sum_v3[num_threads] = { _mm_set1_epi8(0) };
	__m128i addthis_v =  _mm_set1_epi8(scale);

	printf("num_threads %d\n", num_threads);

#pragma omp parallel num_threads(4)
	printf("Number of threads %d\n", omp_get_num_threads());

	for (int run = 0; run<NUMBER_OF_RUNS; run++) {

		gettimeofday(&start_time, NULL);

#pragma omp parallel private(tid) 
		{
			tid = omp_get_thread_num();
#pragma omp for
			for (long test = 0; test < LOCAL_TESTS_PER_RUN; test += 16 * 4) {
				sum_v0[tid] = _mm_add_epi8(addthis_v, sum_v0[tid]);
				sum_v1[tid] = _mm_add_epi8(addthis_v, sum_v1[tid]);
				sum_v2[tid] = _mm_add_epi8(addthis_v, sum_v2[tid]);
				sum_v3[tid] = _mm_add_epi8(addthis_v, sum_v3[tid]);
			}
		}

		gettimeofday(&end_time, NULL);

		double time_in_sec = (end_time.tv_sec + end_time.tv_usec / 1000000.0) - (start_time.tv_sec + start_time.tv_usec / 1000000.0);
		double GOPS = (LOCAL_TESTS_PER_RUN / time_in_sec) / 1000000000;
		sum += reduce(sum_v0, num_threads) + reduce(sum_v1, num_threads) + reduce(sum_v2, num_threads) + reduce(sum_v3, num_threads);
		printf("sum: %d\n", sum);
		printf("Completed %lu adds in %g seconds for %g GOPS.\n", LOCAL_TESTS_PER_RUN, time_in_sec, GOPS);
	}
}

void add_microbenchmark8(char scalar){ ; }
void add_microbenchmark9(char scalar) { ; }


void add_microbenchmark10(char scalar) {
	cl_device_id device;
	cl_context context;
	cl_command_queue queue;
	cl_kernel kernel;
	cl_program program_cl;

	const int total_threads = 96 * 32 * 2;
	char* num_res = (char*)malloc(sizeof(char)*total_threads);
	memset(num_res, 0, total_threads);

	setup_cl(&device, &context, &queue);
	compile_kernel(&device, &context, &kernel, &program_cl);

	cl_mem num = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(char)*total_threads, NULL, NULL);
	clEnqueueWriteBuffer(queue, num, CL_FALSE, 0, sizeof(char)*total_threads, num_res, 0, NULL, NULL);

	int i = 1;
	cl_mem num2 = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(char), NULL, NULL);
	clEnqueueWriteBuffer(queue, num2, CL_FALSE, 0, sizeof(char), &i, 0, NULL, NULL);


	cl_int error;

	error = clSetKernelArg(kernel, 0, sizeof(num), &num);
	if (error != CL_SUCCESS) {
		printf("[ERROR] - Build kernel failed. Error: %d", error);
		exit(-1);
	}

	error = clSetKernelArg(kernel, 1, sizeof(num2), &num2);
	if (error != CL_SUCCESS) {
		printf("[ERROR] - Build kernel failed. Error: %d", error);
		exit(-1);
	}

	size_t global_dimension[] = { total_threads, 0 ,0 };

	//TIME!
	clFinish(queue);

	double start_time = omp_get_wtime();
	clEnqueueNDRangeKernel(queue, kernel, 1, NULL, global_dimension, NULL, 0, NULL, NULL);
	clFinish(queue);
	double total_time = omp_get_wtime() - start_time;
	//TIME AGAIN!!

	double operations = (total_threads * 1000);// / total_time;
	printf("Total time: %lf - Gops: %lf\n", total_time, Gops);

	
	clEnqueueReadBuffer(queue, num, CL_FALSE, 0, sizeof(char)*total_threads, num_res, 0, NULL, NULL);
	clFinish(queue);

	printf("num_res: %d\n", num_res[0]);
	cleanup(&kernel, &program_cl);

	free(num_res);
}