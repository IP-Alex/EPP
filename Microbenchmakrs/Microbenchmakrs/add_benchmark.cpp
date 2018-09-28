//
//  add_benchmark.c
//  Add Microbenchmark Tutorial
//
//  Created by Ricardo Alves on 16/10/14.


#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include "gettimeofdayWin.h"
#include <intrin.h>

#define NUMBER_OF_RUNS 40

#define L_MAX 2147483647L
#define L_meh 1750000000L
#define L_30 1073741824L
#define L_29 536870912L
#define L_25 33554432L
#define TESTS_PER_RUN L_MAX

#define ARRAY_SIZE 2147483640

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
//void add_microbenchmark10(char scale);
//void add_microbenchmark11(char scale);

void pause(void){
	printf("press Enter to exit\n");
	getchar();
}


int main(int argc, const char * argv[])
{
	int benchmark=9;
	printf("Number of runs: %d; Tests per run: %ld\n", NUMBER_OF_RUNS, TESTS_PER_RUN);

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
			add_microbenchmark7(1);
			break;
		case 8:
			printf("== Benchmark 8 == :\n");
			add_microbenchmark8(1);
			break;
		case 9:
			printf("== Benchmark 9 == :\n");
			add_microbenchmark9(1);
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


template<class T>
char reduce(T *sum, int array_size=4, int sse_width = 16) {
	char result = 0;
	for (int i = 0; i < array_size; i++) {
		for (int j = 0; j < sse_width; j++) 
			result += ((char*)(&sum[i]))[j];
	}
	return result;
}

template<class T>
int sum_arr(T* sum_t0, T* sum_t1, T* sum_t2, T* sum_t3, int array_size = 4, int sse_width = 16) {
	char result = 0;
	for (int i = 0; i < array_size; i++) 
		for (int j = 0; j < sse_width; j++)
			result += ((char*)(&sum_t0[i]))[j] + ((char*)(&sum_t1[i]))[j] + ((char*)(&sum_t2[i]))[j] + ((char*)(&sum_t3[i]))[j];

	
	
	return result;
}
/*
	add benchmark 7
 */
void add_microbenchmark7(char scale) {

	struct timeval _start_time, _end_time;
	unsigned long LOCAL_TESTS_PER_RUN = TESTS_PER_RUN * scale;
	char sum = 1;

	__m128i sum_v0 = _mm_set1_epi8(0);
	__m128i sum_v1 = _mm_set1_epi8(0);
	__m128i sum_v2 = _mm_set1_epi8(0);
	__m128i sum_v3 = _mm_set1_epi8(0);
	//__m128i addthis_v = _mm_set1_epi8(scale);

	// cache line 64bytes
	// load 2 lines per fetch = 128bytes
	__m128i sum_v[4][8];
	for (int four = 0; four < 4; four++) {
		for (int fourfour = 0; fourfour < 4 + 4; fourfour++) {
			sum_v[four][fourfour] = _mm_set1_epi8(1);
		}
	}
	__m128i addthis_v[4] = { _mm_set1_epi8(1),_mm_set1_epi8(1) ,_mm_set1_epi8(1) ,_mm_set1_epi8(1) };
	double start_omp_time[4], total_time;

#pragma omp parallel num_threads(4)
	{
		int j = 1;
		printf("%d", j);
	}



		//gettimeofday(&_start_time, NULL);


	for (int run = 0; run < NUMBER_OF_RUNS; run++) {

#pragma omp parallel num_threads(4)
		{
			int tid = omp_get_thread_num();
			start_omp_time[tid] = omp_get_wtime();

			//#pragma omp parallel for  (0 10 20... 40)
#pragma omp for
			for (long test1 = 0; test1 < 200; test1++) {
				for (long test = 0; test < LOCAL_TESTS_PER_RUN; test += 16 * 4) {
					sum_v[tid][0] = _mm_add_epi8(addthis_v[tid], sum_v[tid][0]);
					sum_v[tid][1] = _mm_add_epi8(addthis_v[tid], sum_v[tid][1]);
					sum_v[tid][2] = _mm_add_epi8(addthis_v[tid], sum_v[tid][2]);
					sum_v[tid][3] = _mm_add_epi8(addthis_v[tid], sum_v[tid][3]);
				}
			}
		}


		total_time = omp_get_wtime() - start_omp_time[0];
		//double time_in_sec = (_end_time.tv_sec + _end_time.tv_usec / 1000000.0) - (_start_time.tv_sec + _start_time.tv_usec / 1000000.0);

		double GOPS = 200*((LOCAL_TESTS_PER_RUN / total_time) / 1e9);
		sum += sum_arr(sum_v[0], sum_v[1], sum_v[2], sum_v[3]);

		printf("sum: %d\n", sum);
		printf("Completed %lu adds in %g seconds for %g GOPS.\n", LOCAL_TESTS_PER_RUN, total_time, GOPS);

	}

}


/*
add benchmark 8
*/
void add_microbenchmark8(char scale) {

	struct timeval _start_time, _end_time;
	unsigned long LOCAL_TESTS_PER_RUN = TESTS_PER_RUN * scale;
	char sum = 1;

	// cache line 64bytes
	// load 2 lines per fetch = 128bytes
	__m256 sum_v[4][8];
	for (int four = 0; four < 4; four++) {
		for (int fourfour = 0; fourfour < 4 + 4; fourfour++) {
			sum_v[four][fourfour] = _mm256_set_ps(1,1,1,1,1,1,1,1);
		}
	}
	__m256 addthis_v[4] = { _mm256_set_ps(1,1,1,1,1,1,1,1),_mm256_set_ps(1,1,1,1,1,1,1,1),_mm256_set_ps(1,1,1,1,1,1,1,1),_mm256_set_ps(1,1,1,1,1,1,1,1) };
	double start_omp_time[4], total_time;

#pragma omp parallel num_threads(4)
	{
		int j = 1;
		printf("%d", j);
	}

	for (int run = 0; run < NUMBER_OF_RUNS; run++) {

#pragma omp parallel num_threads(4)
		{
			int tid = omp_get_thread_num();
			start_omp_time[tid] = omp_get_wtime();

#pragma omp for
			for (long test1 = 0; test1 < 200; test1++) {
				for (long test = 0; test < LOCAL_TESTS_PER_RUN; test += 16 * 4) {
					sum_v[tid][0] = _mm256_add_ps(addthis_v[tid], sum_v[tid][0]);
					sum_v[tid][1] = _mm256_add_ps(addthis_v[tid], sum_v[tid][1]);
					sum_v[tid][2] = _mm256_add_ps(addthis_v[tid], sum_v[tid][2]);
					sum_v[tid][3] = _mm256_add_ps(addthis_v[tid], sum_v[tid][3]);
				}
			}
		}


		total_time = omp_get_wtime() - start_omp_time[0];
		//double time_in_sec = (_end_time.tv_sec + _end_time.tv_usec / 1000000.0) - (_start_time.tv_sec + _start_time.tv_usec / 1000000.0);

		double GOPS = 200 * ((LOCAL_TESTS_PER_RUN / total_time) / 1e9);
		sum += sum_arr(sum_v[0], sum_v[1], sum_v[2], sum_v[3]);

		printf("sum: %d\n", sum);
		printf("Completed %lu adds in %g seconds for %g GOPS.\n", LOCAL_TESTS_PER_RUN, total_time, GOPS);

	}

}

/*
add benchmark 9
*/
void add_microbenchmark9(char scale) {
	double start_omp_time[4], total_time;
	
	int* sum_v[4];
	int res[4] = { 0 };

	srand(omp_get_wtime());

	for (int core = 0; core < 4; core++) {
		sum_v[core] = (int*)malloc(sizeof(int)*ARRAY_SIZE);

		for (int i = 0; i < ARRAY_SIZE; i++) {
			sum_v[core][i] = 42 + rand() / RAND_MAX * 100;
		}

	}
#pragma omp parallel num_threads(4)
	{
		int j = 1;
		printf("%d", j);
	}

	for (int run = 0; run < NUMBER_OF_RUNS; run++) {

#pragma omp parallel num_threads(4)
		{
			int tid = omp_get_thread_num();
			start_omp_time[tid] = omp_get_wtime();

#pragma omp for
			for (int test = 0; test < ARRAY_SIZE; test++) 
			{
				res[tid] += sum_v[tid][test];
			}
		}


		total_time = omp_get_wtime() - start_omp_time[0];
		//double time_in_sec = (_end_time.tv_sec + _end_time.tv_usec / 1000000.0) - (_start_time.tv_sec + _start_time.tv_usec / 1000000.0);

		double GBs = (ARRAY_SIZE*4 / total_time) / 1e9;
	
		printf("sum: %d %d %d %d\n", res[0], res[1], res[2], res[3]);
		printf("Completed time: %lf - GBs: %lf.\n", total_time, GBs);

	}

	free(sum_v[0]);
	free(sum_v[1]);
	free(sum_v[2]);
	free(sum_v[3]);

}
