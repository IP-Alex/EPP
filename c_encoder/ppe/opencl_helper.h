//#pragma once
#ifndef opencl_helper_h
#define opencl_helper_h



#include <CL/opencl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <assert.h>

#define blooper_t 128
#define solar_t 2048
#define solar4k_t 4086

#define SIZE blooper_t
#define SIZE_BYTES sizeof(float)*SIZE*SIZE

#define TIME_RUNTIME 1
#define TIME_GPU 2
#define NO_TIMEING 3
#define RUN_MODE TIME_RUNTIME
#define GPU_CORES SIZE*SIZE

typedef enum colorchannel {R=0,G=1,B=2};

void setup_cl(cl_device_id*, cl_context*, cl_command_queue*);
void compile_kernel(cl_device_id*, cl_context*, cl_kernel*, cl_program*, char* function_name);
void cleanup(cl_kernel*, cl_program*);
void kernel_init(cl_mem*, cl_mem*, cl_context*);
//void kernel_calculate(cl_mem*, cl_mem*, cl_context*, cl_command_queue*,
//					cl_kernel*, size_t[], float*, float*, float*, float*, float*, float*);

#endif // !opencl_helper_h