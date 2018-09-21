//#pragma once
#ifndef opencl_helper_h
#define opencl_helper_h



#include <CL/opencl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#define SIZE 2048
#define SIZE_BYTES sizeof(float)*SIZE*SIZE

void setup_cl(cl_device_id*, cl_context*, cl_command_queue*);
void compile_kernel(cl_device_id*, cl_context*, cl_kernel*, cl_program*);
void cleanup(cl_kernel*, cl_program*);
void kernel_init(cl_mem*, cl_mem*, cl_context*);
void kernel_calculate(cl_mem*, cl_mem*, cl_context*, cl_command_queue*,
	cl_kernel*, size_t[], float*, float*);

#endif // !opencl_helper_h