//#pragma once
#ifndef opencl_helper_h
#define opencl_helper_h



#include <CL/opencl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void setup_cl(cl_device_id*, cl_context*, cl_command_queue*);

#endif // !opencl_helper_h