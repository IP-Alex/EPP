#include "opencl_helper.h"

void setup_cl(cl_device_id *device, cl_context *context, cl_command_queue *queue) {
	//Taken from example file given. Opencl_Optimization_Example
	cl_device_type device_type = CL_DEVICE_TYPE_GPU;
	cl_platform_id platforms[10];
	cl_uint number_of_platforms_present;
	cl_int error = clGetPlatformIDs(10, platforms, &number_of_platforms_present);
	if (error != CL_SUCCESS) {
		printf("[ERROR] - Get platform failed. Error: %d", error);
		exit(-1);
	}

	//Get device for our platform
	cl_device_id devices_on_platform[10];
	cl_uint number_of_devices_found;
	for (cl_uint i = 0; i < number_of_platforms_present; i++) {
		char name[1000];
		error = clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 1000, name, NULL);
		if (error != CL_SUCCESS) {
			printf("[ERROR] - Get platform failed. Error: %d", error);
			exit(-1);
		}

		error = clGetDeviceIDs(platforms[i], device_type, 10, devices_on_platform, &number_of_devices_found);
		if (error != CL_DEVICE_NOT_FOUND && error != CL_SUCCESS) {
			printf("[ERROR] - Get platform failed. Error: %d", error);
			exit(-1);
		}

		if (number_of_devices_found > 0)
			break;
	}

	if (number_of_devices_found == 0) {
		printf("Oh no! No devices found.");
		exit(-1);
	}

	*device = devices_on_platform[0];

	*context = clCreateContext(NULL, 1, device, NULL, NULL, &error);
	if (error != CL_SUCCESS) {
		printf("[ERROR] - Get platform failed. Error: %d", error);
		exit(-1);
	}


	*queue = clCreateCommandQueue(*context, *device, CL_QUEUE_PROFILING_ENABLE, &error);
	if (error != CL_SUCCESS) {
		printf("[ERROR] - Get platform failed. Error: %d", error);
		exit(-1);
	}
}