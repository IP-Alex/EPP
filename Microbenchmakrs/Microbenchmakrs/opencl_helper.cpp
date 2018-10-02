#include "opencl_helper.h"

void setup_cl(cl_device_id *device, cl_context *context, cl_command_queue *queue) {
	//Taken from example file given. Opencl_Optimization_Example
	cl_device_type device_type = CL_DEVICE_TYPE_GPU;
	cl_platform_id platforms[10];
	cl_uint number_of_platforms_present;
	cl_int error = clGetPlatformIDs(10, platforms, &number_of_platforms_present);
	if (	error != CL_SUCCESS) {
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

char* load_kernel_file(const char* filename) {
	char* kernel_text;
	FILE* file = fopen(filename, "r");
	struct stat file_info;
	int error_code = stat(filename, &file_info);

	if (error_code != 0) {
		printf("[ERROR] - Get file failed. Error: %s", filename);
		exit(-1);
	}

	

	//Allocate memory for the kernel file
	kernel_text = (char*)malloc(file_info.st_size+1);
	memset(kernel_text, 0, file_info.st_size+1);

	//Read file, exit if unsuccessful
	size_t result = fread(kernel_text, file_info.st_size, 1, file);
	/*if(result != 1) {
		printf("[ERROR] - Get file failed. Error: %s", filename);
		exit(-1);
	}
	*/
	return kernel_text;
}

void compile_kernel(cl_device_id* device, cl_context *context, cl_kernel* kernel, cl_program* program_cl) {
	cl_int error;
	char* program_text = load_kernel_file("kernel.cl");

	*program_cl = clCreateProgramWithSource(*context, 1, (const char**)&program_text, NULL, &error);
	if (error != CL_SUCCESS) {
		printf("[ERROR] - Create context failed. Error: %d", error);
		exit(-1);
	}

	clBuildProgram(*program_cl, 1, device, NULL, NULL, NULL);
	if (error != CL_SUCCESS) {
		printf("[ERROR] - Build kernel failed. Error: %d", error);
		exit(-1);
	}

	*kernel = clCreateKernel(*program_cl, "eightbit", &error);
	if (error != CL_SUCCESS) {
		printf("[ERROR] - Create kernel failed. Error: %d", error);
		exit(-1);
	}

	free(program_text);
}


void cleanup(cl_kernel* kernel, cl_program* program) {
	//clReleaseMemObject(out_buffer);
	//clReleaseMemObject(in_buffer);
	clReleaseKernel(*kernel);
	clReleaseProgram(*program);
}