kernel void convert(global float* in_r, global float* in_g, global float* in_b, 
					global float* out_r, global float* out_g, global float* out_b, 
					global int* workload, global int* image_size) {
	// what’s my id?
	int id = get_global_id(0);
	int threads = get_global_size(0);
	int j = 0;
	for (int i = id; j < *workload && i < *image_size; i += threads) {

		// load in data to local variables
		float R = in_r[i];
		float G = in_g[i];
		float B = in_b[i];

		// compute
		float Y = 0 + ((float)0.299*R) + ((float)0.587*G) + ((float)0.113*B);
		float Cb = 128 - ((float)0.168736*R) - ((float)0.331264*G) + ((float)0.5*B);
		float Cr = 128 + ((float)0.5*R) - ((float)0.418688*G) - ((float)0.081312*B);

		// store data in global results array
		out_r[i] = Y;
		out_g[i] = Cb;
		out_b[i] = Cr;
		j++;
	}
}


kernel void lowpass(global float* in_channel, global float* out_channel,
					global int* workload, global int* width, global int* height) {
	int id = get_global_id(0);
	int threads = get_global_size(0);

	//Blur weights 
	float a = 0.25;
	float b = 0.5;
	float c = 0.25;

	int work_performed = 0;

	for (int index = id; work_performed < *workload && index < (*width) * (*height) - 1; index += threads) {
		int row = index / width;
		int column = index % width;

		if (row > 0 && row < height && column > 0 && column < width)//Would be nice to remove the control statement here if possible
			out_channel[index] = a * in_channel[(row - 1)*width + column] + b * in->data[row*width + column] + c * in->data[(row + 1)*width + column];

		work_performed++;
	}
}