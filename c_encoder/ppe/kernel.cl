kernel void convert(global float* in_r, global float* in_g, global float* in_b, 
					global float* out_r, global float* out_g, global float* out_b, 
					global int* workload, global int* image_size) {
	// what’s my id?
	int id = get_global_id(0);
	int threads = get_global_size(0);
	for (int i = id; i < *image_size; i += threads) {

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
	}
}


kernel void lowpassRow(global float* in_channel, global float* out_channel, global int* width, global int* height) {
	int id = get_global_id(0);
	int threads = get_global_size(0);

	//Blur weights 
	float a = (float)0.25f;
	float b = (float)0.5f;
	float c = (float)0.25f;
	
	/*for (int index = id; index < ((*width) * (*height)); index += threads) {
		int column = index % *width;
		
		if (index >= (*width) && index < (*width) * (*height - 1) && column > 0 && column < *width - 1) {
			out_channel[index] = (float)((a * (float)(in_channel[index - (*width)])) + (b * (float)(in_channel[index])) + (c * (float)(in_channel[index + (*width)])));
		}
	}*/
}

kernel void lowpassColumn(global float* in_channel, global float* out_channel, global int* width, global int* height) {

	int id = get_global_id(0);
	int threads = get_global_size(0);

	//Blur weights 
	float a = (float)0.25f;
	float b = (float)0.5f;
	float c = (float)0.25f;
	int column = 0;

	/*for (int index = id; index < ((*width) * (*height)); index += threads) {
		column = index % *width;

		if (index >= (*width) && index < (*width) * (*height - 1) && column > 0 && column < *width - 1) {
			out_channel[index] = (float)((a * (float)(in_channel[index - 1])) + (b * (float)(in_channel[index])) + (c * (float)(in_channel[index + 1])));
		}
	}*/
}

