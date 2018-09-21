
kernel void convert(global float* in_r, global float* in_g, global float* in_b, global float* out_r, 
					global float* out_g, global float* out_b, global int* workload, global int* image_size) {
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