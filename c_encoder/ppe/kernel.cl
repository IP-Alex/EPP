
kernel void convert(global float* in, global float* out) {
	// what’s my id?
	int id = get_global_id(0);
	
	// load in data to local variables
	float R = in[id];
	float G = in[id];
	float B = in[id];

	// compute
	float Y = 0 + ((float)0.299*R) + ((float)0.587*G) + ((float)0.113*B);
	float Cb = 128 - ((float)0.168736*R) - ((float)0.331264*G) + ((float)0.5*B);
	float Cr = 128 + ((float)0.5*R) - ((float)0.418688*G) - ((float)0.081312*B);

	// store data in global results array
	out[id] = Y;
	out[id] = Cb;
	out[id] = Cr;
}