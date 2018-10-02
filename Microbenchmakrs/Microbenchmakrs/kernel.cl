
kernel void eightbit(global unsigned char *num, global unsigned char *num2) {
	int id = get_global_id(0);
	for (int i = 0; i < 1000; i++) {
		num[id] += *num2;
	}
}