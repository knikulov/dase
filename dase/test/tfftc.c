#include <stdio.h>

#include "dase.h"

#define	N	2048

int
main()
{
	float x[N];
	int i;

	for (i = 0; i < N; i++)
		x[i] = 0.0;
	x[2] = 0.5;
	x[9] = 0.2;
	x[14] = 0.15;

	dase_fft1024(x, DASE_FFT_INV);
	for (i = 0; i < N/2; i++)
		printf("%f\n", x[i*2]);

	return 0;
}

