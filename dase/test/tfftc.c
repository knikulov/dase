#include <stdio.h>

#include "dase.h"

#define	NELEMS(a)	(sizeof(a)/sizeof(a[0]))

int
main()
{
	float x[16];
	int i;

	for (i = 0; i < NELEMS(x); i++)
		x[i] = 0.0;
	x[2] = 0.5;

	dase_fft_common(x, NELEMS(x)/2, 1);
	for (i = 0; i < NELEMS(x)/2; i++)
		printf("%f\n", x[i*2]);

	return 0;
}

