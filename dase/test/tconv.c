#include <stdio.h>

#include "dase.h"

#define	NELEMS(a)	(sizeof(a)/sizeof(a[0]))

extern int dase_convol_simple(const float *x, int nx, const float *h, int nh, float *y);

int
main()
{
	float x[] = {0.0, 0.9, 0.9, 0.0, 0.0, 0.0};
	float h[] = {0.999969482421875, 0.2};
	float y[6];
	int i;

	printf("Simple:\n");
	dase_convol_simple(x, NELEMS(x), h, NELEMS(h), y);

	for (i = 0; i < 6; i++)
		printf("%f\n", y[i]);

	printf("\nWith accumulator:\n");
	dase_convol(x, NELEMS(x), h, NELEMS(h), y);

	for (i = 0; i < 6; i++)
		printf("%f\n", y[i]);

	return 0;
}

