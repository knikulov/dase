#include <stdio.h>

#include "dase.h"

#define	NELEMS(a)	(sizeof(a)/sizeof(a[0]))

int
main()
{
	float x[] = {0.0, 0.9, 0.9, 0.0, 0.0, 0.0};
	float h[] = {0.999969482421875, 0.2};
	float y[6];
	int i;

	dase_convol(x, NELEMS(x), h, NELEMS(h), y);

	for (i = 0; i < 6; i++)
		printf("%f\n", y[i]);

	return 0;
}

