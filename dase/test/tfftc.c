#include <stdio.h>
#include <math.h>

#include "dase.h"

#define	N	32

int
main()
{
	float x[N];
	int i;

	for (i = 0; i < N; i += 2) {
		x[i]   = 0.3 * cosf(M_PI/(i+1));
		x[i+1] = 0.0;
	}

	dase_fft_radix4(x, N/2);
	for (i = 0; i < N; i++)
		printf("%f\n", x[i]);

	return 0;
}

