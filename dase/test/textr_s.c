#include <stdio.h>

#include "dase.h"

int
main()
{
	a64 ac;
	v2q15_union r0, r1;

	r0.b[0] =  -0.6 * 32768.0;
	r0.b[1] =  0.0 * 32768.0;
	r1.b[0] =  0.6 * 32768.0;
	r1.b[1] =  0.0 * 32768.0;

	ac = 0;
	ac = __builtin_mips_mulsaq_s_w_ph(ac, r0.a, r1.a);
	r0.b[0] = __builtin_mips_extr_s_h(ac, 16);

//	printf("%.15f\n", r0.b[0] / 32768.0);
	printf("ac (hex) = %llx\n", ac);
	printf("ac (dec) = %.15f\n", ac / 2147483648.0);

	return 0;
}

