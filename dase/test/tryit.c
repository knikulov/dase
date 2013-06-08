#include <stdio.h>

#include "dase.h"

int
main()
{
	a64 ac0;
	v2q15_union x, y, z;

	x.b[0] = 0x0;
	x.b[1] = 0x4000;
	y.b[0] = 0x8000;
	y.b[1] = 0x4afc;

	printf("x.b[0] = %f\n", x.b[0] / 32768.0);
	printf("x.b[1] = %f\n", x.b[1] / 32768.0);
	printf("y.b[0] = %f\n", y.b[0] / 32768.0);
	printf("y.b[1] = %f\n\n", y.b[1] / 32768.0);

	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, x.a, y.a);
	z.b[0] = __builtin_mips_extr_s_h(ac0, 16);

	printf("%f * %f - %f * %f = %f\n",
		x.b[1] / 32768.0, y.b[1] / 32768.0,
		x.b[0] / 32768.0, y.b[0] / 32768.0,
		z.b[0] / 32768.0);

	return 0;
}

