#include "dase.h"

void
dase_convol(const float *x, int nx, const float *h, int nh, float *y)
{
	v2q15_union r0, r1;
	a64_union ac0;
	int ny;
	int i, j, k;

	ny = nx + nh - 1;
	for (i = 0; i < ny; i++) {

		ac0.a = 0;
		for (j = 0; j < nh; j++) {

			k = i - j;
			if ((0 < k) && (k < nx)) {
				r0.b[0] = h[j] * 32768.0;
				r1.b[0] = x[k] * 32768.0;
				ac0.a = __builtin_mips_maq_s_w_phl(ac0.a,
								    r0.a,
								    r1.a);
			}
		}
		y[i] = ac0.a / 2147483648.0;
	}
}

