#include <stdio.h>

#include "dase.h"

void
dase_convol_simple(const float *x, int nx, const float *h, int nh, float *y)
{
	v2q15_union r0, r1;
	a64 ac0;
	int ny;
	int i, j, k;

	ny = nx + nh - 1;
	for (i = 0; i < ny; i++) {

		ac0 = 0;
		for (j = 0; j < nh; j++) {

			k = i - j;
			if ((0 < k) && (k < nx)) {
				r0.b[0] = h[j] * 32768.0;
				r1.b[0] = x[k] * 32768.0;
				ac0 = __builtin_mips_maq_s_w_phl(ac0,
								  r0.a,
								  r1.a);
			}
		}

		y[i] = ac0 / 2147483648.0;
	}
}

void
dase_convol(const float *x, int nx, const float *h, int nh, float *y)
{
	v2q15_union r0, r1;
	a64 ac0;
	int ny;
	int i, j, k;

	ny = nx + nh - 1;
	for (i = 0; i < ny; i++) {

		ac0 = 0;
		for (j = (nh & 0x1); nh-j >= 2; j += 2) {
			k = i - j;
			if ((0 < k) && (k+1 < nx)) {
				r0.b[0] = h[j]   * 32768.0;
				r0.b[1] = h[j+1] * 32768.0;
				r1.b[0] = x[k]   * 32768.0;
				r1.b[1] = x[k-1] * 32768.0;
				ac0 = __builtin_mips_dpaq_s_w_ph(ac0,
								  r0.a,
								  r1.a);
			}
		}

		if (nh & 0x1) {
			if (i < nx) {
				r0.b[0] = h[0] * 32768.0;
				r1.b[0] = x[i] * 32768.0;
				ac0 = __builtin_mips_maq_s_w_phl(ac0,
								 r0.a,
								 r1.a);
			}
		}

		y[i] = ac0 / 2147483648.0;
	}
}

