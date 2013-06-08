#include <stdio.h>
#include <math.h>

#include "dase.h"

void
dase_fft_common(float *x, int n, int inv)
{
	int frame;
	int i, j;
	float *fp;
	v2q15_union *X;

	fp = x;
	X = (v2q15_union *)x;
	for (i = 0; i < n; i++) {
		X->b[1] = (*fp++) * 32768.0;
		X->b[0] = (*fp++) * 32768.0;
		X++;
	}

	X = (v2q15_union *)x;
	for (i = 0 ; i < n; i++) {
		v2q15_union tmp;
		int c, s;

		s = 0;
		for (c = 1; c < n; c <<= 1) {
			s <<= 1;
			if (i & c)
				s |= 0x1;
		}

		if (s > i) {
			tmp = X[i];
			X[i] = X[s];
			X[s] = tmp;
		}
	}

	for (frame = 1; frame < n; frame <<= 1) {
		v2q15_union w;
		v2q15_union *p, *p2;
		float arg;

		arg = M_PI/frame;
		w.b[1] = 32768.0 * cos(arg);
		w.b[0] = (inv ? 32768.0 : -32768.0) * sin(arg);

		p = X;
		p2 = X + frame;

		for (i = 0; i < n; i += (frame << 1)) {
			v2q15_union s;

			s.b[1] = 32767.0;
			s.b[0] = 0.0;

			for (j = 0; j < frame; j++) {
				a64 ac0;
				v2q15_union tmp, sinv;

				sinv.b[0] = s.b[1];
				sinv.b[1] = s.b[0];

				ac0 = 0;
				ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
				tmp.b[1] = __builtin_mips_extr_s_h(ac0, 16);

				ac0 = 0;
				ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
				tmp.b[0] = __builtin_mips_extr_s_h(ac0, 16);

				p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
				p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);

				ac0 = 0;
				ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
				s.b[1] = __builtin_mips_extr_s_h(ac0, 16);

				ac0 = 0;
				ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
				s.b[0] = __builtin_mips_extr_s_h(ac0, 16);

				p++;
				p2++;
			}
			p  += frame;
			p2 += frame;
		}
	}

	fp = x + 2*n;
	X = (v2q15_union *)x + n;
	for (i = 0; i < n; i++) {
		--X;
		*--fp = X->b[0] / 32768.0;
		*--fp = X->b[1] / 32768.0;
	}
}

