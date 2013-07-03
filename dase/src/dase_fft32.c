#include <stdio.h>
#include <math.h>

#include "dase.h"

void
dase_fft32(float *x, float inv)
{
	a64 ac0;
	v2q15_union w, s, sinv, tmp;
	v2q15_union *X, *p, *p2;
	float *fp;

	fp = x;
	X = (v2q15_union *)x;

	/* conversion to q15 format*/
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;
	X->b[DASE_RE] = (*fp++) * 32768.0;
	X->b[DASE_IM] = (*fp++) * 32768.0;
	X++;

	/* bit reversal sorting */
	X = (v2q15_union *)x;

	tmp = X[1];
	X[1] = X[16];
	X[16] = tmp;
	tmp = X[2];
	X[2] = X[8];
	X[8] = tmp;
	tmp = X[3];
	X[3] = X[24];
	X[24] = tmp;
	tmp = X[5];
	X[5] = X[20];
	X[20] = tmp;
	tmp = X[6];
	X[6] = X[12];
	X[12] = tmp;
	tmp = X[7];
	X[7] = X[28];
	X[28] = tmp;
	tmp = X[9];
	X[9] = X[18];
	X[18] = tmp;
	tmp = X[11];
	X[11] = X[26];
	X[26] = tmp;
	tmp = X[13];
	X[13] = X[22];
	X[22] = tmp;
	tmp = X[15];
	X[15] = X[30];
	X[30] = tmp;
	tmp = X[19];
	X[19] = X[25];
	X[25] = tmp;
	tmp = X[23];
	X[23] = X[29];
	X[29] = tmp;

	/* transform itself */
	/* frame size = 1 */
	w.b[DASE_RE] = -32768;
	w.b[DASE_IM] = inv * 1.05876e-10;

	p = X;
	p2 = X + 1;

	/* frame number = 0 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 1;
	p2 += 1;

	/* frame number = 1 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 1;
	p2 += 1;

	/* frame number = 2 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 1;
	p2 += 1;

	/* frame number = 3 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 1;
	p2 += 1;

	/* frame number = 4 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 1;
	p2 += 1;

	/* frame number = 5 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 1;
	p2 += 1;

	/* frame number = 6 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 1;
	p2 += 1;

	/* frame number = 7 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 1;
	p2 += 1;

	/* frame number = 8 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 1;
	p2 += 1;

	/* frame number = 9 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 1;
	p2 += 1;

	/* frame number = 10 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 1;
	p2 += 1;

	/* frame number = 11 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 1;
	p2 += 1;

	/* frame number = 12 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 1;
	p2 += 1;

	/* frame number = 13 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 1;
	p2 += 1;

	/* frame number = 14 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 1;
	p2 += 1;

	/* frame number = 15 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 1;
	p2 += 1;

	/* frame size = 2 */
	w.b[DASE_RE] = 5.29382e-11;
	w.b[DASE_IM] = inv * 32768;

	p = X;
	p2 = X + 2;

	/* frame number = 0 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 2;
	p2 += 2;

	/* frame number = 1 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 2;
	p2 += 2;

	/* frame number = 2 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 2;
	p2 += 2;

	/* frame number = 3 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 2;
	p2 += 2;

	/* frame number = 4 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 2;
	p2 += 2;

	/* frame number = 5 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 2;
	p2 += 2;

	/* frame number = 6 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 2;
	p2 += 2;

	/* frame number = 7 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 2;
	p2 += 2;

	/* frame size = 4 */
	w.b[DASE_RE] = 23170.5;
	w.b[DASE_IM] = inv * 23170.5;

	p = X;
	p2 = X + 4;

	/* frame number = 0 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 4;
	p2 += 4;

	/* frame number = 1 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 4;
	p2 += 4;

	/* frame number = 2 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 4;
	p2 += 4;

	/* frame number = 3 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 4;
	p2 += 4;

	/* frame size = 8 */
	w.b[DASE_RE] = 30273.7;
	w.b[DASE_IM] = inv * 12539.8;

	p = X;
	p2 = X + 8;

	/* frame number = 0 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 8;
	p2 += 8;

	/* frame number = 1 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 8;
	p2 += 8;

	/* frame size = 16 */
	w.b[DASE_RE] = 32138.4;
	w.b[DASE_IM] = inv * 6392.72;

	p = X;
	p2 = X + 16;

	/* frame number = 0 */
	s.b[DASE_RE] = 32767.0;
	s.b[DASE_IM] = 0.0;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	/* Fourier butterfly */
	sinv.b[DASE_RE] = s.b[DASE_IM];
	sinv.b[DASE_IM] = s.b[DASE_RE];
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, p2->a, s.a);
	tmp.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, p2->a, sinv.a);
	tmp.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);
	p2->a = __builtin_mips_subq_s_ph(p->a, tmp.a);
	p->a  = __builtin_mips_addq_s_ph(p->a,  tmp.a);
	ac0 = 0;
	ac0 = __builtin_mips_mulsaq_s_w_ph(ac0, w.a, s.a);
	s.b[DASE_RE] = __builtin_mips_extr_s_h(ac0, 16);
	ac0 = 0;
	ac0 = __builtin_mips_dpaq_s_w_ph(ac0, w.a, sinv.a);
	s.b[DASE_IM] = __builtin_mips_extr_s_h(ac0, 16);

	p++;
	p2++;

	p  += 16;
	p2 += 16;

	fp = x + 2*32;
	X = (v2q15_union *)x + 32;

	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
}

