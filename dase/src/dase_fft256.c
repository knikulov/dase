#include <stdio.h>
#include <math.h>

#include "dase.h"

void
dase_fft256(float *x, float inv)
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
	X[1] = X[128];
	X[128] = tmp;
	tmp = X[2];
	X[2] = X[64];
	X[64] = tmp;
	tmp = X[3];
	X[3] = X[192];
	X[192] = tmp;
	tmp = X[4];
	X[4] = X[32];
	X[32] = tmp;
	tmp = X[5];
	X[5] = X[160];
	X[160] = tmp;
	tmp = X[6];
	X[6] = X[96];
	X[96] = tmp;
	tmp = X[7];
	X[7] = X[224];
	X[224] = tmp;
	tmp = X[8];
	X[8] = X[16];
	X[16] = tmp;
	tmp = X[9];
	X[9] = X[144];
	X[144] = tmp;
	tmp = X[10];
	X[10] = X[80];
	X[80] = tmp;
	tmp = X[11];
	X[11] = X[208];
	X[208] = tmp;
	tmp = X[12];
	X[12] = X[48];
	X[48] = tmp;
	tmp = X[13];
	X[13] = X[176];
	X[176] = tmp;
	tmp = X[14];
	X[14] = X[112];
	X[112] = tmp;
	tmp = X[15];
	X[15] = X[240];
	X[240] = tmp;
	tmp = X[17];
	X[17] = X[136];
	X[136] = tmp;
	tmp = X[18];
	X[18] = X[72];
	X[72] = tmp;
	tmp = X[19];
	X[19] = X[200];
	X[200] = tmp;
	tmp = X[20];
	X[20] = X[40];
	X[40] = tmp;
	tmp = X[21];
	X[21] = X[168];
	X[168] = tmp;
	tmp = X[22];
	X[22] = X[104];
	X[104] = tmp;
	tmp = X[23];
	X[23] = X[232];
	X[232] = tmp;
	tmp = X[25];
	X[25] = X[152];
	X[152] = tmp;
	tmp = X[26];
	X[26] = X[88];
	X[88] = tmp;
	tmp = X[27];
	X[27] = X[216];
	X[216] = tmp;
	tmp = X[28];
	X[28] = X[56];
	X[56] = tmp;
	tmp = X[29];
	X[29] = X[184];
	X[184] = tmp;
	tmp = X[30];
	X[30] = X[120];
	X[120] = tmp;
	tmp = X[31];
	X[31] = X[248];
	X[248] = tmp;
	tmp = X[33];
	X[33] = X[132];
	X[132] = tmp;
	tmp = X[34];
	X[34] = X[68];
	X[68] = tmp;
	tmp = X[35];
	X[35] = X[196];
	X[196] = tmp;
	tmp = X[37];
	X[37] = X[164];
	X[164] = tmp;
	tmp = X[38];
	X[38] = X[100];
	X[100] = tmp;
	tmp = X[39];
	X[39] = X[228];
	X[228] = tmp;
	tmp = X[41];
	X[41] = X[148];
	X[148] = tmp;
	tmp = X[42];
	X[42] = X[84];
	X[84] = tmp;
	tmp = X[43];
	X[43] = X[212];
	X[212] = tmp;
	tmp = X[44];
	X[44] = X[52];
	X[52] = tmp;
	tmp = X[45];
	X[45] = X[180];
	X[180] = tmp;
	tmp = X[46];
	X[46] = X[116];
	X[116] = tmp;
	tmp = X[47];
	X[47] = X[244];
	X[244] = tmp;
	tmp = X[49];
	X[49] = X[140];
	X[140] = tmp;
	tmp = X[50];
	X[50] = X[76];
	X[76] = tmp;
	tmp = X[51];
	X[51] = X[204];
	X[204] = tmp;
	tmp = X[53];
	X[53] = X[172];
	X[172] = tmp;
	tmp = X[54];
	X[54] = X[108];
	X[108] = tmp;
	tmp = X[55];
	X[55] = X[236];
	X[236] = tmp;
	tmp = X[57];
	X[57] = X[156];
	X[156] = tmp;
	tmp = X[58];
	X[58] = X[92];
	X[92] = tmp;
	tmp = X[59];
	X[59] = X[220];
	X[220] = tmp;
	tmp = X[61];
	X[61] = X[188];
	X[188] = tmp;
	tmp = X[62];
	X[62] = X[124];
	X[124] = tmp;
	tmp = X[63];
	X[63] = X[252];
	X[252] = tmp;
	tmp = X[65];
	X[65] = X[130];
	X[130] = tmp;
	tmp = X[67];
	X[67] = X[194];
	X[194] = tmp;
	tmp = X[69];
	X[69] = X[162];
	X[162] = tmp;
	tmp = X[70];
	X[70] = X[98];
	X[98] = tmp;
	tmp = X[71];
	X[71] = X[226];
	X[226] = tmp;
	tmp = X[73];
	X[73] = X[146];
	X[146] = tmp;
	tmp = X[74];
	X[74] = X[82];
	X[82] = tmp;
	tmp = X[75];
	X[75] = X[210];
	X[210] = tmp;
	tmp = X[77];
	X[77] = X[178];
	X[178] = tmp;
	tmp = X[78];
	X[78] = X[114];
	X[114] = tmp;
	tmp = X[79];
	X[79] = X[242];
	X[242] = tmp;
	tmp = X[81];
	X[81] = X[138];
	X[138] = tmp;
	tmp = X[83];
	X[83] = X[202];
	X[202] = tmp;
	tmp = X[85];
	X[85] = X[170];
	X[170] = tmp;
	tmp = X[86];
	X[86] = X[106];
	X[106] = tmp;
	tmp = X[87];
	X[87] = X[234];
	X[234] = tmp;
	tmp = X[89];
	X[89] = X[154];
	X[154] = tmp;
	tmp = X[91];
	X[91] = X[218];
	X[218] = tmp;
	tmp = X[93];
	X[93] = X[186];
	X[186] = tmp;
	tmp = X[94];
	X[94] = X[122];
	X[122] = tmp;
	tmp = X[95];
	X[95] = X[250];
	X[250] = tmp;
	tmp = X[97];
	X[97] = X[134];
	X[134] = tmp;
	tmp = X[99];
	X[99] = X[198];
	X[198] = tmp;
	tmp = X[101];
	X[101] = X[166];
	X[166] = tmp;
	tmp = X[103];
	X[103] = X[230];
	X[230] = tmp;
	tmp = X[105];
	X[105] = X[150];
	X[150] = tmp;
	tmp = X[107];
	X[107] = X[214];
	X[214] = tmp;
	tmp = X[109];
	X[109] = X[182];
	X[182] = tmp;
	tmp = X[110];
	X[110] = X[118];
	X[118] = tmp;
	tmp = X[111];
	X[111] = X[246];
	X[246] = tmp;
	tmp = X[113];
	X[113] = X[142];
	X[142] = tmp;
	tmp = X[115];
	X[115] = X[206];
	X[206] = tmp;
	tmp = X[117];
	X[117] = X[174];
	X[174] = tmp;
	tmp = X[119];
	X[119] = X[238];
	X[238] = tmp;
	tmp = X[121];
	X[121] = X[158];
	X[158] = tmp;
	tmp = X[123];
	X[123] = X[222];
	X[222] = tmp;
	tmp = X[125];
	X[125] = X[190];
	X[190] = tmp;
	tmp = X[127];
	X[127] = X[254];
	X[254] = tmp;
	tmp = X[131];
	X[131] = X[193];
	X[193] = tmp;
	tmp = X[133];
	X[133] = X[161];
	X[161] = tmp;
	tmp = X[135];
	X[135] = X[225];
	X[225] = tmp;
	tmp = X[137];
	X[137] = X[145];
	X[145] = tmp;
	tmp = X[139];
	X[139] = X[209];
	X[209] = tmp;
	tmp = X[141];
	X[141] = X[177];
	X[177] = tmp;
	tmp = X[143];
	X[143] = X[241];
	X[241] = tmp;
	tmp = X[147];
	X[147] = X[201];
	X[201] = tmp;
	tmp = X[149];
	X[149] = X[169];
	X[169] = tmp;
	tmp = X[151];
	X[151] = X[233];
	X[233] = tmp;
	tmp = X[155];
	X[155] = X[217];
	X[217] = tmp;
	tmp = X[157];
	X[157] = X[185];
	X[185] = tmp;
	tmp = X[159];
	X[159] = X[249];
	X[249] = tmp;
	tmp = X[163];
	X[163] = X[197];
	X[197] = tmp;
	tmp = X[167];
	X[167] = X[229];
	X[229] = tmp;
	tmp = X[171];
	X[171] = X[213];
	X[213] = tmp;
	tmp = X[173];
	X[173] = X[181];
	X[181] = tmp;
	tmp = X[175];
	X[175] = X[245];
	X[245] = tmp;
	tmp = X[179];
	X[179] = X[205];
	X[205] = tmp;
	tmp = X[183];
	X[183] = X[237];
	X[237] = tmp;
	tmp = X[187];
	X[187] = X[221];
	X[221] = tmp;
	tmp = X[191];
	X[191] = X[253];
	X[253] = tmp;
	tmp = X[199];
	X[199] = X[227];
	X[227] = tmp;
	tmp = X[203];
	X[203] = X[211];
	X[211] = tmp;
	tmp = X[207];
	X[207] = X[243];
	X[243] = tmp;
	tmp = X[215];
	X[215] = X[235];
	X[235] = tmp;
	tmp = X[223];
	X[223] = X[251];
	X[251] = tmp;
	tmp = X[239];
	X[239] = X[247];
	X[247] = tmp;

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

	/* frame number = 16 */
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

	/* frame number = 17 */
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

	/* frame number = 18 */
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

	/* frame number = 19 */
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

	/* frame number = 20 */
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

	/* frame number = 21 */
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

	/* frame number = 22 */
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

	/* frame number = 23 */
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

	/* frame number = 24 */
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

	/* frame number = 25 */
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

	/* frame number = 26 */
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

	/* frame number = 27 */
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

	/* frame number = 28 */
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

	/* frame number = 29 */
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

	/* frame number = 30 */
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

	/* frame number = 31 */
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

	/* frame number = 32 */
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

	/* frame number = 33 */
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

	/* frame number = 34 */
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

	/* frame number = 35 */
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

	/* frame number = 36 */
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

	/* frame number = 37 */
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

	/* frame number = 38 */
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

	/* frame number = 39 */
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

	/* frame number = 40 */
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

	/* frame number = 41 */
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

	/* frame number = 42 */
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

	/* frame number = 43 */
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

	/* frame number = 44 */
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

	/* frame number = 45 */
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

	/* frame number = 46 */
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

	/* frame number = 47 */
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

	/* frame number = 48 */
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

	/* frame number = 49 */
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

	/* frame number = 50 */
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

	/* frame number = 51 */
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

	/* frame number = 52 */
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

	/* frame number = 53 */
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

	/* frame number = 54 */
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

	/* frame number = 55 */
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

	/* frame number = 56 */
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

	/* frame number = 57 */
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

	/* frame number = 58 */
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

	/* frame number = 59 */
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

	/* frame number = 60 */
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

	/* frame number = 61 */
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

	/* frame number = 62 */
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

	/* frame number = 63 */
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

	/* frame number = 64 */
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

	/* frame number = 65 */
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

	/* frame number = 66 */
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

	/* frame number = 67 */
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

	/* frame number = 68 */
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

	/* frame number = 69 */
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

	/* frame number = 70 */
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

	/* frame number = 71 */
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

	/* frame number = 72 */
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

	/* frame number = 73 */
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

	/* frame number = 74 */
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

	/* frame number = 75 */
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

	/* frame number = 76 */
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

	/* frame number = 77 */
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

	/* frame number = 78 */
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

	/* frame number = 79 */
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

	/* frame number = 80 */
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

	/* frame number = 81 */
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

	/* frame number = 82 */
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

	/* frame number = 83 */
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

	/* frame number = 84 */
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

	/* frame number = 85 */
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

	/* frame number = 86 */
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

	/* frame number = 87 */
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

	/* frame number = 88 */
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

	/* frame number = 89 */
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

	/* frame number = 90 */
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

	/* frame number = 91 */
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

	/* frame number = 92 */
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

	/* frame number = 93 */
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

	/* frame number = 94 */
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

	/* frame number = 95 */
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

	/* frame number = 96 */
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

	/* frame number = 97 */
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

	/* frame number = 98 */
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

	/* frame number = 99 */
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

	/* frame number = 100 */
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

	/* frame number = 101 */
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

	/* frame number = 102 */
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

	/* frame number = 103 */
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

	/* frame number = 104 */
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

	/* frame number = 105 */
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

	/* frame number = 106 */
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

	/* frame number = 107 */
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

	/* frame number = 108 */
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

	/* frame number = 109 */
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

	/* frame number = 110 */
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

	/* frame number = 111 */
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

	/* frame number = 112 */
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

	/* frame number = 113 */
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

	/* frame number = 114 */
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

	/* frame number = 115 */
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

	/* frame number = 116 */
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

	/* frame number = 117 */
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

	/* frame number = 118 */
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

	/* frame number = 119 */
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

	/* frame number = 120 */
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

	/* frame number = 121 */
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

	/* frame number = 122 */
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

	/* frame number = 123 */
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

	/* frame number = 124 */
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

	/* frame number = 125 */
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

	/* frame number = 126 */
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

	/* frame number = 127 */
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

	/* frame number = 16 */
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

	/* frame number = 17 */
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

	/* frame number = 18 */
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

	/* frame number = 19 */
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

	/* frame number = 20 */
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

	/* frame number = 21 */
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

	/* frame number = 22 */
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

	/* frame number = 23 */
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

	/* frame number = 24 */
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

	/* frame number = 25 */
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

	/* frame number = 26 */
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

	/* frame number = 27 */
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

	/* frame number = 28 */
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

	/* frame number = 29 */
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

	/* frame number = 30 */
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

	/* frame number = 31 */
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

	/* frame number = 32 */
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

	/* frame number = 33 */
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

	/* frame number = 34 */
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

	/* frame number = 35 */
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

	/* frame number = 36 */
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

	/* frame number = 37 */
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

	/* frame number = 38 */
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

	/* frame number = 39 */
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

	/* frame number = 40 */
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

	/* frame number = 41 */
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

	/* frame number = 42 */
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

	/* frame number = 43 */
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

	/* frame number = 44 */
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

	/* frame number = 45 */
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

	/* frame number = 46 */
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

	/* frame number = 47 */
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

	/* frame number = 48 */
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

	/* frame number = 49 */
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

	/* frame number = 50 */
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

	/* frame number = 51 */
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

	/* frame number = 52 */
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

	/* frame number = 53 */
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

	/* frame number = 54 */
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

	/* frame number = 55 */
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

	/* frame number = 56 */
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

	/* frame number = 57 */
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

	/* frame number = 58 */
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

	/* frame number = 59 */
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

	/* frame number = 60 */
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

	/* frame number = 61 */
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

	/* frame number = 62 */
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

	/* frame number = 63 */
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

	/* frame number = 16 */
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

	/* frame number = 17 */
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

	/* frame number = 18 */
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

	/* frame number = 19 */
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

	/* frame number = 20 */
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

	/* frame number = 21 */
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

	/* frame number = 22 */
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

	/* frame number = 23 */
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

	/* frame number = 24 */
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

	/* frame number = 25 */
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

	/* frame number = 26 */
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

	/* frame number = 27 */
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

	/* frame number = 28 */
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

	/* frame number = 29 */
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

	/* frame number = 30 */
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

	/* frame number = 31 */
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

	/* frame size = 32 */
	w.b[DASE_RE] = 32610.2;
	w.b[DASE_IM] = inv * 3211.83;

	p = X;
	p2 = X + 32;

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

	p  += 32;
	p2 += 32;

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

	p  += 32;
	p2 += 32;

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

	p  += 32;
	p2 += 32;

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

	p  += 32;
	p2 += 32;

	/* frame size = 64 */
	w.b[DASE_RE] = 32728.5;
	w.b[DASE_IM] = inv * 1607.85;

	p = X;
	p2 = X + 64;

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

	p  += 64;
	p2 += 64;

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

	p  += 64;
	p2 += 64;

	/* frame size = 128 */
	w.b[DASE_RE] = 32758.1;
	w.b[DASE_IM] = inv * 804.167;

	p = X;
	p2 = X + 128;

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

	p  += 128;
	p2 += 128;

	fp = x + 2*256;
	X = (v2q15_union *)x + 256;

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

