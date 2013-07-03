#include <stdio.h>
#include <math.h>

#include "dase.h"

void
dase_fft128(float *x, float inv)
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

	/* bit reversal sorting */
	X = (v2q15_union *)x;

	tmp = X[1];
	X[1] = X[64];
	X[64] = tmp;
	tmp = X[2];
	X[2] = X[32];
	X[32] = tmp;
	tmp = X[3];
	X[3] = X[96];
	X[96] = tmp;
	tmp = X[4];
	X[4] = X[16];
	X[16] = tmp;
	tmp = X[5];
	X[5] = X[80];
	X[80] = tmp;
	tmp = X[6];
	X[6] = X[48];
	X[48] = tmp;
	tmp = X[7];
	X[7] = X[112];
	X[112] = tmp;
	tmp = X[9];
	X[9] = X[72];
	X[72] = tmp;
	tmp = X[10];
	X[10] = X[40];
	X[40] = tmp;
	tmp = X[11];
	X[11] = X[104];
	X[104] = tmp;
	tmp = X[12];
	X[12] = X[24];
	X[24] = tmp;
	tmp = X[13];
	X[13] = X[88];
	X[88] = tmp;
	tmp = X[14];
	X[14] = X[56];
	X[56] = tmp;
	tmp = X[15];
	X[15] = X[120];
	X[120] = tmp;
	tmp = X[17];
	X[17] = X[68];
	X[68] = tmp;
	tmp = X[18];
	X[18] = X[36];
	X[36] = tmp;
	tmp = X[19];
	X[19] = X[100];
	X[100] = tmp;
	tmp = X[21];
	X[21] = X[84];
	X[84] = tmp;
	tmp = X[22];
	X[22] = X[52];
	X[52] = tmp;
	tmp = X[23];
	X[23] = X[116];
	X[116] = tmp;
	tmp = X[25];
	X[25] = X[76];
	X[76] = tmp;
	tmp = X[26];
	X[26] = X[44];
	X[44] = tmp;
	tmp = X[27];
	X[27] = X[108];
	X[108] = tmp;
	tmp = X[29];
	X[29] = X[92];
	X[92] = tmp;
	tmp = X[30];
	X[30] = X[60];
	X[60] = tmp;
	tmp = X[31];
	X[31] = X[124];
	X[124] = tmp;
	tmp = X[33];
	X[33] = X[66];
	X[66] = tmp;
	tmp = X[35];
	X[35] = X[98];
	X[98] = tmp;
	tmp = X[37];
	X[37] = X[82];
	X[82] = tmp;
	tmp = X[38];
	X[38] = X[50];
	X[50] = tmp;
	tmp = X[39];
	X[39] = X[114];
	X[114] = tmp;
	tmp = X[41];
	X[41] = X[74];
	X[74] = tmp;
	tmp = X[43];
	X[43] = X[106];
	X[106] = tmp;
	tmp = X[45];
	X[45] = X[90];
	X[90] = tmp;
	tmp = X[46];
	X[46] = X[58];
	X[58] = tmp;
	tmp = X[47];
	X[47] = X[122];
	X[122] = tmp;
	tmp = X[49];
	X[49] = X[70];
	X[70] = tmp;
	tmp = X[51];
	X[51] = X[102];
	X[102] = tmp;
	tmp = X[53];
	X[53] = X[86];
	X[86] = tmp;
	tmp = X[55];
	X[55] = X[118];
	X[118] = tmp;
	tmp = X[57];
	X[57] = X[78];
	X[78] = tmp;
	tmp = X[59];
	X[59] = X[110];
	X[110] = tmp;
	tmp = X[61];
	X[61] = X[94];
	X[94] = tmp;
	tmp = X[63];
	X[63] = X[126];
	X[126] = tmp;
	tmp = X[67];
	X[67] = X[97];
	X[97] = tmp;
	tmp = X[69];
	X[69] = X[81];
	X[81] = tmp;
	tmp = X[71];
	X[71] = X[113];
	X[113] = tmp;
	tmp = X[75];
	X[75] = X[105];
	X[105] = tmp;
	tmp = X[77];
	X[77] = X[89];
	X[89] = tmp;
	tmp = X[79];
	X[79] = X[121];
	X[121] = tmp;
	tmp = X[83];
	X[83] = X[101];
	X[101] = tmp;
	tmp = X[87];
	X[87] = X[117];
	X[117] = tmp;
	tmp = X[91];
	X[91] = X[109];
	X[109] = tmp;
	tmp = X[95];
	X[95] = X[125];
	X[125] = tmp;
	tmp = X[103];
	X[103] = X[115];
	X[115] = tmp;
	tmp = X[111];
	X[111] = X[123];
	X[123] = tmp;

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

	fp = x + 2*128;
	X = (v2q15_union *)x + 128;

	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
	--X;
	*--fp = X->b[DASE_IM] / 32768.0;
	*--fp = X->b[DASE_RE] / 32768.0;
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

