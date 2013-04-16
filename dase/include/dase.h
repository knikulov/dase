#ifndef DASE_H_
#define	DASE_H_

#include <stdint.h>

typedef	int16_t q15;
typedef int32_t q31;
typedef int64_t a64;

typedef int16_t v2q15 __attribute__ ((vector_size(4)));

typedef union {
	v2q15 a;
	q15 b[2];
} v2q15_union;

typedef union {
	a64 a;
	q31 b[2];
} a64_union;

void dase_convol(const float *x, int nx, const float *h, int nh, float *y);

#endif

