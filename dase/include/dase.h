#ifndef DASE_H_
#define	DASE_H_

#include <stdint.h>
#include <endian.h>

#if __BYTE_ORDER == __LITTLE_ENDIAN
 #define	DASE_RE	1
 #define	DASE_IM	0
#elif __BYTE_ORDER == __BIG_ENDIAN
 #define	DASE_RE	0
 #define	DASE_IM	1
#endif

typedef	int16_t q15;
typedef int32_t q31;
typedef int64_t a64;

typedef int16_t v2q15 __attribute__ ((vector_size(4)));

typedef union {
	v2q15 a;
	q15 b[2];
} v2q15_union;

void dase_convol(const float *x, int nx, const float *h, int nh, float *y);

void dase_fft_common(float *x, int n, int inv);

#endif

