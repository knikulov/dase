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

enum {
	DASE_FFT_DIR  = -1,
	DASE_FFT_INV  =  1
};

typedef	int16_t q15;
typedef int32_t q31;
typedef int64_t a64;

typedef int16_t v2q15 __attribute__ ((vector_size(4)));

typedef union {
	v2q15 a;
	q15 b[2];
} v2q15_union;

void dase_convol(const float *x, int nx, const float *h, int nh, float *y);

void dase_fft_common(float *x, int n, float inv);
void dase_fft_radix4(float *x, int n);

void dase_fft8(float *x, float inv);
void dase_fft16(float *x, float inv);
void dase_fft32(float *x, float inv);
void dase_fft64(float *x, float inv);
void dase_fft128(float *x, float inv);
void dase_fft256(float *x, float inv);
void dase_fft512(float *x, float inv);
void dase_fft1024(float *x, float inv);
void dase_fft2048(float *x, float inv);

#endif

