ffts =	dase_fft8.c dase_fft16.c dase_fft32.c dase_fft64.c\
	dase_fft128.c dase_fft256.c dase_fft512.c dase_fft1024.c\
	dase_fft2048.c

fftgen = ../tools/fftgen

$(ffts): %.c: $(fftgen)/fftgen.awk
	cd $(fftgen); make
	mv $(fftgen)/dase_fft* .

