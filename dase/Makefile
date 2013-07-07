CC = mips-linux-gnu-gcc

CFLAGS = -I include -Wall -mips32r2 -mdsp -static

headers = dase.h

sources = dase_convol.c dase_fft8.c dase_fft16.c \
	  dase_fft32.c dase_fft64.c dase_fft128.c dase_fft256.c \
	  dase_fft512.c dase_fft1024.c

objects = dase_convol.o dase_fft8.o dase_fft16.o \
	  dase_fft32.o dase_fft64.o dase_fft128.o dase_fft256.o \
	  dase_fft512.o dase_fft1024.o

libfile = libdase.so

vpath %.c src test
vpath %.h include

$(libfile): %.so: $(objects)
	$(CC) -shared -Wl,-soname,$@ -o $@ $^

%: %.c $(objects)
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(objects): %.o: %.c $(headers)
	$(CC) $(CFLAGS) -fPIC -c -o $@ $<

tfftc: tfftc.o dase_fft1024.o
	$(CC) $(CFLAGS) -o $@ $^

.PHONY: clean
clean:
	rm -f $(objects)

