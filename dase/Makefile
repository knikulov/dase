CC = mips-linux-gnu-gcc

CFLAGS = -I include -Wall -mips32r2 -mdsp -static

headers = dase.h

sources = dase_convol.c dase_fft8.c dase_fft16.c dase_fft32.c\
	  dase_fft64.c dase_fft128.c dase_fft256.c dase_fft512.c\
	  dase_fft1024.c dase_fft2048.c

objects = dase_convol.o dase_fft8.o dase_fft16.o dase_fft32.o\
	  dase_fft64.o dase_fft128.o dase_fft256.o dase_fft512.o\
	  dase_fft1024.o dase_fft2048.o

libfile = libdase.so

vpath %.h include
vpath %.c src test

.PHONY: build
build:
	cd src; make
	make $(libfile)

$(libfile): %.so: $(objects)
	$(CC) -shared -Wl,-soname,$@ -o $@ $^

%.o: %.c $(headers)
	$(CC) $(CFLAGS) -fPIC -c -o $@ $<

%: %.c $(objects)
	$(CC) $(CFLAGS) -o $@ $^ -lm

.PHONY: install clean
install:
	install $(libfile) /usr/lib/
	install ./include/dase.h /usr/include/

clean:
	rm -f src/dase_fft*
	rm -f $(objects)

