CC = mips-linux-gnu-gcc

CFLAGS = -I include -Wall -mips32r2 -mdsp -static

headers = dase.h

sources = dase_convol.c dase_fft_common.c

objects = dase_convol.o dase_fft_common.o

vpath %.c src test tools
vpath %.h include

%: %.c $(objects)
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(objects): %.o: %.c $(headers)
	$(CC) $(CFLAGS) -fPIC -c -o $@ $<

.PHONY: clean
clean:
	rm -f $(objects)

