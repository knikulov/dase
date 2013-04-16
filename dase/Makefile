CC = mips-linux-gnu-gcc

CFLAGS = -I include -Wall -mips32r2 -mdsp -static

headers = dase.h

sources = dase_convol.c

objects = dase_convol.o

vpath %.c src test tools
vpath %.h include

%: %.c $(objects)
	$(CC) $(CFLAGS) -o $@ $^

$(objects): %.o: %.c $(headers)
	$(CC) $(CFLAGS) -fPIC -c -o $@ $<

.PHONY: clean
clean:
	rm -f $(objects)

