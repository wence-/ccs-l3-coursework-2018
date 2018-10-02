CFLAGS = -O3 -march=native  -D_GNU_SOURCE
LDFLAGS = -lm
CC = gcc

OBJ = optimised-sparsemm.o basic-sparsemm.o utils.o
HEADER = utils.h

.PHONY: clean help check

all: sparsemm

help:
	@echo "Available targets are"
	@echo "  clean: Remove all build artifacts"
	@echo "  check: Perform a simple test of your optimised routines"
	@echo "  sparsemm: Build the sparse matrix-matrix multiplication binary"

clean:
	-rm -f sparsemm $(OBJ)

check: sparsemm
	./sparsemm CHECK

sparsemm: sparsemm.c $(OBJ)
	$(CC) $(CFLAGS) -o $@ $< $(OBJ) $(LDFLAGS)

%.o: %.c $(HEADER)
	$(CC) $(CFLAGS) -c -o $@ $<
