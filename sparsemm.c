#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "utils.h"

typedef void (*sparse_mm_t)(const COO, const COO, COO*);

void optimised_sparsemm(const COO, const COO, COO*);

int main(int argc, char **argv)
{
    COO A;
    COO B;
    COO C;
    FILE *f;
    if (argc < 4) {
        fprintf(stderr, "Invalid arguments.\n");
        fprintf(stderr, "Usage: %s C A B [...]\n", argv[0]);
        fprintf(stderr, "Where A and B are filenames of matrices to be read.\n");
        fprintf(stderr, "C is the filename of the matrix to be written.\n");
    }

    read_sparse(argv[2], &A);
    read_sparse(argv[3], &B);

    optimised_sparsemm(A, B, &C);

    free_sparse(&A);
    free_sparse(&B);

    f = fopen(argv[1], "w");
    write_sparse(f, C);
    free_sparse(&C);
    fclose(f);
    return 0;
}
