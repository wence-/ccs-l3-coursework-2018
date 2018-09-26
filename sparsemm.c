#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "utils.h"

void optimised_sparsemm(const COO, const COO, COO*);
void optimised_sparsemm_sum(const COO, const COO, const COO,
                            const COO, const COO, const COO,
                            COO *);

int main(int argc, char **argv)
{
    COO O;
    FILE *f;
    if (!(argc == 4 || argc == 8)) {
        fprintf(stderr, "Invalid arguments.\n");
        fprintf(stderr, "Usage: %s O A B\n", argv[0]);
        fprintf(stderr, "  Computes O = A B\n");
        fprintf(stderr, "  Where A and B are filenames of matrices to read.\n");
        fprintf(stderr, "  O is the filename of the matrix to be written.\n\n");
        fprintf(stderr, "Alternate usage: %s O A B C D E F\n", argv[0]);
        fprintf(stderr, "  Computes O = (A + B + C) (D + E + F)\n");
        fprintf(stderr, "  Where A-F are the files names of matrices to read.\n");
        fprintf(stderr, "  O is the filename of the matrix to be written.\n\n");
    }


    if (argc == 4) {
        COO A, B;
        read_sparse(argv[2], &A);
        read_sparse(argv[3], &B);

        optimised_sparsemm(A, B, &O);

        free_sparse(&A);
        free_sparse(&B);
    } else {
        COO A, B, C, D, E, F;
        read_sparse(argv[2], &A);
        read_sparse(argv[3], &B);
        read_sparse(argv[4], &C);
        read_sparse(argv[5], &D);
        read_sparse(argv[6], &E);
        read_sparse(argv[7], &F);

        optimised_sparsemm_sum(A, B, C, D, E, F, &O);

        free_sparse(&A);
        free_sparse(&B);
        free_sparse(&C);
        free_sparse(&D);
        free_sparse(&E);
        free_sparse(&F);
    }

    f = fopen(argv[1], "w");
    write_sparse(f, O);
    free_sparse(&O);
    fclose(f);
    return 0;
}
