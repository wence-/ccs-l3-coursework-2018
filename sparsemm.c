#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "utils.h"

void basic_sparsemm(const COO, const COO, COO*);
void basic_sparsemm_sum(const COO, const COO, const COO,
                            const COO, const COO, const COO,
                            COO *);
void optimised_sparsemm(const COO, const COO, COO*);
void optimised_sparsemm_sum(const COO, const COO, const COO,
                            const COO, const COO, const COO,
                            COO *);

static int check_sparsemm()
{
    COO A, B, Cbasic, Copt;
    double *basic, *opt;
    int i, j, m, n, k;
    int pass = 0;

    m = 20;
    k = 50;
    n = 30;
    random_matrix(m, k, 0.1, &A);
    random_matrix(k, n, 0.2, &B);

    basic_sparsemm(A, B, &Cbasic);
    optimised_sparsemm(A, B, &Copt);

    convert_sparse_to_dense(Cbasic, &basic);
    convert_sparse_to_dense(Copt, &opt);

    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            double diff = fabs(opt[j*m + i] - basic[j*m + i]);
            if (diff != diff || diff > 1e-3) {
                fprintf(stderr, "MM Failed check at entry (%d, %d), basic value %g, opt value %g\n", i, j, basic[j*m + i], opt[j*m + i]);
                pass = 1;
            }
        }
    }

    free(basic);
    free(opt);
    free_sparse(&A);
    free_sparse(&B);
    free_sparse(&Cbasic);
    free_sparse(&Copt);

    return pass;
}

static int check_sparsemm_sum()
{
    COO A, B, C, D, E, F, Obasic, Oopt;
    double *basic, *opt;
    int i, j, m, n, k;
    int pass = 0;

    m = 20;
    k = 50;
    n = 30;
    random_matrix(m, k, 0.1, &A);
    random_matrix(m, k, 0.4, &B);
    random_matrix(m, k, 0.01, &C);
    random_matrix(k, n, 0.2, &D);
    random_matrix(k, n, 0.3, &E);
    random_matrix(k, n, 0.15, &F);
    basic_sparsemm_sum(A, B, C, D, E, F, &Obasic);
    optimised_sparsemm_sum(A, B, C, D, E, F, &Oopt);

    convert_sparse_to_dense(Obasic, &basic);
    convert_sparse_to_dense(Oopt, &opt);

    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            double diff = fabs(opt[j*m + i] - basic[j*m + i]);
            if (diff != diff || diff > 1e-3) {
                fprintf(stderr, "MM-SUM Failed check at entry (%d, %d), basic value %g, opt value %g\n", i, j, basic[j*m + i], opt[j*m + i]);
                pass = 1;
            }
        }
    }

    free(basic);
    free(opt);
    free_sparse(&A);
    free_sparse(&B);
    free_sparse(&C);
    free_sparse(&D);
    free_sparse(&E);
    free_sparse(&F);
    free_sparse(&Obasic);
    free_sparse(&Oopt);

    return pass;
}
int main(int argc, char **argv)
{
    COO O;
    FILE *f;
    if (!(argc == 2 || argc == 4 || argc == 8)) {
        fprintf(stderr, "Invalid arguments.\n");
        fprintf(stderr, "Usage: %s CHECK\n", argv[0]);
        fprintf(stderr, "  Check the implemented routines using randomly generated matrices.\n");
        fprintf(stderr, "Alternate usage: %s O A B\n", argv[0]);
        fprintf(stderr, "  Computes O = A B\n");
        fprintf(stderr, "  Where A and B are filenames of matrices to read.\n");
        fprintf(stderr, "  O is the filename of the matrix to be written.\n\n");
        fprintf(stderr, "Alternate usage: %s O A B C D E F\n", argv[0]);
        fprintf(stderr, "  Computes O = (A + B + C) (D + E + F)\n");
        fprintf(stderr, "  Where A-F are the files names of matrices to read.\n");
        fprintf(stderr, "  O is the filename of the matrix to be written.\n\n");
        return 1;
    }

    if (argc == 2) {
        int pass = 0;
        if (strcmp(argv[1], "CHECK")) {
            fprintf(stderr, "Invalid mode, expecting CHECK, got %s\n", argv[1]);
            return 1;
        }
        pass |= check_sparsemm();
        pass |= check_sparsemm_sum();
        return pass;
    } else if (argc == 4) {
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
