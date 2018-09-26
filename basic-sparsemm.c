#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "utils.h"

/* Compute C = A*B in dense, column major, format. C must be zero on entry. */
static void dgemm(int m, int n, int k, const double *a, const double *b, double *c)
{
    int i, j, p;
    int lda = m;
    int ldb = k;
    int ldc = m;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            /* c[i, j] <-  c[i, j] + a[i, p] * b[p, j]
             * i.e. the dot product of the ith row of a with the jth
             * column of b. */
            for (p = 0; p < k; p++) {
                c[j*ldc + i] = c[j*ldc + i] + a[p*lda + i] * b[j*ldb + p];
            }
        }
    }
}

/* Computes C = A*B by converting A and B to dense column major
 * format, performing the matrix-matrix multiplication, and then
 * converting the result back to sparse. */
void basic_sparsemm(const COO A, const COO B, COO *C)
{
    double *a = NULL;
    double *b = NULL;
    double *c = NULL;
    int m, n, k;
    convert_sparse_to_dense(A, &a);
    convert_sparse_to_dense(B, &b);

    *C = NULL;
    m = A->m;
    k = A->n;
    n = B->n;
    if (k != B->m) {
        fprintf(stderr, "Invalid matrix sizes, got %d x %d and %d x %d\n",
                A->m, A->n, B->m, B->n);
        free(a);
        free(b);
        exit(1);
    }
    c = calloc(m*n, sizeof(double));

    dgemm(m, n, k, a, b, c);
    free(a);
    free(b);
    convert_dense_to_sparse(c, m, n, C);
    free(c);
}
