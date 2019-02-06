#include "utils.h"
#include <stdlib.h>

struct _p_CSR {
  int m, n, NZ;
  int *ai;
  int *aj;
  double *data;
};

typedef struct _p_CSR *CSR;

static void alloc_csr(int m, int n, int NZ, CSR *sparse)
{
  CSR sp = calloc(1, sizeof(*sp));
  sp->m = m;
  sp->n = n;
  sp->NZ = NZ;
  sp->ai = calloc(m + 1, sizeof(*sp->ai));
  sp->aj = calloc(NZ, sizeof(*sp->aj));
  sp->data = calloc(NZ, sizeof(*sp->data));
  *sparse = sp;
}

static void free_csr(CSR *sparse)
{
  CSR sp = *sparse;
  if (!sp) {
    return;
  }
  free(sp->ai);
  free(sp->aj);
  free(sp->data);
  free(sp);
  *sparse = NULL;
}

struct coord_data {
  int i, j;
  double d;
};

static int cmp_csr(const void *a_, const void *b_)
{
  const struct coord_data *a = (const struct coord_data *)a_;
  const struct coord_data *b = (const struct coord_data *)b_;

  if (!(a->i - b->i))
    return a->j - b->j;
  else
    return a->i - b->i;
}

static void convert_coo_to_csr(const COO in, CSR *out)
{
  struct coord_data *tmp = malloc(in->NZ * sizeof(*tmp));
  int i, last_row;

  for (i = 0; i < in->NZ; i++) {
    tmp[i].i = in->coords[i].i;
    tmp[i].j = in->coords[i].j;
    tmp[i].d = in->data[i];
  }

  qsort(tmp, in->NZ, sizeof(*tmp), cmp_csr);

  alloc_csr(in->m, in->n, in->NZ, out);

  for (i = 0; i < in->m + 1; i++) {
    (*out)->ai[i] = 0;
  }
  for (i = 0; i < in->NZ; i++) {
    (*out)->data[i] = tmp[i].d;
    (*out)->aj[i] = tmp[i].j;
    (*out)->ai[tmp[i].i+1]++;
  }
  for (i = 0; i < in->m; i++) {
    (*out)->ai[i+1] += (*out)->ai[i];
  }
  (*out)->ai[in->m] = in->NZ;
}

static void convert_csr_to_coo(const CSR in, COO *out)
{
  int i, j;

  alloc_sparse(in->m, in->n, in->NZ, out);
  for (i = 0; i < in->m; i++) {
    for (j = in->ai[i]; j < in->ai[i+1]; j++) {
      (*out)->coords[j].i = i;
      (*out)->coords[j].j = in->aj[j];
      (*out)->data[j] = in->data[j];
    }
  }
}

static void print_csr(const CSR in)
{
  int i, j;
  printf("%d %d %d\n", in->m, in->n, in->NZ);
  for (i = 0; i < in->m; i++) {
    for (j = in->ai[i]; j < in->ai[i+1]; j++) {
      printf("%d %d %g\n", i, in->aj[j], in->data[j]);
    }
  }
}

/* This doesn't quite do the Gustavson algorithm. Here we're just
 * doing a pre-pass to allocate the correct amount of space. */
static void matmult_symbolic_csr_csr(const CSR A, const CSR B, CSR *out)
{
  CSR C;
  int *seen = malloc(B->n * sizeof(*seen));
  int i;
  alloc_csr(A->m, B->n, 0, &C);
  for (i = 0; i < B->n; i++) {
    seen[i] = -1;
  }
  C->ai[0] = 0;
  for (i = 0; i < A->m; i++) {
    int nnz = 0;
    int j;
    for (j = A->ai[i]; j < A->ai[i+1]; j++) {
      int Acol = A->aj[j];
      int k;
      for (k = B->ai[Acol]; k < B->ai[Acol+1]; k++) {
        int Bcol = B->aj[k];
        /* Have we seen this column against this row already? */
        if (seen[Bcol] != i) {
          /* No, that means that row i in C has another entry. */
          seen[Bcol] = i;
          nnz++;
        }
      }
    }
    C->ai[i+1] = C->ai[i] + nnz;
  }
  /* Allocate the rest of the space. */
  C->aj = malloc(C->ai[C->m] * sizeof(*C->aj));
  C->data = malloc(C->ai[C->m] * sizeof(*C->data));
  C->NZ = C->ai[C->m];
  *out = C;
  free(seen);
}

static void matmult_numeric_csr_csr(const CSR A, const CSR B, CSR C)
{
  double *colsums = calloc(C->n, sizeof(*colsums));
  /* Poor man's hashtable */
  int *nonzero_cols = malloc(C->n * sizeof(*nonzero_cols));
  int ptr = 0;
  int i;

  for (i = 0; i < C->n; i++) {
    nonzero_cols[i] = -1;
  }

  for (i = 0; i < C->m; i++) {
    int len = 0;
    int start = -2;
    int j;
    /* For each row of A */
    for (j = A->ai[i]; j < A->ai[i+1]; j++) {
      int Acol = A->aj[j];
      double val = A->data[j];
      int k;
      /* Select the row of B corresponding the nonzero column entry in
       * A */
      for (k = B->ai[Acol]; k < B->ai[Acol+1]; k++) {
        int Bcol = B->aj[k];
        colsums[Bcol] += val * B->data[k];
        /* "Linked list" of which entries in colsums are nonzero */
        if (nonzero_cols[Bcol] == -1) {
          nonzero_cols[Bcol] = start;
          start = Bcol;
          len++;
        }
      }
    }
    /* Now if any colsums[col] is nonzero, then row i of C has a
     * nonzero in said column. We use a fake "sparse" data structure
     * to traverse this. Could instead use a hash table or similar,
     * but meh. */
    for (j = 0; j < len; j++) {
      int col = start;
      start = nonzero_cols[col];
      /* Zero out for next time round */
      nonzero_cols[col] = -1;
      C->aj[ptr] = col;
      C->data[ptr] = colsums[col];
      colsums[col] = 0.0;
      ptr++;
    }
  }
  free(nonzero_cols);
  free(colsums);
}

/*
 * Computes C = A*B.
 * C should be allocated by this routine.
 */
void optimised_sparsemm(const COO A, const COO B, COO *C)
{
  CSR Ac, Bc, Cc;
  convert_coo_to_csr(A, &Ac);
  convert_coo_to_csr(B, &Bc);
  matmult_symbolic_csr_csr(Ac, Bc, &Cc);
  matmult_numeric_csr_csr(Ac, Bc, Cc);
  free_csr(&Ac);
  free_csr(&Bc);
  convert_csr_to_coo(Cc, C);
  free_csr(&Cc);
}

static void matsum_csr_csr(const CSR A, const CSR B, const CSR C, CSR *out)
{
  CSR O;
  int i, ptr;
  int *seen = malloc(B->n * sizeof(*seen));
  int *which = malloc(B->n * sizeof(*which));
  /* Max possible space (no overlap of nonzero patterns) */
  alloc_csr(A->m, A->n, A->NZ + B->NZ + C->NZ, &O);

  for (i = 0; i < A->n; i++) {
    seen[i] = -1;
    which[i] = -1;
  }
  
#define add_row(in, out) do {                           \
    for (int j = (in)->ai[i]; j < (in)->ai[i+1]; j++) { \
      int col = (in)->aj[j];                            \
      if (seen[col] != i) {                             \
        /* New column */                                \
        seen[col] = i;                                  \
        which[col] = ptr;                               \
        (out)->aj[ptr] = col;                           \
        (out)->data[ptr] = (in)->data[j];               \
        ptr++;                                          \
      } else {                                          \
        (out)->data[which[col]] += (in)->data[j];       \
      }                                                 \
    }                                                   \
  } while (0)

  O->ai[0] = 0;
  ptr = 0;
  for (i = 0; i < A->m; i++) {
    add_row(A, O);
    add_row(B, O);
    add_row(C, O);
    O->ai[i+1] = ptr;
  }
#undef add_row
  free(seen);
  free(which);
  O->NZ = O->ai[O->m];
  *out = O;
}

/* Computes O = (A + B + C) (D + E + F).
 * O should be allocated by this routine.
 */
void optimised_sparsemm_sum(const COO A, const COO B, const COO C,
                            const COO D, const COO E, const COO F,
                            COO *O)
{
  CSR Ac, Bc, Cc, Dc, Ec, Fc, Oc;
  CSR left, right;
  convert_coo_to_csr(A, &Ac);
  convert_coo_to_csr(B, &Bc);
  convert_coo_to_csr(C, &Cc);
  convert_coo_to_csr(D, &Dc);
  convert_coo_to_csr(E, &Ec);
  convert_coo_to_csr(F, &Fc);

  matsum_csr_csr(Ac, Bc, Cc, &left);
  free_csr(&Ac);
  free_csr(&Bc);
  free_csr(&Cc);
  matsum_csr_csr(Dc, Ec, Fc, &right);

  free_csr(&Dc);
  free_csr(&Ec);
  free_csr(&Fc);

  matmult_symbolic_csr_csr(left, right, &Oc);
  matmult_numeric_csr_csr(left, right, Oc);
  free_csr(&left);
  free_csr(&right);
  convert_csr_to_coo(Oc, O);
  free_csr(&Oc);
}
