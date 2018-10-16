#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include "utils.h"


#ifdef _MSC_VER
// -*- C++ -*-
// $Id: drand48.src,v 1.1.1.1.2.1 2004/04/28 06:02:54 garren Exp $
// ---------------------------------------------------------------------------
//
// This code is based on a code extracted from GNU C Library 2.1.3 with
// a purpose to provide drand48() on Windows NT.
//


#define	__LITTLE_ENDIAN	1234
#define	__BIG_ENDIAN	4321

#ifdef __BIG_ENDIAN__
#define __BYTE_ORDER __BIG_ENDIAN     /* powerpc-apple is big-endian. */
#else
#define __BYTE_ORDER __LITTLE_ENDIAN  /* i386 is little-endian.  */
#endif
#define __FLOAT_WORD_ORDER __BYTE_ORDER

#define IEEE754_DOUBLE_BIAS	0x3ff /* Added to exponent.  */

typedef unsigned long long int u_int64_t;

union ieee754_double {
  double d;

    /* This is the IEEE 754 double-precision format.  */
  struct {
#if __BYTE_ORDER == __BIG_ENDIAN
    unsigned int negative:1;
    unsigned int exponent:11;
    /* Together these comprise the mantissa.  */
    unsigned int mantissa0:20;
    unsigned int mantissa1:32;
#endif				/* Big endian.  */
#if __BYTE_ORDER == __LITTLE_ENDIAN
#if __FLOAT_WORD_ORDER == BIG_ENDIAN
    unsigned int mantissa0:20;
    unsigned int exponent:11;
    unsigned int negative:1;
    unsigned int mantissa1:32;
#else
    /* Together these comprise the mantissa.  */
    unsigned int mantissa1:32;
    unsigned int mantissa0:20;
    unsigned int exponent:11;
    unsigned int negative:1;
#endif
#endif				/* Little endian.  */
  } ieee;
};

/* Data structure for communication with thread safe versions.  */
struct drand48_data {
  unsigned short int x[3];     /* Current state.  */
  unsigned short int a[3];     /* Factor in congruential formula.  */
  unsigned short int c;	       /* Additive const. in congruential formula.  */
  unsigned short int old_x[3]; /* Old state.  */
  int init;                    /* Flag for initializing.  */
};

typedef struct drand48_data drand48_data;
/* Global state for non-reentrant functions. */
static drand48_data libc_drand48_data;

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
static int drand48_iterate (unsigned short int xsubi[3], drand48_data *buffer)
{
  u_int64_t X, a, result;

  /* Initialize buffer, if not yet done.  */
  if (!buffer->init) {
#if (USHRT_MAX == 0xffffU)
      buffer->a[2] = 0x5;
      buffer->a[1] = 0xdeec;
      buffer->a[0] = 0xe66d;
#else
      buffer->a[2] = 0x5deecUL;
      buffer->a[1] = 0xe66d0000UL;
      buffer->a[0] = 0;
#endif
      buffer->c = 0xb;
      buffer->init = 1;
  }

  /* Do the real work.  We choose a data type which contains at least
     48 bits.  Because we compute the modulus it does not care how
     many bits really are computed.  */

  if (sizeof (unsigned short int) == 2) {
    X = (u_int64_t)xsubi[2] << 32 | (u_int64_t)xsubi[1] << 16 | xsubi[0];
    a = ((u_int64_t)buffer->a[2] << 32 | (u_int64_t)buffer->a[1] << 16
	 | buffer->a[0]);

    result = X * a + buffer->c;

    xsubi[0] = result & 0xffff;
    xsubi[1] = (result >> 16) & 0xffff;
    xsubi[2] = (result >> 32) & 0xffff;
  }else{
    X = (u_int64_t)xsubi[2] << 16 | xsubi[1] >> 16;
    a = (u_int64_t)buffer->a[2] << 16 | buffer->a[1] >> 16;

    result = X * a + buffer->c;

    xsubi[0] = result >> 16 & 0xffffffffl;
    xsubi[1] = result << 16 & 0xffff0000l;
  }

  return 0;
}


static int erand48_r (unsigned short int xsubi[3],
                      drand48_data *buffer,
                      double *result)
{
  union ieee754_double temp;

  /* Compute next state.  */
  if (drand48_iterate (xsubi, buffer) < 0) return -1;

  /* Construct a positive double with the 48 random bits distributed over
     its fractional part so the resulting FP number is [0.0,1.0).  */

#if USHRT_MAX == 65535
  temp.ieee.negative = 0;
  temp.ieee.exponent = IEEE754_DOUBLE_BIAS;
  temp.ieee.mantissa0 = (xsubi[2] << 4) | (xsubi[1] >> 12);
  temp.ieee.mantissa1 = ((xsubi[1] & 0xfff) << 20) | (xsubi[0] << 4);
#elif USHRT_MAX == 2147483647
  temp.ieee.negative = 0;
  temp.ieee.exponent = IEEE754_DOUBLE_BIAS;
  temp.ieee.mantissa0 = (xsubi[1] << 4) | (xsubi[0] >> 28);
  temp.ieee.mantissa1 = ((xsubi[0] & 0xfffffff) << 4);
#else
# error Unsupported size of short int
#endif

  /* Please note the lower 4 bits of mantissa1 are always 0.  */
  *result = temp.d - 1.0;

  return 0;
}

static double drand48 ()
{
  double result;
  (void) erand48_r (libc_drand48_data.x, &libc_drand48_data, &result);
  return result;
}

// ---------------------------------------------------------------------------
static void srand48 (long int seedval)
{
  (void) srand48_r (seedval, &libc_drand48_data);
}

// ---------------------------------------------------------------------------
static int srand48_r (long int seedval, drand48_data *buffer)
{
  /* The standards say we only have 32 bits.  */
  if (sizeof (long int) > 4)
    seedval &= 0xffffffffl;

#if USHRT_MAX == 0xffffU
  buffer->x[2] = seedval >> 16;
  buffer->x[1] = seedval & 0xffffl;
  buffer->x[0] = 0x330e;

  buffer->a[2] = 0x5;
  buffer->a[1] = 0xdeec;
  buffer->a[0] = 0xe66d;
#else
  buffer->x[2] = seedval;
  buffer->x[1] = 0x330e0000UL;
  buffer->x[0] = 0;

  buffer->a[2] = 0x5deecUL;
  buffer->a[1] = 0xe66d0000UL;
  buffer->a[0] = 0;
#endif
  buffer->c = 0xb;
  buffer->init = 1;

  return 0;
}

// ---------------------------------------------------------------------------
static int seed48_r (unsigned short int seed16v[3], drand48_data *buffer)
{
  /* Save old value at a private place to be used as return value.  */
  memcpy (buffer->old_x, buffer->x, sizeof (buffer->x));

  /* Install new state.  */
#if USHRT_MAX == 0xffffU
  buffer->x[2] = seed16v[2];
  buffer->x[1] = seed16v[1];
  buffer->x[0] = seed16v[0];

  buffer->a[2] = 0x5;
  buffer->a[1] = 0xdeec;
  buffer->a[0] = 0xe66d;
#else
  buffer->x[2] = (seed16v[2] << 16) | seed16v[1];
  buffer->x[1] = seed16v[0] << 16;
  buffer->x[0] = 0;

  buffer->a[2] = 0x5deecUL;
  buffer->a[1] = 0xe66d0000UL;
  buffer->a[0] = 0;
#endif
  buffer->c = 0xb;
  buffer->init = 1;

  return 0;
}
// ---------------------------------------------------------------------------
static unsigned short int * seed48 (unsigned short int seed16v[3])
{
  (void) seed48_r (seed16v, &libc_drand48_data);
  return libc_drand48_data.old_x;
}


#endif  /* _MSC_VER */
/*
 * Allocate a dense matrix
 * m - number of rows
 * n - number of columns
 * dense - newly allocated matrix.
 */
void alloc_dense(int m, int n, double **dense)
{
  *dense = malloc(m*n*sizeof(**dense));
}

/*
 * Free a dense matrix
 * dense - dense matrix, may be NULL
 */
void free_dense(double **dense)
{
    if (!*dense) {
        return;
    }
    free(*dense);
    *dense = NULL;
}

/*
 * Zero a dense matrix
 * m - number of rows
 * n - number of columns
 * dense - matrix to zero.
 */
void zero_dense(int m, int n, double *dense)
{
    int i, j;
    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            dense[j*m + i] = 0;
        }
    }
}

/*
 * Allocate a sparse matrix in coordinate format.
 * m - number of rows
 * n - number of columns
 * NZ - number of nonzeros
 * sparse - newly allocated matrix.
 */
void alloc_sparse(int m, int n, int NZ, COO *sparse)
{
    COO sp = calloc(1, sizeof(struct _p_COO));
    sp->m = m;
    sp->n = n;
    sp->NZ = NZ;
    sp->coords = calloc(NZ, sizeof(struct coord));
    sp->data = calloc(NZ, sizeof(double));
    *sparse = sp;
}

/*
 * Free a sparse matrix.
 * sparse - sparse matrix, may be NULL
 */
void free_sparse(COO *sparse)
{
    COO sp = *sparse;
    if (!sp) {
        return;
    }
    free(sp->coords);
    free(sp->data);
    free(sp);
    *sparse = NULL;
}

/*
 * Convert a sparse matrix to dense format in column major format.
 *
 * sparse - The sparse matrix to convert
 * dense - pointer to output dense matrix (will be allocated)
 */
void convert_sparse_to_dense(const COO sparse, double **dense)
{
    int n;
    int i, j;
    alloc_dense(sparse->m, sparse->n, dense);
    zero_dense(sparse->m, sparse->n, *dense);
    for (n = 0; n < sparse->NZ; n++) {
        i = sparse->coords[n].i;
        j = sparse->coords[n].j;
        (*dense)[j * sparse->m + i] = sparse->data[n];
    }
}

/*
 * Convert a dense matrix in column major format to sparse.
 * Entries with absolute value < 1e-15 are flushed to zero and not
 * stored in the sparse format.
 *
 * dense - the dense array
 * m - number of rows
 * n - number of columns
 * sparse - output sparse matrix (allocated by this routine)
 */
void convert_dense_to_sparse(const double *dense, int m, int n,
                             COO *sparse)
{
    int i, j, NZ;
    COO sp;
    NZ = 0;
    /* Figure out how many nonzeros we're going to have. */
    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            double val = dense[j*m + i];
            if (fabs(val) > 1e-15) {
                NZ++;
            }
        }
    }
    alloc_sparse(m, n, NZ, &sp);

    NZ = 0;
    /* Fill up the sparse matrix */
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            double val = dense[j*m + i];
            if (fabs(val) > 1e-15) {
                sp->coords[NZ].i = i;
                sp->coords[NZ].j = j;
                sp->data[NZ] = val;
                NZ++;
            }
        }
    }
    *sparse = sp;
}

/*
 * Create a random sparse matrix
 *
 * m - number of rows
 * n - number of columns
 * frac - fraction of entries that should be nonzero
 * sparse - newly allocated random matrix.
 */
void random_matrix(int m, int n, double frac, COO *sparse)
{
    int i, j;
    double *d;
    alloc_dense(m, n, &d);
    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            if (drand48() < frac) {
                d[j*m + i] = drand48();
            }
        }
    }
    convert_dense_to_sparse(d, m, n, sparse);
    free_dense(&d);
}

/*
 * Read a sparse matrix from a file.
 *
 * file - The filename to read
 * sparse - The newly read sparse matrix (allocated here)
 */
void read_sparse(const char *file, COO *sparse)
{
    COO sp;
    int i, j, k, m, n, NZ;
    double val;
    int c;
    FILE *f = fopen(file, "r");
    if (!f) {
        fprintf(stderr, "Unable to open %s for reading.\n", file);
        exit(1);
    }
    c = fscanf(f, "%d %d %d\n", &m, &n, &NZ);
    if (c != 3) {
        fprintf(stderr, "File format incorrect on line 1, expecting 3 integers, got %d\n", c);
        fclose(f);
        exit(1);
    }
    if (NZ > m*n) {
        fprintf(stderr, "More nonzeros (%d) than matrix entries (%d x %d)!\n", NZ, m, n);
        fclose(f);
        exit(1);
    }
    alloc_sparse(m, n, NZ, &sp);
    k = 0;
    while ((c = fscanf(f, "%d %d %lg\n", &i, &j, &val)) == 3) {
        if (k >= NZ) {
            fprintf(stderr, "File has nonzero lines than expected (%d)\n", NZ);
            fclose(f);
            free_sparse(&sp);
            exit(1);
        }
        if (i >= m || j >= n) {
            fprintf(stderr, "Entry on line %d incorrect, index (%d, %d) out of bounds for %d x %d matrix\n", k + 2, i, j, m, n);
            fclose(f);
            free_sparse(&sp);
            exit(1);
        }
        sp->coords[k].i = i;
        sp->coords[k].j = j;
        sp->data[k] = val;
        k++;
    }

    if (k != NZ) {
        fprintf(stderr, "File has fewer lines (%d) than expected (%d)\n",
                k, NZ);
        fclose(f);
        free_sparse(&sp);
        exit(1);
    }
    *sparse = sp;
    fclose(f);
}

/*
 * Write a sparse matrix to a file.
 *
 * f - The file handle.
 * sp - The sparse matrix to write.
 */
void write_sparse(FILE *f, COO sp)
{
    int i;
    fprintf(f, "%d %d %d\n", sp->m, sp->n, sp->NZ);
    for (i = 0; i < sp->NZ; i++) {
        fprintf(f, "%d %d %g\n", sp->coords[i].i, sp->coords[i].j, sp->data[i]);
    }
}

/*
 * Print a sparse matrix to stdout
 *
 * sp - The sparse matrix to print.
 */
void print_sparse(COO sp)
{
    write_sparse(stdout, sp);
}
