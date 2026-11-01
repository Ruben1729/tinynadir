#ifndef TINYLINALG_H
#define TINYLINALG_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define TINYLINALG_MAX_DIM 16u

typedef enum {
    TINYLINALG_OK = 0,
    TINYLINALG_INVALID_ARGUMENT,
    TINYLINALG_SINGULAR,
    TINYLINALG_NOT_POSITIVE_DEFINITE,
    TINYLINALG_NO_CONVERGENCE,
    TINYLINALG_ARITHMETIC_ERROR
} tinylinalg_status;

/* All matrices are contiguous row-major arrays. Sizes must be in [1,16].
 * Workspace is caller-owned and must not overlap inputs or outputs. Output
 * arrays are preserved on failure; workspace contents are unspecified.
 * No heap, mutable globals, or platform I/O are used. */
double tinylinalg_dot3(const double a[3], const double b[3]);
double tinylinalg_norm3(const double v[3]);
tinylinalg_status tinylinalg_cross3(const double a[3], const double b[3], double out[3]);
tinylinalg_status tinylinalg_unit3(const double v[3], double out[3]);
tinylinalg_status tinylinalg_multiply(size_t rows, size_t inner, size_t cols,
    const double *a, const double *b, double *out, double *workspace, size_t count);
tinylinalg_status tinylinalg_transpose(size_t rows, size_t cols, const double *a,
    double *out, double *workspace, size_t count);
/* Packed partial-pivoting LU: P*A=L*U, unit diagonal of L is implicit.
 * pivots[i] is the row swap performed at elimination step i. Workspace n*n. */
tinylinalg_status tinylinalg_lu(size_t n, const double *a, double *lu,
    size_t *pivots, double *workspace, size_t count);
/* Workspace n. LU and pivots must come from a successful tinylinalg_lu call. */
tinylinalg_status tinylinalg_lu_solve(size_t n, const double *lu,
    const size_t *pivots, const double *b, double *x, double *workspace, size_t count);
/* Workspace n*n+n. */
tinylinalg_status tinylinalg_solve(size_t n, const double *a, const double *b,
    double *x, double *workspace, size_t count);
/* Lower-triangular factor, A=L*transpose(L); symmetry is checked. Workspace n*n. */
tinylinalg_status tinylinalg_cholesky(size_t n, const double *a, double *lower,
    double *workspace, size_t count);
/* Symmetric Jacobi eigensolver. Eigenvalues ascending; corresponding normalized
 * eigenvectors are columns of vectors. Workspace 2*n*n. At most 64 sweeps. */
tinylinalg_status tinylinalg_eigen_symmetric(size_t n, const double *a,
    double *values, double *vectors, double *workspace, size_t count);

#ifdef __cplusplus
}
#endif
#endif
