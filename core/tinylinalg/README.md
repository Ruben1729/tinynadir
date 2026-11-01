# tinylinalg

Part of [Core](../README.md) in [TinyNadir](../../README.md).

Heap-free C99 linear algebra for dimensions 1–16. Matrices are contiguous and row-major. The library provides vector products and normalization, matrix products/transposes, partial-pivoting LU and solves, Cholesky, and a symmetric Jacobi eigensolver.

Include `tinylinalg/tinylinalg.h` and link `tinylinalg::tinylinalg`. Each routine documents its workspace count in doubles. Workspace must be separate from input/output arrays. An output may alias an input where its shape permits it; failed operations preserve outputs, while workspace is scratch. No allocation, I/O, or mutable global state occurs. The dimension limit bounds matrix work and storage; the eigensolver stops after at most 64 sweeps.

```c
double a[4] = {4, 1, 1, 3}, b[2] = {1, 2}, x[2], work[6];
tinylinalg_status status = tinylinalg_solve(2, a, b, x, work, 6);
```

LU uses partial pivoting and returns a singular status for a zero pivot. The routines do not estimate condition numbers; callers must interpret ill-conditioned systems according to their measurement/model requirements. Cholesky accepts symmetric positive-definite matrices only. Symmetric eigenvalues are sorted ascending and eigenvectors occupy matrix columns; signs and repeated-eigenvalue bases are not unique.

The solve workspace is `n*n+n` doubles, Cholesky/LU `n*n`, and the eigensolver `2*n*n`. Additional bounded stack includes 16 pivots or eigenvalues, loop scalars, and platform math frames. Library code does not store maximum-size matrices in global state.

## Verification

Frozen fixtures were generated with **NumPy 2.5.4 / SciPy 1.18.1**, using seeded independent matrices of dimensions 1, 2, 3, 4, 6, 9, 12, and 16 (24 cases). Results compare solves, Cholesky factors, and eigenvalues, with eigenpair residual/orthogonality checks rather than eigenvector-sign comparisons. Tests also exercise row pivots, singular and indefinite matrices, asymmetry, invalid dimensions, NaNs, overflow, huge-vector normalization, and aliasing.

Tolerances are 2e-12 relative-to-max(1, expected) for solves/Cholesky and 3e-12 for eigenvalues/residuals. Development generators, installed Python packages, and raw comparison outputs are kept in ignored `.validation/tinylinalg/`; retained tests are standalone C and need no Python or external library.

References: [NumPy linear algebra](https://numpy.org/doc/stable/reference/routines.linalg.html), [SciPy linear algebra](https://docs.scipy.org/doc/scipy/reference/linalg.html). New implementation contributions use the root MIT license; see `NOTICE` for provenance.
