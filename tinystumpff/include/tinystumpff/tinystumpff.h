/* SPDX-License-Identifier: MIT */
#ifndef TINYSTUMPFF_H
#define TINYSTUMPFF_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { TINYSTUMPFF_OK, TINYSTUMPFF_INVALID, TINYSTUMPFF_RANGE } tinystumpff_status;
typedef struct { double c[4]; } tinystumpff_values;
/* C_k(z) = sum_n (-z)^n/(2*n+k)!, k=0..3. No allocation.
   For z<0, sqrt(-z) must be <=700. Failure preserves output. */
tinystumpff_status tinystumpff_evaluate(double z, tinystumpff_values *output);
#ifdef __cplusplus
}
#endif
#endif
