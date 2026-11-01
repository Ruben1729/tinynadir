/* SPDX-License-Identifier: MIT */
#include <float.h>

/* C99 compile-time guard; no runtime storage or executable instructions. */
typedef char tinynadir_requires_ieee_binary64_double[
    (sizeof(double) == 8 && FLT_RADIX == 2 && DBL_MANT_DIG == 53 &&
     DBL_MAX_EXP == 1024) ? 1 : -1];
