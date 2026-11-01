# tinystumpff

Heap-free C99 evaluation of the four universal-variable Stumpff functions, `C0` through `C3`. Include `tinystumpff/tinystumpff.h` and link `tinystumpff::tinystumpff`.

```c
tinystumpff_values c;
if (tinystumpff_evaluate(z, &c) == TINYSTUMPFF_OK) {
    /* c.c[2] and c.c[3] are the usual universal-variable C(z), S(z). */
}
```

The definition is `Ck(z) = sum((-z)^n / (2*n+k)!)`. A bounded power series handles small arguments; trigonometric or hyperbolic identities handle the remaining domain. Negative arguments require `sqrt(-z) <= 700` to avoid hyperbolic overflow. Nonfinite inputs and null outputs are rejected. Failure preserves the output.

The implementation is independent of upstream source and uses only caller output, small automatic working storage, and libm; no heap or mutable global state. It can be used by Kepler/Lambert solvers without importing other orbit algorithms. The ordinary build/test is offline.

Validation: 20 frozen cases from mpmath 1.4.1, evaluated with 90-digit trigonometric/hyperbolic functions, cover zero, both sides of the series transition, positive/negative arguments, and a negative argument at the overflow-policy boundary. Maximum observed error scaled by `max(1, abs(reference))` was **4.10e-16** on AppleClang 21. The static tests use a 2e-13 scaled limit and test failure preservation. Raw comparisons and their generator are ignored under `.validation/tinystumpff/`.

References: the Stumpff-function series definition and universal-variable formulations in Vallado, *Fundamentals of Astrodynamics and Applications*; [mpmath functions](https://mpmath.org/doc/current/functions/). The reference library is used only to capture test values. See LICENSE and NOTICE.
