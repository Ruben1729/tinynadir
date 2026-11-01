# tinyelliptic

Part of [Core](../README.md) in [TinyNadir](../../README.md).

Real Legendre elliptic integrals in heap-free C99: complete/incomplete first and second kinds, inverse second integral, and signed ellipse arc length. Include `tinyelliptic/tinyelliptic.h`; link `tinyelliptic::tinyelliptic`.

The parameter is **m=k²**, not modulus k, and this profile accepts m in [0,1]. Amplitudes use radians and are bounded by 1e6 in magnitude. First-kind integrals at m=1 are singular when the path reaches ±pi/2; second-kind integrals remain finite and use the exact degenerate expression. Complex parameters, m<0, m>1, and third-kind integrals are outside this profile.

Carlson symmetric RF/RD duplication plus bounded correction series evaluates the principal incomplete integrals; periodic extension handles other real amplitudes. Each duplication is bounded at 64 iterations. Inversion uses a bracketed Newton method with at most 80 iterations; near m=1 and a quarter turn the inverse is intrinsically ill-conditioned. Exact mathematical quarter-turn endpoints use the complete-integral convention. Output is preserved for invalid, out-of-range, singular, or arithmetic/convergence failures.

`tinyelliptic_arc(a,b,start,end)` integrates the parametric ellipse `x=a*cos(theta),y=b*sin(theta)` with positive semiaxes in any common length unit. It is signed for increasing/decreasing parameter. Period terms are differenced before summation to reduce cancellation; circles use the direct radius-times-angle identity. No heap, mutable global state, I/O, variable-sized arrays, or caller workspace is needed; bounded scalar local scratch and libm stack are used.

## Verification

**SciPy 1.18.1** supplies 210 incomplete integral cases, 11 complete pairs, 49 inverse cases (brentq), and 64 independent numerical ellipse arcs (quad), with NumPy 2.5.4 and seed 776712. Coverage includes negative/multiple-turn amplitudes, zero, tiny m, near-one m, m=1, and quarter-turn boundaries. Exact quarter turns use `ellipk` as the first-kind oracle; `ellipkinc`'s argument-reduction roundoff near m=1 is amplified by the singular endpoint condition and is not used to redefine the complete integral.

Maximum sampled first-integral difference is 4.26e-14. Integral tolerances are 8e-14 times max(1, expected magnitude); inverse angle tolerance is 5e-12 for the sampled noncritical inverse cases, with an additional forward round-trip check. Arc tolerance is 5e-12 times max(1, magnitude). Singular/invalid/range failures and circle identity are also tested. These are numeric implementation comparisons, not a universal inverse forward-error bound near a singular condition.

The oracle generator, raw JSON, and comparison logs remain ignored under `.validation/tinyelliptic/`; retained C fixtures/tests run offline. References: [NIST DLMF duplication algorithms](https://dlmf.nist.gov/19.36), [SciPy elliptic special functions](https://docs.scipy.org/doc/scipy/reference/special.html#elliptic-functions-and-integrals). The implementation is independently written; see `NOTICE` and `LICENSE`.
