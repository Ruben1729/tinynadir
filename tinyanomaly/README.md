# tinyanomaly

Heap-free C99 Kepler-equation solving and true-anomaly conversion for elliptic, hyperbolic, and parabolic conics. Include `tinyanomaly/tinyanomaly.h` and link `tinyanomaly::tinyanomaly`.

`tinyanomaly_from_mean(e, mean, &solution)` returns eccentric anomaly E (`e<1`), hyperbolic anomaly H (`e>1`), or Barker variable D (`e==1`), together with true anomaly. `tinyanomaly_from_true` performs the inverse conversion. Angles are radians. For parabolas the mean argument is dimensionless `D + D^3/3`, not a mean angle. Elliptic inputs are reduced to the principal interval with `remainder(mean, 2*pi)`.

The supported input profile is `0 <= e <= 1000` and `abs(mean) <= 1e12`. True anomaly must be in the physical principal domain; hyperbolic angles beyond the asymptote and the parabolic endpoint are rejected. A safeguarded Newton solve is bounded at 100 iterations. Series evaluation preserves the small cubic term near the parabolic limit. Failure preserves outputs.

There is no heap, shared mutable state, file I/O, or dependence on an orbit propagator. These are anomaly equations, not a complete alternative to SGP4. State propagation belongs to `tinykepler` or another module.

Validation: 80 mean-anomaly cases and 35 inverse conversions were frozen from mpmath 1.4.1 at 90 decimal digits, using its independent bracketing root solver. Coverage includes `e=0`, near-parabolic elliptic/hyperbolic values, `e=1`, large eccentricity, signed/zero times, and invalid domains. Maximum sampled anomaly/true-anomaly errors were **3.56e-15 / 2.23e-16 radians** on AppleClang 21. Static tests also verify failure preservation. Raw oracle output and the generator are ignored under `.validation/tinyanomaly/`; normal tests require only C99 and libm.

References: published elliptic/hyperbolic Kepler and Barker equations, e.g. Vallado, *Fundamentals of Astrodynamics and Applications*; [mpmath root finding](https://mpmath.org/doc/current/calculus/optimization.html). Implementation was written from mathematics without translating AGPL reference source. See LICENSE and NOTICE.
