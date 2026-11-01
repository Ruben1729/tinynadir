# tinygauss

Classical Gauss angles-only initial orbit determination from three inertial optical sight lines and observer positions. This profile computes the eighth-degree geocentric-radius polynomial, reconstructs positive sight-line ranges using the truncated f/g coefficients, and obtains middle-epoch velocity using Gibbs. It is a **Gauss initial estimate**, not an exact three-observation two-body fit or a propagator.

Include `tinygauss/tinygauss.h`; link `tinygauss::tinygauss`. Supply mu [m³/s²], t1<0<t3 [s] relative to the middle epoch, three line-of-sight triples, and three observer-position triples [m] in one inertial frame. Directions are normalized internally. Each time interval is bounded by 1e6 s for arithmetic range, but the truncated model is intended for short arcs; this numeric bound is not an accuracy guarantee.

Multiple positive polynomial roots can represent different candidate orbits. A zero radius hint selects the largest admissible root; a positive hint selects the admissible root nearest that geocentric radius. The result exposes selected radius, ranges, and admissible-root count. All ranges must be positive and the reconstructed positions must admit the Gibbs geometry. No altitude, visibility, sensor uncertainty, or mission-prior constraint is silently imposed. Nearly coplanar/degenerate sight-line matrices are rejected at an absolute unit-column determinant <=1e-12. An application may need stronger quality checks, observation redundancy, or batch differential correction.

The caller workspace is 480 bytes for binary64 doubles: degree-eight polynomial roots/workspace and a 3×3 linear solve workspace. Other constant matrices/vectors add stack. Root solving is bounded by tinypolynomial's 256 iterations. There is no heap, I/O, mutable global state, or external reference runtime. Failure preserves the result; workspace is scratch and must be separate.

## Verification

48 frozen outputs from Vallado's existing Python `astro.iod.angles.gauss` at **CelesTrak/fundamentals-of-astrodynamics commit 49df0479c950917c0a2c5ac31f2db6477bba57f6**, using NumPy 2.5.4/SciPy 1.18.1, compare the initial Gauss/Gibbs estimate. Seed 564712 varies observer location, orbit orientation/eccentricity/radius, and equal/unequal observation timing. Synthetic observations were created independently with SciPy DOP853 and Earth-rotating station coordinates. Inputs preserve the oracle's actual split-Julian-date time differences.

Maximum sampled differences are 2.38e-6 m in position and 1.24e-6 m/s in velocity. Static tolerances are 1e-4 m and 2e-5 m/s. Tests also check positive ranges, explicit radius-root selection, parallel directions, invalid time/mu/NULLs, and failure preservation. Agreement with this approximate reference does not establish orbit truth accuracy; use the result as an initial estimate and validate/refine for the mission.

In the synthetic reference sample, the Gauss initial estimate itself differed from the generating two-body truth by up to approximately 2.77 km and 9.62 m/s. These sampled truncation/initial-estimation differences are distinct from the tiny implementation's agreement with Vallado and are not universal error bounds.

Raw oracle implementation/results/scripts remain ignored under `.validation/tinygauss/`, `.validation/vallado-reference/`, and `.validation/adcs-reference/`; only static numeric inputs/outputs are committed. The runtime was independently written from the classical determinant/polynomial/Gibbs mathematics rather than translated from upstream code. See `NOTICE` and `LICENSE`.
