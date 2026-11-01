# tinylaplace

Classical Laplace angles-only orbit determination using instantaneous line-of-sight curvature and observer kinematics. This is mathematically distinct from Gauss: it eliminates range/range-rate from the line-of-sight acceleration equation and solves an eighth-degree radius polynomial.

Include `tinylaplace/tinylaplace.h`; link `tinylaplace::tinylaplace`. The main API receives a unit direction and its first/second derivatives [1/s,1/s²], observer position/velocity/acceleration [m,m/s,m/s²], and mu [m³/s²], all at the same epoch in one inertial frame. Its derivatives must correspond to the supplied direction; no hidden Earth-rotation model is assumed.

`tinylaplace_three` is a convenience profile for **equally spaced** observations at [-spacing,0,+spacing]. It normalizes the three directions and uses central first/second differences, with the actual middle observer kinematics supplied separately. Finite differentiation can amplify noise and introduce large truncation errors. Unequally spaced observations require externally fitted direction derivatives and the main API; they must not be passed as though equally spaced.

A zero radius hint selects the largest positive-range real root; a positive hint selects the nearest admissible geocentric radius. Root count, selected radius, range, and range-rate are returned. The curvature determinant must exceed 1e-12 times the product of derivative norms. Multiple roots and weak geometry still need application priors, more observations, or refinement. No Earth-altitude or observation-quality prior is silently imposed.

Caller workspace is 384 bytes for binary64 degree-eight complex roots and scratch; other fixed local vectors add stack. Polynomial solving is bounded by 256 iterations, with no heap, I/O, mutable globals, or reference-library dependency. Outputs are preserved on failure.

## Verification

24 frozen, symmetrically spaced cases compare with existing Vallado Python `astro.iod.angles.laplace` at **CelesTrak/fundamentals-of-astrodynamics commit 49df0479c950917c0a2c5ac31f2db6477bba57f6**, using NumPy 2.5.4/SciPy 1.18.1. Seed 564712 varies orbit and station geometry. The oracle's same-station mode is matched by explicit rotating-station velocity/acceleration and actual split-date intervals. Maximum sampled differences are 2.55e-7 m and 1.93e-8 m/s; static tolerances are 1e-4 m and 1e-5 m/s. Invalid timing/NULLs, zero curvature, and result preservation are tested.

These are comparisons of an initial-estimation method, not guarantees of truth accuracy from finite-difference noisy angles. Raw references/generators/results remain ignored under `.validation/tinylaplace/`, `.validation/vallado-reference/`, and `.validation/adcs-reference/`; normal tests use frozen numeric fixtures only. The C implementation was independently written from the classical Laplace equations, not copied/translated from upstream source. See `NOTICE` and `LICENSE`.

The finite-difference Laplace estimates in this sample differed from generating two-body truth by up to approximately 169 km and 458 m/s. This highlights the method's curvature/truncation sensitivity; tiny-to-Vallado agreement is not a truth-accuracy claim. More observations, better derivative fits, or a differential-correction stage are needed when the mission requires them.
