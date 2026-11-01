# tinygravity

Cartesian spherical-harmonic gravity with caller-selected unnormalized coefficients, degree/order 0–32 and a reference radius/gravitational parameter. Solid-harmonic recurrences with analytic forward-mode Cartesian derivatives remain nonsingular at the geographic poles. Includes central gravity when C00=1; potential sign is positive mu/r and acceleration is its gradient. No built-in Earth coefficient archive. Fully normalized geodesy coefficients can be converted explicitly. Runtime workspace is 8*(degree+1)*(degree+2)/2 doubles (384 bytes at degree2; 35,904 bytes at degree32); coefficients are caller-owned and can be const flash data. No heap or mutable globals. Arrays must be disjoint except the documented coefficient in-place conversion.

## Verification

Frozen cases compare an independently written Cartesian-derivative implementation with the pinned Vallado Python Montenbruck gravity oracle, used only outside git. Cases span central gravity and normalized synthetic harmonic fields through degree32; SI acceleration tolerance 3e-11 m/s2. Exact polar J2 and coefficient-normalization tests exercise physical limiting cases. The runtime contains no Vallado source. Oracle repository: https://github.com/CelesTrak/fundamentals-of-astrodynamics/tree/49df0479c950917c0a2c5ac31f2db6477bba57f6 .

Reference scripts, packages and raw outputs are in ignored `.validation/tinygravity/`. Regression fixtures are checked in and run offline without Python or reference libraries.

The reference calls `accel_mont` directly with C00=1. The upstream generic dispatcher adds another central term, so using it with this coefficient convention would double-count central gravity. The frozen central cases and a polar J2 invariant independently check that convention.
