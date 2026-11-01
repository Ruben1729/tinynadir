# tinymeanosculating

Part of [Orbit](../README.md) in [TinyNadir](../../README.md).

Heap-free C99 conversion between an instantaneous Cartesian orbit and a **specified numerical centered-one-orbit average** under central-plus-J2 gravity. The inverse uses bounded fixed-point corrections until applying the same filter recovers the requested mean coordinates. It returns a convergence error when that cannot be achieved within the requested work limits.

This module implements its documented numerical filter. Brouwer–Lyddane, DSST, SGP4 and TLE mean elements have different definitions and cannot be passed to it as though they were interchangeable. In particular, this is not a TLE-to-osculating adapter.

## Definition

At the center epoch, convert the Cartesian state into prograde equinoctial coordinates `a,f,g,h,k,lambda`, where `f=e cos(Omega+omega)`, `g=e sin(Omega+omega)`, `h=tan(i/2) cos(Omega)`, `k=tan(i/2) sin(Omega)`, and `lambda=Omega+omega+M`. This retains useful coordinates for circular/equatorial cases; the retrograde singularity is excluded by the supported profile.

Let `a0` be the initial osculating semimajor axis, `n0=sqrt(mu/a0^3)` and `T0=2pi/n0`. Integrate the central-plus-J2 Cartesian equations forward and backward from the epoch. Sample uniformly over `[-T0/2,+T0/2]` using the caller's even number of intervals and trapezoidal endpoint half-weights. Average `a,f,g,h,k` directly. For longitude, average the locally unwrapped residual `lambda(t)-lambda(0)-n0*t`, add `lambda(0)` and wrap to `[0,2pi)`. The wrap is safe for the selected weak-J2, one-orbit profile; it does not attempt to infer a multi-revolution phase history.

These operations, including the finite sample grid, define the returned mean. Changing the sample count changes the discrete filter. Finite integration steps add numerical error; they do not change the intended equations. This time average is not a canonical perturbation-theory transformation. No universal accuracy equivalence with other mean-element theories is implied.

The inverse initially interprets the supplied coordinates as osculating coordinates, applies the filter, and corrects the six coordinate residuals. It is a local inverse for the weak-perturbation profile. The convergence norm divides the semimajor-axis residual by the requested mean semimajor axis and uses absolute residuals for the dimensionless coordinates and longitude. The caller chooses its tolerance and iteration budget.

## Profile and resources

Inputs use SI and radians; `f,g,h,k` are dimensionless. The gravity symmetry pole is a fixed inertial +z axis. Supported states are elliptic, eccentricity at most 0.3, inclination at most 170 degrees, and osculating pericenter outside the explicit body equatorial radius. Every sampled state's osculating elements must satisfy this profile. `|J2|` is limited to 0.005. The integrator also rejects trajectories that enter the body. This is a weak-J2 conversion model; other harmonics, drag, third bodies and thrust are omitted.

The caller supplies 16–512 even trapezoid intervals, a positive RK4 maximum step, a total integration-step limit up to 1,000,000, an inverse limit of 1–16 averaging passes and an inverse tolerance of 1e-12–1e-4. The total step budget covers every pass of an inverse call. Both functions preserve results and statistics on failure. Workspace contents are unspecified after a call.

Caller workspace is exactly 48 doubles (384 bytes on the required binary64 profile). Additional stack use is bounded and compiler dependent; no sample-sized arrays, heap allocations, mutable globals or runtime data files are used. This conversion integrates at least one orbit and may do several passes: compute it when needed, rather than assuming it is an inexpensive attitude-control-loop operation. Reducing step size or increasing sample count trades runtime for accuracy; choose them against your mission's error budget.

Build with `cmake -S orbit/tinymeanosculating -B build/mean`, `cmake --build build/mean`, and `ctest --test-dir build/mean`. Link `tinymeanosculating::tinymeanosculating`; dependencies `tinycowell`, `tinyelements` and `tinyanomaly` are transitive.

## Verification

16 frozen independent SciPy 1.18.1 DOP853 integrations plus NumPy vector-based element conversions implement the same explicit centered filter using a different integration algorithm. Cases include Earth LEO/MEO, a Mars parameter model, zero and signed J2, circular/equatorial orbits, eccentricity up to 0.25 and inclination up to 2.9 radians. Sample grids contain 64, 128 and 256 intervals.

With maximum RK4 step 2 seconds, the development-host maximum differences were 4.77 micrometers in averaged semimajor axis and 3.72e-13 in other coordinates. Applying the inverse to the independent mean fixtures recovered the original states within 12.3 micrometers position and 5.79e-9 m/s velocity, using at most six passes and 86,528 integration steps. Static regression tolerances allow 1 mm in mean semimajor axis, 1e-10 in other mean coordinates, 2 mm in recovered position and 2e-6 m/s in recovered velocity. These quantify implementation agreement for the selected fixtures, not neglected-force error or accuracy of a different definition of mean elements.

Tests also check filter/inverse consistency, invalid sample counts, exhausted work budgets, nonconvergence, nonfinite input and output preservation. Reference scripts and raw results live only in ignored `.validation/tinymeanosculating/`; checked-in numeric fixtures make tests independent of Python and external libraries.
