# tinyrendezvous

Heap-free C99 maneuver building blocks, with three explicitly selected profiles:

| Function | Meaning | Assumptions |
| --- | --- | --- |
| `tinyrendezvous_phase_circular` | Two tangential burns into/out of a phasing ellipse | Same circular coplanar initial/final orbit; integer target and chaser turns; caller-specified minimum apsis radius |
| `tinyrendezvous_target` | Coast a moving target, solve its intercept, and compute departure/arrival velocity matching impulses | Two-body target, zero-revolution short/long Lambert path; same initial epoch and inertial frame |
| `tinyrendezvous_edelbaum` | Edelbaum characteristic delta-v and ideal transfer duration | Quasi-circular low-thrust radius/inclination change, constant acceleration; no RAAN constraint, eclipses, mass loss, drag, or steering history |

All boundary distances, velocities, gravitational parameters, times and accelerations use SI; angles use radians. Phasing uses **target minus chaser** phase in `[-pi,pi]`. Signed departure/arrival tangential impulses accompany the positive total delta-v. A minimum radius of zero explicitly disables central-body clearance checking. Turn counts are bounded to 100,000 each.

Edelbaum inclination change is limited to 90 degrees and acceleration to 1% of the weaker endpoint gravity. This is an ideal maneuver estimate, not a flight guidance controller. Intercepts inherit the geometry/time limits of `tinykepler` and `tinylambert`; no multi-revolution Lambert search is implied. Failure preserves the output. No heap, files, mutable global state or caller workspace is required.

Build standalone with `cmake -S tinyrendezvous -B build/rendezvous`, `cmake --build build/rendezvous`, and `ctest --test-dir build/rendezvous`. Link `tinyrendezvous::tinyrendezvous`; its Kepler/Lambert dependencies are supplied transitively.

## Verification

Frozen offline C fixtures compare 45 phasing cases with Vallado's `rendezvous_coplanar`, 11 moving-target intercepts with Vallado's Kepler/universal Lambert combination, and 20 Edelbaum cases with poliastro's `change_a_inc` implementation. Vallado reference revision: `49df0479c950917c0a2c5ac31f2db6477bba57f6`; poliastro revision: `21fd7719e89a7d22b4eac63141a60a7f1b01768c`.

Observed maximum delta-v differences on the development host were 4.61e-12 m/s for phasing, 6.70e-10 m/s for Edelbaum, and 1.28e-10 m/s for the transfer departure velocity. Static tolerances are 1e-7 m/s for maneuver estimates, 1e-6 m/s for intercept velocities, and 1 mm for target positions. Every intercept is also propagated to verify arrival position closure. One additional attempted Vallado Lambert reference case reported nonconvergence and returned zero-velocity sentinels; that invalid oracle output was recorded in ignored validation material and was not frozen as a correct answer. This does not imply that all possible Lambert cases were validated.

The C implementation is independently written from orbital mechanics equations and composes this repository's existing propagator and Lambert solver. Reference implementations and generation outputs remain in ignored `.validation/tinyrendezvous/`; tests require only the checked-in numeric fixtures.
