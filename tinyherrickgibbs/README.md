# tinyherrickgibbs

Heap-free C99 Herrick–Gibbs velocity estimation for three closely spaced position observations. Include `tinyherrickgibbs/tinyherrickgibbs.h`, link `tinyherrickgibbs::tinyherrickgibbs`.

Inputs are positions in one inertial frame (metres), mu (metres³/second²), and **positive elapsed differences in seconds** `dt_before=t2-t1`, `dt_after=t3-t2`. Supplying differences directly avoids losing subsecond information by subtracting large Julian dates inside the library. The caller selects maximum adjacent angular separation and coplanarity error, both in (0,pi/2]. Returned velocity belongs to the middle position; geometry diagnostics are radians.

This truncated short-arc formula includes the mu/(12*r³) correction. It is an approximation whose orbital error depends on sampling interval and geometry; matching another implementation does not imply exact orbital velocity. Collinear, excessive-angle/noncoplanar, and invalid-time inputs are rejected. No heap, shared mutable state, or hidden history is used. Failure preserves output.

Validation: **15 frozen Vallado Python `hgibbs` cases** at commit `49df0479c950917c0a2c5ac31f2db6477bba57f6`, with five conic eccentricities and 1–12 second arcs. The exact elapsed differences represented by the reference's binary64 Julian dates are supplied to the C API. Maximum sampled velocity-component difference was **1.02e-9 m/s** on AppleClang 21; the static tolerance is 2e-8 m/s. Geometry, nonpositive/nonfinite timing and failure preservation are also checked. Reference code/raw outputs/generator stay ignored under `.validation/tinyherrickgibbs/`; tests are offline.

Independent implementation from the published Herrick–Gibbs equations; no AGPL source was translated. See LICENSE/NOTICE.
