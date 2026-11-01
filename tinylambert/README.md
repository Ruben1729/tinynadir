# tinylambert

Heap-free C99 **zero-revolution universal-variable Lambert solver**. It returns departure and arrival inertial velocity for given endpoint positions, positive flight time, and mu. Include `tinylambert/tinylambert.h`; link `tinylambert::tinylambert` (only `tinystumpff` is required).

Position is metres, velocity metres/second, time seconds, mu metres³/second². The endpoint states share one inertial frame. `TINYLAMBERT_SHORT_WAY` selects the transfer angle below pi; `LONG_WAY` selects the angle above pi. These geometric choices alone do not mean prograde/retrograde relative to an independently selected reference axis.

This profile has **no Battin or multirevolution branch API**. Collinear/nearly collinear geometry (unit-vector cross magnitude <1e-10) is rejected. Nondimensional flight time is restricted to [1e-6,1e8]. The time equation is nondimensionalized, bracketed, and bisected with a 160-iteration limit; a cancellation-resistant expression handles negative universal arguments. Failure preserves output. No heap, shared mutable state, or hidden gravitational/environment model is used.

Validation compares **18 frozen Vallado universal Lambert cases**, pinned at `49df0479c950917c0a2c5ac31f2db6477bba57f6`, with `nrev=0`, short/long paths, varied endpoint geometry, and 600/1800/5400-second coasts. Maximum sampled velocity-component difference was **8.74e-11 m/s** on AppleClang 21. Static tests use a 2e-6 m/s limit and verify endpoint energy consistency, bounded iteration counts, singularity/time policy, and output preservation. Raw oracle outputs/generator are ignored in `.validation/tinylambert/`; normal tests are offline C99.

Implementation was written independently from the universal-variable Lambert equations; no AGPL source was translated. This is a transfer boundary-value solver, not a replacement orbit propagator or controller. See LICENSE/NOTICE.
