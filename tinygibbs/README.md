# tinygibbs

Heap-free C99 Gibbs initial-orbit velocity determination from three position observations in one inertial frame. Link `tinygibbs::tinygibbs`; include `tinygibbs/tinygibbs.h`. Positions are metres and mu is metres³/second²; the returned velocity is metres/second at the **middle observation**.

`tinygibbs_estimate` also returns adjacent angular separations and signed coplanarity angle. The caller explicitly selects the maximum accepted plane error. Degenerate and nearly collinear geometry is rejected. Scaled positions limit avoidable intermediate overflow; normalized N/D or adjacent cross-product magnitudes below 1e-14 are treated as degenerate. Small angular arcs are ill-conditioned; use Herrick–Gibbs when times and short-arc assumptions are appropriate. This is initial orbit determination, not a propagation model or noise-optimal orbit fit.

There is no allocation, shared mutable state, or observation history stored internally. Failures preserve output. Coplanarity tolerance is in [0,pi/2] radians.

Validation: **15 frozen Vallado Python Gibbs cases**, at commit `49df0479c950917c0a2c5ac31f2db6477bba57f6`, cover circular/elliptic/parabolic/hyperbolic conics and angular spacings .02/.1/.5 radians. Maximum sampled velocity-component difference was **8.85e-8 m/s** on AppleClang 21; the static tolerance is 2e-6 m/s. Additional tests reject collinear/noncoplanar observations and check failure preservation. Raw oracle results and generator are ignored under `.validation/tinygibbs/`; normal tests use frozen arrays only.

This implementation was independently written from the Gibbs N/D/S vector construction. No AGPL source was translated. See LICENSE/NOTICE.
