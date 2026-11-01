# tinyunscented

Heap-free C99 scaled unscented transform for Euclidean vectors of dimension 1..6, up to 13 sigma points. Include `tinyunscented/tinyunscented.h`, link `tinyunscented::tinyunscented`. It is a moment transform, **not a complete UKF**.

`tinyunscented_generate` accepts mean, symmetric PSD covariance, and Merwe parameters alpha>0, beta>=0, n+kappa>0. Points are ordered center, positive columns of the scaled lower Cholesky factor, negative columns. `tinyunscented_statistics` reconstructs signed-weight mean/covariance. `tinyunscented_transform` evaluates a caller nonlinear mapping at each point and computes output moments, at most 13 calls. Arrays are packed row-major and dimensions/counts are explicit.

Covariance factorization uses correlation scaling to preserve quantities expressed in very different units. A nonzero sigma offset that cannot be represented relative to its mean is rejected rather than silently collapsing that uncertainty. Weighted sums use compensation. Arbitrary signed weights need not yield PSD covariance; angle/quaternion/manifold means and process/measurement noise are not implicit. Select scaling/coordinates appropriate to the problem.

The caller owns sigma points and outputs; all internal arrays are bounded. No heap or shared mutable state is used, and failure preserves outputs. Mathematical nonfinite/PSD/model errors return explicit statuses.

Validation compares **nine FilterPy 2.0.0 `MerweScaledSigmaPoints`/`unscented_transform` cases**, spanning dimensions 2/3/6 and alpha .2/.7/1. Maximum sampled transformed mean/covariance component difference was **5.51e-14** on AppleClang 21. Static tests also reconstruct the original moments, exercise 1e20/1e-20 mixed variances, and verify invalid covariance/preservation behavior. Raw oracle outputs and generator are ignored in `.validation/tinyunscented/`; normal tests are offline.

Independent implementation from scaled-unscented equations; no FilterPy code was translated. FilterPy is used only to capture expected numbers. See LICENSE/NOTICE.
