# tinybatchls

Part of [Mission](../README.md) in [TinyNadir](../../README.md).

Heap-free C99 bounded weighted nonlinear least squares using central numerical Jacobians and Levenberg–Marquardt damping. Include `tinybatchls/tinybatchls.h`; link `tinybatchls::tinybatchls`. A caller prediction callback can compose an orbit propagator, observation model, calibration model, or another calculation.

The profile supports **1..8 fitted parameters, at most 64 scalar observations, and at most 100 outer iterations**. Measurement standard deviations are explicit positive scalars (diagonal weighting); whiten correlated data externally. Positive parameter scales define the physical units used by the normalized solve. Numerical-difference steps and convergence tolerances are configurable. No bounds, priors, robust loss, or guaranteed global minimum is implicit. The callback must write every prediction and return nonzero on success.

The caller supplies `tinybatchls_workspace` (**7360 bytes with eight-byte double**) and owns input/output storage. Scratch is modified even on failure; a failed call preserves the result object. At most `1 + (2*n+1)*max_iterations + 2*n` model evaluations are performed. The callback's own memory/work remains the application's responsibility. A pivot threshold detects unresolved rank deficiency; poor scaling/weak observability should be addressed rather than treating the covariance as reliable.

Output covariance is the inverse weighted normal matrix transformed back to physical parameter units, using **known supplied measurement variances**, without reduced-chi-square rescaling. This is local linearized uncertainty at the fitted state.

Validation: **three independent SciPy 1.18.1 LM fits** (linear, exponential decay, logistic) with weighted noisy observations and NumPy covariance from analytical Jacobians. Maximum sampled parameter difference was **4.87e-10**, and covariance difference normalized by diagonal standard-deviation products **2.58e-10**, on AppleClang 21. Static regression limits are 2e-7 for parameter components and 3e-6 normalized covariance; weighted costs, work bounds, rank deficiency, callback errors, and result preservation are checked. Raw comparisons/generator are ignored in `.validation/tinybatchls/`; normal tests are offline C99/libm.

Independent implementation from the weighted Gauss–Newton/Levenberg–Marquardt equations. No SciPy source was translated or linked. See LICENSE/NOTICE.
