# tinycovariance

Part of [Mission](../README.md) in [TinyNadir](../../README.md).

Heap-free C99 first-order uncertainty mapping with dimensions 1..6. Include `tinycovariance/tinycovariance.h`, link `tinycovariance::tinycovariance` (depends on `tinyelements`).

`tinycovariance_transform` computes the congruence `Pout=J*P*J^T` for a supplied Jacobian. `tinycovariance_jacobian` computes bounded central numerical derivatives of a caller model, with explicit physical input steps and optional output periods for angular unwrapping. A callback must fill all output coordinates and return success; it is called at most twice per input dimension. `tinycovariance_rotate_state` handles a fixed known orthogonal frame matrix applied to position and velocity, not a differentiated orbit-dependent frame or angular-rate conversion.

Convenience functions map Cartesian covariance to/from classical `[p,e,i,Omega,omega,nu]` and prograde modified equinoctial `[p,f,g,h,k,L]` elements. **p is semilatus rectum and nu is true anomaly**, unlike interfaces using semimajor axis/mean anomaly. These wrappers use central differences at a cube-root-epsilon scale and explicit angle unwrapping. Classical maps reject e<1e-8 and abs(sin(i))<1e-8; equinoctial maps retain their i=pi singularity. All units match `tinyelements` (SI/radians); covariance entries carry the corresponding product units.

Input covariance must be finite symmetric PSD. The PSD check is performed on the correlation-scaled matrix to avoid discarding small variances merely because coordinate units differ. Numerically indefinite/ill-conditioned transforms can be rejected. Failures preserve outputs. No heap or mutable global state is used. This is local linear uncertainty transformation, not probability integration, an orbit filter, or a guarantee that large nonlinear uncertainties remain Gaussian.

Validation uses **NumPy 2.5.4 congruence/rotation products and six SciPy 1.18.1 numerical Jacobians of pinned Vallado orbit conversions**. Maximum orbital covariance difference normalized by the reference diagonal standard-deviation products was **1.20e-8** on AppleClang 21; the static limit is 3e-6. Inverse covariance round trips, a nonlinear angular Jacobian, and invalid PSD/preservation cases are also checked. Raw references/generators stay ignored in `.validation/tinycovariance/`; normal tests require no Python or reference library.

The C implementation was independently written from covariance and conic equations; no AGPL conversion/Jacobian code was translated. See LICENSE/NOTICE.
