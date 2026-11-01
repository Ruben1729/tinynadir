# tinymekf

A genuine six-error-state multiplicative extended Kalman attitude filter in heap-free C99. It predicts covariance, forms vector innovations/gains, performs Joseph covariance updates, injects quaternion/bias corrections, and resets the local attitude-error covariance. It is not an alias for a complementary filter.

Include `tinymekf/tinymekf.h`; link `tinymekf::tinymekf`. The nominal attitude is scalar-first Hamilton reference-to-body, plus physical body gyro bias [rad/s]. The six-dimensional error is `[delta_theta,delta_bias]` with `C_true=Exp(-cross(delta_theta))*C_estimate`. The caller supplies initial quaternion, bias, and a positive-definite 6×6 covariance; initialization does not invent an uncertainty model.

## Model profile

- Attitude prediction uses the exact held corrected-body-rate exponential from `tinygyro`.
- Covariance uses `Phi=I+F*dt`, with `F=[-cross(omega),-I;0,0]`, and first-order diagonal process covariance `diag(gyro_noise_PSD,bias_walk_PSD)*dt`. Gyro PSD units are rad²/s; bias-walk PSD units are (rad/s)²/s. Higher-order noise-integration/coupling terms are omitted in this profile. Steps must be in (0,1] seconds and short enough for the selected angular rates/noise model.
- Each update consumes one nonzero modeled reference direction and its normalized body observation. The measurement Jacobian is `[cross(predicted_body_direction),0]`. The positive scalar measurement variance is isotropic body-component variance of the normalized vector; it is an approximation of a direction sensor's tangent-plane noise, not a full calibration/correlation model.
- A positive optional NIS gate rejects inconsistent data; zero disables that gate. A proposed physical attitude correction larger than 0.5 rad is also rejected, because this is a small-error filter. Initialize from an appropriate vector/star-tracker attitude first.
- Accepted updates use Joseph covariance, physical quaternion correction `Exp(-delta_theta)`, bias injection, and covariance reset `G=diag(I-0.5*cross(delta_theta),I)`.

Multiple independent vectors can be updated sequentially after each prediction. Do not feed a QUEST solution and its originating raw vectors as independent observations. This profile does not implement a full UKF, arbitrary state dimensions, a full quaternion-measurement update, noise identification, data synchronization, or mode supervision.

State is **344 bytes** and workspace **1,800 bytes** for binary64 doubles; bounded local vectors, call frames, and external math routines add stack. Workspace must not overlap the state/inputs. All state/output diagnostics are preserved on rejection/failure; NIS is written only for an accepted update. Covariances are symmetrized and checked positive-definite through Cholesky. No heap, I/O, shared mutable globals, or external reference runtime exists.

## Verification

The frozen sequence contains **64 FilterPy 2.0.0 covariance predictions and 128 ExtendedKalmanFilter Joseph updates**, with **SciPy 1.18.1 Rotation** handling independent nominal attitude propagation/injection and NumPy 2.5.4 handling reset transforms (seed 792212). It checks predicted and corrected quaternion, bias, full covariance, and NIS at each stage. Maximum sampled quaternion difference is **3.82e-16** modulo sign; maximum covariance-component difference is **1.08e-19**. Limits are 3e-11 for quaternion, 3e-12 for bias, and `5e-17+3e-11*abs(expected)` for covariance.

Tests separately check covariance positive definiteness, single-vector unobservable-axis covariance, zero-rate analytic bias/attitude covariance coupling, physical convergence to an independently propagated truth attitude, NIS rejection, invalid time/covariance, and atomic state/diagnostic preservation. Comparisons verify this stated discretization, not all MEKF variants or flight estimation performance. Raw FilterPy/SciPy oracle scripts, full JSON, package environment, and logs remain ignored under `.validation/tinymekf/` and `.validation/adcs-reference/`; normal C tests use static fixtures only.

References: [Markley, Multiplicative vs. Additive Filtering](https://ntrs.nasa.gov/api/citations/20040037784/downloads/20040037784.pdf), [FilterPy ExtendedKalmanFilter](https://filterpy.readthedocs.io/en/latest/kalman/ExtendedKalmanFilter.html). Independently written equations and code; see `NOTICE` and `LICENSE`.
