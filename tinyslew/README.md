# tinyslew

Shortest-geodesic rest-to-rest attitude references with triangular or trapezoidal angular-speed profiles. Include `tinyslew/tinyslew.h`; link `tinyslew::tinyslew`.

`tinyslew_initialize` normalizes initial/target scalar-first Hamilton reference-to-body quaternions, selects the shortest equivalent relative rotation, and creates a time-minimal profile for positive maximum angular rate [rad/s] and angular acceleration [rad/s²]. Short moves use a triangular speed profile; longer moves cruise at the rate limit. This profile is acceleration-limited, **not jerk-limited**, and assumes zero initial/final rates and no additional pointing/exclusion constraints.

`tinyslew_sample` receives elapsed nonnegative seconds and returns the desired quaternion and **physical angular rate/acceleration in desired body axes**. For this reference-to-body convention, the rate is negative the numeric relative quaternion axis times angular speed. After completion, sampling returns the aligned target representative and zero rate/acceleration. Acceleration jumps at phase boundaries; roundoff can choose either side at a boundary, while attitude/rate remain continuous.

The immutable caller-owned plan contains 17 doubles (136 bytes for binary64), with bounded local scratch. Use plans only after successful initialization; outputs are distinct and preserved on failure. There is no heap, I/O, mutable global state, or clock access. Add guidance, controller tuning, actuator feasibility, and avoidance constraints externally.

## Verification

32 independent analytic rest-to-rest profiles and 256 frozen **SciPy 1.18.1 Slerp** samples (NumPy 2.5.4, seed 231217) compare orientation, rates, acceleration away from discontinuities, and duration. Maximum sampled quaternion difference is 7.22e-16 modulo sign; orientation/rate tolerance is 3e-13. Tests also check rate limits, same-orientation/antipodal endpoints, physical positive-rotation rate sign, invalid limits/time, and preservation. Raw oracle/generator/logs remain ignored under `.validation/tinyslew/` and `.validation/adcs-reference/`; retained C tests run offline. See `NOTICE` and `LICENSE`.
