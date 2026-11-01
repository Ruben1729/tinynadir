# tinybdot

Magnetic B-dot detumbling in heap-free C99. Include `tinybdot/tinybdot.h`; link `tinybdot::tinybdot`. Zero-initialize a caller-owned `tinybdot_state` before first use.

The controller differences consecutive calibrated body magnetic fields (tesla), filters the derivative with backward-Euler coefficient `alpha=dt/(tau+dt)`, then commands `dipole=-gain*filtered_derivative` with per-axis A m² limits. The first sample initializes history and commands zero. Positive sample time, nonnegative gain/time constant/limits, and finite values are required. A NULL saturation mask is allowed. Failed calls preserve state/outputs; state is six doubles plus an initialization flag with ABI padding. Other scratch is bounded stack-local data; no allocation or global state occurs.

B-dot needs no orbit propagator, IGRF model, or attitude estimate. Hardware must provide valid timestamps/measurements and avoid torquer-induced magnetometer contamination. Orbital field variation affects behavior at low rates; gains, filtering, scheduling, and detumbling termination criteria require mission validation.

## Verification

128 frozen samples compare the derivative filter with **SciPy 1.18.1 signal.lfilter**, using NumPy 2.5.4 seeded vector data (949221), followed by independent clipped dipole evaluation. Maximum sampled dipole difference is 5.55e-17 A m² (3e-13 test limit); filtered derivatives match to 3e-18 T/s. Static tests also check first-sample behavior, unfiltered derivative sign, limits, invalid timing, and state/output preservation. Reference JSON/configuration/generator/logs are ignored under `.validation/tinybdot/` and `.validation/adcs-reference/`; offline C tests contain frozen numerical expectations only.

The law follows the B-dot detumbling model discussed in [NASA spacecraft-control design examples](https://ntrs.nasa.gov/api/citations/20110007876/downloads/20110007876.pdf). Reference agreement does not establish closed-loop detumbling performance. See `NOTICE` and `LICENSE`.
