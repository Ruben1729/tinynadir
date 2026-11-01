# tinyephemeris

Part of [Environment](../README.md) in [TinyNadir](../../README.md).

Bounded tabulated-ephemeris interpolation using cubic Hermite position and its derivative for velocity. The caller supplies samples in a single named frame, origin and time scale, preferably mission-relative seconds to retain precision. No automatic JPL archive loading, model downloads or extrapolation. Two to4096 samples; local state storage48 bytes. Interpolation accuracy depends on the supplied truth model and spacing. A selectable alternative to analytical Sun/Moon models or a way to consume uploaded reference trajectories. Link `tinyephemeris::tinyephemeris`.

## Verification

96 independent SciPy CubicHermiteSpline vector trajectories and derivatives are frozen. Tests check endpoints, invalid/sorted data, out-of-range requests and failure preservation. Scaled component tolerance3e-13; velocity additionally allows a1e-9m/s floor for cancellation in synthetic high-curvature position data. Oracle: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicHermiteSpline.html .

Reference scripts, packages and raw outputs are in ignored `.validation/tinyephemeris/`. Regression fixtures are checked in and run offline without Python or reference libraries.
