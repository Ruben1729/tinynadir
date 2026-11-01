# tinyspaceweather

Part of [Environment](../README.md) in [TinyNadir](../../README.md).

Bounded caller-table solar/geomagnetic activity interpolation: F10.7 and its supplied average in solar flux units (1sfu=1e-22W/m2/Hz), and Ap. Time coordinates use an explicit application scale; no network/file readers, defaults, leap handling, or inferred averaging windows. Linear interpolation is the selected profile, not Vallado cubic interpolation or forecast generation. Kp↔Ap uses the standard28 nodes with linear interpolation between thirds. The application must supply any daily/3-hour history needed by an atmosphere model. No heap, mutable globals, workspace or large data archive. Failure preserves output.

## Verification

64 independent NumPy interpolation cases and109 Kp/Ap conversions from standard tabulated nodes are frozen, scaled tolerance3e-14. Includes endpoints, out-of-range/NaN rejection and failure atomicity. Table source: https://kp.gfz.de/en/data ; oracle https://numpy.org/doc/stable/reference/generated/numpy.interp.html .

Reference scripts, packages and raw outputs are in ignored `.validation/tinyspaceweather/`. Regression fixtures are checked in and run offline without Python or reference libraries.
