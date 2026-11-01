# tinyexponential

Piecewise exponential atmospheric density using caller-supplied layers and an explicit altitude validity ceiling. Inputs are metres, kg/m3, and positive scale heights. No default Earth atmosphere is implied, no implicit extrapolation, and discontinuities between supplied layers are retained. Choose a calibrated profile for the mission. No workspace or heap; only loop scalars.

## Verification

96 layered profiles/altitudes are compared with independent NumPy exponential/table selection and frozen as offline cases. Relative tolerance 2e-14; includes each layer boundary and invalid/out-of-range failures. Mathematical reference: exponential scale-height atmosphere, numerical oracle https://numpy.org/doc/stable/reference/generated/numpy.exp.html .

Reference scripts, packages and raw outputs are in ignored `.validation/tinyexponential/`. Regression fixtures are checked in and run offline without Python or reference libraries.
