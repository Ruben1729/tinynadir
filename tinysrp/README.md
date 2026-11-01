# tinysrp

Cannonball solar-radiation-pressure acceleration, using explicit optical reflectivity, area/mass, illumination fraction, reference pressure and distance. The returned force points away from the Sun and scales inversely with distance squared. The caller supplies a Sun position and eclipse model, allowing each dependency to be replaced. SI throughout; no panel attitude effects, self-shadowing, heap, mutable globals or I/O.

## Verification

100 independent NumPy vector force cases freeze expected outputs across geometry, illumination and area/mass. Componentwise relative tolerance 5e-15, plus exact direction/umbra cases. Existing equivalent model description: https://help.agi.com/stkEngineOnUNIX/Content/hpop/hpop.htm ; numerical oracle https://numpy.org/doc/stable/reference/generated/numpy.linalg.norm.html .

Reference scripts, packages and raw outputs are in ignored `.validation/tinysrp/`. Regression fixtures are checked in and run offline without Python or reference libraries.
