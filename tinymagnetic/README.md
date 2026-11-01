# tinymagnetic

Magnetorquer allocation for a body-torque request. Include `tinymagnetic/tinymagnetic.h`; link `tinymagnetic::tinymagnetic`.

The minimum-norm unconstrained dipole is `m=(B cross torque)/|B|^2`. This profile then clips each dipole axis to its supplied nonnegative limit. It returns the command in A m² and the residual `requested - cross(m,B)` in N m. The input magnetic field is expressed in body axes in tesla. **Torque parallel to the instantaneous field cannot be generated.** The clipped solution is not advertised as an optimal constrained allocation; the residual makes lost authority visible.

Zero field is degenerate; invalid inputs/limits or overflow return status. Output buffers must be distinct and are preserved on failure. All scratch is bounded local vectors; no heap, I/O, persistent state, or hardware drivers are included. This allocator is separate from B-dot, which already requests dipole directly.

## Verification

64 frozen **NumPy 2.5.4 SVD pseudoinverse** solutions of the magnetic torque matrix, followed by independently computed component clipping and cross-product residuals (seed 331377), compare dipole to 3e-13 A m² and residual to 3e-17 N m. Maximum sampled dipole difference is 2.5e-16 A m². Tests exercise saturation, unachievable parallel torque, zero field, invalid limits, and output preservation. Raw reference data/scripts/logs remain ignored under `.validation/tinymagnetic/` and `.validation/adcs-reference/`; normal C tests run offline. See `NOTICE` and `LICENSE`.
