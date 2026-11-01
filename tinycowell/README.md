# tinycowell

Cartesian inertial orbit integration with selected RK4 fixed or RKF45 adaptive integration. Central gravity is an explicit nonnegative model parameter; an optional caller callback supplies additional acceleration, allowing gravity/drag/SRP/third-body modules or custom forces. Set central mu=0 when the callback includes it to avoid double-counting. Input states are SI in one inertial frame, epochs are caller-defined seconds, and callbacks must share those conventions. Fixed method uses48 workspace doubles; adaptive54, plus96 bytes bounded local state scratch. No heap, I/O or mutable globals. Outputs preserved on failures; callback side effects are not reversible. No embedded force archives, maneuvers, dense output or event finding.

## Verification

32 independent SciPy DOP853 orbit solutions cover forward/backward central trajectories and linear drag callbacks. Reference rtol2.3e-14; adaptive local atol1e-9/rtol1e-12 is tested against2mm and3micrometres/s global component tolerances; fixed1s RK4 uses3mm/5micrometres/s. Module tests report maximum sampled adaptive differences, and check callback/work-limit failures. Reference: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html .

Reference scripts, packages and raw outputs are in ignored `.validation/tinycowell/`. Regression fixtures are checked in and run offline without Python or reference libraries.
