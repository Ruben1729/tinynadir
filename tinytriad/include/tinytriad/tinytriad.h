#ifndef TINYTRIAD_H
#define TINYTRIAD_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYTRIAD_OK=0,TINYTRIAD_INVALID_ARGUMENT,TINYTRIAD_DEGENERATE,TINYTRIAD_ARITHMETIC_ERROR} tinytriad_status;
/* References and corresponding body observations are direction vectors.
 * They need not be normalized. Both vector pairs must have |sin(angle)|>1e-10.
 * First-vector alignment has priority; second-vector alignment supplies roll.
 * Output DCM is row-major, body=C*reference. Output quaternion is scalar-first
 * Hamilton, in the tinyrotation convention. Outputs are preserved on failure. */
tinytriad_status tinytriad_dcm(const double reference1[3],const double reference2[3],
    const double body1[3],const double body2[3],double dcm[9]);
tinytriad_status tinytriad_solve(const double reference1[3],const double reference2[3],
    const double body1[3],const double body2[3],double quaternion[4]);
#ifdef __cplusplus
}
#endif
#endif
