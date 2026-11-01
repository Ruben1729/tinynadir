#ifndef TINYMEKF_H
#define TINYMEKF_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYMEKF_OK=0,TINYMEKF_INVALID_ARGUMENT,TINYMEKF_INVALID_COVARIANCE,TINYMEKF_DEGENERATE,TINYMEKF_REJECTED,TINYMEKF_ARITHMETIC_ERROR} tinymekf_status;
typedef struct {double q[4],bias[3],covariance[36];} tinymekf_state;
typedef struct {double a[36],b[36],c[36],measurement[18],gain[18],innovation[9],scratch[72];} tinymekf_workspace;
/* q is scalar-first Hamilton reference->body. Error state [dtheta,dbias],
 * C_true=Exp(-[dtheta]cross)*C_estimate; covariance row-major 6x6.
 * Initialize from a normalized/nonzero attitude estimate, bias [rad/s], and
 * symmetric positive-definite covariance. No hidden initial uncertainty. */
tinymekf_status tinymekf_initialize(tinymekf_state *state,const double q[4],const double bias[3],const double covariance[36],tinymekf_workspace *workspace);
/* Held physical body gyro prediction. dt in (0,1] seconds. q uses exact held
 * rate exponential; Phi=I+F*dt, F=[-cross(omega),-I;0,0]. Noise PSD inputs are
 * per-axis gyro [rad²/s] and bias walk [(rad/s)²/s], discretized to Q*dt.
 * This explicit first-order covariance profile requires sufficiently short dt. */
tinymekf_status tinymekf_predict(tinymekf_state *state,const double gyro[3],double dt,
    const double gyro_noise_psd[3],const double bias_walk_psd[3],tinymekf_workspace *workspace);
/* Sequential normalized direction-vector update. measurement_variance >0 is
 * isotropic body-component variance for the normalized vector. H=[cross(bpred),0].
 * gate_nis=0 disables the positive NIS threshold. Corrections exceeding 0.5 rad
 * are rejected to protect the small-error-state profile. Accepted updates use
 * Joseph covariance, quaternion/bias injection, and G=I-0.5*cross(dtheta) reset.
 * Optional nis is written only on success. REJECTED and all failures preserve
 * state and nis. Workspace must be separate from inputs/output state. */
tinymekf_status tinymekf_update_vector(tinymekf_state *state,const double reference_direction[3],const double measured_body_direction[3],
    double measurement_variance,double gate_nis,double *nis,tinymekf_workspace *workspace);
#ifdef __cplusplus
}
#endif
#endif
