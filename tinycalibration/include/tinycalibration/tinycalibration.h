#ifndef TINYCALIBRATION_H
#define TINYCALIBRATION_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYCALIBRATION_OK=0,TINYCALIBRATION_INVALID_ARGUMENT,TINYCALIBRATION_INVALID_ALIGNMENT,TINYCALIBRATION_ARITHMETIC_ERROR} tinycalibration_status;
typedef struct {double bias[3],correction[9],sensor_to_body[9];} tinycalibration_parameters;
/* body=alignment*correction*(raw-bias), with row-major matrices. correction
 * may encode scale/cross-axis/soft-iron corrections; alignment must be a proper
 * rotation. Units follow the supplied raw values and calibration parameters.
 * This applies supplied calibration, it does not fit parameters or own time.
 * Output may alias raw; failures preserve output. */
tinycalibration_status tinycalibration_apply(const tinycalibration_parameters *parameters,const double raw[3],double body[3]);
#ifdef __cplusplus
}
#endif
#endif
