/*
 * File: _coder_L_f_L_f_h_ang_Ccode_api.h
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 25-Jan-2019 20:24:04
 */

#ifndef _CODER_L_F_L_F_H_ANG_CCODE_API_H
#define _CODER_L_F_L_F_H_ANG_CCODE_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_L_f_L_f_h_ang_Ccode_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void L_f_L_f_h_ang_Ccode(real_T yp_dot, real_T psi_dot, real_T epsi,
  real_T ey, real_T s, real_T pos_ob_x, real_T pos_ob_y, real_T vel_ob_x, real_T
  vel_ob_y, real_T acc_ob_x, real_T acc_ob_y, real_T Ds, real_T a_m, real_T
  xp_dot, real_T a, real_T b, real_T cf, real_T cr, real_T m, real_T Iz, real_T
  psi_dot_com, real_T steer, real_T acc, real_T out[8]);
extern void L_f_L_f_h_ang_Ccode_api(const mxArray * const prhs[23], const
  mxArray *plhs[1]);
extern void L_f_L_f_h_ang_Ccode_atexit(void);
extern void L_f_L_f_h_ang_Ccode_initialize(void);
extern void L_f_L_f_h_ang_Ccode_terminate(void);
extern void L_f_L_f_h_ang_Ccode_xil_terminate(void);

#endif

/*
 * File trailer for _coder_L_f_L_f_h_ang_Ccode_api.h
 *
 * [EOF]
 */
