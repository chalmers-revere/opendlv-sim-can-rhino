/*
 * File: _coder_L_f_L_f_h_ang_Ccode_api.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 25-Jan-2019 20:24:04
 */

/* Include Files */
#include "tmwtypes.h"
#include "_coder_L_f_L_f_h_ang_Ccode_api.h"
#include "_coder_L_f_L_f_h_ang_Ccode_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131435U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "L_f_L_f_h_ang_Ccode",               /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *yp_dot,
  const char_T *identifier);
static const mxArray *emlrt_marshallOut(const real_T u[8]);

/* Function Definitions */

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T
 */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = c_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T
 */
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *yp_dot
 *                const char_T *identifier
 * Return Type  : real_T
 */
static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *yp_dot,
  const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(yp_dot), &thisId);
  emlrtDestroyArray(&yp_dot);
  return y;
}

/*
 * Arguments    : const real_T u[8]
 * Return Type  : const mxArray *
 */
static const mxArray *emlrt_marshallOut(const real_T u[8])
{
  const mxArray *y;
  const mxArray *m0;
  static const int32_T iv0[2] = { 0, 0 };

  static const int32_T iv1[2] = { 1, 8 };

  y = NULL;
  m0 = emlrtCreateNumericArray(2, iv0, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m0, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m0, *(int32_T (*)[2])&iv1[0], 2);
  emlrtAssign(&y, m0);
  return y;
}

/*
 * Arguments    : const mxArray * const prhs[23]
 *                const mxArray *plhs[1]
 * Return Type  : void
 */
void L_f_L_f_h_ang_Ccode_api(const mxArray * const prhs[23], const mxArray *
  plhs[1])
{
  real_T (*out)[8];
  real_T yp_dot;
  real_T psi_dot;
  real_T epsi;
  real_T ey;
  real_T s;
  real_T pos_ob_x;
  real_T pos_ob_y;
  real_T vel_ob_x;
  real_T vel_ob_y;
  real_T acc_ob_x;
  real_T acc_ob_y;
  real_T Ds;
  real_T a_m;
  real_T xp_dot;
  real_T a;
  real_T b;
  real_T cf;
  real_T cr;
  real_T m;
  real_T Iz;
  real_T psi_dot_com;
  real_T steer;
  real_T acc;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  out = (real_T (*)[8])mxMalloc(sizeof(real_T [8]));

  /* Marshall function inputs */
  yp_dot = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[0]), "yp_dot");
  psi_dot = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[1]),
    "psi_dot");
  epsi = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[2]), "epsi");
  ey = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[3]), "ey");
  s = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[4]), "s");
  pos_ob_x = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[5]),
    "pos_ob_x");
  pos_ob_y = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[6]),
    "pos_ob_y");
  vel_ob_x = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[7]),
    "vel_ob_x");
  vel_ob_y = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[8]),
    "vel_ob_y");
  acc_ob_x = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[9]),
    "acc_ob_x");
  acc_ob_y = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[10]),
    "acc_ob_y");
  Ds = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[11]), "Ds");
  a_m = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[12]), "a_m");
  xp_dot = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[13]),
    "xp_dot");
  a = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[14]), "a");
  b = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[15]), "b");
  cf = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[16]), "cf");
  cr = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[17]), "cr");
  m = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[18]), "m");
  Iz = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[19]), "Iz");
  psi_dot_com = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[20]),
    "psi_dot_com");
  steer = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[21]), "steer");
  acc = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[22]), "acc");

  /* Invoke the target function */
  L_f_L_f_h_ang_Ccode(yp_dot, psi_dot, epsi, ey, s, pos_ob_x, pos_ob_y, vel_ob_x,
                      vel_ob_y, acc_ob_x, acc_ob_y, Ds, a_m, xp_dot, a, b, cf,
                      cr, m, Iz, psi_dot_com, steer, acc, *out);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*out);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void L_f_L_f_h_ang_Ccode_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  L_f_L_f_h_ang_Ccode_xil_terminate();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void L_f_L_f_h_ang_Ccode_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void L_f_L_f_h_ang_Ccode_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * File trailer for _coder_L_f_L_f_h_ang_Ccode_api.c
 *
 * [EOF]
 */
