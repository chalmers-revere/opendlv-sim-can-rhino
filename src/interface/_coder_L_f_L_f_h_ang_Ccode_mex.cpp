/*
 * File: _coder_L_f_L_f_h_ang_Ccode_mex.cpp
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 25-Jan-2019 20:24:04
 */

/* Include Files */
#include "_coder_L_f_L_f_h_ang_Ccode_api.h"
#include "_coder_L_f_L_f_h_ang_Ccode_mex.h"

/* Function Declarations */
static void L_f_L_f_h_ang_Ccode_mexFunction(int32_T nlhs, mxArray *plhs[1],
  int32_T nrhs, const mxArray *prhs[23]);

/* Function Definitions */

/*
 * Arguments    : int32_T nlhs
 *                const mxArray *plhs[1]
 *                int32_T nrhs
 *                const mxArray *prhs[23]
 * Return Type  : void
 */
static void L_f_L_f_h_ang_Ccode_mexFunction(int32_T nlhs, mxArray *plhs[1],
  int32_T nrhs, const mxArray *prhs[23])
{
  int32_T n;
  const mxArray *inputs[23];
  const mxArray *outputs[1];
  int32_T b_nlhs;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 23) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 23, 4,
                        19, "L_f_L_f_h_ang_Ccode");
  }

  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 19,
                        "L_f_L_f_h_ang_Ccode");
  }

  /* Temporary copy for mex inputs. */
  for (n = 0; n < nrhs; n++) {
    inputs[n] = prhs[n];
  }

  /* Call the function. */
  L_f_L_f_h_ang_Ccode_api(inputs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);

  /* Module termination. */
  L_f_L_f_h_ang_Ccode_terminate();
}

/*
 * Arguments    : int32_T nlhs
 *                const mxArray * const plhs[]
 *                int32_T nrhs
 *                const mxArray * const prhs[]
 * Return Type  : void
 */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(L_f_L_f_h_ang_Ccode_atexit);

  /* Initialize the memory manager. */
  /* Module initialization. */
  L_f_L_f_h_ang_Ccode_initialize();

  /* Dispatch the entry-point. */
  L_f_L_f_h_ang_Ccode_mexFunction(nlhs, plhs, nrhs, prhs);
}

/*
 * Arguments    : void
 * Return Type  : emlrtCTX
 */
emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/*
 * File trailer for _coder_L_f_L_f_h_ang_Ccode_mex.cpp
 *
 * [EOF]
 */
