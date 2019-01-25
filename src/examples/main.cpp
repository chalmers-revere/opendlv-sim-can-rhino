//
// File: main.cpp
//
// MATLAB Coder version            : 3.2
// C/C++ source code generated on  : 25-Jan-2019 20:24:04
//

//***********************************************************************
// This automatically generated example C main file shows how to call
// entry-point functions that MATLAB Coder generated. You must customize
// this file for your application. Do not modify this file directly.
// Instead, make a copy of this file, modify it, and integrate it into
// your development environment.
//
// This file initializes entry-point function arguments to a default
// size and value before calling the entry-point functions. It does
// not store or use any values returned from the entry-point functions.
// If necessary, it does pre-allocate memory for returned values.
// You can use this file as a starting point for a main function that
// you can deploy in your application.
//
// After you copy the file, and before you deploy it, you must make the
// following changes:
// * For variable-size function arguments, change the example sizes to
// the sizes that your application requires.
// * Change the example values of function arguments to the values that
// your application requires.
// * If the entry-point functions return values, store these values or
// otherwise use them as required by your application.
//
//***********************************************************************
// Include Files
#include "rt_nonfinite.h"
#include "L_f_L_f_h_ang_Ccode.h"
#include "main.h"
#include "L_f_L_f_h_ang_Ccode_terminate.h"
#include "L_f_L_f_h_ang_Ccode_initialize.h"

// Function Declarations
static double argInit_real_T();
static void main_L_f_L_f_h_ang_Ccode();

// Function Definitions

//
// Arguments    : void
// Return Type  : double
//
static double argInit_real_T()
{
  return 0.0;
}

//
// Arguments    : void
// Return Type  : void
//
static void main_L_f_L_f_h_ang_Ccode()
{
  double yp_dot;
  double psi_dot;
  double epsi;
  double ey;
  double s;
  double pos_ob_x;
  double pos_ob_y;
  double vel_ob_x;
  double vel_ob_y;
  double acc_ob_x;
  double acc_ob_y;
  double Ds;
  double a_m;
  double out[8];

  // Initialize function 'L_f_L_f_h_ang_Ccode' input arguments.
  yp_dot = argInit_real_T();
  psi_dot = argInit_real_T();
  epsi = argInit_real_T();
  ey = argInit_real_T();
  s = argInit_real_T();
  pos_ob_x = argInit_real_T();
  pos_ob_y = argInit_real_T();
  vel_ob_x = argInit_real_T();
  vel_ob_y = argInit_real_T();
  acc_ob_x = argInit_real_T();
  acc_ob_y = argInit_real_T();
  Ds = argInit_real_T();
  a_m = argInit_real_T();

  // Call the entry-point 'L_f_L_f_h_ang_Ccode'.
  L_f_L_f_h_ang_Ccode(yp_dot, psi_dot, epsi, ey, s, pos_ob_x, pos_ob_y, vel_ob_x,
                      vel_ob_y, acc_ob_x, acc_ob_y, Ds, a_m, argInit_real_T(),
                      argInit_real_T(), argInit_real_T(), argInit_real_T(),
                      argInit_real_T(), argInit_real_T(), argInit_real_T(),
                      argInit_real_T(), argInit_real_T(), argInit_real_T(), out);
}

//
// Arguments    : int argc
//                const char * const argv[]
// Return Type  : int
//
int main(int, const char * const [])
{
  // Initialize the application.
  // You do not need to do this more than one time.
  L_f_L_f_h_ang_Ccode_initialize();

  // Invoke the entry-point functions.
  // You can call entry-point functions multiple times.
  main_L_f_L_f_h_ang_Ccode();

  // Terminate the application.
  // You do not need to do this more than one time.
  L_f_L_f_h_ang_Ccode_terminate();
  return 0;
}

//
// File trailer for main.cpp
//
// [EOF]
//
