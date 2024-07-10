//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: thesis_mpc.h
//
// Code generated for Simulink model 'thesis_mpc'.
//
// Model version                  : 3.0
// Simulink Coder version         : 9.9 (R2023a) 19-Nov-2022
// C/C++ source code generated on : Tue May 16 17:18:32 2023
//
// Target selection: ert.tlc
// Embedded hardware selection: ARM Compatible->ARM Cortex-M
// Code generation objectives:
//    1. Execution efficiency
//    2. RAM efficiency
// Validation result: Not run
//
#ifndef RTW_HEADER_thesis_mpc_h_
#define RTW_HEADER_thesis_mpc_h_
#include "rtwtypes.h"
#include <stddef.h>
#include <cstdlib>
#include <time.h>

// Macros for accessing real-time model data structure
#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_WTmPWsEMvOzNnnAVv5fQNC_
#define DEFINED_TYPEDEF_FOR_struct_WTmPWsEMvOzNnnAVv5fQNC_

struct struct_WTmPWsEMvOzNnnAVv5fQNC
{
  int32_T MaxIterations;
  real_T ConstraintTolerance;
  boolean_T UseWarmStart;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_WHjMt45Sk148iktWsfFxl_
#define DEFINED_TYPEDEF_FOR_struct_WHjMt45Sk148iktWsfFxl_

struct struct_WHjMt45Sk148iktWsfFxl
{
  int32_T MaxIterations;
  real_T ConstraintTolerance;
  real_T OptimalityTolerance;
  real_T ComplementarityTolerance;
  real_T StepTolerance;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_lnQ9KXdSZFplhcBp5LBCc_
#define DEFINED_TYPEDEF_FOR_struct_lnQ9KXdSZFplhcBp5LBCc_

struct struct_lnQ9KXdSZFplhcBp5LBCc
{
  int32_T MaxIterations;
  real_T ConstraintTolerance;
  real_T DiscreteConstraintTolerance;
  boolean_T RoundingAtRootNode;
  int32_T MaxPendingNodes;
};

#endif

extern "C"
{
  static real_T rtGetNaN(void);
  static real32_T rtGetNaNF(void);
}                                      // extern "C"

#define NOT_USING_NONFINITE_LITERALS   1

extern "C"
{
  extern real_T rtInf;
  extern real_T rtMinusInf;
  extern real_T rtNaN;
  extern real32_T rtInfF;
  extern real32_T rtMinusInfF;
  extern real32_T rtNaNF;
  static void rt_InitInfAndNaN(size_t realSize);
  static boolean_T rtIsInf(real_T value);
  static boolean_T rtIsInfF(real32_T value);
  static boolean_T rtIsNaN(real_T value);
  static boolean_T rtIsNaNF(real32_T value);
  struct BigEndianIEEEDouble {
    struct {
      uint32_T wordH;
      uint32_T wordL;
    } words;
  };

  struct LittleEndianIEEEDouble {
    struct {
      uint32_T wordL;
      uint32_T wordH;
    } words;
  };

  struct IEEESingle {
    union {
      real32_T wordLreal;
      uint32_T wordLuint;
    } wordL;
  };
}                                      // extern "C"

extern "C"
{
  static real_T rtGetInf(void);
  static real32_T rtGetInfF(void);
  static real_T rtGetMinusInf(void);
  static real32_T rtGetMinusInfF(void);
}                                      // extern "C"

// Class declaration for model thesis_mpc
class thesis_mpc final
{

  // public data and function members
 public:

 int Ramdon_num = time(0);
  // Block signals and states (default storage) for system '<Root>'
  struct DW {
    real_T last_x_PreviousInput[7];    // '<S2>/last_x'
    real_T last_mv_DSTATE;             // '<S2>/last_mv'
    boolean_T Memory_PreviousInput[68];// '<S2>/Memory'
  };

  // External inputs (root inport signals with default storage)
  struct ExtU {
    real_T uav_r;                      // '<Root>/In1'
    real_T r_cmd;                      // '<Root>/In2'
  };

  // External outputs (root outports fed by signals with default storage)
  struct ExtY {
    real_T Out1;                       // '<Root>/Out1'
  };

  // Real-time Model Data Structure
  struct RT_MODEL {
    const char_T * volatile errorStatus;
  };

  // Copy Constructor
  thesis_mpc(thesis_mpc const&) = delete;

  // Assignment Operator
  thesis_mpc& operator= (thesis_mpc const&) & = delete;

  // Move Constructor
  thesis_mpc(thesis_mpc &&) = delete;

  // Move Assignment Operator
  thesis_mpc& operator= (thesis_mpc &&) = delete;

  // Real-Time Model get method
  thesis_mpc::RT_MODEL * getRTM();

  // External inputs
  ExtU rtU;

  // External outputs
  ExtY rtY;

  // model initialize function
  void initialize();

  // model step function
  void step();

  // Constructor
  thesis_mpc();

  // Destructor
  ~thesis_mpc();

  // private data and function members
 private:
  // Block states
  DW rtDW;

  // private member function(s) for subsystem '<Root>'
  real_T norm(const real_T x[3]);
  real_T maximum(const real_T x[3]);
  real_T xnrm2(int32_T n, const real_T x[9], int32_T ix0);
  void xgemv(int32_T b_m, int32_T n, const real_T b_A[9], int32_T ia0, const
             real_T x[9], int32_T ix0, real_T y[3]);
  void xgerc(int32_T b_m, int32_T n, real_T alpha1, int32_T ix0, const real_T y
             [3], real_T b_A[9], int32_T ia0);
  void KWIKfactor(const real_T b_Ac[204], const int32_T iC[68], int32_T nA,
                  const real_T b_Linv[9], real_T D[9], real_T b_H[9], int32_T n,
                  real_T RLinv[9], real_T *Status);
  void DropConstraint(int32_T kDrop, boolean_T iA[68], int32_T *nA, int32_T iC
                      [68]);
  void qpkwik(const real_T b_Linv[9], const real_T b_Hinv[9], const real_T f[3],
              const real_T b_Ac[204], const real_T b[68], boolean_T iA[68],
              int32_T maxiter, real_T FeasTol, real_T x[3], real_T lambda[68],
              int32_T *status);

  // Real-Time Model
  RT_MODEL rtM;
};

//-
//  These blocks were eliminated from the model due to optimizations:
//
//  Block '<S2>/Constant' : Unused code path elimination
//  Block '<S2>/Floor' : Unused code path elimination
//  Block '<S2>/Floor1' : Unused code path elimination
//  Block '<S3>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S4>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S5>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S6>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S7>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S8>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S9>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S10>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S11>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S12>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S13>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S14>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S15>/Vector Dimension Check' : Unused code path elimination
//  Block '<S16>/Vector Dimension Check' : Unused code path elimination
//  Block '<S17>/Vector Dimension Check' : Unused code path elimination
//  Block '<S18>/Vector Dimension Check' : Unused code path elimination
//  Block '<S19>/Vector Dimension Check' : Unused code path elimination
//  Block '<S20>/Vector Dimension Check' : Unused code path elimination
//  Block '<S2>/Min' : Unused code path elimination
//  Block '<S2>/constant' : Unused code path elimination
//  Block '<S21>/Vector Dimension Check' : Unused code path elimination
//  Block '<S2>/umin_scale2' : Unused code path elimination
//  Block '<S2>/umin_scale3' : Unused code path elimination
//  Block '<S2>/umin_scale5' : Unused code path elimination
//  Block '<S2>/ym_zero' : Unused code path elimination
//  Block '<S1>/m_zero' : Unused code path elimination
//  Block '<S1>/p_zero' : Unused code path elimination
//  Block '<S2>/Reshape' : Reshape block reduction
//  Block '<S2>/Reshape1' : Reshape block reduction
//  Block '<S2>/Reshape2' : Reshape block reduction
//  Block '<S2>/Reshape3' : Reshape block reduction
//  Block '<S2>/Reshape4' : Reshape block reduction
//  Block '<S2>/Reshape5' : Reshape block reduction
//  Block '<S2>/ext.mv_scale' : Eliminated nontunable gain of 1
//  Block '<S2>/ext.mv_scale1' : Eliminated nontunable gain of 1
//  Block '<S2>/umin_scale1' : Eliminated nontunable gain of 1
//  Block '<S2>/umin_scale4' : Eliminated nontunable gain of 1
//  Block '<S2>/ymin_scale1' : Eliminated nontunable gain of 1
//  Block '<S2>/ymin_scale2' : Eliminated nontunable gain of 1


//-
//  The generated code includes comments that allow you to trace directly
//  back to the appropriate location in the model.  The basic format
//  is <system>/block_name, where system is the system number (uniquely
//  assigned by Simulink) and block_name is the name of the block.
//
//  Use the MATLAB hilite_system command to trace the generated code back
//  to the model.  For example,
//
//  hilite_system('<S3>')    - opens system 3
//  hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
//
//  Here is the system hierarchy for this model
//
//  '<Root>' : 'thesis_mpc'
//  '<S1>'   : 'thesis_mpc/MPC Controller'
//  '<S2>'   : 'thesis_mpc/MPC Controller/MPC'
//  '<S3>'   : 'thesis_mpc/MPC Controller/MPC/MPC Matrix Signal Check'
//  '<S4>'   : 'thesis_mpc/MPC Controller/MPC/MPC Matrix Signal Check1'
//  '<S5>'   : 'thesis_mpc/MPC Controller/MPC/MPC Matrix Signal Check2'
//  '<S6>'   : 'thesis_mpc/MPC Controller/MPC/MPC Preview Signal Check'
//  '<S7>'   : 'thesis_mpc/MPC Controller/MPC/MPC Preview Signal Check1'
//  '<S8>'   : 'thesis_mpc/MPC Controller/MPC/MPC Preview Signal Check2'
//  '<S9>'   : 'thesis_mpc/MPC Controller/MPC/MPC Preview Signal Check3'
//  '<S10>'  : 'thesis_mpc/MPC Controller/MPC/MPC Preview Signal Check4'
//  '<S11>'  : 'thesis_mpc/MPC Controller/MPC/MPC Preview Signal Check5'
//  '<S12>'  : 'thesis_mpc/MPC Controller/MPC/MPC Preview Signal Check6'
//  '<S13>'  : 'thesis_mpc/MPC Controller/MPC/MPC Preview Signal Check7'
//  '<S14>'  : 'thesis_mpc/MPC Controller/MPC/MPC Preview Signal Check8'
//  '<S15>'  : 'thesis_mpc/MPC Controller/MPC/MPC Scalar Signal Check'
//  '<S16>'  : 'thesis_mpc/MPC Controller/MPC/MPC Scalar Signal Check1'
//  '<S17>'  : 'thesis_mpc/MPC Controller/MPC/MPC Scalar Signal Check2'
//  '<S18>'  : 'thesis_mpc/MPC Controller/MPC/MPC Vector Signal Check'
//  '<S19>'  : 'thesis_mpc/MPC Controller/MPC/MPC Vector Signal Check1'
//  '<S20>'  : 'thesis_mpc/MPC Controller/MPC/MPC Vector Signal Check6'
//  '<S21>'  : 'thesis_mpc/MPC Controller/MPC/moorx'
//  '<S22>'  : 'thesis_mpc/MPC Controller/MPC/optimizer'
//  '<S23>'  : 'thesis_mpc/MPC Controller/MPC/optimizer/optimizer'

#endif                                 // RTW_HEADER_thesis_mpc_h_

//
// File trailer for generated code.
//
// [EOF]
//
