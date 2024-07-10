//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: thesis_mpc.cpp
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
#include "../include/shao_thesis_mpc/thesis_mpc.h"
#include "../include/shao_thesis_mpc/rtwtypes.h"
#include <cmath>
#include <cstring>
#include <stddef.h>

// Named constants for MATLAB Function: '<S22>/optimizer'
const int32_T degrees{ 3 };

#define NumBitsPerChar                 8U

extern real_T rt_hypotd_snf(real_T u0, real_T u1);
static int32_T div_nde_s32_floor(int32_T numerator, int32_T denominator);
extern "C"
{
  real_T rtInf;
  real_T rtMinusInf;
  real_T rtNaN;
  real32_T rtInfF;
  real32_T rtMinusInfF;
  real32_T rtNaNF;
}

extern "C"
{
  //
  // Initialize rtNaN needed by the generated code.
  // NaN is initialized as non-signaling. Assumes IEEE.
  //
  static real_T rtGetNaN(void)
  {
    size_t bitsPerReal{ sizeof(real_T) * (NumBitsPerChar) };

    real_T nan{ 0.0 };

    if (bitsPerReal == 32U) {
      nan = rtGetNaNF();
    } else {
      union {
        LittleEndianIEEEDouble bitVal;
        real_T fltVal;
      } tmpVal;

      tmpVal.bitVal.words.wordH = 0xFFF80000U;
      tmpVal.bitVal.words.wordL = 0x00000000U;
      nan = tmpVal.fltVal;
    }

    return nan;
  }

  //
  // Initialize rtNaNF needed by the generated code.
  // NaN is initialized as non-signaling. Assumes IEEE.
  //
  static real32_T rtGetNaNF(void)
  {
    IEEESingle nanF{ { 0.0F } };

    nanF.wordL.wordLuint = 0xFFC00000U;
    return nanF.wordL.wordLreal;
  }
}

extern "C"
{
  //
  // Initialize the rtInf, rtMinusInf, and rtNaN needed by the
  // generated code. NaN is initialized as non-signaling. Assumes IEEE.
  //
  static void rt_InitInfAndNaN(size_t realSize)
  {
    (void) (realSize);
    rtNaN = rtGetNaN();
    rtNaNF = rtGetNaNF();
    rtInf = rtGetInf();
    rtInfF = rtGetInfF();
    rtMinusInf = rtGetMinusInf();
    rtMinusInfF = rtGetMinusInfF();
  }

  // Test if value is infinite
  static boolean_T rtIsInf(real_T value)
  {
    return (boolean_T)((value==rtInf || value==rtMinusInf) ? 1U : 0U);
  }

  // Test if single-precision value is infinite
  static boolean_T rtIsInfF(real32_T value)
  {
    return (boolean_T)(((value)==rtInfF || (value)==rtMinusInfF) ? 1U : 0U);
  }

  // Test if value is not a number
  static boolean_T rtIsNaN(real_T value)
  {
    boolean_T result{ (boolean_T) 0 };

    size_t bitsPerReal{ sizeof(real_T) * (NumBitsPerChar) };

    if (bitsPerReal == 32U) {
      result = rtIsNaNF((real32_T)value);
    } else {
      union {
        LittleEndianIEEEDouble bitVal;
        real_T fltVal;
      } tmpVal;

      tmpVal.fltVal = value;
      result = (boolean_T)((tmpVal.bitVal.words.wordH & 0x7FF00000) ==
                           0x7FF00000 &&
                           ( (tmpVal.bitVal.words.wordH & 0x000FFFFF) != 0 ||
                            (tmpVal.bitVal.words.wordL != 0) ));
    }

    return result;
  }

  // Test if single-precision value is not a number
  static boolean_T rtIsNaNF(real32_T value)
  {
    IEEESingle tmp;
    tmp.wordL.wordLreal = value;
    return (boolean_T)( (tmp.wordL.wordLuint & 0x7F800000) == 0x7F800000 &&
                       (tmp.wordL.wordLuint & 0x007FFFFF) != 0 );
  }
}

extern "C"
{
  //
  // Initialize rtInf needed by the generated code.
  // Inf is initialized as non-signaling. Assumes IEEE.
  //
  static real_T rtGetInf(void)
  {
    size_t bitsPerReal{ sizeof(real_T) * (NumBitsPerChar) };

    real_T inf{ 0.0 };

    if (bitsPerReal == 32U) {
      inf = rtGetInfF();
    } else {
      union {
        LittleEndianIEEEDouble bitVal;
        real_T fltVal;
      } tmpVal;

      tmpVal.bitVal.words.wordH = 0x7FF00000U;
      tmpVal.bitVal.words.wordL = 0x00000000U;
      inf = tmpVal.fltVal;
    }

    return inf;
  }

  //
  // Initialize rtInfF needed by the generated code.
  // Inf is initialized as non-signaling. Assumes IEEE.
  //
  static real32_T rtGetInfF(void)
  {
    IEEESingle infF;
    infF.wordL.wordLuint = 0x7F800000U;
    return infF.wordL.wordLreal;
  }

  //
  // Initialize rtMinusInf needed by the generated code.
  // Inf is initialized as non-signaling. Assumes IEEE.
  //
  static real_T rtGetMinusInf(void)
  {
    size_t bitsPerReal{ sizeof(real_T) * (NumBitsPerChar) };

    real_T minf{ 0.0 };

    if (bitsPerReal == 32U) {
      minf = rtGetMinusInfF();
    } else {
      union {
        LittleEndianIEEEDouble bitVal;
        real_T fltVal;
      } tmpVal;

      tmpVal.bitVal.words.wordH = 0xFFF00000U;
      tmpVal.bitVal.words.wordL = 0x00000000U;
      minf = tmpVal.fltVal;
    }

    return minf;
  }

  //
  // Initialize rtMinusInfF needed by the generated code.
  // Inf is initialized as non-signaling. Assumes IEEE.
  //
  static real32_T rtGetMinusInfF(void)
  {
    IEEESingle minfF;
    minfF.wordL.wordLuint = 0xFF800000U;
    return minfF.wordL.wordLreal;
  }
}

static int32_T div_nde_s32_floor(int32_T numerator, int32_T denominator)
{
  return (((numerator < 0) != (denominator < 0)) && (numerator % denominator !=
           0) ? -1 : 0) + numerator / denominator;
}

// Function for MATLAB Function: '<S22>/optimizer'
real_T thesis_mpc::norm(const real_T x[3])
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  scale = 3.3121686421112381E-170;
  absxk = std::abs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }

  absxk = std::abs(x[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  absxk = std::abs(x[2]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  return scale * std::sqrt(y);
}

// Function for MATLAB Function: '<S22>/optimizer'
real_T thesis_mpc::maximum(const real_T x[3])
{
  real_T ex;
  int32_T idx;
  int32_T k;
  if (!std::isnan(x[0])) {
    idx = 1;
  } else {
    boolean_T exitg1;
    idx = 0;
    k = 2;
    exitg1 = false;
    while ((!exitg1) && (k < 4)) {
      if (!std::isnan(x[k - 1])) {
        idx = k;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }

  if (idx == 0) {
    ex = x[0];
  } else {
    ex = x[idx - 1];
    for (k = idx + 1; k < 4; k++) {
      real_T x_0;
      x_0 = x[k - 1];
      if (ex < x_0) {
        ex = x_0;
      }
    }
  }

  return ex;
}

// Function for MATLAB Function: '<S22>/optimizer'
real_T thesis_mpc::xnrm2(int32_T n, const real_T x[9], int32_T ix0)
{
  real_T y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      real_T scale;
      int32_T kend;
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (int32_T k{ix0}; k <= kend; k++) {
        real_T absxk;
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
          real_T t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          real_T t;
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

real_T rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T a;
  real_T b;
  real_T y;
  a = std::abs(u0);
  b = std::abs(u1);
  if (a < b) {
    a /= b;
    y = std::sqrt(a * a + 1.0) * b;
  } else if (a > b) {
    b /= a;
    y = std::sqrt(b * b + 1.0) * a;
  } else if (std::isnan(b)) {
    y = (rtNaN);
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

// Function for MATLAB Function: '<S22>/optimizer'
void thesis_mpc::xgemv(int32_T b_m, int32_T n, const real_T b_A[9], int32_T ia0,
  const real_T x[9], int32_T ix0, real_T y[3])
{
  if ((b_m != 0) && (n != 0)) {
    int32_T b;
    if (n - 1 >= 0) {
      std::memset(&y[0], 0, static_cast<uint32_T>(n) * sizeof(real_T));
    }

    b = (n - 1) * 3 + ia0;
    for (int32_T b_iy{ia0}; b_iy <= b; b_iy += 3) {
      real_T c;
      int32_T d;
      int32_T iyend;
      c = 0.0;
      d = b_iy + b_m;
      for (iyend = b_iy; iyend < d; iyend++) {
        c += x[((ix0 + iyend) - b_iy) - 1] * b_A[iyend - 1];
      }

      iyend = div_nde_s32_floor(b_iy - ia0, 3);
      y[iyend] += c;
    }
  }
}

// Function for MATLAB Function: '<S22>/optimizer'
void thesis_mpc::xgerc(int32_T b_m, int32_T n, real_T alpha1, int32_T ix0, const
  real_T y[3], real_T b_A[9], int32_T ia0)
{
  if (!(alpha1 == 0.0)) {
    int32_T jA;
    jA = ia0;
    for (int32_T j{0}; j < n; j++) {
      real_T temp;
      temp = y[j];
      if (temp != 0.0) {
        int32_T b;
        temp *= alpha1;
        b = b_m + jA;
        for (int32_T ijA{jA}; ijA < b; ijA++) {
          b_A[ijA - 1] += b_A[((ix0 + ijA) - jA) - 1] * temp;
        }
      }

      jA += 3;
    }
  }
}

// Function for MATLAB Function: '<S22>/optimizer'
void thesis_mpc::KWIKfactor(const real_T b_Ac[204], const int32_T iC[68],
  int32_T nA, const real_T b_Linv[9], real_T D[9], real_T b_H[9], int32_T n,
  real_T RLinv[9], real_T *Status)
{
  real_T Q[9];
  real_T R[9];
  real_T TL[9];
  real_T b_A[9];
  real_T tau[3];
  real_T work[3];
  real_T atmp;
  real_T b_A_0;
  real_T beta1;
  int32_T b_coltop;
  int32_T b_lastv;
  int32_T coltop;
  int32_T exitg1;
  int32_T ii;
  int32_T k_i;
  int32_T knt;
  boolean_T exitg2;
  *Status = 1.0;
  std::memset(&RLinv[0], 0, 9U * sizeof(real_T));
  for (k_i = 0; k_i < nA; k_i++) {
    b_lastv = iC[k_i];
    for (b_coltop = 0; b_coltop < 3; b_coltop++) {
      RLinv[b_coltop + 3 * k_i] = (b_Ac[b_lastv - 1] * b_Linv[b_coltop] +
        b_Linv[b_coltop + 3] * b_Ac[b_lastv + 67]) + b_Linv[b_coltop + 6] *
        b_Ac[b_lastv + 135];
    }
  }

  std::memcpy(&b_A[0], &RLinv[0], 9U * sizeof(real_T));
  tau[0] = 0.0;
  work[0] = 0.0;
  tau[1] = 0.0;
  work[1] = 0.0;
  tau[2] = 0.0;
  work[2] = 0.0;
  for (k_i = 0; k_i < 3; k_i++) {
    ii = k_i * 3 + k_i;
    if (k_i + 1 < 3) {
      atmp = b_A[ii];
      b_lastv = ii + 2;
      tau[k_i] = 0.0;
      beta1 = xnrm2(2 - k_i, b_A, ii + 2);
      if (beta1 != 0.0) {
        b_A_0 = b_A[ii];
        beta1 = rt_hypotd_snf(b_A_0, beta1);
        if (b_A_0 >= 0.0) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 1.0020841800044864E-292) {
          knt = 0;
          coltop = (ii - k_i) + 3;
          do {
            knt++;
            for (b_coltop = b_lastv; b_coltop <= coltop; b_coltop++) {
              b_A[b_coltop - 1] *= 9.9792015476736E+291;
            }

            beta1 *= 9.9792015476736E+291;
            atmp *= 9.9792015476736E+291;
          } while ((std::abs(beta1) < 1.0020841800044864E-292) && (knt < 20));

          beta1 = rt_hypotd_snf(atmp, xnrm2(2 - k_i, b_A, ii + 2));
          if (atmp >= 0.0) {
            beta1 = -beta1;
          }

          tau[k_i] = (beta1 - atmp) / beta1;
          atmp = 1.0 / (atmp - beta1);
          for (b_coltop = b_lastv; b_coltop <= coltop; b_coltop++) {
            b_A[b_coltop - 1] *= atmp;
          }

          for (b_lastv = 0; b_lastv < knt; b_lastv++) {
            beta1 *= 1.0020841800044864E-292;
          }

          atmp = beta1;
        } else {
          tau[k_i] = (beta1 - b_A_0) / beta1;
          atmp = 1.0 / (b_A_0 - beta1);
          b_coltop = (ii - k_i) + 3;
          for (knt = b_lastv; knt <= b_coltop; knt++) {
            b_A[knt - 1] *= atmp;
          }

          atmp = beta1;
        }
      }

      b_A[ii] = 1.0;
      if (tau[k_i] != 0.0) {
        b_lastv = 3 - k_i;
        knt = (ii - k_i) + 2;
        while ((b_lastv > 0) && (b_A[knt] == 0.0)) {
          b_lastv--;
          knt--;
        }

        knt = 2 - k_i;
        exitg2 = false;
        while ((!exitg2) && (knt > 0)) {
          b_coltop = ((knt - 1) * 3 + ii) + 3;
          coltop = b_coltop;
          do {
            exitg1 = 0;
            if (coltop + 1 <= b_coltop + b_lastv) {
              if (b_A[coltop] != 0.0) {
                exitg1 = 1;
              } else {
                coltop++;
              }
            } else {
              knt--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        b_lastv = 0;
        knt = 0;
      }

      if (b_lastv > 0) {
        xgemv(b_lastv, knt, b_A, ii + 4, b_A, ii + 1, work);
        xgerc(b_lastv, knt, -tau[k_i], ii + 1, work, b_A, ii + 4);
      }

      b_A[ii] = atmp;
    } else {
      tau[2] = 0.0;
    }
  }

  for (k_i = 0; k_i < 3; k_i++) {
    for (ii = 0; ii <= k_i; ii++) {
      R[ii + 3 * k_i] = b_A[3 * k_i + ii];
    }

    for (ii = k_i + 2; ii < 4; ii++) {
      R[(ii + 3 * k_i) - 1] = 0.0;
    }

    work[k_i] = 0.0;
  }

  for (k_i = 2; k_i >= 0; k_i--) {
    b_lastv = (k_i * 3 + k_i) + 4;
    if (k_i + 1 < 3) {
      b_A[b_lastv - 4] = 1.0;
      if (tau[k_i] != 0.0) {
        knt = 3 - k_i;
        b_coltop = b_lastv - k_i;
        while ((knt > 0) && (b_A[b_coltop - 2] == 0.0)) {
          knt--;
          b_coltop--;
        }

        b_coltop = 2 - k_i;
        exitg2 = false;
        while ((!exitg2) && (b_coltop > 0)) {
          coltop = (b_coltop - 1) * 3 + b_lastv;
          ii = coltop;
          do {
            exitg1 = 0;
            if (ii <= (coltop + knt) - 1) {
              if (b_A[ii - 1] != 0.0) {
                exitg1 = 1;
              } else {
                ii++;
              }
            } else {
              b_coltop--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        knt = 0;
        b_coltop = 0;
      }

      if (knt > 0) {
        xgemv(knt, b_coltop, b_A, b_lastv, b_A, b_lastv - 3, work);
        xgerc(knt, b_coltop, -tau[k_i], b_lastv - 3, work, b_A, b_lastv);
      }

      b_coltop = (b_lastv - k_i) - 1;
      for (knt = b_lastv - 2; knt <= b_coltop; knt++) {
        b_A[knt - 1] *= -tau[k_i];
      }
    }

    b_A[b_lastv - 4] = 1.0 - tau[k_i];
    for (knt = 0; knt < k_i; knt++) {
      b_A[(b_lastv - knt) - 5] = 0.0;
    }
  }

  for (k_i = 0; k_i < 3; k_i++) {
    Q[3 * k_i] = b_A[3 * k_i];
    b_lastv = 3 * k_i + 1;
    Q[b_lastv] = b_A[b_lastv];
    b_lastv = 3 * k_i + 2;
    Q[b_lastv] = b_A[b_lastv];
  }

  k_i = 0;
  do {
    exitg1 = 0;
    if (k_i <= nA - 1) {
      if (std::abs(R[3 * k_i + k_i]) < 1.0E-12) {
        *Status = -2.0;
        exitg1 = 1;
      } else {
        k_i++;
      }
    } else {
      for (k_i = 0; k_i < n; k_i++) {
        for (ii = 0; ii < n; ii++) {
          TL[k_i + 3 * ii] = (b_Linv[3 * k_i + 1] * Q[3 * ii + 1] + b_Linv[3 *
                              k_i] * Q[3 * ii]) + b_Linv[3 * k_i + 2] * Q[3 * ii
            + 2];
        }
      }

      std::memset(&RLinv[0], 0, 9U * sizeof(real_T));
      for (k_i = nA; k_i >= 1; k_i--) {
        b_coltop = (k_i - 1) * 3;
        knt = (k_i + b_coltop) - 1;
        RLinv[knt] = 1.0;
        for (ii = k_i; ii <= nA; ii++) {
          coltop = ((ii - 1) * 3 + k_i) - 1;
          RLinv[coltop] /= R[knt];
        }

        if (k_i > 1) {
          for (ii = 0; ii <= k_i - 2; ii++) {
            for (b_lastv = k_i; b_lastv <= nA; b_lastv++) {
              knt = (b_lastv - 1) * 3;
              coltop = knt + ii;
              RLinv[coltop] -= RLinv[(knt + k_i) - 1] * R[b_coltop + ii];
            }
          }
        }
      }

      for (k_i = 0; k_i < n; k_i++) {
        for (ii = k_i + 1; ii <= n; ii++) {
          b_coltop = (ii - 1) * 3 + k_i;
          b_H[b_coltop] = 0.0;
          for (b_lastv = nA + 1; b_lastv <= n; b_lastv++) {
            knt = (b_lastv - 1) * 3;
            b_H[b_coltop] -= TL[(knt + ii) - 1] * TL[knt + k_i];
          }

          b_H[(ii + 3 * k_i) - 1] = b_H[b_coltop];
        }
      }

      for (k_i = 0; k_i < nA; k_i++) {
        for (ii = 0; ii < n; ii++) {
          b_coltop = 3 * k_i + ii;
          D[b_coltop] = 0.0;
          for (b_lastv = k_i + 1; b_lastv <= nA; b_lastv++) {
            knt = (b_lastv - 1) * 3;
            D[b_coltop] += TL[knt + ii] * RLinv[knt + k_i];
          }
        }
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);
}

// Function for MATLAB Function: '<S22>/optimizer'
void thesis_mpc::DropConstraint(int32_T kDrop, boolean_T iA[68], int32_T *nA,
  int32_T iC[68])
{
  if (kDrop > 0) {
    iA[iC[kDrop - 1] - 1] = false;
    if (kDrop < *nA) {
      for (int32_T i{kDrop}; i < *nA; i++) {
        iC[i - 1] = iC[i];
      }
    }

    iC[*nA - 1] = 0;
    (*nA)--;
  }
}

// Function for MATLAB Function: '<S22>/optimizer'
void thesis_mpc::qpkwik(const real_T b_Linv[9], const real_T b_Hinv[9], const
  real_T f[3], const real_T b_Ac[204], const real_T b[68], boolean_T iA[68],
  int32_T maxiter, real_T FeasTol, real_T x[3], real_T lambda[68], int32_T
  *status)
{
  real_T cTol[68];
  real_T D[9];
  real_T RLinv[9];
  real_T U[9];
  real_T b_H[9];
  real_T Opt[6];
  real_T Rhs[6];
  real_T r[3];
  real_T z[3];
  real_T Xnorm0;
  real_T cMin;
  real_T cVal;
  real_T cVal_tmp;
  real_T rMin;
  real_T t;
  int32_T iC[68];
  int32_T U_tmp;
  int32_T b_exponent;
  int32_T exitg1;
  int32_T exitg3;
  int32_T exponent;
  int32_T i;
  int32_T iC_0;
  int32_T iSave;
  int32_T nA;
  int32_T tmp;
  boolean_T ColdReset;
  boolean_T DualFeasible;
  boolean_T cTolComputed;
  boolean_T exitg2;
  boolean_T exitg4;
  boolean_T guard1{ false };

  boolean_T guard2{ false };

  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.0;
  *status = 1;
  r[0] = 0.0;
  r[1] = 0.0;
  r[2] = 0.0;
  rMin = 0.0;
  cTolComputed = false;
  for (i = 0; i < 68; i++) {
    lambda[i] = 0.0;
    cTol[i] = 1.0;
    iC[i] = 0;
  }

  nA = 0;
  for (tmp = 0; tmp < 68; tmp++) {
    if (iA[tmp]) {
      nA++;
      iC[nA - 1] = tmp + 1;
    }
  }

  guard1 = false;
  if (nA > 0) {
    for (i = 0; i < 6; i++) {
      Opt[i] = 0.0;
    }

    Rhs[0] = f[0];
    Rhs[3] = 0.0;
    Rhs[1] = f[1];
    Rhs[4] = 0.0;
    Rhs[2] = f[2];
    Rhs[5] = 0.0;
    DualFeasible = false;
    tmp = static_cast<int32_T>(std::round(0.3 * static_cast<real_T>(nA)));
    ColdReset = false;
    do {
      exitg3 = 0;
      if ((!DualFeasible) && (nA > 0) && (*status <= maxiter)) {
        KWIKfactor(b_Ac, iC, nA, b_Linv, D, b_H, degrees, RLinv, &Xnorm0);
        if (Xnorm0 < 0.0) {
          if (ColdReset) {
            *status = -2;
            exitg3 = 2;
          } else {
            nA = 0;
            std::memset(&iA[0], 0, 68U * sizeof(boolean_T));
            std::memset(&iC[0], 0, 68U * sizeof(int32_T));
            ColdReset = true;
          }
        } else {
          for (i = 0; i < nA; i++) {
            Rhs[i + 3] = b[iC[i] - 1];
            for (iSave = i + 1; iSave <= nA; iSave++) {
              U_tmp = (3 * i + iSave) - 1;
              U[U_tmp] = 0.0;
              for (iC_0 = 0; iC_0 < nA; iC_0++) {
                U[U_tmp] += RLinv[(3 * iC_0 + iSave) - 1] * RLinv[3 * iC_0 + i];
              }

              U[i + 3 * (iSave - 1)] = U[U_tmp];
            }
          }

          for (i = 0; i < 3; i++) {
            Opt[i] = (b_H[i + 3] * Rhs[1] + b_H[i] * Rhs[0]) + b_H[i + 6] * Rhs
              [2];
            for (iSave = 0; iSave < nA; iSave++) {
              Opt[i] += D[3 * iSave + i] * Rhs[iSave + 3];
            }
          }

          for (i = 0; i < nA; i++) {
            Opt[i + 3] = (D[3 * i + 1] * Rhs[1] + D[3 * i] * Rhs[0]) + D[3 * i +
              2] * Rhs[2];
            for (iSave = 0; iSave < nA; iSave++) {
              Opt[i + 3] += U[3 * iSave + i] * Rhs[iSave + 3];
            }
          }

          Xnorm0 = -1.0E-12;
          i = -1;
          for (iSave = 0; iSave < nA; iSave++) {
            cMin = Opt[iSave + 3];
            lambda[iC[iSave] - 1] = cMin;
            if ((cMin < Xnorm0) && (iSave + 1 <= nA)) {
              i = iSave;
              Xnorm0 = cMin;
            }
          }

          if (i + 1 <= 0) {
            DualFeasible = true;
            x[0] = Opt[0];
            x[1] = Opt[1];
            x[2] = Opt[2];
          } else {
            (*status)++;
            if (tmp <= 5) {
              iC_0 = 5;
            } else {
              iC_0 = tmp;
            }

            if (*status > iC_0) {
              nA = 0;
              std::memset(&iA[0], 0, 68U * sizeof(boolean_T));
              std::memset(&iC[0], 0, 68U * sizeof(int32_T));
              ColdReset = true;
            } else {
              lambda[iC[i] - 1] = 0.0;
              DropConstraint(i + 1, iA, &nA, iC);
            }
          }
        }
      } else {
        if (nA <= 0) {
          std::memset(&lambda[0], 0, 68U * sizeof(real_T));
          Xnorm0 = f[1];
          cMin = f[0];
          cVal = f[2];
          for (tmp = 0; tmp < 3; tmp++) {
            x[tmp] = (-b_Hinv[tmp + 3] * Xnorm0 + -b_Hinv[tmp] * cMin) +
              -b_Hinv[tmp + 6] * cVal;
          }
        }

        exitg3 = 1;
      }
    } while (exitg3 == 0);

    if (exitg3 == 1) {
      guard1 = true;
    }
  } else {
    Xnorm0 = f[1];
    cMin = f[0];
    cVal = f[2];
    for (tmp = 0; tmp < 3; tmp++) {
      x[tmp] = (-b_Hinv[tmp + 3] * Xnorm0 + -b_Hinv[tmp] * cMin) + -b_Hinv[tmp +
        6] * cVal;
    }

    guard1 = true;
  }

  if (guard1) {
    Xnorm0 = norm(x);
    exitg2 = false;
    while ((!exitg2) && (*status <= maxiter)) {
      cMin = -FeasTol;
      tmp = -1;
      for (i = 0; i < 68; i++) {
        if (!cTolComputed) {
          z[0] = std::abs(b_Ac[i] * x[0]);
          z[1] = std::abs(b_Ac[i + 68] * x[1]);
          z[2] = std::abs(b_Ac[i + 136] * x[2]);
          cTol[i] = std::fmax(cTol[i], maximum(z));
        }

        if (!iA[i]) {
          cVal = (((b_Ac[i + 68] * x[1] + b_Ac[i] * x[0]) + b_Ac[i + 136] * x[2])
                  - b[i]) / cTol[i];
          if (cVal < cMin) {
            cMin = cVal;
            tmp = i;
          }
        }
      }

      cTolComputed = true;
      if (tmp + 1 <= 0) {
        exitg2 = true;
      } else if (*status == maxiter) {
        *status = 0;
        exitg2 = true;
      } else {
        do {
          exitg1 = 0;
          if ((tmp + 1 > 0) && (*status <= maxiter)) {
            guard2 = false;
            if (nA == 0) {
              for (iC_0 = 0; iC_0 < 3; iC_0++) {
                z[iC_0] = (b_Hinv[iC_0 + 3] * b_Ac[tmp + 68] + b_Hinv[iC_0] *
                           b_Ac[tmp]) + b_Hinv[iC_0 + 6] * b_Ac[tmp + 136];
              }

              guard2 = true;
            } else {
              KWIKfactor(b_Ac, iC, nA, b_Linv, D, b_H, degrees, RLinv, &cMin);
              if (cMin <= 0.0) {
                *status = -2;
                exitg1 = 1;
              } else {
                for (iC_0 = 0; iC_0 < 9; iC_0++) {
                  U[iC_0] = -b_H[iC_0];
                }

                for (iC_0 = 0; iC_0 < 3; iC_0++) {
                  z[iC_0] = (U[iC_0 + 3] * b_Ac[tmp + 68] + U[iC_0] * b_Ac[tmp])
                    + U[iC_0 + 6] * b_Ac[tmp + 136];
                }

                for (i = 0; i < nA; i++) {
                  r[i] = (D[3 * i + 1] * b_Ac[tmp + 68] + D[3 * i] * b_Ac[tmp])
                    + D[3 * i + 2] * b_Ac[tmp + 136];
                }

                guard2 = true;
              }
            }

            if (guard2) {
              i = 0;
              cMin = 0.0;
              DualFeasible = true;
              ColdReset = true;
              if (nA > 0) {
                iSave = 0;
                exitg4 = false;
                while ((!exitg4) && (iSave <= nA - 1)) {
                  if (r[iSave] >= 1.0E-12) {
                    ColdReset = false;
                    exitg4 = true;
                  } else {
                    iSave++;
                  }
                }
              }

              if ((nA != 0) && (!ColdReset)) {
                for (iSave = 0; iSave < nA; iSave++) {
                  cVal = r[iSave];
                  if (cVal > 1.0E-12) {
                    cVal = lambda[iC[iSave] - 1] / cVal;
                    if ((i == 0) || (cVal < rMin)) {
                      rMin = cVal;
                      i = iSave + 1;
                    }
                  }
                }

                if (i > 0) {
                  cMin = rMin;
                  DualFeasible = false;
                }
              }

              t = b_Ac[tmp + 68];
              cVal_tmp = b_Ac[tmp + 136];
              cVal = (t * z[1] + z[0] * b_Ac[tmp]) + cVal_tmp * z[2];
              if (cVal <= 0.0) {
                cVal = 0.0;
                ColdReset = true;
              } else {
                cVal = (b[tmp] - ((t * x[1] + b_Ac[tmp] * x[0]) + cVal_tmp * x[2]))
                  / cVal;
                ColdReset = false;
              }

              if (DualFeasible && ColdReset) {
                *status = -1;
                exitg1 = 1;
              } else {
                if (ColdReset) {
                  t = cMin;
                } else if (DualFeasible) {
                  t = cVal;
                } else if (cMin < cVal) {
                  t = cMin;
                } else {
                  t = cVal;
                }

                for (iSave = 0; iSave < nA; iSave++) {
                  iC_0 = iC[iSave];
                  lambda[iC_0 - 1] -= t * r[iSave];
                  if ((iC_0 <= 68) && (lambda[iC_0 - 1] < 0.0)) {
                    lambda[iC_0 - 1] = 0.0;
                  }
                }

                lambda[tmp] += t;
                std::frexp(1.0, &exponent);
                if (std::abs(t - cMin) < 2.2204460492503131E-16) {
                  DropConstraint(i, iA, &nA, iC);
                }

                if (!ColdReset) {
                  x[0] += t * z[0];
                  x[1] += t * z[1];
                  x[2] += t * z[2];
                  std::frexp(1.0, &b_exponent);
                  if (std::abs(t - cVal) < 2.2204460492503131E-16) {
                    if (nA == degrees) {
                      *status = -1;
                      exitg1 = 1;
                    } else {
                      nA++;
                      iC[nA - 1] = tmp + 1;
                      i = nA - 1;
                      exitg4 = false;
                      while ((!exitg4) && (i + 1 > 1)) {
                        iC_0 = iC[i - 1];
                        if (iC[i] > iC_0) {
                          exitg4 = true;
                        } else {
                          iSave = iC[i];
                          iC[i] = iC_0;
                          iC[i - 1] = iSave;
                          i--;
                        }
                      }

                      iA[tmp] = true;
                      tmp = -1;
                      (*status)++;
                    }
                  } else {
                    (*status)++;
                  }
                } else {
                  (*status)++;
                }
              }
            }
          } else {
            cMin = norm(x);
            if (std::abs(cMin - Xnorm0) > 0.001) {
              Xnorm0 = cMin;
              for (tmp = 0; tmp < 68; tmp++) {
                cTol[tmp] = std::fmax(std::abs(b[tmp]), 1.0);
              }

              cTolComputed = false;
            }

            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    }
  }
}

// Model step function
void thesis_mpc::step()
{
  static const real_T c_a[476]{ -0.31978899948757172, -0.068365678762982746,
    0.17881357345539545, 0.41167960337812803, 0.62586983456882017,
    0.81903668888730707, 0.98971567748608857, 1.1369723378231549,
    1.2602842049449765, 1.359491976874414, 1.43476953263987, 1.4865977601250475,
    1.5157379405328502, 1.523203732688893, 1.5102318029100819,
    1.4782514436315433, 1.4288536022631566, 1.363759749679756,
    1.2847910012358081, 1.193837876693665, 1.0928310538213519,
    0.9837134355746795, 0.868413813750201, 0.748822373514156,
    0.62676824390809438, 0.50399925988224847, 0.38216406215065185,
    0.26279662269545856, 0.14730324652769147, 0.036952064748414386,
    0.31978899948757172, 0.068365678762982746, -0.17881357345539545,
    -0.41167960337812803, -0.62586983456882017, -0.81903668888730707,
    -0.98971567748608857, -1.1369723378231549, -1.2602842049449765,
    -1.359491976874414, -1.43476953263987, -1.4865977601250475,
    -1.5157379405328502, -1.523203732688893, -1.5102318029100819,
    -1.4782514436315433, -1.4288536022631566, -1.363759749679756,
    -1.2847910012358081, -1.193837876693665, -1.0928310538213519,
    -0.9837134355746795, -0.868413813750201, -0.748822373514156,
    -0.62676824390809438, -0.50399925988224847, -0.38216406215065185,
    -0.26279662269545856, -0.14730324652769147, -0.036952064748414386, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.439210834666996, -0.48802935931175,
    -0.53932010848674383, -0.58455603532415523, -0.62109360892294585,
    -0.648186349035838, -0.66575912232416246, -0.674025291605441,
    -0.67336087354876784, -0.66425685141627788, -0.647295136129868,
    -0.62313165604405463, -0.59248165551289655, -0.5561058184692057,
    -0.51479690455447624, -0.46936690266154085, -0.42063479298415818,
    -0.36941502376280133, -0.3165068019286652, -0.26268428294226254,
    -0.20868772909905989, -0.15521568905547872, -0.10291823495462313,
    -0.052391277629405993, -0.0041719651487191868, 0.0412648443807931,
    0.08350905819941877, 0.1222168071485374, 0.15711093251202335,
    0.18798055692061491, 0.439210834666996, 0.48802935931175,
    0.53932010848674383, 0.58455603532415523, 0.62109360892294585,
    0.648186349035838, 0.66575912232416246, 0.674025291605441,
    0.67336087354876784, 0.66425685141627788, 0.647295136129868,
    0.62313165604405463, 0.59248165551289655, 0.5561058184692057,
    0.51479690455447624, 0.46936690266154085, 0.42063479298415818,
    0.36941502376280133, 0.3165068019286652, 0.26268428294226254,
    0.20868772909905989, 0.15521568905547872, 0.10291823495462313,
    0.052391277629405993, 0.0041719651487191868, -0.0412648443807931,
    -0.08350905819941877, -0.1222168071485374, -0.15711093251202335,
    -0.18798055692061491, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.090079849293012618, -0.00011645920246529963, 0.07735314577705163,
    0.14664814243853705, 0.20879098187526274, 0.26387372413885773,
    0.31177937427191987, 0.35239970721430908, 0.38569843906978013,
    0.41172643024889305, 0.43062161698434026, 0.44260376307207006,
    0.44796721738037115, 0.4470727514187387, 0.44033890382915242,
    0.428233057349559, 0.41126240777348022, 0.38996495884312721,
    0.364900663012939, 0.33664281691863429, 0.3057698096582292,
    0.27285731099184957, 0.238470975224461, 0.20315972494813739,
    0.16744966713292717, 0.13183868241989438, 0.096791717033126182,
    0.0627367956257039, 0.030061762731726172, -0.000888249577261957,
    0.090079849293012618, 0.00011645920246529963, -0.07735314577705163,
    -0.14664814243853705, -0.20879098187526274, -0.26387372413885773,
    -0.31177937427191987, -0.35239970721430908, -0.38569843906978013,
    -0.41172643024889305, -0.43062161698434026, -0.44260376307207006,
    -0.44796721738037115, -0.4470727514187387, -0.44033890382915242,
    -0.428233057349559, -0.41126240777348022, -0.38996495884312721,
    -0.364900663012939, -0.33664281691863429, -0.3057698096582292,
    -0.27285731099184957, -0.238470975224461, -0.20315972494813739,
    -0.16744966713292717, -0.13183868241989438, -0.096791717033126182,
    -0.0627367956257039, -0.030061762731726172, 0.000888249577261957, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.11634721447929182, 0.21919483097646641,
    0.3127544242836352, 0.39633465982224708, 0.46945263927007219,
    0.53180304617785579, 0.58324369164548651, 0.62378495577496751,
    0.65357945480462742, 0.672910966211207, 0.68218247075081728,
    0.6819034219132255, 0.6726764270851503, 0.655183542468182, 0.630172382806373,
    0.59844223883783632, 0.56083038382678951, 0.51819873689375984,
    0.47142103574764554, 0.42137065521489042, 0.36890919096664593,
    0.31487591034471024, 0.26007815444474791, 0.20528275787921815,
    0.15120853514483279, 0.098519865474690771, 0.04782139165856493,
    -0.00034616725947667895, -0.045509104096147474, -0.087262740640512593,
    -0.11634721447929182, -0.21919483097646641, -0.3127544242836352,
    -0.39633465982224708, -0.46945263927007219, -0.53180304617785579,
    -0.58324369164548651, -0.62378495577496751, -0.65357945480462742,
    -0.672910966211207, -0.68218247075081728, -0.6819034219132255,
    -0.6726764270851503, -0.655183542468182, -0.630172382806373,
    -0.59844223883783632, -0.56083038382678951, -0.51819873689375984,
    -0.47142103574764554, -0.42137065521489042, -0.36890919096664593,
    -0.31487591034471024, -0.26007815444474791, -0.20528275787921815,
    -0.15120853514483279, -0.098519865474690771, -0.04782139165856493,
    0.00034616725947667895, 0.045509104096147474, 0.087262740640512593, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0039820934613680374, 0.02076682783683589,
    0.037063921626778393, 0.048921696155707822, 0.055296364304163383,
    0.056081147056626174, 0.051497817523871836, 0.041904397387860184,
    0.027729972236831497, 0.00944839823883413, -0.012435979523164041,
    -0.037398738962822156, -0.06490434676072479, -0.09441436166444557,
    -0.12539503489606779, -0.15732431483435785, -0.18969822437312328,
    -0.22203657528742954, -0.25388799034163445, -0.28483421325274683,
    -0.31449369660469667, -0.34252446759280547, -0.36862628075287523,
    -0.3925420754456877, -0.4140587637341423, -0.43300738134377437,
    -0.4492626405934998, -0.46274192949293452, -0.47340380561026069,
    -0.48124603681803507, -0.0039820934613680374, -0.02076682783683589,
    -0.037063921626778393, -0.048921696155707822, -0.055296364304163383,
    -0.056081147056626174, -0.051497817523871836, -0.041904397387860184,
    -0.027729972236831497, -0.00944839823883413, 0.012435979523164041,
    0.037398738962822156, 0.06490434676072479, 0.09441436166444557,
    0.12539503489606779, 0.15732431483435785, 0.18969822437312328,
    0.22203657528742954, 0.25388799034163445, 0.28483421325274683,
    0.31449369660469667, 0.34252446759280547, 0.36862628075287523,
    0.3925420754456877, 0.4140587637341423, 0.43300738134377437,
    0.4492626405934998, 0.46274192949293452, 0.47340380561026069,
    0.48124603681803507, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0023885274222311755, 0.01251982403342386, 0.022444039183392425,
    0.029775847695277642, 0.033863538476262348, 0.034615149377359371,
    0.032136268135724386, 0.026615723860526408, 0.018286967506234021,
    0.0074125771725973531, -0.0057241253475283412, -0.020825196192422329,
    -0.037583371317974339, -0.055686953694909141, -0.074824378259047808,
    -0.094688455948463077, -0.11498027481676629, -0.13541273359433281,
    -0.15571368674614236, -0.17562868557920067, -0.19492330586726134,
    -0.21338505830588483, -0.23082488369879051, -0.24707824001388437,
    -0.26200579327605861, -0.27549372864213495, -0.28745370190029429,
    -0.29782245503011839, -0.306561122336323, -0.31365425602346708,
    -0.0023885274222311755, -0.01251982403342386, -0.022444039183392425,
    -0.029775847695277642, -0.033863538476262348, -0.034615149377359371,
    -0.032136268135724386, -0.026615723860526408, -0.018286967506234021,
    -0.0074125771725973531, 0.0057241253475283412, 0.020825196192422329,
    0.037583371317974339, 0.055686953694909141, 0.074824378259047808,
    0.094688455948463077, 0.11498027481676629, 0.13541273359433281,
    0.15571368674614236, 0.17562868557920067, 0.19492330586726134,
    0.21338505830588483, 0.23082488369879051, 0.24707824001388437,
    0.26200579327605861, 0.27549372864213495, 0.28745370190029429,
    0.29782245503011839, 0.306561122336323, 0.31365425602346708, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0 };

  static const real_T b_Ac[204]{ 0.0082393983319562145, 0.043241643790057967,
    0.0776020389914959, 0.10308013732709863, 0.11740668642495948,
    0.12024122330596915, 0.11192495647043646, 0.093087329418972228,
    0.064513214663213, 0.0270897107947072, -0.018222606641386398,
    -0.070410680681089924, -0.12842718008483034, -0.191207311425991,
    -0.25768455254736161, -0.32680530976127503, -0.39754242054157474,
    -0.46890741425646021, -0.53996145601447953, -0.60982491755695356,
    -0.677685539447547, -0.74280516892895021, -0.80452507708795107,
    -0.86226987707030656, -0.91555008181220787, -0.96396335495883223,
    -1.0071945222127245, -1.045014422216622, -1.0772776861761704,
    -1.1039195437437468, -0.0082393983319562145, -0.043241643790057967,
    -0.0776020389914959, -0.10308013732709863, -0.11740668642495948,
    -0.12024122330596915, -0.11192495647043646, -0.093087329418972228,
    -0.064513214663213, -0.0270897107947072, 0.018222606641386398,
    0.070410680681089924, 0.12842718008483034, 0.191207311425991,
    0.25768455254736161, 0.32680530976127503, 0.39754242054157474,
    0.46890741425646021, 0.53996145601447953, 0.60982491755695356,
    0.677685539447547, 0.74280516892895021, 0.80452507708795107,
    0.86226987707030656, 0.91555008181220787, 0.96396335495883223,
    1.0071945222127245, 1.045014422216622, 1.0772776861761704,
    1.1039195437437468, -1.0, -1.0, 1.0, 1.0, -1.0, -0.0, 1.0, 0.0, -0.0,
    0.0082393983319562145, 0.043241643790057967, 0.0776020389914959,
    0.10308013732709863, 0.11740668642495948, 0.12024122330596915,
    0.11192495647043646, 0.093087329418972228, 0.064513214663213,
    0.0270897107947072, -0.018222606641386398, -0.070410680681089924,
    -0.12842718008483034, -0.191207311425991, -0.25768455254736161,
    -0.32680530976127503, -0.39754242054157474, -0.46890741425646021,
    -0.53996145601447953, -0.60982491755695356, -0.677685539447547,
    -0.74280516892895021, -0.80452507708795107, -0.86226987707030656,
    -0.91555008181220787, -0.96396335495883223, -1.0071945222127245,
    -1.045014422216622, -1.0772776861761704, 0.0, -0.0082393983319562145,
    -0.043241643790057967, -0.0776020389914959, -0.10308013732709863,
    -0.11740668642495948, -0.12024122330596915, -0.11192495647043646,
    -0.093087329418972228, -0.064513214663213, -0.0270897107947072,
    0.018222606641386398, 0.070410680681089924, 0.12842718008483034,
    0.191207311425991, 0.25768455254736161, 0.32680530976127503,
    0.39754242054157474, 0.46890741425646021, 0.53996145601447953,
    0.60982491755695356, 0.677685539447547, 0.74280516892895021,
    0.80452507708795107, 0.86226987707030656, 0.91555008181220787,
    0.96396335495883223, 1.0071945222127245, 1.045014422216622,
    1.0772776861761704, -0.0, -1.0, 0.0, 1.0, -0.0, -1.0, 0.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T b_Mlim_0[68]{ 1.0471975511966, 1.0471975511966,
    1.0471975511966, 1.0471975511966, 1.0471975511966, 1.0471975511966,
    1.0471975511966, 1.0471975511966, 1.0471975511966, 1.0471975511966,
    1.0471975511966, 1.0471975511966, 1.0471975511966, 1.0471975511966,
    1.0471975511966, 1.0471975511966, 1.0471975511966, 1.0471975511966,
    1.0471975511966, 1.0471975511966, 1.0471975511966, 1.0471975511966,
    1.0471975511966, 1.0471975511966, 1.0471975511966, 1.0471975511966,
    1.0471975511966, 1.0471975511966, 1.0471975511966, 1.0471975511966,
    1.0471975511966, 1.0471975511966, 1.0471975511966, 1.0471975511966,
    1.0471975511966, 1.0471975511966, 1.0471975511966, 1.0471975511966,
    1.0471975511966, 1.0471975511966, 1.0471975511966, 1.0471975511966,
    1.0471975511966, 1.0471975511966, 1.0471975511966, 1.0471975511966,
    1.0471975511966, 1.0471975511966, 1.0471975511966, 1.0471975511966,
    1.0471975511966, 1.0471975511966, 1.0471975511966, 1.0471975511966,
    1.0471975511966, 1.0471975511966, 1.0471975511966, 1.0471975511966,
    1.0471975511966, 1.0471975511966, 0.698131700797732, 0.698131700797732,
    0.698131700797732, 0.698131700797732, 1.0471975511966, 1.0471975511966,
    1.0471975511966, 1.0471975511966 };

  static const real_T d_a[68]{ 0.0082393983319562145, 0.043241643790057967,
    0.0776020389914959, 0.10308013732709863, 0.11740668642495948,
    0.12024122330596915, 0.11192495647043646, 0.093087329418972228,
    0.064513214663213, 0.0270897107947072, -0.018222606641386398,
    -0.070410680681089924, -0.12842718008483034, -0.191207311425991,
    -0.25768455254736161, -0.32680530976127503, -0.39754242054157474,
    -0.46890741425646021, -0.53996145601447953, -0.60982491755695356,
    -0.677685539447547, -0.74280516892895021, -0.80452507708795107,
    -0.86226987707030656, -0.91555008181220787, -0.96396335495883223,
    -1.0071945222127245, -1.045014422216622, -1.0772776861761704,
    -1.1039195437437468, -0.0082393983319562145, -0.043241643790057967,
    -0.0776020389914959, -0.10308013732709863, -0.11740668642495948,
    -0.12024122330596915, -0.11192495647043646, -0.093087329418972228,
    -0.064513214663213, -0.0270897107947072, 0.018222606641386398,
    0.070410680681089924, 0.12842718008483034, 0.191207311425991,
    0.25768455254736161, 0.32680530976127503, 0.39754242054157474,
    0.46890741425646021, 0.53996145601447953, 0.60982491755695356,
    0.677685539447547, 0.74280516892895021, 0.80452507708795107,
    0.86226987707030656, 0.91555008181220787, 0.96396335495883223,
    1.0071945222127245, 1.045014422216622, 1.0772776861761704,
    1.1039195437437468, -1.0, -1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T b_Kr_0[60]{ 0.82393983319562147, 4.3241643790057971,
    7.76020389914959, 10.308013732709863, 11.740668642495947, 12.024122330596915,
    11.192495647043646, 9.3087329418972224, 6.4513214663213, 2.70897107947072,
    -1.8222606641386399, -7.0410680681089923, -12.842718008483034,
    -19.1207311425991, -25.76845525473616, -32.6805309761275,
    -39.754242054157473, -46.89074142564602, -53.996145601447957,
    -60.982491755695357, -67.7685539447547, -74.280516892895022,
    -80.45250770879511, -86.226987707030659, -91.555008181220785,
    -96.396335495883221, -100.71945222127245, -104.5014422216622,
    -107.72776861761703, -110.39195437437468, -0.0, 0.82393983319562147,
    4.3241643790057971, 7.76020389914959, 10.308013732709863, 11.740668642495947,
    12.024122330596915, 11.192495647043646, 9.3087329418972224, 6.4513214663213,
    2.70897107947072, -1.8222606641386399, -7.0410680681089923,
    -12.842718008483034, -19.1207311425991, -25.76845525473616,
    -32.6805309761275, -39.754242054157473, -46.89074142564602,
    -53.996145601447957, -60.982491755695357, -67.7685539447547,
    -74.280516892895022, -80.45250770879511, -86.226987707030659,
    -91.555008181220785, -96.396335495883221, -100.71945222127245,
    -104.5014422216622, -107.72776861761703 };

  static const real_T e_a[49]{ 0.53832173754583346, -0.76417869920176174,
    1.6472151695740147, -0.71774798606464452, 0.0075073117489730017,
    0.0029746638475713636, 0.0, 0.25050651388471762, 1.6967090417446413,
    -1.942455127894442, 1.3955618254410949, -0.010107032280868518,
    -0.0053114380576942842, 0.0, -0.00036887181808894468, 0.451668951918138,
    -0.4448706826053932, 1.3588465532565464, -0.0061024424164344733,
    -0.0017240376034981994, 0.0, -0.12632901196156529, -0.021758959600111479,
    -0.14940665523964092, 1.3389448432481676, -0.00040454936245862245,
    4.58033188390503E-6, 0.0, -0.11822329454014459, -0.30945768789505285,
    0.9015784913267747, -0.65074625082103688, 0.996956642749649,
    -0.00010756532958362895, 0.0, -0.071268494143188407, -0.18665806084912626,
    0.54376736904748513, -0.39254748471046269, 0.0031684550973532627,
    0.99993528989230329, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };

  static const real_T b_Kx[14]{ -829.518231175404, 47.94324089664179,
    -225.55498137463366, -223.51597047253114, 451.34951472329681,
    286.06727354213672, 1144.2772783647592, -697.46613654971, 17.363774773619518,
    -187.9592489948912, -175.83746671227993, 421.5897941586187,
    267.51253657352936, 1033.8853239903847 };

  static const real_T b_Hinv[9]{ 0.03088533347182373, -0.031968424195833407, 0.0,
    -0.031968424195833407, 0.033712857094144989, 0.0, 0.0, 0.0,
    9.9999999999999974E-6 };

  static const real_T b_Linv[9]{ 0.023897248929567894, -0.17410989335882668, 0.0,
    0.0, 0.18361061269475953, 0.0, 0.0, 0.0, 0.003162277660168379 };

  static const real_T a[7]{ 0.5476, 0.4203, 0.2079, -0.0052, 0.0, 0.0, 1.0 };

  static const real_T b_a[7]{ 0.03650822427816084, 0.039227427187085169,
    -0.031787428834962141, -0.05883768397463926, 0.00021786674125768963,
    0.00034981676075266071, 0.093606463238301366 };

  static const real_T f_a[7]{ -0.24614639772968033, -0.64476824672242616,
    1.8782842568026241, -1.3559934948051797, 0.010943492619575232,
    0.0067766169275252208, 0.0 };

  static const real_T g_a[7]{ 0.036873840893038634, 0.025448858519931443,
    0.0072580633878512671, -0.093713151868259681, 0.00031370152338443654,
    0.00030454956090698, 0.093606463238301407 };

  real_T a__1[68];
  real_T b_Mlim[68];
  real_T rseq[30];
  real_T rtb_xest[7];
  real_T xk[7];
  real_T f[3];
  real_T zopt[3];
  real_T b_Kr;
  real_T rtb_last_mv;
  real_T xk_0;
  real_T y_innov;
  int32_T i;
  int32_T i_0;

  // UnitDelay: '<S2>/last_mv'
  rtb_last_mv = rtDW.last_mv_DSTATE;

  // MATLAB Function: '<S22>/optimizer' incorporates:
  //   Inport: '<Root>/In1'
  //   Inport: '<Root>/In2'
  //   Memory: '<S2>/last_x'

  for (i = 0; i < 30; i++) {
    rseq[i] = rtU.r_cmd;
  }

  y_innov = 0.0;
  for (i = 0; i < 7; i++) {
    xk_0 = rtDW.last_x_PreviousInput[i];
    xk[i] = xk_0;
    y_innov += a[i] * xk_0;
  }

  y_innov = rtU.uav_r - y_innov;
  for (i = 0; i < 7; i++) {
    rtb_xest[i] = b_a[i] * y_innov + xk[i];
  }

  f[0] = 0.0;
  f[1] = 0.0;
  f[2] = 0.0;
  for (i = 0; i < 2; i++) {
    xk_0 = 0.0;
    for (i_0 = 0; i_0 < 7; i_0++) {
      xk_0 += b_Kx[7 * i + i_0] * rtb_xest[i_0];
    }

    b_Kr = 0.0;
    for (i_0 = 0; i_0 < 30; i_0++) {
      b_Kr += b_Kr_0[30 * i + i_0] * rseq[i_0];
    }

    f[i] = (-89.60723839903153 * static_cast<real_T>(i) + 1750.072731314694) *
      rtb_last_mv + (xk_0 + b_Kr);
  }

  for (i_0 = 0; i_0 < 68; i_0++) {
    xk_0 = 0.0;
    for (i = 0; i < 7; i++) {
      xk_0 += c_a[68 * i + i_0] * rtb_xest[i];
    }

    b_Mlim[i_0] = -((b_Mlim_0[i_0] + xk_0) + d_a[i_0] * rtb_last_mv);
  }

  // Update for Memory: '<S2>/Memory' incorporates:
  //   MATLAB Function: '<S22>/optimizer'

  qpkwik(b_Linv, b_Hinv, f, b_Ac, b_Mlim, rtDW.Memory_PreviousInput, 284, 1.0E-6,
         zopt, a__1, &i);

  // MATLAB Function: '<S22>/optimizer'
  if ((i < 0) || (i == 0)) {
    zopt[0] = 0.0;
  }

  rtb_last_mv += zopt[0];

  // Outport: '<Root>/Out1' incorporates:
  //   MATLAB Function: '<S22>/optimizer'

  rtY.Out1 = rtb_last_mv;

  // Update for Memory: '<S2>/last_x' incorporates:
  //   MATLAB Function: '<S22>/optimizer'

  for (i_0 = 0; i_0 < 7; i_0++) {
    xk_0 = 0.0;
    for (i = 0; i < 7; i++) {
      xk_0 += e_a[7 * i + i_0] * xk[i];
    }

    rtDW.last_x_PreviousInput[i_0] = (f_a[i_0] * rtb_last_mv + xk_0) + g_a[i_0] *
      y_innov;
  }

  // End of Update for Memory: '<S2>/last_x'

  // Update for UnitDelay: '<S2>/last_mv' incorporates:
  //   MATLAB Function: '<S22>/optimizer'

  rtDW.last_mv_DSTATE = rtb_last_mv;
}

// Model initialize function
void thesis_mpc::initialize()
{
  // Registration code

  // initialize non-finites
  rt_InitInfAndNaN(sizeof(real_T));
}

// Constructor
thesis_mpc::thesis_mpc() :
  rtU(),
  rtY(),
  rtDW(),
  rtM()
{
  // Currently there is no constructor body generated.
}

// Destructor
// Currently there is no destructor body generated.
thesis_mpc::~thesis_mpc() = default;

// Real-Time Model get method
thesis_mpc::RT_MODEL * thesis_mpc::getRTM()
{
  return (&rtM);
}

//
// File trailer for generated code.
//
// [EOF]
//
