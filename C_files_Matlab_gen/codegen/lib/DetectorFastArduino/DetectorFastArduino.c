/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * DetectorFastArduino.c
 *
 * Code generation for function 'DetectorFastArduino'
 *
 */

/* Include files */
#include "DetectorFastArduino.h"
#include "DetectorFastArduino_emxutil.h"
#include "DetectorFastArduino_types.h"
#include "colon.h"
#include "conv.h"
#include "find.h"
#include "fliplr.h"
#include "linsolve.h"
#include "mean.h"
#include "minOrMax.h"
#include "polyfit.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static void b_binary_expand_op(emxArray_real_T *Cir, double maxval,
                               const emxArray_real_T *Gir, double b_maxval,
                               const emxArray_real_T *Cfilter, double c_maxval,
                               const emxArray_real_T *Dfilter, double d_maxval,
                               const emxArray_real_T *Efilter, double e_maxval,
                               const emxArray_real_T *Cir2, double f_maxval);

static void binary_expand_op(emxArray_real_T *Cfilter,
                             const emxArray_real_T *Gir, int i1, int i2, int i3,
                             int i4);

static void c_binary_expand_op(emxArray_real_T *Cir2,
                               const double InterdependencyMatrix[36],
                               const emxArray_real_T *Dir2,
                               const emxArray_real_T *Eir2,
                               const emxArray_real_T *Fir2,
                               const emxArray_real_T *timeSig,
                               const emxArray_real_T *Ffilter);

static void d_binary_expand_op(
    emxArray_real_T *Efilter, const double InterdependencyMatrix[36],
    const emxArray_real_T *Cir2, const emxArray_real_T *Dir2,
    const emxArray_real_T *Eir2, const emxArray_real_T *Fir2,
    const emxArray_real_T *timeSig, const emxArray_real_T *Ffilter);

static int div_s32_floor(int numerator, int denominator);

static void e_binary_expand_op(
    emxArray_real_T *Dfilter, const double InterdependencyMatrix[36],
    const emxArray_real_T *Cir2, const emxArray_real_T *Dir2,
    const emxArray_real_T *Eir2, const emxArray_real_T *Fir2,
    const emxArray_real_T *timeSig, const emxArray_real_T *Ffilter);

static void f_binary_expand_op(
    emxArray_real_T *Cfilter, const double InterdependencyMatrix[36],
    const emxArray_real_T *Cir2, const emxArray_real_T *Dir2,
    const emxArray_real_T *Eir2, const emxArray_real_T *Fir2,
    const emxArray_real_T *timeSig, const emxArray_real_T *Ffilter);

static void
g_binary_expand_op(emxArray_real_T *Gir, const double InterdependencyMatrix[36],
                   const emxArray_real_T *Cir2, const emxArray_real_T *Dir2,
                   const emxArray_real_T *Eir2, const emxArray_real_T *Fir2,
                   const emxArray_real_T *timeSig,
                   const emxArray_real_T *Ffilter);

static void
h_binary_expand_op(emxArray_real_T *Cir, const double InterdependencyMatrix[36],
                   const emxArray_real_T *Cir2, const emxArray_real_T *Dir2,
                   const emxArray_real_T *Eir2, const emxArray_real_T *Fir2,
                   const emxArray_real_T *timeSig,
                   const emxArray_real_T *Ffilter);

static double rt_powd_snf(double u0, double u1);

static double rt_roundd_snf(double u);

/* Function Definitions */
static void b_binary_expand_op(emxArray_real_T *Cir, double maxval,
                               const emxArray_real_T *Gir, double b_maxval,
                               const emxArray_real_T *Cfilter, double c_maxval,
                               const emxArray_real_T *Dfilter, double d_maxval,
                               const emxArray_real_T *Efilter, double e_maxval,
                               const emxArray_real_T *Cir2, double f_maxval)
{
  emxArray_real_T *b_Cir;
  const double *Cfilter_data;
  const double *Cir2_data;
  const double *Dfilter_data;
  const double *Efilter_data;
  const double *Gir_data;
  double d;
  double d1;
  double d2;
  double d3;
  double d4;
  double d5;
  double *Cir_data;
  double *b_Cir_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  int stride_2_1;
  int stride_3_1;
  int stride_4_1;
  int stride_5_1;
  Cir2_data = Cir2->data;
  Efilter_data = Efilter->data;
  Dfilter_data = Dfilter->data;
  Cfilter_data = Cfilter->data;
  Gir_data = Gir->data;
  Cir_data = Cir->data;
  emxInit_real_T(&b_Cir, 2);
  d = 0.5 * maxval;
  d1 = 0.5 * b_maxval;
  d2 = 0.5 * c_maxval;
  d3 = 0.5 * d_maxval;
  d4 = 0.5 * e_maxval;
  d5 = 0.5 * f_maxval;
  i = b_Cir->size[0] * b_Cir->size[1];
  b_Cir->size[0] = 1;
  if (Cir2->size[1] == 1) {
    if (Efilter->size[1] == 1) {
      if (Dfilter->size[1] == 1) {
        if (Cfilter->size[1] == 1) {
          if (Gir->size[1] == 1) {
            b_Cir->size[1] = Cir->size[1];
          } else {
            b_Cir->size[1] = Gir->size[1];
          }
        } else {
          b_Cir->size[1] = Cfilter->size[1];
        }
      } else {
        b_Cir->size[1] = Dfilter->size[1];
      }
    } else {
      b_Cir->size[1] = Efilter->size[1];
    }
  } else {
    b_Cir->size[1] = Cir2->size[1];
  }
  emxEnsureCapacity_real_T(b_Cir, i);
  b_Cir_data = b_Cir->data;
  stride_0_1 = (Cir->size[1] != 1);
  stride_1_1 = (Gir->size[1] != 1);
  stride_2_1 = (Cfilter->size[1] != 1);
  stride_3_1 = (Dfilter->size[1] != 1);
  stride_4_1 = (Efilter->size[1] != 1);
  stride_5_1 = (Cir2->size[1] != 1);
  if (Cir2->size[1] == 1) {
    if (Efilter->size[1] == 1) {
      if (Dfilter->size[1] == 1) {
        if (Cfilter->size[1] == 1) {
          if (Gir->size[1] == 1) {
            loop_ub = Cir->size[1];
          } else {
            loop_ub = Gir->size[1];
          }
        } else {
          loop_ub = Cfilter->size[1];
        }
      } else {
        loop_ub = Dfilter->size[1];
      }
    } else {
      loop_ub = Efilter->size[1];
    }
  } else {
    loop_ub = Cir2->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    b_Cir_data[i] =
        (((((Cir_data[i * stride_0_1] > d) + (Gir_data[i * stride_1_1] > d1)) +
           (Cfilter_data[i * stride_2_1] > d2)) +
          (Dfilter_data[i * stride_3_1] > d3)) +
         (Efilter_data[i * stride_4_1] > d4)) +
        (Cir2_data[i * stride_5_1] > d5);
  }
  c_conv(b_Cir, Cir);
  emxFree_real_T(&b_Cir);
}

static void binary_expand_op(emxArray_real_T *Cfilter,
                             const emxArray_real_T *Gir, int i1, int i2, int i3,
                             int i4)
{
  const double *Gir_data;
  double *Cfilter_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  Gir_data = Gir->data;
  i = Cfilter->size[0] * Cfilter->size[1];
  Cfilter->size[0] = 1;
  if (i4 + 1 == 1) {
    Cfilter->size[1] = (i3 - i2) + 1;
  } else {
    Cfilter->size[1] = i4 + 1;
  }
  emxEnsureCapacity_real_T(Cfilter, i);
  Cfilter_data = Cfilter->data;
  stride_0_1 = ((i3 - i2) + 1 != 1);
  stride_1_1 = (i4 + 1 != 1);
  if (i4 + 1 == 1) {
    loop_ub = (i3 - i2) + 1;
  } else {
    loop_ub = i4 + 1;
  }
  for (i = 0; i < loop_ub; i++) {
    Cfilter_data[i] = (int)(unsigned int)Gir_data[i1 * (i2 + i * stride_0_1)] -
                      (int)(unsigned int)Gir_data[i1 * (i * stride_1_1)];
  }
}

static void c_binary_expand_op(emxArray_real_T *Cir2,
                               const double InterdependencyMatrix[36],
                               const emxArray_real_T *Dir2,
                               const emxArray_real_T *Eir2,
                               const emxArray_real_T *Fir2,
                               const emxArray_real_T *timeSig,
                               const emxArray_real_T *Ffilter)
{
  emxArray_real_T *b_InterdependencyMatrix;
  const double *Dir2_data;
  const double *Eir2_data;
  const double *Ffilter_data;
  const double *Fir2_data;
  const double *timeSig_data;
  double c_InterdependencyMatrix;
  double d_InterdependencyMatrix;
  double e_InterdependencyMatrix;
  double f_InterdependencyMatrix;
  double g_InterdependencyMatrix;
  double h_InterdependencyMatrix;
  double *Cir2_data;
  double *InterdependencyMatrix_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  int stride_2_1;
  int stride_3_1;
  int stride_4_1;
  int stride_5_1;
  Ffilter_data = Ffilter->data;
  timeSig_data = timeSig->data;
  Fir2_data = Fir2->data;
  Eir2_data = Eir2->data;
  Dir2_data = Dir2->data;
  Cir2_data = Cir2->data;
  emxInit_real_T(&b_InterdependencyMatrix, 2);
  c_InterdependencyMatrix = InterdependencyMatrix[5];
  d_InterdependencyMatrix = InterdependencyMatrix[11];
  e_InterdependencyMatrix = InterdependencyMatrix[17];
  f_InterdependencyMatrix = InterdependencyMatrix[23];
  g_InterdependencyMatrix = InterdependencyMatrix[29];
  h_InterdependencyMatrix = InterdependencyMatrix[35];
  i = b_InterdependencyMatrix->size[0] * b_InterdependencyMatrix->size[1];
  b_InterdependencyMatrix->size[0] = 1;
  if (Ffilter->size[1] == 1) {
    if (timeSig->size[1] == 1) {
      if (Fir2->size[1] == 1) {
        if (Eir2->size[1] == 1) {
          if (Dir2->size[1] == 1) {
            b_InterdependencyMatrix->size[1] = Cir2->size[1];
          } else {
            b_InterdependencyMatrix->size[1] = Dir2->size[1];
          }
        } else {
          b_InterdependencyMatrix->size[1] = Eir2->size[1];
        }
      } else {
        b_InterdependencyMatrix->size[1] = Fir2->size[1];
      }
    } else {
      b_InterdependencyMatrix->size[1] = timeSig->size[1];
    }
  } else {
    b_InterdependencyMatrix->size[1] = Ffilter->size[1];
  }
  emxEnsureCapacity_real_T(b_InterdependencyMatrix, i);
  InterdependencyMatrix_data = b_InterdependencyMatrix->data;
  stride_0_1 = (Cir2->size[1] != 1);
  stride_1_1 = (Dir2->size[1] != 1);
  stride_2_1 = (Eir2->size[1] != 1);
  stride_3_1 = (Fir2->size[1] != 1);
  stride_4_1 = (timeSig->size[1] != 1);
  stride_5_1 = (Ffilter->size[1] != 1);
  if (Ffilter->size[1] == 1) {
    if (timeSig->size[1] == 1) {
      if (Fir2->size[1] == 1) {
        if (Eir2->size[1] == 1) {
          if (Dir2->size[1] == 1) {
            loop_ub = Cir2->size[1];
          } else {
            loop_ub = Dir2->size[1];
          }
        } else {
          loop_ub = Eir2->size[1];
        }
      } else {
        loop_ub = Fir2->size[1];
      }
    } else {
      loop_ub = timeSig->size[1];
    }
  } else {
    loop_ub = Ffilter->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    InterdependencyMatrix_data[i] =
        ((((c_InterdependencyMatrix * Cir2_data[i * stride_0_1] +
            d_InterdependencyMatrix * Dir2_data[i * stride_1_1]) +
           e_InterdependencyMatrix * Eir2_data[i * stride_2_1]) +
          f_InterdependencyMatrix * Fir2_data[i * stride_3_1]) +
         g_InterdependencyMatrix * timeSig_data[i * stride_4_1]) +
        h_InterdependencyMatrix * Ffilter_data[i * stride_5_1];
  }
  i = Cir2->size[0] * Cir2->size[1];
  Cir2->size[0] = 1;
  Cir2->size[1] = b_InterdependencyMatrix->size[1];
  emxEnsureCapacity_real_T(Cir2, i);
  Cir2_data = Cir2->data;
  loop_ub = b_InterdependencyMatrix->size[1];
  for (i = 0; i < loop_ub; i++) {
    Cir2_data[i] = InterdependencyMatrix_data[i];
  }
  emxFree_real_T(&b_InterdependencyMatrix);
}

static void d_binary_expand_op(
    emxArray_real_T *Efilter, const double InterdependencyMatrix[36],
    const emxArray_real_T *Cir2, const emxArray_real_T *Dir2,
    const emxArray_real_T *Eir2, const emxArray_real_T *Fir2,
    const emxArray_real_T *timeSig, const emxArray_real_T *Ffilter)
{
  const double *Cir2_data;
  const double *Dir2_data;
  const double *Eir2_data;
  const double *Ffilter_data;
  const double *Fir2_data;
  const double *timeSig_data;
  double b_InterdependencyMatrix;
  double c_InterdependencyMatrix;
  double d_InterdependencyMatrix;
  double e_InterdependencyMatrix;
  double f_InterdependencyMatrix;
  double g_InterdependencyMatrix;
  double *Efilter_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  int stride_2_1;
  int stride_3_1;
  int stride_4_1;
  int stride_5_1;
  Ffilter_data = Ffilter->data;
  timeSig_data = timeSig->data;
  Fir2_data = Fir2->data;
  Eir2_data = Eir2->data;
  Dir2_data = Dir2->data;
  Cir2_data = Cir2->data;
  b_InterdependencyMatrix = InterdependencyMatrix[4];
  c_InterdependencyMatrix = InterdependencyMatrix[10];
  d_InterdependencyMatrix = InterdependencyMatrix[16];
  e_InterdependencyMatrix = InterdependencyMatrix[22];
  f_InterdependencyMatrix = InterdependencyMatrix[28];
  g_InterdependencyMatrix = InterdependencyMatrix[34];
  i = Efilter->size[0] * Efilter->size[1];
  Efilter->size[0] = 1;
  if (Ffilter->size[1] == 1) {
    if (timeSig->size[1] == 1) {
      if (Fir2->size[1] == 1) {
        if (Eir2->size[1] == 1) {
          if (Dir2->size[1] == 1) {
            Efilter->size[1] = Cir2->size[1];
          } else {
            Efilter->size[1] = Dir2->size[1];
          }
        } else {
          Efilter->size[1] = Eir2->size[1];
        }
      } else {
        Efilter->size[1] = Fir2->size[1];
      }
    } else {
      Efilter->size[1] = timeSig->size[1];
    }
  } else {
    Efilter->size[1] = Ffilter->size[1];
  }
  emxEnsureCapacity_real_T(Efilter, i);
  Efilter_data = Efilter->data;
  stride_0_1 = (Cir2->size[1] != 1);
  stride_1_1 = (Dir2->size[1] != 1);
  stride_2_1 = (Eir2->size[1] != 1);
  stride_3_1 = (Fir2->size[1] != 1);
  stride_4_1 = (timeSig->size[1] != 1);
  stride_5_1 = (Ffilter->size[1] != 1);
  if (Ffilter->size[1] == 1) {
    if (timeSig->size[1] == 1) {
      if (Fir2->size[1] == 1) {
        if (Eir2->size[1] == 1) {
          if (Dir2->size[1] == 1) {
            loop_ub = Cir2->size[1];
          } else {
            loop_ub = Dir2->size[1];
          }
        } else {
          loop_ub = Eir2->size[1];
        }
      } else {
        loop_ub = Fir2->size[1];
      }
    } else {
      loop_ub = timeSig->size[1];
    }
  } else {
    loop_ub = Ffilter->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    Efilter_data[i] = ((((b_InterdependencyMatrix * Cir2_data[i * stride_0_1] +
                          c_InterdependencyMatrix * Dir2_data[i * stride_1_1]) +
                         d_InterdependencyMatrix * Eir2_data[i * stride_2_1]) +
                        e_InterdependencyMatrix * Fir2_data[i * stride_3_1]) +
                       f_InterdependencyMatrix * timeSig_data[i * stride_4_1]) +
                      g_InterdependencyMatrix * Ffilter_data[i * stride_5_1];
  }
}

static int div_s32_floor(int numerator, int denominator)
{
  unsigned int absDenominator;
  unsigned int absNumerator;
  int quotient;
  unsigned int tempAbsQuotient;
  boolean_T quotientNeedsNegation;
  if (denominator == 0) {
    if (numerator >= 0) {
      quotient = MAX_int32_T;
    } else {
      quotient = MIN_int32_T;
    }
  } else {
    if (numerator < 0) {
      absNumerator = ~(unsigned int)numerator + 1U;
    } else {
      absNumerator = (unsigned int)numerator;
    }
    if (denominator < 0) {
      absDenominator = ~(unsigned int)denominator + 1U;
    } else {
      absDenominator = (unsigned int)denominator;
    }
    quotientNeedsNegation = ((numerator < 0) != (denominator < 0));
    tempAbsQuotient = absNumerator / absDenominator;
    if (quotientNeedsNegation) {
      absNumerator %= absDenominator;
      if (absNumerator > 0U) {
        tempAbsQuotient++;
      }
      quotient = -(int)tempAbsQuotient;
    } else {
      quotient = (int)tempAbsQuotient;
    }
  }
  return quotient;
}

static void e_binary_expand_op(
    emxArray_real_T *Dfilter, const double InterdependencyMatrix[36],
    const emxArray_real_T *Cir2, const emxArray_real_T *Dir2,
    const emxArray_real_T *Eir2, const emxArray_real_T *Fir2,
    const emxArray_real_T *timeSig, const emxArray_real_T *Ffilter)
{
  const double *Cir2_data;
  const double *Dir2_data;
  const double *Eir2_data;
  const double *Ffilter_data;
  const double *Fir2_data;
  const double *timeSig_data;
  double b_InterdependencyMatrix;
  double c_InterdependencyMatrix;
  double d_InterdependencyMatrix;
  double e_InterdependencyMatrix;
  double f_InterdependencyMatrix;
  double g_InterdependencyMatrix;
  double *Dfilter_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  int stride_2_1;
  int stride_3_1;
  int stride_4_1;
  int stride_5_1;
  Ffilter_data = Ffilter->data;
  timeSig_data = timeSig->data;
  Fir2_data = Fir2->data;
  Eir2_data = Eir2->data;
  Dir2_data = Dir2->data;
  Cir2_data = Cir2->data;
  b_InterdependencyMatrix = InterdependencyMatrix[3];
  c_InterdependencyMatrix = InterdependencyMatrix[9];
  d_InterdependencyMatrix = InterdependencyMatrix[15];
  e_InterdependencyMatrix = InterdependencyMatrix[21];
  f_InterdependencyMatrix = InterdependencyMatrix[27];
  g_InterdependencyMatrix = InterdependencyMatrix[33];
  i = Dfilter->size[0] * Dfilter->size[1];
  Dfilter->size[0] = 1;
  if (Ffilter->size[1] == 1) {
    if (timeSig->size[1] == 1) {
      if (Fir2->size[1] == 1) {
        if (Eir2->size[1] == 1) {
          if (Dir2->size[1] == 1) {
            Dfilter->size[1] = Cir2->size[1];
          } else {
            Dfilter->size[1] = Dir2->size[1];
          }
        } else {
          Dfilter->size[1] = Eir2->size[1];
        }
      } else {
        Dfilter->size[1] = Fir2->size[1];
      }
    } else {
      Dfilter->size[1] = timeSig->size[1];
    }
  } else {
    Dfilter->size[1] = Ffilter->size[1];
  }
  emxEnsureCapacity_real_T(Dfilter, i);
  Dfilter_data = Dfilter->data;
  stride_0_1 = (Cir2->size[1] != 1);
  stride_1_1 = (Dir2->size[1] != 1);
  stride_2_1 = (Eir2->size[1] != 1);
  stride_3_1 = (Fir2->size[1] != 1);
  stride_4_1 = (timeSig->size[1] != 1);
  stride_5_1 = (Ffilter->size[1] != 1);
  if (Ffilter->size[1] == 1) {
    if (timeSig->size[1] == 1) {
      if (Fir2->size[1] == 1) {
        if (Eir2->size[1] == 1) {
          if (Dir2->size[1] == 1) {
            loop_ub = Cir2->size[1];
          } else {
            loop_ub = Dir2->size[1];
          }
        } else {
          loop_ub = Eir2->size[1];
        }
      } else {
        loop_ub = Fir2->size[1];
      }
    } else {
      loop_ub = timeSig->size[1];
    }
  } else {
    loop_ub = Ffilter->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    Dfilter_data[i] = ((((b_InterdependencyMatrix * Cir2_data[i * stride_0_1] +
                          c_InterdependencyMatrix * Dir2_data[i * stride_1_1]) +
                         d_InterdependencyMatrix * Eir2_data[i * stride_2_1]) +
                        e_InterdependencyMatrix * Fir2_data[i * stride_3_1]) +
                       f_InterdependencyMatrix * timeSig_data[i * stride_4_1]) +
                      g_InterdependencyMatrix * Ffilter_data[i * stride_5_1];
  }
}

static void f_binary_expand_op(
    emxArray_real_T *Cfilter, const double InterdependencyMatrix[36],
    const emxArray_real_T *Cir2, const emxArray_real_T *Dir2,
    const emxArray_real_T *Eir2, const emxArray_real_T *Fir2,
    const emxArray_real_T *timeSig, const emxArray_real_T *Ffilter)
{
  const double *Cir2_data;
  const double *Dir2_data;
  const double *Eir2_data;
  const double *Ffilter_data;
  const double *Fir2_data;
  const double *timeSig_data;
  double b_InterdependencyMatrix;
  double c_InterdependencyMatrix;
  double d_InterdependencyMatrix;
  double e_InterdependencyMatrix;
  double f_InterdependencyMatrix;
  double g_InterdependencyMatrix;
  double *Cfilter_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  int stride_2_1;
  int stride_3_1;
  int stride_4_1;
  int stride_5_1;
  Ffilter_data = Ffilter->data;
  timeSig_data = timeSig->data;
  Fir2_data = Fir2->data;
  Eir2_data = Eir2->data;
  Dir2_data = Dir2->data;
  Cir2_data = Cir2->data;
  b_InterdependencyMatrix = InterdependencyMatrix[2];
  c_InterdependencyMatrix = InterdependencyMatrix[8];
  d_InterdependencyMatrix = InterdependencyMatrix[14];
  e_InterdependencyMatrix = InterdependencyMatrix[20];
  f_InterdependencyMatrix = InterdependencyMatrix[26];
  g_InterdependencyMatrix = InterdependencyMatrix[32];
  i = Cfilter->size[0] * Cfilter->size[1];
  Cfilter->size[0] = 1;
  if (Ffilter->size[1] == 1) {
    if (timeSig->size[1] == 1) {
      if (Fir2->size[1] == 1) {
        if (Eir2->size[1] == 1) {
          if (Dir2->size[1] == 1) {
            Cfilter->size[1] = Cir2->size[1];
          } else {
            Cfilter->size[1] = Dir2->size[1];
          }
        } else {
          Cfilter->size[1] = Eir2->size[1];
        }
      } else {
        Cfilter->size[1] = Fir2->size[1];
      }
    } else {
      Cfilter->size[1] = timeSig->size[1];
    }
  } else {
    Cfilter->size[1] = Ffilter->size[1];
  }
  emxEnsureCapacity_real_T(Cfilter, i);
  Cfilter_data = Cfilter->data;
  stride_0_1 = (Cir2->size[1] != 1);
  stride_1_1 = (Dir2->size[1] != 1);
  stride_2_1 = (Eir2->size[1] != 1);
  stride_3_1 = (Fir2->size[1] != 1);
  stride_4_1 = (timeSig->size[1] != 1);
  stride_5_1 = (Ffilter->size[1] != 1);
  if (Ffilter->size[1] == 1) {
    if (timeSig->size[1] == 1) {
      if (Fir2->size[1] == 1) {
        if (Eir2->size[1] == 1) {
          if (Dir2->size[1] == 1) {
            loop_ub = Cir2->size[1];
          } else {
            loop_ub = Dir2->size[1];
          }
        } else {
          loop_ub = Eir2->size[1];
        }
      } else {
        loop_ub = Fir2->size[1];
      }
    } else {
      loop_ub = timeSig->size[1];
    }
  } else {
    loop_ub = Ffilter->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    Cfilter_data[i] = ((((b_InterdependencyMatrix * Cir2_data[i * stride_0_1] +
                          c_InterdependencyMatrix * Dir2_data[i * stride_1_1]) +
                         d_InterdependencyMatrix * Eir2_data[i * stride_2_1]) +
                        e_InterdependencyMatrix * Fir2_data[i * stride_3_1]) +
                       f_InterdependencyMatrix * timeSig_data[i * stride_4_1]) +
                      g_InterdependencyMatrix * Ffilter_data[i * stride_5_1];
  }
}

static void
g_binary_expand_op(emxArray_real_T *Gir, const double InterdependencyMatrix[36],
                   const emxArray_real_T *Cir2, const emxArray_real_T *Dir2,
                   const emxArray_real_T *Eir2, const emxArray_real_T *Fir2,
                   const emxArray_real_T *timeSig,
                   const emxArray_real_T *Ffilter)
{
  const double *Cir2_data;
  const double *Dir2_data;
  const double *Eir2_data;
  const double *Ffilter_data;
  const double *Fir2_data;
  const double *timeSig_data;
  double b_InterdependencyMatrix;
  double c_InterdependencyMatrix;
  double d_InterdependencyMatrix;
  double e_InterdependencyMatrix;
  double f_InterdependencyMatrix;
  double g_InterdependencyMatrix;
  double *Gir_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  int stride_2_1;
  int stride_3_1;
  int stride_4_1;
  int stride_5_1;
  Ffilter_data = Ffilter->data;
  timeSig_data = timeSig->data;
  Fir2_data = Fir2->data;
  Eir2_data = Eir2->data;
  Dir2_data = Dir2->data;
  Cir2_data = Cir2->data;
  b_InterdependencyMatrix = InterdependencyMatrix[1];
  c_InterdependencyMatrix = InterdependencyMatrix[7];
  d_InterdependencyMatrix = InterdependencyMatrix[13];
  e_InterdependencyMatrix = InterdependencyMatrix[19];
  f_InterdependencyMatrix = InterdependencyMatrix[25];
  g_InterdependencyMatrix = InterdependencyMatrix[31];
  i = Gir->size[0] * Gir->size[1];
  Gir->size[0] = 1;
  if (Ffilter->size[1] == 1) {
    if (timeSig->size[1] == 1) {
      if (Fir2->size[1] == 1) {
        if (Eir2->size[1] == 1) {
          if (Dir2->size[1] == 1) {
            Gir->size[1] = Cir2->size[1];
          } else {
            Gir->size[1] = Dir2->size[1];
          }
        } else {
          Gir->size[1] = Eir2->size[1];
        }
      } else {
        Gir->size[1] = Fir2->size[1];
      }
    } else {
      Gir->size[1] = timeSig->size[1];
    }
  } else {
    Gir->size[1] = Ffilter->size[1];
  }
  emxEnsureCapacity_real_T(Gir, i);
  Gir_data = Gir->data;
  stride_0_1 = (Cir2->size[1] != 1);
  stride_1_1 = (Dir2->size[1] != 1);
  stride_2_1 = (Eir2->size[1] != 1);
  stride_3_1 = (Fir2->size[1] != 1);
  stride_4_1 = (timeSig->size[1] != 1);
  stride_5_1 = (Ffilter->size[1] != 1);
  if (Ffilter->size[1] == 1) {
    if (timeSig->size[1] == 1) {
      if (Fir2->size[1] == 1) {
        if (Eir2->size[1] == 1) {
          if (Dir2->size[1] == 1) {
            loop_ub = Cir2->size[1];
          } else {
            loop_ub = Dir2->size[1];
          }
        } else {
          loop_ub = Eir2->size[1];
        }
      } else {
        loop_ub = Fir2->size[1];
      }
    } else {
      loop_ub = timeSig->size[1];
    }
  } else {
    loop_ub = Ffilter->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    Gir_data[i] = ((((b_InterdependencyMatrix * Cir2_data[i * stride_0_1] +
                      c_InterdependencyMatrix * Dir2_data[i * stride_1_1]) +
                     d_InterdependencyMatrix * Eir2_data[i * stride_2_1]) +
                    e_InterdependencyMatrix * Fir2_data[i * stride_3_1]) +
                   f_InterdependencyMatrix * timeSig_data[i * stride_4_1]) +
                  g_InterdependencyMatrix * Ffilter_data[i * stride_5_1];
  }
}

static void
h_binary_expand_op(emxArray_real_T *Cir, const double InterdependencyMatrix[36],
                   const emxArray_real_T *Cir2, const emxArray_real_T *Dir2,
                   const emxArray_real_T *Eir2, const emxArray_real_T *Fir2,
                   const emxArray_real_T *timeSig,
                   const emxArray_real_T *Ffilter)
{
  const double *Cir2_data;
  const double *Dir2_data;
  const double *Eir2_data;
  const double *Ffilter_data;
  const double *Fir2_data;
  const double *timeSig_data;
  double b_InterdependencyMatrix;
  double c_InterdependencyMatrix;
  double d_InterdependencyMatrix;
  double e_InterdependencyMatrix;
  double f_InterdependencyMatrix;
  double g_InterdependencyMatrix;
  double *Cir_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  int stride_2_1;
  int stride_3_1;
  int stride_4_1;
  int stride_5_1;
  Ffilter_data = Ffilter->data;
  timeSig_data = timeSig->data;
  Fir2_data = Fir2->data;
  Eir2_data = Eir2->data;
  Dir2_data = Dir2->data;
  Cir2_data = Cir2->data;
  b_InterdependencyMatrix = InterdependencyMatrix[0];
  c_InterdependencyMatrix = InterdependencyMatrix[6];
  d_InterdependencyMatrix = InterdependencyMatrix[12];
  e_InterdependencyMatrix = InterdependencyMatrix[18];
  f_InterdependencyMatrix = InterdependencyMatrix[24];
  g_InterdependencyMatrix = InterdependencyMatrix[30];
  i = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  if (Ffilter->size[1] == 1) {
    if (timeSig->size[1] == 1) {
      if (Fir2->size[1] == 1) {
        if (Eir2->size[1] == 1) {
          if (Dir2->size[1] == 1) {
            Cir->size[1] = Cir2->size[1];
          } else {
            Cir->size[1] = Dir2->size[1];
          }
        } else {
          Cir->size[1] = Eir2->size[1];
        }
      } else {
        Cir->size[1] = Fir2->size[1];
      }
    } else {
      Cir->size[1] = timeSig->size[1];
    }
  } else {
    Cir->size[1] = Ffilter->size[1];
  }
  emxEnsureCapacity_real_T(Cir, i);
  Cir_data = Cir->data;
  stride_0_1 = (Cir2->size[1] != 1);
  stride_1_1 = (Dir2->size[1] != 1);
  stride_2_1 = (Eir2->size[1] != 1);
  stride_3_1 = (Fir2->size[1] != 1);
  stride_4_1 = (timeSig->size[1] != 1);
  stride_5_1 = (Ffilter->size[1] != 1);
  if (Ffilter->size[1] == 1) {
    if (timeSig->size[1] == 1) {
      if (Fir2->size[1] == 1) {
        if (Eir2->size[1] == 1) {
          if (Dir2->size[1] == 1) {
            loop_ub = Cir2->size[1];
          } else {
            loop_ub = Dir2->size[1];
          }
        } else {
          loop_ub = Eir2->size[1];
        }
      } else {
        loop_ub = Fir2->size[1];
      }
    } else {
      loop_ub = timeSig->size[1];
    }
  } else {
    loop_ub = Ffilter->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    Cir_data[i] = ((((b_InterdependencyMatrix * Cir2_data[i * stride_0_1] +
                      c_InterdependencyMatrix * Dir2_data[i * stride_1_1]) +
                     d_InterdependencyMatrix * Eir2_data[i * stride_2_1]) +
                    e_InterdependencyMatrix * Fir2_data[i * stride_3_1]) +
                   f_InterdependencyMatrix * timeSig_data[i * stride_4_1]) +
                  g_InterdependencyMatrix * Ffilter_data[i * stride_5_1];
  }
}

static double rt_powd_snf(double u0, double u1)
{
  double d;
  double d1;
  double y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d = fabs(u0);
    d1 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d == 1.0) {
        y = 1.0;
      } else if (d > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }
  return y;
}

static double rt_roundd_snf(double u)
{
  double y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }
  return y;
}

double DetectorFastArduino(const double digital[106501], double Octive,
                           double fs_Hz)
{
  emxArray_boolean_T *b_ri;
  emxArray_int32_T *r;
  emxArray_int32_T *r2;
  emxArray_int32_T *r3;
  emxArray_real_T c_Cir_data;
  emxArray_real_T c_peakHistory_data;
  emxArray_real_T *Cfilter;
  emxArray_real_T *Cir;
  emxArray_real_T *Dfilter;
  emxArray_real_T *Dir2;
  emxArray_real_T *Efilter;
  emxArray_real_T *Eir2;
  emxArray_real_T *Ffilter;
  emxArray_real_T *Fir2;
  emxArray_real_T *Gfilter;
  emxArray_real_T *Gir;
  emxArray_real_T *c_Cir2;
  emxArray_real_T *ri;
  emxArray_real_T *timeSig;
  double CC[36];
  double InterdependencyMatrix[36];
  double b_tempo[30];
  double tempo_data[30];
  double b_Cir_data[29];
  double b_peakHistory_data[29];
  double dv[2];
  double C;
  double D;
  double E;
  double F;
  double G;
  double b;
  double d;
  double factor;
  double high_C;
  double tempo;
  double *Cfilter_data;
  double *Cir2_data;
  double *Cir_data;
  double *Dfilter_data;
  double *Dir2_data;
  double *Efilter_data;
  double *Eir2_data;
  double *Ffilter_data;
  double *Fir2_data;
  double *Gfilter_data;
  double *Gir_data;
  double *timeSig_data;
  unsigned int peakHistory[30];
  unsigned int peakHistory_data[30];
  int Cir_size[2];
  int peakHistory_size[2];
  int tempo_size[2];
  int Cir2;
  int b_Cir2;
  int d_Cir2;
  int exitg1;
  int i;
  unsigned int j;
  int k;
  int loop_ub;
  int nx;
  unsigned int qY;
  int *r1;
  signed char tmp_data[30];
  boolean_T *ri_data;
  /*     %% Create matched filters */
  /* tic */
  factor = rt_powd_snf(2.0, Octive);
  C = 16.3516 * factor;
  D = 18.35405 * factor;
  E = 20.60172 * factor;
  F = 21.82676 * factor;
  G = 24.49971 * factor;
  /* A = 27.5    *factor; */
  /* B = 30.86771*factor; */
  high_C = 32.7032 * factor;
  /* timeSig = 0:(1/fs_Hz):1;	%create a time vector */
  /* ------------------------------------------------------------------Parameter
   */
  factor = 1.0 / fs_Hz;
  b = 10.0 / C;
  emxInit_real_T(&timeSig, 2);
  timeSig_data = timeSig->data;
  if (rtIsNaN(factor) || rtIsNaN(b)) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    timeSig->size[1] = 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    timeSig_data[0] = rtNaN;
  } else if ((factor == 0.0) || ((0.0 < b) && (factor < 0.0))) {
    timeSig->size[0] = 1;
    timeSig->size[1] = 0;
  } else if (rtIsInf(b) && (rtIsInf(factor) || (0.0 == b))) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    timeSig->size[1] = 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    timeSig_data[0] = rtNaN;
  } else if (rtIsInf(factor)) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    timeSig->size[1] = 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    timeSig_data[0] = 0.0;
  } else if (floor(factor) == factor) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    loop_ub = (int)floor(b / factor);
    timeSig->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    for (Cir2 = 0; Cir2 <= loop_ub; Cir2++) {
      timeSig_data[Cir2] = factor * (double)Cir2;
    }
  } else {
    eml_float_colon(factor, b, timeSig);
    timeSig_data = timeSig->data;
  }
  emxInit_real_T(&Cfilter, 2);
  Cir2 = Cfilter->size[0] * Cfilter->size[1];
  Cfilter->size[0] = 1;
  Cfilter->size[1] = timeSig->size[1];
  emxEnsureCapacity_real_T(Cfilter, Cir2);
  Cfilter_data = Cfilter->data;
  loop_ub = timeSig->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cfilter_data[Cir2] = C * timeSig_data[Cir2] * 2.0 * 3.1415926535897931;
  }
  nx = Cfilter->size[1];
  for (k = 0; k < nx; k++) {
    Cfilter_data[k] = sin(Cfilter_data[k]);
  }
  fliplr(Cfilter);
  Cfilter_data = Cfilter->data;
  b = 10.0 / D;
  if (rtIsNaN(factor) || rtIsNaN(b)) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    timeSig->size[1] = 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    timeSig_data[0] = rtNaN;
  } else if ((factor == 0.0) || ((0.0 < b) && (factor < 0.0))) {
    timeSig->size[0] = 1;
    timeSig->size[1] = 0;
  } else if (rtIsInf(b) && (rtIsInf(factor) || (0.0 == b))) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    timeSig->size[1] = 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    timeSig_data[0] = rtNaN;
  } else if (rtIsInf(factor)) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    timeSig->size[1] = 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    timeSig_data[0] = 0.0;
  } else if (floor(factor) == factor) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    loop_ub = (int)floor(b / factor);
    timeSig->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    for (Cir2 = 0; Cir2 <= loop_ub; Cir2++) {
      timeSig_data[Cir2] = factor * (double)Cir2;
    }
  } else {
    eml_float_colon(factor, b, timeSig);
    timeSig_data = timeSig->data;
  }
  emxInit_real_T(&Dfilter, 2);
  Cir2 = Dfilter->size[0] * Dfilter->size[1];
  Dfilter->size[0] = 1;
  Dfilter->size[1] = timeSig->size[1];
  emxEnsureCapacity_real_T(Dfilter, Cir2);
  Dfilter_data = Dfilter->data;
  loop_ub = timeSig->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Dfilter_data[Cir2] = D * timeSig_data[Cir2] * 2.0 * 3.1415926535897931;
  }
  nx = Dfilter->size[1];
  for (k = 0; k < nx; k++) {
    Dfilter_data[k] = sin(Dfilter_data[k]);
  }
  fliplr(Dfilter);
  Dfilter_data = Dfilter->data;
  b = 10.0 / E;
  if (rtIsNaN(factor) || rtIsNaN(b)) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    timeSig->size[1] = 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    timeSig_data[0] = rtNaN;
  } else if ((factor == 0.0) || ((0.0 < b) && (factor < 0.0))) {
    timeSig->size[0] = 1;
    timeSig->size[1] = 0;
  } else if (rtIsInf(b) && (rtIsInf(factor) || (0.0 == b))) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    timeSig->size[1] = 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    timeSig_data[0] = rtNaN;
  } else if (rtIsInf(factor)) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    timeSig->size[1] = 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    timeSig_data[0] = 0.0;
  } else if (floor(factor) == factor) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    loop_ub = (int)floor(b / factor);
    timeSig->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    for (Cir2 = 0; Cir2 <= loop_ub; Cir2++) {
      timeSig_data[Cir2] = factor * (double)Cir2;
    }
  } else {
    eml_float_colon(factor, b, timeSig);
    timeSig_data = timeSig->data;
  }
  emxInit_real_T(&Efilter, 2);
  Cir2 = Efilter->size[0] * Efilter->size[1];
  Efilter->size[0] = 1;
  Efilter->size[1] = timeSig->size[1];
  emxEnsureCapacity_real_T(Efilter, Cir2);
  Efilter_data = Efilter->data;
  loop_ub = timeSig->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Efilter_data[Cir2] = E * timeSig_data[Cir2] * 2.0 * 3.1415926535897931;
  }
  nx = Efilter->size[1];
  for (k = 0; k < nx; k++) {
    Efilter_data[k] = sin(Efilter_data[k]);
  }
  fliplr(Efilter);
  Efilter_data = Efilter->data;
  b = 10.0 / F;
  if (rtIsNaN(factor) || rtIsNaN(b)) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    timeSig->size[1] = 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    timeSig_data[0] = rtNaN;
  } else if ((factor == 0.0) || ((0.0 < b) && (factor < 0.0))) {
    timeSig->size[0] = 1;
    timeSig->size[1] = 0;
  } else if (rtIsInf(b) && (rtIsInf(factor) || (0.0 == b))) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    timeSig->size[1] = 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    timeSig_data[0] = rtNaN;
  } else if (rtIsInf(factor)) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    timeSig->size[1] = 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    timeSig_data[0] = 0.0;
  } else if (floor(factor) == factor) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    loop_ub = (int)floor(b / factor);
    timeSig->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    for (Cir2 = 0; Cir2 <= loop_ub; Cir2++) {
      timeSig_data[Cir2] = factor * (double)Cir2;
    }
  } else {
    eml_float_colon(factor, b, timeSig);
    timeSig_data = timeSig->data;
  }
  emxInit_real_T(&Ffilter, 2);
  Cir2 = Ffilter->size[0] * Ffilter->size[1];
  Ffilter->size[0] = 1;
  Ffilter->size[1] = timeSig->size[1];
  emxEnsureCapacity_real_T(Ffilter, Cir2);
  Ffilter_data = Ffilter->data;
  loop_ub = timeSig->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Ffilter_data[Cir2] = F * timeSig_data[Cir2] * 2.0 * 3.1415926535897931;
  }
  nx = Ffilter->size[1];
  for (k = 0; k < nx; k++) {
    Ffilter_data[k] = sin(Ffilter_data[k]);
  }
  fliplr(Ffilter);
  Ffilter_data = Ffilter->data;
  b = 10.0 / G;
  if (rtIsNaN(factor) || rtIsNaN(b)) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    timeSig->size[1] = 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    timeSig_data[0] = rtNaN;
  } else if ((factor == 0.0) || ((0.0 < b) && (factor < 0.0))) {
    timeSig->size[0] = 1;
    timeSig->size[1] = 0;
  } else if (rtIsInf(b) && (rtIsInf(factor) || (0.0 == b))) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    timeSig->size[1] = 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    timeSig_data[0] = rtNaN;
  } else if (rtIsInf(factor)) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    timeSig->size[1] = 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    timeSig_data[0] = 0.0;
  } else if (floor(factor) == factor) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    loop_ub = (int)floor(b / factor);
    timeSig->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    for (Cir2 = 0; Cir2 <= loop_ub; Cir2++) {
      timeSig_data[Cir2] = factor * (double)Cir2;
    }
  } else {
    eml_float_colon(factor, b, timeSig);
    timeSig_data = timeSig->data;
  }
  emxInit_real_T(&Gfilter, 2);
  Cir2 = Gfilter->size[0] * Gfilter->size[1];
  Gfilter->size[0] = 1;
  Gfilter->size[1] = timeSig->size[1];
  emxEnsureCapacity_real_T(Gfilter, Cir2);
  Gfilter_data = Gfilter->data;
  loop_ub = timeSig->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Gfilter_data[Cir2] = G * timeSig_data[Cir2] * 2.0 * 3.1415926535897931;
  }
  nx = Gfilter->size[1];
  for (k = 0; k < nx; k++) {
    Gfilter_data[k] = sin(Gfilter_data[k]);
  }
  fliplr(Gfilter);
  Gfilter_data = Gfilter->data;
  /*      timeSig = 0:factor:numPeriods/A; */
  /*      Afilter = fliplr(sin(A*timeSig*2*pi)); */
  /*      timeSig = 0:factor:numPeriods/B; */
  /*      Bfilter = fliplr(sin(B*timeSig*2*pi)); */
  b = 10.0 / high_C;
  if (rtIsNaN(factor) || rtIsNaN(b)) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    timeSig->size[1] = 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    timeSig_data[0] = rtNaN;
  } else if ((factor == 0.0) || ((0.0 < b) && (factor < 0.0))) {
    timeSig->size[1] = 0;
  } else if (rtIsInf(b) && (rtIsInf(factor) || (0.0 == b))) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    timeSig->size[1] = 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    timeSig_data[0] = rtNaN;
  } else if (rtIsInf(factor)) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    timeSig->size[1] = 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    timeSig_data[0] = 0.0;
  } else if (floor(factor) == factor) {
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    loop_ub = (int)floor(b / factor);
    timeSig->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    for (Cir2 = 0; Cir2 <= loop_ub; Cir2++) {
      timeSig_data[Cir2] = factor * (double)Cir2;
    }
  } else {
    eml_float_colon(factor, b, timeSig);
  }
  Cir2 = timeSig->size[0] * timeSig->size[1];
  timeSig->size[0] = 1;
  emxEnsureCapacity_real_T(timeSig, Cir2);
  timeSig_data = timeSig->data;
  loop_ub = timeSig->size[1] - 1;
  for (Cir2 = 0; Cir2 <= loop_ub; Cir2++) {
    timeSig_data[Cir2] = high_C * timeSig_data[Cir2] * 2.0 * 3.1415926535897931;
  }
  nx = timeSig->size[1];
  for (k = 0; k < nx; k++) {
    timeSig_data[k] = sin(timeSig_data[k]);
  }
  emxInit_real_T(&Cir, 2);
  fliplr(timeSig);
  timeSig_data = timeSig->data;
  /*     %% Calculate frequency interference */
  /*      CC(1,1) = max(conv(Afilter,fliplr(Afilter))); */
  /*      CC(1,2) = max(conv(Afilter,fliplr(Bfilter))); */
  /*      CC(1,3) = max(conv(Afilter,fliplr(Cfilter))); */
  /*      CC(1,4) = max(conv(Afilter,fliplr(Dfilter))); */
  /*      CC(1,5) = max(conv(Afilter,fliplr(Efilter))); */
  /*      CC(1,6) = max(conv(Afilter,fliplr(Ffilter))); */
  /*      CC(1,7) = max(conv(Afilter,fliplr(Gfilter))); */
  /*      CC(1,8) = max(conv(Afilter,fliplr(high_Cfilter))); */
  /*   */
  /*      CC(2,1) = max(conv(Bfilter,fliplr(Afilter))); */
  /*      CC(2,2) = max(conv(Bfilter,fliplr(Bfilter))); */
  /*      CC(2,3) = max(conv(Bfilter,fliplr(Cfilter))); */
  /*      CC(2,4) = max(conv(Bfilter,fliplr(Dfilter))); */
  /*      CC(2,5) = max(conv(Bfilter,fliplr(Efilter))); */
  /*      CC(2,6) = max(conv(Bfilter,fliplr(Ffilter))); */
  /*      CC(2,7) = max(conv(Bfilter,fliplr(Gfilter))); */
  /*      CC(2,8) = max(conv(Bfilter,fliplr(high_Cfilter))); */
  /* CC(3,1) = max(conv(Cfilter,fliplr(Afilter))); */
  /* CC(3,2) = max(conv(Cfilter,fliplr(Bfilter))); */
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Cfilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Cfilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Cfilter_data[Cir2];
  }
  emxInit_real_T(&Gir, 2);
  fliplr(Cir);
  conv(Cfilter, Cir, Gir);
  CC[0] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Dfilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Dfilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Dfilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Cfilter, Cir, Gir);
  CC[6] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Efilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Efilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Efilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Cfilter, Cir, Gir);
  CC[12] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Ffilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Ffilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Ffilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Cfilter, Cir, Gir);
  CC[18] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Gfilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Gfilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Gfilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Cfilter, Cir, Gir);
  CC[24] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = timeSig->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = timeSig->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = timeSig_data[Cir2];
  }
  fliplr(Cir);
  conv(Cfilter, Cir, Gir);
  CC[30] = maximum(Gir);
  /* CC(4,1) = max(conv(Dfilter,fliplr(Afilter))); */
  /* CC(4,2) = max(conv(Dfilter,fliplr(Bfilter))); */
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Cfilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Cfilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Cfilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Dfilter, Cir, Gir);
  CC[1] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Dfilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Dfilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Dfilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Dfilter, Cir, Gir);
  CC[7] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Efilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Efilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Efilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Dfilter, Cir, Gir);
  CC[13] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Ffilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Ffilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Ffilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Dfilter, Cir, Gir);
  CC[19] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Gfilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Gfilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Gfilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Dfilter, Cir, Gir);
  CC[25] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = timeSig->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = timeSig->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = timeSig_data[Cir2];
  }
  fliplr(Cir);
  conv(Dfilter, Cir, Gir);
  CC[31] = maximum(Gir);
  /* CC(5,1) = max(conv(Efilter,fliplr(Afilter))); */
  /* CC(5,2) = max(conv(Efilter,fliplr(Bfilter))); */
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Cfilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Cfilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Cfilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Efilter, Cir, Gir);
  CC[2] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Dfilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Dfilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Dfilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Efilter, Cir, Gir);
  CC[8] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Efilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Efilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Efilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Efilter, Cir, Gir);
  CC[14] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Ffilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Ffilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Ffilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Efilter, Cir, Gir);
  CC[20] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Gfilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Gfilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Gfilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Efilter, Cir, Gir);
  CC[26] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = timeSig->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = timeSig->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = timeSig_data[Cir2];
  }
  fliplr(Cir);
  conv(Efilter, Cir, Gir);
  CC[32] = maximum(Gir);
  /* CC(6,1) = max(conv(Ffilter,fliplr(Afilter))); */
  /* CC(6,2) = max(conv(Ffilter,fliplr(Bfilter))); */
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Cfilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Cfilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Cfilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Ffilter, Cir, Gir);
  CC[3] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Dfilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Dfilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Dfilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Ffilter, Cir, Gir);
  CC[9] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Efilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Efilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Efilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Ffilter, Cir, Gir);
  CC[15] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Ffilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Ffilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Ffilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Ffilter, Cir, Gir);
  CC[21] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Gfilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Gfilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Gfilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Ffilter, Cir, Gir);
  CC[27] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = timeSig->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = timeSig->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = timeSig_data[Cir2];
  }
  fliplr(Cir);
  conv(Ffilter, Cir, Gir);
  CC[33] = maximum(Gir);
  /* CC(7,1) = max(conv(Gfilter,fliplr(Afilter))); */
  /* CC(7,2) = max(conv(Gfilter,fliplr(Bfilter))); */
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Cfilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Cfilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Cfilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Gfilter, Cir, Gir);
  CC[4] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Dfilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Dfilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Dfilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Gfilter, Cir, Gir);
  CC[10] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Efilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Efilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Efilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Gfilter, Cir, Gir);
  CC[16] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Ffilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Ffilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Ffilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Gfilter, Cir, Gir);
  CC[22] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Gfilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Gfilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Gfilter_data[Cir2];
  }
  fliplr(Cir);
  conv(Gfilter, Cir, Gir);
  CC[28] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = timeSig->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = timeSig->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = timeSig_data[Cir2];
  }
  fliplr(Cir);
  conv(Gfilter, Cir, Gir);
  CC[34] = maximum(Gir);
  /* CC(8,1) = max(conv(high_Cfilter,fliplr(Afilter))); */
  /* CC(8,2) = max(conv(high_Cfilter,fliplr(Bfilter))); */
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Cfilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Cfilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Cfilter_data[Cir2];
  }
  fliplr(Cir);
  conv(timeSig, Cir, Gir);
  CC[5] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Dfilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Dfilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Dfilter_data[Cir2];
  }
  fliplr(Cir);
  conv(timeSig, Cir, Gir);
  CC[11] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Efilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Efilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Efilter_data[Cir2];
  }
  fliplr(Cir);
  conv(timeSig, Cir, Gir);
  CC[17] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Ffilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Ffilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Ffilter_data[Cir2];
  }
  fliplr(Cir);
  conv(timeSig, Cir, Gir);
  CC[23] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = Gfilter->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = Gfilter->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = Gfilter_data[Cir2];
  }
  fliplr(Cir);
  conv(timeSig, Cir, Gir);
  CC[29] = maximum(Gir);
  Cir2 = Cir->size[0] * Cir->size[1];
  Cir->size[0] = 1;
  Cir->size[1] = timeSig->size[1];
  emxEnsureCapacity_real_T(Cir, Cir2);
  Cir_data = Cir->data;
  loop_ub = timeSig->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir_data[Cir2] = timeSig_data[Cir2];
  }
  fliplr(Cir);
  conv(timeSig, Cir, Gir);
  CC[35] = maximum(Gir);
  /*      figure() */
  /*      for i = 1:8 */
  /*          subplot(3,3,i) */
  /*          plot(CC(i,:)) */
  /*      end */
  linsolve(CC, InterdependencyMatrix);
  /*     %% Filter dignal for notes      %Do the envelope detection */
  b_conv(digital, Cfilter, Cir);
  Cir_data = Cir->data;
  b_conv(digital, Dfilter, Cfilter);
  Cfilter_data = Cfilter->data;
  b_conv(digital, Efilter, Dfilter);
  Dfilter_data = Dfilter->data;
  b_conv(digital, Ffilter, Efilter);
  Efilter_data = Efilter->data;
  b_conv(digital, Gfilter, Gir);
  Gir_data = Gir->data;
  /* Air = conv(digital,Afilter); */
  /* Bir = conv(digital,Bfilter); */
  b_conv(digital, timeSig, Ffilter);
  Ffilter_data = Ffilter->data;
  /*     %% Filter interdependencies */
  /* ceil(100*fs_ratio); */
  /* Air2 = conv(abs(Air(1:trimLen)),ones(smoothLen,1)*1/smoothLen); */
  /* Bir2 = conv(abs(Bir(1:trimLen)),ones(smoothLen,1)*1/smoothLen); */
  Cir2 = Ffilter->size[1];
  b_Cir2 = Gfilter->size[0] * Gfilter->size[1];
  Gfilter->size[0] = 1;
  Gfilter->size[1] = Ffilter->size[1];
  emxEnsureCapacity_real_T(Gfilter, b_Cir2);
  Gfilter_data = Gfilter->data;
  for (k = 0; k < Cir2; k++) {
    Gfilter_data[k] = fabs(Cir_data[k]);
  }
  emxInit_real_T(&c_Cir2, 2);
  Cir2 = c_Cir2->size[0] * c_Cir2->size[1];
  c_Cir2->size[0] = 1;
  c_Cir2->size[1] = Gfilter->size[1] + 199;
  emxEnsureCapacity_real_T(c_Cir2, Cir2);
  Cir2_data = c_Cir2->data;
  loop_ub = Gfilter->size[1] + 199;
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Cir2_data[Cir2] = 0.0;
  }
  Cir2 = Gfilter->size[1] - 1;
  for (k = 0; k < 200; k++) {
    for (loop_ub = 0; loop_ub <= Cir2; loop_ub++) {
      nx = k + loop_ub;
      Cir2_data[nx] += 0.005 * Gfilter_data[loop_ub];
    }
  }
  Cir2 = Ffilter->size[1];
  b_Cir2 = Gfilter->size[0] * Gfilter->size[1];
  Gfilter->size[0] = 1;
  Gfilter->size[1] = Ffilter->size[1];
  emxEnsureCapacity_real_T(Gfilter, b_Cir2);
  Gfilter_data = Gfilter->data;
  for (k = 0; k < Cir2; k++) {
    Gfilter_data[k] = fabs(Cfilter_data[k]);
  }
  emxInit_real_T(&Dir2, 2);
  Cir2 = Dir2->size[0] * Dir2->size[1];
  Dir2->size[0] = 1;
  Dir2->size[1] = Gfilter->size[1] + 199;
  emxEnsureCapacity_real_T(Dir2, Cir2);
  Dir2_data = Dir2->data;
  loop_ub = Gfilter->size[1] + 199;
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Dir2_data[Cir2] = 0.0;
  }
  Cir2 = Gfilter->size[1] - 1;
  for (k = 0; k < 200; k++) {
    for (loop_ub = 0; loop_ub <= Cir2; loop_ub++) {
      nx = k + loop_ub;
      Dir2_data[nx] += 0.005 * Gfilter_data[loop_ub];
    }
  }
  Cir2 = Ffilter->size[1];
  b_Cir2 = Gfilter->size[0] * Gfilter->size[1];
  Gfilter->size[0] = 1;
  Gfilter->size[1] = Ffilter->size[1];
  emxEnsureCapacity_real_T(Gfilter, b_Cir2);
  Gfilter_data = Gfilter->data;
  for (k = 0; k < Cir2; k++) {
    Gfilter_data[k] = fabs(Dfilter_data[k]);
  }
  emxInit_real_T(&Eir2, 2);
  Cir2 = Eir2->size[0] * Eir2->size[1];
  Eir2->size[0] = 1;
  Eir2->size[1] = Gfilter->size[1] + 199;
  emxEnsureCapacity_real_T(Eir2, Cir2);
  Eir2_data = Eir2->data;
  loop_ub = Gfilter->size[1] + 199;
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Eir2_data[Cir2] = 0.0;
  }
  Cir2 = Gfilter->size[1] - 1;
  for (k = 0; k < 200; k++) {
    for (loop_ub = 0; loop_ub <= Cir2; loop_ub++) {
      nx = k + loop_ub;
      Eir2_data[nx] += 0.005 * Gfilter_data[loop_ub];
    }
  }
  Cir2 = Ffilter->size[1];
  b_Cir2 = Gfilter->size[0] * Gfilter->size[1];
  Gfilter->size[0] = 1;
  Gfilter->size[1] = Ffilter->size[1];
  emxEnsureCapacity_real_T(Gfilter, b_Cir2);
  Gfilter_data = Gfilter->data;
  for (k = 0; k < Cir2; k++) {
    Gfilter_data[k] = fabs(Efilter_data[k]);
  }
  emxInit_real_T(&Fir2, 2);
  Cir2 = Fir2->size[0] * Fir2->size[1];
  Fir2->size[0] = 1;
  Fir2->size[1] = Gfilter->size[1] + 199;
  emxEnsureCapacity_real_T(Fir2, Cir2);
  Fir2_data = Fir2->data;
  loop_ub = Gfilter->size[1] + 199;
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Fir2_data[Cir2] = 0.0;
  }
  Cir2 = Gfilter->size[1] - 1;
  for (k = 0; k < 200; k++) {
    for (loop_ub = 0; loop_ub <= Cir2; loop_ub++) {
      nx = k + loop_ub;
      Fir2_data[nx] += 0.005 * Gfilter_data[loop_ub];
    }
  }
  Cir2 = Ffilter->size[1];
  b_Cir2 = Gfilter->size[0] * Gfilter->size[1];
  Gfilter->size[0] = 1;
  Gfilter->size[1] = Ffilter->size[1];
  emxEnsureCapacity_real_T(Gfilter, b_Cir2);
  Gfilter_data = Gfilter->data;
  for (k = 0; k < Cir2; k++) {
    Gfilter_data[k] = fabs(Gir_data[k]);
  }
  Cir2 = timeSig->size[0] * timeSig->size[1];
  timeSig->size[0] = 1;
  timeSig->size[1] = Gfilter->size[1] + 199;
  emxEnsureCapacity_real_T(timeSig, Cir2);
  timeSig_data = timeSig->data;
  loop_ub = Gfilter->size[1] + 199;
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    timeSig_data[Cir2] = 0.0;
  }
  Cir2 = Gfilter->size[1] - 1;
  for (k = 0; k < 200; k++) {
    for (loop_ub = 0; loop_ub <= Cir2; loop_ub++) {
      nx = k + loop_ub;
      timeSig_data[nx] += 0.005 * Gfilter_data[loop_ub];
    }
  }
  Cir2 = Ffilter->size[1];
  b_Cir2 = Gfilter->size[0] * Gfilter->size[1];
  Gfilter->size[0] = 1;
  Gfilter->size[1] = Ffilter->size[1];
  emxEnsureCapacity_real_T(Gfilter, b_Cir2);
  Gfilter_data = Gfilter->data;
  for (k = 0; k < Cir2; k++) {
    Gfilter_data[k] = fabs(Ffilter_data[k]);
  }
  Cir2 = Ffilter->size[0] * Ffilter->size[1];
  Ffilter->size[0] = 1;
  Ffilter->size[1] = Gfilter->size[1] + 199;
  emxEnsureCapacity_real_T(Ffilter, Cir2);
  Ffilter_data = Ffilter->data;
  loop_ub = Gfilter->size[1] + 199;
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Ffilter_data[Cir2] = 0.0;
  }
  Cir2 = Gfilter->size[1] - 1;
  for (k = 0; k < 200; k++) {
    for (loop_ub = 0; loop_ub <= Cir2; loop_ub++) {
      nx = k + loop_ub;
      Ffilter_data[nx] += 0.005 * Gfilter_data[loop_ub];
    }
  }
  emxFree_real_T(&Gfilter);
  /* figure(); */
  /*      rowNum = 1; */
  /*      Air3 =
   * InterdependencyMatrix(rowNum,1)*Air2+InterdependencyMatrix(rowNum,2)*Bir2+InterdependencyMatrix(rowNum,3)*Cir2+InterdependencyMatrix(rowNum,4)*Dir2+InterdependencyMatrix(rowNum,5)*Eir2+InterdependencyMatrix(rowNum,6)*Fir2+InterdependencyMatrix(rowNum,7)*Gir2+InterdependencyMatrix(rowNum,8)*high_Cir2;
   */
  /*      subplot(3,3,rowNum) */
  /*      plot(Air3) */
  /*      rowNum = 2; */
  /*      Bir3 =
   * InterdependencyMatrix(rowNum,1)*Air2+InterdependencyMatrix(rowNum,2)*Bir2+InterdependencyMatrix(rowNum,3)*Cir2+InterdependencyMatrix(rowNum,4)*Dir2+InterdependencyMatrix(rowNum,5)*Eir2+InterdependencyMatrix(rowNum,6)*Fir2+InterdependencyMatrix(rowNum,7)*Gir2+InterdependencyMatrix(rowNum,8)*high_Cir2;
   */
  /*      subplot(3,3,rowNum) */
  /*      plot(Bir3) */
  if (c_Cir2->size[1] == 1) {
    nx = Dir2->size[1];
  } else {
    nx = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    d_Cir2 = Dir2->size[1];
  } else {
    d_Cir2 = c_Cir2->size[1];
  }
  if (d_Cir2 == 1) {
    d_Cir2 = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    d_Cir2 = Dir2->size[1];
  } else {
    d_Cir2 = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    Cir2 = Dir2->size[1];
  } else {
    Cir2 = c_Cir2->size[1];
  }
  if (Cir2 == 1) {
    Cir2 = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    Cir2 = Dir2->size[1];
  } else {
    Cir2 = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (Cir2 == 1) {
    Cir2 = Fir2->size[1];
  } else if (loop_ub == 1) {
    Cir2 = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    Cir2 = Dir2->size[1];
  } else {
    Cir2 = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    k = Dir2->size[1];
  } else {
    k = c_Cir2->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = Fir2->size[1];
  } else if (k == 1) {
    loop_ub = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    k = Dir2->size[1];
  } else {
    k = c_Cir2->size[1];
  }
  if (k == 1) {
    k = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    k = Dir2->size[1];
  } else {
    k = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    b_Cir2 = Dir2->size[1];
  } else {
    b_Cir2 = c_Cir2->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = timeSig->size[1];
  } else if (k == 1) {
    loop_ub = Fir2->size[1];
  } else if (b_Cir2 == 1) {
    loop_ub = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if ((c_Cir2->size[1] == Dir2->size[1]) && (nx == Eir2->size[1]) &&
      (d_Cir2 == Fir2->size[1]) && (Cir2 == timeSig->size[1]) &&
      (loop_ub == Ffilter->size[1])) {
    Cir2 = Cir->size[0] * Cir->size[1];
    Cir->size[0] = 1;
    Cir->size[1] = c_Cir2->size[1];
    emxEnsureCapacity_real_T(Cir, Cir2);
    Cir_data = Cir->data;
    loop_ub = c_Cir2->size[1];
    for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
      Cir_data[Cir2] = ((((InterdependencyMatrix[0] * Cir2_data[Cir2] +
                           InterdependencyMatrix[6] * Dir2_data[Cir2]) +
                          InterdependencyMatrix[12] * Eir2_data[Cir2]) +
                         InterdependencyMatrix[18] * Fir2_data[Cir2]) +
                        InterdependencyMatrix[24] * timeSig_data[Cir2]) +
                       InterdependencyMatrix[30] * Ffilter_data[Cir2];
    }
  } else {
    h_binary_expand_op(Cir, InterdependencyMatrix, c_Cir2, Dir2, Eir2, Fir2,
                       timeSig, Ffilter);
    Cir_data = Cir->data;
  }
  /* subplot(3,3,rowNum) */
  /* plot(Cir3) */
  if (c_Cir2->size[1] == 1) {
    nx = Dir2->size[1];
  } else {
    nx = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    d_Cir2 = Dir2->size[1];
  } else {
    d_Cir2 = c_Cir2->size[1];
  }
  if (d_Cir2 == 1) {
    d_Cir2 = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    d_Cir2 = Dir2->size[1];
  } else {
    d_Cir2 = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    Cir2 = Dir2->size[1];
  } else {
    Cir2 = c_Cir2->size[1];
  }
  if (Cir2 == 1) {
    Cir2 = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    Cir2 = Dir2->size[1];
  } else {
    Cir2 = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (Cir2 == 1) {
    Cir2 = Fir2->size[1];
  } else if (loop_ub == 1) {
    Cir2 = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    Cir2 = Dir2->size[1];
  } else {
    Cir2 = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    k = Dir2->size[1];
  } else {
    k = c_Cir2->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = Fir2->size[1];
  } else if (k == 1) {
    loop_ub = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    k = Dir2->size[1];
  } else {
    k = c_Cir2->size[1];
  }
  if (k == 1) {
    k = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    k = Dir2->size[1];
  } else {
    k = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    b_Cir2 = Dir2->size[1];
  } else {
    b_Cir2 = c_Cir2->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = timeSig->size[1];
  } else if (k == 1) {
    loop_ub = Fir2->size[1];
  } else if (b_Cir2 == 1) {
    loop_ub = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if ((c_Cir2->size[1] == Dir2->size[1]) && (nx == Eir2->size[1]) &&
      (d_Cir2 == Fir2->size[1]) && (Cir2 == timeSig->size[1]) &&
      (loop_ub == Ffilter->size[1])) {
    Cir2 = Gir->size[0] * Gir->size[1];
    Gir->size[0] = 1;
    Gir->size[1] = c_Cir2->size[1];
    emxEnsureCapacity_real_T(Gir, Cir2);
    Gir_data = Gir->data;
    loop_ub = c_Cir2->size[1];
    for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
      Gir_data[Cir2] = ((((InterdependencyMatrix[1] * Cir2_data[Cir2] +
                           InterdependencyMatrix[7] * Dir2_data[Cir2]) +
                          InterdependencyMatrix[13] * Eir2_data[Cir2]) +
                         InterdependencyMatrix[19] * Fir2_data[Cir2]) +
                        InterdependencyMatrix[25] * timeSig_data[Cir2]) +
                       InterdependencyMatrix[31] * Ffilter_data[Cir2];
    }
  } else {
    g_binary_expand_op(Gir, InterdependencyMatrix, c_Cir2, Dir2, Eir2, Fir2,
                       timeSig, Ffilter);
    Gir_data = Gir->data;
  }
  /* subplot(3,3,rowNum) */
  /* plot(Dir3) */
  if (c_Cir2->size[1] == 1) {
    nx = Dir2->size[1];
  } else {
    nx = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    d_Cir2 = Dir2->size[1];
  } else {
    d_Cir2 = c_Cir2->size[1];
  }
  if (d_Cir2 == 1) {
    d_Cir2 = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    d_Cir2 = Dir2->size[1];
  } else {
    d_Cir2 = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    Cir2 = Dir2->size[1];
  } else {
    Cir2 = c_Cir2->size[1];
  }
  if (Cir2 == 1) {
    Cir2 = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    Cir2 = Dir2->size[1];
  } else {
    Cir2 = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (Cir2 == 1) {
    Cir2 = Fir2->size[1];
  } else if (loop_ub == 1) {
    Cir2 = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    Cir2 = Dir2->size[1];
  } else {
    Cir2 = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    k = Dir2->size[1];
  } else {
    k = c_Cir2->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = Fir2->size[1];
  } else if (k == 1) {
    loop_ub = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    k = Dir2->size[1];
  } else {
    k = c_Cir2->size[1];
  }
  if (k == 1) {
    k = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    k = Dir2->size[1];
  } else {
    k = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    b_Cir2 = Dir2->size[1];
  } else {
    b_Cir2 = c_Cir2->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = timeSig->size[1];
  } else if (k == 1) {
    loop_ub = Fir2->size[1];
  } else if (b_Cir2 == 1) {
    loop_ub = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if ((c_Cir2->size[1] == Dir2->size[1]) && (nx == Eir2->size[1]) &&
      (d_Cir2 == Fir2->size[1]) && (Cir2 == timeSig->size[1]) &&
      (loop_ub == Ffilter->size[1])) {
    Cir2 = Cfilter->size[0] * Cfilter->size[1];
    Cfilter->size[0] = 1;
    Cfilter->size[1] = c_Cir2->size[1];
    emxEnsureCapacity_real_T(Cfilter, Cir2);
    Cfilter_data = Cfilter->data;
    loop_ub = c_Cir2->size[1];
    for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
      Cfilter_data[Cir2] = ((((InterdependencyMatrix[2] * Cir2_data[Cir2] +
                               InterdependencyMatrix[8] * Dir2_data[Cir2]) +
                              InterdependencyMatrix[14] * Eir2_data[Cir2]) +
                             InterdependencyMatrix[20] * Fir2_data[Cir2]) +
                            InterdependencyMatrix[26] * timeSig_data[Cir2]) +
                           InterdependencyMatrix[32] * Ffilter_data[Cir2];
    }
  } else {
    f_binary_expand_op(Cfilter, InterdependencyMatrix, c_Cir2, Dir2, Eir2, Fir2,
                       timeSig, Ffilter);
    Cfilter_data = Cfilter->data;
  }
  /* subplot(3,3,rowNum) */
  /* plot(Eir3) */
  if (c_Cir2->size[1] == 1) {
    nx = Dir2->size[1];
  } else {
    nx = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    d_Cir2 = Dir2->size[1];
  } else {
    d_Cir2 = c_Cir2->size[1];
  }
  if (d_Cir2 == 1) {
    d_Cir2 = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    d_Cir2 = Dir2->size[1];
  } else {
    d_Cir2 = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    Cir2 = Dir2->size[1];
  } else {
    Cir2 = c_Cir2->size[1];
  }
  if (Cir2 == 1) {
    Cir2 = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    Cir2 = Dir2->size[1];
  } else {
    Cir2 = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (Cir2 == 1) {
    Cir2 = Fir2->size[1];
  } else if (loop_ub == 1) {
    Cir2 = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    Cir2 = Dir2->size[1];
  } else {
    Cir2 = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    k = Dir2->size[1];
  } else {
    k = c_Cir2->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = Fir2->size[1];
  } else if (k == 1) {
    loop_ub = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    k = Dir2->size[1];
  } else {
    k = c_Cir2->size[1];
  }
  if (k == 1) {
    k = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    k = Dir2->size[1];
  } else {
    k = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    b_Cir2 = Dir2->size[1];
  } else {
    b_Cir2 = c_Cir2->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = timeSig->size[1];
  } else if (k == 1) {
    loop_ub = Fir2->size[1];
  } else if (b_Cir2 == 1) {
    loop_ub = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if ((c_Cir2->size[1] == Dir2->size[1]) && (nx == Eir2->size[1]) &&
      (d_Cir2 == Fir2->size[1]) && (Cir2 == timeSig->size[1]) &&
      (loop_ub == Ffilter->size[1])) {
    Cir2 = Dfilter->size[0] * Dfilter->size[1];
    Dfilter->size[0] = 1;
    Dfilter->size[1] = c_Cir2->size[1];
    emxEnsureCapacity_real_T(Dfilter, Cir2);
    Dfilter_data = Dfilter->data;
    loop_ub = c_Cir2->size[1];
    for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
      Dfilter_data[Cir2] = ((((InterdependencyMatrix[3] * Cir2_data[Cir2] +
                               InterdependencyMatrix[9] * Dir2_data[Cir2]) +
                              InterdependencyMatrix[15] * Eir2_data[Cir2]) +
                             InterdependencyMatrix[21] * Fir2_data[Cir2]) +
                            InterdependencyMatrix[27] * timeSig_data[Cir2]) +
                           InterdependencyMatrix[33] * Ffilter_data[Cir2];
    }
  } else {
    e_binary_expand_op(Dfilter, InterdependencyMatrix, c_Cir2, Dir2, Eir2, Fir2,
                       timeSig, Ffilter);
    Dfilter_data = Dfilter->data;
  }
  /* subplot(3,3,rowNum) */
  /* plot(Fir3) */
  if (c_Cir2->size[1] == 1) {
    nx = Dir2->size[1];
  } else {
    nx = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    d_Cir2 = Dir2->size[1];
  } else {
    d_Cir2 = c_Cir2->size[1];
  }
  if (d_Cir2 == 1) {
    d_Cir2 = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    d_Cir2 = Dir2->size[1];
  } else {
    d_Cir2 = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    Cir2 = Dir2->size[1];
  } else {
    Cir2 = c_Cir2->size[1];
  }
  if (Cir2 == 1) {
    Cir2 = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    Cir2 = Dir2->size[1];
  } else {
    Cir2 = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (Cir2 == 1) {
    Cir2 = Fir2->size[1];
  } else if (loop_ub == 1) {
    Cir2 = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    Cir2 = Dir2->size[1];
  } else {
    Cir2 = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    k = Dir2->size[1];
  } else {
    k = c_Cir2->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = Fir2->size[1];
  } else if (k == 1) {
    loop_ub = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    k = Dir2->size[1];
  } else {
    k = c_Cir2->size[1];
  }
  if (k == 1) {
    k = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    k = Dir2->size[1];
  } else {
    k = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    b_Cir2 = Dir2->size[1];
  } else {
    b_Cir2 = c_Cir2->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = timeSig->size[1];
  } else if (k == 1) {
    loop_ub = Fir2->size[1];
  } else if (b_Cir2 == 1) {
    loop_ub = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if ((c_Cir2->size[1] == Dir2->size[1]) && (nx == Eir2->size[1]) &&
      (d_Cir2 == Fir2->size[1]) && (Cir2 == timeSig->size[1]) &&
      (loop_ub == Ffilter->size[1])) {
    Cir2 = Efilter->size[0] * Efilter->size[1];
    Efilter->size[0] = 1;
    Efilter->size[1] = c_Cir2->size[1];
    emxEnsureCapacity_real_T(Efilter, Cir2);
    Efilter_data = Efilter->data;
    loop_ub = c_Cir2->size[1];
    for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
      Efilter_data[Cir2] = ((((InterdependencyMatrix[4] * Cir2_data[Cir2] +
                               InterdependencyMatrix[10] * Dir2_data[Cir2]) +
                              InterdependencyMatrix[16] * Eir2_data[Cir2]) +
                             InterdependencyMatrix[22] * Fir2_data[Cir2]) +
                            InterdependencyMatrix[28] * timeSig_data[Cir2]) +
                           InterdependencyMatrix[34] * Ffilter_data[Cir2];
    }
  } else {
    d_binary_expand_op(Efilter, InterdependencyMatrix, c_Cir2, Dir2, Eir2, Fir2,
                       timeSig, Ffilter);
    Efilter_data = Efilter->data;
  }
  /* subplot(3,3,rowNum) */
  /* plot(Gir3) */
  if (c_Cir2->size[1] == 1) {
    nx = Dir2->size[1];
  } else {
    nx = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    d_Cir2 = Dir2->size[1];
  } else {
    d_Cir2 = c_Cir2->size[1];
  }
  if (d_Cir2 == 1) {
    d_Cir2 = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    d_Cir2 = Dir2->size[1];
  } else {
    d_Cir2 = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    Cir2 = Dir2->size[1];
  } else {
    Cir2 = c_Cir2->size[1];
  }
  if (Cir2 == 1) {
    Cir2 = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    Cir2 = Dir2->size[1];
  } else {
    Cir2 = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (Cir2 == 1) {
    Cir2 = Fir2->size[1];
  } else if (loop_ub == 1) {
    Cir2 = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    Cir2 = Dir2->size[1];
  } else {
    Cir2 = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    k = Dir2->size[1];
  } else {
    k = c_Cir2->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = Fir2->size[1];
  } else if (k == 1) {
    loop_ub = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    k = Dir2->size[1];
  } else {
    k = c_Cir2->size[1];
  }
  if (k == 1) {
    k = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    k = Dir2->size[1];
  } else {
    k = c_Cir2->size[1];
  }
  if (c_Cir2->size[1] == 1) {
    b_Cir2 = Dir2->size[1];
  } else {
    b_Cir2 = c_Cir2->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = timeSig->size[1];
  } else if (k == 1) {
    loop_ub = Fir2->size[1];
  } else if (b_Cir2 == 1) {
    loop_ub = Eir2->size[1];
  } else if (c_Cir2->size[1] == 1) {
    loop_ub = Dir2->size[1];
  } else {
    loop_ub = c_Cir2->size[1];
  }
  if ((c_Cir2->size[1] == Dir2->size[1]) && (nx == Eir2->size[1]) &&
      (d_Cir2 == Fir2->size[1]) && (Cir2 == timeSig->size[1]) &&
      (loop_ub == Ffilter->size[1])) {
    loop_ub = c_Cir2->size[1] - 1;
    Cir2 = c_Cir2->size[0] * c_Cir2->size[1];
    c_Cir2->size[0] = 1;
    emxEnsureCapacity_real_T(c_Cir2, Cir2);
    Cir2_data = c_Cir2->data;
    for (Cir2 = 0; Cir2 <= loop_ub; Cir2++) {
      Cir2_data[Cir2] = ((((InterdependencyMatrix[5] * Cir2_data[Cir2] +
                            InterdependencyMatrix[11] * Dir2_data[Cir2]) +
                           InterdependencyMatrix[17] * Eir2_data[Cir2]) +
                          InterdependencyMatrix[23] * Fir2_data[Cir2]) +
                         InterdependencyMatrix[29] * timeSig_data[Cir2]) +
                        InterdependencyMatrix[35] * Ffilter_data[Cir2];
    }
  } else {
    c_binary_expand_op(c_Cir2, InterdependencyMatrix, Dir2, Eir2, Fir2, timeSig,
                       Ffilter);
    Cir2_data = c_Cir2->data;
  }
  emxFree_real_T(&Fir2);
  emxFree_real_T(&Eir2);
  emxFree_real_T(&Dir2);
  /* subplot(3,3,rowNum) */
  /* plot(high_Cir3) */
  /* subplot(3,3,9) */
  /* plot(Cir3+Dir3+Eir3+Fir3+Gir3+high_Cir3) */
  /*      figure() */
  /*      subplot(3,3,3) */
  /*      plot((Cir3>0.5*max(Cir3))) */
  /*      subplot(3,3,4) */
  /*      plot((Dir3>0.5*max(Dir3))) */
  /*      subplot(3,3,5) */
  /*      plot((Eir3>0.5*max(Eir3))) */
  /*      subplot(3,3,6) */
  /*      plot((Fir3>0.5*max(Fir3))) */
  /*      subplot(3,3,7) */
  /*      plot((Gir3>0.5*max(Gir3))) */
  /*      subplot(3,3,8) */
  /*      plot((high_Cir3>0.5*max(high_Cir3))) */
  /*      subplot(3,3,9) */
  factor = maximum(Cir);
  C = maximum(Gir);
  D = maximum(Cfilter);
  E = maximum(Dfilter);
  G = maximum(Efilter);
  b = maximum(c_Cir2);
  /*      plot(pulses) */
  /*      ylim([-1,2]); */
  /*     %% threshold the impulse responses */
  /*      IMPthresh = 1;
   * %------------------------------------------------------------------Parameter
   */
  /*   */
  /*      Air(abs(Air)<IMPthresh) = 0; */
  /*      Bir(abs(Bir)<IMPthresh) = 0; */
  /*      Cir(abs(Cir)<IMPthresh) = 0; */
  /*      Dir(abs(Dir)<IMPthresh) = 0; */
  /*      Eir(abs(Eir)<IMPthresh) = 0; */
  /*      Fir(abs(Fir)<IMPthresh) = 0; */
  /*      Gir(abs(Gir)<IMPthresh) = 0; */
  /*      high_Cir(abs(high_Cir)<IMPthresh) = 0; */
  /*     %% plot note charts */
  /*      maxlength = length(Cir);
   * %------------------------------------------------------------------Parameter
   */
  /*   */
  /*      figure() */
  /*  %     subplot(3,3,1) */
  /*  %     plot(Air) */
  /*  %     subplot(3,3,2) */
  /*  %     plot(Bir) */
  /*      subplot(3,3,3) */
  /*      plot(Cir) */
  /*      subplot(3,3,4) */
  /*      plot(Dir) */
  /*      subplot(3,3,5) */
  /*      plot(Eir) */
  /*      subplot(3,3,6) */
  /*      plot(Fir) */
  /*      subplot(3,3,7) */
  /*      plot(Gir) */
  /*      subplot(3,3,8) */
  /*      plot(high_Cir) */
  /*      subplot(3,3,9) */
  /*      tot_ir = padarray(Cir,[0
   * maxlength-length(Cir)],'post')+padarray(Dir,[0
   * maxlength-length(Dir)],'post')+padarray(Eir,[0
   * maxlength-length(Eir)],'post')+padarray(Fir,[0
   * maxlength-length(Fir)],'post')+padarray(Gir,[0
   * maxlength-length(Gir)],'post')+padarray(high_Cir,[0
   * maxlength-length(high_Cir)],'post'); */
  /*      plot(tot_ir) */
  /* rest_indicator = conv(abs(tot_ir),[.1 .1 .1 .1 .1 .1 .1 .1 .1 .1]);
   * %------------------------------------------------------------------Parameter
   */
  /* rest_indicator(rest_indicator>1)=1; */
  if (Cir->size[1] == 1) {
    nx = Gir->size[1];
  } else {
    nx = Cir->size[1];
  }
  if (Cir->size[1] == 1) {
    d_Cir2 = Gir->size[1];
  } else {
    d_Cir2 = Cir->size[1];
  }
  if (d_Cir2 == 1) {
    d_Cir2 = Cfilter->size[1];
  } else if (Cir->size[1] == 1) {
    d_Cir2 = Gir->size[1];
  } else {
    d_Cir2 = Cir->size[1];
  }
  if (Cir->size[1] == 1) {
    Cir2 = Gir->size[1];
  } else {
    Cir2 = Cir->size[1];
  }
  if (Cir2 == 1) {
    Cir2 = Cfilter->size[1];
  } else if (Cir->size[1] == 1) {
    Cir2 = Gir->size[1];
  } else {
    Cir2 = Cir->size[1];
  }
  if (Cir->size[1] == 1) {
    loop_ub = Gir->size[1];
  } else {
    loop_ub = Cir->size[1];
  }
  if (Cir2 == 1) {
    Cir2 = Dfilter->size[1];
  } else if (loop_ub == 1) {
    Cir2 = Cfilter->size[1];
  } else if (Cir->size[1] == 1) {
    Cir2 = Gir->size[1];
  } else {
    Cir2 = Cir->size[1];
  }
  if (Cir->size[1] == 1) {
    loop_ub = Gir->size[1];
  } else {
    loop_ub = Cir->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = Cfilter->size[1];
  } else if (Cir->size[1] == 1) {
    loop_ub = Gir->size[1];
  } else {
    loop_ub = Cir->size[1];
  }
  if (Cir->size[1] == 1) {
    k = Gir->size[1];
  } else {
    k = Cir->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = Dfilter->size[1];
  } else if (k == 1) {
    loop_ub = Cfilter->size[1];
  } else if (Cir->size[1] == 1) {
    loop_ub = Gir->size[1];
  } else {
    loop_ub = Cir->size[1];
  }
  if (Cir->size[1] == 1) {
    k = Gir->size[1];
  } else {
    k = Cir->size[1];
  }
  if (k == 1) {
    k = Cfilter->size[1];
  } else if (Cir->size[1] == 1) {
    k = Gir->size[1];
  } else {
    k = Cir->size[1];
  }
  if (Cir->size[1] == 1) {
    b_Cir2 = Gir->size[1];
  } else {
    b_Cir2 = Cir->size[1];
  }
  if (loop_ub == 1) {
    loop_ub = Efilter->size[1];
  } else if (k == 1) {
    loop_ub = Dfilter->size[1];
  } else if (b_Cir2 == 1) {
    loop_ub = Cfilter->size[1];
  } else if (Cir->size[1] == 1) {
    loop_ub = Gir->size[1];
  } else {
    loop_ub = Cir->size[1];
  }
  if ((Cir->size[1] == Gir->size[1]) && (nx == Cfilter->size[1]) &&
      (d_Cir2 == Dfilter->size[1]) && (Cir2 == Efilter->size[1]) &&
      (loop_ub == c_Cir2->size[1])) {
    d = 0.5 * factor;
    high_C = 0.5 * C;
    F = 0.5 * D;
    D = 0.5 * E;
    C = 0.5 * G;
    factor = 0.5 * b;
    Cir2 = Ffilter->size[0] * Ffilter->size[1];
    Ffilter->size[0] = 1;
    Ffilter->size[1] = Cir->size[1];
    emxEnsureCapacity_real_T(Ffilter, Cir2);
    Ffilter_data = Ffilter->data;
    loop_ub = Cir->size[1];
    for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
      Ffilter_data[Cir2] =
          (((((Cir_data[Cir2] > d) + (Gir_data[Cir2] > high_C)) +
             (Cfilter_data[Cir2] > F)) +
            (Dfilter_data[Cir2] > D)) +
           (Efilter_data[Cir2] > C)) +
          (Cir2_data[Cir2] > factor);
    }
    c_conv(Ffilter, Cir);
    Cir_data = Cir->data;
  } else {
    b_binary_expand_op(Cir, factor, Gir, C, Cfilter, D, Dfilter, E, Efilter, G,
                       c_Cir2, b);
    Cir_data = Cir->data;
  }
  emxFree_real_T(&c_Cir2);
  emxFree_real_T(&Efilter);
  emxFree_real_T(&Dfilter);
  emxInit_int32_T(&r, 2);
  /* ------------------------------------------------------------------Parameter
   */
  /*  filter out HF noise */
  eml_find(Cir, r);
  r1 = r->data;
  Cir2 = Gir->size[0] * Gir->size[1];
  Gir->size[0] = 1;
  Gir->size[1] = r->size[1];
  emxEnsureCapacity_real_T(Gir, Cir2);
  Gir_data = Gir->data;
  loop_ub = r->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    Gir_data[Cir2] = r1[Cir2];
  }
  emxFree_int32_T(&r);
  /* truelocs = zeros(1,length(locs)); */
  /* truelocs(1) = locs(1); */
  j = 2U;
  Cir2 = Gir->size[1];
  for (i = 0; i <= Cir2 - 2; i++) {
    b_Cir2 = (int)Gir_data[i + 1];
    if ((double)b_Cir2 - (double)(int)Gir_data[i] > 100.0) {
      Gir_data[(int)j - 1] = b_Cir2;
      j++;
    }
  }
  nx = Gir->size[1] - 1;
  k = 0;
  for (i = 0; i <= nx; i++) {
    if ((int)Gir_data[i] > 0) {
      k++;
    }
  }
  loop_ub = 0;
  for (i = 0; i <= nx; i++) {
    Cir2 = (int)Gir_data[i];
    if (Cir2 > 0) {
      Gir_data[loop_ub] = Cir2;
      loop_ub++;
    }
  }
  emxInit_real_T(&ri, 1);
  Cir2 = Gir->size[0] * Gir->size[1];
  Gir->size[0] = 1;
  Gir->size[1] = k;
  emxEnsureCapacity_real_T(Gir, Cir2);
  Gir_data = Gir->data;
  Cir2 = ri->size[0];
  ri->size[0] = Cir->size[1];
  emxEnsureCapacity_real_T(ri, Cir2);
  timeSig_data = ri->data;
  loop_ub = Cir->size[1];
  for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
    timeSig_data[Cir2] = 0.0;
  }
  for (Cir2 = 0; Cir2 < k; Cir2++) {
    d = Gir_data[Cir2];
    timeSig_data[(int)d - 1] = Cir_data[(int)d - 1];
  }
  /*      figure() */
  /*      plot(rest_indicator) */
  /*      ylim([-1,2]); */
  /*      figure() */
  /*      plot(ri) */
  /*      ylim([-1,2]); */
  /*     %% Find the Tempo */
  /* coder.opaque(uint32); */
  j = 1U;
  memset(&peakHistory[0], 0, 30U * sizeof(unsigned int));
  memset(&b_tempo[0], 0, 30U * sizeof(double));
  i = 0;
  /* coder.varsize('upPeak'); */
  /* coder.varsize('upPeak2'); */
  emxInit_boolean_T(&b_ri);
  emxInit_int32_T(&r2, 1);
  do {
    exitg1 = 0;
    if (j > (unsigned int)ri->size[0]) {
      Cir2 = 0;
      b_Cir2 = 0;
    } else {
      Cir2 = (int)j - 1;
      b_Cir2 = ri->size[0];
    }
    loop_ub = b_Cir2 - Cir2;
    b_Cir2 = b_ri->size[0];
    b_ri->size[0] = loop_ub;
    emxEnsureCapacity_boolean_T(b_ri, b_Cir2);
    ri_data = b_ri->data;
    for (b_Cir2 = 0; b_Cir2 < loop_ub; b_Cir2++) {
      ri_data[b_Cir2] = (timeSig_data[Cir2 + b_Cir2] == 1.0);
    }
    b_eml_find(b_ri, r2);
    r1 = r2->data;
    Cir2 = r1[0];
    if (Cir2 < 0) {
      Cir2 = 0;
    }
    qY = j + Cir2;
    if (qY < j) {
      qY = MAX_uint32_T;
    }
    j = qY - j;
    if (j > qY) {
      j = 0U;
    }
    peakHistory[i] = j;
    j = qY;
    if (qY > (unsigned int)ri->size[0]) {
      Cir2 = 0;
      b_Cir2 = 0;
    } else {
      Cir2 = (int)qY - 1;
      b_Cir2 = ri->size[0];
    }
    loop_ub = b_Cir2 - Cir2;
    b_Cir2 = b_ri->size[0];
    b_ri->size[0] = loop_ub;
    emxEnsureCapacity_boolean_T(b_ri, b_Cir2);
    ri_data = b_ri->data;
    for (b_Cir2 = 0; b_Cir2 < loop_ub; b_Cir2++) {
      ri_data[b_Cir2] = (timeSig_data[Cir2 + b_Cir2] == -1.0);
    }
    c_eml_find(b_ri, (int *)&nx, &loop_ub);
    for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
      factor = (double)nx / fs_Hz;
    }
    if (loop_ub != 0) {
      b_tempo[i] = factor;
      i++;
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);
  emxFree_int32_T(&r2);
  emxFree_boolean_T(&b_ri);
  emxFree_real_T(&ri);
  d_Cir2 = 0;
  loop_ub = 0;
  for (i = 0; i < 30; i++) {
    j = peakHistory[i];
    if (j > 0U) {
      d_Cir2++;
      peakHistory_data[loop_ub] = j;
      loop_ub++;
    }
  }
  /*      disp("peak History: ") */
  /*      disp(peakHistory(2:end)) */
  /*      disp("tempo guesses: ") */
  /*      disp(tempo(tempo>0.01)) */
  /* disp("ave tempo guess: ") */
  /* disp(mean(tempo(tempo>0.01))) */
  /*     %% linear regression for slope estimate  */
  /*  if Octive is less than 4 filter the peakHistory signal to reduce high */
  /*  frequency noise from signal edge detection */
  if (Octive < 4.0) {
    /*          ri_spikes = find(rest_indicator); */
    /*          true_spikes = zeros(1,length(ri_spikes)); */
    /*          true_spikes(1) = ri_spikes(1); */
    /*          thresh = 100; */
    /*          j = 1; */
    /*          for i = 2:length(ri_spikes) */
    /*              if((ri_spikes(i)-ri_spikes(i-1))>thresh) */
    /*                  j=j+1; */
    /*                  true_spikes(j) = ri_spikes(i); */
    /*              end */
    /*          end     */
    if (1 > k) {
      Cir2 = 1;
      b_Cir2 = -1;
    } else {
      Cir2 = 2;
      b_Cir2 = k - 1;
    }
    i = div_s32_floor(b_Cir2, Cir2);
    if (2 > i + 1) {
      nx = 0;
      b_Cir2 = 0;
    } else {
      nx = 1;
      b_Cir2 = div_s32_floor(b_Cir2, Cir2) + 1;
    }
    if (1 > i) {
      i = 0;
    }
    loop_ub = b_Cir2 - nx;
    if (loop_ub == i) {
      b_Cir2 = Cfilter->size[0] * Cfilter->size[1];
      Cfilter->size[0] = 1;
      Cfilter->size[1] = loop_ub;
      emxEnsureCapacity_real_T(Cfilter, b_Cir2);
      Cfilter_data = Cfilter->data;
      for (b_Cir2 = 0; b_Cir2 < loop_ub; b_Cir2++) {
        Cfilter_data[b_Cir2] =
            (int)(unsigned int)Gir_data[Cir2 * (nx + b_Cir2)] -
            (int)(unsigned int)Gir_data[Cir2 * b_Cir2];
      }
    } else {
      binary_expand_op(Cfilter, Gir, Cir2, nx, b_Cir2 - 1, i - 1);
      Cfilter_data = Cfilter->data;
    }
    nx = Cfilter->size[1] - 1;
    k = 0;
    for (i = 0; i <= nx; i++) {
      if ((int)Cfilter_data[i] > 0) {
        k++;
      }
    }
    loop_ub = 0;
    for (i = 0; i <= nx; i++) {
      Cir2 = (int)Cfilter_data[i];
      if (Cir2 > 0) {
        Cfilter_data[loop_ub] = Cir2;
        loop_ub++;
      }
    }
    Cir2 = Cfilter->size[0] * Cfilter->size[1];
    Cfilter->size[0] = 1;
    Cfilter->size[1] = k;
    emxEnsureCapacity_real_T(Cfilter, Cir2);
    Cfilter_data = Cfilter->data;
    d_Cir2 = 0;
    loop_ub = 0;
    for (i = 0; i < 30; i++) {
      if (b_tempo[i] > 0.01) {
        d_Cir2++;
        tmp_data[loop_ub] = (signed char)(i + 1);
        loop_ub++;
      }
    }
    Cir2 = Cir->size[0] * Cir->size[1];
    Cir->size[0] = 1;
    Cir->size[1] = k;
    emxEnsureCapacity_real_T(Cir, Cir2);
    Cir_data = Cir->data;
    tempo_size[0] = 1;
    tempo_size[1] = d_Cir2;
    for (Cir2 = 0; Cir2 < d_Cir2; Cir2++) {
      tempo_data[Cir2] = b_tempo[tmp_data[Cir2] - 1];
    }
    d = mean(tempo_data, tempo_size) * fs_Hz;
    for (Cir2 = 0; Cir2 < k; Cir2++) {
      Cir_data[Cir2] = Cfilter_data[Cir2] / d;
    }
    nx = Cir->size[1];
    for (k = 0; k < nx; k++) {
      Cir_data[k] = rt_roundd_snf(Cir_data[k]);
    }
    nx = Cir->size[1] - 1;
    k = 0;
    for (i = 0; i <= nx; i++) {
      if (Cir_data[i] < 5.0) {
        k++;
      }
    }
    emxInit_int32_T(&r3, 2);
    Cir2 = r3->size[0] * r3->size[1];
    r3->size[0] = 1;
    r3->size[1] = k;
    emxEnsureCapacity_int32_T(r3, Cir2);
    r1 = r3->data;
    loop_ub = 0;
    for (i = 0; i <= nx; i++) {
      if (Cir_data[i] < 5.0) {
        r1[loop_ub] = i + 1;
        loop_ub++;
      }
    }
    /* disp("tempo guess from slope: ") */
    /* disp(coefs(1)/fs_Hz) */
    Cir2 = Ffilter->size[0] * Ffilter->size[1];
    Ffilter->size[0] = 1;
    Ffilter->size[1] = r3->size[1];
    emxEnsureCapacity_real_T(Ffilter, Cir2);
    Ffilter_data = Ffilter->data;
    loop_ub = r3->size[1];
    for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
      Ffilter_data[Cir2] = Cir_data[r1[Cir2] - 1];
    }
    Cir2 = timeSig->size[0] * timeSig->size[1];
    timeSig->size[0] = 1;
    timeSig->size[1] = r3->size[1];
    emxEnsureCapacity_real_T(timeSig, Cir2);
    timeSig_data = timeSig->data;
    loop_ub = r3->size[1];
    for (Cir2 = 0; Cir2 < loop_ub; Cir2++) {
      timeSig_data[Cir2] = (int)Cfilter_data[r1[Cir2] - 1];
    }
    emxFree_int32_T(&r3);
    polyfit(Ffilter, timeSig, dv);
    tempo = dv[0] / fs_Hz;
  } else {
    if (2 > d_Cir2) {
      Cir2 = 0;
      b_Cir2 = 0;
    } else {
      Cir2 = 1;
      b_Cir2 = d_Cir2;
    }
    k = 0;
    loop_ub = 0;
    for (i = 0; i < 30; i++) {
      if (b_tempo[i] > 0.01) {
        k++;
        tmp_data[loop_ub] = (signed char)(i + 1);
        loop_ub++;
      }
    }
    tempo_size[0] = 1;
    tempo_size[1] = k;
    for (i = 0; i < k; i++) {
      tempo_data[i] = b_tempo[tmp_data[i] - 1];
    }
    d = mean(tempo_data, tempo_size) * fs_Hz;
    i = Cir->size[0] * Cir->size[1];
    Cir->size[0] = 1;
    loop_ub = b_Cir2 - Cir2;
    Cir->size[1] = loop_ub;
    emxEnsureCapacity_real_T(Cir, i);
    Cir_data = Cir->data;
    for (b_Cir2 = 0; b_Cir2 < loop_ub; b_Cir2++) {
      Cir_data[b_Cir2] = (double)peakHistory_data[Cir2 + b_Cir2] / d;
    }
    nx = Cir->size[1];
    for (k = 0; k < nx; k++) {
      Cir_data[k] = rt_roundd_snf(Cir_data[k]);
    }
    if (1 > Cir->size[1] - 1) {
      loop_ub = 0;
    } else {
      loop_ub = Cir->size[1] - 1;
    }
    if (2 > d_Cir2 - 1) {
      Cir2 = 0;
      b_Cir2 = 0;
    } else {
      Cir2 = 1;
      b_Cir2 = d_Cir2 - 1;
    }
    /* disp("tempo guess from slope: ") */
    /* disp(coefs(1)/fs_Hz) */
    Cir_size[0] = 1;
    Cir_size[1] = loop_ub;
    for (i = 0; i < loop_ub; i++) {
      b_Cir_data[i] = Cir_data[i];
    }
    peakHistory_size[0] = 1;
    loop_ub = b_Cir2 - Cir2;
    peakHistory_size[1] = loop_ub;
    for (b_Cir2 = 0; b_Cir2 < loop_ub; b_Cir2++) {
      b_peakHistory_data[b_Cir2] = peakHistory_data[Cir2 + b_Cir2];
    }
    c_Cir_data.data = &b_Cir_data[0];
    c_Cir_data.size = &Cir_size[0];
    c_Cir_data.allocatedSize = 29;
    c_Cir_data.numDimensions = 2;
    c_Cir_data.canFreeData = false;
    c_peakHistory_data.data = &b_peakHistory_data[0];
    c_peakHistory_data.size = &peakHistory_size[0];
    c_peakHistory_data.allocatedSize = 29;
    c_peakHistory_data.numDimensions = 2;
    c_peakHistory_data.canFreeData = false;
    polyfit(&c_Cir_data, &c_peakHistory_data, dv);
    tempo = dv[0] / fs_Hz;
  }
  emxFree_real_T(&Gir);
  emxFree_real_T(&Cir);
  emxFree_real_T(&Ffilter);
  emxFree_real_T(&Cfilter);
  emxFree_real_T(&timeSig);
  /*     %% Sync the time */
  /* actual_time_offset_s = time_s(time_offset) %print out the actual start time
   * to compare. */
  /* toc */
  return tempo;
}

/* End of code generation (DetectorFastArduino.c) */
