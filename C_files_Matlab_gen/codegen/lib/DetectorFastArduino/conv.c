/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * conv.c
 *
 * Code generation for function 'conv'
 *
 */

/* Include files */
#include "conv.h"
#include "DetectorFastArduino_emxutil.h"
#include "DetectorFastArduino_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void b_conv(const double A[106501], const emxArray_real_T *B,
            emxArray_real_T *C)
{
  const double *B_data;
  double *C_data;
  int b_k;
  int k;
  int nB;
  int nC;
  B_data = B->data;
  nB = B->size[1] - 1;
  if (B->size[1] == 0) {
    nC = 106501;
  } else {
    nC = B->size[1] + 106500;
  }
  k = C->size[0] * C->size[1];
  C->size[0] = 1;
  C->size[1] = nC;
  emxEnsureCapacity_real_T(C, k);
  C_data = C->data;
  for (k = 0; k < nC; k++) {
    C_data[k] = 0.0;
  }
  if (B->size[1] > 0) {
    if (B->size[1] > 106501) {
      for (k = 0; k < 106501; k++) {
        for (b_k = 0; b_k <= nB; b_k++) {
          nC = k + b_k;
          C_data[nC] += A[k] * B_data[b_k];
        }
      }
    } else {
      for (k = 0; k <= nB; k++) {
        for (b_k = 0; b_k < 106501; b_k++) {
          nC = k + b_k;
          C_data[nC] += B_data[k] * A[b_k];
        }
      }
    }
  }
}

void c_conv(const emxArray_real_T *A, emxArray_real_T *C)
{
  const double *A_data;
  double *C_data;
  int b_k;
  int i;
  int k;
  int nC;
  A_data = A->data;
  if (A->size[1] == 0) {
    nC = 2;
  } else {
    nC = A->size[1] + 1;
  }
  i = C->size[0] * C->size[1];
  C->size[0] = 1;
  C->size[1] = nC;
  emxEnsureCapacity_real_T(C, i);
  C_data = C->data;
  for (i = 0; i < nC; i++) {
    C_data[i] = 0.0;
  }
  if (A->size[1] > 0) {
    if (2 > A->size[1]) {
      C_data[0] = -A_data[0];
      C_data[1] += A_data[0];
    } else {
      i = A->size[1] - 1;
      for (k = 0; k < 2; k++) {
        for (b_k = 0; b_k <= i; b_k++) {
          nC = k + b_k;
          C_data[nC] += (2.0 * (((double)k + 1.0) - 1.0) + -1.0) * A_data[b_k];
        }
      }
    }
  }
}

void conv(const emxArray_real_T *A, const emxArray_real_T *B,
          emxArray_real_T *C)
{
  const double *A_data;
  const double *B_data;
  double *C_data;
  int b_k;
  int k;
  int nA;
  int nApnB;
  int nB;
  B_data = B->data;
  A_data = A->data;
  nA = A->size[1] - 1;
  nB = B->size[1] - 1;
  nApnB = A->size[1] + B->size[1];
  if ((A->size[1] != 0) && (B->size[1] != 0)) {
    nApnB--;
  }
  k = C->size[0] * C->size[1];
  C->size[0] = 1;
  C->size[1] = nApnB;
  emxEnsureCapacity_real_T(C, k);
  C_data = C->data;
  for (k = 0; k < nApnB; k++) {
    C_data[k] = 0.0;
  }
  if ((A->size[1] > 0) && (B->size[1] > 0)) {
    if (B->size[1] > A->size[1]) {
      for (k = 0; k <= nA; k++) {
        for (b_k = 0; b_k <= nB; b_k++) {
          nApnB = k + b_k;
          C_data[nApnB] += A_data[k] * B_data[b_k];
        }
      }
    } else {
      for (k = 0; k <= nB; k++) {
        for (b_k = 0; b_k <= nA; b_k++) {
          nApnB = k + b_k;
          C_data[nApnB] += B_data[k] * A_data[b_k];
        }
      }
    }
  }
}

/* End of code generation (conv.c) */
