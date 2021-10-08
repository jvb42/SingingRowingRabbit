/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * conv.h
 *
 * Code generation for function 'conv'
 *
 */

#ifndef CONV_H
#define CONV_H

/* Include files */
#include "DetectorFastArduino_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_conv(const double A[106501], const emxArray_real_T *B,
            emxArray_real_T *C);

void c_conv(const emxArray_real_T *A, emxArray_real_T *C);

void conv(const emxArray_real_T *A, const emxArray_real_T *B,
          emxArray_real_T *C);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (conv.h) */
