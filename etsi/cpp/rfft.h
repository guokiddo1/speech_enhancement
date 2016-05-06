/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: rfft.h
 * PURPOSE: Real valued, in-place split-radix FFT.
 *
 *-------------------------------------------------------------------------------*/
#ifndef _RFFT_H
#define _RFFT_H

#define M_PI        3.14159265358979323846
#define M_SQRT2     1.41421356237309504880

void rfft(float *x, int n, int m);

#endif
