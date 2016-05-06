/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: MelProcExports.h
 * PURPOSE: 1) Mel-filtering of FFT spectrum using the Mel filter bank
 *             frequency windows.
 *          2) Mel inverse DCT applied on the output of Mel filter bank.
 *
 *-------------------------------------------------------------------------------*/
#ifndef _MELPROCEXPORTS_H
#define _MELPROCEXPORTS_H

#define WF_MEL_ORDER        25

/* Structure for FFT window (one triangle in the chained list) */
struct MelFB_Window {
  int StartingPoint;
  int Length;
  float *Data;
  struct MelFB_Window *Next;
};


MelFB_Window *CMelFBAlloc();

void InitMelFBwindows(MelFB_Window * FirstWin, float StFreq, float SmplFreq, int FFTLength, int NumChannels, BOOLEAN normalize);
void 
InitFFTWindows(MelFB_Window * FirstWin,
	       float StFreq,
	       float SmplFreq,
	       int FFTLength,
	       int NumChannels);
void ComputeTriangle(MelFB_Window * FirstWin);

void DoMelFB(float *SigFFT, MelFB_Window * FirstWin);

void ReleaseMelFBwindows(MelFB_Window * FirstWin);

void DoMelIDCT(float *inData, float **melIDCTbasis, int melOrder, int timeLength);
void InitMelIDCTbasis(float **melIDCTbasis, MelFB_Window * FirstWin, short melorder, int sampFreq, int FFTLength);

#endif
