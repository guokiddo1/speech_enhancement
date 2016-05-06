/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: WaveProcExports.h
 * PURPOSE: Apply SNR-dependent Waveform Processing on the de-noised frame.
 *
 *-------------------------------------------------------------------------------*/
#ifndef _WAVEPROCEXPORTS_H
#define _WAVEPROCEXPORTS_H

extern WaveProcStructX *DoWaveProcAlloc(void);
extern void DoWaveProcInit(FEParamsX * This);
extern BOOLEAN DoWaveProc(X_FLOAT32 * Data, FEParamsX * This);
extern void DoWaveProcDelete(WaveProcStructX * WPX);

#endif
