/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: VADExports.h
 * PURPOSE: Logic for the frame dropping VAD. VAD flag is set
 *          according to parameters (SpeechFoundSpec, SpeechFoundMel,
 *          SpeechFoundVar) calculated in NoiseSup.c.
 *
 *-------------------------------------------------------------------------------*/
#ifndef _VADEXPORTS_H
#define _VADEXPORTS_H

extern VADStructX *DoVADAlloc(void);
extern void DoVADInit(FEParamsX * This);
extern BOOLEAN DoVADProc(X_FLOAT32 * FeatureBuffer, FEParamsX * This);
extern BOOLEAN DoVADFlush(X_FLOAT32 * FeatureBuffer, FEParamsX * This);
extern void DoVADDelete(VADStructX * VADX);

#endif
