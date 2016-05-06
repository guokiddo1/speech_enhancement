/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: CompCepsExports.h
 * PURPOSE:   Calculating cepstral features and logE for the de-noised
 *            input frame.
 *
 *-------------------------------------------------------------------------------*/
#ifndef _COMPCEPSEXPORTS_H
#define _COMPCEPSEXPORTS_H

extern CompCepsStructX *DoCompCepsAlloc(void);
extern void DoCompCepsInit(FEParamsX * This);
extern BOOLEAN DoCompCeps(X_FLOAT32 * Data, X_FLOAT32 * Coef, FEParamsX * This);
extern void DoCompCepsDelete(CompCepsStructX * WPX);

#endif
