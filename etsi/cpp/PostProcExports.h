/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: PostProcExports.h
 * PURPOSE: Apply post-processing (blind equalization) on a frame of
 *          cepstral coefficients.
 *
 *-------------------------------------------------------------------------------*/
#ifndef _POSTPROCEXPORTS_H
#define _POSTPROCEXPORTS_H

/* structure allocation */
extern PostProcStructX *DoPostProcAlloc(void);

/* Reset of Post Processing */
extern void DoPostProcInit(FEParamsX * This);

/* Do post proccessing for the current cepstral coefficient frame */
extern BOOLEAN DoPostProc(X_FLOAT32 * Data, FEParamsX * This);

/* structure freeing */
extern void DoPostProcDelete(PostProcStructX * WPX);

#endif
