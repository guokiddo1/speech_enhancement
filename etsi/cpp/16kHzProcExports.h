/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: 16kHzProcExports.h
 * PURPOSE:   For processing of 16 kHz sampled input signal.
 *
 *-------------------------------------------------------------------------------*/
#ifndef _16kHzPROCESSING_H
#define _16kHzPROCESSING_H

/*-------------
 * definitions
 *-------------*/
#define HP16k_MEL_ORDER       5    //
#define HP16k_MEL_USED        3    // must be (HP16k_MEL_ORDER - 2) because both first and last halfbands are not used!
#define USE_INTERP_BAND       0    // use an interpolation band between bands 23 and 24 -> 1=yes, 0=no
#define DO_INTERP             1    // preform interpolation on bands 23 and 24 -> 1=yes, 0=no
#define NB_LP_BANDS_CODING    3    //
#define ENERGY_CORR           1    // correct log energy by using the HP energy
#define PERC_CODED            0.7  //

#define NE16k_FRAMES_THRESH 100    //
#define LAMBDA_NSE16k         0.99 //

/*------------
 * structures
 *------------*/
typedef struct QMF_FIR QMF_FIR;

/*------------------
 * public functions
 *------------------*/
extern QMF_FIR *QMF_FIR_Init();
extern void hq_free(QMF_FIR * fir_ptr);

extern void Do16kProcDelete(DataFor16kProc * pData16k);
extern void Do16kProcInit(FEParamsX * This);
extern DataFor16kProc *Do16kProcAlloc();
extern void 
Do16kProcessing(float *inData,
		DataFor16kProc * pData16k,
		int inDataLength);
extern void 
GetBandsForCoding16k(float *inData,
		     float *bandsForCoding16k,
		     int inDataSize);
extern void 
CodeBands16k(float *fb16k,
	     float *lpBands,
	     float *codeForBands16k);
extern void 
GetBandsForDecoding16k(float *inData,
		       float *bandsForCoding16k,
		       int inDataSize);
extern void 
DecodeBands16k(float *fb16k,
	       float *lpBands,
	       float *codeForBands16k,
	       float *codeWeights);
extern void 
DoSpecSub16k(float *inFB16k,
	     DataFor16kProc * pData16k,
	     long vadCounter16k,
	     int inFB16kSize);
extern void 
MergeSSandCoded(float *FloatBuffer,
		float *hp16kBands,
		int numChannelsWI7,
		X_INT16 * hp16kBandsSize_in,
		DataFor16kProc * pData16k);
extern float 
CorrectEnergy(float LogEnergy,
	      int hp16kBandsSize,
	      float *FloatBuffer,
	      int numChannelsWI7,
	      float preem);

/*----------------------
 * public get functions
 *----------------------*/
extern int Get16k_bufData16kSize(DataFor16kProc * pData16k);
extern int Get16k_hpBandsSize(DataFor16kProc * pData16k);
extern float Get16k_percCoded(DataFor16kProc * pData16k);
extern float 
Get16k_dataHP(DataFor16kProc * pData16k,
	      int i);
extern float *Get16k_p_hpBands(DataFor16kProc * pData16k);
extern float *Get16k_p_bufferCodeForBands16k(DataFor16kProc * pData16k);
extern float *Get16k_p_CodeForBands16k(DataFor16kProc * pData16k);
extern float *Get16k_p_bufferCodeWeights(DataFor16kProc * pData16k);
extern float *Get16k_p_codeWeights(DataFor16kProc * pData16k);
extern float *Get16k_p_bufferData16k(DataFor16kProc * pData16k);
extern float *Get16k_p_BandsForCoding16k(DataFor16kProc * pData16k);
extern MelFB_Window *Get16k_p_FirstWindow16k(DataFor16kProc * pData16k);

#endif
