/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: ParmInterface.h
 * PURPOSE: Processing one input frame in DoAdvProcess():
 *          DoNoiseSup(), DoWaveProc(), DoCompCeps(), DoPostProc(),
 *          and DoVADProc() are called.
 *
 *-------------------------------------------------------------------------------*/
#ifndef _PARMINTERFACE_H
#define _PARMINTERFACE_H

/*-----------------
 * File Inclusions
 *-----------------*/
#include "BufferIn.h"
#include "ParmType.h"

/*------------------------------------
 * Definition of Front-end Parameters
 *------------------------------------*/
#define SAMPLING_FREQ_1    8   // 8kHz
#define SAMPLING_FREQ_2   11   // 11kHz
#define SAMPLING_FREQ_3   16   // 16kHz

#define FRAME_LENGTH     200   // 25ms

#define FRAME_SHIFT       80   // 10ms

#define FFT_LENGTH       256
#define FRAME_BUF_SIZE   241

#define STARTING_FREQ     64.0 // 55.401825

#define NUM_CEP_COEFF     13   //c1...c12 + c0

#define NON_SPEECH_FRAME   0
#define SPEECH_FRAME       1

#define PIx2               6.28318530717958647692


/*------------------------------------------------------
 * Private structures for the modules (noise reduction,
 * SNR waveform processing, cepstrum calculation, ...)
 *------------------------------------------------------*/
typedef struct NoiseSupStructX NoiseSupStructX;
typedef struct WaveProcStructX WaveProcStructX;
typedef struct CompCepsStructX CompCepsStructX;
typedef struct PostProcStructX PostProcStructX;
typedef struct DataFor16kProc DataFor16kProc;
typedef struct MelFB_Window MelFB_Window;
typedef struct VADStructX VADStructX;

/*------------------------------------------------------
 * FEParamsX must (only) contains:
 *  - initialisation parameters
 *    (like FrameLength, FrameShift, FFTLength)
 *  - buffers to transmit data between modules
 *    (like denoisedBuf)
 *  - specific information which must be remanent
 *    between two frames, and which can't be
 *    encapsulated in modules
 *  - if a data is not remanent, they must be transmitted
 *    between modules by a GET ans a SET function
 *  - if a data can be encapsulated in a module, it must
 *    be encapsulated to avoid conflict between modules
 *------------------------------------------------------*/
typedef struct FEParamsX FEParamsX;

struct FEParamsX {
  BOOLEAN Do16kHzProc;
  BOOLEAN Noc0;

  int VAD;
  int CoefNb;
  int FFTLength;
  int FrameShift;
  int FrameLength;
  int FrameCounter;
  int SpeechFoundMel;
  int SpeechFoundVar;
  int SpeechFoundSpec;
  int SpeechFoundVADNS;
  int NbSamplesToRead;
  int SamplingFrequency;
  int offsetDenoisedFrame;

  float StartingFrequency;

  long NonZeroFrameOnset;
  long ZeroFrameCounter;

  BufferIn * denoisedBuf;
  X_FLOAT32 * CurFrame;
  DataFor16kProc * pData16k;

  /*---------------------------
   * module specific expansion
   *---------------------------*/
  NoiseSupStructX * NSX;
  WaveProcStructX *WPX;
  CompCepsStructX *CCX;
  PostProcStructX *PPX;
  VADStructX *VADX;

  /*-------------------------------------
   * methods: use pointer to function to
   *          emulate generic function
   *-------------------------------------*/

  /*-----------------------
   * for noise suppression
   *-----------------------*/
  NoiseSupStructX *(*DoNoiseSupAlloc) (void);
  void (*DoNoiseSupInit) (FEParamsX * that);
      BOOLEAN(*DoNoiseSup) (X_FLOAT32 * InData,
			        X_FLOAT32 * outData,
			        FEParamsX * that);
  void (*DoNoiseSupDelete) (NoiseSupStructX *);

  /*-------------------------
   * for waveform processing
   *-------------------------*/
  WaveProcStructX *(*DoWaveProcAlloc) (void);
  void (*DoWaveProcInit) (FEParamsX * that);
      BOOLEAN(*DoWaveProc) (X_FLOAT32 * Data,
			        FEParamsX * that);
  void (*DoWaveProcDelete) (WaveProcStructX *);

  /*--------------------------
   * for cepstrum calculation
   *--------------------------*/
  CompCepsStructX *(*DoCompCepsAlloc) (void);
  void (*DoCompCepsInit) (FEParamsX * that);
      BOOLEAN(*DoCompCeps) (X_FLOAT32 * Data,
			        X_FLOAT32 * Coef,
			        FEParamsX * that);
  void (*DoCompCepsDelete) (CompCepsStructX *);

  /*---------------------
   * for post processing
   *---------------------*/
  PostProcStructX *(*DoPostProcAlloc) (void);
  void (*DoPostProcInit) (FEParamsX * that);
      BOOLEAN(*DoPostProc) (X_FLOAT32 * Coef,
			        FEParamsX * that);
  void (*DoPostProcDelete) (PostProcStructX *);

  /*----------------------
   * for 16kHz Processing
   *----------------------*/
  DataFor16kProc *(*Do16kProcAlloc) (void);
  void (*Do16kProcInit) (FEParamsX * that);
  void (*Do16kProcessing) (float *inData,
			       DataFor16kProc * pData16k,
			       int inDataLength);
  void (*Do16kProcDelete) (DataFor16kProc * pData16k);

  /*---------
   * for VAD
   *---------*/
  VADStructX *(*DoVADAlloc) (void);
  void (*DoVADInit) (FEParamsX * that);
      BOOLEAN(*DoVADProc) (X_FLOAT32 * FeatureBuffer,
			       FEParamsX * that);
  void (*DoVADDelete) (VADStructX * VADX);

  /*-------------------
   * for miscellaneous
   *-------------------*/
      BOOLEAN(*DoVADFlush) (X_FLOAT32 * FeatureBuffer,
			        FEParamsX * that);

};

/*----------------------------------------------------
 * Global routines : each call the four stages of the
 * front end (noise suppression, wave processing,
 * cepstrum calculation and post processing), and
 * finally call the VAD stage
 *----------------------------------------------------*/

/*------------
 * allocation
 *------------*/
extern FEParamsX *AdvProcessAlloc(int samplingFrequency);

/*----------------
 * initialisation
 *----------------*/
extern void AdvProcessInit(FEParamsX * pComParX);

/*-------------------------------------------------------------
 * Processing "FrameShift" samples
 * Return TRUE if a frame is output, FALSE if no (because some
 * modules may have some latency, hence at the beginning there
 * is no output data for the first input frames
 *-------------------------------------------------------------*/
extern BOOLEAN 
DoAdvProcess(FILE_TYPE * CurrentFrame,
             FILE_TYPE * DenoiseBuffer,
             X_FLOAT32 * FeatureBuffer,
             FEParamsX * pComParX);

/*-------------------------------------------------------------
 * Flushing: (at the end of a record).
 * Return TRUE if a frame is flushed,
 * FALSE if there is no more frame.
 *-------------------------------------------------------------*/
extern BOOLEAN 
FlushAdvProcess(X_FLOAT32 * FeatureBuffer,
		FEParamsX * pComParX);

/*----------
 * deletion
 *----------*/
extern void AdvProcessDelete(FEParamsX **);

#endif
