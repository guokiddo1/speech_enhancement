/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: CompCeps.c
 * PURPOSE:   Calculating cepstral features and logE for the de-noised
 *            input frame.
 *
 *-------------------------------------------------------------------------------*/
/*-----------------
 * File Inclusions
 *-----------------*/
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "ParmInterface.h"
#include "CompCepsExports.h"
#include "16kHzProcExports.h"
#include "MelProcExports.h"
#include "rfft.h"

/*----------------------
 * Constant Definitions
 *----------------------*/
#define PI   3.14159265358979323846
#define PIx2 6.28318530717958647692

/*-------------------------------
 * Global Definitions and Macros
 *-------------------------------*/
#define CC_NUM_CHANNELS_WI8  23
#define CC_PRE_EMPHASIS       0.90
#define CC_ENERGYFLOOR_FB   -10.0
#define CC_ENERGYFLOOR_logE -50.0
#define CC_FRAME_LENGTH     200
#define CC_FFT_LENGTH       256
#define CC_FFT_LENGTH_ORDER   8

struct CompCepsStructX {
  X_INT16 SamplingFrequency;
  X_FLOAT32 StartingFrequency;
  X_FLOAT32 HammingWindow[CC_FRAME_LENGTH];
  X_FLOAT32 * pDCTMatrix;
  MelFB_Window * FirstWindow;
  BOOLEAN Noc0;
  BOOLEAN Do16kHzProc;
  DataFor16kProc * pData16k;
};

/*------------
 * Prototypes
 *------------*/
static void 
WI8CompCeps(CompCepsStructX * CCX, X_FLOAT32 * curFrameFloatBuf,
	    X_FLOAT32 * cepOut, DataFor16kProc * pData16k);
static void InitializeHamming(X_FLOAT32 * win, X_INT16 len);
static void Window(X_FLOAT32 * data, X_FLOAT32 * win, X_INT16 len);
static void DCT(X_FLOAT32 * Data, X_FLOAT32 * Mx, X_INT16 NumCepstralCoeff, X_INT16 NumChannels);
static X_FLOAT32 *InitDCTMatrix(X_INT16 NumCepstralCoeff, X_INT16 NumChannels);

/*-----------
 * Functions
 *-----------*/
/*---------------------------------------------------------------------------
 * FUNCTION NAME: InitializeHamming
 *
 * PURPOSE:       Initializes Hamming window coefficients
 *
 * INPUT:
 *   win          Pointer to window buffer
 *   len          Window length
 *
 * OUTPUT
 *                Hamming window coefficients stored in window buffer pointed
 *                to by *win*
 *
 * RETURN VALUE
 *   none
 *
 *---------------------------------------------------------------------------*/
static void 
InitializeHamming(X_FLOAT32 * win, X_INT16 len)
{
  X_INT16 i;

  for (i = 0; i < len / 2; i++)
    win[i] = 0.54 - 0.46 * cos(PIx2 * (i + 0.5) / len);
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: Window
 *
 * PURPOSE:       Performs windowing on input speech frame (multiplies input
 *                samples by the corresponding window coefficients)
 *
 * INPUT:
 *   data         Pointer to input speech buffer
 *   win          Pointer to window buffer
 *   len          Window (or frame) length
 *
 * OUTPUT
 *                Windowed speech frame stored at the same place as the
 *                original speech samples (pointed by *data*)
 *
 * RETURN VALUE
 *   none
 *
 *---------------------------------------------------------------------------*/
static void 
Window(X_FLOAT32 * data, X_FLOAT32 * win, X_INT16 len)
{
  X_INT16 i;

  for (i = 0; i < len / 2; i++)
    data[i] *= win[i];

  for (i = len / 2; i < len; i++)
    data[i] *= win[len - 1 - i];
}



/*---------------------------------------------------------------------------
 * FUNCTION NAME: InitDCTMatrix
 *
 * PURPOSE:       Initializes matrix for DCT computation (DCT is implemented
 *                as matrix-vector multiplication). The DCT matrix is of size
 *                (NumCepstralCoeff-1)-by-NumChannels. The zeroth cepstral
 *                coefficient is computed separately (needing NumChannels
 *                additions and only one multiplication), so the zeroth row
 *                of DCT matrix corresponds to the first DCT basis vector, the
 *                first one to the second one, and so on up to
 *                NumCepstralCoeff-1.
 *
 * INPUT:
 *   NumCepstralCoeff
 *                Number of cepstral coeffficients
 *   NumChannels  Number of filter bank channels
 *
 * OUTPUT
 *   none
 *
 * RETURN VALUE
 *                Pointer to the initialized DCT matrix
 *
 *---------------------------------------------------------------------------*/
static X_FLOAT32 *
InitDCTMatrix(X_INT16 NumCepstralCoeff, X_INT16 NumChannels)
{
  X_INT16 i, j;
  X_FLOAT32 *Mx;

  /*----------------------------------
   * Allocating memory for DCT-matrix
   *----------------------------------*/
  Mx = (X_FLOAT32 *) malloc(sizeof(X_FLOAT32) * (NumCepstralCoeff - 1) *
			    NumChannels);

  /*--------------------------
   * Computing matrix entries
   *--------------------------*/
  for (i = 1; i < NumCepstralCoeff; i++)
    for (j = 0; j < NumChannels; j++)
      Mx[(i - 1) * NumChannels + j] =
	cos(PI * (X_FLOAT32) i / (X_FLOAT32) NumChannels * ((X_FLOAT32) j + 0.5));
  return Mx;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: DCT
 *
 * PURPOSE:       Computes DCT transformation of filter bank outputs, results
 *                in cepstral coefficients. The DCT transformation is
 *                implemented as matrix-vector multiplication. The zeroth
 *                cepstral coefficient is computed separately and appended.
 *                Final cepstral coefficient order is c1, c2, ...,c12, c0. The
 *                output is stored right after the input values in the memory.
 *                Since the mel filter bank outputs are stored at the beginning
 *                of the FFT magnitude array it shouldn`t cause any problems.
 *                Some memory saving can be done this way.
 *
 * INPUT:
 *   Data         Pointer to input data buffer (filter bank outputs)
 *   Mx           Pointer to DCT matrix
 *   NumCepstralCoeff
 *                Number of cepstral coefficients
 *   NumChannels  Number of filter bank channels
 *
 * OUTPUT
 *                Cepstral coefficients stored after the input filter bank
 *                values pointed to by *Data*
 *
 * RETURN VALUE
 *   none
 *
 *---------------------------------------------------------------------------*/
static void 
DCT(X_FLOAT32 * Data, X_FLOAT32 * Mx, X_INT16 NumCepstralCoeff, X_INT16 NumChannels)
{
  X_INT16 i, j;

  /*-----------------------------------------------------
   * Computing c1..c/NumCepstralCoeff-1/, storing result
   * after the incoming data vector
   *-----------------------------------------------------*/
  for (i = 1; i < NumCepstralCoeff; i++) {
    Data[NumChannels + (i - 1)] = 0.0;
    for (j = 0; j < NumChannels; j++)
      Data[NumChannels + (i - 1)] +=
	Data[j] * Mx[(i - 1) * NumChannels + j];
  }

  /*----------------------------------------------------
   * Computing c0, as the last element of output vector
   *----------------------------------------------------*/
  Data[NumChannels + NumCepstralCoeff - 1] = 0.0;
  for (i = 0; i < NumChannels; i++)
    Data[NumChannels + NumCepstralCoeff - 1] += Data[i];

  return;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoCompCepsAlloc
 *
 * PURPOSE:
 *
 *
 * INPUT:
 *
 *
 * OUTPUT:
 *
 *
 * RETURN VALUE:
 *
 *
 *---------------------------------------------------------------------------*/
extern CompCepsStructX *
DoCompCepsAlloc(void)
{
  CompCepsStructX *CCX = NULL;
  CCX = (CompCepsStructX *) calloc(1, sizeof(CompCepsStructX));
  if (CCX == NULL)
    return NULL;
  return CCX;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoCompCepsInit
 *
 * PURPOSE:
 *
 *
 * INPUT:
 *
 *
 * OUTPUT:
 *
 *
 * RETURN VALUE:
 *
 *
 *---------------------------------------------------------------------------*/
extern void 
DoCompCepsInit(FEParamsX * This)
{
  CompCepsStructX *CCX = This->CCX;

  //16 kHz processing
    CCX->Do16kHzProc = This->Do16kHzProc;
  CCX->pData16k = This->pData16k;
  if (CCX->Do16kHzProc)
    CCX->pDCTMatrix = InitDCTMatrix(NUM_CEP_COEFF, CC_NUM_CHANNELS_WI8 + HP16k_MEL_USED);
  else
    CCX->pDCTMatrix = InitDCTMatrix(NUM_CEP_COEFF, CC_NUM_CHANNELS_WI8);
  CCX->Noc0 = This->Noc0;
  CCX->SamplingFrequency = This->SamplingFrequency;
  CCX->StartingFrequency = This->StartingFrequency;
  CCX->FirstWindow = CMelFBAlloc();
  InitFFTWindows(CCX->FirstWindow, CCX->StartingFrequency,
		 CCX->SamplingFrequency, CC_FFT_LENGTH, CC_NUM_CHANNELS_WI8);
  ComputeTriangle(CCX->FirstWindow);
  InitializeHamming(CCX->HammingWindow, CC_FRAME_LENGTH);
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoCompCeps
 *
 * PURPOSE: Cepstral computation for a frame
 *
 *
 * INPUT:the current frame
 *
 *
 * OUTPUT:NUM_CEP_COEFF + 1 coefficients (including c0) + LogEnergy
 *
 *
 * RETURN VALUE:TRUE  (a frame is always output)
 *
 * WARNING !!!  :   The preceeding sample (index=-1) must be present in Data
 *---------------------------------------------------------------------------*/
extern BOOLEAN 
DoCompCeps(X_FLOAT32 * Data, X_FLOAT32 * Coef, FEParamsX * This)
{
  CompCepsStructX *CCX = This->CCX;
  DataFor16kProc *pData16k = This->pData16k;

  WI8CompCeps(CCX, Data, Coef, pData16k);

  return TRUE;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoCompCepsDelete
 *
 * PURPOSE:
 *
 *
 * INPUT:
 *
 *
 * OUTPUT:
 *
 *
 * RETURN VALUE:
 *
 *
 *---------------------------------------------------------------------------*/
extern void 
DoCompCepsDelete(CompCepsStructX * CCX)
{
  if (CCX->pDCTMatrix != NULL) {
    free(CCX->pDCTMatrix);
  }
  if (CCX != NULL) {
    ReleaseMelFBwindows(CCX->FirstWindow);
    if (CCX->FirstWindow != NULL) {
      free(CCX->FirstWindow);
    }
    free(CCX);
  }
  return;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoCompCepsDelete
 *
 * PURPOSE:
 *
 *
 * INPUT:
 *
 *
 * OUTPUT:
 *
 *
 * RETURN VALUE:
 *
 *
 *---------------------------------------------------------------------------*/
void 
WI8CompCeps(CompCepsStructX * CCX, X_FLOAT32 * curBuffer, X_FLOAT32 * cepOut, DataFor16kProc * pData16k)
{
  X_FLOAT32 LogEnergy;
  X_FLOAT32 EnergyFloor_FB;
  X_FLOAT32 EnergyFloor_logE;
  X_FLOAT32 *FloatBuffer;

  X_INT16 i;

  /*-------------------
   * 16 kHz processing
   *-------------------*/
  X_INT16 hp16kBandsSize;
  X_FLOAT32 hp16kBands[HP16k_MEL_USED];
  X_FLOAT32 *codeWeights;
  X_FLOAT32 *CodeForBands16k;
  X_FLOAT32 percCoded;

  FloatBuffer = (X_FLOAT32 *) malloc(sizeof(X_FLOAT32) * CC_FFT_LENGTH);
  if (FloatBuffer == NULL) {
    fprintf(stderr, "ERROR:   Memory allocation error occured!\r\n");
    exit(0);
  }
  if (CCX->Do16kHzProc) {
    hp16kBandsSize = HP16k_MEL_USED;
    codeWeights = Get16k_p_bufferCodeWeights(CCX->pData16k);
    CodeForBands16k = Get16k_p_bufferCodeForBands16k(CCX->pData16k);
    percCoded = Get16k_percCoded(CCX->pData16k);
  } else {
    hp16kBandsSize = 0;
    codeWeights = NULL;
    CodeForBands16k = NULL;
    percCoded = 0.0;
  }

  /*----------------
   * Initialization
   *----------------*/
  EnergyFloor_FB = (X_FLOAT32) exp((double) CC_ENERGYFLOOR_FB);
  EnergyFloor_logE = (X_FLOAT32) exp((double) CC_ENERGYFLOOR_logE);

  /*------------------
   * logE computation
   *------------------*/
  LogEnergy = 0.0;
  for (i = 0; i < CC_FRAME_LENGTH; i++) {
    LogEnergy += curBuffer[i] * curBuffer[i];
  }

  if (!CCX->Do16kHzProc) {
    if (LogEnergy < EnergyFloor_logE)
      LogEnergy = CC_ENERGYFLOOR_logE;
    else
      LogEnergy = (X_FLOAT32) log((double) LogEnergy);
  }
  /*-------------------------------------------------
   * Pre-emphasis, notice that curBuffer[-1] is used
   *-------------------------------------------------*/
  for (i = 0; i < CC_FRAME_LENGTH; i++) {
    FloatBuffer[i] = curBuffer[i] - CC_PRE_EMPHASIS * curBuffer[i - 1];
  }

  /*-----------
   * Windowing
   *-----------*/
  Window(FloatBuffer, CCX->HammingWindow, CC_FRAME_LENGTH);

  /*--------------
   * Zero padding
   *--------------*/
  for (i = CC_FRAME_LENGTH; i < CC_FFT_LENGTH; i++)
    FloatBuffer[i] = 0.0;

  /*---------------------------------------
   * Real valued, in-place split-radix FFT
   *---------------------------------------*/
  rfft(FloatBuffer, CC_FFT_LENGTH, CC_FFT_LENGTH_ORDER);

  /*----------------
   * Power spectrum
   *----------------*/
  //DC
    FloatBuffer[0] = (double) FloatBuffer[0] * FloatBuffer[0];

  //pi / (N / 2), 2 pi / (N / 2),..., (N / 2 - 1) * pi / (N / 2)
    for (i = 1; i < CC_FFT_LENGTH / 2; i++)
    FloatBuffer[i] = (X_FLOAT32) ((double) FloatBuffer[i] * (double) FloatBuffer[i] +
				  (double) FloatBuffer[CC_FFT_LENGTH - i] * (double) FloatBuffer[CC_FFT_LENGTH - i]);

  //pi / 2
    FloatBuffer[CC_FFT_LENGTH / 2] = (double) FloatBuffer[CC_FFT_LENGTH / 2] * FloatBuffer[CC_FFT_LENGTH / 2];

  /*-------------------
   * 16 kHz processing
   *-------------------*/
  if (CCX->Do16kHzProc) {
    X_FLOAT32 aux[3];

    GetBandsForDecoding16k(FloatBuffer, aux, CC_FFT_LENGTH / 2 + 1);
    for (i = 0; i < 3; i++) {
      if (aux[i] > EnergyFloor_FB)
	aux[i] = (X_FLOAT32) log((double) aux[i]);
      else
	aux[i] = CC_ENERGYFLOOR_FB;
    }
    //temporarily
      codeWeights[0] = 0.1;
    codeWeights[1] = 0.2;
    codeWeights[2] = 0.7;
    DecodeBands16k(hp16kBands, aux, CodeForBands16k, codeWeights);
  }
  /*---------------
   * Mel filtering
   *---------------*/
  DoMelFB(FloatBuffer, CCX->FirstWindow);

  /*-------------------
   * 16 kHz processing
   *-------------------*/
  if (CCX->Do16kHzProc) {
    X_FLOAT32 preemphs_correction = log(1.0 + CC_PRE_EMPHASIS);
    X_FLOAT32 *hpBands;

    hpBands = Get16k_p_hpBands(pData16k);

    /*-----------------------------------------------------------
     * SS HP bands -> including and rough preemphasis correction
     *-----------------------------------------------------------*/
    for (i = CC_NUM_CHANNELS_WI8; i < CC_NUM_CHANNELS_WI8 + hp16kBandsSize; i++)
      FloatBuffer[i] = (1.0 + CC_PRE_EMPHASIS) * hpBands[i - CC_NUM_CHANNELS_WI8];

    /*------------------------------------------------
     * coded HP bands -> rough preemphasis correction
     *------------------------------------------------*/
    for (i = 0; i < hp16kBandsSize; i++)
      hp16kBands[i] += preemphs_correction;
  }
  /*-------------------------------
   * Natural logarithm computation
   *-------------------------------*/
  for (i = 0; i < CC_NUM_CHANNELS_WI8 + hp16kBandsSize; i++)
    if (FloatBuffer[i] < EnergyFloor_FB)
      FloatBuffer[i] = CC_ENERGYFLOOR_FB;
    else
      FloatBuffer[i] = (X_FLOAT32) log((double) FloatBuffer[i]);

  /*-------------------
   * 16 kHz processing
   *-------------------*/
  if (CCX->Do16kHzProc) {
    /*--------------------------------------------
     * mixing coded and SS HP bands in log domain
     *--------------------------------------------*/
    MergeSSandCoded(FloatBuffer, hp16kBands, CC_NUM_CHANNELS_WI8, &hp16kBandsSize, CCX->pData16k);

    LogEnergy = CorrectEnergy(LogEnergy, hp16kBandsSize, FloatBuffer, CC_NUM_CHANNELS_WI8, CC_PRE_EMPHASIS);

    if (LogEnergy < EnergyFloor_logE)
      LogEnergy = CC_ENERGYFLOOR_logE;
    else
      LogEnergy = (X_FLOAT32) log((double) LogEnergy);
  }
  /*---------------------------
   * Discrete Cosine Transform
   *---------------------------*/
  DCT(FloatBuffer, CCX->pDCTMatrix, NUM_CEP_COEFF, CC_NUM_CHANNELS_WI8 + hp16kBandsSize);

  /*--------------------------------------
   * Append logE after c0 or overwrite c0
   *--------------------------------------*/
  FloatBuffer[CC_NUM_CHANNELS_WI8 + hp16kBandsSize + NUM_CEP_COEFF - (CCX->Noc0 ? 1 : 0)] = LogEnergy;

  for (i = 0; i < NUM_CEP_COEFF + 1; i++) {
    cepOut[i] = FloatBuffer[CC_NUM_CHANNELS_WI8 + hp16kBandsSize + i];
  }
  if (FloatBuffer != NULL) {
    free(FloatBuffer);
  }

  return;
}
