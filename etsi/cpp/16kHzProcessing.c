/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: 16kHzProcessing.c
 * PURPOSE:   For processing of 16 kHz sampled input signal.
 *
 *-------------------------------------------------------------------------------*/
/*-----------------
 * File Inclusions
 *-----------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ParmInterface.h"
#include "NoiseSupExports.h"
#include "16kHzProcExports.h"
#include "MelProcExports.h"

/*-------------
 * definitions
 *-------------*/
#define NS_SPEC_ORDER_16K                (X_INT16)(65)
#define NS_HANGOVER_16K                  (X_INT16)(15)
#define NS_MIN_SPEECH_FRAME_HANGOVER_16K (X_INT16)(4)

#define NS_EPS_16K                   (X_FLOAT32)(exp ((double) -10.0))

/*------------
 * structures
 *------------*/
struct QMF_FIR {
  int lengthQMF;
    //number of FIR coefficients
  float *dp;
    //pointer to array with FIR coeff.
  float *hp;
    //pointer to array with FIR coeff.
  float *T;
    //pointer to delay line
};

struct DataFor16kProc {
  //info needed from FEParamsX
  // computation of numFramesInBuffer assumes the same
  // FrameLength and FrameShift for both NSX and CCX
  int FrameLength;
  int FrameShift;
  int numFramesInBuffer;
  int SamplingFrequency;

  BOOLEAN Do16kHzProc;
  float *hpBands;
  int hpBandsSize;
  float CodeForBands16k[HP16k_MEL_USED * NB_LP_BANDS_CODING];
  float bufferCodeForBands16k[3 * HP16k_MEL_USED * NB_LP_BANDS_CODING];	/* 3 = delay between the 1st stage output and MFCC proc */
  float codeWeights[NB_LP_BANDS_CODING];
  float bufferCodeWeights[3 * NB_LP_BANDS_CODING];	/* 3 = delay between the 1st stage output and MFCC proc */

  QMF_FIR *pQMF_Fir;

  float *bufferData16k;
  int bufData16kSize;
  MelFB_Window *FirstWindow16k;
  float noiseSE16k[HP16k_MEL_USED];
  float *HammWin;
  float BandsForCoding16k[9];
  long int vadCounter16k;
  int vad16k;
  int nbSpeechFrames16k;
  int hangOver16k;
  float meanEn16k;

  float percCoded;
  int nb_frame_threshold_nse;
  float lambda_nse;
  float *dataHP;
};


/*---------------------------------------------------------------------------
 * FUNCTION NAME: MergeSSandCoded
 *
 * PURPOSE: Merge spectral subtraction (SS) and coded bands
 *
 * INPUT:
 *   FloatBuffer         Pointer to filter bank bands
 *   hp16kBands          Pointer to hp filter bank bands
 *   numChannelsWI7      Number of filter bank channels
 *   hp16kBandsSize_in   Number of hp filter bank bands
 *   pData16k            Pointer to 16k data structure
 *
 * OUTPUT
 *   Merged SS bands and coded bands at FloatBuffer [numChannelsWI7 ... numChannelsWI7+hp16kBandsSize-1]
 *
 *
 * RETURN VALUE
 *   none
 *
 *---------------------------------------------------------------------------*/
void 
MergeSSandCoded(float *FloatBuffer, float *hp16kBands, int numChannelsWI7, X_INT16 * hp16kBandsSize_in, DataFor16kProc * pData16k)
{
  int i;
  float bandAverage;
  float percCoded = pData16k->percCoded;
  int hp16kBandsSize = *hp16kBandsSize_in;

  for (i = 0; i < hp16kBandsSize; i++)
    FloatBuffer[numChannelsWI7 + i] = percCoded * hp16kBands[i] + (1.0 - percCoded) * FloatBuffer[numChannelsWI7 + i];

  if (hp16kBandsSize > 0) {
    //perform interpolation
      bandAverage = 0.5 * FloatBuffer[numChannelsWI7 - 1] + 0.5 * FloatBuffer[numChannelsWI7];

    FloatBuffer[numChannelsWI7 - 1] = 0.6 * FloatBuffer[numChannelsWI7 - 1] + 0.4 * bandAverage;
    FloatBuffer[numChannelsWI7] = 0.6 * FloatBuffer[numChannelsWI7] + 0.4 * bandAverage;
  }
  *hp16kBandsSize_in = hp16kBandsSize;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: CorrectEnergy
 *
 * PURPOSE: Add hp energy and reflect preemphasis
 *
 * INPUT:
 *  LogEnergy        Input log energy
 *  hp16kBandsSize   Number of hp filter bank bands
 *  FloatBuffer      Pointer to filter bank bands
 *  numChannelsWI7   Number of filter bank channels
 *  preem            Preemphasis coef
 *
 * OUTPUT
 *
 * RETURN VALUE
 *  Returns updated log energy
 *
 *---------------------------------------------------------------------------*/
float 
CorrectEnergy(float LogEnergy, int hp16kBandsSize, float *FloatBuffer, int numChannelsWI7, float preem)
{
  int i;
  float energyHP = 0.0;
  float preemphs_correction = log(1.0 + preem);

  for (i = 0; i < hp16kBandsSize; i++)
    energyHP += exp(FloatBuffer[numChannelsWI7 + i] - preemphs_correction);

  LogEnergy += energyHP;

  return LogEnergy;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: Do16kProcAlloc
 *
 * PURPOSE:
 *
 * INPUT:
 *
 *
 *
 * OUTPUT
 *
 *
 *
 * RETURN VALUE
 *
 *
 *---------------------------------------------------------------------------*/
DataFor16kProc *
Do16kProcAlloc()
{
  DataFor16kProc *newData16k = NULL;
  newData16k = (DataFor16kProc *) calloc(1, sizeof(DataFor16kProc));

  newData16k->FirstWindow16k = CMelFBAlloc();
  if (newData16k == NULL || newData16k->FirstWindow16k == NULL)
    return NULL;
  else
    return newData16k;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: Do16kProcInit
 *
 * PURPOSE: Initialize 16k data structure
 *
 * INPUT:
 *  FEParamsX      Pointer to front-end parameter struct
 *
 *
 * OUTPUT
 *
 *
 *
 * RETURN VALUE
 *
 *
 *---------------------------------------------------------------------------*/
void 
Do16kProcInit(FEParamsX * This)
{
  DataFor16kProc *newData16k = This->pData16k;
  int i;

  newData16k->Do16kHzProc = This->Do16kHzProc;
  newData16k->FrameLength = This->FrameLength;
  newData16k->FrameShift = This->FrameShift;
  newData16k->SamplingFrequency = This->SamplingFrequency;

  newData16k->hpBandsSize = ((int) (newData16k->FrameLength / newData16k->FrameShift) + 1) * HP16k_MEL_USED;
  newData16k->hpBands = (float *) malloc(sizeof(float) * newData16k->hpBandsSize);
  for (i = 0; i < newData16k->hpBandsSize; i++)
    newData16k->hpBands[i] = 0.0;

  //HP - VAD related variables
    newData16k->vadCounter16k = 0;
  newData16k->vad16k = 0;
  newData16k->nbSpeechFrames16k = 0;
  newData16k->hangOver16k = 0;
  newData16k->meanEn16k = 0.0;

  if ((newData16k->FrameLength % newData16k->FrameShift) == 0)
    newData16k->numFramesInBuffer = newData16k->FrameLength / newData16k->FrameShift;
  else
    newData16k->numFramesInBuffer = (newData16k->FrameLength / newData16k->FrameShift) + 1;
  //assuming the same delay for both WF stages
  newData16k->bufData16kSize = 2 * newData16k->numFramesInBuffer * newData16k->FrameShift;
  newData16k->bufferData16k = (float *) malloc(sizeof(float) * newData16k->bufData16kSize);
  for (i = 0; i < newData16k->bufData16kSize; i++)
    newData16k->bufferData16k[i] = 0.0;
  for (i = 0; i < HP16k_MEL_USED; i++)
    newData16k->noiseSE16k[i] = 0.0;

  InitMelFBwindows(newData16k->FirstWindow16k,
		    /* StartingFrequency */ 80.0,
		   (float) newData16k->SamplingFrequency,
		   2 * (NS_SPEC_ORDER_16K - 1),
		   HP16k_MEL_ORDER,
		   0);
  newData16k->HammWin = (float *) malloc(sizeof(float) * newData16k->FrameLength);
  for (i = 0; i < newData16k->FrameLength; i++)
    newData16k->HammWin[i] = 0.54 - 0.46 * cos(PIx2 * (float) i / (float) (newData16k->FrameLength - 1));

  newData16k->pQMF_Fir = QMF_FIR_Init();

  newData16k->percCoded = PERC_CODED;
  newData16k->nb_frame_threshold_nse = NE16k_FRAMES_THRESH;
  newData16k->lambda_nse = LAMBDA_NSE16k;
  newData16k->dataHP = (float *) calloc(sizeof(float), newData16k->FrameShift);
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: Do16kProcDelete
 *
 * PURPOSE:
 *
 * INPUT:
 *
 *
 *
 * OUTPUT
 *
 *
 *
 * RETURN VALUE
 *
 *
 *---------------------------------------------------------------------------*/
void 
Do16kProcDelete(DataFor16kProc * pData16k)
{
  if (pData16k->FirstWindow16k != NULL) {
    ReleaseMelFBwindows(pData16k->FirstWindow16k);
    free(pData16k->FirstWindow16k);
  }
  if (pData16k->bufferData16k != NULL)
    free(pData16k->bufferData16k);

  if (pData16k->hpBands != NULL)
    free(pData16k->hpBands);

  if (pData16k->pQMF_Fir != NULL)
    hq_free(pData16k->pQMF_Fir);

  if (pData16k->HammWin != NULL)
    free(pData16k->HammWin);

  if (pData16k != NULL)
    free(pData16k);
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: GetBandsForDecoding16k
 *
 * PURPOSE: Calculate aux bands for 16k decoding
 *
 * INPUT:
 *  inData                  Pointer to input spectrum data
 *  bandsForCoding16k       Pointer to output aux bands
 *  inDataSize              Size of input spectrum data
 *
 * OUTPUT
 *  Aux bands are stored in bandsForCoding16k
 *
 * RETURN VALUE
 *  none
 *
 *---------------------------------------------------------------------------*/
void 
GetBandsForDecoding16k(float *inData, float *bandsForCoding16k, int inDataSize)
{
  int i;

  for (i = 0; i < 3; i++)
    bandsForCoding16k[i] = 0.0;

  for (i = 66; i < 77; i++)
    bandsForCoding16k[0] += inData[i];

  for (i = 78; i < 97; i++)
    bandsForCoding16k[1] += inData[i];

  for (i = 98; i < 129; i++)
    bandsForCoding16k[2] += inData[i];

  for (i = 0; i < 3; i++)
    bandsForCoding16k[i] /= 2.0;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: DecodeBands16k
 *
 * PURPOSE: Calculate hp bands by using coding/decoding scheme
 *
 * INPUT:
 *  fb16k            Pointer to output decoded bands
 *  lpBands          Aux bands for decoding
 *  codeForBands16k  Pointer to code
 *  codeWeights      Pointer to code weights
 *
 * OUTPUT
 *  Decoded bands are stored in fb16k
 *
 * RETURN VALUE
 *  none
 *
 *---------------------------------------------------------------------------*/
void 
DecodeBands16k(float *fb16k, float *lpBands, float *codeForBands16k, float *codeWeights)
{
  int i, j;

  for (i = 0; i < HP16k_MEL_USED; i++) {
    fb16k[i] = 0.0;
    for (j = 0; j < 3; j++)
      fb16k[i] += (codeWeights[j] * (lpBands[j] - codeForBands16k[3 * i + j]));
  }
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: DoSpecSub16k
 *
 * PURPOSE: Perform spectral subtraction on hp filter bank bands
 *
 * INPUT:
 *  inFB16k         Pointer to 16k filter bank bands
 *  pData16k        Pointer to 16k data struct
 *  vadCounter16k   VAD frame counter
 *  inFB16kSize     Number of 16k filter bank bands
 *
 * OUTPUT
 *  De-noised 16k filter bank bands are stored in inFB16k
 *
 *
 * RETURN VALUE
 *  none
 *
 *---------------------------------------------------------------------------*/
void 
DoSpecSub16k(float *inFB16k, DataFor16kProc * pData16k, long vadCounter16k, int inFB16kSize)
{

  float *noiseFB16k = pData16k->noiseSE16k;
  int *nbSpeechFrames16k = &(pData16k->nbSpeechFrames16k);
  int *hangOver16k = &(pData16k->hangOver16k);
  float *meanEn16k = &(pData16k->meanEn16k);
  int nb_frame_threshold_nse = pData16k->nb_frame_threshold_nse;
  float lambda_nse = pData16k->lambda_nse;

  int i;
  int flagVAD;
  int nbFrame;
  float lambdaNSE;
  float diff, floor;
  float frameEn, meanEn;

  pData16k->vadCounter16k = vadCounter16k;

  nbFrame = vadCounter16k;
  meanEn = *meanEn16k;

  if (nbFrame < nb_frame_threshold_nse)
    //nb_frame_threshold_nse = 100
      lambdaNSE = 1.0 - 1.0 / (X_FLOAT32) nbFrame;
  else
    lambdaNSE = lambda_nse;
  //lambda_nse = 0.99

    // frame energy calc.
    frameEn = 0.0;
  for (i = 0; i < HP16k_MEL_USED; i++)
    frameEn += inFB16k[i];
  if (frameEn > 0.001)
    frameEn = (float) log((double) frameEn);
  else
    frameEn = (float) log((double) 0.001);

  //mean energy calc.
    if (((frameEn - meanEn) < 1.2) || (nbFrame < 10)) {
    if (nbFrame < 10) {
      meanEn += (1 - lambdaNSE) * (frameEn - meanEn);
    } else {
      if (frameEn < meanEn)
	meanEn += (1 - 0.98) * (frameEn - meanEn);
      else
	meanEn += (1 - 0.995) * (frameEn - meanEn);
    }
  }
  if ((frameEn - meanEn) > 2.2) {
    flagVAD = 1;
    (*nbSpeechFrames16k)++;
  } else {
    if (*nbSpeechFrames16k > NS_MIN_SPEECH_FRAME_HANGOVER_16K)
      *hangOver16k = NS_HANGOVER_16K;
    *nbSpeechFrames16k = 0;
    if (*hangOver16k != 0) {
      (*hangOver16k)--;
      flagVAD = 1;
    } else
      flagVAD = 0;
  }

  if (flagVAD == 0) {
    for (i = 0; i < inFB16kSize; i++) {
      if (nbFrame < 10) {
	noiseFB16k[i] = lambdaNSE * noiseFB16k[i] + (1 - lambdaNSE) * inFB16k[i];
      } else {
	if (inFB16k[i] < noiseFB16k[i])
	  noiseFB16k[i] = lambdaNSE * noiseFB16k[i] + (1 - lambdaNSE) * inFB16k[i];
	else
	  noiseFB16k[i] = 0.995 * noiseFB16k[i] + (1 - 0.995) * inFB16k[i];
      }
      if (noiseFB16k[i] < NS_EPS_16K)
	noiseFB16k[i] = NS_EPS_16K;
    }
  }
  for (i = 0; i < inFB16kSize; i++) {
    floor = 0.1 * inFB16k[i];
    diff = inFB16k[i] - 1.5 * noiseFB16k[i];
    if (diff > floor)
      inFB16k[i] = diff;
    else
      inFB16k[i] = floor;
  }

  *meanEn16k = meanEn;

  return;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: GetBandsForCoding16k
 *
 * PURPOSE: Calculate aux bands for 16k coding
 *
 * INPUT:
 *  inData                Pointer to input spectrum data
 *  bandsForCoding16k     Pointer to output aux bands
 *  inDataSize            Size of input spectrum data
 *
 * OUTPUT
 *  Aux bands are stored in bandsForCoding16k
 *
 * RETURN VALUE
 *  none
 *
 *---------------------------------------------------------------------------*/
void 
GetBandsForCoding16k(float *inData, float *bandsForCoding16k, int inDataSize)
{
  int i;

  for (i = 0; i < 3; i++)
    bandsForCoding16k[i] = 0.0;

  for (i = 33; i < 39; i++)
    bandsForCoding16k[0] += inData[i];
  for (i = 39; i < 49; i++)
    bandsForCoding16k[1] += inData[i];
  for (i = 49; i < 65; i++)
    bandsForCoding16k[2] += inData[i];
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: CodeBands16k
 *
 * PURPOSE: Code 16k filter bank bands
 *
 * INPUT:
 *  fb16k             Pointer to 16k filter bank bands
 *  lpBands           Pointer to aux bands for coding
 *  codeForBands16k   Pointer to output code
 *
 * OUTPUT
 *  Output code is stored in codeForBands16k
 *
 * RETURN VALUE
 *  none
 *
 *---------------------------------------------------------------------------*/
void 
CodeBands16k(float *fb16k, float *lpBands, float *codeForBands16k)
{
  int i, j;

  for (i = 0; i < HP16k_MEL_USED; i++)
    for (j = 0; j < 3; j++)
      codeForBands16k[3 * i + j] = lpBands[j] - fb16k[i];
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: DP_HP_filters
 *
 * PURPOSE: Definition of low-pass and high-pass filters
 *
 * INPUT:
 *
 *
 *
 * OUTPUT
 *
 *
 *
 * RETURN VALUE
 *
 *
 *---------------------------------------------------------------------------*/
#define f24 (float)0x00800000
#define LENGTH_QMF 118
void 
DP_HP_filters(float *dp[], float *hp[], long *lengthQMF)
{
  int i, aux;
  static float hp02[LENGTH_QMF];
  static float dp02[LENGTH_QMF] = {
    1584. / f24, 805. / f24, -4192. / f24, -8985. / f24,
    -5987. / f24, 2583. / f24, 4657. / f24, -3035. / f24,
    -7004. / f24, 1542. / f24, 8969. / f24, 567. / f24,
    -10924. / f24, -3757. / f24, 12320. / f24, 7951. / f24,
    -12793. / f24, -13048. / f24, 11923. / f24, 18793. / f24,
    -9331. / f24, -24802. / f24, 4694. / f24, 30570. / f24,
    2233. / f24, -35439. / f24, -11526. / f24, 38680. / f24,
    23114. / f24, -39474. / f24, -36701. / f24, 36999. / f24,
    51797. / f24, -30419. / f24, -67658. / f24, 18962. / f24,
    83318. / f24, -1927. / f24, -97566. / f24, -21284. / f24,
    108971. / f24, 51215. / f24, -115837. / f24, -88430. / f24,
    116130. / f24, 133716. / f24, -107253. / f24, -188497. / f24,
    85497. / f24, 255795. / f24, -44643. / f24, -342699. / f24,
    -28185. / f24, 468096. / f24, 167799. / f24, -696809. / f24,
    -519818. / f24, 1446093. / f24, 3562497. / f24, 3562497. / f24,
    1446093. / f24, -519818. / f24, -696809. / f24, 167799. / f24,
    468096. / f24, -28185. / f24, -342699. / f24, -44643. / f24,
    255795. / f24, 85497. / f24, -188497. / f24, -107253. / f24,
    133716. / f24, 116130. / f24, -88430. / f24, -115837. / f24,
    51215. / f24, 108971. / f24, -21284. / f24, -97566. / f24,
    -1927. / f24, 83318. / f24, 18962. / f24, -67658. / f24,
    -30419. / f24, 51797. / f24, 36999. / f24, -36701. / f24,
    -39474. / f24, 23114. / f24, 38680. / f24, -11526. / f24,
    -35439. / f24, 2233. / f24, 30570. / f24, 4694. / f24,
    -24802. / f24, -9331. / f24, 18793. / f24, 11923. / f24,
    -13048. / f24, -12793. / f24, 7951. / f24, 12320. / f24,
    -3757. / f24, -10924. / f24, 567. / f24, 8969. / f24,
    1542. / f24, -7004. / f24, -3035. / f24, 4657. / f24,
    2583. / f24, -5987. / f24, -8985. / f24, -4192. / f24,
  805. / f24, 1584. / f24};

  aux = 1;
  for (i = 0; i < LENGTH_QMF; i++) {
    hp02[LENGTH_QMF - 1 - i] = (float) aux *dp02[i];
    aux *= (-1);
  }
  *lengthQMF = LENGTH_QMF;
  //store 'number of coefficients'
    * dp = dp02;
  //store pointer to dp02[] - array
    * hp = hp02;
}
#undef f24
#undef LENGTH_QMF

/*---------------------------------------------------------------------------
 * FUNCTION NAME: fir_initialization
 *
 * PURPOSE:
 *
 * INPUT:
 *
 *
 *
 * OUTPUT
 *
 *
 *
 * RETURN VALUE
 *
 *
 *---------------------------------------------------------------------------*/
QMF_FIR *
fir_initialization(int lengthQMF, float dp[], float hp[])
{
  QMF_FIR *ptrFIR; //pointer to the new struct
  int k;

  //.........ALLOCATION OF MEMORY.........
  // Allocate memory for a new struct
  ptrFIR = (QMF_FIR *) malloc(sizeof(QMF_FIR));
  
  //Allocate memory for delay line
  ptrFIR->T = (float *) malloc((lengthQMF - 1) * sizeof(float));
  
  //Allocate memory for impulse response
  ptrFIR->dp = (float *) malloc(lengthQMF * sizeof(float));
  ptrFIR->hp = (float *) malloc(lengthQMF * sizeof(float));

  //.........STORE VARIABLES INTO STATE VARIABLE.........
  // Store number of FIR - coefficients
  ptrFIR->lengthQMF = lengthQMF;
  
  //Fill FIR coefficients into struct
  for (k = 0; k <= ptrFIR->lengthQMF - 1; k++)
    ptrFIR->dp[k] = dp[k];

  //Fill FIR coefficients into struct
  for (k = 0; k <= ptrFIR->lengthQMF - 1; k++)
    ptrFIR->hp[k] = hp[k];

  //Clear Delay Line
  for (k = 0; k < ptrFIR->lengthQMF - 1; k++)
    ptrFIR->T[k] = 0.0;

  return (ptrFIR);
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: QMF_FIR_Init
 *
 * PURPOSE:
 *
 * INPUT:
 *
 *
 *
 * OUTPUT
 *
 *
 *
 * RETURN VALUE
 *
 *
 *---------------------------------------------------------------------------*/
QMF_FIR *
QMF_FIR_Init()
{
  float *dp, *hp; //pointer to array with FIR coeff.
  long lengthQMF; //number of FIR coefficients

  // allocate array for FIR coeff.and fill with coefficients
  DP_HP_filters(&dp, &hp, &lengthQMF);

  //Returns:pointer to QMF_FIR - struct
  return fir_initialization(lengthQMF, dp, hp);
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: Do16kProcessing
 *
 * PURPOSE: Calculate two 0-4kHz and 4-8kHz waveforms from the input 0-8kHz waveform
 *
 * INPUT:
 *  inData         Pointer to 16k waveform
 *  pData16k       Pointer to 16k data struct
 *  inDataLength   Input waveform length
 *
 * OUTPUT
 *  Output 0-4kHz waveform is stored in inData[0 ... inDataLength/2-1]
 *  Output 4-8kHz waveform is stored in pData16k->dataHP[0 ... inDataLength/2-1]
 *
 * RETURN VALUE
 *  none
 *
 *---------------------------------------------------------------------------*/
void 
Do16kProcessing(float *inData, DataFor16kProc * pData16k, int inDataLength)
{
  QMF_FIR *pQMF_Fir = pData16k->pQMF_Fir;

  int i, j, k;
  float *outDataDP, *outDataHP;
  float aux1, aux2;
  int hpSign;

  float *prevData;
  float *filterDP;
  float *filterHP;
  int filterLength;

  prevData = pQMF_Fir->T;
  filterDP = pQMF_Fir->dp;
  filterHP = pQMF_Fir->hp;
  filterLength = pQMF_Fir->lengthQMF;

  outDataDP = (float *) malloc(sizeof(float) * (inDataLength + filterLength - 1));
  outDataHP = (float *) malloc(sizeof(float) * (inDataLength + filterLength - 1));

  for (i = 0; i < filterLength - 1; i++)
    outDataDP[i] = outDataHP[i] = prevData[i];

  for (i = 0; i < inDataLength; i++)
    outDataDP[filterLength - 1 + i] = outDataHP[filterLength - 1 + i] = inData[i];

  //FIR filtering
    hpSign = 1;
  for (i = 0, k = 0; i < inDataLength; i += 2, k++) {
    //computing only every second sample
      aux1 = aux2 = 0.0;
    for (j = 0; j < filterLength; j++) {
      aux1 += outDataDP[i + j] * filterDP[j];
      //LP filtering
	aux2 += outDataHP[i + j] * filterHP[j];
      //HP filtering
    }
    outDataDP[k] = aux1;
    outDataHP[k] = (float) hpSign *aux2;
    //shifting from 4 - 8 kHz to 0 - 4 kHz
      hpSign *= (-1);
  }

  if (inDataLength < filterLength) {
    for (i = 0; i < (filterLength - inDataLength); i++)
      prevData[i] = prevData[inDataLength + i];
    for (i = 0; i < inDataLength - 1; i++)
      prevData[i + (filterLength - inDataLength)] = inData[i];
  } else
    for (i = 0; i < filterLength - 1; i++)
      prevData[i] = inData[i + inDataLength - filterLength + 1];

  for (i = 0; i < inDataLength / 2; i++)
    //LP data to lower half of inData
      inData[i] = outDataDP[i];
  for (i = 0; i < inDataLength / 2; i++)
    //LP data to lower half of inData
      pData16k->dataHP[i] = outDataHP[i];

  free(outDataDP);
  free(outDataHP);
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: hq_free
 *
 * PURPOSE:
 *
 * INPUT:
 *
 *
 *
 * OUTPUT
 *
 *
 *
 * RETURN VALUE
 *
 *
 *---------------------------------------------------------------------------*/
void 
hq_free(QMF_FIR * fir_ptr)
{
  if (fir_ptr->T != NULL)
    free(fir_ptr->T);
  //free state variables
    if (fir_ptr->dp)
    free(fir_ptr->dp);
  //free state impulse response
    if (fir_ptr->hp)
    free(fir_ptr->hp);
  //free state impulse response
    if (fir_ptr)
    free(fir_ptr);
  //free allocated struct
}

/*---------------
 * get functions
 *---------------*/
/*---------------------------------------------------------------------------
 * FUNCTION NAME: Get16k_hpBandsSize
 *
 * PURPOSE: Return number of 16k filter bank bands
 *
 * INPUT:
 *  pData16k       Pointer to 16k data struct
 *
 *
 * OUTPUT
 *
 * RETURN VALUE
 *  Number of 16k filter bank bands (int)
 *
 *---------------------------------------------------------------------------*/
int 
Get16k_hpBandsSize(DataFor16kProc * pData16k)
{
  return pData16k->hpBandsSize;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: Get16k_p_hpBands
 *
 * PURPOSE: Return pointer to 16k filter bank bands
 *
 * INPUT:
 *  pData16k     Pointer to 16k data struct
 *
 *
 * OUTPUT
 *
 *
 *
 * RETURN VALUE
 *  Pointer to 16k filter bank bands (float *)
 *
 *---------------------------------------------------------------------------*/
float *
Get16k_p_hpBands(DataFor16kProc * pData16k)
{
  return pData16k->hpBands;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: Get16k_p_bufferCodeForBands16k
 *
 * PURPOSE: Return pointer to 16k code buffer
 *
 * INPUT:
 *  pData16k     Pointer to 16k data struct
 *
 *
 * OUTPUT
 *
 *
 *
 * RETURN VALUE
 *  Pointer to 16k code buffer (float *)
 *
 *---------------------------------------------------------------------------*/
float *
Get16k_p_bufferCodeForBands16k(DataFor16kProc * pData16k)
{
  return pData16k->bufferCodeForBands16k;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: Get16k_p_CodeForBands16k
 *
 * PURPOSE: Return pointer to 16k code
 *
 * INPUT:
 *  pData16k     Pointer to 16k data struct
 *
 *
 * OUTPUT
 *
 *
 *
 * RETURN VALUE
 *  Pointer to 16k code (float *)
 *
 *---------------------------------------------------------------------------*/

float *
Get16k_p_CodeForBands16k(DataFor16kProc * pData16k)
{
  return pData16k->CodeForBands16k;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: Get16k_p_bufferCodeWeights
 *
 * PURPOSE: Return pointer to 16k buffer code weights
 *
 * INPUT:
 *  pData16k     Pointer to 16k data struct
 *
 *
 * OUTPUT
 *
 *
 *
 * RETURN VALUE
 *  Pointer to 16k buffer code weights (float *)
 *
 *---------------------------------------------------------------------------*/
float *
Get16k_p_bufferCodeWeights(DataFor16kProc * pData16k)
{
  return pData16k->bufferCodeWeights;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: Get16k_p_codeWeights
 *
 * PURPOSE: Return pointer to 16k code weights
 *
 * INPUT:
 *  pData16k     Pointer to 16k data struct
 *
 *
 * OUTPUT
 *
 *
 *
 * RETURN VALUE
 *  Pointer to 16k code weights (float *)
 *
 *---------------------------------------------------------------------------*/
float *
Get16k_p_codeWeights(DataFor16kProc * pData16k)
{
  return pData16k->codeWeights;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: Get16k_percCoded
 *
 * PURPOSE: Return the percentage for coded bands
 *
 * INPUT:
 *  pData16k     Pointer to 16k data struct
 *
 *
 * OUTPUT
 *
 *
 *
 * RETURN VALUE
 *  Percentage for coded bands (float)
 *
 *---------------------------------------------------------------------------*/
float 
Get16k_percCoded(DataFor16kProc * pData16k)
{
  return pData16k->percCoded;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: Get16k_p_bufferData16k
 *
 * PURPOSE: Return pointer to 16k data buffer
 *
 * INPUT:
 *  pData16k     Pointer to 16k data struct
 *
 *
 * OUTPUT
 *
 *
 *
 * RETURN VALUE
 *  Pointer to 16k data buffer (float *)
 *
 *---------------------------------------------------------------------------*/
float *
Get16k_p_bufferData16k(DataFor16kProc * pData16k)
{
  return pData16k->bufferData16k;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: Get16k_bufData16kSize
 *
 * PURPOSE: Return size of 16k data buffer
 *
 * INPUT:
 *  pData16k     Pointer to 16k data struct
 *
 *
 * OUTPUT
 *
 *
 *
 * RETURN VALUE
 *  Size of 16k data buffer (int)
 *
 *---------------------------------------------------------------------------*/
int 
Get16k_bufData16kSize(DataFor16kProc * pData16k)
{
  return pData16k->bufData16kSize;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: Get16k_p_BandsForCoding16k
 *
 * PURPOSE: Return pointer to bands for coding
 *
 * INPUT:
 *  pData16k     Pointer to 16k data struct
 *
 *
 * OUTPUT
 *
 *
 *
 * RETURN VALUE
 *  Pointer to bands for coding (float *)
 *
 *---------------------------------------------------------------------------*/
float *
Get16k_p_BandsForCoding16k(DataFor16kProc * pData16k)
{
  return pData16k->BandsForCoding16k;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: Get16k_p_FirstWindow16k
 *
 * PURPOSE: Return pointer to the first 16k filter bank window
 *
 * INPUT:
 *  pData16k     Pointer to 16k data struct
 *
 *
 * OUTPUT
 *
 *
 *
 * RETURN VALUE
 *  Pointer to the first 16k filter bank window (MelFB_Window *)
 *
 *---------------------------------------------------------------------------*/
MelFB_Window *
Get16k_p_FirstWindow16k(DataFor16kProc * pData16k)
{
  return pData16k->FirstWindow16k;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: Get16k_dataHP
 *
 * PURPOSE: Return i-th value of 4-8kHz waveform
 *
 * INPUT:
 *  pData16k     Pointer to 16k data struct
 *  i            Position of returned value
 *
 * OUTPUT
 *
 *
 *
 * RETURN VALUE
 *  i-th value of 4-8kHz waveform (float)
 *
 *---------------------------------------------------------------------------*/
float 
Get16k_dataHP(DataFor16kProc * pData16k, int i)
{
  return pData16k->dataHP[i];
}
