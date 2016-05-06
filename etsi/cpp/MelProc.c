/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: MelProc.c
 * PURPOSE: 1) Mel-filtering of FFT spectrum using the Mel filter bank
 *             frequency windows.
 *          2) Mel inverse DCT applied on the output of Mel filter bank.
 *
 *-------------------------------------------------------------------------------*/
/*-----------------
 * File Inclusions
 *-----------------*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "ParmInterface.h"
#include "MelProcExports.h"

/*------------------------
 * Definitions and Macros
 *------------------------*/
#define max(a,b) ((a>b)?(a):(b))

/*------------------
 * mel FB functions
 *------------------*/

/*----------------------------------------------------------------------------
 * FUNCTION NAME: CMelFBAlloc
 *
 * PURPOSE:       Memory allocation of first Mel window
 *
 * INPUT:
 *   none
 *
 * OUTPUT:
 *   *FirstWin    Pointer to Mel window
 *
 * RETURN VALUE:
 *   *FirstWin    Pointer to Mel window
 *
 *---------------------------------------------------------------------------*/
MelFB_Window *
CMelFBAlloc()
{
  MelFB_Window *FirstWin = NULL;
  FirstWin = (MelFB_Window *)calloc(1, sizeof(MelFB_Window));

  if (FirstWin == NULL)
    return NULL;
  else
    return FirstWin;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: DoMelFB
 *
 * PURPOSE:       Performs mel filtering on FFT magnitude spectrum using the
 *                filter bank defined by a chained list of FFT window
 *                structures
 *
 * INPUT:
 *   SigFFT       Pointer to signal FFT magnitude spectrum
 *   FirstWin     Pointer to the first channel of the filter bank (first
 *                element in the chained list of FFT window data structures)
 *
 * OUTPUT
 *                Filter bank outputs stored at the beginning of input signal
 *                FFT buffer pointed by *SigFFT*
 *
 * RETURN VALUE
 *   none
 *
 *---------------------------------------------------------------------------*/
void 
DoMelFB(float *SigFFT, MelFB_Window * FirstWin)
{
  MelFB_Window *p1;
  float Sum[WF_MEL_ORDER];
  int i, j;

  p1 = FirstWin;
  j = 0;
  while (p1) {
    Sum[j] = 0.0;
    for (i = 0; i < p1->Length; i++)
      Sum[j] += SigFFT[p1->StartingPoint + i] * p1->Data[i];

    j++;
    p1 = p1->Next;
  }

  for (j--; j >= 0; j--)
    SigFFT[j] = Sum[j];

  return;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: InitMelFBwindows
 *
 * PURPOSE:       Initializes data structure for FFT windows (mel filter bank).
 *                Computes starting point and length of each window, allocates
 *                memory for window coefficients.
 *
 * INPUT:
 *   FirstWin     Pointer to first FFT window structure
 *   StFreq       Starting frequency of mel filter bank
 *   SmplFreq     Sampling frequency
 *   FFTLength    FFT length
 *   NumChannels  Number of channels
 *   normalize    Boolean
 *
 * OUTPUT
 *                Chained list of FFT window data structures.
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
void 
InitMelFBwindows(MelFB_Window * FirstWin, float StFreq, float SmplFreq, int FFTLength, int NumChannels, BOOLEAN normalize)
{
  int i, j, k;
  int *centrFreq;
  float freq;
  float start_mel;
  float fs_per_2_mel;
  float normFBconst;
  MelFB_Window *p1, *p2;

  centrFreq = (int *) malloc(sizeof(int) * NumChannels);

  /* Constants for calculation */
  start_mel = 2595.0 * log10(1.0 + (float) StFreq / 700.0);
  fs_per_2_mel = 2595.0 * log10(1.0 + (SmplFreq / 2) / 700.0);

  for (i = 0; i < NumChannels; i++) {
    freq = 700 * (pow(10, (start_mel + (float) i / (NumChannels - 1) * (fs_per_2_mel - start_mel)) / 2595.0) - 1.0);
    centrFreq[i] = (int) (FFTLength * freq / SmplFreq + 0.5);
  }
  p1 = FirstWin;

  /*----------------------
   * first freq. window 0
   *----------------------*/
  p1->StartingPoint = centrFreq[0];
  p1->Length = centrFreq[1] - centrFreq[0];
  p1->Data = (float *) malloc(sizeof(float) * p1->Length);
  if (p1->Data == NULL) {
    fprintf(stderr, "ERROR:   Memory allocation error occured!\r\n");
    exit(0);
  }
  normFBconst = 0.0;
  for (j = 0; j < p1->Length; j++) {
    p1->Data[j] = 1.0 - (float) j / (float) p1->Length;
    normFBconst += p1->Data[j];
  }
  if (normalize)
    for (j = 0; j < p1->Length; j++)
      p1->Data[j] /= normFBconst;
  p2 = (MelFB_Window *) malloc(sizeof(MelFB_Window));
  if (p2 == NULL) {
    fprintf(stderr, "ERROR:   Memory allocation error occured!\r\n");
    exit(0);
  }
  p1->Next = p2;
  p1 = p2;

  /*----------------------------------
   * freq. windows 1->NumChannels - 2
   *----------------------------------*/
  for (i = 1; i < NumChannels - 1; i++) {
    p1->StartingPoint = centrFreq[i - 1] + 1;
    p1->Length = centrFreq[i + 1] - centrFreq[i - 1] - 1;
    p1->Data = (float *) malloc(sizeof(float) * p1->Length);
    if (p1->Data == NULL) {
      fprintf(stderr, "ERROR:   Memory allocation error occured!\r\n");
      exit(0);
    }
    normFBconst = 0.0;
    for (j = 0; j < (centrFreq[i] - centrFreq[i - 1]); j++) {
      p1->Data[j] = (float) (j + 1) / (float) (centrFreq[i] - centrFreq[i - 1]);
      normFBconst += p1->Data[j];
    }
    for (j = (centrFreq[i] - centrFreq[i - 1]), k = 0; j < p1->Length; j++, k++) {
      p1->Data[j] = 1.0 - (k + 1) / (float) (centrFreq[i + 1] - centrFreq[i]);
      normFBconst += p1->Data[j];
    }
    if (normalize)
      for (j = 0; j < p1->Length; j++)
	p1->Data[j] /= normFBconst;

    p2 = (MelFB_Window *) malloc(sizeof(MelFB_Window));
    if (p2 == NULL) {
      fprintf(stderr, "ERROR:   Memory allocation error occured!\r\n");
      exit(0);
    }
    p1->Next = p2;
    p1 = p2;
  }

  /*-----------------------------------
   * last freq. window NumChannels - 1
   *-----------------------------------*/
  p1->StartingPoint = centrFreq[NumChannels - 2] + 1;
  p1->Length = centrFreq[NumChannels - 1] - centrFreq[NumChannels - 2];
  p1->Data = (float *) malloc(sizeof(float) * p1->Length);
  normFBconst = 0.0;
  for (j = 0; j < p1->Length; j++) {
    p1->Data[j] = (float) (j + 1) / (float) p1->Length;
    normFBconst += p1->Data[j];
  }
  if (normalize)
    for (j = 0; j < p1->Length; j++)
      p1->Data[j] /= normFBconst;

  p1->Next = NULL;

  free(centrFreq);

  return;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: ReleaseMelFBwindows
 *
 * PURPOSE:       Releases memory allocated for FFT windows
 *
 * INPUT:
 *   FirstWin     Pointer to first FFT window structure
 *
 * OUTPUT
 *   none
 *
 * RETURN VALUE
 *   none
 *
 *---------------------------------------------------------------------------*/
void 
ReleaseMelFBwindows(MelFB_Window * FirstWin)
{
  MelFB_Window *p;

  while (FirstWin->Next != NULL) {
    p = FirstWin->Next->Next;
    free(FirstWin->Next->Data);
    free(FirstWin->Next);
    FirstWin->Next = p;
  }

  free(FirstWin->Data);
}

/*--------------------
 * mel IDCT functions
 *--------------------*/
/*----------------------------------------------------------------------------
 * FUNCTION NAME: InitMelIDCTbasis
 *
 * PURPOSE:       Initializes Mel IDCT basis
 *
 * INPUT:
 *   *FirstWin    Pointer to the first Mel filter bank window
 *   melorder     Order of Mel filter bank
 *   sampFreq     Sampling frequency
 *   FFTLength    FFT length
 *
 * OUTPUT:
 *   Inverse DCT basis are stored in **melIDCTbasis
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
void 
InitMelIDCTbasis(float **melIDCTbasis, MelFB_Window * FirstWin, short melorder, int sampFreq, int FFTLength)
{
  MelFB_Window *p1;
  float centrFreq[WF_MEL_ORDER];
  float deltaFreq[WF_MEL_ORDER];
  float startingFreq;
  float sum;
  float linStep;
  int i, j;

  linStep = sampFreq / (float) FFTLength;

  /*-------------------------------------
   * calculating centrFreq and deltaFreq
   *-------------------------------------*/
  p1 = FirstWin;
  for (j = 0; j < melorder; j++) {
    if (j == 0) {		/* freq. window 0 */
      centrFreq[j] = p1->StartingPoint * linStep;
      p1 = p1->Next;
    } else {
      if (j == (melorder - 1)) {/* freq. window melorder- */
	centrFreq[j] = (p1->StartingPoint + p1->Length - 1) * linStep;
      } else {			/* freq. windows 1->(melorder-2) */
	startingFreq = p1->StartingPoint * linStep;
	sum = 0.0;
	centrFreq[j] = 0.0;
	for (i = 0; i < p1->Length; i++) {
	  centrFreq[j] += p1->Data[i] * (startingFreq + i * linStep);
	  sum += p1->Data[i];
	}
	centrFreq[j] /= sum;
	p1 = p1->Next;
      }
    }
  }
  for (j = 0; j < melorder; j++) {	/* first and last deltaFreq is only half */
    if (j == 0)
      deltaFreq[j] = (centrFreq[j + 1] - centrFreq[j]) / sampFreq;
    else {
      if (j == (melorder - 1))
	deltaFreq[j] = (centrFreq[j] - centrFreq[j - 1]) / sampFreq;
      else
	deltaFreq[j] = (centrFreq[j + 1] - centrFreq[j - 1]) / sampFreq;
    }
  }

  /*------
   * IDCT
   *------*/
  for (i = 0; i < melorder; i++)/* time axis */
    for (j = 0; j < melorder; j++)	/* frequency axis */
      melIDCTbasis[i][j] = deltaFreq[j] * cos(PIx2 * i * centrFreq[j] / sampFreq);
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoMelIDCT
 *
 * PURPOSE:       Performs Mel IDCT
 *
 * INPUT:
 *   *inData         Pointer to input Mel filter bank bands
 *   **melIDCTbasis  Double pointer to Mel inverse DCT basis
 *   melOrder        Mel filter bank order
 *   timeLength      Length of output data
 *
 * OUTPUT:
 *   Output impulse response is in *inData
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
void 
DoMelIDCT(float *inData, float **melIDCTbasis, int melOrder, int timeLength)
{
  int t, f;
  float *outData;

  outData = (float *) malloc(sizeof(float) * timeLength);
  if (outData == NULL) {
    fprintf(stderr, "ERROR:   Memory allocation error occured!\r\n");
    exit(0);
  }
  for (t = 0; t < timeLength; t++) {
    outData[t] = 0.0;
    for (f = 0; f < melOrder; f++) {
      outData[t] += inData[f] * melIDCTbasis[t][f];
    }
  }
  for (t = 0; t < timeLength; t++)
    inData[t] = outData[t];

  free(outData);
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: InitFFTWindows
 *
 * PURPOSE:       Initializes data structure for FFT windows (mel filter bank).
 *                Computes starting point and length of each window, allocates
 *                memory for window coefficients.
 *
 * INPUT:
 *   FirstWin     Pointer to first FFT window structure
 *   StFreq       Starting frequency of mel filter bank
 *   SmplFreq     Sampling frequency
 *   FFTLength    FFT length
 *   NumChannels  Number of channels
 *
 * OUTPUT
 *                Chained list of FFT window data structures. NOTE FFT window
 *                coefficients are not computed yet.
 *
 * RETURN VALUE
 *   none
 *
 *---------------------------------------------------------------------------*/
void 
InitFFTWindows(MelFB_Window * FirstWin,
	       float StFreq,
	       float SmplFreq,
	       int FFTLength,
	       int NumChannels)
{
  int i, TmpInt;
  float freq, start_mel, fs_per_2_mel;
  MelFB_Window *p1, *p2;

  /*---------------------------
   * Constants for calculation
   *---------------------------*/
  start_mel = 2595.0 * log10(1.0 + (float) StFreq / 700.0);
  fs_per_2_mel = 2595.0 * log10(1.0 + (SmplFreq / 2) / 700.0);

  p1 = FirstWin;

  for (i = 0; i < NumChannels; i++) {
    /*--------------------------------------------------------
     * Calculating mel-scaled frequency and the corresponding
     * FFT-bin number for the lower edge of the band
     *--------------------------------------------------------*/
    freq = 700 * (pow(10, (start_mel + (float) i / (NumChannels + 1) *
			   (fs_per_2_mel - start_mel)) / 2595.0) - 1.0);
    TmpInt = (int) (FFTLength * freq / SmplFreq + 0.5);

    /*---------
     * Storing
     *---------*/
    p1->StartingPoint = TmpInt;

    /*-----------------------------------------------------------------
     * Calculating mel-scaled frequency for the upper edge of the band
     *-----------------------------------------------------------------*/
    freq = 700 * (pow(10, (start_mel + (float) (i + 2) / (NumChannels + 1)
			   * (fs_per_2_mel - start_mel)) / 2595.0) - 1.0);

    /*---------------------------------------------------------------------
     * Calculating and storing the length of the band in terms of FFT-bins
     *---------------------------------------------------------------------*/
    p1->Length = (int) (FFTLength * freq / SmplFreq + 0.5) - TmpInt + 1;

    /*--------------------------------------
     * Allocating memory for the data field
     *--------------------------------------*/
    p1->Data = (float *) malloc(sizeof(float) * p1->Length);

    /*--------------------------------------------
     * Continuing with the next data structure or
     * close the last structure with NULL
     *--------------------------------------------*/
    if (i < NumChannels - 1) {
      p2 = (MelFB_Window *) malloc(sizeof(MelFB_Window));
      p1->Next = p2;
      p1 = p2;
    } else
      p1->Next = NULL;
  }
  return;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: ComputeTriangle
 *
 * PURPOSE:       Computes and stores FFT window coefficients (triangle points)
 *                into initialized chained list of FFT window structures
 *
 * INPUT:
 *   FirstWin     Pointer to first FFT window structure
 *
 * OUTPUT
 *                Chained list of FFT window data structures with correct
 *                window coefficients
 *
 * RETURN VALUE
 *   none
 *
 *---------------------------------------------------------------------------*/
void 
ComputeTriangle(MelFB_Window * FirstWin)
{
  MelFB_Window *p1;

  int low_part_length, hgh_part_length, TmpInt = 0, i, j;

  p1 = FirstWin;
  j = 0;
  while (p1) {
    low_part_length = p1->Next ?
      p1->Next->StartingPoint - p1->StartingPoint + 1 :
      TmpInt - p1->StartingPoint + 1;

    hgh_part_length = p1->Length - low_part_length + 1;

    /*--------------------------------------
     * Lower frequency part of the triangle
     *--------------------------------------*/
    for (i = 0; i < low_part_length; i++)
      p1->Data[i] = (float) (i + 1) / low_part_length;

    /*---------------------------------------
     * Higher frequency part of the triangle
     *---------------------------------------*/
    for (i = 1; i < hgh_part_length; i++)
      p1->Data[low_part_length + i - 1] =
	(float) (hgh_part_length - i) / hgh_part_length;

    /*------------------------------------------------------
     * Store upper edge (for calculating the last triangle)
     *------------------------------------------------------*/
    TmpInt = p1->StartingPoint + p1->Length - 1;

    /*------------------
     * Next triangle...
     *------------------*/
    p1 = p1->Next;
  }
  return;
}
