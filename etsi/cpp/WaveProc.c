/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: WaveProc.c
 * PURPOSE: Apply SNR-dependent Waveform Processing on the de-noised frame.
 *
 *-------------------------------------------------------------------------------*/
/*-----------------
 * File Inclusions
 *-----------------*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "ParmInterface.h"
#include "WaveProcExports.h"

/*------------------------
 * Definitions and Macros
 *------------------------*/
#define WP_EPS                0.2
#define WP_PERCENT_WIDTH_PAR 80

struct WaveProcStructX {
  X_INT16 FrameLength;
  X_INT16 perCentWidthPar;
  X_FLOAT32 eps;
};

/*------------
 * Prototypes
 *------------*/
static void sort_it(X_INT16 * MaxPos, X_INT16 NumOfMax);
static void TeagerEng(X_FLOAT32 * TeagerWindow, WaveProcStructX * WPX, X_FLOAT32 * Data);
static void GetTeagerFilter(X_FLOAT32 * TeagerWindow, WaveProcStructX * WPX, X_INT16 * MaxPos);
static X_INT16 GetMaximaPositions(X_INT32 * teager, X_INT16 * MaxPos, X_INT16 Length);

/*-----------
 * Functions
 *-----------*/
/*----------------------------------------------------------------------------
 * FUNCTION NAME: sort_it
 *
 * PURPOSE: Sort positions of maxima
 *
 * INPUT:
 *  MaxPos         Pointer to positions of maxima
 *  NumOfMax       Number of maxima
 *
 * OUTPUT:
 *  Sorted maxima positions in MaxPos
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
static void 
sort_it(X_INT16 * MaxPos, X_INT16 NumOfMax)
{
  X_INT16 ii;
  X_INT16 jj;
  X_INT16 jj_aux = 0;
  X_INT16 max_aux[10];

  for (ii = 0; ii < NumOfMax; ii++)
    max_aux[ii] = 0;

  for (ii = 0; ii < NumOfMax; ii++) {
    for (jj = 0; jj < NumOfMax; jj++)
      if (MaxPos[jj] > max_aux[ii]) {
	max_aux[ii] = MaxPos[jj];
	jj_aux = jj;
      }
    MaxPos[jj_aux] = 0;
  }

  for (ii = 0; ii < NumOfMax; ii++)
    MaxPos[ii] = max_aux[NumOfMax - 1 - ii];

  return;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: GetMaximaPositions
 *
 * PURPOSE: Find maxima in energy contour
 *
 * INPUT:
 *  teager         Pointer to Teager energy contour
 *  MaxPos         Pointer to output maxima positions
 *  Length         Length of Teager energy contour
 *
 * OUTPUT:
 *  Maxima positions in MaxPos
 *
 * RETURN VALUE:
 *  Number of maxima (int)
 *
 *---------------------------------------------------------------------------*/
static X_INT16 
GetMaximaPositions(X_INT32 * teager, X_INT16 * MaxPos, X_INT16 Length)
{
  X_INT16 i;
  X_INT16 counterR;
  X_INT16 counterL;
  X_INT16 maxFound = 0;
  X_INT16 maxR[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  X_INT16 maxL[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  X_INT32 teager_max = 0;

  for (i = 0; i < Length; i++) {
    if (teager[i] > teager_max) {
      maxFound = 1;
      teager_max = teager[i];
      maxR[0] = i;
      maxL[0] = i;
    }
  }

  counterR = 0;
  counterL = 0;

  if (maxFound == 1) {
    while (((maxR[counterR] + 25) < Length) && (maxFound == 1)) {
      teager_max = 0;
      maxFound = 0;
      for (i = 25; i < 80; i++)
	if ((maxR[counterR] + i) < Length)
	  if (teager[maxR[counterR] + i] >= teager_max) {
	    maxFound = 1;
	    teager_max = teager[maxR[counterR] + i];
	    maxR[counterR + 1] = maxR[counterR] + i;
	  }
      if (maxFound == 1)
	counterR++;
    }

    maxFound = 1;
    while (((maxL[counterL] - 25) > 0) && (maxFound == 1)) {
      teager_max = 0;
      maxFound = 0;
      for (i = 25; i < 80; i++)
	if ((maxL[counterL] - i) > -1)
	  if (teager[maxL[counterL] - i] >= teager_max) {
	    maxFound = 1;
	    teager_max = teager[maxL[counterL] - i];
	    maxL[counterL + 1] = maxL[counterL] - i;
	  }
      if (maxFound == 1)
	counterL++;
    }
    for (i = 0; i <= counterR; i++)
      MaxPos[i] = maxR[i];
    for (i = (counterR + 1); i <= (counterR + counterL); i++)
      MaxPos[i] = maxL[i - counterR];

    sort_it(MaxPos, counterR + counterL + 1);

    return (counterR + counterL + 1);
  } else {
    return (0);
  }
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: TeagerEng
 *
 * PURPOSE:  Compute Teager energy from input data (time domain)
 *           engT[t] = data[t]*data[t]-data[t-1]*data[t+1]
 *
 * INPUT:
 *  TeagerWindow       Pointer to output Teager energy
 *  WPX                Pointer to waveform processing data stuct
 *  Data               Pointer to input waveform
 *
 * OUTPUT:
 *  Teager energy is stored in TeagerWindow
 *
 * RETURN VALUE:
 *  none
 *
 *---------------------------------------------------------------------------*/
static void 
TeagerEng(X_FLOAT32 * TeagerWindow, WaveProcStructX * WPX, X_FLOAT32 * Data)
{
  X_INT16 i;
  X_INT16 fLength = WPX->FrameLength;

  TeagerWindow[0] = fabs(Data[0] * Data[0] - Data[0] * Data[1]);

  for (i = 1; i < fLength - 1; i++)
    TeagerWindow[i] = fabs(Data[i] * Data[i] - Data[i - 1] * Data[i + 1]);

  TeagerWindow[fLength - 1] = fabs(Data[fLength - 1] * Data[fLength - 1] - Data[fLength - 2] * Data[fLength - 1]);
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: GetTeagerFilter
 *
 * PURPOSE:  Construct Teager energy based windowing function
 *
 *
 * INPUT:
 *  TeagerWindow       Pointer to windowing function
 *  WPX                Pointer to waveform processing data stuct
 *  MaxPos             Pointer to maxima positions
 *
 * OUTPUT:
 *  Teager based windowing function is in TeagerWindow
 *
 * RETURN VALUE:
 *  none
 *
 *---------------------------------------------------------------------------*/
static void 
GetTeagerFilter(X_FLOAT32 * TeagerWindow, WaveProcStructX * WPX, X_INT16 * MaxPos)
{
  X_INT16 i, j;
  X_INT16 io, iom1;
  X_INT16 NoM;
  X_INT16 pWP = WPX->perCentWidthPar;
  X_INT16 FrameLength = WPX->FrameLength;

  X_INT32 teagerInt[9];
  X_INT32 teagerSmoothed[200];

  X_FLOAT32 lowVal;
  X_FLOAT32 highVal;

  /*----------------------------
   * smoothing of Teager energy
   *----------------------------*/
  //JM 25 - 09 - 2002 #define FactIntTeager 1.0
  #define FactIntTeager 0.25

  teagerInt[4] = (int) floor(TeagerWindow[0] * FactIntTeager + 0.5);
  teagerInt[0] = teagerInt[1] = teagerInt[2] = teagerInt[3] = teagerInt[4];
  for (i = 1; i < 5; i++) {
    teagerInt[i + 4] = (int) floor(TeagerWindow[i] * FactIntTeager + 0.5);
  }
  teagerSmoothed[0] = teagerInt[0] + teagerInt[1] + teagerInt[2] + teagerInt[3] + teagerInt[4]
    + teagerInt[5] + teagerInt[6] + teagerInt[7] + teagerInt[8];

  for (i = 1, io = 0; i < FrameLength - 4; i++) {
    X_INT32 val = teagerInt[io];

    teagerInt[io] = (int) floor(TeagerWindow[i + 4] * FactIntTeager + 0.5);
    teagerSmoothed[i] = teagerSmoothed[i - 1] - val + teagerInt[io];
    io++;
    if (io == 9)
      io = 0;
  }

  iom1 = io - 1;
  if (io == 0)
    iom1 = 8;

  for (i = FrameLength - 4; i < FrameLength; i++) {
    X_INT32 val = teagerInt[io];

    teagerInt[io] = teagerInt[iom1];
    teagerSmoothed[i] = teagerSmoothed[i - 1] - val + teagerInt[io];
    iom1 = io;
    io++;
    if (io == 9)
      io = 0;
  }

  /*------------
   *
   *------------*/
  NoM = GetMaximaPositions(teagerSmoothed, MaxPos, FrameLength);

  /*------------
   *
   *------------*/
  lowVal = (1 - WPX->eps) / 2.0;
  highVal = (1 + WPX->eps) / 2.0;

  for (i = 0; i < FrameLength; i++)
    TeagerWindow[i] = lowVal;

  if (NoM > 1) {
    for (i = 0; i < (NoM - 1); i++)
      for (j = (MaxPos[i] - 4); j < ((MaxPos[i] - 4) + ((pWP * (MaxPos[i + 1] - MaxPos[i]) + 99) / 100)); j++)
	if (j >= 0)
	  TeagerWindow[j] = highVal;

    for (j = (MaxPos[NoM - 1] - 4); j < ((MaxPos[NoM - 1] - 4) + ((pWP * (MaxPos[NoM - 1] - MaxPos[NoM - 2]) + 99) / 100)); j++)
      if (j < FrameLength)
	TeagerWindow[j] = highVal;
  }
}

/*------------------------
 * start of Encapsulation
 *------------------------*/
/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoWaveProcAlloc
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
extern WaveProcStructX *
DoWaveProcAlloc(void)
{


  WaveProcStructX *WPX = NULL;
  WPX = (WaveProcStructX *) calloc(1, sizeof(WaveProcStructX));
  return WPX;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoWaveProcInit
 *
 * PURPOSE: Initialize waveform processing data struct
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
DoWaveProcInit(FEParamsX * This)
{
  WaveProcStructX *WPX = This->WPX;

  WPX->eps = WP_EPS;
  WPX->FrameLength = This->FrameLength;
  WPX->perCentWidthPar = WP_PERCENT_WIDTH_PAR;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoWaveProcDelete
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
DoWaveProcDelete(WaveProcStructX * WPX)
{
  if (WPX != NULL)
    free(WPX);
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoWaveProc
 *
 * PURPOSE: Performs waveform processing on input waveform
 *
 *
 * INPUT:
 *  Data      Pointer to input waveform
 *  This      Pointer to front-end data struct
 *
 * OUTPUT:
 *  Processed waveform is in Data
 *
 * RETURN VALUE:
 *   True (boolean)
 *
 *---------------------------------------------------------------------------*/
extern BOOLEAN 
DoWaveProc(X_FLOAT32 * Data, FEParamsX * This)
{
  WaveProcStructX *WPX = This->WPX;

  X_INT16 i;
  X_INT16 MaxPos[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  X_FLOAT32 Win;
  X_FLOAT32 CurWin;
  X_FLOAT32 NxtWin;
  X_FLOAT32 Energy;
  X_FLOAT32 *TeagerWindow;

  /*------------------
   * low energy check
   *------------------*/
  Energy = 0.0;
  for (i = 0; i < WPX->FrameLength; i++)
    Energy += (Data[i] * Data[i]);

  if (Energy >= 100.0) {
    /*-------------------
     * memory allocation
     *-------------------*/
    TeagerWindow = (X_FLOAT32 *) malloc(sizeof(X_FLOAT32) * WPX->FrameLength);
    if (TeagerWindow == NULL) {
      fprintf(stderr, "ERROR:   Memory allocation error occured!\r\n");
      exit(0);
    }
    /*-----------------------
     * compute Teager energy
     *-----------------------*/
    TeagerEng(TeagerWindow, WPX, Data);

    /*-------------------
     * get Teager filter
     *-------------------*/
    GetTeagerFilter(TeagerWindow, WPX, MaxPos);

    /*---------------------
     * apply Teager filter
     *---------------------*/
    CurWin = TeagerWindow[0];

    for (i = 0; i < WPX->FrameLength - 1; i++) {
      NxtWin = TeagerWindow[i + 1];
      Win = CurWin + NxtWin;
      Data[i] *= Win;
      CurWin = NxtWin;
    }

    NxtWin = TeagerWindow[i];
    Win = CurWin + NxtWin;
    Data[i] *= Win;

    free(TeagerWindow);
  }
  return TRUE;
}
