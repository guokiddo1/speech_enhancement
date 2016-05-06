/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: VAD.c
 * PURPOSE: Logic for the frame dropping VAD. VAD flag is set
 *          according to parameters (SpeechFoundSpec, SpeechFoundMel,
 *          SpeechFoundVar) calculated in NoiseSup.c.
 *
 *-------------------------------------------------------------------------------*/
/*-----------------
 * File Inclusions
 *-----------------*/
#include <stdlib.h>
#include <stdio.h>

#include "ParmInterface.h"
#include "VADExports.h"

/*----------------
* VAD Definitions
*-----------------*/
struct VADStructX {
  int Focus;
    //Position of circular buffer
  int HangOver;
    //Hangover length
  int FlushFocus;
    //Position in circular buffer when emptying at end(uses Focus as start / endpoint)
  int H_CountDown;
    //Main hangover countdown
  int V_CountDown;
    //Short hangover countdown
  int DangerWillRobinson;
    //No output flag
  float **OutBuffer;
    //
};

#define PRINTVAD    0
#define PRINTNULLS  0
#define BUFFER_SIZE 7	   // Number of frames in analysis buffer
#define HANGOVERV   23	   // Number of frames by which VAD will stay on after speech
#define H_TRIGGER   4	   // Number of speech frames required to trigger hangover NB Must be <= BUFFER_SIZE!!
#define VANGOVER    5	   // Number of frames for short hangover, triggered by 3 frames of speech
#define V_TRIGGER   3	   // Number of speech frames required to trigger VAD
#define SANGOVER   50	   // Number of frames for safety hangover


/*----------------------------------------------------------------------------
 * FUNCTION NAME: focalpoint
 *
 * PURPOSE: CONTROL VAD BUFFER
 *
 * INPUT:
 *	Focus		current focal point
 *      offset		change in focal point
 *
 * OUTPUT:
 *	None
 *
 *
 * RETURN VALUE:
 *	temp 		new focal point within Buffer boundary
 *
 *
 *---------------------------------------------------------------------------*/
int 
focalpoint(int Focus, int offset)
{
  int temp;

  temp = Focus + offset;
  if (temp > (BUFFER_SIZE - 1))
    temp = temp - BUFFER_SIZE;
  if (temp < 0)
    temp = temp + BUFFER_SIZE;
  return (temp);
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoVADAlloc
 *
 * PURPOSE: Allocation of memory to VAD structure
 *
 *
 * INPUT:
 *	None
 *
 *
 * OUTPUT:
 *	None
 *
 *
 * RETURN VALUE:
 *	VADX		Pointer to VADStructX
 *
 *
 *---------------------------------------------------------------------------*/
extern VADStructX *
DoVADAlloc(void)
{
  VADStructX *VADX = NULL;
  VADX = (VADStructX *)calloc(1, sizeof(VADStructX));
  return VADX;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoVADInit
 *
 * PURPOSE: Initialises VAD structure
 *
 *
 * INPUT:
 *	This		Pointer to front end parameter structure
 *
 *
 * OUTPUT:
 *	This		Pointer to front end parameter structure
 *
 *
 * RETURN VALUE:
 *	None
 *
 *
 *---------------------------------------------------------------------------*/
extern void 
DoVADInit(FEParamsX * This)
{
  int i;
  VADStructX *VADX = This->VADX;

  VADX->Focus = 0;
  VADX->HangOver = HANGOVERV;
  VADX->H_CountDown = 0;
  VADX->V_CountDown = 0;
  VADX->DangerWillRobinson = 1;
  VADX->FlushFocus = -1;

  VADX->OutBuffer = NULL;
  VADX->OutBuffer = (float **)calloc(BUFFER_SIZE, sizeof(VADX->OutBuffer[0]));
  if (VADX->OutBuffer == NULL) {
    fprintf(stderr, "Allocation error Outbuffer\n");
    exit(0);
  }
  for (i = 0; i < BUFFER_SIZE; i++) {
    VADX->OutBuffer[i] = NULL;
    VADX->OutBuffer[i] = (float *)calloc(NUM_CEP_COEFF + 2, sizeof(VADX->OutBuffer[0][0]));
    if (VADX->OutBuffer[i] == NULL) {
      fprintf(stderr, "Allocation error Outbuffer\n");
      exit(0);
    }
  }
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoVADDelete
 *
 * PURPOSE: Free VAD structure
 *
 *
 * INPUT:
 *	VADX		Pointer to VAD structure
 *
 *
 * OUTPUT:
 *	VADX		Pointer to empty VAD structure
 *
 *
 * RETURN VALUE:
 *	None
 *
 *
 *---------------------------------------------------------------------------*/
extern void 
DoVADDelete(VADStructX * VADX)
{
  int i;

  for (i = 0; i < BUFFER_SIZE; i++) {
    if (VADX->OutBuffer[i] != NULL)
      free(VADX->OutBuffer[i]);
  }

  if (VADX->OutBuffer != NULL)
    free(VADX->OutBuffer);

  if (VADX != NULL)
    free(VADX);
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoVADProc
 *
 * PURPOSE: Identifies speech and non-speech frames
 *
 *
 * INPUT:
 *	FeatureBuffer		Buffer of features
 *	This			Pointer to front end parameter structure
 *
 *
 * OUTPUT:
 *	FeatureBuffer		Buffer of features
 *	This			Pointer to front end parameter structure
 *
 *
 * RETURN VALUE:
 *	TRUE 		Once buffer is full
 *	FALSE		For first 10 frames
 *
 *
 *---------------------------------------------------------------------------*/
extern BOOLEAN 
DoVADProc(X_FLOAT32 * FeatureBuffer, FEParamsX * This)
{
  VADStructX *VADX = This->VADX;

  int i;
  int Sum;
  int Trigger;
  int Refocus;

  int Focus = VADX->Focus;
  int HangOver = VADX->HangOver;
  int H_CountDown = VADX->H_CountDown;
  int V_CountDown = VADX->V_CountDown;
  int DangerWillRobinson = VADX->DangerWillRobinson;

  float **OutBuffer = VADX->OutBuffer;

  /*------------------
   *
   *------------------*/
  Focus++;
  if (Focus == BUFFER_SIZE)
    Focus = 0;
  for (i = 0; i < NUM_CEP_COEFF + 1; i++) {
    OutBuffer[Focus][i] = FeatureBuffer[i];
  }

  //!!!VAD Now 3 frames ahead ! !!
    //There is a 3 frame latency between the measures used by the VAD
    // and the feature frame available
    if (This->SpeechFoundSpec || This->SpeechFoundMel || This->SpeechFoundVar || This->SpeechFoundVADNS) {
    OutBuffer[Focus][NUM_CEP_COEFF + 1] = 1;
  } else {
    OutBuffer[Focus][NUM_CEP_COEFF + 1] = 0;
  }

  //+3 due to changes in lookahead
    if (This->FrameCounter > (BUFFER_SIZE + 3)) {
    Sum = 0;
    Trigger = 0;
    for (i = 0; i < BUFFER_SIZE; i++) {
      Refocus = focalpoint(Focus, i + 1);
      if (OutBuffer[Refocus][NUM_CEP_COEFF + 1]) {
	Sum++;
      } else {
	if (Sum > Trigger) {
	  Trigger = Sum;
	}
	Sum = 0;
      }
    }
    if (Sum > Trigger) {
      Trigger = Sum;
    }
    if (Trigger >= H_TRIGGER) {
      H_CountDown = HangOver;
      DangerWillRobinson = 0;
      if (This->FrameCounter <= 35)
	HangOver = SANGOVER;
    }
    if (H_CountDown && Trigger < V_TRIGGER) {
      H_CountDown--;
    }
    if (Trigger >= V_TRIGGER) {
      V_CountDown = VANGOVER;
    }
    if (V_CountDown && Trigger < V_TRIGGER) {
      V_CountDown--;
    }
    Refocus = focalpoint(Focus, 1);

    for (i = 0; i <= NUM_CEP_COEFF + 1; i++)
      FeatureBuffer[i] = OutBuffer[Refocus][i];

    if (V_CountDown || H_CountDown || Trigger >= V_TRIGGER) {
      FeatureBuffer[NUM_CEP_COEFF + 1] = 1.0;
    } else {
      FeatureBuffer[NUM_CEP_COEFF + 1] = 0.0;
      if (PRINTNULLS)
	fprintf(stderr, "%d,", This->FrameCounter - (BUFFER_SIZE + 3));
    }
    This->VAD = (int) FeatureBuffer[NUM_CEP_COEFF + 1];
    VADX->Focus = Focus;
    VADX->HangOver = HangOver;
    VADX->H_CountDown = H_CountDown;
    VADX->V_CountDown = V_CountDown;
    VADX->DangerWillRobinson = DangerWillRobinson;

    return TRUE;
  } else {
    VADX->Focus = Focus;
    VADX->HangOver = HangOver;
    VADX->H_CountDown = H_CountDown;
    VADX->V_CountDown = V_CountDown;
    VADX->DangerWillRobinson = DangerWillRobinson;

    return FALSE;
  }
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoVADFlush
 *
 * PURPOSE: Flush remaining data from VAD Buffer
 *
 *
 * INPUT:
 *	FeatureBuffer		Buffer of features
 *	This			Pointer to front end parameter structure
 *
 *
 * OUTPUT:
 *	FeatureBuffer		Buffer of features
 *	This			Pointer to front end parameter structure
 *
 *
 * RETURN VALUE:
 *	TRUE		When VAD buffer is being flushed
 *	FALSE		When flushing is complete
 *
 *
 *---------------------------------------------------------------------------*/
extern BOOLEAN 
DoVADFlush(X_FLOAT32 * FeatureBuffer, FEParamsX * This)
{
  VADStructX *VADX = This->VADX;

  int i;
  int Sum;
  int Trigger;
  int Refocus;

  int Focus = VADX->Focus;
  int HangOver = VADX->HangOver;
  int H_CountDown = VADX->H_CountDown;
  int V_CountDown = VADX->V_CountDown;
  int DangerWillRobinson = VADX->DangerWillRobinson;

  float **OutBuffer = VADX->OutBuffer;

  //Last frame output from loop
    if (VADX->FlushFocus == -1)
    VADX->FlushFocus = Focus;

  //Move forward one step
    Focus++;
  if (Focus == BUFFER_SIZE)
    Focus = 0;

  if (Focus != VADX->FlushFocus) {
    //Flush remaining data
      // Need to keep the counter incrementing !
      This->FrameCounter++;

    //+3 due to changes in lookahead
      if (This->FrameCounter > (BUFFER_SIZE + 3)) {
      Sum = 0;
      Trigger = 0;
      for (i = 0; i < BUFFER_SIZE; i++) {
	Refocus = focalpoint(Focus, i + 1);
	if (OutBuffer[Refocus][NUM_CEP_COEFF + 1]) {
	  Sum++;
	} else {
	  if (Sum > Trigger) {
	    Trigger = Sum;
	  }
	  Sum = 0;
	}
      }
      if (Sum > Trigger) {
	Trigger = Sum;
      }
      if (Trigger >= H_TRIGGER) {
	H_CountDown = HangOver;
	DangerWillRobinson = 0;
	if (This->FrameCounter <= 35)
	  HangOver = SANGOVER;
      }
      if (H_CountDown && Trigger < V_TRIGGER) {
	H_CountDown--;
      }
      if (Trigger >= V_TRIGGER) {
	V_CountDown = VANGOVER;
      }
      if (V_CountDown && Trigger < V_TRIGGER) {
	V_CountDown--;
      }
      Refocus = focalpoint(Focus, 1);

      for (i = 0; i <= NUM_CEP_COEFF + 1; i++)
	FeatureBuffer[i] = OutBuffer[Refocus][i];

      if (V_CountDown || H_CountDown || Trigger >= V_TRIGGER) {
	FeatureBuffer[NUM_CEP_COEFF + 1] = 1.0;
      } else {
	FeatureBuffer[NUM_CEP_COEFF + 1] = 0.0;
	if (PRINTNULLS)
	  fprintf(stderr, "%d,", This->FrameCounter - (BUFFER_SIZE + 3));
      }
      /*------------------
       * Output result
       *------------------*/
      This->VAD = (int) FeatureBuffer[NUM_CEP_COEFF + 1];
    }
    VADX->Focus = Focus;
    VADX->HangOver = HangOver;
    VADX->H_CountDown = H_CountDown;
    VADX->V_CountDown = V_CountDown;
    VADX->DangerWillRobinson = DangerWillRobinson;

    return TRUE;
  } else {
    return FALSE;
  }
}
