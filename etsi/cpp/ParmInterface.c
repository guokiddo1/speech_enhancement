/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: ParmInterface.c
 * PURPOSE: Processing one input frame in DoAdvProcess():
 *          DoNoiseSup(), DoWaveProc(), DoCompCeps(), DoPostProc(),
 *          and DoVADProc() are called.
 *
 *-------------------------------------------------------------------------------*/
/*-----------------
 * File Inclusions
 *-----------------*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ParmInterface.h"

/*--------------------------------------------
 * The followings include the module-specific
 * and exported methods (or functions).
 *--------------------------------------------*/
#include "NoiseSupExports.h"
#include "WaveProcExports.h"
#include "CompCepsExports.h"
#include "PostProcExports.h"
#include "16kHzProcExports.h"
#include "VADExports.h"

extern FILE* fp_denoised; /* Added by mysellsf  */

/*----------------------------------------------------------------------------
 * FUNCTION NAME: AdvProcessAlloc
 *
 * PURPOSE:       Memory allocation of front end parameter structure
 *
 * INPUT:
 *   SamplingFrequency
 *
 * OUTPUT:
 *   pFEParX      Pointer to front end parameter structure
 *
 * RETURN VALUE:
 *   pFEParX      Pointer to front end parameter structure
 *
 *---------------------------------------------------------------------------*/
FEParamsX *
AdvProcessAlloc(int SamplingFrequency)
{
    FEParamsX *pFEParX = NULL;
    pFEParX = (FEParamsX *) calloc(1, sizeof(FEParamsX));

    pFEParX->DoNoiseSupAlloc = DoNoiseSupAlloc;
    pFEParX->DoNoiseSupInit = DoNoiseSupInit;
    pFEParX->DoNoiseSup = DoNoiseSup;
    pFEParX->DoNoiseSupDelete = DoNoiseSupDelete;

    pFEParX->DoWaveProcAlloc = DoWaveProcAlloc;
    pFEParX->DoWaveProcInit = DoWaveProcInit;
    pFEParX->DoWaveProc = DoWaveProc;
    pFEParX->DoWaveProcDelete = DoWaveProcDelete;

    pFEParX->DoCompCepsAlloc = DoCompCepsAlloc;
    pFEParX->DoCompCepsInit = DoCompCepsInit;
    pFEParX->DoCompCeps = DoCompCeps;
    pFEParX->DoCompCepsDelete = DoCompCepsDelete;

    pFEParX->DoPostProcAlloc = DoPostProcAlloc;
    pFEParX->DoPostProcInit = DoPostProcInit;
    pFEParX->DoPostProc = DoPostProc;
    pFEParX->DoPostProcDelete = DoPostProcDelete;

    pFEParX->DoVADProc = DoVADProc;
    pFEParX->DoVADDelete = DoVADDelete;
    pFEParX->DoVADInit = DoVADInit;
    pFEParX->DoVADAlloc = DoVADAlloc;
    pFEParX->DoVADFlush = DoVADFlush;

    /*----------------------------------------------
     * allocate module-specific expansion structure
     *----------------------------------------------*/
    if (pFEParX->DoNoiseSupAlloc != NULL)
        pFEParX->NSX = pFEParX->DoNoiseSupAlloc();
    if (pFEParX->DoWaveProcAlloc != NULL)
        pFEParX->WPX = pFEParX->DoWaveProcAlloc();
    if (pFEParX->DoCompCepsAlloc != NULL)
        pFEParX->CCX = pFEParX->DoCompCepsAlloc();
    if (pFEParX->DoPostProcAlloc != NULL)
        pFEParX->PPX = pFEParX->DoPostProcAlloc();
    if (pFEParX->DoVADAlloc != NULL)
        pFEParX->VADX = pFEParX->DoVADAlloc();

    pFEParX->SamplingFrequency = SamplingFrequency;

    if (SamplingFrequency == 16000) {
        pFEParX->SamplingFrequency = 8000;
        pFEParX->Do16kHzProc = 1;
        pFEParX->Do16kProcessing = Do16kProcessing;
        pFEParX->Do16kProcDelete = Do16kProcDelete;
        pFEParX->Do16kProcInit = Do16kProcInit;
        pFEParX->Do16kProcAlloc = Do16kProcAlloc;
        if (pFEParX->Do16kProcAlloc != NULL)
            pFEParX->pData16k = pFEParX->Do16kProcAlloc();
    } else {
        pFEParX->Do16kHzProc = 0;
        pFEParX->Do16kProcessing = NULL;
        pFEParX->Do16kProcDelete = NULL;
        pFEParX->Do16kProcInit = NULL;
        pFEParX->Do16kProcAlloc = NULL;
    }

    return pFEParX;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: AdvProcessInit
 *
 * PURPOSE:       Initialization of front end parameter structure
 *
 * INPUT:
 *   pFEParX      Pointer to front end parameter structure
 *
 * OUTPUT:
 *   pFEParX      Pointer to front end parameter structure
 *
 * RETURN VALUE:
 *
 *---------------------------------------------------------------------------*/
void 
AdvProcessInit(FEParamsX * pFEParX)
{
    if (pFEParX->SamplingFrequency == SAMPLING_FREQ_1 * 1000) {
        pFEParX->FrameLength = FRAME_LENGTH;
        pFEParX->FrameShift = FRAME_SHIFT;
        pFEParX->FFTLength = FFT_LENGTH;
        pFEParX->StartingFrequency = STARTING_FREQ;
    } else if (pFEParX->SamplingFrequency == SAMPLING_FREQ_3 * 1000) {
        pFEParX->FrameLength = FRAME_LENGTH;
        pFEParX->FrameShift = FRAME_SHIFT;
        pFEParX->FFTLength = FFT_LENGTH;
        pFEParX->StartingFrequency = STARTING_FREQ;
    } else {
        fprintf(stderr, "ERROR:   Invalid sampling frequency '%d'!\r\n",
                pFEParX->SamplingFrequency);
        exit(0);
    }

    {
        int FrameShift = pFEParX->FrameShift;

        if (pFEParX->Do16kHzProc == 1)
            FrameShift *= 2;
        pFEParX->CurFrame = NULL;
	pFEParX->CurFrame = (X_FLOAT32 *)calloc(1, sizeof(pFEParX->CurFrame[0]) * FrameShift);
        if (pFEParX->CurFrame == NULL) {
            fprintf(stderr, "ERROR:   Memory allocation error occured!\r\n");
            exit(0);
        }
    }
    if (pFEParX->DoNoiseSupInit != NULL)
        pFEParX->DoNoiseSupInit(pFEParX);
    if (pFEParX->DoWaveProcInit != NULL)
        pFEParX->DoWaveProcInit(pFEParX);
    if (pFEParX->DoCompCepsInit != NULL)
        pFEParX->DoCompCepsInit(pFEParX);
    if (pFEParX->DoPostProcInit != NULL)
        pFEParX->DoPostProcInit(pFEParX);
    if (pFEParX->DoVADInit != NULL)
        pFEParX->DoVADInit(pFEParX);
    if (pFEParX->Do16kHzProc == 1) {
        if (pFEParX->Do16kProcInit != NULL)
            pFEParX->Do16kProcInit(pFEParX);
    }
    pFEParX->denoisedBuf =
        BufInAlloc(pFEParX->FrameLength +                         //at least FrameLength
                   (pFEParX->FrameLength % pFEParX->FrameShift) + //integer number times FrameShift
                   1);                                            //for pre-emphasis
    pFEParX->offsetDenoisedFrame = -pFEParX->FrameLength;
    pFEParX->FrameCounter = 0;
    pFEParX->NonZeroFrameOnset = 0;
    pFEParX->ZeroFrameCounter = 0;
    pFEParX->NbSamplesToRead = pFEParX->FrameShift;
    if (pFEParX->Do16kHzProc == 1)
        pFEParX->NbSamplesToRead = 2 * pFEParX->FrameShift;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoAdvProcess
 *
 * PURPOSE:       Processing for a frame
 *
 * INPUT:
 *   CurrentFrame Current Frame to process
 *   pFEParX      Pointer to front end parameter structure
 * OUTPUT:
 *   FeatureBuffer Buffer of features
 *
 * RETURN VALUE:
 *   XTRUE        If a feature frame is output
 *   FALSE        If no frame is output (due to delay)
 *
 *---------------------------------------------------------------------------*/
BOOLEAN 
DoAdvProcess(FILE_TYPE * CurrentFrame,
             FILE_TYPE * DenoiseBuffer,
             X_FLOAT32 * FeatureBuffer,
             FEParamsX * pFEParX)
{
    int i;
    int FrameShift;
    X_FLOAT32 frameBuf[FRAME_BUF_SIZE + HP16k_MEL_USED];
    X_FLOAT32 *CurFrame;
    X_FLOAT32 *ptDenoised;

    float FrameCheck;

    long NonZeroFrameOnset;
    long ZeroFrameCounter;

    ZeroFrameCounter = pFEParX->ZeroFrameCounter;
    NonZeroFrameOnset = pFEParX->NonZeroFrameOnset;
    FrameShift = pFEParX->FrameShift;
    if (pFEParX->Do16kHzProc == 1)
        FrameShift = 2 * pFEParX->FrameShift;

    /*-----------------------------------------------------------------
     * get the pointeur where to input new Frameshift denoised samples
     * shift left data in denoisedBuf by FrameShift
     * and returns pointer to the new available part of
     * denoisedBuf => (denoisedBuf->size - FrameShift + 1)
     *-----------------------------------------------------------------*/
    ptDenoised = BufInShiftToPut(pFEParX->denoisedBuf, pFEParX->FrameShift);

    CurFrame = pFEParX->CurFrame;

    /*-------------------------------------------------
     * convert input samples from X_INT16 to X_FLOAT32
     *-------------------------------------------------*/
    FrameCheck = 0.0;
    for (i = 0; i < FrameShift; i++) {
        CurFrame[i] = (X_FLOAT32) (CurrentFrame[i]);
        FrameCheck += CurrentFrame[i] * CurrentFrame[i];
    }

    if (((int) FrameCheck != 0) || (NonZeroFrameOnset != 0)) {
        pFEParX->NonZeroFrameOnset = 1;

        /*------------------------------------------------------------------------
         * WARNING: 16kHz processing
         * curFrame [0 to (pFEParX->FrameShift - 1)] contains 0-4kHz data
         * pData16k->dataHP [0 to (pFEParX->FrameShift - 1)] contains
         * 4-8kHz data shifted to 0-4kHz freq. range
         *-------------------------------------------------------------------------*/
        if (pFEParX->Do16kHzProc)
            pFEParX->Do16kProcessing(CurFrame, pFEParX->pData16k, FrameShift);

        if (pFEParX->DoNoiseSup(CurFrame, ptDenoised, pFEParX)) {
            // FILE_TYPE s;/* Added by myself */
            for (i=0 ; i<80 ; i++) // #define NS_FRAME_SHIFT (X_INT16)(80)
            {                    
                DenoiseBuffer[i] = (FILE_TYPE) ptDenoised [i];
                // fwrite (&s, sizeof (s), 1, fp_denoised); /*Added by myself */
            }

            /*---------------------------------------------------------------------------
             * the noise suppression module outputs "FrameShift" samples in denoisedBuf.
             * following modules will be processed only when denoisedBuf will be full.
             *---------------------------------------------------------------------------*/
			/*
            if (pFEParX->offsetDenoisedFrame < 0)
                pFEParX->offsetDenoisedFrame += pFEParX->FrameShift;
            if (pFEParX->offsetDenoisedFrame >= 0) {
                BOOLEAN rVal = FALSE;


                BufInGetLast(pFEParX->denoisedBuf, frameBuf, pFEParX->FrameLength + pFEParX->offsetDenoisedFrame + 1);

                /*---------------------
                 * waveform processing
                 *---------------------
                if (pFEParX->DoWaveProc)
                    rVal = pFEParX->DoWaveProc(frameBuf + 1, pFEParX);

                /*----------------------
                 * cepstrum calculation
                 *----------------------
                if (pFEParX->DoCompCeps)
                    rVal = pFEParX->DoCompCeps(frameBuf + 1, FeatureBuffer, pFEParX);

                /*-----------------
                 * post processing
                 *-----------------
                if (pFEParX->DoPostProc)
                    rVal = pFEParX->DoPostProc(FeatureBuffer, pFEParX);

                /*-----
                 * VAD
                 *-----
                if (pFEParX->DoVADProc)
                    return pFEParX->DoVADProc(FeatureBuffer, pFEParX);
                else
                    return rVal;

            } else
                return FALSE;
			*/
        } else
            return FALSE;
    } else {
        pFEParX->ZeroFrameCounter++;

        /*-------------------------
         * Create null MFCC vector
         *-------------------------*/
        for (i = 0; i <= NUM_CEP_COEFF; i++)
            FeatureBuffer[i] = 0.0;

        /*---------------------
         * Update VAD variable
         *---------------------*/
        pFEParX->VAD = NON_SPEECH_FRAME;

        return TRUE;
    }
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: FlushAdvProcess
 *
 * PURPOSE:       Flushing of feature frames
 *
 * INPUT:
 *
 * OUTPUT:
 *   FeatureBuffer Buffer of features
 *
 * RETURN VALUE:
 *   XTRUE        If a feature frame is output
 *   FALSE        If no frame is output (due to delay)
 *
 *---------------------------------------------------------------------------*/
BOOLEAN 
FlushAdvProcess(X_FLOAT32 * FeatureBuffer,
                FEParamsX * pFEParX)
{
    if (pFEParX->DoVADFlush)
        return pFEParX->DoVADFlush(FeatureBuffer, pFEParX);
    else
        return FALSE;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: AdvProcessDelete
 *
 * PURPOSE:       Memory free of front end parameter structure
 *
 * INPUT:
 *   ppFEParX     Pointer to pointer to front end parameter structure
 *
 * OUTPUT:
 *
 * RETURN VALUE:
 *
 *---------------------------------------------------------------------------*/
void 
AdvProcessDelete(FEParamsX ** ppFEParX)
{
    if ((*ppFEParX)->DoNoiseSupDelete != NULL)
        (*ppFEParX)->DoNoiseSupDelete((*ppFEParX)->NSX);
    if ((*ppFEParX)->DoWaveProcDelete != NULL)
        (*ppFEParX)->DoWaveProcDelete((*ppFEParX)->WPX);
    if ((*ppFEParX)->DoCompCepsDelete != NULL)
        (*ppFEParX)->DoCompCepsDelete((*ppFEParX)->CCX);
    if ((*ppFEParX)->DoPostProcDelete != NULL)
        (*ppFEParX)->DoPostProcDelete((*ppFEParX)->PPX);
    if ((*ppFEParX)->Do16kProcDelete != NULL)
        (*ppFEParX)->Do16kProcDelete((*ppFEParX)->pData16k);
    if ((*ppFEParX)->DoVADDelete != NULL)
        (*ppFEParX)->DoVADDelete((*ppFEParX)->VADX);

    BufInFree((*ppFEParX)->denoisedBuf);

    if ((*ppFEParX)->CurFrame != NULL) {
        free((*ppFEParX)->CurFrame);
    }

    free(*ppFEParX);
    *ppFEParX = NULL;
}
