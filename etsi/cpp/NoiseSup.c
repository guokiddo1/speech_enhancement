/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: NoiseSup.c
 * PURPOSE: 1) Apply 2-stage Wiener filter on the input frame.
 *          2) Apply DC offset removal on the output of 2-stage 
 *             Wiener filter.
 *          3) Calculate parameters for the frame dropping VAD (see 
 *             SpeechQSpec(), SpeechQMel(), SpeechQVar()).
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
#include "NoiseSup.h"
#include "16kHzProcExports.h"
#include "NoiseSupExports.h"
#include "MelProcExports.h"
#include "rfft.h"

// extern FILE* fp_denoised; /* Added by mysellsf  */

/*------------------------
 * Definitions and Macros
 *------------------------*/
#define max(a,b) ((a>b)?(a):(b))

typedef struct vad_data_ns VAD_DATA_NS;
typedef struct vad_data_fd VAD_DATA_FD;
typedef struct vad_data_ca VAD_DATA_CA;
typedef struct ns_var      NS_VAR;
typedef struct ns_tmp      NS_TMP;
typedef struct dc_filter   DC_FILTER;
typedef struct gain_fact   GAIN_FACT;
typedef struct buffers     BUFFERS;
typedef struct spectrum    SPECTRUM;

struct dc_filter
{
    X_FLOAT32 lastSampleIn;                  // last input sample of DC offset compensation
    X_FLOAT32 lastDCOut;                     // last output sample of DC offset compensation
} ;

struct vad_data_ns
{
    X_INT32   nbFrame [2];                   // frame number
    X_INT16   flagVAD;                       // VAD flag (1 == SPEECH, 0 == NON SPEECH)
    X_INT16   hangOver;                      // hangover
    X_INT16   nbSpeechFrames;                // nb speech frames (used to set hangover)
    X_FLOAT32 meanEn;                        // mean energy
};

struct vad_data_ca
{
    X_INT16   flagVAD;                       // VAD flag (1 == SPEECH, 0 == NON SPEECH)
    X_INT16   hangOver;                      // hangover
    X_INT16   nbSpeechFrames;                // nb speech frames (used to set hangover)
    X_FLOAT32 meanCA;                        // mean energy
};

struct vad_data_fd
{
    X_FLOAT32 MelMean;			   //
    X_FLOAT32 VarMean;			   //
    X_FLOAT32 AccTest;			   //
    X_FLOAT32 AccTest2;			   //
    X_FLOAT32 SpecMean;			   //
    X_FLOAT32 MelValues[2];		   //
    X_FLOAT32 SpecValues;		           //
    X_FLOAT32 SpeechInVADQ;		   //
    X_FLOAT32 SpeechInVADQ2;		   //
};

struct gain_fact
{
    X_FLOAT32 denEn1 [3];                    // previous denoised frames energies
    X_FLOAT32 lowSNRtrack;                   // low SNR track
    X_FLOAT32 alfaGF;                        // gain factor applied in 2nd stage
};

struct buffers
{
    X_INT32 nbFramesInFirstStage;            // nb frames in first stage
    X_INT32 nbFramesInSecondStage;           // nb frames in second stage
    X_INT32 nbFramesOutSecondStage;          // nb frames out of second stage

    X_FLOAT32 FirstStageInFloatBuffer  [NS_BUFFER_SIZE]; // first stage buffer
    X_FLOAT32 SecondStageInFloatBuffer [NS_BUFFER_SIZE]; // second stage buffer
};

struct spectrum
{
    X_FLOAT32 nSigSE1 [NS_SPEC_ORDER];       // 1st stage noisy signal spectrum estimation
    X_FLOAT32 nSigSE2 [NS_SPEC_ORDER];       // 2nd stage noisy signal spectrum estimation
    X_FLOAT32 noiseSE1 [NS_SPEC_ORDER];      // 1st stage noise spectrum estimation
    X_FLOAT32 noiseSE2 [NS_SPEC_ORDER];      // 2nd stage noise spectrum estimation
    X_FLOAT32 denSigSE1 [NS_SPEC_ORDER];     // 1st stage denoised signal spectrum estimation
    X_FLOAT32 denSigSE2 [NS_SPEC_ORDER];     // 2nd stage denoised signal spectrum estimation

    X_INT16   indexBuffer1;                  // where to enter new PSD, alternatively 0 and 1
    X_INT16   indexBuffer2;                  // where to enter new PSD, alternatively 0 and 1
    X_FLOAT32 PSDMeanBuffer1 [NS_SPEC_ORDER][NS_PSD_MEAN_ORDER]; // 1st stage PSD Mean buffer
    X_FLOAT32 PSDMeanBuffer2 [NS_SPEC_ORDER][NS_PSD_MEAN_ORDER]; // 2nd stage PSD Mean buffer
};

struct ns_var
{
    X_INT16     SampFreq;                    // sampling frequency
    BOOLEAN     Do16kHzProc;                 // do 16 kHz processing ?
    BUFFERS     buffers;                     // signal buffers
    DC_FILTER   prevSamples;                 // previous samples for DC offset removal
    GAIN_FACT   gainFact;                    // gain factorization variables
    VAD_DATA_NS vadNS;                       // VAD for noise suppression data
    VAD_DATA_FD vadFD;                       // VAD for frame dropping data
    VAD_DATA_CA vadCA;                       // VAD for frame dropping data
    SPECTRUM    spectrum;                    // spectrum data
};

struct ns_tmp
{
    MelFB_Window *FirstWindow;               // chained list for Mel filtering
    X_FLOAT32 **melIDCTbasis;                // mel-frequency inverse DCT basis
    X_FLOAT32 IRWindow [NS_FILTER_LENGTH];   // filter impulse response window
    X_FLOAT32 sigWindow [NS_FRAME_LENGTH];   // signal window
    X_FLOAT32 tmpMem [NS_SCRATCH_MEM_SIZE];  // scratch memory
};

struct NoiseSupStructX
{
    NS_VAR    nsVar;                         // non sharable data
    NS_TMP    nsTmp;                         // sharable data
};

/*------------
 * Prototypes
 *------------*/
static void DCOffsetFil (X_FLOAT32 *Data, DC_FILTER *prevSamples, X_INT16 DataLength);
static void DoSigWindowing (X_FLOAT32 *Data, X_FLOAT32 *window, X_INT16 frameLength, X_INT16 FFTLength);
static void FFTtoPSD (X_FLOAT32 *FFTIn, X_FLOAT32 *PSDOut, X_INT16 FFTLength);
static void PSDMean (X_INT16 *indexBuffer, X_FLOAT32 *PSDIn, X_FLOAT32 *PSDOut, X_FLOAT32 *PSDbuffer);
static void ApplyWF (X_FLOAT32 *data, X_FLOAT32 *predata, X_FLOAT32 *filter, X_FLOAT32 *result, X_INT16 frameShift, X_INT16 melOrder);
static void VAD (X_INT16 fstage, VAD_DATA_NS *vadNS, const X_FLOAT32 *newShiftFrame);
static void FilterCalc (X_INT16 fstage, NoiseSupStructX *NSX, X_FLOAT32 *PSDMeaned, X_FLOAT32 *W);
static void DoGainFact (X_INT16 fstage, NoiseSupStructX *NSX, X_FLOAT32 *W);
static void DoFilterWindowing (X_FLOAT32 *filterIRIn, X_FLOAT32 *hanningWindow, X_FLOAT32 *filterIROut);
static X_INT16 SpeechQSpec (FEParamsX *This);
static X_INT16 SpeechQMel (FEParamsX *This);
static X_INT16 SpeechQVar (FEParamsX *This, X_FLOAT32 *W);

/*-----------
 * Functions
 *-----------*/
/*----------------------------------------------------------------------------
 * FUNCTION NAME: DCOffsetFil
 *
 * PURPOSE:       DC offset compensation
 *              
 *
 * INPUT:
 *   *Data        Input frame to filter
 *   *prevSamples Previous samples used in DC offset filtering
 *   DataLength   Number of samples to filter
 *
 * OUTPUT:
 *   *Data        Output filtered frame
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
static void DCOffsetFil (X_FLOAT32 *Data, DC_FILTER *prevSamples, X_INT16 DataLength)
{
    X_INT16 i;
    X_FLOAT32 aux;
    X_FLOAT32 *Prev_x = &(prevSamples->lastSampleIn);
    X_FLOAT32 *Prev_y = &(prevSamples->lastDCOut);

    // y[n] = x[n] - x[n-1] + 0.9990234375 * y[n-1]

    for (i=0 ; i<DataLength ; i++)
    {
        aux = Data[i];
        Data[i] = Data[i] - *Prev_x + 0.9990234375 * *Prev_y;
        *Prev_x = aux;
        *Prev_y = Data[i];
    }
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoSigWindowing
 *
 * PURPOSE:       Signal windowing
 *              
 * INPUT:
 *   *Data        Input frame
 *   *window      Window
 *   frameLength  Frame Length
 *   FFTLength    FFT Length
 *
 * OUTPUT:
 *   *Data        Output frame
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
static void DoSigWindowing (X_FLOAT32 *Data, X_FLOAT32 *window, X_INT16 frameLength, X_INT16 FFTLength)
{
    X_INT16 i;

    // windowing
    for (i=0 ; i<frameLength ; i++)
        Data [i] = Data[i] * window [i];

    // zero padding
    for (i=frameLength ; i<FFTLength ; i++)
        Data [i] = 0.0;

    return;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: FFTtoPSD
 *
 * PURPOSE:       Compute PSD from FFT
 *              
 * INPUT:
 *   *FFTIn       Input Spectrum 
 *   FFTLenth     FFT length
 *
 * OUTPUT:
 *   *PSDOut      Output PSD
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
static void FFTtoPSD (X_FLOAT32 *FFTIn, X_FLOAT32 *PSDOut, X_INT16 FFTLength)
{
    X_INT16 i, j;

    FFTIn[0] = FFTIn[0] * FFTIn[0];
  
    for (i=1,j=FFTLength-1; i<FFTLength/2 ; i++,j--)
    {
        FFTIn[i] = (FFTIn[i] * FFTIn[i] + FFTIn[j] * FFTIn[j]);
    }
  
    FFTIn[i] = FFTIn[i] * FFTIn[i];

    // from 129 to 65 FB
    for (i=0,j=0 ; i<NS_SPEC_ORDER-1 ; i++,j+=2)
    {
        PSDOut[i] = (FFTIn[j] + FFTIn[j+1]) / 2.0;
    }
    PSDOut[i] = FFTIn[j];

    return;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: PSDMean
 *
 * PURPOSE:       Smoothing of the PSD
 *              
 * INPUT:         
 *   *indexBuffer Pointer to buffer index
 *   *PSDIn       Pointer to input PSD
 *   *PSDBuffer   Pointer to PSD buffer
 *
 * OUTPUT:
 *   PSDOut       Output PSD
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
static void PSDMean (X_INT16 *indexBuffer, X_FLOAT32 *PSDIn, X_FLOAT32 *PSDOut, X_FLOAT32 *PSDBuffer)
{
    X_INT16 i, index;

    *indexBuffer = 1 - *indexBuffer;

    for (i=0 ; i<NS_SPEC_ORDER ; i++)
    {
        index = i * NS_PSD_MEAN_ORDER + *indexBuffer;
        PSDBuffer[index] = PSDIn[i];

        index += 1 - 2 * *indexBuffer; 
        PSDOut[i] = (PSDBuffer[index] + PSDIn[i]) / NS_PSD_MEAN_ORDER;
    }
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: ApplyWF
 *
 * PURPOSE:       Convolving input data with denoised filter impulse response
 *              
 * INPUT:
 *   *data        Pointer to input noisy signal
 *   *predata     Pointer to input noisy signal
 *   *filter      Pointer to denoised filter impulse response
 *   frameShift   Frame shift  
 *   melOrder     Impulse response length
 *
 * OUTPUT:
 *   result       Output denoised data 
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
static void ApplyWF (X_FLOAT32 *data, X_FLOAT32 *predata, X_FLOAT32 *filter, X_FLOAT32 *result, X_INT16 frameShift, X_INT16 melOrder)
{
    X_INT16 i, j;

    for (i=0 ; i<frameShift ; i++)
    {
        result[i] = 0.0;
        for (j=-melOrder ; j<=(i<=melOrder ? i : melOrder) ; j++)
        {
            result[i] += (filter[j + melOrder] * data[i-j]);
        }
        for (j=i+1 ; j<=melOrder ; j++)
        {
            result[i] += (filter[j + melOrder] * predata[frameShift - j + i]);
        }
    }
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: VAD
 *
 * PURPOSE:       Voice Ativity Detection
 *              
 * INPUT:
 *   fstage       Stage of two stage noise suppression
 *   *vadNS       Pointer to VAD structure
 *   *newShiftFrame Pointer to new input frame    
 *
 * OUTPUT:        VAD for noise suppression structure updated
 *
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
static void VAD (X_INT16 fstage, VAD_DATA_NS *vadNS, const X_FLOAT32 *newShiftFrame)
{
    X_INT16 i;
    X_INT16 flagVAD;
    X_INT16 hangOver;
    X_INT16 nbSpeechFrames;
    X_INT32 nbFrame;
    X_FLOAT32 meanEn;
    X_FLOAT32 frameEn;
    X_FLOAT32 lambdaLTE;

    nbFrame        = vadNS->nbFrame [fstage];
    meanEn         = vadNS->meanEn;
    flagVAD        = vadNS->flagVAD;
    hangOver       = vadNS->hangOver;
    nbSpeechFrames = vadNS->nbSpeechFrames;

    if (nbFrame < 2147483647) nbFrame++;
    vadNS->nbFrame [fstage] = nbFrame;

    if (fstage == 1) return;

    if (nbFrame < NS_NB_FRAME_THRESHOLD_LTE)
        lambdaLTE = 1 - 1 / (X_FLOAT32) nbFrame;
    else
        lambdaLTE = NS_LAMBDA_LTE_LOWER_E;

    frameEn = 64.0;

    for (i=0 ; i<NS_FRAME_SHIFT ; i++)
        frameEn += newShiftFrame[i] * newShiftFrame[i];

    frameEn = (X_FLOAT32) (0.5 + (log (frameEn / 64.0) / log(2.0)) * 16.0);

    if (((frameEn - meanEn) < NS_SNR_THRESHOLD_UPD_LTE) || (nbFrame < NS_MIN_FRAME))
    {
        if ((frameEn < meanEn) || (nbFrame < NS_MIN_FRAME))
            meanEn += (1 - lambdaLTE) * (frameEn - meanEn);
        else
            meanEn += (1 - NS_LAMBDA_LTE_HIGHER_E) * (frameEn - meanEn);

        if (meanEn < NS_ENERGY_FLOOR)
            meanEn = NS_ENERGY_FLOOR;
    }
    if (nbFrame > 4)
    {
        if ((frameEn - meanEn) > NS_SNR_THRESHOLD_VAD)
        {
            flagVAD = 1;
            nbSpeechFrames++;
        }
        else
        {
            if (nbSpeechFrames > NS_MIN_SPEECH_FRAME_HANGOVER)
                hangOver = NS_HANGOVER;
            nbSpeechFrames = 0;
            if (hangOver != 0)
            {
                hangOver--;
                flagVAD = 1;
            }
            else
                flagVAD = 0;
        }
    }
    vadNS->meanEn         = meanEn;
    vadNS->flagVAD        = flagVAD;
    vadNS->hangOver       = hangOver;
    vadNS->nbSpeechFrames = nbSpeechFrames;

    return;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: FilterCalc
 *
 * PURPOSE:       Computation of noise suppression filter in the frequency domain
 *              
 * INPUT:
 *   fstage       Stage of two stage noise reduction
 *   *NSX         Pointer to noise suppression structure
 *   *PSDMeaned   Pointer to smoothed PSD     
 *
 * OUTPUT:
 *   W            Noise suppression filter
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
static void FilterCalc (X_INT16 fstage, NoiseSupStructX *NSX, X_FLOAT32 *PSDMeaned, X_FLOAT32 *W)
{
    NS_VAR    nsVar = NSX->nsVar;

    X_INT16   i;
    X_INT16   flagVAD;
    X_INT16   nbFrame;

    X_FLOAT32 SNRprio;
    X_FLOAT32 SNRpost;
    X_FLOAT32 lambdaNSE;
    X_FLOAT32 *denSigSE;
    X_FLOAT32 *nSigSE;
    X_FLOAT32 *noiseSE;

    nSigSE    = NSX->nsVar.spectrum.nSigSE1;
    noiseSE   = NSX->nsVar.spectrum.noiseSE1;
    denSigSE  = NSX->nsVar.spectrum.denSigSE1;

    if (fstage == 1)
    {
        nSigSE    = NSX->nsVar.spectrum.nSigSE2;
        noiseSE   = NSX->nsVar.spectrum.noiseSE2;
        denSigSE  = NSX->nsVar.spectrum.denSigSE2;
    }

    nbFrame   = nsVar.vadNS.nbFrame [fstage];
    flagVAD   = nsVar.vadNS.flagVAD;

    /*-------------------------------------------------------
     * Choice of the Noise Estimation according to 2WF stage
     * VAD based NE only at 1st stage of 2WF
     *-------------------------------------------------------*/

    /*--------------------------------
     * non VAD based Noise Estimation
     *--------------------------------*/
    if (fstage == 1)
    {
        // noise estimation in energy
        for (i=0 ; i<NS_SPEC_ORDER ; i++)
            noiseSE[i] *= noiseSE[i];

        if (nbFrame < 11)
        {
            lambdaNSE = 1 - 1 / (X_FLOAT32) nbFrame;
            for (i=0 ; i<NS_SPEC_ORDER ; i++)
            {
                noiseSE[i] = lambdaNSE * noiseSE[i] + (1 - lambdaNSE) * PSDMeaned[i];
            }
        }
        else
        {
            X_FLOAT32 upDate;
            for (i=0 ; i<NS_SPEC_ORDER ; i++)
            {
                upDate = 0.9 + 0.1 * (PSDMeaned [i] / (PSDMeaned [i] + noiseSE [i])) *
                    (1.0 + 1.0 / (1.0 + 0.1 * (PSDMeaned[i] / noiseSE[i])));
                noiseSE [i] *= upDate;
            }
        }

        // store noise estimation values in magnitude
        for (i=0 ; i<NS_SPEC_ORDER ; i++)
        {
            noiseSE [i] = (X_FLOAT32) sqrt (noiseSE[i]);
            if (noiseSE[i] < NS_EPS) noiseSE[i] = NS_EPS;
        }
    }
  
    /*----------------------------------------------------------------------------
     * Spectral estimations of noisy signal are now stored in magnitude in nSigSE
     *----------------------------------------------------------------------------*/
    for (i=0 ; i<NS_SPEC_ORDER ; i++)
    {
        nSigSE[i]    = (X_FLOAT32) sqrt (nSigSE[i]);
        PSDMeaned[i] = (X_FLOAT32) sqrt (PSDMeaned[i]);
    }
  
    /*-------------------------------------------
     * VAD based Noise Estimation (in magnitude)
     *-------------------------------------------*/
    if (fstage == 0)
    {
        if (nbFrame < NS_NB_FRAME_THRESHOLD_NSE)
            lambdaNSE = 1 - 1 / (X_FLOAT32) nbFrame;
        else
            lambdaNSE =  NS_LAMBDA_NSE;

        if (flagVAD == 0)
        {
            for (i=0 ; i<NS_SPEC_ORDER ; i++)
            {
                noiseSE [i] = lambdaNSE * noiseSE [i] + (1 - lambdaNSE) * PSDMeaned[i];
                if (noiseSE[i] < NS_EPS) noiseSE[i] = NS_EPS;
            }
        }
    }

    /*--------------------------------------
     * noise suppression filter calculation
     *--------------------------------------*/
    for (i=0 ; i<NS_SPEC_ORDER ; i++)
    {
        SNRpost = (PSDMeaned[i] / noiseSE[i]) - 1;
        SNRprio = NS_BETA * (denSigSE [i] / noiseSE [i]) + (1 - NS_BETA) * max (0, SNRpost);
        W[i] = SNRprio / (1 + SNRprio);
        SNRprio = W[i] * PSDMeaned[i] / noiseSE[i];
        SNRprio = max (SNRprio, NS_RSB_MIN);
        W[i] = SNRprio / (1 + SNRprio);
        denSigSE [i] = W[i] * nSigSE[i];
    }
  
    return;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoGainFact
 *
 * PURPOSE:       Apply gain factorization to the noise suppression filter
 *              
 * INPUT:
 *   fstage       Stage of the two stage noise reduction
 *   *NSX         Pointer to noise suppression structure     
 *
 * OUTPUT:
 *   W            Noise suppression filter
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
static void DoGainFact (X_INT16 fstage, NoiseSupStructX *NSX, X_FLOAT32 *W)
{
    X_INT16 i;

    X_FLOAT32 averSNR;
    X_FLOAT32 noiseEn;
    X_FLOAT32 lambdaSNR;
    X_FLOAT32 *noiseSE  = NSX->nsVar.spectrum.noiseSE2;
    X_FLOAT32 *denSigSE = NSX->nsVar.spectrum.denSigSE1;

    GAIN_FACT *gainFact = &(NSX->nsVar.gainFact);

    if (fstage == 0)
    {
        gainFact->denEn1 [0] = gainFact->denEn1 [1];
        gainFact->denEn1 [1] = gainFact->denEn1 [2];
        gainFact->denEn1 [2] = 0.0;
        for (i=0 ; i<NS_SPEC_ORDER ; i++) gainFact->denEn1 [2] += denSigSE [i]; // new denEn1
    }
    else
    {
        noiseEn = 0.0;
        for (i=0 ; i<NS_SPEC_ORDER ; i++) noiseEn += noiseSE [i];
        averSNR = (gainFact->denEn1 [0] * gainFact->denEn1 [1] * gainFact->denEn1 [2]) / (noiseEn * noiseEn * noiseEn);

        if (averSNR > 0.00001)
            averSNR = (20 * log10 (averSNR)) / 3.0;
        else
            averSNR = -100.0 / 3.0;

        if ( ((averSNR - gainFact->lowSNRtrack) < 10.0) || (NSX->nsVar.vadNS.nbFrame[fstage] < NS_MIN_FRAME) )
        {
            if (NSX->nsVar.vadNS.nbFrame[fstage] < NS_MIN_FRAME)
                lambdaSNR = 1.0 - 1.0 / (X_FLOAT32) NSX->nsVar.vadNS.nbFrame[fstage];
            else 
            {
                if (averSNR < gainFact->lowSNRtrack)
                    lambdaSNR = 0.95;
                else
                    lambdaSNR = 0.99;
            }
            gainFact->lowSNRtrack += (1.0 - lambdaSNR) * (averSNR - gainFact->lowSNRtrack);
        }

        if (gainFact->denEn1 [2] > 100) // no change if very low signal
        {
            if (averSNR < (gainFact->lowSNRtrack + 3.5))
            {
                gainFact->alfaGF += 0.15;
                if (gainFact->alfaGF > 0.8) gainFact->alfaGF = 0.8;
            }
            else 
            {
                gainFact->alfaGF -= 0.3;
                if (gainFact->alfaGF < 0.1) gainFact->alfaGF = 0.1;
            }
        }

        for (i=0 ; i<WF_MEL_ORDER ; i++)
            W[i] = gainFact->alfaGF * W[i] + (1.0 - gainFact->alfaGF) * 1.0;
    }
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoFilterWindowing
 *
 * PURPOSE:       Filter windowing
 *              
 * INPUT:
 *   *filterIRIn  Pointer to filter impulse response
 *   *hanningWindow Pointer to Hanning window
 *
 * OUTPUT:
 *   filterIROut  Output windowed filter impulse response
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
static void DoFilterWindowing (X_FLOAT32 *filterIRIn, X_FLOAT32 *hanningWindow, X_FLOAT32 *filterIROut)
{
    X_INT16 i, j;

    for (i=NS_HALF_FILTER_LENGTH,j=0 ; i<NS_FILTER_LENGTH ; i++,j++)
    {
        filterIROut [NS_HALF_FILTER_LENGTH+j] = filterIRIn [j] * hanningWindow [i];
        filterIROut [NS_HALF_FILTER_LENGTH-j] = filterIROut [NS_HALF_FILTER_LENGTH+j];
    }
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: SpeechQSpec
 *
 * PURPOSE:       Voice Activity Detection for Frame Dropping
 *              
 * INPUT:
 *   FEParamsX *  Pointer to frontend parameter structure   
 *
 * OUTPUT:        VAD frame dropping structure updated
 *
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
static X_INT16 SpeechQSpec (FEParamsX *This)
{
    NoiseSupStructX *NSX = This->NSX;
    VAD_DATA_FD *vadFD = &(NSX->nsVar.vadFD);
    X_INT16 FrameCounter = This->NSX->nsVar.vadNS.nbFrame[0];
    X_FLOAT32 Acceleration;					

    if (FrameCounter == 1) vadFD->SpecMean = vadFD->SpecValues;

    if (FrameCounter < 15)
    {
        vadFD->AccTest = 1.1 * (vadFD->AccTest * (X_FLOAT32)(FrameCounter - 1) + vadFD->SpecValues) / (X_FLOAT32)(FrameCounter);
        Acceleration = vadFD->SpecValues / vadFD->AccTest;
        if (Acceleration > 2.5)
        {
            vadFD->SpeechInVADQ = 1;	
        }
        if (vadFD->SpeechInVADQ == 0)
        {
            vadFD->SpecMean = max (vadFD->SpecMean, vadFD->SpecValues);
        }
    }
				
    if (vadFD->SpecValues<vadFD->SpecMean*1.5 && vadFD->SpecValues>vadFD->SpecMean*0.75)
    {
        vadFD->SpecMean = (vadFD->SpecMean * 0.8 + vadFD->SpecValues * 0.2);
    }	

    if (vadFD->SpecValues <= vadFD->SpecMean*0.5)
    {
        vadFD->SpecMean = (vadFD->SpecMean * 0.97 + vadFD->SpecValues * 0.03);
    }	

    if (vadFD->SpecValues > vadFD->SpecMean*1.65)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: SpeechQMel
 *
 * PURPOSE:       Voice Activity Detection for Frame Dropping
 *              
 * INPUT:
 *   FEParamsX *  Pointer to frontend parameter structure   
 *
 * OUTPUT:        VAD frame dropping structure updated
 *
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
static X_INT16 SpeechQMel (FEParamsX *This)
{
    NoiseSupStructX *NSX = This->NSX;
    VAD_DATA_FD *vadFD = &(NSX->nsVar.vadFD);
    X_INT16 FrameCounter = This->NSX->nsVar.vadNS.nbFrame[0];
    X_FLOAT32 Acceleration;
    X_FLOAT32 SmoothMel;

    SmoothMel = 0.75 * vadFD->MelValues[1] + 0.25 * vadFD->MelValues[0];
    vadFD->MelValues[0] = vadFD->MelValues[1];

    if (FrameCounter < 15)
    {
        vadFD->MelMean = max (vadFD->MelMean, SmoothMel);
    }

    if (SmoothMel < vadFD->MelMean*1.5 && SmoothMel > vadFD->MelMean*0.75)
    {
        vadFD->MelMean = (vadFD->MelMean * 0.8 + SmoothMel * 0.2);
    }	

    if (SmoothMel <= vadFD->MelMean*0.5)
    {
        vadFD->MelMean = (vadFD->MelMean * 0.97 + SmoothMel * 0.03);
    }	

    Acceleration = SmoothMel / vadFD->MelMean;

    if (SmoothMel > vadFD->MelMean * 3.25)
    {
        return 1;
    }
    else
    {
        return 0;
    }	   	
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: SpeechQVar
 *
 * PURPOSE:       Voice Activity Detection for Frame Dropping
 *              
 * INPUT:
 *   FEParamsX *  Pointer to frontend parameter structure
 *   *W           Pointer to noise suppression filter  
 *
 * OUTPUT:        VAD frame dropping structure updated
 *
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
static X_INT16 SpeechQVar (FEParamsX *This, X_FLOAT32 *W)
{
    NoiseSupStructX *NSX = This->NSX;
    VAD_DATA_FD *vadFD = &(NSX->nsVar.vadFD);
    X_INT16 i; 
    X_INT16 ssize = NS_FFT_LENGTH / 4;
    X_INT16 FrameCounter = This->NSX->nsVar.vadNS.nbFrame[0];
    X_FLOAT32 var = 0.0;
    X_FLOAT32 mean = 0.0;
    X_FLOAT32 specVar = 0.0;
  
    for (i=0 ; i<ssize ; i++)
    {
        mean += W[i];
        var  += W[i] * W[i];
    }
    specVar = (var / ssize)  - mean * mean / (ssize * ssize);

    if (FrameCounter < 15)
    {
        vadFD->VarMean = max (vadFD->VarMean, specVar);
    }

    if (specVar<vadFD->VarMean*1.5 && specVar>vadFD->VarMean*0.85)
    {
        vadFD->VarMean = (vadFD->VarMean * 0.8 + specVar * 0.2);
    }

    if (specVar <= vadFD->VarMean*0.25)
    {
        vadFD->VarMean = (vadFD->VarMean * 0.97 + specVar * 0.03);
    }

    if (specVar > vadFD->VarMean*1.65)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

/*------------------------
 * start of Encapsulation
 *------------------------*/
/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoNoiseSupAlloc
 *
 * PURPOSE:       Memory allocation of noise suppression structure
 *              
 * INPUT:
 *   none     
 *
 * OUTPUT:        Noise suppression structure
 *
 *
 * RETURN VALUE:
 *   NSX          Pointer to noise suppression structure
 *
 *---------------------------------------------------------------------------*/
extern NoiseSupStructX *DoNoiseSupAlloc (void)
{
    NoiseSupStructX *NSX = NULL;
   NSX = (NoiseSupStructX *)calloc (1, sizeof (NoiseSupStructX));

    if (NSX == NULL)
        return NULL;

    return NSX;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoNoiseSupInit
 *
 * PURPOSE:       Initialisation of noise suppression structure
 *              
 * INPUT:
 *   FEParamsX *  Pointer to frontend parameter structure
 *        
 * OUTPUT:        Initialized noise suppression structure
 *
 * RETURN VALUE:
 *   None
 *
 *---------------------------------------------------------------------------*/
extern void DoNoiseSupInit (FEParamsX *This)
{
    X_INT16 i, j;
    NoiseSupStructX *NSX = This->NSX;

    NSX->nsVar.SampFreq    = This->SamplingFrequency;
    NSX->nsVar.Do16kHzProc = This->Do16kHzProc;

    /*---------
     * buffers
     *---------*/
    for (i=0 ; i<NS_BUFFER_SIZE ; i++)
    {
        NSX->nsVar.buffers.FirstStageInFloatBuffer [i]  = 0.0;
        NSX->nsVar.buffers.SecondStageInFloatBuffer [i] = 0.0;
    }

    NSX->nsVar.buffers.nbFramesInFirstStage   = 0;
    NSX->nsVar.buffers.nbFramesInSecondStage  = 0;
    NSX->nsVar.buffers.nbFramesOutSecondStage = 0;

    /*-----------
     * dc_filter
     *-----------*/
    NSX->nsVar.prevSamples.lastSampleIn = 0.0;
    NSX->nsVar.prevSamples.lastDCOut    = 0.0;

    /*-----------
     * gain_fact
     *-----------*/
    for (i=0 ; i<3 ; i++) NSX->nsVar.gainFact.denEn1 [i] = 0.0;
    NSX->nsVar.gainFact.lowSNRtrack = 0.0;
    NSX->nsVar.gainFact.alfaGF = 0.8;

    /*-------------
     * vad_data_ns
     *-------------*/
    NSX->nsVar.vadNS.nbFrame [0]    = NSX->nsVar.vadNS.nbFrame [1] = 0;
    NSX->nsVar.vadNS.hangOver       = 0;
    NSX->nsVar.vadNS.nbSpeechFrames = 0;
    NSX->nsVar.vadNS.meanEn         = 0.0;
    NSX->nsVar.vadNS.flagVAD        = 0;

    NSX->nsVar.vadCA.hangOver       = 0;
    NSX->nsVar.vadCA.nbSpeechFrames = 0;
    NSX->nsVar.vadCA.meanCA         = 0.0;
    NSX->nsVar.vadCA.flagVAD        = 0;

    /*-------------
     * vad_data_fd
     *-------------*/
    NSX->nsVar.vadFD.MelMean = 0.0;
    NSX->nsVar.vadFD.VarMean = 0.0;
    NSX->nsVar.vadFD.AccTest = 0.0;
    NSX->nsVar.vadFD.AccTest2 = 0.0;
    NSX->nsVar.vadFD.SpecMean = 0.0;
    NSX->nsVar.vadFD.MelValues[0] = NSX->nsVar.vadFD.MelValues[1] = 0.0;
    NSX->nsVar.vadFD.SpecValues = 0.0;
    NSX->nsVar.vadFD.SpeechInVADQ = 0.0;
    NSX->nsVar.vadFD.SpeechInVADQ2 = 0.0;

    /*----------
     * spectrum
     *----------*/
    for (i=0 ; i<NS_SPEC_ORDER ; i++)
    {
        NSX->nsVar.spectrum.denSigSE1[i]   = NSX->nsVar.spectrum.denSigSE2[i]   = 0.0;
        NSX->nsVar.spectrum.nSigSE1[i]     = NSX->nsVar.spectrum.nSigSE2[i]     = 0.0;
        NSX->nsVar.spectrum.noiseSE1[i]    = NSX->nsVar.spectrum.noiseSE2[i]    = NS_EPS;
    }

    for (i=0 ; i<NS_SPEC_ORDER; i++)
    {
        for (j=0 ; j<NS_PSD_MEAN_ORDER ; j++)
            NSX->nsVar.spectrum.PSDMeanBuffer1[i][j] = 0.0;
    }

    for (i=0 ; i<NS_SPEC_ORDER ; i++)
    {
        for (j=0 ; j<NS_PSD_MEAN_ORDER ; j++)
            NSX->nsVar.spectrum.PSDMeanBuffer2[i][j] = 0.0;
    }
 
    NSX->nsVar.spectrum.indexBuffer1 = 0;
    NSX->nsVar.spectrum.indexBuffer2 = 0;

    /*--------
     * ns_tmp
     *--------*/
    // Hanning window for spectrum estimation
    for (i=0 ; i<NS_FRAME_LENGTH ; i++)
        NSX->nsTmp.sigWindow[i] = 0.5 - 0.5 * cos ((PIx2 * ((X_FLOAT32) i + 0.5)) / (X_FLOAT32) (NS_FRAME_LENGTH));
  
    // Hanning window for impulse response windowing
    for (i=0 ; i<NS_FILTER_LENGTH ; i++)
        NSX->nsTmp.IRWindow[i] = 0.5 - 0.5 * cos ((PIx2 * ((X_FLOAT32) i + 0.5)) / (X_FLOAT32) (NS_FILTER_LENGTH));

    // mel FB windows
    NSX->nsTmp.FirstWindow = CMelFBAlloc ();
    if (NSX->nsTmp.FirstWindow == NULL)
    {
        fprintf (stderr, "ERROR:   Memory allocation error occured!\r\n");
        exit(0);
    }
    InitMelFBwindows (NSX->nsTmp.FirstWindow, 0.0, (X_FLOAT32)NSX->nsVar.SampFreq, 2*(NS_SPEC_ORDER-1), WF_MEL_ORDER, 1);

    // mel IDCT
    NSX->nsTmp.melIDCTbasis = (X_FLOAT32 **) malloc (sizeof (X_FLOAT32 *) * WF_MEL_ORDER);
    if (NSX->nsTmp.melIDCTbasis == NULL)
    {
        fprintf (stderr, "ERROR:   Memory allocation error occured!\r\n");
        exit(0);
    }

    for (i=0 ; i<WF_MEL_ORDER ; i++)
    {
        NSX->nsTmp.melIDCTbasis[i] = (X_FLOAT32 *) malloc (sizeof (X_FLOAT32) * WF_MEL_ORDER);
        if (NSX->nsTmp.melIDCTbasis[i] == NULL)
        {
            fprintf (stderr, "ERROR:   Memory allocation error occured!\r\n");
            exit(0);
        }
    }

    InitMelIDCTbasis (NSX->nsTmp.melIDCTbasis, NSX->nsTmp.FirstWindow, WF_MEL_ORDER, NSX->nsVar.SampFreq, 2*(NS_SPEC_ORDER-1));
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoNoiseSupDelete
 *
 * PURPOSE:       Memory free of noise suppression structure
 *              
 * INPUT:
 *   *NSX         Pointer to noise suppression structure     
 *
 * OUTPUT:
 *   None
 *
 * RETURN VALUE:
 *   None
 *
 *---------------------------------------------------------------------------*/
extern void DoNoiseSupDelete (NoiseSupStructX *NSX)
{
    X_INT16 i;

    if (NSX != NULL)
    {
        ReleaseMelFBwindows (NSX->nsTmp.FirstWindow);
        free (NSX->nsTmp.FirstWindow);

        for (i=0 ; i<WF_MEL_ORDER ; i++)
            free (NSX->nsTmp.melIDCTbasis[i]);

        free (NSX->nsTmp.melIDCTbasis);

        free (NSX);
    }
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoNoiseSup
 *
 * PURPOSE:       Performs noise suppression on input signal frame
 *              
 * INPUT:
 *   FEParamsX *  Pointer to frontend parameter structure
 *   *InData      Pointer to input signal frame     
 *
 * OUTPUT:
 *   Outdata      Output signal frame
 *
 * RETURN VALUE:
 *   TRUE         If a frame is outputed
 *   FALSE        If no frame outputed (happens at the beginning due to the latency)
 *
 *---------------------------------------------------------------------------*/
extern BOOLEAN DoNoiseSup (X_FLOAT32 *InData, X_FLOAT32 *OutData, FEParamsX *This)
{
    X_INT16 i; 
    X_INT16 fstage;                      // noise suppression filter stage

    NoiseSupStructX *NSX = This->NSX;

    X_FLOAT32 *FirstStageInFloatBuffer  = NSX->nsVar.buffers.FirstStageInFloatBuffer;
    X_FLOAT32 *SecondStageInFloatBuffer = NSX->nsVar.buffers.SecondStageInFloatBuffer;

    X_INT32 *nbFramesInFirstStage     = &(NSX->nsVar.buffers.nbFramesInFirstStage);
    X_INT32 *nbFramesInSecondStage    = &(NSX->nsVar.buffers.nbFramesInSecondStage);
    X_INT32 *nbFramesOutSecondStage   = &(NSX->nsVar.buffers.nbFramesOutSecondStage);

    X_INT16 *indexBuffer = NULL;     //

    X_FLOAT32 *prvFrame = NULL;      //
    X_FLOAT32 *curFrame = NULL;      //
    X_FLOAT32 *nSigSE = NULL;        //
    X_FLOAT32 *PSDMeanBuffer = NULL; //

    X_FLOAT32 *W = NSX->nsTmp.tmpMem + NS_SPEC_ORDER;               // scratch memory
    X_FLOAT32 *filterIR = This->NSX->nsTmp.tmpMem;                  // scratch memory
    X_FLOAT32 *signalIn = This->NSX->nsTmp.tmpMem;                  // scratch memory
    X_FLOAT32 *signalOut = This->NSX->nsTmp.tmpMem + NS_SPEC_ORDER; // scratch memory
    X_FLOAT32 *PSDMeaned = This->NSX->nsTmp.tmpMem;                 // scratch memory

    X_INT16 bufData16kSize;          //

    X_FLOAT32 fb16k [NS_SPEC_ORDER]; //
    X_FLOAT32 *bufferData16k;        //
    X_FLOAT32 *CodeForBands16k;      //
    X_FLOAT32 *tmpWorkSpace16k;      // temporary work space for 16k processing
    X_FLOAT32 *BandsForCoding16k;    //
 
    DataFor16kProc *pData16k;        //
    MelFB_Window *FirstWindow16k;    //

    BOOLEAN dataToProcess;           // data to process ?

    /*------------------
     * 16kHz Processing
     *------------------*/
    if (NSX->nsVar.Do16kHzProc)
    {
        pData16k          = This->pData16k;
        bufferData16k     = Get16k_p_bufferData16k     (pData16k);
        bufData16kSize    = Get16k_bufData16kSize      (pData16k);
        BandsForCoding16k = Get16k_p_BandsForCoding16k (pData16k);
        FirstWindow16k    = Get16k_p_FirstWindow16k    (pData16k);
        CodeForBands16k   = Get16k_p_CodeForBands16k   (pData16k);
        tmpWorkSpace16k = (X_FLOAT32 *) malloc (sizeof (X_FLOAT32) * NS_FFT_LENGTH);
        if (tmpWorkSpace16k == NULL)
        {
            fprintf (stderr, "ERROR:   Memory allocation error occured!\r\n");
            exit(0);
        }

        for (i=0 ; i<NS_FRAME_SHIFT ; i++) 
            bufferData16k [bufData16kSize - NS_FRAME_SHIFT + i] = Get16k_dataHP (pData16k, i);
    }
    else
    {
        pData16k          = NULL;
        bufferData16k     = NULL;
        bufData16kSize    = 0;
        BandsForCoding16k = NULL;
        FirstWindow16k    = NULL;
        CodeForBands16k   = NULL;
        tmpWorkSpace16k   = NULL;
    }

    /*---------------------------------------------------
     * input next shift frame in FirstStageInFloatBuffer
     *---------------------------------------------------*/
    for (i=0 ; i<NS_FRAME_SHIFT ; i++)
        FirstStageInFloatBuffer [NS_DATA_IN_BUFFER + i] = InData[i];
  
    (*nbFramesInFirstStage)++;

    /*------------------------------------
     * Two-Stage Noise Suppression Filter
     *------------------------------------*/
    for (fstage=0 ; fstage<2 ; fstage++) 
    {
        dataToProcess = FALSE;

        /*----------------------------------------------------
         * 1st stage reads data from FirstStageInFloatBuffer
         * 2nd stage reads data from SecondStageInFloatBuffer
         *----------------------------------------------------*/
        if ((fstage == 0) && ((*nbFramesInFirstStage - *nbFramesInSecondStage) > NS_NB_FRAMES_LATENCY))
        {
            /*-------------------
             * process 1st stage
             *-------------------*/
            dataToProcess = TRUE;

            nSigSE = NSX->nsVar.spectrum.nSigSE1;
            PSDMeanBuffer = &(NSX->nsVar.spectrum.PSDMeanBuffer1 [0][0]);
            indexBuffer = &(NSX->nsVar.spectrum.indexBuffer1);

            for (i=0 ; i<NS_FRAME_LENGTH ; i++)
                signalIn[i] = FirstStageInFloatBuffer [NS_ANALYSIS_WINDOW_8K + i];

            if (NSX->nsVar.Do16kHzProc)
            {
                // to get spectrum corresponding to FB in the MFCC block;
                // otherwise there is a difference of 20 samples
                for (i=0 ; i<NS_FRAME_LENGTH ; i++)
                    tmpWorkSpace16k [i] = bufferData16k [NS_ANALYSIS_WINDOW_16K + i];
            }

            prvFrame = &(NSX->nsVar.buffers.FirstStageInFloatBuffer [NS_PRV_FRAME]);
            curFrame = &(NSX->nsVar.buffers.FirstStageInFloatBuffer [NS_CUR_FRAME]);
        }
        else
            if ((fstage == 1) && ((*nbFramesInSecondStage - *nbFramesOutSecondStage) > NS_NB_FRAMES_LATENCY))
            {
                /*-------------------
                 * process 2nd stage
                 *-------------------*/
                dataToProcess = TRUE;

                nSigSE = NSX->nsVar.spectrum.nSigSE2;
                PSDMeanBuffer = &(NSX->nsVar.spectrum.PSDMeanBuffer2 [0][0]);
                indexBuffer = &(NSX->nsVar.spectrum.indexBuffer2);

                for (i=0 ; i<NS_FRAME_LENGTH ; i++)
                    signalIn[i] = SecondStageInFloatBuffer [NS_ANALYSIS_WINDOW_8K + i];

                prvFrame = &(NSX->nsVar.buffers.SecondStageInFloatBuffer [NS_PRV_FRAME]);
                curFrame = &(NSX->nsVar.buffers.SecondStageInFloatBuffer [NS_CUR_FRAME]);
            }

        if (!dataToProcess) continue; // no processing required

        /*-----------------------------------
         * signal windowing and zero padding
         *-----------------------------------*/
        DoSigWindowing (signalIn, NSX->nsTmp.sigWindow, NS_FRAME_LENGTH, NS_FFT_LENGTH);

        /*-----
         * FFT
         *-----*/
        rfft (signalIn, NS_FFT_LENGTH, NS_FFT_ORDER);

        /*-------------------------------------------------------------
         * FFT spectrum (signalIn) -> power spectrum (NSX->nSigSE)
         *-------------------------------------------------------------*/
        FFTtoPSD (signalIn, nSigSE, NS_FFT_LENGTH);

        /*---------------------------------
         * for coding for 16kHz Processing
         *---------------------------------*/
        if (NSX->nsVar.Do16kHzProc && (fstage == 0))
        { 
            X_FLOAT32 aux[3];
      
            // 3 mel-spaced bands 2000-2400-3000-4000 Hz
            // are used for coding 3 bands from 4-8kHz

            GetBandsForCoding16k (nSigSE, aux, NS_SPEC_ORDER);

            for (i=0 ; i<6 ; i++)
                BandsForCoding16k [i] = BandsForCoding16k [i+3];

            for (i=0 ; i<3 ; i++)
            {
                if (aux [i] > NS_SPEC_FLOOR)
                    BandsForCoding16k [6 + i] = (X_FLOAT32) log ((double) aux [i]);
                else
                    BandsForCoding16k [6 + i] = NS_LOG_SPEC_FLOOR;	
            }
        }

        /*---------
         * PSDMean
         *---------*/
        PSDMean (indexBuffer, nSigSE, PSDMeaned, PSDMeanBuffer);

        /*---------------------------
         * VAD for Noise Suppression
         *---------------------------*/
        VAD (fstage, &(NSX->nsVar.vadNS), curFrame); 

        /*--------------------------
         * filter gains calculation
         *--------------------------*/
        FilterCalc (fstage, NSX, PSDMeaned, W);

        /*------------------------
         * VAD for Frame Dropping
         *------------------------*/
        if (fstage == 0)
        {
            This->SpeechFoundVar = SpeechQVar (This, W);
        }

        /*-----------------
         * mel filter bank
         *-----------------*/
        DoMelFB (W, NSX->nsTmp.FirstWindow);

        /*------------------------
         * VAD for Frame Dropping
         *------------------------*/
        if (fstage == 0)
        {
            X_FLOAT32 TempEn = 0.0;

            for(i=0 ; i<WF_MEL_ORDER ; i++) TempEn += (X_FLOAT32)W[i];

            NSX->nsVar.vadFD.SpecValues = TempEn * TempEn - 3.0;
	
            This->SpeechFoundSpec = SpeechQSpec (This);
 
            NSX->nsVar.vadFD.MelValues[1] = (X_FLOAT32)(W[1] + W[2] + W[3]) / 3.0;

            This->SpeechFoundMel = SpeechQMel (This);
        }	 

        /*--------------------
         * gain factorization
         *--------------------*/
        DoGainFact (fstage, NSX, W);

        /*-----------------
         * mel inverse DCT
         *-----------------*/
        DoMelIDCT (W, NSX->nsTmp.melIDCTbasis, WF_MEL_ORDER, WF_MEL_ORDER);
        for (i=1 ; i<WF_MEL_ORDER ; i++) W [2 * WF_MEL_ORDER - 1 - i] = W [i];

        /*------------------
         * filter windowing
         *------------------*/
        DoFilterWindowing (W, NSX->nsTmp.IRWindow, filterIR);

        /*--------------------------------------------------------------
         * apply WF to noisy signal, output samples stored in signalOut
         *--------------------------------------------------------------*/
        ApplyWF (curFrame, prvFrame, filterIR, signalOut, NS_FRAME_SHIFT, NS_HALF_FILTER_LENGTH);

        /*------------------
         * 16kHz Processing 
         *------------------*/
        if (NSX->nsVar.Do16kHzProc && (fstage == 1))
        {
            X_FLOAT32 aux [HP16k_MEL_USED];

            DoSigWindowing (tmpWorkSpace16k, NSX->nsTmp.sigWindow, NS_FRAME_LENGTH, NS_FFT_LENGTH);
            rfft (tmpWorkSpace16k, NS_FFT_LENGTH, NS_FFT_ORDER);
            FFTtoPSD (tmpWorkSpace16k, fb16k, NS_FFT_LENGTH);
            DoMelFB (fb16k, FirstWindow16k);  // WARNING: only fb16k [1:HP16k_MEL_USED] are used
            for (i=0 ; i<HP16k_MEL_USED ; i++) 
                fb16k [i] = fb16k [i+1];

            for (i=0 ; i<HP16k_MEL_USED ; i++)
            {
                if (fb16k [i] > NS_SPEC_FLOOR)
                    aux [i] = log (fb16k [i]);
                else
                    aux [i] = NS_LOG_SPEC_FLOOR;
            }
            CodeBands16k (aux, BandsForCoding16k, CodeForBands16k);
            DoSpecSub16k (fb16k, pData16k, NSX->nsVar.vadNS.nbFrame [1], HP16k_MEL_USED);
        }

        /*---------------------------------------------------
         * update buffer with filtered signal for next stage
         *---------------------------------------------------*/
        if (fstage == 0)
        {
            for (i=0 ; i<NS_FRAME_SHIFT ; i++)
                SecondStageInFloatBuffer[NS_DATA_IN_BUFFER + i] = signalOut[i];

            NSX->nsVar.buffers.nbFramesInSecondStage++;
        }
        else
            if (fstage == 1)
            {
                // FILE_TYPE s;/* Added by */
                
                for (i=0 ; i<NS_FRAME_SHIFT ; i++)
                {                    
                    OutData[i] = signalOut[i];
                    
                    /* /\*Added by myself *\/ */
                    /* s=(FILE_TYPE)OutData[i]; */
                    /* fwrite (&s, sizeof (s), 1, fp_denoised); */
                }

                NSX->nsVar.buffers.nbFramesOutSecondStage++;
            }

        /*------------------------
         * VAD for Frame Dropping
         *------------------------*/
        if (fstage == 0)
        {
            if (NSX->nsVar.vadNS.nbSpeechFrames > 4)
                This->SpeechFoundVADNS = 1;
            else
                This->SpeechFoundVADNS = 0;	 
        }
        This->FrameCounter = This->NSX->nsVar.vadNS.nbFrame[0];
    }

    /*--------------------------------------------
     * shift data in NSX->FirstStageInFloatBuffer
     *--------------------------------------------*/
    if (NSX->nsVar.buffers.nbFramesInFirstStage)
    {
        // copy the content of FirstStageInFloatBuffer [80..319] to [0..239]
        for (i=0 ; i<NS_DATA_IN_BUFFER ; i++)
            FirstStageInFloatBuffer [i] = FirstStageInFloatBuffer [i + NS_FRAME_SHIFT];

        if (NSX->nsVar.Do16kHzProc)
        {
            for (i=0 ; i<(bufData16kSize-NS_FRAME_SHIFT) ; i++)
                bufferData16k [i] = bufferData16k [i + NS_FRAME_SHIFT];
        }
    }

    /*---------------------------------------------
     * shift data in NSX->SecondStageInFloatBuffer
     *---------------------------------------------*/
    if (NSX->nsVar.buffers.nbFramesInSecondStage)
        for (i=0 ; i<NS_DATA_IN_BUFFER ; i++)
            SecondStageInFloatBuffer [i] = SecondStageInFloatBuffer [i + NS_FRAME_SHIFT];

    if (NSX->nsVar.Do16kHzProc)
        free (tmpWorkSpace16k);

    if (NSX->nsVar.buffers.nbFramesOutSecondStage > 0)
    {
        DCOffsetFil (OutData, &(NSX->nsVar.prevSamples), NS_FRAME_SHIFT);

        /*------------------
         * 16kHz Processing 
         *------------------*/
        if (This->Do16kHzProc)
        {
            X_INT16 i;
            X_INT16 hpBandsSize;
            X_FLOAT32 *hpBands, *bufferCodeForBands16k, *CodeForBands16k;
            X_FLOAT32 *bufferCodeWeights, *codeWeights;
      
            hpBandsSize           = Get16k_hpBandsSize (This->pData16k);
            hpBands               = Get16k_p_hpBands (This->pData16k);
            bufferCodeForBands16k = Get16k_p_bufferCodeForBands16k (This->pData16k);
            CodeForBands16k       = Get16k_p_CodeForBands16k (This->pData16k);
            bufferCodeWeights     = Get16k_p_bufferCodeWeights (This->pData16k);
            codeWeights           = Get16k_p_codeWeights (This->pData16k);

            //                          WARNING: 16kHz Processing DM                          
            // CurrentFrame [FrameShift : (FrameShift + HP16k_MEL_USED - 1)] contains HP bands
            for (i=HP16k_MEL_USED ; i<hpBandsSize ; i++)
                hpBands [i - HP16k_MEL_USED] = hpBands [i];

            for (i=0 ; i<HP16k_MEL_USED ; i++)
                hpBands [hpBandsSize - HP16k_MEL_USED + i] = fb16k [i];

            for (i=0 ; i<(2*HP16k_MEL_USED*NB_LP_BANDS_CODING) ; i++)
                bufferCodeForBands16k [i] = bufferCodeForBands16k [i + HP16k_MEL_USED * NB_LP_BANDS_CODING];

            for (i=0 ; i<(HP16k_MEL_USED*NB_LP_BANDS_CODING) ; i++)
                bufferCodeForBands16k [(2 * HP16k_MEL_USED * NB_LP_BANDS_CODING) + i] = CodeForBands16k [i];

            for (i=0 ; i<(2*NB_LP_BANDS_CODING) ; i++)
                bufferCodeWeights [i] = bufferCodeWeights [i + NB_LP_BANDS_CODING];

            for (i=0 ; i<NB_LP_BANDS_CODING ; i++)
                bufferCodeWeights [(2 * NB_LP_BANDS_CODING) + i] = codeWeights [i];
        }
        return TRUE;
    }
    else
        return FALSE;
}
