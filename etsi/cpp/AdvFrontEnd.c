/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: AdvFrontEnd.c
 * PURPOSE:   This file contains the main part of Advanced DSR front-end.
 *            Speech samples are read from input waveform file frame by frame.
 *            Feature extraction is performed for each frame by calling the
 *            function DoAdvProcess () (see ParmInterface.c).
 *            Feature vectors are output to file in HTK format. VAD information
 *            is output to file in ASCII format.
 *            Command line arguments are handled by a command line parsing
 *            function.
 *
 *-------------------------------------------------------------------------------*/

/*-----------------
 * File Inclusions
 *-----------------*/
#include "AdvFrontEnd.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fileio.h"
#include "ParmInterface.h"

/*-------------------------------
 * Global Definitions and Macros
 *-------------------------------*/
#define PRINTMOD 15
#define IEEE_LE 0
#define IEEE_BE 1

/*---------------
 * HTK constants
 *---------------*/
#define HTK_PR_MFCC 6
#define HTK_PR_E    0x40
#define HTK_PR_0    0x2000

/*----------------------------------------------
 * Global Variable Definitions and Declarations
 *----------------------------------------------*/
BOOLEAN QuietMode = FALSE,	/* Supress output to stderr */
    FsSpecified = FALSE,	/* Sampling frequency specified */
    SwapSpecified = FALSE,	/* Byte swap for raw data files specified */
    InputKindSpecified = FALSE,	/* Input file format specified */
    NoOutHeaderSpecified = FALSE,	/* No output HTK header option specified */
    Noc0 = FALSE,		/* No c0 coefficient to output feature vector */
    NologE = FALSE;		/* No logE component to output feature vector */

/* FILE *fp_in = NULL,		/\* Input HTK, NIST or raw data file *\/ */
/*      *fp_denoised=NULL;            /\* Added by myself *\/ */


char InFilename[199],		/* Name of input file */
    DenoisedFilename[199],        /* Added by myself  */  
    InputKind[10] = "RAW";	 /* Input file format */

int SamplingFrequency = 8000,	/* SamplingFrequency */
    NativeByteOrder = IEEE_LE,	/* Native byte ordering */
    InputByteOrder,		/* Default input byte ordering */
    OutParmKind = HTK_PR_MFCC + HTK_PR_E + HTK_PR_0;	/* Output parameter kind MFCC_0_E as default */

#ifndef NS_FRAME_SHIFT
#define NS_FRAME_SHIFT                (X_INT16)(80)
#endif

/*------------
 * Prototypes
 *------------*/
//static int ReadBufWave(FILE * fp_in, short * buf, int nSamples, int Swap);
static int ReadBufWave(short * p_data, long n_start, long n_frame,
                       short * buf, int nSamples, int Swap);

//static int 
//ReadBufWave(FILE * fp_in, short * buf, int nSamples, int Swap)
/* static int */
/* ReadBufWave(FILE * fp_in, short * buf, int nSamples, int Swap) */
static int ReadBufWave(short * p_data, long n_start, long n_frame,
                       short * buf, int nSamples, int Swap)
{
    int i;

    for (i = 0; i < nSamples; i++)
    {
        /* if (fread(&s, sizeof(s), 1, fp_in) != 1) */
        /*     return FALSE; */
        if( n_start + i >= n_frame )
            return FALSE;
        
        /* if (Swap) */
        /*     s = ((s & 0x00ff) << 8) | ((s & 0xff00) >> 8); */

        buf[i] = (short) p_data[ n_start + i ];
    }
    return TRUE;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: main
 *
 * PURPOSE:       Main front-end operations from input speech samples to output
 *                feature vectors. See embedded comments.
 *
 * INPUT:
 *   argc         Number of command line arguments (passed to ParseCommLine)
 *   argv         Array of command line arguments (passed to ParseCommLine)
 *
 * OUTPUT
 *   none
 *
 * RETURN VALUE
 *   TRUE         In case of any errors
 *   FALSE        Otherwise
 *
 *---------------------------------------------------------------------------*/
int etsi_denoise( short* p_data, short* p_denoised, long n_frame )
{
    short *SigBuf;
    X_FLOAT32 FeatureBuffer[NUM_CEP_COEFF + 2];
    FEParamsX *pFEParX;

    int NbSamplesToRead;

    long FrameCounter = 0;
    long SpeechFrameCounter = 0;
    long NonSpeechFrameCounter = 0;
	short *DenoiseBuffer;
	long n_start = 0;
    
    /*----------------*/
    /* Initialization */
    /*----------------*/
    InputByteOrder = NativeByteOrder;

    QuietMode = TRUE;
    SamplingFrequency = 8000;
    SwapSpecified = FALSE;

    /*------------------------------------------
     * Memory allocation for FE data structures
     *------------------------------------------*/
    pFEParX = AdvProcessAlloc(SamplingFrequency);

    /*-------------------------------------------------------
     * Initialization of FE data structures and input buffer
     *-------------------------------------------------------*/
    pFEParX->Noc0 = Noc0;
    AdvProcessInit(pFEParX);

    NbSamplesToRead = pFEParX->NbSamplesToRead;

    DenoiseBuffer = (short*) calloc(1, sizeof( short) * NS_FRAME_SHIFT );
             
    SigBuf = (short int*) calloc(1, sizeof(SigBuf[0]) * NbSamplesToRead) ;
    if (SigBuf == NULL)
    {
        fprintf(stderr, "ERROR:   Memory allocation error occured!\r\n");
        goto _FaultExit;
    }

	/*------------
     * Processing
     *------------*/
    while (ReadBufWave( p_data, n_start, n_frame, SigBuf, NbSamplesToRead,
                       (SwapSpecified || InputByteOrder != NativeByteOrder)))
    {
		int i = 0;
        FrameCounter++;

        if (DoAdvProcess(SigBuf, DenoiseBuffer, FeatureBuffer, pFEParX))
        {
            //if (pFEParX->VAD == SPEECH_FRAME)
            //    SpeechFrameCounter++;
            //else
            //    NonSpeechFrameCounter++;
        }
        for (i=0;i<NbSamplesToRead;i++)
        {
            p_denoised[n_start+i] = DenoiseBuffer[i];
        }
        n_start += NbSamplesToRead;
    }
        
    /*----------------
     * Memory release
     *----------------*/
    AdvProcessDelete(&pFEParX);
    free(SigBuf);

	/*----------------------
     * Display final status
     *----------------------*/
    if (SpeechFrameCounter == 0)
        fprintf(stderr, "NO SPEECH DETECTED !\r\n");
    
    return FALSE;

_FaultExit: 

    return TRUE;
}


int etsi_denoise_synchronization( short* p_data, short* p_denoised, long i_frame )
{
	BOOLEAN Flag;
	short *tmp;
	Flag = TRUE;
	tmp = (short*) malloc(i_frame*sizeof(short));
	Flag = etsi_denoise( p_data, tmp, i_frame );
	if (Flag){
		memset(p_denoised,0,i_frame*sizeof(short));
		memcpy(p_denoised,(tmp+320),i_frame - 320);
	}
	free(tmp);
	return Flag;
}

int etsi_denoise_16k( short* p_data, short* p_denoised, long n_frame )
{
    short *SigBuf;
    X_FLOAT32 FeatureBuffer[NUM_CEP_COEFF + 2];
    FEParamsX *pFEParX;

    int NbSamplesToRead;

    long FrameCounter = 0;
    long SpeechFrameCounter = 0;
    long NonSpeechFrameCounter = 0;
	short *DenoiseBuffer;
	long n_start = 0;
    
    /*----------------*/
    /* Initialization */
    /*----------------*/
    InputByteOrder = NativeByteOrder;

    QuietMode = TRUE;
    SamplingFrequency = 16000;
    SwapSpecified = FALSE;

    /*------------------------------------------
     * Memory allocation for FE data structures
     *------------------------------------------*/
    pFEParX = AdvProcessAlloc(SamplingFrequency);

    /*-------------------------------------------------------
     * Initialization of FE data structures and input buffer
     *-------------------------------------------------------*/
    pFEParX->Noc0 = Noc0;
    AdvProcessInit(pFEParX);

    NbSamplesToRead = pFEParX->NbSamplesToRead;

    DenoiseBuffer = (short*) calloc(1, sizeof( short) * NS_FRAME_SHIFT );
             
    SigBuf = (short int*) calloc(1, sizeof(SigBuf[0]) * NbSamplesToRead) ;
    if (SigBuf == NULL)
    {
        fprintf(stderr, "ERROR:   Memory allocation error occured!\r\n");
        goto _FaultExit;
    }

	/*------------
     * Processing
     *------------*/
    while (ReadBufWave( p_data, n_start, n_frame, SigBuf, NbSamplesToRead,
                       (SwapSpecified || InputByteOrder != NativeByteOrder)))
    {
		int i = 0;
        FrameCounter++;

        if (DoAdvProcess(SigBuf, DenoiseBuffer, FeatureBuffer, pFEParX))
        {
            //if (pFEParX->VAD == SPEECH_FRAME)
              //  SpeechFrameCounter++;
            //else
            //    NonSpeechFrameCounter++;
        }
        for (i=0;i<NbSamplesToRead;i++)
        {
            p_denoised[n_start+i] = DenoiseBuffer[i];
        }
        n_start += NbSamplesToRead;
    }
        
    /*----------------
     * Memory release
     *----------------*/
    AdvProcessDelete(&pFEParX);
    free(SigBuf);

	/*----------------------
     * Display final status
     *----------------------*/
    if (SpeechFrameCounter == 0)
        fprintf(stderr, "NO SPEECH DETECTED !\r\n");
    
    return FALSE;

_FaultExit: 

    return TRUE;
}


int etsi_denoise_16k_synchronization( short* p_data, short* p_denoised, long i_frame )
{
	BOOLEAN Flag;
	short *tmp;
	Flag = TRUE;
	tmp = (short*) malloc(i_frame*sizeof(short));
	Flag = etsi_denoise( p_data, tmp, i_frame );
	if (Flag){
		memset(p_denoised,0,i_frame*sizeof(short));
		memcpy(p_denoised,(tmp+320),i_frame - 320);
	}
	free(tmp);
	return Flag;
}