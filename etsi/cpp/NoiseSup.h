/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: NoiseSup.h
 * PURPOSE: 1) Apply 2-stage Wiener filter on the input frame.
 *          2) Apply DC offset removal on the output of 2-stage
 *             Wiener filter.
 *          3) Calculate parameters for the frame dropping VAD (see
 *             SpeechQSpec(), SpeechQMel(), SpeechQVar()).
 *
 *-------------------------------------------------------------------------------*/
#ifndef _NOISE_SUP_H
#define _NOISE_SUP_H

/*============================================================================
 *                           INCLUDE PUBLIC DEFINITION
 *============================================================================*/
#include "x_default.h"

/*----------------------------------------------------------------------------
 *                             FILTER CALCULATION
 *----------------------------------------------------------------------------*/
#define NS_FILTER_LENGTH              (X_INT16)(17)
#define NS_HALF_FILTER_LENGTH         (X_INT16)(8)
#define NS_NB_FRAME_THRESHOLD_NSE     (X_INT16)(100)

#define NS_EPS                        (X_FLOAT32)(exp ((double) -10.0))
#define NS_BETA                       (X_FLOAT32)(0.98)
#define NS_RSB_MIN                    (X_FLOAT32)(0.079432823)
#define NS_LAMBDA_NSE                 (X_FLOAT32)(0.99)
#define NS_SPEC_FLOOR                 (X_FLOAT32)(exp ((double) -10.0))
#define NS_LOG_SPEC_FLOOR             (X_FLOAT32)(-10.0)

/*----------------------------------------------------------------------------
 *                                     FFT
 *----------------------------------------------------------------------------*/
#define NS_SPEC_ORDER                 (X_INT16)(65)
#define NS_FFT_ORDER                  (X_INT16)(8)
#define NS_FFT_LENGTH                 (X_INT16)(256)
#define NS_FREQUENCY_BINS             (X_INT16)(129)

/*----------------------------------------------------------------------------
 *                                  BUFFERING
 *----------------------------------------------------------------------------*/
#define NS_PRV_FRAME                  (X_INT16)(0)
#define NS_CUR_FRAME                  (X_INT16)(80)
#define NS_FRAME_SHIFT                (X_INT16)(80)
#define NS_BUFFER_SIZE                (X_INT16)(320)
#define NS_FRAME_LENGTH               (X_INT16)(200)
#define NS_DATA_IN_BUFFER             (X_INT16)(240)
#define NS_SCRATCH_MEM_SIZE           (NS_FFT_LENGTH)
#define NS_NB_FRAMES_LATENCY          (X_INT16)(2)
#define NS_ANALYSIS_WINDOW_8K         (X_INT16)(60)
#define NS_ANALYSIS_WINDOW_16K        (X_INT16)(80)
#define NS_NB_FRAMES_IN_BUFFER        (X_INT16)(4)

/*----------------------------------------------------------------------------
 *                                  PSD MEAN
 *----------------------------------------------------------------------------*/
#define NS_PSD_MEAN_ORDER             (X_INT16)(2)

/*----------------------------------------------------------------------------
 *                                     VAD
 *----------------------------------------------------------------------------*/
#define NS_HANGOVER                   (X_INT16)(15)
#define NS_MIN_FRAME                  (X_INT16)(10)
#define NS_SNR_THRESHOLD_VAD          (X_INT16)(15)
#define NS_SNR_THRESHOLD_UPD_LTE      (X_INT16)(20)
#define NS_NB_FRAME_THRESHOLD_LTE     (X_INT16)(10)
#define NS_MIN_SPEECH_FRAME_HANGOVER  (X_INT16)(4)

#define NS_ENERGY_FLOOR               (X_FLOAT32)(80.0)
#define NS_LAMBDA_LTE_LOWER_E         (X_FLOAT32)(0.97)
#define NS_LAMBDA_LTE_HIGHER_E        (X_FLOAT32)(0.99)

#endif				/* _NOISE_SUP_H */

/*============================================================================
 *                                    END
 *============================================================================*/
