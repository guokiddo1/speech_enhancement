
#ifndef _NOISESUPEXPORTS_H
#define _NOISESUPEXPORTS_H

#include "aurora/aurora_include.h"
using namespace aurora;
#define FILE_TYPE   X_INT16
const int BUFFSIZE = 160;//80
struct DENOISEGlobalImpl
{
	int SamplingFrequency;
};

struct esti_denoise_in
{
	X_FLOAT32 *inData;
	int dataNum;
};

struct esti_denoise_out
{
	X_FLOAT32 *outData;
	int *pSpeechFoundVar;
	int *pSpeechFoundSpec;
	int *pSpeechFoundMel;
	int *pSpeechFoundVADNS;
	int *pFrameCounter;
};
#define NUMBER_CHANNEL  25      		        /* maxmimum number of filters ADD_BY_GC_IBM*/
struct mask
{
	float mark[NUMBER_CHANNEL];
};

int32s etsi_denoise_mapping_global_init(PINSTANCE sm_glb_pins, INSTANCE sm_glb_res);
int32s etsi_denoise_mapping_thread_init(PINSTANCE sm_thd_pins, INSTANCE sm_glb_ins);

int32s etsi_denoise_mapping_func_Wiener(INSTANCE sm_glb_ins, INSTANCE sm_thd_ins, INSTANCE in_ins, INSTANCE out_ins, INSTANCE fp_IBM);
int32s etsi_denoise_mapping_func(INSTANCE sm_glb_ins, INSTANCE sm_thd_ins, INSTANCE in_ins, INSTANCE out_ins);

void etsi_denoise_mapping_thread_release(PINSTANCE sm_thd_pins);
void etsi_denoise_mapping_global_release(PINSTANCE sm_glb_pins);

#endif
