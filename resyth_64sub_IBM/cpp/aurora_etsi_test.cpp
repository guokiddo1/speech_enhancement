
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "..\aurora_etsi\NoiseSupExports.h"

#include <fstream>
#include <iomanip>
#include "asdk\Wave.h"
#include "..\wav_extract\HuWang.h"
#include "..\main_part\struct.h"
using namespace std;
using namespace asdk;

void extract_Wiener( char *name_id, FILESNAME opts)
{
	//全局及线程初始化
	INSTANCE sm_glb_ins = NULL;
	etsi_denoise_mapping_global_init(&sm_glb_ins, NULL);
	INSTANCE NSX = NULL;
	etsi_denoise_mapping_thread_init(&NSX, sm_glb_ins);

	short *p_data;
	char list[255],inputname[255],outputname[255],outputWname[255],outputMaskname[255],temp[1024];
	
	CWave::Channel channel=CWave::Channel::LEFT_CHANNEL;
	bool bHead = true;
	CWave::Type type = CWave::Type::PCM;
	int nFs = 16000;
	int nChannel = 1;
	int nMaxSample = -1;
	int len = -1;
	
	
	int count = 0;
	int numM = 0;
	FILE *fp_Wiener =NULL;
	int chan = 0;
	if(chan ==0 )
	{
		CWave currentwav;
		//read wav
		sprintf(temp,"%s%s%s_noisy_%d.wav",opts.outputDictionary,opts.save_subband_noisy_wav_dir,name_id,chan);

		currentwav.Read(temp, channel, bHead, type, nFs, nChannel, nMaxSample);
		p_data = currentwav.GetDataPtr();
		len = currentwav.GetSampleNum();
		//wipe the tail
		numM = len;
		numM -= len%BUFFSIZE;
		float *dateIn = (float*)malloc(sizeof(float)*numM);

		for (int i=0; i<numM; i++)
		{
			dateIn[i] = (float)p_data[i];
		}

		//文件获取输入数据结束，开始准备输入输出结构体
		int buffNum = numM/BUFFSIZE;
		esti_denoise_in ins;
		ins.inData = dateIn;
		ins.dataNum = numM;

		esti_denoise_out deOut;
		deOut.outData = (float*)calloc(numM, sizeof(float)); //要求输出结果提前初始化为零
		deOut.pSpeechFoundMel = (int*)malloc(sizeof(int)*buffNum);
		deOut.pSpeechFoundSpec = (int*)malloc(sizeof(int)*buffNum);
		deOut.pSpeechFoundVADNS = (int*)malloc(sizeof(int)*buffNum);
		deOut.pSpeechFoundVar = (int*)malloc(sizeof(int)*buffNum);
		deOut.pFrameCounter = (int*)malloc(sizeof(int)*buffNum);

					
		//降噪处理
		//etsi_denoise_mapping_func(sm_glb_ins, NSX, (INSTANCE)&ins, (INSTANCE)&deOut);
		sprintf(temp, "%s%s%s_noisy_%d.wiener",opts.outputDictionary, opts.save_subband_noisy_Wiener, name_id, chan);
		fp_Wiener = fopen(temp,"w");
		etsi_denoise_mapping_func_Wiener(sm_glb_ins, NSX, (INSTANCE)&ins, (INSTANCE)&deOut, (INSTANCE) fp_Wiener);//ADD_BY_GC_IBM
		fclose(fp_Wiener);

		
		delete []dateIn;
		delete []deOut.outData;
		delete []deOut.pFrameCounter;
		delete []deOut.pSpeechFoundMel;
		delete []deOut.pSpeechFoundSpec;
		delete []deOut.pSpeechFoundVADNS;
		delete []deOut.pSpeechFoundVar;
	}
	for(int chan=0 ; chan < NUM_CHANNEL ; chan++)
	{
		CWave currentwav;
		//read wav
		sprintf(temp,"%s%s%s_noisy_%d.wav",opts.outputDictionary,opts.save_subband_noisy_wav_dir,name_id,chan);

		currentwav.Read(temp, channel, bHead, type, nFs, nChannel, nMaxSample);
		p_data = currentwav.GetDataPtr();
		len = currentwav.GetSampleNum();
		//wipe the tail
		numM = len;
		numM -= len%BUFFSIZE;
		float *dateIn = (float*)malloc(sizeof(float)*numM);

		for (int i=0; i<numM; i++)
		{
			dateIn[i] = (float)p_data[i];
		}

		//文件获取输入数据结束，开始准备输入输出结构体
		int buffNum = numM/BUFFSIZE;
		esti_denoise_in ins;
		ins.inData = dateIn;
		ins.dataNum = numM;

		esti_denoise_out deOut;
		deOut.outData = (float*)calloc(numM, sizeof(float)); //要求输出结果提前初始化为零
		deOut.pSpeechFoundMel = (int*)malloc(sizeof(int)*buffNum);
		deOut.pSpeechFoundSpec = (int*)malloc(sizeof(int)*buffNum);
		deOut.pSpeechFoundVADNS = (int*)malloc(sizeof(int)*buffNum);
		deOut.pSpeechFoundVar = (int*)malloc(sizeof(int)*buffNum);
		deOut.pFrameCounter = (int*)malloc(sizeof(int)*buffNum);

					
		//降噪处理
		//etsi_denoise_mapping_func(sm_glb_ins, NSX, (INSTANCE)&ins, (INSTANCE)&deOut);
		sprintf(temp, "%s%s%s_noisy_%d.wiener",opts.outputDictionary, opts.save_subband_noisy_Wiener, name_id, chan);
		fp_Wiener = fopen(temp,"w");
		etsi_denoise_mapping_func_Wiener(sm_glb_ins, NSX, (INSTANCE)&ins, (INSTANCE)&deOut, (INSTANCE) fp_Wiener);//ADD_BY_GC_IBM
		fclose(fp_Wiener);

		
		delete []dateIn;
		delete []deOut.outData;
		delete []deOut.pFrameCounter;
		delete []deOut.pSpeechFoundMel;
		delete []deOut.pSpeechFoundSpec;
		delete []deOut.pSpeechFoundVADNS;
		delete []deOut.pSpeechFoundVar;
	}
	//全局及线程结构体销毁
	etsi_denoise_mapping_thread_release(&NSX);
	etsi_denoise_mapping_global_release(&sm_glb_ins);
}