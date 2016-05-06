#include "asdk\Wave.h"
#include "asdk\MFCC.h"
#include "asdk\SpecInfo.h"
#include "asdk\PLP.h"
#include "iostream"
#include "..\main_part\struct.h"
#include "..\wav_extract\HuWang.h"
using namespace std;
using namespace asdk;

void extract_MFCC(char *name_id, FILESNAME opts)
{
	char temp[FILE_LEN];
	int len=0;
	//wav cfg
	CWave::Channel channel=CWave::Channel::LEFT_CHANNEL;
	bool bHead = true;
	CWave::Type type = CWave::Type::PCM;
	int nFs = 16000;
	int nChannel = 1;
	int nMaxSample = -1;
	int plen,nlen = -1;

	//MFCC cfg
	CMFCC::Parm parm1;
    parm1.bHTK = true;
    parm1.nFrameLen = WINDOW;
    parm1.nFrameStep = OFFSET;
    parm1.fLowFreq = 64.0f;
    parm1.fHighFreq = 4000.0f;
    parm1.eE = CMFCC::Parm::RawE;
    parm1.bPower = false;
    parm1.bENorm = true;
    parm1.fEscale = 1.0f;
    parm1.fSilFloor = 50.0f;
    parm1.nFilter = 26;
    parm1.fCepLifter = 22;
    parm1.nCep = 12;
    parm1.bCMS = false;
    parm1.bC0 = false;

	CFeature fea(13, true);
    CMFCC mfcc(parm1);
	

	SpecInfo specpure(WINDOW, OFFSET);
	specpure.m_nFFT = 256;
	CWave noisywav;
	for(int subband =0 ; subband < NUM_CHANNEL; subband++)
	{
		sprintf(temp,"%s%s%s_noisy_%d.wav",opts.outputDictionary ,opts.save_subband_noisy_wav_dir,name_id,subband);
		noisywav.Read(temp, channel, true, type, nFs, nChannel, nMaxSample);
		len = noisywav.GetSampleNum();
		len  -= len%OFFSET;
		mfcc.Extract(noisywav.GetDataPtr(), len, nFs, fea);
		sprintf(temp,"%s%s%s_noisy_%2d.mfcc",opts.outputDictionary, opts.save_subband_noisy_MFCC,name_id,subband);
		fea.WriteHTK(temp);
	}
}
void extract_FBank()
{

}