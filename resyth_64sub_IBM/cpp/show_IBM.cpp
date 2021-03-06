
#include "make_IBM.h"
using namespace asdk;


void makeIBM(char *name_id, FILESNAME opts)
{
	//wav cfg
	CWave::Channel channel=CWave::Channel::LEFT_CHANNEL;
	bool bHead = true;
	CWave::Type type = CWave::Type::PCM;
	int nFs = 16000;
	int nChannel = 1;
	int nMaxSample = -1;
	int plen,nlen = -1;

	//for CHANNEL
	char temp[FILE_LEN];
	FILE *fp_IBM_W = NULL;
	FILE *fp_IRM_W = NULL;
	for(int subband = 0 ; subband < NUM_CHANNEL ; subband++ )
	{
		//load chan pure
		CWave pure;
		CWave noise;
		sprintf(temp,"%s%s%s_%d.wav",opts.outputDictionary,opts.save_subband_pure_wav_dir,name_id,subband);
		pure.Read(temp, channel, true, type, nFs, nChannel, nMaxSample);
		//load chan noise
		sprintf(temp,"%s%s%s_noise_%d.wav",opts.outputDictionary,opts.save_subband_noise_wav_dir,name_id,subband);
		noise.Read(temp, channel, true, type, nFs, nChannel, nMaxSample);
		
		int minlen = pure.GetSampleNum();
		int buffNum = (minlen - WINDOW ) / OFFSET + 1;

		mask *Grp_IBM;
		mask *Grp_IRM;
		Grp_IBM = new mask[buffNum];
		Grp_IRM = new mask[buffNum];
		//IBM
		
		SpecInfo specpure(WINDOW, OFFSET);
		specpure.m_nFFT = 256;
		specpure.GetSpecInfo(pure);

		Gamma_Window *FirstWindow = CGammaAlloc ();//ADD_BY_GC_Gamma
		InitGammawindows (FirstWindow, 80.0, (float)SAMPLING_FREQUENCY, 4*(65-1), (int)NUMBER_CHANNEL, 1);//ADD_BY_GC_Gamma
		for(int i = 0 ; i < buffNum ; i++ )
		{
			DoGamma( specpure.m_pPowerSpec[i], FirstWindow);
		}

		SpecInfo specnoise(WINDOW, OFFSET);
		specnoise.m_nFFT = 256;
		specnoise.GetSpecInfo(noise);

		for(int i = 0 ; i < buffNum ; i++ )
		{
			DoGamma( specnoise.m_pPowerSpec[i], FirstWindow);
		}		
		//calcu IBM & iRM
		for(int i = 0; i < buffNum ; i++)
		{
			for(int j = 0 ; j < NUMBER_CHANNEL ; j++ )
			{
				Grp_IBM[i].mark[j] = (specnoise.m_pPowerSpec[i][j] < specpure.m_pPowerSpec[i][j])?1:0;
				Grp_IRM[i].mark[j] =  specpure.m_pPowerSpec[i][j]/(specnoise.m_pPowerSpec[i][j] + specpure.m_pPowerSpec[i][j] );
			}
		}
		//Write
		sprintf(temp,"%s%s%s_%d.IBM",opts.outputDictionary, opts.save_subband_noisy_IBM_dir,name_id,subband);
		fp_IBM_W = fopen(temp, "w");
		sprintf(temp,"%s%s%s_%d.IRM",opts.outputDictionary, opts.save_subband_noisy_IRM_dir,name_id,subband);
		fp_IRM_W = fopen(temp, "w");
		
		for (int frame=0; frame<buffNum; frame++)
		{
			for(int chan=0; chan<NUMBER_CHANNEL; chan++)
			{
				fprintf(fp_IBM_W, "%f ", Grp_IBM[frame].mark[chan]);
				fprintf(fp_IRM_W, "%f ", Grp_IRM[frame].mark[chan]);
			}
			fprintf(fp_IBM_W,"\n");
			fprintf(fp_IRM_W,"\n");
		}
		fclose(fp_IBM_W);
		fclose(fp_IRM_W);
		delete []Grp_IBM;
		delete []Grp_IRM;
	}
	
}
void make_single_IBM(char *name_id, FILESNAME opts)
{
	//wav cfg
	CWave::Channel channel=CWave::Channel::LEFT_CHANNEL;
	bool bHead = true;
	CWave::Type type = CWave::Type::PCM;
	int nFs = 16000;
	int nChannel = 1;
	int nMaxSample = -1;
	int plen,nlen = -1;

	//for CHANNEL
	char temp[FILE_LEN];
	FILE *fp_IBM_W = NULL;
	FILE *fp_IRM_W = NULL;
	for(int subband = 0 ; subband < NUM_CHANNEL ; subband++ )
	{
		//load chan pure
		CWave pure;
		CWave noise;
		sprintf(temp,"%s%s%s_%d.wav",opts.outputDictionary,opts.save_subband_pure_wav_dir,name_id,subband);
		pure.Read(temp, channel, true, type, nFs, nChannel, nMaxSample);
		//load chan noise
		sprintf(temp,"%s%s%s_noise_%d.wav",opts.outputDictionary,opts.save_subband_noise_wav_dir,name_id,subband);
		noise.Read(temp, channel, true, type, nFs, nChannel, nMaxSample);
		
		int minlen = pure.GetSampleNum();
		int buffNum = (minlen - WINDOW ) / OFFSET + 1;

		int *Grp_IBM;
		int *Grp_IRM;
		Grp_IBM = new int[buffNum];
		//IBM
		
		SpecInfo specpure(WINDOW, OFFSET);
		specpure.m_nFFT = 512;
		specpure.GetSpecInfo(pure);


		SpecInfo specnoise(WINDOW, OFFSET);
		specnoise.m_nFFT = 512;
		specnoise.GetSpecInfo(noise);

		float sum_pure = 0;
		float sum_noise = 0;
		//calcu IBM & iRM
		for(int i = 0; i < buffNum ; i++)
		{
			sum_noise = 0;
			sum_pure = 0;
			for(int j = 0 ; j < NUMBER_CHANNEL ; j++ )
			{
				sum_noise += specnoise.m_pPowerSpec[i][j];
				sum_pure += specpure.m_pPowerSpec[i][j];
			}
			Grp_IBM[i] = (sum_noise < sum_pure)?1:0;
		}
		//Write
		sprintf(temp,"%s%s%s_%d.sIBM",opts.outputDictionary, opts.save_subband_noisy_single_IBM_dir,name_id,subband);
		fp_IBM_W = fopen(temp, "w");
	
		fprintf(fp_IBM_W,"%s_noisy_%d ",name_id, subband);
		for (int frame=0; frame<buffNum; frame++)
		{
			fprintf(fp_IBM_W, "%d ", Grp_IBM[frame]);
			
		}
		fprintf(fp_IBM_W,"\n");
		fclose(fp_IBM_W);
		delete []Grp_IBM;
	}
	
}