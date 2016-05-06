#include <stdlib.h>
#include <iostream>
extern "C"
{
#include "AdvFrontEnd.h"
}
#include "Wave.h"
#include "struct.h"

using namespace std;
using namespace asdk;

void Read_CFG(char *Cfg_File, FILESNAME &opts);
void chop(char *string);


int main( int argc, char *argv[]){
	int run = 0;
	long n_frame = 0;
	short *p_data, *p_denoised;
	char cfg[FILE_LEN], name_id[FILE_LEN], temp[FILE_LEN];
	
	CWave::Channel channel=CWave::LEFT_CHANNEL;
	bool bHead = true;
	CWave::Type type = CWave::PCM;
	int nFs = 16000;
	int nChannel = 1;
	int nMaxSample = -1;
	int len = -1;
	
	//program cfg
	FILESNAME opts;
	strcpy(cfg, argv[1]);
	//read cfg
	Read_CFG(cfg, opts);
	FILE *Log = NULL;
	sprintf(temp, "%s%s", opts.outputDictionary, opts.Log);
	Log = fopen(temp, "a+");
	FILE *fp_purewav = NULL;
	fp_purewav = fopen(opts.purewavlist, "r");

	int count = 0;
	while (fgets(temp, 1024, fp_purewav)){
		strcpy(name_id, temp);
		chop(name_id);
		cout << name_id << endl;
		fprintf(Log, "%s\n ", name_id);
		sprintf(temp, "%s%s%s_noisy.wav", opts.outputDictionary, opts.save_noisy_dir, name_id);

		cout<<temp<<" "<<count<<endl;
		count +=1;

		CWave currentwav;
		currentwav.Read(temp, channel, bHead, type, nFs, nChannel, nMaxSample);
		p_data = currentwav.GetDataPtr();
		len = currentwav.GetSampleNum();
		n_frame = (long)len;
		p_denoised = new short[len];
		run = etsi_denoise( p_data, p_denoised, n_frame );

		CWave denoisedwav(p_denoised, len, nFs);

		sprintf(temp, "%s%s%s_e_resynth.wav", opts.outputDictionary, opts.save_resynth_e_dir, name_id);
		denoisedwav.Write(temp);

		delete [] p_denoised;
	}
	fclose(fp_purewav);
	fclose(Log);

	return 0;
}
void Read_CFG(char *Cfg_File, FILESNAME &opts)
{
	char Char_temp[1024];
	char tmp[FILE_LEN], tmp0[FILE_LEN];

	FILE *fp = NULL;
	fp = fopen(Cfg_File, "r");
	if (fp == NULL)
	{
		printf("Open %s file error!\n", Cfg_File);
		return;
	}
	/*
	**********************INPUT**************************
	*/
	//purewavDictionary=
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.purewavDictionary);
	//purewavlist=
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.purewavlist);
	//numMix=
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %d", tmp0, &opts.numMix);


	/*
	**********************OUTPUT**************************
	*/
	//outputDictionary
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.outputDictionary);

	//save_noisy_dir
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.save_noisy_dir);

	//save_noisy_ebm_dir
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.save_noisy_ebm_dir);

	//save_noisy_sirm_dir
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.save_noisy_sirm_dir);

	//save_resynth_e_dir
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.save_resynth_e_dir);

	//save_resynth_i_dir
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.save_resynth_i_dir);

	//log
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.Log);
	fclose(fp);
	return;
}


void chop(char *string)
{
	int len = strlen(string);
	string[len - 1] = '\0';
}
