/*
	**********************************
	*Writer:
		Cong Guo
	*Function:
		1.mix pure and wav to noisy randomly 
		2.splite wav and noisy into subband
		3.calculate each subband's IBM and IRM
		4.extract the feature of each subband
	*Input_cfg:
		purewav_dir,purewav_list
		noisepath
		output_sum_dir
		//Output wav
		save_noisy_dir
		save_subband_pure_wav_dir
		save_subband_noise_wav_dir
		save_subband_noisy_wav_dir
		//Mask
		save_subband_noisy_IBM_dir
		save_subband_noisy_IRM_dir
		//feature
		save_subband_noisy_MFCC
		save_subband_noisy_ACF
		save_subband_noisy_Wiener
	*Time:
		07-May-2015
	**********************************
*/
#include "iostream"
#include "asdk\Wave.h"
#include "stdlib.h"
#include "struct.h"
#include "string.h"
#include "..\wav_extract\HuWang.h"
#include "..\wav_extract\extractwav.h"
#include <ctime>
#include "cmath"
#include <vector>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;
using namespace asdk;

void Read_CFG(char *Cfg_File, FILESNAME &opts);
void chop(char *string);

int main(int argc, char *argv[])
{
	char cfg[FILE_LEN], temp[FILE_LEN], purewav_name[FILE_LEN], name_id[FILE_LEN];

	//program cfg
	FILESNAME opts;
	strcpy(cfg, argv[1]);
	//read cfg
	Read_CFG(cfg,opts);
	FILE *Log =  NULL;
	sprintf(temp,"%s%s",opts.outputDictionary,opts.Log);
	Log = fopen(temp, "a+");
	//wav cfg
	CWave::Channel channel=CWave::Channel::LEFT_CHANNEL;
	bool bHead = true;
	CWave::Type type = CWave::Type::PCM;
	int nFs = 16000;
	int nChannel = 1;
	int nMaxSample = -1;
	int plen,nlen = -1;
	int FLag = 0;
	int Total_number=0;
	//START¡¡ALL
	srand((unsigned)time(0)); 
	FILE *fp_purewav = NULL;
	fp_purewav = fopen(opts.purewavlist, "r");
	CWave noisywav;
	srand((int)time(0));
	int n = 0;
	char temp1[FILE_LEN];
	FILE *erm;
	sprintf(temp1,"%sresult.txt",opts.outputDictionary);
	cout<<temp1<<endl;
	erm = fopen(temp1,"r");
	char buf[64*12];
	vector<vector<float> > v_erm;
	while(fgets(buf,sizeof(buf),erm))
	{
		int numFrame;
		
		if(strstr(buf,"["))
		{
			FLag++;
			cout<<FLag<<endl;
			//cout<<buf<<endl;
			while(fgets(temp, 1024, fp_purewav))
			{
				strcpy(name_id, temp);
				chop(name_id);
				cout<<name_id<<endl;
				fprintf(Log, "%s\n ", name_id);

				/*****************WAV*********************/
		
				sprintf(purewav_name,"%s%s%s_noisy.wav",opts.outputDictionary,opts.save_noisy_dir,name_id);
				fprintf(Log, "%s\n ", name_id);
				noisywav.Read(purewav_name, channel, true, type, nFs, nChannel, nMaxSample);
				int sigLength;
				sigLength = noisywav.GetSampleNum();
				numFrame =  sigLength/OFFSET;
				cout<<numFrame<<endl;
				for(int i=0; i<numFrame; i++)
				{
					v_erm.push_back( vector<float> (NUM_CHANNEL));
				}
					//cout<<"SS"<<v_erm.size()<<endl;
				break;
			}
		}
		else 
		{
			int i,j;
			for(i=0; i<numFrame;i++)
			{ 
				//cout<<buf<<endl;
				string aaa(buf);
				stringstream ss(aaa);
				for(j=0; j<NUM_CHANNEL; j++)
				{
					ss >> v_erm[i][j];
					//cout<<v_erm.size()<<endl;
					//cout<<v_erm[i].size()<<endl;
					cout<<v_erm[i][j]<<endl;
				}
				cout<<v_erm[i][NUM_CHANNEL-1]<<endl;
				//cout<<v_erm[numFrame-1][NUM_CHANNEL-1]<<endl;
				break;
			}
			//cout<<v_erm[numFrame-1][NUM_CHANNEL-1]<<endl;
			while(v_erm[numFrame-1][NUM_CHANNEL-1]>0)
			{
				cout<<v_erm[numFrame-1][NUM_CHANNEL-1]<<endl;
				fprintf(Log, "resynth\n ");
				resynth(&noisywav, opts, v_erm, name_id);
			}
		}
	}
	fclose(fp_purewav);
	fclose(Log);
	return 1;
}
void Read_CFG(char *Cfg_File, FILESNAME &opts)
{
	char Char_temp[1024];
	char tmp[FILE_LEN], tmp0[FILE_LEN];

	FILE *fp=NULL;
	fp=fopen(Cfg_File, "r");
	if(fp==NULL)
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

	//save_noisy_sibm_dir
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
	string[len-1]='\0';
}
