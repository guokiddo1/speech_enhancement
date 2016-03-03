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
#include "..\wav_extract\extractwav.h"
#include "..\fea_extract\fea_extrac.h"
#include "..\show_IBM\make_IBM.h"
#include <ctime>
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
	int nFs = 8000;
	int nChannel = 1;
	int nMaxSample = -1;
	int plen,nlen = -1;
	
	//read noise
	CWave noisewav1;
	noisewav1.Read(opts.noisepath1, channel, true, type, nFs, nChannel, nMaxSample);
	CWave noisewav2;
	noisewav2.Read(opts.noisepath2, channel, true, type, nFs, nChannel, nMaxSample);
	cout << "read noise" << endl;
	//STARTĦĦALL
	srand((unsigned)time(0)); 
	FILE *fp_purewav = NULL;
	cout << "open purewave_list" << endl;
	fp_purewav = fopen(opts.purewavlist, "r");
	CWave purewav;
	srand((int)time(0));
	float rand_key = 0.0;
	int n1,n2 = 0;
	bool test = false;
	cout << "read purewave_list" << endl;
	if (strcmp(opts.func, "test") != 0)
		test = true;
	while(fgets(temp, 1024, fp_purewav)){
		strcpy(name_id, temp);
		chop(name_id);
		cout<<name_id<<endl;
		fprintf(Log, "%s\n ", name_id);
		/*
		****************WAV*********************
		*/
		//read pure
		sprintf(purewav_name,"%s%s.wav",opts.purewavDictionary,name_id);
		fprintf(Log, "%s\n ", name_id);
		purewav.Read(purewav_name, channel, true, type, nFs, nChannel, nMaxSample);
		//temp noise wav
		int  len_noise = purewav.GetSampleNum();
		rand_key =  rand() % 10001 / 10000.0;
		
		n1 = (int)(rand_key * 6);
		cout << n1 << endl;
		short *point_noise;
		short *p_noisy;
		short *p_noise;
		short *p_pure;

		p_noisy = new short[len_noise];
		point_noise = new short[len_noise];
		p_pure = purewav.GetDataPtr();
		switch (n1)
		{
			case 0:{
					   opts.addnoisedB = 5;
					   n2 = ( int )(rand_key * (noisewav1.GetSampleNum() - purewav.GetSampleNum()));
					   p_noise = noisewav1.GetDataPtr();
					   for (int j = 0; j < len_noise; j++)
					   {
						   point_noise[j] = p_noise[n2 + j];
						   p_noisy[j] = p_pure[j];
					   }
					   break;
			}
			case 1:{
					   opts.addnoisedB = 0;
					   n2 = (int)(rand_key * (noisewav1.GetSampleNum() - purewav.GetSampleNum()));
					   p_noise = noisewav1.GetDataPtr();
					   for (int j = 0; j < len_noise; j++)
					   {
						   point_noise[j] = p_noise[n2 + j];
						   p_noisy[j] = p_pure[j];
					   }
					   break;
			}
			case 2:{
					   opts.addnoisedB = -5;
					   n2 = (int)(rand_key * (noisewav1.GetSampleNum() - purewav.GetSampleNum()));
					   p_noise = noisewav1.GetDataPtr();
					   for (int j = 0; j < len_noise; j++)
					   {
						   point_noise[j] = p_noise[n2 + j];
						   p_noisy[j] = p_pure[j];
					   }
					   break;
			}
			case 3:{
					   opts.addnoisedB = 5;
					   n2 = (int)(rand_key * (noisewav2.GetSampleNum() - purewav.GetSampleNum()));
					   p_noise = noisewav2.GetDataPtr();
					   for (int j = 0; j < len_noise; j++)
					   {
						   point_noise[j] = p_noise[n2 + j];
						   p_noisy[j] = p_pure[j];
					   }
					   break;
			}
			case 4:{
					   opts.addnoisedB = 0;
					   n2 = (int)(rand_key * (noisewav2.GetSampleNum() - purewav.GetSampleNum()));
					   p_noise = noisewav2.GetDataPtr();
					   for (int j = 0; j < len_noise; j++)
					   {
						   point_noise[j] = p_noise[n2 + j];
						   p_noisy[j] = p_pure[j];
					   }
					   break;
			}
			case 5:{
					   opts.addnoisedB = -5;
					   n2 = (int)(rand_key * (noisewav2.GetSampleNum() - purewav.GetSampleNum()));
					   p_noise = noisewav2.GetDataPtr();
					   for (int j = 0; j < len_noise; j++)
					   {
						   point_noise[j] = p_noise[n2 + j];
						   p_noisy[j] = p_pure[j];
					   }
					   break;
			}
			default:
				break;
		}
			
		CWave temp_noisewav(point_noise, len_noise, nFs);
		
		//mix pure + noise = noisy

		CWave noisywav(p_noisy, len_noise,nFs);
		//noisywav.Read(purewav_name, channel, false, type, nFs, nChannel, nMaxSample);
		if (test)
			addnoise(&purewav, &temp_noisewav, &noisywav, name_id, opts);
		else
		{
			sprintf(temp, "%s%s%s_noisy.wav", opts.outputDictionary, opts.save_noisy_dir, name_id);
			noisywav.Write(temp);
		}

		/*
		****************SUBBAND*********************
		*/
		cout<<"subband"<<endl;
		fprintf(Log, "subband\n ");
		char path[FILE_LEN],addtail[FILE_LEN];
		//splite pure into subband_pure
		sprintf(path ,"%s%s",opts.outputDictionary,opts.save_subband_pure_wav_dir);
		strcpy(addtail,"");
		subbband(&purewav, path, addtail, name_id);
		//splite noise into subband_noise
		sprintf(path ,"%s%s",opts.outputDictionary,opts.save_subband_noise_wav_dir);
		strcpy(addtail, "_noise");
		subbband(&temp_noisewav, path, addtail, name_id);
	
		//splite noisy into subband_noisy
		sprintf(path, "%s%s", opts.outputDictionary, opts.save_subband_noisy_wav_dir);
		strcpy(addtail, "_noisy");
		subbband(&noisywav, path, addtail, name_id);

		/*
		****************MASK*********************
		*/
		//for each band do IBM IRM
		//fprintf(Log, "IBM & IRM\n ");
		//cout<<"IBM & IRM"<<endl;
		//makeIBM(name_id, opts);
		if (test)
		{
			fprintf(Log, "single_IBM\n ");
			cout << "single_IBM" << endl;
			make_single_IBM(name_id, opts);
			/*
			****************FEATURE*********************
			*/
			//MFCC
			//extract_MFCC(name_id, opts);
			//Wiener-filter
			//cout<<"Wiener fea"<<endl;
			//fprintf(Log, "Wiener fea\n ");
			//extract_Wiener(name_id, opts);

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
	//func= 
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.func);
	//purewavDictionary= 
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.purewavDictionary);
	//purewavlist= 
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.purewavlist);
	
	//noisepath= 
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.noisepath1);
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.noisepath2);
	//addnoisedB= 5
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %d", tmp0, &opts.addnoisedB);
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

	//save_subband_pure_wav_dir
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.save_subband_pure_wav_dir);

	//save_subband_noise_wav_dir
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.save_subband_noise_wav_dir);

	//save_subband_noisy_wav_dir
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.save_subband_noisy_wav_dir);

	/*
	**********************MASK**************************
	*/
	//save_subband_noisy_IBM_dir
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.save_subband_noisy_IBM_dir);

	//save_subband_noisy_IRM_dir
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.save_subband_noisy_IRM_dir);

	//save_subband_noisy_sIBM_dir
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.save_subband_noisy_single_IBM_dir);

	/*
	**********************FEAT**************************
	*/
	//save_subband_noisy_MFCC
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.save_subband_noisy_MFCC);

	//save_subband_noisy_ACF
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.save_subband_noisy_ACF);

	//save_subband_noisy_Wiener
	fgets(Char_temp, 1024, fp);
	chop(Char_temp);
	sscanf(Char_temp, "%s %s", tmp0, opts.save_subband_noisy_Wiener);

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
