#ifndef _STRUCT_H_
#define _STRUCT_H_


#define FILE_LEN 255
#define NUM_CHANNEL 25
	
typedef struct _FILESNAME
{
	char func[FILE_LEN];
	char purewavDictionary[FILE_LEN];
	char purewavlist[FILE_LEN];
	char noisepath1[FILE_LEN];
	char noisepath2[FILE_LEN];
	int  addnoisedB;
	//Output wav
	char outputDictionary[FILE_LEN];
	char save_noisy_dir[FILE_LEN];
	char save_subband_pure_wav_dir[FILE_LEN];
	char save_subband_noise_wav_dir[FILE_LEN];
	char save_subband_noisy_wav_dir[FILE_LEN];
	//Mask
	char save_subband_noisy_IBM_dir[FILE_LEN];
	char save_subband_noisy_IRM_dir[FILE_LEN];
	char save_subband_noisy_single_IBM_dir[FILE_LEN];
	//feature
	char save_subband_noisy_MFCC[FILE_LEN];
	char save_subband_noisy_ACF[FILE_LEN];
	char save_subband_noisy_Wiener[FILE_LEN];
	//log
	char Log[FILE_LEN];

}FILESNAME;


#endif