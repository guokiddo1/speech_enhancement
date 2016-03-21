#ifndef _STRUCT_H_
#define _STRUCT_H_


#define FILE_LEN 255
#define NUM_CHANNEL 64
	
typedef struct _FILESNAME
{
	char purewavDictionary[FILE_LEN];
	char purewavlist[FILE_LEN];
	int numMix;
	
	//Output wav
	char outputDictionary[FILE_LEN];
	char save_noisy_dir[FILE_LEN];
	char save_noisy_ebm_dir[FILE_LEN];
	char save_noisy_sirm_dir[FILE_LEN];
	char save_resynth_e_dir[FILE_LEN];
	char save_resynth_i_dir[FILE_LEN];
	//log
	char Log[FILE_LEN];

}FILESNAME;


#endif
