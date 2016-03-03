#include "..\main_part\struct.h"
#include "asdk\Wave.h"


void addnoise(asdk::CWave *purewav, asdk::CWave *temp_noisewav, asdk::CWave *noisywav,char* name_id, FILESNAME opts);
void subbband(asdk::CWave *wav, char* path, char* addhead, char* name_id);