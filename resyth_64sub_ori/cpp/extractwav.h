#include "struct.h"
#include "Wave.h"
#include "GMM.h"
#include "fft_core.h"

void resynth(asdk::CWave *wav,FILESNAME opts, const vector<vector<float> > &v_erm, char* name_id);
float Train_thd(int numMix, const vector<vector<float> > &v_erm, int numFrame, int chan);

void fft(float *data_each_frame,  int len, float *spec_mag, complex<double> *spec_angle);
void spec_minus(float *spec_mag, float *noise_spec_mag, int len);
void ifft(float *spec_mag, complex<double> *spec_angle, int len, float *denoise_data);
