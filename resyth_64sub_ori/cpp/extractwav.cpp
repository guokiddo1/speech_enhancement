#include "extractwav.h"
#include "HuWang.h"
#include "cmath"
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;
using namespace asdk;
void resynth(CWave *wav,FILESNAME opts, const vector<vector<float> > &v_erm, char* name_id)
{
	int n=0;
	//wav cfg
	CWave::Channel channel=CWave::LEFT_CHANNEL;
	bool bHead = true;
	CWave::Type type = CWave::PCM;
	int nFs = 16000;
	int nChannel = 1;
	int nMaxSample = -1;
	int plen,nlen = -1;

	gammaTone fChan[NUMBER_CHANNEL];
	float *Input;
	float *gOut[NUMBER_CHANNEL], *hOut[NUMBER_CHANNEL];
	short *s_hOut[NUMBER_CHANNEL];
	float lowerERB, upperERB, spaceERB;
	float cf, phon;
	int chan,sigLength;
	sigLength = wav->GetSampleNum();
	short *p_wav = wav->GetDataPtr();

	Input = new float[sigLength];
	middleEar loudFunc = initMiddleEar();
	for (chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		gOut[chan] = new float[sigLength];
		hOut[chan] = new float[sigLength];
		s_hOut[chan] = new short[sigLength];
	}


	lowerERB = 21.4*log10(MINCF*0.00437 + 1.0);
	upperERB =  21.4*log10(MAXCF*0.00437 + 1.0);
  
	spaceERB = (NUMBER_CHANNEL > 1) ? (upperERB-lowerERB)/(NUMBER_CHANNEL-1) : 0;
	cout<<"gammatone filter init"<<endl;
	for (chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		cf = (pow(10, (lowerERB + chan*spaceERB)/21.4) - 1) / 0.00437 ;
		fChan[chan].cf = cf;
		fChan[chan].bw = 24.7*(cf*0.00437 + 1.0) * BW_CORRECTION;
		
		phon = (loudnessLevelInPhons(cf, loudFunc) - DB);
		fChan[chan].midEarCoeff = pow(10, phon/20); 
	}
	for(int i = 0 ; i < sigLength; i++)
	{
		Input[i] = (float)p_wav[i];
	}
	cout<<"subband to 64 channel"<<endl;
	for(chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		gammaToneFilter(Input, gOut[chan], fChan[chan], sigLength);
		hairCell(gOut[chan], hOut[chan], sigLength);
	}
	char temp[FILE_LEN],temp2[FILE_LEN];
	int frame = 0;
	int numFrame =  (sigLength - WINDOW)/OFFSET +1;

	float *resynth_e,*weight_e,*reverse;//
	resynth_e = new float[sigLength];
	weight_e = new float[sigLength];
	//resynth_i = new float[sigLength];
	//weight_i = new float[sigLength];
	reverse = new float[sigLength];
	for(n=0; n<sigLength; n++)
	{
		resynth_e[n] = 0;
		//resynth_i[n]=0;
	}
	
 	float thd = 0.0;	
	int *noise_frame, *noise_estimate;
	int num_noise_frame = 0;
	noise_frame = new int[numFrame];
	noise_estimate = new int[numFrame];

	float *data_each_frame, *spec_mag, *noise_data_each_frame, *noise_spec_mag;
	complex<double> *spec_angle, *noise_spec_angle;
	float *denoise_data;

	spec_mag = new float[512];
	noise_spec_mag = new float[512];
	spec_angle = new complex<double>[512];
	noise_spec_angle = new complex<double>[512];
	denoise_data = new float[WINDOW];
	for(chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		// Train erm thd for each channal
		thd = Train_thd( opts.numMix, v_erm, numFrame, chan);
		cout<< "thd for " << chan << " channal = "<< thd<<endl;
		
		// noise_frame[] have the noise_frame_idx
		num_noise_frame = 0;
		for(frame = 0 ; frame < numFrame ; frame++)
		{
			if( v_erm[frame][chan] <= thd )
			{
				noise_frame[num_noise_frame] = frame;
				num_noise_frame++;
			}
		}
		//noise estimate for each frame
		int noise_index = 0;
		for(frame = 0 ; frame < numFrame ; frame++)
		{
			if( v_erm[frame][chan] > thd )
			{	if(noise_index == 0)
					noise_estimate[frame] = noise_frame[noise_index];
				else
					noise_estimate[frame] = noise_frame[noise_index-1];
			}
			else
			{
				noise_estimate[frame] = frame;
				noise_index++;
			}
		}
		cout << "noise estimate done:"<<num_noise_frame<<endl;
		//get reverse subband_noisy
		for(n=0; n<sigLength; n++)
			reverse[sigLength - n - 1] = gOut[chan][n] / fChan[chan].midEarCoeff;
		gammaToneFilter(reverse, gOut[chan], fChan[chan], sigLength);
		for(n=0; n<sigLength; n++)
			reverse[sigLength - n - 1] = gOut[chan][n] / fChan[chan].midEarCoeff;
		for(n=0; n<sigLength; n++) 
		{	
			weight_e[n] = 0.0;
		}
		
		for(frame=0; frame < numFrame; frame++)
		{
			if(v_erm[frame][chan]> thd)
			{
				if(frame > 0)
				{
	        			for(n=0; n<OFFSET; n++)
						weight_e[(frame-1)*OFFSET + n] += 0.5 * (1.0 + cos(n * PI/(OFFSET) + PI));
				}
				for(long n=OFFSET; n<WINDOW; n++)
					weight_e[(frame-1)*OFFSET + n] += 0.5 * (1.0 + cos((n-OFFSET) * PI/(OFFSET)));
			}
		}
		for(frame=0; frame < numFrame; frame++)
		{
			if(v_erm[frame][chan]> thd)
			{
				//fft for each frame
				if(num_noise_frame > 0)
				{
					fft( &reverse[ frame * OFFSET],  WINDOW, spec_mag, spec_angle);
					fft( &reverse[ noise_estimate[frame] * OFFSET], WINDOW, noise_spec_mag, noise_spec_angle);
					spec_minus(spec_mag, noise_spec_mag, 512);
					ifft(spec_mag, spec_angle, OFFSET, denoise_data);
				}
				else
				{
					for(long n=0; n<OFFSET; n++)
						denoise_data[n] = reverse[frame * OFFSET + n];
				}				
				
				for(long n=0; n<OFFSET; n++)
					resynth_e[frame*OFFSET + n] += weight_e[frame * OFFSET + n] * denoise_data[n];
			}
			
		}
		delete []gOut[chan];
		delete []hOut[chan];
		delete []s_hOut[chan];
	}
	//save resynth
	short *resynth_s;
	resynth_s = new short[sigLength];
	for(n=0; n<sigLength; n++)
		resynth_s[n] = (short)resynth_e[n];
	CWave resynth_wav_e(resynth_s, sigLength, wav->GetFs());
	sprintf(temp,"%s%s%s_e_resynth.wav",opts.outputDictionary,opts.save_resynth_e_dir, name_id);
	resynth_wav_e.Write(temp);

	delete []noise_frame;
	delete []noise_estimate;
	delete []spec_mag;
	delete []noise_spec_mag;
	delete []spec_angle;
	delete []noise_spec_angle;
	delete []resynth_s;
	delete []resynth_e;
	delete []weight_e;
	delete []reverse;
	delete []Input;
	
}

middleEar initMiddleEar()   
{
	middleEar lF;
	lF.f[0]=20.0;     lF.af[0]=2.347;  lF.bf[0]=0.00561;   lF.tf[0]=74.3;
	lF.f[1]=25.0;     lF.af[1]=2.190;  lF.bf[1]=0.00527;   lF.tf[1]=65.0;
	lF.f[2]=31.5;     lF.af[2]=2.050;  lF.bf[2]=0.00481;   lF.tf[2]=56.3;
	lF.f[3]=40.0;     lF.af[3]=1.879;  lF.bf[3]=0.00404;   lF.tf[3]=48.4;
	lF.f[4]=50.0;     lF.af[4]=1.724;  lF.bf[4]=0.00383;   lF.tf[4]=41.7;
	lF.f[5]=63.0;     lF.af[5]=1.579;  lF.bf[5]=0.00286;   lF.tf[5]=35.5;
	lF.f[6]=80.0;     lF.af[6]=1.512;  lF.bf[6]=0.00259;   lF.tf[6]=29.8;
	lF.f[7]=100.0;    lF.af[7]=1.466;  lF.bf[7]=0.00257;   lF.tf[7]=25.1;
	lF.f[8]=125.0;    lF.af[8]=1.426;  lF.bf[8]=0.00256;   lF.tf[8]=20.7;
	lF.f[9]=160.0;    lF.af[9]=1.394;  lF.bf[9]=0.00255;   lF.tf[9]=16.8;
	lF.f[10]=200.0;   lF.af[10]=1.372; lF.bf[10]=0.00254;  lF.tf[10]=13.8;
	lF.f[11]=250.0;   lF.af[11]=1.344; lF.bf[11]=0.00248;  lF.tf[11]=11.2;
	lF.f[12]=315.0;   lF.af[12]=1.304; lF.bf[12]=0.00229;  lF.tf[12]=8.9;
	lF.f[13]=400.0;   lF.af[13]=1.256; lF.bf[13]=0.00201;  lF.tf[13]=7.2;
	lF.f[14]=500.0;   lF.af[14]=1.203; lF.bf[14]=0.00162;  lF.tf[14]=6.0;
	lF.f[15]=630.0;   lF.af[15]=1.135; lF.bf[15]=0.00111;  lF.tf[15]=5.0;
	lF.f[16]=800.0;   lF.af[16]=1.062; lF.bf[16]=0.00052;  lF.tf[16]=4.4;
	lF.f[17]=1000.0;  lF.af[17]=1.000; lF.bf[17]=0.00000;  lF.tf[17]=4.2;
	lF.f[18]=1250.0;  lF.af[18]=0.967; lF.bf[18]=-0.00039; lF.tf[18]=3.7;
	lF.f[19]=1600.0;  lF.af[19]=0.943; lF.bf[19]=-0.00067; lF.tf[19]=2.6; 
	lF.f[20]=2000.0;  lF.af[20]=0.932; lF.bf[20]=-0.00092; lF.tf[20]=1.0;
	lF.f[21]=2500.0;  lF.af[21]=0.933; lF.bf[21]=-0.00105; lF.tf[21]=-1.2; 
	lF.f[22]=3150.0;  lF.af[22]=0.937; lF.bf[22]=-0.00104; lF.tf[22]=-3.6;
	lF.f[23]=4000.0;  lF.af[23]=0.952; lF.bf[23]=-0.00088; lF.tf[23]=-3.9; 
	lF.f[24]=5000.0;  lF.af[24]=0.974; lF.bf[24]=-0.00055; lF.tf[24]=-1.1;
	lF.f[25]=6300.0;  lF.af[25]=1.027; lF.bf[25]=0.00000;  lF.tf[25]=6.6;
	lF.f[26]=8000.0;  lF.af[26]=1.135; lF.bf[26]=0.00089;  lF.tf[26]=15.3;
	lF.f[27]=10000.0; lF.af[27]=1.266; lF.bf[27]=0.00211;  lF.tf[27]=16.4;
	lF.f[28]=12500.0; lF.af[28]=1.501; lF.bf[28]=0.00488;  lF.tf[28]=11.6;
	return(lF);
}
void gammaToneFilter(float *input, float *output, gammaTone fChan, long sigLength)
{
	float dt, twoPiT, gain, z;
	float f1, f2;
	float x[4], y[4];
	int i;
	long n;
	
	dt= 1/float(SAMPLING_FREQUENCY);

	twoPiT=2 * PI * dt;
	gain = fChan.midEarCoeff * pow(twoPiT * fChan.bw, 4) / 3.0;
	z = exp(-twoPiT * fChan.bw);

	f1 = cos(fChan.cf * twoPiT) * z;
	f2 = sin(fChan.cf * twoPiT) * z;

	for (i=0; i<4; i++)
	{
		fChan.p[i] = 0;
		fChan.q[i] = 0;
	}

	for (n=0; n<sigLength; n++)
	{
		output[n] = fChan.p[3] * gain;
		for (i=0; i<4; i++)
		{
			x[i] = f1*fChan.p[i] - f2*fChan.q[i];
			y[i] = f2*fChan.p[i] + f1*fChan.q[i];
		}

		fChan.p[0] = input[n] * f1 + x[0];
		fChan.q[0] = input[n] * f2 + y[0];
		
		fChan.p[1] = fChan.p[0] + x[1];
		fChan.q[1] = fChan.q[0] + y[1];
		
		fChan.p[2] = fChan.p[1] + x[1] + x[2];
		fChan.q[2] = fChan.q[1] + y[1] + y[2];
		
		fChan.p[3] = fChan.p[2] + x[1] + 2*x[2] + x[3];
		fChan.q[3] = fChan.q[2] + y[1] + 2*y[2] + y[3];
	}
}
void hairCell(float *input, float *output, long sigLength)
{
	float dt, ymdt, xdt, ydt, lplusrdt, rdt, gdt, hdt;
	float kt, c, q, w;
	float replenish, eject, reuptakeandloss, reuptake, reprocess;
	long n;

	dt = 1/float(SAMPLING_FREQUENCY);
	ymdt = MED_Y * MED_M * dt;
	xdt = MED_X * dt;
	ydt = MED_Y * dt;
	lplusrdt = (MED_L + MED_R) * dt;
	rdt = MED_R * dt;
	gdt = MED_G * dt;
	hdt = MED_H; // should be multiplied by dt really

	//initialize variables
	kt = MED_G * MED_A / (MED_A + MED_B);
	c = MED_M * MED_Y * kt / (MED_L * kt + MED_Y * (MED_L + MED_R ));
	q = c * (MED_L + MED_R)/kt;
	w = c * MED_R / MED_X;

	// calculate output
	for (n=0; n<sigLength; n++)
	{
		kt = ((input[n]+MED_A) > 0.0) ? (gdt * (input[n]+MED_A) / (input[n]+MED_A+MED_B)) : 0;

		replenish = (q<MED_M) ? (ymdt - ydt*q) : 0;
        
		eject = kt * q;
		reuptakeandloss = lplusrdt * c;
		reuptake = rdt * c;
		reprocess = xdt * w;

		q = q + replenish - eject + reprocess;
		if (q < 0.0) q = 0.0;

		c = c + eject - reuptakeandloss;
		if (c < 0.0) c = 0.0;
   
        w = w + reuptake - reprocess;
		if (w < 0.0) w = 0.0;
    
		output[n] = hdt * c;
	}
}
float loudnessLevelInPhons(float freq, middleEar lF)
{
	int i;
	float afy, bfy, tfy;
	float ratio;

	if ((freq<20.0) | (freq>12500.0)) {
		fprintf(stderr,"Can't compute a outer/middle ear gain for that frequency\n");
		exit(0);
	}
	
	i = 0;
	while (lF.f[i] < freq) i++;

	ratio = (freq - lF.f[i-1]) / (lF.f[i] - lF.f[i-1]);
	afy = lF.af[i-1] + ratio * (lF.af[i] - lF.af[i-1]);
	bfy	= lF.bf[i-1] + ratio * (lF.bf[i] - lF.bf[i-1]);
	tfy	= lF.tf[i-1] + ratio * (lF.tf[i] - lF.tf[i-1]);
	
	return(4.2 + afy*(DB - tfy) / (1.0 + bfy*(DB - tfy)) );
}
float Train_thd(int numMix, const  vector<vector<float> > &v_erm, int numFrame, int chan)
{
		
  		double *snr_sub;
		snr_sub = new double[numFrame];
		
		double mean , var , thd ;
		mean =0; var = 0;thd = 0;
		double mix_thd , mix_var ;
		mix_thd = 0;mix_var =0;
		double mix_mean = 0.0;
		int idx , frame;
		idx = 0;frame = 0;
		for( frame = 0 ; frame < numFrame; frame++)
		{
			if(v_erm[frame][chan] != 1.0)
				snr_sub[frame] = v_erm[frame][chan]/(1.0 - v_erm[frame][chan]);	
			else
				snr_sub[frame] = 999.9;
			mean += snr_sub[frame];
		}
		mean = mean * 1.0 / numFrame;
		
		for( frame = 0 ; frame < numFrame; frame++)
		{
			var += (snr_sub[frame] - mean) * (snr_sub[frame] - mean);
		}
		var = sqrt( var / numFrame);

		for( frame = 0 ; frame < numFrame; frame++)
		{
			snr_sub[frame] = (snr_sub[frame] - mean)/ var ;
		}
		GMM gmm(1,numMix, frame);     // construct a gmm
		gmm.Read_Feature(snr_sub);
		gmm.Random_Init();
		gmm.Iterate();
		gmm.EM();

		float *gmmmean,*gmmvar;
		float gmmwei = 1;
		int wrong_idx = 0;
		for(int i = 0; i< numMix; i++)
		{
			if(gmm.GetWeight(i) < gmmwei)
			{
				wrong_idx = i;
				gmmwei = gmm.GetWeight(i);
			}
		}
		for(int i = 0; i < numMix; i++ )
		{
			if(i != wrong_idx)
			{
				gmmmean = gmm.GetMean(i);
				gmmvar = gmm.GetVar(i);
				if (gmmmean[0] > mix_mean)
				{
					mix_mean = gmmmean[0];
					mix_var = sqrt(gmmvar[0]);
				}
			}
		}	
		mix_thd = mix_mean -  mix_var;
		thd =  mix_thd  * var  + mean;
		cout<< mean <<"," << var << "," <<mix_mean<<","<<mix_var <<endl;
		delete []snr_sub;
		float thd_erm = thd / (1+thd);
		if( thd_erm > 1)
			return 1;
		else if (thd_erm > 0)
			return thd_erm;
		else
			return 0.0;
}

void fft(float *data_each_frame,  int len, float *spec_mag, complex<double> *spec_angle)
{
	complex<double> *data,*spec;
	data = new complex<double>[512];
	spec = new complex<double>[512];
	for(int i = 0; i < len; i++)
	{
		data[i].real() = data_each_frame[i];
		data[i].imag() = 0;
	}
	for(int i = len; i < 512; i++)
	{
		data[i].real() = 0;
		data[i].imag() = 0;
	}
	discreteFourierFast(data, 512 , spec, ftdFunctionToSpectrum);
	for(int i=0; i < 512 ; i++ )
	{
		spec_mag[i] = (spec[i].real() * spec[i].real()  + spec[i].imag() * spec[i].imag());
		spec_angle[i].real() = spec[i].real() / sqrt(spec_mag[i]);
		spec_angle[i].imag() = spec[i].imag() / sqrt(spec_mag[i]);
	}
	delete []data;
	delete []spec;
}
void spec_minus(float *spec_mag, float *noise_spec_mag, int len)
{
	for(int i = 0 ; i < len ; i++ )
	{
		spec_mag[i] = spec_mag[i] - noise_spec_mag[i];
		if(spec_mag[i] < 0)
		{
			spec_mag[i] = 0;
		}
	}
	
}
void ifft(float *spec_mag, complex<double> *spec_angle, int len, float *denoise_data)
{
	complex<double> *data,*spec;
	data = new complex<double>[512];
	spec = new complex<double>[512];
	for(int i=0; i < 512 ; i++ )
	{
		spec[i].real() = spec_angle[i].real() * sqrt(spec_mag[i]);
		spec[i].imag() = spec_angle[i].imag() * sqrt(spec_mag[i]);
	}
	discreteFourierFast(spec, 512, data, ftdSpectrumToFunction);
	for(int i = 0; i < len; i++)
	{
		denoise_data[i] = data[i].real();
	}
	delete []data;
	delete []spec;
	
}
