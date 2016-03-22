#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "NoiseSupExports.h"

#define PI (3.1415926535897932384626433832795)

#define SAMPLING_FREQUENCY		16000 
#define MAX_SIG_LENGTH			150000


#define NUMBER_CHANNEL  64      		        /* maxmimum number of filters */
#define MINCF			50
#define MAXCF			8000

#define BW_CORRECTION       1.019      			/* ERB bandwidth correction 4th order */

struct gammaTone
{
	float cf, bw;
	float midEarCoeff;
	float p[4];
	float q[4];
};

#define MIDDLE_EAR_SIZE 29
#define DB 60.0

struct middleEar
{
	float f[MIDDLE_EAR_SIZE], af[MIDDLE_EAR_SIZE], bf[MIDDLE_EAR_SIZE], tf[MIDDLE_EAR_SIZE];
};

/* hair cell constants from Meddis 1988 paper */
#define MED_Y 5.05
#define MED_G 2000.0
#define MED_L 2500.0
#define MED_R 6580.0
#define MED_X 66.31
#define MED_A 3.0
#define MED_B 300.0
#define MED_H 48000.0
#define MED_M 1.0

void AudiPeriph();

	
void gammaToneFilter(float *input, float *output, gammaTone fChan, long sigLength);  // gammatone filtering

middleEar initMiddleEar();
// parameters of equal-loudness functions from BS3383,"Normal equal-loudness level
// contours for pure tones under free-field listening conditions", table 1.

float loudnessLevelInPhons(float freq, middleEar loudFunc);
//Uses linear interpolation of the look-up tables to compute the loudness level, 
//in phons, of a pure tone of frequency freq using the reference curve for sound 
// pressure level dB.
//The equation is taken from section 4 of BS3383.

void hairCell(float *input, float *output, long sigLength); // Meddis' auditory nerve transdusction



// Correlogram

#define	MAX_DELAY	(SAMPLING_FREQUENCY/80+1)		/* corresponds to 80 Hz */
#define MIN_DELAY	(SAMPLING_FREQUENCY/500)		/* get ratio from 2 ms delay */
#define WINDOW	    (SAMPLING_FREQUENCY/50)		/* use a window of 20 ms */
#define OFFSET		(SAMPLING_FREQUENCY/100)					/* compute acf every 10 ms */

#define PASSBAND	1000 // Passband for the lowpass filter
#define STOPBAND	1200 // Stopband for the lowpass filter
#define RIPPLE		0.01 // Passband and stopband ripple

struct corrLgm{
	float acf[NUMBER_CHANNEL][MAX_DELAY];
	float cross[NUMBER_CHANNEL];
	float pRatio[NUMBER_CHANNEL];
};

void computeACF(); //compute correlogram

void crossCorr(); // compute cross-channel correlation

void globalPitch(); // find global pitch

void lowPass();

//Initial segregation

#define THETAC	0.985
#define THETAA  50
#define THETAP	0.95
#define THETAT  0.85

#define MAX_NUMBER_SEGMENT	400


struct segment
{
	int sFrame;
	int eFrame;

	int number;
	int relation;

	int *count[2];
};

int segmentation(mask *m, int numFrame, int *Unit[2], int *segMark);

void search(mask *m, int numFrame, int *Unit[2], long n, long &numUnit);

int relationship(segment *Seg, int numSegment);

void relationship2(segment *Seg, int numSegment);

void timeCrn(corrLgm *hc);

void initGroup(mask* grp);

void reGroup(float p, mask* Grp);



// Pitch determination

void pitchDtm(mask* Grp);

void findPitch(int &lFrame, int &rFrame, mask* Grp);

void checkPitch(int *pitch, int *buffer, int numFrame, int lFrame, int rFrame, int &sStreak, int &eStreak);

void leftDevelope(corrLgm *hc, mask *m, int *pitch, int *buffer, int sStreak, int lFrame);
			
void rightDevelope(corrLgm *hc, mask *m, int *pitch, int *buffer, int eStreak, int rFrame);

void linearEstimate(int *pitch, int lFrame, int rFrame);

int sumAuto(float m[NUMBER_CHANNEL], float acg[NUMBER_CHANNEL][MAX_DELAY], int lDelay, int rDelay);
	
int pitchCrn(int d1, int d2);



// AM

#define THETAE					0.20

void computeAM();

void amPatt(float *input, float *output, long sigLength, int *pitch);

void hilbert(float *input, float *output, int N);

void computeAMRate(float *input, int chan, int *pitch, int numFrame, mask *amCrn);

float sinModel(float *sig, float a, float freq);



// Final segregation

#define MIN_SEG_LENGTH		5

void finalSeg(mask * Grp);

void reSegmentation(mask *m);

void develope(mask *m, int v,mask * Grp);



//Resynthesis

void resynthesis(mask *m, FILE *fp);



// tools
void kaiserPara(float delta, float transBw, int &fLength, float &beta);
	
void kaiserLowPass(float *filter, int fLength, float beta, float wn);

void kaiserBandPass(float *filter, int fLength, float beta, float wc, float wn);

float bessi0(float x);

void fft(float *inputR, float *inputI, int N, float direct);
void createIBM(float* inputwav, int len, mask* Grp);
