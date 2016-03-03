#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


//#define PI (3.1415926535897932384626433832795)

#define SAMPLING_FREQUENCY		8000 
#define MAX_SIG_LENGTH			150000


#define NUMBER_CHANNEL  25      		        /* maxmimum number of filters */
#define MINCF			80
#define MAXCF			4000

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






// tools
void kaiserPara(float delta, float transBw, int &fLength, float &beta);
	
void kaiserLowPass(float *filter, int fLength, float beta, float wn);

void kaiserBandPass(float *filter, int fLength, float beta, float wc, float wn);

float bessi0(float x);

void fft(float *inputR, float *inputI, int N, float direct);
//void createIBM(float* inputwav, int len, mask* Grp);