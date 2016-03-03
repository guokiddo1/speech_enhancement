/* HuWang model for monaural speech segregation
   Original code for gammatone filtering and the Meddis model of auditory transduction
   is developed at Univ. of Sheffield and adapted by Guoning Hu and DeLiang Wang, 2003. 

   This software was developed in Microsoft Visual C++ 6.0 */

#include "HuWang.h"

long sigLength;

float Input[MAX_SIG_LENGTH], *gOut[NUMBER_CHANNEL], *hOut[NUMBER_CHANNEL], *hEv[NUMBER_CHANNEL];

gammaTone fChan[NUMBER_CHANNEL];

int numFrame;

corrLgm *corrHc, *corrEv;

int *Pitch;

int *Unit[2];

int segMark[MAX_NUMBER_SEGMENT];

int numSegment;

mask  *AmCrn;

void createIBM(float* inputwav, int len, mask *Grp)
{	
	int chan, frame;

	sigLength = len;
	for(int i = 0 ; i < len ; i++){
		Input[i]=inputwav[i];
	}

	// read input
	numFrame = sigLength/OFFSET;

	// Auditory Periphery
	for (chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		gOut[chan] = new float[sigLength];
		hOut[chan] = new float[sigLength];
	}

	printf(" Auditory Periphery");
	AudiPeriph();

	// Extract features
	printf("\n Correlogram");
	corrHc = new corrLgm[numFrame];
	corrEv = new corrLgm[numFrame];
	Pitch = new int[numFrame];
	for(chan=0; chan<NUMBER_CHANNEL; chan++)
		hEv[chan] = new float[sigLength];
	
	lowPass();	
	computeACF();
	crossCorr();
	globalPitch();
	timeCrn(corrHc);

	for (chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		delete hOut[chan];
		delete hEv[chan];
	}

	// Initial grouping
	printf("\n Initial grouping");

	for(frame=0; frame<numFrame; frame++)
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
			Grp[frame].mark[chan] = ((corrHc[frame].cross[chan]>THETAC) && (corrHc[frame].acf[chan][0]>(THETAA*THETAA))) ? 1 : 0;
			

	Unit[0] = new int[long(numFrame * NUMBER_CHANNEL)];
	Unit[1] = new int[long(numFrame * NUMBER_CHANNEL)];

	initGroup(Grp);


	// Determine Pitch
	printf("\n Pitch determination");
	pitchDtm(Grp);

	// Compute AM pattern;
	AmCrn = new mask[numFrame];
	
	printf("\n Unit labeling");
	computeAM();

	// Final Grouping
	printf("\n Final Grouping");
	finalSeg(Grp);

	for (frame=0; frame<numFrame; frame++)
	{
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
		{
			Grp[frame].mark[chan] = (Grp[frame].mark[chan]>0) ? 1.0:0.0;// MOST_IMPORTANT IBM		
			//fprintf(stdout, "%d ", Grp[frame].mark[chan]);
		}
	}
	delete gOut[chan];
	delete Pitch;
	delete corrHc;
	delete corrEv;
	delete AmCrn;
	delete Unit[0];
	delete Unit[1];
}


//-------------------------------------------------------------------------------------------------

void AudiPeriph()
{	
	float lowerERB, upperERB, spaceERB;
	float cf, phon;
	int chan;

	middleEar loudFunc = initMiddleEar();

	lowerERB = 21.4*log10(MINCF*0.00437 + 1.0);
	upperERB =  21.4*log10(MAXCF*0.00437 + 1.0);
  
	spaceERB = (NUMBER_CHANNEL > 1) ? (upperERB-lowerERB)/(NUMBER_CHANNEL-1) : 0;
	
	for (chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		cf = (pow(10, (lowerERB + chan*spaceERB)/21.4) - 1) / 0.00437 ;
		fChan[chan].cf = cf;
		fChan[chan].bw = 24.7*(cf*0.00437 + 1.0) * BW_CORRECTION;
		
		phon = (loudnessLevelInPhons(cf, loudFunc) - DB);
		fChan[chan].midEarCoeff = pow(10, phon/20); 
	}

	for(chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		gammaToneFilter(Input, gOut[chan], fChan[chan], sigLength);
		hairCell(gOut[chan], hOut[chan], sigLength);
	}
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

//------------------------------------------------------------------------------------------------

void lowPass()
{
	float *filter, beta;
	int chan, fLength, m; 
	long n, tim;
	
	kaiserPara(RIPPLE, float(STOPBAND - PASSBAND) / SAMPLING_FREQUENCY, fLength, beta);
	filter = new float[fLength+1];
	kaiserLowPass(filter, fLength, beta, float(PASSBAND + STOPBAND) / SAMPLING_FREQUENCY);

	int shift=fLength/2;
	for(chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		if( (chan%32)==0 ) printf(".");

		for(n=0; n<sigLength; n++)
		{
			hEv[chan][n] = 0;
		
			for(m=0; m<=fLength; m++)
			{
				tim = n + shift - m;
				if ( (tim >= 0) && (tim < sigLength) ) 
					hEv[chan][n] += hOut[chan][tim] * filter[m];
			}
		}
	}
	
	delete filter;
}

void computeACF()
{
	int numFrame, chan, frame, winsize;
	int delay, step, tim;

	numFrame = sigLength/OFFSET;

	for(chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		if((chan%8) == 0) printf(".");
		
		winsize = int(4 * SAMPLING_FREQUENCY / fChan[chan].cf);
		if(winsize<WINDOW) winsize = WINDOW;

		for(frame=0; frame<numFrame; frame++)
		{
			for(delay=0; delay<MAX_DELAY; delay++)
			{
				corrHc[frame].acf[chan][delay]=0;
				corrEv[frame].acf[chan][delay]=0;
				
				for(step=0; step<winsize; step++)
				{
					tim = (frame+2) * OFFSET - (step+1);
					if ( (tim - delay >= 0)  && (tim<sigLength) )
					{
						corrHc[frame].acf[chan][delay] += hOut[chan][tim] * hOut[chan][tim-delay];
						corrEv[frame].acf[chan][delay] += hEv[chan][tim] * hEv[chan][tim-delay];
					}
				}
				
				corrHc[frame].acf[chan][delay] /= float(winsize);
				corrEv[frame].acf[chan][delay] /= float(winsize);
			}
		}
	}
}

void crossCorr()
{
	int chan, frame, delay;
	float acfHc[NUMBER_CHANNEL][MAX_DELAY], acfEv[NUMBER_CHANNEL][MAX_DELAY], sumHc,sumEv;

	for(frame=0; frame<numFrame; frame++)
	{
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
		{
			sumHc = 0;
			sumEv = 0;
			for(delay=0; delay<MAX_DELAY; delay++) 
			{
				sumHc += corrHc[frame].acf[chan][delay];
				sumEv += corrEv[frame].acf[chan][delay];
			}
				
			sumHc /= float(MAX_DELAY);
			sumEv /= float(MAX_DELAY);

			for(delay=0; delay<MAX_DELAY; delay++) 
			{
				acfHc[chan][delay] = corrHc[frame].acf[chan][delay] - sumHc;
				acfEv[chan][delay] = corrEv[frame].acf[chan][delay] - sumEv;
			}
				
			sumHc = 0;
			sumEv = 0;
			for(delay=0; delay<MAX_DELAY; delay++)
			{
				sumHc += acfHc[chan][delay] * acfHc[chan][delay];
				sumEv += acfEv[chan][delay] * acfEv[chan][delay];
			}
			
			sumHc = sqrt(sumHc / float(MAX_DELAY));
			sumEv = sqrt(sumEv / float(MAX_DELAY));
			if(sumHc !=0 )
			{
				for(delay=0; delay<MAX_DELAY; delay++)
					acfHc[chan][delay] /= sumHc;
			}
			if(sumEv !=0 )
			{
				for(delay=0; delay<MAX_DELAY; delay++)
					acfEv[chan][delay] /= sumEv;
			}
		}
	
		for(chan=0; chan<(NUMBER_CHANNEL-1); chan++)
		{
			corrHc[frame].cross[chan] = 0;
			corrEv[frame].cross[chan] = 0;

			for(delay=0; delay<MAX_DELAY; delay++)
			{
				corrHc[frame].cross[chan] += acfHc[chan][delay] * acfHc[chan+1][delay];
				corrEv[frame].cross[chan] += acfEv[chan][delay] * acfEv[chan+1][delay];
			}
		
			corrHc[frame].cross[chan] /= float(MAX_DELAY);
			corrEv[frame].cross[chan] /= float(MAX_DELAY);
		}
		
		corrHc[frame].cross[NUMBER_CHANNEL-1] = 0;;
		corrEv[frame].cross[NUMBER_CHANNEL-1] = 0;;
	}
}

void globalPitch()
{
	float sumCorr[MAX_DELAY], mp;
	int chan, frame, delay;

	for(frame=0; frame<numFrame; frame++)
	{
		for(delay=0; delay<MAX_DELAY; delay++)
		{
			sumCorr[delay]=0;

			for(chan=0; chan<NUMBER_CHANNEL; chan++)
				sumCorr[delay] += corrHc[frame].acf[chan][delay];
		}

		Pitch[frame] = MIN_DELAY;
		mp = sumCorr[MIN_DELAY];
		for(delay=MIN_DELAY+1; delay<MAX_DELAY; delay++)
		{
			if (sumCorr[delay] > mp){ mp = sumCorr[delay]; Pitch[frame] = delay; }
		}
	}
}
	
//------------------------------------------------------------------------------------------------

int segmentation(mask *m, int numFrame, int *Unit[2], int *segMark)
{
	int frame, chan, numSegment;
	long numUnit, activeUnit;
	
	numSegment = 0;
	numUnit = 0;
	for(frame=1; frame<(numFrame-1); frame++)
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
		{
			if( m[frame].mark[chan] && m[frame-1].mark[chan] && m[frame+1].mark[chan] )
			{		
				Unit[0][numUnit] = chan;
				Unit[1][numUnit] = frame;
				m[frame].mark[chan]=0;

				activeUnit = numUnit;
				numUnit ++;

				while(activeUnit < numUnit)
				{
					search(m, numFrame, Unit, activeUnit, numUnit);
					activeUnit ++;
				}

				segMark[numSegment] = numUnit;
				numSegment ++;
			}
		}

	return(numSegment);
}

int relationship(segment *Seg, int numSegment)
{
	int length, Longest;
	int i, frame;
	int count, n1, n2, n3, n4, number1, number2;
	int s, e;

	length = 0;
	Longest = 0;
	for (i=0; i<numSegment; i++)
	{
		if ((Seg[i].eFrame - Seg[i].sFrame) > length)
		{
			length = Seg[i].eFrame - Seg[i].sFrame;
			Longest = i;
		}
	}

	Seg[Longest].relation = 1;
	count = 0;
	for(frame=Seg[Longest].sFrame; frame<=Seg[Longest].eFrame; frame++)
	{
		n1 = Seg[Longest].count[0][frame];
		n2 = Seg[Longest].count[1][frame];

		if (n1>n2) count ++;
		else count --;
	}	
	if( count <= 0 ) Seg[Longest].relation = -1;

	for(i=0; i<numSegment; i++)
	{
		number1=0;
		number2=0;

		s = Seg[Longest].sFrame;
		if(s < Seg[i].sFrame) s = Seg[i].sFrame;

		e = Seg[Longest].eFrame;
		if(e > Seg[i].eFrame) e = Seg[i].eFrame;

		for(frame=s; frame<=e; frame++)
		{
			n1 = Seg[Longest].count[0][frame];
			n2 = Seg[Longest].count[1][frame];

			n3 = Seg[i].count[0][frame];
			n4 = Seg[i].count[1][frame];

			if((n1 > n2) && (n3 > n4)) number1++;
			else if((n1 < n2) && (n3 < n4)) number1++;
			else number2++;
		}

		if (number1 > (number2*2) ) Seg[i].relation = Seg[Longest].relation;
		else Seg[i].relation = -Seg[Longest].relation;
	}

	return(Longest);
}

void relationship2(segment *Seg, int numSegment)
{
	int i, frame;
	int n1, n2, number1, number2;

	for(i=0; i<numSegment; i++)
	{
		number1 = 0;
		number2 = 0;

		for(frame=Seg[i].sFrame; frame<=Seg[i].eFrame; frame++)
		{
			n1 = Seg[i].count[0][frame];
			n2 = Seg[i].count[1][frame];

			if (n1 >= n2) number1++;
			else number2++;
		}

		if (number1 >= number2) Seg[i].relation = 1;
		else Seg[i].relation = -1;
	}
}

void search(mask *m, int numFrame, int *Unit[2], long n, long &numUnit)
{				
	int chan, frame;

	chan = Unit[0][n];
	frame = Unit[1][n];
	if(frame > 0) 
	{
		if (m[frame-1].mark[chan])
		{
			Unit[0][numUnit] = chan;
			Unit[1][numUnit] = frame-1;

			numUnit ++;
			m[frame-1].mark[chan] = 0;
		}
	}
		
	if(frame < (numFrame-1))
	{
		if(m[frame+1].mark[chan])
		{
			Unit[0][numUnit] = chan;
			Unit[1][numUnit] = frame+1;

			numUnit ++;
			m[frame+1].mark[chan] = 0;
		}
	}

	if(chan > 0)
	{
		if(m[frame].mark[chan-1])
		{
			Unit[0][numUnit] = chan-1;
			Unit[1][numUnit] = frame;

			numUnit ++;
			m[frame].mark[chan-1] = 0;
		}
	}
						
	if(chan < (NUMBER_CHANNEL-1))
	{
		if(m[frame].mark[chan+1])
		{
			Unit[0][numUnit] = chan+1;
			Unit[1][numUnit] = frame;

			numUnit ++;
			m[frame].mark[chan+1]=0;
		}
	}
}
				
void initGroup(mask* Grp)
{
	long n, nUnit, j;
	int i, Longest, frame, chan;

	segment *Seg;
	
	numSegment = segmentation(Grp, numFrame, Unit, segMark);
	Seg = new segment[numSegment];

	nUnit = 0;
	for(i=0; i<numSegment; i++)
	{
		Seg[i].count[0] = new int[numFrame];  
		Seg[i].count[1] = new int[numFrame];

		for(frame=0; frame<numFrame; frame++){ 
			Seg[i].count[0][frame] = 0;  
			Seg[i].count[1][frame] = 0; 
		}

		Seg[i].sFrame = numFrame;  Seg[i].eFrame = 0;

		for(n=nUnit; n<segMark[i]; n++)
		{
			chan = Unit[0][n];  frame = Unit[1][n];
	
			if (Seg[i].sFrame > frame) Seg[i].sFrame = frame;
			if (Seg[i].eFrame < frame) Seg[i].eFrame = frame;

			if (corrHc[frame].pRatio[chan] > THETAP) Seg[i].count[0][frame] ++;
			else Seg[i].count[1][frame] ++;
		}

		Seg[i].number = segMark[i] - nUnit;
		nUnit = segMark[i];
	}

	Longest = relationship(Seg, numSegment);

	for(frame=0; frame<numFrame; frame++)
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
			Grp[frame].mark[chan]=0;

	nUnit = 0;
	for(i=0; i<numSegment; i++)
	{
		for(n=nUnit; n<Seg[i].number+nUnit; n++)
		{
			chan = Unit[0][n]; frame = Unit[1][n];

			if ( (frame >= Seg[Longest].sFrame) && (frame <= Seg[Longest].eFrame) )
				Grp[frame].mark[chan] = Seg[i].relation;
		}
		nUnit = n;
	}

	for(i=0; i<numSegment; i++)
	{
		delete Seg[i].count[0];
		delete Seg[i].count[1];
	}
	delete Seg;
}

void reGroup(float theta, mask* Grp)
{
	int i, frame, chan;
	long n, nUnit, j;

	segment *Seg;
	Seg = new segment[numSegment];
	
	nUnit = 0;
	for(i=0; i<numSegment; i++)
	{
		Seg[i].count[0] = new int[numFrame];  
		Seg[i].count[1] = new int[numFrame];

		for(frame=0; frame<numFrame; frame++){ 
			Seg[i].count[0][frame] = 0;  
			Seg[i].count[1][frame] = 0; 
		}

		Seg[i].sFrame = numFrame;  
		Seg[i].eFrame = 0;

		for(n=nUnit; n<segMark[i]; n++)
		{
			chan = Unit[0][n];  frame = Unit[1][n];
	
			if (Seg[i].sFrame > frame) Seg[i].sFrame = frame;
			if (Seg[i].eFrame < frame) Seg[i].eFrame = frame;

			if (corrHc[frame].pRatio[chan] > theta) Seg[i].count[0][frame] ++;
			else Seg[i].count[1][frame] ++;
		}

		Seg[i].number = segMark[i] - nUnit;
		nUnit = segMark[i];
	}

	relationship2(Seg, numSegment);

	for(frame=0; frame<numFrame; frame++)
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
			Grp[frame].mark[chan]=0;

	nUnit = 0;
	for(i=0; i<numSegment; i++)
	{
		for(j=nUnit; j<Seg[i].number+nUnit; j++)
		{
			chan = Unit[0][j];  frame = Unit[1][j];

			Grp[frame].mark[chan] = Seg[i].relation;
		}
		nUnit = j;
	}

	for(i=0; i<numSegment; i++){
		delete Seg[i].count[0]; 	
		delete Seg[i].count[1];
	}
	delete Seg;
}

void timeCrn(corrLgm *hc)
{
	int frame, chan, delay;
	float mp;

	for(frame=0; frame<numFrame; frame++)
	{
		for(chan=0; chan<NUMBER_CHANNEL; chan++)	
			hc[frame].pRatio[chan] = 0;

		if ( Pitch[frame] >0 )
		{
			for(chan=0; chan<NUMBER_CHANNEL; chan++)	
			{
				mp = hc[frame].acf[chan][Pitch[frame]];
				
				for(delay=MIN_DELAY; delay<MAX_DELAY; delay++)
				{
					if (hc[frame].acf[chan][delay] > mp) mp = hc[frame].acf[chan][delay];	
				}

				hc[frame].pRatio[chan] = hc[frame].acf[chan][Pitch[frame]] / mp;
			}
		}
	}
}

//------------------------------------------------------------------------------------------------

void pitchDtm(mask* Grp)
{
	int lFrame, rFrame;
	int sStreak, eStreak;

	int *buffer;
	buffer = new int[numFrame];
	
	findPitch(lFrame, rFrame, Grp);

	timeCrn(corrHc);

	reGroup(THETAP, Grp);
	
	findPitch(lFrame, rFrame, Grp);

	checkPitch(Pitch, buffer, numFrame, lFrame, rFrame, sStreak, eStreak);

	leftDevelope(corrHc, Grp, Pitch, buffer, sStreak, lFrame);
			
	rightDevelope(corrHc, Grp, Pitch, buffer, eStreak, rFrame);

	linearEstimate(Pitch, lFrame, rFrame);
	
	timeCrn(corrHc);

	reGroup(THETAT, Grp);

	delete buffer;
}

void findPitch(int &lFrame, int &rFrame, mask *Grp)
{
	int frame, chan, delay;
	int n1, n2;
	float mvalue;

	lFrame = numFrame;
	rFrame = 0;
	for(frame=0; frame<numFrame; frame++)
	{
		Pitch[frame]=sumAuto(Grp[frame].mark, corrHc[frame].acf, MIN_DELAY, MAX_DELAY);

		if ((Pitch[frame] > 0) & (lFrame > frame)) lFrame=frame;
		if ((Pitch[frame] > 0) & (rFrame < frame)) rFrame=frame;
			   
		n1 = 0;
		n2 = 0;
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
		{
			mvalue=0;
			for (delay=MIN_DELAY; delay<MAX_DELAY; delay++)
			{
				if(corrHc[frame].acf[chan][delay] > mvalue) mvalue = corrHc[frame].acf[chan][delay];
			}
		
			if(Grp[frame].mark[chan] > 0)
			{
				n1 ++;
				if ((corrHc[frame].acf[chan][Pitch[frame]]/mvalue) > THETAP) n2 ++;
			}
		}
        
		if ((n2*2) <= n1) Pitch[frame]=0;
	}
}

int sumAuto(float m[NUMBER_CHANNEL], float acg[NUMBER_CHANNEL][MAX_DELAY], int lDelay, int rDelay)
{
	int chan, delay, pos;
	float mvalue, sum;
	
	pos=0;
	mvalue=0;
	for (delay=lDelay; delay<rDelay; delay++)
	{	
		float sum=0;
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
		{
			if(m[chan] > 0) sum += acg[chan][delay];
		}
		
		if(sum > mvalue)
		{ 
			mvalue=sum; 
			pos=delay; 
		}
	}

	return(pos);
}
	
int pitchCrn(int d1, int d2)
{
	float lowerP, upperP;

	lowerP = 0.8 * float(d1);
	upperP = 1.2 * float(d1);

	if (d1>d2) upperP = 1.2 * float(d2);
	else lowerP = 0.8 * float(d2);
       
	if( (d1>lowerP) && (d1<upperP) && (d2>lowerP) && (d2<upperP) ) return(1);
	else return(0);
}

void checkPitch(int *pitch, int *buffer, int numFrame, int lFrame, int rFrame, int &sStreak, int &eStreak)
{
	int frame, sFrame, eFrame, tframe;

	for(frame=0; frame<numFrame; frame++)
	{
		buffer[frame] = pitch[frame];
		pitch[frame] = 0;
	}

	sStreak = 0;
	eStreak = 1;

	frame=lFrame;
	sFrame=lFrame;
	eFrame=lFrame;
	
	while(frame<rFrame)
	{
		if( pitchCrn(buffer[frame], buffer[frame+1]) ) eFrame=frame+1;

		else
		{
		    if( (eFrame - sFrame) > (eStreak - sStreak) )
			{
				sStreak = sFrame; eStreak = eFrame;
			}
        
			if( (eFrame-sFrame+1)>2)
			{
				for(tframe=sFrame; tframe<=eFrame; tframe++)
					pitch[tframe] = buffer[tframe];
			}
			
			sFrame=frame+1;
			eFrame=frame+1;
		}
	   
		frame++;
	}

	if( (eFrame - sFrame) > (eStreak - sStreak) )
	{
		sStreak = sFrame; eStreak = eFrame;
	}
            
	if((eFrame-sFrame+1)>2)
	{
		for(tframe=sFrame; tframe<=eFrame; tframe++)
			pitch[tframe] = buffer[tframe];
	}
}

void leftDevelope(corrLgm *hc, mask *m, int *pitch, int *buffer, int sStreak, int lFrame)
{
	int  number, number2;
	float chan_mark[NUMBER_CHANNEL], chan_mark2[NUMBER_CHANNEL];
	int frame, chan, delay;
	float mvalue, k1, k2;
	
	frame = sStreak-1;
	while(frame >= lFrame)
	{
		if(pitch[frame+1])
		{
			for(chan=0; chan<NUMBER_CHANNEL; chan++)
			{
				chan_mark[chan]=0;
				if (chan<40)
				{
					mvalue=0;
					for (delay=MIN_DELAY; delay<MAX_DELAY; delay++)
					{
						if(hc[frame+1].acf[chan][delay] > mvalue) mvalue = hc[frame+1].acf[chan][delay];
					}
		
					if( (m[frame+1].mark[chan] > 0) && ((hc[frame+1].acf[chan][pitch[frame+1]]/mvalue) > THETAP) ) chan_mark[chan]=1;
				}
			}
		}
	 
		if( pitchCrn(buffer[frame], buffer[frame+1]) ) pitch[frame] = buffer[frame];
		else
		{
            number=0; 
			for(chan=0; chan<NUMBER_CHANNEL; chan++)
			{
				chan_mark2[chan]=0;
				if ((m[frame].mark[chan] > 0) && (chan_mark[chan]))
				{
					chan_mark2[chan]=1;
					number++;
				}
			}

			if(number)
			{
				k1 = 0.65 * float(buffer[frame+1] + 1) - 1; 
				if (k1 < MIN_DELAY) k1 = MIN_DELAY; 
				
				k2 = 1.55 * float(buffer[frame+1] + 1) - 1; 
				if (k2 > (MAX_DELAY-1) ) k2 = MAX_DELAY-1;
            
				pitch[frame] = sumAuto(chan_mark2, hc[frame].acf, int(k1), int(k2)); 
                    
				if (pitchCrn(pitch[frame], buffer[frame+1]))
				{
					number2=0;
					for(chan=0; chan<40; chan++)
					{
						mvalue=0;
						for (int delay=MIN_DELAY; delay<MAX_DELAY; delay++)
						{
							if(hc[frame].acf[chan][delay] > mvalue) mvalue = hc[frame].acf[chan][delay];
						}
		
						if( (m[frame].mark[chan] > 0) && ((hc[frame].acf[chan][pitch[frame]]/mvalue) > THETAP) ) number2++;
					}
        
					if ((number2*2) > number) buffer[frame] = pitch[frame];	
					else { pitch[frame] = 0; buffer[frame] = buffer[frame+1]; }
				}
                else { pitch[frame] = 0; buffer[frame] = buffer[frame+1]; }
			}
			else buffer[frame] = buffer[frame+1];
		}
		
		frame--;
	}
}

void rightDevelope(corrLgm *hc, mask *m, int *pitch, int *buffer, int eStreak, int rFrame)
{
	int  number, number2;
	float chan_mark[NUMBER_CHANNEL], chan_mark2[NUMBER_CHANNEL];
	int frame, chan, delay;
	float mvalue, k1, k2;

	frame = eStreak+1;
	while(frame <= rFrame)
	{
		if(pitch[frame-1])
		{
			for(chan=0; chan<NUMBER_CHANNEL; chan++)
			{
				chan_mark[chan]=0;
				if (chan<40)
				{
					mvalue=0;
					for (delay=MIN_DELAY; delay<MAX_DELAY; delay++)
					{
						if(hc[frame-1].acf[chan][delay] > mvalue) mvalue = hc[frame-1].acf[chan][delay];
					}
		
					if( (m[frame-1].mark[chan] > 0) && ((hc[frame-1].acf[chan][pitch[frame-1]]/mvalue) > THETAP) ) chan_mark[chan]=1;
				}
			}
		}
	 
		if( pitchCrn(buffer[frame], buffer[frame-1]) ) pitch[frame] = buffer[frame];
		else
		{
            number=0; 
			for(chan=0; chan<NUMBER_CHANNEL; chan++)
			{
				chan_mark2[chan]=0;
				if ((m[frame].mark[chan]) && (chan_mark[chan]))
				{
					chan_mark2[chan]=1;
					number++;
				}
			}

			if(number)
			{
				k1 = 0.65 * float(buffer[frame-1] + 1) - 1; 
				if (k1 < MIN_DELAY) k1 = MIN_DELAY; 
				
				k2 = 1.55 * float(buffer[frame-1] + 1) - 1; 
				if (k2 > (MAX_DELAY-1) ) k2 = MAX_DELAY-1;
            
				pitch[frame] = sumAuto(chan_mark2, hc[frame].acf, int(k1), int(k2)); 
                    
				if (pitchCrn(pitch[frame], buffer[frame-1]))
				{
					number2=0;
					for(chan=0; chan<40; chan++)
					{
						mvalue=0;
						for (delay=MIN_DELAY; delay<MAX_DELAY; delay++)
						{
							if(hc[frame].acf[chan][delay] > mvalue) mvalue = hc[frame].acf[chan][delay];
						}
		
						if( (m[frame].mark[chan] > 0) && ((hc[frame].acf[chan][pitch[frame]]/mvalue)>THETAP) ) number2++;
					}
        
					if ((number2*2) > number) buffer[frame] = pitch[frame];	
					else { pitch[frame] = 0; buffer[frame] = buffer[frame-1]; }
				}
                else { pitch[frame] = 0; buffer[frame] = buffer[frame-1]; }
			}
			else buffer[frame] = buffer[frame-1];
		}
		
		frame++;
	}
}

void linearEstimate(int *pitch, int lFrame, int rFrame)
{
	int frame, judge, sd, tframe;
	float k, pd;

	frame=0; 
	judge=0; 
	sd=0;
	while(frame <= rFrame)
	{
		if(pitch[frame])
		{
			if(judge)
			{
			    k = float(pitch[frame]-pitch[sd]) / float(frame-sd);
				
				for(tframe=sd+1; tframe<frame; tframe++)
				{
					pd = pitch[sd] + k * float(tframe-sd);
					
					pitch[tframe] = int(pd);
					if ( (pd - pitch[tframe]) >= 0.5) pitch[tframe]++;
				}
                    
				judge=0;    
			}
			sd=frame;
        }
		else if(sd) judge=1;
	   
		frame++;
	}
}

//------------------------------------------------------------------------------------------------

void computeAM()
{
	float *temp1; 
	float *temp2;
	float *inputR, *inputI;
	long n;

	temp1 = new float[sigLength];
	temp2 = new float[sigLength];

	int N = int(ceil(log(sigLength*1.0)/log(2*1.0)));
	long sigL = long(pow(2.0, N));

	inputR = new float[sigL];
	inputI = new float[sigL];

	for(int chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		if ((chan%8) == 0) printf(".");
		for(n=0; n<sigLength; n++)
			temp1[n] = (gOut[chan][n]>0) ? gOut[chan][n] : 0;

		amPatt(temp1, temp2, sigLength, Pitch);

		for(n=0; n<sigL; n++)
		{
			inputR[n] = (n < sigLength) ? temp2[n] : 0;
			inputI[n] = 0;
		}
		hilbert(inputR, inputI, N);

		for(n=0; n<sigLength; n++)
			temp2[n] /= sqrt(inputR[n]*inputR[n]+inputI[n]*inputI[n])/sigL + 1e-30;
		
		computeAMRate(temp2, chan, Pitch, sigLength/OFFSET, AmCrn);
	}

	delete temp1;
	delete temp2;
	delete inputR;
	delete inputI;
}

void amPatt(float *input, float *output, long sigLength, int *pitch)
{
	float *filter;

	for(long n=0; n<sigLength; n++)
		output[n] = 0;

	int numFrame = sigLength/OFFSET;
	for (int frame=2; frame<numFrame; frame+=5)
	{
		float mp = 0;
		int count = 0;
		for (int f=frame-2; f<=frame+2; f++)
		{
			if ( (f>=0) && (f<numFrame) )
			{	
				if (pitch[f]>0){ mp += pitch[f]; count ++; }
			}
		}

		if (count>0)
		{
			mp /= float(count);
			float beta;
			int fLength;

			kaiserPara(0.01, 0.4 / mp, fLength, beta);		
			filter = new float[fLength+1];
			kaiserBandPass(filter, fLength, beta, 2.1 /mp, 0.7 /mp);

			for(long n=(frame-2)*OFFSET; n<(frame+3)*OFFSET; n++)
			{
				if(n<sigLength)
				{
					for(int m=0; m<=fLength; m++)
					{				
						long tim = n + fLength/2 - m;
						if ( (tim >= 0) && (tim < sigLength) ) 
							output[n] += input[tim] * filter[m];
					}
				}
			}
			delete filter;
		}
	}
}

void hilbert(float *inputR, float *inputI, int N)
{
	long n, sigL;
	
	sigL = long(pow(2.0, N));
	
	fft(inputR, inputI, N, 1);
	
	for(n=1; n<sigL/2; n++)
	{
		inputR[n] *= 2;
		inputI[n] *= 2;
	}

	for(n=sigL/2+1; n<sigL; n++)
	{
		inputR[n]=0;
		inputI[n]=0;
	}

	fft(inputR, inputI, N, -1);
}

void computeAMRate(float *input, int chan, int *pitch, int numFrame, mask *amCrn)
{			
	int frame;
	float freq, sEng, error, error2;
	long shift, tim;

	for(frame=0; frame<numFrame; frame++)
	{
		shift = (frame - 1) * OFFSET;
		
		sEng=0;
		if ( (frame <(numFrame-1) ) && (frame>0) )
		{
			for(tim=shift; tim<(shift + WINDOW); tim++)
				sEng += input[tim] * input[tim];
		}

		if ( (pitch[frame] > 0) && (sEng > 0) )
		{
			freq = SAMPLING_FREQUENCY / float(pitch[frame]);

			error = sinModel(input + shift, 1, freq);
			
			amCrn[frame].mark[chan] = ((error/sEng) > THETAE) ? 0 : 1;
		}

		else amCrn[frame].mark[chan] = 0;		
	}
}

float sinModel(float *sig, float a, float freq)
{
	int tim;
	float pa, pb, phase;
	float s, error1, error2;
	
	pa=0, pb=0;
	freq *= 2 * PI / SAMPLING_FREQUENCY;

	for(tim=0; tim<WINDOW; tim++)
	{
		pa += cos(freq * tim) * sig[tim];
		pb += sin(freq * tim) * sig[tim];
	}

	phase = (pa == 0) ? PI / 2 : atan(-pb/pa);

	error1 = 0;
	error2 = 0;
	for (tim=0; tim<WINDOW; tim++)
	{
		s = a * cos(freq*tim + phase);
		error1 += (sig[tim]-s)*(sig[tim]-s);
		error2 += (sig[tim]+s)*(sig[tim]+s);
	}

	if (error1 < error2) return(error1);
	else{ phase += PI;  return(error2); }
}

//------------------------------------------------------------------------------------------------

void finalSeg(mask * Grp)
{
	int frame, chan;
	
	mask *m;

	m = new mask[numFrame];

	// Generate segments for unresolved harmonics
	for(frame=0; frame<numFrame; frame++)
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
			m[frame].mark[chan] = ( (Grp[frame].mark[chan] == 0) && (AmCrn[frame].mark[chan]) && (corrEv[frame].cross[chan] > THETAC) ) ? 1 : 0;
		
	reSegmentation(m);
	
	for(frame=0; frame<numFrame; frame++)
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
			Grp[frame].mark[chan] += m[frame].mark[chan] * 2; 

	// Adjust foreground stream according to unit labels
	for(frame=0; frame<numFrame; frame++)
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
			m[frame].mark[chan] = ( (Grp[frame].mark[chan] == 1) && (corrHc[frame].pRatio[chan] > THETAT) ) ? 1 : 0;
	
	reSegmentation(m);
	
	for(frame=0; frame<numFrame; frame++)
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
		{
			if( (Grp[frame].mark[chan] == 1) && (corrHc[frame].pRatio[chan] > THETAT) )
				Grp[frame].mark[chan] = (m[frame].mark[chan] == 1) ? 1 : -2;
		}

	
	for(frame=0; frame<numFrame; frame++)
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
			m[frame].mark[chan] = ( (Grp[frame].mark[chan] == 1) && (corrHc[frame].pRatio[chan] <= THETAT) ) ? 1 : 0;
	
	reSegmentation(m);
	
	for(frame=0; frame<numFrame; frame++)
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
		{
			if( (Grp[frame].mark[chan] == 1) && (corrHc[frame].pRatio[chan] <= THETAT) )
				Grp[frame].mark[chan] = (m[frame].mark[chan] == 1) ? -1 : -2;
		}

	for(frame=0; frame<numFrame; frame++)
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
			m[frame].mark[chan] = ( (Grp[frame].mark[chan] == -2) ) ? 1 : 0;
	
	develope(m, -1, Grp);
	
	for(frame=0; frame<numFrame; frame++)
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
		{
			if( (Grp[frame].mark[chan] == -2) || (Grp[frame].mark[chan] == 2) )
				Grp[frame].mark[chan] = 1;
		}
	
	for(frame=0; frame<numFrame; frame++)
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
			m[frame].mark[chan] = ( (AmCrn[frame].mark[chan] == 1) && (Grp[frame].mark[chan] == 0) ) ? 1 : 0;

	develope(m, 1, Grp);

	delete m;
}

void reSegmentation(mask *m)
{	
	long nUnit, n;
	int chan, frame, i;

	segment *Seg;

	numSegment = segmentation(m, numFrame, Unit, segMark);
	Seg = new segment[numSegment];

	nUnit = 0;
	for(i=0; i<numSegment; i++)
	{
		Seg[i].sFrame = numFrame;
		Seg[i].eFrame = 0;

		for(n=nUnit; n<segMark[i]; n++)
		{
			chan = Unit[0][n];
			frame = Unit[1][n];
	
			if (Seg[i].sFrame > frame) Seg[i].sFrame = frame;
			if (Seg[i].eFrame < frame) Seg[i].eFrame = frame;
		}

		Seg[i].number = segMark[i] - nUnit;
		nUnit = segMark[i];
	}

	for(frame=0; frame<numFrame; frame++)
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
			m[frame].mark[chan]=0;

	nUnit = 0;
	for(i=0; i<numSegment; i++)
	{
		if ( (Seg[i].eFrame - Seg[i].sFrame + 1) >= MIN_SEG_LENGTH )
		{
			for(n=nUnit; n<Seg[i].number+nUnit; n++)
			{
				chan = Unit[0][n];
				frame = Unit[1][n];

				m[frame].mark[chan] = 1;
			}
		}
		nUnit += Seg[i].number;
	}

	delete Seg;
}

void develope(mask *m, int v,mask * Grp)
{
	int chan, frame, judge;
	
	judge=1;
	while(judge)
	{
		judge=0;

		for(frame=0; frame<numFrame; frame++)
			for(chan=0;chan<NUMBER_CHANNEL; chan++)
			{
				if(Grp[frame].mark[chan] == v)
				{
					if(frame > 0)
					{
						if( (m[frame-1].mark[chan]) && (Grp[frame-1].mark[chan] != v) )
						{
							Grp[frame-1].mark[chan] = v;
							judge++;
						}
					}

					if(frame < (numFrame-1))
					{
						if((m[frame+1].mark[chan]) && (Grp[frame+1].mark[chan] != v))
						{
							Grp[frame+1].mark[chan] = v;
							judge++;
						}
					}

					if(chan > 0)
					{
						if((m[frame].mark[chan-1]) && (Grp[frame].mark[chan-1] != v))
						{
							Grp[frame].mark[chan-1] = v;
							judge++;
						}
					}

					if(chan < (NUMBER_CHANNEL-1))
					{
						if((m[frame].mark[chan+1]) && (Grp[frame].mark[chan+1] != v))
						{
							Grp[frame].mark[chan+1] = v;
							judge++;
						}
					}
				}
			}
	}
}

//------------------------------------------------------------------------------------------------

void resynthesis(mask *m, FILE *fp)
{
	long n; 
	int chan, frame;
	float *resynth, *weight, *reverse;

	resynth = new float[sigLength];
	weight = new float[sigLength];
	reverse = new float[sigLength];

	for(n=0; n<sigLength; n++)
		resynth[n]=0;

	for(chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		for(n=0; n<sigLength; n++)
			reverse[sigLength - n - 1] = gOut[chan][n] / fChan[chan].midEarCoeff;
        
		gammaToneFilter(reverse, gOut[chan], fChan[chan], sigLength);

		for(n=0; n<sigLength; n++)
			reverse[sigLength - n - 1] = gOut[chan][n] / fChan[chan].midEarCoeff;
        
		for(n=0; n<sigLength; n++) 
			weight[n] = 0.0;
        
		for(frame=0; frame<numFrame; frame++)
		{
			if(m[frame].mark[chan]>0)
			{
				if(frame > 0){
        			for(n=0; n<OFFSET; n++)
						weight[(frame-1)*OFFSET + n] += 0.5 * (1.0 + cos(n * PI/(OFFSET) + PI));
				}
				
				for(long n=OFFSET; n<WINDOW; n++)
					weight[(frame-1)*OFFSET + n] += 0.5 * (1.0 + cos((n-OFFSET) * PI/(OFFSET)));
			}
		}

		
		for(n=0; n<sigLength; n++)
			resynth[n] += weight[n] * reverse[n];
	}
	
	for(n=0; n<sigLength; n++)
		fprintf(fp, "%f\n", resynth[n]);

	delete resynth;
	delete weight;
	delete reverse;
}

//-----------------------------------------------------------------------------------------------------

void kaiserPara(float delta, float transBw, int &fLength, float &beta)
{
	float a, len;
	
	a= -20 * log10(delta);

	if (a <= 21) beta = 0;
	else if (a<= 50) beta = 0.5842 * pow(a-21, (float)0.4) + 0.07889 * (a-21);
	else beta = 0.1102 * (a - 8.7);

	len = (a - 7.95) / 14.36 / transBw;
	fLength = int(len);
	if ((len - fLength) < 0.5) fLength++;
	else fLength+=2;

	if (fLength%2 != 0) fLength++;
}
	
void kaiserLowPass(float *filter, int fLength, float beta, float wn)
{
	int tim, step;
	float k, sum;

	for (tim=0; tim<=fLength; tim++)
	{
		k = 2*tim/float(fLength) - 1;
		filter[tim] = bessi0( beta*sqrt(1- k*k)) / bessi0( beta );
	}

	sum=0;
	for (tim=0; tim<=fLength; tim++)
	{
		step = tim - fLength/2;
		if (step !=0) filter[tim] *= sin(wn * PI * step) / PI / step;
		else filter[tim] *= wn;

		sum += filter[tim];
	}
}

void kaiserBandPass(float *filter, int fLength, float beta, float wc, float wn)
{
	int tim;

	kaiserLowPass(filter, fLength, beta, wn);

	for (tim=0; tim<=fLength; tim++)
		filter[tim] *= 2 * cos((tim-fLength/2) * wc * PI);
}

float bessi0(float x)
{
	float ax,ans;
	float y;

	ax = fabs(x);
	if (ax < 3.75)
	{
		y = x/3.75;
		y *= y;
		ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 + y * (0.2659732 + y * (0.360768e-1 + y*0.45813e-2)))));
	}
	else
	{
		y = 3.75 / ax;
		ans = (exp(ax) / sqrt(ax)) * (0.39894228 + y * (0.1328592e-1 + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2 + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1 + y * 0.392377e-2))))))));
	}
	return ans;
}

void fft(float *inputR, float *inputI, int N, float direct)
{
	long sigL, i, j, k, n, period, twoPeriod;
	float tmpR, tmpI, uR, uI, wR, wI;

	sigL = long(pow(2.0, N));

	j = 1;
	for(i=1; i<sigL; i++)
	{
		if(i < j)
		{
			tmpR = inputR[j-1];
			tmpI = inputI[j-1];

			inputR[j-1] = inputR[i-1];
			inputI[j-1] = inputI[i-1];

			inputR[i-1] = tmpR;
			inputI[i-1] = tmpI;
		}

		k = sigL/2;
		while (k < j){ j -=  k;	k /= 2;	}
		j += k;
	}

	for(n=1; n<=N; n++ )
    {  
		twoPeriod = long(pow(2.0, n));
        period = twoPeriod/2;
        uR = 1.0; 
        uI = 0.0; 
        wR = cos( PI/period ); 
        wI = -1.0 * sin( PI/period * direct);

        for(j=0; j<period; j++ ) 
        {  
			for(i=j; i<sigL; i+=twoPeriod)
			{
				tmpR = inputR[i+period]*uR - inputI[i+period]*uI;
                tmpI = inputR[i+period]*uI + inputI[i+period]*uR;
				
				inputR[i+period] = inputR[i] - tmpR; 
				inputI[i+period] = inputI[i] - tmpI; 
				inputR[i] += tmpR ; 
				inputI[i] += tmpI; 
			}
			tmpR = uR*wR - uI*wI; 
			tmpI = uR*wI + uI*wR; 
			uR = tmpR; 
			uI = tmpI; 
		} 
	} 
}
