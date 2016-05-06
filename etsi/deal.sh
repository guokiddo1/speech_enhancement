#!/bin/bash

cd cpp

icpc -O3 -c 16kHzProcessing.c
icpc -O3 -c AdvFrontEnd.c
icpc -O3 -c BufferIn.c
icpc -O3 -c CompCeps.c
icpc -O3 -c MelProc.c
icpc -O3 -c NoiseSup.c
icpc -O3 -c ParmInterface.c
icpc -O3 -c PostProc.c
icpc -O3 -c rfft.c
icpc -O3 -c VAD.c
icpc -O3 -c WaveProc.c

icpc -O3 -c main.cpp -I/home/guocong/tools/asdk_20130701/asdk -L/home/guocong/tools/asdk_20130701/asdk -lasdk 
icpc -O3 -o etsi_denoise main.o 16kHzProcessing.o AdvFrontEnd.o BufferIn.o CompCeps.o MelProc.o NoiseSup.o ParmInterface.o PostProc.o rfft.o VAD.o WaveProc.o -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread -L/home/guocong/tools/asdk_20130701/asdk -lasdk 

#cp Decoder ../

