#!/bin/bash

cd cpp

icpc -O3 -c KMeans.cpp
icpc -O3 -c GMM.cpp
icpc -O3 -c fft_core.cpp
icpc -O3 -c extractwav.cpp -I/home/guocong/tools/asdk_20130701/asdk -L/home/guocong/tools/asdk_20130701/asdk -lasdk
icpc -O3 -c main.cpp -I/home/guocong/tools/asdk_20130701/asdk -L/home/guocong/tools/asdk_20130701/asdk -lasdk 
icpc -O3 -c MelProc.cpp -I/home/guocong/tools/asdk_20130701/asdk -L/home/guocong/tools/asdk_20130701/asdk -lasdk 
icpc -O3 -o resyth_subband main.o extractwav.o MelProc.o KMeans.o GMM.o fft_core.o -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread -L/home/guocong/tools/asdk_20130701/asdk -lasdk 

#cp Decoder ../

