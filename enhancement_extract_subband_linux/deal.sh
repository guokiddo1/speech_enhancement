#!/bin/bash

cd cpp


icpc -O3 -c extractwav.cpp -I/home/guocong/tools/asdk_20130701/asdk -L/home/guocong/tools/asdk_20130701/asdk -lasdk
icpc -O3 -c main.cpp -I/home/guocong/tools/asdk_20130701/asdk -L/home/guocong/tools/asdk_20130701/asdk -lasdk 
icpc -O3 -c MelProc.cpp -I/home/guocong/tools/asdk_20130701/asdk -L/home/guocong/tools/asdk_20130701/asdk -lasdk 
icpc -O3 -c show_IBM.cpp -I/home/guocong/tools/asdk_20130701/asdk -L/home/guocong/tools/asdk_20130701/asdk -lasdk 
icpc -O3 -o enhance_extract_subband main.o extractwav.o show_IBM.o MelProc.o -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread -L/home/guocong/tools/asdk_20130701/asdk -lasdk 

#cp Decoder ../

