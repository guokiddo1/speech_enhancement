############################################################################
# Makefile for building: ASDK
# Author : Shan Yuxiang, Oct. 29, 2008
# Project: ASDK 
# Memo   : General configurations
############################################################################
include ../Makefile.cfg

CXXFLAGS  += -lgomp $(OPENMP)
INCPATH   += 
LIBS      += $(SUBLIBS)  -lasdk $(MKL_LIBS) $(MKL_LIBS) $(MKL_LIBS) -liomp5 -lgomp 

####### Files
SOURCES    = wav_trans.cpp
TARGET     = svf_to_ek

sinclude $(DEPENDENCES)

####### Build rules
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(LINK) $(LFLAGS) -o $(TARGET) ${OBJECTS} $(LIBS)


clean:
	-$(DEL_FILE) $(OBJECTS) $(TARGET) $(DEPENDENCES)

