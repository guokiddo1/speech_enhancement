/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: BufferIn.h
 * PURPOSE:   General purpose function package for buffer operations
 *
 *-------------------------------------------------------------------------------*/
#ifndef _BUFFERIN_H
#define _BUFFERIN_H

#include "ParmType.h"

typedef struct BufferIn BufferIn;

/*-----------------------------------------------------------------------------
Sample 1:
---------

BufferIn *BufI = BufInAlloc (MAX_FRAME_LENGTH);         // Allocation
X_FLOAT32 *ptSig;

// Initialisation (optional)
ptSig = BufInShiftToPut (BufI, frameLength-frameShift); // Where to put samples
ReadSignal (ptSig, frameLength - frameShift);           // Read samples

while (!EOF)
{
  X_FLOAT32 *curFrame;
	
  ptSig = BufInShiftToPut (BufI, frameShift);           // Where to put samples
  ReadSignal (ptSig, frameShift);                       // Read samples
  curFrame = BufInGetPointerToLast (BufI, frameLength); // Get a frame
  Calcul (curFrame);                                    // Process the frame
}
BufI = BufInFree (BufI);                                // Memory free


Sample 2:
---------

BufferIn *BufI = BufInAlloc (MAX_FRAME_LENGTH);    // Allocation of the buffer
X_FLOAT32 SigRead[MAX_FRAME_LENGTH];

   // Initialisation (optional)
ReadSignal (SigRead, frameLength - frameShift);    // Read samples
BufInPut (BufI, SigRead, frameLength - frameShift);// Put samples in the buffer

while (!EOF)
{
  X_FLOAT32 currentFrame[FRAME_LEGNTH];
	
  ReadSignal (SigRead, frameShift);              // Read samples
  BufInPut (BufI, SigRead, frameShift);          // Put samples in the buffer
  BufInGet (BufI, currentFrame, frameLength);    // Get a frame from the buffer
  Process (currentFrame);                        // Process the frame
}
BufI = BufInFree (BufI);                         // Memory free of the buffer
-----------------------------------------------------------------------------*/

struct BufferIn {
  int size;			/* size of the buffer */
  X_FLOAT32 *data;		/* data buffer */
};

#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------------
 * FUNCTION NAME: BufInAlloc
 *
 * PURPOSE:       Memory allocation of a new buffer. The buffer is filled in
 *                with zeroes.
 *
 * INPUT:
 *   BufSize      Buffer size
 *
 * OUTPUT:
 *
 *
 * RETURN VALUE:
 *   NULL         In case of memory allocation error
 *   newBuf       Buffer pointer otherwise
 *
 *---------------------------------------------------------------------------*/
  extern BufferIn *BufInAlloc(int maxSize);

/*----------------------------------------------------------------------------
 * FUNCTION NAME: BufInFree
 *
 * PURPOSE:       Memory free of a buffer.
 *
 * INPUT:
 *   Buf          Pointer to buffer
 *
 * OUTPUT:
 *   none
 *
 * RETURN VALUE:
 *   NULL         In any case, so one can write Buf = BufInFree (Buf);
 *
 *---------------------------------------------------------------------------*/
  extern BufferIn *BufInFree(BufferIn * Buf);

/*----------------------------------------------------------------------------
 * FUNCTION NAME: BufInClear
 *
 * PURPOSE:       Fill in the buffer with zeroes.
 *
 * INPUT:
 *   Buf          Pointer to buffer
 *
 * OUTPUT:
 *                Zeroes in Buf->data[]
 *
 * RETURN VALUE:
 *   NULL         In any case, so one can write Buf = BufInFree (Buf) ;
 *
 *---------------------------------------------------------------------------*/
  extern void BufInClear(BufferIn * Buf);

/*----------------------------------------------------------------------------
 * FUNCTION NAME: BufInPut
 *
 * PURPOSE:       Put data from newData in the buffer. Size of data put in the
 *                buffer is max (newDataSize, Buf->size).
 *
 * INPUT:
 *   Buf          Pointer to buffer
 *   newData      Pointer to new data
 *   newDataSize  Size of new data to put in the buffer
 *
 * OUTPUT:
 *                newData in Buf->data[]
 *
 * RETURN VALUE:
 *   newDataSize  Number of samples of newData put in the buffer
 *
 *---------------------------------------------------------------------------*/
  extern int BufInPut(BufferIn * Buf, X_FLOAT32 * newData, int newDataSize);

/*----------------------------------------------------------------------------
 * FUNCTION NAME: BufInShiftToPut
 *
 * PURPOSE:       Shift data in the buffer to prepare to put new data in the
 *                buffer. Return the location where to put new data.
 *
 * INPUT:
 *   Buf          Pointer to buffer
 *   newDataSize  Size of new data to prepare to put in the buffer
 *
 * OUTPUT:
 *                Shift Buf->data[]
 *
 * RETURN VALUE:
 *   NULL         If newDataSize > Buf->size
 *   Buf->data+i  Otherwise a pointer where to insert new data
 *
 *---------------------------------------------------------------------------*/
  extern X_FLOAT32 *BufInShiftToPut(BufferIn * Buf, int newDataSize);

/*----------------------------------------------------------------------------
 * FUNCTION NAME: BufInGetLast
 *
 * PURPOSE:       Output last nbAskSamples of Buf->data[] to dataOutBuf[]
 *
 * INPUT:
 *   Buf          Buffer pointer
 *   nbAskSamples Number of samples to output, i.e. max (nbAskSamples, Buf->size)
 *
 * OUTPUT:
 *                output Buf->data[] in dataOutBuf[]
 *
 * RETURN VALUE:
 *   nbAskSamples Number of samples written in dataOutBuf
 *
 *---------------------------------------------------------------------------*/
  extern int BufInGetLast(BufferIn * Buf, X_FLOAT32 * dataOutBuf, int nbAskSamples);

/*----------------------------------------------------------------------------
 * FUNCTION NAME: BufInGetPointerToLast
 *
 * PURPOSE:       Return pointer to last nbAskSamples samples of Buf->data[]
 *
 * INPUT:
 *   Buf          Pointer to buffer
 *   nbAskSamples Number of samples to point to
 *
 * OUTPUT:
 *
 * RETURN VALUE:
 *   NULL         If nbAskSamples > Buf->size
 *   pointer      Otherwise pointer to last nbAskSamples samples of Buf->data[]
 *
 *---------------------------------------------------------------------------*/
  extern X_FLOAT32 *BufInGetPointerToLast(BufferIn * Buf, int nbAskSamples);

#ifdef __cplusplus
}
#endif

#endif
