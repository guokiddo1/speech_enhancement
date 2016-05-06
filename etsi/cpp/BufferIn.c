/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: BufferIn.c
 * PURPOSE:   General purpose function package for buffer operations
 *
 *-------------------------------------------------------------------------------*/
/*-----------------
 * File Inclusions
 *-----------------*/
#include <stdlib.h>
#include "BufferIn.h"

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
extern BufferIn *
BufInAlloc(int BufSize)
{
    BufferIn *newBuf = (BufferIn *) calloc(1, sizeof(BufferIn));
  if (newBuf == NULL)
    return NULL;
  newBuf->data = (X_FLOAT32*) calloc(BufSize, sizeof(newBuf->data[0]));
  if (newBuf->data == NULL)
    return BufInFree(newBuf);
  newBuf->size = BufSize;
  return newBuf;
}

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
extern BufferIn *
BufInFree(BufferIn * Buf)
{
  if (Buf != NULL) {
    free(Buf->data);
    free(Buf);
  }
  return NULL;
}

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
extern void 
BufInClear(BufferIn * Buf)
{
  int i;

  for (i = 0; i < Buf->size; i++)
    Buf->data[i] = 0;
}

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
extern int 
BufInPut(BufferIn * Buf, X_FLOAT32 * newData, int newDataSize)
{
  int i;

  if (newDataSize > Buf->size)
    newDataSize = Buf->size;
  /*----------------
   * shift old data
   *----------------*/
  for (i = 0; i < Buf->size - newDataSize; i++) {
    Buf->data[i] = Buf->data[i + newDataSize];
  }
  /*----------------------------
   * put new data in the buffer
   *----------------------------*/
  for (; i < Buf->size; i++) {
    Buf->data[i] = *newData++;
  }
  return newDataSize;
}

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
extern X_FLOAT32 *
BufInShiftToPut(BufferIn * Buf, int newDataSize)
{
  if (newDataSize > Buf->size) {
    return NULL;
  } else {
    int i;
    /*----------------
     * shift old data
     *----------------*/
    for (i = 0; i < Buf->size - newDataSize; i++) {
      Buf->data[i] = Buf->data[i + newDataSize];
    }
    return Buf->data + i;
  }
}

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
extern int 
BufInGetLast(BufferIn * Buf, X_FLOAT32 * dataOutBuf, int nbAskSamples)
{
  int i;

  if (nbAskSamples > Buf->size)
    nbAskSamples = Buf->size;

  /*------------------
   * output last data
   *------------------*/
  for (i = Buf->size - nbAskSamples; i < Buf->size; i++) {
    *dataOutBuf++ = Buf->data[i];
  }
  return nbAskSamples;

}

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
extern X_FLOAT32 *
BufInGetPointerToLast(BufferIn * Buf, int nbAskSamples)
{
  if (nbAskSamples > Buf->size) {
    return NULL;
  } else {
    return Buf->data + (Buf->size - nbAskSamples);
  }
}
