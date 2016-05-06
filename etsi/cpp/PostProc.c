/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: PostProc.c
 * PURPOSE: Apply post-processing (blind equalization) on a frame of
 *          cepstral coefficients.
 *
 *-------------------------------------------------------------------------------*/
/*-----------------
 * File Inclusions
 *-----------------*/
#include <stdlib.h>

#include "ParmInterface.h"
#include "PostProcExports.h"

/*------------------------
 * Definitions and Macros
 *------------------------*/
#define PP_NB_TOPOSTPROC	(NUM_CEP_COEFF-1)
#define PP_INDEX_ENER		(NUM_CEP_COEFF)
#define PP_MIN_ENER_POSTPROC	((X_FLOAT32) 211)
#define PP_DELTA_ENER_POSTPROC	((X_FLOAT32) 64)

struct PostProcStructX {
  X_FLOAT32 weightLMS[PP_NB_TOPOSTPROC];
  BOOLEAN Noc0;
};


static X_FLOAT32 targetLMS[PP_NB_TOPOSTPROC] =
{
  (X_FLOAT32) - 6.618909,
  (X_FLOAT32) 0.198269,
  (X_FLOAT32) - 0.740308,
  (X_FLOAT32) 0.055132,
  (X_FLOAT32) - 0.227086,
  (X_FLOAT32) 0.144280,
  (X_FLOAT32) - 0.112451,
  (X_FLOAT32) - 0.146940,
  (X_FLOAT32) - 0.327466,
  (X_FLOAT32) 0.134571,
  (X_FLOAT32) 0.027884,
  (X_FLOAT32) - 0.114905
};

static X_FLOAT32 lambda = (X_FLOAT32) 0.0087890625;

/*------------------------
 * start of Encapsulation
 *------------------------*/
/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoPostProcAlloc
 *
 * PURPOSE:       Memory allocation of post-processing structure
 *
 * INPUT:
 *   none
 *
 * OUTPUT:        Post-processing structure
 *
 *
 * RETURN VALUE:  Pointer to post-processing structure
 *
 *
 *---------------------------------------------------------------------------*/
extern PostProcStructX *
DoPostProcAlloc(void)
{
  PostProcStructX *newPost = NULL;
  newPost = (PostProcStructX *)calloc(1, sizeof(PostProcStructX));
  return newPost;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoPostProcInit
 *
 * PURPOSE:       Initialisation of post-processing structure
 *
 * INPUT:
 *   FEParamsX *  Pointer to front end parameter structure
 *
 * OUTPUT:        weightLMS in post-processing structure
 *
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
extern void 
DoPostProcInit(FEParamsX * This)
{
  X_INT16 i;
  PostProcStructX *PPX = This->PPX;

  for (i = 0; i < PP_NB_TOPOSTPROC; i++)
    PPX->weightLMS[i] = 0;
  PPX->Noc0 = This->Noc0;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoPostProc
 *
 * PURPOSE:       Cepstral post-processing for a frame
 *
 * INPUT:
 *   *Coef        Pointer to cepstral coefficients
 *   FEParamsX *  Pointer to front end parameter structure
 *
 * OUTPUT:
 *   *Coef        Pointer to cepstral coefficients
 *
 * RETURN VALUE:
 *   TRUE         A frame is always output
 *
 *---------------------------------------------------------------------------*/
extern BOOLEAN 
DoPostProc(X_FLOAT32 * Coef, FEParamsX * This)
{
  X_INT16 i;
  X_FLOAT32 weightingPar;
  PostProcStructX *PPX = This->PPX;

  weightingPar = (Coef[PP_INDEX_ENER - (PPX->Noc0 ? 1 : 0)] * PP_DELTA_ENER_POSTPROC - PP_MIN_ENER_POSTPROC) / PP_DELTA_ENER_POSTPROC;

  if (weightingPar < 0) {
    weightingPar = 0;
  } else if (weightingPar > 1) {
    weightingPar = lambda;
  } else {
    weightingPar *= lambda;
  }

  for (i = 0; i < PP_NB_TOPOSTPROC; i++) {
    X_FLOAT32 val;
    X_FLOAT32 dif = ((Coef[i] - PPX->weightLMS[i]) - targetLMS[i]);
    Coef[i] = Coef[i] - PPX->weightLMS[i];
    val = dif * weightingPar;
    PPX->weightLMS[i] += val;
  }

  return TRUE;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoPostProcDelete
 *
 * PURPOSE:       Memory free of post-processing structure
 *
 * INPUT:
 *   *PPX         Pointer to post-processing structure
 *
 * OUTPUT:
 *   none
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
extern void 
DoPostProcDelete(PostProcStructX * PPX)
{
  if (PPX != NULL) {
    free(PPX);
  }
}
