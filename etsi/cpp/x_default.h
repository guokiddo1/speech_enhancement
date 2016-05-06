/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
#ifndef _X_DEFAUT_H_
#define _X_DEFAUT_H_


/**************** DEFINITION DE SYMBOLE PAR COMPILATEUR ***************************/

#if defined(__BORLANDC__)
#define X_BORL
#if defined(__WIN32__)
#define X_BORL_32
#elif defined(__DPMI32__)
#define X_BORL_DPMI32
#elif defined(__DPMI16__)
#define X_BORL_DPMI16
#else
#define X_BORL_16
#endif
#endif				/* Borland */

#if defined(_MSC_VER)
#if defined(MSDOS)
#define X_MSC_16
#elif defined(WIN32)
#define X_MSC_32
#endif
#endif				/* Microsoft */

#if defined(__ZTC__)
#define X_ZTC_32
#endif				/* Zortech */

#if defined(__WATCOMC__)
#define X_WATC_32
#endif				/* Watcom */

#if defined(unix)
#if defined(__osf__)
#define X_UNIX_64
#else
#define X_UNIX_32
#endif
#endif				/* Unix */

#if defined(_TMS320C30)
#define X_TMS_C30
#endif				/* TMS C30 */
#if defined(__HOS_OS2__)
#define X_OS2_32
#endif


/************************* GESTION DU SYSTEME *************************************/

#if defined(X_BORL) || defined(X_MSC_16) || defined(X_MSC_32) \
|| defined(X_ZTC_32) || defined(X_WATC_32)
#define X_PCDOS
#define X_SYSTEM_STRING     "PCDOS"
/* miscellanous pour Zortech */
#if defined(X_ZTC_32)
#define _dos_getdrive       dos_getdrive
#define _dos_setdrive       dos_setdrive
#define _dos_setfileattr    dos_setfileattr
#endif
#elif defined(X_UNIX_32) || defined(X_UNIX_64)
#define X_UNIX
#define X_SYSTEM_STRING     "UNIX"
#elif defined(X_TMS_C30)
#define X_TMSC30
#define X_SYSTEM_STRING     "TMSC30"
#elif defined(X_OS2_32)
#define X_OS2
#define X_SYSTEM_STRING      "OS/2"
#else
#error Macro de nom du systeme non definies X_SYSTEM
#endif				/* Gestion systeme */


/************** PACKING DES DATA POUR ECRITURE DANS LES FICHIERS ******************/
#if defined(X_UNIX_32) || defined(X_UNIX_64)
#define X_PACKING_DATA  4
#endif				/* Gestion Packing */


/******************************* GESTION cdecl ************************************/

#if defined(X_MSC_16)
#define X_USER  cdecl
#elif defined(X_MSC_32)
#define X_USER  __cdecl
#elif defined(X_BORL_16)
#define X_USER  _USERENTRY
#elif defined(X_BORL_32) || defined(X_BORL_DPMI32) || defined(X_BORL_DPMI16)
#define X_USER  __cdecl
#else
#define X_USER
#endif				/* Gestion cdecl */

/******************************* GESTION X_DLLEXPORT ************************************/

#if defined(X_MSC_32)
#define X_DLLEXPORT  extern __declspec(dllexport)
#define X_DLLEXPORT_C(type)   extern __declspec(dllexport) type X_USER
#define X_DLLEXPORT_H(type)   extern __declspec(dllexport) type X_USER
#elif defined(X_OS2)
#define X_DLLEXPORT_C(type) extern type _Export
#define X_DLLEXPORT_H(type) extern type
#else
#define X_DLLEXPORT X_USER
#define X_DLLEXPORT_C(type) type X_USER
#define X_DLLEXPORT_H(type) type X_USER
#endif				/* Gestion dllexport */



/**************************** DEFINITION DES TYPES *********************************/
/* Definis pour l'instant uniquement en mode standard 32 bits. */

#if defined(X_BORL_32) || defined(X_MSC_32) \
||  defined(X_ZTC_32) || defined(X_WATC_32) \
||  defined(X_UNIX_32) || defined(X_UNIX_64) || defined(X_OS2_32)
#define X_CHAR      char
#define X_INT16     short int
#define X_INT32     int
#define X_FLOAT32   float
#define X_FLOAT64   double
#elif defined(X_MSC_16)
#define X_CHAR      char
#define X_INT16     short int
#define X_INT32     long
#define X_FLOAT32   float
#define X_FLOAT64   double
#elif defined(X_TMS_C30)
#define X_CHAR      char
#define X_INT16     short int
#define X_INT32     long
#define X_FLOAT32   float

#else
#error "Types generaux non definis ou environnement inconnu."
#endif				/* Definition des types. */


/*************************** Types logiques. **************************************/
typedef enum {
xFalse = 0, xTrue = 1} xBoolean;

typedef enum {
xError = -1, xOK = 1} xStatus;

#endif				/* _X_DEFAUT_H_ */
