#ifndef _ESTI_DENOISE_
#define _ESTI_DENOISE_

#pragma warning(disable:4096)
#pragma warning(disable:4244)

#ifdef __cplusplus
extern "C" {
#endif

	int etsi_denoise( short* p_data, short* p_denoised, long i_frame );
	int etsi_denoise_synchronization( short* p_data, short* p_denoised, long i_frame );
	
	int etsi_denoise_16k( short* p_data, short* p_denoised, long i_frame );
	int etsi_denoise_16k_synchronization( short* p_data, short* p_denoised, long i_frame );

#ifdef __cplusplus
}
#endif

#endif
