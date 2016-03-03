// aurora_speech_enhancement.cpp : 定义控制台应用程序的入口点。

#include <Windows.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <time.h>
#include <ipps.h>
#include "aurora/aurora_include.h"
#include "etsi_denoise/AdvFrontEnd.h"

using namespace std;
using namespace aurora;

HANDLE ghMutex;
CRITICAL_SECTION g_section;
int32s g_current = 0;
string src_root_pattern;
string des_root_pattern;
string ext_pattern;

int32s batch_length = 8000 * 60 * 3; // 3 min
int32s etsi_offset_points = 320;

void BufferDenoise(short *src_wave,short *des_wave,
	int32s samples,int32s batch_length, short *work)
{
	int batch_times = samples / batch_length;
	int batch_left = samples % batch_length;

	if (batch_times == 0)
	{
		etsi_denoise(
			src_wave,
			work,
			samples);
		ippsCopy_16s(
			work + etsi_offset_points,
			des_wave,
			samples - etsi_offset_points);
		ippsZero_16s(
			des_wave + samples - etsi_offset_points, 
			etsi_offset_points);
	}
	else
	{
		etsi_denoise(
			src_wave,
			work,
			batch_length);
		ippsCopy_16s(work + etsi_offset_points,
			des_wave,
			batch_length - etsi_offset_points);

		for (int i=1;i<batch_times;i++)
		{
			etsi_denoise(
				src_wave + i * batch_length - etsi_offset_points,
				work, 
				batch_length + etsi_offset_points);
			ippsCopy_16s(
				work + etsi_offset_points, 
				des_wave + i*batch_length - etsi_offset_points, 
				batch_length);
		}

		etsi_denoise(
			src_wave + batch_times * batch_length - etsi_offset_points,
			work, 
			batch_left + etsi_offset_points);
		ippsCopy_16s(
			work + etsi_offset_points,
			des_wave + batch_times * batch_length - etsi_offset_points, 
			batch_left);
		ippsZero_16s(des_wave+samples-etsi_offset_points,
			etsi_offset_points);
	}

	return;
}

void WaveConvert(LPVOID pParam) 
{
	vector<string> *file_point = (vector<string>*) pParam;
	int32s file_count = file_point->size();
	int32s current;


	int8u_vector buf_raw;
	int16s_vector batch_buffer;
	
	Wav buf_wav;
	Wav buf_wav_8k16b;
	Wav buf_denoise;

	float_vector buf_flt_wav;
	float_vector buf_flt_filter;

	int32s max_wave_size = 16000 * 2 * 7200;

	buf_raw.construct(max_wave_size);
	batch_buffer.construct(batch_length + 320);

	buf_wav.construct(max_wave_size);
	buf_wav_8k16b.construct(max_wave_size);
	buf_denoise.construct(max_wave_size);
	
	buf_flt_wav.construct(max_wave_size);
	buf_flt_filter.construct(max_wave_size);

	while(1)
	{
		EnterCriticalSection(&g_section);
		current = g_current;
		g_current++;
		LeaveCriticalSection(&g_section);

		if(g_current > file_count)
		{
			break;
		}

		string *iter = &(*file_point)[current];
		string src_file = *iter;

		printf("== thread %d process %s ..\n", 
			GetCurrentThreadId(),src_file.c_str());

		string des_file = src_file;
		string::size_type pos = des_file.find_first_of(src_root_pattern);
		string::size_type length = src_root_pattern.size();
		des_file.replace(
			des_file.begin() + pos, des_file.begin() + pos + length,
			des_root_pattern.begin(),des_root_pattern.end());

		// check ext
		char tmp_driver[_MAX_DRIVE];
		char tmp_dir[_MAX_DIR];
		char tmp_filename[_MAX_FNAME];
		char tmp_ext[_MAX_EXT];
		_splitpath( src_file.c_str(),tmp_driver,tmp_dir,tmp_filename,tmp_ext);
		if (strcmp(tmp_ext,ext_pattern.c_str()) != 0)
		{
			CopyFileA(src_file.c_str(),des_file.c_str(),FALSE);
			continue;
		}

		if (read_buffer(src_file.c_str(),buf_raw._data,buf_raw.get_mem_len(),buf_raw._len))
		{
			printf("read raw buffer %s failed.\n",src_file.c_str());
			continue;
		}

		if (read_wave(buf_raw._data,buf_raw._len,buf_wav))
		{
			printf("parse wave %s failed.\n",src_file.c_str());
			continue;
		}

		if (speech_coding_standard(buf_wav,buf_wav_8k16b))
		{
			printf("== convert %s wave format failed.\n",src_file.c_str());
			continue;
		}
		if (buf_wav_8k16b._data_size < etsi_offset_points)
		{
			continue;
		}

		
		int32s samples = get_wave_samples(buf_wav_8k16b.wt,buf_wav_8k16b._data_size);
		double energy = 0.0f;
		short *p_data = (short*)buf_wav_8k16b._data;
		for (int t=0;t<samples;t++)
		{
			energy += p_data[t] * p_data[t]; 
		}
		energy /= samples; energy = sqrt(energy);
		double factor = 3000.0f / energy;
		printf("factor %f\n",factor);
		for (int t=0;t<samples;t++)
		{
			p_data[t] = (short)((float)p_data[t] * factor);
		}


		BufferDenoise(
			(short*)buf_wav_8k16b._data,
			(short*)buf_denoise._data,
			get_wave_samples(buf_wav_8k16b.wt,buf_wav_8k16b._data_size),
			batch_buffer.get_mem_len(), 
			batch_buffer._data);

		ippsConvert_16s32f(
			(short*) buf_denoise._data,
			buf_flt_wav._data,
			get_wave_samples(buf_wav_8k16b.wt,buf_wav_8k16b._data_size));

		IppsFIRState_32f *fir_state;
		ippsFIRInitAlloc_32f(
			&fir_state,
			filter_200_to_3600_taps,
			filter_200_to_3600_taps_length,
			NULL);
		ippsFIR_32f(
			buf_flt_wav._data,
			buf_flt_filter._data,
			get_wave_samples(buf_wav_8k16b.wt,buf_wav_8k16b._data_size),
			fir_state);
		ippsFIRFree_32f(fir_state);

		for (int t=0;t<get_wave_samples(buf_wav_8k16b.wt,buf_wav_8k16b._data_size);t++)
		{
			((short*)buf_denoise._data)[t] = (short) buf_flt_filter._data[t]; 
		}
		
		if ( write_wave_sclin8k16b(
			des_file.c_str(),
			(int16s*)buf_denoise._data, 
			get_wave_samples(buf_wav_8k16b.wt,buf_wav_8k16b._data_size),
			buf_wav._data,
			buf_wav.get_mem_size()) )
		{
			printf("== write %s wave failed.\n",des_file.c_str());
			continue;
		}
	}

	return;
}


int main(int argc, char  *argv[])
{
	InitializeCriticalSection(&g_section);

	if (argc != 5)
	{
		printf("%s <thread_num> <src_dir> <des_dir> <ext_format>\n",argv[0]);
		return -1;
	}
	printf("== Aurora ==\n");
	printf("== cmd: %s <%s> <%s> <%s> <%s>\n",
		argv[0],argv[1],argv[2],argv[3],argv[4]);
	printf("== begin process .. \n");

	clock_t start, end;
	start = clock();

	int thread_number = atoi(argv[1]);
	string src_root_dir(argv[2]);
	string des_root_dir(argv[3]);
	string ext_format(argv[4]);

	if (ext_format[0] == '.')
	{
		ext_pattern = ext_format;
	}
	else
	{
		ext_pattern = "." + ext_format;
	}

	// generate list	
	vector<string> dirs;
	walk_dir_to_find_dir(src_root_dir.c_str(),dirs);

	vector<string> files;
	walk_dir_to_find_file(src_root_dir.c_str(),files,NULL);

	// remove last 
	char tmp_src_root[MAX_PATH];
	char tmp_des_root[MAX_PATH];
	strcpy(tmp_src_root,src_root_dir.c_str());
	strcpy(tmp_des_root,des_root_dir.c_str());
	if (tmp_src_root[strlen(tmp_src_root)-1] == '\\' 
		|| tmp_src_root[strlen(tmp_src_root)-1] == '/' ){
			tmp_src_root[strlen(tmp_src_root)-1] = '\0';
	}
	if (tmp_des_root[strlen(tmp_des_root)-1] == '\\' 
		|| tmp_des_root[strlen(tmp_des_root)-1] == '/' ){
			tmp_des_root[strlen(tmp_des_root)-1] = '\0';
	}
	src_root_pattern.assign(tmp_src_root);
	des_root_pattern.assign(tmp_des_root);

	// prepare des dirs
	ghMutex = CreateMutex(NULL,FALSE,L"aurora");
	if (ghMutex == NULL){
		printf("CreateMutex error: %d\n", GetLastError());
		return -1;
	}
	for(vector<string>::iterator iter = dirs.begin();
		iter != dirs.end(); iter ++)
	{
		string src_dir = (*iter);
		string des_dir;

		string::size_type pos = src_dir.find_first_of(src_root_pattern);
		string::size_type length = src_root_pattern.size();
		des_dir = src_dir.replace(
			src_dir.begin() + pos, src_dir.begin() + pos + length,
			des_root_pattern.begin(),des_root_pattern.end());
		// printf("|%s|%s|\n",src_dir.c_str(),des_dir.c_str());
		check_and_make_dir_directory(des_dir.c_str());
	}
	CloseHandle(ghMutex);

	// process, begin multiple process
	g_current = 0;
	HANDLE* hThreads = new HANDLE[thread_number];

	for(int i = 0; i < thread_number; ++i)
	{
		DWORD dwThreadId;
		hThreads[i] = CreateThread
			(
			NULL,				     					//Choose default security
			0,											//Default stack size
			(LPTHREAD_START_ROUTINE) &WaveConvert,	    //Routine to execute
			(LPVOID) &files,						    //Thread parameter
			0,											//Immediately run the thread
			&dwThreadId									//Thread Id
			);

	}
	WaitForMultipleObjects(thread_number, hThreads, TRUE, INFINITE);
	delete [] hThreads;

	end = clock();
	printf("== processed %d files, using %f second.\n", 
		files.size(), ((float)end - start) / CLOCKS_PER_SEC);

	// end multiple process
	DeleteCriticalSection(&g_section);
	return 0;
}

