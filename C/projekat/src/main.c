/** @file main.c
 *
 * 	@brief Source file implementing and applying audio effects algorithms
 */

#include <sys/platform.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cycle_count.h>
#include "compound_signal.h"
#include "filter.h"
#include "equ_filters.h"
#include "adi_initialize.h"

//#define PRINT_OUTPUT
//#define PROFILING

#define PI 3.141592

// TODO: Try to implement parameter configuration with ADSP buttons

// Gain parameters for each band in dB (default values)
float GAIN1 = -12.0;
float GAIN2 = 15.0;
float GAIN3 = 12.0;
float GAIN4 = -9.0;
float GAIN5 = -15.0;

// Parameters for wah-wah filter (default values)
float f_low = 50.0;
float f_high = 1500.0;
float wah_freq = 1800.0;
float damp = 0.1;

// Parameters for flanger (default values)
float delay = 0.001;
float f_osc = 1.0;

// Parameters for tremolo (default values)
float alpha = 0.5;
float f_mod = 1;

// Index for memory allocation
int index = 0;

#ifdef PROFILING
	// Profiling variables
	cycle_t start;
	cycle_t end;
#endif

/*
 * @brief Applies 5 band equalization to the passed audio signal
 * @param[in] signal Pointer to the audio signal data
 * @param[out] output Pointer to the output array initialized to all zeros
 * @param[in] len Size of input and output signal array
 * */
void equalize(float* restrict signal, float* output, int len);

/*
 * @brief Applies the wah-wah effect on the given audio signal
 * @param[in] signal Pointer to the audio signal data
 * @param[out] output Pointer to the output array initialized to all zeros
 * @param[in] len Size of input and output signal array
 * */
void wah_wah(float* restrict signal, float* output, int len);

/*
 * @brief Applies the flanger effect on the given audio signal
 * @param[in] signal Pointer to the audio signal data
 * @param[out] output Pointer to the output array initialized to all zeros
 * @param[in] len Size of input and output signal array
 * */
void flanger(float* restrict signal, float* output, int len);

/*
 * @brief Applies the flanger effect on the given audio signal
 * @param[in] signal Pointer to the audio signal data
 * @param[out] output Pointer to the output array initialized to all zeros
 * @param[in] len Size of input and output signal array
 * */
void tremolo(float* restrict signal, float* output, int len);


// Reserve 500kB of SDRAM memory for signals
#pragma section("seg_sdram")
static char heap_mem[512000];

// UNCOMMENT IF YOU WANT TO USE SRAM MEMORY (WARNING: SIMD FUNCTIONS GIVE WRONG RESULTS)
/*
// Reserve 500kB of SRAM memory for signals
#pragma section("seg_sram")
static char heap_mem[512000];
*/

float dm state[NUM_TAPS + 1];

int main(int argc, char *argv[])
{
	adi_initComponents();
	FILE* fp = NULL;
	int uid = 999;
	float* signal = NULL;
	float* result = NULL;
	index = heap_install(heap_mem, sizeof(heap_mem), uid);
	if (index < 0)
	{
		printf("Heap installation failed\n");
		return -1;
	}
	signal = (float *)heap_malloc(index, LEN*sizeof(float));
	if(signal == NULL)
	{
		printf("Memory allocation failed (signal)\n");
		return -1;
	}
	result = (float *)heap_malloc(index, LEN*sizeof(float));
	if(result == NULL)
	{
			printf("Memory allocation failed (eq_signal)\n");
			return -1;
	}
	for(int i = 0; i < LEN; i++)
		signal[i] = compound_signal[i];
	for(int i = 0; i < LEN; i++)
		result[i] = 0.0;
	#ifdef PROFILING
		START_CYCLE_COUNT(start);
	#endif
	//wah_wah(signal, result, LEN);
	//equalize(signal, result, LEN);
	//flanger(signal, result, LEN);
	//tremolo(signal, result, LEN);
	#ifdef PROFILING
		STOP_CYCLE_COUNT(end, start);
		PRINT_CYCLES("Broj ciklusa: ", end);
	#endif
	#ifdef PRINT_OUTPUT
		fp = fopen("../output_files/tmp_signal.txt", "w");
		if(fp == NULL)
		{
			printf("File opening failed\n");
			return -1;
		}
		for(int i = 0; i < LEN; i++)
			fprintf(fp,"%f\n", result[i]);
		fclose(fp);
	#endif
	heap_free(index, signal);
	heap_free(index, result);
	return 0;
}

/* SLOW IMPLEMENTATION
void equalize(float* restrict signal, float* output, int len)
{
	// Calculate FIR filter delay
	int delay = (NUM_TAPS - 1) / 2;
	float gain1, gain2, gain3, gain4, gain5;
	float* temp = NULL;
	temp = (float *)heap_malloc(index, (len+NUM_TAPS-1)*sizeof(float));
	if(temp == NULL)
	{
		printf("Memory allocation failed (temp)\n");
		return;
	}
	// Convert db gain into unitless gain
	gain1 = pow(10, GAIN1/20);
	gain2 = pow(10, GAIN2/20);
	gain3 = pow(10, GAIN3/20);
	gain4 = pow(10, GAIN4/20);
	gain5 = pow(10, GAIN5/20);
	// Clear temp register for every filtering
	//#pragma SIMD_for
	for(int i = 0; i < (len + NUM_TAPS - 1); i++)
		temp[i] = 0.0;
	convolve(signal, len, FILTER_1, NUM_TAPS, temp);
	//#pragma SIMD_for
	for(int i = 0; i < len; i++)
		output[i] += gain1 * temp[i + delay];
	//#pragma SIMD_for
	for(int i = 0; i < (len + NUM_TAPS - 1); i++)
		temp[i] = 0.0;
	convolve(signal, len, FILTER_2, NUM_TAPS, temp);
	#pragma SIMD_for
	for(int i = 0; i < len; i++)
		output[i] += gain2 * temp[i + delay];
	//#pragma SIMD_for
	for(int i = 0; i < (len + NUM_TAPS - 1); i++)
		temp[i] = 0.0;
	convolve(signal, len, FILTER_3, NUM_TAPS, temp);
	//#pragma SIMD_for
	for(int i = 0; i < len; i++)
		output[i] += gain3 * temp[i + delay];
	//#pragma SIMD_for
	for(int i = 0; i < (len + NUM_TAPS - 1); i++)
		temp[i] = 0.0;
	convolve(signal, len, FILTER_4, NUM_TAPS, temp);
	//#pragma SIMD_for
	for(int i = 0; i < len; i++)
		output[i] += gain4 * temp[i + delay];
	//#pragma SIMD_for
	for(int i = 0; i < (len + NUM_TAPS - 1); i++)
		temp[i] = 0.0;
	convolve(signal, len, FILTER_5, NUM_TAPS, temp);
	//#pragma SIMD_for
	for(int i = 0; i < len; i++)
		output[i] += gain5 * temp[i + delay];
	heap_free(index, temp);
}
*/


void equalize(float* restrict signal, float* output, int len)
{
	float gain1, gain2, gain3, gain4, gain5;
	float* temp = NULL;
	temp = (float *)heap_malloc(index, len*sizeof(float));
	if(temp == NULL)
	{
		printf("Memory allocation failed (temp)\n");
		return;
	}
	// Convert db gain into unitless gain
	gain1 = pow(10, GAIN1/20);
	gain2 = pow(10, GAIN2/20);
	gain3 = pow(10, GAIN3/20);
	gain4 = pow(10, GAIN4/20);
	gain5 = pow(10, GAIN5/20);
	//#pragma SIMD_for
	for(int i = 0; i < (NUM_TAPS + 1); i++)
		state[i] = 0.0;
	//#pragma SIMD_for
	for(int i = 0; i < len; i++)
		temp[i] = 0.0;
	firf(signal, temp, FILTER_1, state, len, NUM_TAPS);
	//#pragma SIMD_for
	for(int i = 0; i < len; i++)
		output[i] += gain1 * temp[i];
	//#pragma SIMD_for
	for(int i = 0; i < (NUM_TAPS + 1); i++)
		state[i] = 0.0;
	//#pragma SIMD_for
	for(int i = 0; i < len; i++)
		temp[i] = 0.0;
	firf(signal, temp, FILTER_2, state, len, NUM_TAPS);
	//#pragma SIMD_for
	for(int i = 0; i < len; i++)
		output[i] += gain2 * temp[i];
	//#pragma SIMD_for
	for(int i = 0; i < (NUM_TAPS + 1); i++)
		state[i] = 0.0;
	//#pragma SIMD_for
	for(int i = 0; i < len; i++)
		temp[i] = 0.0;
	firf(signal, temp, FILTER_3, state, len, NUM_TAPS);
	//#pragma SIMD_for
	for(int i = 0; i < len; i++)
		output[i] += gain3 * temp[i];
	//#pragma SIMD_for
	for(int i = 0; i < (NUM_TAPS + 1); i++)
		state[i] = 0.0;
	//#pragma SIMD_for
	for(int i = 0; i < len; i++)
		temp[i] = 0.0;
	firf(signal, temp, FILTER_4, state, len, NUM_TAPS);
	//#pragma SIMD_for
	for(int i = 0; i < len; i++)
		output[i] += gain4 * temp[i];
	//#pragma SIMD_for
	for(int i = 0; i < (NUM_TAPS + 1); i++)
		state[i] = 0.0;
	//#pragma SIMD_for
	for(int i = 0; i < len; i++)
		temp[i] = 0.0;
	firf(signal, temp, FILTER_5, state, len, NUM_TAPS);
	//#pragma SIMD_for
	for(int i = 0; i < len; i++)
		output[i] += gain5 * temp[i];
	heap_free(index, temp);
}


void wah_wah(float* restrict signal, float* output, int len)
{
	/*
	 * Unlike in python where we used arrays, here we have limited memory so we will
	 * keep only 1 variable for low pass (yl) and high pass (yh) values and place
	 * band pass (yb) values directly into the output :)
	 */
	int i = 1, max_iter;
	float yh, yl, inc, Q1, F1, f1_param;
	float* cutoff_freq = NULL;
	// Generate triangular modulator
	inc = wah_freq / SAMPLE_FREQ;
	max_iter = (int)((f_high - f_low) / inc) + 1;
	/* We will allocate len+2*max_iter elements for cutoff_freq array in order
	 * to avoid if and break statements inside loops which create modulator
	 */
	cutoff_freq = (float *)heap_malloc(index, (len + 2*max_iter)*sizeof(float));
	if(cutoff_freq == NULL)
	{
		printf("Memory allocation failed (cutoff_freq)\n");
		return;
	}
	cutoff_freq[0] = f_low;
	for(; i < max_iter; i++)
	{
		cutoff_freq[i] = cutoff_freq[i-1] + inc;
	}
	while(i < len)
	{
		for(int j = (max_iter - 1); j >= 0; j--)
		{
			cutoff_freq[i] = cutoff_freq[j];
			i++;
		}
		for(int j = 0; j < max_iter; j++)
		{
			cutoff_freq[i] = cutoff_freq[j];
			i++;
		}
	}
	// Implement algorithm
	f1_param = PI / SAMPLE_FREQ; // Calculate once to avoid multiple float divisions in loop
	F1 = 2 * sinf(f1_param * cutoff_freq[0]);
	Q1 = 2 * damp;
	yh = signal[0];
	output[0] = F1 * yh;
	yl = F1 * output[0];
	for(i = 1; i < len; i++)
	{
		yh = signal[i] - yl - Q1 * output[i-1];
		output[i] = F1 * yh + output[i-1];
		yl = F1 * output[i] + yl;
		F1 = 2 * sinf(f1_param * cutoff_freq[i]);
	}
	heap_free(index, cutoff_freq);
}

void flanger(float* restrict signal, float* output, int len)
{
	float* delay_line = NULL;
	float beta, frac, w_osc, R;
	int N, dl_len;
	w_osc = 2 * PI * f_osc / SAMPLE_FREQ;
	R = floorf(delay * SAMPLE_FREQ + 0.5);
	dl_len = 2 * ((int)R + 1);
	delay_line = (float*)heap_malloc(index, dl_len);
	if(delay_line == NULL)
	{
		printf("Memory allocation failed (delay_line)\n");
		return;
	}
	for(int i = 0; i < dl_len; i++)
		delay_line[i] = 0.0;
	for(int i = 0; i < len; i++)
	{
		beta = R * (1 + sinf(w_osc*i));
		N = (int)(floorf(beta));
		frac = beta - N;
		for(int j = (dl_len - 1); j > 0; j--)
			delay_line[j] = delay_line[j-1];
		delay_line[0] = signal[i];
		output[i] = delay_line[N+1] * frac + delay_line[N] * (1 - frac) + signal[i];
	}
	heap_free(index, delay_line);
}

void tremolo(float* restrict signal, float* output, int len)
{
	float* mod_signal = NULL;
	float w_mod = 2 * PI * f_mod / SAMPLE_FREQ;
	// Generate modulating signal with given modulation frequency
	mod_signal = (float*)heap_malloc(index, len);
	if(mod_signal == NULL)
	{
		printf("Memory allocation failed (carrier)\n");
		return;
	}
	//#pragma SIMD_for
	for(int i = 0; i < len; i++)
		mod_signal[i] = sinf(w_mod*i);
	// Generate tremolo output
	//#pragma SIMD_for
	for(int i = 0; i < len; i++)
		output[i] = (1 + alpha * mod_signal[i]) * signal[i];
	heap_free(index, mod_signal);
}
