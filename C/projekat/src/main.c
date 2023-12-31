/** @file main.c
 *
 * 	@brief Source file implementing and applying audio effect algorithms
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "compound_signal.h"
#include "filter.h"
#include "equ_filters.h"

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

// Index for memory allocation
int index = 0;

/*
 * @brief Applies 5 band equalization to the passed audio signal
 * @param[in] signal Pointer to the audio signal data
 * @param[out] output Pointer to the output array initialized to all zeros
 * @param[in] temp Temporary array for storing filtered intermediate results
 * with len + NUM_TAPS - 1 size
 * @param[in] len Size of input and output signal array
 * */
void equalize(float* restrict signal, float* output, float* temp, int len);

/*
 * @brief Applies the wah-wah effect on the given audio signal
 * @param[in] signal Pointer to the audio signal data
 * @param[out] output Pointer to the output array initialized to all zeros
 * @param[in] len Size of input and output signal array
 * */
void wah_wah(float* restrict signal, float* output, int len);

// Reserve ~500kB of SRAM memory for signals
#pragma section("seg_sram")
static char sram_heap[500000];

int main(int argc, char *argv[])
{
	FILE* fp = NULL;
	int uid = 999;
	float* signal = NULL;
	float* result = NULL;
	float* temp = NULL;
	index = heap_install(sram_heap, sizeof(sram_heap), uid);
	if (index < 0)
	{
		printf("Heap installation failed\n");
		return -1;
	}
	signal = (float *)heap_malloc(index, LEN*sizeof(float));
	if(signal == NULL){
		printf("Memory allocation failed (signal)\n");
		return -1;
	}
	result = (float *)heap_malloc(index, LEN*sizeof(float));
	if(result == NULL){
			printf("Memory allocation failed (eq_signal)\n");
			return -1;
	}
	for(int i = 0; i < LEN; i++)
		signal[i] = compound_signal[i];
	for(int i = 0; i < LEN; i++)
		result[i] = 0.0;
	wah_wah(signal, result, LEN);
	/* EQUALIZATOR
	temp = (float *)heap_malloc(index, (LEN+NUM_TAPS-1)*sizeof(float));
	if(temp == NULL){
		printf("Memory allocation failed (temp)\n");
		return -1;
	}
	equalize(signal, result, temp, LEN);
	heap_free(index, temp);
	*/
	fp = fopen("../output_files/wah_signal.txt", "w");
	if(fp == NULL){
		printf("File opening failed\n");
		return -1;
	}
	for(int i = 0; i < LEN; i++)
		fprintf(fp,"%f\n", result[i]);
	fclose(fp);
	heap_free(index, signal);
	heap_free(index, result);
	return 0;
}

void equalize(float* restrict signal, float* output, float* temp, int len)
{
	// Calculate FIR filter delay
	int delay = (NUM_TAPS - 1) / 2;
	// Convert db gain into unitless gain
	float gain1, gain2, gain3, gain4, gain5;
	gain1 = pow(10, GAIN1/20);
	gain2 = pow(10, GAIN2/20);
	gain3 = pow(10, GAIN3/20);
	gain4 = pow(10, GAIN4/20);
	gain5 = pow(10, GAIN5/20);
	// Clear temp register for every filtering
	for(int i = 0; i < (len + NUM_TAPS - 1); i++)
		temp[i] = 0.0;
	convolve(signal, len, FILTER_1, NUM_TAPS, temp);
	for(int i = 0; i < len; i++)
		output[i] += gain1 * temp[i + delay];
	for(int i = 0; i < (len + NUM_TAPS - 1); i++)
		temp[i] = 0.0;
	convolve(signal, len, FILTER_2, NUM_TAPS, temp);
	for(int i = 0; i < len; i++)
		output[i] += gain2 * temp[i + delay];
	for(int i = 0; i < (len + NUM_TAPS - 1); i++)
		temp[i] = 0.0;
	convolve(signal, len, FILTER_3, NUM_TAPS, temp);
	for(int i = 0; i < len; i++)
		output[i] += gain3 * temp[i + delay];
	for(int i = 0; i < (len + NUM_TAPS - 1); i++)
		temp[i] = 0.0;
	convolve(signal, len, FILTER_4, NUM_TAPS, temp);
	for(int i = 0; i < len; i++)
		output[i] += gain4 * temp[i + delay];
	for(int i = 0; i < (len + NUM_TAPS - 1); i++)
		temp[i] = 0.0;
	convolve(signal, len, FILTER_5, NUM_TAPS, temp);
	for(int i = 0; i < len; i++)
		output[i] += gain5 * temp[i + delay];
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
	if(cutoff_freq == NULL){
		printf("Memory allocation failed (cutoff_freq)\n");
		return;
	}
	cutoff_freq[0] = f_low;
	for(; i < max_iter; i++){
		cutoff_freq[i] = cutoff_freq[i-1] + inc;
	}
	while(i < len){
		for(int j = (max_iter - 1); j >= 0; j--){
			cutoff_freq[i] = cutoff_freq[j];
			i++;
		}
		for(int j = 0; j < max_iter; j++){
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
	for(i = 1; i < len; i++){
		yh = signal[i] - yl - Q1 * output[i-1];
		output[i] = F1 * yh + output[i-1];
		yl = F1 * output[i] + yl;
		F1 = 2 * sinf(f1_param * cutoff_freq[i]);
	}
	heap_free(index, cutoff_freq);
}
