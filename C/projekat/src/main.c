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

// TODO: Try to implement parameter configuration with ADSP buttons

// Gain parameters for each band in dB (default values)
float GAIN1 = -12.0;
float GAIN2 = 15.0;
float GAIN3 = 12.0;
float GAIN4 = -9.0;
float GAIN5 = -15.0;

/*
 * @brief Applies 5 band equalization to the passed audio signal
 * @param[in] signal Pointer to the audio signal data
 * @param[out] output Pointer to the output array initialized to all zeros
 * @param[in] temp Temporary array for storing filtered intermediate results
 * with len + NUM_TAPS - 1 size
 * @param[in] len Size of input signal array
 * */
void equalize(float* restrict signal, float* output, float* temp, int len);

// Reserve ~500kB of SRAM memory for signals
#pragma section("seg_sram")
static char sram_heap[500000];

int main(int argc, char *argv[])
{
	FILE* fp = NULL;
	int index, uid = 999;
	float* signal = NULL;
	float* eq_signal = NULL;
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
	eq_signal = (float *)heap_malloc(index, LEN*sizeof(float));
	if(eq_signal == NULL){
			printf("Memory allocation failed (eq_signal)\n");
			return -1;
	}
	temp = (float *)heap_malloc(index, (LEN+NUM_TAPS-1)*sizeof(float));
	if(temp == NULL){
		printf("Memory allocation failed (temp)\n");
		return -1;
	}
	for(int i = 0; i < LEN; i++)
		signal[i] = compound_signal[i];
	for(int i = 0; i < LEN; i++)
		eq_signal[i] = 0.0;
	equalize(signal, eq_signal, temp, LEN);
	fp = fopen("../output_files/equ_signal.txt", "w");
	if(fp == NULL){
		printf("File opening failed\n");
		return -1;
	}
	for(int i = 0; i < LEN; i++)
		fprintf(fp,"%f\n", eq_signal[i]);
	fclose(fp);
	heap_free(index, signal);
	heap_free(index, eq_signal);
	heap_free(index, temp);
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

