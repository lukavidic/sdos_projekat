#include "string.h"
#include "compound_signal.h"
#include "equalizer.h"
#include <stdio.h>
#include <stdlib.h>

// Reserve ~500kB of SRAM memory for signals
#pragma section("seg_sram")
static char sram_heap[500000];

int main(int argc, char *argv[])
{
	printf("START\n");
	FILE* fp = NULL;
	int index, uid = 999;
	float* signal = NULL;
	float* eq_signal = NULL;
	float* temp = NULL;
	printf("PRE HEAP INSTALLATION\n");
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
	printf("POST HEAP INSTALLATION\n");
	for(int i = 0; i < LEN; i++)
		signal[i] = compound_signal[i];
	for(int i = 0; i < LEN; i++)
		eq_signal[i] = 0.0;
	equalize(signal, eq_signal, temp, LEN);
	printf("FILE OPENING\n");
	fp = fopen("../output_files/output.txt", "w");
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
	printf("END\n");
	return 0;
}

