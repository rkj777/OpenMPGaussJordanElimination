/*
* Authors: Eric Smith, Rajan Jassal
* February 26th, 2016
*
* Compile with:
* gcc -Wall -fopenmp Main.c Lab3IO.c -o main
*
*/

#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include "Lab3IO.h"
#include "timer.h"
		
void pivot(double** A, int k, int n);
void swap(double** A, int k, int max_row);

int main(int argc, char* argv[]){
	int thread_count;
	double ** A;
	int n, k, i, j;
	double temp;
	
	double start, end;
	if (argc < 2) {
        printf("Please indicate the number of threads!\n");
        return 1;
	}
	
	thread_count = strtol(argv[1], NULL, 10);
	Lab3LoadInput(&A, &n);

	//Gaussian elimination
	GET_TIME(start);
	
	for(k=0; k < n-1; k++){
		pivot(A,k,n);
#		pragma omp parallel for num_threads(thread_count) private(temp,i,j) 
		for(i=k+1; i < n; i++){
			temp = A[i][k]/A[k][k];
			for(j=k;j<=n;j++){
				A[i][j] = A[i][j] - ( temp * A[k][j]);
			}			
		}
				
	}

	//Jordan elmination	
	for(k=n-1;k>0;k--){
#		pragma omp parallel for num_threads(thread_count) private(i)	
		for(i=0; i<k; i++){
		  A[i][n] = A[i][n] - ((A[i][k]/A[k][k]) * A[k][n]);
		  A[i][k] = 0;
		}
	}

	//Find Solutions
	double* result = malloc(n * sizeof(double));
#	pragma omp parallel for num_threads(thread_count) private(i)	
	for(i=0; i<n;i++){	
		result[i] = A[i][n]/A[i][i];
	}
	
	GET_TIME(end);
	printf("Running-time: %lf\n", end-start);
	Lab3SaveOutput(result,n,end-start);

	free(result);
	DestroyMat(A,n);

	return 0;	
}

void pivot(double** A, int k, int n){
	
	int i;
	int max = abs(A[k][k]);
	int max_row = k;
	//Put in abs
	for(i=k+1; i <= n-1; i++){
		if(max < abs(A[i][k])){
			max = abs(A[i][k]);
			max_row = i;
		}
	
	}
	swap(A,k, max_row);
}

void swap(double** A, int k, int max_row){
	if(max_row != k){
		double* temp = A[max_row];
		A[max_row] = A[k];
		A[k] = temp;
	}
}
