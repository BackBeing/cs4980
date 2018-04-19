#include <iostream>
#include <complex>
#include "mmsystem.h"


using namespace std;


#define M_PI 3.1415926535897932384


//////////////////////////////////////////////////////////
// Sampling frequency = 48000 Hz
// Sampling Resolution = 16 bits   65 536 = 2^16
// 		so for each real time sec of sound my program will 
// send ~2/3 of the data to be transformed, so it will have 
// 17536 * (n - 1) overlapped data for n secs sound.
// I will basically add the related overlapped sound together 
// and /2. 

int main(int argc, char *argv[])
{
	
	int N = 65 536;
	
	///////////////
	// read from some sound file, store in K[ Ntrue ]
	// e.g. K[0] = 16 bits
	int Ntrue = 0 ;
	//TODO: get Ntrue here:// high
	int sampleN = 0;
	while ( Ntrue >= 65 536 ){
		Ntrue -= 65 536;
		sampleN ++;
	}
	
	if( Ntrue > 0 ){
		sampleN ++;
	}
	
	double ** Kgrps; // a group of sample groups
	Kgrps = (double **)malloc( sampleN * sizeof(double *));
	for (int i = 0; i< sampleN; i++){
			Kgrps[i] = (complex<double> *)malloc(N * sizeof(complex<double>));
	}
	
	double * K; // a sample group
	K = (double *)malloc(N * sizeof(double));
	
	//TODO: give K, KK value here:// high
	
	//set up W, which can be used for all sample groups
	//TODO: build fixed W on board:// low
	W = (complex<double> *)malloc(N / 2 * sizeof(complex<double>));
    W[1] = polar(1., -2. * M_PI / N);
    W[0] = 1;
    for(int i = 2; i < N / 2; i++)
		W[i] = pow(W[1], i);
	
	
	// format the Kgrps
	complex<double> ** tKgrps; // a group of sample groups' complex
	tKgrps = (complex<double> **)malloc( sampleN * sizeof(complex<double> *));
	for (int i = 0; i< sampleN; i++){
			tKgrps[i] = (complex<double> *)malloc(N * sizeof(complex<double>));
			for (int j = 0; j< N; j++){
				tKgrps[i][j] = Kgrps[i][j];
			}
	}
	
	
	//TODO: pass tKgrps real, tKgrps imag, W real, W imag, and N(maybe), sampleN to board for FFT convert:// high
	
	//TODO: get tKgrpsNew back from board:// high
	
	
	
	
	//TODO: do high/low filter and sound effect things on tKgrpsNew:// high 
	
	
	
	
	
	//TODO: pass tKgrpsNew real, tKgrpsNew imag, W real, W imag, and N(maybe), sampleN to board for IFFT convert:// high
	
	//TODO: get tKgrpsNewNew back from board:// high
	
	
	
	
	//TODO: rebuild real time sound files on the tKgrpsNewNew
	
	// free
	for (int i = 0; i< sampleN; i++){
			free(Kgrps[i]);
	}
	for (int i = 0; i< sampleN; i++){
			free(tKgrps[i]);
	}
	
	
}