#include <complex>
#include <unistd.h>
#include <stdio.h>
#include <string>
#include <cstring>
#include <stdlib.h>
#include "wave.h"
#define TRUE 1 
#define FALSE 0
#define MAX 65536

using namespace std;


// WAVE header structure

unsigned char buffer4[4];
unsigned char buffer2[2];

char* seconds_to_time(float seconds);
int read_length(int argcc, char **argvv);
int read_array(int argcc, char **argvv, unsigned int & sampleN, unsigned int & channelN, char* ee, unsigned int* Karray);
int read_head(int argcc, char **argvv, char* mm);
int write_array(int argcc, char **argvv, unsigned int Ntruecpy, unsigned int channelN, unsigned int* Karray, char* head, int binc,char* end);
int write_array2(int argcc, char **argvv, unsigned int Ntruecpy, unsigned int channelN, unsigned int* Karray, char* head, int binc,char* end);
int log2(int N);
int check(int n);
int reverse(int N, int n);
void ordina(complex<double>* f1, int N);
void transform(complex<double>* f, int N);
void FFT(complex<double>* f, int N, double d);
void high_low_filter(complex<double>* f, complex<double>* ori, int N, int highfre, int lowfre);


 FILE *ptr;
 char *filename;
 struct HEADER header;

//////////////////////////////////////////////////////////
// Sampling frequency = 44100 Hz						//
// Sampling Resolution = 16 bits   65 536 = 2^16		//
//////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
	
	int N = 65536;
	///////e.g. 2 channel 5 Ntrue 10 data
	/////// 1 2 3 4 5
	///////	6 7 8 9 10
	///////the Karray will be 1 2 3 4 5 6 7 8 9 10
	
	///////////////
	//read from sound file
	//we assume the file will only have two channels!
	//and we will just modify one of the channel.
	unsigned int * Karray;
	
	unsigned int Ntrue = 0;
	unsigned int channelN = 0;

	
	//The wave file is now in an array: Karray!!!!
	//Ntrue is the sample numbers
	char ee [8];
	Karray = (unsigned int *) malloc (read_length(argc, argv) * 2 * sizeof(unsigned int));
	read_array( argc, argv, Ntrue, channelN, ee, Karray );
	char mm [46];

	/**
	for (int i = 0; i< 100; i++){
		printf("%d : %d\n", i, Karray[i]);
	}
	**/
	
	read_head(argc, argv, mm);
	int binc = 2; //bytes in each channel :/assumption
	
	//write to a wave file://
	write_array2( argc, argv, Ntrue, channelN, Karray, mm, binc, ee);	

	int sampleN = 0 ;
	unsigned int Ntruecpy = Ntrue;
	while ( Ntrue >= 65536 ){
		Ntrue -= 65536;
		sampleN ++;
	}
	
	if( Ntrue > 0 ){
		sampleN ++;
	}

/** check point
	printf("sample group numbers (divided by 65536): %d\n", sampleN);
	printf("sample numbers: %d\n", Ntruecpy);
**/
	
	complex<double> ** Kgrps; // a group of sample groups devided by 65536
	Kgrps = (complex<double> **)malloc( sampleN * sizeof(complex<double> *));
	for (int i = 0; i< sampleN; i++){
		Kgrps[i] = (complex<double> *)malloc(N * sizeof(complex<double>));
		
		for(int j = 0; j< N; j++){
			if ( (i * N + j) < Ntruecpy ){
				Kgrps[i][j] = complex<double>(double(1.0 * Karray[ (i * N + j)]),0.0);	
				//printf("expected sample number real: %1f image: 0 \n", real(Kgrps[i][j]));
			}else{
				Kgrps[i][j] = (0.0,0.0);	
			}	
		}
	}

	complex<double> ** Kgrps2; // a group of sample groups devided by 65536
	Kgrps2 = (complex<double> **)malloc( sampleN * sizeof(complex<double> *));
	for (int i = 0; i< sampleN; i++){
		Kgrps2[i] = (complex<double> *)malloc(N * sizeof(complex<double>));
		
		for(int j = 0; j< N; j++){
			if ( (i * N + j) < Ntruecpy ){
				Kgrps2[i][j] = complex<double>(double(1.0 * Karray[ (i * N + j + sampleN)]),0.0);	
				//printf("expected sample number real: %1f image: 0 \n", real(Kgrps[i][j]));
			}else{
				Kgrps2[i][j] = (0.0,0.0);	
			}	
		}
	}
	
/** check point
	for (int i = 0; i< sampleN; i++){
		printf("expected sample number real: %d image: 0 \n", Karray[ (i * N) ]);
		printf("sample number real: %1f  image: %1f \n", real(Kgrps[i][0]), imag(Kgrps[i][0]))	;	
	}
**/	
	
	//set up W, which can be used for all sample groups

	complex<double> * W;
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
	
	complex<double> ** tKgrps2; // a group of sample groups' complex
	tKgrps2 = (complex<double> **)malloc( sampleN * sizeof(complex<double> *));
	for (int i = 0; i< sampleN; i++){
			tKgrps2[i] = (complex<double> *)malloc(N * sizeof(complex<double>));
			for (int j = 0; j< N; j++){
				tKgrps2[i][j] = Kgrps2[i][j];
			}
	}
	printf("Data transformed..\n");
	
	

	//TODO:pass tKgrps real, tKgrps imag, W real, W imag, and N(maybe), sampleN to board for FFT convert:// high
	


	//TODO: get tKgrpsNew back from board:// high

	//here just pass Kgrps to tKgrpsNew
	complex<double> ** tKgrpsNew;
	tKgrpsNew = (complex<double> **)malloc( sampleN * sizeof(complex<double> *));
	for (int i = 0; i< sampleN; i++){
		tKgrpsNew[i] = (complex<double> *)malloc(N * sizeof(complex<double>));
		for (int j = 0; j< N; j++){
			tKgrpsNew[i][j] = Kgrps[i][j];
		}
	}

	for (int i = 0; i< sampleN; i++){
		FFT(tKgrpsNew [i], N, 1);
	}
	
	complex<double> ** tKgrpsNew2;
	tKgrpsNew2 = (complex<double> **)malloc( sampleN * sizeof(complex<double> *));
	for (int i = 0; i< sampleN; i++){
		tKgrpsNew2[i] = (complex<double> *)malloc(N * sizeof(complex<double>));
		for (int j = 0; j< N; j++){
			tKgrpsNew2[i][j] = Kgrps2[i][j];
		}
	}

	for (int i = 0; i< sampleN; i++){
		FFT(tKgrpsNew2 [i], N, 1);
	}
	

	printf("FFT transformed..\n");

	
	


	
	//TODO: do high/low filter and sound effect things on tKgrpsNew:// high 
	/*
	for(int i=0; i<sampleN; i++){	
		high_low_filter(tKgrpsNew[i], Kgrps[i], N, 50000, 20);
	}
	*/
	
	//TODO: pass tKgrpsNew real, tKgrpsNew imag, W real, W imag, and N(maybe), sampleN to board for IFFT convert:// high
	
	//TODO: get tKgrpsNewNew back from board:// high
	
	//here just pass Kgrps to tKgrpsNewNew
	complex<double> ** tKgrpsNewNew;
	tKgrpsNewNew = (complex<double> **)malloc( sampleN * sizeof(complex<double> *));
	for (int i = 0; i< sampleN; i++){
			tKgrpsNewNew[i] = (complex<double> *)malloc(N * sizeof(complex<double>));
			for (int j = 0; j< N; j++){
				tKgrpsNewNew[i][j] = tKgrpsNew[i][j];
			}
	}

	for (int i = 0; i< sampleN; i++){
		FFT(tKgrpsNewNew [i], N, 1);
	}

	complex<double> ** tKgrpsNewNew2;
	tKgrpsNewNew2 = (complex<double> **)malloc( sampleN * sizeof(complex<double> *));
	for (int i = 0; i< sampleN; i++){
			tKgrpsNewNew2[i] = (complex<double> *)malloc(N * sizeof(complex<double>));
			for (int j = 0; j< N; j++){
				tKgrpsNewNew2[i][j] = tKgrpsNew2[i][j];
			}
	}

	for (int i = 0; i< sampleN; i++){
		FFT(tKgrpsNewNew2 [i], N, 1);
	}
	printf("IFFT transformed..\n");
	
	//rebuilding the sound array Karray
	for (int i = 0; i< sampleN; i++){
		for (int j = 0; j< N; j++){
			if ( (i * N + j) < Ntruecpy ){
				Karray[ (i * N + j)] = (int) real(tKgrpsNewNew[i][j]);	
			}else{
				//do nothing
			}	
		}
	}
	for (int i = 0; i< sampleN; i++){
		for (int j = 0; j< N; j++){
			if ( (i * N + j) < Ntruecpy ){
				Karray[ (i * N + j + sampleN)] = (int) real(tKgrpsNewNew2[i][j]);	
			}else{
				//do nothing
			}	
		}
	}
	printf("Sound formatted..\n");
	
	//write to a wave file:// medium
	write_array( argc, argv, Ntruecpy, channelN, Karray, mm, binc, ee);



	// free
	for (int i = 0; i< sampleN; i++){
			free(Kgrps[i]);
	}
	for (int i = 0; i< sampleN; i++){
			free(tKgrps[i]);
	}
	
	
}


int write_array(int argcc, char **argvv, unsigned int Ntruecpy, unsigned int channelN, unsigned int* Karray, char* head, int binc, char* end){
 filename = (char*) malloc(sizeof(char) * 1024);
 if (filename == NULL) {
   printf("Error in malloc\n");
   exit(1);
 }
 

 // get file path
 char cwd[1024];
 if (getcwd(cwd, sizeof(cwd)) != NULL) {
   
	strcpy(filename, cwd);

	// get filename from command line
	if (argcc < 3) {
	  printf("No output wave file specified\n");
	  return 0;
	}
	
	strcat(filename, "/");
	strcat(filename, argvv[2]);
	printf("%s\n", filename);
 }

 // open file
 printf("Writing file..\n");
 ptr = fopen(filename, "wb");
 if (ptr == NULL) {
	printf("Error opening file\n");
	exit(1);
 }

 int write = 0;

 write = fwrite(head, 46, 1, ptr);
//printf("size of head: %d\n",46);
//printf("size of n: %d %d\n",Ntruecpy, channelN);

for (int i =0; i < Ntruecpy; i++) {

	// dump the data read
	unsigned int  xchannels = 0;
	

	for (xchannels = 0; xchannels < channelN; xchannels ++ ) {
		unsigned char data_buffer[binc];
		if (binc == 4) {
			data_buffer[0] = Karray[i + xchannels*Ntruecpy] & 0xFF;
			data_buffer[1] = (Karray[i + xchannels*Ntruecpy] & 0xFF00) >> 8 ;
			data_buffer[2] = (Karray[i + xchannels*Ntruecpy] & 0xFF0000)  >> 16;
			data_buffer[3] = (Karray[i + xchannels*Ntruecpy] & 0xFF000000) >> 24;
			
		}
		else if (binc == 2) {
			data_buffer[0] = Karray[ i + xchannels*Ntruecpy] & 0xFF;
			data_buffer[1] = (Karray[ i + xchannels*Ntruecpy] >> 8) & 0xFF;
		}
		else if (binc == 1) {
			data_buffer[0] = Karray[ i + xchannels*Ntruecpy] ;

		}

		write = fwrite(data_buffer, sizeof(data_buffer), 1, ptr);
	}
	
	//dd if=file ibs=1 skip=? count=?|od -t d1
		
	

} // 	for (i =1; i <= num_samples; i++) {

 write = fwrite(end, 8, 1, ptr);


 printf("Closing file..\n");
 fclose(ptr);

  // cleanup before quitting
 free(filename);

 return 1;
}

int write_array2(int argcc, char **argvv, unsigned int Ntruecpy, unsigned int channelN, unsigned int* Karray, char* head, int binc, char* end){
 filename = (char*) malloc(sizeof(char) * 1024);
 if (filename == NULL) {
   printf("Error in malloc\n");
   exit(1);
 }
 

 // get file path
 char cwd[1024];
 if (getcwd(cwd, sizeof(cwd)) != NULL) {
   
	strcpy(filename, cwd);

	// get filename from command line
	if (argcc < 4) {
	  printf("No output wave file specified\n");
	  return 0;
	}
	
	strcat(filename, "/");
	strcat(filename, argvv[3]);
	printf("%s\n", filename);
 }

 // open file
 printf("Writing file..\n");
 ptr = fopen(filename, "wb");
 if (ptr == NULL) {
	printf("Error opening file\n");
	exit(1);
 }

 int write = 0;

 write = fwrite(head, 46, 1, ptr);
//printf("size of head: %d\n",46);
//printf("size of n: %d %d\n",Ntruecpy, channelN);

for (int i =0; i < Ntruecpy; i++) {

	// dump the data read
	unsigned int  xchannels = 0;
	

	for (xchannels = 0; xchannels < channelN; xchannels ++ ) {
		unsigned char data_buffer[binc];
		if (binc == 4) {
			data_buffer[0] = Karray[i + xchannels*Ntruecpy] & 0xFF;
			data_buffer[1] = (Karray[i + xchannels*Ntruecpy] & 0xFF00) >> 8 ;
			data_buffer[2] = (Karray[i + xchannels*Ntruecpy] & 0xFF0000)  >> 16;
			data_buffer[3] = (Karray[i + xchannels*Ntruecpy] & 0xFF000000) >> 24;
			
		}
		else if (binc == 2) {
			data_buffer[0] = Karray[ i + xchannels*Ntruecpy] & 0xFF;
			data_buffer[1] = (Karray[i + xchannels*Ntruecpy] & 0xFF00) >> 8 ;
		}
		else if (binc == 1) {
			data_buffer[0] = Karray[ i + xchannels*Ntruecpy] ;

		}

		write = fwrite(data_buffer, sizeof(data_buffer), 1, ptr);
	}
	
	//dd if=file ibs=1 skip=? count=?|od -t d1
		
	

} // 	for (i =1; i <= num_samples; i++) {

 write = fwrite(end, 8, 1, ptr);


 printf("Closing file..\n");
 fclose(ptr);

  // cleanup before quitting
 free(filename);

 return 1;
}

int read_head(int argcc, char **argvv, char* mm) {

 filename = (char*) malloc(sizeof(char) * 1024);
 if (filename == NULL) {
   printf("Error in malloc\n");
   exit(1);
 }
 

 // get file path
 char cwd[1024];
 if (getcwd(cwd, sizeof(cwd)) != NULL) {
   
	strcpy(filename, cwd);

	// get filename from command line
	if (argcc < 2) {
	  printf("No wave file specified\n");
	  return 0;
	}
	
	strcat(filename, "/");
	strcat(filename, argvv[1]);
	//printf("%s\n", filename);
 }

 // open file
 //printf("Opening file..\n");
 ptr = fopen(filename, "rb");
 if (ptr == NULL) {
	printf("Error opening file\n");
	exit(1);
 }
 
 int read = 0;
 read = fread(mm, 46, 1, ptr);
 


 return 1;
}

int read_array(int argcc, char **argvv, unsigned int & sampleN, unsigned int & channelN, char * end, unsigned int * Karray) {

 filename = (char*) malloc(sizeof(char) * 1024);
 if (filename == NULL) {
   printf("Error in malloc\n");
   exit(1);
 }
 

 // get file path
 char cwd[1024];
 if (getcwd(cwd, sizeof(cwd)) != NULL) {
   
	strcpy(filename, cwd);

	// get filename from command line
	if (argcc < 2) {
	  printf("No wave file specified\n");
	  return 0;
	}
	
	strcat(filename, "/");
	strcat(filename, argvv[1]);
	printf("%s\n", filename);
 }

 // open file
 printf("Opening file..\n");
 ptr = fopen(filename, "rb");
 if (ptr == NULL) {
	printf("Error opening file\n");
	exit(1);
 }
 
 int read = 0;
 
 // read header parts

 read = fread(header.riff, sizeof(header.riff), 1, ptr);
 //printf("(1-4): %s \n", header.riff); 

 read = fread(buffer4, sizeof(buffer4), 1, ptr);
 //printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);
 
 // convert little endian to big endian 4 byte int
 header.overall_size  = buffer4[0] | 
						(buffer4[1]<<8) | 
						(buffer4[2]<<16) | 
						(buffer4[3]<<24);

 printf("(5-8) Overall size: bytes:%u\n", header.overall_size);

 read = fread(header.wave, sizeof(header.wave), 1, ptr);
 //printf("(9-12) Wave marker: %s\n", header.wave);

 read = fread(header.fmt_chunk_marker, sizeof(header.fmt_chunk_marker), 1, ptr);
 //printf("(13-16) Fmt marker: %s\n", header.fmt_chunk_marker);

 read = fread(buffer4, sizeof(buffer4), 1, ptr);
 //printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);

 // convert little endian to big endian 4 byte integer
 header.length_of_fmt = buffer4[0] |
							(buffer4[1] << 8) |
							(buffer4[2] << 16) |
							(buffer4[3] << 24);
 //printf("(17-20) Length of Fmt header: %u \n", header.length_of_fmt);

 read = fread(buffer2, sizeof(buffer2), 1, ptr); 
 //printf("%u %u \n", buffer2[0], buffer2[1]);
 
 header.format_type = buffer2[0] | (buffer2[1] << 8);
 char format_name[10] = "";
 if (header.format_type == 1)
   strcpy(format_name,"PCM"); 
 else if (header.format_type == 6)
  strcpy(format_name, "A-law");
 else if (header.format_type == 7)
  strcpy(format_name, "Mu-law");

 //printf("(21-22) Format type: %u %s \n", header.format_type, format_name);

 read = fread(buffer2, sizeof(buffer2), 1, ptr);
 //printf("%u %u \n", buffer2[0], buffer2[1]);

 header.channels = buffer2[0] | (buffer2[1] << 8);
 printf("(23-24) Channels: %u \n", header.channels);
 channelN = header.channels;

 read = fread(buffer4, sizeof(buffer4), 1, ptr);
 //printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);

 header.sample_rate = buffer4[0] |
						(buffer4[1] << 8) |
						(buffer4[2] << 16) |
						(buffer4[3] << 24);

 //printf("(25-28) Sample rate: %u\n", header.sample_rate);

 read = fread(buffer4, sizeof(buffer4), 1, ptr);
 //printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);

 header.byterate  = buffer4[0] |
						(buffer4[1] << 8) |
						(buffer4[2] << 16) |
						(buffer4[3] << 24);
 //printf("(29-32) Byte Rate: %u , Bit Rate:%u\n", header.byterate, header.byterate*8);

 read = fread(buffer2, sizeof(buffer2), 1, ptr);
 //printf("%u %u \n", buffer2[0], buffer2[1]);

 header.block_align = buffer2[0] |
					(buffer2[1] << 8);
 //printf("(33-34) Block Alignment: %u \n", header.block_align);

 read = fread(buffer2, sizeof(buffer2), 1, ptr);
 //printf("%u %u \n", buffer2[0], buffer2[1]);

 header.bits_per_sample = buffer2[0] |
					(buffer2[1] << 8);
 printf("(35-36) Bits per sample: %u \n", header.bits_per_sample);

 read = fread(header.data_chunk_header, sizeof(header.data_chunk_header), 1, ptr);
 //printf("(37-40) Data Marker: %s \n", header.data_chunk_header);

 read = fread(buffer4, sizeof(buffer4), 1, ptr);
 //printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);

 header.data_size = buffer4[0] |
				(buffer4[1] << 8) |
				(buffer4[2] << 16) | 
				(buffer4[3] << 24 );
 //printf("(41-44) Size of data chunk: %u \n", header.data_size);
 read = fread(buffer2, sizeof(buffer2), 1, ptr);

 // calculate no.of samples
 long num_samples = 8 * (header.overall_size - 46) / (header.channels * header.bits_per_sample);
 //printf("Number of samples:%lu \n", num_samples);
 sampleN = 8 * (header.overall_size - 46) / (header.channels * header.bits_per_sample);

 long size_of_each_sample = (header.bits_per_sample) / 8;
 //printf("Size of each sample:%ld bytes\n", size_of_each_sample);

 // calculate duration of file
 float duration_in_seconds = (float) header.overall_size / header.byterate;
 //printf("Approx.Duration in seconds=%f\n", duration_in_seconds);
 //printf("Approx.Duration in h:m:s=%s\n", seconds_to_time(duration_in_seconds));

 

 //int* array_of_samples;
 //array_of_samples = (int *)malloc ( num_samples * header.channels * sizeof(int) );

 // read each sample from data chunk if PCM
 if (header.format_type == 1) { // PCM
    /**printf("Dump sample data? Y/N?");
	char c = 'n';
	scanf("%c", &c);**/
	char c = 'y';
	if (c == 'Y' || c == 'y') { 
		long i =0;
		char data_buffer[size_of_each_sample];
		int  size_is_correct = TRUE;

		// make sure that the bytes-per-sample is completely divisible by num.of channels
		long bytes_in_each_channel = (size_of_each_sample / header.channels);
		if ((bytes_in_each_channel  * header.channels) != size_of_each_sample) {
			printf("Error: %ld x %ud <> %ld\n", bytes_in_each_channel, header.channels, size_of_each_sample);
			size_is_correct = FALSE;
		}
 
		if (size_is_correct) { 
					// the valid amplitude range for values based on the bits per sample
			long low_limit = 0l;
			long high_limit = 0l;

			switch (header.bits_per_sample) {
				case 8:
					low_limit = -128;
					high_limit = 127;
					break;
				case 16:
					low_limit = -32768;
					high_limit = 32767;
					break;
				case 32:
					low_limit = -2147483648;
					high_limit = 2147483647;
					break;
			}					

			//printf("\n\n.Valid range for data values : %ld to %ld \n", low_limit, high_limit);			
			
			
			for (i =1; i <= num_samples; i++) {
				//printf("==========Sample %ld / %ld=============\n", i, num_samples);
				
				if (read == 1) {
				
					// dump the data read
					unsigned int  xchannels = 0;
					int data_in_channel = 0;

					for (xchannels = 0; xchannels < header.channels; xchannels ++ ) {
						read = fread(data_buffer, sizeof(data_buffer), 1, ptr);
						// convert data from little endian to big endian based on bytes in each channel sample
						if (bytes_in_each_channel == 4) {
							data_in_channel =	data_buffer[0] | 
												(data_buffer[1]<<8) | 
												(data_buffer[2]<<16) | 
												(data_buffer[3]<<24);

						}
						else if (bytes_in_each_channel == 2) {
							data_in_channel = data_buffer[0] | (data_buffer[1] << 8) ;

						}
						else if (bytes_in_each_channel == 1) {
							data_in_channel = data_buffer[0];

						}

						

						//put data into array_of_samples
						Karray [ (i - 1 + num_samples * xchannels) ] = data_in_channel;
						

						// check if value was in range
						if (data_in_channel < low_limit || data_in_channel > high_limit)
							printf("**value out of range\n");

						
					}

					
				}
				else {
					printf("Error reading file. %d bytes\n", read);
					break;
				}

			} // 	for (i =1; i <= num_samples; i++) {

		} // 	if (size_is_correct) { 

	 } // if (c == 'Y' || c == 'y') { 
 } //  if (header.format_type == 1) { 

read = fread(end, 8, 1, ptr); 

 printf("Closing file..\n");
 fclose(ptr);

  // cleanup before quitting
 free(filename);
 return 1;

}

int read_length(int argcc, char **argvv)
{
	filename = (char*) malloc(sizeof(char) * 1024);
 	if (filename == NULL) {
	   printf("Error in malloc\n");
	   exit(1);
	}
	 

 // get file path
 char cwd[1024];
 if (getcwd(cwd, sizeof(cwd)) != NULL) {
   
	strcpy(filename, cwd);

	// get filename from command line
	if (argcc < 2) {
	  printf("No wave file specified\n");
	  return 0;
	}
	
	strcat(filename, "/");
	strcat(filename, argvv[1]);
	printf("%s\n", filename);
 }

 // open file
 printf("Opening file..\n");
 ptr = fopen(filename, "rb");
 if (ptr == NULL) {
	printf("Error opening file\n");
	exit(1);
 }
 
 int read = 0;
 
 // read header parts

 read = fread(header.riff, sizeof(header.riff), 1, ptr);
 //printf("(1-4): %s \n", header.riff); 

 read = fread(buffer4, sizeof(buffer4), 1, ptr);
 //printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);
 
 // convert little endian to big endian 4 byte int
 header.overall_size  = buffer4[0] | 
						(buffer4[1]<<8) | 
						(buffer4[2]<<16) | 
						(buffer4[3]<<24);

 printf("(5-8) Overall size: bytes:%u\n", header.overall_size);

 read = fread(header.wave, sizeof(header.wave), 1, ptr);
 //printf("(9-12) Wave marker: %s\n", header.wave);

 read = fread(header.fmt_chunk_marker, sizeof(header.fmt_chunk_marker), 1, ptr);
 //printf("(13-16) Fmt marker: %s\n", header.fmt_chunk_marker);

 read = fread(buffer4, sizeof(buffer4), 1, ptr);
 //printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);

 // convert little endian to big endian 4 byte integer
 header.length_of_fmt = buffer4[0] |
							(buffer4[1] << 8) |
							(buffer4[2] << 16) |
							(buffer4[3] << 24);
 //printf("(17-20) Length of Fmt header: %u \n", header.length_of_fmt);

 read = fread(buffer2, sizeof(buffer2), 1, ptr); 
 //printf("%u %u \n", buffer2[0], buffer2[1]);
 
 header.format_type = buffer2[0] | (buffer2[1] << 8);
 char format_name[10] = "";
 if (header.format_type == 1)
   strcpy(format_name,"PCM"); 
 else if (header.format_type == 6)
  strcpy(format_name, "A-law");
 else if (header.format_type == 7)
  strcpy(format_name, "Mu-law");

 //printf("(21-22) Format type: %u %s \n", header.format_type, format_name);

 read = fread(buffer2, sizeof(buffer2), 1, ptr);
 //printf("%u %u \n", buffer2[0], buffer2[1]);

 header.channels = buffer2[0] | (buffer2[1] << 8);
 printf("(23-24) Channels: %u \n", header.channels);


 read = fread(buffer4, sizeof(buffer4), 1, ptr);
 //printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);

 header.sample_rate = buffer4[0] |
						(buffer4[1] << 8) |
						(buffer4[2] << 16) |
						(buffer4[3] << 24);

 //printf("(25-28) Sample rate: %u\n", header.sample_rate);

 read = fread(buffer4, sizeof(buffer4), 1, ptr);
 //printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);

 header.byterate  = buffer4[0] |
						(buffer4[1] << 8) |
						(buffer4[2] << 16) |
						(buffer4[3] << 24);
 //printf("(29-32) Byte Rate: %u , Bit Rate:%u\n", header.byterate, header.byterate*8);

 read = fread(buffer2, sizeof(buffer2), 1, ptr);
 //printf("%u %u \n", buffer2[0], buffer2[1]);

 header.block_align = buffer2[0] |
					(buffer2[1] << 8);
 //printf("(33-34) Block Alignment: %u \n", header.block_align);

 read = fread(buffer2, sizeof(buffer2), 1, ptr);
 //printf("%u %u \n", buffer2[0], buffer2[1]);

 header.bits_per_sample = buffer2[0] |
					(buffer2[1] << 8);
 printf("(35-36) Bits per sample: %u \n", header.bits_per_sample);

 read = fread(header.data_chunk_header, sizeof(header.data_chunk_header), 1, ptr);
 //printf("(37-40) Data Marker: %s \n", header.data_chunk_header);

 read = fread(buffer4, sizeof(buffer4), 1, ptr);
 //printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);

 header.data_size = buffer4[0] |
				(buffer4[1] << 8) |
				(buffer4[2] << 16) | 
				(buffer4[3] << 24 );
 //printf("(41-44) Size of data chunk: %u \n", header.data_size);
 read = fread(buffer2, sizeof(buffer2), 1, ptr);
 
 fclose(ptr);

  // cleanup before quitting
 free(filename);

 // calculate no.of samples
 return 8 * (header.overall_size - 46) / (header.channels * header.bits_per_sample);

}



int log2(int N)    /*function to calculate the log2(.) of int numbers*/
{
  int k = N, i = 0;
  while(k) {
    k >>= 1;
    i++;
  }
  return i - 1;
}

int check(int n)    //checking if the number of element is a power of 2
{
  return n > 0 && (n & (n - 1)) == 0;
}

int reverse(int N, int n)    //calculating revers number
{
  int j, p = 0;
  for(j = 1; j <= log2(N); j++) {
    if(n & (1 << (log2(N) - j)))
      p |= 1 << (j - 1);
  }
  return p;
}

void ordina(complex<double>* f1, int N) //using the reverse order in the array
{
  complex<double> f2[MAX];
  for(int i = 0; i < N; i++)
    f2[i] = f1[reverse(N, i)];
  for(int j = 0; j < N; j++)
    f1[j] = f2[j];
}

void high_low_filter(complex<double>* f, complex<double>* ori, int N, int fre, int lowfre){
	for (int i = 0; i< N; i++){
		if ( (real(f[i])>fre) || (real(f[i])<lowfre) ){
				ori[i] = complex<double> (0, 0);
		}
	}
}

void transform(complex<double>* f, int N) //
{
  ordina(f, N);    //first: reverse order
  complex<double> *W;
  W = (complex<double> *)malloc(N / 2 * sizeof(complex<double>));
  W[1] = polar(1., -2. * M_PI / N);
  W[0] = 1;
  for(int i = 2; i < N / 2; i++)
    W[i] = pow(W[1], i);
  int n = 1;
  int a = N / 2;
  for(int j = 0; j < log2(N); j++) {
    for(int i = 0; i < N; i++) {
      if(!(i & n)) {
        complex<double> temp = f[i];
        complex<double> Temp = W[(i * a) % (n * a)] * f[i + n];
        f[i] = temp + Temp;
        f[i + n] = temp - Temp;
      }
    }
    n *= 2;
    a = a / 2;
  }
}

void FFT(complex<double>* f, int N, double d)
{
  transform(f, N);
  for(int i = 0; i < N; i++)
    f[i] *= d; //multiplying by step
}

