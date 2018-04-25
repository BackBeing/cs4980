#include <complex>
#include <unistd.h>
#include <stdio.h>
#include <string>
#include <cstring>
#include <stdlib.h>
#include "wave.h"
#define TRUE 1 
#define FALSE 0

using namespace std;


// WAVE header structure

unsigned char buffer4[4];
unsigned char buffer2[2];

char* seconds_to_time(float seconds);
int* read_array(int argcc, char **argvv, unsigned int & sampleN, unsigned int & channelN);


 FILE *ptr;
 char *filename;
 struct HEADER header;

//////////////////////////////////////////////////////////
// Sampling frequency = 44100 Hz
// Sampling Resolution = 16 bits   65 536 = 2^16
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
	int * Karray;
	
	unsigned int Ntrue = 0;
	unsigned int channelN = 0;
	
	//The wave file is now in an array: Karray!!!!
	//Ntrue is the sample numbers
	Karray = read_array( argc, argv, Ntrue, channelN );

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
	
/** check point
	for (int i = 0; i< sampleN; i++){
		printf("expected sample number real: %d image: 0 \n", Karray[ (i * N) ]);
		printf("sample number real: %1f  image: %1f \n", real(Kgrps[i][0]), imag(Kgrps[i][0]))	;	
	}
**/	
	
	//set up W, which can be used for all sample groups
	//TODO: build fixed W on board:// low
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
	printf("Data transformed..\n");
	
	//TODO: pass tKgrps real, tKgrps imag, W real, W imag, and N(maybe), sampleN to board for FFT convert:// high
	
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
	
	
	//TODO: do high/low filter and sound effect things on tKgrpsNew:// high 
	
	
	
	
	
	//TODO: pass tKgrpsNew real, tKgrpsNew imag, W real, W imag, and N(maybe), sampleN to board for IFFT convert:// high
	
	//TODO: get tKgrpsNewNew back from board:// high
	
	//here just pass Kgrps to tKgrpsNewNew
	complex<double> ** tKgrpsNewNew;
	tKgrpsNewNew= (complex<double> **)malloc( sampleN * sizeof(complex<double> *));
	for (int i = 0; i< sampleN; i++){
			tKgrpsNewNew[i] = (complex<double> *)malloc(N * sizeof(complex<double>));
			for (int j = 0; j< N; j++){
				tKgrpsNewNew[i][j] = Kgrps[i][j];
			}
	}
	
	
	//rebuilding the sound array Karray
	for (int i = 0; i< sampleN; i++){
		for (int j = 0; j< N; j++){
			if ( (i * N + j) < Ntruecpy ){
				Karray[ (i * N + j)] = (int) real(Kgrps[i][j]);	
			}else{
				//do nothing
			}	
		}
	}
	//TODO: write to a wave file:// medium




	// free
	for (int i = 0; i< sampleN; i++){
			free(Kgrps[i]);
	}
	for (int i = 0; i< sampleN; i++){
			free(tKgrps[i]);
	}
	
	
}

int* read_array(int argcc, char **argvv, unsigned int & sampleN, unsigned int & channelN) {

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

 //printf("(5-8) Overall size: bytes:%u, Kb:%u \n", header.overall_size, header.overall_size/1024);

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
 //printf("(23-24) Channels: %u \n", header.channels);
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
 //printf("(35-36) Bits per sample: %u \n", header.bits_per_sample);

 read = fread(header.data_chunk_header, sizeof(header.data_chunk_header), 1, ptr);
 //printf("(37-40) Data Marker: %s \n", header.data_chunk_header);

 read = fread(buffer4, sizeof(buffer4), 1, ptr);
 //printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);

 header.data_size = buffer4[0] |
				(buffer4[1] << 8) |
				(buffer4[2] << 16) | 
				(buffer4[3] << 24 );
 //printf("(41-44) Size of data chunk: %u \n", header.data_size);


 // calculate no.of samples
 long num_samples = 8 * (header.overall_size - 38) / (header.channels * header.bits_per_sample);
 //printf("Number of samples:%lu \n", num_samples);
 sampleN = 8 * (header.overall_size - 38) / (header.channels * header.bits_per_sample);

 long size_of_each_sample = (header.channels * header.bits_per_sample) / 8;
 //printf("Size of each sample:%ld bytes\n", size_of_each_sample);

 // calculate duration of file
 float duration_in_seconds = (float) header.overall_size / header.byterate;
 //printf("Approx.Duration in seconds=%f\n", duration_in_seconds);
 //printf("Approx.Duration in h:m:s=%s\n", seconds_to_time(duration_in_seconds));

 

 int* array_of_samples;
 array_of_samples = (int *)malloc ( num_samples * header.channels * sizeof(int) );

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
				read = fread(data_buffer, sizeof(data_buffer), 1, ptr);
				if (read == 1) {
				
					// dump the data read
					unsigned int  xchannels = 0;
					int data_in_channel = 0;

					for (xchannels = 0; xchannels < header.channels; xchannels ++ ) {
						
						// convert data from little endian to big endian based on bytes in each channel sample
						if (bytes_in_each_channel == 4) {
							data_in_channel =	data_buffer[0] | 
												(data_buffer[1]<<8) | 
												(data_buffer[2]<<16) | 
												(data_buffer[3]<<24);

						}
						else if (bytes_in_each_channel == 2) {
							data_in_channel = data_buffer[0] |
												(data_buffer[1] << 8);

						}
						else if (bytes_in_each_channel == 1) {
							data_in_channel = data_buffer[0];

						}

						

						//put data into array_of_samples
						array_of_samples [ (i - 1 + num_samples * xchannels) ] = data_in_channel;
						

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

 printf("Closing file..\n");
 fclose(ptr);

  // cleanup before quitting
 free(filename);
 return array_of_samples;

}

