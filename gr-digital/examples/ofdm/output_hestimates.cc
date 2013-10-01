//just output the file that was received on input
//#include <cstdio>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdexcept>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <complex>
#include <iostream>
using namespace std; 

typedef std::complex<float>             gr_complex;

void set_data(const char*);

#define DATA_TYPE gr_complex

#define OCCUPIED_CARRIERS 80
#define DC_CARRIERS 8
#define DATA_CARRIERS (OCCUPIED_CARRIERS-DC_CARRIERS)


DATA_TYPE *d_data; //to store the data stream
long data_length;


int main (int argc, const char* argv[]) {
  if (argc != 2) {
    fprintf (stderr, " usage: %s <data input file> \n",argv[0]);
    exit(1);
  }
  set_data(argv[1]);
  

  return 0; 
}

void set_data(const char* filename) {
  FILE *d_fp;
  
  if((d_fp = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "data file cannot be opened\n");
    assert(false);
  }

  //get file size
  fseek( d_fp, 0L, SEEK_END );
  long endPos = ftell( d_fp );
  fclose(d_fp); 

  d_data = (DATA_TYPE*) malloc(endPos); 

  //re-open file
  if((d_fp = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "data file cannot be opened\n");
    assert(false);
  }
  
  long count = fread_unlocked(d_data, sizeof(DATA_TYPE), endPos/sizeof(DATA_TYPE), d_fp);
  int num_pkts = count/OCCUPIED_CARRIERS;

  int start_dc = DATA_CARRIERS/2;
  int end_dc = start_dc + DC_CARRIERS;

  //printf("start_dc: %d, end_dc: %d, num_pkts: %d\n", start_dc, end_dc, num_pkts); fflush(stdout);

#if 0 
  count = 0;
  for(int i = 0; i < num_pkts; i++) {
	 for(int j = 0; j < OCCUPIED_CARRIERS; j++) {
		if(j >= start_dc && j < end_dc) {
		   j++;
		   continue;
		}
	
		long index = (i*OCCUPIED_CARRIERS) + j;
		float rad = arg(gr_complex(1.0, 0.0)/d_data[index]);		// plot the channel = 1/equalizer
		//float rad = arg(d_data[index]);
		float deg = rad * 180/M_PI;
		if(deg < 0) 
		   deg += 360;

		//assert(deg <= 360);
		cout << count++ << " " << abs(gr_complex(1.0, 0.0)/d_data[index]) << " " << rad << " " << deg << " " << d_data[index] << " " << abs(d_data[index]) << endl;
	 }
	 //count+= 50;
  }

#else
  count = 0;
  const int NUM_SENDERS = 2;
  for(int i = 0; i < num_pkts; i+=NUM_SENDERS) {

     float start_rad[NUM_SENDERS], end_rad[NUM_SENDERS], slope[NUM_SENDERS]; 			// for slope calc within a packet //
     float prev_rad[NUM_SENDERS], diff[NUM_SENDERS];

     for(int j = 0; j < OCCUPIED_CARRIERS; j++) {
	 if(j >= start_dc && j < end_dc) {
            j++;
            continue;
         }
	
	 long index = (i*OCCUPIED_CARRIERS) + j;
	 float rad[NUM_SENDERS], deg[NUM_SENDERS];

	 for(int k = 0; k < NUM_SENDERS; k++) {
	    index += (k * OCCUPIED_CARRIERS); 
	    rad[k] = arg(gr_complex(1.0, 0.0)/d_data[index]); 
	    deg[k] = rad[k] * 180/M_PI;
	    if(deg[k] < 0) 
		deg[k] += 360;

	    /* to calculate the slope within the packet (across subcarriers) */
	    if(j == 0) 	start_rad[k] = rad[k];
	    if(j == OCCUPIED_CARRIERS-1) {
		end_rad[k] = rad[k];
		if(end_rad[k] > start_rad[k])
		   end_rad[k] = -M_PI-(M_PI - end_rad[k]);
		slope[k] = (end_rad[k] - start_rad[k])/((float) DATA_CARRIERS);
	    }

	    if(j  == 0) {
		diff[k] = start_rad[k] - prev_rad[k];
		if(diff[k] < 0) diff[k] += (2*M_PI);
		prev_rad[k] = start_rad[k]; 
	    }
	 }

	cout << count++ << " " << rad[0] << " " << deg[0]  << " " << rad[1] << " " << deg[1] << endl;
     }
     //cout << count++ << " " << slope[0] << " " << slope[1] << endl;

     //cout << count++ << " " << start_rad[0] << " " << start_rad[1] << " " << diff[0] << " " << diff[1] << endl;
     //count += 50;
  }
#endif

  fclose(d_fp);
  d_fp = 0;
}
