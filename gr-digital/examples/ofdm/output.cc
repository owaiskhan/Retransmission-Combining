//just output the file that was received on input
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

DATA_TYPE *d_data; //to store the data stream
int data_length;


int main (int argc, const char* argv[]) {
  if (argc != 2) {
    fprintf (stderr, " usage: %s <data input file>\n",argv[0]);
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
  
  int count = fread_unlocked(d_data, sizeof(DATA_TYPE), endPos/sizeof(DATA_TYPE), d_fp);
  for(int index = 0; index < count; index++) {
    cout << index << " " << d_data[index].real() << " " << d_data[index].imag() << " " << abs(d_data[index]) << '\n';
    //cout << index << " " << arg(d_data[index]) << " " << abs(d_data[index]) << " " << abs(d_data[index]) << '\n';
  }
  data_length = count;
  fclose(d_fp);
  d_fp = 0;
}
