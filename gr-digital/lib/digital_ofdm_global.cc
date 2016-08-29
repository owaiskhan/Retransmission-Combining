
#include "digital_ofdm_global.h"

double current_time(){
      struct timeval tv;
      if ( gettimeofday(&tv, 0) < 0 ){
        std::cout<<"gettimeofday failed in current_time()"<<std::endl;
        exit(0);
      }

      return double(tv.tv_sec) + (double(tv.tv_usec) / 1e6);
}

