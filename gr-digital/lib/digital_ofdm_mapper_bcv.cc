/* -*- c++ -*- */
/*
 * Copyright 2006-2008,2010,2011 Free Software Foundation, Inc.
 * 
 * This file is part of GNU Radio
 * 
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#define REFSNRTX 1
#define LTFPREAMBLE 0


#include <digital_ofdm_mapper_bcv.h>
#include <gr_io_signature.h>
#include <stdexcept>
#include <string.h>
#include <cstdio>

#if REFSNRTX == 1
#include <fstream>
std::ofstream fptr_refsym;
#endif

#if LTFPREAMBLE == 1
int* ltf_preamble;
int d_ltfpreamble_size;
int ltf_counter,ltf_preamble_done;
std::vector<int> d_ltfpreamble_carriers;
#endif


digital_ofdm_mapper_bcv_sptr
digital_make_ofdm_mapper_bcv (const std::vector<gr_complex> &constellation, unsigned int msgq_limit, 
			      unsigned int occupied_carriers, unsigned int fft_length)
{
  return gnuradio::get_initial_sptr(new digital_ofdm_mapper_bcv (constellation, msgq_limit, 
								 occupied_carriers, fft_length));
}

// Consumes 1 packet and produces as many OFDM symbols of fft_length to hold the full packet
digital_ofdm_mapper_bcv::digital_ofdm_mapper_bcv (const std::vector<gr_complex> &constellation, unsigned int msgq_limit, 
						  unsigned int occupied_carriers, unsigned int fft_length)
  : gr_sync_block ("ofdm_mapper_bcv",
		   gr_make_io_signature (0, 0, 0),
		   gr_make_io_signature2 (1, 2, sizeof(gr_complex)*fft_length, sizeof(char))),
    d_constellation(constellation),
    d_msgq(gr_make_msg_queue(msgq_limit)), d_msg_offset(0), d_eof(false),
    d_occupied_carriers(occupied_carriers),
    d_fft_length(fft_length),
    d_bit_offset(0),
    d_pending_flag(0),
    d_resid(0),
    d_nresid(0)
{
  if (!(d_occupied_carriers <= d_fft_length))
    throw std::invalid_argument("digital_ofdm_mapper_bcv: occupied carriers must be <= fft_length");

  assign_subcarriers();


  #if LTFPREAMBLE == 1
    int ltf_preamble_tmp[]= {1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,
                1, -1, -1,1,1, -1,1, -1,1, -1, -1, -1, -1, -1,1,1, -1, -1,1, -1,1,-1,
                1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,
                1, -1, -1,1,1, -1,1, -1,1, -1, -1, -1, -1, -1,1,1, -1, -1,1, -1,1,-1};

    d_ltfpreamble_size = d_data_carriers.size()+d_pilot_carriers.size();
    ltf_preamble = new int[d_ltfpreamble_size];
    memcpy(ltf_preamble,&ltf_preamble_tmp[0],sizeof(int)*d_ltfpreamble_size);

    printf("d_ltfpreamble size:%d\n",d_ltfpreamble_size);
    for(int i=0; i<d_ltfpreamble_size; i++){
        printf("%d,",ltf_preamble[i]);
    }
    printf("\n");
  #endif

  // make sure we stay in the limit currently imposed by the occupied_carriers
  if(d_data_carriers.size() > d_occupied_carriers) {
    throw std::invalid_argument("digital_ofdm_mapper_bcv: subcarriers allocated exceeds size of occupied carriers");
  }
  
  d_nbits = (unsigned long)ceil(log10(float(d_constellation.size())) / log10(2.0));

  #if REFSNRTX == 1
    fptr_refsym.open("ref_symbols_qpsk.txt");
  #endif


}

digital_ofdm_mapper_bcv::~digital_ofdm_mapper_bcv(void)
{
  #if REFSNRTX == 1
    fptr_refsym.close();
  #endif
  #if LTFPREAMBLE == 1
    free(ltf_preamble);
  #endif
}

int digital_ofdm_mapper_bcv::randsym()
{
  return (rand() % d_constellation.size());
}

void
digital_ofdm_mapper_bcv::assign_subcarriers() {
  int dc_tones = 8;
  int num_pilots = 8;
  int pilot_gap = 11;

  int half1_end = (d_occupied_carriers-dc_tones)/2;     //40
  int half2_start = half1_end+dc_tones;                 //48

  int off = (d_fft_length-d_occupied_carriers)/2;	//4
  
  // first half
  for(int i = 0; i < half1_end; i++) {
     if(i%pilot_gap == 0)
        d_pilot_carriers.push_back(i+off);
     else
        d_data_carriers.push_back(i+off);

     #if LTFPREAMBLE == 1
     d_ltfpreamble_carriers.push_back(i+off);
     #endif
  }

  // second half
  for(int i = half2_start, j = 0; i < d_occupied_carriers; i++, j++) {
     if(j%pilot_gap == 0)
        d_pilot_carriers.push_back(i+off);
     else
        d_data_carriers.push_back(i+off);

     #if LTFPREAMBLE == 1
     d_ltfpreamble_carriers.push_back(i+off);
     #endif
  }

  /* debug carriers */
  printf("pilot carriers: \n");
  for(int i = 0; i < d_pilot_carriers.size(); i++) {
     printf("%d ", d_pilot_carriers[i]); fflush(stdout);
  }
  printf("\n");

  printf("data carriers: \n");
  for(int i = 0; i < d_data_carriers.size(); i++) {
     printf("%d ", d_data_carriers[i]); fflush(stdout);
  }
  printf("\n");
}



int
digital_ofdm_mapper_bcv::work(int noutput_items,
			      gr_vector_const_void_star &input_items,
			      gr_vector_void_star &output_items)
{
  gr_complex *out = (gr_complex *)output_items[0];
  
  unsigned int i=0;

  //printf("OFDM BPSK Mapper:  ninput_items: %d   noutput_items: %d\n", ninput_items[0], noutput_items);

  if(d_eof) {
    return -1;
  }
  
  if(!d_msg) {
    d_msg = d_msgq->delete_head();	   // block, waiting for a message
    d_msg_offset = 0;
    d_bit_offset = 0;
    d_pending_flag = 1;			   // new packet, write start of packet flag

    #if LTFPREAMBLE == 1
    ltf_preamble_done = 0;
    ltf_counter = 0;
    #endif
    
    if((d_msg->length() == 0) && (d_msg->type() == 1)) {
      d_msg.reset();
      return -1;		// We're done; no more messages coming.
    }
  }

  char *out_flag = 0;
  if(output_items.size() == 2)
    out_flag = (char *) output_items[1];
  

  // Build a single symbol:
  // Initialize all bins to 0 to set unused carriers
  memset(out, 0, d_fft_length*sizeof(gr_complex));
  
 #if LTFPREAMBLE == 1
 if(ltf_preamble_done == 0){
     generateLTFPreamble(out);

     if (out_flag)
       out_flag[0] = d_pending_flag;
     d_pending_flag = 0;

     ltf_counter += 1;

     if(ltf_counter == 2){
       ltf_preamble_done = 1;
     }
     /*for (int ll=0; ll<noutput_items; ll++){
       printf(" %f + %f j\n",out[ll].real(),out[ll].imag());
     }*/
     return 1;
 }
#endif


  i = 0;
  unsigned char bits = 0;
  while((d_msg_offset < d_msg->length()) && (i < d_data_carriers.size())) {

    // need new data to process
    if(d_bit_offset == 0) {
      d_msgbytes = d_msg->msg()[d_msg_offset];
      //printf("mod message byte: %x\n", d_msgbytes);
    }

    if(d_nresid > 0) {
      // take the residual bits, fill out nbits with info from the new byte, and put them in the symbol
      d_resid |= (((1 << d_nresid)-1) & d_msgbytes) << (d_nbits - d_nresid);
      bits = d_resid;

      out[d_data_carriers[i]] = d_constellation[bits];
      i++;

      #if REFSNRTX == 1
        fptr_refsym<<(int) bits<<",";
      #endif

      d_bit_offset += d_nresid;
      d_nresid = 0;
      d_resid = 0;
      //printf("mod bit(r): %x   resid: %x   nresid: %d    bit_offset: %d\n", 
      //     bits, d_resid, d_nresid, d_bit_offset);
    }
    else {
      if((8 - d_bit_offset) >= d_nbits) {  // test to make sure we can fit nbits
	// take the nbits number of bits at a time from the byte to add to the symbol
	bits = ((1 << d_nbits)-1) & (d_msgbytes >> d_bit_offset);
	d_bit_offset += d_nbits;
	
	out[d_data_carriers[i]] = d_constellation[bits];
	i++;
      #if REFSNRTX == 1
        fptr_refsym<<(int) bits<<",";
      #endif

      }
      else {  // if we can't fit nbits, store them for the next 
	// saves d_nresid bits of this message where d_nresid < d_nbits
	unsigned int extra = 8-d_bit_offset;
	d_resid = ((1 << extra)-1) & (d_msgbytes >> d_bit_offset);
	d_bit_offset += extra;
	d_nresid = d_nbits - extra;
      }
      
    }
            
    if(d_bit_offset == 8) {
      d_bit_offset = 0;
      d_msg_offset++;
    }
  }

  // Ran out of data to put in symbol
  if (d_msg_offset == d_msg->length()) {
    if(d_nresid > 0) {
      d_resid |= 0x00;
      bits = d_resid;
      d_nresid = 0;
      d_resid = 0;
    }

    #if REFSNRTX == 1
     int ref_cntr=0;
     // bpsk
     //int trace_vals[]={1,0,1,0,0,1,1,1,1,0,0,1,1,0,1,1,1,0,0,1,0,0,1,1,0,0,0,1,1,0,1,1,0,0,1,0,1,0,0,0,0,0,0,1,1,0,0,0,1,0,1,1,0,0,1,0,0,0,0,1,0,0,0,1,1,1,1,0,1,1,1,0,0,0,1,0,0,1,1,0,0};
     //qpsk
     int trace_vals[]={3,3,1,0,0,1,0,0,2,2,0,1,1,0,0,3,1,1,3,2,2,0,1,3,1,1,3,2,3,3,1,3,1,1,3,1,2,2,1,3,1,1,2,3,1,2,1,3,1,2,0,3,0,2,1,3,0,1,1,0,0,2,1,0,1,3,0,1,0,2,0,2,3,1,3,1,0,1,2,1,1};
    #endif
    //while(i < d_occupied_carriers) {   // finish filling out the symbol
    while(i < d_data_carriers.size()) {   // finish filling out the symbol

      #if REFSNRTX == 1
        int padbits = trace_vals[ref_cntr];
        ref_cntr++;
        out[d_data_carriers[i]] = d_constellation[padbits];
        fptr_refsym<<(int) padbits<<",";
      #else
        out[d_data_carriers[i]] = d_constellation[randsym()];
      #endif
      i++;
    }

    if (d_msg->type() == 1)	        // type == 1 sets EOF
      d_eof = true;
    d_msg.reset();   			// finished packet, free message
    assert(d_bit_offset == 0);
  }


  #if REFSNRTX == 1
    fptr_refsym<<(int) bits<<"\n";
  #endif


   // now add the pilot symbols //
   double cur_pilot = 1.0;
   for(int i = 0; i < d_pilot_carriers.size(); i++) {
       out[d_pilot_carriers[i]] = gr_complex(cur_pilot, 0.0);
       cur_pilot = -cur_pilot;
   }

  if (out_flag)
    out_flag[0] = d_pending_flag;
  d_pending_flag = 0;

  return 1;  // produced symbol
}

#if LTFPREAMBLE == 1
void
digital_ofdm_mapper_bcv::generateLTFPreamble(gr_complex* out){

    //unsigned int i = 0;

    for(int i=0; i < d_ltfpreamble_carriers.size(); i++){
        out[d_ltfpreamble_carriers[i]] = ltf_preamble[i];
    }

    /*while( i < d_ltf_preamble_carriers.size()){
        out[d_ltfpreamble_carriers[i]] = ltf_preamble[i];
        i++;
    }*/
    //printf("generating ltf preamble.\n");

}
#endif
