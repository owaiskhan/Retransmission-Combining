/* -*- c++ -*- */
/*
 * Copyright 2007,2008,2010,2011 Free Software Foundation, Inc.
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

#include <digital_ofdm_frame_sink.h>
#include <gr_io_signature.h>
#include <gr_expj.h>
#include <gr_math.h>
#include <math.h>
#include <cstdio>
#include <stdexcept>
#include <iostream>
#include <string.h>

#define VERBOSE 0
#define REFSNR 1

#if REFSNR == 1
#include <fstream>
unsigned int ofdm_counter;
std::vector<std::vector<int> > refsym_indices;
std::vector<double> ref_subcarrier_error;
std::ofstream outfile_snr;
#endif

inline void
digital_ofdm_frame_sink::enter_search()
{
  if (VERBOSE)
    fprintf(stderr, "@ enter_search\n");

  d_state = STATE_SYNC_SEARCH;

}
    
inline void
digital_ofdm_frame_sink::enter_have_sync()
{
  if (VERBOSE)
    fprintf(stderr, "@ enter_have_sync\n");

  d_state = STATE_HAVE_SYNC;

  // clear state of demapper
  d_byte_offset = 0;
  d_partial_byte = 0;

  d_header = 0;
  d_headerbytelen_cnt = 0;

  // Resetting PLL
  d_freq = 0.0;
  d_phase = 0.0;
  fill(d_dfe.begin(), d_dfe.end(), gr_complex(1.0,0.0));
}

inline void
digital_ofdm_frame_sink::enter_have_header()
{
  d_state = STATE_HAVE_HEADER;

  // header consists of two 16-bit shorts in network byte order
  // payload length is lower 12 bits
  // whitener offset is upper 4 bits
  d_packetlen = (d_header >> 16) & 0x0fff;
  d_packet_whitener_offset = (d_header >> 28) & 0x000f;
  d_packetlen_cnt = 0;

  if (VERBOSE)
    fprintf(stderr, "@ enter_have_header (payload_len = %d) (offset = %d)\n", 
	    d_packetlen, d_packet_whitener_offset);
}


unsigned char digital_ofdm_frame_sink::slicer(const gr_complex x)
{
  unsigned int table_size = d_sym_value_out.size();
  unsigned int min_index = 0;
  float min_euclid_dist = norm(x - d_sym_position[0]);
  float euclid_dist = 0;
  
  for (unsigned int j = 1; j < table_size; j++){
    euclid_dist = norm(x - d_sym_position[j]);
    if (euclid_dist < min_euclid_dist){
      min_euclid_dist = euclid_dist;
      min_index = j;
    }
  }
  return d_sym_value_out[min_index];
}

/* uses pilots - from rawofdm -- optionally returns if needed */
inline void
digital_ofdm_frame_sink::equalize_interpolate_dfe(gr_complex *in, gr_complex *factor) 
{
  gr_complex phase_error = 0.0;

  float cur_pilot = 1.0;
  for (unsigned int i = 0; i < d_pilot_carriers.size(); i++) {
    gr_complex pilot_sym(cur_pilot, 0.0);
    cur_pilot = -cur_pilot;
    int di = d_pilot_carriers[i];
    phase_error += (in[di] * conj(pilot_sym));
    /*if(d_packetid==6){
        printf("in: (%f, %f), pilot: (%f, %f) in(rad): (%f, %f)\n", in[di].real(), in[di].imag(), pilot_sym.real(), pilot_sym.imag(),abs(in[di]),arg(in[di]));
        fflush(stdout);
    }*/
  } 

  // update phase equalizer
  float angle = arg(phase_error);

  /*if(d_packetid==6){
     std::cout<<"phase error angle: "<<angle<<std::endl;
  }*/
  d_freq = d_freq - d_freq_gain*angle;
  d_phase = d_phase + d_freq - d_phase_gain*angle;
  if (d_phase >= 2*M_PI) d_phase -= 2*M_PI;
  else if (d_phase <0) d_phase += 2*M_PI;

  gr_complex carrier = gr_expj(-angle);

  // update DFE based on pilots
  cur_pilot = 1.0;
  for (unsigned int i = 0; i < d_pilot_carriers.size(); i++) {
    gr_complex pilot_sym(cur_pilot, 0.0);
    cur_pilot = -cur_pilot;
    int di = d_pilot_carriers[i];
    gr_complex sigeq = in[di] * carrier * d_dfe[i];
    // FIX THE FOLLOWING STATEMENT
    if (norm(sigeq)> 0.001)
      d_dfe[i] += d_eq_gain * (pilot_sym/sigeq - d_dfe[i]);
  }

  // equalize all data using interpolated dfe and demap into bytes
  unsigned int pilot_index = 0;
  int pilot_carrier_lo = 0;
  int pilot_carrier_hi = d_pilot_carriers[0];
  gr_complex pilot_dfe_lo = d_dfe[0];
  gr_complex pilot_dfe_hi = d_dfe[0];

  /* alternate implementation of lerp */
  float denom = 1.0/(pilot_carrier_hi - pilot_carrier_lo);					// 1.0/(x_b - x_a)

  for (unsigned int i = 0; i < d_data_carriers.size(); i++) {
      int di = d_data_carriers[i];
      if(di > pilot_carrier_hi) {					// 'di' is beyond the current pilot_hi
	  pilot_index++;						// move to the next pilot

	  pilot_carrier_lo = pilot_carrier_hi;				// the new pilot_lo
	  pilot_dfe_lo = pilot_dfe_hi;					// the new pilot_dfe_lo

	 if (pilot_index < d_pilot_carriers.size()) {
	     pilot_carrier_hi = d_pilot_carriers[pilot_index];
	     pilot_dfe_hi = d_dfe[pilot_index];
	 }
	 else {
	       pilot_carrier_hi = d_occupied_carriers;			// cater to the di's which are beyond the last pilot_carrier
           //pilot_carrier_hi = d_data_carriers.size()+d_pilot_carriers.size()-1;
	 }
	 denom = 1.0/(pilot_carrier_hi - pilot_carrier_lo);
      }

      float alpha = float(pilot_carrier_hi - di) * denom;
      //float alpha = float(di - pilot_carrier_lo) * denom;		// (x - x_a)/(x_b - x_a)
      gr_complex dfe = pilot_dfe_hi + alpha * (pilot_dfe_lo - pilot_dfe_hi);
      //gr_complex dfe = pilot_dfe_lo + alpha * (pilot_dfe_hi - pilot_dfe_lo);	// y = y_a + (y_b - y_a) * (x - x_a)/(x_b - x_a) 
      in[di] *= carrier * dfe;
      /*if(factor != NULL) {
         factor[i] = gr_expj(angle);
	 if(norm(dfe) > 0.0f) 
	    factor[i] /= dfe;
      }*/
  }
}




unsigned int digital_ofdm_frame_sink::demapper(gr_complex *in,
					       unsigned char *out, unsigned int data_flag)
{
  unsigned int i=0, bytes_produced=0;


  equalize_interpolate_dfe(in, NULL);


  while(i < d_data_carriers.size()) {
    if(d_nresid > 0) {
      d_partial_byte |= d_resid;
      d_byte_offset += d_nresid;
      d_nresid = 0;
      d_resid = 0;
    }
    
    while((d_byte_offset < 8) && (i < d_data_carriers.size())) {

      int di = d_data_carriers[i];
      unsigned char bits = slicer(in[di]);

      #if REFSNR == 1
      if(data_flag == 1){
          std::complex<double> ref_value = d_sym_position[refsym_indices[ofdm_counter][i]];
          std::complex<double> recv_value = in[di];
          //std::cout<<"recv symbol: "<<recv_value<<" ref symbol: "<<ref_value<<std::endl;
          double dist_value = std::pow(std::abs(recv_value - ref_value),2);
          //printf("i=%d, ofdm_counter=%d, ref_index=%d dist value=%f\n",i,ofdm_counter,refsym_indices[ofdm_counter][i],dist_value);
 
          ref_subcarrier_error[i]+= dist_value;
      }
      #endif

      i++;

      if((8 - d_byte_offset) >= d_nbits) {
	d_partial_byte |= bits << (d_byte_offset);
	d_byte_offset += d_nbits;
      }
      else {
	d_nresid = d_nbits-(8-d_byte_offset);
	int mask = ((1<<(8-d_byte_offset))-1);
	d_partial_byte |= (bits & mask) << d_byte_offset;
	d_resid = bits >> (8-d_byte_offset);
	d_byte_offset += (d_nbits - d_nresid);
      }
      //printf("demod symbol: %.4f + j%.4f   bits: %x   partial_byte: %x   byte_offset: %d   resid: %x   nresid: %d\n", 
      //     in[i-1].real(), in[i-1].imag(), bits, d_partial_byte, d_byte_offset, d_resid, d_nresid);
    }

    if(d_byte_offset == 8) {
      //printf("demod byte: %x \n\n", d_partial_byte);
      out[bytes_produced++] = d_partial_byte;
      d_byte_offset = 0;
      d_partial_byte = 0;
    }
  }

  #if REFSNR == 1
    if(data_flag==1){ ofdm_counter++;}
  #endif
  //std::cerr << "accum_error " << accum_error << std::endl;
   
  //if(VERBOSE)
  //  std::cerr << angle << "\t" << d_freq << "\t" << d_phase << "\t" << std::endl;
  
  return bytes_produced;
}


digital_ofdm_frame_sink_sptr
digital_make_ofdm_frame_sink(const std::vector<gr_complex> &sym_position, 
			     const std::vector<unsigned char> &sym_value_out,
			     gr_msg_queue_sptr target_queue, unsigned int occupied_carriers,
			     float phase_gain, float freq_gain)
{
  return gnuradio::get_initial_sptr(new digital_ofdm_frame_sink(sym_position, sym_value_out,
								target_queue, occupied_carriers,
								phase_gain, freq_gain));
}


digital_ofdm_frame_sink::digital_ofdm_frame_sink(const std::vector<gr_complex> &sym_position, 
						 const std::vector<unsigned char> &sym_value_out,
						 gr_msg_queue_sptr target_queue, unsigned int occupied_carriers,
						 float phase_gain, float freq_gain)
  : gr_sync_block ("ofdm_frame_sink",
		   gr_make_io_signature2 (2, 2, sizeof(gr_complex)*occupied_carriers, sizeof(char)),
		   gr_make_io_signature (1, 1, sizeof(gr_complex)*occupied_carriers)),
    d_target_queue(target_queue), d_occupied_carriers(occupied_carriers), 
    d_byte_offset(0), d_partial_byte(0),
    d_resid(0), d_nresid(0),d_phase(0),d_freq(0),d_phase_gain(phase_gain),d_freq_gain(freq_gain),
    d_eq_gain(0.05)
{
  
  assign_subcarriers();

  #if REFSNR == 1
    read_refsymbols_file();
    ofdm_counter = 0;
    ref_subcarrier_error.assign(d_data_carriers.size(),0);
    outfile_snr.open("ref_snr_qpsk.txt");
  #endif

  if(d_data_carriers.size() > d_occupied_carriers) {
    throw std::invalid_argument("digital_ofdm_mapper_bcv: subcarriers allocated exceeds size of occupied carriers");
  }

  d_bytes_out = new unsigned char[d_occupied_carriers];
  d_dfe.resize(occupied_carriers);
  fill(d_dfe.begin(), d_dfe.end(), gr_complex(1.0,0.0));

  set_sym_value_out(sym_position, sym_value_out);
  
  enter_search();
}

digital_ofdm_frame_sink::~digital_ofdm_frame_sink ()
{
  delete [] d_bytes_out;

  #if REFSNR == 1
    outfile_snr.close();
  #endif
}

#if REFSNR == 1
void
digital_ofdm_frame_sink::read_refsymbols_file(){

    std::ifstream bit_file("ref_symbols_qpsk.txt");
    std::string line;

    if(bit_file.is_open()){
        while (bit_file.good() ){
            getline(bit_file,line);
            
            if(bit_file.eof()){break;}

            std::vector<int> bit_value(d_data_carriers.size());

            for(int i=0; i<d_data_carriers.size(); i++){
                size_t location = line.find(",");
                std::string value = line.substr(0,location);
                bit_value[i] = std::atoi(value.c_str());
                line = line.substr(location+1);
            }
            refsym_indices.push_back(bit_value);
        }
   }

    bit_file.close();
}
#endif

void
digital_ofdm_frame_sink::assign_subcarriers() {
   int dc_tones = 8;
   int num_pilots = 8;
   int pilot_gap = 11;
 
   int half1_end = (d_occupied_carriers-dc_tones)/2;     //40
   int half2_start = half1_end+dc_tones;                 //48
 
   // first half
   for(int i = 0; i < half1_end; i++) {
      if(i%pilot_gap == 0)
         d_pilot_carriers.push_back(i);
      else
         d_data_carriers.push_back(i);
   }
 
   // second half
   for(int i = half2_start, j = 0; i < d_occupied_carriers; i++, j++) {
      if(j%pilot_gap == 0)
         d_pilot_carriers.push_back(i);
      else
         d_data_carriers.push_back(i);
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



bool
digital_ofdm_frame_sink::set_sym_value_out(const std::vector<gr_complex> &sym_position, 
					   const std::vector<unsigned char> &sym_value_out)
{
  if (sym_position.size() != sym_value_out.size())
    return false;

  if (sym_position.size()<1)
    return false;

  d_sym_position  = sym_position;
  d_sym_value_out = sym_value_out;
  d_nbits = (unsigned long)ceil(log10(float(d_sym_value_out.size())) / log10(2.0));

  return true;
}


int
digital_ofdm_frame_sink::work (int noutput_items,
			       gr_vector_const_void_star &input_items,
			       gr_vector_void_star &output_items)
{
  gr_complex *in = (gr_complex *) input_items[0];
  const char *sig = (const char *) input_items[1];
  unsigned int j = 0;
  unsigned int bytes=0;

  // If the output is connected, send it the derotated symbols
  if(output_items.size() >= 1)
    d_derotated_output = (gr_complex *)output_items[0];
  else
    d_derotated_output = NULL;
  
  if (VERBOSE)
    fprintf(stderr,">>> Entering state machine\n");

  switch(d_state) {
      
  case STATE_SYNC_SEARCH:    // Look for flag indicating beginning of pkt
    if (VERBOSE)
      fprintf(stderr,"SYNC Search, noutput=%d\n", noutput_items);
    
    if (sig[0]) {  // Found it, set up for header decode
      enter_have_sync();
    }
    break;

  case STATE_HAVE_SYNC:
    // only demod after getting the preamble signal; otherwise, the 
    // equalizer taps will screw with the PLL performance
    bytes = demapper(&in[0], d_bytes_out,0);
    
    if (VERBOSE) {
      if(sig[0])
	printf("ERROR -- Found SYNC in HAVE_SYNC\n");
      fprintf(stderr,"Header Search bitcnt=%d, header=0x%08x\n",
	      d_headerbytelen_cnt, d_header);
    }

    j = 0;
    while(j < bytes) {
      d_header = (d_header << 8) | (d_bytes_out[j] & 0xFF);
      j++;
      
      if (++d_headerbytelen_cnt == HEADERBYTELEN) {
	
	if (VERBOSE)
	  fprintf(stderr, "got header: 0x%08x\n", d_header);
	
	// we have a full header, check to see if it has been received properly
	if (header_ok()){
	  enter_have_header();
	  
	  if (VERBOSE)
	    printf("\nPacket Length: %d\n", d_packetlen);
	  
	  while((j < bytes) && (d_packetlen_cnt < d_packetlen)) {
	    d_packet[d_packetlen_cnt++] = d_bytes_out[j++];
	  }
	  
	  if(d_packetlen_cnt == d_packetlen) {
	    gr_message_sptr msg =
	      gr_make_message(0, d_packet_whitener_offset, 0, d_packetlen);
	    memcpy(msg->msg(), d_packet, d_packetlen_cnt);
	    d_target_queue->insert_tail(msg);		// send it
	    msg.reset();  				// free it up
	    
	    enter_search();				
	  }
	}
	else {
	  enter_search();				// bad header
	}
      }
    }
    break;
      
  case STATE_HAVE_HEADER:
    bytes = demapper(&in[0], d_bytes_out,1);

    if (VERBOSE) {
      if(sig[0])
	printf("ERROR -- Found SYNC in HAVE_HEADER at %d, length of %d\n", d_packetlen_cnt, d_packetlen);
      fprintf(stderr,"Packet Build\n");
    }
    
    j = 0;
    while(j < bytes) {
      d_packet[d_packetlen_cnt++] = d_bytes_out[j++];
      
      if (d_packetlen_cnt == d_packetlen){		// packet is filled
	// build a message
	// NOTE: passing header field as arg1 is not scalable
	gr_message_sptr msg =
	  gr_make_message(0, d_packet_whitener_offset, 0, d_packetlen_cnt);
	memcpy(msg->msg(), d_packet, d_packetlen_cnt);
	
	d_target_queue->insert_tail(msg);		// send it
	msg.reset();  				// free it up

    #if REFSNR == 1
        std::vector<double> sub_snr;
        calculate_refsnr(sub_snr);
        ofdm_counter=0;
        std::fill(ref_subcarrier_error.begin(),ref_subcarrier_error.end(),0);

        /*
        for(int lp=0; lp<sub_snr.size(); lp++){
            printf("%2.2f, ",sub_snr[lp]);
        }
        printf("\n");
        */
    #endif
	
	enter_search();
	break;
      }
    }
    break;
    
  default:
    assert(0);
    
  } // switch

  return 1;
}


void
digital_ofdm_frame_sink::calculate_refsnr(std::vector<double>& subcarrier_snr_ref){

    std::vector<double>::iterator it;

    for(it = ref_subcarrier_error.begin(); it != ref_subcarrier_error.end(); it++){

        //FIXME: ref_subcarrier_error should be initialized to zero after calculation
        double noise_value = *it/ofdm_counter;
        double snr_value  = 10.0*std::log10(1/noise_value);
        subcarrier_snr_ref.push_back(snr_value);
        outfile_snr<< snr_value<<",";
    }
    outfile_snr<<"\n";
}
