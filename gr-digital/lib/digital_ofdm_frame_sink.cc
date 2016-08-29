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
#include <sys/time.h>
#include <fstream>


#define VERBOSE 0
#define REFSNR 0
#define LTFPREAMBLE 1
#define PARTIALCOMB 1

// GLOBAL VARIABLES

#if REFSNR == 1
unsigned int ofdm_counter;
std::vector<std::vector<int> > refsym_indices;
std::vector<double> ref_subcarrier_error;
std::ofstream outfile_refsnr;
#endif

#if REFSNR == 1
    #define TRACECHAN 1
#else
    #define TRACECHAN 0
#endif


#if TRACECHAN == 1
unsigned int ofdm_counter_trace;
std::vector<std::complex<double> > ref_trace_chan_per_ofdm;
std::ofstream outfile_chan;
#endif


inline void
conv_bvec_to_chardata(itpp::bvec& vecdata, int valid_bits,unsigned char* data_bytes){

    //printf("valid bits: %d",valid_bits);
    int boffset;
    int ioffset;
    for( int i=0; i<valid_bits; i++){
        boffset = i/8;
        ioffset = i%8;
        int val = (int)vecdata[i];
        //printf("boffset=%d val=%d byte=%x\n",boffset,val,data_bytes[boffset]);
        data_bytes[boffset] |=(val&0x1)<<(8-ioffset-1);
    }

}

inline void
digital_ofdm_frame_sink::enter_search()
{
  if (VERBOSE)
    fprintf(stderr, "@ enter_search\n");

  d_state = STATE_SYNC_SEARCH;

  #if LTFPREAMBLE == 1
    preamble_count = 0;
  #endif
}

#if LTFPREAMBLE == 1
void
digital_ofdm_frame_sink::calculate_ltfpreamble_snr(){

    std::complex<double> h1,h2;
    std::complex<double>* h;
    h = new std::complex<double>[d_data_carriers.size()];
    double noise=0;

    for(int i=0; i < d_data_carriers.size(); i++){
        //estimate channel
        h1 = rcvd_ltfpreambles[0][i]*(double)ltf_preamble_rx[i];
        h2 = rcvd_ltfpreambles[1][i]*(double)ltf_preamble_rx[i];
        h[i] = (h1+h2)/2.0;
        noise += std::pow(std::abs(rcvd_ltfpreambles[0][i] - rcvd_ltfpreambles[1][i]),2);
    } 
    noise/= (2*d_data_carriers.size());

    for(int j=0; j < d_data_carriers.size(); j++){
        //ltf_subcarrier_snr[j] = 10.0*log10(std::pow(std::abs(h[j]),2)/noise); 
        ltf_subcarrier_esno[j] = std::pow(std::abs(h[j]),2)/noise; 
        outfile_ltfsnr<<ltf_subcarrier_esno[j]<<",";
    }
    outfile_ltfsnr<<"\n";

    free(h);
}
#endif
    
inline void
create_constellation_type(int bits_per_mod, double real[], double imag[], std::vector<gr_complex>& constel, std::vector<std::vector<int> >& values){

    int number = pow(2,bits_per_mod);
    for(int i=0; i<number; i++){
        constel.push_back(gr_complex(real[i],imag[i]));
        std::vector<int> bits_per_val;
        for(int b=0; b<bits_per_mod; b++){
            int bit_val = (i&(1<<(bits_per_mod-b-1)))>>(bits_per_mod-b-1);
            bits_per_val.push_back(bit_val);
        }
        values.push_back(bits_per_val);
    }
}

void 
digital_ofdm_frame_sink::create_constellations(){

    int bits_per_mod;
    // bpsk
    bpsk_const.push_back(gr_complex(-1.0,0.0));
    bpsk_const.push_back(gr_complex(1.0,0.0));
    bpsk_values.push_back(std::vector<int>(1,0)); // fill vector size=1 with 0
    bpsk_values.push_back(std::vector<int>(1,1));

    // qpsk
    double real_qpsk[] = {-0.707, -0.707,  0.707, 0.707};
    double imag_qpsk[] = {-0.707,  0.707, -0.707, 0.707};

    bits_per_mod = 2;
    create_constellation_type(bits_per_mod,real_qpsk,imag_qpsk, qpsk_const, qpsk_values);

/*    for(int i=0; i<qpsk_const.size();i++){
        printf("(%f,%f) ",qpsk_const[i].real(),qpsk_const[i].imag());
        for(int j=0; j<qpsk_values[i].size(); j++){printf("%d ",(int)(qpsk_values[i])[j]);}
    printf("\n");
    }
*/
   // qam16
    bits_per_mod = 4; 
    double real_qam16[] = { -0.9487,-0.9487,-0.9487,-0.9487,-0.3162,-0.3162,-0.3162,-0.3162, 0.9487, 0.9487, 0.9487, 0.9487, 0.3162, 0.3162, 0.3162, 0.3162};
    double imag_qam16[] = { -0.9487,-0.3162, 0.9487, 0.3162,-0.9487,-0.3162, 0.9487, 0.3162,-0.9487,-0.3162, 0.9487, 0.3162,-0.9487,-0.3162, 0.9487, 0.3162};
    create_constellation_type(bits_per_mod,real_qam16,imag_qam16, qam16_const, qam16_values);


    // qam64
    bits_per_mod = 6;
    double real_qam64[] = {-1.0801,-1.0801,-1.0801,-1.0801,-1.0801,-1.0801,-1.0801,-1.0801,-0.7715,-0.7715,-0.7715,-0.7715,-0.7715,-0.7715,-0.7715,-0.7715,-0.1543,-0.1543,-0.1543,-0.1543,-0.1543,-0.1543,-0.1543,-0.1543,-0.4629,-0.4629,-0.4629,-0.4629,-0.4629,-0.4629,-0.4629,-0.4629, 1.0801, 1.0801, 1.0801, 1.0801, 1.0801, 1.0801, 1.0801, 1.0801, 0.7715, 0.7715, 0.7715, 0.7715, 0.7715, 0.7715, 0.7715, 0.7715, 0.1543, 0.1543, 0.1543, 0.1543, 0.1543, 0.1543, 0.1543, 0.1543, 0.4629, 0.4629, 0.4629, 0.4629, 0.4629, 0.4629, 0.4629, 0.4629};

    double imag_qam64[] = {-1.0801,-0.7715,-0.1543,-0.4629, 1.0801, 0.7715, 0.1543, 0.4629,-1.0801,-0.7715,-0.1543,-0.4629, 1.0801, 0.7715, 0.1543, 0.4629,-1.0801,-0.7715,-0.1543,-0.4629, 1.0801, 0.7715, 0.1543, 0.4629,-1.0801,-0.7715,-0.1543,-0.4629, 1.0801, 0.7715, 0.1543, 0.4629,-1.0801,-0.7715,-0.1543,-0.4629, 1.0801, 0.7715, 0.1543, 0.4629,-1.0801,-0.7715,-0.1543,-0.4629, 1.0801, 0.7715, 0.1543, 0.4629,-1.0801,-0.7715,-0.1543,-0.4629, 1.0801, 0.7715, 0.1543, 0.4629,-1.0801,-0.7715,-0.1543,-0.4629, 1.0801, 0.7715, 0.1543, 0.4629};

    create_constellation_type(bits_per_mod,real_qam64,imag_qam64, qam64_const,qam64_values);
}

//
void
digital_ofdm_frame_sink::send_mcs_feedback(int mcs){

    unsigned char feedback_data = (unsigned char) mcs;
    printf("mcs=%d\n",mcs);
    socklen_t fromlen = sizeof(struct sockaddr_in);
    int n = sendto(sock,&feedback_data,1,0,(struct sockaddr *)&from,fromlen);
}

//
void 
digital_ofdm_frame_sink::open_feedback_socket(){

    int portno = 51717;
    char buf[1024];

    sock=socket(AF_INET, SOCK_DGRAM, 0);
    if (sock < 0){ fprintf(stderr,"Opening socket"); exit(0);}
    int length = sizeof(server);
    bzero(&server,length);
    server.sin_family=AF_INET;
    server.sin_addr.s_addr=INADDR_ANY;
    server.sin_port=htons(portno);
    if (bind(sock,(struct sockaddr *)&server,length)<0){
        fprintf(stderr,"Error binding\n");
        exit(0);
    }
    socklen_t fromlen = sizeof(struct sockaddr_in);
    int n = recvfrom(sock,buf,1024,0,(struct sockaddr *)&from,&fromlen);
    if (n < 0){ fprintf(stderr,"Error recvfrom\n");exit(0);}
    else{ printf("data received.");}

    n = sendto(sock,"Got your message\n",18,0,(struct sockaddr *)&from,fromlen);

}

/*
*/
inline void
digital_ofdm_frame_sink::enter_have_sync(){

    if (VERBOSE)
        fprintf(stderr, "@ enter_have_sync\n");

    #if LTFPREAMBLE == 1
        d_state = STATE_HAVE_SYNC;
    #else
        d_state = STATE_HAVE_LTFPREAMBLE;
    #endif


    d_state = STATE_HAVE_SYNC;

    // Resetting PLL
    d_freq = 0.0;
    d_phase = 0.0;
    fill(d_dfe.begin(), d_dfe.end(), gr_complex(1.0,0.0));
}

//
inline 
unsigned int dec_convhelper(itpp::bvec& hdr, int offset, int len){

    unsigned int value=0;
    for(int i =0; i<len; i++){
        int tmp =  (int) hdr[offset+i];
        value |= (((unsigned int)tmp)&0x1)<<(len-i-1);
    }
    return value;
}

bool 
digital_ofdm_frame_sink::enter_have_header(){
    
    itpp::bvec hdr_decode;

    fec_decoder.decode(0,hdr_demod,hdr_decode);

    header.pktlength = (unsigned short) dec_convhelper(hdr_decode,0,12);
    header.mcs = (unsigned short) dec_convhelper(hdr_decode,12,4);
    header.pktid = (unsigned short) dec_convhelper(hdr_decode,16,16);
    header.whitener = (unsigned short) dec_convhelper(hdr_decode,32,4);
    header.retx_count = (unsigned short) dec_convhelper(hdr_decode,36,4);
    header.payloadlen = (unsigned short) dec_convhelper(hdr_decode,40,16);
    header.hdr_crc = (unsigned int) dec_convhelper(hdr_decode,56,32);

    printf("pktlen=%u, mcs=%u, pktid=%u, whitener=%u,retx=%u, payloadlen=%u, crc=%u\n",
header.pktlength,header.mcs,header.pktid,header.whitener,header.retx_count,header.payloadlen, header.hdr_crc);

    unsigned char* header_data = (unsigned char*) &header;
    int len_minus_crc = sizeof(PHY_HEADER)-4;                         
    unsigned int hdr_crc = digital_crc32(header_data,len_minus_crc);

    if(header.hdr_crc == hdr_crc){
        return true;
    }
    else{
        return false;
    } 
}



//

void 
digital_ofdm_frame_sink::llrslicer(const gr_complex x, double noise, int bits_per_mod, std::vector<std::vector<gr_complex> >& S0_table,std::vector<std::vector<gr_complex> >& S1_table, std::vector<double>& output){

    unsigned int table_size = S0_table[0].size();

    //printf("noise: %f\n",noise);
    //printf("table_size: %d, bits_per_mod: %d\n",table_size,bits_per_mod);
    for (unsigned int i=0; i < bits_per_mod; i++){
        double sumS0=0,sumS1=0,distS0, distS1;
        for (unsigned int j=0; j < table_size; j++){
            distS0 = norm(x - S0_table[i][j])/noise;
            distS1 = norm(x - S1_table[i][j])/noise;
            //printf("dist: %f, %f\n",distS0,distS1);
            sumS0 += exp(-1*distS0);
            sumS1 += exp(-1*distS1);
        }
        double tmp = log(sumS0/sumS1);
        //printf("%f, ",tmp);
        tmp = std::min(tmp,1000.0);
        output[i] = std::max(tmp,-1000.0);
        //output[i] = log(sumS0/sumS1);
    }
//printf("\n");
}


//
void 
digital_ofdm_frame_sink::bitslicer(const gr_complex x,std::vector<gr_complex>& constellation, std::vector<std::vector<int> >& constel_value, std::vector<int>& output)
{
    unsigned int table_size = constellation.size();
    unsigned int min_index = 0;
    float min_euclid_dist = norm(x - constellation[0]);
    float euclid_dist = 0;
  
    for (unsigned int j = 1; j < table_size; j++){
        euclid_dist = norm(x - constellation[j]);
        if (euclid_dist < min_euclid_dist){
            min_euclid_dist = euclid_dist;
            min_index = j;
        }
    }
    output = constel_value[min_index];
  
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

void
digital_ofdm_frame_sink::get_constellations(int bits_per_mod, std::vector<gr_complex>& constel, std::vector<std::vector<int> >& constel_values){
    if (bits_per_mod == 1){
        constel = bpsk_const;
        constel_values = bpsk_values;
    }
    else if (bits_per_mod == 2){
        constel = qpsk_const;
        constel_values = qpsk_values;
    }
    else if (bits_per_mod == 4){
        constel = qam16_const;
        constel_values = qam16_values;
    }
    else if (bits_per_mod == 6){
        constel = qam64_const;
        constel_values = qam64_values;
    }
    else{
        printf("Incorrect bits_per_mod value.\n");
    }
}

//
void 
digital_ofdm_frame_sink::set_currllr_table(int bits_per_mod){

    if (bits_per_mod == 1){
        currS0_table = S0_table_bpsk; currS1_table = S1_table_bpsk; 
        printf("BPSK set.");
    }
    else if (bits_per_mod == 2){
        currS0_table = S0_table_qpsk; currS1_table = S1_table_qpsk; 
        printf("QPSK set.");
    }
    else if (bits_per_mod == 4){
        currS0_table = S0_table_qam16; currS1_table = S1_table_qam16; 
        printf("QAM16 set.");
    }
    else if (bits_per_mod == 6){
        currS0_table = S0_table_qam64; currS1_table = S1_table_qam64; 
        printf("QAM64 set.");
    }
}


//
void 
digital_ofdm_frame_sink::llr_demapper(gr_complex *in, itpp::vec& out, unsigned int& offset_value, unsigned int bits_per_mod){

    
    //FIXME: should be turned when only llr is used
    //equalize_interpolate_dfe(in, NULL);
    unsigned int i=0;
    std::vector<double> llr_values(bits_per_mod);

    while( i < d_data_carriers.size() ) {

        int di = d_data_carriers[i];
        double noise = 1/ltf_subcarrier_esno[i];
        llrslicer(in[di], noise, bits_per_mod, currS0_table, currS1_table, llr_values);

        for(int b=0; b<bits_per_mod; b++){
            out[offset_value] = llr_values[b];
            offset_value++;
        }
        i++;
    }

}


//
void 
digital_ofdm_frame_sink::bit_demapper(gr_complex *in, itpp::vec& out, unsigned int& offset_value, unsigned int bits_per_mod){

    equalize_interpolate_dfe(in, NULL);
    unsigned int i=0;

    std::vector<int> bit_values;
    std::vector<gr_complex> constel;
    std::vector< std::vector<int> > constel_values;
    get_constellations(bits_per_mod, constel,constel_values);

    while( i < d_data_carriers.size() ) {

      int di = d_data_carriers[i];
      bitslicer(in[di],constel,constel_values,bit_values);


      for(int b=0; b<bits_per_mod; b++){
        out[offset_value] = bit_values[b];
        offset_value++;
      }
      i++;
    }


}

/*
*/
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
    d_phase(0),d_freq(0),d_phase_gain(phase_gain),d_freq_gain(freq_gain),
    d_eq_gain(0.05),
    total_pkts(0),
    succ_pkts(0),
    demodtime_elapsed(0.0),
    fec_decoder(0)
{
  
  assign_subcarriers();

  #if REFSNR == 1
    read_refsymbols_file();
    ofdm_counter = 0;
    ref_subcarrier_error.assign(d_data_carriers.size(),0);
    outfile_refsnr.open("ref_snr_bpsk.txt");
  #endif

  #if TRACECHAN == 1
    ofdm_counter_trace=0;
    outfile_chan.open("trace_channel.txt");
    ref_trace_chan_per_ofdm = std::vector<std::complex<double> > (d_data_carriers.size()); 
  #endif

  if(d_data_carriers.size() > d_occupied_carriers) {
    throw std::invalid_argument("digital_ofdm_mapper_bcv: subcarriers allocated exceeds size of occupied carriers");
  }


  #if LTFPREAMBLE == 1
    int ltf_preamble96[] = {1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1, -1, -1,1,1, -1,1, -1,1, -1, -1, -1, -1, -1,1,1, -1, -1,1, -1,1,-1,0,0,0,0,0,0,0,0,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1, -1, -1,1,1, -1,1, -1,1, -1, -1, -1, -1, -1,1,1, -1, -1,1, -1,1,-1};

    ltf_preamble_rx = new int[d_data_carriers.size()];
    for(int ltf=0; ltf < d_data_carriers.size(); ltf++){
        ltf_preamble_rx[ltf] = ltf_preamble96[d_data_carriers[ltf]];
    }

    int d_ltfpreamble_size = d_data_carriers.size();
    rcvd_ltfpreambles = new std::complex<double>*[2];
    rcvd_ltfpreambles[0] = new std::complex<double>[d_ltfpreamble_size];
    rcvd_ltfpreambles[1] = new std::complex<double>[d_ltfpreamble_size];

    ltf_subcarrier_esno = std::vector<double>(d_data_carriers.size());
    outfile_ltfsnr.open("ltf_snr_bpsk.txt");
    h_estimates = std::vector<double>(d_data_carriers.size(),0);
  #endif

  d_dfe.resize(occupied_carriers);
  fill(d_dfe.begin(), d_dfe.end(), gr_complex(1.0,0.0));

  // hdr variables
  int num_ofdm = ceil(1.0*(sizeof(PHY_HEADER)*8*2 + 12)/d_data_carriers.size());
  coded_hdr_bitlen = num_ofdm*d_data_carriers.size();
  printf("coded_hdr_bitlen=%d\n",coded_hdr_bitlen);
  fflush(stdout);
  d_hdr_offset = 0;
  hdr_demod.set_size(coded_hdr_bitlen);

  create_constellations();
  create_llr_tables();

  // create deinterleaver
  deinterleaver = digital_ofdm_deinterleaver (d_data_carriers.size());

  // rate adaptation
  ra_effsnr = digital_ofdm_effsnr_ra (d_data_carriers.size());
 
  open_feedback_socket();
  enter_search();
}

//
digital_ofdm_frame_sink::~digital_ofdm_frame_sink ()
{
  //delete [] d_bytes_out;
  double dratio = 1.0*succ_pkts/total_pkts;
  printf("total pkts=%d, succ=%d, dratio=%f\n",total_pkts,succ_pkts,dratio);
  fflush(stdout);

  #if REFSNR == 1
    outfile_refsnr.close();
  #endif

  #if TRACECHAN == 1
    outfile_chan.close();
  #endif

  #if LTFPREAMBLE == 1
    free(ltf_preamble_rx);
    free(rcvd_ltfpreambles[0]);
    free(rcvd_ltfpreambles[1]);
    free(rcvd_ltfpreambles);
    outfile_ltfsnr.close();
  #endif
}

#if REFSNR == 1
void
digital_ofdm_frame_sink::read_refsymbols_file(){

    std::ifstream bit_file("ref_symbols_bpsk.txt");
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


/*
*/
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

    if(output_items.size() >= 1){
        d_derotated_output = (gr_complex *)output_items[0];
    }
    else{
        d_derotated_output = NULL;
    }
  
    if (VERBOSE){ fprintf(stderr,">>> Entering state machine\n");}

    switch(d_state) {
      
        case STATE_SYNC_SEARCH:    // Look for flag indicating beginning of pkt
            
            if (VERBOSE)
              fprintf(stderr,"SYNC Search, noutput=%d\n", noutput_items);
   
           // Found it, set up for header decode
            if (sig[0]) { 
                std::vector<gr_tag_t> tags;
                const uint64_t nread = this->nitems_read(0);
                const size_t ninput_items = noutput_items;
                this->get_tags_in_range(tags, 0, nread, nread+ninput_items);

                pmt::pmt_t pmt_vec = tags[0].value;
                size_t vec_size;
                double *elements = pmt_f64vector_writable_elements(pmt_vec, vec_size); 

                for(int i=0; i<d_data_carriers.size();i++){
                    h_estimates[i] = elements[d_data_carriers[i]];
                }
                //for(int i=0; i<vec_size; i++){printf("%f, ",elements[i]);}
                //printf("\n");

                //for(int i=0; i<tags.size();i++){ printf("%f",tags[i]);}
                enter_have_sync();
             }
            break;

        case STATE_HAVE_SYNC:
            #if LTFPREAMBLE==1
                for (int pl = 0; pl < d_data_carriers.size(); pl++){
                    rcvd_ltfpreambles[preamble_count][pl] = in[d_data_carriers[pl]];
                    rcvd_ltfpreambles[preamble_count][pl]*=h_estimates[pl];
                }
                preamble_count += 1;
                if(preamble_count == 2){
                    calculate_ltfpreamble_snr();
                    d_state = STATE_HAVE_LTFPREAMBLE;
                    preamble_count = 0;
                }
                break;
            #endif

        case STATE_HAVE_LTFPREAMBLE:
            // only demod after getting the preamble signal; otherwise, the 
            // equalizer taps will screw with the PLL performance
            bit_demapper(&in[0], hdr_demod, d_hdr_offset,1);


            if (d_hdr_offset==coded_hdr_bitlen) { //start if header-read
                d_hdr_offset = 0;
                d_data_offset= 0;
                llr_data_offset= 0;
                bool hdr_valid = enter_have_header();
                if(hdr_valid){
                    d_state = STATE_HAVE_HEADER;
                    data_demod.set_size(header.payloadlen);
                    llr_demod.set_size(header.payloadlen);

                    int bits_per_mod = MCS2BITSMOD[header.mcs];
                    set_currllr_table(bits_per_mod);
                }            
                else{
    	            enter_search();
                }
            }
            break;

        case STATE_HAVE_HEADER:

            int bits_per_mod = MCS2BITSMOD[header.mcs];

            double t0 = current_time();
            bit_demapper(&in[0], data_demod, d_data_offset,bits_per_mod);
            demodtime_elapsed += (current_time() - t0);
            //llr_demapper(&in[0], llr_demod, llr_data_offset, bits_per_mod);



            if (d_data_offset == data_demod.size()){

                printf("demod time=%f\n",1.0e3*demodtime_elapsed);
                demodtime_elapsed = 0.0;

                //std::cout<<"bit_demapper: "<<data_demod<<std::endl;
                //std::cout<<"llr_demapper: "<<llr_demod<<std::endl;
/*
                for(int i=0; i<llr_demod.size();i++){
                    int diffval = 1;
                    if(llr_demod[i]>0){ diffval =0;}
                    int dval = (int)data_demod[i];
                    printf("%d, ",(diffval-dval));
                }
                printf("\n");
*/              

/*
                #if LTFPREAMBLE == 1
	                enter_have_ltfpreamble(&h_in[0]);
                #endif
*/
                // de-interleave data
                itpp::vec data_deinlv;
                deinterleaver.deinterleave_data(data_demod,data_deinlv,header.mcs);

                itpp::bvec data_decode;
                unsigned int code_rate = MCS2CODERATE[header.mcs];
                fec_decoder.decode(code_rate,data_deinlv,data_decode);
                //fec_decoder.decode(code_rate,data_demod,data_decode);

                unsigned int crc_offset = header.pktlength*8;
                unsigned int pktdata_crc = (unsigned int) dec_convhelper(data_decode,crc_offset,32);
                printf("pkt crc=%x\n",pktdata_crc);

                unsigned int dec_bytes = header.pktlength;
                unsigned char* dec_databytes = new unsigned char[dec_bytes]();
                int valid_bits = header.pktlength*8;
                conv_bvec_to_chardata(data_decode,valid_bits,dec_databytes);

                unsigned int data_crc = digital_crc32(dec_databytes,dec_bytes);

                //printf("data bytes=%d\n",dec_bytes);
                //for(int i=0; i<dec_bytes;i++){ printf("%x",dec_databytes);}
                //printf("\n");
                printf("data crc=%x\n",data_crc);

                //FIXME: Hard coding pkt size for now. Fix it.
                int ra_mcs = ra_effsnr.select_rate(ltf_subcarrier_esno,1000);
               
                if(pktdata_crc==data_crc){
                    printf("pkt succ\n");
                    succ_pkts++;
                    //send_mcs_feedback(std::min((int)(++header.mcs),7));
                }
                else{
                    printf("pkt fail\n");
                    //send_mcs_feedback(std::max((int)(--header.mcs),0));
                }
                total_pkts++;
                send_mcs_feedback(ra_mcs);

                d_data_offset=0;
                llr_data_offset=0;
                free(dec_databytes);
    	        enter_search();// bad header
            }

            break;
                 
//        default:
//            assert(0);
 
    }

    return 1; 
}

#if REFSNR == 1
void
digital_ofdm_frame_sink::calculate_refsnr(std::vector<double>& subcarrier_snr_ref){

    std::vector<double>::iterator it;

    for(it = ref_subcarrier_error.begin(); it != ref_subcarrier_error.end(); it++){

        double noise_value = *it/ofdm_counter;
        double snr_value  = 10.0*std::log10(1/noise_value);
        subcarrier_snr_ref.push_back(snr_value);
        outfile_refsnr<< snr_value<<",";
    }
    outfile_refsnr<<"\n";
}
#endif


inline void 
create_complex_array(int length, double real[], double imag[], std::vector<gr_complex>& out){

    for(int i=0; i<length; i++){
        out.push_back(gr_complex(real[i],imag[i]));
    }
}

//
void 
digital_ofdm_frame_sink::create_llr_tables(){

    // bpsk
    int length = 1;
    double S0_table_bpsk_b0_real[] = {-1};
    double S0_table_bpsk_b0_imag[] = {0};
    double S1_table_bpsk_b0_real[] = {1};
    double S1_table_bpsk_b0_imag[] = {0};
    std::vector<gr_complex> S0_table_bpsk_b0;
    std::vector<gr_complex> S1_table_bpsk_b0;
    create_complex_array(length, S0_table_bpsk_b0_real,S0_table_bpsk_b0_imag, S0_table_bpsk_b0);
    create_complex_array(length, S1_table_bpsk_b0_real,S1_table_bpsk_b0_imag, S1_table_bpsk_b0);
    S0_table_bpsk.push_back(S0_table_bpsk_b0);
    S1_table_bpsk.push_back(S1_table_bpsk_b0);


    // qpsk
    double S0_table_qpsk_b0_real[] = {-0.70711,-0.70711};
    double S0_table_qpsk_b0_imag[] = {0.70711 ,-0.70711};
    double S1_table_qpsk_b0_real[] = {0.70711 , 0.70711};
    double S1_table_qpsk_b0_imag[] = {0.70711 ,-0.70711};
    double S0_table_qpsk_b1_real[] = {-0.70711, 0.70711};
    double S0_table_qpsk_b1_imag[] = {-0.70711,-0.70711};
    double S1_table_qpsk_b1_real[] = {0.70711 ,-0.70711};
    double S1_table_qpsk_b1_imag[] = {0.70711 , 0.70711};

    length = 2;
    std::vector<gr_complex> S0_table_qpsk_b0;
    std::vector<gr_complex> S0_table_qpsk_b1;
    std::vector<gr_complex> S1_table_qpsk_b0;
    std::vector<gr_complex> S1_table_qpsk_b1;
    create_complex_array(length, S0_table_qpsk_b0_real,S0_table_qpsk_b0_imag, S0_table_qpsk_b0);
    create_complex_array(length, S0_table_qpsk_b1_real,S0_table_qpsk_b1_imag, S0_table_qpsk_b1);
    create_complex_array(length, S1_table_qpsk_b0_real,S1_table_qpsk_b0_imag, S1_table_qpsk_b0);
    create_complex_array(length, S1_table_qpsk_b1_real,S1_table_qpsk_b1_imag, S1_table_qpsk_b1);
    S0_table_qpsk.push_back(S0_table_qpsk_b0);
    S0_table_qpsk.push_back(S0_table_qpsk_b1);
    S1_table_qpsk.push_back(S1_table_qpsk_b0);
    S1_table_qpsk.push_back(S1_table_qpsk_b1);



    // qam16
    length = 8;
    double S0_table_qam16_b0_real[] = {-0.94868,-0.94868,-0.94868,-0.94868,-0.31623,-0.31623,-0.31623,-0.31623};
    double S0_table_qam16_b0_imag[] = {-0.94868,-0.31623, 0.94868, 0.31623,-0.94868,-0.31623, 0.94868, 0.31623};
    double S1_table_qam16_b0_real[] = { 0.94868, 0.94868, 0.94868, 0.94868, 0.31623, 0.31623, 0.31623, 0.31623};
    double S1_table_qam16_b0_imag[] = {-0.94868,-0.31623, 0.94868, 0.31623,-0.94868,-0.31623, 0.94868, 0.31623};
    double S0_table_qam16_b1_real[] = {-0.94868,-0.94868,-0.94868,-0.94868, 0.94868, 0.94868, 0.94868, 0.94868};
    double S0_table_qam16_b1_imag[] = {-0.94868,-0.31623, 0.94868, 0.31623,-0.94868,-0.31623, 0.94868, 0.31623};
    double S1_table_qam16_b1_real[] = {-0.31623,-0.31623,-0.31623,-0.31623, 0.31623, 0.31623, 0.31623, 0.31623};
    double S1_table_qam16_b1_imag[] = {-0.94868,-0.31623, 0.94868, 0.31623,-0.94868,-0.31623, 0.94868, 0.31623};
    double S0_table_qam16_b2_real[] = {-0.94868,-0.94868,-0.31623,-0.31623, 0.94868, 0.94868, 0.31623, 0.31623};
    double S0_table_qam16_b2_imag[] = {-0.94868,-0.31623,-0.94868,-0.31623,-0.94868,-0.31623,-0.94868,-0.31623};
    double S1_table_qam16_b2_real[] = {-0.94868,-0.94868,-0.31623,-0.31623, 0.94868, 0.94868, 0.31623, 0.31623};
    double S1_table_qam16_b2_imag[] = { 0.94868, 0.31623, 0.94868, 0.31623, 0.94868, 0.31623, 0.94868, 0.31623};
    double S0_table_qam16_b3_real[] = {-0.94868,-0.94868,-0.31623,-0.31623, 0.94868, 0.94868, 0.31623, 0.31623};
    double S0_table_qam16_b3_imag[] = {-0.94868, 0.94868,-0.94868, 0.94868,-0.94868, 0.94868,-0.94868, 0.94868};
    double S1_table_qam16_b3_real[] = {-0.94868,-0.94868,-0.31623,-0.31623, 0.94868, 0.94868, 0.31623, 0.31623};
    double S1_table_qam16_b3_imag[] = {-0.31623, 0.31623,-0.31623, 0.31623,-0.31623, 0.31623,-0.31623, 0.31623};

    std::vector<gr_complex> S0_table_qam16_b0;
    std::vector<gr_complex> S0_table_qam16_b1;
    std::vector<gr_complex> S0_table_qam16_b2;
    std::vector<gr_complex> S0_table_qam16_b3;
    std::vector<gr_complex> S1_table_qam16_b0;
    std::vector<gr_complex> S1_table_qam16_b1;
    std::vector<gr_complex> S1_table_qam16_b2;
    std::vector<gr_complex> S1_table_qam16_b3;
    create_complex_array(length, S0_table_qam16_b0_real,S0_table_qam16_b0_imag, S0_table_qam16_b0);
    create_complex_array(length, S0_table_qam16_b1_real,S0_table_qam16_b1_imag, S0_table_qam16_b1);
    create_complex_array(length, S0_table_qam16_b2_real,S0_table_qam16_b2_imag, S0_table_qam16_b2);
    create_complex_array(length, S0_table_qam16_b3_real,S0_table_qam16_b3_imag, S0_table_qam16_b3);
    create_complex_array(length, S1_table_qam16_b0_real,S1_table_qam16_b0_imag, S1_table_qam16_b0);
    create_complex_array(length, S1_table_qam16_b1_real,S1_table_qam16_b1_imag, S1_table_qam16_b1);
    create_complex_array(length, S1_table_qam16_b2_real,S1_table_qam16_b2_imag, S1_table_qam16_b2);
    create_complex_array(length, S1_table_qam16_b3_real,S1_table_qam16_b3_imag, S1_table_qam16_b3);
    S0_table_qam16.push_back(S0_table_qam16_b0);
    S0_table_qam16.push_back(S0_table_qam16_b1);
    S0_table_qam16.push_back(S0_table_qam16_b2);
    S0_table_qam16.push_back(S0_table_qam16_b3);
    S1_table_qam16.push_back(S1_table_qam16_b0);
    S1_table_qam16.push_back(S1_table_qam16_b1);
    S1_table_qam16.push_back(S1_table_qam16_b2);
    S1_table_qam16.push_back(S1_table_qam16_b3);



    //qam64
    double S0_table_qam64_b0_real[] = {-1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -0.7715, -0.7715,-0.7715, -0.7715, -0.7715, -0.7715, -0.7715, -0.7715, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629};
    double S0_table_qam64_b0_imag[] = {-1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543, 0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629};
    double S1_table_qam64_b0_real[] = {1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  0.7715,  0.7715, 0.7715,  0.7715,  0.7715,  0.7715,  0.7715,  0.7715,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629};
    double S1_table_qam64_b0_imag[] = {-1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715,-0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629};
    double S0_table_qam64_b1_real[] = {-1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -0.7715, -0.7715,-0.7715, -0.7715, -0.7715, -0.7715, -0.7715, -0.7715,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  0.7715,  0.7715,  0.7715,  0.7715,  0.7715,  0.7715,  0.7715,  0.7715};
    double S0_table_qam64_b1_imag[] = {-1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715,-0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629};
    double S1_table_qam64_b1_real[] = {-0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629};
    double S1_table_qam64_b1_imag[] = {-1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629};
    double S0_table_qam64_b2_real[] = {-1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543};
    double S0_table_qam64_b2_imag[] = {-1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629};
    double S1_table_qam64_b2_real[] = {-0.7715, -0.7715, -0.7715, -0.7715, -0.7715, -0.7715, -0.7715, -0.7715, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629,  0.7715,  0.7715,  0.7715,  0.7715,  0.7715,  0.7715,  0.7715,  0.7715,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629};
    double S1_table_qam64_b2_imag[] = {-1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629};
    double S0_table_qam64_b3_real[] = {-1.0801, -1.0801, -1.0801, -1.0801, -0.7715, -0.7715, -0.7715, -0.7715, -0.1543, -0.1543, -0.1543, -0.1543, -0.4629, -0.4629, -0.4629, -0.4629,  1.0801,  1.0801,  1.0801,  1.0801,  0.7715,  0.7715,  0.7715,  0.7715,  0.1543,  0.1543,  0.1543,  0.1543,  0.4629,  0.4629,  0.4629,  0.4629};
    double S0_table_qam64_b3_imag[] = {-1.0801, -0.7715, -0.1543, -0.4629, -1.0801, -0.7715, -0.1543, -0.4629, -1.0801, -0.7715, -0.1543, -0.4629, -1.0801, -0.7715, -0.1543, -0.4629, -1.0801, -0.7715, -0.1543, -0.4629, -1.0801, -0.7715, -0.1543, -0.4629, -1.0801, -0.7715, -0.1543, -0.4629, -1.0801, -0.7715, -0.1543, -0.4629};
    double S1_table_qam64_b3_real[] = {-1.0801, -1.0801, -1.0801, -1.0801, -0.7715, -0.7715, -0.7715, -0.7715, -0.1543, -0.1543, -0.1543, -0.1543, -0.4629, -0.4629, -0.4629, -0.4629,  1.0801,  1.0801,  1.0801,  1.0801,  0.7715,  0.7715,  0.7715,  0.7715,  0.1543,  0.1543,  0.1543,  0.1543,  0.4629,  0.4629,  0.4629,  0.4629};
    double S1_table_qam64_b3_imag[] = {1.0801,  0.7715,  0.1543,  0.4629,  1.0801,  0.7715,  0.1543,  0.4629,  1.0801,  0.7715,  0.1543,  0.4629,  1.0801,  0.7715,  0.1543,  0.4629,  1.0801,  0.7715,  0.1543,  0.4629,  1.0801,  0.7715,  0.1543,  0.4629,  1.0801,  0.7715,  0.1543,  0.4629,  1.0801,  0.7715,  0.1543,  0.4629};
    double S0_table_qam64_b4_real[] = {-1.0801, -1.0801, -1.0801, -1.0801, -0.7715, -0.7715, -0.7715, -0.7715, -0.1543, -0.1543, -0.1543, -0.1543, -0.4629, -0.4629, -0.4629, -0.4629,  1.0801,  1.0801,  1.0801,  1.0801,  0.7715,  0.7715,  0.7715,  0.7715,  0.1543,  0.1543,  0.1543,  0.1543,  0.4629,  0.4629,  0.4629,  0.4629};
    double S0_table_qam64_b4_imag[] = {-1.0801, -0.7715,  1.0801,  0.7715, -1.0801, -0.7715,  1.0801,  0.7715, -1.0801, -0.7715,  1.0801,  0.7715, -1.0801, -0.7715,  1.0801,  0.7715, -1.0801, -0.7715,  1.0801,  0.7715, -1.0801, -0.7715,  1.0801,  0.7715, -1.0801, -0.7715,  1.0801,  0.7715, -1.0801, -0.7715,  1.0801,  0.7715};
    double S1_table_qam64_b4_real[] = {-1.0801, -1.0801, -1.0801, -1.0801, -0.7715, -0.7715, -0.7715, -0.7715, -0.1543, -0.1543, -0.1543, -0.1543, -0.4629, -0.4629, -0.4629, -0.4629,  1.0801,  1.0801,  1.0801,  1.0801,  0.7715,  0.7715,  0.7715,  0.7715,  0.1543,  0.1543,  0.1543,  0.1543,  0.4629,  0.4629,  0.4629,  0.4629};
    double S1_table_qam64_b4_imag[] = {-0.1543, -0.4629,  0.1543,  0.4629, -0.1543, -0.4629,  0.1543,  0.4629, -0.1543, -0.4629,  0.1543,  0.4629, -0.1543, -0.4629,  0.1543,  0.4629, -0.1543, -0.4629,  0.1543,  0.4629, -0.1543, -0.4629,  0.1543,  0.4629, -0.1543, -0.4629,  0.1543,  0.4629, -0.1543, -0.4629,  0.1543,  0.4629};
    double S0_table_qam64_b5_real[] = {-1.0801, -1.0801, -1.0801, -1.0801, -0.7715, -0.7715, -0.7715, -0.7715, -0.1543, -0.1543, -0.1543, -0.1543, -0.4629, -0.4629, -0.4629, -0.4629,  1.0801,  1.0801,  1.0801,  1.0801,  0.7715,  0.7715,  0.7715,  0.7715,  0.1543,  0.1543,  0.1543,  0.1543,  0.4629,  0.4629,  0.4629,  0.4629};
    double S0_table_qam64_b5_imag[] = {-1.0801, -0.1543,  1.0801,  0.1543, -1.0801, -0.1543,  1.0801,  0.1543, -1.0801, -0.1543,  1.0801,  0.1543, -1.0801, -0.1543,  1.0801,  0.1543, -1.0801, -0.1543,  1.0801,  0.1543, -1.0801, -0.1543,  1.0801,  0.1543, -1.0801, -0.1543,  1.0801,  0.1543, -1.0801, -0.1543,  1.0801,  0.1543};
    double S1_table_qam64_b5_real[] = {-1.0801, -1.0801, -1.0801, -1.0801, -0.7715, -0.7715, -0.7715, -0.7715, -0.1543, -0.1543, -0.1543, -0.1543, -0.4629, -0.4629, -0.4629, -0.4629,  1.0801,  1.0801,  1.0801,  1.0801,  0.7715,  0.7715,  0.7715,  0.7715,  0.1543,  0.1543,  0.1543,  0.1543,  0.4629,  0.4629,  0.4629,  0.4629};
    double S1_table_qam64_b5_imag[] = {-0.7715, -0.4629,  0.7715,  0.4629, -0.7715, -0.4629,  0.7715,  0.4629, -0.7715, -0.4629,  0.7715,  0.4629, -0.7715, -0.4629,  0.7715,  0.4629, -0.7715, -0.4629,  0.7715,  0.4629, -0.7715, -0.4629,  0.7715,  0.4629, -0.7715, -0.4629,  0.7715,  0.4629, -0.7715, -0.4629,  0.7715,  0.4629};

    length = 32;
    std::vector<gr_complex> S0_table_qam64_b0;
    std::vector<gr_complex> S0_table_qam64_b1;
    std::vector<gr_complex> S0_table_qam64_b2;
    std::vector<gr_complex> S0_table_qam64_b3;
    std::vector<gr_complex> S0_table_qam64_b4;
    std::vector<gr_complex> S0_table_qam64_b5;
    std::vector<gr_complex> S1_table_qam64_b0;
    std::vector<gr_complex> S1_table_qam64_b1;
    std::vector<gr_complex> S1_table_qam64_b2;
    std::vector<gr_complex> S1_table_qam64_b3;
    std::vector<gr_complex> S1_table_qam64_b4;
    std::vector<gr_complex> S1_table_qam64_b5;
    create_complex_array(length, S0_table_qam64_b0_real,S0_table_qam64_b0_imag, S0_table_qam64_b0);
    create_complex_array(length, S0_table_qam64_b1_real,S0_table_qam64_b1_imag, S0_table_qam64_b1);
    create_complex_array(length, S0_table_qam64_b2_real,S0_table_qam64_b2_imag, S0_table_qam64_b2);
    create_complex_array(length, S0_table_qam64_b3_real,S0_table_qam64_b3_imag, S0_table_qam64_b3);
    create_complex_array(length, S0_table_qam64_b4_real,S0_table_qam64_b4_imag, S0_table_qam64_b4);
    create_complex_array(length, S0_table_qam64_b5_real,S0_table_qam64_b5_imag, S0_table_qam64_b5);
    create_complex_array(length, S1_table_qam64_b0_real,S1_table_qam64_b0_imag, S1_table_qam64_b0);
    create_complex_array(length, S1_table_qam64_b1_real,S1_table_qam64_b1_imag, S1_table_qam64_b1);
    create_complex_array(length, S1_table_qam64_b2_real,S1_table_qam64_b2_imag, S1_table_qam64_b2);
    create_complex_array(length, S1_table_qam64_b3_real,S1_table_qam64_b3_imag, S1_table_qam64_b3);
    create_complex_array(length, S1_table_qam64_b4_real,S1_table_qam64_b4_imag, S1_table_qam64_b4);
    create_complex_array(length, S1_table_qam64_b5_real,S1_table_qam64_b5_imag, S1_table_qam64_b5);
    S0_table_qam64.push_back(S0_table_qam64_b0);
    S0_table_qam64.push_back(S0_table_qam64_b1);
    S0_table_qam64.push_back(S0_table_qam64_b2);
    S0_table_qam64.push_back(S0_table_qam64_b3);
    S0_table_qam64.push_back(S0_table_qam64_b4);
    S0_table_qam64.push_back(S0_table_qam64_b5);
    S1_table_qam64.push_back(S1_table_qam64_b0);
    S1_table_qam64.push_back(S1_table_qam64_b1);
    S1_table_qam64.push_back(S1_table_qam64_b2);
    S1_table_qam64.push_back(S1_table_qam64_b3);
    S1_table_qam64.push_back(S1_table_qam64_b4);
    S1_table_qam64.push_back(S1_table_qam64_b5);

}
