/* -*- c++ -*- */
/*
 * Copyright 2007,2011 Free Software Foundation, Inc.
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

#ifndef INCLUDED_DIGITAL_OFDM_FRAME_SINK_H
#define INCLUDED_DIGITAL_OFDM_FRAME_SINK_H


#include <digital_api.h>
#include <gr_sync_block.h>
#include <gr_msg_queue.h>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>

#include <itpp/itbase.h>
#include <digital_crc32.h> 
#include "digital_fec_decode.h"
#include "digital_ofdm_deinterleaver.h"
#include "digital_ofdm_effsnr_ra.h"

class digital_ofdm_frame_sink;
typedef boost::shared_ptr<digital_ofdm_frame_sink> digital_ofdm_frame_sink_sptr;

DIGITAL_API digital_ofdm_frame_sink_sptr 
digital_make_ofdm_frame_sink (const std::vector<gr_complex> &sym_position, 
			      const std::vector<unsigned char> &sym_value_out,
			      gr_msg_queue_sptr target_queue, unsigned int occupied_tones,
			      float phase_gain=0.25, float freq_gain=0.25*0.25/4.0);

/*!
 * \brief Takes an OFDM symbol in, demaps it into bits of 0's and 1's, packs
 * them into packets, and sends to to a message queue sink.
 * \ingroup sink_blk
 * \ingroup ofdm_blk
 *
 * NOTE: The mod input parameter simply chooses a pre-defined demapper/slicer. Eventually,
 * we want to be able to pass in a reference to an object to do the demapping and slicing
 * for a given modulation type.
 */
class DIGITAL_API digital_ofdm_frame_sink : public gr_sync_block
{
  friend DIGITAL_API digital_ofdm_frame_sink_sptr 
  digital_make_ofdm_frame_sink (const std::vector<gr_complex> &sym_position, 
				const std::vector<unsigned char> &sym_value_out,
				gr_msg_queue_sptr target_queue, unsigned int occupied_tones,
				float phase_gain, float freq_gain);

 private:
  enum state_t {STATE_SYNC_SEARCH, STATE_HAVE_SYNC, STATE_HAVE_LTFPREAMBLE, STATE_HAVE_HEADER};


  gr_msg_queue_sptr  d_target_queue;		// where to send the packet when received
  state_t            d_state;
  unsigned int       d_occupied_carriers;

  gr_complex * d_derotated_output;  // Pointer to output stream to send deroated symbols out

  std::vector<gr_complex>    d_dfe;

  float d_phase;
  float d_freq;
  float d_phase_gain;
  float d_freq_gain;
  float d_eq_gain;

  std::vector<int> d_data_carriers;
  std::vector<int> d_pilot_carriers;

  // constellation variables
  std::vector<gr_complex> bpsk_const;
  std::vector<gr_complex> qpsk_const;
  std::vector<gr_complex> qam16_const;
  std::vector<gr_complex> qam64_const;

  std::vector<std::vector<int> > bpsk_values;
  std::vector<std::vector<int> > qpsk_values;
  std::vector<std::vector<int> > qam16_values;
  std::vector<std::vector<int> > qam64_values;

  void create_constellations();
  void get_constellations(int bits_per_mod, std::vector<gr_complex>& constel, std::vector<std::vector<int> >& constel_values);
 
  // interleaver
  digital_ofdm_deinterleaver deinterleaver;

  // llr variables
  std::vector<std::vector<gr_complex> > S0_table_bpsk;
  std::vector<std::vector<gr_complex> > S1_table_bpsk;
  std::vector<std::vector<gr_complex> > S0_table_qpsk;
  std::vector<std::vector<gr_complex> > S1_table_qpsk;
  std::vector<std::vector<gr_complex> > S0_table_qam16;
  std::vector<std::vector<gr_complex> > S1_table_qam16;
  std::vector<std::vector<gr_complex> > S0_table_qam64;
  std::vector<std::vector<gr_complex> > S1_table_qam64;

  std::vector<std::vector<gr_complex> > currS0_table;
  std::vector<std::vector<gr_complex> > currS1_table;

  void create_llr_tables();
  void set_currllr_table(int bits_per_mod);
  void llrslicer(const gr_complex x, double noise, int bits_per_mod, std::vector<std::vector<gr_complex> >& S0_table, std::vector<std::vector<gr_complex> >& S1_table, std::vector<double>& output);
  void llr_demapper(gr_complex *in,
			itpp::vec& out_vec, unsigned int& offset, unsigned int bits_per_mod);

  // timing variables
  double demodtime_elapsed;

 protected:
  digital_ofdm_frame_sink(const std::vector<gr_complex> &sym_position, 
			  const std::vector<unsigned char> &sym_value_out,
			  gr_msg_queue_sptr target_queue, unsigned int occupied_tones,
			  float phase_gain, float freq_gain);

  void enter_search();
  void enter_have_sync();

  bool enter_have_header();
  
 void bitslicer(const gr_complex x,std::vector<gr_complex>& constellation, std::vector<std::vector<int> >& constel_value, std::vector<int>& output);
  
  inline void equalize_interpolate_dfe(gr_complex *in, gr_complex *factor);

  void bit_demapper(gr_complex *in,
			itpp::vec& out_vec, unsigned int& offset, unsigned int bits_per_mod);

  void assign_subcarriers();
  void read_refsymbols_file();
  void calculate_refsnr(std::vector<double>& subcarrier_snr_ref);
  void calculate_ltfpreamble_snr();


 public:
  ~digital_ofdm_frame_sink();

  int work(int noutput_items,
	   gr_vector_const_void_star &input_items,
	   gr_vector_void_star &output_items);


  digital_fec_decode fec_decoder;
  typedef struct header_type{ 
    unsigned short pktlength: 12; 
    unsigned short mcs: 4;
    unsigned short pktid: 16;
    unsigned short whitener: 4;
    unsigned short retx_count: 4;
    unsigned short payloadlen: 16;
    unsigned int hdr_crc: 32;
  } __attribute__((packed)) PHY_HEADER;

  PHY_HEADER header;
  unsigned int coded_hdr_bitlen;
  unsigned int coded_data_bitlen;
  unsigned int d_hdr_offset;
  unsigned int d_data_offset;
  unsigned int llr_data_offset;
  itpp::vec hdr_demod;
  itpp::vec data_demod;
  itpp::vec llr_demod;

  // socket info
  void open_feedback_socket();
  void send_mcs_feedback(int mcs);
  int sock;
  struct sockaddr_in server;
  struct sockaddr_in from;

  // rate adaptation
  digital_ofdm_effsnr_ra ra_effsnr;
  int total_pkts;
  int succ_pkts;


  std::complex<double>** rcvd_ltfpreambles;
  std::vector<double> ltf_subcarrier_esno;
  std::vector<double> h_estimates; // used to store h from frame_acq
  int* ltf_preamble_rx;
  int preamble_count;
  std::ofstream outfile_ltfsnr;

};

#endif /* INCLUDED_GR_OFDM_FRAME_SINK_H */
