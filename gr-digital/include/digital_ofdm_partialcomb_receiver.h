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

#ifndef INCLUDED_DIGITAL_OFDM_PARTIALCOMB_RECEIVER_H
#define INCLUDED_DIGITAL_OFDM_PARTIALCOMB_RECEIVER_H


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
#include "digital_ofdm_interleaver.h"
#include "digital_ofdm_deinterleaver.h"
#include "digital_ofdm_global.h"

#include "digital_ofdm_effsnr_ra.h"
#include "digital_ofdm_smart_ra.h"

#define TRACEFILE 1

unsigned int RAMODE=0; // Effsnr = 0, Smart = 1, Soft = 2

class digital_ofdm_partialcomb_receiver;
typedef boost::shared_ptr<digital_ofdm_partialcomb_receiver> digital_ofdm_partialcomb_receiver_sptr;

DIGITAL_API digital_ofdm_partialcomb_receiver_sptr 
digital_make_ofdm_partialcomb_receiver (const std::vector<gr_complex> &sym_position, 
			      const std::vector<unsigned char> &sym_value_out,
			      gr_msg_queue_sptr target_queue, unsigned int occupied_tones,
                  unsigned int ra_mode,unsigned int tx_gain, unsigned int bw_intf, unsigned int run,
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
class DIGITAL_API digital_ofdm_partialcomb_receiver : public gr_sync_block
{
  friend DIGITAL_API digital_ofdm_partialcomb_receiver_sptr 
  digital_make_ofdm_partialcomb_receiver (const std::vector<gr_complex> &sym_position, 
				const std::vector<unsigned char> &sym_value_out,
				gr_msg_queue_sptr target_queue, unsigned int occupied_tones,
                unsigned int ra_mode, unsigned int tx_gain, unsigned int bw_intf,
                unsigned int run,
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
  int Nsubs;
  int firsttx_mcs;
  int d_lastpktsucc_id;

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
  digital_ofdm_interleaver interleaver;
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
  void approxllrslicer(const gr_complex x, double noise, int bits_per_mod, std::vector<std::vector<gr_complex> >& S0_table, std::vector<std::vector<gr_complex> >& S1_table, std::vector<double>& output);
  void llr_demapper(gr_complex *in,
			itpp::vec& out_vec, unsigned int& offset, unsigned int bits_per_mod);

  void llr_demapper_alldata(std::vector<gr_complex>& symbols, std::vector<double>& ref_esno, itpp::vec& out_vec, unsigned int bits_per_mod);
  // combining variables
  //std::vector<int> retx_bits;
  itpp::vec data_comb;
  void combine_data(itpp::vec& recv_data, std::vector<int>& comb_idxs, itpp::vec& output);
  void assign_bits_to_subs(std::vector<std::vector<int> >& bit2submap, int bitmap_size, int mcs, bool inlvFlag);


  //timing variables
  double llrtime_elapsed;


 protected:
  digital_ofdm_partialcomb_receiver(const std::vector<gr_complex> &sym_position, 
			  const std::vector<unsigned char> &sym_value_out,
			  gr_msg_queue_sptr target_queue, unsigned int occupied_tones,
              unsigned int ra_mode, unsigned int tx_gain, unsigned int bw_intf,
              unsigned int run,
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
  void calculate_ltfpreamble_snr();

  int get_codeddata_length(int pktlen, int mcs);

 public:
  ~digital_ofdm_partialcomb_receiver();

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
    unsigned short verifyid: 8;
    unsigned int hdr_crc: 32;
  } __attribute__((packed)) PHY_HEADER;

  typedef struct feedback_type{
    unsigned short pktid;
    unsigned char mcs;
    unsigned char succ;
    unsigned char refid;
    unsigned char retx[80];
  }FeedbackStruct;

  PHY_HEADER header;
  unsigned int coded_hdr_bitlen;
  unsigned int coded_data_bitlen;
  unsigned int d_hdr_offset;
  unsigned int d_data_offset;
  unsigned char d_refid;
  itpp::vec hdr_demod;
  itpp::vec data_demod;

  // socket info
  void open_feedback_socket();
  void send_feedback(int mcs, int succ, int pktid, std::vector<int>& retxbitmap);
  int sock;
  struct sockaddr_in server;
  struct sockaddr_in from;

  // rate adaptation
  digital_ofdm_effsnr_ra ra_effsnr;
  digital_ofdm_smart_ra ra_smart;

  int total_pkts;
  int succ_pkts;
  double total_time;
  int succrcvd_data;

  std::vector<int> retx_subs; //bitmap of subs to retx
  std::vector<std::vector<int> > d_bit2submap;

  // ltf preamble variables
  std::complex<double>** rcvd_ltfpreambles;
  std::vector<std::complex<double> > eqd_ltfpreambles;
  std::vector<double> ltf_subcarrier_esno;
  std::vector<double> h_estimates; // used to store h from frame_acq
  std::vector<int> ltf_preamble_rx;
  int preamble_count;
  std::ofstream outfile_ltfsnr;
  std::ofstream outfile_refsnr;

  // evm variables
  std::vector<double> evm_tables;

  // SNR Calculation
  std::vector<double> prev_refsnr;
  std::vector<double> avg_ltfsnr;
  std::vector<double> Sk_ltfsnr;
  std::vector<double> std_ltfsnr;
  int ltfsnr_count;
  double ewma_alpha;
  double std_beta;

  // ref-snr
  std::vector<gr_complex> retx_refsymbs;
  void calculate_refsnr(std::vector<gr_complex>& data, std::vector<double>& snr, int mcs);
  void calculate_retx_refsnr(std::vector<gr_complex>& data, std::vector<gr_complex>& ref_data, std::vector<double>& snr, int mcs);
  void create_ref_symbols(int mcs, int retx_mcs,std::vector<int>& retxsubs, std::vector<gr_complex>& symbols);
  void calculate_distsnr(std::vector<gr_complex>& symbols, std::vector<double>& esno, std::vector<double>& snr,int mcs);
  void calculate_adjdistsnr(std::vector<gr_complex>& symbols, std::vector<double>& esno, std::vector<double>& snr,int mcs);
  std::vector<gr_complex> rcvd_symbols;
  void calculate_ltfpreamble_distsnr(std::vector<double>& EsNo, std::vector<double>& Snr);
  void calculate_kdistsnr(std::vector<gr_complex>& symbols, std::vector<double>& esno, std::vector<double>& snr, std::vector<double>& mean_snr, std::vector<double>& std_snr,int mcs);

 // SOFT variables
  void soft_combine_data(std::vector<gr_complex>& symbols, std::vector<double>& esno, int retx_count);
  std::vector<gr_complex> soft_rcvd_symbols; // contains combined symbols as sum
  std::vector<gr_complex> soft_combined_data; // contains final combined symbols
  std::vector<double> cum_esno_sqrt;
  std::vector<double> cum_esno;
  int soft_prev_mcs;

 // llr value snr
 void calculate_llrsnr(itpp::vec& llr_values, std::vector<double>& snr, int mcs);

 // time calculation
 double calculate_transmission_time(int mcs, int payload, int txcount, int retx_mod);

 // DEBUG VARIABLES
 std::vector<double> refEffsnrVec;
 std::vector<double> distEffsnrVec;
 std::vector<double> distadjEffsnrVec;
 std::vector<double> ltfEffsnrVec;
 std::vector<double> ltfdistEffsnrVec;
 double calculate_effsnr_from_subsnr(std::vector<double>& esno,digital_ofdm_effsnr_ra& ra, int mcs);

 std::ofstream outfile_effsnr;
 std::ofstream dist_subsnr;
 std::ofstream kdist_subsnr;
 std::ofstream ltf_subsnr;

#if TRACEFILE == 1
 // Trace Info
 std::ofstream traceOutFile;
 std::vector<int> traceMcs;
 std::vector<double> tracePayloadLen;
 std::vector<int> traceSuccPkt;
 std::vector<int> tracePktid;
 std::vector<std::vector<int> > traceRetxSubs;
 std::vector<std::vector<int> > traceErrsPerSub;
 std::vector<int> traceBitErrors;
 std::vector<int> traceTotalBits;
#endif

};
#endif /* INCLUDED_GR_OFDM_PARTIALCOMB_RECEIVER_H */
