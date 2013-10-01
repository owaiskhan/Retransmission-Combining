/* -*- c++ -*- */
/*
 * Copyright 2005 Free Software Foundation, Inc.
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

/*
 * WARNING: This file is automatically generated by
 * generate_gr_rational_resampler_base_XXX.py Any changes made to this
 * file will be overwritten.
 */

#ifndef INCLUDED_GR_RATIONAL_RESAMPLER_BASE_SCC_H
#define	INCLUDED_GR_RATIONAL_RESAMPLER_BASE_SCC_H

#include <gr_core_api.h>
#include <gr_block.h>

class gr_rational_resampler_base_scc;
typedef boost::shared_ptr<gr_rational_resampler_base_scc> gr_rational_resampler_base_scc_sptr;
GR_CORE_API gr_rational_resampler_base_scc_sptr
gr_make_rational_resampler_base_scc (unsigned interpolation,
		     unsigned decimation,
		     const std::vector<gr_complex> &taps);

class gr_fir_scc;

/*!
 * \brief Rational Resampling Polyphase FIR filter with short input, gr_complex output and gr_complex taps
 * \ingroup filter_blk
 */
class GR_CORE_API gr_rational_resampler_base_scc : public gr_block
{
 private:
  unsigned 			d_history;
  unsigned 			d_interpolation, d_decimation;
  unsigned 			d_ctr;
  std::vector<gr_complex>	d_new_taps;
  bool				d_updated;
  std::vector<gr_fir_scc *> d_firs;

  friend GR_CORE_API gr_rational_resampler_base_scc_sptr 
  gr_make_rational_resampler_base_scc (unsigned interpolation, unsigned decimation, const std::vector<gr_complex> &taps);


  /*!
   * Construct a FIR filter with the given taps
   */
  gr_rational_resampler_base_scc (unsigned interpolation, unsigned decimation,
	  const std::vector<gr_complex> &taps);

  void install_taps (const std::vector<gr_complex> &taps);

 public:
  ~gr_rational_resampler_base_scc ();
  unsigned history () const { return d_history; }
  void  set_history (unsigned history) { d_history = history; }

  unsigned interpolation() const { return d_interpolation; }
  unsigned decimation() const { return d_decimation; }

  void set_taps (const std::vector<gr_complex> &taps);

  void forecast (int noutput_items, gr_vector_int &ninput_items_required);
  int  general_work (int noutput_items,
		     gr_vector_int &ninput_items,
		     gr_vector_const_void_star &input_items,
		     gr_vector_void_star &output_items);
};
 

#endif
