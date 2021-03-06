/* -*- c++ -*- */
/*
 * Copyright 2004,2010 Free Software Foundation, Inc.
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
 * WARNING: This file is automatically generated by generate_gr_fir_filter_XXX.py
 * Any changes made to this file will be overwritten.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gr_fir_filter_fff.h>
#include <gr_fir_fff.h>
#include <gr_fir_util.h>
#include <gr_io_signature.h>

gr_fir_filter_fff_sptr gr_make_fir_filter_fff (int decimation, const std::vector<float> &taps)
{
  return gnuradio::get_initial_sptr (new gr_fir_filter_fff (decimation, taps));
}


gr_fir_filter_fff::gr_fir_filter_fff (int decimation, const std::vector<float> &taps)
  : gr_sync_decimator ("fir_filter_fff",
		       gr_make_io_signature (1, 1, sizeof (float)),
		       gr_make_io_signature (1, 1, sizeof (float)),
		       decimation),
    d_updated (false)
{
  d_fir = gr_fir_util::create_gr_fir_fff (taps);
  set_history (d_fir->ntaps ());
}

gr_fir_filter_fff::~gr_fir_filter_fff ()
{
  delete d_fir;
}

void
gr_fir_filter_fff::set_taps (const std::vector<float> &taps)
{
  d_new_taps = taps;
  d_updated = true;
}

std::vector<float>
gr_fir_filter_fff::taps () const
{
  return d_new_taps;
}

int
gr_fir_filter_fff::work (int noutput_items,
		   gr_vector_const_void_star &input_items,
		   gr_vector_void_star &output_items)
{
  float *in = (float *) input_items[0];
  float *out = (float *) output_items[0];

  if (d_updated) {
    d_fir->set_taps (d_new_taps);
    set_history (d_fir->ntaps ());
    d_updated = false;
    return 0;		     // history requirements may have changed.
  }

  if (decimation() == 1)
    d_fir->filterN (out, in, noutput_items);

  else
    d_fir->filterNdec (out, in, noutput_items, decimation());

  return noutput_items;
}
