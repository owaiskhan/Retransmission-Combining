/* -*- c++ -*- */
/*
 * Copyright 2012 Free Software Foundation, Inc.
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

#include <gr_multiply_cc.h>
#include <gr_io_signature.h>
#include <volk/volk.h>

gr_multiply_cc_sptr
gr_make_multiply_cc (size_t vlen)
{
  return gnuradio::get_initial_sptr(new gr_multiply_cc (vlen));
}

gr_multiply_cc::gr_multiply_cc (size_t vlen)
  : gr_sync_block ("gr_multiply_cc",
		   gr_make_io_signature (1, -1, sizeof (gr_complex)*vlen),
		   gr_make_io_signature (1, 1, sizeof (gr_complex)*vlen)),
    d_vlen(vlen)
{
 const int alignment_multiple =
   volk_get_alignment() / sizeof(gr_complex);
 set_alignment(alignment_multiple);
}

int
gr_multiply_cc::work (int noutput_items,
		      gr_vector_const_void_star &input_items,
		      gr_vector_void_star &output_items)
{
  gr_complex *out = (gr_complex *) output_items[0];
  int noi = d_vlen*noutput_items;

  memcpy(out, input_items[0], noi*sizeof(gr_complex));
  if(is_unaligned()) {
    for(size_t i = 1; i < input_items.size(); i++)
      volk_32fc_x2_multiply_32fc_u(out, out, (gr_complex*)input_items[i], noi);
  }
  else {
    for(size_t i = 1; i < input_items.size(); i++)
      volk_32fc_x2_multiply_32fc_a(out, out, (gr_complex*)input_items[i], noi);
  }
  return noutput_items;
}



