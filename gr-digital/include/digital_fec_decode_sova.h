/*  	
 * Title:  Error Control Decoding at receiver for 802.11n
 * Created By: Robert Daniels
 * Creation Date: 1/02/2007
 * 
 * Description: Header for "digital_fec_decode.cc"
 *
 * Revision: v0.00 - 1/02/2007 - Initial Release
 * 
 * Copyright (C) 2009  The University of Texas at Austin.
 * 
 * This file is part of Hydra: A wireless multihop testbed.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */ 


#ifndef Digital_Fec_Decode_11N_H
#define Digital_Fec_Decode_11N_H

#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <vector>
#include <assert.h>

#define NSC 1

class digital_fec_decode
{
	public:
	digital_fec_decode(unsigned int code_rate, unsigned int sova); // for init, apurv++ added sova
        std::string decode(unsigned int code_type, std::string code_input);
	std::string decode_ppr(unsigned int code_type, unsigned int data_bytes, 
			       const std::vector<float>& in_llrs); 	// apurv++ 
        void decode_internal(itpp::vec& code_input, itpp::bvec& code_output);
	std::vector<float> get_out_llrs();

	private:
		itpp::bmat puncture_matrix;
		itpp::ivec generator;

	/* apurv++ start */
	unsigned int d_sova;
	float d_code_rate;
	itpp::SISO d_siso;
	itpp::vec d_apriori_data, d_extrinsic_coded, d_extrinsic_data;
	itpp::bvec d_puncture_ind;
	std::vector<float> d_out_llrs;
	/* apurv++ end */
};
#endif
