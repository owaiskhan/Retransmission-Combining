/*   	
 * Title:  Bit level error control coding for 802.11n
 * Created By: Robert Daniels
 * Creation Date: 10/08/2006
 * 
 * Description: Header for "channel_encoding_11n.h"
 * 
 * Revision: v0.00 - 10/08/2006 - Initial Release
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


#ifndef Channel_Encoding_11N
#define Channel_Encoding_11N

#include <itpp/itbase.h>

class digital_fec_encode
{
	public:
		digital_fec_encode(unsigned int code_rate); // for init

        void set_puncture_matrix(unsigned int code_rate);

        std::string encode(unsigned int code_rate, const std::string code_input);
        void encode_internal(itpp::bvec& code_input, itpp::bvec& code_output);
		// forward error correction encoding (21.3.10.5 and 21.3.10.6)
	private:
		itpp::bmat puncture_matrix; // for setting up BCC puncturing
		itpp::ivec generator; // for BCC generator polynomials
};
#endif
