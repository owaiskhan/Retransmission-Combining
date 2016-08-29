/*  	
 * Title:  Bit level error control coding for 802.11n
 * Created By: Robert Daniels
 * Creation Date: 10/08/2006
 * 
 * Description: Allow coding at bit level (convolutional or LDPC)
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

#include <itpp/itcomm.h>

#include "digital_fec_encode.h"

using namespace itpp;

digital_fec_encode::digital_fec_encode(unsigned int code_rate){

    set_puncture_matrix(code_rate);
}

void digital_fec_encode::set_puncture_matrix(unsigned int code_rate){

	switch (code_rate) // defining appropriate puncture matrix (17.3.5.6)
	{
		case 0: // 1/2
			puncture_matrix = "1;1";
			break;
		case 1: // 2/3
			puncture_matrix = "1 1;1 0";
			break;
		case 2: // 3/4
			puncture_matrix = "1 1 0;1 0 1";
			break;
		case 3: // 5/6
			puncture_matrix= "1 1 0 1 0;1 0 1 0 1";
			break;
		default:
			puncture_matrix = "1;1";
	}

}

std::string digital_fec_encode::encode(unsigned int code_rate, const std::string code_input){


    set_puncture_matrix(code_rate);

    const char* data = code_input.data(); 
    unsigned int data_size = code_input.size();

    unsigned int size_byte = 8;
    bvec uncoded_data(data_size*size_byte);

    // convert data from string to bvec
    for (unsigned int i=0; i<data_size; i++){
        for (unsigned int j=0; j< size_byte; j++){
            uncoded_data[j+(i*size_byte)] = (data[i]&(0x1<<(size_byte-j-1)))>>(size_byte-j-1);
        }
    }

    bvec code_output;
    encode_internal(uncoded_data, code_output);
    //std::cout<<"uncoded data len: "<<uncoded_data.length()<<std::endl;
    //std::cout<<"uncoded data: "<<uncoded_data<<std::endl;
    //std::cout<<"coded data len: "<<code_output.length()<<std::endl;
    //std::cout<<"coded data: "<<code_output<<std::endl;

    unsigned int num_output_bytes = (unsigned int) code_output.length()/size_byte;
    unsigned int rem_bits = code_output.length()%size_byte;
    unsigned int extra = (rem_bits==0) ? 0 : 1;

    // convert coded data from bvec to string
    char* output = new char[num_output_bytes+extra];
    memset(output,0,num_output_bytes+extra);

    int i=0, j=0;
    for (int k=0; k < code_output.length(); k++){
        i = k/size_byte;
        j = j%8;
        //std::cout<<"before: "<<int(output[i])<<std::endl;
        if(code_output[k]==1){
            output[i] |= 0x1<<(7-j);
        } 
        //std::cout<<"after: "<<int(output[i])<<std::endl;
        j++;
    }
/*
    for(int i=0; i<num_output_bytes+extra;i++){
        printf("%x",output[i]);
    }
    printf("\n");
*/
    //std::string coded_string = output;
    std::string coded_string;
    coded_string.assign(output,num_output_bytes+extra);
 
return coded_string;
}

void digital_fec_encode::encode_internal(bvec& code_input, bvec& code_output){

	Punctured_Convolutional_Code code; // IT++ punctured BCC class
	Convolutional_Code states; // IT++ BCC class
	
	generator.set_size(2,false); // create array
    generator(0) = 0133; // define BCC generators
  	generator(1) = 0171;

	code.set_generator_polynomials(generator, 7); // as name suggests
	code.set_puncture_matrix(puncture_matrix); // as name suggests
	states.set_start_state(0); // start at 0 state
	
    //std::cout<<"code_input: "<<code_input<<std::endl;
	code.encode_tail(code_input, code_output); // encode the stream
    #if 0
    //std::cout<<"code_output: "<<code_output<<std::endl;
    vec rx_data = to_vec(code_output);
    rx_data = -2.0*rx_data + 1.0;
    std::cout<<"rx_data: "<<rx_data<<std::endl;
    #endif
    #if 0
    bvec decoded_output;
    code.decode_tail(rx_data, decoded_output); // viterbi decoder
    //std::cout<<"decoded_output: "<<decoded_output<<std::endl;
    std::cout<<"bit diff: "<<(decoded_output - code_input)<<std::endl;
    #endif

}
