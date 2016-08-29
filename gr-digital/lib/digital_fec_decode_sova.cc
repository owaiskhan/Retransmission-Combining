/*  	
 * Title:  Bit level error control decoding for 802.11n
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

#include "digital_fec_decode.h"

using namespace itpp;

#include <itpp/itcomm.h>

digital_fec_decode::digital_fec_decode(unsigned int code_rate, unsigned int sova)	// apurv++ sova
{
	switch (code_rate) // defining appropriate puncture matrix (17.3.5.6)
	{
		case 0: // 1/2
			puncture_matrix = "1;1";
			d_code_rate = 0.5;
			d_puncture_ind.set_size(2);
			break;
		case 1: // 2/3
			puncture_matrix = "1 1;1 0";
			d_code_rate = 2.0/3.0;
			d_puncture_ind.set_size(4);
			break;
		case 2: // 3/4
			puncture_matrix = "1 1 0;1 0 1";
			d_code_rate = 3.0/4.0;
			d_puncture_ind.set_size(6);
			break;
		case 3: // 5/6
			puncture_matrix= "1 1 0 1 0;1 0 1 0 1";
			d_code_rate = 5.0/6.0;
			d_puncture_ind.set_size(10);
			break;
		default:
			puncture_matrix = "1;1";
			d_code_rate = 0.5;
			d_puncture_ind.set_size(2);
			break;
	}

	/* apurv++ start for sova */
	d_sova = sova;
	if(d_sova) {
	   printf("digital_fec_decode:: sova\n"); fflush(stdout);
	   generator.set_size(2,false); // create array
	   generator(0) = 0133; // define BCC generators
	   generator(1) = 0171;

	   int cl = 7;
	   std::string map_metric="maxlogMAP";
	   d_siso.set_generators(generator, cl);
	   d_siso.set_map_metric(map_metric); 
	
	   /* get the column vector of mat */
	   for(int i=0; i<puncture_matrix.cols(); i++) {
	       d_puncture_ind.set_subvector((i*2), puncture_matrix.get_col(i));
	   }
	}
	/* apurv++ end */
}

std::string digital_fec_decode::decode(unsigned int code_type, std::string code_input){

    const char* data = code_input.data(); 
    unsigned int data_size = code_input.size();

    unsigned int size_byte = 8;
    vec coded_data(data_size*size_byte);

    // convert data from string to bvec
    for (unsigned int i=0; i<data_size; i++){
        for (unsigned int j=0; j< size_byte; j++){
            coded_data[j+(i*size_byte)] = (data[i]&(0x1<<(size_byte-j-1)))>>(size_byte-j-1);
        }
    }

    bvec code_output;
    decode_internal(coded_data, code_output);
    //std::cout<<"coded data len: "<<coded_data.length()<<std::endl;
    //std::cout<<"uncoded data len: "<<code_output.length()<<std::endl;

    //std::cout<<code_output<<std::endl;

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
        if(code_output[k]==1){
            output[i] |= 0x1<<(7-j);
        } 
        j++;
    }

    std::string uncoded_string;// = output;
    uncoded_string.assign(output,num_output_bytes+extra);    

return uncoded_string;

}

void digital_fec_decode::decode_internal(vec& code_input, bvec& code_output){
    Punctured_Convolutional_Code code; // IT++ punctured BCC class
    Convolutional_Code states; // IT++ BCC class
		
    generator.set_size(2,false); // create array
    generator(0) = 0133; // define BCC generators
    generator(1) = 0171;

    code.set_generator_polynomials(generator, 7); // as name suggests
    code.set_puncture_matrix(puncture_matrix); // as name suggests
    states.set_start_state(0); // start at 0 state

    code_input = -2.0*(code_input) + 1.0;

    code.decode_tail(code_input, code_output); // viterbi decoder
}

/* apurv++ start */
std::string 
digital_fec_decode::decode_ppr(unsigned int code_type, unsigned int data_bytes, const std::vector<float>& in_llrs) {

    d_out_llrs.clear();
    int size_byte = 8;
    int input_size = in_llrs.size();
    //int total_punc_size = ceil((input_size/2)*(1.0/d_code_rate));
    int total_punc_size = ceil((input_size*d_code_rate)*2);
    int input_data_size = (data_bytes*8)+6;			// tail bits: 6

    printf("decode_ppr: input_size: %d, size_after_punc: %d\n", input_size, total_punc_size);
    printf("decode_ppr: input_data_size: %d, data_bytes: %d\n", input_data_size, data_bytes); fflush(stdout); 

    vec intrinsic_coded;
    intrinsic_coded.set_size(total_punc_size);
    d_apriori_data.set_size(input_data_size);

    d_apriori_data.zeros();
    intrinsic_coded.zeros();

    int punc_bits = d_puncture_ind.size();
    int ind = 0;
    for(int i = 0; i < total_punc_size; i++) {
#if 1
	if(d_puncture_ind[i%punc_bits] == 1)
	   intrinsic_coded[i] = in_llrs[ind++];
#else
	if(in_llrs[ind++] < 0) 
	   intrinsic_coded[i] = -1;
	else
	   intrinsic_coded[i] = 1;
#endif
    }
    assert(ind == input_size);

#ifdef NSC
    d_siso.nsc(d_extrinsic_coded, d_extrinsic_data, intrinsic_coded, d_apriori_data, false);
    int output_bits = d_extrinsic_data.size();
    printf("output_coded_bits(NSC): %d, output_data_bits: %d\n", d_extrinsic_coded.size(), d_extrinsic_data.size());
#else

    Punctured_Convolutional_Code code; // IT++ punctured BCC class
    code.set_generator_polynomials(generator, 7); // as name suggests
    code.set_puncture_matrix(puncture_matrix); // as name suggests

    bvec code_output;
    code.decode_tail(intrinsic_coded, code_output); // viterbi decoder
	/*
    printf("input size: %d, output_size: %d\n", intrinsic_coded.size(), code_output.length());
    std::cout<<"input:: "<<intrinsic_coded<<std::endl;
    std::cout<<"output:: "<<code_output<<std::endl;
    fflush(stdout);
	*/

    int output_bits = code_output.length();
    printf("output_bits: %d\n", output_bits);
#endif
  
    unsigned int num_output_bytes = (unsigned int) output_bits/size_byte;
    unsigned int rem_bits = output_bits % size_byte;
    unsigned int extra = (rem_bits == 0) ? 0 : 1;

    // convert coded data from bvec to string
    char* output = new char[num_output_bytes+extra];
    memset(output, 0, num_output_bytes+extra);
    int i = 0, j = 0;
    printf("out_llrs: "); fflush(stdout);
    for (int k = 0; k < output_bits; k++) {
        i = k/size_byte;
        j = j % 8;
#ifdef NSC
	d_out_llrs.push_back(d_extrinsic_data[k]);
	if(d_extrinsic_data[k] < 0.0)				// bit = 1 iff llr > 0.0
#else
	if(code_output[k]==1)
#endif
	{
            output[i] |= 0x1 << (7 - j);
        }
	if(k < 20) { printf("%.2f ", d_extrinsic_data[k]); fflush(stdout); }
        j++;
    }
    printf("\n"); fflush(stdout);
    std::string uncoded_string;
    uncoded_string.assign(output, num_output_bytes + extra);

    return uncoded_string;
}

std::vector<float>
digital_fec_decode::get_out_llrs() {
    return d_out_llrs; 
}
/* apurv++ end */
