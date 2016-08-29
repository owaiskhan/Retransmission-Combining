/*  	
 * 
 */ 

#include "digital_fec_decode.h"

using namespace itpp;

#include <itpp/itcomm.h>

digital_fec_decode::digital_fec_decode(unsigned int code_rate){

    set_puncture_matrix(code_rate);
}



void 
digital_fec_decode::set_puncture_matrix(unsigned int code_rate){

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
			break;
	}
}

/*
*/

void 
digital_fec_decode::decode(unsigned int code_rate, vec& code_input, bvec& code_output){


    set_puncture_matrix(code_rate);

    Punctured_Convolutional_Code code; // IT++ punctured BCC class
    Convolutional_Code states; // IT++ BCC class
		
    generator.set_size(2,false); // create array
    generator(0) = 0133; // define BCC generators
    generator(1) = 0171;

    //code.set_truncation_length(100);
    code.set_generator_polynomials(generator, 7); // as name suggests
    code.set_puncture_matrix(puncture_matrix); // as name suggests
    states.set_start_state(0); // start at 0 state

    //code_input = -2.0*(code_input) + 1.0;

    code.decode_tail(code_input, code_output); // viterbi decoder
}

