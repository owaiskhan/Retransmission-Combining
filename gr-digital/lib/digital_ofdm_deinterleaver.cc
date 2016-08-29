/*
 * Title: Data De-Interleaver
 * Created By: Owais Khan
 * Creation Date: 3/3/2015
 * 
 * Description: Data De-Interleaver for Nsubs = 80
 *
*/ 

#include "digital_ofdm_deinterleaver.h"

digital_ofdm_deinterleaver::digital_ofdm_deinterleaver() { }

digital_ofdm_deinterleaver::digital_ofdm_deinterleaver(int nSubs):Nsubs(nSubs){

    Ncol = 16;
    calculate_deinterleave_patterns();
}


void digital_ofdm_deinterleaver::calculate_deinterleave_patterns(){

    int numMCS = 8;
    for(int mcs=0; mcs<numMCS; mcs++){

        Modulation mod = MCS2MODULATION[mcs];
        int bits_per_mod = BITSPERMOD[mod];

        int Ncbps = Nsubs*bits_per_mod;
        int Nrow = Ncbps/16;
        int s = std::max(bits_per_mod/2,1);

        int perm1, k;
        std::vector<int> deinlv_pattern;
        for(int perm2 = 0; perm2 < Ncbps; perm2++){
            perm1 = s*floor(perm2/s) + ( (perm2 + (int)floor(Ncol*perm2/Ncbps) )%s);
            k = (Ncol*perm1) - (Ncbps-1)*floor(perm1/Nrow);

            deinlv_pattern.push_back(k);
        }
        deinterleave_patterns.push_back(deinlv_pattern);
    }

}

void
digital_ofdm_deinterleaver::deinterleave_data(itpp::vec& data, itpp::vec& deinlv_data, int mcs){

    
    Modulation mod = MCS2MODULATION[mcs];
    int bits_per_mod = BITSPERMOD[mod];
    int Ncbps = Nsubs*bits_per_mod;

    int bits_per_ofdm = Nsubs*bits_per_mod;

    if((data.size()%bits_per_ofdm)!=0){
        printf("Data should be ofdm aligned.");
        exit(0);
    }

    deinlv_data.set_size(data.size()); 

    std::vector<int> deinlv_pattern;

    deinlv_pattern = deinterleave_patterns[mcs];

    int ofdm = 0, index = 0;
    for(int i=0; i<data.size(); i++){
        ofdm = floor(i/Ncbps);
        index = i%Ncbps;
        deinlv_data[ ofdm*(deinlv_pattern.size()) + deinlv_pattern[index] ] = data[i];
    }

}
