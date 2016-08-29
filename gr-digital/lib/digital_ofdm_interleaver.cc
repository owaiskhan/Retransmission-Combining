/*
 * Title: Data Interleaver
 * Created By: Owais Khan
 * Creation Date: 3/3/2015
 * 
 * Description: Data Interleaver for Nsubs = 80
 *
*/ 

#include "digital_ofdm_interleaver.h"

digital_ofdm_interleaver::digital_ofdm_interleaver(){ }

digital_ofdm_interleaver::digital_ofdm_interleaver(int nSubs):Nsubs(nSubs){

    Ncol = 16;
    calculate_interleave_patterns();
}

int 
digital_ofdm_interleaver::get_subcarrier_index(int index, int mcs){

    std::vector<int> inlv_pattern;
    inlv_pattern = interleave_patterns[mcs];


    Modulation mod = MCS2MODULATION[mcs];
    int bits_per_mod = BITSPERMOD[mod];
    int Ncbps = Nsubs*bits_per_mod;

    int ofdm_loc = index%Ncbps;
    int inlv_loc = inlv_pattern[ofdm_loc];
    int sub_num = floor(inlv_loc/bits_per_mod);
    return sub_num;

}

void digital_ofdm_interleaver::calculate_interleave_patterns(){

    int numMCS = 8;
    for(int mcs=0; mcs<numMCS; mcs++){

        Modulation mod = MCS2MODULATION[mcs];
        int bits_per_mod = BITSPERMOD[mod];

        int Ncbps = Nsubs*bits_per_mod;
        int Nrow = Ncbps/16;
        int s = std::max(bits_per_mod/2,1);

        int perm1, perm2;
        std::vector<int> inlv_pattern;
        for(int k = 0; k < Ncbps; k++){
            perm1 = Nrow*(k%Ncol) + floor(k/Ncol);
            perm2 = s*floor(perm1/s) + ( (perm1 + Ncbps - (int)floor(Ncol*perm1/Ncbps) )%s);

            inlv_pattern.push_back(perm2);
        }
        interleave_patterns.push_back(inlv_pattern);
    }

}

void
digital_ofdm_interleaver::interleave_data(std::vector<int>& data, std::vector<int>& inlv_data, int mcs){

    Modulation mod = MCS2MODULATION[mcs];
    int bits_per_mod = BITSPERMOD[mod];
    int Ncbps = Nsubs*bits_per_mod;
    int Nrow = Ncbps/16;

    int bits_per_ofdm = Nsubs*bits_per_mod;

    if((data.size()%bits_per_ofdm)!=0){
        printf("Data should be ofdm aligned.");
        exit(0);
    }

    inlv_data = std::vector<int> (data.size()); 

    std::vector<int> inlv_pattern;

    inlv_pattern = interleave_patterns[mcs];

    //printf("mcs=%d, inlv-pattern size:%d\n",mcs,inlv_pattern.size());

    int ofdm = 0, index = 0;
    for(int i=0; i<data.size(); i++){
        ofdm = floor(i/Ncbps);
        index = i%Ncbps;
        //int tmp = ofdm*inlv_pattern.size() + inlv_pattern[index];
        //printf("ofdm=%d, index=%d, inlv_index:%d\n",ofdm,index,tmp);
        inlv_data[ ofdm*(inlv_pattern.size()) + inlv_pattern[index] ] = data[i];
    }

}
void
digital_ofdm_interleaver::interleave_data(itpp::bvec& data, itpp::bvec& inlv_data, int mcs){

    
    Modulation mod = MCS2MODULATION[mcs];
    int bits_per_mod = BITSPERMOD[mod];
    int Ncbps = Nsubs*bits_per_mod;
    int Nrow = Ncbps/16;

    int bits_per_ofdm = Nsubs*bits_per_mod;

    if((data.size()%bits_per_ofdm)!=0){
        printf("Data should be ofdm aligned.");
        exit(0);
    }

    inlv_data.set_size(data.size()); 

    std::vector<int> inlv_pattern;

    inlv_pattern = interleave_patterns[mcs];

    //printf("mcs=%d, inlv-pattern size:%d\n",mcs,inlv_pattern.size());

    int ofdm = 0, index = 0;
    for(int i=0; i<data.size(); i++){
        ofdm = floor(i/Ncbps);
        index = i%Ncbps;
        //int tmp = ofdm*inlv_pattern.size() + inlv_pattern[index];
        //printf("ofdm=%d, index=%d, inlv_index:%d\n",ofdm,index,tmp);
        inlv_data[ ofdm*(inlv_pattern.size()) + inlv_pattern[index] ] = data[i];
    }

}
