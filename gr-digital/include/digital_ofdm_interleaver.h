#ifndef INCLUDED_DIGITAL_OFDM_INTERLEAVER_H
#define INCLUDED_DIGITAL_OFDM_INTERLEAVER_H 

#include <itpp/itbase.h>

#include "digital_ofdm_global.h"

class digital_ofdm_interleaver{

    private:
        int Nsubs;
        int Ncol;
        std::vector< std::vector<int> > interleave_patterns;

    public:
        digital_ofdm_interleaver();
        digital_ofdm_interleaver(int nSubs);
        void calculate_interleave_patterns();
        void interleave_data(std::vector<int>& data, std::vector<int>& inlv_data, int mcs);
        void interleave_data(itpp::bvec& data, itpp::bvec& inlv_data, int mcs);
        int get_subcarrier_index(int index, int mcs);

};


#endif /* INCLUDED_DIGITAL_OFDM_INTERLEAVER_H */
