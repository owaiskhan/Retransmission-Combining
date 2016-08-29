#ifndef INCLUDED_DIGITAL_OFDM_DEINTERLEAVER_H
#define INCLUDED_DIGITAL_OFDM_DEINTERLEAVER_H 

#include <itpp/itbase.h>

#include "digital_ofdm_global.h"

class digital_ofdm_deinterleaver{

    private:
        int Nsubs;
        int Ncol;
        std::vector< std::vector<int> > deinterleave_patterns;

    public:
        digital_ofdm_deinterleaver();
        digital_ofdm_deinterleaver(int nSubs);
        void calculate_deinterleave_patterns();
        void deinterleave_data(itpp::vec& data, itpp::vec& deinlv_data, int mcs);

};


#endif /* INCLUDED_DIGITAL_OFDM_DEINTERLEAVER_H */
