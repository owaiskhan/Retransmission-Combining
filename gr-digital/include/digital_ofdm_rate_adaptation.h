#ifndef INCLUDED_DIGITAL_OFDM_RATE_ADAPTATION_H
#define INCLUDED_DIGITAL_OFDM_RATE_ADAPTATION_H

#include <itpp/itbase.h>

#include "digital_ofdm_global.h"

/*
*/
template <typename T>
void print_vec(std::vector<T> vec){

    for(int i=0; i<vec.size(); i++){
        std::cout<<vec[i]<<", ";
    }
    std::cout<<"\n";
}


class digital_ofdm_rate_adaptation{

    public:
        digital_ofdm_rate_adaptation();
        digital_ofdm_rate_adaptation(int nSubs);
        double dratio_helper(double esno, int mcs, int pktsize_idx);
        double calculate_delivery_ratio(double esno, int mcs, int pktSize);
        void read_effsnr2dr_file(double data[][201],const char* filename);
        int bits_per_modulation(Modulation mod);

        double calculate_effesno(double uber, Modulation mod);
        double qfuncInv(double Qx);

    protected:
        unsigned int Nsubs, selected_mcs;
        double tOfdm, tpreamble, MACsize, PHYsize, ACKsize;
#if WIFITIME == 1
        double tSIFS, tDIFS, tSlottime;
#endif

        int numRates;
        double drRate[8][10][201];

};

#endif /* INCLUDED_DIGITAL_OFDM_RATEADAPTATION_H */
