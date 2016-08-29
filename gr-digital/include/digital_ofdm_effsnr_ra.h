#ifndef INCLUDED_DIGITAL_OFDM_EFFSNR_RA_H
#define INCLUDED_DIGITAL_OFDM_EFFSNR_RA_H

#include <itpp/itbase.h>
#include "digital_ofdm_rate_adaptation.h"

class digital_ofdm_effsnr_ra: public digital_ofdm_rate_adaptation{

    public:
        digital_ofdm_effsnr_ra();
        digital_ofdm_effsnr_ra(int nSubs);
        void calculate_uber(std::vector<double>& uber, std::vector<double>& esno, Modulation mod);
        double calculate_delivery_ratio_esno(std::vector<double>& esno, int mcs, int pktSize);
        int select_rate(std::vector<double>& esno, int pktSize);
        double calculate_transmission_time(int mcs, int payload);

};
#endif /* INCLUDED_DIGITAL_OFDM_EFFSNR_RA_H */
