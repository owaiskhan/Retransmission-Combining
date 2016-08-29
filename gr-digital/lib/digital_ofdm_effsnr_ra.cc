/*
 * Title: Rate adaptation
 * Created By: Owais Khan                                            
 * Creation Date: 3/2/2015
 *      
 * Description: Base Class for Rate adaptaiton                       
 *  
*/                                                                   
                                                                     
#include "digital_ofdm_effsnr_ra.h"                            

digital_ofdm_effsnr_ra::digital_ofdm_effsnr_ra(){ }

digital_ofdm_effsnr_ra::digital_ofdm_effsnr_ra(int nSubs):
            digital_ofdm_rate_adaptation(nSubs){

}


/*
*/
double
digital_ofdm_effsnr_ra::calculate_transmission_time(int mcs, int payload){

    int num_subs = Nsubs;
    double cbits_per_symbol = CBITSPERMCS[mcs];
    double numOfdmSymbols = ceil(payload/(cbits_per_symbol*num_subs));

    double tpayload = tOfdm*numOfdmSymbols;

    double base_bits_per_symbol = 0.5;
    double tphyheader = tOfdm*ceil(PHYsize/(base_bits_per_symbol*num_subs));
    double tmacheader = tOfdm*ceil(MACsize/(cbits_per_symbol*num_subs));

    double numAckOfdmSymbols = ceil(ACKsize/(cbits_per_symbol*num_subs));
    double tack = tpreamble + tphyheader + tOfdm*numAckOfdmSymbols;
#if WIFITIME == 1
    double txtime = tpreamble + tSIFS + tDIFS + tpayload + tphyheader + tmacheader + tack;
#else
    double txtime = tpreamble + tpayload + tphyheader + tmacheader + tack;
#endif

    return txtime;
}


/*
*/
void
digital_ofdm_effsnr_ra::calculate_uber(std::vector<double>& uber, std::vector<double>& esno, Modulation mod){


    double const1[] ={0.5, 1, 5.0, 21.0};
    double const2[] ={1, 1, 0.75, 0.5833};
    uber.resize(esno.size());

    for(int i=0; i<esno.size(); i++){
        uber[i] = const2[mod]*itpp::Qfunc(sqrt(esno[i]/const1[mod]));
    }
}


/*
*/
double
digital_ofdm_effsnr_ra::calculate_delivery_ratio_esno(std::vector<double>& esno, int mcs, int pktSize){

    Modulation mod = MCS2MODULATION[mcs];

    std::vector<double> uber;
    calculate_uber(uber, esno, mod);

    double avg_ber=0;
    for(int i=0; i<uber.size(); i++){
        avg_ber+=uber[i];
    }
    avg_ber = avg_ber/esno.size();

    double effesno = calculate_effesno(avg_ber, mod);

    double dr = calculate_delivery_ratio(effesno, mcs, pktSize);

    return dr;

}


/*
*/
int 
digital_ofdm_effsnr_ra::select_rate(std::vector<double>& esno, int pktSize){
  
    
    double* tput = new double[numRates];
    double max_tput=-9999.99;
    int max_mcs=-1;
    for (int crate_index=numRates-1; crate_index >=0; crate_index--){

        int curr_mcs = crate_index;

        // calculate delivery ratio
        double dratio = calculate_delivery_ratio_esno(esno,curr_mcs, pktSize);

        // calculate transmission time
        int payload = pktSize*8;
        double txtime = calculate_transmission_time(curr_mcs, payload);

        tput[crate_index] = dratio*payload/txtime;
        //printf("mcs=%d, dratio=%f, tput=%f\n",curr_mcs,dratio,tput[crate_index]);

        if(max_tput<tput[crate_index]){
            max_tput = tput[crate_index];
            max_mcs = crate_index;
        }
        
    }

    free(tput);
    if(max_mcs == -1){
        printf("Invalid mcs value.\n");
        exit(0);
    }
    return max_mcs;

}
