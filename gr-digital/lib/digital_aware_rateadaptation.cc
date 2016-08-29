/*
 * Title: Aware rate adaptation
 * Created By: Owais Khan
 * Creation Date: 01/25/2015
 * 
 * Description: Rate adaptaiton which consider re-transmission till depth-k
 *
*/ 


#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <map>
#include <itpp/base/sort.h>
#include <assert.h>
#include <sys/time.h>
#include "digital_aware_rateadaptation.h"


using namespace itpp;


digital_aware_ra::digital_aware_ra(int nsubs, int depth, int listlen):selected_mcs(0),
curr_pktid(-1),tpreamble(160.0e-6),tSIFS(16.0e-6),tSlottime(9.0e-6),
tOfdm(80.0e-6),Nsubs(nsubs),maxDepth(depth),maxList(listlen),
tx_count(0){

    retxmaxDepth = maxDepth;
    MACsize = 28*8;
    ACKsize = 14*8;
    PHYsize = 10*8;
    tDIFS = tSIFS + 2*tSlottime;

} 


int
digital_aware_ra::select_rate(const std::vector<float>& esno_per_bit, const std::vector<float>& subcarrier_esno){


}
