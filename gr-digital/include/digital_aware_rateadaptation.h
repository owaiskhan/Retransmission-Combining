/*  	
 * Title: Aware rate adaptation
 * Created By: Owais Khan
 * Creation Date: 01/25/2053
 * 
 * Description: Rate adaptaiton which consider re-transmission till depth-k
 *
*/ 


#ifndef Digital_Aware_RA_11N_H
#define Digital_Aware_RA_11N_H

#include <itpp/itbase.h>

struct Node{
    
    itpp::vec comb_esno_bit;
    itpp::ivec bit2submap;
    itpp::vec txrates;
    double delivery_ratio;
    double tx_time;
    double tput;
};

class digital_aware_ra
{
	public:
		digital_aware_ra(int nsubs, int depth, int listlen); // for init

        int select_rate(const std::vector<float> &esno_per_bit,const std::vector<float> &subcarrier_esno);



    private:
        std::vector<Node> ActiveList;
        unsigned int maxList, maxDepth, Nsubs, selected_mcs;
        double tOfdm, tpreamble, tSIFS, tDIFS, tSlottime, MACsize, PHYsize, ACKsize;

};
#endif
