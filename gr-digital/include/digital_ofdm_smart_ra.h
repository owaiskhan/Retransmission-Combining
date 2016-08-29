#ifndef INCLUDED_DIGITAL_OFDM_SMART_RA_H
#define INCLUDED_DIGITAL_OFDM_SMART_RA_H

#include <itpp/itbase.h>
#include "digital_ofdm_rate_adaptation.h"
#include "digital_ofdm_interleaver.h"

class digital_ofdm_smart_ra: public digital_ofdm_rate_adaptation{

    private:
        // Struct definitions
        struct TxInfo{
            int pktSize;
            int txCount;
            // should be setup in constructor
            int nSubs;
            int maxlistLen;
            int searchDepth;
            int maxRetx;
        };

        struct Node{
            std::vector<double> comb_esno_bit;
            std::vector<std::vector<int> > bit2submap;
            std::vector<int> txmcs;
            std::vector< std::vector<int> > tx_bits;
            std::vector< std::vector<int> > tx_subs;
            std::vector<double> comb_uber;
            double deliveryRatio;
            double txTime;
            double tput;
        };


    digital_ofdm_interleaver interleaver;

    // variable declarations
    TxInfo d_txinfo;
    Node d_activeNode;
    int d_pktId;
    std::vector<std::vector<double> > d_EsNo;
    std::vector<Node> nodeList;
    // function declarations
    public:
        // Constructors
        digital_ofdm_smart_ra();
        digital_ofdm_smart_ra(int nSubs,int maxdepth,int maxlist,int maxretx);

        // Calculation functions
        void calculate_uber(std::vector<double>& uber, std::vector<double>& esno, std::vector<int>& tx_bits, Modulation mod);

        double calculate_delivery_ratio_esno(Node& node, std::vector<double>& esno_bitmap, std::vector<int>& tx_bits, int mcs, TxInfo& info);

        double calculate_transmission_time(int mcs, int payload, TxInfo& info);
        double calculate_retransmission_time(int payload, int mod2, int mcs, TxInfo& info);


        // First Transmission
        void first_transmission(Node& first_tx_node, std::vector<double>& EsNo, int first_tx_mcs, TxInfo& info);

        void find_rate_firsttransmission(Node& bestNode, std::vector<std::vector<double> >& EsNo, TxInfo& info);


        // Re-Transmission

        void retransmission(Node& node, double base_time, int depth, int retx_mod, int retx_bits, std::vector<int>& retransmitted_bits,std::vector<int>& bit_indices, std::vector<double>& EsNo, TxInfo& info);

        void linear_search(std::vector<Node>& tmpList, Node curr_node, int depth, int retx_mod, std::vector<int>& sub_indices, std::vector<int>& bit_indices, std::vector<double>& EsNo, TxInfo& info);


        void search_retransmission_combination(std::vector<Node>& tmpList, Node& curr_node, int depth, std::vector<double>& EsNo,TxInfo& info);


        void find_rate_retransmission(Node& bestNode, std::vector<std::vector<double> >& EsNo, TxInfo& info);


        void select_rate_internal(Node& active_node, std::vector<std::vector<double> >& subcarrier_esno, TxInfo& info);

        int select_rate(std::vector<double>& subsEsNo, itpp::vec& combllr, std::vector<int>& retx_subs, 
                        int txcount, int mcs, int pktSize, int pktId,int succ);


    // These were helper functions in standalone implementation
        void pick_best_rate(Node& best_node, std::vector<Node>& nodeList);
        void update_node_list(std::vector<Node>& list, Node& node, TxInfo& info);
        int get_rate_bitmap_len(int mcs, int pktlen, TxInfo& info);
        
        void assign_bits_to_subs(std::vector<std::vector<int> >& bit2submap, int bitmap_size, int mcs, TxInfo& info, bool inlvFlag);

        void assign_esno_per_bit(std::vector<double>& esno_bitmap, std::vector<double>& EsNo, int tx_bits_size, Modulation mod, TxInfo& info);

        void assign_esno_per_bit_retx(std::vector<double>& esno_bitmap, std::vector<double>& EsNo, int tx_bits_size, Modulation mod1, Modulation mod2, TxInfo& info);

        void remap_esno_values(std::vector<double>& remapped_esno, std::vector<double>& EsNo, Modulation mod1, Modulation mod2, TxInfo& info);


        // new helper functions
        void clear_activenode(Node& node);
        void setup_newtransmission(std::vector<double>& subsEsNo, TxInfo& info, int pktId);
        void update_combesno(Node& node, itpp::vec& combllr);
        void reboot_transmission(Node& node, TxInfo& info, int mcs);





/*                                                                          
*/  
struct sortCompareIndices
{
    sortCompareIndices(const std::vector<double>& data): m_data(data) {}
    bool operator()(int left, int right) const { return m_data[left] < m_data[right]; }
    const std::vector<double>& m_data;
};



};
#endif /* INCLUDED_DIGITAL_OFDM_SMART_RA_H */
