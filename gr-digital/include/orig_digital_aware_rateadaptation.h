/*  	
 * Title: Aware rate adaptation
 * Created By: Owais Khan
 * Creation Date: 08/08/2013
 * 
 * Description: Rate adaptaiton which consider re-transmission till depth-k
 *
*/ 


#ifndef Digital_Aware_RA_11N_H
#define Digital_Aware_RA_11N_H

#include <itpp/itbase.h>

enum Modulation {BPSK, QPSK, QAM16, QAM64};
enum CodeRate {Oneby2, Twoby3, Threeby4, Fiveby6};

double MCS2RATE[] = {6.5, 13, 19.5, 26, 39, 52, 58.5, 65.0};
int CBITSPERMCS[] = {0.5,1.0,1.5,2.0,3.0,4.0,4.5,5.0};


class digital_aware_ra
{
	public:
		digital_aware_ra(int nsubs, int depth, int listlen); // for init

        void read_cberfile();
        void read_fecfiles();
        void read_compressionfiles();
        inline void read_compfile(const char* filename, Modulation mod);
        void read_file(const char* filename, CodeRate crate);
        int select_rate(const std::vector<float> &recvbit_esno,const std::vector<float> &subcarrier_esno, int pkt_len, int pktid, int pkt_mcs, bool pkt_succ);
        void select_retransmission_rate(itpp::vec& bit_esno, itpp::vec& sub_esno, int pktlen);
        void firstpacket_retransmission(itpp::vec &esno_values, int firstmcs, int pktlen);
        void select_firsttransmission_rate(itpp::vec &esno_values, int pktlen);
        void first_transmission(itpp::Array<itpp::ivec> &mcslist, itpp::vec &tput_list,
            itpp::Array<itpp::vec> &list_esno_map, itpp::Array<itpp::Array<itpp::ivec> > &bitmap_list,
            itpp::Array<itpp::Array<itpp::ivec> > &subs_list,
            itpp::vec &txtime, itpp::vec &esno, int pktlen);


void second_transmission(itpp::Array<itpp::ivec> &curr_mod_list, itpp::Array<itpp::Array<itpp::ivec> > &curr_bitmap,
         itpp::Array<itpp::Array<itpp::ivec> > &curr_subs,
         itpp::Array<itpp::vec> &curr_llr_map, itpp::vec &curr_txtime_list,
         itpp::vec &curr_tput_list, const itpp::ivec &prev_mods, const itpp::vec &prev_llr_map, const itpp::Array<itpp::ivec>  &prev_bitmap,
         const itpp::Array<itpp::ivec>  &prev_subs, 
         double prev_txtime, double prev_tput, itpp::vec &sub_esno, int mcs1, int pktlen, int l, int depth, bool first_flag);

        void calculate_uber(itpp::vec &esno, itpp::vec &uber, Modulation mod);
        double calculate_avgber(itpp::vec &uber, itpp::vec &num_per_ber);
        double calculate_fecber(double ber, int mcs);
        void get_rate_bitmap(int mcs, itpp::ivec &bitmap,int pktlen);
        void assign_llr_per_bit(itpp::vec &esno_values, itpp::ivec &bitmap, itpp::vec &llr_bitmap, Modulation mod, int baselen);
        inline void generate_interleave_pattern(itpp::ivec &inlv_indices, Modulation mod);
        void assign_bits_to_subs(itpp::ivec& bitmap, int mcs);
        double sort_subcarriers(itpp::vec& llr_values, itpp::ivec& sorted_sub_indices,int firsttx_mcs);


 
        void assign_llr_per_bit_retx(itpp::vec &llr_values, itpp::ivec &bitmap, itpp::vec &llr_bitmap, Modulation mod, int baselen);
 
        double calculate_deliveryratio(itpp::vec &esno_map, Modulation mod, int mcs, int pktlen);
        void find_unique_values(itpp::vec &values, itpp::vec &unique_values, itpp::vec &val_count);

        double calculate_transmissiontime(int payload, Modulation mod, int mcs, bool ack_flag);
        double calculate_feedback_overhead(int ofdm_num, Modulation mod);
        double calculate_additional_time(int payload, double dr, itpp::vec& esno);

        void remap_esno_values(itpp::vec &esno_values, itpp::vec &dest_esno, Modulation first_mod, Modulation dest_mod);

        void find_sorted_snr_groups(itpp::vec &esno_values, itpp::vec &unique_values, itpp::Array<itpp::ivec> &bits_per_group);
        void quantize_llr_values(itpp::vec& llr_values);

        void calculate_cber_for_esno(itpp::vec& subcarrier_esno);

        std::vector<int> get_curr_bitmap(){
            itpp::ivec tmp_bitmap = selected_bitmap(tx_count-1);
            std::vector<int> bitmap(tmp_bitmap.size());
            for(int i=0; i< tmp_bitmap.size();i++){
                bitmap[i] = tmp_bitmap[i];
            }
        return bitmap;   
        }

        std::vector<int> get_selected_subcarriers(){
            itpp::ivec tmp_subs = selected_subs(tx_count-1);
            std::vector<int> subs(tmp_subs.size());
            for(int i=0; i< tmp_subs.size();i++){
                subs[i] = tmp_subs[i];
            }
        return subs;   
        }


	private:
        int maxDepth, maxList;
        int retxmaxDepth;
        unsigned int selected_mcs;
        int curr_pktid;
        Modulation MCS2MODULATION[8];
        CodeRate MCS2CODERATE[8];
        int RETXMOD2MCS[4];

        //itpp::ivec bucket_map;


        itpp::Array<itpp::vec> UBER;
        itpp::Array<itpp::vec> CBER;
        itpp::Array<itpp::ivec> FeedbackBits;
        itpp::Array<itpp::ivec> CompressBits;
        itpp::Array<itpp::vec> LlrsperTx;

        itpp::Array<itpp::Array<itpp::ivec> > Bits2Submap;

        itpp::vec cber_curr_esno;
    
        double tOfdm, tpreamble, tSIFS, tDIFS, tSlottime, MACsize, PHYsize, ACKsize;
        int Nsubs;

        int tx_count;

        itpp::ivec selected_rate_list;
        itpp::Array<itpp::ivec> selected_bitmap;
        itpp::Array<itpp::ivec> selected_subs;
        double selected_tput;
        double feedback_mcs;
        double feedback_compressed_size;


    // some debugging variables
    itpp::vec debug_dr_list;
};
#endif
