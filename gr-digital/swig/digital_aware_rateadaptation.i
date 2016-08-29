class digital_aware_ra
{
	public:
		digital_aware_ra(int nsubs, int depth, int listlen); // for init
        int select_rate(const std::vector<float> &recvbit_esno,const std::vector<float> &subcarrier_esno, int pkt_len, int pktid, int pkt_mcs, bool pkt_succ);
        std::vector<int> get_curr_bitmap();
        std::vector<int> get_selected_subcarriers();
};

