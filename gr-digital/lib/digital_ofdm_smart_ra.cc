/*
 * Title: Smart Rate adaptation
 * Created By: Owais Khan                                            
 * Creation Date: 3/10/2015
 *      
 * Description: Partial Combining Aware Rate adaptation 
 *
*/

#include <algorithm>
#include <sys/time.h>
#include "digital_ofdm_smart_ra.h"

digital_ofdm_smart_ra::digital_ofdm_smart_ra(){ }

digital_ofdm_smart_ra::digital_ofdm_smart_ra(int nSubs,int maxdepth,int maxlist,int maxretx):
            digital_ofdm_rate_adaptation(nSubs){

    d_txinfo.pktSize = 0;
    d_txinfo.txCount = 0;
    d_txinfo.nSubs = Nsubs;
    d_txinfo.maxlistLen = maxlist;
    d_txinfo.searchDepth = maxdepth;
    d_txinfo.maxRetx = maxretx;

    interleaver = digital_ofdm_interleaver (Nsubs);
}

/*
*/
void 
digital_ofdm_smart_ra::pick_best_rate(Node& best_node, std::vector<Node>& nodeList){

    double max_tput = -9e99;
    int best_idx = -1;
    for (int i=0; i< nodeList.size(); i++){
        if( nodeList[i].tput > max_tput) {
            best_idx = i; 
            max_tput = nodeList[i].tput;
        }
    }

    if (best_idx == -1){
        std::cout<<"No node selected. This should never happen."<<std::endl;
        exit(0);
    }
    best_node = nodeList[best_idx];
}


/*
*/
void 
digital_ofdm_smart_ra::update_node_list(std::vector<Node>& list, Node& node, TxInfo& info){


    if (list.size() < info.maxlistLen){
        list.push_back(node);
    }
    else{
        // find worst node
        double min_tput = list[0].tput;
        int min_idx = 0;
        for(int i=1; i<list.size(); i++){
            if (list[i].tput < min_tput){
                min_tput  = list[i].tput;
                min_idx = i;
            }
        }

        if (node.tput > min_tput){
            list[min_idx] = node;
        }
    }

}



/*
*/
int 
digital_ofdm_smart_ra::get_rate_bitmap_len(int mcs, int pktlen, TxInfo& info){

    double const1[] = {2.0,2.0,4.0,2.0,4.0,3.0,4.0,6.0};
    double const2[] = {1.0,1.0,3.0,1.0,3.0,2.0,3.0,5.0};

    int bitmap_size = ceil(((pktlen*8)+6)*(const1[mcs]/const2[mcs]));

    Modulation mod = MCS2MODULATION[mcs];
    int Nbpsc = BITSPERMOD[mod];
    int Ncbps = Nbpsc*info.nSubs;
    int bit_pad = Ncbps - (bitmap_size%Ncbps);
    if (bit_pad == Ncbps){ bit_pad = 0; }
   
    bitmap_size +=bit_pad;

    return bitmap_size;
}

/*
*/
void 
digital_ofdm_smart_ra::remap_esno_values(std::vector<double>& remapped_esno, std::vector<double>& EsNo, Modulation mod1, Modulation mod2, TxInfo& info){

    double const0_table[][4] = { {1.0, 0.5, 0.5, 0.5}, {2.0, 1.0, 1.0, 1.0}, {5.0, 5.0, 1.0, 5.0}, {21.0, 21.0, 21.0, 1.0} };
    double const1_table[][4] = { {1.0, 1.0, 0.75, 0.5833}, {1.0, 1.0, 0.75, 0.5833}, {1.333, 1.333, 1.0, 1.333}, {1.7143, 1.7143, 1.7143, 1.0} };
    double const2_table[][4] = { {1.0, 1.0, 1.0, 1.0}, {1.0, 1.0, 1.0, 1.0}, {1.0, 1.0, 1.0, 0.5833}, {1.0, 1.0, 0.75, 1.0} };
    double const3_table[][4] = { {1.0, 1.0, 5.0, 21.0}, {0.5, 1.0, 5.0, 21.0}, {0.5, 1.0, 1.0, 21.0}, {0.5, 1.0, 5.0, 1.0} };

    double const0 = const0_table[mod1][mod2];
    double const1 = const1_table[mod1][mod2];
    double const2 = const2_table[mod1][mod2];
    double const3 = const3_table[mod1][mod2];

    for(int i=0; i<EsNo.size(); i++){
        double remap_esno = const0*pow( qfuncInv((const1)*(const2)*itpp::Qfunc(sqrt(EsNo[i]/const3) ) ), 2);
        if (remap_esno > 1.0e10){
            remap_esno = 300.0;
        }
        remapped_esno.push_back(remap_esno);
    }

}


/*
*/
void 
digital_ofdm_smart_ra::assign_bits_to_subs(std::vector<std::vector<int> >& bit2submap, int bitmap_size, int mcs, TxInfo& info, bool inlvFlag){

    std::vector<int> inlvBitmap;
    std::vector<int> bitmap;

    Modulation mod = MCS2MODULATION[mcs];
    int Nbpsc = bits_per_modulation(mod);
    int Ncbps = Nbpsc*info.nSubs;
    int bit_pad = Ncbps - (bitmap_size%Ncbps);
    if (bit_pad == Ncbps){ bit_pad = 0; }

    for(int i=0; i<bitmap_size+bit_pad; i++){ bitmap.push_back(i);}

    if (inlvFlag==true){
        interleaver.interleave_data(bitmap,inlvBitmap, mcs);
    }
    else{
        inlvBitmap = bitmap;
    }

    int sub_value = 0;
    int cntr = 0;

    bit2submap = std::vector<std::vector<int> >(info.nSubs);

    for (int i=0; i<inlvBitmap.size();i++){
        sub_value = floor(cntr/Nbpsc);
        if (sub_value >= info.nSubs){
            sub_value = 0;
            cntr = 0;
        }
        bit2submap[sub_value].push_back(inlvBitmap[i]);
        cntr++;
    }
}   // end of assign bit to subs

/*
*/
void 
digital_ofdm_smart_ra::assign_esno_per_bit(std::vector<double>& esno_bitmap, std::vector<double>& EsNo, int tx_bits_size, Modulation mod, TxInfo& info){

    std::vector<int> bitmap;
    std::vector<int> inlvBitmap;

    int Nbpsc = bits_per_modulation(mod);
    int Ncbps = Nbpsc*info.nSubs;
    int bit_pad = Ncbps - (tx_bits_size%Ncbps);
    if (bit_pad == Ncbps){ bit_pad = 0; }

    for(int i=0; i<tx_bits_size+bit_pad; i++){ bitmap.push_back(i);}

    //FIXME: Changing modulation to mcs to fit interleaver implementation. Need to check if this will not cause any problems
    int mcs = MOD2MCS[mod];
    interleaver.interleave_data(bitmap,inlvBitmap, mcs);

    int sub_value = 0;
    int cntr = 0;

    esno_bitmap.resize(inlvBitmap.size());
    for (int i=0; i<inlvBitmap.size(); i++){

        sub_value = floor(cntr/Nbpsc);
        if (sub_value >= info.nSubs){
            sub_value = 0;
            cntr = 0;
        }

        esno_bitmap[inlvBitmap[i]] = EsNo[sub_value];
        cntr++;
    }
}

void digital_ofdm_smart_ra::assign_esno_per_bit_retx(std::vector<double>& esno_bitmap, std::vector<double>& EsNo, int tx_bits_size, Modulation mod1, Modulation mod2, TxInfo& info){

    int Nbpsc = bits_per_modulation(mod2);
    int Ncbps = Nbpsc*info.nSubs;

    std::vector<double> remapped_esno;

    if (mod1 == mod2){
        remapped_esno = EsNo;
    }
    else{
        remap_esno_values(remapped_esno, EsNo, mod1, mod2,info);
    }


    int sub_value = 0;
    int cntr = 0;
    for(int i=0; i<tx_bits_size; i++){
        sub_value = floor(cntr/Nbpsc);
        if (sub_value >= info.nSubs){
            sub_value = 0;
            cntr = 0;
        }
        esno_bitmap.push_back(remapped_esno[sub_value]);
        cntr++;

    }
}



/*
*/
double 
digital_ofdm_smart_ra::calculate_transmission_time(int mcs, int payload, TxInfo& info){

    int num_subs = info.nSubs;
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
double 
digital_ofdm_smart_ra::calculate_retransmission_time(int payload, int mod2, int mcs, TxInfo& info){

    int num_subs = info.nSubs;
    Modulation mod = MCS2MODULATION[mcs];
    int bits_per_symbol = BITSPERMOD[mod2];
    double numOfdmSymbols = ceil(payload/(bits_per_symbol*num_subs));

    double tpayload = tOfdm*numOfdmSymbols;

    double cbits_per_symbol = CBITSPERMCS[mcs];
    double base_bits_per_symbol = 0.5;
    double tphyheader = tOfdm*ceil(PHYsize/(base_bits_per_symbol*num_subs));
    double tmacheader = tOfdm*ceil(MACsize/(cbits_per_symbol*num_subs));

//    double numNackOfdmSymbols = ceil(NACKsize/(cbits_per_symbol*num_subs));
//    double tnack = tpreamble + tphyheader + tOfdm*numNackOfdmSymbols;
    double tnack = tOfdm;
#if WIFITIME == 1
    double txtime = tpreamble + tSIFS + tDIFS + tpayload + tphyheader + tmacheader + tnack;
#else
    double txtime = tpreamble + tpayload + tphyheader + tmacheader + tnack;
#endif

    return txtime;

}

/*
*/ 
void
digital_ofdm_smart_ra::calculate_uber(std::vector<double>& uber, std::vector<double>& esno, std::vector<int>& tx_bits, Modulation mod){
    
    double const1[] ={0.5, 1, 5.0, 21.0};
    double const2[] ={1, 1, 0.75, 0.5833};
 
    uber.resize(tx_bits.size());

    for(int i=0; i<tx_bits.size(); i++){
        uber[i] = const2[mod]*itpp::Qfunc(sqrt(esno[tx_bits[i] ]/const1[mod]));
    }
}


/*
*/
double 
digital_ofdm_smart_ra::calculate_delivery_ratio_esno(Node& node , std::vector<double>& esno_bitmap, std::vector<int>& tx_bits, int mcs, TxInfo& info){

    Modulation mod = MCS2MODULATION[mcs];

    std::vector<double> uber;
    calculate_uber(uber, esno_bitmap,tx_bits, mod);
    node.comb_uber = uber;

    double avg_ber=0;
    for(int i=0; i<uber.size(); i++){
        avg_ber+=uber[i];
    }
    avg_ber = avg_ber/tx_bits.size();
    
    double effesno = calculate_effesno(avg_ber, mod);

    double dr = calculate_delivery_ratio(effesno, mcs, info.pktSize);

    return dr;
    
}

/*
*/
void 
digital_ofdm_smart_ra::retransmission(Node& node, double base_time, int depth, int retx_mod, int retx_bits,std::vector<int>& retransmitted_bits, std::vector<int>& bit_indices, std::vector<double>& EsNo, TxInfo& info){

    int mcs = node.txmcs[0];
    Modulation mod2 = Modulation (retx_mod);
    Modulation mod = MCS2MODULATION[mcs];
    std::vector<double> esno_bitmap;
    assign_esno_per_bit_retx(esno_bitmap, EsNo, retransmitted_bits.size(), mod, mod2, info);

    for(int i=0; i<esno_bitmap.size(); i++){
        node.comb_esno_bit[retransmitted_bits[i]] += esno_bitmap[i];
    }

    std::vector<double> uber;
    calculate_uber(uber, node.comb_esno_bit,retransmitted_bits, mod);

    for(int i = 0; i<uber.size(); i++){
        node.comb_uber[retransmitted_bits[i]] = uber[i];
    }

    if(info.txCount==1){
    //print_vec(node.comb_esno_bit);
    //    print_vec(esno_bitmap);
    }

    double avg_ber=0;
    for(int i=0; i<node.comb_uber.size(); i++){
        avg_ber+=node.comb_uber[i];
    }
    avg_ber = avg_ber/node.comb_uber.size();
    
    double effesno = calculate_effesno(avg_ber, mod);

    node.deliveryRatio = calculate_delivery_ratio(effesno, mcs, info.pktSize);

    // calculate transmission time
    double txtime = calculate_retransmission_time(retx_bits, mod2, mcs, info);
    node.txTime = base_time + txtime;



    int payload = info.pktSize*8;
    node.tput = node.deliveryRatio*payload/node.txTime;
    node.txmcs[depth] = MOD2MCS[mod2];
    node.tx_bits[depth].assign(bit_indices.begin(),bit_indices.begin()+retx_bits-1);

//    std::cout<<"mcs: "<<mcs<<" tx_count: "<<info.txCount<<" mod2: "<<mod2<<" bits: "<<node.tx_bits[depth].size()
//             <<" avg-ber: "<<avg_ber<<" txtime: "<<node.txTime<<" dr: "<<node.deliveryRatio<<" tput: "<<node.tput<<std::endl;
}

/*
*
*/
void
digital_ofdm_smart_ra::first_transmission(Node& first_tx_node, std::vector<double>& EsNo, int mcs, TxInfo& info){

    Modulation mod = MCS2MODULATION[mcs];

    // get rate bitmap
    std::vector<int> tx_bits;
    std::vector<int> tx_subs;
    //get_rate_bitmap(tx_bits, mcs, info.pktSize);
    int tx_bits_len = get_rate_bitmap_len(mcs,info.pktSize, info);
    for(int i=0; i<tx_bits_len; i++){ tx_bits.push_back(i); }
    for(int i=0; i<Nsubs; i++) { tx_subs.push_back(i); }

    // assign bits to subcarriers
    std::vector<std::vector<int> > bit2submap;
    assign_bits_to_subs(bit2submap, tx_bits_len, mcs, info, true);

    // assign esno per bit
    std::vector<double> esno_bitmap;
    assign_esno_per_bit(esno_bitmap, EsNo, tx_bits_len, mod, info);

    // calculate delivery ratio
    double dratio = calculate_delivery_ratio_esno(first_tx_node, esno_bitmap, tx_bits, mcs, info);
    //std::cout<<"dratio: "<<dratio<<std::endl;

    // calculate transmission time
    int payload = info.pktSize*8;
    double txtime = calculate_transmission_time(mcs, payload, info);

    double tput = dratio*payload/txtime;

    // finalize node and return
    first_tx_node.tx_bits = std::vector<std::vector<int> > (info.maxRetx+1);
    first_tx_node.tx_bits[0] = tx_bits;
    first_tx_node.tx_subs = std::vector<std::vector<int> > (info.maxRetx+1);
    first_tx_node.tx_subs[0] = tx_subs;
    first_tx_node.txmcs = std::vector<int> (info.maxRetx+1,-1);
    first_tx_node.txmcs[0] = mcs;
    first_tx_node.bit2submap = bit2submap;
    first_tx_node.comb_esno_bit = esno_bitmap;
    first_tx_node.deliveryRatio = dratio;
    first_tx_node.txTime = txtime;
    first_tx_node.tput = tput;

}

void 
digital_ofdm_smart_ra::linear_search(std::vector<Node>& tmpList, Node curr_node, int depth, int retx_mod, std::vector<int>& sub_indices, std::vector<int>& bit_indices, std::vector<double>& EsNo, TxInfo& info){

    int max_reached = 0;
    std::vector<double> tput_list;
    int retx_bits = 0;

    double max_tput=-9e99;
    int dec_cntr=0;
    double base_time = curr_node.txTime;
    for( int i =0; i<info.nSubs; i++){

        retx_bits += curr_node.bit2submap[sub_indices[i]].size();
        retransmission(curr_node,base_time,depth,retx_mod,retx_bits,curr_node.bit2submap[sub_indices[i]],bit_indices, EsNo, info);
        (curr_node.tx_subs[depth]).push_back(sub_indices[i]);
        update_node_list(tmpList, curr_node, info);

        // check for max
        if (curr_node.tput < max_tput){ 
            dec_cntr++; 
        }
        else{
            max_tput = curr_node.tput;
        }
        if (dec_cntr == 2) { break; }
    }

}


/*
*/
void 
digital_ofdm_smart_ra::search_retransmission_combination(std::vector<Node>& tmpList, Node& curr_node, int depth, std::vector<double>& EsNo,TxInfo& info){

    Modulation mod = MCS2MODULATION[curr_node.txmcs[0]];

    // find effsnr per subcarrier
    std::vector<double> uber_per_bit;
    
    calculate_uber(uber_per_bit,curr_node.comb_esno_bit,curr_node.tx_bits[0],mod);
    curr_node.comb_uber = uber_per_bit;

    std::vector<double> effEsNoPerSub(info.nSubs);
    std::vector<int> sub_indices;
    std::vector<int> bit_indices;

    for(int i=0; i<info.nSubs; i++){
        double avg_ber = 0;
        std::vector<int> sub_bitmap = curr_node.bit2submap[i];

        for(int j=0; j < sub_bitmap.size(); j++){
            avg_ber += uber_per_bit[sub_bitmap[j]];
        }

        avg_ber = avg_ber/sub_bitmap.size();
        effEsNoPerSub[i] = calculate_effesno(avg_ber, mod);
        sub_indices.push_back(i);
    }

    std::sort(sub_indices.begin(), sub_indices.end(), sortCompareIndices(effEsNoPerSub));

    std::vector<int>* vec_ptr;
    for(int i=0; i<info.nSubs; i++){
        vec_ptr = &curr_node.bit2submap[sub_indices[i]];
        for(int j=0; j < vec_ptr->size(); j++){
            bit_indices.push_back(vec_ptr->at(j));
        }
    }

    for(int mod_idx = 0; mod_idx < 4; mod_idx++){

        linear_search(tmpList, curr_node, depth, mod_idx, sub_indices, bit_indices, EsNo, info);

    }

}


/*
*
*/
void 
digital_ofdm_smart_ra::find_rate_firsttransmission(Node& bestNode, std::vector<std::vector<double> >& EsNo, TxInfo& info){

    std::vector<double> firstTxEsNo = EsNo[0];
    // first transmission rate selection
    for (int crate_index=numRates-1; crate_index >=0; crate_index--){
//    for (int crate_index=numRates-1; crate_index >=numRates-1; crate_index--){

//        int first_tx_mcs = 4;
        int first_tx_mcs = crate_index;
        Node first_tx_node;
        first_transmission(first_tx_node, firstTxEsNo, first_tx_mcs, info); 
        update_node_list(nodeList, first_tx_node, info);

        if(first_tx_node.deliveryRatio > 0.99){ break; }
    }

    // second transmission
    std::vector<Node> tmpList(nodeList);
    Node curr_node;
    for (int depth = 1; depth < info.searchDepth; depth++){

        for( int node = 0; node < nodeList.size(); node++){

            curr_node = nodeList[node];
            search_retransmission_combination(tmpList, curr_node, depth, EsNo[depth],info);

        }
        nodeList = tmpList;
    }

    // pick best rate
    pick_best_rate(bestNode, nodeList);

    print_vec(bestNode.txmcs);
    //std::cout<<"tput:"<<bestNode.tput<<std::endl;
    nodeList.clear();
}

/*
*/
void 
digital_ofdm_smart_ra::find_rate_retransmission(Node& bestNode, std::vector<std::vector<double> >& EsNo, TxInfo& info){

    int tx_count = info.txCount;
    nodeList.push_back(bestNode);
    std::vector<Node> tmpList(nodeList);

    Node curr_node;
    //for(int i=0;i<100;i++){printf("%f, ",bestNode.comb_esno_bit[i]);}
    //printf("\n");
    //std::cout<<" tx_count: "<<tx_count<<std::endl;
    for (int depth = tx_count; depth < info.searchDepth +tx_count; depth++){
        for( int node = 0; node < nodeList.size(); node++){
            // std::cout<<"depth: "<<depth<<" searchDepth: "<<info.searchDepth<<std::endl;
            curr_node = nodeList[node];
            search_retransmission_combination(tmpList, curr_node, depth, EsNo[depth],info);
        }
        nodeList = tmpList;
    }
    // pick best rate
    pick_best_rate(bestNode, nodeList);

    print_vec(bestNode.txmcs);
    //std::cout<<"tput:"<<bestNode.tput<<std::endl;
    nodeList.clear();
}

/*
*/
void 
digital_ofdm_smart_ra::reboot_transmission(Node& node, TxInfo& info, int mcs){

    Modulation mod = MCS2MODULATION[mcs];

    // get rate bitmap
    std::vector<int> tx_bits;
    std::vector<int> tx_subs;
    //get_rate_bitmap(tx_bits, mcs, info.pktSize);
    int tx_bits_len = get_rate_bitmap_len(mcs,info.pktSize, info);
    for(int i=0; i<tx_bits_len; i++){ tx_bits.push_back(i); }
    for(int i=0; i<Nsubs; i++) { tx_subs.push_back(i); }

    // assign bits to subcarriers
    std::vector<std::vector<int> > bit2submap;
    assign_bits_to_subs(bit2submap, tx_bits_len, mcs, info, true);

    // calculate transmission time
    int payload = info.pktSize*8;
    double txtime = calculate_transmission_time(mcs, payload, info);

    // finalize node and return
    node.tx_bits = std::vector<std::vector<int> > (info.maxRetx+1);
    node.tx_bits[0] = tx_bits;
    node.tx_subs = std::vector<std::vector<int> > (info.maxRetx+1);
    node.tx_subs[0] = tx_subs;
    node.txmcs = std::vector<int> (info.maxRetx+1,-1);
    node.txmcs[0] = mcs;
    node.bit2submap = bit2submap;
    node.deliveryRatio = 0;
    node.txTime = txtime;
    node.tput = 0;
    node.comb_esno_bit.clear(); // Updated when update_combesno() called.
    node.comb_uber.clear();     // Same

}


/*
*/
void
digital_ofdm_smart_ra::select_rate_internal(Node& active_node, std::vector<std::vector<double> >& EsNo, TxInfo& info){

    double t0 = current_time();
    if (info.txCount == 0){
        find_rate_firsttransmission(active_node, EsNo, info);
    }
    else{
        find_rate_retransmission(active_node, EsNo, info);
    }

    double time_elapsed = current_time() - t0;
    std::cout<<"time elapsed: "<<time_elapsed*1.0e3<<"ms"<<std::endl;
}


/*
*/
int
digital_ofdm_smart_ra::select_rate(std::vector<double>& subsEsNo, itpp::vec& combllr, std::vector<int>& retx_subs,int txcount, int mcs, int pktSize, int pktId, int succ){


    if(txcount>=(int)d_txinfo.maxRetx){
        printf("This is awkward. This shouldn't happen.\n");
        printf("The maxRetx numbers should match.\n");
        printf("I am not dealing with this right now.\n");
        printf("Calling exit(0).\n");
        exit(0);
    }

    d_txinfo.pktSize = pktSize;
    d_activeNode.deliveryRatio = 0;
    d_activeNode.tput = 0;

    // Case 1: succ==1
    // Calculate for a new transmission

    // Case 2: succ==0
        // a: txcount=0
        //  New transmission, first tx unsuccessful
        //  If mcs is not expected mcs, recalculate

        // b: txcount > 0
        //  if txcount is expected txcount count then contiune
        //  if unexpected bootstrap

    if (succ == 1){
        setup_newtransmission(subsEsNo, d_txinfo, pktId);
        d_pktId = pktId+1;
        d_txinfo.txCount = 0;
    }
    else{

        // check to see if the rxed packet is expected
        if((d_pktId != pktId)&&(txcount!=0)){
            printf("OK. First tx was dropped.\n");
            printf("Will need to bootstrap the rate search.\n");
            printf("For now. Do not allow it.\n");
            printf("Calling exit(0).\n");
            exit(0);
            setup_newtransmission(subsEsNo, d_txinfo, pktId);

        }
        else{

            if(txcount == 0){

                if(d_pktId != pktId){
                    d_pktId = pktId;
                }

                if(d_activeNode.txmcs[0]!=mcs){
                    printf("Transmission mcs does not match expected mcs.\n");
                    printf("PktSize: %d\n",d_txinfo.pktSize);
                    reboot_transmission(d_activeNode, d_txinfo, mcs);
                    printf("Values rebooted.\n");
                }
            }
            else if(d_txinfo.txCount!=(txcount+1)){
                // combllr conversion remains the same
                // need to fill up d_EsNo for the missing packets
                // For now use the current EsNo for missing packets
                for(int i=d_EsNo.size(); i<(txcount+d_txinfo.searchDepth);i++){
                    d_EsNo.push_back(subsEsNo);
                }
            }

            d_txinfo.txCount = txcount+1;
            update_combesno(d_activeNode, combllr);

            for(int i=0; i<d_txinfo.searchDepth;i++){
                d_EsNo[txcount+i] = subsEsNo;
            }
            d_EsNo.push_back(subsEsNo);

        }
    } // end else succ=1


    // call select rate
    select_rate_internal(d_activeNode, d_EsNo, d_txinfo);

    // set retx subs
    retx_subs.resize(d_txinfo.nSubs);
    for(int i=0; i<d_txinfo.nSubs; i++){
        retx_subs[i]=0;
    }

    for(int i=0; i<d_activeNode.tx_subs[d_txinfo.txCount].size(); i++){
        retx_subs[d_activeNode.tx_subs[d_txinfo.txCount][i] ] = 1;
    }
    int ra_mcs = d_activeNode.txmcs[d_txinfo.txCount];
    return ra_mcs;

}

/*
*/
void 
digital_ofdm_smart_ra::update_combesno(Node& node, itpp::vec& combllr){

    // convert llr values to EsNo values
    std::vector<double> effesno_per_bit;
    double pe,effesno;
    // FIXME: double check if first mcs value should be used.
    // Also above FIXME, what they are different.

    Modulation mod = MCS2MODULATION[node.txmcs[0]];
    for(int i=0; i<combllr.size(); i++){
        pe = 1.0/(1+exp(abs(combllr[i])));
        effesno = calculate_effesno(pe, mod);
        effesno_per_bit.push_back(effesno);
    }

    // Update activeNode.comb_esno_bit
    node.comb_esno_bit = effesno_per_bit;
    std::vector<double> uber;
    calculate_uber(uber, node.comb_esno_bit, node.tx_bits[0], mod);
    node.comb_uber = uber;

    //print_vec(effesno_per_bit);
    //print_vec(uber);

}

/*
*/
void 
digital_ofdm_smart_ra::setup_newtransmission(std::vector<double>& subsEsNo, TxInfo& info, int pktId){

    
    d_EsNo = std::vector<std::vector<double> >(0);
    // set up esno values and start rate selection process
    for(int i=0; i<info.searchDepth; i++){
        d_EsNo.push_back(subsEsNo);
    }

    clear_activenode(d_activeNode);


}

//
void 
digital_ofdm_smart_ra::clear_activenode(Node& node){

    node.deliveryRatio = 0;
    node.tput = 0;
    node.comb_esno_bit.clear();
    node.bit2submap.clear();
    node.txmcs.clear();
    node.tx_bits.clear();
    node.tx_subs.clear();
    node.comb_uber.clear();
}
