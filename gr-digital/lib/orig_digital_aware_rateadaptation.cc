/*
 * Title: Aware rate adaptation
 * Created By: Owais Khan
 * Creation Date: 08/08/2013
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

#define FECBER 1

using namespace itpp;

digital_aware_ra::digital_aware_ra(int nsubs, int depth, int listlen):selected_mcs(0),
curr_pktid(-1),tpreamble(160.0e-6),tSIFS(16.0e-6),tSlottime(9.0e-6),
tOfdm(80.0e-6),Nsubs(nsubs),maxDepth(depth),maxList(listlen),
tx_count(0){

    retxmaxDepth = maxDepth;
    MACsize = 28*8;
    ACKsize = 14*8;
    //PHYsize = 22;
    PHYsize = 10*8;
    tDIFS = tSIFS + 2*tSlottime;

    
    #if FECBER == 1
    UBER.set_size(4);
    CBER.set_size(4);
    #else
    UBER.set_size(8);
    CBER.set_size(8);
    #endif

    Bits2Submap.set_size(8);

    cber_curr_esno.set_size(8);
    FeedbackBits.set_size(4);
    CompressBits.set_size(4);

    MCS2MODULATION = {BPSK, QPSK, QPSK, QAM16, QAM16, QAM64, QAM64, QAM64};
    MCS2CODERATE = {Oneby2, Oneby2, Threeby4, Oneby2, Threeby4, Twoby3, Threeby4, Fiveby6};
 
    RETXMOD2MCS = {0,1,3,5};

    //bucket_map.set_size(20000);
    
    #if FECBER == 1
    read_fecfiles();
    #else
    read_cberfile();
    #endif

    read_compressionfiles();
}

inline
double current_time(){
      struct timeval tv;
      if ( gettimeofday(&tv, 0) < 0 ){
        std::cout<<"gettimeofday failed in current_time()"<<std::endl;
        exit(0);
      }
    
      return double(tv.tv_sec) + double(tv.tv_usec) / 1e6;
}

inline Modulation int2modulation(int val){
    
    Modulation mod;
    switch(val){
        case 0: mod = BPSK; break;
        case 1: mod = QPSK; break;
        case 2: mod = QAM16; break;
        case 3: mod = QAM64; break;
        default: std::cout<<"passed int value doesnot match modulation."<<std::endl;
    }
return mod;
}


inline int rate2mcs(std::string rate){

    int mcs=-1;

    if (strcmp(rate.data(),"65") == 0){mcs = 7;}
    else if(strcmp(rate.data(),"58.5") == 0){ mcs = 6;}
    else if(strcmp(rate.data(),"52") == 0){ mcs = 5;}
    else if(strcmp(rate.data(),"39") == 0){ mcs = 4;}
    else if(strcmp(rate.data(),"26") == 0){ mcs = 3;}
    else if(strcmp(rate.data(),"19.5") == 0){ mcs = 2;}
    else if(strcmp(rate.data(),"13") == 0){ mcs = 1;}
    else if(strcmp(rate.data(),"6.5")== 0){ mcs = 0;}
    else{
        std::cout<<"rate doesnot match any mcs."<<std::endl;
        exit(0);
    }
    return mcs;
}

inline int bits_per_modulation(Modulation mod){

    int bits = 0;

    switch(mod){
        case BPSK : bits = 1; break;
        case QPSK : bits = 2; break;
        case QAM16: bits = 4; break;
        case QAM64: bits = 6; break;
        default: std::cout<<"Incorrect modulation."<<std::endl;
    }

return bits;
}

inline void 
digital_aware_ra::read_compfile(const char* filename, Modulation mod){

    std::ifstream compfile(filename);
    std::string line;
   
    itpp::ivec feedback_bits_vec;
    itpp::ivec compress_vec;

    if(compfile.is_open()){
        while (compfile.good() ){
            getline(compfile,line);
            
            if(compfile.eof()){break;}

            size_t location = line.find(", ");
            std::string uber_str = line.substr(0,location);
            double uber = std::atoi(uber_str.c_str());

            std::string line2 = line.substr(location+1);
            double cber = std::atoi(line2.c_str());

            feedback_bits_vec.ins(feedback_bits_vec.size(),uber);
            compress_vec.ins(compress_vec.size(),cber);

        }
        FeedbackBits(mod) = feedback_bits_vec;
        CompressBits(mod) = compress_vec;
    }
    else{
        std::cout<<"Failed to open "<<filename<<std::endl;
        exit(0);
    }

    compfile.close();
}


void digital_aware_ra::read_compressionfiles(){

    read_compfile("bpskcompress_data.dat",BPSK);
    read_compfile("qpskcompress_data.dat",QPSK);
    read_compfile("qam16compress_data.dat",QAM16);
    read_compfile("qam64compress_data.dat",QAM64);
    
}




void 
digital_aware_ra::read_file(const char* filename, CodeRate crate){

    std::ifstream cberfile(filename);
    std::string line;
   
    itpp::vec uber_vec;
    itpp::vec cber_vec;

    if(cberfile.is_open()){
        while (cberfile.good() ){
            getline(cberfile,line);
            
            if(cberfile.eof()){break;}

            size_t location = line.find(", ");
            std::string uber_str = line.substr(0,location);
            double uber = std::atof(uber_str.c_str());

            std::string line2 = line.substr(location+1);
            double cber = std::atof(line2.c_str());

            uber_vec.ins(uber_vec.size(),uber);
            cber_vec.ins(cber_vec.size(),cber);

        }
        UBER(crate) = uber_vec;
        CBER(crate) = cber_vec;
    }
    else{
        std::cout<<"Failed to open "<<filename<<std::endl;
        exit(0);
    }

    cberfile.close();
}

void digital_aware_ra::read_fecfiles(){

    read_file("rate1by2.data", Oneby2);
    read_file("rate2by3.data", Twoby3);
    read_file("rate3by4.data", Threeby4);
    
}

void digital_aware_ra::read_cberfile(){

    std::ifstream cberfile("cberdata11n.dat");
    std::string line;
   
    int prev_mcs = 7; //assuming file starts from rate 65
    itpp::vec uber_vec;
    itpp::vec cber_vec;

    if(cberfile.is_open()){
        while (cberfile.good() ){
            getline(cberfile,line);
            
            if(cberfile.eof()){break;}

            size_t location = line.find(" ");
            std::string rate = line.substr(0,location);
            int mcs = rate2mcs(rate);

            if (mcs!=prev_mcs){
                UBER(prev_mcs) = uber_vec;
                CBER(prev_mcs) = cber_vec;
                uber_vec.set_size(0);
                cber_vec.set_size(0);
                prev_mcs = mcs;
            }

            std::string line2 = line.substr(location+1);
            size_t location2 = line2.find(" ");
            //double snr = std::atof(line2.substr(0,location2).c_str());
    
            std::string line3 = line2.substr(location2+1);
            size_t location3 = line3.find(" ");
            double uber = std::atof(line3.substr(0,location3).c_str());
    
            std::string line4 = line3.substr(location3+1);
            double cber = std::atof(line4.c_str());

            uber_vec.ins(uber_vec.size(),uber);
            cber_vec.ins(cber_vec.size(),cber);

        }
        UBER(0) = uber_vec;
        CBER(0) = cber_vec;
    }

    cberfile.close();
}


void
digital_aware_ra::calculate_uber(itpp::vec &esno, itpp::vec &uber, Modulation mod){


    switch(mod){

        case BPSK:
            uber = itpp::Qfunc(itpp::sqrt(2*esno));
            break;
        case QPSK:
            uber = itpp::Qfunc(itpp::sqrt(esno));
            break;

        case QAM16:
            //uber = (1/4.0)*itpp::Qfunc(itpp::sqrt((esno/15.0)));
            uber = (3/4.0)*itpp::Qfunc(itpp::sqrt((esno/5.0)));
            break;

        case QAM64:
            //uber = (7/26.0)*itpp::Qfunc(itpp::sqrt((esno/63.0)));
            uber = (7/12.0)*itpp::Qfunc(itpp::sqrt((esno/21.0)));
            break;
        default:
            std::cout<<"Incorrect modulation selection."<<std::endl;

    }
}

double
digital_aware_ra::calculate_avgber(itpp::vec &uber, itpp::vec &num_per_ber){
    int total_bits = itpp::sum(num_per_ber);
    double total_sum = itpp::elem_mult_sum(uber,num_per_ber);
    double avgber = total_sum/total_bits;
return avgber;
}

#if FECBER == 1
double
digital_aware_ra::calculate_fecber(double ber, int mcs){

    CodeRate code_rate = MCS2CODERATE[mcs]; 
    itpp::vec ubertable = UBER(code_rate);
    itpp::vec cbertable = CBER(code_rate);

    int i=0;
    for(; i<ubertable.size(); i++){
        if(ber>=ubertable[i]){
            break;
        }
    }
    int index = std::min(i,ubertable.size()-1);
    double cber =  cbertable[index];

return cber;
}
#else
double
digital_aware_ra::calculate_fecber(double ber, int mcs){

    itpp::vec ubertable = UBER(mcs);
    itpp::vec cbertable = CBER(mcs);

    int i=0;
    for(; i<ubertable.size(); i++){
        if(ber>=ubertable[i]){
            break;
        }
    }
    int index = std::min(i,ubertable.size()-1);
    double cber =  cbertable[index];

return cber;
}
#endif


void
digital_aware_ra::get_rate_bitmap(int mcs, itpp::ivec &bitmap,int pktlen){

    CodeRate coderate = MCS2CODERATE[mcs];

    int baselen = ((pktlen*8)+6)*2;
    
    itpp::ivec basemap(baselen);
    for(int i=0; i<baselen; i++){ basemap[i] = i; }

    itpp::bvec puncture;
    int bitmap_len;

    switch(coderate){

        case Oneby2:
            puncture = "1 1 1 1";
            bitmap_len = baselen;
            break;

        case Twoby3:
            puncture = "1 1 1 0";
            bitmap_len = ceil(((pktlen*8)+6)*(3.0/2));
            break;

        case Threeby4:
            puncture = "1 1 1 0 0 1";
            bitmap_len = ceil(((pktlen*8)+6)*(4.0/3));
            break;

        case Fiveby6:
            puncture = "1 1 1 0 0 1 1 0 0 1";
            bitmap_len = ceil(((pktlen*8)+6)*(6.0/5));
            break;

        default:
            std::cout<<"mcs doesnot match any code rate."<<std::endl;

    }
    int psize = puncture.size();
    bitmap.set_size(bitmap_len);
    bitmap.clear();

    int j=0;
    for(int i=0; i<basemap.size(); i++){
        if(puncture[i%psize]==1){
            bitmap[j] = i;
            j++;
            //bitmap.ins(bitmap.size(),i);
        }
    }
    //std::cout<<"grb -> basemap size: "<<basemap.size()<<std::endl;
    //std::cout<<"grb -> bitmap size: "<<bitmap.size()<<std::endl;

}

void
digital_aware_ra::assign_llr_per_bit_retx(itpp::vec &llr_values, itpp::ivec &bitmap, itpp::vec &llr_bitmap, Modulation mod, int baselen){

    llr_bitmap.set_size(baselen);
    llr_bitmap.zeros();

    //std::cout<<"baselen: "<<baselen<<std::endl;
    //std::cout<<"esno values size: "<<esno_values.size()<<std::endl; 
    //std::cout<<"esno bitmap size: "<<esno_bitmap.size()<<std::endl; 
    //std::cout<<"bitmap size: "<<bitmap.size()<<" largest index: "<<bitmap(bitmap.size()-1)<<std::endl;
    for(int index=0; index < bitmap.size(); index++){
            llr_bitmap[bitmap[index]] = llr_values[bitmap[index]];
    }
}

inline void
get_ber_values(itpp::vec& esno, Modulation mod, int bits_per_symbol, itpp::mat& ber_values){

    ber_values.set_size(esno.size(),bits_per_symbol);

    itpp::vec uber1;
    itpp::vec uber2;
    itpp::vec uber3;

    switch(mod){
        case BPSK:
            uber1 = itpp::Qfunc(itpp::sqrt(2*esno));
            ber_values.set_col(0,uber1);
            break;
        case QPSK:
            uber1 = itpp::Qfunc(itpp::sqrt(esno));
            ber_values.set_col(0,uber1);
            ber_values.set_col(1,uber1);
            break;
        case QAM16:
            uber1 = 0.5*itpp::Qfunc(itpp::sqrt((esno/5.0)));
            uber2 = itpp::Qfunc(itpp::sqrt((esno/5.0)));
            ber_values.set_col(0,uber1);
            ber_values.set_col(1,uber2);
            ber_values.set_col(2,uber1);
            ber_values.set_col(3,uber2);
            break;
        case QAM64:
            uber1 = 0.25*itpp::Qfunc(itpp::sqrt((esno/21.0)));
            uber2 = 0.5*itpp::Qfunc(itpp::sqrt((esno/21.0)));
            uber3 = itpp::Qfunc(itpp::sqrt((esno/21.0)));
            ber_values.set_col(0,uber1);
            ber_values.set_col(1,uber2);
            ber_values.set_col(2,uber3);
            ber_values.set_col(3,uber1);
            ber_values.set_col(4,uber2);
            ber_values.set_col(5,uber3);
            break;
        default:
            std::cout<<"Incorrect modulation selection."<<std::endl;
    }
}

inline void
digital_aware_ra::generate_interleave_pattern(itpp::ivec &inlv_indices, Modulation mod){
    
    int Nbpsc = bits_per_modulation(mod);
    int Ncbps = Nbpsc*Nsubs;
    int s = std::max(Nbpsc/2,1);
 
    inlv_indices.set_size(Ncbps);

    int inlv_index1;
    for(int k=0; k < Ncbps; k++){
        inlv_index1 =(int)( (Ncbps/16)*(k%16)+floor(k/16) );
        inlv_indices[k] = s*floor(inlv_index1/s)+(inlv_index1+Ncbps - (int)floor(16*inlv_index1/Ncbps))%s;
    }
}

void
digital_aware_ra::assign_llr_per_bit(itpp::vec &esno_values, itpp::ivec &bitmap, itpp::vec &llr_bitmap, Modulation mod, int baselen){

    int bits_per_symbol = bits_per_modulation(mod);
    //int num_of_bits = bitmap.size();
    //int numofsymbols = (int) ceil(num_of_bits/bits_per_symbol);
    int num_subs = esno_values.size();
    itpp::ivec order_array;
    switch(mod){
        case BPSK: order_array = "7 6 5 4 3 2 1 0"; break;
        case QPSK: order_array = "6 7 4 5 2 3 0 1"; break;
        case QAM16:order_array = "4 5 6 7 0 1 2 3"; break;
        case QAM64:order_array = "10 11 0 1 2 3 4 5 14 15 16 17 6 7 8 9 18 19 20 21 22 23 12 13"; break;
        default: std::cout<<"assign esno. Invalid mod."<<std::endl;
    }

    itpp::ivec inlv_indices;
    generate_interleave_pattern(inlv_indices,mod);

    itpp::mat ber_values;
    // Optimization not needed to be called every time
    get_ber_values(esno_values, mod, bits_per_symbol,ber_values);

    llr_bitmap.set_size(baselen);
    llr_bitmap.zeros();
    int offset = order_array.size();
    int k=0, sub_index=0, bit_loc=0;

    int bit_index_in_sub=0;
    double est_ber;
    int inlv_index,inlv_size, ofdm_num;

    inlv_size = inlv_indices.size();

    for(int i=0; i < bitmap.size(); i++){
        ofdm_num = (int) i/inlv_size;
        inlv_index = inlv_indices[i%inlv_size] + (ofdm_num*inlv_size);
        k = (int) inlv_index/offset;
        bit_loc = (order_array[inlv_index%offset]) + (k*offset);
        sub_index = (bit_loc/bits_per_symbol)%num_subs;
        bit_index_in_sub = bit_loc%bits_per_symbol;
        //std::cout<<i<<" inlv_index: "<<inlv_index<<" k: "<<k<<" subindex: "<<sub_index<<" bit loc: "<<bit_loc<<" bitmap val: "<<bitmap[i]<<std::endl; 

        // Figure out the number of bits per snr and map them accordingly
        //bit_counter = i%bits_per_symbol; 
        est_ber = ber_values(sub_index,bit_index_in_sub);

        llr_bitmap[bitmap[i]] = log((1-est_ber)/est_ber);
    }
    //std::cout<<"llr_bitmap:"<<llr_bitmap<<std::endl;
}

void
digital_aware_ra::assign_bits_to_subs(itpp::ivec& bitmap, int mcs){

    itpp::Array<itpp::ivec> bits2sub_mcs(Nsubs);

    for(int i=0; i< Nsubs; i++){
        bits2sub_mcs(i) = itpp::ivec(0); 
    }

    Modulation mod = MCS2MODULATION[mcs];
    int bits_per_symbol = bits_per_modulation(mod);

    int num_subs = Nsubs;
    itpp::ivec order_array;
    switch(mod){
        case BPSK: order_array = "7 6 5 4 3 2 1 0"; break;
        case QPSK: order_array = "6 7 4 5 2 3 0 1"; break;
        case QAM16:order_array = "4 5 6 7 0 1 2 3"; break;
        case QAM64:order_array = "10 11 0 1 2 3 4 5 14 15 16 17 6 7 8 9 18 19 20 21 22 23 12 13"; break;
        default: std::cout<<"assign esno. Invalid mod."<<std::endl;
    }

    itpp::ivec inlv_indices;
    generate_interleave_pattern(inlv_indices,mod);

    int offset = order_array.size();
    int k=0, sub_index=0, bit_loc=0;

    //int bit_index_in_sub=0;
    double est_ber;
    int inlv_index,inlv_size, ofdm_num;

    inlv_size = inlv_indices.size();

    for(int i=0; i < bitmap.size(); i++){
        ofdm_num = (int) i/inlv_size;
        inlv_index = inlv_indices[i%inlv_size] + (ofdm_num*inlv_size);
        k = (int) inlv_index/offset;
        bit_loc = (order_array[inlv_index%offset]) + (k*offset);
        sub_index = (bit_loc/bits_per_symbol)%num_subs;
        //bit_index_in_sub = bit_loc%bits_per_symbol;
        //std::cout<<i<<" inlv_index: "<<inlv_index<<" k: "<<k<<" subindex: "<<sub_index<<" bit loc: "<<bit_loc<<" bitmap val: "<<bitmap[i]<<std::endl; 
        bits2sub_mcs(sub_index).ins(bits2sub_mcs(sub_index).size(),bitmap[i]);
    }

    Bits2Submap(mcs) = bits2sub_mcs;

    //std::cout<<"Bits2Submap: "<<Bits2Submap(mcs)<<std::endl;

}



#if 1
void digital_aware_ra::find_unique_values(itpp::vec &values, itpp::vec &unique_values, itpp::vec &val_count){

    std::map<double,int> esno_map;
    std::map<double,int>::iterator it_find;

    //std::cout<<"fuv -> values: "<<values<<std::endl; 
    int curr_cnt = 0;
    for( int i=0; i<values.size(); i++){
        it_find = esno_map.find(values[i]);
        if (it_find != esno_map.end()){
            curr_cnt = (*it_find).second;
            curr_cnt++;
            (*it_find).second = curr_cnt;

        }
        else{
            if(values[i]!=0){
                esno_map[values[i]] = 1;
            }
        }
    }
    int map_size = esno_map.size();
    unique_values.set_size(map_size);
    val_count.set_size(map_size);

    int i =0 ;
    for(std::map<double,int>::iterator it=esno_map.begin(); it!=esno_map.end(); it++){
       double key = (*it).first;
       unique_values(i) = key;
       val_count(i) = (*it).second;
       i++;
       //unique_values.ins(unique_values.size(),key); 
       //val_count.ins(val_count.size(),(*it).second);
       //std::cout<<"fuv -> key: "<<key<<" count: "<<(*it).second<<std::endl;
    }


}
#else
void digital_aware_ra::find_unique_values(itpp::vec &values, itpp::vec &unique_values, itpp::vec &val_count){

    std::set<double> currset;
    for (int i=0; i<values.size(); i++){
        currset.insert(values[i]);
    }

    for(std::set<double>::iterator it=currset.begin(); it != currset.end(); it++){
        unique_values.ins(unique_values.size(),*it);
    }

    if( unique_values[0]==0){
        unique_values.del(0);
    }

    val_count.set_size(unique_values.size());
    for(int i=0; i <unique_values.size(); i++){
        val_count[i] = currset.count(unique_values[i]);
    }
    /*
    for(int u=0; u < unique_values.size();u++){
        int counter = 0;
        for(int i=0; i < values.size(); i++){
            if(unique_values[u] == values[i]){
                counter++;
            }

        }
        val_count[u] = counter;
    }*/
}
#endif

double time_unique = 0;
double time_uber = 0;
double time_fec = 0;
double 
digital_aware_ra::calculate_deliveryratio(itpp::vec &llr_map, Modulation mod, int mcs, int pktlen){

    double t0 = current_time();
    //find_unique_values(comb_esno, unique_esno, num_per_uber);
    double t1 = current_time();
    time_unique+=(t1-t0);

    t0 = current_time();
    itpp::vec abs_llrs = itpp::abs(llr_map);
    //itpp::vec abs_llrs = itpp::abs(comb_llr);
    double uber=0;
    double llr_value=0;
    int uber_count=0;
    for(int i=0; i<abs_llrs.size(); i++){
        llr_value = abs_llrs[i];

        if(llr_value>0){
            uber+= 1/(1+exp(llr_value));
            uber_count++;
        } 
    }
    //std::cout<<"dr-> uber: "<<uber<<std::endl;
    t1 = current_time();
    time_uber+=(t1-t0);

    //double avgber = itpp::sum(uber)/uber.size();
    double avgber = uber/uber_count;
    //std::cout<<"dr-> avgber: "<<avgber<<std::endl;

    // calculate delivery ratio
    t0 = current_time();
    double cber = calculate_fecber(avgber, mcs);
    //std::cout<<"dr-> fecber: "<<cber<<std::endl;
    t1 = current_time();
    time_fec+=(t1-t0);

    double dr = pow((1 - cber),pktlen*8);

    std::cout<<"dr-> uber: "<<uber<<" avg ber: "<<avgber<<" cber: "<<cber<<" dr: "<<dr<<std::endl;
return dr;
}

double
digital_aware_ra::calculate_additional_time(int payload, double dr, itpp::vec& esno){

    double mcs2coderate[] = {2,2,1.33,2,1.33,1.5,1.33,1.2};
    int txbits = ceil(payload*(1-dr));

    double* avg_ber = new double[4];
    for(int i=0; i < 4; i++){
        Modulation mod = int2modulation(i);
        itpp::vec uber;
        calculate_uber(esno, uber, mod);
        avg_ber[i] = itpp::sum(uber)/uber.size();
    }

    double fec_ber, retx_dr;
    double th_dr = 0.99;
    int chosen_mcs = 0;
    for(int mcs = 6; mcs >=0; mcs--){
        Modulation mod = MCS2MODULATION[mcs];
        fec_ber = calculate_fecber(avg_ber[mod],mcs);
        retx_dr = pow((1-fec_ber),txbits);
        if(retx_dr > th_dr){
            chosen_mcs = mcs;
        }
    }

    Modulation chosen_mod = MCS2MODULATION[chosen_mcs];

    int total_payload = ceil(txbits*mcs2coderate[chosen_mcs]);

    double bits_per_symbol = bits_per_modulation(chosen_mod);
    int payload_with_mac_hdr = total_payload + MACsize;
    double numOfdmSymbols = ceil(payload_with_mac_hdr/(bits_per_symbol*Nsubs));
    double numAckOfdmSymbols = ceil(ACKsize/(bits_per_symbol*Nsubs));
    double tpayload = tOfdm*numOfdmSymbols;
    double tack = tpreamble + tOfdm*numAckOfdmSymbols;
    double tphyheader = tOfdm*ceil(PHYsize/(bits_per_symbol*Nsubs));

    double txtime;
    txtime = tpreamble + tpayload + tphyheader +tack;

return txtime;
}

double
digital_aware_ra::calculate_feedback_overhead(int ofdm_num, Modulation mod){

    double cbits_per_symbol, cber, dr, txtime, tput;
    double selected_txtime, maxtput = -100000;

    itpp::ivec feedback = FeedbackBits(mod);
    itpp::ivec compress = CompressBits(mod);

    int feedback_bits = feedback[ofdm_num];
    int compress_bits = compress[ofdm_num];

    std::cout<<"compress bits: "<<compress_bits<<std::endl;

    for(int mcs=6; mcs >=0; mcs--){

        cbits_per_symbol = CBITSPERMCS[mcs];
        cber = cber_curr_esno[mcs];
        dr = pow((1-cber),compress_bits);
        txtime = tOfdm*(compress_bits/(cbits_per_symbol*Nsubs));

        tput = compress_bits*dr/txtime;
        if(tput> maxtput){
            maxtput = tput;
            selected_txtime = txtime;
        }
    }
    std::cout<<"selected txtime: "<<selected_txtime<<std::endl;
return selected_txtime;
}

double
digital_aware_ra::calculate_transmissiontime(int pktlen, Modulation mod, int mcs, bool ack_flag){

    //double data_rate = 1.0e6*MCS2RATE[mcs];
    //double base_rate = 1.0e6*MCS2RATE[0];
    double bits_per_symbol = bits_per_modulation(mod);
    double cbits_per_symbol = CBITSPERMCS[mcs];
    double numOfdmSymbols = ceil(pktlen/(bits_per_symbol*Nsubs));
    double tpayload = tOfdm*numOfdmSymbols;

    double tphyheader = tOfdm*ceil( (PHYsize*2)/Nsubs);

    double numAckOfdmSymbols = ceil(ACKsize/(cbits_per_symbol*Nsubs));
    double tack = tpreamble + tphyheader + tOfdm*numAckOfdmSymbols;
    //double tmacheader = tOfdm*ceil(MACsize/(cbits_per_symbol*Nsubs));

    //double txtime = tpreamble + tSIFS + tDIFS + tpayload + (MACsize/data_rate) + (PHYsize/base_rate) + (ACKsize/data_rate);
    //double txtime = tpayload;
    //double txtime = tpreamble + tSIFS + tDIFS + tpayload;
    //double txtime = tpreamble + tSIFS + tDIFS + tpayload + (PHYsize/base_rate);
    //double txtime = tpreamble + tSIFS + tDIFS + tpayload + (MACsize/data_rate) + (PHYsize/base_rate);
    double txtime;
    txtime = tpreamble + tpayload + tphyheader + tack;
    //txtime = tpreamble + tpayload + tphyheader +tmacheader + tack;
    /*if (ack_flag==true){
        txtime = tpreamble + tpayload + tphyheader + tack;
    }
    else{
        txtime = tpreamble + tSIFS + tDIFS + tpayload + tphyheader;
    }*/

return txtime;
}

void 
digital_aware_ra::select_retransmission_rate(itpp::vec& bit_llrs, itpp::vec& sub_esno, int pktlen){

    // first assign llrs for transmission already happened
    int baselen = ((pktlen*8)+6)*2;
    itpp::vec full_llr_bitmap(baselen);
    itpp::Array<itpp::ivec> bitmap_list;
    itpp::Array<itpp::ivec> subs_list;
    itpp::ivec prev_mod;
    itpp::Array<itpp::ivec> prev_list;
    itpp::Array<itpp::Array<itpp::ivec> > prev_bitmap_list;
    itpp::Array<itpp::Array<itpp::ivec> > prev_subs_list;
    itpp::Array<itpp::vec> prev_llr_map;
    itpp::vec prev_txtime_list(1);
    itpp::vec prev_tput_list(1);
    debug_dr_list.set_size(1);

    prev_bitmap_list.set_size(1);
    prev_subs_list.set_size(1);
    prev_llr_map.set_size(1);
    bitmap_list.set_size(tx_count);
    subs_list.set_size(tx_count);
    full_llr_bitmap.zeros();

    double txtime = 0;
    int mcs_firsttx = selected_rate_list(0);

    // setting up data of prev transmissions
    for(int i =0; i<tx_count; i++){
        itpp::ivec bitmap = selected_bitmap(i);
        itpp::vec llrs = LlrsperTx(i);
        itpp::vec llr_bitmap;
        Modulation mod = (i == 0) ? MCS2MODULATION[selected_rate_list(i)]:int2modulation(selected_rate_list(i));
        prev_mod.ins(prev_mod.size(),selected_rate_list(i));
        assign_llr_per_bit_retx(llrs, bitmap, llr_bitmap, mod, baselen);
        //full_llr_bitmap.ins_row(i,llr_bitmap);
        full_llr_bitmap+=llr_bitmap;

        bool ack_flag = (i==0)?true:false;
        txtime += calculate_transmissiontime(bitmap.size(), mod, mcs_firsttx,ack_flag);
      
        bitmap_list(i) = bitmap;
        subs_list(i) = selected_subs(i);
        std::cout<<"subs list:"<<subs_list(i)<<std::endl;
        /* 
        std::cout<<"srr -> i = "<<i<<std::endl;
        std::cout<<"srr -> bitmap: "<<bitmap<<std::endl;
        std::cout<<"srr -> Modulation mod: "<<mod<<std::endl;
        std::cout<<"srr -> baselen: "<<baselen<<std::endl;
        std::cout<<"srr -> tx time: "<<txtime<<std::endl;
        std::cout<<"srr -> llr bitmap: "<<llr_bitmap<<std::endl;
        std::cout<<"srr -> llr: "<<llrs<<std::endl;
        std::cout<<"srr -> full llr bitmap: "<<full_llr_bitmap<<std::endl;
        */
    }

    //std::cout<<"tx count: "<<tx_count<<std::endl;
    //std::cout<<"curr pktid: "<<curr_pktid<<std::endl;
    Modulation first_mod = MCS2MODULATION[mcs_firsttx];
    //double dr = calculate_deliveryratio(full_llr_bitmap, first_mod, mcs_firsttx, pktlen);
    //double tput = dr*pktlen*8/txtime;
   
    prev_list.set_size(1);
    prev_list(0) = prev_mod;
    prev_llr_map(0) = full_llr_bitmap;   
    prev_bitmap_list(0) = bitmap_list;
    prev_subs_list(0) = subs_list;
    prev_txtime_list(0) = txtime;
    prev_tput_list(0) = 0;//tput;
    debug_dr_list(0) = 0;//tput;


    itpp::Array<itpp::ivec> curr_mod_list = prev_list;
    itpp::Array<itpp::Array<itpp::ivec> > curr_bitmap_list = prev_bitmap_list;
    itpp::Array<itpp::Array<itpp::ivec> > curr_subs_list = prev_subs_list;
    itpp::Array<itpp::vec> curr_llr_map = prev_llr_map;
    itpp::vec curr_txtime_list = prev_txtime_list;
    itpp::vec curr_tput_list = prev_tput_list;

    //std::cout<<"srr:curr txtime list ->"<<curr_txtime_list<<std::endl;

    for(int depth=tx_count; depth < retxmaxDepth; depth++){

        // all curr maintain the maxList vatlues for present depth
        // prev is for the previous depth
        for(int l=0; l < prev_list.size(); l++){

            itpp::ivec prev_mods = prev_list(l);
            int mcs1 = prev_mods[0]; // we want mcs of 1st transmission
            double prev_time = prev_txtime_list[l];
            double prev_tput = prev_tput_list[l];
            itpp::vec prev_llr = prev_llr_map(l);
            itpp::Array<itpp::ivec> prev_bitmap = prev_bitmap_list(l);
            itpp::Array<itpp::ivec> prev_subs = prev_subs_list(l);
            //std::cout<<"depth = "<<depth<<" l = "<<l<<std::endl;

            double t0 = current_time();
            second_transmission(curr_mod_list, curr_bitmap_list,curr_subs_list, curr_llr_map, curr_txtime_list,
                curr_tput_list, prev_mods, prev_llr, prev_bitmap, prev_subs, prev_time, prev_tput, sub_esno, mcs1, pktlen, l, depth, false);
            double t1 = current_time();
            //std::cout<<"second transmission time:"<<(t1-t0)<<std::endl;
            //fflush(stdout);
        }
        /* 
        std::cout<<"curr_mod_list: "<<curr_mod_list<<std::endl;
        std::cout<<"curr_tput_list: "<<curr_tput_list<<std::endl;
        std::cout<<"curr_txtime_list: "<<curr_txtime_list<<std::endl;
        for(int i=0; i<curr_llr_map.size(); i++){
            std::cout<<"curr_llr_map("<<i<<"): "<<curr_llr_map(i).get(0,depth-1,0,20)<<std::endl;
        }
        */
        // at the end of a depth curr becomes previous
        // llr_map and bitmap are cumulative for depth hence mats
        prev_list = curr_mod_list;
        prev_llr_map = curr_llr_map;
        prev_bitmap_list = curr_bitmap_list;
        prev_subs_list = curr_subs_list;
        prev_tput_list = curr_tput_list;
        prev_txtime_list = curr_txtime_list;
        //std::cout<<"prev bitmap "<<prev_bitmap_list<<std::endl;

    }
    
    std::cout<<"curr_mod_list: "<<curr_mod_list<<std::endl;
    std::cout<<"curr_tput_list: "<<curr_tput_list<<std::endl;
    std::cout<<"curr_txtime_list: "<<curr_txtime_list<<std::endl;
    std::cout<<"debug_dr_list: "<<debug_dr_list<<std::endl;
    for (int bl=0; bl < curr_bitmap_list.size();bl++){
        std::cout<<"bl "<<bl<<std::endl;
        itpp::Array<itpp::ivec> tmp_blist = curr_bitmap_list(bl);
        for(int int_l =0 ;int_l < tmp_blist.size(); int_l++){
            itpp::ivec tmp_int_l = tmp_blist(int_l);
            std::cout<<"bitmap size "<<tmp_int_l.size()<<std::endl;
        }
    }
    
    //std::cout<<"curr_bitmap_list: "<<curr_bitmap_list<<std::endl;

    int final_max_index = itpp::max_index(curr_tput_list);
    selected_rate_list =  curr_mod_list(final_max_index);
    selected_bitmap = curr_bitmap_list(final_max_index);
    selected_subs = curr_subs_list(final_max_index);
    selected_tput = curr_tput_list(final_max_index);
    // check for really low tput
    if (selected_tput < 1.0e3){
        std::cout<<"Poor tput fix called."<<std::endl;
        int sel_rate = selected_rate_list[0];
        selected_rate_list.set_size(1);
        selected_rate_list[0] = sel_rate;
        itpp::ivec new_bitmap, new_subs(Nsubs);
        //get_rate_bitmap(sel_rate,new_bitmap, baselen);
        get_rate_bitmap(sel_rate,new_bitmap, pktlen);
        selected_bitmap(tx_count-1) = new_bitmap;
        for(int s_tmp=0; s_tmp < Nsubs; s_tmp++){new_subs(s_tmp)=s_tmp;} 
        selected_subs(tx_count-1) = new_subs;
    }
    
    //std::cout<<"selected bitmap"<<selected_bitmap.size()<<std::endl;

}

void
digital_aware_ra::first_transmission(itpp::Array<itpp::ivec> &mcslist, itpp::vec &tput_list,
    itpp::Array<itpp::vec> &list_llr_map, itpp::Array< itpp::Array<itpp::ivec> > &bitmap_list, 
    itpp::Array< itpp::Array<itpp::ivec> > &subs_list,
    itpp::vec &txtime, itpp::vec &sub_esno, int pktlen){

    double time_ratebitmap=0;
    double time_assignllr=0;
    double time_list=0;
    double time_dr=0;
    double t0,t1;

    int maxMCS = 6;
    int minMCS = 0;

    int baselen = ((pktlen*8)+6)*2;

    itpp::ivec subs_firsttx(Nsubs);
    for(int i=0;i<Nsubs;i++){subs_firsttx(i)=i;}

    for(int mcs = maxMCS; mcs>=minMCS; mcs--){

        itpp::ivec bitmap;
        itpp::vec llr_bitmap;

        Modulation mod = MCS2MODULATION[mcs];

        //std::cout<<"baselen "<<baselen<<std::endl;
        t0 = current_time();
        //get_rate_bitmap(mcs, bitmap, baselen);
        get_rate_bitmap(mcs, bitmap, pktlen);
        t1 = current_time();
        time_ratebitmap+=(t1-t0);
        //std::cout<<"bitmap "<<bitmap.size()<<std::endl;
        t0 = current_time();
        assign_llr_per_bit(sub_esno, bitmap, llr_bitmap, mod, baselen);

        assign_bits_to_subs(bitmap, mcs);
        t1 = current_time();
        time_assignllr+=(t1-t0);
        //std::cout<<"llr bitmap: "<<llr_bitmap<<std::endl;
       

        t0 = current_time();
        double dr = calculate_deliveryratio(llr_bitmap,mod,mcs,pktlen);
        t1 = current_time();
        time_dr+=(t1-t0);

        double tx_time = calculate_transmissiontime(bitmap.size(), mod, mcs, true);
        double tput = pktlen*8*dr/tx_time;

        //std::cout<<"number bits: "<<bitmap.size()<<" mcs: "<<mcs<<std::endl;
        //std::cout<<"ft -> dr: "<<dr<<" tx_time: "<<tx_time<<" tput: "<<tput<<std::endl;
        t0 = current_time();
        if(mcslist.size()<maxList){

            mcslist.set_size(mcslist.size()+1,true);
            itpp::ivec tmpmcs(1);
            tmpmcs[0] = mcs;
            mcslist(mcslist.size()-1) = tmpmcs;

            txtime.ins(txtime.size(),tx_time);
            tput_list.ins(tput_list.size(),tput);
            debug_dr_list.ins(debug_dr_list.size(),dr);

            list_llr_map.set_size(mcslist.size(),true);
            list_llr_map(mcslist.size()-1) = llr_bitmap;

            bitmap_list.set_size(mcslist.size(),true);
            itpp::Array<itpp::ivec> bitmap_array(1);
            bitmap_array(0) = bitmap;
            bitmap_list(mcslist.size()-1) = bitmap_array;

            subs_list.set_size(mcslist.size(),true);
            itpp::Array<itpp::ivec> subs_array(1);
            subs_array(0) = subs_firsttx;
            subs_list(mcslist.size()-1) = subs_array;
 
        }
        else{
            // check if current tput is larger than previous values
            double min_tput = itpp::min(tput_list);
            if(min_tput < tput){
                int min_index = itpp::min_index(tput_list);
                itpp::ivec tmpmcs(1);
                tmpmcs[0] = mcs;
                mcslist(min_index) = tmpmcs;
                txtime[min_index] = tx_time;
                tput_list[min_index] = tput;
                debug_dr_list[min_index] = dr;

                list_llr_map(min_index) = llr_bitmap;

                itpp::Array<itpp::ivec> bitmap_array(1);
                bitmap_array(0) = bitmap;
                bitmap_list(min_index) = bitmap_array;

                itpp::Array<itpp::ivec> subs_array(1);
                subs_array(0) = subs_firsttx;
                subs_list(min_index) = subs_array;
            }

        }
        t1 = current_time();
        time_list+=(t1-t0);

        if(dr>0.99){
            break;
        }
    }
    /*
    std::cout<<" first transmission -->"<<std::endl; 
    std::cout<<"ft-> time rate bitmap: "<<time_ratebitmap<<std::endl;
    std::cout<<"ft-> time assign llr: "<<time_assignllr<<std::endl;
    std::cout<<"ft-> time list:"<<time_list<<std::endl;
    std::cout<<"ft-> time dr:"<<time_dr<<std::endl;
    */
}

void 
digital_aware_ra::remap_esno_values(itpp::vec &esno_values, itpp::vec &dest_esno, Modulation mod1, Modulation destmod){

    double scale_factor;

    switch(mod1){
        case BPSK:
            switch(destmod){
                case BPSK: scale_factor = 1.0; break;
                case QPSK: scale_factor = 0.5; break;
                case QAM16: scale_factor = 0.1; break;
                case QAM64: scale_factor = 0.02381; break;
                default: std::cout<<"Incorrect destmod selection."<<std::endl;
            }
            break;
        case QPSK:
            switch(destmod){
                case BPSK: scale_factor = 2.0; break;
                case QPSK: scale_factor = 1.0; break;
                case QAM16: scale_factor = 0.2; break;
                case QAM64: scale_factor = 0.04761; break;
                default: std::cout<<"Incorrect destmod selection."<<std::endl;
            }
            break;
        case QAM16:
            switch(destmod){
                case BPSK: scale_factor = 10.0; break;
                case QPSK: scale_factor = 5.0; break;
                case QAM16: scale_factor = 1.0; break;
                case QAM64: scale_factor = 0.2381; break;
                default: std::cout<<"Incorrect destmod selection."<<std::endl;
            }
            break;
        case QAM64:
            switch(destmod){
                case BPSK: scale_factor = 42.0; break;
                case QPSK: scale_factor = 21.0; break;
                case QAM16: scale_factor = 4.2; break;
                case QAM64: scale_factor = 1.0; break;
                default: std::cout<<"Incorrect destmod selection."<<std::endl;
            }
            break;
        default:
            std::cout<<"Incorrect first mod selection"<<std::endl;
    }

    dest_esno.set_size(esno_values.size());
    for(int i=0; i<esno_values.size(); i++){
        dest_esno[i] = scale_factor*esno_values[i];
    }
}

void 
digital_aware_ra::find_sorted_snr_groups(itpp::vec &esno_values,
    itpp::vec &unique_values, itpp::Array<itpp::ivec> &bits_per_group){

    itpp::vec val_count;
    //std::cout<<"fssg -> esno values: "<<esno_values<<std::endl;
    find_unique_values(esno_values, unique_values, val_count);
    //std::cout<<"fssg -> val count: "<<val_count<<std::endl;

    bits_per_group.set_size(unique_values.size());

    for(int i = 0; i <unique_values.size(); i++){

        itpp::ivec bitindex_vec;
        for(int j=0; j < esno_values.size(); j++){
            if (unique_values[i] == esno_values[j]){
                bitindex_vec.ins(bitindex_vec.size(),j);
            }
        }
        bits_per_group(i) = bitindex_vec;
    }

}


inline bool
check_for_max(itpp::vec &tput_list, int offset){

    int max_index = itpp::max_index(tput_list);
    bool max_reached;

    if (max_index+offset < tput_list.size()){
        max_reached = true;
    }
    else{
        max_reached = false;
    }
return max_reached;
}

void
digital_aware_ra::quantize_llr_values(itpp::vec& llr_values){

    // 1. convert to snr
    // 2. round to 1dB
    // 3. convert back to esno

    //itpp::vec snr = itpp::dB(llr_values);
    itpp::vec rounded_values = itpp::round(llr_values);
    llr_values = itpp::inv_dB(rounded_values);

}

double
digital_aware_ra::sort_subcarriers(itpp::vec& llr_values, itpp::ivec& sorted_sub_indices,int firsttx_mcs){

    itpp::vec sub_llr_values;
    itpp::ivec sub_bit_indices;
    itpp::Array<itpp::ivec> bits_per_sub_mod = Bits2Submap(firsttx_mcs);

    itpp::vec avg_llr_per_sub(Nsubs);

    double avg_ber, total_num_bits=0;
    for(int sub = 0; sub <Nsubs; sub++){
        sub_bit_indices =  bits_per_sub_mod(sub);   
        sub_llr_values = llr_values.get(sub_bit_indices);
        avg_ber = itpp::sum(1.0/(1+itpp::exp(sub_llr_values)))/sub_bit_indices.size();
        avg_llr_per_sub(sub) = log((1-avg_ber)/avg_ber);

        total_num_bits+= sub_bit_indices.size();
    }
    itpp::SORTING_METHOD method = itpp::INTROSORT;
    itpp::Sort<double> sort_llrs = Sort<double>(method);
    sorted_sub_indices = sort_llrs.sort_index(0,avg_llr_per_sub.size()-1,avg_llr_per_sub);

return total_num_bits;
}

void 
digital_aware_ra::second_transmission(itpp::Array<itpp::ivec> &curr_mod_list, itpp::Array<itpp::Array<itpp::ivec> > &curr_bitmap,
        itpp::Array<itpp::Array<itpp::ivec> > &curr_subs,
        itpp::Array<itpp::vec> &curr_llr_map, itpp::vec &curr_txtime_list,
        itpp::vec &curr_tput_list, const itpp::ivec &prev_mods, const itpp::vec &prev_llr_map, const itpp::Array<itpp::ivec>  &prev_bitmap, 
        const itpp::Array<itpp::ivec>  &prev_subs, double prev_txtime,
        double prev_tput, itpp::vec &sub_esno, int mcs1, int pktlen, int l, int depth, bool first_flag){


    double t0,t1;
    Modulation first_mod = MCS2MODULATION[mcs1];
    int dr_counter=0;

    double time_comb=0;
    t0 = current_time();
    itpp::vec comb_llr = prev_llr_map;
    itpp::vec comb_llr_abs = itpp::abs(comb_llr);
    int baselen = comb_llr.size();
    t1 = current_time();
    time_comb += (t1-t0);
 

    double time_snrgroups=0;
    t0 = current_time();

    //sort subcarriers based on average llr values
    itpp::ivec sorted_sub_indices;
    double total_num_bits = sort_subcarriers(comb_llr_abs, sorted_sub_indices, mcs1);


    //setup bit indices in ascending order for llr assignment
    itpp::ivec relevant_indices(total_num_bits);
    itpp::Array<itpp::ivec> relevant_submap = Bits2Submap(mcs1);

    int curr_index=0;
    itpp::ivec relevant_indices_per_sub;
    itpp::ivec bits_per_sub(Nsubs);
    for(int i=0; i<Nsubs;i++){
        relevant_indices_per_sub = relevant_submap(sorted_sub_indices(i));
        relevant_indices.set_subvector(curr_index,relevant_indices_per_sub);
        curr_index+=relevant_indices_per_sub.size();
        bits_per_sub(sorted_sub_indices(i)) = relevant_indices_per_sub.size();
    }

    std::cout<<"mcs 1: "<<mcs1<<std::endl;
    t1 = current_time();
    time_snrgroups+=(t1-t0);

    double time_assignllr =0;
    double time_dr = 0, time_drlist=0;
    double time_list = 0;
    double time_checkmax = 0;

    for(int i=0; i <4; i++){

        Modulation curr_mod = int2modulation(i);
        itpp::ivec bitmap;

        itpp::vec group_tput;
        int num_of_groups = Nsubs;

        itpp::vec llr_bitmap_complete;
        t0 = current_time();
        assign_llr_per_bit(sub_esno, relevant_indices, llr_bitmap_complete, curr_mod, baselen);
        t1 = current_time();
        time_assignllr+=(t1-t0);
        
        //std::cout<<"llr_bitmap_complete: "<<llr_bitmap_complete<<std::endl;
        t0 = current_time();
        itpp::vec uber_per_mod(comb_llr_abs.size());
        itpp::vec full_llr_bitmap = prev_llr_map;
        uber_per_mod.zeros();

        double llr_value;
        int uber_count=0;
        for(int cc=0; cc<comb_llr_abs.size(); cc++){
                llr_value = comb_llr_abs[cc];
                if(llr_value>0){
                    uber_per_mod[cc] = 1/(1+exp(llr_value));
                    uber_count++;
                } 
            }
        t1 = current_time();
        time_uber+=(t1-t0);

        int index1=0,index2=-1; 
        itpp::vec llr_bitmap(comb_llr_abs.size());
        llr_bitmap.zeros();
        double group_uber = 0;
        int curr_sub = 0;
        for(int group=0; group < Nsubs; group++){

            t0 = current_time();
            curr_sub = sorted_sub_indices(group);
            index1 = index2+1;
            index2 = index1+ bits_per_sub(curr_sub)-1;
            bitmap = relevant_indices.get(index1,index2);

            double llr_value, old_llr_value, comb_value;
            for(int v=0; v < bitmap.size(); v++){
                old_llr_value = comb_llr_abs[bitmap[v]];
                llr_value = llr_bitmap_complete[bitmap[v]];
                comb_value =old_llr_value + llr_value;
                llr_bitmap(bitmap[v]) = llr_value;
                full_llr_bitmap(bitmap[v]) = comb_value;
                uber_per_mod[bitmap[v]] = 1/(1+exp(comb_value));
            }

            group_uber = itpp::sum(uber_per_mod);
            double avgber = group_uber/uber_count;
            double cber = calculate_fecber(avgber, mcs1);
            double dr = pow((1 - cber),pktlen*8);

            if(curr_pktid==13){
                std::cout<<"max uber per mod: "<<max(uber_per_mod)<<std::endl;
                std::cout<<"avgber: "<<avgber<<" cber: "<<cber<<" dr: "<<dr<<std::endl;
            }

            t1 = current_time();
            time_dr+=(t1-t0);

            t0 = current_time();
            //itpp::mat full_llr_bitmap = prev_llr_map;
            //full_llr_bitmap.ins_row(full_llr_bitmap.rows(), llr_bitmap);
            //full_llr_bitmap.set_row(full_llr_bitmap.rows()-1,llr_bitmap);
            //std::cout<<"full llr bitmap: "<<full_llr_bitmap<<std::endl; 

            //double dr = calculate_deliveryratio(full_llr_bitmap, first_mod, mcs1, pktlen);
            //std::cout<<"group ber ->"<<group_uber<<" index2 "<<index2<<"avg ber ->"<<avgber<<std::endl;
            //std::cout<<"dr -> "<<dr<<" tmp dr -> "<<dr_tmp<<std::endl;
            dr_counter++;
            t1 = current_time();
            time_drlist+=(t1-t0);

            t0 = current_time();
            int bits_txed = index2;
            double tx_time = calculate_transmissiontime(bits_txed, curr_mod, mcs1, false);
            // FIXME:Feedback time needs to be calculated more accurately
            double feedback_time = 2*tOfdm;
            tx_time+=feedback_time;
            double tput = pktlen*8*dr/(prev_txtime+tx_time);
            group_tput.ins(group_tput.size(),tput);

            if(curr_mod_list.size()<maxList){

                curr_mod_list.set_size(curr_mod_list.size()+1,true);
                itpp::ivec tmp_mod = prev_mods;
                tmp_mod.ins(tmp_mod.size(),i);
                curr_mod_list(curr_mod_list.size()-1) = tmp_mod;

                curr_bitmap.set_size(curr_mod_list.size(),true);
                itpp::Array<itpp::ivec> bitmap_array = prev_bitmap;
                bitmap_array.set_size(bitmap_array.size()+1,true);
                bitmap_array(bitmap_array.size()-1) = relevant_indices.get(0,index2);
                curr_bitmap(curr_mod_list.size()-1) = bitmap_array;

                curr_subs.set_size(curr_mod_list.size(),true);
                itpp::Array<itpp::ivec> subs_array = prev_subs;
                subs_array.set_size(subs_array.size()+1,true);
                if(group==0){
                    itpp::ivec tmp_sub_vec(1);
                    tmp_sub_vec(0)= sorted_sub_indices(group);
                    subs_array(subs_array.size()-1) = tmp_sub_vec; 
                }
                else{
                    subs_array(subs_array.size()-1) = sorted_sub_indices.get(0,group);
                }
                curr_subs(curr_mod_list.size()-1) = subs_array;


                curr_tput_list.ins(curr_tput_list.size(),tput);
                curr_txtime_list.ins(curr_txtime_list.size(),prev_txtime+tx_time);
                debug_dr_list.ins(debug_dr_list.size(), dr);

                curr_llr_map.set_size(curr_mod_list.size(),true);
                curr_llr_map(curr_mod_list.size()-1) = full_llr_bitmap;
            }
            else{
                // check if current tput is larger than previous values
                double min_tput = itpp::min(curr_tput_list);
                if(min_tput < tput){
                    int min_index = itpp::min_index(curr_tput_list);

                    itpp::ivec tmp_mod = prev_mods;
                    tmp_mod.ins(tmp_mod.size(),i);
                    curr_mod_list(min_index) = tmp_mod;

                    itpp::Array<itpp::ivec> bitmap_array = prev_bitmap;
                    bitmap_array.set_size(bitmap_array.size()+1,true);
                    bitmap_array(bitmap_array.size()-1) = relevant_indices.get(0,index2);
                    curr_bitmap(min_index) = bitmap_array;

                    itpp::Array<itpp::ivec> subs_array = prev_subs;
                    subs_array.set_size(subs_array.size()+1,true);
                    subs_array(subs_array.size()-1) = sorted_sub_indices.get(0,group);
                    curr_subs(min_index) = subs_array;


                    curr_txtime_list[min_index] = prev_txtime+tx_time;
                    curr_tput_list[min_index] = tput;
                    debug_dr_list[min_index] = dr;
                    curr_llr_map(min_index) = full_llr_bitmap;
                }
            }
            t1 = current_time();
            time_list+=(t1-t0);
            
            //std::cout<<"tput_list: "<<curr_tput_list<<std::endl;
            //std::cout<<"txtime_list: "<<curr_txtime_list<<std::endl;
            //std::cout<<"crr mod_list: "<<curr_mod_list<<std::endl;
            
            t0 = current_time();
            bool max_reached = check_for_max(group_tput, 4);
            t1 = current_time();
            time_checkmax +=(t1-t0);

            if(max_reached){
                break;
            }

            if(dr > 0.99){
                break;
            }

        }
    }
    /* 
    std::cout<<"second_transmission timing info-->"<<std::endl;
    std::cout<<"time comb: "<<time_comb<<std::endl;
    std::cout<<"time check max: "<<time_checkmax<<std::endl;
    std::cout<<"time uber: "<<time_uber<<std::endl;
    std::cout<<"time llr: "<<time_assignllr<<std::endl;
    std::cout<<"time dr: "<<time_dr<<std::endl;
    std::cout<<"time dr list: "<<time_drlist<<std::endl;
    std::cout<<"time list: "<<time_list<<std::endl;
    std::cout<<"time snrgroups: "<<time_snrgroups<<std::endl;
    std::cout<<"dr counter: "<<dr_counter<<std::endl;
    fflush(stdout);
    */
    time_unique=0; time_uber=0; time_fec=0;
    dr_counter=0;
}


void 
digital_aware_ra::select_firsttransmission_rate(itpp::vec &sub_esno, int pktlen){

    itpp::Array<itpp::ivec> prev_list;
    itpp::vec prev_tput_list;
    itpp::Array<itpp::vec> prev_llr_map;
    itpp::Array< itpp::Array<itpp::ivec> > prev_bitmap_list;
    itpp::Array< itpp::Array<itpp::ivec> > prev_subs_list;
    itpp::vec prev_txtime_list;

    double t0 = current_time();
    first_transmission(prev_list, prev_tput_list, prev_llr_map, prev_bitmap_list, prev_subs_list,prev_txtime_list, sub_esno, pktlen);
    double t1 = current_time();

    itpp::Array<itpp::ivec> curr_mod_list = prev_list;
    itpp::Array<itpp::Array<itpp::ivec> > curr_bitmap_list = prev_bitmap_list;
    itpp::Array<itpp::Array<itpp::ivec> > curr_subs_list = prev_subs_list;
    itpp::Array<itpp::vec> curr_llr_map = prev_llr_map;
    itpp::vec curr_txtime_list = prev_txtime_list;
    itpp::vec curr_tput_list = prev_tput_list;


    for(int depth=1; depth < maxDepth; depth++){

        // all curr maintain the maxList vatlues for present depth
        // prev is for the previous depth
        for(int l=0; l < prev_list.size(); l++){

            itpp::ivec prev_mods = prev_list(l);
            int mcs1 = prev_mods[0]; // we want mcs of 1st transmission
            double prev_time = prev_txtime_list[l];
            double prev_tput = prev_tput_list[l];
            itpp::vec prev_llr = prev_llr_map(l);
            itpp::Array<itpp::ivec> prev_bitmap = prev_bitmap_list(l);
            itpp::Array<itpp::ivec> prev_subs = prev_subs_list(l);
            //std::cout<<"depth = "<<depth<<" l = "<<l<<" mcs: "<<mcs1<<std::endl;
            t0 = current_time();
            second_transmission(curr_mod_list, curr_bitmap_list, curr_subs_list, curr_llr_map, curr_txtime_list,
                curr_tput_list, prev_mods, prev_llr, prev_bitmap, prev_subs, prev_time, prev_tput, sub_esno, mcs1, pktlen, l, depth, true);
            t1 = current_time();
            //std::cout<<"first-tx second transmission time:"<<(t1-t0)<<std::endl;
        }
         
         
        // at the end of a depth curr becomes previous
        // llr_map and bitmap are cumulative for depth hence mats
        prev_list = curr_mod_list;
        prev_llr_map = curr_llr_map;
        prev_bitmap_list = curr_bitmap_list;
        prev_subs_list = curr_subs_list;
        prev_tput_list = curr_tput_list;
        prev_txtime_list = curr_txtime_list;
        //std::cout<<"prev bitmap "<<prev_bitmap_list<<std::endl;

    }
    
    std::cout<<"curr_mod_list: "<<curr_mod_list<<std::endl;
    std::cout<<"curr_tput_list: "<<curr_tput_list<<std::endl;
    std::cout<<"debug_dr_list: "<<debug_dr_list<<std::endl;
    
    int final_max_index = itpp::max_index(curr_tput_list);
    selected_rate_list =  curr_mod_list(final_max_index);
    selected_bitmap = curr_bitmap_list(final_max_index);
    selected_subs = curr_subs_list(final_max_index);
     
    for(int z=0; z<selected_bitmap.size(); z++){
        std::cout<<"len tx "<<z<<" -> "<<selected_bitmap(z).size()<<std::endl;
    }
    

}

void 
digital_aware_ra::firstpacket_retransmission(itpp::vec &llr_values, int firstmcs, int pktlen){

    LlrsperTx.set_size(1);
    LlrsperTx(0) = llr_values;

    selected_rate_list.set_size(1);
    selected_rate_list(0) = firstmcs;
    
    int baselen = ((pktlen*8)+6)*2;
    itpp::ivec bitmap;
    get_rate_bitmap(firstmcs, bitmap, pktlen);
    //get_rate_bitmap(firstmcs, bitmap, baselen);

    selected_bitmap.set_size(1);
    selected_bitmap(0) = bitmap;

    itpp::ivec subs(Nsubs);
    for(int i=0; i<Nsubs; i++){subs(i) = i;} 
    selected_subs.set_size(1);
    selected_subs(0) = subs;


    /*
    std::cout<<"EsnperTx: "<<EsNoperTx<<std::endl;
    std::cout<<"selected_rate_list: "<<selected_rate_list<<std::endl;
    std::cout<<"selected_bitmap: "<<selected_bitmap<<std::endl;
    */ 
}

void 
digital_aware_ra::calculate_cber_for_esno(itpp::vec& subcarrier_esno){

    for(int mcs = 0; mcs<7; mcs++){
        itpp::vec uber;
        Modulation mod = MCS2MODULATION[mcs];
        calculate_uber(subcarrier_esno, uber,mod);
        double avgber = itpp::sum(uber)/subcarrier_esno.size();
        double cber = calculate_fecber(avgber, mcs);
        cber_curr_esno[mcs] = cber;
    }
}

int
digital_aware_ra::select_rate(const std::vector<float> &recvbit_llrs, const std::vector<float> &subcarrier_esno, int pktlen, int pktid, int pkt_mcs, bool pkt_succ){


    itpp::vec bit_llrs(recvbit_llrs.size());
    int llr_sign;
    for(int i=0; i<recvbit_llrs.size(); i++){ 
        // 0.01 is added to differentiate it from punctured bits
        // llr sign is to make sure that the number does not accidently become zero
        llr_sign = (recvbit_llrs[i]>0)?1:-1;
        bit_llrs[i]=recvbit_llrs[i]+ (llr_sign*0.01); 
    }

    itpp::vec sub_esno(subcarrier_esno.size());
    for(int i=0; i<subcarrier_esno.size(); i++){ 
        if(subcarrier_esno[i]<0.5){
            sub_esno[i]=subcarrier_esno[i]; 
        }
        else{
            sub_esno[i] = roundf(subcarrier_esno[i]); 
            //sub_esno[i] = subcarrier_esno[i]; 
        }
    }

    calculate_cber_for_esno(sub_esno);

    //clear debug varaibles
    debug_dr_list.set_size(0);

    if(tx_count > 0){

        LlrsperTx.set_size(tx_count,true);
        LlrsperTx(LlrsperTx.size()-1) = bit_llrs;
    }
    
    if (pkt_succ == true){
        retxmaxDepth = maxDepth;
        tx_count = 0;
        LlrsperTx.set_size(0);
        // if pkt succ is true bit esno is subcarrier esno
        assert(sub_esno.size()==Nsubs);
        double t0 = current_time();
        select_firsttransmission_rate(sub_esno, pktlen);
        double t1 = current_time();
        std::cout<<"select first transmission time: "<<(t1-t0)<<std::endl;
    } 
    else{
        if(curr_pktid!=pktid){
            tx_count=1;
            retxmaxDepth = maxDepth;
            firstpacket_retransmission(bit_llrs, pkt_mcs, pktlen);
        }
        retxmaxDepth += 1;
        double t0 = current_time();
        select_retransmission_rate(bit_llrs, sub_esno,pktlen);
        double t1 = current_time();
        std::cout<<"select re-transmission time: "<<(t1-t0)<<std::endl;
    }

    if (tx_count == 0){
        selected_mcs = selected_rate_list(tx_count);
    }
    else{
        selected_mcs = RETXMOD2MCS[selected_rate_list(tx_count)];
    }
    tx_count += 1;
    curr_pktid = pktid;


return selected_mcs;
}

