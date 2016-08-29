/* -*- c++ -*- */
/*
 * Copyright 2007,2008,2010,2011 Free Software Foundation, Inc.
 * 
 * This file is part of GNU Radio
 * 
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <digital_ofdm_partialcomb_receiver.h>
#include <gr_io_signature.h>
#include <gr_expj.h>
#include <gr_math.h>
#include <math.h>
#include <cstdio>
#include <stdexcept>
#include <iostream>
#include <string.h>
#include <fstream>
#include <algorithm>

#include "digital_ofdm_refsymbols.h"
#include "digital_ofdm_evm_tables.h"
#include "digital_ofdm_coded_data.h"

#define VERBOSE 0
#define REFSNR 1
#define LTFPREAMBLE 1
#define REFCODEDDATA 1
#define REFORLTF 2 // REF=0, LTF=1, K-DIST=2
#define OFFLINE 0


inline void
conv_bvec_to_chardata(itpp::bvec& vecdata, int valid_bits,unsigned char* data_bytes){

    //printf("valid bits: %d",valid_bits);
    int boffset;
    int ioffset;
    for( int i=0; i<valid_bits; i++){
        boffset = i/8;
        ioffset = i%8;
        int val = (int)vecdata[i];
        //printf("boffset=%d val=%d byte=%x\n",boffset,val,data_bytes[boffset]);
        data_bytes[boffset] |=(val&0x1)<<(8-ioffset-1);
    }

}

inline void
digital_ofdm_partialcomb_receiver::enter_search()
{
  if (VERBOSE)
    fprintf(stderr, "@ enter_search\n");

  d_state = STATE_SYNC_SEARCH;

  #if LTFPREAMBLE == 1
    preamble_count = 0;
  #endif
}


/*
*/
int
digital_ofdm_partialcomb_receiver::get_codeddata_length(int pktlen, int mcs){

    double const1[] = {2.0,2.0,4.0,2.0,4.0,3.0,4.0,6.0};
    double const2[] = {1.0,1.0,3.0,1.0,3.0,2.0,3.0,5.0};

    int bitmap_size = ceil(((pktlen*8)+6)*(const1[mcs]/const2[mcs]));
    return bitmap_size;

}

/*
*/

double 
digital_ofdm_partialcomb_receiver::calculate_effsnr_from_subsnr(std::vector<double>& esno, digital_ofdm_effsnr_ra& effra, int mcs){

    Modulation effmod = MCS2MODULATION[mcs];
    std::vector<double> refuber;
    effra.calculate_uber(refuber,esno,effmod);
    double eff_avguber = 0;
    for(size_t i=0; i<refuber.size();i++){eff_avguber+=refuber[i];}
    eff_avguber/=refuber.size();
    double ref_effesno = effra.calculate_effesno(eff_avguber,effmod);
    double ref_effsnr = 10.0*log10(ref_effesno);

    return ref_effsnr;
}

#if REFSNR == 1
void 
digital_ofdm_partialcomb_receiver::calculate_refsnr(std::vector<gr_complex>& symbols, std::vector<double>& snr,int mcs){

    int symb_size=0;
    std::vector<std::complex<double> >* refsymbs;
    if(mcs == 0){ symb_size = symbsMCS0; refsymbs = &refsymbsMCS0; }
    else if(mcs==1){ symb_size = symbsMCS1; refsymbs = &refsymbsMCS1; }
    else if(mcs==2){ symb_size = symbsMCS2; refsymbs = &refsymbsMCS2; }
    else if(mcs==3){ symb_size = symbsMCS3; refsymbs = &refsymbsMCS3; }
    else if(mcs==4){ symb_size = symbsMCS4; refsymbs = &refsymbsMCS4; }
    else if(mcs==5){ symb_size = symbsMCS5; refsymbs = &refsymbsMCS5; }
    else if(mcs==6){ symb_size = symbsMCS6; refsymbs = &refsymbsMCS6; }
    else if(mcs==7){ symb_size = symbsMCS7; refsymbs = &refsymbsMCS7; }
    
    if(refsymbs->size()!=symbols.size()){
        printf("Symbol sizes not equal. Ref=%d, Rcvd=%d\n",symb_size,symbols.size());
    }
    std::vector<double> dist_values(Nsubs,0.0);
    std::vector<int> count(Nsubs,0);
    for(size_t i=0; i<symbols.size(); i++){
        //std::complex<double> tmp = std::complex<double>(symbols[i]);
        double tmpreal = (double) symbols[i].real();
        double tmpimag = (double) symbols[i].imag();
        std::complex<double> tmp = std::complex<double>(tmpreal,tmpimag);
        //printf("tmp=(%f,%f)\n",tmp.real(),tmp.imag());
        //printf("refsymbs=(%f,%f)\n",(*refsymbs)[i].real(),(*refsymbs)[i].imag());
        double diff = std::abs(tmp-(*refsymbs)[i]);
        double value = (diff*diff);
        int sub = i%Nsubs;
        dist_values[sub]+= value;
        count[sub]+=1;
    }
    for(int i=0; i<Nsubs; i++){ 
        double evm = dist_values[i]/count[i];
        double esno_val = 1.0/evm;
        double snr_value = 10.0*log10(esno_val);
        outfile_refsnr<<snr_value<<",";
        snr.push_back(snr_value);
    }
    outfile_refsnr<<"\n";
}
/*
*/
void 
digital_ofdm_partialcomb_receiver::calculate_retx_refsnr(std::vector<gr_complex>& symbols, std::vector<gr_complex>& ref_symbols, std::vector<double>& snr,int mcs){

    if(ref_symbols.size()!=symbols.size()){
        printf("Symbol sizes not equal. Ref=%d, Rcvd=%d\n",ref_symbols.size(),symbols.size());
    }

    std::vector<double> dist_values(Nsubs,0.0);
    std::vector<int> count(Nsubs,0);
    for(size_t i=0; i<symbols.size(); i++){
        double value = std::norm(ref_symbols[i]-symbols[i]);

        int sub = i%Nsubs;
        dist_values[sub]+= value;
        count[sub]+=1;
    }
    for(int i=0; i<Nsubs; i++){ 
        double evm = dist_values[i]/count[i];
        double esno_val = 1.0/evm;
        double snr_value = 10.0*log10(esno_val);
        outfile_refsnr<<snr_value<<",";
        snr.push_back(snr_value);
    }
    outfile_refsnr<<"\n";
}

void 
digital_ofdm_partialcomb_receiver::create_ref_symbols(int first_mcs, int retx_mcs, std::vector<int>& subs, std::vector<gr_complex>& symbols){

    itpp::bvec output_data;
    std::vector<int>& txRefCodedData = TxCodedData[first_mcs];
    // assignbit2submap() ensures all subs have same bits
    // printf("subs size=%d, submapsize = %d\n",subs.size(),d_bit2submap[0].size());
    int valid_subs = 0;
    for(size_t i=0; i<subs.size(); i++){ if(subs[i]==1){valid_subs++;} }

    int output_size = valid_subs*d_bit2submap[0].size();
    output_data.set_size(output_size);
    // printf("output size=%d\n",output_size);

    // printf("txRefCodedDataSize %d\n",txRefCodedData.size());

    int bit_count = 0;
    for(size_t i=0; i<subs.size();i++){
        if(subs[i]==1){
            std::vector<int>* bitmapPtr = &d_bit2submap[i];
            for(size_t j=0; j<bitmapPtr->size(); j++){
                if(bitmapPtr->at(j)>=txRefCodedData.size()){
                printf("%d,",bitmapPtr->at(j));
                }
                output_data[bit_count] = txRefCodedData[bitmapPtr->at(j)];
                bit_count++;
            }
        }
    }
    //printf("\n");

    unsigned int bits_per_mod = MCS2BITSMOD[retx_mcs];
    int bits_per_ofdm = bits_per_mod*d_data_carriers.size();
    int missing_bits = bit_count%bits_per_ofdm;
    if(missing_bits > 0){
        unsigned int bit_pad = bits_per_ofdm - missing_bits;
        int new_size = bit_count+bit_pad;
        output_data.set_size(new_size,true);
        for (unsigned int i=0;i<bit_pad;i++){
            //output_data[bit_count+i] = (rand()%2);
            output_data[bit_count+i] = fixedrand[i];
        }
    }
    //for(size_t i=0; i<output_data.size();i++){printf("%d,",(int)output_data[i]);}
    //printf("\n");

    std::vector<gr_complex> constellation;
    std::vector< std::vector<int> > constel_values;
    get_constellations(bits_per_mod, constellation,constel_values);

    int num_symbols = output_data.size()/bits_per_mod;
    printf("2nd: output_datasize=%d, num_symbols=%d\n",output_data.size(),num_symbols);
    int bit_offset = 0;
    for(int i=0; i<num_symbols; i++){

        unsigned int bits = itpp::bin2dec(output_data(bit_offset,bit_offset+bits_per_mod-1));
        bit_offset += bits_per_mod;
        symbols.push_back(constellation[bits]);
        //printf("(%.2f,%.2f),",symbols[i].real(),symbols[i].imag());
    }
    //printf("\n");
}
#endif


/*
*/
inline double 
map_evm_to_snr(double evm, Modulation mod){

    double* evm_table;
    if(mod==BPSK){ evm_table = &evm_values_array_bpsk[0]; }
    else if(mod==QPSK){ evm_table = &evm_values_array_qpsk[0];}
    else if(mod==QAM16){ evm_table = &evm_values_array_qam16[0];}
    else if(mod==QAM64){ evm_table = &evm_values_array_qam64[0];}

    int low = 0;
    int high = num_evm-1;
    int mid = 0;
    while(1){
        mid = (low+high)/2;

        if((low==mid) || (high==mid)){ break;}

        if(evm_table[mid]==evm){ break;}
        else if(evm_table[mid]<evm){ high = mid;}
        else{ low = mid;}
    }

    double snr = -10 + (0.2*mid);
    return snr;
}

void 
digital_ofdm_partialcomb_receiver::calculate_adjdistsnr(std::vector<gr_complex>& symbols, std::vector<double>& esno, std::vector<double>& snr,int mcs){
    Modulation mod = MCS2MODULATION[mcs];
    int bits_per_mod = BITSPERMOD[mod];

    std::vector<gr_complex> constellation;
    std::vector< std::vector<int> > constel_values;
    get_constellations(bits_per_mod, constellation,constel_values);
    unsigned int table_size = constellation.size();

    std::vector<double> dist_values(Nsubs,0.0);
    std::vector<int> count(Nsubs,0);


    for(size_t i=0; i<symbols.size(); i++){

        gr_complex x = symbols[i];
        float min_euclid_dist = norm(x - constellation[0]);
        float euclid_dist = 0;
  
        for (unsigned int j = 1; j < table_size; j++){
            euclid_dist = norm(x - constellation[j]);
            if (euclid_dist < min_euclid_dist){
                min_euclid_dist = euclid_dist;
            }
        }
        //double value = 1.0/min_euclid_dist;
        int sub = i%Nsubs;
        //esno_values[sub]+= value;
        dist_values[sub]+= min_euclid_dist;
        count[sub]+=1;
    }
    for(int i=0; i<Nsubs; i++){ 
        double evm_value = std::sqrt(dist_values[i]/count[i]);
        double snr_value = map_evm_to_snr(evm_value,mod);
        double esno_value = std::pow(10,snr_value/10.0);
        esno.push_back(esno_value);
        snr.push_back(snr_value);
    }

}
/*
*/
void 
digital_ofdm_partialcomb_receiver::calculate_distsnr(std::vector<gr_complex>& symbols, std::vector<double>& esno, std::vector<double>& snr,int mcs){

    Modulation mod = MCS2MODULATION[mcs];
    int bits_per_mod = BITSPERMOD[mod];

    std::vector<gr_complex> constellation;
    std::vector< std::vector<int> > constel_values;
    get_constellations(bits_per_mod, constellation,constel_values);
    unsigned int table_size = constellation.size();

    std::vector<double> dist_values(Nsubs,0.0);
    std::vector<int> count(Nsubs,0);

    for(size_t i=0; i<symbols.size(); i++){

        //unsigned int min_index = 0;
        gr_complex x = symbols[i];
        float min_euclid_dist = norm(x - constellation[0]);
        float euclid_dist = 0;
  
        for (unsigned int j = 1; j < table_size; j++){
            euclid_dist = norm(x - constellation[j]);
            if (euclid_dist < min_euclid_dist){
                min_euclid_dist = euclid_dist;
            }
        }
        //double value = 1.0/min_euclid_dist;
        int sub = i%Nsubs;
        //esno_values[sub]+= value;
        dist_values[sub]+= min_euclid_dist;
        count[sub]+=1;
    }
    for(int i=0; i<Nsubs; i++){ 
        double evm_value = dist_values[i]/count[i];
        double esno_value = 1.0/evm_value;
        double snr_value = 10.0*log10(esno_value);
        esno.push_back(esno_value);
        snr.push_back(snr_value);
    }
}

/*
*/
void 
digital_ofdm_partialcomb_receiver::calculate_kdistsnr(std::vector<gr_complex>& symbols, std::vector<double>& esno, std::vector<double>& snr, std::vector<double>& mean_snr, std::vector<double>& std_snr,int mcs){

    Modulation mod = MCS2MODULATION[mcs];
    int bits_per_mod = BITSPERMOD[mod];

    std::vector<gr_complex> constellation;
    std::vector< std::vector<int> > constel_values;
    get_constellations(bits_per_mod, constellation,constel_values);
    unsigned int table_size = constellation.size();

    std::vector<std::vector<double> > all_dist_values(Nsubs);
    std::vector<int> count(Nsubs,0);

    for(int i=0; i<Nsubs; i++){
        all_dist_values[i] = std::vector<double>(table_size,0.0);
    }
        

    for(size_t i=0; i<symbols.size(); i++){

        std::vector<double> symb_dist;

        //unsigned int min_index = 0;
        gr_complex x = symbols[i];
  
        for (unsigned int j = 0; j < table_size; j++){
            symb_dist.push_back(norm(x - constellation[j]));
        }

        std::sort(symb_dist.begin(),symb_dist.end());

        int sub = i%Nsubs;
        for(unsigned int k=0; k<table_size; k++){
            all_dist_values[sub][k]+= symb_dist[k];
        }
        count[sub]+=1;
    }


    std::vector<double> snr_value(table_size,0.0);
    std::vector<double> prob_occ(table_size,0.0);

    for(int sub=0; sub< Nsubs; sub++){
        double prob_sum=0;
        for(unsigned int k=0; k<table_size; k++){
            double evm_value = all_dist_values[sub][k]/count[sub];
            double esno_value = 1.0/evm_value;
            snr_value[k] = 10.0*log10(esno_value);
            prob_occ[k] = std::exp(-1.0*(abs(snr_value[k] - mean_snr[sub]))/std_snr[sub]);
            prob_sum+=prob_occ[k];
        }
    
        double prob_sanity=0;
        double final_sub_snr=0;
        for(unsigned int i=0; i<table_size; i++){
            prob_occ[i]=prob_occ[i]/prob_sum;
            final_sub_snr+= (prob_occ[i]*snr_value[i]);
            prob_sanity+=prob_occ[i];
        }
        snr.push_back(final_sub_snr);
        double subesno = std::pow(10,final_sub_snr/10.0);
        esno.push_back(subesno);

    }
}


/*
*/
void 
digital_ofdm_partialcomb_receiver::calculate_llrsnr(itpp::vec& llr_values, std::vector<double>& snr, int mcs){

    Modulation mod = MCS2MODULATION[mcs];
    int bits_per_mod = BITSPERMOD[mod];

    std::vector<double> uber(Nsubs,0.0);
    std::vector<int> count(Nsubs,0);
    for(int i=0; i<llr_values.size(); i++){
        double uber_value = 1.0/(1.0 + std::exp(std::abs(llr_values[i])));
        int down_idx = i%(bits_per_mod*Nsubs);
        int sub = floor(down_idx/bits_per_mod);
        uber[sub]+=uber_value;
        count[sub]+=1;
        //printf("(%d %f %f), ",sub,llr_values[i],uber_value);
    }
    //printf("\n");

    digital_ofdm_rate_adaptation ra;
    for(int i=0; i<Nsubs; i++){
        double mean_uber = uber[i]/count[i];
        //printf("(%f, %f, %d) ",uber[i],mean_uber,count[i]);
        double snr_value = 10.0*log10(ra.calculate_effesno(mean_uber,mod));
        if (std::isinf(snr_value)){snr_value=30.0;}
        snr.push_back(snr_value);
    }
    //printf("\n");
}


#if LTFPREAMBLE == 1
void 
digital_ofdm_partialcomb_receiver::calculate_ltfpreamble_distsnr(std::vector<double>& EsNo, std::vector<double>& Snr){

    std::vector<double> dist_val(Nsubs,0.0);
    double ltf_value;
    unsigned int offset = 23;
    for(size_t i=0; i<eqd_ltfpreambles.size(); i++){
        ltf_value = -2*fixedrand[i+offset] + 1; 
        dist_val[(i%Nsubs)] +=std::norm(eqd_ltfpreambles[i] - ltf_value);
    }

    for(int i=0; i<Nsubs; i++){
        double avg_dist = dist_val[i]/2.0;
        double esno_val = 1/avg_dist; 
        double snr_val = 10.0*log10(esno_val);
        EsNo.push_back(esno_val);
        Snr.push_back(snr_val); 
    }

}

#if 0
void 
digital_ofdm_partialcomb_receiver::calculate_ltfpreamble_distsnr(std::vector<double>& EsNo, std::vector<double>& Snr){

    std::vector<double> dist_val(Nsubs,0.0);

    
    for(int i=0; i<eqd_ltfpreambles.size(); i++){
        dist_val[(i%Nsubs)] += norm(eqd_ltfpreambles[i] - (double)ltf_preamble_rx[(i%Nsubs)]);
    }

    for(int i=0; i<Nsubs; i++){
        double avg_dist = dist_val[i]/2.0;
        double esno_val = 1/avg_dist; 
        double snr_val = 10.0*log10(esno_val);
        EsNo.push_back(esno_val);
        Snr.push_back(snr_val); 
    }

}
#endif

void
digital_ofdm_partialcomb_receiver::calculate_ltfpreamble_snr(){

    std::complex<double> h1,h2;
    std::complex<double>* h;
    h = new std::complex<double>[d_data_carriers.size()];
    double noise=0;

    for(size_t i=0; i < d_data_carriers.size(); i++){
        //estimate channel
        h1 = rcvd_ltfpreambles[0][i]*(double)ltf_preamble_rx[i];
        h2 = rcvd_ltfpreambles[1][i]*(double)ltf_preamble_rx[i];
        h[i] = (h1+h2)/2.0;
        noise += std::pow(std::abs(rcvd_ltfpreambles[0][i] - rcvd_ltfpreambles[1][i]),2);
    } 
    noise/= (2*d_data_carriers.size());

    for(size_t j=0; j < d_data_carriers.size(); j++){
        //ltf_subcarrier_snr[j] = 10.0*log10(std::pow(std::abs(h[j]),2)/noise); 
        ltf_subcarrier_esno[j] = std::pow(std::abs(h[j]),2)/noise; 
        outfile_ltfsnr<<ltf_subcarrier_esno[j]<<",";
    }
    outfile_ltfsnr<<"\n";

    free(h);
    std::cout<<"noise: "<<noise<<std::endl;
    std::vector<double> tmpsnr;
    for(size_t i=0; i<d_data_carriers.size();i++){
        double tmp = 10.0*log10(ltf_subcarrier_esno[i]);
        tmpsnr.push_back(tmp);
    }
//    std::cout<<"ltf snr: "<<std::endl;
//    for(int i=0; i<tmpsnr.size();i++){printf("%f, ",tmpsnr[i]);}
//    printf("\n");
}
#endif
    
inline void
create_constellation_type(int bits_per_mod, double real[], double imag[], std::vector<gr_complex>& constel, std::vector<std::vector<int> >& values){

    int number = pow(2,bits_per_mod);
    for(int i=0; i<number; i++){
        constel.push_back(gr_complex(real[i],imag[i]));
        std::vector<int> bits_per_val;
        for(int b=0; b<bits_per_mod; b++){
            int bit_val = (i&(1<<(bits_per_mod-b-1)))>>(bits_per_mod-b-1);
            bits_per_val.push_back(bit_val);
        }
        values.push_back(bits_per_val);
    }
}

void 
digital_ofdm_partialcomb_receiver::create_constellations(){

    int bits_per_mod;
    // bpsk
    bpsk_const.push_back(gr_complex(-1.0,0.0));
    bpsk_const.push_back(gr_complex(1.0,0.0));
    bpsk_values.push_back(std::vector<int>(1,0)); // fill vector size=1 with 0
    bpsk_values.push_back(std::vector<int>(1,1));

    // qpsk
    double real_qpsk[] = {-0.707, -0.707,  0.707, 0.707};
    double imag_qpsk[] = {-0.707,  0.707, -0.707, 0.707};

    bits_per_mod = 2;
    create_constellation_type(bits_per_mod,real_qpsk,imag_qpsk, qpsk_const, qpsk_values);

/*    for(int i=0; i<qpsk_const.size();i++){
        printf("(%f,%f) ",qpsk_const[i].real(),qpsk_const[i].imag());
        for(int j=0; j<qpsk_values[i].size(); j++){printf("%d ",(int)(qpsk_values[i])[j]);}
    printf("\n");
    }
*/
   // qam16
    bits_per_mod = 4; 
    double real_qam16[] = { -0.9487,-0.9487,-0.9487,-0.9487,-0.3162,-0.3162,-0.3162,-0.3162, 0.9487, 0.9487, 0.9487, 0.9487, 0.3162, 0.3162, 0.3162, 0.3162};
    double imag_qam16[] = { -0.9487,-0.3162, 0.9487, 0.3162,-0.9487,-0.3162, 0.9487, 0.3162,-0.9487,-0.3162, 0.9487, 0.3162,-0.9487,-0.3162, 0.9487, 0.3162};
    create_constellation_type(bits_per_mod,real_qam16,imag_qam16, qam16_const, qam16_values);


    // qam64
    bits_per_mod = 6;
    double real_qam64[] = {-1.0801,-1.0801,-1.0801,-1.0801,-1.0801,-1.0801,-1.0801,-1.0801,-0.7715,-0.7715,-0.7715,-0.7715,-0.7715,-0.7715,-0.7715,-0.7715,-0.1543,-0.1543,-0.1543,-0.1543,-0.1543,-0.1543,-0.1543,-0.1543,-0.4629,-0.4629,-0.4629,-0.4629,-0.4629,-0.4629,-0.4629,-0.4629, 1.0801, 1.0801, 1.0801, 1.0801, 1.0801, 1.0801, 1.0801, 1.0801, 0.7715, 0.7715, 0.7715, 0.7715, 0.7715, 0.7715, 0.7715, 0.7715, 0.1543, 0.1543, 0.1543, 0.1543, 0.1543, 0.1543, 0.1543, 0.1543, 0.4629, 0.4629, 0.4629, 0.4629, 0.4629, 0.4629, 0.4629, 0.4629};

    double imag_qam64[] = {-1.0801,-0.7715,-0.1543,-0.4629, 1.0801, 0.7715, 0.1543, 0.4629,-1.0801,-0.7715,-0.1543,-0.4629, 1.0801, 0.7715, 0.1543, 0.4629,-1.0801,-0.7715,-0.1543,-0.4629, 1.0801, 0.7715, 0.1543, 0.4629,-1.0801,-0.7715,-0.1543,-0.4629, 1.0801, 0.7715, 0.1543, 0.4629,-1.0801,-0.7715,-0.1543,-0.4629, 1.0801, 0.7715, 0.1543, 0.4629,-1.0801,-0.7715,-0.1543,-0.4629, 1.0801, 0.7715, 0.1543, 0.4629,-1.0801,-0.7715,-0.1543,-0.4629, 1.0801, 0.7715, 0.1543, 0.4629,-1.0801,-0.7715,-0.1543,-0.4629, 1.0801, 0.7715, 0.1543, 0.4629};

    create_constellation_type(bits_per_mod,real_qam64,imag_qam64, qam64_const,qam64_values);
}

//
void
digital_ofdm_partialcomb_receiver::send_feedback(int mcs, int succ, int pktid, std::vector<int> & retxbitmap){


    d_refid = rand()%128;

    FeedbackStruct info;
    info.pktid = pktid;
    info.mcs = mcs;
    info.succ = succ;
    info.refid = d_refid;
    for(int i=0; i<Nsubs; i++){ info.retx[i] = retxbitmap[i]; }

    printf("pktid=%d, mcs=%d, succ=%d, refid=%d\n",pktid, mcs,succ,d_refid);
    for(int i=0; i<Nsubs;i++){ printf("%d, ",info.retx[i]); } printf("\n");
    socklen_t fromlen = sizeof(struct sockaddr_in);
    int n = sendto(sock,&info,sizeof(FeedbackStruct),0,(struct sockaddr *)&from,fromlen);

#if 0
    unsigned char buffer[100];
    int size = retxbitmap.size() + 3;
    buffer[0] = (unsigned char) mcs;
    buffer[1] = (unsigned char) succ;
    buffer[2] = d_refid;
    
    for(size_t i=0; i<retxbitmap.size();i++){

        buffer[i+3] = retxbitmap[i];
    }
    printf("pktid=%d, mcs=%d, succ=%d, refid=%d\n",pktid, mcs,succ,d_refid);
    socklen_t fromlen = sizeof(struct sockaddr_in);
    int n = sendto(sock,&buffer,size,0,(struct sockaddr *)&from,fromlen);
#endif
}

//
void 
digital_ofdm_partialcomb_receiver::open_feedback_socket(){

    int portno = 51717;
    char buf[1024];

    sock=socket(AF_INET, SOCK_DGRAM, 0);
    if (sock < 0){ fprintf(stderr,"Opening socket"); exit(0);}
    int length = sizeof(server);
    bzero(&server,length);
    server.sin_family=AF_INET;
    server.sin_addr.s_addr=INADDR_ANY;
    server.sin_port=htons(portno);
    if (bind(sock,(struct sockaddr *)&server,length)<0){
        fprintf(stderr,"Error binding\n");
        exit(0);
    }
    socklen_t fromlen = sizeof(struct sockaddr_in);
    int n = recvfrom(sock,buf,1024,0,(struct sockaddr *)&from,&fromlen);
    if (n < 0){ fprintf(stderr,"Error recvfrom\n");exit(0);}
    else{ printf("data received.");}

    n = sendto(sock,"Got your message\n",18,0,(struct sockaddr *)&from,fromlen);

}

/*
*/
void
digital_ofdm_partialcomb_receiver::assign_bits_to_subs(std::vector<std::vector<int> >& bit2submap, int bitmap_size, int mcs, bool inlvFlag){

    std::vector<int> inlvBitmap;
    std::vector<int> bitmap;

    Modulation mod = MCS2MODULATION[mcs];
    int Nbpsc = BITSPERMOD[mod];
    int Ncbps = Nbpsc*Nsubs;
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

    bit2submap = std::vector<std::vector<int> >(Nsubs);

    for (size_t i=0; i<inlvBitmap.size();i++){
        sub_value = floor(cntr/Nbpsc); // divison by bits-per-sub
        if (sub_value >= Nsubs){
            sub_value = 0;
            cntr = 0;
        }
        bit2submap[sub_value].push_back(inlvBitmap[i]);
        cntr++;
    }

}   // end of assign bit to subs



//
void 
digital_ofdm_partialcomb_receiver::combine_data(itpp::vec& recv_data, std::vector<int>& comb_subs, itpp::vec& output){

    std::vector<int> comb_idxs;

    int count = 0;
    for(size_t i=0; i<comb_subs.size();i++){
        std::vector<int>* bitmapPtr = &d_bit2submap[i];
        if (comb_subs[i]==1){
            for(size_t j=0; j<bitmapPtr->size(); j++){
                //printf("%d,",bitmapPtr->at(j));
                output[bitmapPtr->at(j)] += recv_data[count];
                count++;
            }
        }
    }
    //printf("\n");
}


/*
*/
inline void
digital_ofdm_partialcomb_receiver::enter_have_sync(){

    if (VERBOSE)
        fprintf(stderr, "@ enter_have_sync\n");

    #if LTFPREAMBLE == 1
        d_state = STATE_HAVE_SYNC;
    #else
        d_state = STATE_HAVE_LTFPREAMBLE;
    #endif


    d_state = STATE_HAVE_SYNC;

    // Resetting PLL
    d_freq = 0.0;
    d_phase = 0.0;
    fill(d_dfe.begin(), d_dfe.end(), gr_complex(1.0,0.0));
}

//
inline 
unsigned int dec_convhelper(itpp::bvec& hdr, int offset, int len){

    unsigned int value=0;
    for(int i =0; i<len; i++){
        int tmp =  (int) hdr[offset+i];
        value |= (((unsigned int)tmp)&0x1)<<(len-i-1);
    }
    return value;
}

bool 
digital_ofdm_partialcomb_receiver::enter_have_header(){
    
    itpp::bvec hdr_decode;

    hdr_demod = -2.0*(hdr_demod) + 1.0;
    fec_decoder.decode(0,hdr_demod,hdr_decode);

    header.pktlength = (unsigned short) dec_convhelper(hdr_decode,0,12);
    header.mcs = (unsigned short) dec_convhelper(hdr_decode,12,4);
    header.pktid = (unsigned short) dec_convhelper(hdr_decode,16,16);
    header.whitener = (unsigned short) dec_convhelper(hdr_decode,32,4);
    header.retx_count = (unsigned short) dec_convhelper(hdr_decode,36,4);
    header.payloadlen = (unsigned short) dec_convhelper(hdr_decode,40,16);
    header.verifyid = (unsigned short) dec_convhelper(hdr_decode,56,8);
    header.hdr_crc = (unsigned int) dec_convhelper(hdr_decode,64,32);

    printf("pktlen=%u, mcs=%u, pktid=%u, whitener=%u,retx=%u, payloadlen=%u, verifyid=%u, crc=%u\n",
header.pktlength,header.mcs,header.pktid,header.whitener,header.retx_count,header.payloadlen, header.verifyid, header.hdr_crc);

    unsigned char* header_data = (unsigned char*) &header;
    int len_minus_crc = sizeof(PHY_HEADER)-4;
    unsigned int hdr_crc = digital_crc32(header_data,len_minus_crc);
    printf("header crc=%u\n",hdr_crc);

    if(header.hdr_crc == hdr_crc){
        return true;
    }
    else{
        return false;
    } 
}

//
void 
digital_ofdm_partialcomb_receiver::approxllrslicer(const gr_complex x, double noise, int bits_per_mod, std::vector<std::vector<gr_complex> >& S0_table,std::vector<std::vector<gr_complex> >& S1_table, std::vector<double>& output){

    unsigned int table_size = S0_table[0].size();

    //printf("noise: %f\n",noise);
    //printf("table_size: %d, bits_per_mod: %d\n",table_size,bits_per_mod);
    double dist0, dist1;
    for (int i=0; i < bits_per_mod; i++){
        double min_dist0 = norm(x - S0_table[i][0]);
        double min_dist1 = norm(x - S1_table[i][0]);

        for (unsigned int j=1; j < table_size; j++){
            dist0 = norm(x - S0_table[i][j]);
            dist1 = norm(x - S1_table[i][j]);
            if(dist0 < min_dist0){ min_dist0 = dist0; }
            if(dist1 < min_dist1){ min_dist1 = dist1; }
        }

        output[i] = -1*(min_dist0 - min_dist1)/noise;
    }
}

//
void 
digital_ofdm_partialcomb_receiver::llrslicer(const gr_complex x, double noise, int bits_per_mod, std::vector<std::vector<gr_complex> >& S0_table,std::vector<std::vector<gr_complex> >& S1_table, std::vector<double>& output){

    unsigned int table_size = S0_table[0].size();

    //printf("noise: %f\n",noise);
    //printf("table_size: %d, bits_per_mod: %d\n",table_size,bits_per_mod);
    for (int i=0; i < bits_per_mod; i++){
        double sumS0=0,sumS1=0,distS0, distS1;
        for (unsigned int j=0; j < table_size; j++){
            distS0 = norm(x - S0_table[i][j])/noise;
            distS1 = norm(x - S1_table[i][j])/noise;
            //printf("dist: %f, %f\n",distS0,distS1);
            sumS0 += exp(-1*distS0);
            sumS1 += exp(-1*distS1);
        }
        double tmp = log(sumS0/sumS1);
        //printf("%f, ",tmp);
        tmp = std::min(tmp,1000.0);
        output[i] = std::max(tmp,-1000.0);
        //output[i] = log(sumS0/sumS1);
    }
//printf("\n");
}


//
void 
digital_ofdm_partialcomb_receiver::bitslicer(const gr_complex x,std::vector<gr_complex>& constellation, std::vector<std::vector<int> >& constel_value, std::vector<int>& output)
{
    unsigned int table_size = constellation.size();
    unsigned int min_index = 0;
    float min_euclid_dist = norm(x - constellation[0]);
    float euclid_dist = 0;
  
    for (unsigned int j = 1; j < table_size; j++){
        euclid_dist = norm(x - constellation[j]);
        if (euclid_dist < min_euclid_dist){
            min_euclid_dist = euclid_dist;
            min_index = j;
        }
    }
    output = constel_value[min_index];
  
}


/* uses pilots - from rawofdm -- optionally returns if needed */
inline void
digital_ofdm_partialcomb_receiver::equalize_interpolate_dfe(gr_complex *in, gr_complex *factor) 
{
  gr_complex phase_error = 0.0;

  float cur_pilot = 1.0;
  for (unsigned int i = 0; i < d_pilot_carriers.size(); i++) {
    gr_complex pilot_sym(cur_pilot, 0.0);
    cur_pilot = -cur_pilot;
    int di = d_pilot_carriers[i];
    phase_error += (in[di] * conj(pilot_sym));
    /*if(d_packetid==6){
        printf("in: (%f, %f), pilot: (%f, %f) in(rad): (%f, %f)\n", in[di].real(), in[di].imag(), pilot_sym.real(), pilot_sym.imag(),abs(in[di]),arg(in[di]));
        fflush(stdout);
    }*/
  } 

  // update phase equalizer
  float angle = arg(phase_error);

  /*if(d_packetid==6){
     std::cout<<"phase error angle: "<<angle<<std::endl;
  }*/
  d_freq = d_freq - d_freq_gain*angle;
  d_phase = d_phase + d_freq - d_phase_gain*angle;
  if (d_phase >= 2*M_PI) d_phase -= 2*M_PI;
  else if (d_phase <0) d_phase += 2*M_PI;

  gr_complex carrier = gr_expj(-angle);

  // update DFE based on pilots
  cur_pilot = 1.0;
  for (unsigned int i = 0; i < d_pilot_carriers.size(); i++) {
    gr_complex pilot_sym(cur_pilot, 0.0);
    cur_pilot = -cur_pilot;
    int di = d_pilot_carriers[i];
    gr_complex sigeq = in[di] * carrier * d_dfe[i];
    // FIX THE FOLLOWING STATEMENT
    if (norm(sigeq)> 0.001)
      d_dfe[i] += d_eq_gain * (pilot_sym/sigeq - d_dfe[i]);
  }

  // equalize all data using interpolated dfe and demap into bytes
  unsigned int pilot_index = 0;
  int pilot_carrier_lo = 0;
  int pilot_carrier_hi = d_pilot_carriers[0];
  gr_complex pilot_dfe_lo = d_dfe[0];
  gr_complex pilot_dfe_hi = d_dfe[0];

  /* alternate implementation of lerp */
  float denom = 1.0/(pilot_carrier_hi - pilot_carrier_lo);					// 1.0/(x_b - x_a)

  for (unsigned int i = 0; i < d_data_carriers.size(); i++) {
      int di = d_data_carriers[i];
      if(di > pilot_carrier_hi) {					// 'di' is beyond the current pilot_hi
	  pilot_index++;						// move to the next pilot

	  pilot_carrier_lo = pilot_carrier_hi;				// the new pilot_lo
	  pilot_dfe_lo = pilot_dfe_hi;					// the new pilot_dfe_lo

	 if (pilot_index < d_pilot_carriers.size()) {
	     pilot_carrier_hi = d_pilot_carriers[pilot_index];
	     pilot_dfe_hi = d_dfe[pilot_index];
	 }
	 else {
	       pilot_carrier_hi = d_occupied_carriers;			// cater to the di's which are beyond the last pilot_carrier
           //pilot_carrier_hi = d_data_carriers.size()+d_pilot_carriers.size()-1;
	 }
	 denom = 1.0/(pilot_carrier_hi - pilot_carrier_lo);
      }

      float alpha = float(pilot_carrier_hi - di) * denom;
      //float alpha = float(di - pilot_carrier_lo) * denom;		// (x - x_a)/(x_b - x_a)
      gr_complex dfe = pilot_dfe_hi + alpha * (pilot_dfe_lo - pilot_dfe_hi);
      //gr_complex dfe = pilot_dfe_lo + alpha * (pilot_dfe_hi - pilot_dfe_lo);	// y = y_a + (y_b - y_a) * (x - x_a)/(x_b - x_a) 
      in[di] *= carrier * dfe;
      /*if(factor != NULL) {
         factor[i] = gr_expj(angle);
	 if(norm(dfe) > 0.0f) 
	    factor[i] /= dfe;
      }*/
  }
}

void
digital_ofdm_partialcomb_receiver::get_constellations(int bits_per_mod, std::vector<gr_complex>& constel, std::vector<std::vector<int> >& constel_values){
    if (bits_per_mod == 1){
        constel = bpsk_const;
        constel_values = bpsk_values;
    }
    else if (bits_per_mod == 2){
        constel = qpsk_const;
        constel_values = qpsk_values;
    }
    else if (bits_per_mod == 4){
        constel = qam16_const;
        constel_values = qam16_values;
    }
    else if (bits_per_mod == 6){
        constel = qam64_const;
        constel_values = qam64_values;
    }
    else{
        printf("Incorrect bits_per_mod value.\n");
    }
}

//
void 
digital_ofdm_partialcomb_receiver::set_currllr_table(int bits_per_mod){

    if (bits_per_mod == 1){
        currS0_table = S0_table_bpsk; currS1_table = S1_table_bpsk; 
        //printf("BPSK set.");
    }
    else if (bits_per_mod == 2){
        currS0_table = S0_table_qpsk; currS1_table = S1_table_qpsk; 
        //printf("QPSK set.");
    }
    else if (bits_per_mod == 4){
        currS0_table = S0_table_qam16; currS1_table = S1_table_qam16; 
        //printf("QAM16 set.");
    }
    else if (bits_per_mod == 6){
        currS0_table = S0_table_qam64; currS1_table = S1_table_qam64; 
        //printf("QAM64 set.");
    }
}

/*
*/
void
digital_ofdm_partialcomb_receiver::soft_combine_data(std::vector<gr_complex>& symbols, std::vector<double>& esno, int retx_count){

    if (retx_count == 0){
        soft_rcvd_symbols.clear();
        soft_combined_data.clear();
        cum_esno.clear();
        cum_esno_sqrt.clear();
        int idx = 0;
        for (size_t k=0; k<symbols.size(); k+=1){
            idx = k%esno.size();
            float sqrt_esno = std::sqrt(esno[idx]);
            gr_complex scaled_symb = sqrt_esno*symbols[k];
            soft_rcvd_symbols.push_back(scaled_symb);
            soft_combined_data.push_back(symbols[k]);
        }
        for (size_t k=0; k<esno.size(); k+=1){ 
            cum_esno_sqrt.push_back(std::sqrt(esno[k])); 
            cum_esno.push_back(esno[k]); 
        }
    }
    else{
        for (size_t k=0; k<esno.size(); k+=1){
            cum_esno_sqrt[k] += std::sqrt(esno[k]);
            cum_esno[k] = cum_esno_sqrt[k]*cum_esno_sqrt[k];
        }

        int idx = 0;
        for (size_t k=0; k<symbols.size(); k+=1){
            idx = k%(esno.size());
            float sqrt_esno = std::sqrt(esno[idx]);
            soft_rcvd_symbols[k] += (sqrt_esno*symbols[k]);
            gr_complex tmp = gr_complex(soft_rcvd_symbols[k]);
            soft_combined_data[k] = gr_complex(tmp.real()/cum_esno_sqrt[idx],tmp.imag()/cum_esno_sqrt[idx]);
        }
    }
}

/*
*/
void 
digital_ofdm_partialcomb_receiver::llr_demapper_alldata(std::vector<gr_complex>& symbols, std::vector<double>& ref_esno, itpp::vec& out, unsigned int bits_per_mod){
    
    //int num_symbols = symbols.size()/Nsubs;
    std::vector<double> llr_values(bits_per_mod);
    int offset_value = 0;
    for(size_t i=0; i<symbols.size(); i++){

        int sub = i%Nsubs;
        double noise = 1/ref_esno[sub];
        llrslicer(symbols[i], noise, bits_per_mod, currS0_table, currS1_table, llr_values);
        for(unsigned int b=0; b<bits_per_mod; b++){
            out[offset_value] = llr_values[b];
            offset_value++;
        }


    }

}

//
void 
digital_ofdm_partialcomb_receiver::llr_demapper(gr_complex *in, itpp::vec& out, unsigned int& offset_value, unsigned int bits_per_mod){

    
    equalize_interpolate_dfe(in, NULL);
    unsigned int i=0;
    std::vector<double> llr_values(bits_per_mod);

    while( i < d_data_carriers.size() ) {

        int di = d_data_carriers[i];
        double noise = 1/ltf_subcarrier_esno[i];
        llrslicer(in[di], noise, bits_per_mod, currS0_table, currS1_table, llr_values);
        //approxllrslicer(in[di], noise, bits_per_mod, currS0_table, currS1_table, llr_values);

        for(unsigned int b=0; b<bits_per_mod; b++){
            out[offset_value] = llr_values[b];
            offset_value++;
        }
        i++;
    }

}


//
void 
digital_ofdm_partialcomb_receiver::bit_demapper(gr_complex *in, itpp::vec& out, unsigned int& offset_value, unsigned int bits_per_mod){

    equalize_interpolate_dfe(in, NULL);
    unsigned int i=0;

    std::vector<int> bit_values;
    std::vector<gr_complex> constel;
    std::vector< std::vector<int> > constel_values;
    get_constellations(bits_per_mod, constel,constel_values);

    while( i < d_data_carriers.size() ) {

      int di = d_data_carriers[i];
      bitslicer(in[di],constel,constel_values,bit_values);


      for(unsigned int b=0; b<bits_per_mod; b++){
        out[offset_value] = bit_values[b];
        offset_value++;
      }
      i++;
    }


}

/*
*/
digital_ofdm_partialcomb_receiver_sptr
digital_make_ofdm_partialcomb_receiver(const std::vector<gr_complex> &sym_position, 
			     const std::vector<unsigned char> &sym_value_out,
			     gr_msg_queue_sptr target_queue, unsigned int occupied_carriers,
                 unsigned int ra_mode,
                 unsigned int tx_gain,
                 unsigned int bw_intf,
                 unsigned int run,
			     float phase_gain, float freq_gain)
{
  return gnuradio::get_initial_sptr(new digital_ofdm_partialcomb_receiver(sym_position, sym_value_out,
								target_queue, occupied_carriers, ra_mode, tx_gain, bw_intf,run,
								phase_gain, freq_gain));
}

digital_ofdm_partialcomb_receiver::digital_ofdm_partialcomb_receiver(const std::vector<gr_complex> &sym_position, 
						 const std::vector<unsigned char> &sym_value_out,
						 gr_msg_queue_sptr target_queue, unsigned int occupied_carriers,
                         unsigned int ra_mode, unsigned int tx_gain, unsigned int bw_intf,unsigned int run,
						 float phase_gain, float freq_gain)
  : gr_sync_block ("ofdm_frame_sink",
		   gr_make_io_signature2 (2, 2, sizeof(gr_complex)*occupied_carriers, sizeof(char)),
		   gr_make_io_signature (1, 1, sizeof(gr_complex)*occupied_carriers)),
    d_target_queue(target_queue), d_occupied_carriers(occupied_carriers), 
    d_phase(0),d_freq(0),d_phase_gain(phase_gain),d_freq_gain(freq_gain),
    d_eq_gain(0.05),
    total_pkts(0), succ_pkts(0), total_time(0.0),succrcvd_data(0),
    firsttx_mcs(0),
    llrtime_elapsed(0.0),
    fec_decoder(0),
    ewma_alpha(0.125),std_beta(0.25),
    d_lastpktsucc_id(-1){
  
  if (RAMODE > 2){
    std::cout<<"Only ra modes=0,1,2 supported.\n";
    assert(0);
  }
  else{ RAMODE = ra_mode;}

  printf("RAMODE = %d\n",RAMODE);

  assign_subcarriers();

  Nsubs = d_data_carriers.size();

  if(d_data_carriers.size() > d_occupied_carriers) {
    throw std::invalid_argument("digital_ofdm_mapper_bcv: subcarriers allocated exceeds size of occupied carriers");
  }


  #if LTFPREAMBLE == 1
//    int ltf_preamble96[] = {1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1, -1, -1,1,1, -1,1, -1,1, -1, -1, -1, -1, -1,1,1, -1, -1,1, -1,1,-1,0,0,0,0,0,0,0,0,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1, -1, -1,1,1, -1,1, -1,1, -1, -1, -1, -1, -1,1,1, -1, -1,1, -1,1,-1};

    int ltf_preamble_values[] ={1,-1,-1,1,1,-1,1,-1,1,-1,1,-1,1,-1,-1,1,1,-1,1,-1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,1,1,-1,1,-1,-1,1,-1,1,1,-1,-1,1,1,-1,1,-1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1};

    for(size_t ltf=0; ltf < d_data_carriers.size(); ltf++){
        ltf_preamble_rx.push_back(ltf_preamble_values[ltf]);
        printf("%d,",ltf_preamble_rx[ltf]);
    }
    printf("\n");

    int d_ltfpreamble_size = d_data_carriers.size();
    rcvd_ltfpreambles = new std::complex<double>*[2];
    rcvd_ltfpreambles[0] = new std::complex<double>[d_ltfpreamble_size];
    rcvd_ltfpreambles[1] = new std::complex<double>[d_ltfpreamble_size];

    ltf_subcarrier_esno = std::vector<double>(d_data_carriers.size());
    outfile_ltfsnr.open("ltf_snr_bpsk.txt");
    outfile_refsnr.open("SubSnrFiles/ref_snr.txt");
    outfile_effsnr.open("EffSnrFiles/ref_snr.txt");
    dist_subsnr.open("SubSnrFiles/dist_snr.txt");
    kdist_subsnr.open("SubSnrFiles/kdist_snr.txt");
    ltf_subsnr.open("SubSnrFiles/ltf_snr.txt");
    h_estimates = std::vector<double>(d_data_carriers.size(),0);
  #endif

#if TRACEFILE == 1
    std::string ra_folder;

    if(RAMODE == 0){ ra_folder = "EffSNR/";}
    else if(RAMODE == 1){ ra_folder = "Smart/";}
    else if(RAMODE == 2) { ra_folder = "Soft/";}
    std::string trace_file = "TraceFolder/CameraReady/";
    trace_file+=ra_folder;
    trace_file+=std::string("trace_bwintf");
    std::stringstream outbwintf; outbwintf<<bw_intf;
    trace_file+=std::string(outbwintf.str());
    trace_file+=std::string("_txgain");
    std::stringstream outtxgain; outtxgain<<tx_gain;
    trace_file+=std::string(outtxgain.str());
    trace_file+=std::string("_run");
    //trace_file+=std::string("trace_tmp_run");
    std::stringstream outrun; outrun<<run;
    trace_file+=std::string(outrun.str());
    trace_file+=std::string(".txt");
    traceOutFile.open(trace_file.c_str());
#endif

  d_dfe.resize(occupied_carriers);
  fill(d_dfe.begin(), d_dfe.end(), gr_complex(1.0,0.0));

  // hdr variables
  int num_ofdm = ceil(1.0*(sizeof(PHY_HEADER)*8*2 + 12)/d_data_carriers.size());
  coded_hdr_bitlen = num_ofdm*d_data_carriers.size();
  printf("coded_hdr_bitlen=%d\n",coded_hdr_bitlen);
  fflush(stdout);
  d_hdr_offset = 0;
  hdr_demod.set_size(coded_hdr_bitlen);

  create_constellations();
  create_llr_tables();

#if REFSNR == 1
  for(int i=0; i<8; i++){
    //create_refvector(i);
    create_refvector1000Bytes(i);
  }
#endif

  ltfsnr_count = 0;
  Sk_ltfsnr  = std::vector<double>(Nsubs,0.0);
  avg_ltfsnr = std::vector<double>(Nsubs,0.0);
  std_ltfsnr = std::vector<double>(Nsubs,0.0);

#if REFCODEDDATA == 1
    for(int i=0; i<8; i++){
        create_refcodeddata(i,1000);
    }
#endif

  // create deinterleaver
  interleaver = digital_ofdm_interleaver (d_data_carriers.size());
  deinterleaver = digital_ofdm_deinterleaver (d_data_carriers.size());

  // rate adaptation
  int maxDepth = 2;
  int maxList = 4;
  int maxRetx = 7;
  if( (RAMODE == 0) || (RAMODE == 2) ){
      ra_effsnr = digital_ofdm_effsnr_ra (d_data_carriers.size());
  }
  else{
      ra_smart = digital_ofdm_smart_ra (d_data_carriers.size(),maxDepth,maxList,maxRetx);
  }

 
#if OFFLINE == 0
  open_feedback_socket();
#endif
  enter_search();
}

//
digital_ofdm_partialcomb_receiver::~digital_ofdm_partialcomb_receiver(){

  double dratio = 1.0*succ_pkts/total_pkts;
  double tput = 1.0e-3*succrcvd_data/total_time;
  printf("total pkts=%d, succ=%d, dratio=%f, tput=%f kbps\n",total_pkts,succ_pkts,dratio,tput);
  fflush(stdout);

  #if LTFPREAMBLE == 1
    free(rcvd_ltfpreambles[0]);
    free(rcvd_ltfpreambles[1]);
    free(rcvd_ltfpreambles);
    outfile_ltfsnr.close();
    outfile_refsnr.close();
  #endif

#if TRACEFILE == 1
    for(size_t i=0; i<traceMcs.size(); i++){

        traceOutFile<<"{ ";
        traceOutFile<<" 'mcs': "<<traceMcs[i]<<", ";
        traceOutFile<<" 'payload': "<<tracePayloadLen[i]<<", ";
        traceOutFile<<" 'succpkt': "<<traceSuccPkt[i]<<", ";
        traceOutFile<<" 'pktid': "<<tracePktid[i]<<", ";
        traceOutFile<<" 'biterrors': "<<traceBitErrors[i]<<", ";
        traceOutFile<<" 'totalbits': "<<traceTotalBits[i]<<", ";
        traceOutFile<<" 'retx': [ ";
        for(size_t k=0; k<Nsubs; k++){
            traceOutFile<<traceRetxSubs[i][k]<<", ";
        }
        traceOutFile<<"], ";
        traceOutFile<<" 'biterrsub': [ ";
        for(size_t k=0; k<Nsubs; k++){
            traceOutFile<<traceErrsPerSub[i][k]<<", ";
        }
        traceOutFile<<"] }\n";

    }
    traceOutFile<<"{ 'Dratio': "<<dratio<<", 'TotalTime': "<<total_time<<", 'TotalPkts': "<<total_pkts<<", 'Tput': "<<tput<<" }\n";
    
    traceOutFile.close();

#endif

#if 0
    outfile_effsnr<<"RefSNR: ";
    for(int i=0; i<refEffsnrVec.size(); i++){ outfile_effsnr<<refEffsnrVec[i]<<","; }
    outfile_effsnr<<std::endl;

    outfile_effsnr<<"DistSNR: ";
    for(int i=0; i<distEffsnrVec.size(); i++){ outfile_effsnr<<distEffsnrVec[i]<<","; }
    outfile_effsnr<<std::endl;

    outfile_effsnr<<"DistAdjSNR: ";
    for(int i=0; i<distadjEffsnrVec.size(); i++){ outfile_effsnr<<distadjEffsnrVec[i]<<","; }
    outfile_effsnr<<std::endl;

    outfile_effsnr<<"LtfSNR: ";
    for(int i=0; i<ltfEffsnrVec.size(); i++){ outfile_effsnr<<ltfEffsnrVec[i]<<","; }
    outfile_effsnr<<std::endl;

    outfile_effsnr<<"LtfDistSNR: ";
    for(int i=0; i<ltfdistEffsnrVec.size(); i++){ outfile_effsnr<<ltfdistEffsnrVec[i]<<","; }
    outfile_effsnr<<std::endl;
#endif
    outfile_effsnr.close();
    dist_subsnr.close();
    kdist_subsnr.close();
    ltf_subsnr.close();
}

/*
*/
double 
digital_ofdm_partialcomb_receiver::calculate_transmission_time(int mcs, int payload, int txcount, int mcs2){

    int MACsize = 28*8;
    int ACKsize = 14*8;
    int PHYsize = 22;
    int num_subs = Nsubs;
#if WIFITIME == 1
    double tSlottime = 9.0e-6;
    double tSIFS = 16.0e-6;
    double tpreamble = 24.0e-6;
    double tDIFS = tSIFS + 2*tSlottime;
    double tOfdm = 4.0e-6;
#else
    double tOfdm = 128.0e-6;
    double tpreamble = 6*tOfdm;
#endif
    double cbits_per_symbol = CBITSPERMCS[mcs];
    double numOfdmSymbols = 0.0;
    if( (RAMODE == 0) || (RAMODE == 2) ){
        int bits_per_symbol = BITSPERMOD[ MCS2MODULATION[mcs] ];
        numOfdmSymbols = ceil(payload/(bits_per_symbol*num_subs));
    }
    else{
        if (txcount == 0){
            int bits_per_symbol = BITSPERMOD[ MCS2MODULATION[mcs] ];
            numOfdmSymbols = ceil(payload/(bits_per_symbol*num_subs));
        }
        else{
            Modulation mod2 = MCS2MODULATION[mcs2];
            int bits_per_symbol = BITSPERMOD[mod2];
            numOfdmSymbols = ceil(payload/(bits_per_symbol*num_subs));
        }
    }
    double tpayload = tOfdm*numOfdmSymbols;

    double base_bits_per_symbol = 0.5;
    double tphyheader = tOfdm*ceil(PHYsize/(base_bits_per_symbol*num_subs));
    double tmacheader = tOfdm*ceil(MACsize/(cbits_per_symbol*num_subs));

    double numAckOfdmSymbols = ceil(ACKsize/(cbits_per_symbol*num_subs));
    double tack;
    if( (RAMODE == 0) || (RAMODE == 2) ){
            tack = tpreamble + tphyheader + tOfdm*numAckOfdmSymbols;
    }
    else{
        if(txcount == 0){
            tack = tpreamble + tphyheader + tOfdm*numAckOfdmSymbols;
        }
        else{ tack = tOfdm; }
    }
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
digital_ofdm_partialcomb_receiver::assign_subcarriers() {
   int dc_tones = 8;
   //int num_pilots = 8;
   int pilot_gap = 11;
 
   int half1_end = (d_occupied_carriers-dc_tones)/2;     //40
   int half2_start = half1_end+dc_tones;                 //48
 
   // first half
   for(int i = 0; i < half1_end; i++) {
      if(i%pilot_gap == 0)
         d_pilot_carriers.push_back(i);
      else
         d_data_carriers.push_back(i);
   }
 
   // second half
   for(unsigned int i = half2_start, j = 0; i < d_occupied_carriers; i++, j++) {
      if(j%pilot_gap == 0)
         d_pilot_carriers.push_back(i);
      else
         d_data_carriers.push_back(i);
   }
 
   /* debug carriers */
   printf("pilot carriers: \n");
   for(size_t i = 0; i < d_pilot_carriers.size(); i++) {
      printf("%d ", d_pilot_carriers[i]); fflush(stdout);
   }
   printf("\n");
 
   printf("data carriers: \n");
   for(size_t i = 0; i < d_data_carriers.size(); i++) {
      printf("%d ", d_data_carriers[i]); fflush(stdout);
   }
   printf("\n");
 }


/*
*/
int
digital_ofdm_partialcomb_receiver::work (int noutput_items,
			       gr_vector_const_void_star &input_items,
			       gr_vector_void_star &output_items)
{
    gr_complex *in = (gr_complex *) input_items[0];
    const char *sig = (const char *) input_items[1];
    // If the output is connected, send it the derotated symbols

    if(output_items.size() >= 1){
        d_derotated_output = (gr_complex *)output_items[0];
    }
    else{
        d_derotated_output = NULL;
    }
  
    if (VERBOSE){ fprintf(stderr,">>> Entering state machine\n");}


    switch(d_state) {
      
        case STATE_SYNC_SEARCH:    // Look for flag indicating beginning of pkt
            
            if (VERBOSE)
              fprintf(stderr,"SYNC Search, noutput=%d\n", noutput_items);
   
            // Found it, set up for header decode
            if (sig[0]) { 
                std::vector<gr_tag_t> tags;
                const uint64_t nread = this->nitems_read(0);
                const size_t ninput_items = noutput_items;
                this->get_tags_in_range(tags, 0, nread, nread+ninput_items);

                pmt::pmt_t pmt_vec = tags[0].value;
                size_t vec_size;
                double *elements = pmt_f64vector_writable_elements(pmt_vec, vec_size); 

                for(size_t i=0; i<d_data_carriers.size();i++){
                    h_estimates[i] = elements[d_data_carriers[i]];
                }
                //for(int i=0; i<vec_size; i++){printf("%f, ",elements[i]);}
                //printf("\n");

                //for(int i=0; i<tags.size();i++){ printf("%f",tags[i]);}
                enter_have_sync();
            }
            break;

        case STATE_HAVE_SYNC:
            #if LTFPREAMBLE==1
                equalize_interpolate_dfe(in, NULL);
                if(preamble_count == 0){ eqd_ltfpreambles.clear(); }
                for (size_t pl = 0; pl < d_data_carriers.size(); pl++){
                    eqd_ltfpreambles.push_back(in[d_data_carriers[pl]]);
                    rcvd_ltfpreambles[preamble_count][pl] = in[d_data_carriers[pl]];
                    rcvd_ltfpreambles[preamble_count][pl]*=h_estimates[pl];
                }
                preamble_count += 1;
                if(preamble_count == 2){
                    calculate_ltfpreamble_snr();
                    d_state = STATE_HAVE_LTFPREAMBLE;
                    preamble_count = 0;
                }
                break;
            #endif

        case STATE_HAVE_LTFPREAMBLE:
            // only demod after getting the preamble signal; otherwise, the 
            // equalizer taps will screw with the PLL performance
            bit_demapper(&in[0], hdr_demod, d_hdr_offset,1);


            if (d_hdr_offset==coded_hdr_bitlen) { //start if header-read
                d_hdr_offset = 0;
                d_data_offset= 0;
                bool hdr_valid = enter_have_header();

                if(hdr_valid){

                    if(d_lastpktsucc_id >= header.pktid){
                        printf("Dropping packets because already received successfully.\n");
                        enter_search();
                        break;
                    }

                    d_state = STATE_HAVE_HEADER;

                    data_demod.set_size(header.payloadlen);
                    //if(header.retx_count == 0){
                        //data_comb.set_size(header.payloadlen);
                        //retx_bits.resize(header.payloadlen);
                        // used by comined data on first transmission
                        /*for(int i=0; i<header.payloadlen; i++){
                            data_comb[i] = 0;
                            retx_bits[i] = i;
                        }*/
                    //}

                    int bits_per_mod = MCS2BITSMOD[header.mcs];
                    set_currllr_table(bits_per_mod);
                }            
                else{
    	            enter_search();
                }
            }
            break;

        case STATE_HAVE_HEADER:

            int bits_per_mod = MCS2BITSMOD[header.mcs];
            //bit_demapper(&in[0], data_demod, d_data_offset,bits_per_mod);
            double t0 = current_time();


            //llr_demapper(&in[0], data_demod, d_data_offset, bits_per_mod);
            equalize_interpolate_dfe(in,NULL);
            d_data_offset += (Nsubs*bits_per_mod);

           
            // llr demapper equalizes in data
            for(size_t i=0; i<d_data_carriers.size();i++){ rcvd_symbols.push_back(in[d_data_carriers[i]]); }

            llrtime_elapsed += (current_time() - t0);


            if (d_data_offset == data_demod.size()){

                //std::cout<<"bit_demapper: "<<data_demod<<std::endl;
                //std::cout<<"llr_demapper: "<<llr_demod<<std::endl;
                // printf("llr time=%f\n",1.0e3*llrtime_elapsed);
                llrtime_elapsed = 0.0;

                // calculate reference symbol snr
                std::vector<double> refsnr;
                std::vector<double> ref_esno;
#if REFSNR == 1
                if( (RAMODE == 0) || (RAMODE ==2) ) {
                    calculate_refsnr(rcvd_symbols, refsnr, header.mcs);
                }
                else{
                    if(header.retx_count == 0){
                        calculate_refsnr(rcvd_symbols, refsnr, header.mcs);
                    }
                    else{
                        if(d_refid != header.verifyid){
#if OFFLINE == 0
                            refsnr = prev_refsnr;
                            for(size_t i=0; i<refsnr.size();i++){
                                ref_esno.push_back(std::pow(10.0,refsnr[i]/10.0));
                            }
#else
                            create_ref_symbols(firsttx_mcs,header.mcs,retx_subs,retx_refsymbs);
                            calculate_retx_refsnr(rcvd_symbols, retx_refsymbs,refsnr, header.mcs);
                            retx_refsymbs.resize(0);
#endif
              
                        }
                        else{
                            create_ref_symbols(firsttx_mcs,header.mcs,retx_subs,retx_refsymbs);
                            calculate_retx_refsnr(rcvd_symbols, retx_refsymbs,refsnr, header.mcs);
                            prev_refsnr = refsnr;
                            retx_refsymbs.resize(0);
                        }
                    }   
                }


                // printf("ref-snr\n");
                // for(size_t i=0; i<refsnr.size();i++){printf("%.1f, ",refsnr[i]);}
                // printf("\n");

                for(size_t i=0; i<refsnr.size();i++){
                    double tmp = std::pow(10.0,refsnr[i]/10);
                    ref_esno.push_back(tmp);
                }
                // printf("\n");

#endif
                // printf("ltf-dist\n");
                std::vector<double> ltfdistesno;
                std::vector<double> ltf_snr;
                calculate_ltfpreamble_distsnr(ltfdistesno,ltf_snr);

                // EWMA calculation
                if(ltfsnr_count==0){
                    for(size_t i=0; i<ltf_snr.size(); i++){
                        avg_ltfsnr[i] = ltf_snr[i];
                        std_ltfsnr[i] = 0;
                    }
                }
                else{
                    for(size_t i=0; i<ltf_snr.size(); i++){
                            avg_ltfsnr[i] = (ewma_alpha*ltf_snr[i]) + (1-ewma_alpha)*avg_ltfsnr[i];
                            std_ltfsnr[i] = (1-std_beta)*std_ltfsnr[i] + std_beta*(abs(ltf_snr[i] - avg_ltfsnr[i]));
                        //printf("(%.2f,%.2f), ",avg_ltfsnr[i],std_ltfsnr[i]);
                    }
                    //printf("\n");
                }
               
#if 0
                // mean and std-dev reference (Donald Knuth)
                // http://www.johndcook.com/blog/standard_deviation/
                if(ltfsnr_count==0){
                    for(int i=0; i<ltf_snr.size(); i++){
                        avg_ltfsnr[i] = ltf_snr[i];
                        Sk_ltfsnr[i] = 0;
                    }
                }
                else{
                    for(int i=0; i<ltf_snr.size(); i++){
                        double Mk_1 = avg_ltfsnr[i];
                        avg_ltfsnr[i] = Mk_1+(ltf_snr[i]-Mk_1)/(ltfsnr_count+1);
                        Sk_ltfsnr[i] = Sk_ltfsnr[i]+(ltf_snr[i]- Mk_1)*(ltf_snr[i] - Mk_1);
                        std_ltfsnr[i] = std::sqrt(Sk_ltfsnr[i]/ltfsnr_count);
                    }
                }
 #endif
                ltfsnr_count+=1;

                for(size_t i=0; i<ltf_snr.size(); i++){ 
                //    printf("%.1f, ",ltf_snr[i]);
                    ltf_subsnr<<ltf_snr[i]<<", ";
                }
                ltf_subsnr<<std::endl;
                // printf("\n");

            

#if 0
                std::vector<double> distadjesno;
                std::vector<double> distadjsnr;
                calculate_adjdistsnr(rcvd_symbols, distadjesno, distadjsnr,header.mcs);
                printf("dist-adj-snr\n");
                for(int i=0; i<distadjsnr.size();i++){printf("%f, ",distadjsnr[i]);}
                printf("\n");

                printf("ltf-dist\n");
                std::vector<double> ltfdistesno;
                std::vector<double> ltf_snr;
                calculate_ltfpreamble_distsnr(ltfdistesno,ltf_snr);

                for(int i=0; i<ltf_snr.size(); i++){ printf("%f, ",ltf_snr[i]);}
                printf("\n");
                double diff_dist=0,diff_distadj=0;
                for(int i=0; i<refsnr.size(); i++){
                    diff_dist += refsnr[i] - distsnr[i];
                    diff_distadj+= refsnr[i] - distadjsnr[i];
                }
                diff_dist/=refsnr.size();
                diff_distadj/=refsnr.size();

                for(int i=0; i<ltf_subcarrier_esno.size();i++){printf("%f, ",(10.0*log10(ltf_subcarrier_esno[i])));}
                printf("\n");
                printf("diff: dist=%f, distadj=%f\n",diff_dist,diff_distadj);

                // calculate effsnr
                digital_ofdm_effsnr_ra eff_ra;
                double ref_effsnr = calculate_effsnr_from_subsnr(ref_esno, eff_ra, header.mcs);
                double dist_effsnr = calculate_effsnr_from_subsnr(distesno, eff_ra, header.mcs);
                double distadj_effsnr = calculate_effsnr_from_subsnr(distadjesno, eff_ra, header.mcs);
                double ltf_effsnr = calculate_effsnr_from_subsnr(ltf_subcarrier_esno, eff_ra, header.mcs);
                double ltfdist_effsnr = calculate_effsnr_from_subsnr(ltfdistesno, eff_ra, header.mcs);
                refEffsnrVec.push_back(ref_effsnr);
                distEffsnrVec.push_back(dist_effsnr);
                distadjEffsnrVec.push_back(distadj_effsnr);
                ltfEffsnrVec.push_back(ltf_effsnr);
                ltfdistEffsnrVec.push_back(ltfdist_effsnr);
                
                //printf("effsnr-> ref=%f, dist=%f\n",ref_effsnr,dist_effsnr);
               

                //llr_demapper_alldata(rcvd_symbols,ref_esno,data_demod,bits_per_mod);
                //llr_demapper_alldata(rcvd_symbols,ltf_subcarrier_esno,data_demod,bits_per_mod);
                // calculate snr from llr
                /*
                std::vector<double> llr_snr;
                calculate_llrsnr(data_demod, llr_snr, header.mcs);
                std::vector<double> llr_esno;
                for(int i=0;i<llr_snr.size();i++){
                    double tmp = std::pow(10.0,llr_snr[i]/10.0);
                    llr_esno.push_back(tmp);
                }
                */
                //printf("llr-snr\n");
                //for(int i=0; i<llr_snr.size();i++){ printf("%f, ",llr_snr[i]); }
               //printf("\n");

                //printf("llr\n");
                //for(int i=0; i<refsnr.size();i++){printf("%f, ",(refsnr[i]-llr_snr[i]));} printf("\n");

                //printf("ltf\n");
                //for(int i=0; i<refsnr.size();i++){printf("%f, ",(refsnr[i]-(10.0*log10(ltf_subcarrier_esno[i]))));} printf("\n");
                
#endif
                // dist snr
                std::vector<double> distesno;
                std::vector<double> distsnr;
                calculate_distsnr(rcvd_symbols, distesno, distsnr,header.mcs);
                // printf("dist-snr\n");
                for(size_t i=0; i<distsnr.size();i++){
                //    printf("%.1f, ",distsnr[i]);
                    dist_subsnr<<distsnr[i]<<", ";
                }
                dist_subsnr<<std::endl;
                // printf("\n");


                // k-dist snr
                std::vector<double> kdistesno;
                std::vector<double> kdistsnr;
                if(total_pkts<5){
                    kdistesno = ltfdistesno;
                    kdistsnr = ltf_snr;
                }
                if(total_pkts>=5){
                    calculate_kdistsnr(rcvd_symbols, kdistesno, kdistsnr,avg_ltfsnr,std_ltfsnr,header.mcs);

                }

                // printf("k-dist-snr\n");
                for(size_t i=0; i<kdistsnr.size();i++){
                //    printf("%.1f, ",kdistsnr[i]);
                    kdist_subsnr<<kdistsnr[i]<<", ";
                }
                kdist_subsnr<<std::endl;
                // printf("\n");


                // For SOFT combine symbols
                if (RAMODE == 2){
#if REFORLTF == 0                
                    soft_combine_data(rcvd_symbols,ref_esno,header.retx_count);
#elif REFORLTF == 1
                    soft_combine_data(rcvd_symbols,ltfdistesno,header.retx_count);
#else
                    soft_combine_data(rcvd_symbols,kdistesno,header.retx_count);
#endif
                    llr_demapper_alldata(soft_combined_data,cum_esno,data_demod,bits_per_mod);
                }
                else{
#if REFORLTF == 0                
                    llr_demapper_alldata(rcvd_symbols,ref_esno,data_demod,bits_per_mod);
#elif REFORLTF == 1
                    llr_demapper_alldata(rcvd_symbols,ltfdistesno,data_demod,bits_per_mod);
#else
                    llr_demapper_alldata(rcvd_symbols,kdistesno,data_demod,bits_per_mod);
#endif
                }
                rcvd_symbols.clear();

                
                // de-interleave data
                if( (RAMODE == 0) || (RAMODE == 2) ) {
                    deinterleaver.deinterleave_data(data_demod,data_comb,header.mcs);
                }
                else{
                    if(header.retx_count == 0){
                        firsttx_mcs = header.mcs;
                        deinterleaver.deinterleave_data(data_demod,data_comb,header.mcs);
                        assign_bits_to_subs(d_bit2submap, header.payloadlen, header.mcs,true);
                    }
                    else{
                        //data_comb = itpp::vec(data_demod);
                        combine_data(data_demod, retx_subs, data_comb);
                    }
                }

                retx_subs.resize(0);
                itpp::vec fec_input;
                int pktlen_with_crc = header.pktlength+4;
                int coded_data_len = 0;
                if ( (RAMODE == 0) || (RAMODE == 2)){
                   coded_data_len = get_codeddata_length(pktlen_with_crc,header.mcs);
                }
                else{
                    // printf("first_txmcs =%d\n",firsttx_mcs);
                    coded_data_len = get_codeddata_length(pktlen_with_crc,firsttx_mcs);
                }
                // printf("coded data length=%d\n",coded_data_len);
                fec_input.set_size(coded_data_len);
                //fec_input.set_size(data_comb.size());
                //for(int i=0; i<data_comb.size(); i++){fec_input[i] = (data_comb[i]>0)?0:1;}
                for(int i=0; i<fec_input.size(); i++){fec_input[i] = (data_comb[i]>0)?0:1;}
                // printf("fec input size: %d\n",fec_input.size());
                //for(int i=0; i<fec_input.size(); i++){printf("%d, ",(int)fec_input[i]);}
                //printf("\n");
                //for(int i=0; i<data_comb.size(); i++){printf("%d, ",fec_input[i]);}
#if REFCODEDDATA  == 1
                int diffcomb=0,difftx=0,total_bits;
                std::vector<int> difflocs;
                std::vector<int> diffsubs(Nsubs,0);
                int mcs_to_use;
                if( (RAMODE == 0) || (RAMODE == 2) ){ mcs_to_use = header.mcs;}
                else{ mcs_to_use = firsttx_mcs; }

                std::vector<int>& txRefCodedData = TxCodedData[mcs_to_use];
                // for (size_t j=0; j<txRefCodedData.size(); j+=1) { printf("%d, ",txRefCodedData[j]);}
                // printf("\n");
                total_bits = fec_input.size();
                for(size_t i=0; i<fec_input.size();i++){
                    /*if(header.retx_count>0){
                        int tmp = (data_demod[i] > 0) ? 0 : 1;
                        difftx+=   abs(tmp -txCodedPacket[i]);
                    }*/
//                    printf("(%d,%d) ",fec_input[i],txRefCodedData[i]);
                    int tmp= abs(fec_input[i]-txRefCodedData[i]);
                    diffcomb+=tmp;
                    if (tmp==1){
                        difflocs.push_back(i);
                    if( (RAMODE == 0) || (RAMODE == 2) ){
                            diffsubs[interleaver.get_subcarrier_index(i,header.mcs)]+=1;
                    }
                    else{
                            diffsubs[interleaver.get_subcarrier_index(i,firsttx_mcs)]+=1;
                    }
                        //diffsubs.push_back(interleaver.get_subcarrier_index(i,header.mcs));
                    }
                }
                printf("difftx: %d, diffcomb: %d\n",difftx,diffcomb);
                // printf("difflocs= ");
                // for(size_t i=0; i<difflocs.size();i++){ printf("%d,",difflocs[i]); }
                // printf("\n");
                // for(size_t i=0; i<diffsubs.size();i++){ printf("%d,",diffsubs[i]); }
                // printf("\n");
                difflocs.clear();
    #if TRACEFILE == 0
                diffsubs.clear();
    #endif
#endif

                itpp::bvec data_decode;
                unsigned int code_rate;
                if( (RAMODE == 0)|| (RAMODE == 2) ){ code_rate = MCS2CODERATE[header.mcs];}
                else{ code_rate = MCS2CODERATE[firsttx_mcs]; }

                fec_input = -2.0*(fec_input) + 1.0;
                fec_decoder.decode(code_rate,fec_input,data_decode);
                //fec_decoder.decode(code_rate,data_comb,data_decode);
                //fec_decoder.decode(code_rate,data_deinlv,data_decode);

                unsigned int crc_offset = header.pktlength*8;
                unsigned int pktdata_crc = (unsigned int) dec_convhelper(data_decode,crc_offset,32);

                //unsigned int pktdata_crc = 0x19c9a0d9;
                printf("pkt crc=%x\n",pktdata_crc);

                unsigned int dec_bytes = header.pktlength;
                unsigned char* dec_databytes = new unsigned char[dec_bytes]();
                int valid_bits = header.pktlength*8;
                conv_bvec_to_chardata(data_decode,valid_bits,dec_databytes);

                unsigned int data_crc = digital_crc32(dec_databytes,dec_bytes);

                //for(int i=0; i<dec_bytes;i++){ printf("%x",dec_databytes);}
                //printf("\n");
                printf("data crc=%x\n",data_crc);

                //FIXME: Hard coding pkt size for now. Fix it.
                /*int ra_mcs=0;
                for (int i=0; i<d_data_carriers.size();i++){
                    retx_subs.push_back(1);
                }*/

                int succ = (pktdata_crc==data_crc);

                if( (RAMODE == 0) || (RAMODE == 2) ){
                    total_time+=calculate_transmission_time(header.mcs, header.payloadlen, header.retx_count, -1);
                }
                else{
                    if(header.retx_count == 0){
                        total_time+=calculate_transmission_time(header.mcs, header.payloadlen, header.retx_count, -1);
                    }
                    else{
                        total_time+=calculate_transmission_time(header.mcs, header.payloadlen, header.retx_count,firsttx_mcs);
                    }
                }
                succrcvd_data+=(succ*header.pktlength*8);

                int ra_mcs;
                if(RAMODE == 0){
#if REFORLTF == 0
                    ra_mcs = ra_effsnr.select_rate(ref_esno,header.pktlength);
#elif REFORLTF == 1
                    ra_mcs = ra_effsnr.select_rate(ltfdistesno,header.pktlength);
#else
                    ra_mcs = ra_effsnr.select_rate(kdistesno,header.pktlength);
#endif
                //printf("ltf-snr ra\n");
                //int ra_mcs = ra_effsnr.select_rate(ltf_subcarrier_esno,header.pktlength);

                //int ra_mcs = ra_effsnr.select_rate(llr_esno,header.pktlength);
                //if(ra_mcs==7){ ra_mcs=5;}
                retx_subs.resize(Nsubs);
                for(int i=0; i<Nsubs; i++){ retx_subs[i]=1;}
                }
                else if (RAMODE == 2){
                    if (succ == 1){
#if REFORLTF == 0
                        ra_mcs = ra_effsnr.select_rate(ref_esno,header.pktlength);
#elif REFORLTF == 1
                        ra_mcs = ra_effsnr.select_rate(ltfdistesno,header.pktlength);
#else
                        ra_mcs = ra_effsnr.select_rate(kdistesno,header.pktlength);
#endif
                        //soft_prev_mcs = ra_mcs;
                    }
                    else{
                        ra_mcs = header.mcs;
                    }
                    retx_subs.resize(Nsubs);
                    for(int i=0; i<Nsubs; i++){ retx_subs[i]=1;}
                }
                else{
                    if(header.retx_count==0){
#if REFORLTF == 0
                        ra_mcs = ra_smart.select_rate(ref_esno,data_comb,retx_subs,header.retx_count,header.mcs,header.pktlength,header.pktid,succ);
#elif REFORLTF == 1
                        ra_mcs = ra_smart.select_rate(ltfdistesno,data_comb,retx_subs,header.retx_count,header.mcs,header.pktlength,header.pktid,succ);
#else
                        ra_mcs = ra_smart.select_rate(kdistesno,data_comb,retx_subs,header.retx_count,header.mcs,header.pktlength,header.pktid,succ);
#endif
                    }
                    else{
#if REFORLTF == 0
                        ra_mcs = ra_smart.select_rate(ref_esno,data_comb,retx_subs,header.retx_count,header.mcs,header.pktlength,header.pktid,succ);
#elif REFORLTF == 1
                        ra_mcs = ra_smart.select_rate(ltfdistesno,data_comb,retx_subs,header.retx_count,header.mcs,header.pktlength,header.pktid,succ);
#else
                        ra_mcs = ra_smart.select_rate(kdistesno,data_comb,retx_subs,header.retx_count,header.mcs,header.pktlength,header.pktid,succ);
#endif
                    }
                }


                int retxsub_count = 0;
                for(size_t i=0; i<retx_subs.size(); i++){
                    if(retx_subs[i]){retxsub_count++;}
                }
                printf("retx subs count=%d\n",retxsub_count);
                if(succ){
                    printf("pkt succ\n");
                    succ_pkts++;
                    d_lastpktsucc_id = header.pktid;
                    printf("d_lastpktsuccid=%d\n",d_lastpktsucc_id);
                }
                else{
                    printf("pkt fail\n");
                }
                total_pkts++;

#if TRACEFILE == 1
                traceMcs.push_back(header.mcs);
                tracePayloadLen.push_back(header.payloadlen);
                traceSuccPkt.push_back(succ);
                tracePktid.push_back(header.pktid);
                traceRetxSubs.push_back(retx_subs);
    #if REFCODEDDATA == 1
                traceBitErrors.push_back(diffcomb);
                traceTotalBits.push_back(total_bits);
                traceErrsPerSub.push_back(diffsubs);
                diffsubs.clear();
    #endif
#endif


#if OFFLINE == 0
                send_feedback(ra_mcs,succ,header.pktid,retx_subs);
#endif

                d_data_offset=0;
                free(dec_databytes);
    	        enter_search();// bad header
            }

            break;
                 
//        default:
//            assert(0);
 
    }

    return 1; 
}


inline void 
create_complex_array(int length, double real[], double imag[], std::vector<gr_complex>& out){

    for(int i=0; i<length; i++){
        out.push_back(gr_complex(real[i],imag[i]));
    }
}

//
void 
digital_ofdm_partialcomb_receiver::create_llr_tables(){

    // bpsk
    int length = 1;
    double S0_table_bpsk_b0_real[] = {-1};
    double S0_table_bpsk_b0_imag[] = {0};
    double S1_table_bpsk_b0_real[] = {1};
    double S1_table_bpsk_b0_imag[] = {0};
    std::vector<gr_complex> S0_table_bpsk_b0;
    std::vector<gr_complex> S1_table_bpsk_b0;
    create_complex_array(length, S0_table_bpsk_b0_real,S0_table_bpsk_b0_imag, S0_table_bpsk_b0);
    create_complex_array(length, S1_table_bpsk_b0_real,S1_table_bpsk_b0_imag, S1_table_bpsk_b0);
    S0_table_bpsk.push_back(S0_table_bpsk_b0);
    S1_table_bpsk.push_back(S1_table_bpsk_b0);


    // qpsk
    double S0_table_qpsk_b0_real[] = {-0.70711,-0.70711};
    double S0_table_qpsk_b0_imag[] = {0.70711 ,-0.70711};
    double S1_table_qpsk_b0_real[] = {0.70711 , 0.70711};
    double S1_table_qpsk_b0_imag[] = {0.70711 ,-0.70711};
    double S0_table_qpsk_b1_real[] = {-0.70711, 0.70711};
    double S0_table_qpsk_b1_imag[] = {-0.70711,-0.70711};
    double S1_table_qpsk_b1_real[] = {0.70711 ,-0.70711};
    double S1_table_qpsk_b1_imag[] = {0.70711 , 0.70711};

    length = 2;
    std::vector<gr_complex> S0_table_qpsk_b0;
    std::vector<gr_complex> S0_table_qpsk_b1;
    std::vector<gr_complex> S1_table_qpsk_b0;
    std::vector<gr_complex> S1_table_qpsk_b1;
    create_complex_array(length, S0_table_qpsk_b0_real,S0_table_qpsk_b0_imag, S0_table_qpsk_b0);
    create_complex_array(length, S0_table_qpsk_b1_real,S0_table_qpsk_b1_imag, S0_table_qpsk_b1);
    create_complex_array(length, S1_table_qpsk_b0_real,S1_table_qpsk_b0_imag, S1_table_qpsk_b0);
    create_complex_array(length, S1_table_qpsk_b1_real,S1_table_qpsk_b1_imag, S1_table_qpsk_b1);
    S0_table_qpsk.push_back(S0_table_qpsk_b0);
    S0_table_qpsk.push_back(S0_table_qpsk_b1);
    S1_table_qpsk.push_back(S1_table_qpsk_b0);
    S1_table_qpsk.push_back(S1_table_qpsk_b1);



    // qam16
    length = 8;
    double S0_table_qam16_b0_real[] = {-0.94868,-0.94868,-0.94868,-0.94868,-0.31623,-0.31623,-0.31623,-0.31623};
    double S0_table_qam16_b0_imag[] = {-0.94868,-0.31623, 0.94868, 0.31623,-0.94868,-0.31623, 0.94868, 0.31623};
    double S1_table_qam16_b0_real[] = { 0.94868, 0.94868, 0.94868, 0.94868, 0.31623, 0.31623, 0.31623, 0.31623};
    double S1_table_qam16_b0_imag[] = {-0.94868,-0.31623, 0.94868, 0.31623,-0.94868,-0.31623, 0.94868, 0.31623};
    double S0_table_qam16_b1_real[] = {-0.94868,-0.94868,-0.94868,-0.94868, 0.94868, 0.94868, 0.94868, 0.94868};
    double S0_table_qam16_b1_imag[] = {-0.94868,-0.31623, 0.94868, 0.31623,-0.94868,-0.31623, 0.94868, 0.31623};
    double S1_table_qam16_b1_real[] = {-0.31623,-0.31623,-0.31623,-0.31623, 0.31623, 0.31623, 0.31623, 0.31623};
    double S1_table_qam16_b1_imag[] = {-0.94868,-0.31623, 0.94868, 0.31623,-0.94868,-0.31623, 0.94868, 0.31623};
    double S0_table_qam16_b2_real[] = {-0.94868,-0.94868,-0.31623,-0.31623, 0.94868, 0.94868, 0.31623, 0.31623};
    double S0_table_qam16_b2_imag[] = {-0.94868,-0.31623,-0.94868,-0.31623,-0.94868,-0.31623,-0.94868,-0.31623};
    double S1_table_qam16_b2_real[] = {-0.94868,-0.94868,-0.31623,-0.31623, 0.94868, 0.94868, 0.31623, 0.31623};
    double S1_table_qam16_b2_imag[] = { 0.94868, 0.31623, 0.94868, 0.31623, 0.94868, 0.31623, 0.94868, 0.31623};
    double S0_table_qam16_b3_real[] = {-0.94868,-0.94868,-0.31623,-0.31623, 0.94868, 0.94868, 0.31623, 0.31623};
    double S0_table_qam16_b3_imag[] = {-0.94868, 0.94868,-0.94868, 0.94868,-0.94868, 0.94868,-0.94868, 0.94868};
    double S1_table_qam16_b3_real[] = {-0.94868,-0.94868,-0.31623,-0.31623, 0.94868, 0.94868, 0.31623, 0.31623};
    double S1_table_qam16_b3_imag[] = {-0.31623, 0.31623,-0.31623, 0.31623,-0.31623, 0.31623,-0.31623, 0.31623};

    std::vector<gr_complex> S0_table_qam16_b0;
    std::vector<gr_complex> S0_table_qam16_b1;
    std::vector<gr_complex> S0_table_qam16_b2;
    std::vector<gr_complex> S0_table_qam16_b3;
    std::vector<gr_complex> S1_table_qam16_b0;
    std::vector<gr_complex> S1_table_qam16_b1;
    std::vector<gr_complex> S1_table_qam16_b2;
    std::vector<gr_complex> S1_table_qam16_b3;
    create_complex_array(length, S0_table_qam16_b0_real,S0_table_qam16_b0_imag, S0_table_qam16_b0);
    create_complex_array(length, S0_table_qam16_b1_real,S0_table_qam16_b1_imag, S0_table_qam16_b1);
    create_complex_array(length, S0_table_qam16_b2_real,S0_table_qam16_b2_imag, S0_table_qam16_b2);
    create_complex_array(length, S0_table_qam16_b3_real,S0_table_qam16_b3_imag, S0_table_qam16_b3);
    create_complex_array(length, S1_table_qam16_b0_real,S1_table_qam16_b0_imag, S1_table_qam16_b0);
    create_complex_array(length, S1_table_qam16_b1_real,S1_table_qam16_b1_imag, S1_table_qam16_b1);
    create_complex_array(length, S1_table_qam16_b2_real,S1_table_qam16_b2_imag, S1_table_qam16_b2);
    create_complex_array(length, S1_table_qam16_b3_real,S1_table_qam16_b3_imag, S1_table_qam16_b3);
    S0_table_qam16.push_back(S0_table_qam16_b0);
    S0_table_qam16.push_back(S0_table_qam16_b1);
    S0_table_qam16.push_back(S0_table_qam16_b2);
    S0_table_qam16.push_back(S0_table_qam16_b3);
    S1_table_qam16.push_back(S1_table_qam16_b0);
    S1_table_qam16.push_back(S1_table_qam16_b1);
    S1_table_qam16.push_back(S1_table_qam16_b2);
    S1_table_qam16.push_back(S1_table_qam16_b3);



    //qam64
    double S0_table_qam64_b0_real[] = {-1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -0.7715, -0.7715,-0.7715, -0.7715, -0.7715, -0.7715, -0.7715, -0.7715, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629};
    double S0_table_qam64_b0_imag[] = {-1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543, 0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629};
    double S1_table_qam64_b0_real[] = {1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  0.7715,  0.7715, 0.7715,  0.7715,  0.7715,  0.7715,  0.7715,  0.7715,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629};
    double S1_table_qam64_b0_imag[] = {-1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715,-0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629};
    double S0_table_qam64_b1_real[] = {-1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -0.7715, -0.7715,-0.7715, -0.7715, -0.7715, -0.7715, -0.7715, -0.7715,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  0.7715,  0.7715,  0.7715,  0.7715,  0.7715,  0.7715,  0.7715,  0.7715};
    double S0_table_qam64_b1_imag[] = {-1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715,-0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629};
    double S1_table_qam64_b1_real[] = {-0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629};
    double S1_table_qam64_b1_imag[] = {-1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629};
    double S0_table_qam64_b2_real[] = {-1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  1.0801,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543,  0.1543};
    double S0_table_qam64_b2_imag[] = {-1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629};
    double S1_table_qam64_b2_real[] = {-0.7715, -0.7715, -0.7715, -0.7715, -0.7715, -0.7715, -0.7715, -0.7715, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629, -0.4629,  0.7715,  0.7715,  0.7715,  0.7715,  0.7715,  0.7715,  0.7715,  0.7715,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629,  0.4629};
    double S1_table_qam64_b2_imag[] = {-1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629, -1.0801, -0.7715, -0.1543, -0.4629,  1.0801,  0.7715,  0.1543,  0.4629};
    double S0_table_qam64_b3_real[] = {-1.0801, -1.0801, -1.0801, -1.0801, -0.7715, -0.7715, -0.7715, -0.7715, -0.1543, -0.1543, -0.1543, -0.1543, -0.4629, -0.4629, -0.4629, -0.4629,  1.0801,  1.0801,  1.0801,  1.0801,  0.7715,  0.7715,  0.7715,  0.7715,  0.1543,  0.1543,  0.1543,  0.1543,  0.4629,  0.4629,  0.4629,  0.4629};
    double S0_table_qam64_b3_imag[] = {-1.0801, -0.7715, -0.1543, -0.4629, -1.0801, -0.7715, -0.1543, -0.4629, -1.0801, -0.7715, -0.1543, -0.4629, -1.0801, -0.7715, -0.1543, -0.4629, -1.0801, -0.7715, -0.1543, -0.4629, -1.0801, -0.7715, -0.1543, -0.4629, -1.0801, -0.7715, -0.1543, -0.4629, -1.0801, -0.7715, -0.1543, -0.4629};
    double S1_table_qam64_b3_real[] = {-1.0801, -1.0801, -1.0801, -1.0801, -0.7715, -0.7715, -0.7715, -0.7715, -0.1543, -0.1543, -0.1543, -0.1543, -0.4629, -0.4629, -0.4629, -0.4629,  1.0801,  1.0801,  1.0801,  1.0801,  0.7715,  0.7715,  0.7715,  0.7715,  0.1543,  0.1543,  0.1543,  0.1543,  0.4629,  0.4629,  0.4629,  0.4629};
    double S1_table_qam64_b3_imag[] = {1.0801,  0.7715,  0.1543,  0.4629,  1.0801,  0.7715,  0.1543,  0.4629,  1.0801,  0.7715,  0.1543,  0.4629,  1.0801,  0.7715,  0.1543,  0.4629,  1.0801,  0.7715,  0.1543,  0.4629,  1.0801,  0.7715,  0.1543,  0.4629,  1.0801,  0.7715,  0.1543,  0.4629,  1.0801,  0.7715,  0.1543,  0.4629};
    double S0_table_qam64_b4_real[] = {-1.0801, -1.0801, -1.0801, -1.0801, -0.7715, -0.7715, -0.7715, -0.7715, -0.1543, -0.1543, -0.1543, -0.1543, -0.4629, -0.4629, -0.4629, -0.4629,  1.0801,  1.0801,  1.0801,  1.0801,  0.7715,  0.7715,  0.7715,  0.7715,  0.1543,  0.1543,  0.1543,  0.1543,  0.4629,  0.4629,  0.4629,  0.4629};
    double S0_table_qam64_b4_imag[] = {-1.0801, -0.7715,  1.0801,  0.7715, -1.0801, -0.7715,  1.0801,  0.7715, -1.0801, -0.7715,  1.0801,  0.7715, -1.0801, -0.7715,  1.0801,  0.7715, -1.0801, -0.7715,  1.0801,  0.7715, -1.0801, -0.7715,  1.0801,  0.7715, -1.0801, -0.7715,  1.0801,  0.7715, -1.0801, -0.7715,  1.0801,  0.7715};
    double S1_table_qam64_b4_real[] = {-1.0801, -1.0801, -1.0801, -1.0801, -0.7715, -0.7715, -0.7715, -0.7715, -0.1543, -0.1543, -0.1543, -0.1543, -0.4629, -0.4629, -0.4629, -0.4629,  1.0801,  1.0801,  1.0801,  1.0801,  0.7715,  0.7715,  0.7715,  0.7715,  0.1543,  0.1543,  0.1543,  0.1543,  0.4629,  0.4629,  0.4629,  0.4629};
    double S1_table_qam64_b4_imag[] = {-0.1543, -0.4629,  0.1543,  0.4629, -0.1543, -0.4629,  0.1543,  0.4629, -0.1543, -0.4629,  0.1543,  0.4629, -0.1543, -0.4629,  0.1543,  0.4629, -0.1543, -0.4629,  0.1543,  0.4629, -0.1543, -0.4629,  0.1543,  0.4629, -0.1543, -0.4629,  0.1543,  0.4629, -0.1543, -0.4629,  0.1543,  0.4629};
    double S0_table_qam64_b5_real[] = {-1.0801, -1.0801, -1.0801, -1.0801, -0.7715, -0.7715, -0.7715, -0.7715, -0.1543, -0.1543, -0.1543, -0.1543, -0.4629, -0.4629, -0.4629, -0.4629,  1.0801,  1.0801,  1.0801,  1.0801,  0.7715,  0.7715,  0.7715,  0.7715,  0.1543,  0.1543,  0.1543,  0.1543,  0.4629,  0.4629,  0.4629,  0.4629};
    double S0_table_qam64_b5_imag[] = {-1.0801, -0.1543,  1.0801,  0.1543, -1.0801, -0.1543,  1.0801,  0.1543, -1.0801, -0.1543,  1.0801,  0.1543, -1.0801, -0.1543,  1.0801,  0.1543, -1.0801, -0.1543,  1.0801,  0.1543, -1.0801, -0.1543,  1.0801,  0.1543, -1.0801, -0.1543,  1.0801,  0.1543, -1.0801, -0.1543,  1.0801,  0.1543};
    double S1_table_qam64_b5_real[] = {-1.0801, -1.0801, -1.0801, -1.0801, -0.7715, -0.7715, -0.7715, -0.7715, -0.1543, -0.1543, -0.1543, -0.1543, -0.4629, -0.4629, -0.4629, -0.4629,  1.0801,  1.0801,  1.0801,  1.0801,  0.7715,  0.7715,  0.7715,  0.7715,  0.1543,  0.1543,  0.1543,  0.1543,  0.4629,  0.4629,  0.4629,  0.4629};
    double S1_table_qam64_b5_imag[] = {-0.7715, -0.4629,  0.7715,  0.4629, -0.7715, -0.4629,  0.7715,  0.4629, -0.7715, -0.4629,  0.7715,  0.4629, -0.7715, -0.4629,  0.7715,  0.4629, -0.7715, -0.4629,  0.7715,  0.4629, -0.7715, -0.4629,  0.7715,  0.4629, -0.7715, -0.4629,  0.7715,  0.4629, -0.7715, -0.4629,  0.7715,  0.4629};

    length = 32;
    std::vector<gr_complex> S0_table_qam64_b0;
    std::vector<gr_complex> S0_table_qam64_b1;
    std::vector<gr_complex> S0_table_qam64_b2;
    std::vector<gr_complex> S0_table_qam64_b3;
    std::vector<gr_complex> S0_table_qam64_b4;
    std::vector<gr_complex> S0_table_qam64_b5;
    std::vector<gr_complex> S1_table_qam64_b0;
    std::vector<gr_complex> S1_table_qam64_b1;
    std::vector<gr_complex> S1_table_qam64_b2;
    std::vector<gr_complex> S1_table_qam64_b3;
    std::vector<gr_complex> S1_table_qam64_b4;
    std::vector<gr_complex> S1_table_qam64_b5;
    create_complex_array(length, S0_table_qam64_b0_real,S0_table_qam64_b0_imag, S0_table_qam64_b0);
    create_complex_array(length, S0_table_qam64_b1_real,S0_table_qam64_b1_imag, S0_table_qam64_b1);
    create_complex_array(length, S0_table_qam64_b2_real,S0_table_qam64_b2_imag, S0_table_qam64_b2);
    create_complex_array(length, S0_table_qam64_b3_real,S0_table_qam64_b3_imag, S0_table_qam64_b3);
    create_complex_array(length, S0_table_qam64_b4_real,S0_table_qam64_b4_imag, S0_table_qam64_b4);
    create_complex_array(length, S0_table_qam64_b5_real,S0_table_qam64_b5_imag, S0_table_qam64_b5);
    create_complex_array(length, S1_table_qam64_b0_real,S1_table_qam64_b0_imag, S1_table_qam64_b0);
    create_complex_array(length, S1_table_qam64_b1_real,S1_table_qam64_b1_imag, S1_table_qam64_b1);
    create_complex_array(length, S1_table_qam64_b2_real,S1_table_qam64_b2_imag, S1_table_qam64_b2);
    create_complex_array(length, S1_table_qam64_b3_real,S1_table_qam64_b3_imag, S1_table_qam64_b3);
    create_complex_array(length, S1_table_qam64_b4_real,S1_table_qam64_b4_imag, S1_table_qam64_b4);
    create_complex_array(length, S1_table_qam64_b5_real,S1_table_qam64_b5_imag, S1_table_qam64_b5);
    S0_table_qam64.push_back(S0_table_qam64_b0);
    S0_table_qam64.push_back(S0_table_qam64_b1);
    S0_table_qam64.push_back(S0_table_qam64_b2);
    S0_table_qam64.push_back(S0_table_qam64_b3);
    S0_table_qam64.push_back(S0_table_qam64_b4);
    S0_table_qam64.push_back(S0_table_qam64_b5);
    S1_table_qam64.push_back(S1_table_qam64_b0);
    S1_table_qam64.push_back(S1_table_qam64_b1);
    S1_table_qam64.push_back(S1_table_qam64_b2);
    S1_table_qam64.push_back(S1_table_qam64_b3);
    S1_table_qam64.push_back(S1_table_qam64_b4);
    S1_table_qam64.push_back(S1_table_qam64_b5);

}
