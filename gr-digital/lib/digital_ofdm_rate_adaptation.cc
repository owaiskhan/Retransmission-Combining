/*
 * Title: Rate adaptation
 * Created By: Owais Khan
 * Creation Date: 3/2/2015
 * 
 * Description: Base Class for Rate adaptaiton
 *
*/ 

#include "digital_ofdm_rate_adaptation.h"

digital_ofdm_rate_adaptation::digital_ofdm_rate_adaptation(){ }

#if WIFITIME == 1
digital_ofdm_rate_adaptation::digital_ofdm_rate_adaptation(int nSubs):selected_mcs(0),tpreamble(24.0e-6),tSIFS(16.0e-6),tSlottime(9.0e-6),
tOfdm(4.0e-6){
#else
digital_ofdm_rate_adaptation::digital_ofdm_rate_adaptation(int nSubs):selected_mcs(0),tOfdm(128.0e-6){

    tpreamble = 6*tOfdm;
#endif

    MACsize = 28*8;
    ACKsize = 14*8;
    PHYsize = 22;
#if WIFITIME == 1
    tDIFS = tSIFS + 2*tSlottime;
#endif

    numRates = 8;
    Nsubs = nSubs;

    read_effsnr2dr_file(drRate[0],"./FecTables/effsnr_rate6_5.dat");
    read_effsnr2dr_file(drRate[1],"./FecTables/effsnr_rate13.dat");
    read_effsnr2dr_file(drRate[2],"./FecTables/effsnr_rate19_5.dat");
    read_effsnr2dr_file(drRate[3],"./FecTables/effsnr_rate26.dat");
    read_effsnr2dr_file(drRate[4],"./FecTables/effsnr_rate39.dat");
    read_effsnr2dr_file(drRate[5],"./FecTables/effsnr_rate52.dat");
    read_effsnr2dr_file(drRate[6],"./FecTables/effsnr_rate58_5.dat");
    read_effsnr2dr_file(drRate[7],"./FecTables/effsnr_rate65.dat");

}

/*
*/
double
digital_ofdm_rate_adaptation::dratio_helper(double esno, int mcs, int pktsize_idx){

    double dratio=0;

    double snr = 10*log10(esno);
    snr = std::min(snr,30.0);
    snr = std::max(snr,-10.0);

    int snr_sc10 = (int) itpp::round(snr*10.0);
    double snr_rnd = snr_sc10/10.0;
    int index = (int) ceil( (snr_rnd+10)/0.2);


    if( (snr_sc10%2) == 0 ){
        dratio = drRate[mcs][pktsize_idx][index];
    }
    else{

        double dratio1 = drRate[mcs][pktsize_idx][index-1];
        double dratio2 = drRate[mcs][pktsize_idx][index];
        double snr1 = (index-1)*0.2 -10;
        double snr2 = (index)*0.2 - 10;
        double esno1 = pow(10,snr1/10);
        double esno2 = pow(10,snr2/10);
        double slope = (dratio2 - dratio1)/(esno2 - esno1);
        double intercept = dratio1 - (slope*esno1);
        dratio = (slope*esno) + intercept;
    }

    return dratio;
}

/*
*/
double 
digital_ofdm_rate_adaptation::calculate_delivery_ratio(double esno, int mcs, int pktSize){

    double dratio = 0;

    int supPktSizes[] = {25,50,100,200,400,600,800,1000,2000,4000};
    int numSizes = 10;
    int lower_idx = -1;
    int upper_idx = -1;


    int pktsize_idx= -1;
    for(int i=0; i<numSizes; i++){

        if(pktSize<=supPktSizes[numSizes-i-1]){
            lower_idx = numSizes-i-1;
        }

        if(pktSize>=supPktSizes[i]){
            upper_idx = i;
        }
    }

    if(lower_idx == upper_idx){
        pktsize_idx = lower_idx;
        dratio = dratio_helper(esno, mcs, pktsize_idx);
    }
    else{
        if(lower_idx == -1){
            upper_idx = numSizes-1;
            lower_idx = numSizes-2;
        }
        else if(upper_idx == -1){
            upper_idx = 1;
            lower_idx = 0;
        }
        double dratio1 = dratio_helper(esno, mcs, lower_idx);
        double dratio2 = dratio_helper(esno, mcs, upper_idx);
        int pktsize1 = supPktSizes[lower_idx];        
        int pktsize2 = supPktSizes[upper_idx];        
        double slope = (dratio2-dratio1)/(pktsize2 - pktsize1);
        double intercept = dratio1 - (slope*pktsize1);
        dratio = (slope*pktSize) + intercept;

    }
    return dratio;
}

/*
*/
#if 0
double 
digital_ofdm_rate_adaptation::calculate_delivery_ratio(double esno, int mcs, int pktSize){


    double dratio = 0;

    int supPktSizes[] = {25,50,100,200,400,600,800,1000,2000,4000};

    int pktsize_idx= -1;
    for(int i=0; i<10; i++){
        if(pktSize<=supPktSizes[i]){
            pktsize_idx = i;
            break;
        }
    /*    if(pktSize==supPktSizes[i]){
            pktsize_idx = i;
            break;
        }*/
    }
    if(pktsize_idx == -1){
        std::cout<<"Invalid pktsize."<<std::endl;
        exit(0);
    }

    double snr = 10*log10(esno);
    snr = std::min(snr,30.0);
    snr = std::max(snr,-10.0);

    int snr_sc10 = (int) itpp::round(snr*10.0);
    double snr_rnd = snr_sc10/10.0;
    int index = (int) ceil( (snr_rnd+10)/0.2);


    if( (snr_sc10%2) == 0 ){
        dratio = drRate[mcs][pktsize_idx][index];
    }
    else{

        double dratio1 = drRate[mcs][pktsize_idx][index-1];
        double dratio2 = drRate[mcs][pktsize_idx][index];
        double snr1 = (index-1)*0.2 -10;
        double snr2 = (index)*0.2 - 10;
        double esno1 = pow(10,snr1/10);
        double esno2 = pow(10,snr2/10);
        double slope = (dratio2 - dratio1)/(esno2 - esno1);
        double intercept = dratio1 - (slope*esno1);
        dratio = (slope*esno) + intercept;
    }

    return dratio;
}
#endif

/*
*/
void
digital_ofdm_rate_adaptation::read_effsnr2dr_file(double data[][201],const char* filename){

    std::ifstream drfile(filename);
    std::string line;

    int Rows = 10;
    int Cols = 201;


    if(drfile.is_open()){
        for (int r=0; r <Rows; r++){
            getline(drfile,line);
            for(int c=0; c<Cols; c++){

                size_t location = line.find(", ");
                std::string uber_str = line.substr(0,location);
                double dr = std::atof(uber_str.c_str());
                data[r][c] = dr;
                line = line.substr(location+1,line.length());
            }
        }
    }
    else{
        std::cout<<"Failed to open "<<filename<<std::endl;
        exit(0);
    }
}

/*
*/
int digital_ofdm_rate_adaptation::bits_per_modulation(Modulation mod){

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

double 
digital_ofdm_rate_adaptation::qfuncInv(double Qx){
return sqrt(2.0)*itpp::erfinv(1-(2*Qx));
}


/*
*/
double
digital_ofdm_rate_adaptation::calculate_effesno(double uber, Modulation mod){

    double const1[] ={0.5, 1, 5.0, 21.0};
    double const2[] ={1, 1, 1.333, 1.7143};

    double effesno = const1[mod]*pow(qfuncInv(const2[mod]*uber),2);
    return effesno;
}

