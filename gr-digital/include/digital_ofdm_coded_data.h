#include<iostream>

#ifndef INCLUDED_DIGITAL_OFDM_CODED_DATA_H
#define INCLUDED_DIGITAL_OFDM_CODED_DATA_H


int txCodedPacketMcs0Size=1820;
int txCodedPacketMcs2Size=1214;
int txCodedPacketMcs5Size=1365;
int txCodedPacketMcs7Size=1092;

std::vector<std::vector<int> > TxCodedData;

/*
*/
void read_codeddatafile(std::string filename, std::vector<int >& coded_bits){

    std::ifstream datafile(filename.c_str());
    std::string line;

    if(datafile.is_open()){
            getline(datafile,line);
    }
    else{
        std::cout<<"Failed to open "<<filename<<std::endl;
        exit(0);
    }

    while(1){

        size_t location = line.find(", ");

        if (location == std::string::npos){ break; }
        std::string value_str = line.substr(0,location);

        int value = std::atoi(value_str.c_str());

        coded_bits.push_back(value);
        line = line.substr(location+1,line.length());
    }
}

/*
*/
void create_refcodeddata(int mcs, int pktSize){

    std::vector<int> txcodeddata;


    std::string filename = "RefCodedDataFiles/coded_data_mcs";
    std::stringstream outmcs; outmcs<<mcs;
    filename+=std::string(outmcs.str());
    filename+=std::string("_pktsize");
    std::stringstream outpktsize; outpktsize<<pktSize;
    filename+=std::string(outpktsize.str());
    filename+=std::string(".txt");

    read_codeddatafile(filename, txcodeddata);

    TxCodedData.push_back(txcodeddata);
}

#endif
