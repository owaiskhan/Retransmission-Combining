#ifndef INCLUDED_DIGITAL_OFDM_GLOBAL_H
#define INCLUDED_DIGITAL_OFDM_GLOBAL_H 


#include <itpp/itbase.h>
#include <sys/time.h>

#define WIFITIME 0

// global variables
const int BITSPERMOD[] = {1,2,4,6};
const double RATELIST[] = {6.5, 13, 19.5, 26, 39, 52, 58.5, 65.0};
const double CBITSPERMCS[] = {0.5,1.0,1.5,2.0,3.0,4.0,4.5,5.0};
const int MOD2MCS[] = {0, 2, 4, 7};
const int MCS2BITSMOD[8] = {1,2,2,4,4,6,6,6};


enum Modulation {BPSK, QPSK, QAM16, QAM64};
enum CodeRate {Oneby2, Twoby3, Threeby4, Fiveby6};

const Modulation MCS2MODULATION[] = {BPSK, QPSK, QPSK, QAM16, QAM16, QAM64, QAM64, QAM64};
const CodeRate MCS2CODERATE[] = {Oneby2, Oneby2, Threeby4, Oneby2, Threeby4, Twoby3, Threeby4, Fiveby6};


double current_time();

#endif /* INCLUDED_DIGITAL_OFDM_GLOBAL_H */
