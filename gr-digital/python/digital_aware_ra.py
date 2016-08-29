#! /usr/bin/env python

from gnuradio.digital import digital_swig
import math
import ast

#import os
#print os.getpid()
#raw_input("Press enter to continue")

class digitalrateAdaptationAware:

    def __init__(self,nsubs,maxdepth,maxnodes):

        self.ra = digital_swig.digital_aware_ra(nsubs,maxdepth, maxnodes)
        self.mcs = None
        self.bitmap = None
        self.selected_subs = None


    def select_rate(self,bit_esno,sub_esno, pktlen, pktid, pkt_mcs, pkt_succ):
        #print "aware ra"
        self.mcs = self.ra.select_rate(bit_esno,sub_esno,pktlen,pktid, pkt_mcs, pkt_succ)
        #print "mcs,",self.mcs
        # This must be called after select_rate
        self.bitmap = self.ra.get_curr_bitmap()
        self.selected_subs = self.ra.get_selected_subcarriers()
        return [self.mcs, self.selected_subs, self.bitmap]

def main():

    fptr = open('snr_trace.dat')
    snrs=[]
    for line in fptr:
         data = ast.literal_eval(line)
         snrs.append(data['subesno'])
    fptr.close()

    #subsnr = [11.58, 8.31, 11.97, 15.47, 11.27, 11.37, 12.53, 10.41, 10.03, 10.91, 14.66, 12.28, 14.47, 12.22, 11.22, 9.10, 7.07, 10.35, 14.02, 16.77, 11.80, 11.06, 15.05, 15.42, 15.58, 14.69, 18.30, 15.45, 13.53, 10.01, 11.83, 15.33, 15.11, 17.13, 14.16, 16.39, 13.95, 14.77, 13.32, 14.97, 13.24, 14.15, 15.00, 14.92, 14.41, 12.48, 12.40, 10.05, 15.31, 11.53, 4.77, 10.81, 13.63, 13.15, 13.49, 17.53, 14.32, 14.45, 9.79, 12.98, 15.09, 12.04, 12.13, 11.21, 13.16, 10.76, 12.20, 8.43, 14.06, 11.02, 13.51, 15.58, 11.01, 12.37, 12.13, 12.79, 12.07, 12.80, 12.70, 12.49]


    #subsnr2= [-5.13479,-5.13479,-0.196516,-0.196516,9.91671,9.91671,14.7012,17.2324,17.2324,17.9288,17.9288,17.1015,17.1015,13.6124,13.6124,8.59217,8.59217,12.7017,12.7017,15.2582,16.9461,16.9461,15.8608,15.8608,13.3993,11.9659,8.32414,8.32414,11.3746,11.3746,13.7202,13.7202,16.1267,15.396,15.396,13.4984,13.4984,9.7291,9.7291,8.75166,8.75166,10.5816,10.5816,13.6901,13.6901,15.5558,16.1874,16.1874,16.4789,16.4789,16.2529,16.2547]

    #subcarrier_esno = [10.0**(s/10.0) for s in subsnr]
    #subcarrier_esno2 = [10.0**(s/10.0) for s in subsnr2]
    subsnr = [12.31, 11.89, 14.18, 12.08, 13.87, 12.65, 12.85, 12.58, 12.55, 10.43, 13.46, 14.07, 11.97, 12.40, 12.68, 11.56, 11.01, 13.20, 14.94, 15.09, 14.12, 14.31, 15.05, 12.92, 13.76, 10.29, 12.61, 13.49, 14.36, 14.47, 14.85, 14.35, 13.63, 15.45, 15.73, 12.77, 11.75, 13.44, 14.34, 14.04, 11.66, 9.35, 13.51, 12.61, 14.78, 14.09, 14.18, 11.12, 13.18, 14.01, 14.08, 15.16, 11.79, 11.41, 11.54, 14.14, 14.31, 13.97, 11.32, 11.68, 13.73, 12.14, 11.98, 12.35, 11.68, 12.24, 12.43, 11.66, 13.12, 13.98, 11.98, 12.72, 11.36, 5.62, 0.56, 2.55, 4.96, 6.19, 7.14, 6.52]

    #subsnr = [18.35, 17.91, 15.50, 15.22, 17.73, 17.86, 17.83, 18.64, 15.23, 13.31, 17.17, 19.06, 18.69, 17.37, 18.52, 18.17, 18.84, 18.62, 16.23, 15.95, 20.29, 18.28, 17.74, 19.30, 19.00, 19.63, 18.01, 17.60, 18.00, 15.69, 10.29, 15.72, 13.78, 17.71, 17.62, 18.15, 19.10, 17.92, 16.56, 16.40, 17.30, 17.89, 13.24, 17.38, 15.82, 18.10, 17.87, 16.82, 17.48, 17.94, 15.21, 14.95, 14.52, 17.01, 18.27, 17.27, 18.02, 17.37, 18.35, 20.36, 19.13, 17.35, 18.08, 16.70, 17.23, 17.94, 16.11, 17.02, 17.40, 18.13, 16.74, 16.94, 15.52, 16.93, 18.55, 18.51, 15.81, 17.38, 11.72, 17.36]

    subsnr2 = [8.12, 7.19, 10.36, 11.48, 10.89, 9.55, 9.90, 9.48, 14.76, 11.61, 10.38, 13.38, 14.09, 14.37, 13.13, 14.22, 14.13, 14.13, 13.62, 15.40, 14.84, 13.71, 13.07, 12.66, 12.29, 13.42, 13.56, 13.32, 14.17, 14.72, 10.86, 13.45, 15.31, 14.38, 13.61, 13.07, 15.05, 15.00, 11.21, 12.64, 12.90, 12.77, 13.12, 12.63, 12.60, 13.68, 13.34, 13.74, 15.50, 13.60, 9.12, 13.49, 11.64, 11.03, 11.12, 12.86, 11.61, 14.34, 12.20, 14.74, 12.00, 11.16, 10.17, 10.42, 11.76, 14.59, 13.73, 12.69, 12.68, 13.09, 13.44, 12.46, 9.82, 11.50, 12.17, 11.52, 12.11, 11.75, 11.84, 12.30]

    maxdepth = 2
    maxnodes = 4
    pktmcs = 0
    Nsubs = 80
    pktlen = 1000
    CodeFactor = {'1/2':2.0, '2/3':3.0/2, '3/4':4.0/3, '5/6':6.0/5}
    MCS2CODERATE = ['1/2','1/2','3/4','1/2','3/4','2/3','3/4','5/6']
    INTF = 0
    intf_subs = range(5,10)

    code_rate = MCS2CODERATE[pktmcs]


    ra = digitalrateAdaptationAware(Nsubs, maxdepth, maxnodes)
    pkt_succ = [True,False]

    snrs = [subsnr,subsnr2]
    for pktnum,snr in zip(range(0,len(snrs)),snrs):

        subcarrier_esno = snr

        if INTF==1:
            for isub in intf_subs:
                subcarrier_esno[isub] = subcarrier_esno[isub]/2.0

        avg_snr = sum([10.0*math.log10(esno) for esno in subcarrier_esno])/len(subcarrier_esno)
        print 'pktnum ->'+str(pktnum)+' snr(db) ->'+str(avg_snr)
        #print subcarrier_esno
        bit_esno_base = get_bit_esno_base(pktmcs,pktlen,subcarrier_esno)
        ra.select_rate(bit_esno_base, subcarrier_esno, pktlen,0, pktmcs ,pkt_succ[pktnum])
        print ra.mcs
        #print ra.bitmap
        print len(ra.bitmap)

#    #ra.select_rate(subcarrier_esno2, 200,0, pktmcs,False)
    #bit_esno_base = get_bit_esno_base(ra.mcs,pktlen,subcarrier_esno2)
    #pkt_succ = False
    #ra.select_rate(bit_esno_base, subcarrier_esno, pktlen,0, ra.mcs ,pkt_succ)
    #print ra.mcs
    #print ra.bitmap
    #print len(ra.bitmap)

    #ra.select_rate(subcarrier_esno2, 200,0,pktmcs,False)
    #print ra.mcs
    #print ra.bitmap
    #print len(ra.bitmap)

def get_bit_esno_base(pktmcs,pktlen,subcarrier_esno):

    CodeFactor = {'1/2':2.0, '2/3':3.0/2, '3/4':4.0/3, '5/6':6.0/5}
    MCS2CODERATE = ['1/2','1/2','3/4','1/2','3/4','2/3','3/4','5/6']

    bitmap = get_rate_bitmap(pktmcs,pktlen)
    baselen = int(math.ceil(((pktlen*8)+6)*CodeFactor['1/2']))
    bit_esno = assign_esno_per_bit(subcarrier_esno,'BPSK',bitmap,baselen)
    bit_esno_base = [0]*baselen
    for index in range(0,len(bit_esno)):
        bit_esno_base[bitmap[index]] = bit_esno[index]
    return bit_esno_base
 
def get_rate_bitmap(mcs, pktsize):

    CodeFactor = {'1/2':2.0, '2/3':3.0/2, '3/4':4.0/3, '5/6':6.0/5}
    MCS2CODERATE = ['1/2','1/2','3/4','1/2','3/4','2/3','3/4','5/6']

    code_rate =  MCS2CODERATE[mcs]
    codedpktsize = math.ceil(( (pktsize*8) +6)*CodeFactor[code_rate])

    if code_rate == '1/2':
        basecodedpktsize = int(math.ceil(( (pktsize*8) +6)*CodeFactor['1/2']))
        base_map = range(0,basecodedpktsize)
        return base_map
    elif code_rate == '2/3':
        puncture = [1,1,1,0]

    elif code_rate == '3/4':
        puncture = [1,1,1,0,0,1]

    elif code_rate == '5/6':
        puncture = [1,1,1,0,0,1,1,0,0,1]

    counter = 0
    bit_counter = 0
    rate_bitmap = []
    while 1:
        pvalue = puncture[counter%len(puncture)]
        if pvalue == 1:
            rate_bitmap.append(counter) 
            bit_counter+=1
        counter +=1
        if bit_counter >=codedpktsize:
            break;

    return rate_bitmap
 
def assign_esno_per_bit(esno,mod,bitmap,baselen):
    
    bits_in_modulation ={'BPSK':1, 'QPSK':2, 'QAM16':4, 'QAM64':6}
    # first assign each bit to a subcarrier
    bits_per_symbol = bits_in_modulation[mod]
    numofbits = len(bitmap)
    numofsymbols = int(math.ceil(1.0*numofbits/bits_per_symbol))
    num_of_subs = len(esno)

    subassignment=[]
    for index in range(0,numofsymbols):
        sub_num = index%num_of_subs
        subassignment[index*bits_per_symbol:(index+1)*bits_per_symbol] = [sub_num]*bits_per_symbol
        
    esno_per_bit = [esno[sub] for sub in subassignment]
    esno_map=[0]*baselen
    for index,val in zip(bitmap, esno_per_bit): esno_map[index] = val

    return esno_map



if __name__ == "__main__":
    main()

