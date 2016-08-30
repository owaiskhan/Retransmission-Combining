#! /usr/bin/env python

import os
import sys

# --rate-adapt ,Effsnr=0, Smart=1
def main(args):

    #cmd = 'python benchmark_rx.py --rx-freq 5007.996M --rx-gain 15 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 1.5e-6 -W1.0M'

    # First five runs -experiment LOS
    #cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 12 --fft-length 128 --occupied-tones 96 --cp-length 32 -m bpsk -v --threshold 6.0e-6 -W 1.2M'

    # Second set of experiments Non-LOS
    #cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m bpsk -v --threshold 3.0e-6 -W 1.2M'
    #cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 5e-5 -W 1.0M'

    # For tx-gain = 3
    #cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 15 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 3.0e-7 -W 1.0M'
    # For tx-gain = 5
    #cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 4.0e-6 -W 1.0M'
    # For tx-gain = 10
    #cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 8e-6 -W 1.0M'

    # args[0] : rate-adaptation EffSnr=0, Smart=1, Soft=2
    # args[1] : tx-gain used at the transmitter (Used for threshold setting)
    # args[2] : bw-intf based on interference overlap
    # args[3] : run number
    tx_gain = int(args[1])
    if(tx_gain==1): 
        cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 2.0e-6 -W 1.0M --rate-adapt '+args[0]+' --tx-gain '+args[1]+' --bw-intf '+args[2]+' --run '+args[3]
    elif(tx_gain==5): 
        cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 4.5e-6 -W 1.0M --rate-adapt '+args[0]+' --tx-gain '+args[1]+' --bw-intf '+args[2]+' --run '+args[3]
    elif(tx_gain==7): 
        cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 6.0e-6 -W 1.0M --rate-adapt '+args[0]+' --tx-gain '+args[1]+' --bw-intf '+args[2]+' --run '+args[3]
    elif (tx_gain == 10):
        cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 1.1e-5 -W 1.0M --rate-adapt '+args[0]+' --tx-gain '+args[1]+' --bw-intf '+args[2]+' --run '+args[3]
    elif (tx_gain == 12):
        cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 1.3e-5 -W 1.0M --rate-adapt '+args[0]+' --tx-gain '+args[1]+' --bw-intf '+args[2]+' --run '+args[3]
    elif (tx_gain == 15):
        cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 2.0e-5 -W 1.0M --rate-adapt '+args[0]+' --tx-gain '+args[1]+' --bw-intf '+args[2]+' --run '+args[3]
    elif (tx_gain == 17):
        cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 3.0e-5 -W 1.0M --rate-adapt '+args[0]+' --tx-gain '+args[1]+' --bw-intf '+args[2]+' --run '+args[3]
    elif(tx_gain==20): 
        cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 4.0e-5 -W 1.0M --rate-adapt '+args[0]+' --tx-gain '+args[1]+' --bw-intf '+args[2]+' --run '+args[3]
    else:
        print "tx-gain="+str(tx_gain)+" not supported."
        return

    # For tx-gain = 10 tx-amplitued=0.2
    #cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 1.0e-5 -W 1.0M    --rate-adapt '+args[0]


    # For tx-gain = 20 tx-amplitude=0.2
    #cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 9.0e-5 -W 1.0M    --rate-adapt '+args[0]
    #cmd = 'python benchmark_rx.py --rx-freq 2490.028M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 1.0e-6'
    #cmd = 'python benchmark_rx.py -v -m bpsk --rx-gain 25 --fft-length 64 --occupied-tones 48 --cp-length 16 --from-file sample_output_local.dat --log'
    os.system(cmd)


if __name__ == "__main__":
    main(sys.argv[1:])
