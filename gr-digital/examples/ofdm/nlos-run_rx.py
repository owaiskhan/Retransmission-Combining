#! /usr/bin/env python

import os

def main():

    #cmd = 'python benchmark_rx.py --rx-freq 5007.996M --rx-gain 15 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 1.5e-6 -W1.0M'

    # First five runs -experiment LOS
    #cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 12 --fft-length 128 --occupied-tones 96 --cp-length 32 -m bpsk -v --threshold 6.0e-6 -W 1.2M'

    # Second set of experiments Non-LOS
    #cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m bpsk -v --threshold 3.0e-6 -W 1.2M'
    #cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 5e-5 -W 1.0M'

    # For tx-gain = 3
    #cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 15 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 3.0e-7 -W 1.0M'
    # For tx-gain = 5
    #cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 2.0e-6 -W 1.0M'
    # For tx-gain = 10
    #cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 8e-6 -W 1.0M'

    # For tx-gain = 15
    cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 5.0e-7 -W 1.0M'
    # For tx-gain = 20
    #cmd = 'python benchmark_rx.py --rx-freq 2490.002M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 2.2e-5 -W 1.0M'


    #cmd = 'python benchmark_rx.py --rx-freq 2490.028M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 1.0e-6'
    #cmd = 'python benchmark_rx.py -v -m bpsk --rx-gain 25 --fft-length 64 --occupied-tones 48 --cp-length 16 --from-file sample_output_local.dat --log'
    os.system(cmd)


if __name__ == "__main__":
    main()
