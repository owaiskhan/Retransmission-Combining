#! /usr/bin/env python

import os

def main():

    cmd = 'python benchmark_rx.py --rx-freq 2490.028M --rx-gain 20 --fft-length 128 --occupied-tones 96 --cp-length 32 -m qpsk -v --threshold 8.0e-6'
    #cmd = 'python benchmark_rx.py -v -m bpsk --rx-gain 25 --fft-length 64 --occupied-tones 48 --cp-length 16 --from-file sample_output_local.dat --log'
    os.system(cmd)


if __name__ == "__main__":
    main()
