#! /usr/bin/env python

import os

def main():

    cmd = 'python benchmark_tx.py --tx-freq 2490M --tx-gain 15 --tx-amplitude 0.4 --size 200 --fft-length 128 --occupied-tones 96 --cp-length 32 -M 0.2 -m qpsk -v'
    #cmd = 'python rate_adaptation_tx.py --tx-freq 2490M --tx-gain 5 --tx-amplitude 0.4 --size 1000 --fft-length 128 --occupied-tones 96 --cp-length 32 -M 0.2 -v --scheme aware -W 2.0M'

    #cmd = 'python rate_adaptation_tx.py --tx-gain 4 --tx-amplitude 0.4 --size 10 --fft-length 128 --occupied-tones 96 --cp-length 32 -M 0.04 -v --scheme baseline -W 2.0M --to-file sample_output.dat --log'

    #cmd = 'python rate_adaptation_tx.py --tx-gain 15 --tx-amplitude 0.4 --size 200 --fft-length 128 --occupied-tones 96 --cp-length 32 -M 0.02 -v --to-file sample_output.dat --log'

    #cmd = 'python benchmark_tx.py --tx-freq 2490M --tx-gain 15 --tx-amplitude 0.4 --size 200 --fft-length 64 --occupied-tones 52 --cp-length 16 -M 0.02 --bandwidth 200k -v'
    #cmd = 'python benchmark_tx.py --tx-freq 2.4G --tx-gain 25 --tx-amplitude 0.3 --size 200 --fft-length 64 --occupied-tones 48 --cp-length 16 -M 0.2 -v'
    #cmd = 'python rate_adaptation_tx.py --tx-gain 15 --tx-amplitude 0.4 --size 200 --fft-length 128 --occupied-tones 96 --cp-length 32 --to-file sample_output.dat -M 0.002 --log'
    #cmd = 'python rate_adaptation_tx.py --tx-gain 25 --tx-amplitude 0.4 --size 200 --fft-length 64 --occupied-tones 56 --cp-length 16 --to-file sample_output.dat -M 0.0002 --log'
    os.system(cmd)


if __name__ == "__main__":
    main()
