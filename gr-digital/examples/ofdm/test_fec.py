#! /usr/bin/env python

import struct
from gnuradio.digital import fec

def main():

    mcs = 7
    uncoded_pkt = 200*chr(3 & 0xff) #"Hello"*10
    coded_pkt = fec.fec_encode(uncoded_pkt, mcs)
    print "uncoded pkt: ",uncoded_pkt.encode("hex")

    print "coded pkt: ",coded_pkt.encode("hex")
    print "length ->",len(coded_pkt)
    decoded_pkt = fec.fec_decode(coded_pkt[0:], mcs)
    print "decoded pkt: ",decoded_pkt.encode("hex")
    print "length ->",len(decoded_pkt)

if __name__ == "__main__":
    main()
