#! /usr/bin/env python

from scipy import io

def main():

    #fptr = open("ltf_snr_bpsk.txt")
    fptr = open("ref_snr.txt")

    esno=[]
    for line in fptr:
        tmp = [float(val) for val in line.strip().split(',')[0:-1]]
        esno.append(tmp)
    #io.savemat('ltf_esno.mat',{'esno':esno})
    io.savemat('ref_esno.mat',{'esno':esno})


if __name__ == "__main__":
    main()
