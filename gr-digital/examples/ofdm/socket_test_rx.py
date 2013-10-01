#!/usr/bin/env python

import socket
import struct

def main():

    server_socket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    server_socket.bind(("",5309))
    server_socket.settimeout(1)

    flag = False
    while 1:
        try:
            data,address = server_socket.recvfrom(1000)
            print data.encode("hex")
            #type = struct.unpack_from('!3s',data,0)
            #print type
            decoded_data = struct.unpack_from('>3sHBB',data,0)
            print decoded_data

            if decoded_data[0] == 'NAK':
                bitmap_len = struct.unpack_from('>H',data,7)
                data = data[9:]
                bitmap = struct.unpack('>%sI'%bitmap_len[0],data)
                print bitmap


            print "decoded data ",decoded_data
            flag = True
            if flag == True:
                break;
        except socket.timeout:
            print "timeout"

#        if flag == True:
#            break;

if __name__ == "__main__":
    main()
