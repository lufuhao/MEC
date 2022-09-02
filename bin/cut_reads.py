#! /usr/bin/env python3

import os,sys

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'

def cut_reads(file01,file02):
    n =0
    outf = open(file02,"w")
    if is_gz_file(file01):
            print ("Info: detected gzipped file: "+str(file01))
            inf=gzip.open(file01, 'rt')
        else:
            print ("Info: detected flat input: "+str(file01))
            inf=open(file01, 'r')
    if re.match('.*\.[gG]{1,1}[zZ]{1,1}$', str(file02)):
        print ("Info: detected gzipped output: "+str(file02))
        outf=gzip.open(file02, mode='wt', compresslevel=9)
    else:
        print ("Info: detected flat output: "+file02)
        outf=open(file02, 'w')
    for line in inf:
        l = line.rstrip()
        n+=1
        if n%4==1:
            outf.write(l+"\n")
        elif n%4==2 or n%4==0:
            if len(l)>=100:
                outf.write(l[0:100]+"\n")
            else:
                outf.write(l+"\n")
        elif n%4==3:
            outf.write("+\n")
        else:
            print ("Error: this error should not happen at line number: "+ str(n))
            sys.exit(100)
            
    inf.close()
    outf.close()

    return 0

def main():
    # read command line arguments
    import argparse
    parser = argparse.ArgumentParser(description='Tool to cut the paired-end reads')
    parser.add_argument('r1', metavar='read1', type=str, help='Input: Paired R1 left read')
    parser.add_argument('r2', metavar='read2', type=str, help='Input: Paired R2 left read')
    parser.add_argument('o1', metavar='out1', type=str, help='Output: Paired R1 left read')
    parser.add_argument('o2', metavar='out2', type=str, help='Output: Paired R2 left read')
    args = parser.parse_args()

    cut_reads01(args.r1,args.o1)
    cut_reads02(args.r2,args.o2)    
    
if __name__ == "__main__":
    main()
