#! /usr/bin/python

import sys
import os
from multiprocessing import Pool

snp = sys.argv[1]
counts = sys.argv[2]
cis = sys.argv[3]
out = sys.argv[4]

chr_num = list( range(1,25) )
chr_num = [str(num) for num in chr_num]

def f(chr_num):
    cmd = 'fastQTL.static --vcf {} --bed {} --permute 1000 --window {} --region {}:1-9999999999999 --out {}'.format(snp, counts, cis, chr_num, out + chr_num)
    os.system(cmd)
    return cmd
    
if __name__ == '__main__':
    # run all the 
    p = Pool(25)
    print(p.map(f, chr_num))
    
    # concatinate all the files
    names = [out + num for num in chr_num]
    cmd = "zcat " + " ".join( names ) + " > " + out
    os.system(cmd)
    
    # remove the old files
    cmd = "rm " + " ".join( names )
    os.system(cmd)
    #print(cmd)
