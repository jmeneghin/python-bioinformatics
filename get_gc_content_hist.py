#!/usr/bin/python3
##########################################################################################
### Get GC Content                                                                     ###
### Usage: get_gc_content_hist.py -i <fasta input file> -o <tab delimited output file> ###
### This program reads a fasta file and returns a tab delimited file with:             ###
### column 1 = header (every thing between > and new line)                             ###
### column 2 = %gc content for the fasta entry                                         ###
### column 3 = total count of letters (upper or lower case)                            ###
### column 4 = total count of Gs (or gs)                                               ###
### column 5 = total count of Cs (or cs)                                               ###
### column 6 = total count of As (or as)                                               ###
### column 7 = total count of Ts (or ts)                                               ###
###                                                                                    ###
### Jennifer Meneghin                                                                  ###
### 01/20/2020                                                                         ###
###                                                                                    ###
### Update:                                                                            ###
### This script now also displays a histogram of % GC Content and a line graph         ###
### of nucleotide counts                                                               ###
###                                                                                    ###
### Jennifer Meneghin                                                                  ###
### 02/11/2020                                                                         ###
###                                                                                    ###
##########################################################################################

import re, sys, getopt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#################
#Process Sequence
#################
def process_seq(seq):
    seq = seq.upper()
    acount = len(re.findall("A",seq))
    ccount = len(re.findall("C",seq))
    gcount = len(re.findall("G",seq))
    tcount = len(re.findall("T",seq))
    gccount = ccount + gcount
    totalcount = len(re.findall("[A-Z]",seq))
    if totalcount > 0:
        gccontent = (100 * gccount) / totalcount
    else:
        gccontent = 0
    return [gccontent, totalcount, gcount, ccount, acount, tcount]

###############
#The Main Event
###############
def main(argv):

    #---------------------------
    #Read command line arguments
    #---------------------------
    in_file = ""
    out_file = "gc_out.txt"
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print("\nNot a valid argument")
        print("get_gc_content.py -i <fasta format input file> -o <tab delimited output file>\n")
        sys.exit(2)        
    for opt, arg in opts:
        if opt == "-h":
            print("\nget_gc_content.py -i <fasta format input file> -o <tab delimited output file>\n")
            sys.exit()
        elif opt in ("-i", "--ifile"):
            in_file = arg
        elif opt in ("-o", "--ofile"):
            out_file = arg

    #----------------------------------
    #Open Files for reading and writing
    #----------------------------------
    try:
        IN = open(in_file,"r")
    except FileNotFoundError:
        print("\nFASTA Input File could not be found")
        print("get_gc_content.py -i <fasta format input file> -o <tab delimited output file>\n")
        sys.exit(2)

    #---------------------------
    #Read and Process FASTA file
    #---------------------------
    seq = ""
    lines = IN.readlines()
    gcresults = {}
    for line in lines:
        if re.search("^>",line):
            if len(seq) > 0:
                gcresults[hid] = process_seq(seq)
                seq = ""
            hid = line[1:len(line)-1]
        else:
            seq = seq + line
    gcresults[hid] = process_seq(seq)
    df = pd.DataFrame(gcresults,index=["%GC Content","Total Count","G Count","C Count","A Count","T Count"])
    df.T.to_csv(out_file,sep='\t')
    IN.close()

    #-------------------------------
    #Create and Write Histogram Plot
    #-------------------------------
    plt.hist(df.loc['%GC Content'])
    plt.title('GC Content Histogram')
    plt.show()

    plt.plot(df.loc['A Count'], label='A')
    plt.plot(df.loc['C Count'], label='C')
    plt.plot(df.loc['G Count'], label='G')
    plt.plot(df.loc['T Count'], label='T')
    plt.title('Nucleotide Counts')
    plt.legend(loc='upper right', bbox_to_anchor=(1.1,1.0))
    plt.show()
    
if __name__ == "__main__":
    main(sys.argv[1:])

