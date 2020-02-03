#!/usr/bin/python3
#####################################################################################
### Get GC Content                                                                ###
### Usage: get_gc_content.pl -i <fasta input file> -o <tab delimited output file> ###
### This program reads a fasta file and returns a tab delimited file with:        ###
### column 1 = header (every thing between > and new line)                        ###
### column 2 = %gc content for the fasta entry                                    ###
### column 3 = total count of letters (upper or lower case)                       ###
### column 4 = total count of Gs (or gs)                                          ###
### column 5 = total count of Cs (or cs)                                          ###
### column 6 = total count of As (or as)                                          ###
### column 7 = total count of Ts (or ts)                                          ###
###                                                                               ###
### Jennifer Meneghin                                                             ###
### 01/20/2020                                                                    ###
###                                                                               ###
#####################################################################################

import re, sys, getopt

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
    OUT = open(out_file,"w")


    #---------------------------
    #Read and Process FASTA file
    #---------------------------
    OUT.write("ID\t%GC Content\tTotal Count\tG Count\tC Count\tA Count\tT Count\n")
    seq = ""
    lines = IN.readlines()
    for line in lines:
        if re.search("^>",line):
            if len(seq) > 0:
                result = process_seq(seq)
                OUT.write("%f\t%d\t%d\t%d\t%d\t%d\n" % (result[0],result[1],result[2],result[3],result[4],result[5]))
                seq = ""
            len_id = len(line)
            chopped_id = line[1:len_id-1] #removing  the ">" and the "\n" from the ends of the header line
            OUT.write(chopped_id + "\t")
        else:
            seq = seq + line
    #process last record's sequence.
    result = process_seq(seq)
    OUT.write("%f\t%d\t%d\t%d\t%d\t%d\n" % (result[0],result[1],result[2],result[3],result[4],result[5]))

    #-----------
    #Close Files
    #-----------
    IN.close()
    OUT.close()

if __name__ == "__main__":
    main(sys.argv[1:])

