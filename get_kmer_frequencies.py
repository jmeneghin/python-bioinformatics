#!/usr/bin/python3
###################################################################################
### Get kmer Frequencies                                                        ###
### Usage: get_kmer_frequencies.pl -i <fasta file> -p <prefix> -k <kmer length> ###
### This program takes a fasta file, k and a prefix as it's parameters.         ###
###                                                                             ###
### It returns a tab delimited file of kmer counts.                             ###
### Revers complimented kmers are summed                                        ###
###                                                                             ###
### Jennifer Meneghin                                                           ###
### Original Perl Code: 05/11/2015                                              ###
### 01/28/2020                                                                  ###
###################################################################################

import re, sys, getopt

def usage ():
    usage = "\nGet kmer Frequencies\n"
    usage = usage + "\nUsage: get_kmer_frequencies.pl -i <fasta file> -p <prefix> -k <kmer length>\n"
    usage = usage + "\nThis program takes a fasta file, k and prefix as it's parameters.\n\n"
    usage = usage + "It returns a tab delimited file (prefix_kmers.txt) of kmer counts. (columns = records, rows = kmer counts.)\n\n"
    usage = usage + "Jennifer Meneghin\n"
    usage = usage + "Janurary 28, 2020\n\n"
    return usage

def rc_seq(mykmer, k):
    rcmykmer = ""
    i = int(k)-1
    while i >= 0:
        if mykmer[i] == "A":
            rcmykmer = rcmykmer + "T"
        elif mykmer[i] == "C":
            rcmykmer = rcmykmer + "G"
        elif mykmer[i] == "G":
            rcmykmer = rcmykmer + "C"
        elif mykmer[i] == "T":
            rcmykmer = rcmykmer + "A"
        else:
            rcmykmer = rcmykmer + mykmer[i]
        i = i - 1
    return rcmykmer

def process_it(knucs, seq, k, header):
    i = 0
    seq = seq + seq[0:int(k)-1] #Tack the first kmer onto the end to get the wrap around kmers
    end = len(seq) - int(k) + 1
    while i < end:
        thiskmer = seq[i:i+int(k)]
        i+=1
        rckmer = rc_seq(thiskmer,k)
        if thiskmer <= rckmer:
            key = header + "\t" + thiskmer
        else:
            key = header + "\t" + rckmer
        if key in knucs:
            knucs[key] = knucs[key] + 1
        else:
            knucs[key] = 1
    return knucs

def main(argv):
    #---------------------------
    #Read command line arguments
    #---------------------------
    in_file = ""
    k = 4
    prefix = "My"
    try:
        opts, args = getopt.getopt(argv,"hi:k:p:",["ifile=","kint=","pstr="])
    except getopt.GetoptError:
        print("\nNot a valid argument or value")
        print("get_kmer_frequencies.pl -i <fasta file> -p <prefix> -k <kmer length>\n")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print(usage())
            sys.exit()
        elif opt in ("-i", "--ifile"):
            in_file = arg
        elif opt in ("-k", "--kint"):
            k = arg
        elif opt in ("-p", "--pstr"):
            prefix = arg

    #----------------------------------
    #Open Files for reading and writing
    #----------------------------------
    try:
        IN = open(in_file,"r")
    except FileNotFoundError:
        print("\nFASTA Input File could not be found")
        print("get_kmer_frequencies.pl -i <fasta file> -p <prefix> -k <kmer length>\n")
        sys.exit(2)

    out_file = prefix + "_kmers.txt";
    OUT = open(out_file,"w")
    print("\nParameters:\nfasta file = "+in_file)
    print("k = "+str(k))
    print("prefix = "+str(prefix))
    print("output file= "+str(out_file)+"\n")

    #--------------
    #The main event
    #--------------
    seq = ""
    pc = 0
    linecount = 0
    knucs = {}
    
    print("Reading FASTA File...")
    lines = IN.readlines()
    for line in lines:
        line=line.rstrip()
        if re.search("^>",line):
            if len(seq) > 0:
                knucs  = process_it(knucs, seq, k, header)
                seq = ""
            ll = len(line)
            header = line[1:ll] #remove the > from the front of the line
            pc+=1
            if pc % 100 == 0:
                print("record count = "+str(pc))
        else:
            seq = seq + line.upper()
        linecount+=1
        if linecount % 10000 == 0:
            print("line count = "+str(linecount))
    knucs = process_it(knucs, seq, k, header)
    
    print("Sorting and Counting...")
    kmers = {}
    records = {}
    for key, value in sorted(knucs.items()):
        parts = re.split("\t", key)
        record = parts[0]
        kmer = parts[1]
        if kmer in kmers:
            kmers[kmer] = kmers[kmer] + 1
        else:
            kmers[kmer] = 1
        if record in records:
            records[record] = records[record] + 1
        else:
            records[record] = 1

    print("Printing...")
    OUT.write(str(k)+"-mer")
    for i, value in sorted(records.items()):
        OUT.write("\t"+str(prefix)+"_"+str(i))
    OUT.write("\n")
    for i, value in sorted(kmers.items()):
        OUT.write(i)
        for j, value2 in sorted(records.items()):
            key = j + "\t" + i
            if key in knucs:
                OUT.write("\t"+str(knucs[key]))
            else:
                OUT.write("\t0")
        OUT.write("\n")
            
    IN.close()
    OUT.close()
    
if __name__ == "__main__":
    main(sys.argv[1:])
