#!/usr/bin/python3
####################################################################################
# This script removes duplicate records from a "short" format BLAST output file, ###
# and keeps only the "best" records                                              ###
# (sorts by smallest e-value and then biggest percent identity)                  ###
#                                                                                ###                
# Usage: blast_best.pl -i <input file> -o <output file>                          ###
#                                                                                ###
# Jennifer Meneghin                                                              ###
# Original Perl version: 08/14/2007                                              ###
# Python version: 01/20/2020                                                     ###
####################################################################################

import re, sys, getopt

###############
#The Main Event
###############
def main(argv):

    #---------------------------
    #Read command line arguments
    #---------------------------
    in_file = ""
    out_file = "best.blast"
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print("\nNot a valid argument or value")
        print("blast_best.py -i <BLAST input file> -o <BLAST output file>\n")
        sys.exit(2)        
    for opt, arg in opts:
        if opt == "-h":
            print("\nblast_best.py -i <BLAST input file> -o <BLAST output file>\n")
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
        print("\nBLAST Input File could not be found")
        print("blast_best.py -i <BLAST input file> -o <BLAST output file>\n")
        sys.exit(2)
    OUT = open(out_file,"w")

    
    #Everything looks good. Print the parameters we've found.
    print("\nParameters:\ninput file = "+in_file+"\noutput file = "+out_file+"\n")
    
    #--------------------------------------
    #Read and Process Short Form BLAST File
    #--------------------------------------
    counter = 0
    total_counter = 0
    lines = IN.readlines()
    dedupe = {}
    for line in lines:
        if not re.search("^#",line):
            total_counter+=1
            fields = re.split("\t",line)
            if len(fields) >= 12:
                this_id = fields[0]
                percent_identity = fields[2]
                evalue = fields[10]
            else:
                print("Got a bad line:"+str(line))
                sys.exit(2)               
            if this_id in dedupe:
                oldvalue = re.split("\t",dedupe[this_id])
                list_percent_identity = oldvalue[2]
                list_evalue = oldvalue[10]
                if float(evalue) < float(list_evalue):
                    dedupe[this_id] = line
                elif float(evalue) == float(list_evalue):
                    if percent_identity > list_percent_identity:
                        dedupe[this_id] = line
            else:
              dedupe[this_id] = line
              counter+=1
    print("Total # records = "+str(total_counter)+"\nBest only # records = "+str(counter))
    print("Writing to output file...")
    for key,value in sorted(dedupe.items()):
        OUT.write(str(value))

    #-----------
    #Close Files
    #-----------
    IN.close()
    OUT.close()

if __name__ == "__main__":
    main(sys.argv[1:])

