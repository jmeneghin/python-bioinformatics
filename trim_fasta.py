#!/usr/bin/python3
#########################################
### Jennifer Meneghin                 ###
### Original Perl Version: 08/10/2010 ###
### Python Version: 01/25/2020        ###
#########################################

import re, sys, getopt

def usage ():
    usage = "Usage: trim_fasta.pl -i fasta_file -m integer -p integer -c 1\n\n"
    usage = usage + "Parameters:\n"
    usage = usage + "-i fasta_file:\tThe fasta file to trim.\n"
    usage = usage + "-m integer:\tThe minimum length allowed for a fasta record. (Optional. Default is 50.)\n"
    usage = usage + "-p integer:\tThe minimum length of a Poly(A/T) sequence to consider. (Optional. Default is 8.)\n"
    usage = usage + "-c 1\t\t1 = Remove sequences with Poly(A/T) in the middle of the sequence, 0 = Don't remove.\n"
    usage = usage + "\t\t(Optional. Default is 0 = Don't remove.)\n\n"
    usage = usage + "This script trims any Poly(A/T) tails from the beginning and end of each sequence,\n"
    usage = usage + "then removes the sequence if it's length is less than the minimum number provided (-m integer).\n"
    usage = usage + "If -c is set to 1, it also removes any sequences with a Poly(A/T) sequence in the middle of the sequence (assumed to be chimera).\n\n"
    usage = usage + "A Poly(A/T) sequence is defined as a Poly(A) or a Poly(T) sequence.\n"
    usage = usage + "A Poly(A) sequence is defined as (-p integer) or more consecutive As, Xs or Ns.\n"
    usage = usage + "A Poly(T) sequence is defined as (-p integer) or more consecutive Ts, Xs or Ns.\n"
    usage = usage + "Please run homopolymer_count.pl to help you decide the most appropriate length of Poly(A/T) sequence for your dataset.\n\n"
    usage = usage + "It returns two new fasta files: one with the removed sequences (fasta_file.removed), and one with the salvaged sequences (fasta_file.salvaged).\n\n"
    usage = usage + "Jennifer Meneghin\n"
    usage = usage + "01/27/2020\n"
    return usage

def trim_it(header, string, stringa, stringt):
    if re.search("^"+stringa, string):
        length = len(string)
        count = 0
        for i in string:
            if not(i == "A") and not(i == "N") and not(i == "X"):
                string = string[count:length]
                break
            count+=1
    if re.search("^"+stringt, string):
        length = len(string)
        count = 0
        for i in string:
            if not(i == "T") and not(i == "N") and not(i == "X"):
                string = string[count:length]
                break
            count+=1
    if re.search(stringa+"$",string):
        length = len(string)
        while length > 0:
            if not(string[length-1] == "A") and not(string[length-1] == "N") and not(string[length-1] == "X"):
                string = string[0:length]
                break
            length=length-1
    if re.search(stringt+"$",string):
        length = len(string)
        while length > 0:
            if not(string[length-1] == "T") and not(string[length-1] == "N") and not(string[length-1] == "X"):
                string = string[0:length]
                break
            length=length-1
    return string

def main(argv):

    #---------------------------
    #Read command line arguments
    #---------------------------
    in_file = ""
    minimum_sequence_length = 50
    poly_length = 8
    remove_chimera_flag = 0
    try:
        opts, args = getopt.getopt(argv,"hi:m:p:c:",["ifile=","mint=","pint=","cint="])
    except getopt.GetoptError:
        print("\nNot a valid argument or value")
        print("trim_fasta.pl -i <fasta input file> -m <minimum sequence length -p <minimum polyt(A/T) sequence length> -c <chimera flag 1=yes or 0=no>\n")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print(usage())
            print("trim_fasta.pl -i <fasta input file> -m <minimum sequence length -p <minimum polyt(A/T) sequence length> -c <chimera flag 1=yes or 0=no>\n")
            sys.exit()
        elif opt in ("-i", "--ifile"):
            in_file = arg
        elif opt in ("-m", "--mint"):
            minimum_sequence_length = arg
        elif opt in ("-p", "--pint"):
            poly_length = arg
        elif opt in ("-c", "--cint"):
            remove_chimera_flag = arg

    remove_chimera_flag = int(remove_chimera_flag)
    if not remove_chimera_flag == 0:
        if not remove_chimera_flag == 1:
            #print("bad flag")
            print("\nNot a valid chimera flag: -c must be 1 or 0. You entered: "+str(remove_chimera_flag))
            print("trim_fasta.pl -i <fasta input file> -m <minimum sequence length -p <minimum polyt(A/T) sequence length> -c <chimera flag 1=yes or 0=no>\n")
            sys.exit(2)
        
    #----------------------------------
    #Open Files for reading and writing
    #----------------------------------
    try:
        IN = open(in_file,"r")
    except FileNotFoundError:
        print("\nFASTA Input File could not be found")
        print("trim_fasta.pl -i <fasta input file> -m <minimum sequence length -p <minimum polyt(A/T) sequence length> -c <chimera flag 1=yes or 0=no>\n")
        sys.exit(2)

    salvaged_file = in_file + ".salvaged"
    removed_file = in_file + ".removed"
    SALVAGED = open(salvaged_file,"w")
    REMOVED = open(removed_file,"w")
    print("\nParameters:\nfasta file = "+in_file)
    print("minimum sequence length = "+str(minimum_sequence_length))
    print("minimum Poly(A/T) length = "+str(poly_length))
    print("remove chimera flag = "+str(remove_chimera_flag)+"\n")


    #--------------
    #The main event
    #--------------
    header = ""
    string = ""
    stringa = "(A|X|N){"+str(poly_length)+"}"
    stringt = "(T|X|N){"+str(poly_length)+"}"
    count_removed = 0
    count_chimera = 0
    count_too_short = 0
    count_salvaged = 0
    count_total = 0
    count_trimmed = 0

    lines = IN.readlines()
    for line in lines:
        if re.search("^>",line):
            if len(header)>0 and len(string)>0:
                count_total+=1
                if re.search(stringa,string) or re.search(stringt,string):
                    newstring = trim_it(header, string, stringa, stringt)
                    if len(newstring) < len(string):
                        count_trimmed += 1
                        string = newstring
                    if len(string) < int(minimum_sequence_length):
                        REMOVED.write(header)
                        REMOVED.write(string+"\n")
                        print("Removed too short: "+str(header))
                        count_removed+=1
                        count_too_short+=1
                    elif re.search(stringa,string) or re.search(stringt,string):
                        if remove_chimera_flag == 1:
                            REMOVED.write(header)
                            REMOVED.write(string+"\n")
                            print("Removed as chimera (Poly(A/T) in middle): "+str(header))
                            count_removed+=1
                            count_chimera+=1
                        else:
                            SALVAGED.write(header)
                            SALVAGED.write(string+"\n")
                            count_salvaged+=1
                            count_chimera+=1
                    else:
                        SALVAGED.write(header)
                        SALVAGED.write(string+"\n")
                        count_salvaged+=1
                else:
                    SALVAGED.write(header)
                    SALVAGED.write(string+"\n")
                    count_salvaged+=1
            string = ""
            header = line
        else:
            line = line.rstrip()
            string = string + line
    #ENDFOR
    if len(header)>0 and len(string)>0:
        count_total+=1
        if re.search(stringa,line) or re.search(stringt,line):
            newstring = trim_it(header, string, stringa, stringt)
            if len(newstring) < len(string):
                count_trimmed += 1
                string = newstring
            if len(string) < int(minimum_sequence_length):
                REMOVED.write(header)
                REMOVED.write(string+"\n")
                print("Removed too short: "+str(header))
                count_removed+=1
                count_too_short+=1
            elif re.search(stringa,string) or re.search(stringt,string):
                if remove_chimera_flag == 1:
                    REMOVED.write(header)
                    REMOVED.write(string+"\n")
                    print("Removed as chimera (Poly(A/T) in middle): "+str(header))
                    count_removed+=1
                    count_chimera+=1
                else:
                    SALVAGED.write(header)
                    SALVAGED.write(string+"\n")
                    count_salvaged+=1
                    count_chimera+=1
            else:
                SALVAGED.write(header)
                SALVAGED.write(string+"\n")
                count_salvaged+=1
        else:
            SALVAGED.write(header)
            SALVAGED.write(string+"\n")
            count_salvaged+=1

    IN.close()
    SALVAGED.close()
    REMOVED.close()

    print("REMOVED = "+str(count_removed))
    print("REMOVED TOO SHORT = "+str(count_too_short))
    if remove_chimera_flag == 1:
        print("REMOVED CHIMERA = "+str(count_chimera))
    else:
        print("FOUND CHIMERA (not removed) = "+str(count_chimera))
    print("SALVAGED = "+str(count_salvaged))
    print("TRIMMED = "+str(count_trimmed))
    print("TOTAL = "+str(count_total))


    
if __name__ == "__main__":
    main(sys.argv[1:])
