#!/usr/bin/python3
#######################################################################################
### This script reads a fasta file and returns a fasta file containing the records  ###
### found with more than the minimum number of consecutive nucleotides, along with  ###
### nucleotide and count found (e.g. A=24, C=8).                                    ###
###                                                                                 ###
### Additionally, this script returns the number of fasta entries with at least N   ###
### consecutive Zs, where N = is the range from the integer provided to the point   ###
### that there are no more homopolymers, and Z = A, G, C, T, X, and N.              ###
###                                                                                 ###
### It also returns the number of entries in the fasta file, and the average        ###
### sequence length. The output is a tab delimmited matrix written to standard out. ###
### To write this output to a file, run as follows:                                 ###
### homopolymer_count.pl -i fasta_file -m an_integer >my_file.txt                   ###
###                                                                                 ###                
### Usage: homopolymer_count.pl -i <fasta input file> -m N -o <fasta output file>   ###
###                                                                                 ###
### Jennifer Meneghin                                                               ###
### Original Perl version: 02/04/2009                                               ###
### Python version: 01/21/2020                                                      ###
#######################################################################################

import re, sys, getopt

############
# Usage (-h)
############
def usage():
    usage = "\nThis script reads a fasta file and returns a fasta file containing the records\n"
    usage = usage + "found with more than the minimum number of consecutive nucleotides, along with\n"
    usage = usage + "nucleotide and count found (e.g. A=24, C=8).\n\n"
    usage = usage + "Additionally, this script returns the number of fasta entries with at least N\n"
    usage = usage + "consecutive Zs, where N = is the range from the integer provided to the point\n"
    usage = usage + "that there are no more homopolymers, and Z = A, G, C, T, X, and N.\n\n"
    usage = usage + "It also returns the number of entries in the fasta file, and the average\n"
    usage = usage + "sequence length. The output is a tab delimmited matrix written to standard out.\n"
    usage = usage + "To write this output to a file, run as follows:\n"
    usage = usage + "homopolymer_count.pl -i fasta_file -m an_integer >my_file.txt\n\n"
    usage = usage + "Usage: homopolymer_count.pl -i <fasta input file> -m N -o <fasta output file>\n\n"
    usage = usage + "Jennifer Meneghin\n"
    usage = usage + "Original Perl version: 02/04/2009\n"
    usage = usage + "Python version: 01/24/2020\n"
    return usage

###############
#The Main Event
###############
def main(argv):
    
    #---------------------------
    #Read command line arguments
    #---------------------------
    in_file = ""
    minimum = 5
    out_file = "homopolymers.fasta"
    try:
        opts, args = getopt.getopt(argv,"hi:o:m:",["ifile=","ofile=","mint="])
    except getopt.GetoptError:
        print("\nNot a valid argument or value")
        print("homopolymer_count.pl -i <fasta input file> -m N -o <fasta output file>\n")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print(usage())
            print("homopolymer_count.pl -i <fasta input file> -m N -o <fasta output file>\n")
            sys.exit()
        elif opt in ("-i", "--ifile"):
            in_file = arg
        elif opt in ("-o", "--ofile"):
            out_file = arg
        elif opt in ("-m", "--mint"):
            minimum = arg

    #----------------------------------
    #Open Files for reading and writing
    #----------------------------------
    try:
        IN = open(in_file,"r")
    except FileNotFoundError:
        print("\nFASTA Input File could not be found")
        print("homopolymer_count.pl -i <fasta input file> -m N -o <fasta output file>\n")
        sys.exit(2)
    OUT = open(out_file,"w")
    print("\nParameters:\ninput file = "+in_file+"\noutput file = "+out_file+"\nMinimum Number to Count = "+str(minimum))

    #-----------------------------------------
    #Read the Input FASTA and Print The Matrix
    #-----------------------------------------
    headers = {}
    records = {}
    homopolymer_count = 9999
    i = minimum
    print("Poly Seq Length\tPolyA\tPolyT\tPolyG\tPolyC\tPolyN\tPolyX")
    lines = IN.readlines()
    while homopolymer_count > 0:
        total_seq_length = 0
        count = 0
        counta = 0
        countt = 0
        countg = 0
        countc = 0
        countx = 0
        countn = 0
        seq = ""
        header = ""
        for line in lines:
            if re.search("^>",line):
                if len(header)>0:
                    count+=1
                    records[header] = seq
                    seqa = "A{"+str(i)+"}"              
                    seqt = "T{"+str(i)+"}"              
                    seqg = "G{"+str(i)+"}"
                    seqc = "C{"+str(i)+"}"
                    seqx = "X{"+str(i)+"}"
                    seqn = "N{"+str(i)+"}"
                    if re.search(seqa, seq):
                        counta+=1
                        if header in headers:
                            headers[header] = headers[header] + "\tA=" + str(i)
                        else:
                            headers[header] = "A=" + str(i)
                    if re.search(seqa, seq):
                        countt+=1
                        if header in headers:
                            headers[header] = headers[header] + "\tT=" + str(i)
                        else:
                            headers[header] = "T=" + str(i)
                    if re.search(seqg, seq):
                        countg+=1
                        if header in headers:
                            headers[header] = headers[header] + "\tG=" + str(i)
                        else:
                            headers[header] = "G=" + str(i)
                    if re.search(seqc, seq):
                        countc+=1
                        if header in headers:
                            headers[header] = headers[header] + "\tC=" + str(i)
                        else:
                            headers[header] = "C=" + str(i)
                    if re.search(seqx, seq):
                        countx+=1
                        if header in headers:
                            headers[header] = headers[header] + "\tX=" + str(i)
                        else:
                            headers[header] = "X=" + str(i)
                    if re.search(seqn, seq):
                        countn+=1
                        if header in headers:
                            headers[header] = headers[header] + "\tN=" + str(i)
                        else:
                            headers[header] = "N=" + str(i)
                    homopolymer_count = counta + countt + countg + countc + countx + countn
                    total_seq_length = total_seq_length + len(seq)
                #ENDIF
                header = line
                header = header.rstrip()
                len_header = len(header)
                header = header[1:len_header] #removing  the ">" from the front of the header line
                seq = ""
            else:
                seq = seq + line
                seq = seq.rstrip()
                seq = seq.upper()
        #ENDFOR
        count+=1
        records[header] = seq
        seqa = "A{"+str(i)+"}"
        seqt = "T{"+str(i)+"}"
        seqg = "G{"+str(i)+"}"
        seqc = "C{"+str(i)+"}"
        seqx = "X{"+str(i)+"}"
        seqn = "N{"+str(i)+"}"
        if re.search(seqa, seq):
            counta+=1
            if header in headers:
                headers[header] = headers[header] + "\tA=" + str(i)
            else:
                headers[header] = "A=" + str(i)
        if re.search(seqt, seq):
            countt+=1
            if header in headers:
                headers[header] = headers[header] + "\tT=" + str(i)
            else:
                headers[header] = "T=" + str(i)
        if re.search(seqg, seq):
            countg+=1
            if header in headers:
                headers[header] = headers[header] + "\tG=" + str(i)
            else:
                headers[header] = "G=" + str(i)
        if re.search(seqc, seq):
            countc+=1
            if header in headers:
                headers[header] = headers[header] + "\tC=" + str(i)
            else:
                headers[header] = "C=" + str(i)
        if re.search(seqx, seq):
            countx+=1
            if header in headers:
                headers[header] = headers[header] + "\tX=" + str(i)
            else:
                headers[header] = "X=" + str(i)
        if re.search(seqn, seq):
            countn+=1
            if header in headers:
                headers[header] = headers[header] + "\tN=" + str(i)
            else:
                headers[header] = "N=" + str(i)
        homopolymer_count = counta + countt + countg + countc + countx + countn
        total_seq_length = total_seq_length + len(seq)
        average_seq_length = total_seq_length / count
        total_seq_length = 0
        print(str(i)+"\t"+str(counta)+"\t"+str(countt)+"\t"+str(countg)+"\t"+str(countc)+"\t"+str(countn)+"\t"+str(countx)) #runs once per homopolymer length
        i=int(i)+1
    #ENDWHILE
    print("Number of Fasta Entries = " +str (count))
    print("Average Sequence Length = " +str (average_seq_length))

    #-------------------------------------
    #Write the Results to the Output FASTA
    #-------------------------------------
    #for key, value in sorted(headers.items()):
    for key, value in sorted(records.items()):
        key = key.rstrip()
        if key in headers:
            if re.search("\t",headers[key]):
                acount = 0
                tcount = 0
                gcount = 0
                ccount = 0
                xcount = 0
                ncount = 0
                polycounts = re.split("\t", headers[key])
                for j in polycounts:
                    parts = re.split("=", j)
                    if (parts[0] == "A") and (int(parts[1]) > int(acount)):
                        acount = parts[1]
                    elif (parts[0] == "T") and (int(parts[1]) > int(tcount)):
                        tcount = parts[1]
                    elif (parts[0] == "G") and (int(parts[1]) > int(gcount)):
                        gcount = parts[1]
                    elif (parts[0] == "C") and (int(parts[1]) > int(ccount)):
                        ccount = parts[1]
                    elif (parts[0] == "X") and (int(parts[1]) > int(xcount)):
                        xcount = parts[1]
                    elif (parts[0] == "N") and (int(parts[1]) > int(ncount)):
                        ncount = parts[1]
                    else:
                        print("ERROR!")
                        print(usage())
                        sys.exit(2)
                OUT.write(">"+str(key)+" ")
                if int(acount) > 0:
                    OUT.write("A="+str(acount))
                    if (int(tcount) > 0 or int(gcount) > 0 or int(ccount) > 0 or int(xcount) > 0 or int(ncount) > 0):
                        OUT.write(", ")
                if int(tcount) > 0:
                    OUT.write("T="+str(tcount))
                    if (int(gcount) > 0 or int(ccount) > 0 or int(xcount) > 0 or int(ncount) > 0):
                        OUT.write(", ")
                if int(gcount) > 0:
                    OUT.write("G="+str(gcount))
                    if (int(ccount) > 0 or int(xcount) > 0 or int(ncount) > 0):
                        OUT.write(", ")
                if int(ccount) > 0:
                    OUT.write("C="+str(ccount))
                    if (int(xcount) > 0 or int(ncount) > 0):
                        OUT.write(", ")
                if int(xcount) > 0:
                    OUT.write("X="+str(xcount))
                    if int(ncount) > 0:
                        OUT.write(", ")
                if int(ncount) > 0:
                    OUT.write("N="+str(ncount))
                OUT.write("\n")
                OUT.write(records[key]+"\n")
            else:
                OUT.write(">"+str(key)+" "+str(headers[key])+"\n")
                OUT.write(value+"\n")
            #ENDIF
        else:
            OUT.write(">"+str(key)+"\n")
            OUT.write(value+"\n")
            #ENDIF
    #ENDFOR
    #-----------
    #Close Files
    #-----------
    IN.close()
    OUT.close()

if __name__ == "__main__":
    main(sys.argv[1:])

