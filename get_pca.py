#!/usr/bin/python3
#########################################################################
### Get kmer Frequencies                                              ###
### Usage: get_pca.py -i <tab delim file>                             ###
### This program takes a delimited file as its only parameter, with:  ###
### Column 1 = ID, rest of columns are data points                    ###
###                                                                   ###
### It returns a PCA plot of the two principle components of the data ###
###                                                                   ###
### Jennifer Meneghin                                                 ###
### 02/13/2020                                                        ###
#########################################################################

import sys, getopt
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.decomposition import PCA

def usage ():
    usage = "\nGet PCA\n"
    usage = usage + "\nUsage: get_pca.py -i <tab delim file> -o <optional output filename>\n"
    usage = usage + "\nThis program reads a tab delimited file where,\n"
    usage = usage + "Column 1 = Group, and the rest of the columns are data points.\nFirst row should contain headers (not data).\n"
    usage = usage + "Colors will be duplicated if there are more than eight groups.\n\n"
    usage = usage + "-o out_file_name is optional (default = pca_image), .png will be appended to the name\n\n"
    usage = usage + "It returns a PCA plot (in .png file) of the first two principle components of the data\n\n"
    usage = usage + "Jennifer Meneghin\n"
    usage = usage + "February 13, 2020\n\n"
    return usage

def main(argv):
    #---------------------------
    #Read command line arguments
    #---------------------------
    in_file = ""
    out_file = "pca_image"
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print("\nNot a valid argument or value")
        print(usage())
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print(usage())
            sys.exit()
        elif opt in ("-i", "--ifile"):
            in_file = arg
        elif opt in ("-o", "--ofile"):
            out_file = arg

    #---------------------
    #Open File for reading
    #---------------------
    try:
        df = pd.read_csv(in_file,sep='\t')
        print("\nParameters:\ndata file = "+in_file)
        print("output file = "+out_file+".png")
    except FileNotFoundError:
        print("\nPlease enter a tab delimited file.")
        print(usage())
        sys.exit(2)

    #--------------------------------------------------------------------
    #Get the First Column from the File and Put in Dictionary with Colors
    #--------------------------------------------------------------------
    colors = ['red','green','orange','blue','yellow','purple','pink','turquoise']
    groups = {}
    bigarray = df.to_numpy()
    i = 0
    for item in bigarray[:,0]:
        if not(item in groups):
            groups[item] = colors[i]
            i+=1
            if i >= len(colors):
                print("Warning: More than "+str(len(colors))+" groups -- colors will be duplicated")
                i = 0
    print("Number of groups found = "+str(len(groups)))

    #-------------
    #Calculate PCA
    #-------------
    pca = PCA(n_components=2)
    bignums = bigarray[:,1:]
    pca.fit(bignums)
    bnt = pca.transform(bignums)
    df_pca = pd.DataFrame(bnt.T, index=['PC1','PC2'])
    
    #----------------------------------
    #Display PCA (first two components)
    #----------------------------------
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    my_patches = []
    i = 0
    while i < len(df_pca.T):
        ax1.scatter(df_pca[i].loc['PC1'],df_pca[i].loc['PC2'],c=groups[bigarray[i,0]])
        i +=1
    for item in groups:
        my_patches.append(mpatches.Patch(color=groups[item], label=item))
    lgd = ax1.legend(handles=my_patches,loc='center right', bbox_to_anchor=(1.5,0.5))
    pcatotal = PCA().fit(bignums)
    pepc1 = pcatotal.explained_variance_ratio_[0]
    pepc2 = pcatotal.explained_variance_ratio_[1]
    petotal = (pepc1 + pepc2)*100
    ax1.set_title('Two Dimensional PCA of '+in_file+'\n(Variance Explained = '+str("%.2f" % petotal)+"%)")
    ax1.set_xlabel("PC 1")
    ax1.set_ylabel("PC 2")
    
    #------------------
    #Write Plot to File
    #------------------
    fig.savefig(out_file+".png", dpi=300, format='png', bbox_extra_artists=[lgd], bbox_inches='tight')

if __name__ == "__main__":
    main(sys.argv[1:])
