#!/usr/bin/python3
################################
### Cellular Automata Viewer ###
### Jennifer Meneghin        ###
### 02/18/2020               ###
################################

import sys, getopt
import pandas as pd
import matplotlib.pyplot as plt

def usage ():
    usage = "\nCellular Autotmata Viewer\n"
    usage = usage + "\nUsage: ca_viewer.py -n num_gens -c CA\n"
    usage = usage + "\nThis program takes a 128 bit CA, and requests user input for bit strings.\n"
    usage = usage + "It then runs the CA on the bit string for num_gens generations.\n\n"
    usage = usage + "The CA must be 128 bits long exactly.\nThe input bit strings must be at least 7 bits long.\n"
    usage = usage + "If num_gens is not provided the CA will run for 100 generations.\n"
    usage = usage + "At any prompt, <enter> or Q<enter> will quit the program.\n\n"
    usage = usage + "Jennifer Meneghin\n"
    usage = usage + "February 18, 2020\n\n"
    return usage

def main(argv):
    #---------------------------
    #Read command line arguments
    #---------------------------
    ca_string = "00000000001111111111000000000011111111110000000000111111111100000000001111111111000000000011111111110000000000111111111111111111"
    num_gens = 100
    try:
        opts, args = getopt.getopt(argv,"hi:n:",["istring=","nint="])
    except getopt.GetoptError:
        print("\nNot a valid argument or value")
        print(usage())
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print(usage())
            sys.exit()
        elif opt in ("-i", "--istring"):
            ca_string = arg
        elif opt in ("-n", "--nint"):
            num_gens = arg
    print("Running CA = "+ca_string)
    print("for "+str(num_gens)+" generations...")

    #-------------------------------------
    #Keep going until user decides to quit
    #-------------------------------------
    input_string = ""
    while not(input_string == "Q"):

        #------------------------------------
        # Get Bit String at least 7 bits long
        #------------------------------------
        input_string = input("Please enter a bit string at least 7 bits long or Q to quit: ")
        if input_string == "Q":
            print("\n\nThank you for using the CA Viewer.\n\n")
            sys.exit(0)
        if len(input_string) < 7:
            print("\n\nInput String must be at least 7 bits. You only entered "+str(len(input_string))+" bits.\n\n")
            sys.exit(2)
        try:
            checker = int(input_string,2)
        except ValueError:
            print("\n\nInput String must be a bit string -- 0s and 1s only please.\n\n")
            sys.exit(2)

        #--------------------------------------------------------------------------------------------------------
        #This program runs an 128 bit CA ... therefore neighborhoods are size 7 (because 2**7=128)
        #So, tack last three bits to beginning and first three bits to end .. then start at [3] and go to [len-3]
        #--------------------------------------------------------------------------------------------------------
        total_results = []
        this_result = []
        for i in input_string:
            this_result.append(int(i))
        total_results.append(this_result)
        lis = len(input_string)
        current_string = input_string[lis-3:lis] + input_string + input_string[0:3]
        for i in range(1,int(num_gens)):
            j = 3
            end = len(current_string)-3
            this_result = []               #as array of integers for plotting 2D grid
            temp_string = ""               #as string for converting from binary to int
            while j < end:
                neighborhood = str(current_string[j-3:j+4])
                ca_index = int(neighborhood,2)                       #convert from binary string to int ... nifty
                this_result.append(int(ca_string[ca_index]))         #add ca result to array
                temp_string = temp_string + str(ca_string[ca_index]) #add ca result to string
                j+=1
            lis = len(temp_string)
            current_string = temp_string[lis-3:lis] + temp_string + temp_string[0:3] #new string for next generation
            total_results.append(this_result)                                        #new row in grid of integers for display

        #----------
        # Plot Grid
        #----------
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ndf = pd.DataFrame(total_results, dtype=int)
        ax.imshow(ndf)
        ax.set_title('Results of 2D CA for '+str(num_gens)+' generations')
        plt.show()

        #-----------
        # Save Plot?
        #-----------
        value = input("Would you like to save this file (Y/N)? ")
        if value == "Y" or value == "y" or value == "Yes" or value == "yes":
            png_file = input("Please enter a file name to save to (.png will be appended to it): ")
            png_file = png_file + ".png"
            fig.savefig(png_file, dpi=300, format='png')

if __name__ == "__main__":
    main(sys.argv[1:])
