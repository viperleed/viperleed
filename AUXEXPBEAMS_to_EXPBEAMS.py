# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 17:45 2019

@author: Florian Kraushofer

Reads an AUXEXPBEAMS file and writes the contents in EXPBEAMS.csv format.
"""

import tleedmlib as tl


###############################################
#                  MAIN                       #
###############################################

def main():
    # print some info
    print("This utility reads an AUXEXPBEAMS file (ie, beams formatted as "
          "TensErLEED experimental input) and writes the contents in the csv "
          "formatting applied in tleedmap for THEOBEAMS.csv and EXPBEAMS.csv "
          "files.\n")
    
    # read the AUXEXPBEAMS file  
    filename = ""
    while filename == "":
        filename = input("Enter AUXEXPBEAMS file name (default:"
                         " 'AUXEXPBEAMS'): ")
        if filename == "":
            filename = "AUXEXPBEAMS"
        try:
            beams = tl.readAUXEXPBEAMS(filename, interactive=True)
        except FileNotFoundError:
            print("File "+filename+" not found.")
            filename = ""
        except:
            print("Exception while reading AUXEXPBEAMS file")
            return 1
    
    if len(beams) == 0:
        print("Error reading AUXEXPBEAMS file.")
        return 1
    # print some info
    print("Found file with "+str(len(beams))+" beams.\n")
    
    # relabel beams
    mw = max([beam.lwidth for beam in beams])
    for beam in beams:
        beam.lwidth = mw
        beam.getLabel()
    
    # get output file name
    filename = ""
    filename = input("Enter output file name (default: 'EXPBEAMS.csv'): ")
    if filename == "":
        filename = "EXPBEAMS.csv"
    
    # write new file
    try:
        tl.writeOUTBEAMS(beams, filename)
        print("Wrote output as "+filename)
    except:
        print("Error writing new file "+filename)
        return 1
    return 0
        

if __name__ == "__main__":
    main()
