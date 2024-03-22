"""ViPErLEED utility to convert AUXEXPBEAMS to EXPBEAMS.

Reads an AUXEXPBEAMS file and writes the contents in EXPBEAMS.csv format.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2019-12-16'
__license__ = 'GPLv3+'

from viperleed.calc.files.beams import readAUXEXPBEAMS, writeOUTBEAMS


def main():
    # print some info
    print("This utility reads an AUXEXPBEAMS file (i.e., beams formatted as "
          "TensErLEED experimental input) and writes the contents in the "
          "csv formatting applied in viperleed.calc for THEOBEAMS.csv and "
          "EXPBEAMS.csv files.\n")

    # read the AUXEXPBEAMS file
    filename = ""
    while filename == "":
        filename = input("Enter AUXEXPBEAMS file name (default:"
                         " 'AUXEXPBEAMS'): ")
        if filename == "":
            filename = "AUXEXPBEAMS"
        try:
            beams = readAUXEXPBEAMS(filename, interactive=True)
        except FileNotFoundError:
            print("File "+filename+" not found.")
            filename = ""
        except Exception:
            print("Exception while reading AUXEXPBEAMS file")
            return 1

    if len(beams) == 0:
        print("Error reading AUXEXPBEAMS file.")
        return 1
    # print some info
    print("Found file with "+str(len(beams))+" beams.\n")

    # relabel beams
    mw = max([beam.getLabel()[1] for beam in beams])
    for beam in beams:
        beam.label = beam.getLabel(lwidth=mw)[0]

    # get output file name
    filename = ""
    filename = input("Enter output file name (default: 'EXPBEAMS.csv'): ")
    if filename == "":
        filename = "EXPBEAMS.csv"

    # write new file
    try:
        writeOUTBEAMS(beams, filename)
        print("Wrote output as "+filename)
    except Exception:
        print("Error writing new file "+filename)
        return 1
    return 0


if __name__ == "__main__":
    main()
