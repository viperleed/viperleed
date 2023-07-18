#!/usr/bin/env python
"""
Quickly checks if atomic density files produce correct integrated charge, since EEASISSS was raising a suspicious error.
"""

import numpy as np
import os

def main():
    filelist = os.listdir('./')
    for filename in filelist:
        if filename == "charge_density_check.py":
            continue
        with open(filename) as file:
            file.readline() #skip 1
            z_line = file.readline().split()
            z = float(z_line[0])

            n = []
            rho = []
            # integrate
            for line in file:
                split = line.split()
                n.append(float(split[0]))
                rho.append(float(split[1]))

        z_integrated = np.trapz(rho, x = n)
        print('Z: '+ str(z) + '\tResult: ' + str(z_integrated) + '\tDiff: ' + str(z-z_integrated)+ '\t\tPercent: ' + str((z_integrated/z -1)*100))


if __name__ == "__main__":
    main()
