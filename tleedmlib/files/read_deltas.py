
from re import I
import sys
import numpy as np
import fortranformat as ff
import matplotlib.pyplot as plt
import scipy
import os

def read_delta_file(filename, n_E):
    #TODO:
    # - get rid of HeaderBlocks ang give them propper names
    # - better readability!
    # - can we improve performance?
    # - comments in Englisch
    # - comments that explain why things are done a certain way
        
    #Daten- und Hilfslisten des jeweiligen Files
    HeaderBlock1=[]
    HeaderBlock2=[]
    HeaderBlock3=[]
    HeaderBlock4=[]
    HeaderBlock5=[]
    HeaderBlock6=[]
    Parameter_list=[]
    E_array= np.full(n_E, fill_value = np.NaN)
    VPI_array= np.full(n_E, fill_value = np.NaN)
    VV_array= np.full(n_E, fill_value = np.NaN)
    listdummy=[]
    listdummy2=[]

    # we need two fortran format readers
    ff_reader_6E13_7 = ff.FortranRecordReader("6E13.7")
    ff_reader_10F7_4 = ff.FortranRecordReader("10F7.4")

    #Einlesen der Daten eines Files
    with open(filename, mode = "r") as file:
        content = file.readlines()
    # make into an iterator
    file_lines = iter(content)

    #1.Block of Header - only 1 line
    line=next(file_lines)
    HeaderBlock1=ff_reader_6E13_7.read(line)
    trar1=[0,0]
    trar2=[0,0]
    theta, phi, trar1[0], trar1[1], trar2[0], trar2[1] = HeaderBlock1
    trar1 = np.array(trar1)
    trar2 = np.array(trar2)

    #2.Block of Header - also only 1 line
    line=next(file_lines)
    z2_elemente=line.split()
    for part in z2_elemente:
        HeaderBlock2.append(int(part))
    Int0, n_atoms, NCSteps = HeaderBlock2

    #3.Block of Header
    while(len(HeaderBlock3)<Int0):
        line=next(file_lines)
        Delta_Amp=line.split()
        for part in Delta_Amp:
            listdummy.append(part)
        # AMI: this is unreadable, make nicer
        if(len(listdummy)==2*Int0):
            for j in range(2*Int0):
                if(j%2==0):
                    e1=float(listdummy[j])
                elif(j%2==1):
                    e2=float(listdummy[j])
                    Element=[e1,e2]
                    HeaderBlock3.append(Element)
            listdummy.clear()

    #4.Block of Header
    while(len(HeaderBlock4)<n_atoms):
        line=next(file_lines)
        Cundisp=line.split()
        for part in Cundisp:
            listdummy.append(part)
        # same as above
        if(len(listdummy)>=3*n_atoms):
            for j in range(3*n_atoms):
                if(j%3==0):
                    e1=float(listdummy[j])
                elif(j%3==1):
                    e2=float(listdummy[j])
                elif(j%3==2):
                    e3=float(listdummy[j])
                    Element=[e1,e2,e3]
                    HeaderBlock4.append(Element)
            listdummy.clear()

    #5.Block of Header
    while(len(HeaderBlock5)<1):
        line=next(file_lines)
        listdummy=ff_reader_10F7_4.read(line)
        for part in listdummy:
            listdummy2.append(part)
        if(len(listdummy2)>=3*n_atoms*NCSteps):
        # d=1
            #zeile6block=i
            matrix=np.zeros([NCSteps,n_atoms,3])
            for j in range(NCSteps):
                for k in range(n_atoms):
                    for l in range(3):
                        matrix[j,k,l]=listdummy2[n_atoms*3*j+3*k+l]
            HeaderBlock5.append(matrix)
            listdummy.clear()
            listdummy2.clear()

    #6.Block of Header
    while(len(HeaderBlock6)<NCSteps):
        line=next(file_lines)
        listdummy=ff_reader_10F7_4.read(line)
        for part in listdummy:
            listdummy2.append(part)
        if(len(listdummy2)>=NCSteps):
            for j in range(NCSteps):
                HeaderBlock6.append(listdummy2[j])
            listdummy.clear()
            listdummy2.clear()

    # Initialize arrays for reference and delta amplitudes
    amplitudes_ref = (np.full(shape = [n_E, Int0], 
        fill_value = np.nan, dtype = complex))
    amplitudes_del = (np.full(shape = [n_E, NCSteps, Int0], 
        fill_value = np.nan, dtype = complex))

    # maybe working arrays for transfer into amplitude arrays ?

    #End of the Header - Start of Reading in the Delta Data
    for e_index in range(n_E): # Energy loop

        # Energy, VPI and VV header
        line=next(file_lines)
        listdummy=ff_reader_6E13_7.read(line)
        for part in listdummy:
            if(part is not None):
                Parameter_list.append(part)
        E, VPI, VV = Parameter_list
        E_array[e_index] = E
        VPI_array[e_index] = VPI
        VV_array[e_index] = VV
        listdummy.clear()
        Parameter_list.clear()

        # @Tobi: these conditions below are not nice.
        # think about:
        # do you need to have the if inside the loop?
        # do you need the if at all?
        #
        # for performance and readability:
        # move loops over j  / j&k into an extra function

        # Reference amplitudes
        while(len(listdummy2)<2*Int0):
            line=next(file_lines)
            listdummy=ff_reader_6E13_7.read(line)
            for part in listdummy:
                listdummy2.append(part)
            if(len(listdummy2)>=2*Int0):
                for j in range(Int0):
                    amplitudes_ref[e_index, j] = complex(listdummy2[2*j], listdummy2[2*j+1])
        listdummy.clear()
        listdummy2.clear()

        # Delta amplitudes
        while(len(listdummy2)<2*Int0*NCSteps):
            line=next(file_lines)
            listdummy=ff_reader_6E13_7.read(line)
            for part in listdummy:
                listdummy2.append(part)
            if(len(listdummy2)>=2*Int0*NCSteps):
                for j in range(NCSteps):
                    for k in range(Int0):
                        amplitudes_del[e_index, j, k] = complex(listdummy2[2*Int0*j+k*2], listdummy2[2*Int0*j+k*2+1]) 
        listdummy.clear()
        listdummy2.clear()

    # think about what is sensible to return here
    return ((phi, theta),
        (trar1, trar2),
        (Int0, n_atoms, NCSteps),
        (E_array, VPI_array, VV_array),
        HeaderBlock3,
        HeaderBlock4,
        HeaderBlock5,
        HeaderBlock6,
        amplitudes_ref,
        amplitudes_del)
