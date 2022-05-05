# -*- coding: utf-8 -*-
"""

@author: Tobias Hable, Alexander M. Imre

Reads in delta files.
"""
from re import I
import sys
import numpy as np
import fortranformat as ff
import matplotlib.pyplot as plt
import scipy
import os

def read_delta_file(filename, n_E):
    """This function reads in one file of data and stores the data in arrays, which can be used in later functions (ex.: GetInt)
    
    Parameters
    ----------
    filename : string
    The filename describes which file you want to read in (path from the function location to the file)
    
    n_E : int
    Number of different energy levels. This should be a known factor for your files
    
    
    Returns
    -------
    (phi, theta): tuple of float
    Angles of how the beam hits the sample
    
    (trar1, trar2) : tuple of ndarray
    Vectors of the normal and the reciprocal unit cell of the sample
    
    int0 : int
    Number of beams that are reflected by the sample 
    
    n_atoms: int
    TBD
    
    nc_steps: ndarray
    Number of permutations between direction deltas and vibration deltas
    
    E_array : ndarray
    Array that contains all the energies of the file
    
    VPI_array : ndarray
    Imaginary part of the inner potential of the surface
    
    VV_array : ndarray
    Real part of the inner potential of the surface
    
    Beam_places : ndarray
    Array with the order of beams
    
    Cundisp : ndarray
    TBD always 0
    
    CDisp : ndarry
    Geometric displacement of given delta
    
    amplitudes_ref : ndarray
    Array that contains all values of the reference amplitudes
    
    amplitudes_del : ndarray
    Array that contains all values of the delta amplitudes
    """
        
    #Lists and Arrays needed; E,VPI,VV are numpy arrays that get returned, the others are lists that help saving the other data
    HeaderBlock1=[]
    HeaderBlock2=[]
    E_array= np.full(n_E, fill_value = np.NaN)
    VPI_array= np.full(n_E, fill_value = np.NaN)
    VV_array= np.full(n_E, fill_value = np.NaN)
    listdummy=[]
    listdummy2=[]

    # we need two fortran format readers
    ff_reader_6E13_7 = ff.FortranRecordReader("6E13.7")
    ff_reader_10F7_4 = ff.FortranRecordReader("10F7.4")

    #Reading in the data of a file
    with open(filename, mode = "r") as file:
        content = file.readlines()
    # make into an iterator
    file_lines = iter(content)

    #1.Block of Header - only 1 line - theta, phi, trar1, trar2 variables
    line=next(file_lines)
    HeaderBlock1=ff_reader_6E13_7.read(line)
    trar1=[0,0]
    trar2=[0,0]
    theta, phi, trar1[0], trar1[1], trar2[0], trar2[1] = HeaderBlock1
    trar1 = np.array(trar1)
    trar2 = np.array(trar2)

    #2.Block of Header - also only 1 line - int0, n_atoms, nc_steps variables
    line=next(file_lines)
    z2_elemente=line.split()
    for part in z2_elemente:
        HeaderBlock2.append(int(part))
    int0, n_atoms, nc_steps = HeaderBlock2

    #3.Block of Header - Positions of the beams
    while(len(listdummy2)<2*int0):
        line=next(file_lines)
        listdummy=line.split()
        for part in listdummy:
            listdummy2.append(part)
        #dieses if wieder redundant?
        if(len(listdummy2)>=2*int0):
            Beam_places=(np.full(shape=[int0,2],fill_value=np.NaN))
            for i in range(int0):
                for j in range(2):
                    Beam_places[i,j]=listdummy2[2*i+j]
    listdummy.clear()
    listdummy2.clear()

    #4.Block of Header - Cundisp
    while(len(listdummy2)<3*n_atoms):
        line=next(file_lines)
        listdummy=line.split()
        for part in listdummy:
            listdummy2.append(part)
        #if redundant?
        if(len(listdummy2)>=3*n_atoms):
            Cundisp=(np.full(shape=[n_atoms,3],fill_value=np.NaN))
            for i in range(n_atoms):
                for j in range(3):
                    Cundisp[i,j]=listdummy2[3*i+j]
    listdummy.clear()
    listdummy2.clear()

    #5.Block of Header - CDisp
    while(len(listdummy2)<3*n_atoms*nc_steps):
        line=next(file_lines)
        listdummy=ff_reader_10F7_4.read(line)
        for part in listdummy:
            listdummy2.append(part)
            #maybe getting rid of the if again
        if(len(listdummy2)>=3*n_atoms*nc_steps):
            CDisp=(np.full(shape=[nc_steps,n_atoms,3],fill_value=np.NaN))
            for j in range(nc_steps): #some syntax error here
                for k in range(n_atoms):
                    for l in range(3):
                        CDisp[j,k,l]=listdummy2[n_atoms*3*j+3*k+l]
    listdummy.clear()
    listdummy2.clear()

    #6.Block of Header - Aid
    while(len(listdummy2)<nc_steps):
        line=next(file_lines)
        listdummy=ff_reader_10F7_4.read(line)
        for part in listdummy:
            listdummy2.append(part)
        if(len(listdummy2)>=nc_steps):
            Aid=(np.full(shape=[nc_steps],fill_value=np.NaN))
            for i in range(nc_steps):
                Aid[i]=listdummy2[i]
    listdummy.clear()
    listdummy2.clear()

    # Initialize arrays for reference and delta amplitudes
    amplitudes_ref = (np.full(shape = [n_E, int0], 
        fill_value = np.nan, dtype = complex))
    amplitudes_del = (np.full(shape = [n_E, nc_steps, int0], 
        fill_value = np.nan, dtype = complex))

    # maybe working arrays for transfer into amplitude arrays ?

    #End of the Header - Start of Reading in the Delta Data
    for e_index in range(n_E): # Energy loop

        # Energy, VPI and VV header
        line=next(file_lines)
        listdummy=ff_reader_6E13_7.read(line)
        for part in listdummy:
            if(part is not None):
                listdummy2.append(part)
        E, VPI, VV = listdummy2
        #transversing  Hartree to eV
        E=E*27.211396641308                            
        E_array[e_index] = E
        VPI_array[e_index] = VPI
        VV_array[e_index] = VV
        listdummy.clear()
        listdummy2.clear()

        # Reference amplitudes
        while(len(listdummy2)<2*int0):
            line=next(file_lines)
            listdummy=ff_reader_6E13_7.read(line)
            for part in listdummy:
                listdummy2.append(part)
        for j in range(int0):
            amplitudes_ref[e_index, j] = complex(listdummy2[2*j], listdummy2[2*j+1])
        listdummy.clear()
        listdummy2.clear()       
        
        # Delta amplitudes
        while(len(listdummy2)<2*int0*nc_steps):
            line=next(file_lines)
            listdummy=ff_reader_6E13_7.read(line)
            for part in listdummy:
                listdummy2.append(part)
        for j in range(nc_steps):
            for k in range(int0):
                amplitudes_del[e_index, j, k] = complex(listdummy2[2*int0*j+k*2], listdummy2[2*int0*j+k*2+1]) 
        listdummy.clear()
        listdummy2.clear()
        
    return (
        (phi, theta),
        (trar1, trar2),
        (int0, n_atoms, nc_steps),
        Beam_places,
        Cundisp,
        CDisp,
        (E_array, VPI_array, VV_array),
        amplitudes_ref,
        amplitudes_del)


#ganz dringend die -0.0000001 weggeben, nur da wegen ganzzahligen Indizes bei denen x2 zu groß werden würde
def GetInt(HeaderBlock1_g, HeaderBlock2_g, HeaderBlock3_g, HeaderBlock5_g, E_list_g, VPI_list_g, VV_list_g, ReferenzAmplituden_g, DeltaAmplituden_g, \
          NCSurf, delta_steps, filename_list):
    # TODO:
    #  - all inputs and outputs should be numpy arrays. This is VERY important for performance
    # inputs should not contain filenames or dicts
    # for case of multiple input arrays, we should use touples, keyword arguments etc
    # - output ATSAS needs to be in a very specific format - talk to me about what it should look like
    # - improve readability: add comments and refactor stuff into short functions where possible
    # - comments in Englisch

    #Conc wird probably auch als input parameter genommen
    CXDisp=1000
    XDisp=0
    counter=0
    Conc=1
    for name in filename_list:
        HeaderBlock5=HeaderBlock5_g[name]
        if(NCSurf[counter]==1):
            Int0, NAtoms, NCSteps = HeaderBlock2_g[name]  
            for j in range(NAtoms):
                fill=delta_steps[counter]-0.0000001     #probably besserer Name als fill finden
                x1=int(fill//1)
                x2=x1+1
                x=fill-x1
                y1=HeaderBlock5[0][x1,j,0]
                y2=HeaderBlock5[0][x2,j,0]
                CDisp=x*(y2-y1)+y1
                XDisp=Conc*CDisp
                if(XDisp<CXDisp):
                    CXDisp=XDisp
        counter=counter+1

    
    #unter der Annahme, dass HB1, HB3, E, VPI, VV, RA Listen immer gleich sind für alle Files
    for name in filename_list:
        HeaderBlock1=HeaderBlock1_g[name]
        #HeaderBlock2=HeaderBlock2_g[name]
        HeaderBlock3=HeaderBlock3_g[name]
        #HeaderBlock5=HeaderBlock5_g[name]
        E_list=E_list_g[name]
        VPI_list=VPI_list_g[name]
        VV_list=VV_list_g[name]
        ReferenzAmplituden=ReferenzAmplituden_g[name]
        #DeltaAmplituden=DeltaAmplituden_g[name]
        
        Trar1=[0,0]
        Trar2=[0,0]
        theat, phi, Trar1[0], Trar1[1], Trar2[0], Trar2[1] = HeaderBlock1
        break # @Tobi: what does this loop do? It break at first iteration..

        #Int0, NAtoms, NCSteps = HeaderBlock2    

    ATSAS_matrix=np.zeros([len(E_list), Int0])
        
              
    #Loop über die Energien
    for Energy in range(len(E_list)):
        #Definieren von Variablen, die in der jeweiligen Energie gleichbleiben
        E=E_list[Energy]
        VV=VV_list[Energy]
        VPI=VPI_list[Energy]

        AK=np.sqrt(max(2*E-2*VV,0))
        C=AK*np.cos(theat)
        BK2=AK*np.sin(theat)*np.cos(phi)
        BK3=AK*np.sin(theat)*np.sin(phi)
        BKZ=np.sqrt(complex(2*E-BK2**2-BK3**2,-2*VPI))

        #Loop über die Beams
        for Beam in range(Int0):
            #Variablen pro Beam
            h=HeaderBlock3[Beam][0]
            k=HeaderBlock3[Beam][1]
            AK2=BK2+h*Trar1[0]+k*Trar2[0]
            AK3=BK3+h*Trar1[1]+k*Trar2[1]
            AK=2*E-AK2**2-AK3**2
            AKZ=complex(AK,-2*VPI)
            APERP=AK-2*VV

            #Hier kommen später noch die Delta Amplituden ebenfalls ins Spiel, wird gerade gemacht
            #Herausfinden welche NCStep deltas man nimmt mit delta_step matrix
            DelAct=ReferenzAmplituden[Energy][Beam]
            counter=0
            for name in filename_list:
                #noch genau schauen wie genau gewollt
                Int0, NAtoms, NCSteps = HeaderBlock2_g[name]
                
                #Interpolation der Float delta_step Werte
                #iwie random dass x1 extra als int definiert werden muss, damit sich python nicht aufregt
                j=delta_steps[counter]-0.0000001
                x1=int(j//1)
                x2=x1+1
                x=j-x1
                y1=DeltaAmplituden_g[name][Energy][x1,Beam]
                y2=DeltaAmplituden_g[name][Energy][x2,Beam]
                Interpolation=x*(y2-y1)+y1

                
                DelAct=DelAct+Interpolation
                counter=counter+1


            if(APERP>0):
                A=np.sqrt(APERP)
                PRE=(BKZ+AKZ)*CXDisp
                PRE=np.exp(complex(0,-1)*PRE)
                amp_abs=abs(DelAct)
                if(amp_abs>10e10):
                    ATSAS=10e20
                else:
                    RPRE=abs(PRE)
                    ATSAS=amp_abs**2*RPRE**2*A/C
            else:
                ATSAS=0

            ATSAS_matrix[Energy,Beam]=ATSAS
            #Atsas hier in ein Dictionary mit dem Key name abspeichern; nope worked so immernoch
    
    return ATSAS_matrix