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
from tqdm import tqdm
from numba import njit, prange

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
    
    E_kin_array : ndarray
    Array that contains the kinetic energies of the elastically scatted
    electrons inside the crystal. (Not the incidence energy!)
    
    VPI_array : ndarray
    Imaginary part of the inner potential of the surface
    
    VV_array : ndarray
    Real part of the inner potential of the surface
    
    beam_indices : ndarray
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
    E_kin_array= np.full(n_E, fill_value = np.NaN)
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
            beam_indices=(np.full(shape=[int0,2],fill_value=np.NaN))
            for i in range(int0):
                for j in range(2):
                    beam_indices[i,j]=listdummy2[2*i+j]
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
        E_kin, VPI, VV = listdummy2
        #transversing  Hartree to eV
        E_kin=hartree_to_eV(E_kin)
        E_kin_array[e_index] = E_kin
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
        beam_indices,
        Cundisp,
        CDisp,
        (E_kin_array, VPI_array, VV_array),
        amplitudes_ref,
        amplitudes_del)


@njit(fastmath=True, parallel=True, nogil=True)
def GetInt(phi, theta, trar1, trar2, Beam_variables, beam_indices, ph_CDisp, E_array, VPI_array, VV_array, amplitudes_ref, amplitudes_del, \
          n_files, NCSurf, delta_steps,number_z_steps):
    """This function reads in the values of the function Transform and uses them to get the ATSAS_matrix
    
    Parameters
    ----------
    n_files: int
    Number of files
    
    phi, theta:
    Angles of how the beam hits the sample
    
    trar1, trar2:
    Vectors of the normal and the reciprocal unit cell of the sample
    
    Beam_variables:
    The variables int0, n_atoms, nc_steps for each file stored in an array
    
    beam_indices:
    Array with the order of beams
    
    ph_CDisp:
    Geometric displacements of the atom
    
    E_array:
    Array that contains all the energies of the file
    
    VPI_array:
    Imaginary part of the inner potential of the surface
    
    VV_array:
    Real part of the inner potential of the surface
    
    amplitudes_ref:
    Array that contains all values of the reference amplitudes
    
    amplitudes_del:
    Array that contains all values of the delta amplitudes
    
    filename_list:
    List of the filenames that contain the data
    
    NCSurf:
    List of 0 and 1 to decide which file takes part in creating the CXDisp
    
    delta_steps:
    List of numbers that decide which geometric displacement this atom has
    
    number_z_steps:
    Total number of different delta_z values
    
    
    Returns
    ----------
    ATSAS_matrix:
    Array that contains the intensities of the different beams with each energy
    """
    

    #Conc will probably taken as Input parameter aswell
    CXDisp=1000
    XDisp=0
    Conc=1
    for i in range(n_files):
        if(NCSurf[i]==1):
            int0, n_atoms, nc_steps = Beam_variables[i,:]
            int0=int(int0)
            n_atoms=int(n_atoms)
            nc_steps=int(nc_steps)
            for j in range(n_atoms):
                fill=delta_steps[i]    #probably besserer Name als fill finden
                x1=int(fill//1)
                #necessary so that the x2 doesnt go out of the index range
                if(x1==number_z_steps-1):
                    x1=x1-1
                x2=x1+1
                x=fill-x1
                y1=ph_CDisp[i,x1,j,0]
                y2=ph_CDisp[i,x2,j,0]
                CDisp=x*(y2-y1)+y1
                XDisp=Conc*CDisp
                if(XDisp<CXDisp):
                    CXDisp=XDisp

    ATSAS_matrix=np.zeros((len(E_array), int0))       
              
    #Loop über die Energien
    for e_index in prange(len(E_array)):
        #Definieren von Variablen, die in der jeweiligen Energie gleichbleiben
        E=E_array[e_index]
        VV=VV_array[e_index]
        VPI=VPI_array[e_index]

        AK=np.sqrt(max(2*E-2*VV,0))
        C=AK*np.cos(theta)
        BK2=AK*np.sin(theta)*np.cos(phi)
        BK3=AK*np.sin(theta)*np.sin(phi)
        BKZ=np.sqrt(complex(2*E-BK2**2-BK3**2,-2*VPI))

        #Loop über die Beams
        for b_index in range(int0):
            #Variablen per Beam
            h=beam_indices[b_index,0]
            k=beam_indices[b_index,1]
            AK2=BK2+h*trar1[0]+k*trar2[0]
            AK3=BK3+h*trar1[1]+k*trar2[1]
            AK=2*E-AK2**2-AK3**2
            AKZ=complex(AK,-2*VPI)
            APERP=AK-2*VV

            #Herausfinden welche NCStep deltas man nimmt mit delta_step matrix
            DelAct=amplitudes_ref[e_index,b_index]
            for i in range(n_files):
                #noch genau schauen wie genau gewollt
                int0, n_atoms, nc_steps = Beam_variables[i]
                int0=int(int0)
                n_atoms=int(n_atoms)
                nc_steps=int(nc_steps)
                
                #Interpolation of float delta_step values
                #iwie random dass x1 extra als int definiert werden muss, damit sich python nicht aufregt
                fill=delta_steps[i]
                x1=int(fill//1)
                #necessary so that the x2 doesnt go out of the index range
                if(x1==number_z_steps-1):
                    x1=x1-1
                x2=x1+1
                x=fill-x1
                y1=amplitudes_del[i,e_index,x1,b_index]
                y2=amplitudes_del[i,e_index,x2,b_index]
                Interpolation=x*(y2-y1)+y1

                DelAct=DelAct+Interpolation



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

            ATSAS_matrix[e_index,b_index]=ATSAS

    
    return ATSAS_matrix




@njit(fastmath=True,parallel=True, nogil=True)
def GetInt2D(phi, theta, trar1, trar2, Beam_variables, beam_indices, ph_CDisp, E_array, VPI_array, VV_array, amplitudes_ref, amplitudes_del, \
          n_files, NCSurf, delta_steps,number_z_steps,number_vib_steps):
    
    """This function reads in the values of the function Transform and uses them to get the ATSAS_matrix
    
    Parameters
    ----------
    phi, theta:
    Angles of how the beam hits the sample
    
    trar1, trar2:
    Vectors of the normal and the reciprocal unit cell of the sample
    
    Beam_variables:
    The variables int0, n_atoms, nc_steps for each file stored in an array
    
    beam_indices:
    Array with the order of beams
    
    ph_CDisp:
    Geometric displacements of the atom
    
    E_array:
    Array that contains all the energies of the file
    
    VPI_array:
    Imaginary part of the inner potential of the surface
    
    VV_array:
    Real part of the inner potential of the surface
    
    amplitudes_ref:
    Array that contains all values of the reference amplitudes
    
    amplitudes_del:
    Array that contains all values of the delta amplitudes
    
    n_files: int
    Number of files
    
    NCSurf:
    List of 0 and 1 to decide which file takes part in creating the CXDisp
    
    delta_step:
    List of numbers that decide which geometric displacement this atom has
    
    number_z_steps:
    Total number of different delta_z values
    
    number_vib_steps:
    Total number of different delta_vib values
    
    Returns
    ----------
    ATSAS_matrix:
    Array that contains the intensities of the different beams with each energy
    """

    #Conc will probably taken in as an Input parameter
    CXDisp=1000
    XDisp=0
    Conc=1
    for i in range(n_files):
        if(NCSurf[i]==1):
            int0, n_atoms, nc_steps = Beam_variables[i,:]
            int0=int(int0)
            n_atoms=int(n_atoms)
            nc_steps=int(nc_steps)
            for j in range(n_atoms):
                fill=delta_steps[i,0]  #probably besserer Name als fill finden #2er numpy array machen
                x1=int(fill//1)
                #necessary so that the x2 doesnt go out of the index range
                if(x1==number_z_steps-1):
                    x1=x1-1
                x2=x1+1
                x=fill-x1
                y1=ph_CDisp[i,x1,j,0]
                y2=ph_CDisp[i,x2,j,0]
                CDisp=x*(y2-y1)+y1
                XDisp=Conc*CDisp
                if(XDisp<CXDisp):
                    CXDisp=XDisp

    ATSAS_matrix=np.zeros((len(E_array), int0))       
              
    #Loop over the energies
    for e_index in prange(len(E_array)):
        #defining variables that stay the same for their respective energy
        E=E_array[e_index]
        VV=VV_array[e_index]
        VPI=VPI_array[e_index]

        AK=np.sqrt(max(2*E-2*VV,0))
        C=AK*np.cos(theta)
        BK2=AK*np.sin(theta)*np.cos(phi)
        BK3=AK*np.sin(theta)*np.sin(phi)
        BKZ=np.sqrt(complex(2*E-BK2**2-BK3**2,-2*VPI))

        #Loop over the beams
        for b_index in range(int0):
            #variables that differ with each beam
            h=beam_indices[b_index,0]
            k=beam_indices[b_index,1]
            AK2=BK2+h*trar1[0]+k*trar2[0]
            AK3=BK3+h*trar1[1]+k*trar2[1]
            AK=2*E-AK2**2-AK3**2
            AKZ=complex(AK,-2*VPI)
            APERP=AK-2*VV

            #finding out which NCStep you take with the delta_step matrix
            DelAct=amplitudes_ref[e_index,b_index]
            for i in range(n_files):
                int0, n_atoms, nc_steps = Beam_variables[i]
                int0=int(int0)
                n_atoms=int(n_atoms)
                nc_steps=int(nc_steps)
                
                #Interpolation of the float delta_step values
                #iwie random dass x1 extra als int definiert werden muss, damit sich python nicht aufregt

                z_step=delta_steps[i,0]
                vib_step=delta_steps[i,1]
                
                z1=int(z_step//1)
                #necessary so that the z2 doesnt go out of the index range
                if(z1==number_z_steps-1):
                    z1=z1-1
                z2=z1+1
                v1=int(vib_step//1)
                #necessary so that the v2 doesnt go out of the index range
                if(v1==number_vib_steps-1):
                    v1=v1-1
                v2=v1+1
                
                z=z_step-z1
                v=vib_step-v1
                
               
                #could be an error here, need to look closer at this
                #maybe change the delta_step intervalls to see how the peaks behave
                value1=amplitudes_del[i,e_index,v1*11+z1,b_index]
                value2=amplitudes_del[i,e_index,v1*11+z2,b_index]
                ip1=z*(value2-value1)+value1
               
                value3=amplitudes_del[i,e_index,v2*11+z1,b_index]
                value4=amplitudes_del[i,e_index,v2*11+z2,b_index]
                ip2=z*(value4-value3)+value3
                
                ip_final=v*(ip2-ip1)+ip1
                        
                DelAct=DelAct+ip_final

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

            ATSAS_matrix[e_index,b_index]=ATSAS

    
    return ATSAS_matrix



def hartree_to_eV(E_hartree):
    return E_hartree*27.211396



def Transform(n_E,directory):
    """This function transforms the read in data to a form, where it can be read in by the GetInt function

    Parameters
    ----------
    n_E:
    Number of different energies of one file
    
    directory:
    Relative path to the file that gets read in
    
    
    Returns 
    ----------
    phi, theta:
    Angles of how the beam hits the sample
    
    trar1, trar2:
    Vectors of the normal and the reciprocal unit cell of the sample
    
    Beam_variables:
    The variables int0, n_atoms, nc_steps for each file stored in an array
    
    beam_indices:
    Array with the order of beams
    
    CDisp:
    Geometric displacements of the atom
    
    E_array:
    Array that contains all the energies of the file
    
    VPI_array:
    Imaginary part of the inner potential of the surface
    
    VV_array:
    Real part of the inner potential of the surface
    
    amplitudes_ref:
    Array that contains all values of the reference amplitudes
    
    amplitudes_del:
    Array that contains all values of the delta amplitudes
    
    filename_list:
    List of the filenames that contain the data
    """
    
    filename_list=[]
    data_list_all={}
    for filenames in os.walk(directory):
        filename_list.append(filenames) 
    filename_list=filename_list[0][2]
    
    for name in tqdm(filename_list):
        filename=directory+name
        data_list_all[name]=read_delta_file(filename,n_E)

    
    #constant variables
    phi,theta=data_list_all[filename_list[0]][0]
    trar1,trar2=data_list_all[filename_list[0]][1]
    beam_indices=data_list_all[filename_list[0]][3]    
    E_array,VPI_array,VV_array=data_list_all[filename_list[0]][6]
    amplitudes_ref=data_list_all[filename_list[0]][7]
    n_files = len(filename_list)
    
    #int0, n_atoms, nc_steps in an array(nc_steps can change, n_atoms maybe too)
    counter=0
    Beam_variables=np.full(shape=[n_files,3],fill_value=np.NaN)
    for name in filename_list:
        Beam_variables[counter,:]=data_list_all[name][2]
        counter=counter+1
     
    int0=int(np.max(Beam_variables[:,0]))
    n_atoms_max=int(np.max(Beam_variables[:,1]))
    nc_steps_max=int(np.max(Beam_variables[:,2]))   
    
    #saving the changing data in arrays
    counter=0
    CDisp=np.full(shape=(n_files,nc_steps_max,n_atoms_max,3),fill_value=np.NaN)
    amplitudes_del=np.full(shape=(n_files,n_E,nc_steps_max,int0),fill_value=np.NaN, dtype=complex)
    for name in filename_list:
        int0,n_atoms,nc_steps=Beam_variables[counter]
        int0=int(int0)
        n_atoms=int(n_atoms)
        nc_steps=int(nc_steps)
        CDisp[counter,0:nc_steps,0:n_atoms,:]=data_list_all[name][5]
        amplitudes_del[counter,:,0:nc_steps,:]=data_list_all[name][8]
        counter=counter+1
    
    return (
            phi, theta,
            trar1, trar2,
            Beam_variables, 
            beam_indices,
            CDisp, 
            E_array, VPI_array, VV_array,
            amplitudes_ref,
            amplitudes_del, 
            filename_list)



def PlotMaker(phi,theta,trar1,trar2,Beam_variables,beam_indices,CDisp,E_array,VPI_array,VV_array,amplitudes_ref,amplitudes_del, \
                          n_files,NCSurf,delta_step_einzeln, phi_g,theta_g,trar1_g,trar2_g,Beam_variables_g,beam_indices_g, \
                          CDisp_g,E_array_g,VPI_array_g,VV_array_g,amplitudes_ref_g,amplitudes_del_g,n_files_g,delta_step_g, \
                          int0,number_z_steps,number_vib_steps):  
   
    #This function determines the difference between a file, where 2 contributors to the intensities are considered, and 2 files, where only one of these 
    #contributors each is considered. That way we can determine, if it is viable to use linear combination of the 3 dimensions + the vib files to describe
    #each delta of the atom.
    #
    #INPUT:
    #
    #phi, theta, phi_g, theta_g:
    #Angles of how the beam hits the sample
    #
    #trar1, trar2, trar1_g, trar2_g:
    #Vectors of the normal and the reciprocal unit cell of the sample
    #
    #Beam_variables, Beam_variables_g:
    #The variables int0, n_atoms, nc_steps for each file stored in an array
    #
    #beam_indices, beam_indices_g:
    #Array with the order of beams
    #
    #CDisp, CDisp_g:
    #Geometric displacements of the atom
    #
    #E_array, E_array_g:
    #Array that contains all the energies of the file
    #
    #VPI_array, VPI_array_g:
    #Imaginary part of the inner potential of the surface
    #
    #VV_array, VV_array_g:
    #Real part of the inner potential of the surface
    #
    #amplitudes_ref, amplitudes_ref_g:
    #Array that contains all values of the reference amplitudes
    #
    #amplitudes_del, amplitudes_del_g:
    #Array that contains all values of the delta amplitudes
    #
    #filename_list, filename_list_g:
    #List of the filenames that contain the data
    #
    #NCSurf:
    #List of 0 and 1 to decide which file takes part in creating the CXDisp
    #
    #delta_step, delta_step_g:
    #List of numbers that decide which geometric displacement this atom has
    #
    #int0:
    #number of beams
    #
    #number_z_steps:
    #Total number of different delta_z values
    #
    #number_vib_steps:
    #Total number of different delta_vib values
    #
    #OUTPUT:
    #ATSAS_matrix:
    #Array that contains the intensities of the different beams with each energy
    
    matrix_differenz=np.full(shape=(len(E_array),int0),fill_value=np.NaN)
    array=np.full(shape=(101,101),fill_value=np.NaN)
    for i in tqdm(range(101)): #gerade weird dringend drüber nachdenken, 0.1 schritte ist dann diese range sinnvoll?
        for j in range(101):
            delta_step_einzeln=np.array([j/10,5.,5.,5.,5.,i/10,5.,5.,5.,5.])

            delta_step_g=np.full(shape=(5,2),fill_value=5.)
            delta_step_g[0,0]=j/10
            delta_step_g[0,1]=i/10

            matrix_z_vib_einzeln=GetInt(phi,theta,trar1,trar2,Beam_variables,beam_indices,CDisp,E_array,VPI_array,VV_array,amplitudes_ref,amplitudes_del, \
                                      n_files,NCSurf,delta_step_einzeln,number_z_steps)

            matrix_z_vib_g=GetInt2D(phi_g,theta_g,trar1_g,trar2_g,Beam_variables_g,beam_indices_g,CDisp_g,E_array_g,VPI_array_g,VV_array_g, \
                                  amplitudes_ref_g,amplitudes_del_g, n_files_g,NCSurf,delta_step_g,number_z_steps,number_vib_steps)

            for e_index in range(len(E_array)):
                for b_index in range(int0):
                    matrix_differenz[e_index,b_index]=matrix_z_vib_g[e_index,b_index]-matrix_z_vib_einzeln[e_index,b_index]


            Integral=0
            for b in range(int0):
                Integral=Integral+abs(np.trapz(matrix_differenz[:,b]))

            array[j,i]=Integral

    return array



def PlotMaker1D(phi,theta,trar1,trar2,Beam_variables,beam_indices,CDisp,E_array,VPI_array,VV_array,amplitudes_ref,amplitudes_del, \
                n_files,NCSurf,phi_g,theta_g,trar1_g,trar2_g,Beam_variables_g,beam_indices_g,CDisp_g, \
                E_array_g,VPI_array_g,VV_array_g, amplitudes_ref_g,amplitudes_del_g,n_files_g,int0,number_z_steps):    
    
    array4=np.full(shape=(101),fill_value=np.NaN)
    matrix_differenz=np.full(shape=(len(E_array),int0),fill_value=np.NaN)
    
    for i in tqdm(range(101)):

        delta_step_einzeln=np.full(shape=(40),fill_value=5.)
        delta_step_einzeln[0]=i*1/(10*np.sqrt(2))+1.464
        delta_step_einzeln[1]=i*1/(10*np.sqrt(2))+1.464

        delta_step_dig=np.full(shape=(20),fill_value=5.)
        delta_step_dig[0]=i/10


        matrix_z_vib_einzeln=GetInt(phi,theta,trar1,trar2,Beam_variables,beam_indices,CDisp,E_array,VPI_array,VV_array,amplitudes_ref,amplitudes_del, \
                      n_files,NCSurf,delta_step_einzeln,number_z_steps)

        matrix_z_vib_g=GetInt(phi_g,theta_g,trar1_g,trar2_g,Beam_variables_g,beam_indices_g,CDisp_g,E_array_g,VPI_array_g,VV_array_g, \
                              amplitudes_ref_g,amplitudes_del_g,n_files_g,NCSurf,delta_step_dig,number_z_steps)

        for e_index in range(len(E_array)):
            for b_index in range(int0):
                matrix_differenz[e_index,b_index]=matrix_z_vib_g[e_index,b_index]-matrix_z_vib_einzeln[e_index,b_index]


        Integral=0
        for b in range(int0):
              Integral=Integral+abs(np.trapz(matrix_differenz[:,b]))
        
        array4[i]=Integral

    return array4 