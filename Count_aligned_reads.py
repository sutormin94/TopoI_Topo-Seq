###############################################
##Dmitry Sutormin, 2021##
##Topo-Seq analysis##

#Script takes WIG files as input (N3E, N5E) and counts the number of aligned reads.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import random
from Bio import SeqIO

#Path to the directory with input files.
Path_to_input_files='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\WIG_aligned\\'


#######
#Parses WIG file with N3/5E values.
#Computes a total number of Ends.
#######

def wig_parsing(wigfile):
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    NE_values=[]
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(float(line[0]))
    wigin.close()
    NE_values=np.array(NE_values)
    print(f'Total genome length: {len(NE_values)}, total number of aligned read pairs: {sum(NE_values)}')
    return NE_values


#######
#Wrapper function.
#######

def wrapper(input_files_path):
    
    #List WIG files to be masked.
    files_list=os.listdir(input_files_path)
    
    #Read, mask, write WIG files.
    for file in files_list:
        print(file)
        WIG_data=wig_parsing(input_files_path+file)
    
    return

wrapper(Path_to_input_files)