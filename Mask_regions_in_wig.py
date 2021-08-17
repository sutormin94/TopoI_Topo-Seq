###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#Script takes WIG files as input and masks specified regions (by default substitutes by 0).
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
Path_to_input_files='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\WIG_NE_strand_specific\\'

#Path to the directory with output files.
Path_to_output_files='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\WIG_NE_strand_specific_masked\\'

#Path to the file with intervals to be masked (bed).
Path_to_masked_regions='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\Scripts_TopoI_Topo-seq\Additional_genome_features\\Regions_to_be_masked.broadPeak'


#######
#Checks if directory exists and if not creates.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return


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
    return NE_values


#######
#Write .wig file.
#######

def write_wig(ar, fileout_path, name):
    fileout=open(fileout_path, 'w')
    fileout.write(f'track type=wiggle_0 name="{name}" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom=NC_007779.1_w3110_Mu start=1 step=1\n')
    for point in ar:
        fileout.write(f'{point}\n')
    fileout.close()
    return


#######
#Opens and reads BED file with deletions coordinates.
#Example:
#GenomeID\tStart\tEnd
#NC_007779.1_w3110_Mu\t274500\t372148
#######

def deletions_info(del_path):
    del_ar=[]
    filein=open(del_path, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        del_ar.append([int(line[1]), int(line[2])])
    filein.close()
    return del_ar


#######
#Mask regions of a WIG file.
#######

def mask_array(NE_values, del_ar, sub):
    
    NE_values_masked=NE_values
    
    for region in del_ar:
        maska=[sub]*(region[1]-region[0])
        print(len(maska), len(NE_values[region[0]:region[1]]))
        NE_values_masked[region[0]:region[1]]=maska
    
    return NE_values_masked


#######
#Wrapper function.
#######

def wrapper(input_files_path, output_files_path, masked_path, sub):
    
    #Read masked regions data.
    masked_regions_ar=deletions_info(masked_path)
    #List WIG files to be masked.
    files_list=os.listdir(input_files_path)
    #Create output path.
    Dir_check_create(output_files_path)
    
    #Read, mask, write WIG files.
    for file in files_list:
        print(file)
        WIG_data=wig_parsing(input_files_path+file)
        print(len(WIG_data))
        WIG_data_masked=mask_array(WIG_data, masked_regions_ar, sub)
        print(len(WIG_data_masked))
        write_wig(WIG_data_masked, output_files_path+file, file.split('.')[0])
    
    return

wrapper(Path_to_input_files, Path_to_output_files, Path_to_masked_regions, 0)