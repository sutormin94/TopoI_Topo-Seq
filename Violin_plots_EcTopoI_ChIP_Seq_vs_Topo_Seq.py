###############################################
##Dmitry Sutormin, 2021##
##TopoI Topo-Seq analysis##

#The script tests sets of genomic intervals (Peaks, TUs, BIMEs-1, BIMEs-2, IHF sites, Fis sites, H-NS sites, MatP sites, etc.)
#for the enrichment with some continously distributed character (RNApol fold enrichment, score, GC%, etc.) (t-test).
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.stats import pearsonr
from scipy.stats import binom

#######
#Variables to be defined.
#######

#Path to the working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\\"

#Input: Intervals (e.g. EcTopoI peaks) with additional info (FE, GC%, width, etc.).
path_to_intervals_data_dict={'EcTopoI_Topo_all'    : PWD + 'EcTopoI_ChIP_Seq_vs_Topo_Seq\TopoI_Ara_TCSs_called_15.BroadPeak',
                             'EcTopoI_Topo_shared' : PWD + 'EcTopoI_ChIP_Seq_vs_Topo_Seq\EcTopoI_ChIP_shared_with_EcTopoI_Topo.narrowPeak',
                             'EcTopoI_Topo_unique' : PWD + 'EcTopoI_ChIP_Seq_vs_Topo_Seq\EcTopoI_Topo_unique_from_EcTopoI_ChIP.narrowPeak',
                             }

path_to_wig_files_dict={'RpoC_ChIP_FE' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\Borukhov_RpoC_Pol_Sofi_LB_FE.wig'}

#Output: path to the dir to store output
Outputpath=PWD + "EcTopoI_ChIP_Seq_vs_Topo_Seq\\"
if not os.path.exists(Outputpath):
    os.makedirs(Outputpath)
    
    

#######
#Read peaks data.
#######

def read_peaks(path_dict, enlarge):
    peaks_dict={}
    for name, inpath in path_dict.items():
        print(f'Now is processing {name} peaks')
        filein=open(inpath, 'r')
        peaks_ar=[]        
        if enlarge==False:
            for line in filein:
                line=line.rstrip().split('\t')
                start=int(line[1])
                end=int(line[2])
                peaks_ar.append([start, end])
        elif enlarge==True:
            for line in filein:
                line=line.rstrip().split('\t')
                start=int(line[1])-100
                end=int(line[2])+100
                peaks_ar.append([start, end])  
        filein.close() 
        
        peaks_dict[name]=peaks_ar

    return peaks_dict


#######
#Read wig data.
#######

def read_wig(path_dict):
    wig_dict={}
    for name, inpath in path_dict.items():
        print(f'Now is processing {name} wig')
        filein=open(inpath, 'r')
        wig_ar=[]    
        for line in filein:
            line=line.rstrip().split(' ')
            if line[0] not in ['track', 'fixedStep']:
                wig_ar.append(float(line[0])) 
        wig_dict[name]=wig_ar
    
    return wig_dict

#######
#Mark peaks.
#######

def mark_peaks_return_subsets(peaks_dict, wig_dict):
    deletions_ar=[[274500, 372148], [793800, 807500], [1199000, 1214000]]
    sets_dict={}
    for wig_name, wig in wig_dict.items():
        for peaks_set_name, peaks_set in peaks_dict.items():
            print(f'{wig_name} masked with {peaks_set_name} peaks')
            mask=[0]*4647454
            for peak in peaks_set:
                mask[peak[0]:(peak[1]+1)]=[1]*((peak[1]+1)-peak[0])
            
            wig_in_peaks=[]
            wig_out_peaks=[]
            for i in range(len(mask)):
                deletion_check=0
                for deletion in deletions_ar:
                    if (i>deletion[0]) & (i<deletion[1]):
                        deletion_check=1
                if deletion_check==0:
                    if wig[i]>0: 
                        if mask[i]==0:
                            wig_out_peaks.append(wig[i])
                        else:
                            wig_in_peaks.append(wig[i])
                    
            sets_dict[f'{wig_name}_wig_in_{peaks_set_name}_peaks']=wig_in_peaks
            sets_dict[f'{wig_name}_wig_out_{peaks_set_name}_peaks']=wig_out_peaks
            
    return sets_dict

#######
#Draw violin plots.
#######

def set_axis_style(ax, labels, pos):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(pos)
    ax.set_xticklabels(labels, size=6)
    ax.set_xlim(0.25, max(pos)+0.75)
    return

def draw_violins(dataset1, signal_name, interval_names, pwd):
    
    #Positions of violins.
    pos1=[1, 2, 3, 4]

    #Draw violin-plots.
    fig=plt.figure(figsize=(6,6), dpi=100)
    plt1=fig.add_subplot(1,1,1) 
    violins=plt1.violinplot(dataset1, positions=pos1, widths=0.77, showmeans=True, showmedians=True, points=2000)
    for i in range(len(violins['bodies'])):
        violins['bodies'][i].set_facecolor('#ff7762')
        violins['bodies'][i].set_edgecolor('black')
        violins['bodies'][i].set_alpha(1)
    
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
      
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    
    vbars=violins['cbars']
    vbars.set_linewidth(1)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    
    labels=interval_names
    set_axis_style(plt1, labels, pos1)
    
    yticknames1=np.arange(0, 10, 1)
    plt1.set_yticks(yticknames1, minor=False)
    plt1.set_yticklabels(yticknames1)
    plt1.set_ylabel(f'{signal_name} FE', size=15)
    plt1.set_ylim(-0.5, 10)
    plt.setp(plt1.set_yticklabels(yticknames1), rotation=0, fontsize=15)   
    #plt1.set_yscale('log')
    
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(dataset1[0]),2)}', xy=(0.5, 2), xycoords='data', size=12, rotation=90)
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(dataset1[1]),2)}', xy=(1.5, 2), xycoords='data', size=12, rotation=90)
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(dataset1[2]),2)}', xy=(2.5, 2), xycoords='data', size=12, rotation=90)
    plt1.annotate(r"   $\overline{X}$"+f'={round(np.mean(dataset1[3]),2)}', xy=(3.5, 2), xycoords='data', size=12, rotation=90)    
    
    #Welch t-test.
    for i in range(len(dataset1)):
        for j in range(len(dataset1)):
            if j>i:
                Intervals_stat=stats.ttest_ind(dataset1[i], dataset1[j], equal_var=False)
                print(f'Sets in comparison: {interval_names[i]} vs {interval_names[j]}')
                print(f'\nT-test FE means\nMean1={round(np.mean(dataset1[i]),2)} Mean2={round(np.mean(dataset1[j]),2)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')

    plt.show()
    plt.tight_layout()
    plt.savefig(f'{pwd}\{signal_name}_FE_and_TopoI_Topo_Seq_TCSs.png', dpi=400, figsize=(4.5,6))   
    
    return


#######
#Wrapper function.
#######

def violin_wrapper(path_to_intervals_data_dict, path_to_wig_files_dict, enlarge, Outpath):

    Peaks_dictionary=read_peaks(path_to_intervals_data_dict, enlarge)
    WIG_dictionary=read_wig(path_to_wig_files_dict)
    Sets_dictionary=mark_peaks_return_subsets(Peaks_dictionary, WIG_dictionary)
    
    #Dataset. Specify names of datasets to be compared.
    signal_name='RpoC_ChIP_FE'
    interval_names=list(Peaks_dictionary.keys())
    interval_names.append('Not_EcTopoI_TCSs')
    dataset1=[Sets_dictionary[f'{signal_name}_wig_in_{interval_names[0]}_peaks'], Sets_dictionary[f'{signal_name}_wig_in_{interval_names[1]}_peaks'], Sets_dictionary[f'{signal_name}_wig_in_{interval_names[2]}_peaks'], Sets_dictionary[f'{signal_name}_wig_out_{interval_names[0]}_peaks']]
    
    draw_violins(dataset1, signal_name, interval_names, Outpath)
    
    return

violin_wrapper(path_to_intervals_data_dict, path_to_wig_files_dict, True, Outputpath)