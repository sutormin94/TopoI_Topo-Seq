###############################################
##Dmitry Sutormin, 2020##
##ChIP-Seq data analysis##

#Tests the significance of the overlap between two ChIP-Seq peak sets by Monte-Carlo simulation.
###############################################

#######
#Packages to be imported.
#######

import random as rd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import scipy
from scipy.stats import norm


#######
#Import data.
#######

#Path to working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\\"

#Dictionary of datasets.
Input_dict={'EcTopoI_ChIP_Seq' : PWD + "Peak_calling\Reproducible_peaks\\TopoA_noCTD_noRif_rep123_thr_3_nm_0.001_peaks.narrowPeak",
            'EcTopoI_Topo_Seq' : PWD + "EcTopoI_vs_RNApol_Gyrase\RNApol_peaks\RpoC_Borukhov\RNApol_peaks_threshold_3.BroadPeak"}

#Pre-calculated data.
Simulated_data_dict={'EcTopoI_vs_EcRpoC' : PWD + "EcTopoI_overlap_EcRpoC_simulation\EcTopoI_vs_EcRpoC_peaks_overlay_10000_it_Real_delcor.txt",
                     'EcTopoI' : PWD + "EcTopoI_overlap_EcRpoC_simulation\EcTopoI_peaks_in_overlay_10000_it_Real_delcor.txt",
                     'EcRpoC' : PWD + "EcTopoI_overlap_EcRpoC_simulation\EcRpoC_peaks_in_overlay_10000_it_Real_delcor.txt"}

#Genome length.
Genome_len=4647454 #E. coli DY330 4647454; M. tuberculosis 4411532; M. smegmatis 6988302
Deletion_correction=126348 #Correcion on deletions cumulative length: E. coli DY330 126348; M. tuberculosis 0; M. smegmatis 0;

#Output data path.
Outpath=PWD + "EcTopoI_overlap_EcRpoC_simulation\\"

#################
### Read input sets. Calculate bin width and number of peaks. Get number of real overlaying peaks.
#################

def read_peak_file(in_path):
    peak_ar=[]
    filein=open(in_path, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        peak_ar.append([int(line[1]), int(line[2])])
    filein.close()
    return peak_ar

def read_data_bin_width(input_data):
    bin_info_dict={}
    Real_peaks_dict={}
    for name, datapath in input_data.items():
        peak_ar=read_peak_file(datapath)
        peak_width_ar=[]
        for peak in peak_ar:
            peak_width_ar.append(peak[1]-peak[0])
        bin_info_dict[name]=[len(peak_width_ar), int(np.mean(peak_width_ar))]
        Real_peaks_dict[name]=peak_ar
        print(f'{name} peaks number {len(peak_width_ar)} and mean width {int(np.mean(peak_width_ar))}')
    return bin_info_dict, Real_peaks_dict

#################
### Get number of real overlaying peaks.
#################

#######
#Indicate where peaks occures by addition of 1 to these positions to genome-length array.
#######

def Indicate_where_peaks(genome_ar, peaks_ar):
    for peak in peaks_ar:
        for i in range (peak[1]-peak[0]):
            genome_ar[peak[0]+i]+=1
    return genome_ar

#######
#Find overlapping peaks in genome-length array.
#######  

def Find_rep_peaks(genome_ar, thr):
    peak=0
    rep_peaks_ar=[]
    for i in range(len(genome_ar)):
        
        #Rules how to enter a peak.
        if genome_ar[i]<thr and peak==0: #We are not in peak.
            continue
        elif genome_ar[i]>=thr and peak==0: #We are at left peak border.
            peak=1
            current_peak=[i]
            continue
        
        #Rules how go deeper into the peak.
        elif genome_ar[i]>=thr and peak==1: #We are within a peak.
            continue
        elif genome_ar[i]<thr and peak==1: #We are at the right peak border.
            peak=0
            current_peak.append(i)
            rep_peaks_ar.append(current_peak)
            continue
    return rep_peaks_ar

#######
#Identifies reproducible peaks with threshold (number of samples in which a peak should be present) given.
#######  

def overlap_call(Real_peaks_dict, thr, genome_len):
    #Create template genome-long array.
    genome_ar=[0]*genome_len
    #Indicates peaks.
    for name, peaks_ar in Real_peaks_dict.items():
        genome_ar=Indicate_where_peaks(genome_ar, peaks_ar)
    #Identify reproducible peaks.
    Rep_peaks_array=Find_rep_peaks(genome_ar, thr)
    return Rep_peaks_array

#######
#Count peaks shared between 2 samples. Is non-symmetric, because different possible overlappings between peaks of the two sets.
#######  

def Find_rep_peaks_pairwise(rep_data, keys_list1, keys_list2, peaks_shared_ar):
    peaks_ar1=rep_data[keys_list1] 
    peaks_ar2=rep_data[keys_list2]
    #Number of peaks shared.
    num_peak12_shared=0
    num_peak21_shared=0
    
    #Count the number of peaks of sample 1 shared with sample 2.
    for peak1 in peaks_ar1:
        shared=0
        for shared_peak in peaks_shared_ar:
            shared_peak_med=int((shared_peak[0]+shared_peak[1])/2)
            if peak1[1]>=shared_peak_med>=peak1[0]:
                shared=1
        if shared==1:
            num_peak12_shared+=1
                
    #Count the number of peaks of sample 2 shared with sample 1.
    for peak2 in peaks_ar2:
        shared=0
        for shared_peak in peaks_shared_ar:
            shared_peak_med=int((shared_peak[0]+shared_peak[1])/2)
            if peak2[1]>=shared_peak_med>=peak2[0]:
                shared=1
        if shared==1:
            num_peak21_shared+=1            
    print("Number of peaks of sample ", keys_list1, " shared with peaks of sample ", keys_list2, ": ", num_peak12_shared)
    print("Number of peaks of sample ", keys_list2, " shared with peaks of sample ", keys_list1, ": ", num_peak21_shared)
    return num_peak12_shared, num_peak21_shared 


#################
### Model peaks random placement.
#################

def model_genome_place_peaks(genome_len, peak_type, bin_info_dict):
    Generic_tracks_dict={}
    names_list=[]
    Sim_peaks_dict={}
    
    #Place peaks.
    if peak_type=='Generic':
        for name, peak_info in bin_info_dict.items():
            #Generate genome-wide array.
            genome_model=np.array([0]*genome_len)
            #Keep coordinates of simulated peaks.
            sim_peaks_ar=[]
            for i in range(peak_info[0]):
                test_place=0
                while test_place!=3:
                    #Choose a random position of a peak start.
                    peak_start_position=rd.randrange(0, genome_len-1, 1)
                    if (peak_start_position+peak_info[1])<genome_len: #Check that peak is within the range of a genome.
                        if (genome_model[peak_start_position]==1) or (genome_model[peak_start_position+peak_info[1]]==1):
                            test_place=2 #Region is already occupied by a peak or out of the genome range.
                        elif (genome_model[peak_start_position]==0) and (genome_model[peak_start_position+peak_info[1]]==0):
                            test_place=3 #Region is free for peak placement.
                    else:
                        test_place=1
                #Place a peak.
                genome_model[peak_start_position:peak_start_position+peak_info[1]]=np.array([1]*peak_info[1])
                sim_peaks_ar.append([peak_start_position, peak_start_position+peak_info[1]])
                
            Generic_tracks_dict[name]=genome_model
            names_list.append(name)
            Sim_peaks_dict[name]=sim_peaks_ar
    
    elif peak_type=='Real':
        for name, peak_info in bin_info_dict.items():
            #Generate genome-wide array.
            genome_model=np.array([0]*genome_len)
            #Keep coordinates of simulated peaks.
            sim_peaks_ar=[]            
            for i in peak_info:
                peak_width=i[1]-i[0]
                test_place=0
                while test_place!=3:
                    #Choose a random position of a peak start.
                    peak_start_position=rd.randrange(0, genome_len-1, 1)
                    if (peak_start_position+peak_width)<genome_len: #Check that peak is within the range of a genome.
                        if np.sum(genome_model[peak_start_position:peak_start_position+peak_width+1])>0:
                            test_place=2 #Region is already occupied by a peak or out of the genome range.
                        else:
                            test_place=3 #Region is free for peak placement.
                    else:
                        test_place=1
                #Place a peak.
                genome_model[peak_start_position:peak_start_position+peak_width]=np.array([1]*peak_width)
                sim_peaks_ar.append([peak_start_position, peak_start_position+peak_width])
    
            Generic_tracks_dict[name]=genome_model
            names_list.append(name)     
            Sim_peaks_dict[name]=sim_peaks_ar
    
    #Overlay two simulated genome tracks.
    Track_overlay=Generic_tracks_dict[names_list[0]] + Generic_tracks_dict[names_list[1]] #Numpy arrays much faster for operations with long arrays.
    Track_overlay_list=Track_overlay.tolist()
    Generic_tracks_dict[names_list[0]+'_vs_'+names_list[1]]=Track_overlay_list
    
    names_list.append(names_list[0]+'_vs_'+names_list[1])
    
    return Generic_tracks_dict, Track_overlay_list, names_list, Sim_peaks_dict


#################
### Simulate experiments. Add statistics.
#################

def run_simulation(input_data, genome_len, deletion_correction, peak_type, Number_of_simulations, output_path):
    
    #Read data, get number of peaks and mean peak width.
    bin_info_dict, Real_peaks_dict=read_data_bin_width(input_data)
    
    #Run Monte-Carlo simulations.
    #Keep the number of overlaying peaks.
    Overlay_peaks_num=[]
    #Keep the number of peaks of a set1 and set2 in overlay.
    Sets_overlaying_peaks=[[], []]
    for i in range(Number_of_simulations):
        print(f'Number of iteration: {i}')
        #Place peaks randomly.
        if peak_type=="Real":
            Generic_tracks_dict, Track_overlay_list, names_list, Sim_peaks_dict=model_genome_place_peaks(genome_len-deletion_correction, peak_type, Real_peaks_dict) #!Genome with corrected length!
        elif peak_type=="Generic":
            Generic_tracks_dict, Track_overlay_list, names_list, Sim_peaks_dict=model_genome_place_peaks(genome_len-deletion_correction, peak_type, bin_info_dict) #!Genome with corrected length!
        #Identify ovelaying peaks.
        Overlay_peaks=Find_rep_peaks(Track_overlay_list, 2)
        Overlay_peaks_num.append(len(Overlay_peaks))
        #Number of peaks of the set 1 in overlay and number of peaks of the set 2 in overlay.
        Num_peaks12_shared_sim, Num_peaks21_shared_sim=Find_rep_peaks_pairwise(Sim_peaks_dict, names_list[0], names_list[1], Overlay_peaks) 
        Sets_overlaying_peaks[0].append(Num_peaks12_shared_sim)
        Sets_overlaying_peaks[1].append(Num_peaks21_shared_sim)
    
    print(Overlay_peaks_num)
    print(Sets_overlaying_peaks[0])
    print(Sets_overlaying_peaks[1])
    
    ##Number of overlaying peaks.   
    #Get real data.
    Real_overlay_peaks=overlap_call(Real_peaks_dict, 2, genome_len)
    Real_overlay_peaks_num=len(Real_overlay_peaks)
    
    #Write simulated values for further processing.
    fileout=open(f'{output_path}{names_list[2]}_peaks_overlay_{Number_of_simulations}_it_{peak_type}_delcor.txt', 'w+')
    fileout.write(str(Overlay_peaks_num))
    fileout.write(str(Real_overlay_peaks_num))
    fileout.close()
    
    ##Numbers of set1 and set2 peaks in overlay.  
    #Get real data.
    Num_peaks12_shared_real, Num_peaks21_shared_real=Find_rep_peaks_pairwise(Real_peaks_dict, names_list[0], names_list[1], Real_overlay_peaks) 
    
    #Write simulated values for further processing.
    fileout=open(f'{output_path}{names_list[0]}_peaks_in_overlay_{Number_of_simulations}_it_{peak_type}_delcor.txt', 'w+')
    fileout.write(str(Sets_overlaying_peaks[0]))
    fileout.write(str(Num_peaks12_shared_real))
    fileout.close()   
    
    fileout=open(f'{output_path}{names_list[1]}_peaks_in_overlay_{Number_of_simulations}_it_{peak_type}_delcor.txt', 'w+')
    fileout.write(str(Sets_overlaying_peaks[1]))
    fileout.write(str(Num_peaks21_shared_real))
    fileout.close()  
    
    #Assemble dict of calculated data.
    simulated_data_dict={names_list[2] : [Real_overlay_peaks_num, Overlay_peaks_num],
                         names_list[0] : [Num_peaks12_shared_real, Sets_overlaying_peaks[0]],
                         names_list[1] : [Num_peaks21_shared_real, Sets_overlaying_peaks[1]]}
    
    return simulated_data_dict


#################
### Analyse data. Add statistics.
#################

def read_precalc_data(simulated_data_inpath):
    #Read pre-simulated data.
    filein=open(simulated_data_inpath, 'r')
    for line in filein:
        line=line.rstrip().split(']')
        Real_overlay_peaks_num=int(line[1])
        Overlay_peaks_num=list(map(int, line[0].lstrip('[').split(', '))) #Stolen from https://stackoverflow.com/questions/3371269/call-int-function-on-every-list-element
    
    return Real_overlay_peaks_num, Overlay_peaks_num


def analyse_simulated_data(input_data, genome_len, deletion_correction, peak_type, Number_of_simulations, analysis_regime, simulated_data_dict_in, output_path):
    
    #Read data, get number of peaks and mean peak width.
    bin_info_dict, Real_peaks_dict=read_data_bin_width(input_data)
    
    #Simulate data de-novo or read data with pre-calculations.
    if analysis_regime=='Simulate':
        simulated_data_dict=run_simulation(input_data, genome_len, deletion_correction, peak_type, Number_of_simulations, output_path)
        
    elif analysis_regime=='Read':
        simulated_data_dict={}
        for name, simulated_data_inpath in simulated_data_dict_in.items():
            Real_overlay_peaks_num, Sim_overlay_peaks_num=read_precalc_data(simulated_data_inpath)
            print(name, Real_overlay_peaks_num, Sim_overlay_peaks_num)
            simulated_data_dict[name]=[Real_overlay_peaks_num, Sim_overlay_peaks_num]
    
    #Get disctionary keys.
    names_list=list(simulated_data_dict.keys())
    
    ## Number of overlaying peaks.    
    #Empirical average and variance are computed. Stolen from https://stackoverflow.com/questions/7805552/fitting-a-histogram-with-python
    Overlay_peaks_num=simulated_data_dict[names_list[0]][1]
    peaks_num_av=np.mean(Overlay_peaks_num)
    peaks_num_var=np.var(Overlay_peaks_num)
    #From that, we know the shape of the fitted Gaussian.
    pdf_x=np.linspace(np.min(Overlay_peaks_num)-100,np.max(Overlay_peaks_num)+100,200)
    pdf_y=norm.pdf(pdf_x, peaks_num_av, np.sqrt(peaks_num_var))
    
    #Calculate p-value.
    Real_overlay_peaks_num=simulated_data_dict[names_list[0]][0]
    p_value=1-norm.cdf(Real_overlay_peaks_num, peaks_num_av, np.sqrt(peaks_num_var))
    print(f'Simulated distribution of {names_list[0]} overlaying peaks: norm({np.round(peaks_num_av,1)}, {np.round(peaks_num_var,1)})')
    print(f'Overlaying peaks detected={Real_overlay_peaks_num}; p-value={p_value}')
    
    #Plot distribution of the number of overlaying peaks.
    plt.hist(Overlay_peaks_num, bins=int(np.sqrt(Number_of_simulations)/4), density=True, facecolor='#ff9088', edgecolor='black', linewidth=0.1)
    plt.axvline(x=Real_overlay_peaks_num, color='r', linestyle='dashed', linewidth=2)
    plt.plot(pdf_x, pdf_y, 'b--', linewidth=2)   
    plt.annotate(f'Iterations: {Number_of_simulations}', xy=(peaks_num_av+30, 0.038), xycoords='data', size=12)
    plt.annotate(f'{names_list[1]}: {bin_info_dict[names_list[1]][0]} peaks {bin_info_dict[names_list[1]][1]} bp', xy=(peaks_num_av+30, 0.035), xycoords='data', size=12)
    plt.annotate(f'{names_list[2]}: {bin_info_dict[names_list[2]][0]} peaks {bin_info_dict[names_list[2]][1]} bp', xy=(peaks_num_av+30, 0.032), xycoords='data', size=12)
    plt.annotate(f'$\mu$={np.round(peaks_num_av,1)}', xy=(peaks_num_av+30, 0.029), xycoords='data', size=12)
    plt.annotate(f'$\sigma^{2}$={np.round(peaks_num_var,1)}', xy=(peaks_num_av+30, 0.026), xycoords='data', size=12)
    plt.annotate(f'Real peaks in overlay: {Real_overlay_peaks_num}', xy=(Real_overlay_peaks_num-13, 0.004), xycoords='data', size=12, rotation=90)
    plt.annotate(f'p-value: {"{:.1e}".format(p_value)}', xy=(Real_overlay_peaks_num+5, 0.004), xycoords='data', size=12, rotation=90)    
    plt.ylabel('Probability', size=13)
    plt.xlabel('Number of peaks in overlay', size=13)
    plt.xlim([np.min(Overlay_peaks_num)-30, Real_overlay_peaks_num+30])
    plt.show()
    plt.savefig(f'{output_path}{names_list[0]}_peaks_overlay_{Number_of_simulations}_it_{peak_type}_delcor.png', dpi=300, figsize=(5, 3))
    plt.close()
    
    
    ## Numbers of set1 and set2 peaks in overlay.   
    #Empirical average and variance are computed. Stolen from https://stackoverflow.com/questions/7805552/fitting-a-histogram-with-python
    Peaks12_in_overlay=simulated_data_dict[names_list[1]][1]
    Peaks21_in_overlay=simulated_data_dict[names_list[2]][1]
    #Empirical average and variance are computed. Stolen from https://stackoverflow.com/questions/7805552/fitting-a-histogram-with-python
    peaks12_num_av=np.mean(Peaks12_in_overlay)
    peaks12_num_var=np.var(Peaks12_in_overlay)
    peaks21_num_av=np.mean(Peaks21_in_overlay)
    peaks21_num_var=np.var(Peaks21_in_overlay)
    #From that, we know the shape of the fitted Gaussian.
    pdf_x12=np.linspace(np.min(Peaks12_in_overlay)-100,np.max(Peaks12_in_overlay)+100,200)
    pdf_y12=norm.pdf(pdf_x12, peaks12_num_av, np.sqrt(peaks12_num_var))
    pdf_x21=np.linspace(np.min(Peaks21_in_overlay)-20,np.max(Peaks21_in_overlay)+20,100)  
    pdf_y21=norm.pdf(pdf_x21, peaks21_num_av, np.sqrt(peaks21_num_var))
    
    #Calculate p-value.
    Real_peaks12_in_overlay_num=simulated_data_dict[names_list[1]][0]
    Real_peaks21_in_overlay_num=simulated_data_dict[names_list[2]][0]
    p_value12=1-norm.cdf(Real_peaks12_in_overlay_num, peaks12_num_av, np.sqrt(peaks12_num_var))
    p_value21=1-norm.cdf(Real_peaks21_in_overlay_num, peaks21_num_av, np.sqrt(peaks21_num_var))
    print(f'Simulated distribution of {names_list[1]} peaks in overlay: norm({np.round(peaks12_num_av,1)}, {np.round(peaks12_num_var,1)})')
    print(f'{names_list[1]} overlaying peaks detected={Real_peaks12_in_overlay_num}; p-value={p_value12}')
    print(f'Simulated distribution of {names_list[2]} peaks in overlay: norm({np.round(peaks21_num_av,1)}, {np.round(peaks21_num_var,1)})')
    print(f'{names_list[2]} overlaying peaks detected={Real_peaks21_in_overlay_num}; p-value={p_value21}')    
    
    #Plot distribution of the number of peaks of set1 and set2 in overlay.
    plt.hist(Peaks12_in_overlay, bins=int(np.sqrt(Number_of_simulations)/4), density=True, facecolor='#ff9088', edgecolor='black', linewidth=0.1)
    plt.axvline(x=Real_peaks12_in_overlay_num, color='r', linestyle='dashed', linewidth=2)
    plt.plot(pdf_x12, pdf_y12, 'b--', linewidth=2)   
    plt.annotate(f'Iterations: {Number_of_simulations}', xy=(peaks12_num_av-45, 0.12), xycoords='data', size=12)
    plt.annotate(f'{names_list[1]}: {bin_info_dict[names_list[1]][0]} peaks {bin_info_dict[names_list[1]][1]} bp', xy=(peaks12_num_av-45, 0.11), xycoords='data', size=11)
    plt.annotate(f'{names_list[2]}: {bin_info_dict[names_list[2]][0]} peaks {bin_info_dict[names_list[2]][1]} bp', xy=(peaks12_num_av-45, 0.10), xycoords='data', size=11)
    plt.annotate(f'$\mu$={np.round(peaks12_num_av,1)}', xy=(peaks12_num_av-45, 0.09), xycoords='data', size=12)
    plt.annotate(f'$\sigma^{2}$={np.round(peaks12_num_var,1)}', xy=(peaks12_num_av-45, 0.08), xycoords='data', size=12) 
    plt.annotate(f'{names_list[1]} real peaks in overlay: {Real_peaks12_in_overlay_num}', xy=(Real_peaks12_in_overlay_num-3, 0.02), xycoords='data', size=12, rotation=90)
    plt.annotate(f'p-value: {"{:.1e}".format(p_value12)}', xy=(Real_peaks12_in_overlay_num+1, 0.02), xycoords='data', size=12, rotation=90)    
    plt.ylabel('Probability', size=13)
    plt.xlabel(f'Number of {names_list[1]} peaks in overlay', size=13)
    plt.xlim([np.min(Peaks12_in_overlay)-30, Real_peaks12_in_overlay_num+10])
    plt.show()
    plt.savefig(f'{output_path}{names_list[1]}_peaks_in_overlay_{Number_of_simulations}_it_{peak_type}_delcor.png', dpi=300, figsize=(5, 3))    
    plt.close()
    
    #Plot distribution of the number of peaks of set1 and set2 in overlay.
    plt.hist(Peaks21_in_overlay, bins=int(np.sqrt(Number_of_simulations)/5), density=True, facecolor='#ff9088', edgecolor='black', linewidth=0.1)
    plt.axvline(x=Real_peaks21_in_overlay_num, color='r', linestyle='dashed', linewidth=2)
    plt.plot(pdf_x21, pdf_y21, 'b--', linewidth=2)   
    plt.annotate(f'Iterations: {Number_of_simulations}', xy=(peaks21_num_av+30, 0.033), xycoords='data', size=12)
    plt.annotate(f'{names_list[1]}: {bin_info_dict[names_list[1]][0]} peaks {bin_info_dict[names_list[1]][1]} bp', xy=(peaks21_num_av+30, 0.030), xycoords='data', size=11)
    plt.annotate(f'{names_list[2]}: {bin_info_dict[names_list[2]][0]} peaks {bin_info_dict[names_list[2]][1]} bp', xy=(peaks21_num_av+30, 0.027), xycoords='data', size=11)
    plt.annotate(f'$\mu$={np.round(peaks21_num_av,1)}', xy=(peaks21_num_av+30, 0.024), xycoords='data', size=12)
    plt.annotate(f'$\sigma^{2}$={np.round(peaks21_num_var,1)}', xy=(peaks21_num_av+30, 0.021), xycoords='data', size=12)
    plt.annotate(f'{names_list[2]} real peaks in overlay: {Real_peaks21_in_overlay_num}', xy=(Real_peaks21_in_overlay_num-13, 0.004), xycoords='data', size=12, rotation=90)
    plt.annotate(f'p-value: {"{:.1e}".format(p_value21)}', xy=(Real_peaks21_in_overlay_num+5, 0.004), xycoords='data', size=12, rotation=90)    
    plt.ylabel('Probability', size=13)
    plt.xlabel(f'Number of {names_list[2]} peaks in overlay', size=13)
    plt.xlim([np.min(Peaks21_in_overlay)-30, Real_peaks21_in_overlay_num+30])
    plt.show()
    plt.savefig(f'{output_path}{names_list[2]}_peaks_in_overlay_{Number_of_simulations}_it_{peak_type}_delcor.png', dpi=300, figsize=(5, 3))    
    plt.close()    
        
    return


#analyse_simulated_data(Input_dict, Genome_len, Deletion_correction, 'Real', 10000, 'Read', Simulated_data_dict, Outpath)