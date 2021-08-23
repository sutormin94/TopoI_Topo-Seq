###############################################
##Dmitry Sutormin, 2021##
##TopoI Topo-Seq analysis##

#The script takes TCSs from BroadPeak file. Extracts 
# and returns sequences under the TCSs, then plots sequence motifs.
#Also it writes sequences and motif to files.
###############################################

#######
#Packages to be imported.
#######

import os
import numpy as np
from scipy.stats import binom_test
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from Bio import SeqIO
from Bio import SeqUtils
from Bio.SeqUtils import GC as GC_count
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import AlignIO, motifs
from Bio.Align import AlignInfo
import weblogo
from weblogo import LogoOptions, LogoFormat, eps_formatter, read_seq_data, LogoData, png_formatter, pdf_formatter

#######
#Variables to be defined.
#######



print('Variables to be defined:')

#PWD
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\\"

#Input data - WIG.
TCSs_set_name='TopoI_Ara'
path_to_wig_files={'F' : PWD + "WIG_NE_strand_specific_masked_scaled_av_masked_accB_subtract_mock\TopoI_Ara_N3E_F_masked_scaled_av_123_mock_subtr.wig",    'R' : PWD + "WIG_NE_strand_specific_masked_scaled_av_masked_accB_subtract_mock\TopoI_Ara_N3E_R_masked_scaled_av_123_mock_subtr.wig"}
TCSs_set_name_1='TopoI_no_Ara'
path_to_wig_files_1={'F' : PWD + "WIG_NE_strand_specific_masked_scaled_av_masked_accB_subtract_mock\TopoI_N3E_F_masked_scaled_av_123_mock_subtr.wig", 'R' : PWD + "WIG_NE_strand_specific_masked_scaled_av_masked_accB_subtract_mock\TopoI_N3E_R_masked_scaled_av_123_mock_subtr.wig"}

#Input data - TSCs coordinates, BroadPeak.
All_TCSs_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\EcTopoI_ChIP_Seq_vs_Topo_Seq\\TopoI_Ara_TCSs_called_15.BroadPeak"
Unique_TCSs_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\EcTopoI_ChIP_Seq_vs_Topo_Seq\\EcTopoI_ChIP_shared_with_EcTopoI_Topo.narrowPeak"

#Path to the genome FASTA.
Genome_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\Scripts_TopoI_Topo-seq\TopoI_Topo-Seq\Additional_genome_features\E_coli_w3110_G_Mu.fasta"

#Genome ID.
Genome_ID='NC_007779.1_w3110_Mu'

#Path for the output.
Output_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\EcTopoI_ChIP_Seq_vs_Topo_Seq\\TEST_Shared_TCSs_signal_motif\\"
if not os.path.exists(Output_path):
    os.makedirs(Output_path)
    
    
    
#######
#Parses WIG file with N3/5E values.
#Computes a total number of Ends.
#######

def wig_parsing(wigfile):
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    NE_values=[]
    Total_NE=0
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(float(line[0]))
            Total_NE+=float(line[0])
    print('Total number of ends: ' + str(Total_NE))
    wigin.close()
    return NE_values, Total_NE
    
#######
#Opens and reads BroadPeak file with TCSs coordinates.
#######

def read_TCSs(TCSs_path):
    TCSs_dict={}
    filein=open(TCSs_path, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        TCSs_dict[int(line[1])]=[int(line[2]), line[3]]
    filein.close()
    return TCSs_dict

#######
#Return strandness of TCSs, classify TCSs by strand.
#######

def get_TCSs_strandness(Unique_TCSs_dict, All_TCSs_dict):
    #Get strandness of unique TCSs.
    TCSs_F_ar=[]
    TCSs_R_ar=[]
    for TCSs_coord, TCSs_info in Unique_TCSs_dict.items():
        if TCSs_coord in All_TCSs_dict:
            strand=All_TCSs_dict[TCSs_coord][1].split('_')[2]
            if strand=='F':
                TCSs_F_ar.append(TCSs_coord-1)
            elif strand=='R':
                TCSs_R_ar.append(TCSs_coord)
        else:
            print(f'{TCSs_coord} not found!')
    return TCSs_F_ar, TCSs_R_ar

#######
#Genome sequence parsing.
#######

def genome_seq(genome_path):
    genome=open(genome_path, 'r')
    for record in SeqIO.parse(genome, "fasta"):
        genome_sequence=str(record.seq)
    genome.close()
    print('Whole genome average GC: ' + str(SeqUtils.GC(genome_sequence)))
    print('Whole genome length: ' + str(len(genome_sequence)))        
    return genome_sequence

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
#Detect TCSs.
#######

def get_TCSs_signal(wig_data, TCSs_ar, win_range_local):
    
    TCSs_list=[]
    TCSs_metasignal=np.array([0.0]*sum(win_range_local))
    for i in TCSs_ar:
        TCSs_list.append([i, wig_data[i]])
        cleavage_region_wig=wig_data[i-win_range_local[0]:i+win_range_local[1]]
        TCSs_metasignal+=np.array(cleavage_region_wig)
        
    return TCSs_list, TCSs_metasignal


#######
#Returns list of DNA seqs under the GCSs. Seqs have len=win_width.
#Writes sequences under the GCSs to file.
#######

def return_seqs(TCS_coords, Seq_type, win_range, win_range_local, genomefa, filepath_full_len, filepath_nbp_LOGO): 
    fileout=open(filepath_full_len, 'w')
    fileout_nbp_LOGO=open(filepath_nbp_LOGO, 'w')
    seqs={}
    seqs_nbp={}
    for i in range(len(TCS_coords)):
        if (TCS_coords[i][0]-win_range[0]-1>0) and (TCS_coords[i][0]+win_range[1]-1)<len(genomefa):
            seq=genomefa[int(TCS_coords[i][0]- win_range[0] - 1):int(TCS_coords[i][0]+ win_range[1] - 1)]
            seq_nbp_LOGO=genomefa[int(TCS_coords[i][0]-win_range_local[0]):int(TCS_coords[i][0]+win_range_local[1])]
            seqs[str(TCS_coords[i][0])+'_'+Seq_type]=seq
            seqs_nbp[str(TCS_coords[i][0])+'_'+Seq_type]=seq_nbp_LOGO
            fileout.write('>'+str(TCS_coords[i][0])+'\n'+str(seq)+'\n')
            fileout_nbp_LOGO.write('>'+str(TCS_coords[i][0])+'\n'+str(seq_nbp_LOGO)+'\n')
    fileout.close()
    fileout_nbp_LOGO.close()
    print('Number of sequences (TCSs) analysing: ' + str(len(seqs)))
    return seqs, seqs_nbp

#######
#Prepare list of reverse-complement sequences.
#######

def RC_merge_R(F_sequences_dict, R_sequences_dict, filepath_full_len):

    fileout=open(filepath_full_len, 'w')        

    R_rc_sequences_dict={}
    F_R_rc_sequences_dict={}

    for seq_name, seq in R_sequences_dict.items():
        seq=Seq(seq)
        seq_rc=seq.reverse_complement()   
        R_rc_sequences_dict[seq_name]=str(seq_rc)
        F_R_rc_sequences_dict[seq_name]=str(seq_rc)
        fileout.write('>'+str(seq_name)+'\n'+str(str(seq_rc))+'\n')

    for seq_name, seq in F_sequences_dict.items():
        F_R_rc_sequences_dict[seq_name]=seq
        fileout.write('>'+str(seq_name)+'\n'+str(seq)+'\n')

    fileout.close()

    return R_rc_sequences_dict, F_R_rc_sequences_dict

#######
#Converts dictionary to list.
#######

def dict_to_list(indict):

    outlist=[]

    for name, data in indict.items():
        outlist.append(data)

    return outlist

#######
#PFM construction.
#Scans sequences stack by columns, counts the number of particular letters.
#Returns a range of PFMs - "positional frequencies matrixes" .
#######

def make_PFM(seqs_list):
    matrix=[]
    template=seqs_list[0]
    for i in range(len(template)):
        column=[0, 0, 0, 0]
        for j in range(len(seqs_list)):
            #print(seqs_list[j])
            if seqs_list[j][i] == str('A'):
                column[0] = column[0] + 1
            elif seqs_list[j][i] == str('T'):
                column[1] = column[1] + 1
            elif seqs_list[j][i] == str('G'):
                column[2] = column[2] + 1
            elif seqs_list[j][i] == str('C'):
                column[3] = column[3] + 1
        matrix.append(column)

    #Returns a range of PFMs.
    GC_percent = []
    GT_percent = []
    CT_percent = []
    A_percent = []
    T_percent = []
    G_percent = []
    C_percent = []
    for i in range(len(matrix)):
        GC = float((int(matrix[i][2]) + int(matrix[i][3]))) / (
                    int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        GT = float((int(matrix[i][1]) + int(matrix[i][2]))) / (
                    int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        CT = float((int(matrix[i][1]) + int(matrix[i][3]))) / (
                    int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        A = float((int(matrix[i][0]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        T = float((int(matrix[i][1]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        G = float((int(matrix[i][2]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        C = float((int(matrix[i][3]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        GC_percent.append(GC)
        GT_percent.append(GT)
        CT_percent.append(CT)
        A_percent.append(A)
        T_percent.append(T)
        G_percent.append(G)
        C_percent.append(C)
    return {'Num_seqs': len(seqs_list), 'A': A_percent, 'T': T_percent, 'G': G_percent, 'C': C_percent, 'CT': CT_percent, 'GT': GT_percent, 'GC': GC_percent}

#######
#Writes PFM data to file.
#######

def write_motif(ar, filepath, coord_shift):
    fileout=open(filepath, 'w')
    fileout.write("#X\tY\n")
    for i in range(len(ar)):
        fileout.write(str((-coord_shift/2)+1+i) + '\t' + str(ar[i])+'\n')
    fileout.close()
    return

#######
#Plotting the average cleavage signal (N3E) over the cleavage sites.
#######

def Plotting_N3E_metasignal(metasignal, num_seq, title, write_out, win_width):

    #Plotting   
    x_axis=[]
    for i in range(len(metasignal)):
        x_axis.append(-(win_width/2)+i) 

    print('len(x_axis)=' + str(len(x_axis)))
    #xticks_list=list(range(-11,0,2))+list(range(1,10,2))
    xticks_list=[-11,'',-9,'',-7,'',-5,'',-3,'',-1,1,'',3,'',5,'',7,'',9]

    plt.figure(dpi=100, figsize=(7, 3))
    plot1 = plt.subplot()          
    plot1.bar(x_axis, metasignal, width=0.8, align='center', color='cyan', edgecolor='#3f7570', linewidth=1.5)
    plot1.annotate(f'F strand: {num_seq[0]} sites\nR strand: {num_seq[1]} sites', xycoords='axes fraction', xytext=(0.65, 0.8), xy=(0.65, 0.8), color='k', weight="bold", size=12)   
    plot1.set_xlabel('Position, nt', size=17)
    plot1.set_ylabel('TopoI Topo-Seq N3E', size=17)                 
    plot1.tick_params(axis='both', direction='out', bottom='on', top=False, left='on', right=False, labelsize=15)
    plot1.set_xticks(x_axis)   
    plot1.set_xticklabels(xticks_list)    
    plot1.set_xlim(-win_width/2-0.5, win_width/2-0.5)
    plot1.set_ylim(0, 25)
    plt.tight_layout()
    plt.savefig(write_out, dpi=400, figsize=(7, 3)) 
    plt.close()

    return

#######
#Plotting the motif with statistic.
#Matrix type - type of the PFM to plot.
#######

def Plotting_stat(GC_PFM, num_seq, title, matrix_type, genome_sequence, write_out, win_width):
    #GC statistics module
    #Counts average GC% over the whole genome
    GC_genome=GC_count(genome_sequence)/100
    print('GC% of the reference genome: ' + str(GC_genome))

    #Counts GC% p-value in the particular pwm column.
    #Returns p-value array and auxiliary Zero array for plotting.
    alignment_thick=num_seq
    pValue=[]
    Zero=[]
    for i in range(len(GC_PFM)):
        pValue.append(binom_test(float(GC_PFM[i]) * alignment_thick, n=alignment_thick, p=GC_genome))
        Zero.append(1)
    #Plotting   
    x_axis=[]
    for i in range(len(GC_PFM)):
        x_axis.append(-(win_width/2)+1+i)      
    print('len(x_axis)=' + str(len(x_axis)))
    ax_range = [-win_width/2, win_width/2, 0.35, 0.9]
    plt.figure(dpi=100, figsize=(16, 6))
    plt.suptitle(str(title), fontsize=20)
    plot1 = plt.subplot()
    plot1.set_xticks([0], minor=True)
    plot1.xaxis.grid(True, which='minor', linewidth=0.5, linestyle='--', alpha=1)            
    #GC% pwm plotting
    plot1.plot(x_axis, GC_PFM, color='green', linewidth=1)
    plot1.plot(x_axis, GC_PFM, 'go', markersize=3)
    plot1.axis(ax_range) 
    plot1.annotate(matrix_type+'%', xytext=(65, 0.65), xy=(40, 0.85), color='green', weight="bold", size=15)
    txt=plot1.annotate('p-value', xytext=(80, 0.60), xy=(-105, 0.64), color='cyan', weight="bold", size=15)
    txt.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='black')])   
    plot1.set_xlabel('Position, nt', size=17)
    plot1.set_ylabel(matrix_type+'%', size=17)                 
    #Set axis parameters
    plot1.tick_params(axis='both', direction='in', bottom='on', top='on', left='on', right='on')
    plot1.set_ylim(0.0, 1.0)
    plot1.set_xlim(-win_width/2, win_width/2)
    plot1.set_xticks(np.concatenate((np.arange(-(win_width/2)+5, (win_width/2)+2, 10), [0, 3, -63, -17, 20, 66])))
    #p-value plotting
    plot2=plot1.twinx()
    plot2.plot(x_axis, pValue, 'k', linewidth=0.5, alpha=0.6)
    plot2.fill_between(x_axis, pValue, Zero, color='cyan', alpha=0.2)
    plot2.set_yticks(np.arange(0, 1.01, 0.01), minor=False)
    plot2.set_yscale('log')
    plot2.set_yticks([0.005], minor=True)
    plot2.yaxis.grid(True, which='minor', linewidth=1, linestyle='--', alpha=0.8)
    plot2.annotate('Confidence level = 0.005', xytext=(45, 0.0025), xy=(40, 0.8), color='black', size=15)
    plot2.set_ylim(0.000000001, 1.0)
    plot2.set_xlim(-win_width/2, win_width/2)
    plot2.set_ylabel('p-value, logarithmic scale', size=17)        
    #plt.show()
    plt.savefig(write_out, dpi=400, figsize=(16, 6)) 
    plt.close()
    return

#######
#Takes fasta file containing sequences under the cleavage sites.
#Returns consensus sequence and plots motif logo.
#######

def get_consensus_plot_motif(output_path, datapath):

    #Read mfa file, identify consensus sequences.
    alignment=AlignIO.read(output_path, "fasta")
    alignment_summary=AlignInfo.SummaryInfo(alignment)
    consensus=alignment_summary.dumb_consensus(threshold=0.35,  ambiguous='X')
    print('Consensus sequence:\n' + consensus)
    consensus_rc=consensus.reverse_complement()
    print('Reverse-complement consensus sequence:\n' + consensus_rc)
    print('Done!')

    #Read mfa file, draw motif. + strand.
    MFA_data=open(output_path)
    MFA_seqs=read_seq_data(MFA_data)
    logodata=LogoData.from_seqs(MFA_seqs)
    logooptions=LogoOptions(yaxis_scale=1.5, pad_right=True, stacks_per_line=20)
    logooptions.show_errorbars=False
    logoformat=LogoFormat(logodata, logooptions)
    pdf=weblogo.logo_formatter.pdf_formatter(logodata, logoformat)
    logout=open(datapath + ".pdf", 'wb')
    logout.write(pdf)
    logout.close()  

    return consensus

#######
#Wraps all the functions together.
#######

def wrap_function(dict_of_wigs_path, genome_input_path, All_TCSs_path, Unique_TCSs_path, TCSs_set_name, Genome_ID, output_path):

    win_width=100
    win_range_F=[(win_width/2)-2, (win_width/2)+2]
    win_range_R=[(win_width/2)-3, (win_width/2)+1]
    win_width_local=20
    win_range_local_F=[10, 10]
    win_range_local_R=[9, 11]
    PFM_type='GC'
    plot_title='TopoI cleavage motif obtained with Topo-Seq'

    #Reads input data in wig files.
    wig_F_data=wig_parsing(dict_of_wigs_path['F'])[0]
    wig_R_data=wig_parsing(dict_of_wigs_path['R'])[0]
    dict_of_wigs={'F' : wig_F_data, 'R' : wig_R_data}

    #Extract genomic sequence.        
    genome_sequence=genome_seq(genome_input_path)

    #Read TCSs - TopoI cleavage sites.
    All_TCSs_dict=read_TCSs(All_TCSs_path)
    Unique_TCSs_dict=read_TCSs(Unique_TCSs_path)
    
    #Get strandness of unique TCSs.
    TCSs_F_ar, TCSs_R_ar=get_TCSs_strandness(Unique_TCSs_dict, All_TCSs_dict)
    
    #Return signal at TCSs.
    TCSs_F, TCSs_metasignal_F=get_TCSs_signal(dict_of_wigs['F'], TCSs_F_ar, win_range_local_F)
    TCSs_R, TCSs_metasignal_R=get_TCSs_signal(dict_of_wigs['R'], TCSs_R_ar, win_range_local_R)
    TCSs_metasignal_F_R_rc=TCSs_metasignal_F+np.array(TCSs_metasignal_R.tolist()[::-1])
    TCSs_metasignal_F_R_rc_scaled=TCSs_metasignal_F_R_rc/(len(TCSs_F)+len(TCSs_R))

    #Plot N3E signal around cleavage sites.
    Plotting_N3E_metasignal(TCSs_metasignal_F_R_rc_scaled, [len(TCSs_F), len(TCSs_R)], f'TopoI Topo-Seq N3E statistic for {TCSs_set_name}', f'{output_path}_TopoI_Topo_Seq_unique_TCSs_cleavage_signal_F_R_rc_{TCSs_set_name}.svg', win_width_local)

    #Return sequences, plot motif.                
    #Get sequences.
    F_sequences_dict, F_sequences_dict_nbp=return_seqs(TCSs_F, 'F', win_range_F, win_range_local_F, genome_sequence, f'{output_path}{TCSs_set_name}_sequences_under_unique_TCSs_{TCSs_set_name}_full_F.fasta', f'{output_path}{TCSs_set_name}_sequences_under_unique_TCSs_{sum(win_range_local_F)}bp_{TCSs_set_name}_LOGO_F.fasta')
    R_sequences_dict, R_sequences_dict_nbp=return_seqs(TCSs_R, 'R', win_range_R, win_range_local_R, genome_sequence, f'{output_path}{TCSs_set_name}_sequences_under_unique_TCSs_{TCSs_set_name}_full_R.fasta', f'{output_path}{TCSs_set_name}_sequences_under_unique_TCSs_{sum(win_range_local_R)}bp_{TCSs_set_name}_LOGO_R.fasta')

    #Reverse complement R sequences, merge F and R, write merged fasta.
    R_rc_sequences_dict, F_R_rc_Sequences_dict=RC_merge_R(F_sequences_dict, R_sequences_dict, f'{output_path}{TCSs_set_name}_sequences_under_unique_TCSs_{TCSs_set_name}_full_F_R_rc.fasta')
    R_rc_sequences_dict_nbp, F_R_rc_Sequences_dict_nbp=RC_merge_R(F_sequences_dict_nbp, R_sequences_dict_nbp, f'{output_path}{TCSs_set_name}_sequences_under_unique_TCSs_{sum(win_range_local_F)}bp_{TCSs_set_name}_LOGO_F_R_rc.fasta')

    #Convert dict to list.
    F_sequences_list=dict_to_list(F_sequences_dict)
    R_sequences_list=dict_to_list(R_sequences_dict)
    F_R_rc_Sequences_list=dict_to_list(F_R_rc_Sequences_dict)

    #Compute PWM.
    F_PFMs=make_PFM(F_sequences_list)
    R_PFMs=make_PFM(R_sequences_list)
    All_PFMs=make_PFM(F_R_rc_Sequences_list)

    #Get consensus sequence, plot motif using weblogo.
    get_consensus_plot_motif(f'{output_path}{TCSs_set_name}_sequences_under_unique_TCSs_{sum(win_range_local_F)}bp_{TCSs_set_name}_LOGO_F_R_rc.fasta', f'{output_path}{TCSs_set_name}_sequences_under_unique_TCSs_{sum(win_range_local_F)}bp_{TCSs_set_name}_LOGO_F_R_rc')  

    #Write motif.
    write_motif(F_PFMs[PFM_type],   f'{output_path}{TCSs_set_name}_unique_TCSs_GC_pfm_F.txt', win_width)
    write_motif(R_PFMs[PFM_type],   f'{output_path}{TCSs_set_name}_unique_TCSs_GC_pfm_R.txt', win_width)
    write_motif(All_PFMs[PFM_type], f'{output_path}{TCSs_set_name}_unique_TCSs_GC_pfm_F_R_rc.txt', win_width)

    #Plot nucleotides frequences around cleavage sites, add statistics.
    Plotting_stat(F_PFMs[PFM_type],   F_PFMs['Num_seqs'],   f'TopoI motif F_statistic for {TCSs_set_name}', PFM_type, genome_sequence, f'{output_path}{PFM_type}_TopoI_motif_unique_TCSs_statistic_F_{TCSs_set_name}.png', win_width)
    Plotting_stat(R_PFMs[PFM_type],   R_PFMs['Num_seqs'],   f'TopoI motif R_statistic for {TCSs_set_name}', PFM_type, genome_sequence, f'{output_path}{PFM_type}_TopoI_motif_unique_TCSs_statistic_R_{TCSs_set_name}.png', win_width)
    Plotting_stat(All_PFMs[PFM_type], All_PFMs['Num_seqs'], f'TopoI motif R_statistic for {TCSs_set_name}', PFM_type, genome_sequence, f'{output_path}{PFM_type}_TopoI_motif_unique_TCSs_statistic_F_R_rc_{TCSs_set_name}.png', win_width)

    return

wrap_function(path_to_wig_files, Genome_path, All_TCSs_path, Unique_TCSs_path, TCSs_set_name, Genome_ID, Output_path)


print('Script ended its work succesfully!')