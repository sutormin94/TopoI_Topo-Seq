###############################################
##Dmitry Sutormin, 2021##
##TopoI Topo-Seq analysis##

#The script detects TCSs in input WIG. Return broadPeak file with TCSs coordinates, extracts 
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
path_to_wig_files={'TopoI_Ara': {'F' : PWD + "WIG_NE_strand_specific_masked_FE_smoothed_av\TopoI_Ara_F_N3E_FE_av_123.wig",    'R' : PWD + "WIG_NE_strand_specific_masked_FE_smoothed_av\TopoI_Ara_R_N3E_FE_av_123.wig"},
                   'TopoI':     {'F' : PWD + "WIG_NE_strand_specific_masked_FE_smoothed_av\TopoI_no_Ara_F_N3E_FE_av_123.wig", 'R' : PWD + "WIG_NE_strand_specific_masked_FE_smoothed_av\TopoI_no_Ara_R_N3E_FE_av_123.wig"}}

path_to_wig_files_1={'TopoI_Ara_IP':      {'F' : PWD + "WIG_NE_strand_specific_masked\DSu_14_S14_edt_N3E_F.wig", 'R' : PWD + "WIG_NE_strand_specific_masked\DSu_14_S14_edt_N3E_R.wig"},
                     'TopoI_Ara_mock':    {'F' : PWD + "WIG_NE_strand_specific_masked\DSu_13_S13_edt_N3E_F.wig", 'R' : PWD + "WIG_NE_strand_specific_masked\DSu_13_S13_edt_N3E_R.wig"},
                     'TopoI_no_Ara_IP':   {'F' : PWD + "WIG_NE_strand_specific_masked\DSu_16_S16_edt_N3E_F.wig", 'R' : PWD + "WIG_NE_strand_specific_masked\DSu_16_S16_edt_N3E_R.wig"},
                     'TopoI_no_Ara_mock': {'F' : PWD + "WIG_NE_strand_specific_masked\DSu_15_S15_edt_N3E_F.wig", 'R' : PWD + "WIG_NE_strand_specific_masked\DSu_15_S15_edt_N3E_R.wig"}}
path_to_wig_files_2={'TopoI_Ara_IP':      {'F' : PWD + "WIG_NE_strand_specific_masked\DSu_18_S18_edt_N3E_F.wig", 'R' : PWD + "WIG_NE_strand_specific_masked\DSu_18_S18_edt_N3E_R.wig"},
                     'TopoI_Ara_mock':    {'F' : PWD + "WIG_NE_strand_specific_masked\DSu_17_S17_edt_N3E_F.wig", 'R' : PWD + "WIG_NE_strand_specific_masked\DSu_17_S17_edt_N3E_R.wig"},
                     'TopoI_no_Ara_IP':   {'F' : PWD + "WIG_NE_strand_specific_masked\DSu_20_S20_edt_N3E_F.wig", 'R' : PWD + "WIG_NE_strand_specific_masked\DSu_20_S20_edt_N3E_R.wig"},
                     'TopoI_no_Ara_mock': {'F' : PWD + "WIG_NE_strand_specific_masked\DSu_19_S19_edt_N3E_F.wig", 'R' : PWD + "WIG_NE_strand_specific_masked\DSu_19_S19_edt_N3E_R.wig"}}
path_to_wig_files_3={'TopoI_Ara_IP':      {'F' : PWD + "WIG_NE_strand_specific_masked\DSu_22_S22_edt_N3E_F.wig", 'R' : PWD + "WIG_NE_strand_specific_masked\DSu_22_S22_edt_N3E_R.wig"},
                     'TopoI_Ara_mock':    {'F' : PWD + "WIG_NE_strand_specific_masked\DSu_21_S21_edt_N3E_F.wig", 'R' : PWD + "WIG_NE_strand_specific_masked\DSu_21_S21_edt_N3E_R.wig"},
                     'TopoI_no_Ara_IP':   {'F' : PWD + "WIG_NE_strand_specific_masked\DSu_24_S24_edt_N3E_F.wig", 'R' : PWD + "WIG_NE_strand_specific_masked\DSu_24_S24_edt_N3E_R.wig"},
                     'TopoI_no_Ara_mock': {'F' : PWD + "WIG_NE_strand_specific_masked\DSu_23_S23_edt_N3E_F.wig", 'R' : PWD + "WIG_NE_strand_specific_masked\DSu_23_S23_edt_N3E_R.wig"}}



#Path to the genome FASTA.
Genome_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\Scripts_TopoI_Topo-seq\Additional_genome_features\E_coli_w3110_G_Mu.fasta"

#Signal threshold.
Threshold=10

#Genome ID.
Genome_ID='NC_007779.1_w3110_Mu'

#Path for the output.
Output_path=PWD + "TCS_motifs\\Replics_123_av_Thresholds\\"
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

def detect_TCSs(wig_data, threshold, win_range_local):
        TCSs_list=[]
        TCSs_metasignal=np.array([0.0]*sum(win_range_local))
        for i in range(len(wig_data)):
                if (i>win_range_local[0]) & (i<len(wig_data)-win_range_local[1]):
                        if wig_data[i]>threshold:
                                TCSs_list.append([i, wig_data[i]])
                                cleavage_region_wig=wig_data[i-win_range_local[0]:i+win_range_local[1]]
                                TCSs_metasignal+=np.array(cleavage_region_wig)
        return TCSs_list, TCSs_metasignal


#######
#Write TCSs coordinates in BroadPeak format.
#######

def write_TCSs_coords(TCSs_list_F, TCSs_list_R, Genome_ID, path_out):
        #Converts coordinates to 1-start system.
        TCSs_out=open(path_out, 'w')
        for i in range(len(TCSs_list_F)):
                TCSs_out.write(f'{Genome_ID}\t{TCSs_list_F[i][0]}\t{TCSs_list_F[i][0]+1}\tTCS_{i}_F\t10\t.\t{TCSs_list_F[i][1]}\t-1\t-1\n')  
        for j in range(len(TCSs_list_R)):
                TCSs_out.write(f'{Genome_ID}\t{TCSs_list_R[j][0]}\t{TCSs_list_R[j][0]+1}\tTCS_{j}_R\t10\t.\t{TCSs_list_R[j][1]}\t-1\t-1\n') 
        
        TCSs_out.close()        
        return

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

def wrap_function(dict_of_wigs_path, genome_input_path, threshold, Genome_ID, output_path):
        
        win_width=100
        win_range_F=[(win_width/2)-2, (win_width/2)+2]
        win_range_R=[(win_width/2)-3, (win_width/2)+1]
        win_width_local=20
        win_range_local_F=[10, 10]
        win_range_local_R=[9, 11]
        PFM_type='GC'
        plot_title='TopoI cleavage motif obtained with Topo-Seq'
        
        #Reads input data in wig files.
        dict_of_wigs={}
        for pair_name, pair_dict in dict_of_wigs_path.items():
                wig_F_data=wig_parsing(pair_dict['F'])[0]
                wig_R_data=wig_parsing(pair_dict['R'])[0]
                dict_of_wigs[pair_name]={'F' : wig_F_data, 'R' : wig_R_data}
        
        #Extract genomic sequence.        
        genome_sequence=genome_seq(genome_input_path)
        
        #Detect TCSs - TopoI cleavage sites.
        TCSs_dict={}
        for pair_name, pair_dict in dict_of_wigs.items():
                TCSs_F, TCSs_metasignal_F=detect_TCSs(pair_dict['F'], threshold, win_range_local_F)
                TCSs_R, TCSs_metasignal_R=detect_TCSs(pair_dict['R'], threshold, win_range_local_R)
                TCSs_metasignal_F_R_rc=TCSs_metasignal_F+np.array(TCSs_metasignal_R.tolist()[::-1])
                TCSs_metasignal_F_R_rc_scaled=TCSs_metasignal_F_R_rc/(len(TCSs_F)+len(TCSs_R))
                
                #Write TCSs coordinates into BroadPeak file.
                write_TCSs_coords(TCSs_F, TCSs_R, Genome_ID, f'{output_path}{pair_name}_TCSs_called_{threshold}.BroadPeak')
                
                #Plot N3E signal around cleavage sites.
                Plotting_N3E_metasignal(TCSs_metasignal_F_R_rc_scaled, [len(TCSs_F), len(TCSs_R)], 'TopoI Topo-Seq N3E statistic for ' + pair_name, output_path+'_TopoI_Topo_Seq_trusted_TCSs_cleavage_signal_'+str(threshold)+'_F_R_rc_'+str(pair_name)+'.svg', win_width_local)
                
                TCSs_dict[pair_name]={'F' : TCSs_F, 'R' : TCSs_R}
        
        #Return sequences, plot motif.
        dict_of_PFMs={}
        for pair_name, pair_dict in TCSs_dict.items():
                TCSs_F=pair_dict['F']
                TCSs_R=pair_dict['R']                
                
                #Get sequences.
                F_sequences_dict, F_sequences_dict_nbp=return_seqs(TCSs_F, 'F', win_range_F, win_range_local_F, genome_sequence, output_path+str(pair_name)+'_sequences_under_TCSs_full_'+str(threshold)+'_F.fasta', output_path+str(pair_name)+'_sequences_under_TCSs_'+str(sum(win_range_local_F))+'bp_LOGO_'+str(threshold)+'_F.fasta')
                R_sequences_dict, R_sequences_dict_nbp=return_seqs(TCSs_R, 'R', win_range_R, win_range_local_R, genome_sequence, output_path+str(pair_name)+'_sequences_under_TCSs_full_'+str(threshold)+'_R.fasta', output_path+str(pair_name)+'_sequences_under_TCSs_'+str(sum(win_range_local_R))+'bp_LOGO_'+str(threshold)+'_R.fasta')
                
                #Reverse complement R sequences, merge F and R, write merged fasta.
                R_rc_sequences_dict, F_R_rc_Sequences_dict=RC_merge_R(F_sequences_dict, R_sequences_dict, output_path+str(pair_name)+'_sequences_under_TCSs_full_'+str(threshold)+'_F_R_rc.fasta')
                R_rc_sequences_dict_nbp, F_R_rc_Sequences_dict_nbp=RC_merge_R(F_sequences_dict_nbp, R_sequences_dict_nbp, output_path+str(pair_name)+'_sequences_under_TCSs_'+str(sum(win_range_local_F))+'bp_LOGO_'+str(threshold)+'_F_R_rc.fasta')
                
                
                #Convert dict to list.
                F_sequences_list=dict_to_list(F_sequences_dict)
                R_sequences_list=dict_to_list(R_sequences_dict)
                F_R_rc_Sequences_list=dict_to_list(F_R_rc_Sequences_dict)
                
                #Compute PWM.
                F_PFMs=make_PFM(F_sequences_list)
                R_PFMs=make_PFM(R_sequences_list)
                All_PFMs=make_PFM(F_R_rc_Sequences_list)
                
                #Get consensus sequence, plot motif using weblogo.
                get_consensus_plot_motif(output_path+str(pair_name)+'_sequences_under_TCSs_'+str(sum(win_range_local_F))+'bp_LOGO_'+str(threshold)+'_F_R_rc.fasta', output_path+str(pair_name)+'_sequences_under_TCSs_'+str(sum(win_range_local_F))+'bp_LOGO_'+str(threshold)+'_F_R_rc')  
                
                #Write motif.
                write_motif(F_PFMs[PFM_type], output_path+str(pair_name)+'_GC_pfm_F.txt', win_width)
                write_motif(R_PFMs[PFM_type], output_path+str(pair_name)+'_GC_pfm_R.txt', win_width)
                write_motif(All_PFMs[PFM_type], output_path+str(pair_name)+'_GC_pfm_F_R_rc.txt', win_width)
                
                dict_of_PFMs[pair_name]={'F' : F_PFMs[PFM_type], 'R' : R_PFMs[PFM_type], 'FR' : All_PFMs[PFM_type]}
                
                #Plot nucleotides frequences around cleavage sites, add statistics.
                Plotting_stat(F_PFMs[PFM_type], F_PFMs['Num_seqs'], 'TopoI motif F_statistic for ' + pair_name, PFM_type, genome_sequence, output_path+PFM_type+'_TopoI_motif_trusted_TCSs_statistic_'+str(threshold)+'_F_'+str(pair_name)+'.png', win_width)
                Plotting_stat(R_PFMs[PFM_type], R_PFMs['Num_seqs'], 'TopoI motif R_statistic for ' + pair_name, PFM_type, genome_sequence, output_path+PFM_type+'_TopoI_motif_trusted_TCSs_statistic_'+str(threshold)+'_R_'+str(pair_name)+'.png', win_width)
                Plotting_stat(All_PFMs[PFM_type], All_PFMs['Num_seqs'], 'TopoI motif R_statistic for ' + pair_name, PFM_type, genome_sequence, output_path+PFM_type+'_TopoI_motif_trusted_TCSs_statistic_'+str(threshold)+'_F_R_rc_'+str(pair_name)+'.png', win_width)
                
        return

wrap_function(path_to_wig_files, Genome_path, Threshold, Genome_ID, Output_path)



#Path to the file with intervals to be masked (bed).
Path_to_masked_regions='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\Scripts_TopoI_Topo-seq\Additional_genome_features\\Regions_to_be_masked.broadPeak'

#Signal confidence interval.
Confidence=0.01

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
#Returns nearby NE value if current position falls into deleted region of the genome.
#######

def get_value(i, ends, deletions):
        if i<0: #coordinate is out of the left genome border (start)
                j=len(ends)+i
        elif i>=len(ends): #coordinate is out of the right genome border (end)
                j=i-len(ends)
        else: #coordinate is within the genome borders
                check_in_del=0
                for dl in deletions: #check if coordinate falls into deletion
                        if dl[1]>=i>=dl[0]:
                                j=dl[1]-dl[0]+i+1
                                check_in_del=1
                if check_in_del==0:
                        j=i
        return ends[j]

#######
#Returns smoothed N3/5E tracks.
#Smoothing using sliding window (default - 200000 nt).
#######

def Smoothing(ends, deletions):
        smoothed=[]
        #Calculating the value for the first genome position
        mean=0.0
        window=100000
        window_float=float(window)
        for i in range(-window, window):
                mean=mean + get_value(i, ends, deletions)
        mean=mean/(2*window_float)
        smoothed.append(mean)
        #Calculating values for the part of the genome remains
        for i in range(1, len(ends)):
                mean=mean + (get_value(i+window, ends, deletions) - get_value(i-window, ends, deletions))/(2*window_float)
                smoothed.append(mean)
        return smoothed

#######
#Returns Ara+IP+/smoothed(Ara+IP-) and Ara-IP+/smoothed(Ara-IP-) tracks ready for TCSs calling.
#######

def norm_smooth_devide(ex_file_path, cont_file_path, un_ex_file_path, un_cont_file_path, deletions):
        
        #WIG parsing, total NE counting (for further normalization on reads number)
        treated_experiment=wig_parsing(ex_file_path) #+Ara+IP
        treated_control=wig_parsing(cont_file_path) #+Ara-IP
        untreated_experiment=wig_parsing(un_ex_file_path) #-Ara+IP
        untreated_control=wig_parsing(un_cont_file_path) #-Ara-IP
        #Normalization on the reads number
        #Adds pseudocounts to avoid zero values
        Min_total_NE=min(treated_experiment[1], treated_control[1], 
                     untreated_experiment[1], untreated_control[1])
        print('Min_total_NE: ' + str(Min_total_NE))
        treated_experiment_norm=[1.0 * (x + 1) * Min_total_NE/treated_experiment[1] for x in treated_experiment[0]] #+Ara+IP norm
        treated_control_norm=[1.0 * (x + 1) * Min_total_NE/treated_control[1] for x in treated_control[0]] #+Ara-IP norm
        untreated_experiment_norm=[1.0 * (x + 1) * Min_total_NE/untreated_experiment[1] for x in untreated_experiment[0]] #-Ara+IP norm
        untreated_control_norm=[1.0 * (x + 1) * Min_total_NE/untreated_control[1] for x in untreated_control[0]] #-Ara-IP norm
        #Control samples smoothing: +Ara-IP and -Ara-IP
        treated_control_norm_sm=Smoothing(treated_control_norm, deletions) #+Ara-IP norm sm 
        untreated_control_norm_sm=Smoothing(untreated_control_norm, deletions) #-Ara-IP norm sm
        #Pairwise division: +Ara+IP/+Ara-IP and -Ara+IP/-Ara-IP
        ends_divide_Ara=[] #+Ara+IP/+Ara-IP
        ends_divide_no_Ara=[] #-Ara+IP/-Ara-IP
        for i in range (len(treated_experiment_norm)):
                if treated_experiment_norm[i]!=0 and treated_control_norm_sm[i]!=0:
                        ends_divide_Ara.append(treated_experiment_norm[i]/treated_control_norm_sm[i])
                else:
                        ends_divide_Ara.append(0)
                if untreated_experiment_norm[i]!=0 and untreated_control_norm_sm[i]!=0:
                        ends_divide_no_Ara.append(untreated_experiment_norm[i]/untreated_control_norm_sm[i])
                else:
                        ends_divide_no_Ara.append(0) 
        return ends_divide_Ara, ends_divide_no_Ara, treated_control_norm_sm, untreated_control_norm_sm

#######
#Plots the enrichment signal over the genome: +Ara+IP/smoothed(+Ara-IP) and -Ara+IP/smoothed(-Ara-IP)
#######

def plot_enrichment_signal(fname, track_ar, deletions, path_out):
        
        ends_divide_Ara=track_ar[0]
        ends_divide_no_Ara=track_ar[1]
        treated_control_norm_sm=track_ar[2]
        untreated_control_norm_sm=track_ar[3]
        
        #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
        #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
        mpl.rcParams['agg.path.chunksize']=10000
        #Scaling smoothed tracks to make them visible on the plot.
        max_element=max(ends_divide_Ara+ends_divide_no_Ara) #Max N3E value of experimental tracks
        max_element_tc_sm=max(treated_control_norm_sm)
        max_element_utc_sm=max(untreated_control_norm_sm)
        treated_control_norm_sm=[(max_element/2)*x/max_element_tc_sm for x in treated_control_norm_sm]
        untreated_control_norm_sm=[(max_element/2)*x/max_element_utc_sm for x in untreated_control_norm_sm]
        #Regions to be masked (e.g. deletions).  
        mask_array=[]
        for k in range(len(ends_divide_Ara)):
                check_in_del=0
                for dl in deletions:
                        if dl[1]>=k>=dl[0]:
                                mask_array.append(True)
                                check_in_del=1
                if check_in_del==0:
                        mask_array.append(False)
        Ara_exp=np.ma.masked_array(ends_divide_Ara, mask=mask_array)
        No_Ara_exp=np.ma.masked_array(ends_divide_no_Ara, mask=mask_array)
        Ara_cont_sm=np.ma.masked_array(treated_control_norm_sm, mask=mask_array)
        No_Ara_cont_sm=np.ma.masked_array(untreated_control_norm_sm, mask=mask_array)
        #Plotting the distribution of the signal over the genome.
        xcoord=np.arange(0,4647999)
        plt.figure(figsize=(16, 8), dpi=100)
        plt.suptitle(fname, fontsize=20)
        plot1=plt.subplot() 
        plot1.plot(xcoord, Ara_exp, '-', label='+Ara+IP/smoothed(+Ara-IP)', color='blue', linewidth=1)
        plot1.plot(xcoord, No_Ara_exp, '-', label='-Ara+IP/smoothed(-Ara-IP)', color='orange', linewidth=1)
        plot1.plot(xcoord, Ara_cont_sm, '-', label='smoothed(+Ara-IP)', color='#5bbdff', linewidth=3)
        plot1.plot(xcoord, No_Ara_cont_sm, '-', label='smoothed(-Ara-IP)', color='#ed781f', linewidth=3)    
        plot1.set_xlabel('Genome position, nt', size=17)
        plot1.set_ylabel('Signal enrichment', size=17)
        plot1.legend(loc='upper right')
        plt.show()
        plt.savefig(f'{path_out}{fname}_signal_enrichment.png', dpi=300, figsize=(16, 8))
        plt.close()
        return

#######
#Audic & Claverie statistics: borders of the confidential intervals (p-value=0.05, two-tailed test).
#From Audic & Claverie, 1997
#######

def AC_stat(x, confidence):
        x+=-1
        #Confidential intervals borders (from Audic & Claverie, 1997).
        if confidence==0.05:
                AU_test=[5,7,9,11,12,14,16,17,19,20,22,23,24,26,27,28,30,31,32,34,35]
                AU_test20=20*1.75
                AU_test25=25*1.64
                AU_test30=30*1.60
                AU_test40=40*1.50
                AU_test50=50*1.44
                AU_test75=75*1.36   
                AU_test100=100*1.30
        elif confidence==0.01:
                AU_test=[7,10,12,14,16,18,19,21,23,24,26,27,29,30,32,33,35,36,38,39,40]
                AU_test20=20*2
                AU_test25=25*1.88
                AU_test30=30*1.80
                AU_test40=40*1.68
                AU_test50=50*1.60
                AU_test75=75*1.48  
                AU_test100=100*1.40     
        #Estimation of a confidential interval higher border according to the value given - x.
        #Returns the interval border.
        if x<len(AU_test):
                int_border=AU_test[int(x)]
        elif 25>x>=20:
                int_border=AU_test20
        elif 30>x>=25:
                int_border=AU_test25
        elif 40>x>=30:
                int_border=AU_test30
        elif 50>x>=40:
                int_border=AU_test40
        elif 75>x>=50:
                int_border=AU_test50
        elif 100>x>=75:
                int_border=AU_test75
        else:
                int_border=AU_test100
        return int_border

#######
#Detect TCSs.
#######

def detect_TCSs_AC(wig_data_norm_Ara, wig_data_norm_no_Ara, confidence, win_range_local):
        
        TCSs_list=[]
        TCSs_metasignal=np.array([0.0]*sum(win_range_local))        
        for i in range(len(wig_data_norm_Ara)):
                if wig_data_norm_Ara[i]>AC_stat(wig_data_norm_no_Ara[i], confidence):
                        TCSs_list.append([i, wig_data_norm_Ara[i]])                
                
                        cleavage_region_wig=wig_data_norm_Ara[i-win_range_local[0]:i+win_range_local[1]]
                        TCSs_metasignal+=np.array(cleavage_region_wig)
                        
        print('Number of TCSs just found: ' + str(len(TCSs_list)))
        thr=25000
        if (len(TCSs_list)>thr):
                print('Number of TCSs is extremely high! The threshold is ' + str(thr) + '.\nJust warning...') 
                
        return TCSs_list, TCSs_metasignal


#######
#Wraps all the functions together.
#######

def wrap_function_AC(dict_of_wigs_path, del_path, fname, genome_input_path, confidence, Genome_ID, output_path):
        
        win_width=100
        win_range_F=[(win_width/2)-2, (win_width/2)+2]
        win_range_R=[(win_width/2)-3, (win_width/2)+1]
        win_width_local=20
        win_range_local_F=[10, 10]
        win_range_local_R=[9, 11]
        PFM_type='GC'
        plot_title='TopoI cleavage motif obtained with Topo-Seq'
        
        #Read regions to be ommitted (e.g., deletions).
        deletions=deletions_info(del_path)
                
        #Read and prepare tracks for signal plotting and TCSs calling.
        Prepared_tracks_F=norm_smooth_devide(dict_of_wigs_path['TopoI_Ara_IP']['F'], dict_of_wigs_path['TopoI_Ara_mock']['F'], dict_of_wigs_path['TopoI_no_Ara_IP']['F'], dict_of_wigs_path['TopoI_no_Ara_mock']['F'], deletions)
        Prepared_tracks_R=norm_smooth_devide(dict_of_wigs_path['TopoI_Ara_IP']['R'], dict_of_wigs_path['TopoI_Ara_mock']['R'], dict_of_wigs_path['TopoI_no_Ara_IP']['R'], dict_of_wigs_path['TopoI_no_Ara_mock']['R'], deletions)
        
        #Plot signal over the genome.
        plot_enrichment_signal(f'{fname}_F', Prepared_tracks_F, deletions, output_path)
        plot_enrichment_signal(f'{fname}_R', Prepared_tracks_R, deletions, output_path)
        
        #Detect TCSs - TopoI cleavage sites.
        TCSs_F, TCSs_metasignal_F=detect_TCSs_AC(Prepared_tracks_F[0], Prepared_tracks_F[1], confidence, win_range_local_F)
        TCSs_R, TCSs_metasignal_R=detect_TCSs_AC(Prepared_tracks_R[0], Prepared_tracks_R[1], confidence, win_range_local_R)
        
        #Cleavage signal in the vicinity of TCSs.
        TCSs_metasignal_F_R_rc=TCSs_metasignal_F+np.array(TCSs_metasignal_R.tolist()[::-1])
        TCSs_metasignal_F_R_rc_scaled=TCSs_metasignal_F_R_rc/(len(TCSs_F)+len(TCSs_R))
        
        #Write TCSs coordinates into BroadPeak file.
        write_TCSs_coords(TCSs_F, TCSs_R, Genome_ID, f'{output_path}{fname}_TCSs_called_AC_{confidence}.BroadPeak')  
        
        #Plot N3E signal around cleavage sites.
        Plotting_N3E_metasignal(TCSs_metasignal_F_R_rc_scaled, [len(TCSs_F), len(TCSs_R)], f'TopoI Topo-Seq N3E statistic for {fname}', f'{output_path}_TopoI_Topo_Seq_TCSs_cleavage_signal_AC_{confidence}_F_R_rc_{fname}.png', win_width_local)
        
        #Return sequences under the TCSs, plot motif.
        #Extract genomic sequence.        
        genome_sequence=genome_seq(genome_input_path)  
        
        #Get sequences.
        F_sequences_dict, F_sequences_dict_nbp=return_seqs(TCSs_F, 'F', win_range_F, win_range_local_F, genome_sequence, f'{output_path}{fname}_sequences_under_TCSs_full_AC_{confidence}_F.fasta', f'{output_path}{fname}_sequences_under_TCSs_{sum(win_range_local_F)}bp_LOGO_AC_{confidence}_F.fasta')
        R_sequences_dict, R_sequences_dict_nbp=return_seqs(TCSs_R, 'R', win_range_R, win_range_local_R, genome_sequence, f'{output_path}{fname}_sequences_under_TCSs_full_AC_{confidence}_R.fasta', f'{output_path}{fname}_sequences_under_TCSs_{sum(win_range_local_R)}bp_LOGO_AC_{confidence}_R.fasta')
                
        #Reverse complement R sequences, merge F and R, write merged fasta.
        R_rc_sequences_dict, F_R_rc_Sequences_dict=RC_merge_R(F_sequences_dict, R_sequences_dict, f'{output_path}{fname}_sequences_under_TCSs_full_AC_{confidence}_F_R_rc.fasta')
        R_rc_sequences_dict_nbp, F_R_rc_Sequences_dict_nbp=RC_merge_R(F_sequences_dict_nbp, R_sequences_dict_nbp, f'{output_path}{fname}_sequences_under_TCSs_{sum(win_range_local_F)}bp_LOGO_AC_{confidence}_F_R_rc.fasta')
                
        #Convert dict to list.
        F_sequences_list=dict_to_list(F_sequences_dict)
        R_sequences_list=dict_to_list(R_sequences_dict)
        F_R_rc_Sequences_list=dict_to_list(F_R_rc_Sequences_dict)
        
        #Compute PWM.
        F_PFMs=make_PFM(F_sequences_list)
        R_PFMs=make_PFM(R_sequences_list)
        All_PFMs=make_PFM(F_R_rc_Sequences_list)
        
        #Get consensus sequence, plot motif using weblogo.
        get_consensus_plot_motif(f'{output_path}{fname}_sequences_under_TCSs_{sum(win_range_local_F)}bp_LOGO_AC_{confidence}_F_R_rc.fasta', f'{output_path}{fname}_sequences_under_TCSs_{sum(win_range_local_F)}bp_LOGO_AC_{confidence}_F_R_rc')  
        
        #Write motif.
        write_motif(F_PFMs[PFM_type],   f'{output_path}{fname}_GC_pfm_F.txt', win_width)
        write_motif(R_PFMs[PFM_type],   f'{output_path}{fname}_GC_pfm_R.txt', win_width)
        write_motif(All_PFMs[PFM_type], f'{output_path}{fname}_GC_pfm_F_R_rc.txt', win_width)
        
        #Plot nucleotides frequences around cleavage sites, add statistics.
        Plotting_stat(F_PFMs[PFM_type],   F_PFMs['Num_seqs'],   f'TopoI motif F_statistic for {fname}',   PFM_type, genome_sequence, f'{output_path}{PFM_type}_TopoI_motif_TCSs_statistic_AC_{confidence}_F_{fname}.png', win_width)
        Plotting_stat(R_PFMs[PFM_type],   R_PFMs['Num_seqs'],   f'TopoI motif R_statistic for {fname}',   PFM_type, genome_sequence, f'{output_path}{PFM_type}_TopoI_motif_TCSs_statistic_AC_{confidence}_R_{fname}.png', win_width)
        Plotting_stat(All_PFMs[PFM_type], All_PFMs['Num_seqs'], f'TopoI motif F_R_statistic for {fname}', PFM_type, genome_sequence, f'{output_path}{PFM_type}_TopoI_motif_TCSs_statistic_AC_{confidence}_F_R_rc_{fname}.png', win_width)
                
        return

#wrap_function_AC(path_to_wig_files_1, Path_to_masked_regions, 'TopoI_Topo_Seq_1', Genome_path, Confidence, Genome_ID, Output_path)



#######
#Wraps all the functions together.
#######

def wrap_function_Threshold(dict_of_wigs_path, del_path, fname, genome_input_path, threshold, Genome_ID, output_path, pwd):
        
        win_width=100
        win_range_F=[(win_width/2)-2, (win_width/2)+2]
        win_range_R=[(win_width/2)-3, (win_width/2)+1]
        win_width_local=20
        win_range_local_F=[10, 10]
        win_range_local_R=[9, 11]
        PFM_type='GC'
        plot_title='TopoI cleavage motif obtained with Topo-Seq'
        
        #Read regions to be ommitted (e.g., deletions).
        deletions=deletions_info(del_path)
                
        #Read and prepare tracks for signal plotting and TCSs calling.
        Prepared_tracks_F=norm_smooth_devide(dict_of_wigs_path['TopoI_Ara_IP']['F'], dict_of_wigs_path['TopoI_Ara_mock']['F'], dict_of_wigs_path['TopoI_no_Ara_IP']['F'], dict_of_wigs_path['TopoI_no_Ara_mock']['F'], deletions)
        Prepared_tracks_R=norm_smooth_devide(dict_of_wigs_path['TopoI_Ara_IP']['R'], dict_of_wigs_path['TopoI_Ara_mock']['R'], dict_of_wigs_path['TopoI_no_Ara_IP']['R'], dict_of_wigs_path['TopoI_no_Ara_mock']['R'], deletions)
        
        #Plot signal over the genome.
        plot_enrichment_signal(f'{fname}_F', Prepared_tracks_F, deletions, output_path)
        plot_enrichment_signal(f'{fname}_R', Prepared_tracks_R, deletions, output_path)
        
        #Write resultant wig files.
        write_wig(Prepared_tracks_F[0], f'{pwd}\WIG_NE_strand_specific_masked_FE_smoothed\\{fname}_Ara_IP_div_by_smoothed_Ara_mock_F.wig', f'{fname}_Ara_IP_div_by_smoothed_Ara_mock_F')
        write_wig(Prepared_tracks_R[0], f'{pwd}\WIG_NE_strand_specific_masked_FE_smoothed\\{fname}_Ara_IP_div_by_smoothed_Ara_mock_R.wig', f'{fname}_Ara_IP_div_by_smoothed_Ara_mock_R')
        write_wig(Prepared_tracks_F[1], f'{pwd}\WIG_NE_strand_specific_masked_FE_smoothed\\{fname}_no_Ara_IP_div_by_smoothed_no_Ara_mock_F.wig', f'{fname}_no_Ara_IP_div_by_smoothed_no_Ara_mock_F')
        write_wig(Prepared_tracks_R[1], f'{pwd}\WIG_NE_strand_specific_masked_FE_smoothed\\{fname}_no_Ara_IP_div_by_smoothed_no_Ara_mock_R.wig', f'{fname}_no_Ara_IP_div_by_smoothed_no_Ara_mock_R')
        
        #Detect TCSs - TopoI cleavage sites.
        #+Ara tracks.
        TCSs_Ara_F, TCSs_Ara_metasignal_F=detect_TCSs(Prepared_tracks_F[0], threshold, win_range_local_F)
        TCSs_Ara_R, TCSs_Ara_metasignal_R=detect_TCSs(Prepared_tracks_R[0], threshold, win_range_local_R)
        #-Ara tracks.
        TCSs_no_Ara_F, TCSs_no_Ara_metasignal_F=detect_TCSs(Prepared_tracks_F[1], threshold, win_range_local_F)
        TCSs_no_Ara_R, TCSs_no_Ara_metasignal_R=detect_TCSs(Prepared_tracks_R[1], threshold, win_range_local_R)        
        
        
        #Cleavage signal in the vicinity of TCSs.
        TCSs_Ara_metasignal_F_R_rc=TCSs_Ara_metasignal_F+np.array(TCSs_Ara_metasignal_R.tolist()[::-1])
        TCSs_Ara_metasignal_F_R_rc_scaled=TCSs_Ara_metasignal_F_R_rc/(len(TCSs_Ara_F)+len(TCSs_Ara_R))
        TCSs_no_Ara_metasignal_F_R_rc=TCSs_no_Ara_metasignal_F+np.array(TCSs_no_Ara_metasignal_R.tolist()[::-1])
        TCSs_no_Ara_metasignal_F_R_rc_scaled=TCSs_no_Ara_metasignal_F_R_rc/(len(TCSs_no_Ara_F)+len(TCSs_no_Ara_R))        
        
        #Write TCSs coordinates into BroadPeak file.
        write_TCSs_coords(TCSs_Ara_F, TCSs_Ara_R, Genome_ID, f'{output_path}{fname}_Ara_TCSs_called_thr_{threshold}.BroadPeak') 
        write_TCSs_coords(TCSs_no_Ara_F, TCSs_no_Ara_R, Genome_ID, f'{output_path}{fname}_no_Ara_TCSs_called_thr_{threshold}.BroadPeak')  
        
        #Plot N3E signal around cleavage sites.
        Plotting_N3E_metasignal(TCSs_Ara_metasignal_F_R_rc_scaled,    [len(TCSs_Ara_F), len(TCSs_Ara_R)],       f'TopoI Topo-Seq N3E statistic for {fname} Ara',    f'{output_path}_TopoI_Topo_Seq_TCSs_cleavage_signal_Ara_thr_{threshold}_F_R_rc_{fname}.png', win_width_local)
        Plotting_N3E_metasignal(TCSs_no_Ara_metasignal_F_R_rc_scaled, [len(TCSs_no_Ara_F), len(TCSs_no_Ara_R)], f'TopoI Topo-Seq N3E statistic for {fname} no Ara', f'{output_path}_TopoI_Topo_Seq_TCSs_cleavage_signal_no_Ara_thr_{threshold}_F_R_rc_{fname}.png', win_width_local)        
        
        #Return sequences under the TCSs, plot motif.
        #Extract genomic sequence.        
        genome_sequence=genome_seq(genome_input_path)  
        
        #Get sequences.
        F_sequences_dict_Ara, F_sequences_dict_nbp_Ara=return_seqs(TCSs_Ara_F, 'F', win_range_F, win_range_local_F, genome_sequence, f'{output_path}{fname}_sequences_under_TCSs_Ara_full_thr_{threshold}_F.fasta', f'{output_path}{fname}_sequences_under_TCSs_Ara_{sum(win_range_local_F)}bp_LOGO_thr_{threshold}_F.fasta')
        R_sequences_dict_Ara, R_sequences_dict_nbp_Ara=return_seqs(TCSs_Ara_R, 'R', win_range_R, win_range_local_R, genome_sequence, f'{output_path}{fname}_sequences_under_TCSs_Ara_full_thr_{threshold}_R.fasta', f'{output_path}{fname}_sequences_under_TCSs_Ara_{sum(win_range_local_R)}bp_LOGO_thr_{threshold}_R.fasta')
        F_sequences_dict_no_Ara, F_sequences_dict_nbp_no_Ara=return_seqs(TCSs_no_Ara_F, 'F', win_range_F, win_range_local_F, genome_sequence, f'{output_path}{fname}_sequences_under_TCSs_no_Ara_full_thr_{threshold}_F.fasta', f'{output_path}{fname}_sequences_under_TCSs_no_Ara_{sum(win_range_local_F)}bp_LOGO_thr_{threshold}_F.fasta')
        R_sequences_dict_no_Ara, R_sequences_dict_nbp_no_Ara=return_seqs(TCSs_no_Ara_R, 'R', win_range_R, win_range_local_R, genome_sequence, f'{output_path}{fname}_sequences_under_TCSs_no_Ara_full_thr_{threshold}_R.fasta', f'{output_path}{fname}_sequences_under_TCSs_no_Ara_{sum(win_range_local_R)}bp_LOGO_thr_{threshold}_R.fasta')
         
        #Reverse complement R sequences, merge F and R, write merged fasta.
        R_rc_sequences_dict_Ara, F_R_rc_Sequences_dict_Ara=RC_merge_R(F_sequences_dict_Ara, R_sequences_dict_Ara, f'{output_path}{fname}_sequences_under_TCSs_Ara_full_thr_{threshold}_F_R_rc.fasta')
        R_rc_sequences_dict_nbp_Ara, F_R_rc_Sequences_dict_nbp_Ara=RC_merge_R(F_sequences_dict_nbp_Ara, R_sequences_dict_nbp_Ara, f'{output_path}{fname}_sequences_under_TCSs_Ara_{sum(win_range_local_F)}bp_LOGO_thr_{threshold}_F_R_rc.fasta')
        R_rc_sequences_dict_no_Ara, F_R_rc_Sequences_dict_no_Ara=RC_merge_R(F_sequences_dict_no_Ara, R_sequences_dict_no_Ara, f'{output_path}{fname}_sequences_under_TCSs_no_Ara_full_thr_{threshold}_F_R_rc.fasta')
        R_rc_sequences_dict_nbp_no_Ara, F_R_rc_Sequences_dict_nbp_no_Ara=RC_merge_R(F_sequences_dict_nbp_no_Ara, R_sequences_dict_nbp_no_Ara, f'{output_path}{fname}_sequences_under_TCSs_no_Ara_{sum(win_range_local_F)}bp_LOGO_thr_{threshold}_F_R_rc.fasta')
                
        #Convert dict to list.
        F_sequences_list_Ara=dict_to_list(F_sequences_dict_Ara)
        R_sequences_list_Ara=dict_to_list(R_sequences_dict_Ara)
        F_R_rc_Sequences_list_Ara=dict_to_list(F_R_rc_Sequences_dict_Ara)
        F_sequences_list_no_Ara=dict_to_list(F_sequences_dict_no_Ara)
        R_sequences_list_no_Ara=dict_to_list(R_sequences_dict_no_Ara)
        F_R_rc_Sequences_list_no_Ara=dict_to_list(F_R_rc_Sequences_dict_no_Ara)        
        
        #Compute PWM.
        F_PFMs_Ara=make_PFM(F_sequences_list_Ara)
        R_PFMs_Ara=make_PFM(R_sequences_list_Ara)
        All_PFMs_Ara=make_PFM(F_R_rc_Sequences_list_Ara)
        F_PFMs_no_Ara=make_PFM(F_sequences_list_no_Ara)
        R_PFMs_no_Ara=make_PFM(R_sequences_list_no_Ara)
        All_PFMs_no_Ara=make_PFM(F_R_rc_Sequences_list_no_Ara)        
        
        #Get consensus sequence, plot motif using weblogo.
        get_consensus_plot_motif(f'{output_path}{fname}_sequences_under_TCSs_Ara_{sum(win_range_local_F)}bp_LOGO_thr_{threshold}_F_R_rc.fasta', f'{output_path}{fname}_sequences_under_TCSs_Ara_{sum(win_range_local_F)}bp_LOGO_thr_{threshold}_F_R_rc')  
        get_consensus_plot_motif(f'{output_path}{fname}_sequences_under_TCSs_no_Ara_{sum(win_range_local_F)}bp_LOGO_thr_{threshold}_F_R_rc.fasta', f'{output_path}{fname}_sequences_under_TCSs_no_Ara_{sum(win_range_local_F)}bp_LOGO_thr_{threshold}_F_R_rc')  
        
        #Write motif.
        write_motif(F_PFMs_Ara[PFM_type],   f'{output_path}{fname}_Ara_thr_{threshold}_GC_pfm_F.txt', win_width)
        write_motif(R_PFMs_Ara[PFM_type],   f'{output_path}{fname}_Ara_thr_{threshold}_GC_pfm_R.txt', win_width)
        write_motif(All_PFMs_Ara[PFM_type], f'{output_path}{fname}_Ara_thr_{threshold}_GC_pfm_F_R_rc.txt', win_width)
        write_motif(F_PFMs_no_Ara[PFM_type],   f'{output_path}{fname}_no_Ara_thr_{threshold}_GC_pfm_F.txt', win_width)
        write_motif(R_PFMs_no_Ara[PFM_type],   f'{output_path}{fname}_no_Ara_thr_{threshold}_GC_pfm_R.txt', win_width)
        write_motif(All_PFMs_no_Ara[PFM_type], f'{output_path}{fname}_no_Ara_thr_{threshold}_GC_pfm_F_R_rc.txt', win_width)        
        
        #Plot nucleotides frequences around cleavage sites, add statistics.
        Plotting_stat(F_PFMs_Ara[PFM_type],   F_PFMs_Ara['Num_seqs'],   f'TopoI motif F_statistic for {fname} Ara',   PFM_type, genome_sequence, f'{output_path}{PFM_type}_TopoI_motif_TCSs_Ara_statistic_thr_{threshold}_F_{fname}.png', win_width)
        Plotting_stat(R_PFMs_Ara[PFM_type],   R_PFMs_Ara['Num_seqs'],   f'TopoI motif R_statistic for {fname} Ara',   PFM_type, genome_sequence, f'{output_path}{PFM_type}_TopoI_motif_TCSs_Ara_statistic_thr_{threshold}_R_{fname}.png', win_width)
        Plotting_stat(All_PFMs_Ara[PFM_type], All_PFMs_Ara['Num_seqs'], f'TopoI motif F_R_statistic for {fname} Ara', PFM_type, genome_sequence, f'{output_path}{PFM_type}_TopoI_motif_TCSs_Ara_statistic_thr_{threshold}_F_R_rc_{fname}.png', win_width)
        Plotting_stat(F_PFMs_no_Ara[PFM_type],   F_PFMs_no_Ara['Num_seqs'],   f'TopoI motif F_statistic for {fname} no Ara',   PFM_type, genome_sequence, f'{output_path}{PFM_type}_TopoI_motif_TCSs_no_Ara_statistic_thr_{threshold}_F_{fname}.png', win_width)
        Plotting_stat(R_PFMs_no_Ara[PFM_type],   R_PFMs_no_Ara['Num_seqs'],   f'TopoI motif R_statistic for {fname} no Ara',   PFM_type, genome_sequence, f'{output_path}{PFM_type}_TopoI_motif_TCSs_no_Ara_statistic_thr_{threshold}_R_{fname}.png', win_width)
        Plotting_stat(All_PFMs_no_Ara[PFM_type], All_PFMs_no_Ara['Num_seqs'], f'TopoI motif F_R_statistic for {fname} no Ara', PFM_type, genome_sequence, f'{output_path}{PFM_type}_TopoI_motif_TCSs_no_Ara_statistic_thr_{threshold}_F_R_rc_{fname}.png', win_width)
             
        return

#wrap_function_Threshold(path_to_wig_files_3, Path_to_masked_regions, 'TopoI_Topo_Seq_1', Genome_path, Threshold, Genome_ID, Output_path, PWD)

print('Script ended its work succesfully!')