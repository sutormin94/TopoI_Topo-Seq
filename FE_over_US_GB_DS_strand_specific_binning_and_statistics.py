###############################################
##Dmitry Sutormin, 2021##
##TopoA Topo-Seq analysis##

#Script computes Fold Enrichment (FE) over upstream (US)
#and downstream (DS) regions of transcription units (TUs)
#and over TUs bodies.
#Script is for strand-specific data: signal for F and R strands is provided separately.
#Script is dedicated for noisy data. To handle the data it performes data binning for smoothing.
#Also it keeps signal data for all TUs to contstruct confidential interval.
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
from scipy import stats
from scipy.interpolate import CubicSpline


#Path to the directory with input files.
PWD='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\\'
#Path to TUs groups file.
TUs_groups_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_RNA-Seq\Expression_data\DY330_transcripts\Representative_transcripts\\"
#Path to the input annotation, type of annotation and name of TUs set.
##1##
Path_to_annotation_1=TUs_groups_path + 'DY330_RNA-Seq_transcripts_representative_EP_del_cor.txt'
Type_of_annot_1='broadPeak'             
Genes_set_name_1='All_TUs_1672'    
##2##                                   
Path_to_annotation_2=TUs_groups_path + 'DY330_RNA-Seq_transcripts_representative_no_tRNA_rRNA_EP_del_cor.txt'
Type_of_annot_2='broadPeak'             
Genes_set_name_2='All_TUs_no_tRNA_rRNA_1623'    
##3##                                   
Path_to_annotation_3=TUs_groups_path + 'DY330_RNA-Seq_transcripts_representative_no_tRNA_rRNA_ompX_EP_del_cor_HETU_200.txt'
Type_of_annot_3='broadPeak'             
Genes_set_name_3='HETU_no_ompX_200'         
##4##                                   
Path_to_annotation_4=TUs_groups_path + 'DY330_RNA-Seq_transcripts_representative_no_tRNA_rRNA_appY_ybiI_EP_del_cor_LETU_200.txt'
Type_of_annot_4='broadPeak'             
Genes_set_name_4='LETU_no_appY_ybiI_200'             
##5##                                   
Path_to_annotation_5=TUs_groups_path + 'DY330_RNA-Seq_transcripts_EP_del_cor_rRNA_7.txt'
Type_of_annot_5='broadPeak'             
Genes_set_name_5='rRNA_7'               
##6##                                   
Path_to_annotation_6=TUs_groups_path + 'DY330_RNA-Seq_transcripts_EP_del_cor_tRNA_49.txt'
Type_of_annot_6='broadPeak'             
Genes_set_name_6='tRNA_49'              
##7##                                   
Path_to_annotation_7=TUs_groups_path + 'DY330_RNA-Seq_transcripts_representative_NO_DPS_EP_del_cor.txt'
Type_of_annot_7='broadPeak'             
Genes_set_name_7='All_TUs_no_dps_1660'    
##8##                                   
Path_to_annotation_8=TUs_groups_path + 'DY330_RNA-Seq_transcripts_representative_NO_DPS_NO_RFA_no_tRNA_rRNA_EP_del_cor_HETU_200.txt'
Type_of_annot_8='broadPeak'             
Genes_set_name_8='HETU_no_dps_rfa_200'         
##9##                                   
Path_to_annotation_9=TUs_groups_path + 'DY330_RNA-Seq_transcripts_representative_NO_DPS_no_tRNA_rRNA_EP_del_cor_LETU_200.txt'
Type_of_annot_9='broadPeak'             
Genes_set_name_9='LETU_no_dps_200'  
##10##                                   
Path_to_annotation_10=TUs_groups_path + 'DY330_RNA-Seq_transcripts_representative_NO_DPS_EP_Long_Active_TUs_del_cor.txt'
Type_of_annot_10='broadPeak'             
Genes_set_name_10='LATU_no_dps'  

#Path to the file with regions to be omitted (e.g. deletions).
Deletions_inpath='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\Scripts_TopoI_Topo-seq\TopoI_Topo-Seq\Additional_genome_features\\Deletions_w3110_G_Mu_SGS.broadPeak'
#Width of US, DS regions.
Win_width=15000
#Length of GB.
Length=5000
#Bin width.
Bin_width=200

#Dictionary of pathes to input data.
Dict_of_wigs_path={'TopoI_Ara_FE' :     {'F' : PWD + 'FE_cov_depth_masked_av\TopoI_Ara_f_FE_av_123.wig', 'R' : PWD + 'FE_cov_depth_masked_av\TopoI_Ara_r_FE_av_123.wig'},
                   'TopoI_FE' :         {'F' : PWD + 'FE_cov_depth_masked_av\TopoI_f_FE_av_123.wig',     'R' : PWD + 'FE_cov_depth_masked_av\TopoI_r_FE_av_123.wig'}
                   }
Dict_of_wigs_path_1={'TopoI_Ara_FE_subtr' :     {'F' : PWD + 'FE_cov_depth_masked_av_subtr_no_Ara\TopoI_Ara_f_FE_av_123_subtr_mock.wig', 'R' : PWD + 'FE_cov_depth_masked_av_subtr_no_Ara\TopoI_Ara_r_FE_av_123_subtr_mock.wig'},
                     }
Dict_of_wigs_ss_path_2={'TopoI_Ara_N3E_subtr_mock_subtr_no_Ara' :     {'F' : PWD + 'WIG_NE_strand_specific_masked_scaled_av_masked_accB_subtract_mock_subtract_no_Ara\TopoI_Ara_N3E_F_masked_scaled_av_123_subtr_mock_subtr_no_Ara.wig', 'R' : PWD + 'WIG_NE_strand_specific_masked_scaled_av_masked_accB_subtract_mock_subtract_no_Ara\TopoI_Ara_N3E_R_masked_scaled_av_123_subtr_mock_subtr_no_Ara.wig'},
                     }
Dict_of_wigs_path_2={'TopoI_no_CTD_no_Rif_FE' :  "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\Sutormin_TopA_ChIP_CTD_minus_Rif_minus_FE_av_346.wig"
                     }
Dict_of_wigs_path_3={'TopoI_Ara_FE_sm' :              {'F' : PWD + 'WIG_NE_strand_specific_masked_FE_smoothed_av\TopoI_Ara_F_N3E_FE_av_123.wig',    'R' : PWD + 'WIG_NE_strand_specific_masked_FE_smoothed_av\TopoI_Ara_R_N3E_FE_av_123.wig'},
                     'TopoI_no_Ara_FE_sm' :           {'F' : PWD + 'WIG_NE_strand_specific_masked_FE_smoothed_av\TopoI_no_Ara_F_N3E_FE_av_123.wig', 'R' : PWD + 'WIG_NE_strand_specific_masked_FE_smoothed_av\TopoI_no_Ara_R_N3E_FE_av_123.wig'},
                     'TopoI_Ara_FE_sm_FE_no_Ara' :    {'F' : PWD + 'WIG_NE_strand_specific_masked_FE_smoothed_av_FE_no_Ara\TopoI_Ara_F_N3E_FE_av_123_div_by_TopoI_no_Ara_F_N3E_FE_av_123.wig', 'R' : PWD + 'WIG_NE_strand_specific_masked_FE_smoothed_av_FE_no_Ara\TopoI_Ara_R_N3E_FE_av_123_div_by_TopoI_no_Ara_R_N3E_FE_av_123.wig'},
                     'TopoI_Ara_FE_sm_FE_no_Ara_sc' : {'F' : PWD + 'WIG_NE_strand_specific_masked_FE_smoothed_av_FE_no_Ara\TopoI_Ara_F_N3E_FE_av_123_div_by_TopoI_no_Ara_F_N3E_FE_av_123_scaled.wig', 'R' : PWD + 'WIG_NE_strand_specific_masked_FE_smoothed_av_FE_no_Ara\TopoI_Ara_R_N3E_FE_av_123_div_by_TopoI_no_Ara_R_N3E_FE_av_123_scaled.wig'}
                     }

#######
#Checks if directory exists and if not creates.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return

#Path to the output directory.
Out_path='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\Metagene_analysis\Transcripts_strand_specific_binned_stat_for_Fig_5\\'

#Output path.
def create_out_dirs(out_path, genes_set_name):
    Dir_check_create(out_path)
    Dir_check_create(out_path+'\Figures\Plots\\'+genes_set_name)
    Dir_check_create(out_path+'\Figures\Histograms\\'+genes_set_name)
    Dir_check_create(out_path+'\Signal_of_TUs_tab\\'+genes_set_name)
    Dir_check_create(out_path+'\Signal_of_TUs_wig\\'+genes_set_name)    
    return

#create_out_dirs(Out_path, Genes_set_name_1)
#create_out_dirs(Out_path, Genes_set_name_2)
#create_out_dirs(Out_path, Genes_set_name_3)
#create_out_dirs(Out_path, Genes_set_name_4)
#create_out_dirs(Out_path, Genes_set_name_2)
create_out_dirs(Out_path, Genes_set_name_5)
#create_out_dirs(Out_path, Genes_set_name_6)
create_out_dirs(Out_path, Genes_set_name_7)
create_out_dirs(Out_path, Genes_set_name_8)
create_out_dirs(Out_path, Genes_set_name_9)
#create_out_dirs(Out_path, Genes_set_name_10)


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
#Parsing gff file and preparing gene annotation.
#######

def parse_gff_annotation(gff_inpath, deletions_inpath):
    #Parsing deletions
    deletions=deletions_info(deletions_inpath)
    
    filein=open(gff_inpath, 'r')
    genes_annotation={'Gene': {},
                      'rRNA': {},
                      'tRNA': {},
                      'ncRNA': {}
                      }
    data_source={}
    for line in filein:
        line=line.rstrip()
        if line[0]!='#':
            line=line.split('\t')
            #What occurs in the annotation:
            if line[1] not in data_source:
                data_source[line[1]]={line[2]: 1}
            else:
                if line[2] not in data_source[line[1]]:
                    data_source[line[1]][line[2]]=1
                else:
                    data_source[line[1]][line[2]]+=1
            #Classify genes:
            #Protein coding genes.
            if line[1]=='ena' and line[2]=='gene':
                gene_start=int(line[3])
                gene_end=int(line[4])
                gene_strand=line[6]
                annot=line[8].split(';')
                if annot[1][0]=="N":
                    if annot[2].split('=')[1]=='protein_coding':
                        gene_name=annot[1].split('=')[1]
                        genes_annotation['Gene'][gene_name]=[gene_start, gene_end, gene_strand]
            #rRNA genes.
            elif line[1]=='ena' and line[2]=='rRNA_gene':
                gene_start=int(line[3])
                gene_end=int(line[4])
                gene_strand=line[6]
                annot=line[8].split(';')
                if annot[1][0]=="N":
                    if annot[2].split('=')[1]=='rRNA':
                        gene_name=annot[1].split('=')[1]
                        genes_annotation['rRNA'][gene_name]=[gene_start, gene_end, gene_strand] 
            #tRNA genes.
            elif line[1]=='ena' and line[2]=='tRNA_gene':
                gene_start=int(line[3])
                gene_end=int(line[4])
                gene_strand=line[6]
                annot=line[8].split(';')
                if annot[1][0]=="N":
                    if annot[2].split('=')[1]=='tRNA':
                        gene_name=annot[1].split('=')[1]
                        genes_annotation['tRNA'][gene_name]=[gene_start, gene_end, gene_strand]
            #Other non-coding RNAs.
            elif line[1]=='Rfam' and line[2]=='ncRNA_gene':
                gene_start=int(line[3])
                gene_end=int(line[4])
                gene_strand=line[6]
                annot=line[8].split(';')
                if annot[1][0]=="N":
                    if annot[2].split('=')[1]=='ncRNA':
                        gene_name=annot[1].split('=')[1]
                        genes_annotation['ncRNA'][gene_name]=[gene_start, gene_end, gene_strand]
    filein.close()            
    return genes_annotation, data_source


#######
#Reads annotation of particular set of genes .tab BroadPeak-like (determined on a basis of expression level).
#######

def parse_expression_annotation(annot_inpath):
    genes_annotation={}
    filein=open(annot_inpath, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        if line[0] not in ['GeneID', 'OperonID', 'TU_ID']:
            TU_name=line[1].lstrip('"').rstrip(';"')
            TU_start=int(line[2])
            TU_end=int(line[3])
            TU_strand=line[4]
            TU_expression=float(line[5].replace(',','.'))
            genes_annotation[TU_name]=[TU_start, TU_end, TU_strand, TU_expression]
    filein.close()            
    return genes_annotation


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
#Returns binned tracks.
#######

def Binning(data_ar, bin_width):
    Binned=[]
    #Calculating number of bins.
    N_bins=int(len(data_ar)/bin_width)
    for i in range(N_bins):
        bin_value=np.mean(data_ar[i*bin_width:(i+1)*bin_width])
        Binned.append(bin_value)
    if N_bins!=0:
        if len(data_ar)>((i+1)*bin_width):
            Binned.append(np.mean(data_ar[(i+1)*bin_width:]))
    elif N_bins==0:
        Binned.append(np.mean(data_ar))
    return Binned

#######
#Scale regions (gene bodies) to equal length: make long shorter and short longer.
#######

def scale_gene_body(ar, length):
    scaled=[]
    if len(ar)>length: #array should be shrinked
        #Determines positions to be taken (other positions will be discarded).
        positions_to_take=[]
        while len(positions_to_take)!=length:
            position=random.randint(0,len(ar)-1)
            if position not in positions_to_take:
                positions_to_take.append(position)
            else:
                continue
        positions_to_take.sort()
        for pos in positions_to_take:
            scaled.append(ar[pos])
    elif len(ar)<length:
        #Determine positions to be duplicated (other positions will be discarded).
        scaled=ar
        for i in range(length-len(ar)):
            position=random.randint(0,len(scaled))
            if position==0:
                scaled=scaled[:position+1]+scaled[position:position+1]+scaled[position+1:]
            else:
                scaled=scaled[:position]+scaled[position-1:position]+scaled[position:]        
    elif len(ar)==length:
        scaled=ar

    return scaled

#######
#Write .tab file with FE info for genes US, GB, and DS.
#######

def write_genes_FE(dict1, dict2, dict3, FE_track_name, path_out):
    fileout=open(path_out, 'w')
    fileout.write(f'Gene_name\tStart\tEnd\tStrand\t{FE_track_name}_FE_US\t{FE_track_name}_FE_GB\t{FE_track_name}_FE_DS\n')
    for k, v in dict1.items():
        fileout.write(f'{k}\t{v[1]}\t{v[2]}\t{v[3]}\t{v[0]}\t{dict2[k][0]}\t{dict3[k][0]}\n')
    fileout.close()
    return

#######
#Convert dictionary to array, discard keys.
#######

def dict_to_ar(dictionary):
    ar=[]
    for k,v in dictionary.items():
        ar.append(v[0]) 
    return ar


#########
##Makes histogram for FE over TUs: US, GB, DS.
#########

def plot_FE_dist_UDB(ar0, name0, ar1, name1, ar2, name2, pathout):
    #Plot distribution of FE values.
    
    mean_FE0=round(np.mean(ar0),2)
    print(f'Mean FE in {name0}={mean_FE0}')
    fig=plt.figure(figsize=(15, 3), dpi=100)
    bins0=np.arange(min(ar0+ar1+ar2), max(ar0+ar1+ar2), 0.25)
    plot0=plt.subplot2grid((1,3),(0,0), rowspan=1, colspan=1)
    plot0.hist(ar0, bins0, color='#ff878b', edgecolor='black', alpha=0.8, label=f'{name0}')
    plot0.annotate(f'Mean FE={mean_FE0}', xy=(0.60, 0.8), xycoords='axes fraction', size=15)
    plot0.set_yscale('log')
    plot0.set_xlabel('Fold enrichment', size=17)
    plot0.set_ylabel('Number of TUs', size=17)
    plot0.set_title(name0, size=18)  
    #plot0.legend(fontsize=22)
       
    mean_FE1=round(np.mean(ar1),2)
    print(f'Mean FE in {name1}={mean_FE1}')
    bins1=np.arange(min(ar0+ar1+ar2), max(ar0+ar1+ar2), 0.25)
    plot1=plt.subplot2grid((1,3),(0,1), rowspan=1, colspan=1)     
    plot1.hist(ar1, bins1, color='#ffce91', edgecolor='black', alpha=0.5, label=f'{name1}')
    plot1.annotate(f'Mean FE={mean_FE1}', xy=(0.60, 0.8), xycoords='axes fraction', size=15)
    plot1.set_yscale('log')
    plot1.set_xlabel('Fold enrichment', size=17)
    plot1.set_ylabel('Number of TUs', size=17)
    plot1.set_title(name1, size=18) 
    #plot1.legend(fontsize=22)
    
    mean_FE2=round(np.mean(ar2),2)
    print(f'Mean FE in {name2}={mean_FE2}')
    bins2=np.arange(min(ar0+ar1+ar2), max(ar0+ar1+ar2), 0.25)
    plot2=plt.subplot2grid((1,3),(0,2), rowspan=1, colspan=1) 
    plot2.hist(ar2, bins2, color='#7FCE79', edgecolor='black', alpha=0.5, label=f'{name2}')
    plot2.annotate(f'Mean FE={mean_FE2}', xy=(0.60, 0.8), xycoords='axes fraction', size=15)
    plot2.set_yscale('log')
    plot2.set_xlabel('Fold enrichment', size=17)
    plot2.set_ylabel('Number of TUs', size=17)
    plot2.set_title(name2, size=18)  
    #plot2.legend(fontsize=22)    
    
    plt.tight_layout()
    plt.show()
    plt.savefig(pathout, dpi=300, figsize=(15, 3))
    plt.close() 
    return


#######
#Computes standard error of mean.
#######

def compute_standard_error(ar):
    std_err=np.std(ar)/np.sqrt(len(ar))
    return std_err


#######
#Returns FE or Ded FE over the set of genes (US, GB, DS) - for each gene separately.
#For strand-specific data.
#Compares signal at US, DS, TU start, and TU end. Makes a barplot.
#######

def genes_and_FE_ss(gene_annotation, genes_set_name, FE_track_pair, FE_track_pair_name, out_path, deletions_inpath, win_width, length, bin_width):
    
    #Strand-specific tracks.
    FE_track_F=FE_track_pair['F']
    FE_track_R=FE_track_pair['R']
    
    #Parsing deletions
    deletions=deletions_info(deletions_inpath) 

    #Calculate FE over genes.
    gene_US_F=[]
    gene_DS_F=[]
    gene_B_F=[]
    gene_US_F_mean_dict={}
    gene_DS_F_mean_dict={}
    gene_B_F_mean_dict={}
    gene_US_R=[]
    gene_DS_R=[]
    gene_B_R=[]
    gene_US_R_mean_dict={}
    gene_DS_R_mean_dict={}
    gene_B_R_mean_dict={}   
    
    #Number of genes.
    Num_genes=len(gene_annotation)
    
    for gene_name, gene_info in gene_annotation.items():
        delited=0
        
        for deletion in deletions:
            if deletion[1]>=gene_info[0]>=deletion[0] or deletion[1]>=gene_info[1]>=deletion[0]:
                delited=1
                
        if delited==0:
            start=gene_info[0]
            end=gene_info[1]
            glen=len(FE_track_F)
            
            if gene_info[2]=='+':
                if start<win_width:
                    gene_US_F.append(FE_track_F[glen-(win_width-start):] + FE_track_F[:start])
                    gene_US_F_mean_dict[gene_name]=[np.mean(FE_track_F[glen-(win_width-start):] + FE_track_F[:start]), start, end, gene_info[2]]
                    gene_US_R.append(FE_track_R[glen-(win_width-start):] + FE_track_R[:start])
                    gene_US_R_mean_dict[gene_name]=[np.mean(FE_track_R[glen-(win_width-start):] + FE_track_R[:start]), start, end, gene_info[2]]                   
                else:
                    gene_US_F.append(FE_track_F[start-win_width:start])
                    gene_US_F_mean_dict[gene_name]=[np.mean(FE_track_F[start-win_width:start]), start, end, gene_info[2]]
                    gene_US_R.append(FE_track_R[start-win_width:start])
                    gene_US_R_mean_dict[gene_name]=[np.mean(FE_track_R[start-win_width:start]), start, end, gene_info[2]]                    
                if end+win_width>glen:
                    gene_DS_F.append(FE_track_F[end:] + FE_track_F[:end+win_width-glen])
                    gene_DS_F_mean_dict[gene_name]=[np.mean(FE_track_F[end:] + FE_track_F[:end+win_width-glen]), start, end, gene_info[2]]
                    gene_DS_R.append(FE_track_R[end:] + FE_track_R[:end+win_width-glen])
                    gene_DS_R_mean_dict[gene_name]=[np.mean(FE_track_R[end:] + FE_track_R[:end+win_width-glen]), start, end, gene_info[2]]                    
                else:
                    gene_DS_F.append(FE_track_F[end:end+win_width])
                    gene_DS_F_mean_dict[gene_name]=[np.mean(FE_track_F[end:end+win_width]), start, end, gene_info[2]]
                    gene_DS_R.append(FE_track_R[end:end+win_width])
                    gene_DS_R_mean_dict[gene_name]=[np.mean(FE_track_R[end:end+win_width]), start, end, gene_info[2]] 
                gene_B_F.append(FE_track_F[start:end])
                gene_B_F_mean_dict[gene_name]=[np.mean(FE_track_F[start:end]), start, end, gene_info[2]]
                gene_B_R.append(FE_track_R[start:end])
                gene_B_R_mean_dict[gene_name]=[np.mean(FE_track_R[start:end]), start, end, gene_info[2]]                
                
            elif gene_info[2]=='-':
                if start<win_width:
                    gene_DS_F.append(FE_track_R[:start][::-1] + FE_track_R[glen-(win_width-start):][::-1])                   #Take care of strand-specificity.
                    gene_DS_F_mean_dict[gene_name]=[np.mean(FE_track_R[:start][::-1] + FE_track_R[glen-(win_width-start):][::-1]), start, end, gene_info[2]] #Take care of strand-specificity.
                    gene_DS_R.append(FE_track_F[:start][::-1] + FE_track_F[glen-(win_width-start):][::-1])                   #Take care of strand-specificity.
                    gene_DS_R_mean_dict[gene_name]=[np.mean(FE_track_F[:start][::-1] + FE_track_F[glen-(win_width-start):][::-1]), start, end, gene_info[2]] #Take care of strand-specificity.                    
                else:
                    gene_DS_F.append(FE_track_R[start-win_width:start][::-1])                                                #Take care of strand-specificity.
                    gene_DS_F_mean_dict[gene_name]=[np.mean(FE_track_R[start-win_width:start][::-1]), start, end, gene_info[2]] #Take care of strand-specificity.
                    gene_DS_R.append(FE_track_F[start-win_width:start][::-1])                                                #Take care of strand-specificity.
                    gene_DS_R_mean_dict[gene_name]=[np.mean(FE_track_F[start-win_width:start][::-1]), start, end, gene_info[2]] #Take care of strand-specificity.                    
                if end+win_width>glen:
                    gene_US_F.append(FE_track_R[:end+win_width-glen][::-1] + FE_track_R[end:][::-1])                         #Take care of strand-specificity.
                    gene_US_F_mean_dict[gene_name]=[np.mean(FE_track_R[:end+win_width-glen][::-1] + FE_track_R[end:][::-1]), start, end, gene_info[2]] #Take care of strand-specificity.
                    gene_US_R.append(FE_track_F[:end+win_width-glen][::-1] + FE_track_F[end:][::-1])                         #Take care of strand-specificity.
                    gene_US_R_mean_dict[gene_name]=[np.mean(FE_track_F[:end+win_width-glen][::-1] + FE_track_F[end:][::-1]), start, end, gene_info[2]] #Take care of strand-specificity.                    
                else:
                    gene_US_F.append(FE_track_R[end:end+win_width][::-1])                                                    #Take care of strand-specificity.
                    gene_US_F_mean_dict[gene_name]=[np.mean(FE_track_R[end:end+win_width][::-1]), start, end, gene_info[2]]     #Take care of strand-specificity.
                    gene_US_R.append(FE_track_F[end:end+win_width][::-1])                                                    #Take care of strand-specificity.
                    gene_US_R_mean_dict[gene_name]=[np.mean(FE_track_F[end:end+win_width][::-1]), start, end, gene_info[2]]     #Take care of strand-specificity.                                                                                 #Take care of strand-specificity.
                gene_B_F.append(FE_track_R[start:end][::-1])
                gene_B_F_mean_dict[gene_name]=[np.mean(FE_track_R[start:end][::-1]), start, end, gene_info[2]]                                                              #Take care of strand-specificity.
                gene_B_R.append(FE_track_F[start:end][::-1])
                gene_B_R_mean_dict[gene_name]=[np.mean(FE_track_F[start:end][::-1]), start, end, gene_info[2]]
                
    #Data binning.
    print(len(gene_US_F), len(gene_B_F), len(gene_DS_F), len(gene_US_R), len(gene_B_R), len(gene_DS_R))
    
    gene_US_F_binned=[]
    gene_DS_F_binned=[]
    gene_B_F_binned=[]
    gene_US_R_binned=[]
    gene_DS_R_binned=[]
    gene_B_R_binned=[]    
    
    for i in range(len(gene_US_F)):
        gene_US_F_binned.append(Binning(gene_US_F[i], bin_width))
        gene_DS_F_binned.append(Binning(gene_DS_F[i], bin_width))
        gene_B_F_binned.append(Binning(gene_B_F[i], bin_width))
        gene_US_R_binned.append(Binning(gene_US_R[i], bin_width))
        gene_DS_R_binned.append(Binning(gene_DS_R[i], bin_width))
        gene_B_R_binned.append(Binning(gene_B_R[i], bin_width))
    
    #Scale GB length.
    print(f'GB F scaling in progress, it takes some time...')
    length_binned=int(length/bin_width)
    gene_B_F_binned_sc=[]
    for gene in gene_B_F_binned:
        gene_B_F_binned_sc.append(scale_gene_body(gene, length_binned))

    gene_B_R_binned_sc=[]
    print(f'GB R scaling in progress, it takes some time...')
    for gene in gene_B_R_binned:
        gene_B_R_binned_sc.append(scale_gene_body(gene, length_binned))
       
    #Calculate mean, std.
    win_width_binned=int(win_width/bin_width)
    
    gene_US_F_binned_mean=[]
    gene_DS_F_binned_mean=[]
    gene_US_R_binned_mean=[]
    gene_DS_R_binned_mean=[]
    gene_US_F_binned_std=[]
    gene_DS_F_binned_std=[]
    gene_US_R_binned_std=[]
    gene_DS_R_binned_std=[]  
  
    for i in range(win_width_binned):
        gene_US_F_bin_ar=[]
        gene_DS_F_bin_ar=[]
        gene_US_R_bin_ar=[]
        gene_DS_R_bin_ar=[]
       
        for j in range(len(gene_US_F_binned)):
            gene_US_F_bin_ar.append(gene_US_F_binned[j][i])
            gene_DS_F_bin_ar.append(gene_DS_F_binned[j][i])
            gene_US_R_bin_ar.append(gene_US_R_binned[j][i])
            gene_DS_R_bin_ar.append(gene_DS_R_binned[j][i])
        
        gene_US_F_binned_mean.append(np.mean(gene_US_F_bin_ar))
        gene_US_F_binned_std.append(np.std(gene_US_F_bin_ar))
        gene_DS_F_binned_mean.append(np.mean(gene_DS_F_bin_ar))
        gene_DS_F_binned_std.append(np.std(gene_DS_F_bin_ar))    
        gene_US_R_binned_mean.append(np.mean(gene_US_R_bin_ar))
        gene_US_R_binned_std.append(np.std(gene_US_R_bin_ar))    
        gene_DS_R_binned_mean.append(np.mean(gene_DS_R_bin_ar))
        gene_DS_R_binned_std.append(np.std(gene_DS_R_bin_ar))
        
    gene_B_F_binned_mean=[]
    gene_B_R_binned_mean=[]
    gene_B_F_binned_std=[]
    gene_B_R_binned_std=[]   
    
    for i in range(length_binned):
        gene_B_F_bin_ar=[]
        gene_B_R_bin_ar=[]
        
        for j in range(len(gene_B_F_binned_sc)):
            gene_B_F_bin_ar.append(gene_B_F_binned_sc[j][i])
            gene_B_R_bin_ar.append(gene_B_R_binned_sc[j][i])
        
        gene_B_F_binned_mean.append(np.mean(gene_B_F_bin_ar))
        gene_B_F_binned_std.append(np.std(gene_B_F_bin_ar))
        gene_B_R_binned_mean.append(np.mean(gene_B_R_bin_ar))
        gene_B_R_binned_std.append(np.std(gene_B_R_bin_ar))  
  

    #Write wig-like file with FE over US, GB, DS.
    print(f'Writing FE over TU, GB, DS...')
    gene_F_binned_mean=np.concatenate((gene_US_F_binned_mean, gene_B_F_binned_mean, gene_DS_F_binned_mean), axis=None)
    gene_R_binned_mean=np.concatenate((gene_US_R_binned_mean, gene_B_R_binned_mean, gene_DS_R_binned_mean), axis=None)
    gene_F_binned_std=np.concatenate((gene_US_F_binned_std,  gene_B_F_binned_std,  gene_DS_F_binned_std),  axis=None)
    gene_R_binned_std=np.concatenate((gene_US_R_binned_std,  gene_B_R_binned_std,  gene_DS_R_binned_std),  axis=None)
    write_wig(gene_F_binned_mean, f'{out_path}\Signal_of_TUs_wig\\{genes_set_name}\Mean_signal_{FE_track_pair_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp_F.wig', f'{win_width}_{length}_{bin_width}')
    write_wig(gene_R_binned_mean, f'{out_path}\Signal_of_TUs_wig\\{genes_set_name}\Mean_signal_{FE_track_pair_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp_R.wig', f'{win_width}_{length}_{bin_width}')
    write_wig(gene_F_binned_std, f'{out_path}\Signal_of_TUs_wig\\{genes_set_name}\STD_signal_{FE_track_pair_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp_F.wig',  f'{win_width}_{length}_{bin_width}')
    write_wig(gene_R_binned_std, f'{out_path}\Signal_of_TUs_wig\\{genes_set_name}\STD_signal_{FE_track_pair_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp_R.wig',  f'{win_width}_{length}_{bin_width}')


    #Compute confidential interval borders (+/- 1*SEM).
    Upper_conf_interval_F=np.array(gene_F_binned_mean)+(np.array(gene_F_binned_std)/np.sqrt(Num_genes))
    Lower_conf_interval_F=np.array(gene_F_binned_mean)-(np.array(gene_F_binned_std)/np.sqrt(Num_genes))
    Upper_conf_interval_R=np.array(gene_R_binned_mean)+(np.array(gene_R_binned_std)/np.sqrt(Num_genes))
    Lower_conf_interval_R=np.array(gene_R_binned_mean)-(np.array(gene_R_binned_std)/np.sqrt(Num_genes))     
    
    #Plot FE over US, GB, DS. 
    print(f'Plotting FE over TU, GB, DS...')
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)  
    
    positions_bn=np.arange(-win_width+(bin_width/2), win_width+length-(bin_width/2)+1, bin_width)
    
    #plot1.plot(positions_bn,  gene_F_binned_mean, linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label='Coding strand') 
    #plot1.plot(positions_bn,  Upper_conf_interval_F,  linestyle='-', color='#B6B8BD', linewidth=0.8, alpha=0.6)
    #plot1.plot(positions_bn,  Lower_conf_interval_F,  linestyle='-', color='#B6B8BD', linewidth=0.8, alpha=0.6)
    #plot1.fill_between(positions_bn, Lower_conf_interval_F, Upper_conf_interval_F, facecolor='#43c287', alpha=0.4, interpolate=True) 
    #Make interpolation with splines.
    cs_F=CubicSpline(positions_bn, gene_F_binned_mean)
    #regressionLineOrder=20
    #regressionLine_F=np.polyfit(positions_bn, gene_F_binned_mean, regressionLineOrder)
    #p_F=np.poly1d(regressionLine_F)    
    x_range=np.arange(np.min(positions_bn), np.max(positions_bn), 20)
    plot1.plot(x_range, cs_F(x_range), linestyle='--', color='#B6B8BD', linewidth=2, alpha=1, label='Cubic Spline F')
    
    #plot1.plot(positions_bn,  gene_R_binned_mean, linestyle='-', color='#333738', linewidth=2.5, alpha=1, label='Template strand')
    #plot1.plot(positions_bn,  Upper_conf_interval_R,  linestyle='-', color='#333738', linewidth=0.8, alpha=0.6)   
    #plot1.plot(positions_bn,  Lower_conf_interval_R,  linestyle='-', color='#333738', linewidth=0.8, alpha=0.6)
    #plot1.fill_between(positions_bn, Lower_conf_interval_R, Upper_conf_interval_R, facecolor='#7ce0ff', alpha=0.4, interpolate=True) 
    #Make interpolation with splines.
    cs_R=CubicSpline(positions_bn, gene_R_binned_mean)
    #regressionLine_R=np.polyfit(positions_bn, gene_R_binned_mean, regressionLineOrder)
    #p_R=np.poly1d(regressionLine_R)    
    plot1.plot(x_range, cs_R(x_range), linestyle='--', color='#333738', linewidth=2, alpha=1, label='Cubic Spline R')    
    
    #plot1.set_ylim(-0.6, 0.5) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1)
    plot1.set_xticks([0, length], minor='True')
    plot1.set_yticks([1], minor='True') 
    plot1.axhline(0, color='black', linestyle='--', alpha=0.5, linewidth=0.5, zorder=20) 
    plot1.axvline(0, color='black', linestyle='--', alpha=0.5, linewidth=0.5, zorder=21)
    plot1.axvline(length, color='black', linestyle='--', alpha=0.5, linewidth=0.5, zorder=22)       
    plot1.legend(fontsize=12, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel('TopoI N3E', size=20)  
    plot1.set_xlim(-5000, 10000)
    plt.savefig(f'{out_path}\Figures\Plots\\{genes_set_name}\\{FE_track_pair_name}_over_{genes_set_name}_{win_width}bp_bin_width_{bin_width}bp.png', dpi=400, figsize=(10, 6))   
    plt.close()      
    
    
    #Compare signal in US, GB, DS. Do statistics.
    print(len(gene_US_F_binned), len(gene_B_F_binned_sc), len(gene_DS_F_binned), len(gene_US_R_binned), len(gene_B_R_binned_sc), len(gene_DS_R_binned))   
    us_ds_distance=5000
    bplot_win_width=int(us_ds_distance/bin_width)
    tub_start_end_distance=2000
    tub_bplot_win_width=int(tub_start_end_distance/bin_width)
    
    US_F_means_ar=[]
    DS_F_means_ar=[]
    US_R_means_ar=[]
    DS_R_means_ar=[] 
    GBs_F_means_ar=[]
    GBs_R_means_ar=[] 
    GBe_F_means_ar=[]
    GBe_R_means_ar=[]     
    
    for i in range(len(gene_US_F_binned)):
        US_F_means_ar.append(np.mean(gene_US_F_binned[i][-bplot_win_width-4:-4]))
        DS_F_means_ar.append(np.mean(gene_DS_F_binned[i][4:bplot_win_width+4]))
        US_R_means_ar.append(np.mean(gene_US_R_binned[i][-bplot_win_width-4:-4]))
        DS_R_means_ar.append(np.mean(gene_DS_R_binned[i][4:bplot_win_width+4]))   
        GBs_F_means_ar.append(np.mean(gene_B_F_binned_sc[i][1:tub_bplot_win_width]))
        GBs_R_means_ar.append(np.mean(gene_B_R_binned_sc[i][1:tub_bplot_win_width]))  
        GBe_F_means_ar.append(np.mean(gene_B_F_binned_sc[i][-tub_bplot_win_width-1:-1]))
        GBe_R_means_ar.append(np.mean(gene_B_R_binned_sc[i][-tub_bplot_win_width-1:-1]))  
        
        
        
    #Compare US, TSS, TUB.
    fig, plot_av=plt.subplots(1,1,figsize=(4,3), dpi=100)

    Data_ar=[US_F_means_ar, GBs_F_means_ar, GBe_F_means_ar, DS_F_means_ar, US_R_means_ar, GBs_R_means_ar, GBe_R_means_ar, DS_R_means_ar]
    Mean_data_ar=[np.mean(x) for x in Data_ar]
    StEr_data_ar=[compute_standard_error(x) for x in Data_ar]
    print(Mean_data_ar)
    print(np.mean(US_F_means_ar), np.std(US_F_means_ar), np.std(US_F_means_ar)/np.sqrt(len(US_F_means_ar)))
    
    Conditions=['US', 'TU\nstart\nreg', 'TU\nend\nreg', 'DS', 'US', 'TU\nstart\nreg', 'TU\nend\nreg', 'DS']
    
    color_list=['#6a65c7', '#6a65c7', '#6a65c7', '#6a65c7', '#c96458', '#c96458', '#c96458', '#c96458']
    X_coords=[1,2,3,4,6,7,8,9]
    X_coords_m=[1,2,3,4,6,7,8,9]
      
    Bars=plot_av.bar(X_coords, Mean_data_ar, yerr=StEr_data_ar, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.9, color=color_list, edgecolor='k', linewidth=1)
    plot_av.set_ylabel('EcTopoI N3E', size=16)
    plot_av.set_xticks(X_coords_m)
    plot_av.set_xticklabels(Conditions, rotation=0, size=9)     
    plt.legend((Bars[0],Bars[4]), ('Coding strand', 'Template strand'), frameon=False, loc='best', markerscale=1, handlelength=0.7, handletextpad=0.3)  
    plt.tight_layout()
    plt.show()
    plt.savefig(f'{out_path}\Figures\Plots\\{genes_set_name}\\US_GB_DS_mean_{FE_track_pair_name}_over_{genes_set_name}_{int(bplot_win_width*bin_width)}bp_bin_width_{bin_width}bp.svg', dpi=300, size=(8,3), transparent=True)   
    plt.close() 
    
    #Welch t-test.
    #For F strand.
    Data_ar_F=Data_ar[:4]
    for i in range(len(Data_ar_F)):
        for j in range(len(Data_ar_F)):
            if j>i:
                Intervals_stat=stats.ttest_ind(Data_ar_F[i], Data_ar_F[j], equal_var=False)
                print(f'Sample size: {len(Data_ar_F[i])}, Sample size: {len(Data_ar_F[j])}')
                print(f'\nT-test N3E, F strand Mean1={round(np.mean(Data_ar_F[i]),3)}; Mean2={round(np.mean(Data_ar_F[j]),3)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n') 
    
    #For R strand.           
    Data_ar_R=Data_ar[4:]
    for i in range(len(Data_ar_R)):
        for j in range(len(Data_ar_R)):
            if j>i:
                Intervals_stat=stats.ttest_ind(Data_ar_R[i], Data_ar_R[j], equal_var=False)
                print(f'Sample size: {len(Data_ar_R[i])}, Sample size: {len(Data_ar_R[j])}')
                print(f'\nT-test N3E, R strand Mean1={round(np.mean(Data_ar_R[i]),3)}; Mean2={round(np.mean(Data_ar_R[j]),3)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n') 
    
    return gene_US_F_binned, gene_B_F_binned_sc, gene_DS_F_binned, gene_US_R_binned, gene_B_R_binned_sc, gene_DS_R_binned


#######
#Returns FE or Ded FE over the set of genes (US, GB, DS) - for each gene separately. For non-strand-specific data.
#######

def genes_and_FE_nss(gene_annotation, genes_set_name, FE_track, FE_track_name, out_path, deletions_inpath, win_width, length, bin_width):
    
    #Parsing deletions
    deletions=deletions_info(deletions_inpath) 

    #Calculate FE over genes.
    gene_US=[]
    gene_DS=[]
    gene_B=[]
    gene_US_mean_dict={}
    gene_DS_mean_dict={}
    gene_B_mean_dict={} 
    
    #Number of genes.
    Num_genes=len(gene_annotation)
    
    for gene_name, gene_info in gene_annotation.items():
        delited=0
        for deletion in deletions:
            if deletion[1]>=gene_info[0]>=deletion[0] or deletion[1]>=gene_info[1]>=deletion[0]:
                delited=1
        if delited==0:
            start=gene_info[0]
            end=gene_info[1]
            glen=len(FE_track)
            if gene_info[2]=='+':
                if start<win_width:
                    gene_US.append(FE_track[glen-(win_width-start):] + FE_track[:start])
                    gene_US_mean_dict[gene_name]=[np.mean(FE_track[glen-(win_width-start):] + FE_track[:start]), start, end, gene_info[2]]                   
                else:
                    gene_US.append(FE_track[start-win_width:start])
                    gene_US_mean_dict[gene_name]=[np.mean(FE_track[start-win_width:start]), start, end, gene_info[2]]                    
                if end+win_width>glen:
                    gene_DS.append(FE_track[end:] + FE_track[:end+win_width-glen])
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[end:] + FE_track[:end+win_width-glen]), start, end, gene_info[2]]                    
                else:
                    gene_DS.append(FE_track[end:end+win_width])
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[end:end+win_width]), start, end, gene_info[2]]
                
                gene_B.append(FE_track[start:end])
                gene_B_mean_dict[gene_name]=[np.mean(FE_track[start:end]), start, end, gene_info[2]]               
                
            elif gene_info[2]=='-':
                if start<win_width:
                    gene_DS.append(FE_track[:start][::-1] + FE_track[glen-(win_width-start):][::-1])                  
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[:start][::-1] + FE_track[glen-(win_width-start):][::-1]), start, end, gene_info[2]]
                else:
                    gene_DS.append(FE_track[start-win_width:start][::-1])                                               
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[start-win_width:start][::-1]), start, end, gene_info[2]]                  
                if end+win_width>glen:
                    gene_US.append(FE_track[:end+win_width-glen][::-1] + FE_track[end:][::-1])                         
                    gene_US_mean_dict[gene_name]=[np.mean(FE_track[:end+win_width-glen][::-1] + FE_track[end:][::-1]), start, end, gene_info[2]]               
                else:
                    gene_US.append(FE_track[end:end+win_width][::-1])                                                   
                    gene_US_mean_dict[gene_name]=[np.mean(FE_track[end:end+win_width][::-1]), start, end, gene_info[2]]                                                                                   
                
                gene_B.append(FE_track[start:end][::-1])
                gene_B_mean_dict[gene_name]=[np.mean(FE_track[start:end][::-1]), start, end, gene_info[2]]                                                         
                
    #Data binning.
    print(len(gene_US), len(gene_B), len(gene_DS))
    
    gene_US_binned=[]
    gene_DS_binned=[]
    gene_B_binned=[]   
    
    for i in range(len(gene_US)):
        gene_US_binned.append(Binning(gene_US[i], bin_width))
        gene_DS_binned.append(Binning(gene_DS[i], bin_width))
        gene_B_binned.append(Binning(gene_B[i], bin_width))
    
    #Scale GB length.
    print(f'GB F scaling in progress, it takes some time...')
    length_binned=int(length/bin_width)
    gene_B_binned_sc=[]
    for gene in gene_B_binned:
        gene_B_binned_sc.append(scale_gene_body(gene, length_binned))
       
    #Calculate mean, std.
    win_width_binned=int(win_width/bin_width)
    
    gene_US_binned_mean=[]
    gene_DS_binned_mean=[]
    gene_US_binned_std=[]
    gene_DS_binned_std=[]
  
    for i in range(win_width_binned):
        gene_US_bin_ar=[]
        gene_DS_bin_ar=[]
       
        for j in range(len(gene_US_binned)):
            gene_US_bin_ar.append(gene_US_binned[j][i])
            gene_DS_bin_ar.append(gene_DS_binned[j][i])
        
        gene_US_binned_mean.append(np.mean(gene_US_bin_ar))
        gene_US_binned_std.append(np.std(gene_US_bin_ar))
        gene_DS_binned_mean.append(np.mean(gene_DS_bin_ar))
        gene_DS_binned_std.append(np.std(gene_DS_bin_ar))    
        
    gene_B_binned_mean=[]
    gene_B_binned_std=[]   
    
    for i in range(length_binned):
        gene_B_bin_ar=[]
        
        for j in range(len(gene_B_binned_sc)):
            gene_B_bin_ar.append(gene_B_binned_sc[j][i])
        
        gene_B_binned_mean.append(np.mean(gene_B_bin_ar))
        gene_B_binned_std.append(np.std(gene_B_bin_ar))  

    #Write wig-like file with FE over US, GB, DS.
    print(f'Writing FE over TU, GB, DS...')
    gene_binned_mean=np.concatenate((gene_US_binned_mean, gene_B_binned_mean, gene_DS_binned_mean), axis=None)
    gene_binned_std=np.concatenate((gene_US_binned_std,  gene_B_binned_std,  gene_DS_binned_std),  axis=None)
    write_wig(gene_binned_mean, f'{out_path}\Signal_of_TUs_wig\\{genes_set_name}\Mean_signal_{FE_track_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp.wig', f'{win_width}_{length}_{bin_width}')
    write_wig(gene_binned_std,  f'{out_path}\Signal_of_TUs_wig\\{genes_set_name}\STD_signal_{FE_track_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp.wig',  f'{win_width}_{length}_{bin_width}')


    #Plot FE over US, GB, DS. 
    Upper_conf_interval=np.array(gene_binned_mean)+(np.array(gene_binned_std)/np.sqrt(Num_genes))
    Lower_conf_interval=np.array(gene_binned_mean)-(np.array(gene_binned_std)/np.sqrt(Num_genes))   
    
    print(f'Plotting FE over TU, GB, DS...')
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)  
    
    positions_bn=np.arange(-win_width+(bin_width/2), win_width+length-(bin_width/2)+1, bin_width)
    
    plot1.plot(positions_bn,  gene_binned_mean, linestyle='-', color='#c44733', linewidth=2.5, alpha=1, label='Coding strand') 
    plot1.plot(positions_bn,  Upper_conf_interval,  linestyle='-', color='#c44733', linewidth=0.8, alpha=0.6)
    plot1.plot(positions_bn,  Lower_conf_interval,  linestyle='-', color='#c44733', linewidth=0.8, alpha=0.6)
    plot1.fill_between(positions_bn, Lower_conf_interval, Upper_conf_interval, facecolor='#dca0ff', alpha=0.4, interpolate=True)       
    
    #plot1.set_ylim(-0.6, 0.5) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1)
    plot1.set_xticks([0, length], minor='True')
    plot1.set_yticks([1], minor='True') 
    plot1.axhline(1, color='black', linestyle='--', alpha=0.5, linewidth=0.5, zorder=20) 
    plot1.axvline(0, color='black', linestyle='--', alpha=0.5, linewidth=0.5, zorder=21)
    plot1.axvline(length, color='black', linestyle='--', alpha=0.5, linewidth=0.5, zorder=22)       
    plot1.legend(fontsize=12, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel('TopoI FE', size=20)    
    plt.savefig(f'{out_path}\Figures\Plots\\{genes_set_name}\\{FE_track_name}_over_{genes_set_name}_{win_width}bp_bin_width_{bin_width}bp.png', dpi=400, figsize=(10, 6))   
    plt.close()      
  
    return gene_US_binned, gene_B_binned_sc, gene_DS_binned


#######
#Normalizes tracks and divides strand-specific tracks by a non-strand-specific track.
#######

def genes_FE_ss_vs_nss(gene_annotation, genes_set_name, ss_track_pair_data, ss_track_pair_name, track_data, track_name, out_path, deletions_inpath, win_width, length, bin_width):
    
    #Strand-specific tracks.
    gene_US_F_binned, gene_B_F_binned_sc, gene_DS_F_binned, gene_US_R_binned, gene_B_R_binned_sc, gene_DS_R_binned=ss_track_pair_data
    #Non strand-specific tracks.
    gene_US_binned, gene_B_binned_sc, gene_DS_binned=track_data
    
    #Normalization parameters: mean and std (for 200bp binning window width).
    ss_mean_F=0.0032  #200:0.0032;500:0.0032 
    ss_std_F=0.121    #200:0.121 ;500:0.089
    ss_mean_R=0.0027  #200:0.0027;500:0.0027
    ss_std_R=0.121    #200:0.121 ;500:0.089
    nss_mean=1.26     #200:1.26  ;500:1.26
    nss_std=2.20      #200:2.20  ;500:2.34
    
    #Number of genes.
    Num_genes=len(gene_annotation)
    
    #Normalize and divide tracks, compute mean and std.
    #For US and DS regions.
    win_width_binned=int(win_width/bin_width)
    
    gene_US_F_binned_div_mean=[]
    gene_DS_F_binned_div_mean=[]
    gene_US_R_binned_div_mean=[]
    gene_DS_R_binned_div_mean=[]
    gene_US_F_binned_div_std=[]
    gene_DS_F_binned_div_std=[]
    gene_US_R_binned_div_std=[]
    gene_DS_R_binned_div_std=[]    
    
    for i in range(win_width_binned):
        gene_US_F_bin_div_ar=[]
        gene_DS_F_bin_div_ar=[]
        gene_US_R_bin_div_ar=[]
        gene_DS_R_bin_div_ar=[]        
       
        for j in range(len(gene_US_F_binned)):  
            gene_US_F_bin_div_ar.append((((gene_US_F_binned[j][i]-ss_mean_F)/ss_std_F)+1)/(((gene_US_binned[j][i]-nss_mean)/nss_std)+1))
            gene_US_R_bin_div_ar.append((((gene_US_R_binned[j][i]-ss_mean_R)/ss_std_R)+1)/(((gene_US_binned[j][i]-nss_mean)/nss_std)+1))
            gene_DS_F_bin_div_ar.append((((gene_DS_F_binned[j][i]-ss_mean_F)/ss_std_F)+1)/(((gene_DS_binned[j][i]-nss_mean)/nss_std)+1))
            gene_DS_R_bin_div_ar.append((((gene_DS_R_binned[j][i]-ss_mean_R)/ss_std_R)+1)/(((gene_DS_binned[j][i]-nss_mean)/nss_std)+1))   
            
        gene_US_F_binned_div_mean.append(np.mean(gene_US_F_bin_div_ar))
        gene_US_F_binned_div_std.append(np.std(gene_US_F_bin_div_ar))
        gene_DS_F_binned_div_mean.append(np.mean(gene_DS_F_bin_div_ar))
        gene_DS_F_binned_div_std.append(np.std(gene_DS_F_bin_div_ar))    
        gene_US_R_binned_div_mean.append(np.mean(gene_US_R_bin_div_ar))
        gene_US_R_binned_div_std.append(np.std(gene_US_R_bin_div_ar))    
        gene_DS_R_binned_div_mean.append(np.mean(gene_DS_R_bin_div_ar))
        gene_DS_R_binned_div_std.append(np.std(gene_DS_R_bin_div_ar))    
    
    #For gene body regions.     
    length_binned=int(length/bin_width)
    
    gene_B_F_binned_div_mean=[]
    gene_B_R_binned_div_mean=[]
    gene_B_F_binned_div_std=[]
    gene_B_R_binned_div_std=[]   
    
    for i in range(length_binned):
        gene_B_F_bin_div_ar=[]
        gene_B_R_bin_div_ar=[]
        
        for j in range(len(gene_B_F_binned_sc)):
            gene_B_F_bin_div_ar.append((((gene_B_F_binned_sc[j][i]-ss_mean_F)/ss_std_F)+1)/(((gene_B_binned_sc[j][i]-nss_mean)/nss_std)+1))
            gene_B_R_bin_div_ar.append((((gene_B_R_binned_sc[j][i]-ss_mean_R)/ss_std_R)+1)/(((gene_B_binned_sc[j][i]-nss_mean)/nss_std)+1))
        
        gene_B_F_binned_div_mean.append(np.mean(gene_B_F_bin_div_ar))
        gene_B_F_binned_div_std.append(np.std(gene_B_F_bin_div_ar))
        gene_B_R_binned_div_mean.append(np.mean(gene_B_R_bin_div_ar))
        gene_B_R_binned_div_std.append(np.std(gene_B_R_bin_div_ar)) 
    
       
    #Normalize tracks, calculate mean, std.
    #For US and DS regions.
    gene_US_F_binned_mean=[]
    gene_DS_F_binned_mean=[]
    gene_US_R_binned_mean=[]
    gene_DS_R_binned_mean=[]
    gene_US_F_binned_std=[]
    gene_DS_F_binned_std=[]
    gene_US_R_binned_std=[]
    gene_DS_R_binned_std=[]  
    gene_US_binned_mean=[]
    gene_US_binned_std=[]
    gene_DS_binned_mean=[]
    gene_DS_binned_std=[]
  
    for i in range(win_width_binned):
        gene_US_F_bin_ar=[]
        gene_DS_F_bin_ar=[]
        gene_US_R_bin_ar=[]
        gene_DS_R_bin_ar=[]
        gene_US_bin_ar=[]
        gene_DS_bin_ar=[]        
       
        for j in range(len(gene_US_F_binned)):
            gene_US_F_bin_ar.append(((gene_US_F_binned[j][i]-ss_mean_F)/ss_std_F)+1)
            gene_DS_F_bin_ar.append(((gene_DS_F_binned[j][i]-ss_mean_F)/ss_std_F)+1)
            gene_US_R_bin_ar.append(((gene_US_R_binned[j][i]-ss_mean_R)/ss_std_R)+1)
            gene_DS_R_bin_ar.append(((gene_DS_R_binned[j][i]-ss_mean_R)/ss_std_R)+1)
            gene_US_bin_ar.append(((gene_US_binned[j][i]-nss_mean)/nss_std)+1)
            gene_DS_bin_ar.append(((gene_DS_binned[j][i]-nss_mean)/nss_std)+1)            
        
        gene_US_F_binned_mean.append(np.mean(gene_US_F_bin_ar))
        gene_US_F_binned_std.append(np.std(gene_US_F_bin_ar))
        gene_DS_F_binned_mean.append(np.mean(gene_DS_F_bin_ar))
        gene_DS_F_binned_std.append(np.std(gene_DS_F_bin_ar))    
        gene_US_R_binned_mean.append(np.mean(gene_US_R_bin_ar))
        gene_US_R_binned_std.append(np.std(gene_US_R_bin_ar))    
        gene_DS_R_binned_mean.append(np.mean(gene_DS_R_bin_ar))
        gene_DS_R_binned_std.append(np.std(gene_DS_R_bin_ar))
        gene_US_binned_mean.append(np.mean(gene_US_bin_ar))
        gene_US_binned_std.append(np.std(gene_US_bin_ar))    
        gene_DS_binned_mean.append(np.mean(gene_DS_bin_ar))
        gene_DS_binned_std.append(np.std(gene_DS_bin_ar))        
    
    #For gene body regions.    
    gene_B_F_binned_mean=[]
    gene_B_R_binned_mean=[]
    gene_B_F_binned_std=[]
    gene_B_R_binned_std=[] 
    gene_B_binned_mean=[]
    gene_B_binned_std=[]    
    
    for i in range(length_binned):
        gene_B_F_bin_ar=[]
        gene_B_R_bin_ar=[]
        gene_B_bin_ar=[]
        
        for j in range(len(gene_B_F_binned_sc)):
            gene_B_F_bin_ar.append(((gene_B_F_binned_sc[j][i]-ss_mean_F)/ss_std_F)+1)
            gene_B_R_bin_ar.append(((gene_B_R_binned_sc[j][i]-ss_mean_R)/ss_std_R)+1)
            gene_B_bin_ar.append(((gene_B_binned_sc[j][i]-nss_mean)/nss_std)+1)
        
        gene_B_F_binned_mean.append(np.mean(gene_B_F_bin_ar))
        gene_B_F_binned_std.append(np.std(gene_B_F_bin_ar))
        gene_B_R_binned_mean.append(np.mean(gene_B_R_bin_ar))
        gene_B_R_binned_std.append(np.std(gene_B_R_bin_ar))  
        gene_B_binned_mean.append(np.mean(gene_B_bin_ar))
        gene_B_binned_std.append(np.std(gene_B_bin_ar))          
  

    #Write wig-like file with FE over US, GB, DS for normalized divided tracks.
    print(f'Writing FE over TU, GB, DS... (normalized divided)')
    gene_F_binned_div_mean=np.concatenate((gene_US_F_binned_div_mean, gene_B_F_binned_div_mean, gene_DS_F_binned_div_mean), axis=None)
    gene_R_binned_div_mean=np.concatenate((gene_US_R_binned_div_mean, gene_B_R_binned_div_mean, gene_DS_R_binned_div_mean), axis=None)
    gene_F_binned_div_std=np.concatenate((gene_US_F_binned_div_std,  gene_B_F_binned_div_std,  gene_DS_F_binned_div_std),  axis=None)
    gene_R_binned_div_std=np.concatenate((gene_US_R_binned_div_std,  gene_B_R_binned_div_std,  gene_DS_R_binned_div_std),  axis=None)
    write_wig(gene_F_binned_div_mean, f'{out_path}\Signal_of_TUs_wig\\{genes_set_name}\Divided_mean_normalized_signal_{ss_track_pair_name}_by_{track_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp_F.wig', f'{win_width}_{length}_{bin_width}')
    write_wig(gene_R_binned_div_mean, f'{out_path}\Signal_of_TUs_wig\\{genes_set_name}\Divided_mean_normalized_signal_{ss_track_pair_name}_by_{track_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp_R.wig', f'{win_width}_{length}_{bin_width}')
    write_wig(gene_F_binned_div_std,  f'{out_path}\Signal_of_TUs_wig\\{genes_set_name}\Divided_STD_normalized_signal_{ss_track_pair_name}_by_{track_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp_F.wig',  f'{win_width}_{length}_{bin_width}')
    write_wig(gene_R_binned_div_std,  f'{out_path}\Signal_of_TUs_wig\\{genes_set_name}\Divided_STD_normalized_signal_{ss_track_pair_name}_by_{track_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp_R.wig',  f'{win_width}_{length}_{bin_width}')
     
    #Write wig-like file with FE over US, GB, DS for normalized tracks.
    print(f'Writing FE over TU, GB, DS...(normalized)')    
    gene_F_binned_mean=np.concatenate((gene_US_F_binned_mean, gene_B_F_binned_mean, gene_DS_F_binned_mean), axis=None)
    gene_R_binned_mean=np.concatenate((gene_US_R_binned_mean, gene_B_R_binned_mean, gene_DS_R_binned_mean), axis=None)
    gene_binned_mean=np.concatenate((gene_US_binned_mean, gene_B_binned_mean, gene_DS_binned_mean), axis=None)
    gene_F_binned_std=np.concatenate((gene_US_F_binned_std,  gene_B_F_binned_std,  gene_DS_F_binned_std),  axis=None)
    gene_R_binned_std=np.concatenate((gene_US_R_binned_std,  gene_B_R_binned_std,  gene_DS_R_binned_std),  axis=None)
    gene_binned_std=np.concatenate((gene_US_binned_std,  gene_B_binned_std,  gene_DS_binned_std),  axis=None)
    write_wig(gene_F_binned_mean, f'{out_path}\Signal_of_TUs_wig\\{genes_set_name}\Mean_normalized_signal_{ss_track_pair_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp_F.wig', f'{win_width}_{length}_{bin_width}')
    write_wig(gene_R_binned_mean, f'{out_path}\Signal_of_TUs_wig\\{genes_set_name}\Mean_normalized_signal_{ss_track_pair_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp_R.wig', f'{win_width}_{length}_{bin_width}')
    write_wig(gene_binned_mean,   f'{out_path}\Signal_of_TUs_wig\\{genes_set_name}\Mean_normalized_signal_{track_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp.wig',   f'{win_width}_{length}_{bin_width}')    
    write_wig(gene_F_binned_std,  f'{out_path}\Signal_of_TUs_wig\\{genes_set_name}\STD_normalized_signal_{ss_track_pair_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp_F.wig',  f'{win_width}_{length}_{bin_width}')
    write_wig(gene_R_binned_std,  f'{out_path}\Signal_of_TUs_wig\\{genes_set_name}\STD_normalized_signal_{ss_track_pair_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp_R.wig',  f'{win_width}_{length}_{bin_width}')
    write_wig(gene_binned_std,    f'{out_path}\Signal_of_TUs_wig\\{genes_set_name}\STD_normalized_signal_{track_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp.wig',    f'{win_width}_{length}_{bin_width}')
       

    #Plot FE over US, GB, DS. 
    Upper_conf_interval_F=np.array(gene_F_binned_mean)+(np.array(gene_F_binned_std)/np.sqrt(Num_genes))
    Lower_conf_interval_F=np.array(gene_F_binned_mean)-(np.array(gene_F_binned_std)/np.sqrt(Num_genes))
    Upper_conf_interval_R=np.array(gene_R_binned_mean)+(np.array(gene_R_binned_std)/np.sqrt(Num_genes))
    Lower_conf_interval_R=np.array(gene_R_binned_mean)-(np.array(gene_R_binned_std)/np.sqrt(Num_genes))
    
    Upper_conf_interval=np.array(gene_binned_mean)+(np.array(gene_binned_std)/np.sqrt(Num_genes))
    Lower_conf_interval=np.array(gene_binned_mean)-(np.array(gene_binned_std)/np.sqrt(Num_genes))  
    
    Upper_conf_interval_div_F=np.array(gene_F_binned_div_mean)+(np.array(gene_F_binned_div_std)/np.sqrt(Num_genes))
    Lower_conf_interval_div_F=np.array(gene_F_binned_div_mean)-(np.array(gene_F_binned_div_std)/np.sqrt(Num_genes))
    Upper_conf_interval_div_R=np.array(gene_R_binned_div_mean)+(np.array(gene_R_binned_div_std)/np.sqrt(Num_genes))
    Lower_conf_interval_div_R=np.array(gene_R_binned_div_mean)-(np.array(gene_R_binned_div_std)/np.sqrt(Num_genes))    
    
    print(f'Plotting FE over TU, GB, DS...')
    plt.figure(figsize=(7.5, 4.5), dpi=100)
    plot1=plt.subplot(111)  
    
    positions_bn=np.arange(-win_width+(bin_width/2), win_width+length-(bin_width/2)+1, bin_width)
    #F strand normalized.
    plot1.plot(positions_bn,  gene_F_binned_mean, linestyle='-', color='#6a65c7', linewidth=1.5, alpha=1, label='Cleavage (N3E): coding strand') 
    plot1.plot(positions_bn,  Upper_conf_interval_F,  linestyle='-', color='#6a65c7', linewidth=0.8, alpha=0.2)
    plot1.plot(positions_bn,  Lower_conf_interval_F,  linestyle='-', color='#6a65c7', linewidth=0.8, alpha=0.2)
    plot1.fill_between(positions_bn, Lower_conf_interval_F, Upper_conf_interval_F, facecolor='#6a65c7', alpha=0.1, interpolate=True)       
    #R strand normalized.
    plot1.plot(positions_bn,  gene_R_binned_mean, linestyle='-', color='#c96458', linewidth=1.5, alpha=1, label='Cleavage (N3E): template strand')
    plot1.plot(positions_bn,  Upper_conf_interval_R,  linestyle='-', color='#c96458', linewidth=0.8, alpha=0.2)   
    plot1.plot(positions_bn,  Lower_conf_interval_R,  linestyle='-', color='#c96458', linewidth=0.8, alpha=0.2)
    plot1.fill_between(positions_bn, Lower_conf_interval_R, Upper_conf_interval_R, facecolor='#c96458', alpha=0.1, interpolate=True) 
    #No strand normalized.
    plot1.plot(positions_bn,  gene_binned_mean, linestyle='-', color='#333738', linewidth=1.5, alpha=1, label='Binding (FE)')
    plot1.plot(positions_bn,  Upper_conf_interval,  linestyle='-', color='#333738', linewidth=0.8, alpha=0.2)   
    plot1.plot(positions_bn,  Lower_conf_interval,  linestyle='-', color='#333738', linewidth=0.8, alpha=0.2)
    plot1.fill_between(positions_bn, Lower_conf_interval, Upper_conf_interval, facecolor='#333738', alpha=0.05, interpolate=True)  
    #F strand normalized divided.
    #plot1.plot(positions_bn,  gene_F_binned_div_mean, linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label='Cleavage/Binding: coding strand') 
    #plot1.plot(positions_bn,  Upper_conf_interval_div_F,  linestyle='-', color='#B6B8BD', linewidth=0.8, alpha=0.6)
    #plot1.plot(positions_bn,  Lower_conf_interval_div_F,  linestyle='-', color='#B6B8BD', linewidth=0.8, alpha=0.6)
    #plot1.fill_between(positions_bn, Lower_conf_interval_div_F, Upper_conf_interval_div_F, facecolor='#43c287', alpha=0.4, interpolate=True)       
    #F strand normalized divided.
    #plot1.plot(positions_bn,  gene_R_binned_div_mean, linestyle='-', color='#333738', linewidth=2.5, alpha=1, label='Cleavage/Binding: template strand')
    #plot1.plot(positions_bn,  Upper_conf_interval_div_R,  linestyle='-', color='#333738', linewidth=0.8, alpha=0.6)   
    #plot1.plot(positions_bn,  Lower_conf_interval_div_R,  linestyle='-', color='#333738', linewidth=0.8, alpha=0.6)
    #plot1.fill_between(positions_bn, Lower_conf_interval_div_R, Upper_conf_interval_div_R, facecolor='#7ce0ff', alpha=0.4, interpolate=True)    
    
    #plot1.set_ylim(-0.6, 0.5) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='$TU_{start}$'
    ticks_lables[ticks.index(length)]='$TU_{end}$'
    ticks_lables1=ticks_lables[:ticks_lables.index('$TU_{end}$')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1)
    plot1.tick_params(axis='x', which='major', labelsize=15, pad=1)
    plot1.set_xticks([0, length], minor='True')
    plot1.set_yticks([1], minor='True') 
    plot1.axhline(1, color='black', linestyle='--', alpha=0.5, linewidth=0.5, zorder=20) 
    plot1.axvline(0, color='black', linestyle='--', alpha=0.5, linewidth=0.5, zorder=21)
    plot1.axvline(length, color='black', linestyle='--', alpha=0.5, linewidth=0.5, zorder=22) 
    plot1.tick_params(axis='y', which='major', labelsize=15, pad=2)
    plot1.legend(fontsize=10, frameon=False, markerscale=5, handlelength=0.5)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel('TopoI N3E (Topo-Seq), norm\nTopoI FE (ChIP-Seq), norm', size=20)    
    plt.tight_layout()
    plt.savefig(f'{out_path}\Figures\Plots\\{genes_set_name}\\Cleavage_vs_binding_{ss_track_pair_name}_by_{track_name}_over_{genes_set_name}_{win_width}bp_bin_width_{bin_width}bp.png', dpi=400, figsize=(7.5, 4.5), transparent=True)   
    plt.savefig(f'{out_path}\Figures\Plots\\{genes_set_name}\\Cleavage_vs_binding_{ss_track_pair_name}_by_{track_name}_over_{genes_set_name}_{win_width}bp_bin_width_{bin_width}bp.svg', dpi=400, figsize=(7.5, 4.5), transparent=True)       
    plt.close()      
   
    return

#######
#Wrapper: reads data and gene annotation, computes signal over US, GB, DS, plots and writes the output.
#######


def Wrapper_signal_over_TUs(dict_of_wigs_ss_path, dict_of_wigs_path, path_to_annotation, type_of_annot, genes_set_name, deletions_inpath, win_width, length, bin_width, out_path):
    
    #Reads annotation.
    print(f'Now working with {path_to_annotation}')
    if type_of_annot=='gff':
        gene_annotation=parse_gff_annotation(path_to_annotation, deletions_inpath)[0]['Gene']
    elif type_of_annot=='broadPeak':
        gene_annotation=parse_expression_annotation(path_to_annotation)    
    
    #Reads strand-specific input data in wig files.
    dict_of_wigs_ss={}
    for pair_name, pair_dict in dict_of_wigs_ss_path.items():
        wig_F_data=wig_parsing(pair_dict['F'])
        wig_R_data=wig_parsing(pair_dict['R'])
        dict_of_wigs_ss[pair_name]={'F' : wig_F_data, 'R' : wig_R_data}
    
    #Calculate and plot signal over TUs for strand-specific data.
    ss_binned_data_dict={}
    for FE_track_pair_name, FE_track_pair in dict_of_wigs_ss.items():
        ss_binned_data_dict[FE_track_pair_name]=genes_and_FE_ss(gene_annotation, genes_set_name, FE_track_pair, FE_track_pair_name, out_path, deletions_inpath, win_width, length, bin_width)
          
    #Reads input data in wig files (not strand-specific).
    dict_of_wigs={}
    for wig_name, wig_data in dict_of_wigs_path.items():
        dict_of_wigs[wig_name]=wig_parsing(wig_data)
        
    #Calculate and plot signal over TUs (not strand-specific).
    binned_data_dict={}
    for FE_track_name, FE_track in dict_of_wigs.items():
        binned_data_dict[FE_track_name]=genes_and_FE_nss(gene_annotation, genes_set_name, FE_track, FE_track_name, out_path, deletions_inpath, win_width, length, bin_width)
    
    #Normalize tracks, divide strand-specific tracks by not strand-specific.
    genes_FE_ss_vs_nss(gene_annotation, genes_set_name, ss_binned_data_dict['TopoI_Ara_N3E_subtr_mock_subtr_no_Ara'], 'TopoI_Ara_N3E_subtr_mock_subtr_no_Ara', binned_data_dict['TopoI_no_CTD_no_Rif_FE'], 'TopoI_no_CTD_no_Rif_FE', out_path, deletions_inpath, win_width, length, bin_width)
    
    return

Wrapper_signal_over_TUs(Dict_of_wigs_ss_path_2, Dict_of_wigs_path_2, Path_to_annotation_5,  Type_of_annot_5,  Genes_set_name_5,  Deletions_inpath, Win_width, Length, Bin_width, Out_path)
#Wrapper_signal_over_TUs(Dict_of_wigs_ss_path_2, Dict_of_wigs_path_2, Path_to_annotation_6,  Type_of_annot_6,  Genes_set_name_6,  Deletions_inpath, Win_width, Length, Bin_width, Out_path)
Wrapper_signal_over_TUs(Dict_of_wigs_ss_path_2, Dict_of_wigs_path_2, Path_to_annotation_7,  Type_of_annot_7,  Genes_set_name_7,  Deletions_inpath, Win_width, Length, Bin_width, Out_path)
Wrapper_signal_over_TUs(Dict_of_wigs_ss_path_2, Dict_of_wigs_path_2, Path_to_annotation_8,  Type_of_annot_8,  Genes_set_name_8,  Deletions_inpath, Win_width, Length, Bin_width, Out_path)
Wrapper_signal_over_TUs(Dict_of_wigs_ss_path_2, Dict_of_wigs_path_2, Path_to_annotation_9,  Type_of_annot_9,  Genes_set_name_9,  Deletions_inpath, Win_width, Length, Bin_width, Out_path)
#Wrapper_signal_over_TUs(Dict_of_wigs_ss_path_2, Dict_of_wigs_path_2, Path_to_annotation_10, Type_of_annot_10, Genes_set_name_10, Deletions_inpath, Win_width, Length, Bin_width, Out_path)