###############################################
##Dmitry Sutormin, 2021##
##TopoA Topo-Seq analysis##

#Script computes Fold Enrichment (FE) over upstream (US)
#and downstream (DS) regions of transcription units (TUs)
#and over TUs bodies.
#Script is for strand-specific data.
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
Deletions_inpath='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\Scripts_TopoI_Topo-seq\Additional_genome_features\\Deletions_w3110_G_Mu_SGS.broadPeak'
#Width of US, DS regions.
Win_width=15000
#Length of GB.
Length=5000

#Dictionary of pathes to input data.
Dict_of_wigs_path={'TopoI_Ara_FE' :     {'F' : PWD + 'FE_cov_depth_masked_av\TopoI_Ara_f_FE_av_123.wig', 'R' : PWD + 'FE_cov_depth_masked_av\TopoI_Ara_r_FE_av_123.wig'},
                   'TopoI_FE' :         {'F' : PWD + 'FE_cov_depth_masked_av\TopoI_f_FE_av_123.wig',     'R' : PWD + 'FE_cov_depth_masked_av\TopoI_r_FE_av_123.wig'}
                   }
Dict_of_wigs_path_1={'TopoI_Ara_FE_subtr' :     {'F' : PWD + 'FE_cov_depth_masked_av_subtr_no_Ara\TopoI_Ara_f_FE_av_123_subtr_mock.wig', 'R' : PWD + 'FE_cov_depth_masked_av_subtr_no_Ara\TopoI_Ara_r_FE_av_123_subtr_mock.wig'},
                     }
Dict_of_wigs_path_2={'TopoI_Ara_N3E_subtr_mock_subtr_no_Ara' :     {'F' : PWD + 'WIG_NE_strand_specific_masked_scaled_av_masked_accB_subtract_mock_subtract_no_Ara\TopoI_Ara_N3E_F_masked_scaled_av_123_subtr_mock_subtr_no_Ara.wig', 'R' : PWD + 'WIG_NE_strand_specific_masked_scaled_av_masked_accB_subtract_mock_subtract_no_Ara\TopoI_Ara_N3E_R_masked_scaled_av_123_subtr_mock_subtr_no_Ara.wig'},
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
Out_path='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\Metagene_analysis\Transcripts_strand_specific_N3E_FE\\'

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
#create_out_dirs(Out_path, Genes_set_name_5)
#create_out_dirs(Out_path, Genes_set_name_6)
#create_out_dirs(Out_path, Genes_set_name_7)
#create_out_dirs(Out_path, Genes_set_name_8)
#create_out_dirs(Out_path, Genes_set_name_9)
create_out_dirs(Out_path, Genes_set_name_10)


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
#Returns FE or Ded FE over the set of genes (US, GB, DS) - for each gene separately.
#######

def genes_and_FE(gene_annotation, genes_set_name, FE_track_pair, FE_track_pair_name, out_path, deletions_inpath, win_width, length):
    
    #Strand-specific tracks.
    FE_track_F=FE_track_pair['F']
    FE_track_R=FE_track_pair['R']
    
    #Parsing deletions
    deletions=deletions_info(deletions_inpath) 

    #Calculate FE over genes.
    gene_US_F=np.array([0.0]*win_width)
    gene_DS_F=np.array([0.0]*win_width)
    gene_B_F=[]
    gene_US_F_mean_dict={}
    gene_DS_F_mean_dict={}
    gene_B_F_mean_dict={}
    gene_US_R=np.array([0.0]*win_width)
    gene_DS_R=np.array([0.0]*win_width)
    gene_B_R=[]
    gene_US_R_mean_dict={}
    gene_DS_R_mean_dict={}
    gene_B_R_mean_dict={}   
    
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
                    gene_US_F+=np.array(FE_track_F[glen-(win_width-start):] + FE_track_F[:start])
                    gene_US_F_mean_dict[gene_name]=[np.mean(FE_track_F[glen-(win_width-start):] + FE_track_F[:start]), start, end, gene_info[2]]
                    gene_US_R+=np.array(FE_track_R[glen-(win_width-start):] + FE_track_R[:start])
                    gene_US_R_mean_dict[gene_name]=[np.mean(FE_track_R[glen-(win_width-start):] + FE_track_R[:start]), start, end, gene_info[2]]                   
                else:
                    gene_US_F+=np.array(FE_track_F[start-win_width:start])
                    gene_US_F_mean_dict[gene_name]=[np.mean(FE_track_F[start-win_width:start]), start, end, gene_info[2]]
                    gene_US_R+=np.array(FE_track_R[start-win_width:start])
                    gene_US_R_mean_dict[gene_name]=[np.mean(FE_track_R[start-win_width:start]), start, end, gene_info[2]]                    
                if end+win_width>glen:
                    gene_DS_F+=np.array(FE_track_F[end:] + FE_track_F[:end+win_width-glen])
                    gene_DS_F_mean_dict[gene_name]=[np.mean(FE_track_F[end:] + FE_track_F[:end+win_width-glen]), start, end, gene_info[2]]
                    gene_DS_R+=np.array(FE_track_R[end:] + FE_track_R[:end+win_width-glen])
                    gene_DS_R_mean_dict[gene_name]=[np.mean(FE_track_R[end:] + FE_track_R[:end+win_width-glen]), start, end, gene_info[2]]                    
                else:
                    gene_DS_F+=np.array(FE_track_F[end:end+win_width])
                    gene_DS_F_mean_dict[gene_name]=[np.mean(FE_track_F[end:end+win_width]), start, end, gene_info[2]]
                    gene_DS_R+=np.array(FE_track_R[end:end+win_width])
                    gene_DS_R_mean_dict[gene_name]=[np.mean(FE_track_R[end:end+win_width]), start, end, gene_info[2]] 
                gene_B_F.append(FE_track_F[start:end])
                gene_B_F_mean_dict[gene_name]=[np.mean(FE_track_F[start:end]), start, end, gene_info[2]]
                gene_B_R.append(FE_track_R[start:end])
                gene_B_R_mean_dict[gene_name]=[np.mean(FE_track_R[start:end]), start, end, gene_info[2]]                
                
            elif gene_info[2]=='-':
                if start<win_width:
                    gene_DS_F+=np.array(FE_track_R[:start][::-1] + FE_track_R[glen-(win_width-start):][::-1])                   #Take care of strand-specificity.
                    gene_DS_F_mean_dict[gene_name]=[np.mean(FE_track_R[:start][::-1] + FE_track_R[glen-(win_width-start):][::-1]), start, end, gene_info[2]] #Take care of strand-specificity.
                    gene_DS_R+=np.array(FE_track_F[:start][::-1] + FE_track_F[glen-(win_width-start):][::-1])                   #Take care of strand-specificity.
                    gene_DS_R_mean_dict[gene_name]=[np.mean(FE_track_F[:start][::-1] + FE_track_F[glen-(win_width-start):][::-1]), start, end, gene_info[2]] #Take care of strand-specificity.                    
                else:
                    gene_DS_F+=np.array(FE_track_R[start-win_width:start][::-1])                                                #Take care of strand-specificity.
                    gene_DS_F_mean_dict[gene_name]=[np.mean(FE_track_R[start-win_width:start][::-1]), start, end, gene_info[2]] #Take care of strand-specificity.
                    gene_DS_R+=np.array(FE_track_F[start-win_width:start][::-1])                                                #Take care of strand-specificity.
                    gene_DS_R_mean_dict[gene_name]=[np.mean(FE_track_F[start-win_width:start][::-1]), start, end, gene_info[2]] #Take care of strand-specificity.                    
                if end+win_width>glen:
                    gene_US_F+=np.array(FE_track_R[:end+win_width-glen][::-1] + FE_track_R[end:][::-1])                         #Take care of strand-specificity.
                    gene_US_F_mean_dict[gene_name]=[np.mean(FE_track_R[:end+win_width-glen][::-1] + FE_track_R[end:][::-1]), start, end, gene_info[2]] #Take care of strand-specificity.
                    gene_US_R+=np.array(FE_track_F[:end+win_width-glen][::-1] + FE_track_F[end:][::-1])                         #Take care of strand-specificity.
                    gene_US_R_mean_dict[gene_name]=[np.mean(FE_track_F[:end+win_width-glen][::-1] + FE_track_F[end:][::-1]), start, end, gene_info[2]] #Take care of strand-specificity.                    
                else:
                    gene_US_F+=np.array(FE_track_R[end:end+win_width][::-1])                                                    #Take care of strand-specificity.
                    gene_US_F_mean_dict[gene_name]=[np.mean(FE_track_R[end:end+win_width][::-1]), start, end, gene_info[2]]     #Take care of strand-specificity.
                    gene_US_R+=np.array(FE_track_F[end:end+win_width][::-1])                                                    #Take care of strand-specificity.
                    gene_US_R_mean_dict[gene_name]=[np.mean(FE_track_F[end:end+win_width][::-1]), start, end, gene_info[2]]     #Take care of strand-specificity.                    
                gene_B_F_ar=np.array(FE_track_R[start:end][::-1])                                                               #Take care of strand-specificity.
                gene_B_F.append(gene_B_F_ar.tolist())
                gene_B_F_mean_dict[gene_name]=[np.mean(gene_B_F_ar), start, end, gene_info[2]]
                gene_B_R_ar=np.array(FE_track_F[start:end][::-1])                                                               #Take care of strand-specificity.
                gene_B_R.append(gene_B_R_ar.tolist())
                gene_B_R_mean_dict[gene_name]=[np.mean(gene_B_R_ar), start, end, gene_info[2]]
                
    Num_genes=len(gene_annotation)
    gene_US_F=gene_US_F/Num_genes
    gene_DS_F=gene_DS_F/Num_genes
    gene_US_R=gene_US_R/Num_genes
    gene_DS_R=gene_DS_R/Num_genes    
    print(f'FE over TUs computed...')

    #Scale GB length.
    print(f'GB F scaling in progress, it takes some time...')
    gene_body_F=np.array([0.0]*length)
    for gene in gene_B_F:
        gene_scaled=np.array(scale_gene_body(gene, length))
        gene_body_F+=gene_scaled
    gene_B_F=gene_body_F/Num_genes
    
    print(f'GB R scaling in progress, it takes some time...')
    gene_body_R=np.array([0.0]*length)
    for gene in gene_B_R:
        gene_scaled=np.array(scale_gene_body(gene, length))
        gene_body_R+=gene_scaled
    gene_B_R=gene_body_R/Num_genes    

    #Write wig-like file with FE over US, GB, DS.
    print(f'Writing FE over TU, GB, DS...')
    write_wig(np.concatenate((gene_US_F, gene_B_F, gene_DS_F), axis=None), f'{out_path}\Signal_of_TUs_wig\\{genes_set_name}\Signal_{FE_track_pair_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_F.wig', f'{win_width}_{length}')
    write_wig(np.concatenate((gene_US_R, gene_B_R, gene_DS_R), axis=None), f'{out_path}\Signal_of_TUs_wig\\{genes_set_name}\Signal_{FE_track_pair_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_R.wig', f'{win_width}_{length}')

    #Plot FE over US, GB, DS. 
    print(f'Plotting FE over TU, GB, DS...')
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)  
    positions_DS=np.arange(1+length,win_width+1+length,1)
    positions_US=np.arange(-1-win_width,-1,1)  
    positions_Body=np.arange(0,length,1)
    plot1.plot(positions_US,  gene_US_F, linestyle='-', color='#D17E7E', linewidth=1, label='F') 
    plot1.plot(positions_DS,  gene_DS_F, linestyle='-', color='#D17E7E', linewidth=1)
    plot1.plot(positions_Body, gene_B_F, linestyle='-', color='#D17E7E', linewidth=2)
    plot1.plot(positions_US,  gene_US_R, linestyle='-', color='#7996d4', linewidth=1, label='R') 
    plot1.plot(positions_DS,  gene_DS_R, linestyle='-', color='#7996d4', linewidth=1)
    plot1.plot(positions_Body, gene_B_R, linestyle='-', color='#7996d4', linewidth=2)    
    #plot1.set_ylim(-0.6, 0.5) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1)
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=0.5)
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=0.5)       
    plot1.legend(fontsize=12)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel('TopoI N3E scaled', size=20)
    plot1.set_title(f'{FE_track_pair_name} over the {genes_set_name}s', size=20)     
    plt.savefig(f'{out_path}\Figures\Plots\\{genes_set_name}\\{FE_track_pair_name}_over_{genes_set_name}_{win_width}bp.png', dpi=400, figsize=(10, 6))   
    plt.close()      

    #Make ar from dict.
    gene_US_F_mean=dict_to_ar(gene_US_F_mean_dict)
    gene_DS_F_mean=dict_to_ar(gene_DS_F_mean_dict)
    gene_B_F_mean=dict_to_ar(gene_B_F_mean_dict)
    gene_US_R_mean=dict_to_ar(gene_US_R_mean_dict)
    gene_DS_R_mean=dict_to_ar(gene_DS_R_mean_dict)
    gene_B_R_mean=dict_to_ar(gene_B_R_mean_dict)    

    #Write table contains FE for US, GB, DS of TUs in a set.
    print(f'Writing FE for TUs\' TU, GB, DS...')
    write_genes_FE(gene_US_F_mean_dict, gene_B_F_mean_dict, gene_DS_F_mean_dict, FE_track_pair_name, f'{out_path}\Signal_of_TUs_tab\\{genes_set_name}\{FE_track_pair_name}_over_{genes_set_name}_{win_width}bp_F.txt')
    write_genes_FE(gene_US_R_mean_dict, gene_B_R_mean_dict, gene_DS_R_mean_dict, FE_track_pair_name, f'{out_path}\Signal_of_TUs_tab\\{genes_set_name}\{FE_track_pair_name}_over_{genes_set_name}_{win_width}bp_R.txt')    

    #Plot distribution of mean TUs' FEs.
    print(f'Plotting FE distribution over TU, GB, DS...')
    plot_FE_dist_UDB(gene_US_F_mean, f'{FE_track_pair_name} US', gene_B_F_mean, f'{FE_track_pair_name} TU body', gene_DS_F_mean, f'{FE_track_pair_name} DS', f'{out_path}\Figures\Histograms\\{genes_set_name}\Signal_distribution_{FE_track_pair_name}_over_{genes_set_name}_{win_width}bp_F.png')
    print(len(gene_US_F_mean), len(gene_DS_F_mean), len(gene_B_F_mean))
    plot_FE_dist_UDB(gene_US_R_mean, f'{FE_track_pair_name} US', gene_B_R_mean, f'{FE_track_pair_name} TU body', gene_DS_R_mean, f'{FE_track_pair_name} DS', f'{out_path}\Figures\Histograms\\{genes_set_name}\Signal_distribution_{FE_track_pair_name}_over_{genes_set_name}_{win_width}bp_R.png')
    print(len(gene_US_R_mean), len(gene_DS_R_mean), len(gene_B_R_mean))    
    return gene_US_F, gene_DS_F, gene_B_F, gene_US_R, gene_DS_R, gene_B_R

#######
#Wrapper: reads data and gene annotation, computes signal over US, GB, DS, plots and writes the output.
#######


def Wrapper_signal_over_TUs(dict_of_wigs_path, path_to_annotation, type_of_annot, genes_set_name, deletions_inpath, win_width, length, out_path):
    
    #Reads input data in wig files.
    dict_of_wigs={}
    for pair_name, pair_dict in dict_of_wigs_path.items():
        wig_F_data=wig_parsing(pair_dict['F'])
        wig_R_data=wig_parsing(pair_dict['R'])
        dict_of_wigs[pair_name]={'F' : wig_F_data, 'R' : wig_R_data}
    
    #Reads annotation.
    print(f'Now working with {path_to_annotation}')
    if type_of_annot=='gff':
        gene_annotation=parse_gff_annotation(path_to_annotation, deletions_inpath)[0]['Gene']
    elif type_of_annot=='broadPeak':
        gene_annotation=parse_expression_annotation(path_to_annotation)
    
    #Calculate and plot signal over TUs.
    for FE_track_pair_name, FE_track_pair in dict_of_wigs.items():
        genes_and_FE(gene_annotation, genes_set_name, FE_track_pair, FE_track_pair_name, out_path, deletions_inpath, win_width, length)
    return

#Wrapper_signal_over_TUs(Dict_of_wigs_path_3, Path_to_annotation_2, Type_of_annot_2, Genes_set_name_2, Deletions_inpath, Win_width, Length, Out_path)
#Wrapper_signal_over_TUs(Dict_of_wigs_path_3, Path_to_annotation_5, Type_of_annot_5, Genes_set_name_5, Deletions_inpath, Win_width, Length, Out_path)
#Wrapper_signal_over_TUs(Dict_of_wigs_path_3, Path_to_annotation_6, Type_of_annot_6, Genes_set_name_6, Deletions_inpath, Win_width, Length, Out_path)
#Wrapper_signal_over_TUs(Dict_of_wigs_path_3, Path_to_annotation_7, Type_of_annot_7, Genes_set_name_7, Deletions_inpath, Win_width, Length, Out_path)
#Wrapper_signal_over_TUs(Dict_of_wigs_path_3, Path_to_annotation_8, Type_of_annot_8, Genes_set_name_8, Deletions_inpath, Win_width, Length, Out_path)
#Wrapper_signal_over_TUs(Dict_of_wigs_path_3, Path_to_annotation_9, Type_of_annot_9, Genes_set_name_9, Deletions_inpath, Win_width, Length, Out_path)

Wrapper_signal_over_TUs(Dict_of_wigs_path_2, Path_to_annotation_10, Type_of_annot_10, Genes_set_name_10, Deletions_inpath, Win_width, Length, Out_path)