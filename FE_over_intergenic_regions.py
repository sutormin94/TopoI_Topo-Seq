###############################################
##Dmitry Sutormin, 2021##
##TopoI ChIP-Seq analysis##

#Script computes Fold Enrichment (FE) over intergenic regions of transcription units (IGRs).

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
import pandas as pd


#Path to the directory with features files.
Path_to_input_files='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\\'
#Path to intergenic region file.
IGR_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\E_coli_TF_sites\\"
#Path to the intergenic regions annotation, type of annotation and name of IGR set.
##1##
Path_to_annotation_1=IGR_path + 'DY330_intergenic_regions_filtrated_50_1000bp_no_dps.txt'
Type_of_annot_1='broadPeak'             
Genes_set_name_1='IGR_50_1000_all_no_dps'    

#Path to the file with regions to be omitted (e.g. deletions).
Deletions_inpath='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\\Deletions_w3110_G_Mu_SGS.broadPeak'
#Width of US, DS regions.
Win_width=200
#Length of GB.
Length=100

#Dictionary of pathes to input data.
Dict_of_wigs_path={'RpoC_Borukhov' : Path_to_input_files+'Borukhov_RpoC_Pol_Sofi_LB_FE.wig',
                   'RpoB_Kahramanoglou' : Path_to_input_files+'Kahramanoglou_RpoB_IP_ME.wig',
                   'RpoS_Peano' : Path_to_input_files+'Peano_RpoS_FE_av.wig',
                   'RpoS_Seo' : Path_to_input_files+'Seo_RpoS_FE_Rep1.wig',
                   'RpoD_Myers' : Path_to_input_files + 'Myers_RpoD_FE_av.wig',
                   'HNS_Kahramanoglou' : Path_to_input_files+'Kahramanoglou_HNS_IP_ME.wig',
                   'Fis_Kahramanoglou' : Path_to_input_files+'Kahramanoglou_Fis_IP_ME.wig',
                   'MukB_Nolivos' : Path_to_input_files+'Nolivos_MukB_IP_av.wig',
                   'MatP_Nolivos' : Path_to_input_files+'Nolivos_MatP_IP_av.wig',
                   'TopA_CTD_minus_Rif_minus_av_1_2_3' : Path_to_input_files + 'Sutormin_TopA_ChIP_CTD_minus_Rif_minus_FE_av_123.wig',
                   'TopA_CTD_minus_Rif_plus_av_1_2_3' : Path_to_input_files + 'Sutormin_TopA_ChIP_CTD_minus_Rif_plus_FE_av_123.wig',
                   'TopA_CTD_plus_Rif_minus_av_1_2_3' : Path_to_input_files + 'Sutormin_TopA_ChIP_CTD_plus_Rif_minus_FE_av.wig',
                   'TopA_CTD_plus_Rif_plus_av_2_3' :   Path_to_input_files + 'Sutormin_TopA_ChIP_CTD_plus_Rif_plus_FE_av.wig',                   
                   'Gyrase_Cfx' : Path_to_input_files+'Sutormin_Gyrase_Cfx_10mkM_FE_av.wig',
                   'Gyrase_Cfx_Rif' : Path_to_input_files+'Sutormin_Gyrase_RifCfx_122mkM_10mkM_FE_av.wig',
                   'TopoIV_Cfx' : Path_to_input_files+'Sutormin_TopoIV_Cfx_FE_av.wig',
                   'GC' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_Gyrase_Topo-Seq\Scripts\Gyrase_Topo-seq\Additional_genome_features\E_coli_w3110_Mu_GC_133bp.wig',
                   'RNA_Seq' : Path_to_input_files+'Sutormin_RNA_Seq_Exponential_av.wig',
                   'EcTopoI_motif_plus' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_scanning\SARUS_new\EcTopoI_noCTD_Rif_motif_1_w3110_Mu_scanned_plus.wig',
                   'EcTopoI_motif_minus' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_scanning\SARUS_new\EcTopoI_noCTD_Rif_motif_1_w3110_Mu_scanned_minus.wig',
                   'EcTopoI_motif_both' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_scanning\SARUS_new\EcTopoI_noCTD_Rif_motif_1_w3110_Mu_scanned_both.wig'
                   }
                   

Dict_of_wigs_path_1={'GC' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_Gyrase_Topo-Seq\Scripts\Gyrase_Topo-seq\Additional_genome_features\E_coli_w3110_Mu_GC_133bp.wig',
                     'TopA_CTD_minus_Rif_minus_av_1_2_3' : Path_to_input_files + 'Sutormin_TopA_ChIP_CTD_minus_Rif_minus_FE_av_123.wig',}

Dict_of_wigs_path_2={'TopA_CTD_minus_Rif_minus_av_3_4_6' : Path_to_input_files + 'Sutormin_TopA_ChIP_CTD_minus_Rif_minus_FE_av_346.wig',}

Dict_of_wigs_path_3={'EcTopoI_Ara_N3E_F' : Path_to_input_files + 'TopoI_Ara_N3E_F_masked_scaled_av_123_mock_subtr.wig',
                     'EcTopoI_Ara_N3E_R' : Path_to_input_files + 'TopoI_Ara_N3E_R_masked_scaled_av_123_mock_subtr.wig',
                     'EcTopoI_N3E_F' :     Path_to_input_files + 'TopoI_N3E_F_masked_scaled_av_123_mock_subtr.wig',
                     'EcTopoI_N3E_R' :     Path_to_input_files + 'TopoI_N3E_R_masked_scaled_av_123_mock_subtr.wig',}

Dict_of_wigs_path_4={'EcTopoI_Ara_N3E_F' : Path_to_input_files + 'TopoI_Ara_N3E_F_masked_scaled_av_123_mock_subtr.wig',
                     'EcTopoI_Ara_N3E_R' : Path_to_input_files + 'TopoI_Ara_N3E_R_masked_scaled_av_123_mock_subtr.wig',}

Dict_of_wigs_path_5={'TopoI_Ara_N3E_subtr_mock_subtr_no_Ara_F'  : Path_to_input_files + 'WIG_NE_strand_specific_masked_scaled_av_masked_accB_subtract_mock_subtract_no_Ara\TopoI_Ara_N3E_F_masked_scaled_av_123_subtr_mock_subtr_no_Ara.wig', 
                     'TopoI_Ara_N3E_subtr_mock_subtr_no_Ara_R'  : Path_to_input_files + 'WIG_NE_strand_specific_masked_scaled_av_masked_accB_subtract_mock_subtract_no_Ara\TopoI_Ara_N3E_R_masked_scaled_av_123_subtr_mock_subtr_no_Ara.wig',
                     'TopoI_Ara_N3E_subtr_mock_subtr_no_Ara_FR' : Path_to_input_files + 'WIG_NE_strand_specific_masked_scaled_av_masked_accB_subtract_mock_subtract_no_Ara\TopoI_Ara_N3E_FR_masked_scaled_av_123_subtr_mock_subtr_no_Ara.wig',                     
                     }


#######
#Checks if directory exists and if not creates.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return

#Path to the output directory.
Out_path='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\Metagene_analysis\IGRs_Ara_st\\'

#Output path.
def create_out_dirs(out_path, genes_set_name):
    Dir_check_create(out_path)
    Dir_check_create(out_path+'\Figures\Plots\\'+genes_set_name)
    Dir_check_create(out_path+'\Figures\Histograms\\'+genes_set_name)
    Dir_check_create(out_path+'\Signal_of_TUs_tab\\'+genes_set_name)
    Dir_check_create(out_path+'\Signal_of_TUs_wig\\'+genes_set_name)    
    return

create_out_dirs(Out_path, Genes_set_name_1)


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
    IGR_type_dict={'--': {}, '++' : {}, '-+' : {}, '+-' : {}}
    for line in filein:
        line=line.rstrip().split('\t')
        if line[0] not in ['GeneID', 'OperonID', 'TU_ID']:
            IGR_name=line[0]+'_'+line[6]
            IGR_start=int(line[2])
            IGR_end=int(line[7])
            IGR_type=line[3]+line[9]
            IGR_type_dict[IGR_type][IGR_name]=[IGR_start, IGR_end, IGR_type]
    filein.close()            
    return IGR_type_dict


#######
#Reads annotation of particular set of intergenic regions in xlsx format, select sets by flanking genes orientation.
#######

def read_classify_IGR(IGR_inpath):
    IGR_filein=pd.read_excel(IGR_inpath, index_col=0, header=0)
    IGR_type_dict={'--': {}, '++' : {}, '-+' : {}, '+-' : {}}
    for IGR_type in IGR_type_dict.keys():
        IGR_dict={}
        IGR_set=IGR_filein[(IGR_filein['G1_strand']==IGR_type[0]) & (IGR_filein['G2_strand']==IGR_type[1])][['G1_name', 'G2_name', 'G1_end', 'G2_start', 'G1_strand', 'G2_strand']]
        for index, row in IGR_set.iterrows():
            IGR_dict[row['G1_name']+'_'+row['G2_name']]=[row['G1_end'], row['G2_start'], IGR_type]
            
        IGR_type_dict[IGR_type]=IGR_dict
    
    return IGR_type_dict


#######
#Combine intergenic annotation irrespectively to flanking genes orientation.
#######

def combine_intergenic_regions(IGR_type_dict):
    IGR_type_dict_update=IGR_type_dict
    
    IGR_all_set={}
    for IGR_type, IGR_sets in IGR_type_dict.items():
        IGR_all_set={**IGR_all_set, **IGR_sets}         #A method taken from https://stackoverflow.com/questions/38987/how-do-i-merge-two-dictionaries-in-a-single-expression-in-python-taking-union-o
    
    IGR_type_dict_update['all']=IGR_all_set
    
    return IGR_type_dict_update


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
    if len(ar)>length: #array should be shrinked.
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
    elif len(ar)<length: #array should be stretched out.
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
    fileout.write(f'IGR_name\tStart\tEnd\tStrand\t{FE_track_name}_FE_US\t{FE_track_name}_FE_IGR\t{FE_track_name}_FE_DS\n')
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
#Makes histogram for FE over TUs: US, GB, DS.
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

def genes_and_FE(IGR_type_dict, genes_set_name, FE_track, FE_track_name, out_path, deletions_inpath, win_width, length):
    
    #Parsing deletions
    deletions=deletions_info(deletions_inpath) 
    
    for IGR_type, IGR_dict in IGR_type_dict.items():
        print(f'Now are processed {IGR_type} IGRs')

        #Calculate FE over genes.
        gene_US=np.array([0.0]*win_width)
        gene_DS=np.array([0.0]*win_width)
        IGR=[]
        gene_US_mean_dict={}
        gene_DS_mean_dict={}
        IGR_mean_dict={}
        
        for gene_name, gene_info in IGR_dict.items():
            delited=0
            for deletion in deletions:
                if deletion[1]>=gene_info[0]>=deletion[0] or deletion[1]>=gene_info[1]>=deletion[0]:
                    delited=1
            if delited==0:
                start=gene_info[0]
                end=gene_info[1]
                glen=len(FE_track)
                
                if start<win_width:
                    gene_US+=np.array(FE_track[glen-(win_width-start):] + FE_track[:start])
                    gene_US_mean_dict[gene_name]=[np.mean(FE_track[glen-(win_width-start):] + FE_track[:start]), start, end, gene_info[2]]
                else:
                    gene_US+=np.array(FE_track[start-win_width:start])
                    gene_US_mean_dict[gene_name]=[np.mean(FE_track[start-win_width:start]), start, end, gene_info[2]]
                if end+win_width>glen:
                    gene_DS+=np.array(FE_track[end:] + FE_track[:end+win_width-glen])
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[end:] + FE_track[:end+win_width-glen]), start, end, gene_info[2]]
                else:
                    gene_DS+=np.array(FE_track[end:end+win_width])
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[end:end+win_width]), start, end, gene_info[2]]
                if start<end:
                    IGR.append(FE_track[start:end])
                    IGR_mean_dict[gene_name]=[np.mean(FE_track[start:end]), start, end, gene_info[2]]
                else:
                    IGR.append(FE_track[start:]+FE_track[:end])
                    IGR_mean_dict[gene_name]=[np.mean(FE_track[start:]+FE_track[:end]), start, end, gene_info[2]]                    
    
        Num_genes=len(IGR_dict)
        gene_US=gene_US/Num_genes
        gene_DS=gene_DS/Num_genes
        print(f'FE over IGRs computed...')
    
        #Scale GB length.
        print(f'IGR scaling in progress, it takes some time...')
        gene_body=np.array([0.0]*length)
        for gene in IGR:
            gene_scaled=np.array(scale_gene_body(gene, length))
            print(len(gene_body), len(gene_scaled), len(gene))
            gene_body+=gene_scaled
        IGR=gene_body/Num_genes
    
        #Write wig-like file with FE over US, GB, DS.
        print(f'Writing FE over IGR and US, DS genes...')
        write_wig(np.concatenate((gene_US, IGR, gene_DS), axis=None), f'{out_path}\Signal_of_TUs_wig\\{genes_set_name}\Signal_{FE_track_name}_over_{IGR_type}{genes_set_name}_width_{win_width}bp_gb_{length}bp.wig', f'{win_width}_{length}')
    
        #Plot FE over US, GB, DS. 
        print(f'Plotting FE over IGR and US, DS genes...')
        plt.figure(figsize=(10, 6), dpi=100)
        plot1=plt.subplot(111)  
        positions_DS=np.arange(1+length,win_width+1+length,1)
        positions_US=np.arange(-1-win_width,-1,1)  
        positions_Body=np.arange(0,length,1)
        plot1.plot(positions_US, gene_US, linestyle='-', color='#D17E7E', linewidth=1, label='Rep12') 
        plot1.plot(positions_DS, gene_DS, linestyle='-', color='#D17E7E', linewidth=1)
        plot1.plot(positions_Body, IGR, linestyle='-', color='#D17E7E', linewidth=2)
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
        plot1.set_ylabel('TopoA fold enrichment', size=20)
        plot1.set_title(f'{FE_track_name} over the {IGR_type}{genes_set_name}s', size=20)     
        plt.savefig(f'{out_path}\Figures\Plots\\{genes_set_name}\\{FE_track_name}_over_{IGR_type}{genes_set_name}_{win_width}bp.png', dpi=400, figsize=(10, 6))   
        plt.close()      
    
        #Make ar from dict.
        gene_US_mean=dict_to_ar(gene_US_mean_dict)
        gene_DS_mean=dict_to_ar(gene_DS_mean_dict)
        IGR_mean=dict_to_ar(IGR_mean_dict)
    
        #Write table contains FE for US, GB, DS of TUs in a set.
        print(f'Writing FE for TUs\' IGR and US, DS genes...')
        write_genes_FE(gene_US_mean_dict, IGR_mean_dict, gene_DS_mean_dict, FE_track_name, f'{out_path}\Signal_of_TUs_tab\\{genes_set_name}\{FE_track_name}_over_{IGR_type}{genes_set_name}_{win_width}bp.txt')
    
        #Plot distribution of mean TUs' FEs.
        print(f'Plotting FE distribution over IGR and US, DS genes...')
        plot_FE_dist_UDB(gene_US_mean, f'{IGR_type}{FE_track_name} US', IGR_mean, f'{IGR_type}{FE_track_name} IGR', gene_DS_mean, f'{IGR_type}{FE_track_name} DS', f'{out_path}\Figures\Histograms\\{genes_set_name}\Signal_distribution_{FE_track_name}_over_{IGR_type}{genes_set_name}_{win_width}bp.png')
        print(len(gene_US_mean), len(gene_DS_mean), len(IGR_mean))
    return


#######
#Wrapper: reads data and gene annotation, computes signal over US, GB, DS, plots and writes the output.
#######


def Wrapper_signal_over_TUs(dict_of_wigs_path, path_to_annotation, type_of_annot, genes_set_name, deletions_inpath, win_width, length, out_path):
    #Reads input data in wig files.
    dict_of_wigs={}
    for name, data_path in dict_of_wigs_path.items():
        dict_of_wigs[name]=wig_parsing(data_path)
    
    #Reads annotation.
    print(f'Now working with {path_to_annotation}')
    if type_of_annot=='gff':
        gene_annotation=parse_gff_annotation(path_to_annotation, deletions_inpath)[0]['Gene']
    elif type_of_annot=='broadPeak':
        IGR_type_dict=parse_expression_annotation(path_to_annotation)
    elif type_of_annot=='excel':
        IGR_type_dict=read_classify_IGR(path_to_annotation)
        
    #Combine all intergenic regions.
    IGR_type_dict=combine_intergenic_regions(IGR_type_dict)
    
    #Calculate and plot signal over TUs.
    for FE_track_name, FE_track in dict_of_wigs.items():
        genes_and_FE(IGR_type_dict, genes_set_name, FE_track, FE_track_name, out_path, deletions_inpath, win_width, length)
    return

Wrapper_signal_over_TUs(Dict_of_wigs_path_5, Path_to_annotation_1, Type_of_annot_1, Genes_set_name_1, Deletions_inpath, Win_width, Length, Out_path)
