###############################################
##Dmitry Sutormin, 2021##
##Topo-Seq data analysis##

# Classify TCSs by localization in IGRs or if in TU than by strand orientation which is cleaved. 
###############################################

#######
#Packages to be imported.
#######

import random as rd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.stats import binom


#######
#Import data.
#######

#TCSs input.
PWD_peaks="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\TCS_motifs\\TopoI_Ara_TCSs_called_15.BroadPeak"

#TUs annotation.
TUs_groups_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_RNA-Seq\Expression_data\DY330_transcripts\Representative_transcripts\DY330_RNA-Seq_transcripts_representative_EP_del_cor.txt"


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
#Opens and reads BED or narrowPeak files.
#######

def deletions_info(del_path):
    del_ar=[]
    filein=open(del_path, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        del_ar.append([int(line[1]), int(line[2]), line[3].split('_')[2]])
    filein.close()
    return del_ar

#######
#Localize TCSs.
#######

def loc_TCSs(tcss_inpath, tus_inpath):
    #Read TUs data.
    TUs_annot=parse_expression_annotation(tus_inpath)
    #Read TCSs data.
    TCSs_data=deletions_info(tcss_inpath)
    
    #TCSs classification.
    TCSs_classif={'IGR' : 0, 'Coding' : 0, 'Template' : 0}
    
    for TCS in TCSs_data:
        if_IGR=0
        TUs_ar=[]
        for TU_name, TU_data in TUs_annot.items():
            if (TCS[0]>=TU_data[0]) and (TCS[0]<=TU_data[1]):
                
                TUs_ar.append(TU_data)
                if_IGR=1
                
                TCS_strand=TCS[2]
                if TCS_strand=='F':
                    TCS_strand='+'
                elif TCS_strand=='R':
                    TCS_strand='-'
                else:
                    print(TCS_strand)
                    
                TU_strand=TU_data[2]
                                       
        if if_IGR==0:
            TCSs_classif['IGR']+=1
        
        if len(TUs_ar)==1:
            if (TU_strand=="+") and (TCS_strand=="+"):
                TCSs_classif['Coding']+=1
            elif (TU_strand=="-") and (TCS_strand=="-"):
                TCSs_classif['Coding']+=1
            elif (TU_strand=="+") and (TCS_strand=="-"):
                TCSs_classif['Template']+=1
            elif (TU_strand=="-") and (TCS_strand=="+"):
                TCSs_classif['Template']+=1 
        
        elif len(TUs_ar)>1:
            TU_strand_ar=[]
            for TU_data in TUs_ar:
                TU_strand_ar.append(TU_data[2])
                
            if len(list(set(TU_strand_ar)))>1:
                print(TUs_ar)
                continue
            elif len(list(set(TU_strand_ar)))==1:
                TU_strand=list(set(TU_strand_ar))[0]
                
                if (TU_strand=="+") and (TCS_strand=="+"):
                    TCSs_classif['Coding']+=1
                elif (TU_strand=="-") and (TCS_strand=="-"):
                    TCSs_classif['Coding']+=1
                elif (TU_strand=="+") and (TCS_strand=="-"):
                    TCSs_classif['Template']+=1
                elif (TU_strand=="-") and (TCS_strand=="+"):
                    TCSs_classif['Template']+=1                 
    
    Genome_len=4647454
    Intergenic_length=617293
    Deletions_length=126348

    Genes_ratio=(Genome_len-Intergenic_length)/float(Genome_len)
    Igenes_ratio=1-Genes_ratio
    
    Genes_length_cor=int((Genome_len-Deletions_length)*Genes_ratio)
    Intergenic_length_cor=int((Genome_len-Deletions_length)*Igenes_ratio)
    
    Genes_ratio_cor=float(Genes_length_cor)/(Genes_length_cor+Intergenic_length_cor)
    Igenes_ratio_cor=float(Intergenic_length_cor)/(Genes_length_cor+Intergenic_length_cor)
    
    print(f'Binom test p-value for the number of TCSs in intergenic regions: {1-binom.cdf(TCSs_classif["IGR"], len(TCSs_data), Igenes_ratio_cor)}')
    print(f'Enrichment of TCSs in intergenic regions: {TCSs_classif["IGR"]/float(len(TCSs_data)*Igenes_ratio_cor)}')
    
    print(f'Binom test p-value for the number of TCSs on coding strand: {binom.cdf(TCSs_classif["Coding"], len(TCSs_data), Genes_ratio_cor/2)}')
    print(f'Enrichment of TCSs on coding strand: {TCSs_classif["Coding"]/float(len(TCSs_data)*Genes_ratio_cor/2)}')    
    
    print(f'Binom test p-value for the number of TCSs on template strand: {binom.cdf(TCSs_classif["Template"], len(TCSs_data), Genes_ratio_cor/2)}')
    print(f'Enrichment of TCSs on template strand: {TCSs_classif["Template"]/float(len(TCSs_data)*Genes_ratio_cor/2)}')     
        
    print(TCSs_classif)
    
    return

loc_TCSs(PWD_peaks, TUs_groups_path)