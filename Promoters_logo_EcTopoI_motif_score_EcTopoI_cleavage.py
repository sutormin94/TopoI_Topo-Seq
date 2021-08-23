###############################################
##Dmitry Sutormin, 2021##
##TopoI ChIP-Seq analysis##

#Script extracts upstream regions of transcription units and alignes them by transcription start site.
#Then frequency of nucleotides is calculated and logo is constructed. Mean position-wise score of EcTopoI binding motif is
#calculated and aligned with the logo.
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
from Bio.Seq import Seq
from Bio import AlignIO
import weblogo
from weblogo import LogoOptions, LogoFormat, eps_formatter, read_seq_data, LogoData, png_formatter, pdf_formatter
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline


#Path to the working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\\"

#Path to genome sequence, fasta.
Genome_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\E_coli_w3110_G_Mu.fasta"

#Path to files with regions to be omitted, bed.
Deletions_inpath='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\Deletions_w3110_G_Mu_SGS.broadPeak'

#Output path.
Outpath=PWD + "Promoter_structure\TopoI_cleavage_over_promoters_st_mock_st_no_Ara\\"


#######
#Checks if directory exists and if not creates.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return

#Output path.
def create_out_dirs(out_path, genes_set_name):
    Dir_check_create(out_path)
    Dir_check_create(out_path+'\Figures\Plots_TSS_TES\\'+genes_set_name)
    Dir_check_create(out_path+'\Signal_of_TSS_TES_tab\\'+genes_set_name)
    Dir_check_create(out_path+'\Signal_of_TSS_TES_wig\\'+genes_set_name)    
    Dir_check_create(out_path+'\Sequences_of_TSS_TES_fasta\\') 
    return


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
#Parses WIG file.
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
#Opens and reads FASTA file with reference genome
#######

def read_genome(genome_path):
    #Opens FASTA file with the reference genome.
    genome=open(genome_path, 'r')
    for record in SeqIO.parse(genome, "fasta"):
        genomefa=str(record.seq)
    return genomefa


#######
#Opens and reads BED file with deletions coordinates.
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
#Write .tab file with FE info for genes US, GB, and DS.
#######

def write_genes_FE(dict1, dict2, FE_track_name, path_out):
    fileout=open(path_out, 'w')
    fileout.write(f'Gene_name\tStart\tEnd\tStrand\t{FE_track_name}_FE_US\t{FE_track_name}_FE_DS\n')
    for k, v in dict1.items():
        fileout.write(f'{k}\t{v[1]}\t{v[2]}\t{v[3]}\t{v[0]}\t{dict2[k][0]}\n')
    fileout.close()
    return


#######
#Write .fasta file with sequences of TUs US and DS regions.
#######

def write_sequences_dict(seq_dict, win_width, length, outpath):
    fileout=open(outpath, 'w')
    for TU_name, sequence_info in seq_dict.items():
        fileout.write(f'>{TU_name}_{sequence_info[1]}_{sequence_info[2]}_{sequence_info[3]}_{win_width}_{length}\n{sequence_info[0]}\n')
    fileout.close()
    return


#######
#Creates logo for US and DS regions.
#######

def Create_logo(alig_inpath, out_path):
    MFA_data=open(alig_inpath)
    MFA_seqs=read_seq_data(MFA_data)
    logodata=LogoData.from_seqs(MFA_seqs)
    logooptions=LogoOptions(yaxis_scale=0.9, pad_right=True, stacks_per_line=90)
    logooptions.show_errorbars=False
    logoformat=LogoFormat(logodata, logooptions)
    pdf=weblogo.logo_formatter.pdf_formatter(logodata, logoformat)
    logout=open(out_path, 'wb')
    logout.write(pdf)
    logout.close()
    return



#######################
## Part 1.
## Classify promoter sequences form RegulonDB by promoter type.
## Scan the secuences with a motif (e.g., EcTopoI motif) using Sarus.
## Plot strand-specific motif score (the result of the scanning procedure) alongside with promoters LOGO.
#######################

#######
#Reads promoters annotation from RegulonDB. Classify promoters by type of a sigma-factor.
#Aligns selected promoters and returns logo. Sends the promoters to SARUS and returns constructs metagene plot on SARUS output.
#######

def Select_promoters_run_sarus_plot(rdb_promoters_inpath, TSS_position, outpath):
    
    #Read regulonDB file.
    promoters_dict={'Strong' : {'Sigma70' : [], 'Sigma54' : [], 'Sigma32' : [], 'Sigma38' : [], 'Sigma28' : [], 'Sigma24' : [], 'mixed' : [], 'unknown' : []}, 'Weak' : {'Sigma70' : [], 'Sigma54' : [], 'Sigma32' : [], 'Sigma38' : [], 'Sigma28' : [], 'Sigma24' : [], 'mixed' : [], 'unknown' : []}}
    promoters_infile=open(rdb_promoters_inpath, 'r')
    for line in promoters_infile:
        line=line.rstrip().split('\t')
        if line[0] not in ['RegulonDB ID']:
            name=line[1]
            sigma=line[4]
            if ', ' in sigma:
                sigma='mixed'
            sequence=line[5]
            confidence=line[7]
                
            promoters_dict[confidence][sigma].append(sequence)
    
    #Write mfa files for SARUS.
    #Run SARUS.
    for conf_name, conf_set in promoters_dict.items():
        for sigma_name, sigma_set in conf_set.items():
            fileout=open(f'{outpath}\{conf_name}_{sigma_name}_promoters.fasta', 'w')
            for i in range(len(sigma_set)):
                fileout.write(f'>{i}_{conf_name}_{sigma_name}\n{sigma_set[i]}\n')
            fileout.close()
            #Tips are taken from https://datatofish.com/command-prompt-python/
            os.system(f'cmd /c "java -jar C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Dists\sarus-2.0.2.jar {outpath}\{conf_name}_{sigma_name}_promoters.fasta C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_identification\EcTopoI_noCTD_Rif_rep123_nm_0.001_peaks_motifs_3.txt -100 > {outpath}\{conf_name}_{sigma_name}_sarus_scan_noCTD_Rif_3.txt"')
            
    
    #Read Sarus output.
    promoter_score_dict={}
    for conf_name, conf_set in promoters_dict.items():
        promoter_score_dict[conf_name]={}
        for sigma_name, sigma_set in conf_set.items():
            sarus_output=f'{outpath}\{conf_name}_{sigma_name}_sarus_scan_noCTD_Rif_3.txt'
            filein=open(sarus_output, 'r')
            
            init=0
            dict_plus_minus={}
            for line in filein:
                if line[0]=='>':
        
                    if init==0:
                        ar_plus=[]
                        ar_minus=[]  
                    elif init==1:
                        seq_len=len(ar_plus)
                        dict_plus_minus[seq_id]=[ar_plus, ar_minus]
                        ar_plus=[]
                        ar_minus=[]                
                        
                    line=line.lstrip('>').rstrip().split(' ')
                    seq_id=line[0] 
                    init=1
                    
                else:
                    line=line.rstrip().split('\t')
                    score=float(line[0])
                    position=int(line[1])
                    strand=line[2]
                    if strand=='+':
                        ar_plus.append(score)
                    elif strand=='-':
                        ar_minus.append(score)
        
            filein.close()
            print(f'{sarus_output} was parsed succesfully')
            
            Num_proms=len(dict_plus_minus)
    
            #Calculate signal over TSS/TES.    
            #Calculate mean, std.
            promoter_F_mean=[]
            promoter_R_mean=[]
            promoter_F_std=[]
            promoter_R_std=[]
            
            for i in range(seq_len):
                promoter_F_ar=[]
                promoter_R_ar=[]
                
                for gene_name, scan_info in dict_plus_minus.items():
                    promoter_F_ar.append(scan_info[0][i])
                    promoter_R_ar.append(scan_info[1][i])
                
                promoter_F_mean.append(np.mean(promoter_F_ar))
                promoter_F_std.append(np.std(promoter_F_ar))
                promoter_R_mean.append(np.mean(promoter_R_ar))
                promoter_R_std.append(np.std(promoter_R_ar))    
                
            #Keep data.
            promoter_score_dict[conf_name][sigma_name]={'F mean' : promoter_F_mean, 'R mean' : promoter_R_mean, 'F std' : promoter_F_std, 'R std' : promoter_R_std, 'Num proms' : Num_proms}
            
            #Compute confidential interval borders (+/- 1*SEM).
            Upper_conf_interval_F=np.array(promoter_F_mean)+(np.array(promoter_F_std)/np.sqrt(Num_proms))
            Lower_conf_interval_F=np.array(promoter_F_mean)-(np.array(promoter_F_std)/np.sqrt(Num_proms))
            Upper_conf_interval_R=np.array(promoter_R_mean)+(np.array(promoter_R_std)/np.sqrt(Num_proms))
            Lower_conf_interval_R=np.array(promoter_R_mean)-(np.array(promoter_R_std)/np.sqrt(Num_proms))      
            
            #Plot Sarus output as metagene plot.
            #Plot signal over anchor. 
            print(f'Plotting signal over anchor...')
            plt.figure(figsize=(13.1, 3), dpi=100)
            plot1=plt.subplot(111)  
            positions=np.arange(-TSS_position+1, seq_len+(-TSS_position+1), 1)  
            
            #F strand.
            plot1.plot(positions,  promoter_F_mean, linestyle='--', color='#6a65c7', linewidth=1.5, alpha=1, label='Motif score: coding strand') 
            plot1.plot(positions,  Upper_conf_interval_F,  linestyle='--', color='#6a65c7', linewidth=0.8, alpha=0.2)
            plot1.plot(positions,  Lower_conf_interval_F,  linestyle='--', color='#6a65c7', linewidth=0.8, alpha=0.2)
            plot1.fill_between(positions, Lower_conf_interval_F, Upper_conf_interval_F, facecolor='#6a65c7', alpha=0.1, interpolate=True)       
            #R strand.
            plot1.plot(positions,  promoter_R_mean, linestyle='--', color='#c96458', linewidth=1.5, alpha=1, label='Motif score: template strand')
            plot1.plot(positions,  Upper_conf_interval_R,  linestyle='--', color='#c96458', linewidth=0.8, alpha=0.2)   
            plot1.plot(positions,  Lower_conf_interval_R,  linestyle='--', color='#c96458', linewidth=0.8, alpha=0.2)
            plot1.fill_between(positions, Lower_conf_interval_R, Upper_conf_interval_R, facecolor='#c96458', alpha=0.1, interpolate=True)             
            #Mark TSS.
            plot1.axvline(0, color='black', linestyle='--', alpha=0.5, linewidth=1, zorder=20)
            
            ticks=np.arange(-TSS_position+1, seq_len+(-TSS_position+1), 10).tolist()
            plot1.set_xticks(ticks)
            ticks_lables=np.arange(-TSS_position+1, seq_len+(-TSS_position+1), 10).tolist()
            ticks_lables[ticks.index(0)]='TSS'
            plot1.set_xticklabels(ticks_lables)
            plot1.legend(fontsize=12, frameon=False)    
            plot1.set_xlabel('Distance, bp', size=20)
            plot1.set_ylabel(f'{conf_name} {sigma_name} ' + '$\it{E. coli}$ promoters score' + f'({Num_proms})', size=10) 
            plt.xlim([-TSS_position, seq_len+(-TSS_position+1)+1])
            plt.tight_layout()
            plt.savefig(f'{outpath}\{conf_name}_{sigma_name}_sarus_scan_noCTD_Rif_3.png', dpi=400, figsize=(13.1, 3))   
            plt.close()         
            
            #Create LOGO for a set of promoters.
            Create_logo(f'{outpath}\{conf_name}_{sigma_name}_promoters.fasta', f'{outpath}\{conf_name}_{sigma_name}_promoters_LOGO.pdf')
    
    return promoter_score_dict

#TSS position in promoters annotation. 1-based.
TSS_position=61

## Takes promoters and classifies them by sigma factor type; Runs SARUS with TopoI motif; Takes SARUS output to create metagene score-plot for promoters of different types separately; Creates LOGO for promoter consensus sequence.
Promoter_score_dict=Select_promoters_run_sarus_plot("C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Promoter_structure\RegulonDB_promoters.txt", TSS_position, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Promoter_structure\Try3_regulonDB_promoters")



#######################
## Part 2.
## Classify promoter sequences form RegulonDB by promoter type.
## Blast the promoter sequences in a reference genome to get real coordinates.
## Extract by the coordinates cleavage signal form wig files (strand-specifically).
## Plot strand-specific cleavage signal over aligned promoters.
#######################

#######
#Reads RegulonDB promoters set and classifies promoters by sigma factor. BLAST promoter sequences in reference genome.
#######

def Select_promoters_blast_sequences(rdb_promoters_inpath, ref_genome_path, outpath):
    
    #Read regulonDB file.
    promoters_dict={'Strong' : {'Sigma70' : [], 'Sigma54' : [], 'Sigma32' : [], 'Sigma38' : [], 'Sigma28' : [], 'Sigma24' : [], 'mixed' : [], 'unknown' : []}, 'Weak' : {'Sigma70' : [], 'Sigma54' : [], 'Sigma32' : [], 'Sigma38' : [], 'Sigma28' : [], 'Sigma24' : [], 'mixed' : [], 'unknown' : []}}
    promoters_infile=open(rdb_promoters_inpath, 'r')
    for line in promoters_infile:
        line=line.rstrip().split('\t')
        if line[0] not in ['RegulonDB ID']:
            name=line[1]
            sigma=line[4]
            if ', ' in sigma:
                sigma='mixed'
            sequence=line[5]
            confidence=line[7]
                
            promoters_dict[confidence][sigma].append(sequence)
            
    #Create BLAST database from the reference genome.
    Make_ref_genome_db=NcbimakeblastdbCommandline(dbtype="nucl", input_file=ref_genome_path)    
    print('Making blast database: ' + str(Make_ref_genome_db))
    Make_ref_genome_db()
    
    #Iterate over sets of promoters, write promoter sequences and BLAST them.
    for conf_name, conf_set in promoters_dict.items():
        for sigma_name, sigma_set in conf_set.items():
            
            #Write mfa files for BLAST.
            fileout=open(f'{outpath}\{conf_name}_{sigma_name}_promoters.fasta', 'w')
            for i in range(len(sigma_set)):
                fileout.write(f'>{i}_{conf_name}_{sigma_name}\n{sigma_set[i]}\n')
            fileout.close()
            
            #Run BLAST for promoters.
            Promoters_blast=NcbiblastnCommandline(query=f'{outpath}\{conf_name}_{sigma_name}_promoters.fasta', db=ref_genome_path, out=f'{outpath}\{conf_name}_{sigma_name}_promoters_blast_result.txt', outfmt=6)  
            print('Run blast of promoters sequences: ' + str(Promoters_blast))
            Promoters_blast()  
            
    return promoters_dict


#######
#Read promoters BLAST output.
#######

def read_blast_output_combine(promoters_dict, outpath):
    
    #Keep coordinates of promoters.
    promoters_coords_dict={}
    
    for conf_name, conf_set in promoters_dict.items():
        promoters_coords_dict[conf_name]={}
        for sigma_name, sigma_set in conf_set.items():  
            promoters_coords_dict[conf_name][sigma_name]={}
    
            #Read promoter sequences BLAST results.
            Pblast_output=open(f'{outpath}\{conf_name}_{sigma_name}_promoters_blast_result.txt', 'r')
            
            Promoter_name_list=[]
            i=0
            for line in Pblast_output:
                line=line.rstrip().split('\t')
                
                if int(line[8]) < int(line[9]):
                    start=int(line[8])-1 #Correction for 0-based enumeration of python (BLAST is 1-based)
                    end=int(line[9])     #Correction for slice method in python (end of the interval is not included)
                    strand="+"
                else:
                    start=int(line[9])-1 #Correction for 0-based enumeration of python (BLAST is 1-based)
                    end=int(line[8])     #Correction for slice method in python (end of the interval is not included)
                    strand="-"
                    
                P_name=line[0]
                if P_name in Promoter_name_list:
                    i+=1
                    if line[2]=='100.000' and line[3]=='81':
                        promoters_coords_dict[conf_name][sigma_name][f'{P_name}_{i}']=[start, end, strand]
                else:
                    i=1
                    if line[2]=='100.000' and line[3]=='81':
                        promoters_coords_dict[conf_name][sigma_name][P_name]=[start, end, strand]                    
                                   
    return promoters_coords_dict


#######
#Convert dictionary to array, discard keys.
#######

def dict_to_ar(dictionary):
    ar=[]
    for k,v in dictionary.items():
        ar.append(v[0]) 
    return ar


#######
#Write .tab file with signal info for promoter regions.
#######

def write_promoter_signal(dict1, FE_track_name, path_out):
    fileout=open(path_out, 'w')
    fileout.write(f'Promoter_name\tStart\tEnd\tStrand\t{FE_track_name}_N3E\n')
    for k, v in dict1.items():
        fileout.write(f'{k}\t{v[1]}\t{v[2]}\t{v[3]}\t{v[0]}\n')
    fileout.close()
    return


#######
#Output path.
#######

def create_out_dirs_N3E(out_path):
    Dir_check_create(out_path)
    Dir_check_create(f'{out_path}\\N3E_of_promoters_tab\\')
    Dir_check_create(f'{out_path}\\N3E_of_promoters_wig\\')    
    return


#######
#Returns TopoI strand-specific cleavage signal for promoters.
#######

def promoters_N3E_ss(promoter_annotation, promoter_set_name, FE_track_pair, FE_track_pair_name, out_path, deletions_inpath):
    
    #Standard length of a promoter region in RegulonDB, bp.
    win_width=81
    
    #Number of promoters.
    Num_proms=len(promoter_annotation)
    
    #Create sub-directories for output files.
    create_out_dirs_N3E(out_path)
    
    #Strand-specific tracks.
    FE_track_F=FE_track_pair['F']
    FE_track_R=FE_track_pair['R']
    
    #Parsing deletions
    deletions=deletions_info(deletions_inpath) 

    #Calculate FE over genes.
    promoter_F=[]
    promoter_R=[]
    promoter_F_mean_dict={}
    promoter_R_mean_dict={}
    
    for promoter_name, promoter_info in promoter_annotation.items():
        delited=0
        
        for deletion in deletions:
            if deletion[1]>=promoter_info[0]>=deletion[0] or deletion[1]>=promoter_info[1]>=deletion[0]:
                delited=1
                
        if delited==0:
            start=promoter_info[0]
            end=promoter_info[1]
            strand=promoter_info[2]
            glen=len(FE_track_F)
            
            if strand=='+': 
                promoter_F.append(FE_track_F[start:end])
                promoter_F_mean_dict[promoter_name]=[np.mean(FE_track_F[start:end]), start, end, strand]
                promoter_R.append(FE_track_R[start:end])
                promoter_R_mean_dict[promoter_name]=[np.mean(FE_track_R[start:end]), start, end, strand]                
                
            elif strand=='-':                                         
                promoter_F.append(FE_track_R[start:end][::-1])                                                  #Take care of strand-specificity.                                     
                promoter_F_mean_dict[promoter_name]=[np.mean(FE_track_R[start:end][::-1]), start, end, strand]                          
                promoter_R.append(FE_track_F[start:end][::-1])                                                  #Take care of strand-specificity.
                promoter_R_mean_dict[promoter_name]=[np.mean(FE_track_F[start:end][::-1]), start, end, strand]
                
    #Calculate mean, std.
    promoter_F_mean=[]
    promoter_R_mean=[]
    promoter_F_std=[]
    promoter_R_std=[]
    
    for i in range(win_width):
        promoter_F_ar=[]
        promoter_R_ar=[]
       
        for j in range(len(promoter_F)):
            promoter_F_ar.append(promoter_F[j][i])
            promoter_R_ar.append(promoter_R[j][i])
        
        promoter_F_mean.append(np.mean(promoter_F_ar))
        promoter_F_std.append(np.std(promoter_F_ar))
        promoter_R_mean.append(np.mean(promoter_R_ar))
        promoter_R_std.append(np.std(promoter_R_ar))    
        
    #Write wig-like file with signal over promoters.
    print(f'Writing signal over promoters...')
    write_wig(promoter_F_mean, f'{out_path}\\N3E_of_promoters_wig\\Mean_signal_{FE_track_pair_name}_N3E_over_{promoter_set_name}_F.wig', f'{promoter_set_name}_{FE_track_pair_name}')
    write_wig(promoter_R_mean, f'{out_path}\\N3E_of_promoters_wig\\Mean_signal_{FE_track_pair_name}_N3E_over_{promoter_set_name}_R.wig', f'{promoter_set_name}_{FE_track_pair_name}')
    write_wig(promoter_F_std,  f'{out_path}\\N3E_of_promoters_wig\\STD_signal_{FE_track_pair_name}_N3E_over_{promoter_set_name}_F.wig', f'{promoter_set_name}_{FE_track_pair_name}')
    write_wig(promoter_R_std,  f'{out_path}\\N3E_of_promoters_wig\\STD_signal_{FE_track_pair_name}_N3E_over_{promoter_set_name}_R.wig', f'{promoter_set_name}_{FE_track_pair_name}')    
        
        
    #Compute confidential interval borders (+/- 1*SEM).
    Upper_conf_interval_F=np.array(promoter_F_mean)+(np.array(promoter_F_std)/np.sqrt(Num_proms))
    Lower_conf_interval_F=np.array(promoter_F_mean)-(np.array(promoter_F_std)/np.sqrt(Num_proms))
    Upper_conf_interval_R=np.array(promoter_R_mean)+(np.array(promoter_R_std)/np.sqrt(Num_proms))
    Lower_conf_interval_R=np.array(promoter_R_mean)-(np.array(promoter_R_std)/np.sqrt(Num_proms))      

    #Plot signal over promoter. 
    print(f'Plotting signal over promoter...')
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)  
    positions=np.arange(1,win_width+1,1) 
    #F strand.
    plot1.plot(positions,  promoter_F_mean, linestyle='-', color='#6a65c7', linewidth=1.5, alpha=1, label='Cleavage (N3E): coding strand') 
    plot1.plot(positions,  Upper_conf_interval_F,  linestyle='-', color='#6a65c7', linewidth=0.8, alpha=0.2)
    plot1.plot(positions,  Lower_conf_interval_F,  linestyle='-', color='#6a65c7', linewidth=0.8, alpha=0.2)
    plot1.fill_between(positions, Lower_conf_interval_F, Upper_conf_interval_F, facecolor='#6a65c7', alpha=0.1, interpolate=True)       
    #R strand.
    plot1.plot(positions,  promoter_R_mean, linestyle='-', color='#c96458', linewidth=1.5, alpha=1, label='Cleavage (N3E): template strand')
    plot1.plot(positions,  Upper_conf_interval_R,  linestyle='-', color='#c96458', linewidth=0.8, alpha=0.2)   
    plot1.plot(positions,  Lower_conf_interval_R,  linestyle='-', color='#c96458', linewidth=0.8, alpha=0.2)
    plot1.fill_between(positions, Lower_conf_interval_R, Upper_conf_interval_R, facecolor='#c96458', alpha=0.1, interpolate=True)     
    
    ticks=np.arange(0,win_width+1,10).tolist()
    ticks+=[1]
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    plt.axvline(x=60)
    plot1.legend(fontsize=12, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'{promoter_set_name} ' + '$\it{E. coli}$ TopoI N3E ' + f'({Num_proms})', size=10) 
    plt.xlim([0, win_width+1])
    plt.savefig(f'{out_path}\Signal_{FE_track_pair_name}_{promoter_set_name}_promoter_TopoI_N3E.png', dpi=400, figsize=(13.1, 3))   
    plt.close()            

    #Write table contains mean signal for promoters.
    print(f'Writing mean signal for promoters...')
    write_promoter_signal(promoter_F_mean_dict, FE_track_pair_name, f'{out_path}\\N3E_of_promoters_tab\\{FE_track_pair_name}_N3E_over_{promoter_set_name}_F.txt')
    write_promoter_signal(promoter_R_mean_dict, FE_track_pair_name, f'{out_path}\\N3E_of_promoters_tab\\{FE_track_pair_name}_N3E_over_{promoter_set_name}_R.txt')        
    
    Promoter_cleavage_signal_ss={'F mean' : promoter_F_mean, 'R mean' : promoter_R_mean, 'F std' : promoter_F_std, 'R std' : promoter_R_std, 'Num proms' : Num_proms}
    return Promoter_cleavage_signal_ss


#######
#Wrapper: classify promoters, extract and BLAST promoter sequences in reference genome.
#Reads BLAST output and collects coordinates of the promoters.
#######

def wrapper_TopoI_cleavage_over_promoters(rdb_promoters_inpath, ref_genome_path, deletions_inpath, dict_of_wigs_path, outpath):
    
    #Read, classify promoters, BLAST sequences in reference genome.
    promoters_dict=Select_promoters_blast_sequences(rdb_promoters_inpath, ref_genome_path, outpath)
    
    #Read BLAST output.
    promoters_coords_dict=read_blast_output_combine(promoters_dict, outpath)
    
    #Reads input data in wig files.
    dict_of_wigs={}
    for pair_name, pair_dict in dict_of_wigs_path.items():
        wig_F_data=wig_parsing(pair_dict['F'])
        wig_R_data=wig_parsing(pair_dict['R'])
        dict_of_wigs[pair_name]={'F' : wig_F_data, 'R' : wig_R_data}  
    
    #Make metagene plot.
    promoter_cleavage_dict={}
    for wig_pair_name, wig_pair_data in dict_of_wigs.items():
        promoter_cleavage_dict[wig_pair_name]={}
        for conf_name, conf_set in promoters_coords_dict.items():
            promoter_cleavage_dict[wig_pair_name][conf_name]={}
            for sigma_name, sigma_set in conf_set.items(): 
                print(f'Now on {conf_name} {sigma_name} promoters applying {wig_pair_name} data...')
                Promoter_cleavage_signal_ss=promoters_N3E_ss(sigma_set, f'{conf_name}_{sigma_name}', wig_pair_data, wig_pair_name, outpath, deletions_inpath)
                promoter_cleavage_dict[wig_pair_name][conf_name][sigma_name]=Promoter_cleavage_signal_ss
    
    return promoter_cleavage_dict

#Path to RegulonDB promoters.
Path_to_RegulonDB_promoters_db="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Promoter_structure\RegulonDB_promoters.txt"

Dict_of_wig_pairs_path={'TopoI_N3E_st_mock_st_no_Ara' :     {'F' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\WIG_NE_strand_specific_masked_scaled_av_masked_accB_subtract_mock_subtract_no_Ara\TopoI_Ara_N3E_F_masked_scaled_av_123_subtr_mock_subtr_no_Ara.wig', 'R' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\WIG_NE_strand_specific_masked_scaled_av_masked_accB_subtract_mock_subtract_no_Ara\TopoI_Ara_N3E_R_masked_scaled_av_123_subtr_mock_subtr_no_Ara.wig'},
                          }

Promoter_cleavage_dict=wrapper_TopoI_cleavage_over_promoters(Path_to_RegulonDB_promoters_db, Genome_path, Deletions_inpath, Dict_of_wig_pairs_path, Outpath)



#######################
## Part 3.
## Plot score and cleavage data together.
#######################

def plot_score_and_cleavage(promoter_score_dict, promoter_cleavage_dict, TSS_position, out_path):
    
    #Iterate over datasets, groups of promoters.
    for cleavage_pair_name, cleavage_pair_data in promoter_cleavage_dict.items():
        for conf_level, conf_data in cleavage_pair_data.items():
            for sigma_name, sigma_data in conf_data.items():
                
                #Extract score and cleavage data.
                Score_data=promoter_score_dict[conf_level][sigma_name]
                Cleavage_data=promoter_cleavage_dict[cleavage_pair_name][conf_level][sigma_name]
                
                #Number of promoters.
                Num_proms=Score_data['Num proms']               
                
                #Compute confidential interval borders (+/- 1*SEM) for score.
                Score_Upper_conf_interval_F=np.array(Score_data['F mean'])+(np.array(Score_data['F std'])/np.sqrt(Num_proms))
                Score_Lower_conf_interval_F=np.array(Score_data['F mean'])-(np.array(Score_data['F std'])/np.sqrt(Num_proms))
                Score_Upper_conf_interval_R=np.array(Score_data['R mean'])+(np.array(Score_data['R std'])/np.sqrt(Num_proms))
                Score_Lower_conf_interval_R=np.array(Score_data['R mean'])-(np.array(Score_data['R std'])/np.sqrt(Num_proms))  
                
                #Length of the scanned part of the sequences.
                seq_len=len(Score_data['F mean'])
            
                print(f'Plotting signal over promoter...')
                plt.figure(figsize=(13.1, 3), dpi=100)
                plot1=plt.subplot(111)  
                #Plot score over promoter. 
                positions=np.arange(-TSS_position+1, seq_len+(-TSS_position+1), 1)  
                #F strand score.
                plot1.plot(positions,  Score_data['F mean'], linestyle='--', color='#6a65c7', linewidth=1, alpha=1, label='Score: coding strand') 
                plot1.plot(positions,  Score_Upper_conf_interval_F,  linestyle='--', color='#6a65c7', linewidth=0.5, alpha=0.2)
                plot1.plot(positions,  Score_Lower_conf_interval_F,  linestyle='--', color='#6a65c7', linewidth=0.5, alpha=0.2)
                plot1.fill_between(positions, Score_Lower_conf_interval_F, Score_Upper_conf_interval_F, facecolor='#6a65c7', alpha=0.1, interpolate=True)       
                #R strand score.
                plot1.plot(positions,  Score_data['R mean'], linestyle='--', color='#c96458', linewidth=1, alpha=1, label='Score: template strand')
                plot1.plot(positions,  Score_Upper_conf_interval_R,  linestyle='--', color='#c96458', linewidth=0.5, alpha=0.2)   
                plot1.plot(positions,  Score_Lower_conf_interval_R,  linestyle='--', color='#c96458', linewidth=0.5, alpha=0.2)
                plot1.fill_between(positions, Score_Lower_conf_interval_R, Score_Upper_conf_interval_R, facecolor='#c96458', alpha=0.1, interpolate=True)     
                #Mark TSS.
                plot1.axvline(0, color='black', linestyle='--', alpha=0.5, linewidth=1, zorder=20)
                
                ticks=np.arange(-TSS_position+1, seq_len+(-TSS_position+1), 10).tolist()
                plot1.set_xticks(ticks)
                ticks_lables=np.arange(-TSS_position+1, seq_len+(-TSS_position+1), 10).tolist()
                ticks_lables[ticks.index(0)]='TSS'
                plot1.set_xticklabels(ticks_lables)
                plot1.legend(fontsize=12, frameon=False, loc='upper left', handlelength=0.5, handletextpad=0.5)    
                plot1.set_xlabel('Distance, bp', size=20)
                plot1.set_ylabel(f'{sigma_name} promoters EcTopoI motif score ({Num_proms})', size=9) 
                plot1.set_xlim([-TSS_position, seq_len+(-TSS_position+1)+1])
                
                
                #Plot cleavage signal over promoter. 
                plot2=plot1.twinx()
                
                #Compute confidential interval borders (+/- 1*SEM) for cleavage.
                Cleavage_Upper_conf_interval_F=np.array(Cleavage_data['F mean'])+(np.array(Cleavage_data['F std'])/np.sqrt(Num_proms))
                Cleavage_Lower_conf_interval_F=np.array(Cleavage_data['F mean'])-(np.array(Cleavage_data['F std'])/np.sqrt(Num_proms))
                Cleavage_Upper_conf_interval_R=np.array(Cleavage_data['R mean'])+(np.array(Cleavage_data['R std'])/np.sqrt(Num_proms))
                Cleavage_Lower_conf_interval_R=np.array(Cleavage_data['R mean'])-(np.array(Cleavage_data['R std'])/np.sqrt(Num_proms))                 
                
                #Make datasets of the same dimensions.
                Cleavage_Upper_conf_interval_F=Cleavage_Upper_conf_interval_F.tolist()[:len(Score_data['F mean'])]
                Cleavage_Lower_conf_interval_F=Cleavage_Lower_conf_interval_F.tolist()[:len(Score_data['F mean'])]
                Cleavage_Upper_conf_interval_R=Cleavage_Upper_conf_interval_R.tolist()[:len(Score_data['F mean'])]
                Cleavage_Lower_conf_interval_R=Cleavage_Lower_conf_interval_R.tolist()[:len(Score_data['F mean'])]   
                Cleavage_data_F_mean=Cleavage_data['F mean'][:len(Score_data['F mean'])]
                Cleavage_data_R_mean=Cleavage_data['R mean'][:len(Score_data['F mean'])]
                
                #F strand cleavage.
                plot2.plot(positions,  Cleavage_data_F_mean, linestyle='-', color='#4b478a', linewidth=2, alpha=1, label='Cleavage: coding strand') 
                plot2.plot(positions,  Cleavage_Upper_conf_interval_F,  linestyle='-', color='#4b478a', linewidth=0.8, alpha=0.2)
                plot2.plot(positions,  Cleavage_Lower_conf_interval_F,  linestyle='-', color='#4b478a', linewidth=0.8, alpha=0.2)
                plot2.fill_between(positions, Cleavage_Lower_conf_interval_F, Cleavage_Upper_conf_interval_F, facecolor='#4b478a', alpha=0.2, interpolate=True)       
                #R strand cleavage.
                plot2.plot(positions,  Cleavage_data_R_mean, linestyle='-', color='#a34f44', linewidth=2, alpha=1, label='Cleavage: template strand')
                plot2.plot(positions,  Cleavage_Upper_conf_interval_R,  linestyle='-', color='#a34f44', linewidth=0.8, alpha=0.2)   
                plot2.plot(positions,  Cleavage_Lower_conf_interval_R,  linestyle='-', color='#a34f44', linewidth=0.8, alpha=0.2)
                plot2.fill_between(positions, Cleavage_Lower_conf_interval_R, Cleavage_Upper_conf_interval_R, facecolor='#a34f44', alpha=0.2, interpolate=True)     
                #Mark TSS.
                plot2.axvline(0, color='black', linestyle='--', alpha=0.5, linewidth=1, zorder=20)
                
                ticks=np.arange(-TSS_position+1, seq_len+(-TSS_position+1), 10).tolist()
                plot2.set_xticks(ticks)
                ticks_lables=np.arange(-TSS_position+1, seq_len+(-TSS_position+1), 10).tolist()
                ticks_lables[ticks.index(0)]='TSS'
                plot2.set_xticklabels(ticks_lables)
                plot2.legend(fontsize=12, frameon=False, loc='lower left', handlelength=0.5, handletextpad=0.5)    
                plot2.set_xlabel('Distance, bp', size=20)
                plot2.set_ylabel(f'{sigma_name} promoters EcTopoI N3E ({Num_proms})', size=9) 
                plot2.set_xlim([-TSS_position, seq_len+(-TSS_position+1)+1])                
                
                plt.tight_layout()
                plt.savefig(f'{out_path}\Signal_{cleavage_pair_name}_{conf_level}_{sigma_name}_promoter_EcTopoI_score_and_N3E.png', dpi=400, figsize=(13.1, 3))
                plt.savefig(f'{out_path}\Signal_{cleavage_pair_name}_{conf_level}_{sigma_name}_promoter_EcTopoI_score_and_N3E.svg', dpi=400, figsize=(13.1, 3))
                plt.close()                            
    
    
    return

plot_score_and_cleavage(Promoter_score_dict, Promoter_cleavage_dict, TSS_position, 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Promoter_structure\TopoI_cleavage_st_mock_st_no_Ara_and_score_over_promoters\\')