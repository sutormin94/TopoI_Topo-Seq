###############################################
##Dmitry Sutormin, 2020##
##DRIP-Seq analysis##

####
#The only purpose - to subtract "+" and "-" strands of DRIP-Seq experiment (or any other strand-specific experiments).
####

###############################################

import numpy as np

#Path to the working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoI-ChIP-Seq\TopA_ChIP-Seq\EcTopoI_G116S_M320V_Topo-Seq\\"
#Path to the files with data to be scaled.
IN_path_dict={'1' :  PWD +        "WIG_masked\DSu_13_S13_edt_N3E_F.wig",
              '2' :  PWD +        "WIG_masked\DSu_13_S13_edt_N3E_R.wig",  
              '3' :  PWD +        "WIG_masked\DSu_14_S14_edt_N3E_F.wig",
              '4' :  PWD +        "WIG_masked\DSu_14_S14_edt_N3E_R.wig", 
              '5' :  PWD +        "WIG_masked\DSu_15_S15_edt_N3E_F.wig",
              '6' :  PWD +        "WIG_masked\DSu_15_S15_edt_N3E_R.wig",  
              '7' :  PWD +        "WIG_masked\DSu_16_S16_edt_N3E_F.wig",
              '8' :  PWD +        "WIG_masked\DSu_16_S16_edt_N3E_R.wig",
              '9' :  PWD +        "WIG_masked\DSu_17_S17_edt_N3E_F.wig",
              '10' : PWD +        "WIG_masked\DSu_17_S17_edt_N3E_R.wig",  
              '11' : PWD +        "WIG_masked\DSu_18_S18_edt_N3E_F.wig",
              '12' : PWD +        "WIG_masked\DSu_18_S18_edt_N3E_R.wig",
              '13' : PWD +        "WIG_masked\DSu_19_S19_edt_N3E_F.wig",
              '14' : PWD +        "WIG_masked\DSu_19_S19_edt_N3E_R.wig",  
              '15' : PWD +        "WIG_masked\DSu_20_S20_edt_N3E_F.wig",
              '16' : PWD +        "WIG_masked\DSu_20_S20_edt_N3E_R.wig",
              '17' : PWD +        "WIG_masked\DSu_21_S21_edt_N3E_F.wig",
              '18' : PWD +        "WIG_masked\DSu_21_S21_edt_N3E_R.wig",  
              '19' : PWD +        "WIG_masked\DSu_22_S22_edt_N3E_F.wig",
              '20' : PWD +        "WIG_masked\DSu_22_S22_edt_N3E_R.wig",  
              '21' : PWD +        "WIG_masked\DSu_23_S23_edt_N3E_F.wig",
              '22' : PWD +        "WIG_masked\DSu_23_S23_edt_N3E_R.wig",  
              '23' : PWD +        "WIG_masked\DSu_24_S24_edt_N3E_F.wig",
              '24' : PWD +        "WIG_masked\DSu_24_S24_edt_N3E_R.wig",               
              }

#Dict with scaling coefficients.
Scaling_dict={'1' :  0.737154157, '2' :  0.737154157, '3' :  0.76483886,  '4' :  0.76483886, 
              '5' :  0.970584457, '6' :  0.970584457, '7' :  1.0,         '8' :  1.0,  
              '9' :  0.674768937, '10' : 0.674768937, '11' : 0.703533009, '12' : 0.703533009, 
              '13' : 0.589865485, '14' : 0.589865485, '15' : 0.816500128, '16' : 0.816500128, 
              '17' : 0.700313237, '18' : 0.700313237, '19' : 0.360237308, '20' : 0.360237308, 
              '21' : 0.505132503, '22' : 0.505132503, '23' : 0.839625826, '24' : 0.839625826} 


#ID or short description of the track (will be the name of a track in IGV).
name_dict={'1' :  "1_TopoI_Ara_Mock_DSu_13_edt_N3E_F_scaled",
           '2' :  "1_TopoI_Ara_Mock_DSu_13_edt_N3E_R_scaled",  
           '3' :  "1_TopoI_Ara_IP_DSu_14_edt_N3E_F_scaled",
           '4' :  "1_TopoI_Ara_IP_DSu_14_edt_N3E_R_scaled", 
           '5' :  "1_TopoI_Mock_DSu_15_edt_N3E_F_scaled",
           '6' :  "1_TopoI_Mock_DSu_15_edt_N3E_R_scaled",  
           '7' :  "1_TopoI_IP_DSu_16_edt_N3E_F_scaled",
           '8' :  "1_TopoI_IP_DSu_16_edt_N3E_R_scaled",
           '9' :  "2_TopoI_Ara_Mock_DSu_17_edt_N3E_F_scaled",
           '10' : "2_TopoI_Ara_Mock_DSu_17_edt_N3E_R_scaled",  
           '11' : "2_TopoI_Ara_IP_DSu_18_edt_N3E_F_scaled",
           '12' : "2_TopoI_Ara_IP_DSu_18_edt_N3E_R_scaled",
           '13' : "2_TopoI_Mock_DSu_19_edt_N3E_F_scaled",
           '14' : "2_TopoI_Mock_DSu_19_edt_N3E_R_scaled",  
           '15' : "2_TopoI_IP_DSu_20_edt_N3E_F_scaled",
           '16' : "2_TopoI_IP_DSu_20_edt_N3E_R_scaled",
           '17' : "3_TopoI_Ara_Mock_DSu_21_edt_N3E_F_scaled",
           '18' : "3_TopoI_Ara_Mock_DSu_21_edt_N3E_R_scaled",  
           '19' : "3_TopoI_Ara_IP_DSu_22_edt_N3E_F_scaled",
           '20' : "3_TopoI_Ara_IP_DSu_22_edt_N3E_R_scaled",  
           '21' : "3_TopoI_Mock_DSu_23_edt_N3E_F_scaled",
           '22' : "3_TopoI_Mock_DSu_23_edt_N3E_R_scaled",  
           '23' : "3_TopoI_IP_DSu_24_edt_N3E_F_scaled",
           '24' : "3_TopoI_IP_DSu_24_edt_N3E_R_scaled",            
           }

#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name_manual=''
#Mode for Chromosome name writing: 0 - auto detection from bed file provided, 1 - manualy provided by user in Chromosome_name variable.
Auto_or_manual=int(0)

#Path to the output directory.
PWD_out=PWD + "WIG_masked_scaled\\"
#Output path to the final file (fold enrichment).
FE_file_path_dict={'1' :  PWD_out + "1_TopoI_Ara_Mock_DSu_13_edt_N3E_F_scaled.wig",
                   '2' :  PWD_out + "1_TopoI_Ara_Mock_DSu_13_edt_N3E_R_scaled.wig",  
                   '3' :  PWD_out + "1_TopoI_Ara_IP_DSu_14_edt_N3E_F_scaled.wig",
                   '4' :  PWD_out + "1_TopoI_Ara_IP_DSu_14_edt_N3E_R_scaled.wig", 
                   '5' :  PWD_out + "1_TopoI_Mock_DSu_15_edt_N3E_F_scaled.wig",
                   '6' :  PWD_out + "1_TopoI_Mock_DSu_15_edt_N3E_R_scaled.wig",  
                   '7' :  PWD_out + "1_TopoI_IP_DSu_16_edt_N3E_F_scaled.wig",
                   '8' :  PWD_out + "1_TopoI_IP_DSu_16_edt_N3E_R_scaled.wig",
                   '9' :  PWD_out + "2_TopoI_Ara_Mock_DSu_17_edt_N3E_F_scaled.wig",
                   '10' : PWD_out + "2_TopoI_Ara_Mock_DSu_17_edt_N3E_R_scaled.wig",  
                   '11' : PWD_out + "2_TopoI_Ara_IP_DSu_18_edt_N3E_F_scaled.wig",
                   '12' : PWD_out + "2_TopoI_Ara_IP_DSu_18_edt_N3E_R_scaled.wig",
                   '13' : PWD_out + "2_TopoI_Mock_DSu_19_edt_N3E_F_scaled.wig",
                   '14' : PWD_out + "2_TopoI_Mock_DSu_19_edt_N3E_R_scaled.wig",  
                   '15' : PWD_out + "2_TopoI_IP_DSu_20_edt_N3E_F_scaled.wig",
                   '16' : PWD_out + "2_TopoI_IP_DSu_20_edt_N3E_R_scaled.wig",
                   '17' : PWD_out + "3_TopoI_Ara_Mock_DSu_21_edt_N3E_F_scaled.wig",
                   '18' : PWD_out + "3_TopoI_Ara_Mock_DSu_21_edt_N3E_R_scaled.wig",  
                   '19' : PWD_out + "3_TopoI_Ara_IP_DSu_22_edt_N3E_F_scaled.wig",
                   '20' : PWD_out + "3_TopoI_Ara_IP_DSu_22_edt_N3E_R_scaled.wig",  
                   '21' : PWD_out + "3_TopoI_Mock_DSu_23_edt_N3E_F_scaled.wig",
                   '22' : PWD_out + "3_TopoI_Mock_DSu_23_edt_N3E_R_scaled.wig",  
                   '23' : PWD_out + "3_TopoI_IP_DSu_24_edt_N3E_F_scaled.wig",
                   '24' : PWD_out + "3_TopoI_IP_DSu_24_edt_N3E_R_scaled.wig",                   
                   }


#######
#Parses WIG file.
#######

def wig_parsing(wigfile):
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    Dict_of_chromosomes_data={}
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0]=='fixedStep':
            chrom_name=line[1].split('=')[1]
            Dict_of_chromosomes_data[chrom_name]=[]
        if line[0] not in ['track', 'fixedStep']:
            Dict_of_chromosomes_data[chrom_name].append(float(line[0]))
    wigin.close()
    
    for Chromosome_name, data in Dict_of_chromosomes_data.items():
        data_array=np.array(data)
        data_mean=np.mean(data_array)
        print(f'Mean coverage of {Chromosome_name}: {data_mean}')
        Dict_of_chromosomes_data[Chromosome_name]=data_array
    return Dict_of_chromosomes_data


def read_files(input_dict):
    Data_dict={}
    for name, path in input_dict.items():
        Data_dict[name]=wig_parsing(path)
        print(f'Progress: {name}/{len(input_dict)}')
    return Data_dict

IN_dict=read_files(IN_path_dict)

def scale_write(IN_dict, Scale_dict, name_dict, Auto_or_manual, Chromosome_name_manual, FE_file_path_dict):
    for sample_name, sample_data in IN_dict.items():
        print(f'Now is processing: {sample_name}')
        print(f'Progress: {sample_name}/{len(IN_dict)}')
        FE_for_out=open(FE_file_path_dict[sample_name]+'.wig', 'w')
        #Write file with scaled data.
        for Chromosome_name, data in sample_data.items():
            print(f'Average covarage before scaling: {np.mean(data)}')
            For_data_scaled=np.array(data)*Scale_dict[sample_name]
            print(f'Average covarage of after scaling: {np.mean(For_data_scaled)}')            
            if Auto_or_manual==0:
                FE_for_out.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
            elif Auto_or_manual==1:
                FE_for_out.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name_manual+' start=1 step=1\n')
            for i in range(len(data)):
                FE_for_out.write(str(For_data_scaled[i])+'\n')
        FE_for_out.close()      
    return

scale_write(IN_dict, Scaling_dict, name_dict, Auto_or_manual, Chromosome_name_manual, FE_file_path_dict)
