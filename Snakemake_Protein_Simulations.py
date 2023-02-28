###### CONDA ENVIRONMENT for Slim Tree: conda install -c conda-forge -c bioconda biopython slim==3.6 matplotlib-base pandas numpy java-jdk c-compiler snakemake biopython
###### Requires GitClone SlimTree


import os
import os.path
from os import listdir
from os.path import isfile, join, isdir
from Bio import SeqIO
################################################################################################################################################################################################################################################################################################################################################
################################################################### USEFULL FUNCTIONS



################################################################################################################ SET UP ################################################################################################################################################################################################################################

SAMPLES=[]
SAMPLES_TO_NEWICK={}
SAMPLES_TO_NUMBER_OF_PROTEINS={}

NEWICK_LIST_FILE=open('Newick_Trees.txt','r')




for line in NEWICK_LIST_FILE:
    line=line.strip().split()
    SAMPLE=line[0]
    SAMPLES.append(SAMPLE)
    SAMPLES_TO_NEWICK[SAMPLE]=line[1]
    SAMPLES_TO_NUMBER_OF_PROTEINS[SAMPLE]=int(line[2])

SAMPLES=list(set(SAMPLES))









############################################################################################################################################################################################################################################################



###### One Rule to Rule them ALL

rule all:
    input:
        'Newick_Trees.txt',
        Slim_Tree_Folder=directory('slim-tree'),
        Newick_Files=expand('Newick_Files/{SAMPLE}/Newick_File.txt',SAMPLE=SAMPLES),
        Slim_Tree_Output_Files=expand('Newick_Files/{SAMPLE}/Simulations_Complete',SAMPLE=SAMPLES),
        PaleoProPhyler_Input_File=expand('Newick_Files/{SAMPLE}/PaleoProPhyler_Input/{SAMPLE}_Dataset.fasta',SAMPLE=SAMPLES),
        Datasets_File='Newick_Files/For_PaleoProPhyler/Datasets.txt',
        Copied_Dataset=expand('Newick_Files/For_PaleoProPhyler/{SAMPLE}_Dataset.fasta',SAMPLE=SAMPLES)













######################################################################################################################
#### Step 0 - Make sure SlimTree exists in folder


#### Step 0
rule Git_Clone_SlimTree_And_Set_Up:
    input:
        'Newick_Trees.txt' 
    output:
        Slim_Tree_Folder=directory('slim-tree'),
        BackUpFiles=directory('backupFiles'),
    threads:1
    run:
        shell('rm -rf slim-tree') # remove if pre-existing repo
        shell(F"git clone -b dev https://github.com/dekoning-lab/slim-tree.git")
        if os.path.isdir('backupFiles')==False:#check if backupFiles folder exists, otherwise bug will occur!
            shell(F"mkdir backupFiles")
    











#######################################################################################################################
##### Step 1 - Create a txt file for each newick



##### Step 1
rule Create_Newick_TXT_Files:
    input:
        'Newick_Trees.txt',
        Slim_Tree_Folder=directory('slim-tree'),
        BackUpFiles=directory('backupFiles'),
    output:
        Newick_Files=expand('Newick_Files/{sample}/Newick_File.txt',sample=SAMPLES)
    threads:1
    run:
        if os.path.isdir('Newick_Files/')==False:
            shell('mkdir Newick_Files')
            
            
        for SAMPLE in SAMPLES: 
        
            if os.path.isdir(F'Newick_Files/{SAMPLE}')==False:  #### Create sample folder
                shell(F'mkdir Newick_Files/{SAMPLE}')
            if os.path.isdir(F'Newick_Files/{SAMPLE}')==True:  #### If it exist, remove previous fodler, create it anew
                shell(F'rm -rf Newick_Files/{SAMPLE}')
                shell(F'mkdir Newick_Files/{SAMPLE}')
            
            
            shell(F"echo '{SAMPLES_TO_NEWICK[SAMPLE]}' > Newick_Files/{SAMPLE}/Newick_File.txt") ### Create Newick txt inside the sample folder
















######################################################################################################################
#### Step 2 - Generate the Trees ###



#### Step 2
rule Generate_The_Trees_Through_SlimTree:
    input:
        Newick_File='Newick_Files/{sample}/Newick_File.txt'
    output:
        Slim_Tree_Output_Files='Newick_Files/{sample}/Simulations_Complete'
    threads:4
    run:
    
        SAMPLE=wildcards.sample
        
        ### For each simulated gene
        for GENE in range(0,SAMPLES_TO_NUMBER_OF_PROTEINS[SAMPLE]): ### Loop for each gene for this phylogeny
            
            if os.path.isdir(F'Newick_Files/{SAMPLE}/GENE_{GENE}')==True:
                shell(F'rm -rf Newick_Files/{SAMPLE}/GENE_{GENE}')
                
            shell(F'mkdir Newick_Files/{SAMPLE}/GENE_{GENE}')
            shell(F'python3 slim-tree/SLiMTree.py -g 1000 -i {input.Newick_File} -T SLim-Tree -k 1')
            
            shell(F'mv  Newick_Files/{SAMPLE}/*fasta Newick_Files/{SAMPLE}/GENE_{GENE}/')
        
        shell(F'touch Newick_Files/{SAMPLE}/Simulations_Complete')

















######################################################################################################################
#### Step 3 - Prepare Slim Tree Output for PaleoProPhyler



#### Step 3
rule Format_Data_For_PaleoProPhyler:
    input:
        'Newick_Files/{sample}/Simulations_Complete'
    output:
        'Newick_Files/{sample}/PaleoProPhyler_Input/{sample}_Dataset.fasta'
    threads:1
    run:
        SAMPLE=wildcards.sample
        
        #### Set up folder
        if os.path.isdir(F'Newick_Files/{SAMPLE}/PaleoProPhyler_Input/')==True:
            shell(F'rm -rf Newick_Files/{SAMPLE}/PaleoProPhyler_Input/')
        shell(F'mkdir Newick_Files/{SAMPLE}/PaleoProPhyler_Input/')  


        ### Cycle through Simulated Genes of Phylogeny
        for GENE in range(0,SAMPLES_TO_NUMBER_OF_PROTEINS[SAMPLE]):
        
            ### Read fasta through Biopython
            FASTA_FILE=SeqIO.parse(open(F'Newick_Files/{SAMPLE}/GENE_{GENE}/Newick_File_aa.fasta'),'fasta')
            OUTPUT_FASTA=open(F'Newick_Files/{SAMPLE}/PaleoProPhyler_Input/{SAMPLE}_Dataset.fasta','a')
            
            
            #### Go through Generated fasta file, rename things and transfer to new file
            counter=0
            for fasta in FASTA_FILE:
                name, sequence = fasta.id, str(fasta.seq)
                name = F'>sample{counter}_GENE{GENE}'
                counter+=1
                
                OUTPUT_FASTA.write(name+'\n')
                OUTPUT_FASTA.write(sequence+'\n')
                
                
                
             











######################################################################################################################
#### Step 4 - Final Output Prep form all Runs



#### Step 4
rule Create_Datasets:
    input:
        expand('Newick_Files/{SAMPLE}/PaleoProPhyler_Input/{SAMPLE}_Dataset.fasta',SAMPLE=SAMPLES)
    output:
        'Newick_Files/For_PaleoProPhyler/Datasets.txt'
    threads:1
    run:
        ### Prep Final Directory
        if os.path.isdir(F'Newick_Files/For_PaleoProPhyler/')==True:
            shell(F'rm -rf Newick_Files/For_PaleoProPhyler/')
        shell(F'mkdir Newick_Files/For_PaleoProPhyler/')
        
        ### Create Datasets.txt file with headers for PPP
        shell('''echo "Dataset Ancient_Samples" >  Newick_Files/For_PaleoProPhyler/Datasets.txt''')




######################################################################################################################
#### Step 5 - Final Output Prep form all Runs



#### Step 5
rule Final_Data_Prep:
    input:
        Datasets_File='Newick_Files/For_PaleoProPhyler/Datasets.txt',
        Dataset='Newick_Files/{sample}/PaleoProPhyler_Input/{sample}_Dataset.fasta'
    output:
        Copied_Dataset='Newick_Files/For_PaleoProPhyler/{sample}_Dataset.fasta'
    threads:1
    run:
    
        #### For each generated dataset copy the fasta file over and create an entry in the Dataset.txt file
        
        shell(F'cp {input.Dataset} {output.Copied_Dataset}')
        shell(F'''echo "{wildcards.sample}_Dataset.fasta sample0" >>  {input.Datasets_File}''')




