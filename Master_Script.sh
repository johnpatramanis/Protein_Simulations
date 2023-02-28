#!/bin/bash
#### Pre install both conda environments using conda -f


eval "$(conda shell.bash hook)"
conda activate SlimTree

snakemake -F -j10

cp Newick_Files/For_PaleoProPhyler/Datasets.txt Dataset_Analysis/Workspace
cp Newick_Files/For_PaleoProPhyler/*_Dataset.fasta Dataset_Analysis/Workspace/1_OG_Dataset/

conda deactivate





eval "$(conda shell.bash hook)"
conda activate Analyser
cd Dataset_Analysis
snakemake -F -j24
conda deactivate

