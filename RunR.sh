bin/bash
#BSUB -n 16
#BSUB -R "rusage[mem=30GB]"
#BSUB -W 40:30
#BSUB -J Amplicon
#BSUB -o out.%J
#BSUB -e err.%J
#BSUB -R "span[hosts=1]"
#BSUB -q gage

module load conda
conda init
conda activate /usr/local/usrapps/gage/njkorth/myR

#Run from home directory
cd ..
#Rscript ./Scripts/MakeMeta_Quality.R

#look at where to trim before running the next 2 scripts:

#Rscript ./Scripts/ASVPick.R
#Rscript ./Scripts/AssignTax.R

#generate tree using mafft and FastTree
module load mafft
mafft --auto --thread -1 Output/dada2_nochim.fa > dada2_nochim_mafft_msa.fa
fasttree -gtr -nt < dada2_nochim_mafft_msa.fa > dada2_nochim.tree

#Can run the ASVprocessing script here but I prefer to run it in R studio to see figure output in real time.
