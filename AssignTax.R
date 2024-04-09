library("tidyverse")
library("dada2")
library("gridExtra")
no_of_cores=16

metadata_df<-read.csv("AmpliconTest.csv")

load("./intermediate/dada2_lean.rda")

ls()

#the silva database is availiable here: zenodo.org/records/4587955

silva<-"RawData/silva_nr99_v138.1_wSpecies_train_set.fa.gz"
taxa <- assignTaxonomy(seqtab.nochim, silva, multithread=no_of_cores, verbose=TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Name the sequence file
sequence_outfile <- "dada2_nochim.fa"
print(paste0('Writing FASTA file ', sequence_outfile, ' with ', ncol(seqtab.nochim), ' sequences.'))

file.remove(sequence_outfile) # Make sure the file does not exist

for (i in 1:ncol(seqtab.nochim)){
    cat(paste0(">ASV_", i), file = sequence_outfile, sep="\n", append = TRUE)
    cat(colnames(seqtab.nochim)[i], file = sequence_outfile, sep="\n", append= TRUE)
}

system("head dada2_nochim.fa", intern=TRUE) # use system() to run bash command in R

#For now skip multiple sequence alignment and tree construction:
# We need to make the external packages MAFFT and FastTree available for our use;
# this is what the next three lines do.
#source(file.path(Sys.getenv("LMOD_PKG"), "init/R"))
#module("load", "mafft")
#module("load", "fasttree")

# Multiple sequence alignment with MAFFT
#system("mafft --auto --thread -1 dada2_nochim.fa > dada2_nochim_mafft_msa.fa", intern=TRUE)

# Phylogenetic tree reconstruction with FastTree
#system("fasttree -gtr -nt < dada2_nochim_mafft_msa.fa > dada2_nochim.tree", intern=TRUE)

list.files(pattern = "*dada2*")

colnames(seqtab.nochim) <- paste0("ASV_", seq(1:ncol(seqtab.nochim)))
rownames(taxa) <- paste0("ASV_", seq(1:nrow(taxa)))

seqtab.nochim[1:3,1:3]
taxa[1:3,]
metadata_df[1:3,]

save(metadata_df, seqtab.nochim,taxa, file = "./intermediate/dada2_seqtab_taxa.rda")

