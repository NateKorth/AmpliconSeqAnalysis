library("tidyverse")
library("dada2")
library("gridExtra")
no_of_cores=16

#this file should contain the names of forward and reverse reads:
metadata_df<-read.csv("AmpliconTest.csv")

#forward reads
fnFs<-paste("./RawData/",metadata_df$R1,sep="")
#reverse reads
fnRs<-paste("./RawData/",metadata_df$R2,sep="")

# We will use a separate folder called "Filtered" in our working directory
# directory where the filtered fastq files will be stored
sample.names<-metadata_df$ID
filt_path <- paste0(getwd(), '/Filtered') # don't change this
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))

any(duplicated(c(filtFs, filtRs)))

FORWARD_TRUNC <- 270 # determine from quality plots
REVERSE_TRUNC <- 235 # determine from quality plots

#see dada 2 documentation, can trim based on quality or other  
#In this case, remove primers with trimLeft call:forward primer: 17bases, reverse primer: 21 bases

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(FORWARD_TRUNC,REVERSE_TRUNC),
                     trimLeft=c(17, 21), maxEE=c(5,8),
                     multithread=no_of_cores,
                     matchIDs=TRUE, compress=TRUE,
                     verbose=TRUE,n=1e6)

derepFs <- derepFastq(filtFs, n = 5e06, verbose = TRUE)
derepRs <- derepFastq(filtRs, n = 5e06, verbose = TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

errF <- learnErrors(filtFs, verbose=TRUE, multithread=no_of_cores,nbases=1e9)
errR <- learnErrors(filtRs, verbose=TRUE, multithread=no_of_cores,nbases=1e9)
errPlotF <- plotErrors(errF, nominalQ = TRUE) + labs(title = "Forward")
errPlotR <- plotErrors(errR, nominalQ = TRUE) + labs(title = "Reverse")
g1<-grid.arrange(errPlotF, errPlotR, nrow = 1)
ggsave("errPlot.pdf",plot=g1)

dadaFs <- dada(derepFs, err=errF, pool=TRUE, multithread=no_of_cores, 
               verbose=TRUE)

dadaRs <- dada(derepRs, err=errR, pool=TRUE, multithread=no_of_cores, 
               verbose=TRUE)

cat("dada-class: object describing DADA2 denoising results", sep="\n")
cat(paste(length(dadaFs[[1]]$denoised), 
          "sequence variants were inferred from", 
          length(derepFs[[1]]$uniques), 
          "input unique sequences."), sep="\n")
cat(paste("Key parameters: OMEGA_A =", dadaFs[[1]]$opts$OMEGA_A, 
          "OMEGA_C =", 
          dadaFs[[1]]$opts$OMEGA_C, "BAND_SIZE =", 
          dadaFs[[1]]$opts$BAND_SIZE), sep="\n")

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=no_of_cores, verbose=TRUE)

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x)) # getUniques() gets abundance of unique sequences
# Calculate the number of reads at different steps
track <- cbind(out, 
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "Filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

track.long <- as.data.frame(track,
                            row.names = row.names(track)) %>%
  gather(., key = steps, value = counts, input:nonchim, factor_key = TRUE)

 g2<-ggplot(track.long, aes(x = steps, y = counts, color = steps)) +
    theme_classic() +
    geom_boxplot() +
    geom_jitter(shape = 16, position = position_jitter(0.3)) +
    scale_y_continuous(labels = scales::comma)
ggsave("FilteredReads.pdf",plot=g2)

dir.create("./intermediate/", showWarnings = FALSE)
save(metadata_df, seqtab.nochim, file = "./intermediate/dada2_lean.rda")
write.csv(seqtab.nochim,"seqTable1.csv")
